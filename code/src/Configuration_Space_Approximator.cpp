#include "Configuration_Space_Approximator.h"

#include "CageEscapePlanner.h"
#include "Filtration.h"
#include "ShapeFactory.hpp"
#include "Util.h"

#include <algorithm>
#include <glog/logging.h>
//#include <omp.h>
#include <unistd.h>

//#define DISPLAY_AS
#define DISPLAY_RATE 10000

Configuration_Space_Approximator::Configuration_Space_Approximator(CfApproxConfig config, CGAL::Geomview_stream& gv)
  : connectivity_checker_(NULL), config_(config), gprm_norm_(0.0f, 10.0f), gprm_distribution_(rng_, gprm_norm_), gv_(gv)
{
  compute_cf_padding();
  add_points_at_infinity();
}

Configuration_Space_Approximator::~Configuration_Space_Approximator()
{
  collision_triangulation_.clear();
  free_space_triangulation_.clear();
  if (connectivity_checker_)
    delete connectivity_checker_;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Configuration_Space_Approximator::Clear_Stored_Data
 *  Description:  This class clears all of the data associated with the Planner.  After
 *                this function is called, the planner will be reset to initialized state.
 * =====================================================================================
 */
void Configuration_Space_Approximator::reset()
{
  collision_triangulation_.clear();
  free_space_triangulation_.clear();
  points_at_infinity_.clear();
  if (connectivity_checker_)
    delete connectivity_checker_;

  compute_cf_padding();
  add_points_at_infinity();
}

void Configuration_Space_Approximator::compute_cf_padding()
{
  float x_dim = config_.x_scale * (config_.x_max - config_.x_min);
  float y_dim = config_.y_scale * (config_.x_max - config_.x_min);
  float theta_dim = config_.theta_scale * (config_.theta_max - config_.theta_min);
  cf_padding_ = 2 * std::max<float>(theta_dim, std::max<float>(x_dim, y_dim));
}

void Configuration_Space_Approximator::pose_to_point(Pose2D pose, float radius_sq, Weighted_point& point)
{
  point = Weighted_point(Bare_point(config_.x_scale * pose.x, 
                                    config_.y_scale * pose.y,
                                    config_.theta_scale * pose.theta),
                         radius_sq);
}

void Configuration_Space_Approximator::point_to_pose(Weighted_point point, Pose2D& pose)
{
  pose.x = (float)(point.point().x() / config_.x_scale);
  pose.y = (float)(point.point().y() / config_.y_scale);
  pose.theta = (float)(point.point().z() / config_.theta_scale);
}

void Configuration_Space_Approximator::point_3_to_pose(Point_3 point, Pose2D& pose)
{
  pose.x = CGAL::to_double(point.x()) / config_.x_scale;
  pose.y = CGAL::to_double(point.y()) / config_.y_scale;
  pose.theta = CGAL::to_double(point.z()) / config_.theta_scale;
}

bool Configuration_Space_Approximator::add_to_triangulation(Pose2D pose, float radius_sq, Triangulation_3& triangulation, CGAL::Color c)
{
  // add all identifications of the point to the triangulation
  Pose2D pose_id;
  pose_id.x = pose.x;
  pose_id.y = pose.y;

  if (pose.x > max_x_val_) {
    max_x_val_ = pose.x;
  }

  for (int i = -config_.num_rots-1; i <= config_.num_rots; i++) {
    // get an identification of theta
    Weighted_point new_point;
    pose_id.theta = pose.theta + 2 * M_PI * i;
    pose_to_point(pose_id, radius_sq, new_point);

    //#pragma omp critical
    {
      triangulation.insert(new_point);
    }
    // visualization (still with macros :/)
  }    
#ifdef DISPLAY_SPHERES
  gv_ << c;
  gv_ << K::Sphere_3(Bare_point(pose.x, pose.y, pose.theta), radius_sq);
#endif
}

bool Configuration_Space_Approximator::add_to_collision_triangulation(Pose2D pose, float radius_sq)
{
  add_to_triangulation(pose, radius_sq, collision_triangulation_, CGAL::PURPLE);
}

bool Configuration_Space_Approximator::add_to_free_space_triangulation(Pose2D pose, float radius_sq)
{
  add_to_triangulation(pose, radius_sq, free_space_triangulation_, CGAL::ORANGE);
}

// Stores all of the boundary points for all equivalent values of theta
void Configuration_Space_Approximator::add_points_at_infinity()
{
  points_at_infinity_.clear();
  Weighted_point new_point;

  // check collision for free space points (since we need them in the triangulation
  float theta_rot_min = 2 * M_PI * (-config_.num_rots-1);
  float theta_rot_max = 2 * M_PI * (config_.num_rots+1);

  float x_infty_min = config_.x_scale * config_.x_min - cf_padding_;
  float x_infty_max = config_.x_scale * config_.x_max + cf_padding_;
  float y_infty_min = config_.y_scale * config_.y_min - cf_padding_;
  float y_infty_max = config_.y_scale * config_.y_max + cf_padding_;
  float theta_infty_min = config_.theta_scale * theta_rot_min;// - cf_padding_;
  float theta_infty_max = config_.theta_scale * theta_rot_max;// + cf_padding_;

  Pose2D pose_infty;
  for (float x = x_infty_min; x <= x_infty_max; x += x_infty_max - x_infty_min) { 
    for (float y = y_infty_min; y <= y_infty_max; y += y_infty_max - y_infty_min) {
      for (float theta = theta_infty_min; theta <= theta_infty_max; theta += theta_infty_max - theta_infty_min) {

        // add free space point at "infinty"
        new_point = Weighted_point(Bare_point(x, y, theta),
                                   3*cf_padding_*cf_padding_);
        points_at_infinity_.push_back(new_point);
      }
    }
  }
}

// samples a number of random poses and adds to the triangulation
unsigned int Configuration_Space_Approximator::sample_poses(Pose2D pose_orig, unsigned int num_samples)
{
  unsigned int num_collisions = 0;
  Pose2D pose;

  //  #pragma omp parallel for
  LOG(INFO) << num_samples;
  for (unsigned int i = 0; i < num_samples; i++) {

    if (i % DISPLAY_RATE == 0) {
      LOG(INFO) << "Sample " << i;
      //      LOG(INFO) << "Tri " << collision_triangulation_.number_of_vertices();
    }

    // sample a pose in collision using REJECTION sampling
    MeshCollisionResult coll_result = sample_random_pose(pose); 
    float nu = 0.95f; // scaling factor to account for libccd tolerance while preserving nonzero pen depth
    float distance = coll_result.distance;
    distance = nu * distance;

    if (coll_result.collision) {
      if (distance > 0.0f)
        add_to_collision_triangulation(pose, distance*distance);
        num_collisions++;
    } else {
      //add_to_free_space_triangulation(pose, distance*distance);
    }
  }
  return num_collisions;
}

MeshCollisionResult Configuration_Space_Approximator::sample_random_pose(Pose2D& pose)
{
  //  uniform_random_pose(pose);
  // gaussian_prm_pose(pose);
  return uniform_collision_pose(pose);
  // float eps = 0.0f;
  // float toss = (float)rand()/(float)RAND_MAX;
  // if (toss < eps) {
  //   uniform_collision_pose(pose);
  // }
  // else {
  //   gaussian_prm_pose(pose);
  // }
}

void Configuration_Space_Approximator::uniform_random_pose(Pose2D& pose)
{
  pose.x = ((float)rand()/(float)RAND_MAX) * (config_.x_max - config_.x_min) + config_.x_min;
  pose.y = ((float)rand()/(float)RAND_MAX) * (config_.y_max - config_.y_min) + config_.y_min;
  pose.theta = ((float)rand()/(float)RAND_MAX) * (config_.theta_max - config_.theta_min) + config_.theta_min;; //((2 * M_PI); // don't use bounds, always uniform at random
}

MeshCollisionResult Configuration_Space_Approximator::uniform_collision_pose(Pose2D& pose)
{ 
  float eps = 1e-3;
  MeshCollisionResult coll_result;
  coll_result.collision = false;
  coll_result.distance = 0.0f;

  while (!coll_result.collision && coll_result.distance < eps) {
    uniform_random_pose(pose);
    coll_result = check_object_obstacle_intersection(pose);
  }
  //Return that this pose is in collision...
  return coll_result;
}

void Configuration_Space_Approximator::gaussian_prm_pose(Pose2D& pose)
{
  Pose2D other_pose;
  MeshCollisionResult coll_result;
  MeshCollisionResult other_coll_result;
  float dist;
  coll_result.collision = false;
  other_coll_result.collision = false;

  // gaussian sampling until triggered
  while ((coll_result.collision && other_coll_result.collision) || 
         (!coll_result.collision && !other_coll_result.collision)) {
    uniform_random_pose(other_pose);

    dist = abs(gprm_distribution_());
    float theta_pos = ((float)rand()/(float)RAND_MAX) * (2 * M_PI);
    //    float theta_rot = ((float)rand()/(float)RAND_MAX) * (2 * M_PI);
    pose.x = other_pose.x + dist * cos(theta_pos);
    pose.y = other_pose.y + dist * sin(theta_pos);
    pose.theta = other_pose.theta;//theta_rot;

    coll_result = check_object_obstacle_intersection(pose);
    other_coll_result = check_object_obstacle_intersection(other_pose);
  }           

  // always keep the pose in collision
  if (!coll_result.collision && other_coll_result.collision) {
    pose = other_pose;
  }
}

void Configuration_Space_Approximator::level_set_random_pose(Pose2D& pose)
{
  // float level_set = 2.0f;
  // std::vector<Bare_point> pose_particles;
  // std::vector<float> particle_weights;

  // boost::mt19337 gen;
  // boost::random::discrete_distribution<float> d(particle_weights);
  // int sample_ind = d(gen);
  
  // sample normal random variable around point
  

}

bool Configuration_Space_Approximator::check_collisions(Pose2D pose)
{
  MeshCollisionResult coll_result = check_object_obstacle_intersection(pose);
  //distance = std::max<float>(distance - 0.01f, 0.0f);

  // cap the distance by the known distance to the origin to correct for CCD tolerance issues
  // float theta = pose.theta;
  // if (theta > M_PI) {
  //   theta = -2 * M_PI + theta;
  // }
  // float dist_to_origin = sqrt((pose.x * pose.x) + (pose.y * pose.y) + config_.theta_scale * config_.theta_scale * (theta * theta));
  // if (distance > dist_to_origin) {
    // LOG(INFO) << "Pose " << pose.x << " " << pose.y << " " << pose.theta; 
    // LOG(INFO) << "Dist " << distance;
    // LOG(INFO) << "Dist to origin " << dist_to_origin;
    // distance = dist_to_origin - 1e-3f;
  //  }
  // if (fabs(pose.x - -0.830232) < 1e-5)
  //   while(true);

  // add to appropriate triangulation
  return coll_result.collision;
}


//Check collisions without adding to a triangulation....
bool Configuration_Space_Approximator::check_collisions_safe(Pose2D pose)
{
  // check for collisions
  MeshCollisionResult coll_result = check_object_obstacle_intersection(pose, false);
  if (coll_result.collision) {
    return true;
  } else {
    return false;
  }
}

MeshCollisionResult Configuration_Space_Approximator::check_object_obstacle_intersection(Pose2D pose, bool print, bool upper)
{
  //Timer timeout;
  //timeout.start();

  float pen_depth = 0.0f;
  float dist = FLT_MAX;
  float time = 0.0f;
  bool collision = false;

  if (object_ != NULL && obstacles_.size() > 0) {
    Eigen::Matrix4f pose_mat = CreatePose(pose.x, pose.y, pose.theta);
    object_->SetPose(pose_mat, false, false);

    if (print) {
      std::cout << "Set pose " << std::endl << pose_mat << std::endl;
    }

    // computate penetration depth for each potential obstacle
    for (unsigned int j = 0; j < obstacles_.size(); j++) {

      // TODO: change back to lower bound
      MeshCollisionResult c;
      if (!upper)
        c = Mesh::LowerBoundCollision(object_, obstacles_[j], false);
      else
        c = Mesh::UpperBoundCollision(object_, obstacles_[j]);

      time += c.time;
      if (c.collision) {
        collision = true;

        if (c.distance > pen_depth) {
          pen_depth = c.distance;
        }
      }
      else {
        if (c.distance < dist) {
          dist = c.distance;
        }
      }
    }
  }
  // fill in collision result
  MeshCollisionResult result;
  result.collision = collision;
  result.time = time;

  //timeout.stop();
  //if (timeout.time() > 2.0f){
    //std::cout << "COLL CHECK TIME :" << timeout.time() << std::endl; 
    //std::cout << "IN COLLISION: " << collision << std::endl;
    //std::cout << "PEN DEPTH: " << pen_depth << std::endl; 
    //visualize_pose(pose);
    //wait(5);
  //}
  
  if (!collision) {
    result.distance = dist;
  }
  else {
    result.distance = pen_depth;
  }
  return result;

}

void Configuration_Space_Approximator::visualize_pose(Pose2D pose) 
{
  gv_.clear();
  gv_ << CGAL::WHITE;
  Eigen::Matrix4f object_pose = CreatePose(pose.x, pose.y, pose.theta);
  Eigen::Matrix4f temp_pose = object_->Pose();
  object_->SetPose(object_pose, true); 
  object_->RenderToGeomview(gv_);
  for (unsigned int j = 0; j < obstacles_.size(); j++) { 
    obstacles_[j]->RenderToGeomview(gv_);
  }
  wait(1);
  object_->SetPose(temp_pose, true); 
}



float Configuration_Space_Approximator::void_volume(std::vector<ASE::Simplex> void_simplices) {

  float volume = 0.0f;
  for (int i = 0; i < void_simplices.size(); i++) {
    ASE::Simplex cur_simplex = void_simplices[i];
    Bare_point cur_centroid = cur_simplex.centroid();
    Pose2D cur_pose;
    point_to_pose(cur_centroid, cur_pose);         
    if (check_collisions_safe(cur_pose)) {
      continue;
    }
    // get vertices and compute volume
    Tetrahedron_3 tetra = Convert_Simplex_To_Tetrahedron(cur_simplex);    
    volume = volume + tetra.volume();
  }

  return volume;

}

//Check push reacability via a simple line search
bool Configuration_Space_Approximator::check_push_reach(Pose2D pose, Eigen::Vector2d force_direction, float disc, float backup)
{
  Eigen::VectorXd cur_vector(2);
  cur_vector(0) = pose.x;
  cur_vector(1) = pose.y;
  Eigen::Vector2d inc_vector(force_direction(0), force_direction(1));
  inc_vector = -(backup/disc)*inc_vector;
  float theta = pose.theta;

  for (int i = 0; i < disc; i++) {
    Pose2D test_pose;
    test_pose.x = (float) cur_vector[0];
    test_pose.y = (float) cur_vector[1];
    test_pose.theta = theta;
    if (check_collisions_safe(test_pose)) {
      return false;
    }
    cur_vector = cur_vector + inc_vector; 
  }

  return true;
}

float Configuration_Space_Approximator::void_opening_area(std::vector<ASE::Simplex> void_simplices, Bare_point birth_vertex, AlphaCustomPotential potential_func, Eigen::Vector2d force_vector,
                                        Potential birth_energy) 
{
  float total_area = 0.0f;

  for (int i = 0; i < void_simplices.size(); i++) {
    ASE::Simplex cur_simplex = void_simplices[i];

    float projected_area = simplex_projected_area(cur_simplex, birth_vertex, potential_func, force_vector, birth_energy);
    total_area = total_area + projected_area;   
  }

  return total_area;
}

//What area does a simplex on the boundary of a void's CC cast onto the opening of the void
float Configuration_Space_Approximator::simplex_projected_area(ASE::Simplex cur_simplex, Bare_point birth_vertex, AlphaCustomPotential potential_func, Eigen::Vector2d force_vector, Potential birth_energy)
{
  //Get point vecotrs for corners of a triangle completely abve void opening
  std::vector< Eigen::Vector3d > triangle_points;
  for (int i = 0; i < 4; i ++) {
    ASE::Vertex cur_vertex = cur_simplex.vertex(i);
    if (potential_func.potential(cur_vertex) <= birth_energy && cur_simplex.potential() != -FLT_MAX) {
      
      Eigen::Vector3d point_vector(cur_vertex.point().x(), cur_vertex.point().y(), cur_vertex.point().z());
        
      triangle_points.push_back(point_vector);
    }
  }

  if (triangle_points.size() < 3) {
    return 0.0f;
  }

  Eigen::Vector3d known_point(birth_vertex.x(), birth_vertex.y(), birth_vertex.z());

  //Force direction is the normal vector of the plane represetning the void opening
  Eigen::Vector3d normal_vector(force_vector(0), force_vector(1), 0.0f);

  //Some math to actually get projected points onto plabe
  double d = -known_point.dot(normal_vector);

  for (int i = 0; i < 4; i++) {
    Eigen::Vector3d point_vector = triangle_points[i];
    double backup = normal_vector.dot(point_vector) + d;
    point_vector = point_vector - backup*normal_vector;
    triangle_points[i] = point_vector;
  }

  return triangle_area(triangle_points);
}



// get the min energy to separate the initial pose from the "escape" space
std::vector< std::vector<synthesis_info> > Configuration_Space_Approximator::synthesize_grasps(int num_searches, Pose2D pose, 
			    EscapeEnergyConfig energy_config, EscapeEnergyResult& result, float angle_sweep, float angle_disc, bool check_reachability, float max_push_force)
{ 
  LOG(INFO) << "Beginning of escape energy";

  // currently must reset the triangulation
  reset();
  max_x_val_ = 0.0f;

  // params for checking reachability
  float disc = 200;
  float backup = 100;
  float min_energy_thresh = 0.5f;

  // check collisions for initial pose
  check_collisions(pose);
  Weighted_point start_point;
  pose_to_point(pose, 0, start_point);

  // add points at infinity to triangulation (in order for the convex hull to contain the workspace)
  for (unsigned int i = 0; i < points_at_infinity_.size(); i++) {
    collision_triangulation_.insert(Weighted_point(points_at_infinity_[i].point(), -points_at_infinity_[i].weight()));
    free_space_triangulation_.insert(points_at_infinity_[i]);
  }

  LOG(INFO) << "Inserted points at infinity";

  // sample a bunch of random points inside the rectangloid to begin the search
  timer_.start();
  unsigned int num_collisions = sample_poses(pose, energy_config.num_samples);
  timer_.stop();
  result.sample_time = timer_.time();
  float sample_time = timer_.time();
  timer_.reset();
  LOG(INFO) << "Sampled poses with collision ratio " << (float)num_collisions / (float) energy_config.num_samples;
  LOG(INFO) << "Sampling took " << result.sample_time;
  LOG(INFO) << "Max X Val was " << max_x_val_;

  // create alpha shape for collision space and use it to generate a filtration over gravity and alpha
  timer_.start();
  float default_alpha = 0.0f;
  Alpha_shape_3 alpha_shape(collision_triangulation_, default_alpha, Alpha_shape_3::GENERAL);
  std::vector<CGAL::Object> filtered_simplices;
  std::vector<K::FT> alpha_values;
  alpha_shape.filtration_with_alpha_values(CGAL::Dispatch_output_iterator<CGAL::cpp11::tuple<CGAL::Object, K::FT>,
                                                                          CGAL::cpp11::tuple<std::back_insert_iterator<std::vector<CGAL::Object> >,
                                                                                             std::back_insert_iterator<std::vector<K::FT> > > >(std::back_inserter(filtered_simplices),
                                                                                                                                                std::back_inserter(alpha_values))); 
  LOG(INFO) << "Created alpha shape " << filtered_simplices.size();
  timer_.stop();
  LOG(INFO) << "Cage alpha shape constructed in: " << timer_.time() << " seconds.";
  result.triangulation_time = timer_.time();
  float triangulation_time = timer_.time();
  timer_.reset();
  
  std::vector< std::vector<synthesis_info> > all_poses; 

  // sweep through push directions
  angle_sweep = std::min<float>(angle_sweep, M_PI); // cap angle for sweep
  for (float angle_offset = -angle_sweep; angle_offset <= angle_sweep; angle_offset += angle_disc) { 
    LOG(INFO) << "Checking angle offset " << angle_offset;
    
    float pi = M_PI;
    float pi_neg = -M_PI;
   
    //If there's overlap, don't do the same angle twice!
    if (angle_offset == pi_neg && angle_sweep == pi) {
      continue;
    }

    Timer persistence_timer;
    persistence_timer.start();
    
    // create potential function
    LinearPotentialFunction push_force(start_point.point(), max_push_force, config_.x_scale,
                                       config_.y_scale, config_.theta_scale, angle_offset);
    AlphaCustomPotential potential_func(alpha_shape, push_force, cf_padding_*cf_padding_); 
    LOG(INFO) << "Created potential func";

    // create filtration based on alpha shape and gravity potential function
    Manifold_Sweep_Filtration_3 collision_filtration(alpha_shape, filtered_simplices, potential_func);
    LOG(INFO) << "Created filtration";
    persistence_timer.stop();
    float filtration_time = persistence_timer.time();

    // create connectivity checker for parsing void
    LOG(INFO) << "Create connectivity checker";
    bool cache_connected_components = false;
    float p = 0.0f;
    connectivity_checker_ = new Connectivity_Checker(alpha_shape, collision_filtration,
                                                     gv_, p, cache_connected_components);

    LOG(INFO) << "Running persistence";
    persistence_timer.start();
    ASE::Vertex birth_point;
    std::vector< Persistence_info > delta_pairs = collision_filtration.run_persistence(points_at_infinity_, -FLT_MAX,
                                                                                       config_.theta_scale,
                                                                                       birth_point,
                                                                                       min_energy_thresh);
    persistence_timer.stop();
    float persistence_time = persistence_timer.time();

    //Just useful info for debugging purposes. Not needed for core algo..
    Pose2D best_pose;
    std::vector<ASE::Simplex> best_void_simplices; 
    Bare_point best_centroid;
    Persistence_info optimal_info;
    Potential best_delta = 0.0f;

    //Do the region growing algorithm for each void.
    std::vector<synthesis_info> extracted_poses; 
    Potential seen_void_potential = 0.0;
    int num_searched = 0;

    std::ofstream of("chosen_persistance_pairs.csv", std::ios::app);
    for (int i = delta_pairs.size() - 1; i >= 0 && num_searched < num_searches; i--) {
      Potential void_potential = delta_pairs[i].get_potential();
      //LOG(INFO) << "Checking void " << i << " with potential " << void_potential;
      
      //The voids are sorted by delta, so we know when we're done...
      if (void_potential != seen_void_potential) {
        Timer void_timer;
        void_timer.start();

        // walk through the void, looking for highest energy void in collision 
        ASE::Simplex death_simplex = delta_pairs[i].get_death_simplex();    
        ASE::Simplex birth_simplex = delta_pairs[i].get_birth_simplex();

        // get birth pose for debugging
        Bare_point birth_vertex;
        Potential birth_energy = -FLT_MAX;
        ASE::Vertex v;
        for (int a = 0; a < 3; a++) {
          v = birth_simplex.vertex(a);
          Potential v_energy = push_force.potential(v);
          if (v_energy > birth_energy) {
            birth_vertex = v.point();
            birth_energy = v_energy;
          }
          //std::cout << "VISUALIZE POSE!" << std::endl;
	  //Pose2D cur_pose;
	  //point_to_pose(v, cur_pose);
	  //std::cout << cur_pose.x << " " << cur_pose.y << " " << cur_pose.theta << std::endl;
	  //visualize_pose(cur_pose);
	  //wait(10);
	}
        Pose2D birth_pose;
        point_to_pose(birth_vertex, birth_pose);         

        // get simplices in void and sort
        //LOG(INFO) << "Getting connected component";
        std::vector<ASE::Simplex> void_simplices =                      \
          connectivity_checker_->connected_component(death_simplex.centroid(), 
                                                     birth_simplex.potential());
        std::sort(void_simplices.begin(), void_simplices.end());
        
        // check collisions for each simplex
        bool void_searched = true;
        bool found_coll_free_pose = false;
        for (int j = void_simplices.size() - 1; j >= 0 && !found_coll_free_pose; j--) {
          //LOG(INFO) << "Checking simplex " << j;

          // get next pose in void
          ASE::Simplex cur_simplex = void_simplices[j];
          Bare_point cur_centroid = cur_simplex.centroid();
          Pose2D cur_pose;
          point_to_pose(cur_centroid, cur_pose);         
          bool cur_collision = check_collisions_safe(cur_pose);   
          ASE::Vertex cur_vertex = Weighted_point(cur_centroid, 0.0f);
          Potential cur_potential = potential_func.potential(cur_vertex);
          Potential cur_delta = abs(cur_potential - birth_simplex.potential());

          // check validity
          // 1. centroid within void (external in current filtration)
          // 2. not in collision
          // 3. reachable by linear push, IF we care about that (see check_reachability)
          // 4. within 0 to 2pi angle bounds
          if (cur_potential > birth_simplex.potential() && cur_delta > min_energy_thresh &&
              !cur_collision &&
              (!check_reachability ||  check_push_reach(cur_pose, push_force.get_force_vector(), disc, backup)) &&
              cur_pose.theta >= 0.0f && cur_pose.theta <= 2*M_PI) {

            // set best pose (for debugging only)
            if (cur_delta > best_delta) {
              best_delta = cur_delta;
              best_pose = cur_pose;
              best_void_simplices = void_simplices;
              best_centroid = cur_centroid;
              optimal_info = delta_pairs[i];
            }          

            void_timer.stop();

            LOG(INFO) << "Birth pose " << birth_pose.x << " " << birth_pose.y << " " << birth_pose.theta;
            LOG(INFO) << "EBC pose " << cur_pose.x << " " << cur_pose.y << " " << cur_pose.theta;
            LOG(INFO) << "EBC delta " << cur_delta;
            LOG(INFO) << "EBC energy " << cur_simplex.potential();

            // MeshCollisionResult lower_res = check_object_obstacle_intersection(pose, false);
            // MeshCollisionResult upper_res = check_object_obstacle_intersection(pose, true);
            // LOG(INFO) << "Lower collision: " << lower_res.collision;
            // LOG(INFO) << "Upper collision: " << upper_res.collision;
            // visualize_pose(cur_pose);
            // sleep(10);

            ASE::Vertex cur_vertex;
            Potential min_energy = FLT_MAX;
            for (int a = 0; a < 3; a++) {
              v = cur_simplex.vertex(a);
              Potential v_energy = push_force.potential(v);
              if (v_energy < min_energy) {
                cur_vertex = v.point();
                min_energy = v_energy;
              }
            }
            LOG(INFO) << "Cur energy " << min_energy;
            LOG(INFO) << "Cur vertex " << cur_vertex;

            
	    
            //DEBUG
	    //gv_.clear();
            //gv_ << CGAL::YELLOW;
	    //gv_ << Kernel::Sphere_3(Point_3(cur_centroid.x(), cur_centroid.y(), cur_centroid.z()), 0.05f);
	    //gv_ << CGAL::RED;
	    //gv_ << Convert_Simplex_To_Triangle(birth_simplex);
	    //alpha_shape.set_alpha(0.0f);
	    //gv_ << CGAL::GREEN;  
	    //for (int x = 0; x < void_simplices.size(); x++){
            //  for (int a = 0; a < 3; a++) {
            //    ASE::Vertex v  = void_simplices[x].vertex(a);
	    //    Pose2D cur_pose;
	    //    point_to_pose(v, cur_pose);
	    //    std::cout << "CUR POSE :" << " " << cur_pose.x << " " << cur_pose.y << " " << cur_pose.theta << std::endl;
	    //  }
	    //  gv_ << Convert_Simplex_To_Tetrahedron(void_simplices[x]);
	    //}
	    //wait(60);
            //gv_ << CGAL::DEEPBLUE;
	    ////gv_ << alpha_shape; 
            //gv_ << Convert_Simplex_To_Triangle(birth_simplex); 
	    //wait(1000);
	    
	    
	    // extract void info
            float void_time = void_timer.time();
            visualize_pose(cur_pose); 
            //DEBUG
	    float v_volume = void_volume(void_simplices); 
            float v_area = void_opening_area(void_simplices, birth_vertex, potential_func, push_force.get_force_vector(), birth_simplex.potential()); 
	    float norm_delta = cur_delta / (max_push_force);
            synthesis_info extracted(cur_pose,
                                     birth_pose,
                                     cur_delta,
                                     norm_delta,
                                     v_volume,
				     v_area,
                                     sample_time,
                                     triangulation_time,
                                     filtration_time,
                                     persistence_time,
                                     void_time,
                                     angle_offset); 
            extracted_poses.push_back(extracted); 
            
            // set flag to exit loop
            found_coll_free_pose = true;
            seen_void_potential = void_potential;
            num_searched++;

	    //Second damp of pair info ONLY for the chosen poses///
	    of << delta_pairs[i].birth_index_ << ", " << delta_pairs[i].death_index_ << "\n";

          }
        }
      }
    } 

    ASE::Simplex death_simplex = optimal_info.get_death_simplex();
    ASE::Simplex birth_simplex = optimal_info.get_birth_simplex();
    Potential birth_energy = birth_simplex.potential();
    
    if (extracted_poses.size() > 0) {
      LOG(INFO) << "FOUND ENERGY-BOUNDED CAGES";
      LOG(INFO) << "Here is the best pose: " << best_pose.x << " " << best_pose.y << " " << best_pose.theta;
      LOG(INFO) << "Here is the largest delta: " << best_delta;
    }
    else {
      LOG(INFO) << "DID NOT FIND ANY ENERGY-BOUNDED CAGES";
    }

    // sort in order of increasing potential
    std::sort(extracted_poses.begin(), extracted_poses.end());
    all_poses.push_back(extracted_poses); 
    LOG(INFO) << "SORTED";
   
    //for (int i = 0; i < extracted_poses.size(); i++) {
    //  Pose2D cur_pose;
    //  cur_pose = extracted_poses[i].get_pose();
    //  Potential cur_delta = extracted_poses[i].get_delta();
    //  std::cout << "Here is the cur pose: " << cur_pose.x << " " << cur_pose.y << " " << cur_pose.theta << std::endl;
    //  std::cout << "Here is the cur delta: " << cur_delta << std::endl; 
    //}


    //Cell_handle death_cell;
    //if (!CGAL::assign(death_cell, death_simplex.object())) {
    //  LOG(INFO) << "Could not extract cell handle from death simplex";
    //  exit(0);
    //}

    ////Raw_triangulation_3 cc_tri =                                        \
    //    //  connectivity_checker_->connected_component(death_cell, birth_energy);

    //// optional display alpha shape
    std::vector<Cell_handle> incident_cells;
    Locate_type cs_lt, cg_lt, cb_lt;
    int cs_li, cs_lj;
    //int cg_li, cg_lj;
    //int cb_li, cb_lj;
    Cell_handle start_cell, birth_cell;
    start_cell = alpha_shape.locate(start_point.point(), cs_lt, cs_li, cs_lj); //find start cell
    ////birth_cell = alpha_shape.locate(birth_point.point(), cb_lt, cb_li, cb_lj); //find start cell //TODO: possible bug here to due with the locate query returning something 
    ////death_cell = alpha_shape.locate(death_simplex.point(), cg_lt, cg_li, cg_lj); //find start cell //TODO: possible bug here to due with the locate query returning something 

    if (cs_lt == Triangulation_3::VERTEX) {
      Vertex_handle vh = start_cell->vertex(cs_li);
      alpha_shape.incident_cells(vh, std::back_inserter(incident_cells));
      start_cell = incident_cells[0];
    }

    of.close();

#ifdef DISPLAY_AS
    gv_.clear();
    alpha_shape.set_alpha(0.0f);
    gv_ << CGAL::DEEPBLUE;
    gv_ << alpha_shape; 
    wait(1000);
    gv_ << CGAL::GREEN;
    gv_ << Convert_Cell_To_Tetrahedron(start_cell);

    gv_ << CGAL::YELLOW;
    gv_ << Convert_Simplex_To_Tetrahedron(death_simplex);

    gv_ << CGAL::RED;
    gv_ << Convert_Simplex_To_Triangle(birth_simplex);

    // float y = birth_energy / (-9.81f * object_->Mass());
    // Bare_point p0(y, 1000, 1000); 
    // Bare_point p1(y, -1000, 1000); 
    // Bare_point p2(y, 1000, -1000); 
    // Tri_3 t(p0, p1, p2);
    // gv_ << CGAL::PURPLE;
    // gv_ << t;

    // Bare_point p3(y, -1000, 1000); 
    // Bare_point p4(y, 1000, -1000); 
    // Bare_point p5(y, -1000, -1000); 
    // Tri_3 t2(p3, p4, p5);
    // gv_ << CGAL::PURPLE;
    // gv_ << t2;
    while(true);
#endif
    
  }
  return all_poses;
}


