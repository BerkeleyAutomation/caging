#include "Configuration_Space_Approximator.h"

#include "CageEscapePlanner.h"
#include "Filtration.h"
#include "ShapeFactory.hpp"
#include "Util.h"

#include <algorithm>
#include <glog/logging.h>
//#include <omp.h>
#include <unistd.h>

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
  pose.x = point.point().x() / config_.x_scale;
  pose.y = point.point().y() / config_.y_scale;
  pose.theta = point.point().z() / config_.theta_scale;
}

bool Configuration_Space_Approximator::add_to_triangulation(Pose2D pose, float radius_sq, Triangulation_3& triangulation, CGAL::Color c)
{
  // add all identifications of the point to the triangulation
  Pose2D pose_id;
  pose_id.x = pose.x;
  pose_id.y = pose.y;

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
  for (unsigned int i = 0; i < num_samples; i++) {

    if (i % 10000 == 0) {
      LOG(INFO) << "Sample " << i;
      //      LOG(INFO) << "Tri " << collision_triangulation_.number_of_vertices();
    }

    // sample a pose uniformly at random
    sample_random_pose(pose);
    
    bool collision = check_collisions(pose);
    if (collision) {
      num_collisions++;
    }
  }
  return num_collisions;
}

void Configuration_Space_Approximator::sample_random_pose(Pose2D& pose)
{
  //  uniform_random_pose(pose);
  // gaussian_prm_pose(pose);
  uniform_collision_pose(pose);
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

void Configuration_Space_Approximator::uniform_collision_pose(Pose2D& pose)
{
  float eps = 1e-3;
  MeshCollisionResult coll_result;
  coll_result.collision = false;
  coll_result.distance = 0.0f;

  while (!coll_result.collision && coll_result.distance < eps) {
    uniform_random_pose(pose);
    coll_result = check_object_obstacle_intersection(pose);
  }
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
  // check for collisions
  float eps = 0.25f; // tolerance of libccd according to config (but not actually)
  float nu = 0.95f; // scaling factor to account for libccd tolerance while preserving nonzero pen depth
  MeshCollisionResult coll_result = check_object_obstacle_intersection(pose);
  float distance = coll_result.distance;
  distance = nu * distance;

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
  if (coll_result.collision) {
    if (distance > 0.0f)
      add_to_collision_triangulation(pose, distance*distance);
    return true;
  }
  add_to_free_space_triangulation(pose, distance*distance);
  return false;
}

MeshCollisionResult Configuration_Space_Approximator::check_object_obstacle_intersection(Pose2D pose)
{
  float pen_depth = 0.0f;
  float dist = FLT_MAX;
  float time = 0.0f;
  bool collision = false;

  bool print = false;
  // if (fabs(pose.x - -0.830232) < 1e-5)
  //   print = true;

  if (object_ != NULL && obstacles_.size() > 0) {
    Eigen::Matrix4f pose_mat = CreatePose(pose.x, pose.y, pose.theta);
    object_->SetPose(pose_mat, false, print);

    // computate penetration depth for each potential obstacle
    for (unsigned int j = 0; j < obstacles_.size(); j++) {

      // TODO: change back to lower bound
      MeshCollisionResult c = Mesh::LowerBoundCollision(object_, obstacles_[j], print);
      //      MeshCollisionResult c2 = Mesh::UpperBoundCollision(object_, obstacles_[j]);

      time += c.time;
      if (c.collision) {
        collision = true;
        // LOG(INFO) << "Lower bd " << c.distance;
        // LOG(INFO) << "Upper bd " << c2.distance;

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
  // if (fabs(pose.x - -0.830232) < 1e-5)
  //   while(true);
  // fill in collision result
  MeshCollisionResult result;
  result.collision = collision;
  result.time = time;

  if (!collision) {
    result.distance = dist;
  }
  else {
    result.distance = pen_depth;
  }
  return result;
}

// get the min energy to separate the initial pose from the "escape" space
bool Configuration_Space_Approximator::min_escape_energy(Pose2D pose, EscapeEnergyConfig energy_config, EscapeEnergyResult& result)
{ 
  LOG(INFO) << "Beginning of escape energy";
  // currently must reset the triangulation
  reset();

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
  timer_.reset();
  LOG(INFO) << "Sampled poses with collision ratio " << (float)num_collisions / (float) energy_config.num_samples;
  LOG(INFO) << "Sampling took " << result.sample_time;

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
  LOG(INFO) << "Created alpha shape " << filtered_simplices.size() ;

  // create filtration based on alpha shape and gravity potential function
  GravityPotential gp(start_point.point(), object_->Mass(), config_.x_scale, config_.y_scale, config_.theta_scale);
  //  DistancePotential gp(start_point.point(), config_.x_scale, config_.y_scale, config_.theta_scale);
  AlphaCustomPotential potential_func(alpha_shape, gp, cf_padding_*cf_padding_); 
  LOG(INFO) << "Created potential func";
  Manifold_Sweep_Filtration_3 collision_filtration(alpha_shape, filtered_simplices, potential_func);
  LOG(INFO) << "Created filtration";
  connectivity_checker_ = new Connectivity_Checker(alpha_shape, collision_filtration, gv_);
  timer_.stop();
  LOG(INFO) << "Cage alpha shape constructed in: " << timer_.time() << " seconds.";
  result.triangulation_time = timer_.time();
  timer_.reset();

  // binary search for potential value
  bool caged = false;
  bool connected = false;
  float potential_low = energy_config.energy_min;
  float potential_high = energy_config.energy_max;
  float potential_mid;
  Weighted_point end_point;
  timer_.start();
  while ((potential_high - potential_low) > energy_config.energy_res) {
    // check for connection between the initial point and every other point
    connected = false;
    potential_mid = potential_high / 2.0f + potential_low / 2.0f;
    LOG(INFO) << "Checking potential " << potential_mid << " from bounds " << potential_low << " " << potential_high;
    for (unsigned int i = 0; i < points_at_infinity_.size() && !connected; i++) {
      end_point = points_at_infinity_[i];    
      bool out = connectivity_checker_->connected(start_point.point(), end_point.point(), potential_mid);
      connected |= out;
    }

    // update search range based on outcome
    if (connected) {
      potential_low = potential_mid;
    }
    else {
      potential_high = potential_mid;
    }
  }

  if (potential_high == energy_config.energy_max) {
    LOG(INFO) << "Could not find a disconnection. Upper potential bound may have been too low";
  }
  timer_.stop();
  result.iter_time = timer_.time();
  timer_.reset();

  // same for free space
  bool use_free = false;
  if (use_free) {
    LOG(INFO) << "Using free alpha shape";
    timer_.start();
    Alpha_shape_3 free_alpha_shape(free_space_triangulation_, default_alpha, Alpha_shape_3::GENERAL);
    std::vector<CGAL::Object> free_filtered_simplices;
    std::vector<K::FT> free_alpha_values;
    free_alpha_shape.filtration_with_alpha_values(CGAL::Dispatch_output_iterator<CGAL::cpp11::tuple<CGAL::Object, K::FT>,
                                                                                 CGAL::cpp11::tuple<std::back_insert_iterator<std::vector<CGAL::Object> >,
                                                                                                    std::back_insert_iterator<std::vector<K::FT> > > >(std::back_inserter(free_filtered_simplices),
                                                                                                                                                       std::back_inserter(free_alpha_values))); 
    LOG(INFO) << "Created free alpha shape " << free_filtered_simplices.size();
 
    // create filtration based on alpha shape and gravity potential function
    //  DistancePotential gp(start_point.point(), config_.x_scale, config_.y_scale, config_.theta_scale);
    AlphaCustomPotential free_potential_func(free_alpha_shape, gp, cf_padding_*cf_padding_); 
    LOG(INFO) << "Created potential func";
    Manifold_Sweep_Filtration_3 free_collision_filtration(free_alpha_shape, free_filtered_simplices, free_potential_func);
    LOG(INFO) << "Created filtration";
    connectivity_checker_ = new Connectivity_Checker(free_alpha_shape, free_collision_filtration, gv_);
    timer_.stop();
    LOG(INFO) << "Free space alpha shape constructed in: " << timer_.time() << " seconds.";
    timer_.reset();

    // binary search for potential value
    potential_low = energy_config.energy_min;
    potential_high = energy_config.energy_max;
    timer_.start();
    while ((potential_high - potential_low) > energy_config.energy_res) {
      // check for connection between the initial point and every other point
      connected = false;
      potential_mid = potential_high / 2.0f + potential_low / 2.0f;
      LOG(INFO) << "Checking potential " << potential_mid << " from bounds " << potential_low << " " << potential_high;
      for (unsigned int i = 0; i < points_at_infinity_.size() && !connected; i++) {
        end_point = points_at_infinity_[i];
        bool out = connectivity_checker_->connected(start_point.point(), end_point.point(), potential_mid);
        connected |= out;
      }

      // update search range based on outcome
      if (connected) {
        potential_high = potential_mid;
      }
      else {
        potential_low = potential_mid;
      }
    }
    LOG(INFO) << "Upper potential " << potential_mid;

#ifdef DISPLAY_AS
    std::vector<Cell_handle> incident_cells;
    Locate_type cs_lt, cg_lt;
    int cs_li, cs_lj;
    int cg_li, cg_lj;
    Cell_handle start_cell, end_cell;
    start_cell = free_alpha_shape.locate(start_point.point(), cs_lt, cs_li, cs_lj); //find start cell //TODO: possible bug here to due with the locate query returning something 
    end_cell = free_alpha_shape.locate(end_point.point(), cg_lt, cg_li, cg_lj);


    if (cs_lt == Triangulation_3::VERTEX) {
      Vertex_handle vh = start_cell->vertex(cs_li);
      free_alpha_shape.incident_cells(vh, std::back_inserter(incident_cells));
      start_cell = incident_cells[0];
      if (free_alpha_shape.classify(start_cell) == Alpha_shape_3::EXTERIOR) {
        LOG(INFO) << "Start cell exterior";      
      }
    }

    free_alpha_shape.set_alpha(0.0f);
    gv_.clear();
    gv_ << CGAL::DEEPBLUE;
    gv_ << free_alpha_shape;
    gv_ << CGAL::PURPLE;

    float y = potential_mid / (-9.81f * object_->Mass());
    Bare_point p0(y, 1000, 1000); 
    Bare_point p1(y, -1000, 1000); 
    Bare_point p2(y, 500, -1000); 
    Tri_3 t(p0, p1, p2);
    //  gv_ << t;

    gv_ << CGAL::GREEN;
    gv_ << Convert_Cell_To_Tetrahedron(start_cell);

    incident_cells.clear();
    if (cg_lt == Triangulation_3::VERTEX) {
      Vertex_handle vh = end_cell->vertex(cg_li);
      free_alpha_shape.incident_cells(vh, std::back_inserter(incident_cells));
      end_cell = incident_cells[0];
      gv_ << CGAL::RED;
      gv_ << Convert_Cell_To_Tetrahedron(end_cell);      
      //    LOG(INFO) << "End v";
    }
#endif
    while(true);
  }

  // fill in result
  result.energy = potential_mid; //std::min<float>(potential_mid, 0.0f); // assumed that the initial configuration is zero energy
  result.normalized_energy = potential_mid / (GravityPotential::GRAVITY_ACCEL * object_->Mass());
  result.path_exists = (potential_low > energy_config.energy_min); 
  result.coll_ratio = (float)num_collisions / (float)energy_config.num_samples;
  //  connectivity_checker_->connected_component_volume(start_point.point(), result.volume);

  // optional display alpha shape
#ifdef DISPLAY_AS
  std::vector<Cell_handle> incident_cells;
  Locate_type cs_lt, cg_lt;
  int cs_li, cs_lj;
  int cg_li, cg_lj;
  Cell_handle start_cell, end_cell;
  start_cell = alpha_shape.locate(start_point.point(), cs_lt, cs_li, cs_lj); //find start cell //TODO: possible bug here to due with the locate query returning something 
  end_cell = alpha_shape.locate(end_point.point(), cg_lt, cg_li, cg_lj);

  if (cs_lt == Triangulation_3::VERTEX) {
    Vertex_handle vh = start_cell->vertex(cs_li);
    alpha_shape.incident_cells(vh, std::back_inserter(incident_cells));
    start_cell = incident_cells[0];
    LOG(INFO) << "Start v";
  }

  alpha_shape.set_alpha(0.0f);
  gv_.clear();
  gv_ << CGAL::DEEPBLUE;
  gv_ << alpha_shape;
  gv_ << CGAL::PURPLE;

  float y = potential_mid / (-9.81f * object_->Mass());
  Bare_point p0(y, 1000, 1000); 
  Bare_point p1(y, -1000, 1000); 
  Bare_point p2(y, 500, -1000); 
  Tri_3 t(p0, p1, p2);
  //  gv_ << t;

  gv_ << CGAL::GREEN;
  gv_ << Convert_Cell_To_Tetrahedron(start_cell);

  incident_cells.clear();
  if (cg_lt == Triangulation_3::VERTEX) {
    Vertex_handle vh = end_cell->vertex(cg_li);
    alpha_shape.incident_cells(vh, std::back_inserter(incident_cells));
    end_cell = incident_cells[0];
    gv_ << CGAL::RED;
    gv_ << Convert_Cell_To_Tetrahedron(end_cell);      
    //    LOG(INFO) << "End v";
  }
#endif
}
