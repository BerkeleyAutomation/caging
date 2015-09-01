#include "Triangulation_Planner.h"

#include "CageEscapePlanner.h"
#include "Filtration.h"
#include "ShapeFactory.hpp"
#include "Util.h"

#include <unistd.h>

#define GOAL_DIST_THRESH 2.5f
#define GRAVITY_ACCEL 9.81f
#define MAX_AREA 1000.0f

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Triangulation_Planner::Triangulation_Planner
 *  Description:  This function implements the constructor for a Triangulation Planner
 *                class.  
 * =====================================================================================
 */
    
Triangulation_Planner::Triangulation_Planner(float xscale, float yscale, float theta_scale, 
                                             float xmin, float xmax, float ymin, float ymax, float theta_min, float theta_max, int num_rots,
                                             CGAL::Geomview_stream& gv,
                                             float alpha_max, float alpha_res, float level_set)
  : x_scale_(xscale), y_scale_(yscale), theta_scale_(theta_scale),
    x_min_(xmin), x_max_(xmax), y_min_(ymin), y_max_(ymax), 
    theta_min_(theta_min), theta_max_(theta_max), num_rots_(num_rots),
    gv_(gv), alpha_max_(alpha_max), alpha_res_(alpha_res), level_set_(level_set)
{
  cage_prover_ = NULL;
  path_prover_ = NULL;
}

Triangulation_Planner::~Triangulation_Planner()
{
  collision_triangulation_.clear();
  delete cage_prover_;
  delete path_prover_;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Triangulation_Planner::Clear_Stored_Data
 *  Description:  This class clears all of the data associated with the Planner.  After
 *                this function is called, the planner will be reset to initialized state.
 * =====================================================================================
 */
void
Triangulation_Planner::Clear_Stored_Data()
{
  collision_triangulation_.clear();
  cell_classification_.clear();
  point_in_collision_.clear();
  delete cage_prover_;
  delete path_prover_;
}

// Stores all of the boundary points for all equivalent values of theta
// TODO: check this when rescaling
void
Triangulation_Planner::Add_Boundary_Points()
{
  Weighted_point new_point;
  float theta_id;
  float boundary_alpha = 0.0f;
  boundary_points_.clear();
  free_space_points_.clear();

  // check collision for free space points (since we need them in the triangulation
  Check_Collisions(x_min_, y_min_, theta_min_);

  for (int i = -num_rots_-1; i <= num_rots_; i++) { 
    theta_id = 2 * M_PI * i;

    // add free space point outside of max alpha radius
    new_point = Weighted_point(Bare_point(x_min_, y_min_,
                                          theta_scale_ * (theta_min_ + theta_id)),
                               boundary_alpha);
    free_space_points_.push_back(new_point);

    // add max alpha to boundary points to ensure query point remain exterior to all valid alpha shape queries
    // TODO: find more elegant solution
    new_point = Weighted_point(Bare_point(x_min_ - alpha_max_, y_min_ - alpha_max_,
                                          theta_scale_ * (theta_min_ + theta_id) - alpha_max_),
                               boundary_alpha);
    boundary_points_.push_back(new_point);

    new_point = Weighted_point(Bare_point(x_max_ + alpha_max_, y_min_ - alpha_max_,
                                          theta_scale_ * (theta_min_ + theta_id) - alpha_max_),
                               boundary_alpha);
    boundary_points_.push_back(new_point);

    new_point = Weighted_point(Bare_point(x_min_ - alpha_max_, y_max_ + alpha_max_,
                                          theta_scale_ * (theta_min_ + theta_id) - alpha_max_),
                               boundary_alpha);
    boundary_points_.push_back(new_point);

    new_point = Weighted_point(Bare_point(x_min_ - alpha_max_, y_min_ - alpha_max_,
                                          theta_scale_ * (theta_max_ + theta_id) + alpha_max_),
                               boundary_alpha);
    boundary_points_.push_back(new_point);

    new_point = Weighted_point(Bare_point(x_max_ + alpha_max_, y_min_ - alpha_max_,
                                          theta_scale_ * (theta_max_ + theta_id) + alpha_max_),
                               boundary_alpha);
    boundary_points_.push_back(new_point);

    new_point = Weighted_point(Bare_point(x_max_ + alpha_max_, y_max_ + alpha_max_,
                                          theta_scale_ * (theta_min_ + theta_id) - alpha_max_),
                               boundary_alpha);
    boundary_points_.push_back(new_point);

    new_point = Weighted_point(Bare_point(x_min_ - alpha_max_, y_max_ + alpha_max_,
                                          theta_scale_ * (theta_max_ + theta_id) + alpha_max_),
                               boundary_alpha);
    boundary_points_.push_back(new_point);

    new_point = Weighted_point(Bare_point(x_max_ + alpha_max_, y_max_ + alpha_max_,
                                          theta_scale_ * (theta_max_ + theta_id) + alpha_max_),
                               boundary_alpha);
    boundary_points_.push_back(new_point);
  }
}

// Sets the scaling of theta for the alpha shape
// NOTE: This should be the maximum moment arm of the shapeo
void
Triangulation_Planner::Set_Theta_Scaling(float scale)
{
  theta_scale_ = scale;
}

float
Triangulation_Planner::Potential_Energy(float x1, float y1, float theta1,
                                        float x2, float y2, float theta2,
                                        float& distance_sq)
{
  if (potential_ == ENERGY_GRAVITY) {
    float dh = x2 - x1;
    float energy_diff = object_->Mass() * GRAVITY_ACCEL * dh;
    float x_boundary = x1 + energy_thresh_ / (object_->Mass() * GRAVITY_ACCEL);
    distance_sq = x2 - x_boundary;
    distance_sq = distance_sq * distance_sq;
    //    std::cout << energy_diff << std::endl;
    return energy_diff;
  }
}

bool
Triangulation_Planner::High_Potential(float x, float y, float theta,
                                      float x_orig, float y_orig, float theta_orig)
{
  float distance_sq;
  float potential_diff = Potential_Energy(x_orig, y_orig, theta_orig, x, y, theta, distance_sq);
  if (potential_diff > energy_thresh_) {
    float theta_id;
    for (int i = -num_rots_-1; i <= num_rots_; i++) {

      theta_id = theta + 2 * M_PI * i;
      float theta_scaled = theta_scale_ * theta_id; // scale the angle by te max moment arm
      Weighted_point new_point(Bare_point(x, y, theta_scaled), distance_sq);

      // insert into collision space triangulation
      collision_triangulation_.insert(new_point);

#ifdef DISPLAY_SPHERES
      gv_ << CGAL::GREEN;
      gv_ << K::Sphere_3(Bare_point(x, y, theta_scaled), distance_sq);
#endif
    }
    return true;
  }
  return false;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Triangulation_Planner::Sample_Random_Configuration
 *  Description:  This function chooses a random configuration and checks whether or not
 *                it places the robot in collision with the obstacles.  It then calculates
 *                either penetration depth or separation depth and stores this data in a map.
 * =====================================================================================
 */
bool
Triangulation_Planner::Sample_Random_Pose(float x_orig, float y_orig, float theta_orig)
{
  // sample a pose uniformly at random
  float x, y, theta, theta_id;
  x = ((float)rand()/(float)RAND_MAX) * (x_max_ - x_min_) + x_min_;
  y = ((float)rand()/(float)RAND_MAX) * (y_max_ - y_min_) + y_min_;
  theta = ((float)rand()/(float)RAND_MAX) * (theta_max_ - theta_min_) + theta_min_;

  // first check if the energy of the configuration is too high
  if (High_Potential(x, y, theta, x_orig, y_orig, theta_orig)) {
    return true;
  }

  // otherwise check for collisions
  return Check_Collisions(x, y, theta);
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Triangulation_Planner::Sample_Configuration
 *  Description:  This function takes in a sampled configuration and inserts it into the
 *                triangulation structure after det ermining appropriate distance.
 *  * =====================================================================================
 */
bool
Triangulation_Planner::Check_Collisions(float tx, float ty, float theta)
{
  Eigen::Matrix4f pose = CreatePose(tx, ty, theta);
  float distance_squared_to_valid_region = 0.0f;  
  if (tx < x_min_) {
    distance_squared_to_valid_region += (x_min_ - tx) * (x_min_ - tx);
  } 
  else if (tx > x_max_) {
    distance_squared_to_valid_region += (tx - x_max_) * (tx - x_max_);
  } 

  if (ty < y_min_) {
    distance_squared_to_valid_region += (y_min_ - ty) * (y_min_ - ty);
  } 
  else if (ty > y_max_) {
    distance_squared_to_valid_region += (ty - y_max_) * (ty - y_max_);
  }

  // get actual boundaries based on number of allowed rotations
  float actual_theta_min = -(num_rots_+1) * 2 * M_PI  + theta_min_;
  float actual_theta_max = num_rots_ * 2 * M_PI + theta_max_;

  if (theta < actual_theta_min) {
    distance_squared_to_valid_region += theta_scale_ * theta_scale_ * \
      (actual_theta_min - theta) * (actual_theta_min - theta);
  } 
  else if (theta > actual_theta_max) {
    distance_squared_to_valid_region += theta_scale_ * theta_scale_ * \
      (theta - actual_theta_max) * (theta - actual_theta_max);
  } 

  // check for collisions and get distance to nearest collision object
  MeshCollisionResult result = Check_Intersection(tx, ty, theta);
  bool collision = result.collision;
  double distance = result.distance;
  double distance_sq = distance * distance;

  // take max of joint limit collision and obstacle collision
  if (collision && sqrt(distance_squared_to_valid_region) > distance) {
    distance = sqrt(distance_squared_to_valid_region);
    distance_sq = distance * distance;
  }

  // float pen_depth_thresh = 1.0f;
  // if (distance_sq > pen_depth_thresh)
  //   distance_sq = distance_sq - pen_depth_thresh;
  distance_sq = sqrt(distance_sq);

  // add all identifications of the point to space triangulation
  float theta_id;
  for (int i = -num_rots_-1; i <= num_rots_; i++) {

    theta_id = theta + 2 * M_PI * i;
    float theta_scaled = theta_scale_ * theta_id; // scale the angle by te max moment arm
    Weighted_point new_point(Bare_point(tx, ty, theta_scaled), distance_sq);

    if (collision) {
      // insert into collision space triangulation
      collision_triangulation_.insert(new_point);

#ifdef DISPLAY_SPHERES
      gv_ << CGAL::PURPLE;
      gv_ << K::Sphere_3(Bare_point(tx, ty, theta_scaled), distance_sq);
#endif
    }
    else { 
      // insert into the free space triangulation
      free_space_triangulation_.insert(new_point);

#ifdef DISPLAY_SPHERES
      gv_ << CGAL::ORANGE;
      gv_ << K::Sphere_3(Bare_point(tx, ty, theta_scaled), distance_sq);
#endif
    }
  }
  return collision;
}

MeshCollisionResult
Triangulation_Planner::Check_Intersection(float tx, float ty, float theta)
{
  float pen_depth = 0.0f;
  float dist = FLT_MAX;
  bool collision = false;
  float time = 0.0f;

  if (object_ != NULL && obstacles_.size() > 0) {
    Eigen::Matrix4f pose = CreatePose(tx, ty, theta);
    object_->SetPose(pose);

    // test computation of penetration depth
    for (unsigned int j = 0; j < obstacles_.size(); j++) {
      MeshCollisionResult c = Mesh::UpperBoundCollision(object_, obstacles_[j]);
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

  //  std::cout << "Collision : " << collision << " distance " << dist << std::endl;
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

float
Triangulation_Planner::Penetration_Depth(const vectord& pose)
{
  // check dimensions
  if (pose.size() != 3) {
    std::cout << "Warning: Only supports 3-dimensional poses but queried with " << pose.size() << std::endl;
    return 0;
  }

  // check collision and return depth
  float tx = pose(0);
  float ty = pose(1);
  float theta = pose(2);
  MeshCollisionResult coll_result = Check_Intersection(tx, ty, theta);

  // compute penetration depth
  bool collision = coll_result.collision;
  float distance = coll_result.distance;
  float pen_depth = distance;
  if (!collision) {
    pen_depth = -pen_depth; // negative if not in collision
  }
  else {
    num_collisions_++;
  }
  //  std::cout << "Collision : "<< collision << " " << pen_depth << std::endl;

  // add all identifications of theta to the structure
  float theta_id;
  for (int i = -num_rots_-1; i <= num_rots_; i++) {
    theta_id = theta + 2 * M_PI * i;
    float theta_scaled = theta_scale_ * theta_id; // scale the angle by te max moment arm
    Weighted_point new_point(Bare_point(tx, ty, theta_scaled), distance);

    if (collision) {
      // insert into collision space triangulation
      collision_triangulation_.insert(new_point);

#ifdef DISPLAY_SPHERES
      gv_ << CGAL::PURPLE;
      gv_ << K::Sphere_3(Bare_point(tx, ty, theta_scaled), distance);
#endif
    }
    else { 
      // insert into the free space triangulation
      free_space_triangulation_.insert(new_point);

#ifdef DISPLAY_FREE_SPHERES
      gv_ << CGAL::ORANGE;
      gv_ << K::Sphere_3(Bare_point(tx, ty, theta_scaled), distance);
#endif
    }
  }  
  return pen_depth;
}

void
Triangulation_Planner::Add_Pose(const vectord& pose)
{
  // check dimensions
  if (pose.size() != 3) {
    std::cout << "Warning: Only supports 3-dimensional poses but queried with " << pose.size() << std::endl;
    return;
  }
  // do nothing until interface is changed
}

void
Triangulation_Planner::Adaptive_Sample(unsigned int num_samples)
{
  // set level set probing params
  int num_init_samples = 150;
  bopt_params par = initialize_parameters_to_default();
  par.n_iterations = num_samples - num_init_samples;
  par.random_seed = 0;
  par.verbose_level = 1;
  par.noise = 1e-10; // no actual noise but used to ensure invertability
  par.n_init_samples = num_init_samples;
  par.epsilon = 0.5;

  // set gp kernel priors
  set_kernel(&par,"kSEISO");
  par.kernel.hp_mean[0] = 0.0f;
  par.kernel.hp_std[0] = 10.0f;
  par.kernel.n_hp = 1;

  int d = 3;
  boost::function<double (const vectord& x)> f = boost::bind(&Triangulation_Planner::Penetration_Depth, this, _1);
  boost::function<void (const vectord& x)> cb = boost::bind(&Triangulation_Planner::Add_Pose, this, _1);
  LevelSetProber ols(par, level_set_, f, d, cb);

  vectord lower_bound(d);
  vectord upper_bound(d);
  lower_bound(0) = x_min_;
  lower_bound(1) = y_min_;
  lower_bound(2) = theta_min_;
  upper_bound(0) = x_max_;
  upper_bound(1) = y_max_;
  upper_bound(2) = theta_max_;
  ols.setBoundingBox(lower_bound, upper_bound);

  num_collisions_ = 0; // logging for coll ratio
  vectord result(d);
  ols.optimize(result);
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Triangulation_Planner::Find_Escape_Path
 *  Description:  This function computes a path to any of the exterion vertices and returns that no path exists if it
 *                detects that none exists.
 * =====================================================================================
 */
void
Triangulation_Planner::Find_Escape_Path(float tx, float ty, float theta, PathResult& path_result,
                                        unsigned int num_samples, PotentialType potential, float energy_thresh,
                                        bool prove_exists)
{ 
  Eigen::Matrix4f pose = CreatePose(tx, ty, theta);
  potential_ = potential;
  energy_thresh_ = energy_thresh;

  // set initial random samples
  int size = collision_triangulation_.number_of_vertices();

  Add_Boundary_Points();
  Check_Collisions(tx, ty, theta); // add samples of the query pose
  num_collisions_ = 0;
  if (size == 0) {
    // add boundary points to triangulation (in order for the convex hull to contain the workspace)
    for (unsigned int i = 0; i < boundary_points_.size(); i++) {
      collision_triangulation_.insert(boundary_points_[i]);
    }

    // sample a bunch of random points inside the rectangloid to begin the search
    timer_.start();
    //    Adaptive_Sample(num_samples);

    for(int i = 0; i < num_samples; i++) {
      bool collision = Sample_Random_Pose(tx, ty, theta);
      if (collision) {
        num_collisions_++;
      }
    }

    timer_.stop();
    path_result.sample_time = timer_.time();
    timer_.reset();
  }

  path_result.coll_ratio = (float)num_collisions_ / (float)num_samples;
  float free_ratio = 1.0f - path_result.coll_ratio;

  // create alpha shape for collision space
  delete cage_prover_;
  timer_.start();
  cage_prover_ = new Alpha_shapes_disconnection(collision_triangulation_);
  timer_.stop();
  //  std::cout << "Cage alpha shape constructed in: " << timer_.time() << " seconds." << std::endl;
  path_result.alpha_time = timer_.time();
  timer_.reset();


  // create new filtration
  std::vector<CGAL::Object> filtered_simplices;
  std::vector<K::FT> alpha_values;
  cage_prover_->alpha_shape_.filtration_with_alpha_values(CGAL::Dispatch_output_iterator<CGAL::cpp11::tuple<CGAL::Object, K::FT>,
                                                                                         CGAL::cpp11::tuple<std::back_insert_iterator<std::vector<CGAL::Object> >,
                                                                                                            std::back_insert_iterator<std::vector<K::FT> > > >(std::back_inserter(filtered_simplices),
                                                                                                                                                               std::back_inserter(alpha_values))); 
  std::cout << "Num simplices " << filtered_simplices.size() << std::endl;
  timer_.start();
  Alpha_Shape_Filtration_3 filtration(cage_prover_->alpha_shape_);
  filtration.save_boundary_matrix("alpha_filtration.bin");
  timer_.stop();
  std::cout << "Boundary matrix computation took " << timer_.time() << " seconds" << std::endl;
  timer_.reset();
  while(true);

  // iterate through 
  Triangulation_3::Finite_vertices_iterator it;
  int ind = 0;
  for (it = cage_prover_->alpha_shape_.finite_vertices_begin(); it != cage_prover_->alpha_shape_.finite_vertices_end(); it++) {
    //std::cout << it->point()[0] << std::endl;
    ind++;
  }

  Triangulation_3::Finite_edges_iterator edge_it;
  for (edge_it = cage_prover_->alpha_shape_.finite_edges_begin(); edge_it != cage_prover_->alpha_shape_.finite_edges_end(); edge_it++) {
    //    std::cout << Triangulation_3::segment(*edge_it) << std::endl;
    ind++;
  }

  while(true);

  // create alpha shape for free space
  if (prove_exists) {
    delete path_prover_;
    timer_.start();
    path_prover_ = new Alpha_shapes_connection(free_space_triangulation_);
    timer_.stop();
    //    std::cout << "Path alpha shape constructed in: " << timer_.time() << " seconds." << std::endl;
    timer_.reset();
  }

  // check for a path to any of the points that are guaranteed to be in free space
  float query_theta_scaled = theta_scale_ * theta;
  Bare_point start_point = Bare_point(tx, ty, query_theta_scaled);
  Bare_point end_point;
  bool path_exists = true; // default to path exists
  Cell_handle start_cell, end_cell;
  start_cell = collision_triangulation_.locate(start_point); //find start cell //TODO: possible bug here to due with the locate query returning something 

  // debugging only!
  // float max_area = MAX_AREA;
  // float conv_hull_area = 0.0f;
  // float conv_hull_volume = 0.0f;
  // std::vector<DifferentialTriangle> separator_tris;
  // cage_prover_->Escape_Boundary(start_point, alpha_max_, max_area, separator_tris, conv_hull_area, conv_hull_volume, gv_);

  // // plan paths to all DTs
  // PathPlanningParams params;
  // params.start_tx = tx;
  // params.start_ty = ty;
  // params.start_theta = theta;
  // params.goal_tx = x_min_ + 0.1;
  // params.goal_ty = y_min_ + 0.1;
  // params.goal_theta = 0;

  // params.min_tx = x_min_;
  // params.min_ty = y_min_;
  // params.min_theta = 2 * M_PI * (-num_rots_ - 1);

  // params.max_tx = x_max_;
  // params.max_ty = y_max_;
  // params.max_theta = 2 * M_PI * (num_rots_ + 1);

  // params.tx_scale = 1.0f;
  // params.ty_scale = 1.0f;
  // params.theta_scale = 1.0f / theta_scale_;

  // params.goal_dist_thresh = GOAL_DIST_THRESH;

  // CageEscapePlanner cep(object_, obstacles_);
  // float space_volume = (alpha_max_ + x_scale_ * x_max_ - (x_scale_ * x_min_ + alpha_max_)) * \
  //   (alpha_max_ + y_scale_ * y_max_ - (y_scale_ * y_min_ - alpha_max_)) *                     \
  //   (alpha_max_ + theta_scale_ * params.max_theta - (theta_scale_ * params.min_theta - alpha_max_));
  // float volume_ratio = conv_hull_volume / space_volume;
  // float timeout = 1e-2;

  // timer_.start();
  // float escape_energy = cep.IntegrateEscapePathEnergy(params, separator_tris, timeout);
  // escape_energy = escape_energy / conv_hull_area;
  // timer_.stop();
  // path_result.energy = escape_energy;
  // path_result.conv_volume_ratio = volume_ratio;
  // path_result.energy_time = timer_.time();


  // iterate over "interesting" values of alpha and check the cage condition
  Alpha_iterator alpha_it;
  float alpha;
  int k = 0;

  // BELOW: extra code to iterate over all interesting alpha
  timer_.start();
  for (alpha = 0.0f; path_exists && alpha < alpha_max_; alpha += alpha_res_) {
    cage_prover_->Set_Alpha(alpha);

    for (unsigned int i = 0; i < free_space_points_.size() && path_exists; i++) {
      // locate the start and end points (this is for path planning, so probably will be removed)
      end_point = free_space_points_[i];
      end_cell = collision_triangulation_.locate(end_point); //find end cell   //other than a cell.  also, likely has a bug when the cell returned is a 
    
      // check for disconnection between the two points
      path_exists = cage_prover_->Query_disconnection(start_point, end_point);
    }
    k++;
  }
  alpha = alpha - alpha_res_; // undo the add from the final loop
  timer_.stop();
  path_result.iter_time = timer_.time();
  timer_.reset();

  cage_prover_->Component_Volume(start_point, path_result.volume);

  // optional display alpha shape
#ifdef DISPLAY_AS
  cage_prover_->Set_Alpha(0.0f);
  gv_.clear();
  gv_ << CGAL::DEEPBLUE;
  gv_ << cage_prover_->alpha_shape_;
  gv_ << CGAL::PURPLE;
  
  gv_ << CGAL::GREEN;
  gv_ << Convert_Cell_To_Tetrahedron(start_cell);
  gv_ << CGAL::RED;
  gv_ << Convert_Cell_To_Tetrahedron(end_cell);      
#endif
  
  // path set to exist if alpha is nonzero
  path_exists = (alpha > 0);
  path_result.alpha = alpha;
  path_result.exists = path_exists;   // if path does not exist, then it truly must not exist
  start_cell = free_space_triangulation_.locate(start_point);

  // attempt to prove existence if specified
  // TODO: this doesn't work right now, need to track interior cells ONLY
  bool free_space_path_exists = false;
  if (prove_exists) {
    for (unsigned int i = 0; i < free_space_points_.size() && !free_space_path_exists; i++) {
      // locate the start and end points (this is for path planning, so probably will be removed)
      end_point = free_space_points_[i];
      end_cell = free_space_triangulation_.locate(end_point); //find end cell   //other than a cell.  also, likely has a bug when the cell returned is a 
    
      // check for disconnection between the two points
      free_space_path_exists = path_prover_->Query_connection(start_point, end_point);

#ifdef DISPLAY_AS
      if (i == 0) {
        gv_.clear();
        gv_ << CGAL::DEEPBLUE;
        gv_ << path_prover_->alpha_shape_;
        gv_ << CGAL::PURPLE;
      
        gv_ << CGAL::GREEN;
        gv_ << Convert_Cell_To_Tetrahedron(start_cell);
        gv_ << CGAL::RED;
        gv_ << Convert_Cell_To_Tetrahedron(end_cell);      
      }
#endif
    }
    //    std::cout << "Free space path exists is " << free_space_path_exists << std::endl;
    path_result.exists = free_space_path_exists;
  }
}

Alpha_status
Triangulation_Planner::classify_cell(Cell_handle& s)
{
  if(!cell_map_.is_defined(s)) {
    double alpha = Compute_squared_radius(collision_triangulation_, s);
    cell_map_[s] =  alpha;
    s->set_alpha(alpha);
  }
  
  if (collision_triangulation_.is_infinite(s)) return EXTERIOR;
  return (s->get_alpha() <=  0.0) ? INTERIOR : EXTERIOR;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Triangulation_Planner::Sample_Tetrahedron_Configurations
 *  Description:  This function takes in a Cell_path, assumed to be composed of all mixed
 *                cells, and then subdivides in each one by sampling a random point within
 *                each cell.
 * =====================================================================================
 */
void
Triangulation_Planner::Sample_Tetrahedron_Configurations(Cell_path& path)
{
    std::list<Configuration> sampled_points;
    // sample all the points first so that the cells are still in the triangulation
    for(Cell_path::iterator iter = path.begin(); iter != path.end(); iter++)
    {
      sampled_points.push_back(Random_Point_In_Tetrahedron(*iter));
    }
    // then insert all of them into the triangulation

    // TODO: make the below work with poses
    // for(std::list<Configuration>::iterator iter = sampled_points.begin(); iter != sampled_points.end(); iter++)
    // {
    //   Sample_Configuration(*iter);
    // }
}		/* -----  end of function Triangulation_Planner::Sample_Tetrahedron_Configurations  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Triangulation_Planner::Random_Point_In_Tetrahedron
 *  Description:  Generates a uniform point in a tetrahedron using barycentric coordinates.
 *                Method from Generating Random Points in a Tetrahedron by C. Rocchini,
 *                P. Cignoni
 * =====================================================================================
 */
Configuration
Triangulation_Planner::Random_Point_In_Tetrahedron (Cell_handle in_cell)
{
    double s, t, u;
    Configuration conf;
    //generate 3 random values between 0 and 1
    s = ((double)rand()/(double)RAND_MAX);
    t = ((double)rand()/(double)RAND_MAX);
    u = ((double)rand()/(double)RAND_MAX);


    double new_s, new_t, new_u;
    //fold the sampled point into a tetrahedron from a cube
    //resulting doubles are barycentric coordinates for a tetrahedron
    if(s+t+u > 1) {
      if(t+u > 1) {
        new_s = s;
        new_t = 1-u;
        new_u = 1-s-t;
      }
      else {
        new_s = 1-t-u;
        new_t = t;
        new_u = s+t+u-1;
      }
    }
    else {
      new_s = s;
      new_t = t;
      new_u = u;
    }

    Bare_point base(in_cell->vertex(0)->point());
    Bare_point a(in_cell->vertex(1)->point());
    Bare_point b(in_cell->vertex(2)->point());
    Bare_point c(in_cell->vertex(3)->point());

    Bare_point end = base + (new_s*(a-base) + new_t*(b-base) + new_u*(c-base));
    conf.push_back(end[0]);
    conf.push_back(end[1]);
    conf.push_back(end[2]);
    return conf;
}

Configuration Triangulation_Planner::Centroid(Cell_handle in_cell) {
  Configuration conf;
  //prep return
  double array[] = {0.0, 0.0, 0.0};
  //calculate centroid
  for(int i = 0; i < 4; i++)
    {
      for(int j = 0; j < 3; j++)
        {
          array[j] += .25*in_cell->vertex(i)->point()[j];
        }
    }
  for(int i = 0; i < 3; i++)
    {
      conf.push_back(array[i]);
    }
  return conf;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Convert_Cell_To_Tetrahedron
 *  Description:  Takes in an input cell_handle and returns the tetrahedron for gv
 * =====================================================================================
 */
Tetrahedron
Convert_Cell_To_Tetrahedron ( Cell_handle in_cell )
{
  return Tetrahedron(in_cell->vertex(0)->point(), in_cell->vertex(1)->point(),
                     in_cell->vertex(2)->point(), in_cell->vertex(3)->point());
}		/* -----  end of function Convert_Cell_To_Tetrahedron  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Compute_squared_radius
 *  Description:  This function takes in a cell and computes the orthogonal circumsphere
 *                of it.  I haven't dealt at all with degeneracy. TODO: FIX THAT!
 * =====================================================================================
 */
double
Compute_squared_radius (Triangulation_3& dt, Cell_handle& cell_in)
{return dt.geom_traits().compute_squared_radius_smallest_orthogonal_sphere_3_object()(
                                                                                      dt.point(cell_in,0), dt.point(cell_in,1),
                                                                                      dt.point(cell_in,2), dt.point(cell_in,3));
}		/* -----  end of function Compute_squared_radius  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Compute_squared_radius
 *  Description:  This function takes in a cell and computes the orthogonal circumsphere
 *                of it.  I haven't dealt at all with degeneracy. TODO: FIX THAT!
 * =====================================================================================
 */
double
Compute_squared_radius_facet (Triangulation_3& dt, Cell_handle& cell_in, int neighbor)
{   int i[3];
  int j = 0;
  for(int m = 0; m < 4; m++)
    {
      if(neighbor != m)
        {
          i[j++] = m;
        }
    }
  return dt.geom_traits().compute_squared_radius_smallest_orthogonal_sphere_3_object()(
                                                                                       dt.point(cell_in,i[0]), dt.point(cell_in,i[1]),
                                                                                       dt.point(cell_in,i[2]));
}		/* -----  end of function Compute_squared_radius  ----- */
