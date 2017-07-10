#pragma once
#define NOMINMAX

#ifndef DISPLAY_TETRA
//#define DISPLAY_TETRA
#endif

#ifndef DISPLAY_SPHERES
//#define DISPLAY_SPHERES
#endif

#ifndef DISPLAY_FREE_SPHERES
//#define DISPLAY_FREE_SPHERES
#endif

#ifndef DISPLAY_AS
//#define DISPLAY_AS
#endif

#include "Typedef.h"
#include "Connectivity_Checker.h"
#include "Mesh.hpp"

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>


struct Pose2D {
  Pose2D() {}
  Pose2D(float xp, float yp, float thetap) {x=xp; y=yp; theta=thetap;}
  Pose2D(Pose2D& other) {x=other.x; y=other.y; theta=other.theta;}

  float x;
  float y;
  float theta;
};

struct PoseParticle {
  Bare_point p;
  float weight;
};

struct CfApproxConfig {
  float x_scale;
  float y_scale;
  float theta_scale;
  float x_min;
  float x_max;
  float y_min;
  float y_max;
  float theta_min;
  float theta_max;
  int num_rots;
};

// config for computing the "escape energy"
struct EscapeEnergyConfig {
  /* EscapeEnergyConfig(unsigned int n, PotentialType pt = ENERGY_GRAVITY, bool pe = false) : num_samples(n), potential_type(pt), prove_exists(pe) {}  */ 
  float energy_min; // min energy to search over
  float energy_max; // max energy to search over
  float energy_res;
  unsigned int num_samples;
  PotentialType potential_type;
  bool prove_exists;
};

struct EscapeEnergyResult {
  bool  path_exists;            // whether or not a path exists from start to goal
  float energy;            // uper bound on energy required to escape cage
  float normalized_energy; // energy normalized by mass and moment arm 
  float volume;            // volume of the cage interior
  float coll_ratio;        // number of samples in collision
  float sample_time;       // time to sample poses
  float triangulation_time;// time to construct triangulation and filtration
  float iter_time;         // time to iterate through alpha values
};

//Struct that is used by the synthesis algorithm to return important information on the best persisence pair.
struct synthesis_info {
  //synthesis_info() { 
  //} 

  synthesis_info(Pose2D selected_pose, Pose2D birth_pose,
                 Potential pose_delta, float norm_delta, float v_volume, float v_area,
                 float sample_time, float triangulation_time,
                 float filtration_time,
  		 float persistence_time, float void_time, float rotation){
    pose_ = selected_pose;
    birth_pose_ = birth_pose;
    delta_ = pose_delta;
    norm_delta_ = norm_delta;
    v_volume_ = v_volume;
    v_area_ = v_area;
    sample_time_ = sample_time;
    triangulation_time_ = triangulation_time;
    filtration_time_ = filtration_time;
    persistence_time_ = persistence_time;
    void_time_ = void_time;
    rotation_ = rotation;
  }

  synthesis_info(const synthesis_info& obj) {
    pose_ = obj.pose_;
    birth_pose_ = obj.birth_pose_;
    delta_ = obj.delta_;
    norm_delta_ = obj.norm_delta_;
    v_volume_ = obj.v_volume_;
    v_area_ = obj.v_area_;
    sample_time_ = obj.sample_time_;
    triangulation_time_ = obj.triangulation_time_;
    filtration_time_ = obj.filtration_time_;
    persistence_time_ = obj.persistence_time_;
    void_time_ = obj.void_time_;
    rotation_ = obj.rotation_;
  }
  
  Pose2D get_pose() {
    return pose_;
  }

  Pose2D get_birth_pose() {
    return birth_pose_;
  }
  
  Potential get_delta() const {
    return delta_;
  }

  Potential get_norm_delta() const {
    return norm_delta_;
  
}
  float get_void_volume() {
    return v_volume_;
  }

  float get_sample_time() {
    return sample_time_;
  }

  float get_triangulation_time() {
    return triangulation_time_;
  }

  float get_filtration_time() {
    return filtration_time_;
  }

  float get_persistence_time() {
    return persistence_time_;
  }

  float get_void_time() {
    return void_time_;
  }

  float get_rotation() {
    return rotation_;
  }
  
  bool operator<(const synthesis_info& other) const {
    return delta_ < other.get_delta();
  } 

public: 
  Pose2D pose_;
  Pose2D birth_pose_;
  Potential delta_;
  Potential norm_delta_;
  float v_volume_;
  float v_area_;
  float sample_time_;
  float triangulation_time_;
  float filtration_time_;
  float persistence_time_;
  float void_time_;
  float rotation_;
};

// approximates the configuration space of a gripper and object using sampling, and can additionally check the
// connectivity of the underlying configuration space
// essentially a lean version of the old and bloated and horrible triangulation planner
class Configuration_Space_Approximator
{
 public:
  Configuration_Space_Approximator(CfApproxConfig config, CGAL::Geomview_stream& gv);
  ~Configuration_Space_Approximator();

  //Utility function to visualize a given pose.
  void visualize_pose(Pose2D pose);

  //Check that a grasp is actually reachable in the context of push caging....
  bool check_push_reach(Pose2D pose, Eigen::Vector2d force_direction, float disc, float backup);

   
  float void_opening_area(std::vector<ASE::Simplex> void_simplices, Bare_point birth_vertex, AlphaCustomPotential potential_func, Eigen::Vector2d force_vector, Potential birth_energy); 
  float simplex_projected_area(ASE::Simplex cur_simplex, Bare_point birth_vertex, AlphaCustomPotential potential_func, Eigen::Vector2d force_vector, Potential birth_energy);
  
  // main functionality
  std::vector< std::vector<synthesis_info> > synthesize_grasps(int num_searches, Pose2D pose, EscapeEnergyConfig config, 
  							     EscapeEnergyResult& result, float angle_sweep, float angle_disc, bool check_reachability = true, float max_push_force = 1.0f);

  // set the object and obstacles
  void set_object(Mesh* object) {object_ = object;}
  void set_obstacles(std::vector<Mesh*> obstacles) {obstacles_ = obstacles;}

 public:
  // some setters for convenienve
  void set_x_scale(float scale) { config_.x_scale = scale; reset(); }
  void set_y_scale(float scale) { config_.y_scale = scale; reset(); }
  void set_theta_scale(float scale) {config_.theta_scale = scale; reset(); }

  // reset the object
  void reset();


 private:
  // adds points at infinity
  void compute_cf_padding();

  // conversions between points and poses
  void pose_to_point(Pose2D pose, float radius_sq, Weighted_point& point);
  void point_to_pose(Weighted_point point, Pose2D& pose);
  void point_3_to_pose(Point_3 point, Pose2D& pose);

  // adding to triangulations
  bool add_to_triangulation(Pose2D pose, float radius_sq, Triangulation_3& triangulation, CGAL::Color c);
  bool add_to_collision_triangulation(Pose2D pose, float radius_sq);
  bool add_to_free_space_triangulation(Pose2D pose, float radius_sq);

  // adds the infinite points
  void add_points_at_infinity();

  // samples random poses and checks collisions
  unsigned int sample_poses(Pose2D pose_orig, unsigned int num_samples);
  MeshCollisionResult sample_random_pose(Pose2D& pose);
  void uniform_random_pose(Pose2D& pose);
  MeshCollisionResult uniform_collision_pose(Pose2D& pose);
  void gaussian_prm_pose(Pose2D& pose);
  void level_set_random_pose(Pose2D& pose);

  bool check_collisions(Pose2D pose);
  //Check collisions without adding to a triangulation....
  bool check_collisions_safe(Pose2D pose);
  MeshCollisionResult check_object_obstacle_intersection(Pose2D pose, bool print = false, bool upper=false);
  float void_volume(std::vector<ASE::Simplex> void_simplices);
 
 private:
  // core data members
  Triangulation_3             collision_triangulation_;  // triangulation of points in collision
  Triangulation_3             free_space_triangulation_; // triangulation of points out of collision
  std::vector<Weighted_point> points_at_infinity_;       // vector of points on the boundary outside of the bounded region of conf space
  Connectivity_Checker*       connectivity_checker_;     // checks the connectivity of a filtration and base triangulation

  // sampling
  boost::mt19937 rng_;
  boost::normal_distribution<float> gprm_norm_;
  boost::variate_generator<boost::mt19937, boost::normal_distribution<float> > gprm_distribution_;  

  // debugging
  Timer                       timer_;
  CGAL::Geomview_stream&      gv_;

  // configuration parameters for the class (see above)
  CfApproxConfig config_;
  float cf_padding_; // amount of excess around configuration space for points "at infinity"
  int num_collisions_;

  // meshes for the collision space
  Mesh* object_;
  std::vector<Mesh*> obstacles_;
  float max_x_val_;
};
