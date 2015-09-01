#pragma once

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

#include "Connectivity_Checker.h"
#include "Mesh.hpp"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Union_find.h>
#include <CGAL/Timer.h>
#include <CGAL/utility.h>
#include <CGAL/IO/Geomview_stream.h>

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

#include <algorithm>
#include <iostream>
#include <list>
#include <vector>
#include <queue>
#include <math.h>

typedef std::vector<double>                                 Configuration;
typedef std::list< Configuration >                          Path;
typedef std::pair< Path, bool >                             Path_with_exist;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Regular_triangulation_euclidean_traits_3<K>   Gt;
typedef CGAL::Timer                                         Timer;

typedef CGAL::Alpha_shape_vertex_base_3<Gt>         Vb;
typedef CGAL::Alpha_shape_cell_base_3<Gt>           Fb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb> Tds;
typedef CGAL::Regular_triangulation_3<Gt,Tds>       Triangulation_3;
typedef CGAL::Triangulation_3<K>       Simple_Triangulation_3;

typedef Triangulation_3::Cell_handle                        Cell_handle;
typedef Triangulation_3::Vertex_handle                      Vertex_handle;
typedef Triangulation_3::Facet                              Facet;
typedef Triangulation_3::Edge                               Edge;
typedef Gt::Weighted_point                                  Weighted_point;
typedef Gt::Bare_point                                      Bare_point;

typedef K::Tetrahedron_3                                    Tetrahedron;

typedef std::list< Cell_handle >                            Cell_path;
typedef std::vector<Cell_handle>                            Cell_vector;
typedef std::pair< Cell_path, Alpha_Classification >        Cell_path_with_exist;

typedef CGAL::Unique_hash_map<Vertex_handle, bool>         Point_map;
typedef CGAL::Unique_hash_map<Cell_handle, Alpha_Classification>  Cell_class_map;
typedef CGAL::Unique_hash_map<Cell_handle, double>         Alpha_cell_map;
typedef CGAL::Unique_hash_map<Facet, double>               Alpha_facet_map;

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

// approximates the configuration space of a gripper and object using sampling, and can additionally check the
// connectivity of the underlying configuration space
// essentially a lean version of the old and bloated and horrible triangulation planner
class Configuration_Space_Approximator
{
 public:
  Configuration_Space_Approximator(CfApproxConfig config, CGAL::Geomview_stream& gv);
  ~Configuration_Space_Approximator();

  // main functionality
  bool min_escape_energy(Pose2D pose, EscapeEnergyConfig config, EscapeEnergyResult& result); //find a path between configurations

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

  // adding to triangulations
  bool add_to_triangulation(Pose2D pose, float radius_sq, Triangulation_3& triangulation, CGAL::Color c);
  bool add_to_collision_triangulation(Pose2D pose, float radius_sq);
  bool add_to_free_space_triangulation(Pose2D pose, float radius_sq);

  // adds the infinite points
  void add_points_at_infinity();

  // samples random poses and checks collisions
  unsigned int sample_poses(Pose2D pose_orig, unsigned int num_samples);
  void sample_random_pose(Pose2D& pose);
  void uniform_random_pose(Pose2D& pose);
  void uniform_collision_pose(Pose2D& pose);
  void gaussian_prm_pose(Pose2D& pose);
  void level_set_random_pose(Pose2D& pose);

  bool check_collisions(Pose2D pose);
  MeshCollisionResult check_object_obstacle_intersection(Pose2D pose);

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

};
