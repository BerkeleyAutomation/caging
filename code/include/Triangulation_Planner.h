/*
 * =====================================================================================
 *
 *       Filename:  Triangulation_Planner.h
 *
 *    Description:  This file contains the code for a planner that plans using Alpha Shapes
 *                  and maintains a triangulation of the sampled point set.
 *
 *        Version:  1.0
 *        Created:  07/26/2011 09:32:52 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Zoe McCarthy (zm), ZoeMcCarthy12@gmail.com
 *        Company:  University of Illinois at Urbana-Champaign
 *
 * =====================================================================================
 */


#ifndef TRIANGULATION_PLANNER_H
#define TRIANGULATION_PLANNER_H

#ifndef DEBUG
#define DEBUG
#endif

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
// #define DISPLAY_AS
#endif

//#include "Planner.h"
#include "CollisionChecker.h"
#include "ConfigurationMapper.h"
#include "Alpha_shapes_disconnection.h"
#include "Alpha_shapes_connection.h"
#include "Mesh.hpp"
#include "LevelSetProber.hpp"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Union_find.h>
#include <CGAL/Timer.h>
#include <CGAL/utility.h>
#include <CGAL/IO/Geomview_stream.h>

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


enum Alpha_Classification {DEFAULT, EMPTY, FULL, MIXED};

enum Alpha_status { EXTERIOR, INTERIOR };

enum PotentialType {ENERGY_GRAVITY};

typedef std::list< Cell_handle >                            Cell_path;
typedef std::vector<Cell_handle>                            Cell_vector;
typedef std::pair< Cell_path, Alpha_Classification >        Cell_path_with_exist;

typedef CGAL::Unique_hash_map<Vertex_handle, bool>         Point_map;
typedef CGAL::Unique_hash_map<Cell_handle, Alpha_Classification>  Cell_class_map;
typedef CGAL::Unique_hash_map<Cell_handle, double>         Alpha_cell_map;
typedef CGAL::Unique_hash_map<Facet, double>               Alpha_facet_map;


Tetrahedron Convert_Cell_To_Tetrahedron ( Cell_handle in_cell );
double Compute_squared_radius (Triangulation_3& dt, Cell_handle& cell_in);
double Compute_squared_radius_facet (Triangulation_3& dt, Cell_handle& cell_in, int neighbor);

struct PathResult {
  bool  exists;            // whether or not a path exists out of the cage
  float alpha;             // alpha value at which the object becomes caged
  float volume;            // volume of the cage interior
  float energy;            // uper bound on energy required to escape cage
  float conv_volume_ratio; // ratio of cage convex hull volume to universe volume (for energy re-normalization)
  float coll_ratio;        // number of samples in collision
  float sample_time;       // time to sample poses
  float alpha_time;        // time to construct alpha shape
  float iter_time;         // time to iterate through alpha values
  float energy_time;       // time to compute cage "energy"
};

/*
 * =====================================================================================
 *        Class:  Triangulation_Planner
 *  Description:  This class is an implementation of the Planner class.  It is a PRM planner
 *                that utilizes penetration depth and separation depth calculations while 
 *                sampling in order to guide sampling and inform of path existence.
 * =====================================================================================
 */
class Triangulation_Planner
{
  friend class LevelSetProber;
  friend class Path_Tree;

 public:
  /* ====================  LIFECYCLE     ======================================= */
  Triangulation_Planner (float xscale, float yscale, float theta_scale, 
                         float xmin, float xmax, float ymin, float ymax,
                         float theta_min, float theta_max, int num_rots,
                         CGAL::Geomview_stream& gv,
                         float alpha_max = 1.0f, float alpha_res = 1e-3, float level_set = 0.0f);          /* constructor */
  ~Triangulation_Planner();
  void Clear_Stored_Data();

  void Set_Object(Mesh* object) {object_ = object;}
  void Set_Obstacles(std::vector<Mesh*> obstacles) {obstacles_ = obstacles;}

  // main functionality
  void Find_Escape_Path(float tx, float ty, float theta, PathResult& path_result, unsigned int num_samples = 1e4,
                        PotentialType potential = ENERGY_GRAVITY, float energy_thresh = 1000.0f, bool prove_exists = false); //find a path between configurations

  // todo: make private
  Point_map  point_in_collision_; 
  CGAL::Geomview_stream&      gv_;

 public:
  void Set_Theta_Scaling(float scale);

  // computes the potential energy diff btw x2 and x1
  float Potential_Energy(float x1, float y1, float theta1,
                         float x2, float y2, float theta2, float& distance_sq);
  bool High_Potential(float x, float y, float z,
                      float x_orig, float y_orig, float theta_orig);

 private:
  // private functions to find paths
  void Add_Boundary_Points();

  //  Cell_vector Find_Path_Between_And_Classify(Bare_point qs, Bare_point qg);

  // Classifies the cell `f' of the underlying Delaunay
  // tetrahedralization with respect to `A'.
  // s->radius == alpha => f interior
  Alpha_status classify_cell(Cell_handle& s);
  
  // samples a random pose within bounds
  bool Sample_Random_Pose(float x_orig, float y_orig, float theta_orig);
  bool Check_Collisions(float tx, float ty, float theta);
  MeshCollisionResult Check_Intersection(float tx, float ty, float theta);

  // functions for level set probing
  void Adaptive_Sample(unsigned int num_samples);
  float Penetration_Depth(const vectord& pose);   // negative if not in collision
  void Add_Pose(const vectord& pose);

 private:
  Configuration Random_Point_In_Tetrahedron(Cell_handle in_cell);
  Configuration Centroid(Cell_handle in_cell);
  void Sample_Tetrahedron_Configurations(Cell_path& in_Path);

 private:
  /* ====================  DATA MEMBERS  ======================================= */
  Triangulation_3             collision_triangulation_;
  Triangulation_3             free_space_triangulation_;
  Alpha_shapes_disconnection* cage_prover_;
  Alpha_shapes_connection*    path_prover_;
  Timer                       timer_;

  // currently only allows for linear scaling
  float x_scale_;
  float y_scale_;
  float theta_scale_;

  // min/max are in PLAN_CONF units
  float x_min_;
  float x_max_;
  float y_min_;
  float y_max_;
  float theta_min_;
  float theta_max_;

  // params for stepping over continuous alpha values
  float alpha_max_;
  float alpha_res_;

  std::vector<Weighted_point> boundary_points_;
  std::vector<Weighted_point> free_space_points_;
  int num_rots_;
  float level_set_;  // not sure what this is
  int num_collisions_;

  PotentialType potential_;
  float energy_thresh_;

  // alpha shape factors
  Cell_class_map cell_classification_; 
  Alpha_facet_map facet_map_;
  Alpha_cell_map  cell_map_;
  Mesh* object_;
  std::vector<Mesh*> obstacles_;

}; /* -----  end of class Triangulation_Planner  ----- */


#endif // TRIANGULATION_PLANNER_H
