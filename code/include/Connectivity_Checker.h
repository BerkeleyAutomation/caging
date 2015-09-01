#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Union_find.h>
#include <CGAL/Timer.h>

#include <iostream>
#include <list>
#include <map>

#include "Filtration.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Regular_triangulation_euclidean_traits_3<K> Gt;

typedef CGAL::Timer                         Timer;

typedef CGAL::Alpha_shape_vertex_base_3<Gt>         Vb;
typedef CGAL::Alpha_shape_cell_base_3<Gt>           Fb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb> Tds;
typedef CGAL::Regular_triangulation_3<Gt,Tds>       Triangulation_3;
typedef CGAL::Alpha_shape_3<Triangulation_3>        Alpha_shape_3;
typedef CGAL::Tetrahedron_3<K>                      Tetrahedron_3;
typedef CGAL::Triangle_3<K>                         Tri_3;
typedef CGAL::Point_3<K>                            Pt_3;
typedef CGAL::Direction_3<K>                        Dir_3;

typedef Triangulation_3::Locate_type    Locate_type;

typedef Alpha_shape_3::Cell_handle          Cell_handle;
typedef Alpha_shape_3::Vertex_handle        Vertex_handle;
typedef Alpha_shape_3::Facet                Facet;
typedef Alpha_shape_3::Edge                 Edge;
typedef Alpha_shape_3::Alpha_iterator       Alpha_iterator;
typedef Gt::Weighted_point                  Weighted_point;
typedef Gt::Bare_point                      Bare_point;

typedef Cell_handle                         T;
typedef CGAL::Union_find<T>                 Disjoint_set;  
typedef Disjoint_set::handle                handle;

typedef CGAL::Unique_hash_map<Cell_handle, Cell_handle>     Cell_map;
typedef std::vector<Cell_handle>                            Cell_vector;

////// ********************    CLASS DEFINITION    ********************/

// this class will allow querying of two points to test 
// whether or not they are disconnected
// it uses a filtration to determine which areas of the space are inside / outside of the shape
class Connectivity_Checker
{
 public:
  // construct the disjoint set that contains the connectivity information
  Connectivity_Checker(Triangulation_3& raw_triangulation, Filtration_3& filtration, CGAL::Geomview_stream& gv, Potential potential_thresh = 0.0f);

  // query for connections / disconnections between two points
  bool connected(Bare_point qs, Bare_point qg);
  bool connected(Bare_point qs, Bare_point qg, Potential potential);
  bool disconnected(Bare_point qs, Bare_point qg);
  bool disconnected(Bare_point qs, Bare_point qg, Potential potential);

  // get the volume of the connected component containing the point q
  bool connected_component_volume(Bare_point q, float& volume);

  // set a new potential threshold
  inline Potential get_potential() { return potential_thresh_; } 
  void set_potential(Potential p);

 private:
  // for resetting the internal connectivity structures
  void reset_classification();
  void compute_connected_components();

  // return the number of sets
  int num_connected_components() {
    return disjoint_set_.number_of_sets();
  }

 private:
  Triangulation_3& triangulation_;
  Filtration_3& filtration_;
  Potential potential_thresh_;
  std::vector<Cell_handle> external_cells_;

  CGAL::Geomview_stream& gv_;

  Timer timer_;
  Disjoint_set disjoint_set_;
  std::map<Cell_handle, handle> map_; //map converts cell_handles from location queries to ds structure handles
};
