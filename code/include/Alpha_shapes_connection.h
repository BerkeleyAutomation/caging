/*********************
Zoe McCarthy, May 2011

//defines  Alpha_shapes_disconnection class.
//supports an alpha shape for weighted points (spheres) in R^3
//has an auxillary structure for querying two points and determining
//if they are in the same connected component of the free space

// SAME as alpha_shapes_disconnection but for interior cells
********************/

#ifndef ALPHA_CON_H
#define ALPHA_CON_H

////// ********************    PREAMBLE    ********************/

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Union_find.h>
#include <CGAL/Timer.h>

#include "Alpha_shapes_disconnection.h"

#include <iostream>
#include <list>
#include <map>

#ifndef DEBUG
#define DEBUG
#endif


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Regular_triangulation_euclidean_traits_3<K> Gt;

typedef CGAL::Timer                         Timer;

typedef CGAL::Alpha_shape_vertex_base_3<Gt>         Vb;
typedef CGAL::Alpha_shape_cell_base_3<Gt>           Fb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb> Tds;
typedef CGAL::Regular_triangulation_3<Gt,Tds>       Triangulation_3;
typedef CGAL::Alpha_shape_3<Triangulation_3>        Alpha_shape_3;

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
// its underlying structure is a disjoint set of 
// connected components of cells that are not in the alpha shape
class Alpha_shapes_connection
{
 public:
  // construct the disjoint set that contains the connectivity information
  Alpha_shapes_connection(Triangulation_3 raw_triangulation, float alpha = 0.0f);

  Cell_vector Find_Path_Between(Bare_point qs, Bare_point qg);
  bool Query_connection(Bare_point qs, Bare_point qg);

  void Reset_Classification();
  void Set_Alpha(float alpha);
  void Classify_Facets();

  // iterators for looking at cages
  Alpha_iterator Alpha_Value_Begin();
  Alpha_iterator Alpha_Value_End();

  // return the number of sets
  int Num_sets() {
    return disjoint_set_.number_of_sets();
  }

 public:
  float alpha_;
  Alpha_shape_3 alpha_shape_;
  std::list<Cell_handle> interior_cells_;

 private:
  Timer timer_;
  Disjoint_set disjoint_set_;
  std::map<Cell_handle, handle> map_; //map converts cell_handles from location queries to ds structure handles

};


#endif // ALPHA_CON_H
