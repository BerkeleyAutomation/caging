/*********************
Zoe McCarthy, May 2011

//defines  Alpha_shapes_disconnection class.
//supports an alpha shape for weighted points (spheres) in R^3
//has an auxillary structure for querying two points and determining
//if they are in the same connected component of the free space
********************/

#ifndef ALPHA_H
#define ALPHA_H

////// ********************    PREAMBLE    ********************/

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Union_find.h>
#include <CGAL/Timer.h>

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



// function to generate some kind of spheres of spheres
std::list<Weighted_point> generate_sphere_of_spheres(double x, double y, double z, double large_radius,
                                                     double small_radius, int stacks, int slices);

// function to make random spheres
std::list<Weighted_point> generate_random_spheres(double x_min, double y_min, double z_min, 
                                                  double x_max, double y_max, double z_max, 
                                                  double rad_min, double rad_max, int number);

struct DifferentialTriangle
{
  Bare_point center;
  Dir_3 dir;
  Tri_3 tri;
  float area;
};

////// ********************    CLASS DEFINITION    ********************/

// this class will allow querying of two points to test 
// whether or not they are disconnected
// its underlying structure is a disjoint set of 
// connected components of cells that are not in the alpha shape
class Alpha_shapes_disconnection
{
 public:
  // construct the disjoint set that contains the connectivity information
  Alpha_shapes_disconnection(Triangulation_3 raw_triangulation, float alpha = 0.0f);

  Cell_vector Find_Path_Between(Bare_point qs, Bare_point qg);
  bool Query_disconnection(Bare_point qs, Bare_point qg);

  // get the volume of the connected component containing the point q
  bool Component_Volume(Bare_point q, float& volume);

  // discretize the boundary of the free region of the zero-shape containing q and the area outside of the "alpha"-shape 
  // convex hull alpha should be related to the size of the bounded region without infinity points
  bool Escape_Boundary(Bare_point q, float conv_alpha, float max_area, std::vector<DifferentialTriangle>& separator_tris,
                       float& conv_hull_area, float& conv_hull_volume, CGAL::Geomview_stream& gv);

  bool Add_Tri_Centers(Tri_3 tri, float max_area, std::vector<DifferentialTriangle>& tri_vec);


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
  std::list<Cell_handle> external_cells_;

 private:
  Timer timer_;
  Disjoint_set disjoint_set_;
  std::map<Cell_handle, handle> map_; //map converts cell_handles from location queries to ds structure handles

};


#endif // ALPHA_H
