#ifndef TYPEDEF_H
#define TYPEDEF_H

#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/IO/Polyhedron_geomview_ostream.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Triangulation_geomview_ostream_2.h>
#include <CGAL/IO/Triangulation_geomview_ostream_3.h>
#include <CGAL/Nef_3/SNC_indexed_items.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Timer.h>
#include <CGAL/Triangulation_euclidean_traits_xy_3.h>
#include <CGAL/Union_find.h>
#include <CGAL/convex_decomposition_3.h> 
#include <CGAL/intersections.h>
#include <CGAL/number_utils.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/utility.h>
#include <algorithm>
#include <cmath>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <math.h>
#include <opencv2/opencv.hpp>
#include <queue>
#include <set>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <phat/boundary_matrix.h>

using namespace std;

enum Alpha_Classification {DEFAULT, EMPTY, FULL, MIXED};
enum Alpha_status { EXTERIOR, INTERIOR };
enum PotentialType {ENERGY_GRAVITY};

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Regular_triangulation_euclidean_traits_3<K> Gt;
typedef CGAL::Alpha_shape_cell_base_3<Gt>           Fb;
typedef CGAL::Alpha_shape_vertex_base_3<Gt>         Vb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb> Tds;
typedef CGAL::Regular_triangulation_3<Gt,Tds>       Triangulation_3;
typedef CGAL::Triangulation_3<Gt, Tds>               Raw_triangulation_3;              
typedef CGAL::Alpha_shape_3<Triangulation_3>        Alpha_shape_3;
typedef Alpha_shape_3::Alpha_iterator       Alpha_iterator;
typedef Alpha_shape_3::Cell_handle          Cell_handle;
typedef Alpha_shape_3::Edge                 Edge;
typedef Alpha_shape_3::Facet                Facet;
typedef Alpha_shape_3::Vertex_handle        Vertex_handle;
typedef Cell_handle                         T;
typedef std::vector<Cell_handle>                            Cell_vector;
typedef std::vector<double>                                 Configuration;
typedef std::list< Cell_handle >                            Cell_path;
typedef std::list< Configuration >                          Path;
typedef std::pair< Cell_path, Alpha_Classification >        Cell_path_with_exist;
typedef std::pair< Path, bool >                             Path_with_exist;
//typedef std::vector<Cell_handle>               Raw_triangulation_3;              

typedef CGAL::Aff_transformation_3<Kernel> CGAL_Aff_Transform;
typedef CGAL::Direction_3<K>                        Dir_3;
typedef CGAL::Nef_polyhedron_3<Kernel, CGAL::SNC_indexed_items> Nef_polyhedron_3;
typedef CGAL::Point_3<K>                            Pt_3;
typedef CGAL::Polygon_2<Kernel> Polygon;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
typedef CGAL::Tetrahedron_3<K>                      Tetrahedron_3;
typedef CGAL::Timer                                         Timer;
typedef CGAL::Triangle_3<K>                         Tri_3;
typedef CGAL::Triangulation_3<K>       Simple_Triangulation_3;
typedef CGAL::Union_find<T>                 Disjoint_set;  
typedef CGAL::Unique_hash_map<Cell_handle, Alpha_Classification>  Cell_class_map;
typedef CGAL::Unique_hash_map<Cell_handle, Alpha_Classification>  Cell_map;
typedef CGAL::Unique_hash_map<Cell_handle, double>         Alpha_cell_map;
typedef CGAL::Unique_hash_map<Facet, double>               Alpha_facet_map;
typedef CGAL::Unique_hash_map<Vertex_handle, bool>         Point_map;

typedef Disjoint_set::handle                handle;
typedef Gt::Bare_point                      Bare_point;
typedef Gt::Weighted_point                  Weighted_point;
typedef Kernel::Point_2 			            Point_2;
typedef Kernel::Point_3 			            Point_3;
typedef Kernel::Point_3                                     Point;
typedef Kernel::Tetrahedron_3 				    Tetrahedron;
typedef Kernel::Triangle_3 				    Triangle;
typedef Kernel::Vector_3                                     Vector;
typedef Nef_polyhedron_3::Volume_const_iterator Volume_const_iterator;
typedef Polyhedron::HalfedgeDS 				    HalfedgeDS;
typedef Polyhedron::Halfedge_handle        		    Halfedge_handle;
typedef Polyhedron_3::Edge_iterator                            Edge_iterator;
typedef Polyhedron_3::Facet_iterator                           Facet_iterator;
typedef Polyhedron_3::Halfedge_around_facet_circulator         HF_circulator;
typedef Polyhedron_3::Halfedge_around_vertex_const_circulator  HV_circulator;
typedef Polyhedron_3::Halfedge_handle                          Halfedge_handle;
typedef Polyhedron_3::Vertex                                   Vertex;
typedef Polyhedron_3::Vertex_iterator                          Vertex_iterator;
typedef Triangulation_3::Cell_circulator 		    Cell_circulator;
typedef Triangulation_3::Cell_handle                        Cell_handle;
typedef Triangulation_3::Edge                               Edge;
typedef Triangulation_3::Facet                              Facet;
typedef Triangulation_3::Locate_type    		    Locate_type;
typedef Triangulation_3::Vertex_handle                      Vertex_handle;
typedef double 						    Potential;
typedef phat::index  					    Index;

#endif // TYPEDEF_H
