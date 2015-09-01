#ifndef SHAPE_FACT_H
#define SHAPE_FACT_H

#include <list>
#include <set>
#include <vector>
#include <math.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

#include <SOLID/SOLID.h>
#include <SOLID/MT_Scalar.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Segment_3<Kernel> LineSegment;
typedef CGAL::Polygon_2<Kernel> Polygon;
typedef Polyhedron::HalfedgeDS HalfedgeDS;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Triangle_3 Triangle_3;
typedef Polyhedron::Point_iterator Point_iterator;
typedef CGAL::Aff_transformation_3<Kernel> CGAL_Aff_Transform;

struct BoxConfig
{
  float height;
  float width;
  float cx;
  float cy;
  float theta;
};

// equilateral only!
struct TriangleConfig
{
  float scale;
  float cx;
  float cy;
  float theta;
};

struct ComponentConfig
{
  float scale;
  float cx;
  float cy;
  float theta; // in multiples of pi
  std::string sdf_filename;
  std::string obj_filename;
};

struct CompositeObjectConfig
{
  std::vector<BoxConfig> boxes;
  std::vector<TriangleConfig> tris;
};

typedef std::vector<ComponentConfig> MultiObjectConfig;

struct WorkspacePoly
{
  Polyhedron polyhedron;
  CGAL_Aff_Transform tf;  
  Point_2 centroid;
};

// helper class to extrude polygons
template <class HDS>
class PolyhedronBuilder : public CGAL::Modifier_base<HDS> {
 public:
  PolyhedronBuilder(Polygon in, double in_extrude_width);
  void operator()(HDS& hds);

 private:
  Polygon polygon_;
  float extrude_width_2_;
};

// main class to construct polygons
class ShapeFactory
{
 public:
  static WorkspacePoly CreateBox(BoxConfig box, float extrude_width);
  static WorkspacePoly CreateTriangle(TriangleConfig tri, float extrude_width);
  static LineSegment CreateLine(DT_Vector3 p1, DT_Vector3 p2, float extrude_width);
  static CGAL_Aff_Transform CreateRigidTf(float x, float y, float theta);
  static Polyhedron ExtrudePolygon(Polygon in, float extrude_width);
};

// IMPLEMENTATION
template<typename HDS>
PolyhedronBuilder<HDS>::PolyhedronBuilder(Polygon p, double extrude_width)
  : polygon_(p),
    extrude_width_2_(extrude_width / 2)
{
}


template<typename HDS>
void PolyhedronBuilder<HDS>::operator()(HDS& hds)
{
  // Postcondition: `hds' is a valid polyhedral surface.
  CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
  int size = polygon_.size();

  //there are 2*size vertices and 2+size faces
  B.begin_surface( 2*size, 2+size);

  // extremely ugly, artifact of old code
  typedef typename HDS::Vertex   Vertex;
  typedef typename Vertex::Point Point;

  //build two sets of vertices, one extruded in negative direction, one in positive
  for(std::vector<Point_2>::iterator iter = polygon_.vertices_begin(); iter != polygon_.vertices_end(); iter++) {
    B.add_vertex(Point( iter->x(), iter->y(), -extrude_width_2_));
  }
  for(std::vector<Point_2>::iterator iter = polygon_.vertices_begin(); iter != polygon_.vertices_end(); iter++) {
    B.add_vertex(Point( iter->x(), iter->y(), extrude_width_2_));
  }

  //get the side faces of the polyhedron.  there are size of them
  for(int i=0; i < size; i++) {
    B.begin_facet();
    B.add_vertex_to_facet(i);
    B.add_vertex_to_facet((i + 1)%size);
    B.add_vertex_to_facet((i + 1)%size + size);
    B.add_vertex_to_facet(i + size);
    B.end_facet();
  }
        
  //get the top and bottom polygon faces of the polyhedron. there are two
  B.begin_facet();
  for(int i=size-1; i >= 0; i--) {
    B.add_vertex_to_facet(i);
  }
  B.end_facet();

  B.begin_facet();
  for(int i=size; i < 2*size; i++) {
    B.add_vertex_to_facet(i);
  }
  B.end_facet();
  B.end_surface();
}


#endif // SHAPE_FACT_H
