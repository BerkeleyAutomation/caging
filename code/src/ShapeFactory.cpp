#include "ShapeFactory.hpp"

WorkspacePoly ShapeFactory::CreateBox(BoxConfig box, float extrude_width)
{
  float width_2 = box.width / 2;
  float height_2 = box.height / 2;

  Point_2 points[] = { Point_2( width_2, -height_2),
                       Point_2(-width_2, -height_2),
                       Point_2(-width_2,  height_2),
                       Point_2( width_2,  height_2)};
  Polygon poly(points, points+4);
  Polyhedron ph(ExtrudePolygon(poly, extrude_width));

  CGAL_Aff_Transform tf = CreateRigidTf(box.cx, box.cy, box.theta);

  WorkspacePoly pg;
  pg.polyhedron = ph;
  pg.tf = tf;
  pg.centroid = Point_2(0.0f, 0.0f);
  return pg;
}

WorkspacePoly
ShapeFactory::CreateTriangle(TriangleConfig tri, float extrude_width)
{
  float scale_2 = tri.scale / 2; 
  float height = scale_2 * sqrt(3);
  float height_2 = height / 2;
  Point_2 points[] = { Point_2(       0,  height_2),
                       Point_2(-scale_2, -height_2),
                       Point_2( scale_2, -height_2)};
  Polygon poly(points, points+3);
  Polyhedron ph(ExtrudePolygon(poly, extrude_width));
  CGAL_Aff_Transform tf = CreateRigidTf(tri.cx, tri.cy, tri.theta);

  WorkspacePoly pg;
  pg.polyhedron = ph;
  pg.tf = tf;
  pg.centroid = Point_2(0.0f, 0.0f);
  return pg;
}		/* -----  end of function triangle  ----- */

LineSegment ShapeFactory::CreateLine(DT_Vector3 p1, DT_Vector3 p2, float extrude_width)
{
  Point_3 point1 = Point_3(p1[0], p1[1], extrude_width / 2);
  Point_3 point2 = Point_3(p2[0], p2[1], extrude_width / 2);
  LineSegment ls(point1, point2);
  return ls;
}

CGAL_Aff_Transform
ShapeFactory::CreateRigidTf(float x, float y, float theta)
{
  double c = cos(theta);
  double s = sin(theta);
  CGAL_Aff_Transform ret_trans(c,-s,0,x,
                               s,c,0,y,
                               0,0,1,0);
  return ret_trans;
}

//actually build the polyhedron
Polyhedron
ShapeFactory::ExtrudePolygon(Polygon in, float extrude_width)
{
  Polyhedron ph;
  PolyhedronBuilder<HalfedgeDS> polybuilder(in, extrude_width);
  ph.delegate(polybuilder);
  return ph;
}

