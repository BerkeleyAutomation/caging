/*
 * =====================================================================================
 *
 *       Filename:  obstaclefiles.h
 *
 *    Description:  This file contains prototypes for obstacles for the PRM planner.
 *
 *        Version:  1.0
 *        Created:  07/22/2011 06:49:08 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Zoe McCarthy (zm), ZoeMcCarthy12@gmail.com
 *        Company:  University of Illinois at Urbana-Champaign
 *
 * =====================================================================================
 */

#ifndef OBSTACLE_FILES_H
#define OBSTACLE_FILES_H

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
#include "PRM.h"

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Polygon_2<Kernel> Polygon;
typedef Polyhedron::HalfedgeDS HalfedgeDS;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Aff_transformation_3<Kernel> CGAL_Aff_Transform;

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Extrude_Polygon
 *  Description:  This function takes a polygon and extrudes it in 3 dimensions so that it
 *                is a polyhedron. There is a class, Build_polyhedron that does the work.
 * =====================================================================================
 */
template <class HDS>
class Build_polyhedron : public CGAL::Modifier_base<HDS> {
public:
    Build_polyhedron(Polygon in, double in_extrude_width) {to_extrude = in; extrude_width_2 = in_extrude_width/2.0;}
    void operator()( HDS& hds) {
        // Postcondition: `hds' is a valid polyhedral surface.
        CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
        int size = to_extrude.size();

        //there are 2*size vertices and 2+size faces
        B.begin_surface( 2*size, 2+size);
        typedef typename HDS::Vertex   Vertex;
        typedef typename Vertex::Point Point;

        //build two sets of vertices, one extruded in negative direction, one in positive
        for(std::vector<Point_2>::iterator iter = to_extrude.vertices_begin(); iter != to_extrude.vertices_end(); iter++)
        {
            B.add_vertex(Point( iter->x(), iter->y(), -extrude_width_2));
        }
        for(std::vector<Point_2>::iterator iter = to_extrude.vertices_begin(); iter != to_extrude.vertices_end(); iter++)
        {
            B.add_vertex(Point( iter->x(), iter->y(), extrude_width_2));
        }

        //get the side faces of the polyhedron.  there are size of them
        for(int i=0; i < size; i++)
        {
            B.begin_facet();
            B.add_vertex_to_facet(i);
            B.add_vertex_to_facet((i + 1)%size);
            B.add_vertex_to_facet((i + 1)%size + size);
            B.add_vertex_to_facet(i + size);
            B.end_facet();
        }
        
        //get the top and bottom polygon faces of the polyhedron. there are two
        B.begin_facet();
        for(int i=size-1; i >= 0; i--)
        {
            B.add_vertex_to_facet(i);
        }
        B.end_facet();

        B.begin_facet();
        for(int i=size; i < 2*size; i++)
        {
            B.add_vertex_to_facet(i);
        }
        B.end_facet();
        B.end_surface();
    }
private:
    Polygon to_extrude;
    double extrude_width_2;
};


//actually build the polyhedron
    Polyhedron
Extrude_Polygon (Polygon in, double extrusion_width)
{
    Polyhedron P;
    Build_polyhedron<HalfedgeDS> polybuilder(in, extrusion_width);
    P.delegate(polybuilder);
    return P;
}		/* -----  end of function Extrude_Polygon  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Rotate_Translate_Plane
 *  Description:  This function produces the CGAL affine transformation for an xy affine
 *                transformation.
 * =====================================================================================
 */
    CGAL_Aff_Transform
Rotate_Translate_Plane (double x, double y, double theta)
{
    double c = cos(theta);
    double s = sin(theta);
    CGAL_Aff_Transform ret_trans(c,-s,0,x,
                                 s,c,0,y,
                                 0,0,1,0);
    return ret_trans;
}		/* -----  end of function Rotate_Translate_Plane  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Create_Triangle
 *  Description:  This function creates a right triangle triangle with scaling factor input
 * =====================================================================================
 */
    Polygon
Create_Triangle (double scale)
{
    Point_2 points[] = {Point_2(0,0),Point_2(scale,0),Point_2(scale,scale)};
    Polygon poly(points, points+3);
    return poly;
}		/* -----  end of function Create_Triangle  ----- */




    /* 
     * ===  FUNCTION  ======================================================================
     *         Name:  triangle
     *  Description:  This class loads a triangle into the obstacle for a PRM
     * =====================================================================================
     */
    void
triangle(PRM& in_PRM, double scale, double x, double y, double theta, double extrude_width)
{
    Polygon pg(Create_Triangle(scale));
    Polyhedron ph(Extrude_Polygon(pg, extrude_width));
    in_PRM.Load_Obstacle(ph,Rotate_Translate_Plane(x,y,theta));
}		/* -----  end of function triangle  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  rectangle
 *  Description:  This class loads a rectangle into the obstacle for a PRM
 * =====================================================================================
 */
    void
rectangle(PRM& in_PRM, double width, double height, double x, double y, double theta, double extrude_width)
{
    Point_2 points[] = {Point_2(0,0),Point_2(width,0),Point_2(width,height), Point_2(0, height)};
    Polygon poly(points, points+4);
    Polyhedron ph(Extrude_Polygon(poly, extrude_width));
    in_PRM.Load_Obstacle(ph,Rotate_Translate_Plane(x,y,theta));
}		/* -----  end of function rectangle  ----- */
#endif //OBSTACLE_FILES_H
