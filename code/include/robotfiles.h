/*
 * =====================================================================================
 *
 *       Filename:  robotfiles.h
 *
 *    Description:  This file contains functions for loading robots and obstacles into
 *                  the PRM planner.
 *
 *        Version:  1.0
 *        Created:  07/22/2011 01:16:46 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Zoe McCarthy (zm), ZoeMcCarthy12@gmail.com
 *        Company:  University of Illinois at Urbana-Champaign
 *
 * =====================================================================================
 */


#ifndef ROBOT_FILES_H
#define ROBOT_FILES_H

#define DEBUG

#ifdef DEBUG
#include <stdio.h>
using namespace std;
#endif

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Aff_transformation_3.h>
#include <iostream>
#include "PRM.h"
#include "ConfigurationMapper.h"
#include <vector>
#include <list>
#include <math.h>


typedef std::vector<double>                Configuration;

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Aff_transformation_3<Kernel> CGAL_Aff_Transform;
typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;
typedef Polyhedron::Halfedge_handle        Halfedge_handle;

template <class Poly>
typename Poly::Halfedge_handle make_rectangle_3( Poly& P, double width_div_2, double joint_length, double height_div_2) {
    // appends a cube of size [0,1]^3 to the polyhedron P.
    CGAL_precondition( P.is_valid());
    typedef typename Poly::Point_3         Point;
    typedef typename Poly::Halfedge_handle Halfedge_handle;
    Halfedge_handle h = P.make_tetrahedron( Point( joint_length, -width_div_2, -height_div_2),
                                            Point( 0, -width_div_2, height_div_2),
                                            Point( 0, -width_div_2, -height_div_2),
                                            Point( 0, width_div_2, -height_div_2));
    Halfedge_handle g = h->next()->opposite()->next();             // Fig. (a)
    P.split_edge( h->next());
    P.split_edge( g->next());
    P.split_edge( g);                                              // Fig. (b)
    h->next()->vertex()->point()     = Point( joint_length, -width_div_2, height_div_2);
    g->next()->vertex()->point()     = Point( 0, width_div_2, height_div_2);
    g->opposite()->vertex()->point() = Point( joint_length, width_div_2, -height_div_2);            // Fig. (c)
    Halfedge_handle f = P.split_facet( g->next(),
                                       g->next()->next()->next()); // Fig. (d)
    Halfedge_handle e = P.split_edge( f);
    e->vertex()->point() = Point( joint_length, width_div_2, height_div_2);                        // Fig. (e)
    P.split_facet( e, f->next()->next());                          // Fig. (f)
    CGAL_postcondition( P.is_valid());
    return h;
}


/*
 * =====================================================================================
 *        Class:  LinkMapper
 *  Description:  This class is the ConfigurationMapper class for the 3 link robot.
 * =====================================================================================
 */
class LinkMapper: public ConfigurationMapper
{
    public:

        /* ====================  LIFECYCLE     ======================================= */
        LinkMapper (double in_width, double height, double in_joint_length0, double in_joint_length1, double in_joint_length2);                             /* constructor      */

        /* ====================  ACCESSORS     ======================================= */

        /* ====================  MUTATORS      ======================================= */

        /* ====================  OPERATORS     ======================================= */

        std::list<CGAL_Aff_Transform> Display_Transformations(Configuration& configuration);

        double** Collision_Transformations(Configuration& configuration);

    protected:
        /* ====================  DATA MEMBERS  ======================================= */

    private:
        /* ====================  DATA MEMBERS  ======================================= */
        double width, height;
        double joint_length0, joint_length1, joint_length2;

}; /* -----  end of class LinkMapper  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  revolute_3_link
 *  Description:  This function takes in a PRM planner and places a 3 link revolute
 *                joint robot into its robot.  It also provides the necessary ConfigurationMapper
 *                for the planner to function.
 * =====================================================================================
 */
    bool
revolute_3_link (PRM& in_PRM, double width, double height, double joint_length1, double joint_length2, double joint_length3, double angle_max,
        double extra_factor)
{
    double width_div_2 = width/2.0;
    double height_div_2 = height/2.0;

    //collision and geometry
    Polyhedron P1;
    Halfedge_handle h1 = make_rectangle_3(P1, width_div_2, joint_length1, height_div_2);
    in_PRM.Load_Robot(P1);
    Polyhedron P2;
    Halfedge_handle h2 = make_rectangle_3(P2, width_div_2, joint_length2, height_div_2);
    in_PRM.Load_Robot(P2);
    Polyhedron P3;
    Halfedge_handle h3 = make_rectangle_3(P3, width_div_2, joint_length3, height_div_2);
    in_PRM.Load_Robot(P3);

    //configuration mapper
    LinkMapper* link_map = new LinkMapper(width, height, joint_length1, joint_length2, joint_length3);
    in_PRM.Load_Configuration_Mapper(link_map);

    //planner
    double x_scale, y_scale, z_scale;
    x_scale = sqrt(pow(width_div_2,2) + pow(joint_length1 + joint_length2 + joint_length3,2));
    std::cout << "x_scale : " << x_scale << std::endl;
    y_scale = sqrt(pow(width_div_2,2) + pow(joint_length2 + joint_length3,2));
    std::cout << "y_scale : " << y_scale << std::endl;
    z_scale = sqrt(pow(width_div_2,2) + pow(joint_length3,2));
    std::cout << "z_scale : " << z_scale << std::endl;
    in_PRM.Set_Triangulation_Planner(x_scale, y_scale, z_scale, -angle_max, angle_max,
            -angle_max, angle_max, -angle_max, angle_max, extra_factor);
    return true;
}		/* -----  end of function revolute_3_link  ----- */






/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  LinkMapper::3LinkMapper
 *  Description:  This function implements the constructor for the 3 link mapper
 * =====================================================================================
 */
    
LinkMapper::LinkMapper (double in_width, double in_height, double in_joint_length0, double in_joint_length1, double in_joint_length2)
{
    width = in_width;
    height = in_height;
    joint_length0 = in_joint_length0;
    joint_length1 = in_joint_length1;
    joint_length2 = in_joint_length2;
}		/* -----  end of function LinkMapper::LinkMapper  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  LinkMapper::Display_Transformations
 *  Description:  This method implements the mapping to display transformations for a 
 *                3 link robot arm. The robot is constrained to the x-y plane currently
 *                and all rotations are in this plane.
 * =====================================================================================
 */
std::list<CGAL_Aff_Transform>
LinkMapper::Display_Transformations (Configuration& configuration)
{
    double c0, s0, c1, s1, c2, s2;
    c0 = cos(configuration[0]);
    s0 = sin(configuration[0]);
    c1 = cos(configuration[1]);
    s1 = sin(configuration[1]);
    c2 = cos(configuration[2]);
    s2 = sin(configuration[2]);
    CGAL_Aff_Transform transform_1(c0,-s0,0,0,
                                    s0,c0,0,0,
                                    0,0,1,0,1);
    std::list<CGAL_Aff_Transform> ret_list;
    ret_list.push_back(transform_1);
    CGAL_Aff_Transform transform_2(c1,-s1,0,joint_length0,
                                    s1,c1,0,0,
                                    0,0,1,0,1);

    ret_list.push_back(transform_1*transform_2);
    CGAL_Aff_Transform transform_3(c2,-s2,0,joint_length1,
                                    s2,c2,0,0,
                                    0,0,1,0,1);

    ret_list.push_back(transform_1*transform_2*transform_3);
    

    return ret_list;
}		/* -----  end of function LinkMapper::Display_Transformations  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  LinkMapper::Collision_Transformations
 *  Description:  This method implements the mapping to collision transformations for a 
 *                3 link robot arm. The robot is constrained to the x-y plane currently
 *                and all rotations are in this plane.
 * =====================================================================================
 */
double**
LinkMapper::Collision_Transformations (Configuration& configuration)
{
    double c0, s0, c1, s1, c2, s2;
    c0 = cos(configuration[0]);
    s0 = sin(configuration[0]);
    c1 = cos(configuration[1]);
    s1 = sin(configuration[1]);
    c2 = cos(configuration[2]);
    s2 = sin(configuration[2]);

    // do common math operations for speed
    double c01 = c0*c1 - s0*s1;
    double s01 = c0*s1 + c1*s0;
    double c012 = c01*c2 - s01*s2;
    double s012 = c01*s2 + c2*s01;
    double c0j0 = c0*joint_length0;
    double s0j0 = s0*joint_length0;

    // write down affine transformation matrices without doing matrix multiplications
    // matrices stored in column-major format as per specification from SOLID
    //double T0[] = { c0, -s0, 0, 0,  s0, c0, 0, 0,  0, 0, 1, 0,  0, 0, 0, 1 };
    //double T1[] = { c01, -s01, 0, c0j0,  s01, c01, 0, s0j0,  0, 0, 1, 0, 0, 0, 0, 1 };
    //double T2[] = { c012, -s012, 0, c0j0 + c01*joint_length1,  s012, c012, 0, s0j0 + s01*joint_length1,  0, 0, 1, 0,  0, 0, 0, 1 };
    
    double** ret_list = new double*[3];
    ret_list[0] = new double[16]; 
    ret_list[1] = new double[16];
    ret_list[2] = new double[16];

    ret_list[0][0] = c0;
    ret_list[0][1] = s0;
    ret_list[0][2] = 0;
    ret_list[0][3] = 0;
    ret_list[0][4] = -s0;
    ret_list[0][5] = c0;
    ret_list[0][6] = 0;
    ret_list[0][7] = 0;
    ret_list[0][8] = 0;
    ret_list[0][9] = 0;
    ret_list[0][10] = 1;
    ret_list[0][11] = 0;
    ret_list[0][12] = 0;
    ret_list[0][13] = 0;
    ret_list[0][14] = 0;
    ret_list[0][15] = 1;
    
    ret_list[1][0] = c01;
    ret_list[1][1] = s01;
    ret_list[1][2] = 0;
    ret_list[1][3] = 0;
    ret_list[1][4] = -s01;
    ret_list[1][5] = c01;
    ret_list[1][6] = 0;
    ret_list[1][7] = 0;
    ret_list[1][8] = 0;
    ret_list[1][9] = 0;
    ret_list[1][10] = 1;
    ret_list[1][11] = 0;
    ret_list[1][12] = c0j0;
    ret_list[1][13] = s0j0;
    ret_list[1][14] = 0;
    ret_list[1][15] = 1;

    ret_list[2][0] = c012;
    ret_list[2][1] = s012;
    ret_list[2][2] = 0;
    ret_list[2][3] = 0;
    ret_list[2][4] = -s012;
    ret_list[2][5] = c012;
    ret_list[2][6] = 0;
    ret_list[2][7] = 0;
    ret_list[2][8] = 0;
    ret_list[2][9] = 0;
    ret_list[2][10] = 1;
    ret_list[2][11] = 0;
    ret_list[2][12] = c0j0 + c01*joint_length1;
    ret_list[2][13] = s0j0 + s01*joint_length1;
    ret_list[2][14] = 0;
    ret_list[2][15] = 1;

#ifdef DEBUG
/*      cout << "c0: " << c0 << endl << "s0: " << s0 << endl;
    for(int i = 0; i < 16; i++)
    {
        cout << "T0[" << i << "]: " << ret_list[0][i] << endl;
    }
    for(int i = 0; i < 16; i++)
    {
        cout << "T1[" << i << "]: " << ret_list[1][i] << endl;
    }
    for(int i = 0; i < 16; i++)
    {
        cout << "T2[" << i << "]: " << ret_list[2][i] << endl;
    }*/
#endif


    return ret_list;
}		/* -----  end of function LinkMapper::Display_Transformations  ----- */

#endif //ROBOT_FILES_H
