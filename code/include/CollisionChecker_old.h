/*
 * =====================================================================================
 *
 *       Filename:  CollisionChecker.h
 *
 *    Description:  This is a header file for a collision checker class.  It will work 
 *                  with the Robot class to check whether given configurations collide
 *                  and what their separation/penetration distance is
 *
 *        Version:  1.0
 *        Created:  07/14/2011 07:51:14 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Zoe McCarthy (zm), ZoeMcCarthy12@gmail.com
 *        Company:  University of Illinois at Urbana-Champaign
 *
 * =====================================================================================
 */
#ifndef COLLCHECK_H
#define COLLCHECK_H

//#define DEBUG //enable bounding box printouts

#include "WorkspaceObject.h"
#include <list>
#include <set>
#include <vector>
#include <math.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/number_utils.h>
#include <CGAL/utility.h>
#include <SOLID/SOLID.h>
#include <SOLID/MT_Scalar.h>
#include <iostream>
#include <algorithm>
#include <cmath>



#include <stdio.h>
using namespace std;

typedef CGAL::Exact_predicates_exact_constructions_kernel    Kernel;
typedef Kernel::Vector_3                                     Vector;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Polyhedron_3<Kernel>                           Polyhedron;
typedef CGAL::Aff_transformation_3<Kernel>                   CGAL_Aff_Transform;

typedef Polyhedron::Vertex                                   Vertex;
typedef Polyhedron::Vertex_iterator                          Vertex_iterator;
typedef Polyhedron::Halfedge_handle                          Halfedge_handle;
typedef Polyhedron::Edge_iterator                            Edge_iterator;
typedef Polyhedron::Facet_iterator                           Facet_iterator;
typedef Polyhedron::Halfedge_around_vertex_const_circulator  HV_circulator;
typedef Polyhedron::Halfedge_around_facet_circulator         HF_circulator;

typedef CGAL::Aff_transformation_3<Kernel>                   CGAL_Aff_Transform;

/*
 * =====================================================================================
 *        Class:  Object
 *  Description:  This class encapsulates the SOLID collision checking environment for
 *                a 
 * =====================================================================================
 */
class Object
{
    public:
        /* ====================  LIFECYCLE     ======================================= */
	    Object() {}
        Object(DT_ShapeHandle* shape, MT_Scalar margin = 0.0f)
          : m_shape(shape),
            m_object(DT_CreateObject(this, *shape))
        {
            DT_SetMargin(m_object, margin);
        }

        ~Object() //removed delete m_shape due to segfault error on transferring object
        {
            DT_DestroyObject(m_object);
        }

        /* ====================  ACCESSORS     ======================================= */
        DT_ObjectHandle getHandle() const { return m_object; }

        /* ====================  MUTATORS      ======================================= */
        void setShape(DT_ShapeHandle shape)
        {
            *m_shape = shape;
        }

        void Set_Transformation(double* transformation)
        {
            DT_SetMatrixd(m_object, transformation);
        }

        void Set_Response_Class(DT_RespTableHandle respTable, DT_ResponseClass responseClass)
        {
            DT_SetResponseClass(respTable, m_object, responseClass);
        }
        /* ====================  OPERATORS     ======================================= */

        DT_Scalar GetClosestPair(Object* other_obj, DT_Vector3 point1, DT_Vector3 point2)
        {
            return DT_GetClosestPair(m_object, other_obj->m_object, point1, point2);
        }

        DT_ObjectHandle   m_object;
    protected:
        /* ====================  DATA MEMBERS  ======================================= */

    private:
        /* ====================  DATA MEMBERS  ======================================= */
        DT_ShapeHandle*   m_shape;


}; /* -----  end of class Object  ----- */

//global variables for the callback function

//callback function
DT_Bool Collision(void * client_data, void *obj1, void *obj2, const DT_CollData *coll_data);

//helper functions
double Distance_Squared(DT_Vector3 p1, DT_Vector3 p2)
{
    double ret_dist = 0.0;
    for(int i = 0; i < 3; i++)
    {
        ret_dist += (p1[i]-p2[i])*(p1[i]-p2[i]);
    }
    return ret_dist;
}

double Magnitude_Squared(const DT_Vector3 p)
{
    double ret_dist = 0.0;
    for(int i = 0; i < 3; i++)
    {
        ret_dist += p[i]*p[i];
    }
    return ret_dist;
}

double Magnitude_Squared(DT_Vector3 p)
{
    double ret_dist = 0.0;
    for(int i = 0; i < 3; i++)
    {
        ret_dist += p[i]*p[i];
    }
    return ret_dist;
}

/*
 * =====================================================================================
 *        Class:  CollisionChecker
 *  Description:  This class is the structure that allows for collision checking of a sampled
 *                configuration for a given robot and obstacles.
 * =====================================================================================
 */
class CollisionChecker
{
    public:
        /* ====================  LIFECYCLE     ======================================= */
        CollisionChecker ();                             /* constructor */

        /* ====================  ACCESSORS     ======================================= */

        /* ====================  MUTATORS      ======================================= */

        void Load_Robot(WorkspaceObject& robot);
        
        void Load_Obstacle(WorkspaceObject& obstacle, CGAL_Aff_Transform transform);

        void Load_Obstacle(WorkspaceObject& obstacle);
        /* ====================  OPERATORS     ======================================= */
        //check intersection given a list of SOLID affine transformations for the robot's 
        //configuration
        CGAL::Triple<bool, DT_Vector3, DT_Vector3> Check_Intersection(double** m);


        /* ====================  DATA MEMBERS     ======================================= */
        bool               collision_detected;
        DT_Vector3         global_point_1, global_point_2;
        double             distance_squared;
    protected:
        /* ====================  DATA MEMBERS  ======================================= */

    private:
        /* ====================  DATA MEMBERS  ======================================= */
        DT_SceneHandle     g_scene;
        DT_ResponseClass   g_robotClass;
        DT_ResponseClass   g_obstacleClass;
        DT_RespTableHandle g_respTable;

        std::list< std::list< Object* > > robot;
        std::list< Object* >            obstacles;


        /* ====================  HELPER FUNCTIONS  ======================================= */

        DT_ShapeHandle Convert_CGAL_Polyhedron(Polyhedron& input);

        double* Convert_CGAL_DT_Transform(const CGAL_Aff_Transform& transform);

        void Set_Robot_Configuration(double** m);

        void Get_Closest_Pairs();

        int Get_Penetration_Depth();
        

}; /* -----  end of class CollisionChecker  ----- */





/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Convert_CGAL_Polyhedron
 *  Description:  Converts a CGAL polyhedron into a SOLID polyhedron.  Note that in 
 *                SOLID, all polyhedra are convex and thus the input CGAL polyhedron
 *                is assumed to be convex.  This is ensured by performing a convex
 *                decomposition on the workspace objects by CGAL. Note that this function
 *                doesn't currently clean up its own memory and nothing does.
 * =====================================================================================
 */ //TODO: add memory management for verts
    DT_ShapeHandle
CollisionChecker::Convert_CGAL_Polyhedron(Polyhedron& input)
{
    DT_Vector3* verts = new DT_Vector3[input.size_of_vertices()];
    int i = 0;
    Vertex_iterator iter = input.vertices_begin();
    Vertex_iterator iter_end = input.vertices_end();
    for (; iter != iter_end; iter++)
    {
        for(int j = 0; j<3; j++)
        {
            verts[i][j] = (float)( CGAL::to_double( iter->point()[j] ) );
            
#ifdef DEBUG
            //cout << "verts[" << i << "][" << j << "] = " << verts[i][j] << endl;
#endif
        }
        i++;
    }
    DT_VertexBaseHandle base = DT_NewVertexBase(verts[0], sizeof(DT_Vector3));

    DT_ShapeHandle ret_shape = DT_NewComplexShape(base);
    DT_VertexRange(0,input.size_of_vertices());
    DT_EndComplexShape();
    //delete verts; //might need to get rid of this line.  not sure
    return ret_shape;
}		/* -----  end of function CollisionChecker::Convert_CGAL_Polyhedron  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  CollisionChecker
 *  Description:  Initializes all of the tables and classes and responses for a collision
 *                scene.
 * =====================================================================================
 */

 CollisionChecker::CollisionChecker ()
{
    g_scene = DT_CreateScene();
    g_respTable = DT_CreateRespTable();
    g_robotClass = DT_GenResponseClass(g_respTable);
    g_obstacleClass = DT_GenResponseClass(g_respTable);
    collision_detected = false;
    for(int i = 0; i < 3; i++)
    {
        global_point_1[i] = 0.0;
        global_point_2[i] = 100.0;
    }
    distance_squared = 10000.0;
    DT_AddPairResponse(g_respTable, g_robotClass, g_obstacleClass, Collision, DT_DEPTH_RESPONSE, this);
    cout << "CollisionChecker initialized" << endl;


}		/* -----  end of function CollisionChecker  ----- */






        /* 
         * ===  FUNCTION  ======================================================================
         *         Name:  CollisionChecker::Load_Obstacle
         *  Description:  Loads in a workspace object representing an obstacle, and a transform,
         *                and converts it to a SOLID representation, registered in the collision
         *                functions.
         * =====================================================================================
         */
    void
CollisionChecker::Load_Obstacle(WorkspaceObject& obstacle, CGAL_Aff_Transform transform)
{
    DT_ShapeHandle* cur_shape;
    double* SOLID_transform;
    for(std::list<Polyhedron>::iterator iter = obstacle.convex_pieces.begin(); iter != obstacle.convex_pieces.end(); iter++)
    {
        cur_shape = new DT_ShapeHandle(Convert_CGAL_Polyhedron(*iter));
        Object* cur_object = new Object(cur_shape);
        obstacles.push_back(cur_object);
        SOLID_transform = Convert_CGAL_DT_Transform(transform);
        obstacles.back()->Set_Transformation(SOLID_transform);
        delete [] SOLID_transform;
        obstacles.back()->Set_Response_Class(g_respTable, g_obstacleClass);
#ifdef DEBUG
        /*    DT_Vector3 min, max;
        DT_GetBBox(obstacles.back()->m_object, min, max);
        cout << "The min point is: " << endl;
        for(int i = 0; i < 3 ; i++)
        {
            cout << min[i] << endl;
        }
        cout << "The max point is: " << endl;
        for(int i = 0; i < 3 ; i++)
        {
            cout << max[i] << endl;
        }
        */

#endif
    }
    
}		/* -----  end of function CollisionChecker::Load_Obstacle  ----- */


        /* 
         * ===  FUNCTION  ======================================================================
         *         Name:  CollisionChecker::Load_Obstacle
         *  Description:  Loads in a workspace object representing an obstacle
         *                and converts it to a SOLID representation, registered in the collision
         *                functions.
         * =====================================================================================
         */
    void
CollisionChecker::Load_Obstacle(WorkspaceObject& obstacle)
{
    DT_ShapeHandle* cur_shape;
    for(std::list<Polyhedron>::iterator iter = obstacle.convex_pieces.begin(); iter != obstacle.convex_pieces.end(); iter++)
    {
        cur_shape = new DT_ShapeHandle(Convert_CGAL_Polyhedron(*iter));
        Object* cur_object = new Object(cur_shape);
        obstacles.push_back(cur_object);
        obstacles.back()->Set_Response_Class(g_respTable, g_obstacleClass);
    }
    
}		/* -----  end of function CollisionChecker::Load_Obstacle  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  CollisionChecker::Load_Robot
 *  Description:  Loads in a list of workspace objects representing a robot.  Registers
 *                each part of the robot in the collision functions.
 * =====================================================================================
 */
    void
CollisionChecker::Load_Robot (WorkspaceObject& in_robot)
{
    DT_ShapeHandle* cur_shape;
    std::list< Object* > current_robot_section;
    for(std::list<Polyhedron>::iterator iter = in_robot.convex_pieces.begin(); iter != in_robot.convex_pieces.end(); iter++)
    {
        cur_shape = new DT_ShapeHandle(Convert_CGAL_Polyhedron(*iter));
        Object* cur_object = new Object(cur_shape);
        current_robot_section.push_back(cur_object);
        current_robot_section.back()->Set_Response_Class(g_respTable, g_robotClass);
    }
    robot.push_back(current_robot_section);
    
}		/* -----  end of function CollisionChecker::Load_Robot  ----- */



        /* 
         * ===  FUNCTION  ======================================================================
         *         Name:  CollisionChecker::Convert_CGAL_DT_Transform
         *  Description:  Converts the representation of an affine transformation from a CGAL
         *                transformation to a SOLID transformation.  Array needs to be deleted
         *                or a memory leak will occur.
         * =====================================================================================
         */
        double*
CollisionChecker::Convert_CGAL_DT_Transform (const CGAL_Aff_Transform& transform)
{
    //convert from matrix to column-major array
    double* ret_array = new double[16];

    //first 3 rows
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            ret_array[i+4*j] = CGAL::to_double(transform.m(i,j));

        }
    }

    //last row
    for(int j = 0; j< 3; j++)
    {
        ret_array[3+4*j] = 0;
    }
    ret_array[15] = 1;

    return ret_array;
}		/* -----  end of function CollisionChecker::Convert_CGAL_DT_Transform  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  CollisionChecker::Set_Robot_Configuration
 *  Description:  Iterates through the robot's different pieces and sets the right transformation for
 *                each part.
 * =====================================================================================
 */
     void
CollisionChecker::Set_Robot_Configuration (double** m)
{
    std::list< std::list< Object* > >::iterator robot_piece_iterator = robot.begin();
    for(int j = 0; j < 3; j++)
    {
        for(std::list<Object*>::iterator convex_piece_iter = robot_piece_iterator->begin(); convex_piece_iter != robot_piece_iterator->end(); convex_piece_iter++)
        {
            (*convex_piece_iter)->Set_Transformation(m[j]);
#ifdef DEBUG
            /*  DT_Vector3 min, max;
            DT_GetBBox((*convex_piece_iter)->m_object, min, max);
            cout << "the min point is: " << endl;
            for(int i = 0; i < 3 ; i++)
            {
                cout << min[i] << endl;
            }
            cout << "the max point is: " << endl;
            for(int i = 0; i < 3 ; i++)
            {
                cout << max[i] << endl;
            }
             

            for(int i = 0; i < 16; i++)
            {
                cout << i << ": " << m[j][i] << endl;
            }*/

            
#endif
        }
        delete [] m[j];
        robot_piece_iterator++;
    }
    delete [] m;
}		/* -----  end of static function CollisionChecker::Set_Robot_Configuration  ----- */



        /* 
         * ===  FUNCTION  ======================================================================
         *         Name:  CollisionChecker::Check_Intersection
         *  Description:  This is the main API function for the collisionChecker class.  It allows
         *                the PRM class to query whether or not the robot intersects with the obstacles,
         *                and if so, provides a penetration depth.  If not, it provides a separation
         *                distance.
         * =====================================================================================
         */
    CGAL::Triple<bool, DT_Vector3, DT_Vector3>
CollisionChecker::Check_Intersection (double** m)
{
#ifdef DEBUG
/*      for(int j = 0; j < 3; j++)
    {
        cout << "Received in Check_Intersection: " << endl;
        for(int i = 0; i < 16; i++)
        {
            cout << i << ": " << m[j][i] << endl;
        }
    }*/
#endif
    Set_Robot_Configuration(m);

    //set up distances
    collision_detected = false;
    for(int i = 0; i < 3; i++)
    {
        global_point_1[i] = 0.0;
        global_point_2[i] = 0.0;
    }
    distance_squared = 0.0;

    //check for collision
    int num_collisions = Get_Penetration_Depth();


    //set up return
    CGAL::Triple<bool, DT_Vector3, DT_Vector3> triple;
    triple.first = collision_detected;

    //iterate to find separation distance. TODO: fix very ugly O(n^2) implementation
    if(!collision_detected)
    {
        distance_squared = 100000.0;
        Get_Closest_Pairs();
    }
    for(int i = 0; i < 3; i++)
    {
        triple.second[i] = global_point_1[i];
        triple.third[i] = global_point_2[i];
    }

#ifdef DEBUG
    cout << "Number of collisions: " << num_collisions << endl;
    if(num_collisions > 0)
    {
        cout << "Penetration depth squared: " << distance_squared << endl;
    }
    else
    {
        cout << "Separation distance squared: " << distance_squared << endl;
    }
    cout << "Point1: " << global_point_1[0] << " " << global_point_1[1] << " " << global_point_1[2] << endl;
    cout << "Point2: " << global_point_2[0] << " " << global_point_2[1] << " " << global_point_2[2] << endl;
#endif


    return triple;
}		/* -----  end of function CollisionChecker::Check_Intersection  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  CollisionChecker::Collision
 *  Description:  The response callback function for collisions between the robot
 *                and the obstacles. 
 * =====================================================================================
 */
    DT_Bool
Collision (void *client_data,  
			   void *object1,
			   void *object2,
			   const DT_CollData *coll_data)
{
    CollisionChecker* this_checker = (CollisionChecker*)(client_data);
    this_checker->collision_detected = true;
    double new_distance_squared = Magnitude_Squared(coll_data->normal);
    if(new_distance_squared < this_checker->distance_squared)
    {
        for(int i = 0; i < 3; i++)
        {
            this_checker->global_point_1[i] = coll_data->point1[i];
            this_checker->global_point_2[i] = coll_data->point2[i];
        }
        this_checker->distance_squared = new_distance_squared;
    }

	return DT_CONTINUE;
}		/* -----  end of function CollisionChecker::Collision  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  CollisionChecker::Get_Closest_Pairs
 *  Description:  Run after no intersection is detected.  This will return the closest 
 *                pair of points between the obstacles and the robot.
 * =====================================================================================
 */
    void
CollisionChecker::Get_Closest_Pairs ()
{
    DT_Vector3 point_1, point_2;
    double new_distance_squared;
    for( std::list< std::list< Object*> >::iterator robot_piece_iter = robot.begin(); robot_piece_iter != robot.end(); robot_piece_iter++)
    {
        for( std::list<Object*>::iterator convex_piece_iter = robot_piece_iter->begin(); convex_piece_iter != robot_piece_iter->end(); convex_piece_iter++)
        {
            for(std::list<Object*>::iterator obstacle_iter = obstacles.begin(); obstacle_iter != obstacles.end();
                    obstacle_iter++)
            {
                (*convex_piece_iter)->GetClosestPair(*obstacle_iter, point_1, point_2);
                new_distance_squared = Distance_Squared(point_1, point_2);

                //find minimum separation distance
                if(new_distance_squared < distance_squared)
                {
                    for(int i = 0; i < 3; i++)
                    {
                        global_point_1[i] = point_1[i];
                        global_point_2[i] = point_2[i];
                    }
                    distance_squared = new_distance_squared;
                }
            }
        }
    }
}		/* -----  end of function CollisionChecker::Get_Closest_Pairs  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  CollisionChecker::Get_Penetration_Depth
 *  Description:  Run to detect a collision.  This simply checks each pair of obstacle/robot
 *                for collision and then returns the minimum penetration depth.
 * =====================================================================================
 */
    int
CollisionChecker::Get_Penetration_Depth ()
{
    DT_Vector3 point_1, point_2;
    double new_distance_squared;
    bool intersection = false;
    int num_coll = 0;
    for( std::list< std::list< Object*> >::iterator robot_piece_iter = robot.begin(); robot_piece_iter != robot.end(); robot_piece_iter++)
    {
        for( std::list<Object*>::iterator convex_piece_iter = robot_piece_iter->begin(); convex_piece_iter != robot_piece_iter->end(); convex_piece_iter++)
        {
            for(std::list<Object*>::iterator obstacle_iter = obstacles.begin(); obstacle_iter != obstacles.end();
                    obstacle_iter++)
            {


                intersection = DT_GetPenDepth((*convex_piece_iter)->m_object, (*obstacle_iter)->m_object, point_1, point_2);
                if(intersection)
                {
                    num_coll++;
                    collision_detected = true;
                    new_distance_squared = Distance_Squared(point_1, point_2);

                    //find maximum penetration depth
                    if(new_distance_squared > distance_squared)
                    {
                        for(int i = 0; i < 3; i++)
                        {
                            global_point_1[i] = point_1[i];
                            global_point_2[i] = point_2[i];
                        }
                        distance_squared = new_distance_squared;
                    }
                }
            }
        }
    }
    return num_coll;
}		/* -----  end of function CollisionChecker::Get_Penetration_Depth  ----- */

#endif // COLLCHECK_H

