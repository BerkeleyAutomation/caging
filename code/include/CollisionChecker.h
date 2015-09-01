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
#include "Util.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <list>
#include <math.h>
#include <set>
#include <vector>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/number_utils.h>
#include <CGAL/utility.h>
#include <SOLID/SOLID.h>
#include <SOLID/MT_Scalar.h>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

#include <stdio.h>
using namespace std;

typedef CGAL::Exact_predicates_exact_constructions_kernel    Kernel;
typedef Kernel::Vector_3                                     Vector;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Aff_transformation_3<Kernel>                   CGAL_Aff_Transform;

typedef Polyhedron_3::Vertex                                   Vertex;
typedef Polyhedron_3::Vertex_iterator                          Vertex_iterator;
typedef Polyhedron_3::Halfedge_handle                          Halfedge_handle;
typedef Polyhedron_3::Edge_iterator                            Edge_iterator;
typedef Polyhedron_3::Facet_iterator                           Facet_iterator;
typedef Polyhedron_3::Halfedge_around_vertex_const_circulator  HV_circulator;
typedef Polyhedron_3::Halfedge_around_facet_circulator         HF_circulator;

typedef CGAL::Aff_transformation_3<Kernel>                   CGAL_Aff_Transform;

/*
 * =====================================================================================
 *        Class:  CollisionObject
 *  Description:  This class encapsulates the SOLID collision checking environment for objects
 * =====================================================================================
 */
class CollisionObject
{
 public:
  /* ====================  LIFECYCLE     ======================================= */
  CollisionObject() {}
  CollisionObject(DT_ShapeHandle* shape, MT_Scalar margin = 0.0f);
  ~CollisionObject();
  /* ====================  ACCESSORS     ======================================= */
  DT_ObjectHandle getHandle();

  void setShape(DT_ShapeHandle shape);
  void Set_Transformation(Eigen::Matrix4f transform);
  void Set_Response_Class(DT_RespTableHandle respTable, DT_ResponseClass responseClass);
  DT_Scalar GetClosestPair(CollisionObject* other_obj, DT_Vector3 point1, DT_Vector3 point2);

 private:
  DT_ShapeHandle* shape_;
  DT_ObjectHandle object_;
};

// utility functions
// callback function
DT_Bool Collision(void * client_data, void *obj1, void *obj2, const DT_CollData *coll_data);

/*
 * =====================================================================================
 *        Class:  CollisionChecker
 *  Description:  This class is the structure that allows for collision checking of a sampled
 *                pose for an object given a gripper configuration
 * =====================================================================================
 */
class CollisionChecker
{
 public:
  /* ====================  LIFECYCLE     ======================================= */
  CollisionChecker ();                             /* constructor */

  // load an object or gripper into the collision checker
  void Load_Object(WorkspaceObject& object);
  void Load_Object(WorkspaceObject& object, CGAL_Aff_Transform transform);
  void Load_Gripper_Finger(WorkspaceObject& gripper_finger);
  void Load_Gripper_Finger(WorkspaceObject& gripper_finger, CGAL_Aff_Transform transform);

  // check intersection for a list of object poses
  // TODO: probably change to a different way of specifying poses
  CGAL::Triple<bool, DT_Vector3, DT_Vector3> Check_Intersection(Eigen::Matrix4f pose);

  // Set the object pose (typically only called for collision checking, but also can be set for vis)
  void Set_Object_Pose(Eigen::Matrix4f pose);

 private:
  // conversions between CGAL, SOLID
  DT_ShapeHandle Convert_CGAL_Polyhedron(Polyhedron_3& input);
  Eigen::Matrix4f Convert_CGAL_DT_Transform(const CGAL_Aff_Transform& transform);

  // depth calcs
  void Get_Closest_Pairs();
  int Get_Penetration_Depth();
       
 public:
  // collision vars
  bool               collision_detected_;
  DT_Vector3         global_point_1_, global_point_2_;
  double             distance_sq_;

 private:
  // solid handlers
  DT_SceneHandle     scene_;
  DT_ResponseClass   object_class_;
  DT_ResponseClass   gripper_class_;
  DT_RespTableHandle resp_table_;

  // collision objects to handle object and gripper 
  std::vector< CollisionObject* > objects_;
  std::vector< CollisionObject* > gripper_fingers_;
};

#endif // COLLCHECK_H

