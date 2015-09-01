#include "CollisionChecker.h"

// CollisionObect: functions to encapculate tranforms on SOLID collision objects
CollisionObject::CollisionObject(DT_ShapeHandle* shape, MT_Scalar margin)
  : shape_(shape),
    object_(DT_CreateObject(this, *shape))
{
  DT_SetMargin(object_, margin);
}

CollisionObject::~CollisionObject() {
  DT_DestroyObject(object_);
}

// CollisionObject Accessors and Setters
DT_ObjectHandle CollisionObject::getHandle()
{ 
  return object_;
}

void CollisionObject::setShape(DT_ShapeHandle shape)
{
  *shape_ = shape;
}

void CollisionObject::Set_Transformation(Eigen::Matrix4f transformation)
{
  DT_SetMatrixf(object_, (float*)transformation.data());
}

void CollisionObject::Set_Response_Class(DT_RespTableHandle respTable, DT_ResponseClass responseClass) {
  DT_SetResponseClass(respTable, object_, responseClass);
}

DT_Scalar CollisionObject::GetClosestPair(CollisionObject* other_obj, DT_Vector3 point1, DT_Vector3 point2) {
  return DT_GetClosestPair(object_, other_obj->object_, point1, point2);
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Collision
 *  Description:  The response callback function for collisions between the robot
 *                and the obstacles. 
 * =====================================================================================
 */
DT_Bool Collision(void *client_data, void *object1, void *object2, const DT_CollData *coll_data)
{
  CollisionChecker* this_checker = (CollisionChecker*)(client_data);
  this_checker->collision_detected_ = true;
  double new_distance_squared = Magnitude_Squared(coll_data->normal);
  if (new_distance_squared < this_checker->distance_sq_) {
    for (int i = 0; i < 3; i++) {
      this_checker->global_point_1_[i] = coll_data->point1[i];
      this_checker->global_point_2_[i] = coll_data->point2[i];
    }
    this_checker->distance_sq_ = new_distance_squared;
  }
  return DT_CONTINUE;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  CollisionChecker
 *  Description:  Initializes all of the tables and classes and responses for a collision
 *                scene.
 * =====================================================================================
 */

CollisionChecker::CollisionChecker ()
{
  scene_ = DT_CreateScene();
  resp_table_ = DT_CreateRespTable();
  object_class_ = DT_GenResponseClass(resp_table_);
  gripper_class_ = DT_GenResponseClass(resp_table_);
  collision_detected_ = false;
  for(int i = 0; i < 3; i++) {
    global_point_1_[i] = 0.0;
    global_point_2_[i] = 100.0;
  }
  distance_sq_ = 10000.0;
  DT_AddPairResponse(resp_table_, object_class_, gripper_class_, Collision, DT_DEPTH_RESPONSE, this);
  std::cout << "CollisionChecker initialized" << std::endl;
}

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
CollisionChecker::Convert_CGAL_Polyhedron(Polyhedron_3& input)
{
  DT_Vector3* verts = new DT_Vector3[input.size_of_vertices()];
  Vertex_iterator iter = input.vertices_begin();
  Vertex_iterator iter_end = input.vertices_end();
  int i = 0;
  for (; iter != iter_end; iter++, i++) {
    for(int j = 0; j < 3; j++) {
      verts[i][j] = (float)( CGAL::to_double( iter->point()[j] ) );
    }
  }
  DT_VertexBaseHandle base = DT_NewVertexBase(verts[0], sizeof(DT_Vector3));
  DT_ShapeHandle ret_shape = DT_NewComplexShape(base);
  DT_VertexRange(0, input.size_of_vertices());
  DT_EndComplexShape();
  //delete verts; //might need to get rid of this line.  not sure
  return ret_shape;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  CollisionChecker::Load_Object
 *  Description:  Loads in a workspace object representing an object to grasp, and a transform,
 *                and converts it to a SOLID representation, registered in the collision
 *                functions.
 * =====================================================================================
 */
void
CollisionChecker::Load_Object(WorkspaceObject& object)
{
  DT_ShapeHandle* cur_shape;
  for(std::list<Polyhedron_3>::iterator iter = object.convex_pieces.begin(); iter != object.convex_pieces.end(); iter++) {
    cur_shape = new DT_ShapeHandle(Convert_CGAL_Polyhedron(*iter));
    CollisionObject* cur_object = new CollisionObject(cur_shape);
    objects_.push_back(cur_object);
    objects_.back()->Set_Response_Class(resp_table_, object_class_);
  } 
}

void
CollisionChecker::Load_Object(WorkspaceObject& object, CGAL_Aff_Transform transform)
{
  DT_ShapeHandle* cur_shape;
  Eigen::Matrix4f pose = Convert_CGAL_DT_Transform(transform);
  for(std::list<Polyhedron_3>::iterator iter = object.convex_pieces.begin(); iter != object.convex_pieces.end(); iter++) {
    cur_shape = new DT_ShapeHandle(Convert_CGAL_Polyhedron(*iter));
    CollisionObject* cur_object = new CollisionObject(cur_shape);
    objects_.push_back(cur_object);
    objects_.back()->Set_Transformation(pose);
    objects_.back()->Set_Response_Class(resp_table_, object_class_);
  }  
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  CollisionChecker::Load_Gripper
 *  Description:  Loads in a list of workspace objects representing a robot gripper.  Registers
 *                each part of the gripper in the collision functions.
 * =====================================================================================
 */
void
CollisionChecker::Load_Gripper_Finger(WorkspaceObject& gripper_finger)
{
  DT_ShapeHandle* cur_shape;
  std::list<Polyhedron_3>::iterator iter;
  for(iter = gripper_finger.convex_pieces.begin();
      iter != gripper_finger.convex_pieces.end(); iter++) {
    cur_shape = new DT_ShapeHandle(Convert_CGAL_Polyhedron(*iter));
    CollisionObject* cur_object = new CollisionObject(cur_shape);
    gripper_fingers_.push_back(cur_object); 
    gripper_fingers_.back()->Set_Response_Class(resp_table_, gripper_class_);
  }
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  CollisionChecker::Load_Gripper
 *  Description:  Loads in a list of workspace objects representing a robot gripper.  Registers
 *                each part of the gripper in the collision functions.
 * =====================================================================================
 */
void
CollisionChecker::Load_Gripper_Finger(WorkspaceObject& gripper_finger, CGAL_Aff_Transform transform)
{
  DT_ShapeHandle* cur_shape;
  Eigen::Matrix4f pose = Convert_CGAL_DT_Transform(transform);
  std::list<Polyhedron_3>::iterator iter;
  for(iter = gripper_finger.convex_pieces.begin();
      iter != gripper_finger.convex_pieces.end(); iter++) {
    cur_shape = new DT_ShapeHandle(Convert_CGAL_Polyhedron(*iter));
    CollisionObject* cur_object = new CollisionObject(cur_shape);
    gripper_fingers_.push_back(cur_object); 
    gripper_fingers_.back()->Set_Transformation(pose);
    gripper_fingers_.back()->Set_Response_Class(resp_table_, gripper_class_);
  }
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  CollisionChecker::Convert_CGAL_DT_Transform
 *  Description:  Converts the representation of an affine transformation from a CGAL
 *                transformation to an Eigen Matrix.
 * =====================================================================================
 */
Eigen::Matrix4f
CollisionChecker::Convert_CGAL_DT_Transform(const CGAL_Aff_Transform& transform)
{
  //convert from matrix to column-major array
  Eigen::Matrix4f pose;

  //first 3 rows
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 4; j++) {
      pose(i, j) = (float)CGAL::to_double(transform.m(i,j));  
    }
  }

  //last row
  for(int j = 0; j < 3; j++) {
    pose(3, j) = 0;
  }
  pose(3,3) = 1;
  return pose;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  CollisionChecker::Set_Object_Poses
 *  Description:  Iterates through the object pieces and sets each to the corresponding pose in
 *                the corresponding N x 16 pose array
 * =====================================================================================
 */
// Sets all pieces of the objec to the same pose
void CollisionChecker::Set_Object_Pose(Eigen::Matrix4f pose)
{
  for(unsigned int j = 0; j < objects_.size(); j++) {
    objects_[j]->Set_Transformation(pose);
  }
}

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
CollisionChecker::Check_Intersection(Eigen::Matrix4f pose)
{
  Set_Object_Pose(pose);

  //set up distances
  collision_detected_ = false;
  for(int i = 0; i < 3; i++) {
    global_point_1_[i] = 0.0;
    global_point_2_[i] = 0.0;
  }
  distance_sq_ = 0.0;

  //check for collision
  int num_collisions = Get_Penetration_Depth();

  //set up return
  CGAL::Triple<bool, DT_Vector3, DT_Vector3> triple;
  triple.first = collision_detected_;

  //iterate to find separation distance. TODO: fix very ugly O(n^2) implementation
  if (!collision_detected_) {
    distance_sq_ = 100000.0;
    Get_Closest_Pairs();
  }
  for (int i = 0; i < 3; i++) {
    triple.second[i] = global_point_1_[i];
    triple.third[i] = global_point_2_[i];
  }

#ifdef DEBUG
  cout << "Number of collisions: " << num_collisions << endl;
  if (num_collisions > 0) {
    cout << "Penetration depth squared: " << distance_squared << endl;
  }
  else {
    cout << "Separation distance squared: " << distance_squared << endl;
  }
  cout << "Point1: " << global_point_1[0] << " " << global_point_1[1] << " " << global_point_1[2] << endl;
  cout << "Point2: " << global_point_2[0] << " " << global_point_2[1] << " " << global_point_2[2] << endl;
#endif
  return triple;
}

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

  // loop through robot gripper fingers and objects to look for closest pairs of points
  for (unsigned int i = 0; i < gripper_fingers_.size(); i++) {
    for (unsigned int j = 0; j < objects_.size(); j++) {

      // find closest points
      gripper_fingers_[i]->GetClosestPair(objects_[j], point_1, point_2);
      new_distance_squared = Distance_Squared(point_1, point_2);
        
      // update minimum separation distance
      if(new_distance_squared < distance_sq_) {
        for(int i = 0; i < 3; i++) {
          global_point_1_[i] = point_1[i];
          global_point_2_[i] = point_2[i];
        }
        distance_sq_ = new_distance_squared;
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

  for (unsigned int i = 0; i < gripper_fingers_.size(); i++) {
    for (unsigned int j = 0; j < objects_.size(); j++) {
      
      // find penetration depth
      intersection = DT_GetPenDepth(gripper_fingers_[i]->getHandle(), objects_[j]->getHandle(), point_1, point_2);

      if (intersection) {
        // update num collisions
        num_coll++;
        collision_detected_ = true;
        new_distance_squared = Distance_Squared(point_1, point_2);
        
        // update maximum penetration depth
        if (new_distance_squared > distance_sq_) {
          for (int i = 0; i < 3; i++) {
            global_point_1_[i] = point_1[i];
            global_point_2_[i] = point_2[i];
          }
          distance_sq_ = new_distance_squared;
        }
      }
    }
  }
  return num_coll;
}

