#include "WorkspaceObject.h"

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  PRM::operator>>
 *  Description:  
 * =====================================================================================
 */
    

bool
WorkspaceObject::Load_Polyhedron (std::istream& in)  /* loads in a new polygon for the object */
{

  in >> object;
  transformed_object = object;
  
  Nef_polyhedron_3 N(object);



  CGAL::convex_decomposition_3(N);
  
  // the first volume is the outer volume, which is 
  // ignored in the decomposition
  Volume_const_iterator ci = ++N.volumes_begin();
  for( ; ci != N.volumes_end(); ++ci) {
    if(ci->mark()) {
      Polyhedron_3 P;
      N.convert_inner_shell_to_polyhedron(ci->shells_begin(), P);
      convex_pieces.push_back(P);
      transformed_convex_pieces.push_back(P);
    }
  }
  return true;
}//end function operator>>


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  WorkspaceObject::WorkspaceObject
 *  Description:  
 * =====================================================================================
 */
WorkspaceObject::WorkspaceObject ()
{
}		/* -----  end of function WorkspaceObject::WorkspaceObject  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  WorkspaceObject::Load_Object
 *  Description:  
 * =====================================================================================
 */
bool
WorkspaceObject::Load_Polyhedron (const Polyhedron_3& P, Point_2 C)
{
  object = P;
  transformed_object = object;
  centroid = C;

  Nef_polyhedron_3 N(object);

  CGAL::convex_decomposition_3(N);
  
  // the first volume is the outer volume, which is 
  // ignored in the decomposition
  Volume_const_iterator ci = ++N.volumes_begin();
  for( ; ci != N.volumes_end(); ++ci) {
    if(ci->mark()) {
      Polyhedron_3 P;
      N.convert_inner_shell_to_polyhedron(ci->shells_begin(), P);
      convex_pieces.push_back(P);
      transformed_convex_pieces.push_back(P);
    }
  }
  return true;

}		/* -----  end of function WorkspaceObject::Load_Object  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  WorkspaceObject::Transform
 *  Description:  This function transforms the transformed_object polyhedron using the
 *                given affine transformation.
 * =====================================================================================
 */
    void
WorkspaceObject::Transform (CGAL_Aff_Transform transform)
{
  std::transform(object.points_begin(), object.points_end(), transformed_object.points_begin(), transform);
  std::list<Polyhedron_3>::iterator trans_iter = transformed_convex_pieces.begin();
  for(std::list<Polyhedron_3>::iterator iter = convex_pieces.begin(); iter != convex_pieces.end(); iter++) {
    std::transform(iter->points_begin(), iter->points_end(), trans_iter->points_begin(), transform);
    trans_iter++;
  }
}		/* -----  end of function WorkspaceObject::Transform  ----- */

// Computes the maximum moment arm of the objects
float
WorkspaceObject::Moment_Arm()
{
  float max_moment_arm = 0.0f;
  float dist = 0.0f;
  float x = 0.0f;
  float y = 0.0f;
  float cx = CGAL::to_double(centroid.x());
  float cy = CGAL::to_double(centroid.y());

  for(Point_iterator iter = object.points_begin(); iter != object.points_end(); iter++) {
    x = CGAL::to_double((*iter)[0]);
    y = CGAL::to_double((*iter)[1]);
    dist = (x - cx) * (x - cx) + (y - cy)  * (y - cy);

    if (dist > max_moment_arm) {
      max_moment_arm = dist;
    }
  }
  return max_moment_arm;
}
