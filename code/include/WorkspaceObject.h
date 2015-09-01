/*
 * =====================================================================================
 *
 *       Filename:  WorkspaceObject.h
 *
 *    Description:  This is a class for a WorkspaceObject.  It is a polyhedron as well as a convex
 *                  decomposition of the polyhedron
 *
 *        Version:  1.0
 *        Created:  07/14/2011 09:15:03 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Zoe McCarthy (zm), ZoeMcCarthy12@gmail.com
 *        Company:  University of Illinois at Urbana-Champaign
 *
 * =====================================================================================
 */
#ifndef WORKSPACE_OBJ_H
#define WORKSPACE_OBJ_H

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/Nef_3/SNC_indexed_items.h>
#include <CGAL/convex_decomposition_3.h> 
#include "Polyhedron_iostream.h"
#include "ShapeFactory.hpp"
#include <CGAL/Aff_transformation_3.h>
#include <list>
#include <vector>
#include <algorithm>



typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
typedef CGAL::Aff_transformation_3<Kernel> CGAL_Aff_Transform;
typedef CGAL::Nef_polyhedron_3<Kernel, CGAL::SNC_indexed_items> Nef_polyhedron_3;
typedef Nef_polyhedron_3::Volume_const_iterator Volume_const_iterator;

/*
 * =====================================================================================
 *        Class:  WorkspaceObject
 *  Description:  This class serve as the 3d model for a workspace object
 *                such as a robot, or obstacle.  It contains a polyhedron as well as the
 *                convex decomposition of that object.
 * =====================================================================================
 */
class WorkspaceObject
{
    public:
        /* ====================  LIFECYCLE     ======================================= */
        WorkspaceObject ();                             /* constructor */

        /* ====================  ACCESSORS     ======================================= */


        /* ====================  MUTATORS      ======================================= */
        bool Load_Polyhedron (std::istream& in);  /* loads in a new polygon for the object */

        bool Load_Polyhedron (const Polyhedron_3& P, Point_2 C = Point_2(0.0f, 0.0f));  /* loads in a new polygon for the object */

        void Transform(CGAL_Aff_Transform transform); //transforms the transformed_object

        /* ====================  OPERATORS     ======================================= */
        float Moment_Arm();
        
 public:
        Polyhedron_3 Object() {return object;}
        Polyhedron_3 TransformedObject() {return transformed_object;}
        
        /* ====================  DATA MEMBERS     ======================================= */
        Polyhedron_3            transformed_object;
        std::list<Polyhedron_3> transformed_convex_pieces; 
        std::list<Polyhedron_3> convex_pieces; 
    protected:
        /* ====================  DATA MEMBERS  ======================================= */

    private:
        /* ====================  DATA MEMBERS  ======================================= */
        Polyhedron_3            object;
        Point_2 centroid;


}; /* -----  end of class WorkspaceObject  ----- */

#endif // WORKSPACE_OBJ_H
