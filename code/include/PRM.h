/*
 * =====================================================================================
 *
 *       Filename:  PRM.h
 *
 *    Description:  header file for the PRM class that is the main class used in the PRM
 *                  planner
 *
 *        Version:  1.0
 *        Created:  07/14/2011 06:57:29 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Zoe McCarthy (zm), ZoeMcCarthy12@gmail.com
 *        Company:  University of Illinois at Urbana-Champaign
 *
 * =====================================================================================
 */


#ifndef PRM_H
#define PRM_H

#include <CGAL/Cartesian.h>
#include <CGAL/Aff_transformation_3.h>
#include <iostream>
#include <CGAL/Timer.h>
#include "CollisionChecker.h"
#include "WorkspaceObject.h"
#include "ConfigurationMapper.h"
#include "Planner.h"
#include "Triangulation_Planner.h"

#include <list>
#include <vector>
#include <math.h>

#include <fstream>
#include <unistd.h> // for sleep()


#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/Triangulation_geomview_ostream_3.h>
#include <CGAL/IO/Polyhedron_geomview_ostream.h>
#include <CGAL/utility.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
typedef CGAL::Aff_transformation_3<Kernel> CGAL_Aff_Transform;


typedef CGAL::Timer                         Timer;

typedef std::vector<double>                                 Configuration;
typedef std::list< Configuration >                          Path;
typedef std::pair< Path, bool >                             Path_with_exist;



//typedef CGAL::Unique_hash_map<Vertex_handle, bool>         Point_map;

/*#include <CGAL/Triangulation_euclidean_traits_xy_3.h>
#include <CGAL/intersections.h>

#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Regular_triangulation_euclidean_traits_3<K> Gt;


typedef CGAL::Alpha_shape_vertex_base_3<Gt>         Vb;
typedef CGAL::Alpha_shape_cell_base_3<Gt>           Fb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb> Tds;
typedef CGAL::Regular_triangulation_3<Gt,Tds>       Triangulation_3;
typedef CGAL::Alpha_shape_3<Triangulation_3>        Alpha_shape_3;

typedef Alpha_shape_3::Cell_handle          Cell_handle;
typedef Alpha_shape_3::Vertex_handle        Vertex_handle;
typedef Alpha_shape_3::Facet                Facet;
typedef Alpha_shape_3::Edge                 Edge;
typedef Gt::Weighted_point                  Weighted_point;
typedef Gt::Bare_point                      Bare_point;
*/



/*
 * =====================================================================================
 *        Class:  PRM
 *  Description:  this class contains the machinery to implement a Probablistic RoadMap
 *                planner.  It will contain a robot model, obstacles, a Regular Triangulation
 *                for sampled configurations,  collision checking machinery, and an alpha shape
 *                it will support queries of path construction/existence and also display
 *                queries
 *
 * =====================================================================================
 */
class PRM
{


    public:
        /* ====================  LIFECYCLE     ======================================= */
        PRM (double bounding);                             /* constructor */


        ~PRM()                                             /* destructor */
        {
            delete conf_mapper;
            delete planner;
        }

        /* ====================  ACCESSORS     ======================================= */

        void Draw_Spheres();                             /* draw the spheres in the triangulation */

        bool Draw_Robot(Configuration configuration);                             /* draw the scene for a given parameter */


        void Clear_Geomview(){gv.clear(); draw_count = 0;}

        Path_with_exist Pad_For_Animation(Path_with_exist& in_path);                       /*  pads the path for animation */

//        void Draw_Alpha_Shape();                             /* draw the scene for a given parameter */

        /* ====================  MUTATORS      ======================================= */
        bool Load_Robot(std::istream& in);

        bool Load_Obstacle(std::istream& in);

        bool Load_Robot(const Polyhedron_3& P);

        bool Load_Obstacle (const Polyhedron_3& P);

        bool Load_Obstacle(const Polyhedron_3& P, CGAL_Aff_Transform transform);

        void Load_Configuration_Mapper (ConfigurationMapper* in){conf_mapper = in;}

        void Set_Triangulation_Planner(double in_xscale, double in_yscale, double in_zscale, 
        double in_xmin, double in_xmax, double in_ymin, double in_ymax, double in_zmin, double in_zmax,
        double extra_factor)
        {   if (planner != NULL)
            {
                delete planner;
            }
            planner = new Triangulation_Planner(in_xscale, in_yscale, in_zscale,
            in_xmin, in_xmax, in_ymin, in_ymax, in_zmin, in_zmax, 
            coll_checker, conf_mapper, gv, extra_factor);} //set up a triangulation_planner


        Path_with_exist Find_Path(Configuration& start, Configuration& end);

        void Clear_Planner(){if(planner != NULL) planner->Clear_Stored_Data();}


        /* ====================  OPERATORS     ======================================= */

        CGAL::Triple<bool, DT_Vector3, DT_Vector3> Check_Collision(Configuration& configuration);

    protected:
        /* ====================  DATA MEMBERS  ======================================= */

    private:
        /* ====================  DATA MEMBERS  ======================================= */
        std::list<WorkspaceObject>  robot;
        std::list<WorkspaceObject>  obstacles;
        CollisionChecker            coll_checker;
        ConfigurationMapper*        conf_mapper;
        Timer                       timer;
        CGAL::Geomview_stream       gv;
        Planner*                    planner;
        int                         draw_count;

}; /* -----  end of class PRM  ----- */


//START IMPLEMENTATION

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  PRM::PRM
 *  Description:  
 * =====================================================================================
 */
    
PRM::PRM (double bounding)
    :gv(CGAL::Bbox_3(-bounding, -bounding, -bounding, bounding, bounding, bounding))
{ 
  draw_count = 0;
  gv.set_line_width(2);
  // gv.set_trace(true);
  gv.set_bg_color(CGAL::Color(0, 200, 200));
  gv.set_wired(false);
  conf_mapper = NULL;
  planner = NULL;
}		/* -----  end of function PRM::PRM  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Draw
 *  Description:  This function draws the robot and obstacle in the geomview stream. 
 * =====================================================================================
 */
    bool
PRM::Draw_Robot(Configuration configuration)
{


  std::list<WorkspaceObject>::iterator iter;
  std::list<CGAL_Aff_Transform> transform_list;
  transform_list = conf_mapper->Display_Transformations(configuration);
  std::list<CGAL_Aff_Transform>::iterator transf_iter = transform_list.begin();

 
  for (iter = robot.begin(); iter != robot.end(); iter++) {
      iter->Transform(*transf_iter);
      if(draw_count == 0)
      {
          gv << CGAL::GREEN;
      }
      else
      {
          gv << CGAL::PURPLE;
      }
      gv << iter->transformed_object;
      transf_iter++;
  }

  for (iter = obstacles.begin(); iter != obstacles.end(); iter++) {
      gv << CGAL::RED;
      gv << iter->transformed_object;
  }

  CGAL::Triple<bool, DT_Vector3, DT_Vector3> coll_data(Check_Collision(configuration));

  //gv << CGAL::BLUE;
  //gv << Kernel::Triangle_3(Kernel::Point_3(coll_data.second[0], coll_data.second[1], coll_data.second[2] - 10.0), Kernel::Point_3(coll_data.second[0], coll_data.second[1], coll_data.second[2] + 10.0), Kernel::Point_3(coll_data.third[0], coll_data.third[1], coll_data.third[2] + 10.0));
  //gv << Kernel::Triangle_3(Kernel::Point_3(coll_data.second[0], coll_data.second[1], coll_data.second[2] - 10.0), Kernel::Point_3(coll_data.third[0], coll_data.third[1], coll_data.third[2] - 10.0), Kernel::Point_3(coll_data.third[0], coll_data.third[1], coll_data.third[2] + 10.0));
  draw_count++;

  return coll_data.first;
  

  } 		/* -----  end of function Draw  ----- */




/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Draw_Spheres
 *  Description:  This function draws the sampled spheres for the triangulation structure.
 * =====================================================================================
 */
    void
PRM::Draw_Spheres()
{

  gv.clear();


 
  /* for (iter = planner->point_in_collision.begin(); iter != planner->point_in_collision.end(); iter++) 
  {
      if(*
      gv << CGAL::GREEN;
      gv << iter->transformed_object;
      transf_iter++;
  }*/

  

  } 		/* -----  end of function Draw_Spheres  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  PRM::Load_Robot
 *  Description:  Pushes a new link onto the robot. Reads from an input stream in OFF format
 * =====================================================================================
 */
    bool
PRM::Load_Robot(std::istream& in) 
{
    WorkspaceObject object;
    bool ret = object.Load_Polyhedron(in);
    robot.push_back(object);
    coll_checker.Load_Robot(object);
    return ret;
}		/* -----  end of function PRM::Load_Robot  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  PRM::Load_Robot
 *  Description:  Pushes a new link onto the robot. Reads from a polyhedron model
 * =====================================================================================
 */
    bool
PRM::Load_Robot(const Polyhedron_3& P) 
{
    WorkspaceObject object;
    bool ret = object.Load_Polyhedron(P);
    robot.push_back(object);
    coll_checker.Load_Robot(object);
    return ret;
}		/* -----  end of function PRM::Load_Robot  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  PRM::Load_Obstacle
 *  Description:  Pushes a obstacle into the obstacles set. Reads from an input stream in OFF format
 * =====================================================================================
 */
    bool
PRM::Load_Obstacle(std::istream& in) 
{
    WorkspaceObject object;
    bool ret = object.Load_Polyhedron(in);
    obstacles.push_back(object);
    coll_checker.Load_Obstacle(object);
    return ret;
}		/* -----  end of function PRM::Load_Obstacle  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  PRM::Load_Obstacle
 *  Description:  Pushes a obstacle into the obstacles set. Reads from a polyhedron model
 * =====================================================================================
 */
    bool
PRM::Load_Obstacle(const Polyhedron_3& P) 
{
    WorkspaceObject object;
    bool ret = object.Load_Polyhedron(P);
    obstacles.push_back(object);
    coll_checker.Load_Obstacle(object);
    return ret;
}		/* -----  end of function PRM::Load_Robot  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  PRM::Load_Obstacle
 *  Description:  Pushes a obstacle into the obstacles set. Reads from a polyhedron model
 *                Also transforms the object before pushing it.
 * =====================================================================================
 */
    bool
PRM::Load_Obstacle(const Polyhedron_3& P, CGAL_Aff_Transform transform)
{
    WorkspaceObject object;
    bool ret = object.Load_Polyhedron(P);
    object.Transform(transform);
    obstacles.push_back(object);
    coll_checker.Load_Obstacle(object, transform);
    return ret;
}		/* -----  end of function PRM::Load_Robot  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  PRM::Check_Collision
 *  Description:  Links the CollisionChecker and ConfigurationMapper
 * =====================================================================================
 */
CGAL::Triple<bool, DT_Vector3, DT_Vector3>
PRM::Check_Collision (Configuration& configuration)
{
    return coll_checker.Check_Intersection(conf_mapper->Collision_Transformations(configuration));
}		/* -----  end of function PRM::Check_Collision  ----- */

;
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  PRM::Find_Path
 *  Description:  Wrapper around the planner's find path function
 * =====================================================================================
 */
    Path_with_exist
PRM::Find_Path (Configuration& start, Configuration& end)
{
    Path_with_exist ret_path;
    timer.reset();
    timer.start();
    ret_path = planner->Find_Path(start,end);
    timer.stop();
    std::cout << "Path Planning time: " << timer.time() << std::endl;
    timer.reset();
    return ret_path;
}		/* -----  end of function PRM::Find_Path  ----- */




Configuration Average_Configuration(Configuration p1, Configuration p2)
{
    Configuration ret_conf;
    for(int i = 0; i < 3; i++)
    {
        ret_conf.push_back((p1[i]+p2[i])/2.0);
    }
    return ret_conf;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  PRM::Pad_For_Animation
 *  Description:  Adds in extra 
 * =====================================================================================
 */


    Path_with_exist
PRM::Pad_For_Animation (Path_with_exist& in_path)
{
    Path_with_exist ret_path;
    ret_path.second = in_path.second;
    ret_path.first.push_back(in_path.first.front());
    Path::iterator iter = in_path.first.begin();
    iter++;
    int count = 0;
    for(; iter != in_path.first.end(); iter++)
    {
        while(Distance_Squared(*iter, ret_path.first.back()) > .5)
        {
            ret_path.first.push_back(Average_Configuration(*iter, ret_path.first.back()));
            count++;
        }
        ret_path.first.push_back(*iter);
    }
    std::cout << "Number of pads: " << count << std::endl;
    return ret_path;

}		/* -----  end of function PRM::Pad_For_Animation  ----- */


#endif //PRM_H




