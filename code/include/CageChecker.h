/*
 * CageChecker.h
 * 
 * Header file for the CageChecker class that is used to check a given gripper
 *     
 * Author:  Jeff Mahler, jmahler@berkeley.edu
 *          Zoe McCarthy (zm), ZoeMcCarthy12@gmail.com
 * Affiliation:  University of California Berkeley
 *
 */
#ifndef CAGING_H
#define CAGING_H

#include <CGAL/Cartesian.h>
#include <CGAL/Aff_transformation_3.h>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

#include <iostream>
#include <CGAL/Timer.h>
#include "ShapeFactory.hpp"
#include "CollisionChecker.h"
#include "WorkspaceObject.h"
#include "ConfigurationMapper.h"
#include "Mesh.hpp"
#include "CageEscapePlanner.h"
#include "Configuration_Space_Approximator.h" 
#include "StaticCageSimulator.h"
#include "Util.h"

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

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

typedef CGAL::Timer                         Timer;

typedef std::vector<double>                                 Configuration;
typedef std::list< Configuration >                          Path;
typedef std::pair< Path, bool >                             Path_with_exist;

struct CageResult
{
  EscapeEnergyResult energy_result;
  float total_time;
};

/*
 * =====================================================================================
 *        Class:  CageChecker
 *  Description:  This class contains the machinery to check a gripper and object configuration
 *                for caging
 * =====================================================================================
 */
class CageChecker
{
 public:
  CageChecker(CfApproxConfig cf_config, Mesh* object, std::vector<Mesh*> gripper_fingers);
  CageChecker(CGAL::Geomview_stream* gv, CfApproxConfig cf_config, Mesh* object, std::vector<Mesh*> gripper_fingers);
  ~CageChecker();                

 public:
  // Visualization
  void Draw_Object();
  void Draw_Gripper();
  void Render_Object(cv::Mat& image);
  void Render_Gripper(cv::Mat& image);
  void Render_Object_Tris(cv::Mat& image);
  void Render_Gripper_Tris(cv::Mat& image);
  void Draw_Collision(CGAL::Triple<bool, DT_Vector3, DT_Vector3> coll);
  void Clear_Geomview() { gv_->clear();}
  
 public:
  // adds an object to grasp
  bool Set_Object_Pose(float tx, float ty, float theta);
  void Set_Gripper_Pose(float tx, float ty, float theta);

  // adds gripper "finger"
  // in the future this might be a single cage
  void Rotate_World(float theta);

  // setup configuration space approximator
  void Set_Configuration_Space_Approximator();
  
 public:
  // Check cages
  CageResult Check_Cage(EscapeEnergyConfig energy_config);
  CageResult Check_Cage(Pose2D pose, EscapeEnergyConfig energy_config);
  PathPlanningResult Find_Escape_Path(float tx, float ty, float theta, float energy_thresh, float timeout = 15.0f, float range = 1.0f);
  SimulationResult Ratio_Escape_Paths(float tx, float ty, float theta, int num_samples = 100, float potential_thresh = FLT_MAX, float timeout = 1.0f);
  void Clear_Planner();

 private:
  std::vector<Mesh*> gripper_fingers_;
  Mesh*              object_;
  Configuration_Space_Approximator* cf_approximator_;
  CfApproxConfig              cf_config_;
  CGAL::Geomview_stream*      gv_;
  Timer                       timer_;
};

#endif //CAGE_H




