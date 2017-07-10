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

#include "ShapeFactory.hpp"
#include "WorkspaceObject.h"
#include "ConfigurationMapper.h"
#include "Mesh.hpp"
#include "CageEscapePlanner.h"
#include "Configuration_Space_Approximator.h" 
#include "StaticCageSimulator.h"
#include "Util.h"
#include "Typedef.h"

#include <unistd.h> // for sleep()

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
  //Synthesize Grasps
  std::vector< std::vector<synthesis_info> > synthesize_grasps(EscapeEnergyConfig energy_config, int num_searches = 10, 
						float angle_sweep = 0, float angle_disc = M_PI/4, bool check_reachability = true, float max_push_force = 1.0f);
  std::vector< std::vector<synthesis_info> > synthesize_grasps(Pose2D pose, EscapeEnergyConfig energy_config, int num_searches = 10, 
						float angle_sweep = 0, float angle_disc = M_PI/4, bool check_reachability = true, float max_push_force = 1.0f);
  
  // Check cages
  CageResult Check_Cage(EscapeEnergyConfig energy_config);
  CageResult Check_Cage(Pose2D pose, EscapeEnergyConfig energy_config); 
  PathPlanningResult Find_Escape_Path(float tx, float ty, float theta, float energy_thresh, float timeout = 15.0f, float range = 1.0f);
  PathPlanningResult Upper_Bound_Escape_Energy(float tx, float ty, float theta, float energy_thresh, float energy_angle, float timeout = 15.0f, float range = 1.0f, float max_push_force = 1.0f);
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




