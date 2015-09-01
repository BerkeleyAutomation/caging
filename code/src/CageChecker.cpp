#include "CageChecker.h"
#include "Grasp.hpp"
#include "StaticCageSimulator.h"

#include <glog/logging.h>

#define GOAL_DIST_THRESH 2.5f
#define IMAGE_SCALE 15.0f

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  CageChecker::CageChecker
 *  Description:  
 * =====================================================================================
 */

CageChecker::CageChecker (CfApproxConfig cf_config, Mesh* object, std::vector<Mesh*> gripper_fingers)
  : cf_approximator_(NULL),    
    object_(object),
    gripper_fingers_(gripper_fingers),
    cf_config_(cf_config), 
    gv_(NULL)
{ 
  Set_Configuration_Space_Approximator();
}		/* -----  end of function CageChecker::CageChecker  ----- */
    
CageChecker::CageChecker (CGAL::Geomview_stream* gv, CfApproxConfig cf_config, Mesh* object, std::vector<Mesh*> gripper_fingers)
  : cf_approximator_(NULL),    
    object_(object),
    gripper_fingers_(gripper_fingers),
    cf_config_(cf_config),
    gv_(gv)
{ 
  Set_Configuration_Space_Approximator();
}

CageChecker::~CageChecker()
{
  delete cf_approximator_;
  // delete object_;

  // for (unsigned int i = 0; i < gripper_fingers_.size(); i++) {
  //   delete gripper_fingers_[i];
  // }
}

void CageChecker::Draw_Object()
{
  (*gv_) << CGAL::GREEN;
  object_->RenderToGeomview(*gv_);
}

void CageChecker::Draw_Gripper()
{
  (*gv_) << CGAL::RED;

  for (unsigned int i = 0; i < gripper_fingers_.size(); i++) {
    gripper_fingers_[i]->RenderToGeomview(*gv_);
  }
}

void CageChecker::Draw_Collision(CGAL::Triple<bool, DT_Vector3, DT_Vector3> coll)
{
  // gv_ << CGAL::BLUE;
  // LineSegment line_seg = ShapeFactory::CreateLine(coll.second, coll.third, extrude_width_); 
  // gv_ << line_seg;
}

void CageChecker::Render_Object(cv::Mat& image)
{
  object_->RenderToImage(image, 'g');
}

void CageChecker::Render_Gripper(cv::Mat& image)
{
  for (unsigned int i = 0; i < gripper_fingers_.size(); i++) {
    gripper_fingers_[i]->RenderToImage(image, 'r');
  }
}

void CageChecker::Render_Object_Tris(cv::Mat& image)
{
  bool draw_com = true;
  object_->RenderTrisToImage(image, IMAGE_SCALE, 'b', draw_com);
}

void CageChecker::Render_Gripper_Tris(cv::Mat& image)
{
  for (unsigned int i = 0; i < gripper_fingers_.size(); i++) {
    gripper_fingers_[i]->RenderTrisToImage(image, IMAGE_SCALE, 'r');
  }
}

bool CageChecker::Set_Object_Pose(float tx, float ty, float theta)
{
  Eigen::Matrix4f pose = CreatePose(tx, ty, theta);
  // std::cout << "SetObj" << std::endl;
  // std::cout << pose << std::endl;
  object_->SetPose(pose, true); // TODO: make this actualy work
}

void CageChecker::Set_Gripper_Pose(float tx, float ty, float theta)
{
  for (unsigned int j = 0; j < gripper_fingers_.size(); j++) {
    Eigen::Matrix4f pose = CreatePose(tx, ty, theta);
    // std::cout << "SetGrip" << std::endl;
    // std::cout << pose << std::endl;
    gripper_fingers_[j]->SetPose(pose, true);
  }
}

// rotate everythin about hte origin by theta
void CageChecker::Rotate_World(float theta)
{
  Eigen::Matrix4f rot_pose = CreatePose(0, 0, theta);
  Eigen::Matrix4f obj_pose = object_->Pose();
  object_->SetPose(rot_pose * obj_pose);

  for (unsigned int j = 0; j < gripper_fingers_.size(); j++) {
    Eigen::Matrix4f gripper_pose = gripper_fingers_[j]->Pose();
    // std::cout << "Rotate" << std::endl;
    // std::cout << rot_pose*gripper_pose << std::endl;
    gripper_fingers_[j]->SetPose(rot_pose * gripper_pose);
  }
}

// Sets the params of the triangulation planner_
void CageChecker::Set_Configuration_Space_Approximator()
{
  if (cf_approximator_ != NULL) {
    delete cf_approximator_;
  }
  
  cf_approximator_ = new Configuration_Space_Approximator(cf_config_, (*gv_));
  cf_approximator_->set_object(object_);
  cf_approximator_->set_obstacles(gripper_fingers_);
}

void CageChecker::Clear_Planner()
{
  if(cf_approximator_ != NULL)
    cf_approximator_->reset();
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  CageChecker::Find_Path
 *  Description:  Wrapper around the planner's find path function
 * =====================================================================================
 */
CageResult
CageChecker::Check_Cage(EscapeEnergyConfig energy_config)
{
  // in future call the function with the default value
  Pose2D p;
  Pose3DToParams2D(object_->Pose(), p.x, p.y, p.theta);
  return Check_Cage(p, energy_config);
}

CageResult
CageChecker::Check_Cage(Pose2D pose, EscapeEnergyConfig energy_config)
{
  timer_.reset();
  timer_.start();
  EscapeEnergyResult energy_result;
  energy_config.energy_min = GravityPotential::GRAVITY_ACCEL * object_->Mass() * cf_config_.x_min;
  energy_config.energy_max = GravityPotential::GRAVITY_ACCEL * object_->Mass() * cf_config_.x_max;

  LOG(INFO) << "Min " << energy_config.energy_min;
  LOG(INFO) << "Max " << energy_config.energy_max;

  cf_approximator_->min_escape_energy(pose, energy_config, energy_result);
  timer_.stop();
  LOG(INFO) << "Cage check time: " << timer_.time();

  bool caged = !energy_result.path_exists;
  CageResult cage_result;
  cage_result.energy_result = energy_result;
  cage_result.total_time = timer_.time();

  timer_.reset();
  return cage_result;
}

PathPlanningResult
CageChecker::Find_Escape_Path(float tx, float ty, float theta, float energy_thresh, float timeout, float range)
{
  float mult = 5.0f;
  PathPlanningParams params;
  params.start_tx = tx;
  params.start_ty = ty;
  params.start_theta = theta;
  params.goal_tx = mult*cf_config_.x_min + 0.1;
  params.goal_ty = 0.0f;
  params.goal_theta = 0;

  params.min_tx = mult*cf_config_.x_min;
  params.min_ty = mult*cf_config_.y_min;;
  params.min_theta = mult * M_PI * (-cf_config_.num_rots - 1);

  params.max_tx = mult*cf_config_.x_max;
  params.max_ty = mult*cf_config_.y_max;
  params.max_theta = mult * M_PI * cf_config_.num_rots;

  params.goal_dist_thresh = GOAL_DIST_THRESH;
  params.planner_range = range;

  CageEscapePlanner cep(object_, gripper_fingers_, energy_thresh);
  return cep.FindEscapePath(params, timeout);
}

SimulationResult
CageChecker::Ratio_Escape_Paths(float tx, float ty, float theta, int num_samples, float potential_thresh, float timeout)
{
  PathPlanningParams params;
  params.start_tx = tx;
  params.start_ty = ty;
  params.start_theta = theta;
  params.goal_tx = cf_config_.x_min + 0.1;
  params.goal_ty = cf_config_.y_min + 0.1;
  params.goal_theta = 0;

  params.min_tx = cf_config_.x_min;
  params.min_ty = cf_config_.y_min;
  params.min_theta = 2 * M_PI * (cf_config_.num_rots - 1);

  params.max_tx = cf_config_.x_max;
  params.max_ty = cf_config_.y_max;
  params.max_theta = 2 * M_PI * cf_config_.num_rots;

  bool use_gui = false;
  StaticCageSimulator scs(object_, gripper_fingers_);
  return scs.RatioEscapesBox2D(num_samples, use_gui, potential_thresh);
}

