#ifndef CAGE_ESCAPE_PLANNER_H
#define CAGE_ESCAPE_PLANNER_H

#include <ompl/base/objectives/PathLengthOptimizationObjective.h>
#include <ompl/base/objectives/MinimaxObjective.h>
#include <ompl/base/ProblemDefinition.h>
#include <ompl/base/StateSpace.h>
#include <ompl/base/StateValidityChecker.h>

#include "Typedef.h"
#include "Mesh.hpp"
#include "Util.h"
#include "PotentialFunction.h"

struct PathPlanningParams {
  float start_tx, start_ty, start_theta;
  float goal_tx, goal_ty, goal_theta;
  float min_tx, max_tx;
  float min_ty, max_ty;
  float tx_scale, ty_scale, theta_scale;
  float min_theta, max_theta;
  float goal_dist_thresh;
  float planner_range;
  float energy_angle;
};

struct PathPlanningResult {
  bool path_exists;
  std::vector<Eigen::Vector3f> states;
  float max_energy;
  float normalized_max_energy;
  int max_state;
};

class PlanningCollisionChecker : public ompl::base::StateValidityChecker
{
 public:
  PlanningCollisionChecker(const ompl::base::SpaceInformationPtr& si, Mesh* object, std::vector<Mesh*> obstacles);

 public:
  // checks whether or not the object at | state | collides with the obstacles 
  bool isValid(const ompl::base::State* state) const;

private:
  Mesh* object_;
  std::vector<Mesh*> obstacles_;
};

class PlanningCollisionAndEnergyChecker : public ompl::base::StateValidityChecker
{
 public:
  PlanningCollisionAndEnergyChecker(const ompl::base::SpaceInformationPtr& si, Mesh* object, std::vector<Mesh*> obstacles,
                                    float energy_thresh, PotentialFunction* potential_func, ompl::base::ScopedState<> start);

 public:
  // checks whether or not the object at | state | collides with the obstacles 
  bool isValid(const ompl::base::State* state) const;

private:
  Mesh* object_;
  std::vector<Mesh*> obstacles_;
  float energy_thresh_;
  float tx_start_;
  float ty_start_;
  float theta_start_;
  PotentialFunction* potential_func_;
};

class WorkOptimizationObjective : public ompl::base::PathLengthOptimizationObjective
{
 public:
  WorkOptimizationObjective(const ompl::base::SpaceInformationPtr& si, float mass, float moment_of_intertia, float cost_offset = 1e4);

 public:
  ompl::base::Cost motionCost(const ompl::base::State* s1, const ompl::base::State* s2) const;
  ompl::base::Cost motionCostHeuristic(const ompl::base::State* s1, const ompl::base::State* s2) const;

 private:
  float mass_;
  float moment_of_inertia_;
  float cost_offset_;
};

class MinimaxGravityOptimizationObjective : public ompl::base::MinimaxObjective
{
 public:
  MinimaxGravityOptimizationObjective(const ompl::base::SpaceInformationPtr& si, float mass, float moment_of_inertia, float h0,
                                      PotentialFunction* potential_func, float cost_offset = 1e4);

 public:
  ompl::base::Cost stateCost(const ompl::base::State* s) const;

 private:
  float mass_;
  float moment_of_inertia_;
  float cost_offset_;
  float h0_; // initial object height
  PotentialFunction* potential_func_;
};

class CageEscapePlanner
{
 public:
  CageEscapePlanner(Mesh* object, std::vector<Mesh*> obstacles, float energy_thresh);
  ~CageEscapePlanner() {}

 public:
  PathPlanningResult FindEscapePath(PathPlanningParams params, float timeout = 1.0f, float max_push_force = 1.0f);
  PathPlanningResult UpperBoundEscapeEnergy(PathPlanningParams params, float timeout = 1.0f, float max_push_force = 1.0f);
  float IntegrateEscapePathEnergy(PathPlanningParams params, std::vector<DifferentialTriangle> separator_tris,
                                  float timeout = 1e-4f);
 
 private:
  Mesh* object_;
  std::vector<Mesh*> obstacles_;
  float energy_thresh_;
};

#endif // CAGE_ESCAPE_PLANNER_H
