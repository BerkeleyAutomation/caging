#ifndef STATIC_CAGE_SIMULATOR_H
#define STATIC_CAGE_SIMULATOR_H

#include <ompl/base/ProblemDefinition.h>
#include <ompl/base/StateSpace.h>
#include <ompl/base/StateValidityChecker.h>

#include <eigen3/Eigen/Core>

#include "CageEscapePlanner.h"
#include "Mesh.hpp"
#include "PotentialFunction.h"

struct SimulationResult
{
  float cage_rate;
  float raw_cage_rate;
};

class AllValidChecker : public ompl::base::StateValidityChecker
{
 public:
  AllValidChecker(const ompl::base::SpaceInformationPtr& si) : ompl::base::StateValidityChecker(si) {}

 public:
  // checks whether or not the object at | state | collides with the obstacles 
  bool isValid(const ompl::base::State* state) const { return true; }
};

class StaticCageSimulator
{
 public:
  StaticCageSimulator(Mesh* object, std::vector<Mesh*> obstacles);
  ~StaticCageSimulator() {}

 public:
  float RatioPathsInCollision(int num_paths, PathPlanningParams params, float timeout = 1.0f);
  float RatioPathsFound(int num_paths, PathPlanningParams params, float timeout = 0.1f);
  SimulationResult RatioEscapesBox2D(unsigned int num_trials, bool use_gui = true, Potential potential_thresh = FLT_MAX);
 
 private:
  Mesh* object_;
  std::vector<Mesh*> obstacles_;  
};

#endif // STATIC_CAGE_SIMULATOR_H
