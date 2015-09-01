#include "CageEscapePlanner.h"

#include "PotentialFunction.h"

#include <glog/logging.h>

#include <ompl/base/Planner.h>
#include <ompl/base/objectives/PathLengthOptimizationObjective.h>
#include <ompl/base/spaces/RealVectorStateSpace.h>
#include <ompl/geometric/planners/prm/PRM.h>
#include <ompl/geometric/planners/prm/PRMstar.h>
#include <ompl/geometric/planners/rrt/RRTConnect.h>
#include <ompl/geometric/planners/rrt/RRTstar.h>
#include <ompl/geometric/PathGeometric.h>
#include <ompl/geometric/SimpleSetup.h>

#define BETA 1e-2 // 5e-1

PlanningCollisionChecker::PlanningCollisionChecker(const ompl::base::SpaceInformationPtr& si, Mesh* object, std::vector<Mesh*> obstacles)
  : ompl::base::StateValidityChecker(si),
    object_(object),
    obstacles_(obstacles)
{
}

bool PlanningCollisionChecker::isValid(const ompl::base::State* state) const
{
  // We know we're working with a RealVectorStateSpace in this
  // example, so we downcast state into the specific type.
  const ompl::base::RealVectorStateSpace::StateType* state3D =
    state->as<ompl::base::RealVectorStateSpace::StateType>();

  // Extract the robot's (x,y) position from its state
  double tx = state3D->values[0];
  double ty = state3D->values[1];
  double theta = state3D->values[2];
  bool collision = false;

  if (object_ != NULL && obstacles_.size() > 0) {
    Eigen::Matrix4f pose = CreatePose(tx, ty, theta);
    object_->SetPose(pose);

    // test computation of penetration depth
    for (unsigned int j = 0; j < obstacles_.size(); j++) {
      MeshCollisionResult c = Mesh::UpperBoundCollision(object_, obstacles_[j]);
      if (c.collision) {
        collision = true;
      }
    }
  }

  return !collision;
}

PlanningCollisionAndEnergyChecker::PlanningCollisionAndEnergyChecker(const ompl::base::SpaceInformationPtr& si, Mesh* object, std::vector<Mesh*> obstacles,
                                                   float energy_thresh, ompl::base::ScopedState<> start)
  : ompl::base::StateValidityChecker(si),
    object_(object),
    obstacles_(obstacles),
    energy_thresh_(energy_thresh)
{
  tx_start_ = start->as<ompl::base::RealVectorStateSpace::StateType>()->values[0];
  ty_start_ = start->as<ompl::base::RealVectorStateSpace::StateType>()->values[1];
  theta_start_ = start->as<ompl::base::RealVectorStateSpace::StateType>()->values[2];
}

bool PlanningCollisionAndEnergyChecker::isValid(const ompl::base::State* state) const
{
  // We know we're working with a RealVectorStateSpace in this
  // example, so we downcast state into the specific type.
  const ompl::base::RealVectorStateSpace::StateType* state3D =
    state->as<ompl::base::RealVectorStateSpace::StateType>();

  // Extract the robot's (x,y) position from its state
  double tx = state3D->values[0];
  double ty = state3D->values[1];
  double theta = state3D->values[2];
  bool collision = false;

  if (object_ != NULL && obstacles_.size() > 0) {
    Eigen::Matrix4f pose = CreatePose(tx, ty, theta);
    object_->SetPose(pose);

    // test computation of penetration depth
    for (unsigned int j = 0; j < obstacles_.size(); j++) {
      MeshCollisionResult c = Mesh::UpperBoundCollision(object_, obstacles_[j]);
      if (c.collision) {
        collision = true;
      }
    }
  }

  // check energy
  float energy = object_->Mass() * GravityPotential::GRAVITY_ACCEL * (tx - tx_start_);
  bool exceeded_energy = (energy > energy_thresh_);

  return !collision && !exceeded_energy;
}

WorkOptimizationObjective::WorkOptimizationObjective(const ompl::base::SpaceInformationPtr& si, float mass,
                                                     float moment_of_inertia, float cost_offset)
  : PathLengthOptimizationObjective(si),
    mass_(mass),
    moment_of_inertia_(moment_of_inertia),
    cost_offset_(cost_offset)
{
}

ompl::base::Cost WorkOptimizationObjective::motionCost(const ompl::base::State* s1, const ompl::base::State* s2) const

{
  // convert states for usability
  const ompl::base::RealVectorStateSpace::StateType* s1_3d =
    s1->as<ompl::base::RealVectorStateSpace::StateType>();
  const ompl::base::RealVectorStateSpace::StateType* s2_3d =
    s2->as<ompl::base::RealVectorStateSpace::StateType>();

  Eigen::Vector2f s1_cart;
  s1_cart << s1_3d->values[0], s1_3d->values[1];
  Eigen::Vector2f s2_cart;
  s2_cart << s2_3d->values[0], s2_3d->values[1];

  // calculate work from forces and torques, assuming unit acceleration
  float v_cart = fabs(s2_3d->values[1] - s1_3d->values[1]); // assume work along x axis only   //(s1_cart - s2_cart).norm();
  float work_cart = 1.0f * mass_ * v_cart;

  float w_angular = fabs(s1_3d->values[2] - s2_3d->values[2]);
  float work_angular = 1.0f * moment_of_inertia_ * w_angular; 

  // calculate gravity cost from the origin
  float dh = (s2_3d->values[0]);// - s1_3d->values[0]);
  float gravity_cost = mass_ * GravityPotential::GRAVITY_ACCEL * dh;

  // add cost offset to prevent negative edge weights
  // std::cout << "Work cart: " << work_cart << std::endl;
  // std::cout << "Work angular: " << work_angular << std::endl;
  // std::cout << "Work grav: " << gravity_cost << std::endl;

  return ompl::base::Cost(/*work_cart + work_angular + */gravity_cost + cost_offset_);
}

ompl::base::Cost WorkOptimizationObjective::motionCostHeuristic(const ompl::base::State* s1, const ompl::base::State* s2) const
{
  return motionCost(s1, s2);
}

CageEscapePlanner::CageEscapePlanner(Mesh* object, std::vector<Mesh*> obstacles, float energy_thresh)
  : object_(object),
    obstacles_(obstacles),
    energy_thresh_(energy_thresh)
{
}

PathPlanningResult CageEscapePlanner::FindEscapePath(PathPlanningParams params, float timeout)
{
  const int dim = 3;
  ompl::base::StateSpacePtr space(new ompl::base::RealVectorStateSpace(dim));

  // set the bounds of space
  ompl::base::RealVectorBounds bounds(dim);
  bounds.setLow(0, params.min_tx);
  bounds.setLow(1, params.min_ty);
  bounds.setLow(2, params.min_theta);
  bounds.setHigh(0, params.max_tx);
  bounds.setHigh(1, params.max_ty);
  bounds.setHigh(2, params.max_theta);
  space->as<ompl::base::RealVectorStateSpace>()->setBounds(bounds);

  // Construct a space information instance for this state space
  ompl::base::SpaceInformationPtr si(new ompl::base::SpaceInformation(space));

  // set start state
  ompl::base::ScopedState<> start(space);
  start->as<ompl::base::RealVectorStateSpace::StateType>()->values[0] = params.start_tx;
  start->as<ompl::base::RealVectorStateSpace::StateType>()->values[1] = params.start_ty;
  start->as<ompl::base::RealVectorStateSpace::StateType>()->values[2] = params.start_theta;
  
  // set goal state
  ompl::base::ScopedState<> goal(space);
  goal->as<ompl::base::RealVectorStateSpace::StateType>()->values[0] = params.goal_tx;
  goal->as<ompl::base::RealVectorStateSpace::StateType>()->values[1] = params.goal_ty;
  goal->as<ompl::base::RealVectorStateSpace::StateType>()->values[2] = params.goal_theta;

  std::cout << "Start " << params.start_tx << " " << params.start_ty << " " << params.start_theta << std::endl;
  std::cout << "Goal " << params.goal_tx << " " << params.goal_ty << " " << params.goal_theta << std::endl;

  // get center of mass for cost
  float mass = object_->Mass();
  float moment_of_inertia = object_->MomentOfInertia();
  float cost_offset = mass * GravityPotential::GRAVITY_ACCEL * (params.max_ty - params.min_ty); // max gravity penalty to keep nonegative
  float delta = 0.0f * mass * GravityPotential::GRAVITY_ACCEL; // unit increase in y coordinate
  std::cout << "Planning with energy thresh " << energy_thresh_ + delta << " and range " << params.planner_range << std::endl;

  // state validity checker
  si->setStateValidityChecker(ompl::base::StateValidityCheckerPtr(new PlanningCollisionAndEnergyChecker(si, object_, obstacles_,
                                                                                               energy_thresh_ + delta,
                                                                                               start)));
  si->setup();

  // create a problem instance
  ompl::base::ProblemDefinitionPtr pdef(new ompl::base::ProblemDefinition(si));
  pdef->setStartAndGoalStates(start, goal);
  // pdef->setOptimizationObjective(ompl::base::OptimizationObjectivePtr(new WorkOptimizationObjective(si, mass, moment_of_inertia, cost_offset)));

  // construct our optimizing planner using the RRTstar algorithm.
  ompl::geometric::RRTstar* rrt = new ompl::geometric::RRTstar(si);
  rrt->setRange(params.planner_range);
  ompl::base::PlannerPtr optimizingPlanner(rrt);//new ompl::geometric::RRTstar(si));

  // set the problem instance for our planner to solve
  optimizingPlanner->setProblemDefinition(pdef);//->getProblemDefinition());
  optimizingPlanner->setup();

  // attempt to solve the planning problem within configured time
  ompl::base::PlannerStatus solved = optimizingPlanner->solve(timeout);

  if (solved == ompl::base::PlannerStatus::APPROXIMATE_SOLUTION) {
    LOG(INFO) << "APPROX SOLUTION!";    
  }
  else if (solved == ompl::base::PlannerStatus::EXACT_SOLUTION) {
    LOG(INFO) << "EXACT SOLUTION!";
  }
  else {
    LOG(INFO) << "NO SOLUTION";
  }

  PathPlanningResult result;
  result.path_exists = false;

  if (!solved)
    std::cout << "NOT SOLVED" << std::endl;
  
  if (solved) {
    boost::shared_ptr<ompl::geometric::PathGeometric> solution_path = \
      boost::static_pointer_cast<ompl::geometric::PathGeometric>(pdef->getSolutionPath());

    float dist_from_goal = pdef->getSolutionDifference();
    if (dist_from_goal < params.goal_dist_thresh) {
      result.path_exists = true;

      // get the resulting path and max energy achieved
      float max_path_energy = 0.0f;
      int max_state = 0;
      std::vector<ompl::base::State*> soln_states = solution_path->getStates();
      for (unsigned int i = 0; i < soln_states.size(); i++) {
        // store state
        Eigen::Vector3f state;
        state << soln_states[i]->as<ompl::base::RealVectorStateSpace::StateType>()->values[0],
          soln_states[i]->as<ompl::base::RealVectorStateSpace::StateType>()->values[1],
          soln_states[i]->as<ompl::base::RealVectorStateSpace::StateType>()->values[2];
        result.states.push_back(state);

        // compute energy
        float energy = mass * GravityPotential::GRAVITY_ACCEL * soln_states[i]->as<ompl::base::RealVectorStateSpace::StateType>()->values[0];
        if (energy > max_path_energy) {
          max_path_energy = energy;
          max_state = i;
        }
      }
      result.max_energy = max_path_energy;
      result.normalized_max_energy = max_path_energy / (mass * GravityPotential::GRAVITY_ACCEL);
      result.max_state = max_state;
    }
  }

  return result;
}

float CageEscapePlanner::IntegrateEscapePathEnergy(PathPlanningParams params, std::vector<DifferentialTriangle> separator_tris,
                                                   float timeout)
{
  const int dim = 3;
  ompl::msg::setLogLevel(ompl::msg::LOG_ERROR);
  ompl::base::StateSpacePtr space(new ompl::base::RealVectorStateSpace(dim));

  // set the bounds of space
  ompl::base::RealVectorBounds bounds(dim);
  bounds.setLow(0, params.min_tx);
  bounds.setLow(1, params.min_ty);
  bounds.setLow(2, params.min_theta);
  bounds.setHigh(0, params.max_tx);
  bounds.setHigh(1, params.max_ty);
  bounds.setHigh(2, params.max_theta);
  space->as<ompl::base::RealVectorStateSpace>()->setBounds(bounds);

  std::vector<Eigen::Matrix4f> p = obstacles_[0]->PosesWorldFrame();
  float tx, ty, theta;
  for (unsigned int i = 0; i < p.size(); i++) {
    Pose3DToParams2D(p[i], tx, ty, theta);
  }

  // Construct a space information instance for this state space
  ompl::base::SpaceInformationPtr si(new ompl::base::SpaceInformation(space));

  // state validity checker
  PlanningCollisionChecker checker(si, object_, obstacles_);
  si->setStateValidityChecker(ompl::base::StateValidityCheckerPtr(new PlanningCollisionChecker(si, object_, obstacles_)));
  si->setup();

  // set start state
  ompl::base::ScopedState<> start(space);
  start->as<ompl::base::RealVectorStateSpace::StateType>()->values[0] = params.start_tx;
  start->as<ompl::base::RealVectorStateSpace::StateType>()->values[1] = params.start_ty;
  start->as<ompl::base::RealVectorStateSpace::StateType>()->values[2] = params.start_theta;

  // construct our optimizing planner using the RRTstar algorithm.
  ompl::geometric::PRMstar* optimizingPlanner = new ompl::geometric::PRMstar(si);
  float mass = object_->Mass();
  float moment_of_inertia = object_->MomentOfInertia();
  float cost_offset = mass * GravityPotential::GRAVITY_ACCEL * (params.max_ty - params.min_ty); // max gravity penalty
    
  // loop through the separator tris, plan to each one, and summ result
  float gibbs_prob = 0.0f;
  for (unsigned int i = 0; i < separator_tris.size(); i++) {
    // get next goal tri
    DifferentialTriangle tri = separator_tris[i];

    if (i % 100 == 0) {
      //      std::cout << "Planning for tri: " << i << std::endl;
      // std::cout << "Goal " << params.tx_scale * tri.center.x() << " " << \
      //   params.ty_scale * tri.center.y() << " " << params.theta_scale * tri.center.z() << std::endl;
    }

    // set goal state
    ompl::base::ScopedState<> goal(space);
    // goal->as<ompl::base::RealVectorStateSpace::StateType>()->values[0] = params.start_tx - 20; // params.tx_scale * tri.center.x();
    // goal->as<ompl::base::RealVectorStateSpace::StateType>()->values[1] = params.start_ty; //params.ty_scale * tri.center.y();
    // goal->as<ompl::base::RealVectorStateSpace::StateType>()->values[2] = params.start_theta; //params.theta_scale * tri.center.z();
 
    goal->as<ompl::base::RealVectorStateSpace::StateType>()->values[0] = params.tx_scale * tri.center.x();
    goal->as<ompl::base::RealVectorStateSpace::StateType>()->values[1] = params.ty_scale * tri.center.y();
    goal->as<ompl::base::RealVectorStateSpace::StateType>()->values[2] = params.theta_scale * tri.center.z();

    // check goal validity
    if (checker.isValid(goal.get())) {
      // set new problem definition
      ompl::base::ProblemDefinitionPtr pdef(new ompl::base::ProblemDefinition(si));
      pdef->setStartAndGoalStates(start, goal);
      pdef->setOptimizationObjective(ompl::base::OptimizationObjectivePtr(new WorkOptimizationObjective(si, mass, moment_of_inertia, cost_offset)));

      // set the problem instance for our planner to solve
      optimizingPlanner->clearQuery();
      optimizingPlanner->setProblemDefinition(pdef);//->getProblemDefinition());
      optimizingPlanner->setup();

      // attempt to solve the planning problem within configured time
      ompl::base::PlannerStatus solved = optimizingPlanner->solve(ompl::base::timedPlannerTerminationCondition(timeout));//timeout);

      // integrate if solved
      if (solved) {
        // check that problem was actually solved
        float dist_from_goal = pdef->getSolutionDifference();
        if (fabs(dist_from_goal) < params.goal_dist_thresh) { // && pdef->hasOptimizedSolution()) {
          // get solution path
          boost::shared_ptr<ompl::geometric::PathGeometric> solution_path = \
            boost::static_pointer_cast<ompl::geometric::PathGeometric>(pdef->getSolutionPath());

          // get "direction" of solution path
          std::vector<ompl::base::State*> soln_states = solution_path->getStates();
          unsigned int final_ind = soln_states.size() - 1;
          Eigen::Vector3f final_state;
          final_state << soln_states[final_ind]->as<ompl::base::RealVectorStateSpace::StateType>()->values[0],
            soln_states[final_ind]->as<ompl::base::RealVectorStateSpace::StateType>()->values[1],
            soln_states[final_ind]->as<ompl::base::RealVectorStateSpace::StateType>()->values[2];

          unsigned int sec_last_ind = soln_states.size() - 2;
          Eigen::Vector3f second_to_last_state;
          second_to_last_state << soln_states[sec_last_ind]->as<ompl::base::RealVectorStateSpace::StateType>()->values[0],
            soln_states[sec_last_ind]->as<ompl::base::RealVectorStateSpace::StateType>()->values[1],
            soln_states[sec_last_ind]->as<ompl::base::RealVectorStateSpace::StateType>()->values[2];

          // add new path cost
          float path_energy = solution_path->cost(pdef->getOptimizationObjective()).value();

          // remove "offset" energy but leave one multiple of offset per path to keep nonegative
          path_energy = path_energy - (soln_states.size() - 1) * cost_offset;

          // check path direction
          Eigen::Vector3f path_dir = final_state - second_to_last_state;

          if (path_dir.norm() > 0) {
            path_dir.normalize();
            Eigen::Vector3f area_norm;
            area_norm << tri.dir.dx(), tri.dir.dy(), tri.dir.dz();
            area_norm.normalize();

            // calculate alignment (negative due to default tri orientation)
            float alignment = -path_dir.dot(area_norm);

            if (alignment > 0) {
              gibbs_prob = gibbs_prob + exp(-BETA * path_energy) * tri.area;
            }
          }
        } 
      }
    }
  }

  return gibbs_prob;
}
