#include "StaticCageSimulator.h"

#include <ompl/base/Planner.h>
#include <ompl/base/spaces/RealVectorStateSpace.h>
#include <ompl/geometric/planners/rrt/RRT.h>
#include <ompl/geometric/planners/rrt/RRTConnect.h>
#include <ompl/geometric/planners/rrt/RRTstar.h>
#include <ompl/geometric/SimpleSetup.h>

#include <Box2D/Box2D.h>
#include "Simulation/Render.h"
#include "Simulation/Test.h"
#include "Simulation/ConvergenceTest.h"
#include "glui/glui.h"

// global params for simulation (has to be global to use OGL... :/)
namespace Box2DSimulation {
  ConvergenceTest* convergence_test;
  int main_window;
  int width = 640;
  int height = 480;
  int frame_period = 16;
  Settings settings;
  bool sim_finished;
  float ratio_cages;
  GLUI* glui;
};

static void Resize(int32 w, int32 h)
{
  int tx, ty, tw, th;
  GLUI_Master.get_viewport_area(&tx, &ty, &tw, &th);
  glViewport(tx, ty, tw, th);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  float32 ratio = float32(tw) / float32(th);
  b2Vec2 extents(ratio * 25.0f, 25.0f);
  b2Vec2 lower = Box2DSimulation::settings.viewCenter - extents;
  b2Vec2 upper = Box2DSimulation::settings.viewCenter + extents;
  gluOrtho2D(lower.x, upper.x, lower.y, upper.y);
}

// This is used to control the frame rate (60Hz).
static void SimTimer(int)
{
  glutSetWindow(Box2DSimulation::main_window);
  glutPostRedisplay();
  glutTimerFunc(Box2DSimulation::frame_period, SimTimer, 0);
}

static void ConvergenceSimulationLoop()
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  Box2DSimulation::convergence_test->DrawTitle("Convergence Test");

  //  std::cout << "Step" << std::endl;
  Box2DSimulation::convergence_test->Step(&Box2DSimulation::settings);

  if (Box2DSimulation::convergence_test->IsFinished()) { // hack (need to figure out how to see states)
    Box2DSimulation::sim_finished = true;
  }

  glutSwapBuffers();
}

StaticCageSimulator::StaticCageSimulator(Mesh* object, std::vector<Mesh*> obstacles)
  : object_(object),
    obstacles_(obstacles)
{
}

float StaticCageSimulator::RatioPathsInCollision(int num_paths, PathPlanningParams params, float timeout)
{
  const int dim = 3;
  ompl::base::StateSpacePtr space(new ompl::base::RealVectorStateSpace(dim));

  // set the bounds of space
  float tol = 0.1f;
  ompl::base::RealVectorBounds bounds(dim);
  bounds.setLow(0, params.min_tx - tol);
  bounds.setLow(1, params.min_ty - tol);
  bounds.setLow(2, params.min_theta - tol);
  bounds.setHigh(0, params.max_tx + tol);
  bounds.setHigh(1, params.max_ty + tol);
  bounds.setHigh(2, params.max_theta + tol);
  space->as<ompl::base::RealVectorStateSpace>()->setBounds(bounds);

  // Construct a space information instance for this state space
  ompl::base::SpaceInformationPtr si(new ompl::base::SpaceInformation(space));

  // state validity checker
  PlanningCollisionChecker collision_checker(si, object_, obstacles_);
  si->setStateValidityChecker(ompl::base::StateValidityCheckerPtr(new AllValidChecker(si)));
  si->setup();

  // set start state
  ompl::base::ScopedState<> start(space);
  start->as<ompl::base::RealVectorStateSpace::StateType>()->values[0] = params.start_tx;
  start->as<ompl::base::RealVectorStateSpace::StateType>()->values[1] = params.start_ty;
  start->as<ompl::base::RealVectorStateSpace::StateType>()->values[2] = params.start_theta;

  // try out random paths (as chosen using randomly sampled goals and RRT)
  int num_in_collision = 0;
  float min_dim[dim] = {params.min_tx, params.min_ty, params.min_theta};
  float max_dim[dim] = {params.max_tx, params.max_ty, params.max_theta};
  float goal_state[dim];

  for (int i = 0; i < num_paths; i++) { 
    std::cout << "Iteration " << i << std::endl;
  
    // sample facet of boundary
    float rand_sample = dim * ((float)rand()/((float)RAND_MAX + 0.1f));
    int rand_index = (int)rand_sample;
    int use_max = (int)(2 * ((float)rand()/((float)RAND_MAX + 0.1f)));
    
    // sample relative position on face
    int ind_i = (rand_index+1) % dim;
    int ind_j = (rand_index+2) % dim;
    goal_state[ind_i] = ((float)rand()/(float)RAND_MAX) * (max_dim[ind_i] - min_dim[ind_i]) + min_dim[ind_i];
    goal_state[ind_j] = ((float)rand()/(float)RAND_MAX) * (max_dim[ind_j] - min_dim[ind_j]) + min_dim[ind_j];
    if (use_max) {
      goal_state[rand_index] = max_dim[rand_index];
    }
    else {
      goal_state[rand_index] = min_dim[rand_index];
    }

    // set goal state
    ompl::base::ScopedState<> goal(space);
    goal->as<ompl::base::RealVectorStateSpace::StateType>()->values[0] = goal_state[0];
    goal->as<ompl::base::RealVectorStateSpace::StateType>()->values[1] = goal_state[1];
    goal->as<ompl::base::RealVectorStateSpace::StateType>()->values[2] = goal_state[2];
    std::cout << "Goal " << goal_state[0] << " " << goal_state[1] << " " << goal_state[2] << std::endl;

    // create a problem instance
    ompl::geometric::SimpleSetupPtr pdef(new ompl::geometric::SimpleSetup(si));
    pdef->setStartAndGoalStates(start, goal);

    // construct our optimizing planner using the RRTstar algorithm.
    ompl::geometric::RRT* rrt = new ompl::geometric::RRT(si);
    rrt->setRange(0.1f);
    ompl::base::PlannerPtr optimizingPlanner(rrt);

    // set the problem instance for our planner to solve
    optimizingPlanner->setProblemDefinition(pdef->getProblemDefinition());
    optimizingPlanner->setup();

    // attempt to solve the planning problem within configured time
    ompl::base::PlannerStatus solved = optimizingPlanner->solve(timeout);
    std::vector<ompl::base::State*> soln_states = pdef->getSolutionPath().getStates();

    if (solved) {      
      // try out the path, check for collisions
      bool path_in_collision = false; 
      for (unsigned int i = 0; i < soln_states.size() && !path_in_collision; i++) {
        if (!collision_checker.isValid(soln_states[i])) {
          std::cout << "Collision at state " << i << " of " << soln_states.size() << std::endl;
          path_in_collision = true;
          num_in_collision++;

          Eigen::Vector3f state;
          state << soln_states[i]->as<ompl::base::RealVectorStateSpace::StateType>()->values[0],
            soln_states[i]->as<ompl::base::RealVectorStateSpace::StateType>()->values[1],
            soln_states[i]->as<ompl::base::RealVectorStateSpace::StateType>()->values[2];
          std::cout << "State " << state.transpose() << std::endl;
        }
      }

      if (!path_in_collision) {
        std::cout << "Path exists!" << std::endl;
      }
    }
  }

  float ratio_caged = (float)num_in_collision / (float)num_paths;
  return ratio_caged;
}

float
StaticCageSimulator::RatioPathsFound(int num_paths, PathPlanningParams params, float timeout)
{
  const int dim = 3;
  ompl::base::StateSpacePtr space(new ompl::base::RealVectorStateSpace(dim));

  // set the bounds of space
  float tol = 0.1f;
  ompl::base::RealVectorBounds bounds(dim);
  bounds.setLow(0, params.min_tx - tol);
  bounds.setLow(1, params.min_ty - tol);
  bounds.setLow(2, params.min_theta - tol);
  bounds.setHigh(0, params.max_tx + tol);
  bounds.setHigh(1, params.max_ty + tol);
  bounds.setHigh(2, params.max_theta + tol);
  space->as<ompl::base::RealVectorStateSpace>()->setBounds(bounds);

  // Construct a space information instance for this state space
  ompl::base::SpaceInformationPtr si(new ompl::base::SpaceInformation(space));
  si->setStateValidityChecker(ompl::base::StateValidityCheckerPtr(new PlanningCollisionChecker(si, object_, obstacles_)));
  si->setup();

  // set start state
  ompl::base::ScopedState<> start(space);
  start->as<ompl::base::RealVectorStateSpace::StateType>()->values[0] = params.start_tx;
  start->as<ompl::base::RealVectorStateSpace::StateType>()->values[1] = params.start_ty;
  start->as<ompl::base::RealVectorStateSpace::StateType>()->values[2] = params.start_theta;

  // try out random paths (as chosen using randomly sampled goals and RRT)
  int num_in_collision = 0;
  float min_dim[dim] = {params.min_tx, params.min_ty, params.min_theta};
  float max_dim[dim] = {params.max_tx, params.max_ty, params.max_theta};
  float goal_state[dim];

  for (int i = 0; i < num_paths; i++) { 
    std::cout << "Iteration " << i << std::endl;
  
    // sample facet of boundary
    float rand_sample = dim * ((float)rand()/((float)RAND_MAX + 0.1f));
    int rand_index = (int)rand_sample;
    int use_max = (int)(2 * ((float)rand()/((float)RAND_MAX + 0.1f)));
    
    // sample relative position on face
    int ind_i = (rand_index+1) % dim;
    int ind_j = (rand_index+2) % dim;
    goal_state[ind_i] = ((float)rand()/(float)RAND_MAX) * (max_dim[ind_i] - min_dim[ind_i]) + min_dim[ind_i];
    goal_state[ind_j] = ((float)rand()/(float)RAND_MAX) * (max_dim[ind_j] - min_dim[ind_j]) + min_dim[ind_j];
    if (use_max) {
      goal_state[rand_index] = max_dim[rand_index];
    }
    else {
      goal_state[rand_index] = min_dim[rand_index];
    }

    // set goal state
    ompl::base::ScopedState<> goal(space);
    goal->as<ompl::base::RealVectorStateSpace::StateType>()->values[0] = goal_state[0];
    goal->as<ompl::base::RealVectorStateSpace::StateType>()->values[1] = goal_state[1];
    goal->as<ompl::base::RealVectorStateSpace::StateType>()->values[2] = goal_state[2];
    std::cout << "Goal " << goal_state[0] << " " << goal_state[1] << " " << goal_state[2] << std::endl;

    // create a problem instance
    ompl::geometric::SimpleSetupPtr pdef(new ompl::geometric::SimpleSetup(si));
    pdef->setStartAndGoalStates(start, goal);

    // construct our optimizing planner using the RRTstar algorithm.
    ompl::geometric::RRT* rrt = new ompl::geometric::RRT(si);
    rrt->setRange(0.005f);
    ompl::base::PlannerPtr optimizingPlanner(rrt);

    // set the problem instance for our planner to solve
    optimizingPlanner->setProblemDefinition(pdef->getProblemDefinition());
    optimizingPlanner->setup();

    // attempt to solve the planning problem within configured time
    ompl::base::PlannerStatus solved = optimizingPlanner->solve(timeout);

    if (solved == ompl::base::PlannerStatus::APPROXIMATE_SOLUTION) {
      std::cout << "Approx solution" << std::endl;      
    }
    else if (solved == ompl::base::PlannerStatus::EXACT_SOLUTION) {
      std::cout << "Exact solution" << std::endl;      
    }

    // check if solved in the allotted time
    if (!solved) {     
      std::cout << "Path does not exist" << std::endl;
      num_in_collision++;
    }
  }

  float ratio_caged = (float)num_in_collision / (float)num_paths;
  return ratio_caged;
}

SimulationResult
StaticCageSimulator::RatioEscapesBox2D(unsigned int num_trials, bool use_gui, Potential potential_thresh)
{
  // create test for given params
  Box2DSimulation::convergence_test = new ConvergenceTest(num_trials, object_, obstacles_, potential_thresh);
  Box2DSimulation::sim_finished = false;

  if (use_gui) {
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
    glutInitWindowSize(Box2DSimulation::width, Box2DSimulation::height);

    // add title
    char title[32];
    sprintf(title, "Box2D Version %d.%d.%d", b2_version.major, b2_version.minor, b2_version.revision);
    Box2DSimulation::main_window = glutCreateWindow(title);
    glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_GLUTMAINLOOP_RETURNS);

    // add functions
    glutDisplayFunc(ConvergenceSimulationLoop);
    GLUI_Master.set_glutReshapeFunc(Resize);  
    Box2DSimulation::glui = GLUI_Master.create_glui_subwindow(Box2DSimulation::main_window, GLUI_SUBWINDOW_RIGHT);

    // Use a timer to control the frame rate.
    glutTimerFunc(Box2DSimulation::frame_period, SimTimer, 0);
    glutMainLoop();

    // block until finished
    while (!Box2DSimulation::sim_finished) {
      usleep(100000);
    }
  }
  else {

    // non-gui main loop
    while (!Box2DSimulation::sim_finished) {
      Box2DSimulation::convergence_test->Step(&Box2DSimulation::settings);
      Box2DSimulation::sim_finished = Box2DSimulation::convergence_test->IsFinished();
    }
  }

  // grab the cage ratio and quit
  std::vector<float> cage_rates = Box2DSimulation::convergence_test->CageRates();
  std::vector<float> raw_cage_rates = Box2DSimulation::convergence_test->RawCageRates();
  delete Box2DSimulation::convergence_test;

  SimulationResult result;
  result.cage_rate = cage_rates[cage_rates.size()-1];
  result.raw_cage_rate = raw_cage_rates[cage_rates.size()-1];
  return result;
}
