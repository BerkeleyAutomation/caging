/*
* Copyright (c) 2006-2009 Erin Catto http://www.box2d.org
*
* This software is provided 'as-is', without any express or implied
* warranty.  In no event will the authors be held liable for any damages
* arising from the use of this software.
* Permission is granted to anyone to use this software for any purpose,
* including commercial applications, and to alter it and redistribute it
* freely, subject to the following restrictions:
* 1. The origin of this software must not be misrepresented; you must not
* claim that you wrote the original software. If you use this software
* in a product, an acknowledgment in the product documentation would be
* appreciated but is not required.
* 2. Altered source versions must be plainly marked as such, and must not be
* misrepresented as being the original software.
* 3. This notice may not be removed or altered from any source distribution.
*/

#ifndef CONVERGENCE_TEST_H
#define CONVERGENCE_TEST_H

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

#include <string>
#include <vector>

#include <Box2D/Dynamics/b2Separator.h>

#include "Mesh.hpp"
#include "PotentialFunction.h"
#include "Simulation/Test.h"
#include "Simulation/Render.h"

enum Box2DConvergenceState
{
  START_TRIAL = 0,
  CUR_TRIAL,
  LOG_RESULTS,
  TRIALS_COMPLETE
};

class ConvergenceTest : public Test
{
 public:
  ConvergenceTest(unsigned int num_iters, Mesh* object, std::vector<Mesh*> obstacles, Potential potential_thresh = FLT_MAX);

  // NOTE: constructor not up-to-date with potential thresholding. Do not use
  ConvergenceTest(unsigned int num_iters, Mesh* object, std::vector<Mesh*> obstacles,
                  b2Vec2 gravity_vec, float dir_scale, float force_scale, float cage_tol,
                  float object_density, float object_friction, float num_timesteps);
  
 public:
  void SetGlobalVars();
  void InitBox2D();
  void Step(Settings* settings);

  // returns results of simulation
  bool IsFinished();
  std::vector<float> CageRates();
  std::vector<float> RawCageRates();

 private:
  b2Vec2 RandomForce();
  float RandomTorque();
  int CheckValid(std::vector<b2Vec2> vertices);

 private:
  // box2d storage of items
  b2Body* object_body_;
  std::vector<b2Body*> obstacle_bodies_;
  b2Vec2  gravity_vec_;

  // object state for simulator
  b2Vec2 object_cur_dir_;
  float object_cur_angle_;
  float object_cur_torque_;
  std::vector<float> cage_rates_;
  std::vector<float> raw_cage_rates_;
  unsigned int num_escapes_; // total failed cages with potential thresh
  unsigned int num_raw_escapes_; // total number of escapes, ignoring thresholding 

  // mesh obejcts in workspace
  Box2DConvergenceState state_; // current state
  Mesh* object_;
  std::vector<Mesh*> obstacles_;

  // initial pose of object
  b2Vec2 box2d_center_pos_;
  b2Vec2 object_start_pos_;
  float object_start_angle_;

  // random other params
  float obstacle_radius_; // max radius of obstacles
  float dir_perturb_scale_; // scale of direction perturbations
  float force_sigma_; // scale of force relative to circle size
  float torque_sigma_; // scale of torque relative to circle size
  boost::mt19937 rng_;
  boost::normal_distribution<float> force_norm_;
  boost::normal_distribution<float> torque_norm_;
  boost::variate_generator<boost::mt19937, boost::normal_distribution<float> > force_distribution_;
  boost::variate_generator<boost::mt19937, boost::normal_distribution<float> > torque_distribution_;
  b2Vec2 prev_force_perturb_;
  float prev_torque_perturb_;

  float object_mass_;
  float object_inertia_;
  Potential max_potential_; // the max potential reached on one simulation
  Potential potential_thresh_; // the threshold of permissible potentials

  float cage_tol_; // tolerance of caging
  float object_density_; // mass density of object
  float object_friction_; // friction coef of grippers and object

  // variables for logging trial progress, transitions
  unsigned int cur_iter_; // current num iters
  unsigned int num_iters_; // total iters to run
  unsigned int cur_timestep_;
  unsigned int num_timesteps_;
  bool finished_;
  bool trial_in_progress_;
};

#endif
