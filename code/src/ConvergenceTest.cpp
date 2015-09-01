#include "Simulation/ConvergenceTest.h"

#include "Simulation/Box2DShapeFactory.hpp"
#include "PotentialFunction.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <unistd.h>

#include <glog/logging.h>

#define STD_MASS 0.03f

// global defaults
b2Vec2 g_gravityVec = b2Vec2(0.0f, -GravityPotential::GRAVITY_ACCEL);
float32 g_forceSigma = 0.25f;
float32 g_torqueSigma = M_PI / 256.0f;
float32 g_cagedTol = 1e-1f; // tolerance of caging

float32 g_dirScale = M_PI / 8;
int g_allowedTimesteps = 1000;

ConvergenceTest::ConvergenceTest(unsigned int num_iters, Mesh* object, std::vector<Mesh*> obstacles, Potential potential_thresh)
  : gravity_vec_(g_gravityVec),
    num_escapes_(0),
    num_raw_escapes_(0),
    state_(START_TRIAL),
    object_(object),
    obstacles_(obstacles),
    obstacle_radius_(0.0f),
    dir_perturb_scale_(g_dirScale),
    force_sigma_(g_forceSigma),
    torque_sigma_(g_torqueSigma),
    force_norm_(0.0f, force_sigma_),
    torque_norm_(0.0f, torque_sigma_),
    force_distribution_(rng_, force_norm_),
    torque_distribution_(rng_, torque_norm_),
    potential_thresh_(potential_thresh),
    cage_tol_(g_cagedTol),
    object_density_(Box2DShapeFactory::default_density),
    object_friction_(Box2DShapeFactory::default_friction),
    num_timesteps_(g_allowedTimesteps),
    cur_iter_(0),
    num_iters_(num_iters),
    finished_(false),
    trial_in_progress_(false)
{
  SetGlobalVars();
  InitBox2D();
}

ConvergenceTest::ConvergenceTest(unsigned int num_iters, Mesh* object, std::vector<Mesh*> obstacles,
                                 b2Vec2 gravity_vec, float dir_scale, float force_scale, float cage_tol,
                                 float object_density, float object_friction, float num_timesteps)
  : gravity_vec_(gravity_vec),
    num_escapes_(0),
    num_raw_escapes_(0),
    state_(START_TRIAL),
    object_(object),
    obstacles_(obstacles),
    obstacle_radius_(0.0f),
    dir_perturb_scale_(dir_scale),
    force_sigma_(force_scale),
    torque_sigma_(force_scale),
    force_norm_(0.0f, force_sigma_),
    torque_norm_(0.0f, torque_sigma_),
    force_distribution_(rng_, force_norm_),
    torque_distribution_(rng_, torque_norm_),
    potential_thresh_(FLT_MAX),
    cage_tol_(cage_tol),
    object_density_(object_density),
    object_friction_(object_friction),
    num_timesteps_(num_timesteps),
    cur_iter_(0),
    num_iters_(num_iters),
    finished_(false),
    trial_in_progress_(false)
{
  SetGlobalVars();
  InitBox2D();
}



// Set some class vars to global variables (to change once config set up)
void ConvergenceTest::SetGlobalVars()
{
  cage_rates_.resize(num_iters_);
  raw_cage_rates_.resize(num_iters_);
  // force_distribution_ = std::normal_distribution<float>(0.0f, force_sigma_);
  // torque_distribution_ = std::normal_distribution<float>(0.0f, torque_sigma_);
}

void ConvergenceTest::InitBox2D()
{
  LOG(INFO) << "Using potential thresh " << potential_thresh_;

  // center vectors
  box2d_center_pos_.Set(0, 20);

  // set gravity
  m_world->SetGravity(gravity_vec_);

  // TODO replace object creation with box2d shape factory stuff

  // get object initial pose
  std::vector<Eigen::Matrix4f> object_poses_3d = object_->PosesWorldFrame();    
  float object_tx, object_ty, object_theta;
  Pose3DToParams2D(object_poses_3d[0], object_tx, object_ty, object_theta);

  object_start_pos_.Set(object_ty, object_tx);
  object_start_pos_ = object_start_pos_ + box2d_center_pos_;
  object_start_angle_ = object_theta;

  // create dynamic circle
  b2BodyDef body_def;
  body_def.type = b2_dynamicBody;
  body_def.position.Set(object_tx, object_ty);
  body_def.angle = object_theta;
  object_body_ = m_world->CreateBody(&body_def);

  // add vertices from object, obstacles
  b2FixtureDef object_fixture_def;
  b2Separator object_sep;
  std::vector<std::vector<fcl::Vec3f> > object_vertices_3d = object_->BoundaryVertices();
  std::vector<b2Vec2> object_vertices_2d;
  if (object_vertices_3d.size() == 0) {
    LOG(ERROR) << "Error: No object vertices";;
    return;
  }

  b2Vec2 object_com(0.0f, 0.0f);
  int num_v = 0;
  for (unsigned int i = 0; i < object_vertices_3d[0].size(); i++) {
    fcl::Vec3f v = object_vertices_3d[0][i];;

    // only add extruded fronts
    if (v[2] > 0) {
      object_vertices_2d.push_back(b2Vec2(v[0], v[1]));
      object_com += b2Vec2(v[0], v[1]);
      num_v++;
    }
  }
  object_com = (1.0f / num_v) * object_com;

  // convex decomposition of vertices for simulation
  float scale = 1.0f; // some random param for separator that seems to work
  object_sep.Separate(object_body_, &object_fixture_def, &object_vertices_2d, scale);

  // set circle density and friction
  object_fixture_def.density = object_density_;
  object_fixture_def.friction = object_friction_;

  // add the shape to the body.
  object_body_->CreateFixture(&object_fixture_def);

  float obstacle_tx, obstacle_ty, obstacle_theta; 
  b2Separator obstacle_sep;
  obstacle_radius_ = 0.0f;

  // create obstacles
  unsigned int obstacle_index = 0;
  for (unsigned int i = 0; i < obstacles_.size(); i++) {

    std::vector<std::vector<fcl::Vec3f> > obstacle_vertices_3d = obstacles_[i]->BoundaryVertices();
    std::vector<Eigen::Matrix4f> obstacle_poses_3d = obstacles_[i]->PosesWorldFrame();

    for (unsigned int j = 0; j < obstacles_[i]->NumComponents(); j++) {

      // setup body defs
      Pose3DToParams2D(obstacle_poses_3d[j], obstacle_tx, obstacle_ty, obstacle_theta);
      b2BodyDef obstacle_body_def;
      b2FixtureDef obstacle_fixture_def;
      b2Separator obstacle_sep;
      obstacle_body_def.type = b2_staticBody;
      obstacle_body_def.position.Set(obstacle_ty, obstacle_tx); // x and y flipped in box2d
      obstacle_body_def.position += box2d_center_pos_;
      obstacle_body_def.angle = obstacle_theta;

      obstacle_bodies_.push_back(m_world->CreateBody(&obstacle_body_def));

      // update max radius
      float dist_from_start = sqrt((obstacle_tx - object_start_pos_.x) * (obstacle_tx - object_start_pos_.x) + 
                                   (obstacle_ty - object_start_pos_.y) * (obstacle_ty - object_start_pos_.y));
      if (dist_from_start + 2*object_->MaxMomentArm()  > obstacle_radius_) {
        obstacle_radius_ = dist_from_start + 2*object_->MaxMomentArm();
      }

      // add vertices from object, obstacles
      std::vector<b2Vec2> obstacle_vertices_2d;
      for (unsigned int k = 0; k < obstacle_vertices_3d[j].size(); k++) {
        fcl::Vec3f v = obstacle_vertices_3d[j][k];
        // only add extruded fronts
        if (v[2] > 0) {
          obstacle_vertices_2d.push_back(b2Vec2(v[0], v[1]));
        }
      }

      // convex decomposition of vertices for simulation
      obstacle_sep.Separate(obstacle_bodies_[obstacle_index], &obstacle_fixture_def, &obstacle_vertices_2d, scale);

      // set circle density and friction
      obstacle_fixture_def.density = Box2DShapeFactory::default_density;
      obstacle_fixture_def.friction = Box2DShapeFactory::default_friction;

      // add the shape to the body.
      obstacle_bodies_[obstacle_index]->CreateFixture(&obstacle_fixture_def);
      obstacle_index++;
    }
  }

  // set object mass at true center (b2 separator messes it up badly)
  b2MassData mass_data;
  object_body_->GetMassData(&mass_data);
  mass_data.center = object_com;
  object_body_->SetMassData(&mass_data);
  object_mass_ = mass_data.mass;
  object_inertia_ = mass_data.I;

  // convert from euclidean and rotational acceration sigma to actual force sigma
  force_sigma_ = force_sigma_ * mass_data.mass;
  torque_sigma_ = torque_sigma_ * mass_data.I;
  force_norm_ = boost::normal_distribution<float>(0.0f, force_sigma_);
  torque_norm_ = boost::normal_distribution<float>(0.0f, torque_sigma_);
  force_distribution_ = boost::variate_generator<boost::mt19937, boost::normal_distribution<float> >(rng_, force_norm_);
  torque_distribution_ = boost::variate_generator<boost::mt19937, boost::normal_distribution<float> >(rng_, torque_norm_); 
}

b2Vec2 ConvergenceTest::RandomForce()
{
  b2Vec2 random_force(0.0f, 0.0f);
  random_force.x = force_distribution_();
  random_force.x += prev_force_perturb_.x;
  random_force.y = force_distribution_();
  random_force.y += prev_force_perturb_.y;
  prev_force_perturb_ = random_force;
  // std::cout << "Force " << random_force.x << " " << random_force.y << std::endl;
  return random_force;
}

float ConvergenceTest::RandomTorque()
{
  float random_torque = torque_distribution_();
  random_torque += prev_torque_perturb_;
  prev_torque_perturb_ = random_torque;
  // std::cout << "Torque " << random_torque << std::endl;
  return random_torque;
}

void ConvergenceTest::Step(Settings* settings)
{
  Test::Step(settings);

  if (state_ == START_TRIAL) {
    object_body_->SetTransform(object_start_pos_, object_start_angle_); // reset circle position
    object_body_->SetLinearVelocity(b2Vec2(0.0f, 0.0f)); // reset circle position
    object_body_->SetAngularVelocity(0.0f); // reset circle position

    max_potential_ = 0.0f;
    prev_force_perturb_ = b2Vec2(0.0f, 0.0f);
    prev_torque_perturb_ = 0.0f;
    // choose initial direction randomly
    // object_cur_angle_ = 2.0f * M_PI * ((float)rand() / RAND_MAX);
    // object_cur_dir_.x = cos(object_cur_angle_);
    // object_cur_dir_.y = sin(object_cur_angle_);
    // object_cur_torque_ = 2.0f * torque_scale_ * (((float)rand() / RAND_MAX) - 0.5f);

    cur_timestep_ = 1;
    trial_in_progress_ = true;
    state_ = CUR_TRIAL;
  }
  else if (state_ == CUR_TRIAL) {
    // calculate new force (must be smooth)
    // float32 perturb_angle = dir_perturb_scale_ * (((float)rand() / RAND_MAX) - 0.5f);
    // object_cur_angle_ = object_cur_angle_ + perturb_angle;
    // object_cur_dir_.x = cos(object_cur_angle_);
    // object_cur_dir_.y = sin(object_cur_angle_);

    // check whether or not the current position exceeds the energy potential
    b2Vec2 object_pos = object_body_->GetPosition();
    float object_angle = object_body_->GetAngle();

    // TODO:: change to use a single potential
    float cur_gravity_potential = GravityPotential::GRAVITY_ACCEL * object_mass_ * (object_pos.y - object_start_pos_.y);
    max_potential_ = std::max<float>(cur_gravity_potential, max_potential_);

    // calculate new force, torque
    b2Vec2 perturb_force = RandomForce();
    float32 perturb_torque = RandomTorque();

    // apply forces, torques
    object_body_->ApplyForceToCenter(perturb_force, true);
    object_body_->ApplyTorque(perturb_torque, true);

    float32 dist_from_start = (object_pos.x - object_start_pos_.x) * (object_pos.x - object_start_pos_.x) + \
      (object_pos.y - object_start_pos_.y) * (object_pos.y - object_start_pos_.y);

    // check for cage (also works when square object)
    if (dist_from_start > obstacle_radius_ * obstacle_radius_ + cage_tol_) {
      if (max_potential_ <= potential_thresh_) {
        num_escapes_++;
      }
      num_raw_escapes_++;
      trial_in_progress_ = false;
    }

    // check for timesteps over
    if (cur_timestep_ > num_timesteps_) {
      trial_in_progress_ = false;
    }
    cur_timestep_++;
      
    // trial over, state transition
    if (!trial_in_progress_) {

      cage_rates_[cur_iter_] = 1.0f - (float)num_escapes_ / (float)(cur_iter_+1);
      raw_cage_rates_[cur_iter_] = 1.0f - (float)num_raw_escapes_ / (float)(cur_iter_+1);
      if (cur_iter_ % 100 == 0)
        LOG(INFO) << "Trial " << cur_iter_ << " complete. Thresholded cage rate = " << cage_rates_[cur_iter_] << ". Raw Cage Rate = " << raw_cage_rates_[cur_iter_];
      cur_iter_++; // increase trial count

      if (cur_iter_ >= num_iters_) {
        state_ = TRIALS_COMPLETE;
      }
      else {
        state_ = START_TRIAL; // start next trial
      }
    } 
  }
  else if (state_ == TRIALS_COMPLETE) {
    finished_ = true;
  }
}

bool ConvergenceTest::IsFinished()
{
  return finished_;
}

std::vector<float> ConvergenceTest::CageRates()
{
  return cage_rates_;
}

std::vector<float> ConvergenceTest::RawCageRates()
{
  return raw_cage_rates_;
}

int ConvergenceTest::CheckValid(std::vector<b2Vec2> vertices)
{
  b2Separator sep;
  return sep.Validate(vertices);
}
