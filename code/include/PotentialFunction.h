#pragma once

#include "Typedef.h"
#include "Simplex.h"

#define GRAVITY_ACCEL 9.81f

// temporary ASE lab simplex "classes"

class PotentialFunction
{
 public:
  // Evaluate a potential energy function over the simplex defined by the vertices supplied in | vertices |
  virtual Potential evaluate(ASE::Simplex simplex) = 0;
  virtual Potential potential(ASE::Vertex vertex) = 0;
  Potential operator()(ASE::Simplex simplex) {
    return evaluate(simplex);
  }

};

class LinearPotentialFunction : public PotentialFunction
{
 public:
  LinearPotentialFunction(ASE::Vertex reference_pose, double push_force,
                          float x_scale, float y_scale, float theta_scale,
                          float theta_offset = 0.0f);
  ~LinearPotentialFunction() {}

 public:
  Potential evaluate(ASE::Simplex simplex);
  Potential potential(ASE::Vertex vertex);
  Eigen::Vector2d get_force_vector();

 private:
  ASE::Vertex reference_pose_;
  double push_force_;
  float x_scale_;
  float y_scale_;
  float theta_scale_;
  Eigen::Vector2d force_vector_;
};


class DistancePotential : public PotentialFunction
{
 public:
  DistancePotential(ASE::Vertex reference_pose,
                    float x_scale, float y_scale, float theta_scale);
  ~DistancePotential() {}

 public:
  Potential evaluate(ASE::Simplex simplex);
  Potential potential(ASE::Vertex vertex);

 private:
  ASE::Vertex reference_pose_;
  float x_scale_;
  float y_scale_;
  float theta_scale_;
};

// uses the alpha shape with alpha value 0 as the 1st simplex in the filtration with -infty potential, then uses gravity from thereon out
class AlphaCustomPotential : public PotentialFunction
{
 public:
  AlphaCustomPotential(Alpha_shape_3& alpha_shape_, PotentialFunction& potential_fn,
                        float infinity_radius, float alpha_thresh = 0.0f);
  ~AlphaCustomPotential() {}

 public:
  Potential evaluate(ASE::Simplex simplex);
  Potential potential(ASE::Vertex vertex);

 private:
  Alpha_shape_3& alpha_shape_;
  float alpha_thresh_;
  float infinity_radius_;
  PotentialFunction& pf_;  
};
