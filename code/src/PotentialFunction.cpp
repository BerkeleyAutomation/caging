#include "PotentialFunction.h"
#include "Util.h"

bool compare(Bare_point p0, Bare_point p1) {
  if (fabs(p0.x() - p1.x()) < 1e-2 && 
      fabs(p0.y() - p1.y()) < 1e-2 && 
      fabs(p0.z() - p1.z()) < 1e-2) {
    return true;
  }
  return false;
}

LinearPotentialFunction::LinearPotentialFunction(ASE::Vertex reference_pose, double push_force,
                                                 float x_scale, float y_scale, float theta_scale,
                                                 float theta_offset)
  : reference_pose_(reference_pose),
    push_force_(push_force),
    x_scale_(x_scale),
    y_scale_(y_scale),
    theta_scale_(theta_scale)
{
  Eigen::Vector2d unit_down(-1, 0);
  force_vector_ = rotate_vector(unit_down, theta_offset);
}


Eigen::Vector2d LinearPotentialFunction::get_force_vector()
{
  return force_vector_;
}

Potential LinearPotentialFunction::evaluate(ASE::Simplex simplex)
{
  // take min gravitational potential over the vertices
  Potential min_energy = FLT_MAX;
  Potential p;
  for (unsigned int i = 0; i <= simplex.dim(); i++) {
    p = potential(simplex.vertex(i));
    if (p < min_energy) {
      min_energy = p;
    }
  }
  return min_energy;
}


Potential LinearPotentialFunction::potential(ASE::Vertex vertex)
{
  Eigen::Vector2d vertex_vector(vertex.x(), vertex.y());
  Eigen::Vector2d reference_vector(reference_pose_.x(), reference_pose_.y());
  Eigen::Vector2d difference_vector = reference_vector -  vertex_vector;
  double height = difference_vector.dot(force_vector_);
  return push_force_ * height / y_scale_;
}


DistancePotential::DistancePotential(ASE::Vertex reference_pose,
                                   float x_scale, float y_scale, float theta_scale)
  : reference_pose_(reference_pose),
    x_scale_(x_scale),
    y_scale_(y_scale),
    theta_scale_(theta_scale)
{
}

Potential DistancePotential::evaluate(ASE::Simplex simplex)
{
  // take min gravitational potential over the vertices
  Potential min_energy = FLT_MAX;
  Potential p;
  for (unsigned int i = 0; i <= simplex.dim(); i++) {
    p = potential(simplex.vertex(i));
    if (p < min_energy) {
      min_energy = p;
    }
  }
  return min_energy;
}

Potential DistancePotential::potential(ASE::Vertex vertex)
{
  float x_diff = (vertex.point().x() - reference_pose_.point().x()) / x_scale_;
  float y_diff = (vertex.point().y() - reference_pose_.point().y()) / y_scale_;
  float theta_diff = (vertex.point().z() - reference_pose_.point().z()) / theta_scale_;
  return sqrt(x_diff * x_diff + y_diff * y_diff + theta_diff * theta_diff);
}

AlphaCustomPotential::AlphaCustomPotential(Alpha_shape_3& alpha_shape, PotentialFunction& potential_fn,
                                           float infinity_radius, float alpha_thresh)
  : alpha_shape_(alpha_shape),
    alpha_thresh_(alpha_thresh),
    infinity_radius_(infinity_radius),
    pf_(potential_fn)
{
}

Potential AlphaCustomPotential::evaluate(ASE::Simplex simplex)
{
  Vertex_handle vh;
  Edge e;
  Facet f;
  Cell_handle ch;

  // assign negative infty potential if interior to the alpha shape
  if (CGAL::assign(vh, simplex.object())) {
    if (alpha_shape_.classify(vh, alpha_thresh_) != Alpha_shape_3::EXTERIOR) {
      return -FLT_MAX;
    } 
  }
  else if (CGAL::assign(e, simplex.object())) {
    if (alpha_shape_.classify(e, alpha_thresh_) != Alpha_shape_3::EXTERIOR) {
      return -FLT_MAX;
    } 
  }
  else if (CGAL::assign(f, simplex.object())) {
    if (alpha_shape_.classify(f, alpha_thresh_) != Alpha_shape_3::EXTERIOR) {
      return -FLT_MAX;
    } 
  }
  else if (CGAL::assign(ch, simplex.object())) {
    if (alpha_shape_.classify(ch, alpha_thresh_) != Alpha_shape_3::EXTERIOR) {
      return -FLT_MAX;
    } 
  }
  
  // filter out simplices with vertices at infinity
  for (unsigned int i = 0; i <= simplex.dim(); i++) {
    float rad = simplex.vertex(i).weight();
    if (rad >= infinity_radius_) {
      return FLT_MAX;
    }
  }

  // otherwise return the negative gravity potential
  return -pf_(simplex);
}

Potential AlphaCustomPotential::potential(ASE::Vertex vertex)
{
  return -pf_.potential(vertex);
}
