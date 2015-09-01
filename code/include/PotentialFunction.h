#pragma once

#include <map>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Union_find.h>
#include <CGAL/Timer.h>
#include <CGAL/utility.h>
#include <CGAL/IO/Geomview_stream.h>

#include "Triangulation_Planner.h"

typedef double Potential;

// temporary ASE lab simplex "classes"
namespace ASE {
  typedef Weighted_point Vertex;

  class Simplex {
   public:
    Simplex(CGAL::Object object, Potential potential = 0) {
      object_ = object;
      potential_ = potential;

      Edge e;
      Facet f;
      Cell_handle c;

      // assign dimension
      dim_ = 0;
      if (CGAL::assign(e, object_)) {
        dim_ = 1;
      }
      else if (CGAL::assign(f, object_)) {
        dim_ = 2;
      }
      else if (CGAL::assign(c, object_)) {
        dim_ = 3;
      }
    }

    unsigned int dim() const {
      return dim_;
    }

    Vertex vertex(unsigned int index) {
      Vertex_handle v;
      Edge e;
      Facet f;
      Cell_handle c;
      
      if (CGAL::assign(v, object_)) {
        return v->point();
      }
      else if (CGAL::assign(e, object_)) {
        if (index == 0) {
          return e.first->vertex(e.second)->point();
        }
        return e.first->vertex(e.third)->point();
      }
      else if (CGAL::assign(f, object_)) {
        unsigned int cell_index = index;
        if (index >= f.second) {
          cell_index = index+1;
        }
        return f.first->vertex(cell_index)->point();
      }
      else if (CGAL::assign(c, object_)) {
        return c->vertex(index)->point();
      }
      else {
        CGAL_assertion(false);
      }
    }

    CGAL::Object object() {
      return object_;
    }

    Potential potential() const {
      return potential_;
    }

    void set_potential(Potential p) {
      potential_ = p;
    }

    bool operator<(const Simplex& other) const {
      // for equality order by dimension
      if (potential_ == other.potential()) {
        return (dim_ < other.dim());
      }
      return (potential_ < other.potential());
    }
    
   private:
    unsigned int dim_;
    /* std::vector<Vertex> vertices_; */
    CGAL::Object object_; // not super-happy about this
    Potential potential_;
  };
};

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

class GravityPotential : public PotentialFunction
{
 public:
  static double GRAVITY_ACCEL; // in meters per second

 public:
  GravityPotential(ASE::Vertex reference_pose, double mass,
                   float x_scale, float y_scale, float theta_scale);
  ~GravityPotential() {}

 public:
  Potential evaluate(ASE::Simplex simplex);
  Potential potential(ASE::Vertex vertex);

 private:
  ASE::Vertex reference_pose_;
  double mass_;
  float x_scale_;
  float y_scale_;
  float theta_scale_;
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
