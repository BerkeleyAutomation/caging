#ifndef SIMPLEX_H
#define SIMPLEX_H

#include "Typedef.h"

namespace ASE {
  typedef Weighted_point Vertex;

  class Simplex {
   public:
    Simplex() {
    }
 
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

    Bare_point centroid() {
      float centroid_x  = 0;
      float centroid_y  = 0;
      float centroid_z  = 0;

      for (int i = 0; i <= dim_; i++) {
        centroid_x = centroid_x + vertex(i).point().x();
        centroid_y = centroid_y + vertex(i).point().y();
        centroid_z = centroid_z + vertex(i).point().z();
      }

      centroid_x = centroid_x/(dim_ + 1);
      centroid_y = centroid_y/(dim_ + 1);
      centroid_z = centroid_z/(dim_ + 1);
       
      return Bare_point(centroid_x, centroid_y, centroid_z); 
    }

    void print_vertices() {
     std::cout << "Dimension: " << dim_ << std::endl; 
     for (int i = 0; i <= dim_; i++) {
        Vertex cur_vertex = vertex(i);
        std::cout << "Vertex " << i << ": " << std::endl;
        std::cout << cur_vertex.point() << std::endl;
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

#endif // SIMPLEX_H
