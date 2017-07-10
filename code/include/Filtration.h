#ifndef FILTRATION_H
#define FILTRATION_H

#include "PotentialFunction.h"
#include "Typedef.h"

namespace ASE {
  enum Filtration_Cell_Class {
    INTERIOR, // interior to current subcomplex
    EXTERIOR, // exterior to current subcomplex
    INVALID   // not part of the subcomplex or not a cell 
  };
}

// vertex and cell info borrowed from phat/addons/alpha_3.cpp
struct Vertex_info_3 {
  Vertex_info_3() {
    index_ = -1;
  }
  
  bool has_index() {
    return (index_ >= 0);
  }
  
  Index index() {
    CGAL_assertion(has_index());
    return index_;
  }

  void set_index(Index I) {
    index_ = I;
  }
  
private: 
  Index index_;
};


//Struct that is used by the run_persistance method to return important information on the best persisence pair.
struct Persistence_info {
  Persistence_info() { 
  } 

  Persistence_info(Potential potential, ASE::Simplex birth_simplex, ASE::Simplex death_simplex, int birth_index, int death_index) {
    potential_ = potential;
    birth_simplex_ = birth_simplex;
    death_simplex_ = death_simplex; 
    birth_index_ = birth_index;
    death_index_ = death_index;
  }
  
  Potential get_potential() const {
    return potential_;
  }
  
  ASE::Simplex get_birth_simplex() {
    return birth_simplex_;
  }
  
  ASE::Simplex get_death_simplex() {
    return death_simplex_;
  }

  bool operator<(const Persistence_info& other) const {
    return potential_ < other.get_potential();
  } 

public: 
  Potential potential_;
  ASE::Simplex birth_simplex_;
  ASE::Simplex death_simplex_;
  int birth_index_;
  int death_index_;
};


// Struct to store info about the indices of cells for 3D triangulations. Necessary because CGAL is very lean
struct Cell_info_3 {
  Cell_info_3() {
    for(std::size_t i = 0; i < 6; i++) {
      edge_index_[i] = -1;
    }
    for(std::size_t i = 0; i < 4; i++) {
      facet_index_[i] = -1;
    }
  }

  int edge_conv(int i, int j) {
    if(i>j) std::swap(i,j);
    if(i==0 && j==1) return 0;
    if(i==0 && j==2) return 1;
    if(i==0 && j==3) return 2;
    if(i==1 && j==2) return 3;
    if(i==1 && j==3) return 4;
    if(i==2 && j==3) return 5;
  }


  bool has_edge_index(int i, int j) {
    return (edge_index_[edge_conv(i,j)] >= 0);
  }

  Index edge_index(int i, int j) {
    CGAL_assertion(has_edge_index(i,j));
    int k = edge_conv(i,j);
    return edge_index_[k];
  }

  bool has_facet_index(int i) {
    CGAL_assertion(i>=0 && i<4);
    return (facet_index_[i] >= 0);
  }

  Index facet_index(int i) {
    CGAL_assertion(has_facet_index(i));
    return facet_index_[i];
  }

  void set_edge_index(int i, int j, Index I) {
    edge_index_[edge_conv(i,j)] = I;
  }

  void set_facet_index(int i, Index I) {
    facet_index_[i] = I;
  }

private:
  // stores the indices of edges and stuff
  Index edge_index_[6];
  Index facet_index_[4];
};

class Filtration_3
{
 protected:
  Filtration_3(Potential potential = 0);

 public:
  // Puts the subcomplex below potential into the vector | complex |
  void subcomplex(std::vector<Cell_handle>& complex);
  void subcomplex(std::vector<Cell_handle>& complex, Potential potential);
  // Puts the simplices exterior to the subcomplex below potential into the vector | complex |
  void subcomplex_exterior(std::vector<Cell_handle>& complex);
  void subcomplex_exterior(std::vector<Cell_handle>& complex, Potential potential, bool print = false);
  void subcomplex_exterior_simplices(std::vector<ASE::Simplex>& complex);
  // Sets the potential in the filtration
  void set_potential(Potential potential);
  Potential get_potential();
  Potential simplex_potential(Index i);
  bool simplex(Cell_handle ch, ASE::Simplex& s);

  // Classifies the cell as being inside / outside the current filtered complex
  ASE::Filtration_Cell_Class classify(Cell_handle ch, bool print = false); 
  ASE::Filtration_Cell_Class classify(Cell_handle ch1,Cell_handle ch2, bool print = false);

  // Saves the current boundary matrix
  void save_boundary_matrix(std::string filename = "filtration.bin");

  // Render the triangles of the subcomplex to the geomview stream
  void render_subcomplex(CGAL::Geomview_stream& gv, float potential_thresh);
  void render_subcomplex_exterior(CGAL::Geomview_stream& gv);
  
  // Return the ideal simplex to place the object through Persistent Homology
  std::vector< Persistence_info > run_persistence(std::vector<Weighted_point> points_at_infinity, Potential cutoff,
                                                  float theta_scale, ASE::Vertex& birth_point,
                                                  float min_energy_thresh = 0.0f);
  bool touches_infinity(ASE::Simplex given_simplex, std::vector<Weighted_point> points_at_infinity);

 protected:
  // Updates the classification cutoff index based on the current value
  void update_classification();
  // compute a boundary matrix
  bool create_boundary_matrix(const Triangulation_3& triangulation);

  // helper functions
  void set_index_of_edge(const Triangulation_3& T, const Edge& e, Index I);
  void set_index_of_facet(const Triangulation_3& T, const Facet& f, Index I);

 protected:
  std::vector<ASE::Simplex> simplices_; // list of simplices
  std::map<Cell_handle, Index> indices_; // indices of cells
  Potential potential_; // current potential threshold for filtering
  Index cutoff_index_; // index under which all values are in the complex and all others are not

  std::map<Vertex_handle, Vertex_info_3> vertex_info_map_; // maps vertices to indices for boundary matrix
  std::map<Cell_handle, Cell_info_3> cell_info_map_;       // maps cell handles to edge and face indices for boundary matrix
  phat::boundary_matrix<phat::vector_vector> boundary_matrix_; // boundary matrix for persistent homology
};

// class to encapsulate filtrations generated by CGAL (e.g. alpha shapes)
class Alpha_Shape_Filtration_3 : public Filtration_3
{
 public:
  Alpha_Shape_Filtration_3(const Alpha_shape_3& alpha_shape);
  ~Alpha_Shape_Filtration_3() {}
  
 private:
  // get the filtration and convert to simplices
  void init(const Alpha_shape_3& alpha_shape);
  // assess the potential of each simplex and sort
  void filter(std::vector<CGAL::Object> unordered_simplices, std::vector<Potential> potentials);
};

// class to encapsulate filtrations generated by CGAL (e.g. alpha shapes)
class Manifold_Sweep_Filtration_3 : public Filtration_3
{
 public:
  Manifold_Sweep_Filtration_3(const Triangulation_3& triangulation, std::vector<CGAL::Object> unordered_simplices, PotentialFunction& p_func);
  ~Manifold_Sweep_Filtration_3() {}
  
 private:
  // get the filtration and convert to simplices
  void init(const Triangulation_3& triangulation, std::vector<CGAL::Object> unordered_simplices, PotentialFunction& p_func);
  // compute a boundary matrix given a list of in order and their values (assumed to be in order)
  bool filter(std::vector<CGAL::Object> unordered_simplices, PotentialFunction& p_func);
};

#endif // FILTRATION_H
