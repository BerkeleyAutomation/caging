#include "Filtration.h"

#include <glog/logging.h>

Filtration_3::Filtration_3(Potential potential)
  : potential_(potential)
{
}

void Filtration_3::set_potential(Potential potential)
{
  potential_ = potential;
  update_classification();
}

Potential Filtration_3::get_potential()
{
  return potential_;
}

Potential Filtration_3::simplex_potential(Index i)
{
  CGAL_assertion(i >= 0 && i < simplices_.size());
  return simplices_[i+cutoff_index_].potential();
}

void Filtration_3::update_classification()
{
  Index i;
  bool found_cutoff = false;
  for (i = 0; i < simplices_.size() && !found_cutoff; i++) {
    // LOG(INFO) << "Simplex " << i << " " << simplices_[i].potential() << " cutoff " << potential_;
    if (simplices_[i].potential() >= potential_) {
      found_cutoff = true;
    }
  }
  cutoff_index_ = i-1;
  //  LOG(INFO) << "UPDATED " << cutoff_index_;
}

ASE::Filtration_Cell_Class Filtration_3::classify(Cell_handle ch, bool print)
{
  if (indices_.find(ch) == indices_.end()) {
    if (print)
      LOG(INFO) << "invalid";
    return ASE::INVALID;
  }
  if (indices_[ch] < cutoff_index_) {
    return ASE::INTERIOR;
  }
  return ASE::EXTERIOR;
}

void Filtration_3::subcomplex(std::vector<Cell_handle>& complex)
{
  subcomplex(complex, potential_);
}

void Filtration_3::subcomplex(std::vector<Cell_handle>& complex, Potential potential)
{
  complex.clear();
  Cell_handle ch;
  for (unsigned int i = 0; i < cutoff_index_; i++) {
    if (CGAL::assign(ch, simplices_[i].object())) {
      complex.push_back(ch);
    }
  }
}

void Filtration_3::subcomplex_exterior(std::vector<Cell_handle>& complex)
{
  subcomplex_exterior(complex, potential_);
}

void Filtration_3::subcomplex_exterior(std::vector<Cell_handle>& complex, Potential potential, bool print)
{
  complex.clear();
  Cell_handle ch;
  for (unsigned int i = cutoff_index_; i < simplices_.size(); i++) {
    if (CGAL::assign(ch, simplices_[i].object())) {
      complex.push_back(ch);
    }
  }
}

void Filtration_3::render_subcomplex(CGAL::Geomview_stream& gv, float potential_thresh)
{
  Cell_handle ch;
  unsigned int count = 0;
  for (Index i = 0; i < cutoff_index_ && count < 100; i++) {
    if (simplices_[i].potential() > potential_thresh && CGAL::assign(ch, simplices_[i].object())) {
      gv << Convert_Cell_To_Tetrahedron(ch);
      count++;
    }
  }
}

void Filtration_3::render_subcomplex_exterior(CGAL::Geomview_stream& gv)
{
  Cell_handle ch;
  unsigned int count = 0;
  for (Index i = simplices_.size()-1; i >= 0 && count < 1000; i--) {
    if (CGAL::assign(ch, simplices_[i].object())) {
      gv << Convert_Cell_To_Tetrahedron(ch);
      count++;
    }
  }
}

bool Filtration_3::create_boundary_matrix(const Triangulation_3& triangulation)
{
  boundary_matrix_.set_num_cols(simplices_.size());

  for (Index i = 0; i < simplices_.size(); i++) {
    // clear column of boundary matrix
    boundary_matrix_.clear(i);
    phat::column col;
    ASE::Simplex simplex = simplices_[i];

    // handle vertex
    if (simplex.object().is<Vertex_handle>()) {
      boundary_matrix_.set_dim(i, 0);

      Vertex_handle v = CGAL::object_cast<Vertex_handle>(simplex.object());
      vertex_info_map_[v].set_index(i);
      boundary_matrix_.set_col(i, col);
    }
    // handle edge
    else if (simplex.object().is<Edge>()) {
      boundary_matrix_.set_dim(i, 1);

      Edge e = CGAL::object_cast<Edge>(simplex.object());
      Vertex_handle v1 = e.first->vertex(e.second);
      Vertex_handle v2 = e.first->vertex(e.third);
      
      Index i1 = vertex_info_map_[v1].index();
      Index i2 = vertex_info_map_[v2].index();

      if (i1 > i2) {
        std::swap(v1, v2);
        std::swap(i1, i2);
      }

      if (i <= i1 || i <= i2) {
        LOG(INFO) << "ERROR";
        LOG(INFO) << "Edge index: " << i << " referencing vertices " << i1 << " and " << i2;
      }
      col.push_back(i1);
      col.push_back(i2);
      boundary_matrix_.set_col(i, col);
      set_index_of_edge(triangulation, e, i);
    }
    // handle face
    else if (simplex.object().is<Facet>()) {
      boundary_matrix_.set_dim(i, 2);
      
      Facet f = CGAL::object_cast<Facet>(simplex.object());
      Index i1 = cell_info_map_[f.first].edge_index( (f.second+1)%4, (f.second+2)%4 );
      col.push_back(i1);
      Index i2 = cell_info_map_[f.first].edge_index( (f.second+1)%4, (f.second+3)%4 );
      col.push_back(i2);
      Index i3 = cell_info_map_[f.first].edge_index( (f.second+2)%4, (f.second+3)%4 );
      col.push_back(i3);

      if (i <= i1 || i <= i2 || i <= i3) {
        LOG(INFO) << "ERROR";
        LOG(INFO) << "Facet index: " << i << " referencing edges " << i1 << ", " << i2 << " and " << i3;
      }

      std::sort(col.begin(),col.end());
      boundary_matrix_.set_col(i, col);
      set_index_of_facet(triangulation, f, i);
    }
    // handle cell
    else if (simplex.object().is<Cell_handle>()) {
      boundary_matrix_.set_dim(i, 3);
      
      Cell_handle c = CGAL::object_cast<Cell_handle>(simplex.object());
      if (!cell_info_map_[c].has_facet_index(0)) {
        std::cout << i << " " << 0 << std::endl; 
      }
      else if (!cell_info_map_[c].has_facet_index(1)) {
        std::cout << i << " " << 1 << std::endl; 
      }
      else if (!cell_info_map_[c].has_facet_index(2)) {
        std::cout << i << " " << 2 << std::endl; 
      }
      else if (!cell_info_map_[c].has_facet_index(3)) {
        std::cout << i << " " << 3 << std::endl; 
      }

      col.push_back(cell_info_map_[c].facet_index(0));
      col.push_back(cell_info_map_[c].facet_index(1));
      col.push_back(cell_info_map_[c].facet_index(2));
      col.push_back(cell_info_map_[c].facet_index(3));
      std::sort(col.begin(), col.end());
      boundary_matrix_.set_col(i, col);
    }
  }
  return true;
}

void Filtration_3::set_index_of_edge(const Triangulation_3& T, const Edge& e, Index I)
{
  Vertex_handle v1 = e.first->vertex(e.second);
  Vertex_handle v2 = e.first->vertex(e.third);

  Cell_circulator ch = T.incident_cells(e);
  Cell_circulator ch_start = ch;
  int count = 0;
  do {
    cell_info_map_[ch].set_edge_index(ch->index(v1), ch->index(v2), I);
    ch++;
    count++;
  } while (ch != ch_start);
}

void Filtration_3::set_index_of_facet(const Triangulation_3& T, const Facet& f, Index I)
{
  cell_info_map_[f.first].set_facet_index(f.second, I);
  Facet mf = T.mirror_facet(f);
  cell_info_map_[mf.first].set_facet_index(mf.second, I);    
}

void Filtration_3::save_boundary_matrix(std::string filename)
{
  boundary_matrix_.save_binary(filename.c_str());
}

Alpha_Shape_Filtration_3::Alpha_Shape_Filtration_3(const Alpha_shape_3& alpha_shape)
  : Filtration_3(0)
{
  init(alpha_shape);
}

void Alpha_Shape_Filtration_3::init(const Alpha_shape_3& alpha_shape)
{
  std::vector<CGAL::Object> filtered_simplices;
  std::vector<K::FT> alpha_values;
  alpha_shape.filtration_with_alpha_values(CGAL::Dispatch_output_iterator<CGAL::cpp11::tuple<CGAL::Object, K::FT>,
                                                                          CGAL::cpp11::tuple<std::back_insert_iterator<std::vector<CGAL::Object> >,
                                                                                              std::back_insert_iterator<std::vector<K::FT> > > >(std::back_inserter(filtered_simplices),
                                                                                                                                                 std::back_inserter(alpha_values))); 
  filter(filtered_simplices, alpha_values);
  Filtration_3::create_boundary_matrix(alpha_shape);
}

// assess the potential of each simplex
void Alpha_Shape_Filtration_3::filter(std::vector<CGAL::Object> unordered_simplices, std::vector<Potential> potentials)
{
  simplices_.clear();

  // create simplex list
  for (Index i = 0; i < unordered_simplices.size(); i++) {
    simplices_.push_back(ASE::Simplex(unordered_simplices[i], potentials[i]));
  }  

  // sort simplices by potential
  std::sort(simplices_.begin(), simplices_.end());

  // update the indices
  Cell_handle ch;
  for (Index i = 0; i < simplices_.size(); i++) {
    if (CGAL::assign(ch, simplices_[i].object())) {
      indices_[ch] = i;
    }    
  }

  // update the classification
  Filtration_3::update_classification();
}

Manifold_Sweep_Filtration_3::Manifold_Sweep_Filtration_3(const Triangulation_3& triangulation, std::vector<CGAL::Object> unordered_simplices, PotentialFunction& p_func)
  : Filtration_3(0)
{
  init(triangulation, unordered_simplices, p_func);
}

void Manifold_Sweep_Filtration_3::init(const Triangulation_3& triangulation, std::vector<CGAL::Object> unordered_simplices, PotentialFunction& p_func)
{
  filter(unordered_simplices, p_func);
  Filtration_3::create_boundary_matrix(triangulation);
}

bool Manifold_Sweep_Filtration_3::filter(std::vector<CGAL::Object> unordered_simplices, PotentialFunction& p_func)
{
  simplices_.clear();

  // create simplex list
  LOG(INFO) << "Filtering from " << unordered_simplices.size() << " simplices";
  for (Index i = 0; i < unordered_simplices.size(); i++) {
    ASE::Simplex s(unordered_simplices[i]);
    Potential p = p_func(s);
    s.set_potential(p);
    simplices_.push_back(s);
  }  

  // sort simplices by potential
  LOG(INFO) << "Sorting simplices by potential";
  std::sort(simplices_.begin(), simplices_.end());

  // update the indices
  LOG(INFO) << "Setting index map";
  Cell_handle ch;
  for (Index i = 0; i < simplices_.size(); i++) {
    if (CGAL::assign(ch, simplices_[i].object())) {
      indices_[ch] = i;
    }    
  }

  // update the classification
  LOG(INFO) << "Updating simplex classification";
  Filtration_3::update_classification();
}

