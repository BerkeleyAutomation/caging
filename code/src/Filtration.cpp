#include "Filtration.h"
#include "Configuration_Space_Approximator.h"
#include "Util.h"

#include <glog/logging.h>
#include "phat/compute_persistence_pairs.h"

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

bool Filtration_3::simplex(Cell_handle ch, ASE::Simplex& s)
{
  if (indices_.find(ch) == indices_.end()) {
    LOG(ERROR) << "Could not find cell handle";
    return false;
  }
  s = simplices_[indices_[ch]];
  return true;
}

void Filtration_3::update_classification()
{
  Index i;
  bool found_cutoff = false;
  for (i = 0; i < simplices_.size() && !found_cutoff; i++) {
    // LOG(INFO) << "Simplex " << i << " " << simplices_[i].potential() << " cutoff " << potential_;
    if (simplices_[i].potential() > potential_) {
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


ASE::Filtration_Cell_Class Filtration_3::classify(Cell_handle ch1, Cell_handle ch2, bool print)
{
  if (indices_.find(ch1) == indices_.end() || indices_.find(ch2) == indices_.end()) {
    if (print)
      LOG(INFO) << "invalid";
    return ASE::INVALID;
  }

  int cur_index;
  bool match_found = false;
  for (int i = 0; i < 4 && !match_found; i++) {
    cur_index = cell_info_map_[ch1].facet_index(i);
    for (int j = 0; j < 4 && !match_found; j++) {
      int other_index = cell_info_map_[ch2].facet_index(j);
      if (cur_index == other_index) {
       match_found = true;
      }
    }
  }

  if (cur_index < cutoff_index_) {
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

void Filtration_3::subcomplex_exterior_simplices(std::vector<ASE::Simplex>& complex)
{
  complex.clear();
  Cell_handle ch;
  for (unsigned int i = cutoff_index_; i < simplices_.size(); i++) {
    if (CGAL::assign(ch, simplices_[i].object())) {
      complex.push_back(simplices_[i]);
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

//Synthesis method for deciding the bext simplex in which to place the object.
std::vector< Persistence_info > Filtration_3::run_persistence(std::vector<Weighted_point> points_at_infinity,
                                                              Potential cutoff,
                                                              float theta_scale,
                                                              ASE::Vertex& birth_point,
                                                              float min_energy_thresh)
{
  //Save pair data out to a scv!
  std::ofstream of("persistance_pairs.csv", std::ios::app);  	
  
  set_potential(cutoff);
  
  // define the object to hold the resulting persistence pairs
  phat::persistence_pairs pairs;

  // choose an algorithm (choice affects performance) and compute the persistence pair
  // (modifies boundary_matrix)
  phat::compute_persistence_pairs< phat::twist_reduction >( pairs, boundary_matrix_ );
 
  //Vector of pairs of energy deltas and the corresponding death simplex.
  std::vector< Persistence_info > delta_pairs;
      
  for (phat::index i = 0; i < pairs.get_num_pairs(); i++)
  {
    phat::index birth_index = pairs.get_pair(i).first;
    phat::index death_index = pairs.get_pair(i).second;
    Potential birth_energy = simplices_[birth_index].potential();
    Potential death_energy = simplices_[death_index].potential();
    Potential energy_delta = abs(birth_energy - death_energy);
    ASE::Simplex birth_simplex = simplices_[birth_index];
    ASE::Simplex death_simplex = simplices_[death_index];
          
    if ( (birth_index < cutoff_index_ && death_index < cutoff_index_) 
         || energy_delta < min_energy_thresh
         || touches_infinity(simplices_[birth_index], points_at_infinity) || touches_infinity(simplices_[death_index], points_at_infinity) 
         || simplices_[death_index].dim() != 3 )
    {
      continue;
    }

    if (energy_delta == 0) {
      continue;
    }

    Persistence_info temp_info(energy_delta, birth_simplex, death_simplex, birth_index, death_index);
    delta_pairs.push_back(temp_info);  

    of << birth_index << ", " << death_index <<  "\n";
  
    ////DEBUG Sanity check that the birth simplex is two dimensional....
    //if (simplices_[birth_index].dim() != 2)
    //{
    //  std::cout << "SHADY SHIT!" << std::endl;
    //  std::cout << "ENERGY LEVEL" << energy_delta << std::endl;
    //}

    //DEBUG Check potential
    //if (energy_delta > 50.0){
    //  std::cout << "FILTRATION POTENTIAL: " << energy_delta << std::endl;
    //  std::cout << "DEATH INDEX " << birth_index << std::endl;
    //  std::cout << "BIRTH INDEX " << death_index << std::endl;
    //  std::cout << "CUTOFF INDEX " << cutoff_index_ << std::endl;
    //}
  }

  //while(1){}

  //We want the death index (simplex) corresponding to the the larget energy delta!
  std::sort(delta_pairs.begin(), delta_pairs.end());

 
  ////Just some good 'ol debug to print the largest energy_delta and other information....
  //for (unsigned int i = delta_pairs.size()-1; i >= delta_pairs.size()-20; i--)
  //  std::cout << "Delta for " << delta_pairs.size()-1-i << "th highest energy: " << delta_pairs[i].get_potential() << std::endl;    
  //std::cout << "Largest Delta: " << delta_pairs.back().get_potential() << std::endl;
  //ASE::Simplex birth_simplex = delta_pairs.back().get_birth_simplex();
  //ASE::Simplex death_simplex = delta_pairs.back().get_death_simplex();
 

  //ASE::Vertex death_vertex;
  //Potential highest_potential = -FLT_MAX; 
  //std::cout << "Type: " << death_simplex.dim() << std::endl;
  //for (int i = 0; i <= death_simplex.dim(); i++) 
  //{
  //  ASE::Vertex cur_vertex = death_simplex.vertex(i); 
  //  Potential cur_potential = free_potential_func.potential(cur_vertex);
  //  if (cur_potential > highest_potential)
  //  {
  //    highest_potential = cur_potential;
  //    death_vertex = cur_vertex;
  //  }
  //} 

  //highest_potential = -FLT_MAX; 
  //for (int i = 0; i <= birth_simplex.dim(); i++) 
  //{
  //  ASE::Vertex cur_vertex = birth_simplex.vertex(i); 
  //  Potential cur_potential = free_potential_func.potential(cur_vertex);
  //  if (cur_potential > highest_potential)
  //  {
  //    highest_potential = cur_potential;
  //    birth_point = cur_vertex;
  //  }
  //}

  ////Point_3 centroid = death_simplex.centroid();
  ////ASE::Vertex best_pose(Bare_point((float) CGAL::to_double(centroid.x()),
  //  //                               (float) CGAL::to_double(centroid.y()),
  //    //                             fmod(((float) CGAL::to_double(centroid.z()) / theta_scale),(2*M_PI))),
  //      //                1.0);
  ////std::cout << "best pose: " << best_pose << std::endl;
  ////std::cout << "highest_potential: " << highest_potential << std::endl;  

  of.close();

  return delta_pairs;   
}

bool Filtration_3::touches_infinity(ASE::Simplex given_simplex, std::vector<Weighted_point> points_at_infinity)
{
  for (int i = 0; i <= given_simplex.dim(); i++)
  {
    ASE::Vertex cur_vertex = given_simplex.vertex(i); 
    for (int j = 0; j < points_at_infinity.size(); j++)
    {   
      if (points_at_infinity[j].point() == cur_vertex.point())
      {
        return true;
      } 
    }  
  }
  return false; 
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

