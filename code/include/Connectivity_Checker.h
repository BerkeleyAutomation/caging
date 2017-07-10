#pragma once

#include "Filtration.h"
#include "Typedef.h"

////// ********************    CLASS DEFINITION    ********************/

// this class will allow querying of two points to test 
// whether or not they are disconnected
// it uses a filtration to determine which areas of the space are inside / outside of the shape
class Connectivity_Checker
{
 public:
  // construct the disjoint set that contains the connectivity information
  Connectivity_Checker(Triangulation_3& raw_triangulation, Filtration_3& filtration,
                       CGAL::Geomview_stream& gv, Potential potential_thresh = 0.0f,
                       bool cache_connected_components = true);

  // query for connections / disconnections between two points
  bool connected(Bare_point qs, Bare_point qg);
  bool connected(Bare_point qs, Bare_point qg, Potential potential);
  bool disconnected(Bare_point qs, Bare_point qg);
  bool disconnected(Bare_point qs, Bare_point qg, Potential potential);
  bool set_goal_start_cells(Bare_point qs, Bare_point qg, Potential potential, 
                            std::vector<Cell_handle>& start_cells,
  	                    std::vector<Cell_handle>& goal_cells);
  
  //DFS algorithm to show paths between points...			    
  bool render_cc_path(Bare_point qs, Bare_point qg, Potential potential);
  //Recursive DFS implementation
  bool Cell_DFS(Cell_handle cur_cell, std::map<Cell_handle, Cell_handle>& closed_set, 
  Cell_handle& found_goal, std::vector<Cell_handle>& goal_cells);
  
  // return a triangulation of the connected component containing the query cell
  std::vector<ASE::Simplex> connected_component(Bare_point bp);
  std::vector<ASE::Simplex> connected_component(Bare_point bp, Potential potential);

  // get the volume of the connected component containing the point q
  bool connected_component_volume(Bare_point q, float& volume);

  // set a new potential threshold
  inline Potential get_potential() { return potential_thresh_; } 
  void set_potential(Potential p);

 private:
  // for resetting the internal connectivity structures
  void reset_classification();
  void compute_connected_components();

  // return the number of sets
  int num_connected_components() {
    return disjoint_set_.number_of_sets();
  }

 private:
  Triangulation_3& triangulation_;
  Filtration_3& filtration_;
  Potential potential_thresh_;
  std::vector<Cell_handle> external_cells_;
  bool cache_connected_components_;

  CGAL::Geomview_stream& gv_;

  Timer timer_;
  Disjoint_set disjoint_set_;
  std::map<Cell_handle, handle> map_; //map converts cell_handles from location queries to ds structure handles
};
