#include "Connectivity_Checker.h"

#include "Util.h"

#include <glog/logging.h>
#include <queue>

Connectivity_Checker::Connectivity_Checker(Triangulation_3& raw_triangulation,
                                           Filtration_3& filtration,
                                           CGAL::Geomview_stream& gv,
                                           Potential potential_thresh,
                                           bool cache_connected_components)
  : triangulation_(raw_triangulation),
    filtration_(filtration),
    potential_thresh_(potential_thresh),
    cache_connected_components_(cache_connected_components),
    gv_(gv)
{
  set_potential(potential_thresh_);
}

void Connectivity_Checker::set_potential(Potential p)
{
  potential_thresh_ = p;
  filtration_.set_potential(potential_thresh_);
  //LOG(INFO) << "Set filtration potential to " << potential_thresh_;
  if (cache_connected_components_) {
    reset_classification();
    //LOG(INFO) << "Reset classification";
  }
}

void Connectivity_Checker::reset_classification()
{
  external_cells_.clear();
  disjoint_set_.clear();
  map_.clear();
  compute_connected_components();
}

void Connectivity_Checker::compute_connected_components()
{
  // get cells exterior to current subcomplex
  LOG(INFO) << "Getting the simplices exterior to the current filtered complex";
  filtration_.subcomplex_exterior(external_cells_);

  timer_.start();
  // initialize disjoint set structure with all of the finite external cells
  LOG(INFO) << "Computing connected comps";
  Cell_handle current_cell;
  unsigned int j = 0;
  for (std::vector<Cell_handle>::iterator it = external_cells_.begin(); it != external_cells_.end(); it++) {
    current_cell = *it;
    if (triangulation_.is_infinite(current_cell)) {
      LOG(INFO) << "Infinite cell found in external cells iterator.";
      continue;
    }

    // if (j < 1000) {
    //   gv_ << CGAL::PURPLE;
    //   gv_ << Convert_Cell_To_Tetrahedron(current_cell);
    // }

    map_.insert(std::pair<Cell_handle, handle>(current_cell, disjoint_set_.push_back(current_cell)));
    j++;
  }

  // iterate on the neighbors of each exterior cell
  LOG(INFO) << "Merging sets";
  Cell_handle current_neighbor;
  for (std::vector<Cell_handle>::iterator it = external_cells_.begin(); it != external_cells_.end(); it++) {
    // check for infinite cells
    current_cell = *it;
    if (triangulation_.is_infinite(current_cell)) {
      LOG(INFO) << "Infinite cell found in external cells iterator.";
      continue;
    }

    // join all cells that are both in the exterior set
    // do not join cells that are near the edge
    for (int neighbor = 0; neighbor <= 3; neighbor++) { // iterate over all 4 neighbors
      current_neighbor = current_cell->neighbor(neighbor); // get current neighbor

      // LOG(INFO) << "Is Vertex " << triangulation_.is_vertex(current_neighbor);
      // LOG(INFO) << "Is Edge " << triangulation_.is_edge(current_neighbor);
      // LOG(INFO) << "Is Facet " << triangulation_.is_facet(current_neighbor);
      // LOG(INFO) << "Is Cell " << triangulation_.is_cell(current_neighbor);
      
      // check if neighbor is an exterior cell
      if (!triangulation_.is_infinite(current_neighbor) && filtration_.classify(current_neighbor) == ASE::EXTERIOR &&
          filtration_.classify(current_cell, current_neighbor) == ASE::EXTERIOR) {

        // if it is, join the set with the infinite set
        disjoint_set_.unify_sets(map_[current_cell], map_[current_neighbor]);
      }  
    } 
  }

  // stop timer
  timer_.stop();
  //  LOG(INFO) << "Disjoint set construction time is: " << timer_.time();
  timer_.reset();
}


std::vector<ASE::Simplex> Connectivity_Checker::connected_component(Bare_point bp)
{
  return connected_component(bp, potential_thresh_);
}

std::vector<ASE::Simplex> Connectivity_Checker::connected_component(Bare_point bp, Potential potential)
{  
  
  Locate_type cs_lt;
  int cs_li, cs_lj;
  Cell_handle query_cell;
  query_cell = triangulation_.locate(bp, cs_lt, cs_li, cs_lj); // find start cell

  std::vector<ASE::Simplex>  cc_triangulation;

  // set potential
  if (potential != potential_thresh_) {
    set_potential(potential);
  }

  // check valid initial cell
  if (triangulation_.is_infinite(query_cell)) {
    LOG(INFO) << "Query cell is infinite";
    return cc_triangulation;
  }

  if (filtration_.classify(query_cell) != ASE::EXTERIOR) {
    LOG(INFO) << "Query cell is not in exterior of current subcomplex";
    return cc_triangulation;
  }

  // check connectivity of exterior_cells using region growing (BFS)
  Cell_handle cur_cell;
  Cell_handle cur_neighbor;
  std::set<Cell_handle> checked_cells;
  std::queue<Cell_handle> cells_to_process;
  cells_to_process.push(query_cell);
  checked_cells.insert(query_cell);
  while (cells_to_process.size() > 0) {
    // get front cell
    cur_cell = cells_to_process.front();
    cells_to_process.pop();

    // add to triangulation if valid
    if (!triangulation_.is_infinite(cur_cell) &&
        filtration_.classify(cur_cell) == ASE::EXTERIOR) {
      ASE::Simplex s;
      filtration_.simplex(cur_cell, s);
      cc_triangulation.push_back(s);
    }

    // add all valid neighbors to check set
    for (int neighbor = 0; neighbor <= 3; neighbor++) {
      cur_neighbor = cur_cell->neighbor(neighbor); // get cur neighbor

      // add neighbor to query set if separating simplex is valid
      if (filtration_.classify(cur_cell, cur_neighbor) == ASE::EXTERIOR &&
          checked_cells.find(cur_neighbor) == checked_cells.end() ) {
        cells_to_process.push(cur_neighbor);
        checked_cells.insert(cur_neighbor);
      }
    }
  }

  // add points to triangulation
  return cc_triangulation;
}


bool Connectivity_Checker::render_cc_path(Bare_point qs, Bare_point qg, Potential potential)
{ 
  std::vector<Cell_handle> start_cells;
  std::vector<Cell_handle> goal_cells;

  if (!set_goal_start_cells(qs, qg, potential, start_cells, goal_cells)) {
    return false;
  }
  
  std::map<Cell_handle, Cell_handle> closed_set; //map converts cell_handles from location queries to ds structure handles
  Cell_handle found_goal;
  int i;
  for (i = 0; i < start_cells.size(); i++) {
    
    closed_set.clear();
    //Bit of a hack, but we never follow the original cell's backpointer...should be OK!
    closed_set.insert(std::pair<Cell_handle, Cell_handle>(start_cells[i], start_cells[i]));
    if (Cell_DFS(start_cells[i], closed_set, found_goal, goal_cells)) { 
      Cell_handle plotted_cell = found_goal;
      //Follow "backpinters" incooreperated into closed set to visualize path...
      gv_ << CGAL::GREEN;
      while (plotted_cell != start_cells[i]) {
        gv_ << Convert_Cell_To_Tetrahedron(plotted_cell);
        plotted_cell = closed_set[plotted_cell];
      }
      gv_ << CGAL::RED;
      gv_ << Convert_Cell_To_Tetrahedron(plotted_cell);
      
      return true;
    }
  }

  return false;  
}



bool Connectivity_Checker::Cell_DFS(Cell_handle cur_cell, std::map<Cell_handle, Cell_handle>& closed_set, 
Cell_handle& found_goal, std::vector<Cell_handle>& goal_cells)
{
  for (int i = 0; i < goal_cells.size(); i++) {
    if (goal_cells[i] == cur_cell) {
      found_goal = goal_cells[i];
      return true;
    }
  }

  
  for (int neighbor = 0; neighbor <= 3; neighbor++) { // iterate over all 4 neighbors
    Cell_handle current_neighbor = cur_cell->neighbor(neighbor); // get current neighbor
  
    // check if neighbor is an exterior cell
    if (!triangulation_.is_infinite(current_neighbor) && filtration_.classify(current_neighbor) == ASE::EXTERIOR
        && filtration_.classify(cur_cell, current_neighbor) == ASE::EXTERIOR 
        && closed_set.count(current_neighbor) == 0) {
      closed_set.insert(std::pair<Cell_handle, Cell_handle>(current_neighbor, cur_cell));
      if (Cell_DFS(current_neighbor, closed_set, found_goal, goal_cells)) {
        return true;
      }
    }
  }

  return false;
} 



bool Connectivity_Checker::set_goal_start_cells(Bare_point qs, Bare_point qg, Potential potential, 
 						std::vector<Cell_handle>& start_cells,
  						std::vector<Cell_handle>& goal_cells)
{
  // set potential
  if (potential != potential_thresh_) {
    set_potential(potential);
  }

  //  timer_.start();

  LOG(INFO) << "Starting" << std::endl;

  // locate the start and goal cells
  Locate_type cs_lt, cg_lt;
  int cs_li, cs_lj;
  int cg_li, cg_lj;
  Cell_handle cs, cg; // start and goal cells
  cs = triangulation_.locate(qs, cs_lt, cs_li, cs_lj); // find start cell
  cg = triangulation_.locate(qg, cg_lt, cg_li, cg_lj); // find end cell

  // if either cell is infinite, query the disjoint_set_ with the infinite cell set
  if (triangulation_.is_infinite(cs)) {
    LOG(INFO) << "Start is infinite";
    return false;
  }
  if (triangulation_.is_infinite(cg)) {
    LOG(INFO) << "End is infinite";
    return false;
  }

  LOG(INFO) << "Init" << std::endl;

  std::vector<Cell_handle> incident_cells;

  // search through all possible pairs of neighboring cells to the start point if located on a vertex
  //  LOG(INFO) << "Start point " << qs;
  LOG(INFO) << "Start loc" << std::endl;
  if (cs_lt == Triangulation_3::VERTEX) {
    Vertex_handle vh = cs->vertex(cs_li);
    triangulation_.incident_cells(vh, std::back_inserter(incident_cells));
    LOG(INFO) << "Start vertex ";

    for (unsigned int i = 0; i < incident_cells.size(); i++) {
      if (triangulation_.is_cell(incident_cells[i]) && filtration_.classify(incident_cells[i]) == ASE::EXTERIOR) {
        start_cells.push_back(incident_cells[i]);
      }    
    }
    if (start_cells.size() == 0) {
      return false; // cannot be connected when the start cells are all interior
    }
  }
  // otherwise set the goal cell as the single cell that it falls into
  else if (cs_lt == Triangulation_3::CELL) {
    // LOG(INFO) << "Start cell ";
    ASE::Filtration_Cell_Class start_fc = filtration_.classify(cs);
    if (start_fc == ASE::INTERIOR) {
      //      LOG(INFO) << "Start cell interior";
      return false; // cannot be connected when the start cell is interior
    }
    else if (start_fc == ASE::EXTERIOR) {
      start_cells.push_back(cs);
    }
    else {
      LOG(INFO) << "Error: Start point is on invalid cell"; 
      return false;      
    }
  }
  else {
    LOG(INFO) << "Error: Start point is on an edge or facet"; 
    return false;
  }

  incident_cells.clear();
  LOG(INFO) << "Goal point " << qg;
  // search through all possible pairs of neighboring cells to the goal point if located on a vertex
  if (cg_lt == Triangulation_3::VERTEX) {
    Vertex_handle vh = cg->vertex(cg_li);
    triangulation_.incident_cells(vh, std::back_inserter(incident_cells));
    // LOG(INFO) << "Goal vertex ";

    for (unsigned int i = 0; i < incident_cells.size(); i++) {
      if (triangulation_.is_cell(incident_cells[i]) && filtration_.classify(incident_cells[i]) == ASE::EXTERIOR) {
        goal_cells.push_back(incident_cells[i]);
      }    
    }
    if (goal_cells.size() == 0) {
      return false; // cannot be connected when the goal cells are all interior
    }
  }
  // otherwise set the goal cell as the single cell that it falls into
  else if (cg_lt == Triangulation_3::CELL) {
    // LOG(INFO) << "Goal cell";
    ASE::Filtration_Cell_Class goal_fc = filtration_.classify(cg);

    if (goal_fc == ASE::INTERIOR) {
      return false; // cannot be connected when the goal cell is interior      
    }
    else if (goal_fc == ASE::EXTERIOR) {
      //      LOG(INFO) << "Goal cell exterior";
      goal_cells.push_back(cg);
    }
    else {
      LOG(INFO) << "Error: Goal point is on invalid cell"; 
      return false;      
    }
  }
  else {
    LOG(INFO) << "Error: Goal point is on an edge or facet"; 
    return false;
  }

  return true;
}


bool Connectivity_Checker::connected(Bare_point qs, Bare_point qg, Potential potential)
{
  if (!cache_connected_components_) {
    reset_classification();
  }

  std::vector<Cell_handle> start_cells;
  std::vector<Cell_handle> goal_cells;

  if (!set_goal_start_cells(qs, qg, potential, start_cells, goal_cells)) {
    return false;
  }
  
  // check all connectivity of all possible goals against all possible starts
  bool same_set = false;
  LOG(INFO) << "Same sets" << std::endl;
  for (unsigned int i = 0; i < !same_set && goal_cells.size(); i++) {
    for (unsigned int j = 0; j < !same_set && start_cells.size(); j++) {
      LOG(INFO) << "Checking " << i << " " << j << std::endl;
      gv_ << CGAL::DEEPBLUE;
      //gv_ << Convert_Cell_To_Tetrahedron(goal_cells[i]);
      gv_ << CGAL::WHITE;
      //gv_ << Convert_Cell_To_Tetrahedron(start_cells[j]);
      same_set = disjoint_set_.same_set(map_[goal_cells[i]], map_[start_cells[j]]);      
    }
  }

  //  LOG(INFO) << "Connected " << same_set;

  // check the time
  //  timer_.stop();
  //  LOG(INFO) << "Disconnection query time is: " << timer_.time();
  //  timer_.reset();
  return same_set; //return if they are the same set and thus connected
}

bool Connectivity_Checker::connected(Bare_point qs, Bare_point qg)
{
  return connected(qs, qg, potential_thresh_);
}

bool Connectivity_Checker::disconnected(Bare_point qs, Bare_point qg, Potential potential)
{
  return !connected(qs, qg, potential);
}

bool Connectivity_Checker::disconnected(Bare_point qs, Bare_point qg)
{
  return !connected(qs, qg, potential_thresh_);
}
bool Connectivity_Checker::connected_component_volume(Bare_point q, float& volume)
{
  // locate start cell
  Locate_type c_lt;
  int c_li, c_lj;
  Cell_handle c;
  c = triangulation_.locate(q, c_lt, c_li, c_lj);

  // if cell is infinite, query the disjoint_set_ with the infinite cell set
  if (triangulation_.is_infinite(c)) {
    LOG(INFO) << "Cell is infinite";
    return false;
  }

  // add potential neighboring initial cells
  std::vector<Cell_handle> valid_cells;
  std::vector<Cell_handle> incident_cells;
  if (c_lt == Triangulation_3::VERTEX) {
    Vertex_handle vh = c->vertex(c_li);
    triangulation_.incident_cells(vh, std::back_inserter(incident_cells));

    for (unsigned int i = 0; i < incident_cells.size(); i++) {
      if (triangulation_.is_cell(incident_cells[i]) && filtration_.classify(incident_cells[i]) == ASE::EXTERIOR) {
        valid_cells.push_back(incident_cells[i]);
      }    
    }
  }
  else if (c_lt == Triangulation_3::CELL) {
    if (filtration_.classify(c) == ASE::EXTERIOR) {
      valid_cells.push_back(c);
    }
  }


  // check all cells connected to start
  volume = 0.0f;
  filtration_.subcomplex_exterior(external_cells_);

  // initialize disjoint set structure with all of the finite external cells
  Cell_handle current_cell;
  for (std::vector<Cell_handle>::iterator it = external_cells_.begin(); it != external_cells_.end(); it++) {
    current_cell = *it;
    if (triangulation_.is_infinite(current_cell)) {
      LOG(INFO) << "Infinite cell found in external cells iterator.";
      continue;
    }

    // check if same set
    bool found_same_set = false;
    for (unsigned int i = 0; i < valid_cells.size() && !found_same_set; i++) {
      if (disjoint_set_.same_set(map_[valid_cells[i]], map_[current_cell])) {
        found_same_set = true;

        // get vertices and compute volume
        Tetrahedron_3 tetra(triangulation_.point(current_cell, 0), triangulation_.point(current_cell, 1),
                            triangulation_.point(current_cell, 2), triangulation_.point(current_cell, 3));
        volume = volume + tetra.volume();
      }      
    }
  }
  return true;
}
