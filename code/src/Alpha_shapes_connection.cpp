#include "Alpha_shapes_connection.h"

// alpha shapes are disconnected
Alpha_shapes_connection::Alpha_shapes_connection(Triangulation_3 raw_triangulation, float alpha)
  : alpha_(alpha),
    alpha_shape_(raw_triangulation, alpha, Alpha_shape_3::GENERAL)  //construct the alpha shape, O(n^2)
{
  // Alpha_iterator alpha_it;
  // int i = 0;
  // for (alpha_it = alpha_shape_.alpha_begin(); alpha_it != alpha_shape_.alpha_end(); alpha_it++, i++) {
  //   std::cout << "Interesting alpha " << i << " " << (*alpha_it) << std::endl;
  // }

  Classify_Facets();
}  

void
Alpha_shapes_connection::Reset_Classification()
{
  interior_cells_.clear();
  disjoint_set_.clear();
  map_.clear();
}

// shitty interface to change the alpha value and recompute
// can be far more elegant if the alpha values step down (alpha = infty to alpha = 0)
void
Alpha_shapes_connection::Set_Alpha(float alpha)
{
  Reset_Classification();
  alpha_ = alpha;
  alpha_shape_.set_alpha(alpha);
  Classify_Facets();
}

Alpha_iterator
Alpha_shapes_connection::Alpha_Value_Begin()
{
  return alpha_shape_.alpha_begin();
}

Alpha_iterator
Alpha_shapes_connection::Alpha_Value_End()
{
  return alpha_shape_.alpha_end();
}

void
Alpha_shapes_connection::Classify_Facets()
{
  // get finite external cells
  alpha_shape_.get_alpha_shape_cells(std::back_inserter(interior_cells_), 
                                     Alpha_shape_3::INTERIOR);
  timer_.start();

  // initialize disjoint set structure with all of the finite external cells
  Cell_handle current_cell;
  for (std::list<Cell_handle>::iterator it = interior_cells_.begin(); it != interior_cells_.end(); it++) {
    current_cell = *it;
    if (alpha_shape_.is_infinite(current_cell)) {
      std::cout << "Infinite cell found in interior cells iterator." << std::endl;
      continue;
    }
    map_.insert(std::pair<Cell_handle, handle>(current_cell, disjoint_set_.push_back(current_cell)));
  }

  // iterate on the neighbors of each exterior cell
  Cell_handle current_neighbor;
  for (std::list<Cell_handle>::iterator it = interior_cells_.begin(); it != interior_cells_.end(); it++) {

    // check for infinite cells
    current_cell = *it;
    if (alpha_shape_.is_infinite(current_cell)) {
      std::cout << "Infinite cell found in interior cells iterator." << std::endl;
      continue;
    }

    // join all cells that are both in the exterior set
    // do not join cells that are near the edge
    // debugging process?
    for (int neighbor = 0; neighbor <= 3; neighbor++) { // iterate over all 4 neighbors
      current_neighbor = current_cell->neighbor(neighbor); // get current neighbor

      if (!alpha_shape_.is_infinite(current_neighbor) && alpha_shape_.classify(current_neighbor) == Alpha_shape_3::INTERIOR) {  // check if neighbor is an exterior cell

        if (alpha_shape_.classify(current_cell, neighbor) == Alpha_shape_3::INTERIOR) { //check if the joining face is exterior

          // if it is, join the set with the infinite set
          disjoint_set_.unify_sets(map_[current_cell], map_[current_neighbor]);
        }
      }  
    } 
  }

  // stop timer
  timer_.stop();
  //  std::cout << "Disjoint set construction time is: " << timer_.time() << std::endl;
  timer_.reset();
}

Cell_vector
Alpha_shapes_connection::Find_Path_Between(Bare_point qs, Bare_point qg)
{
  Cell_handle cs, cg; // start and goal cells
  cs = alpha_shape_.locate(qs); // find start cell //TODO: possible bug here to due with the locate query returning something 
  cg = alpha_shape_.locate(qg); // find end cell   //other than a cell.  also, likely has a bug when the cell returned is a 

  Cell_vector cell_path;

  // if either cell is infinite, query the disjoint_set_ with the infinite cell set
  if(alpha_shape_.is_infinite(cs)) {
    std::cout << "Start is infinite" << std::endl;
    return cell_path;
  }
  if(alpha_shape_.is_infinite(cg)) {
    std::cout << "End is infinite" << std::endl;
    return cell_path;
  }

  Cell_vector open_cells;
  open_cells.push_back(cs);

  Cell_handle current_cell, current_neighbor;
  Cell_map came_from;
  came_from[cs] = NULL;

  // iterate through open cells
  while (!open_cells.empty()) {

    // check open cells
    current_cell = open_cells.back();
    open_cells.pop_back();

    // skip infinte cells
    if(alpha_shape_.is_infinite(current_cell)) {
      std::cout << "Infinite cell found in external cells iterator." << std::endl;
      continue;
    }

    //do not join cells that are near the edge
    //debugging process?
    for(int neighbor = 0; neighbor <=3; neighbor++) { // iterate over all 4 neighbors

      current_neighbor = current_cell->neighbor(neighbor); //get current neighbor

      if(!came_from.is_defined(current_neighbor)) {

        if (!alpha_shape_.is_infinite(current_neighbor) && alpha_shape_.classify(current_neighbor) == Alpha_shape_3::INTERIOR) {

          if(alpha_shape_.classify(current_cell, neighbor) == Alpha_shape_3::INTERIOR) {
              //if it is, join the set with the infinite set
            came_from[current_neighbor] = current_cell;
            open_cells.push_back(current_neighbor);
            if(current_neighbor == cg)  { // end found
              break;
            }
          }  
        }  
      }
    } 

    // found the end
    if(current_neighbor == cg) { // end found
      break;
    }
  }

  // check if the end was found
  if(current_neighbor == cg) {//end found
     current_cell = current_neighbor;

     while(current_cell != NULL) { // backtrace until the start cell
       cell_path.push_back(current_cell);
       current_cell = came_from[current_cell];
     }
  }
  
  return cell_path;
}

bool
Alpha_shapes_connection::Query_connection(Bare_point qs, Bare_point qg)
{
  //  timer_.start();
  Locate_type cs_lt, cg_lt;
  int cs_li, cs_lj;
  int cg_li, cg_lj;
  Cell_handle cs, cg; //start and goal cells
  cs = alpha_shape_.locate(qs, cs_lt, cs_li, cs_lj); //find start cell
  cg = alpha_shape_.locate(qg, cg_lt, cg_li, cg_lj); //find end cell

  // if either cell is infinite, query the disjoint_set_ with the infinite cell set
  if (alpha_shape_.is_infinite(cs)) {
    std::cout << "Start is infinite" << std::endl;
    return false;
  }
  if (alpha_shape_.is_infinite(cg)) {
    std::cout << "End is infinite" << std::endl;
    return false;
  }

  std::vector<Cell_handle> start_cells;
  std::vector<Cell_handle> goal_cells;
  std::vector<Cell_handle> incident_cells;

  // search through all possible pairs of neighboring cells to the start point if located on a vertex
  if (cs_lt == Triangulation_3::VERTEX) {
    Vertex_handle vh = cs->vertex(cs_li);
    alpha_shape_.incident_cells(vh, std::back_inserter(incident_cells));

    for (unsigned int i = 0; i < incident_cells.size(); i++) {
      if (alpha_shape_.is_cell(incident_cells[i]) && alpha_shape_.classify(incident_cells[i]) == Alpha_shape_3::INTERIOR) {
        start_cells.push_back(incident_cells[i]);
      }    
    }
  }
  else if (cs_lt == Triangulation_3::CELL) {
    start_cells.push_back(cs);
  }
  else {
    std::cout << "Error: Start point is on an edge or facet" << std::endl; 
    return false;
  }

  // search through all possible pairs of neighboring cells to the goal point if located on a vertex
  incident_cells.clear();
  if (cg_lt == Triangulation_3::VERTEX) {
    Vertex_handle vh = cg->vertex(cg_li);
    alpha_shape_.incident_cells(vh, std::back_inserter(incident_cells));

    for (unsigned int i = 0; i < incident_cells.size(); i++) {
      if (alpha_shape_.is_cell(incident_cells[i]) && alpha_shape_.classify(incident_cells[i]) == Alpha_shape_3::INTERIOR) {
        goal_cells.push_back(incident_cells[i]);
      }    
    }
  }
  else if (cg_lt == Triangulation_3::CELL) {
    goal_cells.push_back(cg);
  }
  else {
    std::cout << "Error: Goal point is on an edge or facet" << std::endl; 
    return false;
  }

  // std::cout << "Possible incident start " << start_cells.size() << std::endl;
  // std::cout << "Possible incident goal " << goal_cells.size() << std::endl;

  // check all possible goals against all possible starts
  bool same_set = false;
  for (unsigned int i = 0; i < !same_set && goal_cells.size(); i++) {
    for (unsigned int j = 0; j < !same_set && start_cells.size(); j++) {
      same_set = disjoint_set_.same_set(map_[goal_cells[i]], map_[start_cells[j]]);      
    }
  }

  // check the time
  //  timer_.stop();
  //  std::cout << "Disconnection query time is: " << timer_.time() << std::endl;
  //  timer_.reset();
  return same_set; //return if they are the same set and thus connected
}

