#include "Alpha_shapes_disconnection.h"

// alpha shapes are disconnected
Alpha_shapes_disconnection::Alpha_shapes_disconnection(Triangulation_3 raw_triangulation, float alpha)
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
Alpha_shapes_disconnection::Reset_Classification()
{
  external_cells_.clear();
  disjoint_set_.clear();
  map_.clear();
}

// shitty interface to change the alpha value and recompute
// can be far more elegant if the alpha values step down (alpha = infty to alpha = 0)
void
Alpha_shapes_disconnection::Set_Alpha(float alpha)
{
  Reset_Classification();
  alpha_ = alpha;
  alpha_shape_.set_alpha(alpha);
  Classify_Facets();
}

Alpha_iterator
Alpha_shapes_disconnection::Alpha_Value_Begin()
{
  return alpha_shape_.alpha_begin();
}

Alpha_iterator
Alpha_shapes_disconnection::Alpha_Value_End()
{
  return alpha_shape_.alpha_end();
}

void
Alpha_shapes_disconnection::Classify_Facets()
{
  // get finite external cells
  alpha_shape_.get_alpha_shape_cells(std::back_inserter(external_cells_), 
                                     Alpha_shape_3::EXTERIOR);

  timer_.start();

  // initialize disjoint set structure with all of the finite external cells
  //  std::cout << "Init ds struct" << std::endl;
  Cell_handle current_cell;
  for (std::list<Cell_handle>::iterator it = external_cells_.begin(); it != external_cells_.end(); it++) {
    current_cell = *it;
    if (alpha_shape_.is_infinite(current_cell)) {
      std::cout << "Infinite cell found in external cells iterator." << std::endl;
      continue;
    }
    map_.insert(std::pair<Cell_handle, handle>(current_cell, disjoint_set_.push_back(current_cell)));
  }

  //  std::cout << "It on neighbors" << std::endl;
  // iterate on the neighbors of each exterior cell
  Cell_handle current_neighbor;
  for (std::list<Cell_handle>::iterator it = external_cells_.begin(); it != external_cells_.end(); it++) {

    // check for infinite cells
    current_cell = *it;
    if (alpha_shape_.is_infinite(current_cell)) {
      std::cout << "Infinite cell found in external cells iterator." << std::endl;
      continue;
    }

    // join all cells that are both in the exterior set
    // do not join cells that are near the edge
    // debugging process?
    for (int neighbor = 0; neighbor <= 3; neighbor++) { // iterate over all 4 neighbors
      current_neighbor = current_cell->neighbor(neighbor); // get current neighbor

      if (!alpha_shape_.is_infinite(current_neighbor) && alpha_shape_.classify(current_neighbor) == Alpha_shape_3::EXTERIOR) {  // check if neighbor is an exterior cell

        if (alpha_shape_.classify(current_cell, neighbor) == Alpha_shape_3::EXTERIOR) { //check if the joining face is exterior

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
Alpha_shapes_disconnection::Find_Path_Between(Bare_point qs, Bare_point qg)
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

        if (!alpha_shape_.is_infinite(current_neighbor) && alpha_shape_.classify(current_neighbor) == Alpha_shape_3::EXTERIOR) { // check if neighbor is an exterior cell

          if(alpha_shape_.classify(current_cell, neighbor) == Alpha_shape_3::EXTERIOR) { // check if the joining face is exterior
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
Alpha_shapes_disconnection::Query_disconnection(Bare_point qs, Bare_point qg)
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
      if (alpha_shape_.is_cell(incident_cells[i]) && alpha_shape_.classify(incident_cells[i]) == Alpha_shape_3::EXTERIOR) {
        start_cells.push_back(incident_cells[i]);
      }    
    }
  }
  else if (cs_lt == Triangulation_3::CELL) {
    if (alpha_shape_.classify(cs) == Alpha_shape_3::EXTERIOR) {
      start_cells.push_back(cs);
    }
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
      if (alpha_shape_.is_cell(incident_cells[i]) && alpha_shape_.classify(incident_cells[i]) == Alpha_shape_3::EXTERIOR) {
        goal_cells.push_back(incident_cells[i]);
      }    
    }
  }
  else if (cg_lt == Triangulation_3::CELL) {
    if (alpha_shape_.classify(cg) == Alpha_shape_3::EXTERIOR) {
      goal_cells.push_back(cg);
    }
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

bool
Alpha_shapes_disconnection::Component_Volume(Bare_point q, float& volume)
{
  Locate_type c_lt;
  int c_li, c_lj;
  Cell_handle c; //start and goal cells
  c = alpha_shape_.locate(q, c_lt, c_li, c_lj); //find start cell

  // if either cell is infinite, query the disjoint_set_ with the infinite cell set
  if (alpha_shape_.is_infinite(c)) {
    std::cout << "Cell is infinite" << std::endl;
    return false;
  }

  // add potential neighboring initial cells
  std::vector<Cell_handle> valid_cells;
  std::vector<Cell_handle> incident_cells;
  if (c_lt == Triangulation_3::VERTEX) {
    Vertex_handle vh = c->vertex(c_li);
    alpha_shape_.incident_cells(vh, std::back_inserter(incident_cells));

    for (unsigned int i = 0; i < incident_cells.size(); i++) {
      if (alpha_shape_.is_cell(incident_cells[i]) && alpha_shape_.classify(incident_cells[i]) == Alpha_shape_3::EXTERIOR) {
        valid_cells.push_back(incident_cells[i]);
      }    
    }
  }
  else if (c_lt == Triangulation_3::CELL) {
    if (alpha_shape_.classify(c) == Alpha_shape_3::EXTERIOR) {
      valid_cells.push_back(c);
    }
  }


  // check all cells connected to start
  volume = 0.0f;
  alpha_shape_.get_alpha_shape_cells(std::back_inserter(external_cells_), 
                                     Alpha_shape_3::EXTERIOR);

  // initialize disjoint set structure with all of the finite external cells
  Cell_handle current_cell;
  for (std::list<Cell_handle>::iterator it = external_cells_.begin(); it != external_cells_.end(); it++) {
    current_cell = *it;
    if (alpha_shape_.is_infinite(current_cell)) {
      std::cout << "Infinite cell found in external cells iterator." << std::endl;
      continue;
    }

    // check if same set
    bool found_same_set = false;
    for (unsigned int i = 0; i < valid_cells.size() && !found_same_set; i++) {
      if (disjoint_set_.same_set(map_[valid_cells[i]], map_[current_cell])) {
        found_same_set = true;

        // get vertices and compute volume
        Tetrahedron_3 tetra(alpha_shape_.point(current_cell, 0), alpha_shape_.point(current_cell, 1),
                            alpha_shape_.point(current_cell, 2), alpha_shape_.point(current_cell, 3));
        volume = volume + tetra.volume();
      }      
    }
  }
  return true;
}

bool Alpha_shapes_disconnection::Escape_Boundary(Bare_point q, float conv_alpha, float max_area,
                                                 std::vector<DifferentialTriangle>& separator_tris,
                                                 float& conv_hull_area, float& conv_hull_volume, CGAL::Geomview_stream& gv)
{
  // 1. add cells exterior to 0 shape but not exterior to the conv alpha shape to a union find struct
  //  std::cout << "Adding cells to map" << std::endl;
  std::map<Cell_handle, handle> interior_map; //map converts cell_handles from location queries to ds structure handles
  Disjoint_set interior_disjoint_set;
  Set_Alpha(0.0f);
  alpha_shape_.get_alpha_shape_cells(std::back_inserter(external_cells_), 
                                     Alpha_shape_3::EXTERIOR);

  Cell_handle current_cell;
  for (std::list<Cell_handle>::iterator it = external_cells_.begin(); it != external_cells_.end(); it++) {
    current_cell = *it;
    if (alpha_shape_.is_infinite(current_cell)) {
      std::cout << "Infinite cell found in external cells iterator." << std::endl;
      continue;
    }
    if (alpha_shape_.classify(current_cell) == Alpha_shape_3::EXTERIOR && 
        alpha_shape_.classify(current_cell, conv_alpha) != Alpha_shape_3::EXTERIOR) {
      interior_map.insert(std::pair<Cell_handle, handle>(current_cell, interior_disjoint_set.push_back(current_cell)));
    }
  }

  // 2. union cells to find connected components
  //  std::cout << "Finding connected components" << std::endl;
  Cell_handle current_neighbor;
  for (std::list<Cell_handle>::iterator it = external_cells_.begin(); it != external_cells_.end(); it++) {

    // check for infinite cells
    current_cell = *it;
    if (alpha_shape_.is_infinite(current_cell)) {
      std::cout << "Infinite cell found in external cells iterator." << std::endl;
      continue;
    }

    // check if cell is in the union find structure
    if (alpha_shape_.classify(current_cell, conv_alpha) != Alpha_shape_3::EXTERIOR) {

      // join all cells that are both in the exterior set
      // do not join cells that are near the edge
      for (int neighbor = 0; neighbor <= 3; neighbor++) { // iterate over all 4 neighbors
        current_neighbor = current_cell->neighbor(neighbor); // get current neighbor

        // check if neighboring cell satisfies criterion
        if (!alpha_shape_.is_infinite(current_neighbor) && alpha_shape_.classify(current_neighbor) == Alpha_shape_3::EXTERIOR && 
            alpha_shape_.classify(current_neighbor, conv_alpha) != Alpha_shape_3::EXTERIOR) {

          if (alpha_shape_.classify(current_cell, neighbor) == Alpha_shape_3::EXTERIOR) { //check if the joining face is exterior

            // if it is, join the set with the infinite set
            interior_disjoint_set.unify_sets(interior_map[current_cell], interior_map[current_neighbor]);
          }
        }  
      } 
    }
  }

  // 3. find cells containing q
  //  std::cout << "Finding cell(s) containing q" << std::endl;
  Locate_type c_lt;
  int c_li, c_lj;
  Cell_handle c; //start and goal cells
  c = alpha_shape_.locate(q, c_lt, c_li, c_lj); //find start cell

  // if either cell is infinite, query the disjoint_set_ with the infinite cell set
  if (alpha_shape_.is_infinite(c)) {
    std::cout << "Cell is infinite" << std::endl;
    return false;
  }

  // add potential neighboring initial cells
  std::vector<Cell_handle> valid_cells;
  std::vector<Cell_handle> incident_cells;
  if (c_lt == Triangulation_3::VERTEX) {
    Vertex_handle vh = c->vertex(c_li);
    alpha_shape_.incident_cells(vh, std::back_inserter(incident_cells));

    for (unsigned int i = 0; i < incident_cells.size(); i++) {
      if (alpha_shape_.is_cell(incident_cells[i]) && alpha_shape_.classify(incident_cells[i]) == Alpha_shape_3::EXTERIOR && 
          alpha_shape_.classify(incident_cells[i], conv_alpha) != Alpha_shape_3::EXTERIOR) {
        valid_cells.push_back(incident_cells[i]);
      }    
    }
  }
  else if (c_lt == Triangulation_3::CELL) {
    if (alpha_shape_.classify(c) == Alpha_shape_3::EXTERIOR && alpha_shape_.classify(c, conv_alpha) != Alpha_shape_3::EXTERIOR) {
      valid_cells.push_back(c);
    }
  }

  // std::cout << "Num valid cells " << valid_cells.size() << std::endl;
  // std::cout << "CONV ALPHA: " << conv_alpha << std::endl;
  //  alpha_shape_.set_alpha(conv_alpha);
  // gv.clear();
  // gv << CGAL::DEEPBLUE;
  // gv << alpha_shape_;

  // //  gv << alpha_shape_.tetrahedron(c);
  // for (unsigned int i = 0; i < valid_cells.size(); i++) {
  //   gv << CGAL::GREEN;
  //   gv << alpha_shape_.tetrahedron(valid_cells[i]);
  // }

  // 4. get boundary tris by finding cells in same connected component with neighbors in exterior
  //  std::cout << "Getting boundary tris" << std::endl;
  std::vector<Tri_3> boundary_tris;
  conv_hull_area = 0.0f;
  conv_hull_volume = 0.0f;
  for (std::list<Cell_handle>::iterator it = external_cells_.begin(); it != external_cells_.end(); it++) {
    current_cell = *it;
    if (alpha_shape_.is_infinite(current_cell)) {
      std::cout << "Infinite cell found in external cells iterator." << std::endl;
      continue;
    }

    // check if cell is in the valid set
    if (alpha_shape_.classify(current_cell) == Alpha_shape_3::EXTERIOR &&
        alpha_shape_.classify(current_cell, conv_alpha) != Alpha_shape_3::EXTERIOR) {

      for (unsigned int i = 0; i < valid_cells.size(); i++) {

        // check if cell in same connected component as containing cell
        if (interior_disjoint_set.same_set(interior_map[valid_cells[i]], interior_map[current_cell])) {

          // check neighbors to get boundary cells
          for (int neighbor = 0; neighbor <= 3; neighbor++) {
            current_neighbor = current_cell->neighbor(neighbor);

            // check if neighboring cell satisfies criterion
            if (!alpha_shape_.is_infinite(current_neighbor) &&
                alpha_shape_.classify(current_neighbor, conv_alpha) == Alpha_shape_3::EXTERIOR) {

              // add separating triangle
              Tri_3 tri = alpha_shape_.triangle(current_cell, neighbor);
              boundary_tris.push_back(tri);
            }  
          } 
        }      
      }
    }

    // add to area of convex hull boundary
    if (alpha_shape_.classify(current_cell, conv_alpha) != Alpha_shape_3::EXTERIOR) {
      Tetrahedron_3 tetra = alpha_shape_.tetrahedron(current_cell);
      conv_hull_volume += tetra.volume();

      for (int neighbor = 0; neighbor <= 3; neighbor++) {
        current_neighbor = current_cell->neighbor(neighbor);

        // check if neighboring cell satisfies criterion
        if (!alpha_shape_.is_infinite(current_neighbor) &&
            alpha_shape_.classify(current_neighbor, conv_alpha) == Alpha_shape_3::EXTERIOR) {

          // add separating triangle
          Tri_3 tri = alpha_shape_.triangle(current_cell, neighbor);
          conv_hull_area += sqrt(tri.squared_area());
        }  
      } 
    }
  }
  
  // 5. decompose each triangle into smaller components
  //  std::cout << "Boundary tris size: " << boundary_tris.size() << std::endl;
  for (unsigned int i = 0; i < boundary_tris.size(); i++) {
    Add_Tri_Centers(boundary_tris[i], max_area, separator_tris);
  }
  
  //  std::cout << "Num separators: " << separator_tris.size() << std::endl;
  // for (unsigned int j = 0; j < separator_tris.size(); j++) {
  //   gv << CGAL::RED;
  //   gv << separator_tris[j].tri;
  // }

  //  std::cout << "Done " << std::endl;
  return true;
}

bool Alpha_shapes_disconnection::Add_Tri_Centers(Tri_3 tri, float max_area, std::vector<DifferentialTriangle>& tri_vec)
{
  float sq_area = tri.squared_area();
  if (sq_area < max_area * max_area) {

    Bare_point tri_center(0.333f * (tri[0].x() + tri[1].x() + tri[2].x()),
                          0.333f * (tri[0].y() + tri[1].y() + tri[2].y()),
                          0.333f * (tri[0].z() + tri[1].z() + tri[2].z()));
    DifferentialTriangle diff_tri;
    diff_tri.center = tri_center;
    diff_tri.dir = tri.supporting_plane().orthogonal_direction();
    diff_tri.tri = tri;
    diff_tri.area = sqrt(sq_area);
    
    tri_vec.push_back(diff_tri);
    return true;
  }

  // points on the 
  Pt_3 midpoints[3];
  for (unsigned int i = 0; i < 3; i++) {
    midpoints[i] = Pt_3(0.5f * (tri[(i + 1) % 3].x() + tri[i].x()),
                        0.5f * (tri[(i + 1) % 3].y() + tri[i].y()),
                        0.5f * (tri[(i + 1) % 3].z() + tri[i].z()));
  }
  
  // decompose triangle into 4 pieces and recurse
  Tri_3 t0(tri[0], midpoints[0], midpoints[2]);
  Tri_3 t1(midpoints[0], tri[1], midpoints[1]);
  Tri_3 t2(midpoints[1], tri[2], midpoints[2]);
  Tri_3 t3(midpoints[0], midpoints[1], midpoints[2]);

  Add_Tri_Centers(t0, max_area, tri_vec);
  Add_Tri_Centers(t1, max_area, tri_vec);
  Add_Tri_Centers(t2, max_area, tri_vec);
  Add_Tri_Centers(t3, max_area, tri_vec);

  return false;
}

std::list<Weighted_point> generate_sphere_of_spheres(double x, double y, double z, double large_radius,
                                                     double small_radius, int stacks, int slices)
{
  std::list<Weighted_point> ret_list; //initialize

  // bottom point first
  Weighted_point current_append(Bare_point(x, y, (z - large_radius)), small_radius); 
  ret_list.push_back(current_append);

  //random variables needed
  double pi = 3.14159;
  double dphi = pi/stacks;
  double phi = -pi/2;
  double dtheta = 2*pi/slices;
  double theta, radial, height;

  for (int i = 0; i < stacks; i++) {
    phi += dphi;
    radial = large_radius*cos(phi);
    height = large_radius*sin(phi);
    //std::cout << "Height: " << height << "Phi: " << phi << "\n";
  
    theta = 0.0;
    for (int j = 0; j < slices; j++) {
      current_append = Weighted_point(Bare_point(x + radial*cos(theta), 
                                                 y + radial*sin(theta), z + height), small_radius);
      ret_list.push_back(current_append);
      theta += dtheta;
    }
  }

  // finish with top point 
  current_append = Weighted_point(Bare_point(x, y, (z + large_radius)), small_radius);
  ret_list.push_back(current_append);
  return ret_list;
}

// function to make random spheres
std::list<Weighted_point> generate_random_spheres(double x_min, double y_min, double z_min, 
                                                  double x_max, double y_max, double z_max, 
                                                  double rad_min, double rad_max, int number)
{
  std::list<Weighted_point> ret_list; //initialize

  for(int i = 0; i < number; i++) {
    srand((unsigned)time(NULL));
    ret_list.push_back(Weighted_point(Bare_point(((double)rand()/(double)RAND_MAX)*(x_max-x_min)+x_min,
                                                 ((double)rand()/(double)RAND_MAX)*(y_max-y_min)+y_min,
                                                 ((double)rand()/(double)RAND_MAX)*(z_max-z_min)+z_min), 
                                      ((double)rand()/(double)RAND_MAX)*(rad_max-rad_min)+rad_min));
  }
  return ret_list;
}

