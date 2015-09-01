/*********************
Zoe McCarthy, May 2011

//defines  Alpha_shapes_disconnection class.
//supports an alpha shape for weighted points (spheres) in R^3
//has an auxillary structure for querying two points and determining
//if they are in the same connected component of the free space
********************/

#ifndef ALPHA_H
#define ALPHA_H

////// ********************    PREAMBLE    ********************/

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Union_find.h>
#include <list>
#include <CGAL/Timer.h>
#include <iostream>

#ifndef DEBUG
#define DEBUG
#endif


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Regular_triangulation_euclidean_traits_3<K> Gt;

typedef CGAL::Timer                         Timer;

typedef CGAL::Alpha_shape_vertex_base_3<Gt>         Vb;
typedef CGAL::Alpha_shape_cell_base_3<Gt>           Fb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb> Tds;
typedef CGAL::Regular_triangulation_3<Gt,Tds>       Triangulation_3;
typedef CGAL::Alpha_shape_3<Triangulation_3>        Alpha_shape_3;

typedef Alpha_shape_3::Cell_handle          Cell_handle;
typedef Alpha_shape_3::Vertex_handle        Vertex_handle;
typedef Alpha_shape_3::Facet                Facet;
typedef Alpha_shape_3::Edge                 Edge;
typedef Gt::Weighted_point                  Weighted_point;
typedef Gt::Bare_point                      Bare_point;

typedef Cell_handle                         T;
typedef CGAL::Union_find<T>                 Disjoint_set;  
typedef Disjoint_set::handle                handle;

typedef std::vector<Cell_handle>            Cell_vector;



////// ********************    CLASS DEFINITION    ********************/



//this class will allow querying of two points to test 
//whether or not they are disconnected
//its underlying structure is a disjoint set of 
//connected components of cells that are not in the alpha shape
class Alpha_shapes_disconnection
{
public:

Alpha_shape_3 as;
std::list<Cell_handle> external_cells;


//construct the disjoint set that contains the connectivity information
Alpha_shapes_disconnection(Triangulation_3 RT) : as(RT, 0, 
			Alpha_shape_3::GENERAL)  //construct the alpha shape, O(n^2)
{
    bool debug_print = true;


  //get finite external cells
  as.get_alpha_shape_cells(std::back_inserter(external_cells), 
				Alpha_shape_3::EXTERIOR);



  if(debug_print)
    timer.start();
  Cell_handle current_cell;
  //initialize disjoint set structure with all of the finite external cells
  for(std::list<Cell_handle>::iterator it = external_cells.begin(); it != external_cells.end(); it++)
  {
    current_cell = *it;
    if(as.is_infinite(current_cell))
    {
        std::cout << "Infinite cell found in external cells iterator." << std::endl;
        continue;
    }
    map.insert(std::pair<Cell_handle, handle> (current_cell, ds.push_back(current_cell) ) );
  }


  Cell_handle current_neighbor;
  //iterate on the neighbors of each exterior cell
  for(std::list<Cell_handle>::iterator it = external_cells.begin(); it != external_cells.end(); it++)
  {
    current_cell = *it;
    if(as.is_infinite(current_cell))
    {
        std::cout << "Infinite cell found in external cells iterator." << std::endl;
        continue;
    }
    //std::cout << "Alpha for exterior cell is: " << current_cell->get_alpha() << std::endl;
    //do not join cells that are near the edge
    //debugging process?
    for(int neighbor = 0; neighbor <=3; neighbor++) //iterate over all 4 neighbors
    {
      current_neighbor = current_cell->neighbor(neighbor); //get current neighbor
      if(as.is_infinite(current_neighbor))  //if infinite
      {
        //std::cout << "Infinite cell skipped in disjoint set initialization" << std::endl;
        //do nothing if infinite
      }
      else if(as.classify(current_neighbor) == Alpha_shape_3::EXTERIOR)  // check if neighbor is an exterior cell
      {
        if(as.classify(current_cell, neighbor) == Alpha_shape_3::EXTERIOR) //check if the joining face is exterior
        {
          //if it is, join the set with the infinite set
          ds.unify_sets(map[current_cell], map[current_neighbor]);
        }

      }  
    } 

  }
  if(debug_print){
    timer.stop();
    std::cout << "Disjoint set construction time is: " << timer.time() << std::endl;
    timer.reset();}


}  




Cell_vector Find_Path_Between(Bare_point qs, Bare_point qg)
{
  Cell_handle cs, cg; //start and goal cells

  cs = as.locate(qs); //find start cell //TODO: possible bug here to due with the locate query returning something 
  cg = as.locate(qg); //find end cell   //other than a cell.  also, likely has a bug when the cell returned is a 

  Cell_vector  cell_path;

  //if either cell is infinite, query the ds with the infinite cell set
  if(as.is_infinite(cs))
  {
      std::cout << "Start is infinite" << std::endl;
      return cell_path;
  }


  if(as.is_infinite(cg))
  {
      std::cout << "End is infinite" << std::endl;
      return cell_path;
  }


  Cell_vector open_cells;

  open_cells.push_back(cs);

  Cell_handle current_cell, current_neighbor;

  typedef CGAL::Unique_hash_map<Cell_handle, Cell_handle>         Cell_map;
  Cell_map came_from;

  came_from[cs] = NULL;

  while(!open_cells.empty())
  {
    current_cell = open_cells.back();
    open_cells.pop_back();

    if(as.is_infinite(current_cell))
    {
        std::cout << "Infinite cell found in external cells iterator." << std::endl;
        continue;
    }
    //do not join cells that are near the edge
    //debugging process?
    for(int neighbor = 0; neighbor <=3; neighbor++) //iterate over all 4 neighbors
    {
      current_neighbor = current_cell->neighbor(neighbor); //get current neighbor
      if(!came_from.is_defined(current_neighbor))
      {
          if(as.is_infinite(current_neighbor))  //if infinite
          {
            //std::cout << "Infinite cell skipped in disjoint set initialization" << std::endl;
            //do nothing if infinite
          }
          else if(as.classify(current_neighbor) == Alpha_shape_3::EXTERIOR)  // check if neighbor is an exterior cell
          {
            if(as.classify(current_cell, neighbor) == Alpha_shape_3::EXTERIOR) //check if the joining face is exterior
            {
              //if it is, join the set with the infinite set
              came_from[current_neighbor] = current_cell;
              open_cells.push_back(current_neighbor);
              if(current_neighbor == cg)  //end found
              {
                  break;
              }
            }

          }  

      }
    } 
    if(current_neighbor == cg) //end found
    {
        break;
    }

  }

  if(current_neighbor == cg) //end found
  {
     current_cell = current_neighbor;
     while(current_cell != NULL) //backtrace until the start cell
     {
         cell_path.push_back(current_cell);
         current_cell = came_from[current_cell];
     }

  }
  
  return cell_path;


}




bool Query_disconnection(Bare_point qs, Bare_point qg)
{
  bool debug_print = true;
  if(debug_print)
    timer.start();
  Cell_handle cs, cg; //start and goal cells

  cs = as.locate(qs); //find start cell //TODO: possible bug here to due with the locate query returning something 
  cg = as.locate(qg); //find end cell   //other than a cell.  also, likely has a bug when the cell returned is a 
                                        // 0,1, or 2 simplex with adj cells that are both interior and ext


  //if either cell is infinite, query the ds with the infinite cell set
  if(as.is_infinite(cs))
  {
      std::cout << "Start is infinite" << std::endl;
      return false;
  }


  if(as.is_infinite(cg))
  {
      std::cout << "End is infinite" << std::endl;
      return false;
  }


  bool same_set(ds.same_set(map[cs],map[cg]));
  if(debug_print){
    timer.stop();
    std::cout << "Disconnection query time is: " << timer.time() << std::endl;
    timer.reset();}
  return same_set; //return if they are the same set and thus connected



}


int Num_sets()
{

  return ds.number_of_sets();

}





private:


Timer timer;
Disjoint_set ds;
std::map<Cell_handle, handle> map; //map converts cell_handles from location queries to ds structure handles


};


#endif // ALPHA_H
