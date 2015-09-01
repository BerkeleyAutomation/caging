// A star search for a path given a triangulation
// Considered deprecated for the purposes of the caging project

#include "Planner.h"
#include "CollisionChecker.h"
#include "ConfigurationMapper.h"
#include "Alpha_shapes_disconnection.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Union_find.h>
#include <CGAL/Timer.h>
#include <CGAL/utility.h>
#include <CGAL/IO/Geomview_stream.h>

#include <algorithm>
#include <iostream>
#include <list>
#include <vector>
#include <queue>
#include <math.h>

typedef std::vector<double>                                 Configuration;
typedef std::list< Configuration >                          Path;
typedef std::pair< Path, bool >                             Path_with_exist;




typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Regular_triangulation_euclidean_traits_3<K>   Gt;
typedef CGAL::Timer                                         Timer;

typedef CGAL::Alpha_shape_vertex_base_3<Gt>         Vb;
typedef CGAL::Alpha_shape_cell_base_3<Gt>           Fb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb> Tds;
typedef CGAL::Regular_triangulation_3<Gt,Tds>       Triangulation_3;

typedef Triangulation_3::Cell_handle                        Cell_handle;
typedef Triangulation_3::Vertex_handle                      Vertex_handle;
typedef Triangulation_3::Facet                              Facet;
typedef Triangulation_3::Edge                               Edge;
typedef Gt::Weighted_point                                  Weighted_point;
typedef Gt::Bare_point                                      Bare_point;

typedef K::Tetrahedron_3                                    Tetrahedron;


enum Alpha_Classification {DEFAULT, EMPTY, FULL, MIXED};

enum Alpha_status { EXTERIOR, INTERIOR };


typedef std::list< Cell_handle >                            Cell_path;
typedef std::vector<Cell_handle>                            Cell_vector;
typedef std::pair< Cell_path, Alpha_Classification >        Cell_path_with_exist;

typedef CGAL::Unique_hash_map<Vertex_handle, bool>         Point_map;
typedef CGAL::Unique_hash_map<Cell_handle, Alpha_Classification>  Cell_map;
typedef CGAL::Unique_hash_map<Cell_handle, double>         Alpha_cell_map;
typedef CGAL::Unique_hash_map<Facet, double>               Alpha_facet_map;

class Path_Tree
{
 public:
  //constructor
 Path_Tree(Cell_handle in_start, Cell_handle in_end, Triangulation_Planner* in_planner):
  start(in_start), end(in_end), planner(*in_planner), unchecked_empty(), unchecked_mixed()
    {
    }

  //destructor
  ~Path_Tree()
    {
      Clear();
    }



  //main function to call. This initiates search for a path and returns if none is found.
  //implements A* with heuristic euclidean distance between centroids of cells
  //modified to give preference to empty cells first
  Cell_path_with_exist Search_For_Path()
  {
#ifdef DISPLAY_TETRA
    planner.gv.clear();
#endif 
#ifdef DISPLAY_SPHERES
    //planner.gv.clear();
#endif 
    Path_Tree_Node* current_node;
    Path_Tree_Node* end_node = NULL;
    Cell_path_with_exist return_path;

    //std::cout << "Default value = to DEFAULT? " << ((planner.cell_classification.default_value() == DEFAULT) ? 1 : 0) << std::endl;

    //classify start
    root = Classify(start, NULL);
    g_score[root] = 0.0;
    double dist_start_end = Distance_Between(start, end);
    h_score[root] = dist_start_end;
    f_score[root] = dist_start_end;
                    
    if(planner.cell_classification[start] == FULL)
      {
        std::cout << "Error.  Start Cell is full." << std::endl;
      }
    else if(planner.cell_classification[start] == EMPTY)
      {
        unchecked_empty.push(Node_with_cost(root, dist_start_end));
      }
    else if(planner.cell_classification[start] == MIXED)
      {
        unchecked_mixed.push(Node_with_cost(root, dist_start_end));
      }

    //check for a wholly empty path first
    while((!unchecked_empty.empty()) && (end_node == NULL))
      {
        current_node = unchecked_empty.top().first;
        unchecked_empty.pop();
        end_node = Add_Children_Nodes(current_node);
        closed_set[current_node->Give_Cell()] = current_node;
      }
    while((!unchecked_mixed.empty()) && (end_node == NULL))
      {
        while((!unchecked_empty.empty()) && (end_node == NULL))
          {
            current_node = unchecked_empty.top().first;
            unchecked_empty.pop();
            end_node = Add_Children_Nodes(current_node);
            closed_set[current_node->Give_Cell()] = current_node;
          }
        if(end_node == NULL)
          {
            current_node = unchecked_mixed.top().first;
            unchecked_mixed.pop();
            end_node = Add_Children_Nodes(current_node);
            closed_set[current_node->Give_Cell()] = current_node;
          }
      }
                    
    if(end_node != NULL)//end found. mixed or free path?
      {
        bool mixed = false;
        current_node = end_node;
        Cell_handle current_cell;

        //follow the leaf of the tree up to the root
        while(current_node != NULL)
          {
            current_cell = current_node->Give_Cell();
            if(!mixed)//all previous cells in the path have been empty.
              {
                if(planner.cell_classification[current_cell] == MIXED)
                  {
                    return_path.first.clear(); //we only want mixed in order to subdivide
                    mixed = true;
                  }
                return_path.first.push_front(current_cell);
              }
            else if(planner.cell_classification[current_cell] == MIXED)
              {
                return_path.first.push_front(current_cell);
              }
            current_node = current_node->Give_Parent();
          }
        if(mixed)
          {
            return_path.second = MIXED;
          }
        else
          {
            return_path.second = EMPTY;
          }
      }
                        
    else //no path found. none exists.
      {
        return_path.second = FULL;
      }
                    
    return return_path;
  }


                


 private:
  Cell_handle start, end;
  Triangulation_Planner& planner;


  class Path_Tree_Node
  {
  public:

  Path_Tree_Node(Cell_handle in_cell, Path_Tree_Node* in_parent):cell(in_cell), parent(in_parent)
    { 
      if(parent != NULL)
        {
          parent->Add_Child(this);
        }
    }

    Path_Tree_Node* Give_Parent()
    {
      return parent;
    }

    Cell_handle Give_Cell()
    {
      return cell;
    }

    std::list<Path_Tree_Node*>& Give_Children()
      {
        return children;
      }

    void Add_Child(Path_Tree_Node* child)
    {
      children.push_back(child);
    }

    void Remove_Child(Path_Tree_Node* child)
    {
      for(std::list<Path_Tree_Node*>::iterator iter = children.begin(); iter != children.end(); iter++)
        {
          if ((*iter) == child)
            {
              children.erase(iter);
              break;
            }
        }
    }

    void Set_New_Parent(Path_Tree_Node* new_parent)
    {
      if(new_parent != parent)
        {
          if(parent != NULL)
            {
              parent->Remove_Child(this);
            }
          if(new_parent != NULL)
            {
              new_parent->Add_Child(this);
            }
          parent = new_parent;
        }
    }

  private:
    Cell_handle cell;
    Path_Tree_Node* parent;
    std::list<Path_Tree_Node*> children;
  };

  Path_Tree_Node* root;
  typedef CGAL::Unique_hash_map<Cell_handle, Path_Tree_Node*>         Node_map;
  Node_map                                                            closed_set, set_of_nodes;

  typedef CGAL::Unique_hash_map<Path_Tree_Node*, double>                  Cost_map;
  Cost_map                                                            g_score, h_score, f_score;

  typedef std::pair<Path_Tree_Node*, double>               Node_with_cost;


  class Node_comparison
  {
  public:
    bool operator() (const Node_with_cost& lhs, const Node_with_cost& rhs) const
    {
      return (lhs.second > rhs.second);
    }
  };

  typedef std::priority_queue< Node_with_cost, std::vector<Node_with_cost>, Node_comparison > Node_Vector;
  Node_Vector unchecked_empty;
  Node_Vector unchecked_mixed;

  //this function allows the tree to be grown for a given root.  it analyzes
  //the neighbors cells and classifies them according to their alpha value
  //and whether they are in collision or not to determine if they are empty,
  //mixed, or full.  returns non-NULL only if the goal is reached
  Path_Tree_Node* Add_Children_Nodes(Path_Tree_Node* cur_root)
  {
    for(int j = 0; j < 4; j++)
      {
        Cell_handle new_cell = cur_root->Give_Cell()->neighbor(j);
        //check if its been visited
        if(!closed_set.is_defined(new_cell))
          {
            //reject infinite cells
            if(planner.underlying_triangulation.is_infinite(new_cell))
              {
                break;
              }

            Path_Tree_Node* new_node = NULL;
            double tentative_g_score = g_score[cur_root] + Distance_Between(cur_root->Give_Cell(), new_cell);
            bool tentative_is_better;
            Alpha_Classification new_classification(planner.cell_classification[new_cell]);
            Node_with_cost node_with_cost;

            //check if its been adjacent to before
            if(!set_of_nodes.is_defined(new_cell))
              {
                tentative_is_better = true;
                //check if its been classified
                if(new_classification == DEFAULT) //not been classified
                  {
                    new_node = Classify(new_cell, cur_root);
                    new_classification = planner.cell_classification[new_cell];
                    if(new_node != NULL)//FULL if still NULL
                      {
                        if(new_node->Give_Cell() == end) // end found
                          {
                            return new_node;
                          }
                      }
                  }
                else //been classified
                  {
                    if(new_classification == MIXED)
                      {
                        new_node = new Path_Tree_Node(new_cell, cur_root);
                        set_of_nodes[new_cell] = new_node;
                        if(new_cell == end)
                          {
                            return new_node;
                          }
                      }
                    else if(new_classification == EMPTY)
                      {
                        new_node = new Path_Tree_Node(new_cell, cur_root);
                        set_of_nodes[new_cell] = new_node;
#ifdef DISPLAY_TETRA
                        //planner.gv << CGAL::PURPLE;
                        //planner.gv << Convert_Cell_To_Tetrahedron(new_cell);
#endif
                        if(new_cell == end)
                          {
                            return new_node;
                          }
                      }
                    else if(new_classification == FULL)//full
                      {
                        set_of_nodes[new_cell] = NULL;
#ifdef DISPLAY_TETRA
                        planner.gv << CGAL::DEEPBLUE;
                        planner.gv << Convert_Cell_To_Tetrahedron(new_cell);
#endif
                        //do not add into tree structure
                        if(new_cell == end)
                          {
                            new_node = new Path_Tree_Node(new_cell, cur_root);
                            std::cout << "Error.  End cell is full." << std::endl;
                            return new_node;
                          }
                      }
                  }
              }
            else
              {
                new_node = set_of_nodes[new_cell];
                if (new_classification != FULL and tentative_g_score < g_score[new_node])
                  {
                    tentative_is_better = true;
                  }
                else
                  {
                    tentative_is_better = false;
                  }
              }

            if (tentative_is_better and new_classification != FULL)
              {
                new_node->Set_New_Parent(cur_root);
                double tentative_h_score = Distance_Between(new_cell, end);
                double tentative_f_score = tentative_g_score + tentative_h_score;
                g_score[new_node] = tentative_g_score;
                h_score[new_node] = tentative_h_score;
                f_score[new_node] = tentative_f_score;
                Node_with_cost new_node_with_cost(new_node, tentative_f_score);
                if( new_classification == EMPTY)
                  {
                    unchecked_empty.push(new_node_with_cost);
                  }
                else if( new_classification == MIXED)
                  {
                    unchecked_mixed.push(new_node_with_cost);
                  }


              }
          }
      }
    return NULL; //signal the end has not been found
  }


  //this function is a cost estimate that is based around the distance between centroids
  double Distance_Between(Cell_handle cell_start, Cell_handle cell_end)
  {
    return sqrt( Distance_Squared(planner.Centroid(cell_start), planner.Centroid(cell_end)) );
  }


  //this function classifies a cell, creates a new node for it, 
  //and returns a pointer to the resulting node
  Path_Tree_Node* Classify(Cell_handle new_cell, Path_Tree_Node* cur_root)
  {
    Path_Tree_Node* new_node = NULL;
    if(!planner.underlying_triangulation.is_infinite(new_cell))
      {
        double radius = Compute_squared_radius(planner.underlying_triangulation, new_cell);
        if(radius > 1000.0)
          {
            //std::cout << "Large radius, something may be wrong: " << std::endl;
            //for(int i = 0; i < 4; i++)
            //{
            //    std::cout << "Point" << i << " =: ";
            //    for(int j = 0; j < 3; j++)
            //    {
            //        std::cout << planner.underlying_triangulation.point(new_cell, i)[j] << ", ";
            //    }
            //    std::cout << endl;
            //}
          }
        //std::cout << "Radius of cell: " << radius << std::endl;
        if(radius > 0) // not in alpha-shape
          {
            planner.cell_classification[new_cell] = MIXED;
            new_node = new Path_Tree_Node(new_cell, cur_root);
            set_of_nodes[new_cell] = new_node;
            //std::cout << "Mixed cell found." << std::endl;
          }
        else //in alpha shape. full or empty?
          {
            //sampling one of the four vertices is enough since it is in the alpha shape
            Vertex_handle sample_vertex_in_cell0 = new_cell->vertex(0);
            Vertex_handle sample_vertex_in_cell1 = new_cell->vertex(1);
            Vertex_handle sample_vertex_in_cell2 = new_cell->vertex(2);
            Vertex_handle sample_vertex_in_cell3 = new_cell->vertex(3);
            bool collision0 = planner.point_in_collision[sample_vertex_in_cell0];
            bool collision1 = planner.point_in_collision[sample_vertex_in_cell1];
            bool collision2 = planner.point_in_collision[sample_vertex_in_cell2];
            bool collision3 = planner.point_in_collision[sample_vertex_in_cell3];
            //make sure that all vertices agree or I'm FUCKED
            if((collision0 == collision1) && (collision1 == collision2) && (collision2 == collision3) && (collision0 == true) )
              ;//std::cout << "all vertices agree and are in collision" << std::endl;
            else if((collision0 == collision1) && (collision1 == collision2) && (collision2 == collision3) && (collision0 == false) )
              ;//std::cout << "all vertices agree and are not in collision" << std::endl;
            else
              {
                std::cout << "all vertices do not agree" << std::endl;
              }

            if(collision0) // full
              {
                planner.cell_classification[new_cell] = FULL;
                set_of_nodes[new_cell] = NULL;
                //do not add into tree structure
                //std::cout << "Full cell found." << std::endl;
#ifdef DISPLAY_TETRA
                planner.gv << CGAL::DEEPBLUE;
                planner.gv << Convert_Cell_To_Tetrahedron(new_cell);
#endif
                if(new_cell == end)
                  {
                    new_node = new Path_Tree_Node(new_cell, cur_root);
                    std::cout << "Error.  End cell is full." << std::endl;
                  }
              }
            else // empty
              {
                planner.cell_classification[new_cell] = EMPTY;
                new_node = new Path_Tree_Node(new_cell, cur_root);
                //std::cout << "Empty cell found." << std::endl;
                set_of_nodes[new_cell] = new_node;
#ifdef DISPLAY_TETRA
                //planner.gv << CGAL::PURPLE;
                //planner.gv << Convert_Cell_To_Tetrahedron(new_cell);
#endif
              }
          }
      }
    return new_node;
  }

  //clears the entire tree when called on the root
  void Clear(Path_Tree_Node* clear_root)
  {
    if(!clear_root->Give_Children().empty())
      {
        for(std::list<Path_Tree_Node*>::iterator iter = clear_root->Give_Children().begin(); 
            iter != clear_root->Give_Children().end(); iter++)
          {
            Clear(*iter);
          }
      }
    delete clear_root;
  }
                    
  //clears the entire tree
  void Clear()
  {
    Clear(root);
  }



};
