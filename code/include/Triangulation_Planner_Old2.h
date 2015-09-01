/*
 * =====================================================================================
 *
 *       Filename:  Triangulation_Planner.h
 *
 *    Description:  This file contains the code for a planner that plans using Alpha Shapes
 *                  and maintains a triangulation of the sampled point set.
 *
 *        Version:  1.0
 *        Created:  07/26/2011 09:32:52 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Zoe McCarthy (zm), ZoeMcCarthy12@gmail.com
 *        Company:  University of Illinois at Urbana-Champaign
 *
 * =====================================================================================
 */


#ifndef TRIANGULATION_PLANNER_H
#define TRIANGULATION_PLANNER_H

#ifndef DEBUG
#define DEBUG
#endif

#ifndef DISPLAY_TETRA
//#define DISPLAY_TETRA
#endif

#ifndef DISPLAY_SPHERES
//#define DISPLAY_SPHERES
#endif

#ifndef DISPLAY_AS
#define DISPLAY_AS
#endif

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





//helper functions


double Distance_Squared(Configuration p1, Configuration p2)
{
    double ret_dist = 0.0;
    for(int i = 0; i < 3; i++)
    {
        ret_dist += (p1[i]-p2[i])*(p1[i]-p2[i]);
    }
    return ret_dist;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Convert_Cell_To_Tetrahedron
 *  Description:  Takes in an input cell_handle and returns the tetrahedron for gv
 * =====================================================================================
 */
    Tetrahedron
Convert_Cell_To_Tetrahedron ( Cell_handle in_cell )
{
    return Tetrahedron(in_cell->vertex(0)->point(), in_cell->vertex(1)->point(),
            in_cell->vertex(2)->point(), in_cell->vertex(3)->point());
}		/* -----  end of function Convert_Cell_To_Tetrahedron  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Compute_squared_radius
 *  Description:  This function takes in a cell and computes the orthogonal circumsphere
 *                of it.  I haven't dealt at all with degeneracy. TODO: FIX THAT!
 * =====================================================================================
 */
    double
Compute_squared_radius (Triangulation_3& dt, Cell_handle& cell_in)
{return dt.geom_traits().compute_squared_radius_smallest_orthogonal_sphere_3_object()(
	  dt.point(cell_in,0), dt.point(cell_in,1),
	  dt.point(cell_in,2), dt.point(cell_in,3));
}		/* -----  end of function Compute_squared_radius  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Compute_squared_radius
 *  Description:  This function takes in a cell and computes the orthogonal circumsphere
 *                of it.  I haven't dealt at all with degeneracy. TODO: FIX THAT!
 * =====================================================================================
 */
    double
Compute_squared_radius_facet (Triangulation_3& dt, Cell_handle& cell_in, int neighbor)
{   int i[3];
    int j = 0;
    for(int m = 0; m < 4; m++)
    {
        if(neighbor != m)
        {
            i[j++] = m;
        }
    }
    return dt.geom_traits().compute_squared_radius_smallest_orthogonal_sphere_3_object()(
	  dt.point(cell_in,i[0]), dt.point(cell_in,i[1]),
	  dt.point(cell_in,i[2]));
}		/* -----  end of function Compute_squared_radius  ----- */



/*
 * =====================================================================================
 *        Class:  Triangulation_Planner
 *  Description:  This class is an implementation of the Planner class.  It is a PRM planner
 *                that utilizes penetration depth and separation depth calculations while 
 *                sampling in order to guide sampling and inform of path existence.
 * =====================================================================================
 */
class Triangulation_Planner : public Planner
{
    public:
        /* ====================  LIFECYCLE     ======================================= */
        Triangulation_Planner (double in_xscale, double in_yscale, double in_zscale, 
        double in_xmin, double in_xmax, double in_ymin, double in_ymax, double in_zmin, double in_zmax,
        CollisionChecker& in_coll_checker, ConfigurationMapper* in_conf_mapper, CGAL::Geomview_stream& in_gv,
        double in_extra_factor);          /* constructor */
        ~Triangulation_Planner(){underlying_triangulation.clear();}

        /* ====================  ACCESSORS     ======================================= */

        /* ====================  MUTATORS      ======================================= */
        void Clear_Stored_Data();

        /* ====================  OPERATORS     ======================================= */
        Path_with_exist Find_Path(Configuration& c1, Configuration& c2); //find a path between configurations



        /* ====================  DATA MEMBERS  ======================================= */
        Point_map  point_in_collision; 

        CGAL::Geomview_stream&      gv;
    protected:
        /* ====================  DATA MEMBERS  ======================================= */
        Triangulation_3             underlying_triangulation;
        Alpha_shapes_disconnection* underlying_prover;
        Timer                       timer;



        //currently only allows for linear scaling
        //scaling is: PLAN_CONF = _scale*PRM_CONF
        double x_scale, y_scale, z_scale;

        // min/max are in PLAN_CONF units
        double x_min, x_max, y_min, y_max, z_min, z_max;

        double extra_factor;

        Cell_map cell_classification; 

        Alpha_facet_map             facet_map;
        Alpha_cell_map              cell_map;

        CollisionChecker&            coll_checker;

        ConfigurationMapper*        conf_mapper;
        
        /* ====================  MEMBER FUNCTIONS  =================================== */
        Cell_vector Find_Path_Between_And_Classify(Bare_point qs, Bare_point qg);

        Alpha_status  classify_cell(Cell_handle& s)
        // Classifies the cell `f' of the underlying Delaunay
        // tetrahedralization with respect to `A'.
          // s->radius == alpha => f interior
        {
          if(!cell_map.is_defined(s))
          {
              double alpha = Compute_squared_radius(underlying_triangulation, s);
              cell_map[s] =  alpha;
              s->set_alpha(alpha);
          }

          if (underlying_triangulation.is_infinite(s)) return EXTERIOR;
          return (s->get_alpha() <=  0.0) ? INTERIOR : EXTERIOR;
        }



        /*  Alpha_status  classify_facet(Cell_handle& s, int neighbor)
        // Classifies the cell `f' of the underlying Delaunay
        // tetrahedralization with respect to `A'.
          // s->radius == alpha => f interior
        {
          
          if(!facet_map.is_defined(Facet(s,neighbor)))
          {
              double alpha = Compute_squared_radius_facet(underlying_triangulation, s, neighbor);
              facet_map[Facet(s,neighbor)] = alpha;
              
          }

          if (underlying_triangulation.is_infinite(s)) return EXTERIOR;
          return (facet_map[Facet(s,neighbor)] <=  0.0) ? INTERIOR : EXTERIOR;
        }*/

        void Sample_Random_Configuration();


        void Sample_Configuration(Configuration& c);

        CGAL::Triple<bool, DT_Vector3, DT_Vector3> Check_Collision(Configuration& configuration);

        void Sample_Tetrahedron_Configurations(Cell_path& in_Path);

        Configuration Convert_Config_To_PRM(Configuration& in_conf)
        {
            Configuration::iterator iter = in_conf.begin();
            Configuration ret_conf;
            ret_conf.push_back((*iter)/x_scale);
            iter++;
            ret_conf.push_back((*iter)/y_scale);
            iter++;
            ret_conf.push_back((*iter)/z_scale);
            return ret_conf;
        }


        Configuration Convert_Config_From_PRM(Configuration& in_conf)
        {
            Configuration::iterator iter = in_conf.begin();
            Configuration ret_conf;
            ret_conf.push_back((*iter)*x_scale);
            iter++;
            ret_conf.push_back((*iter)*y_scale);
            iter++;
            ret_conf.push_back((*iter)*z_scale);
            return ret_conf;
        }


        Configuration Random_Point_In_Tetrahedron(Cell_handle in_cell);

        Configuration Centroid(Cell_handle in_cell)
        {
            Configuration conf;
            //prep return
            double array[] = {0.0, 0.0, 0.0};
            //calculate centroid
            for(int i = 0; i < 4; i++)
            {
                for(int j = 0; j < 3; j++)
                {
                    array[j] += .25*in_cell->vertex(i)->point()[j];
                }
            }
            for(int i = 0; i < 3; i++)
            {
                conf.push_back(array[i]);
            }
            return conf;
        }






        friend class Path_Tree;

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

    private:
        /* ====================  DATA MEMBERS  ======================================= */




}; /* -----  end of class Triangulation_Planner  ----- */




/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Triangulation_Planner::Triangulation_Planner
 *  Description:  This function implements the constructor for a Triangulation Planner
 *                class.  
 * =====================================================================================
 */
    
Triangulation_Planner::Triangulation_Planner (double in_xscale, double in_yscale, double in_zscale, 
        double in_xmin, double in_xmax, double in_ymin, double in_ymax, double in_zmin, double in_zmax,
        CollisionChecker& in_coll_checker, ConfigurationMapper* in_conf_mapper, CGAL::Geomview_stream& in_gv,
        double in_extra_factor)
        : x_scale(in_xscale), y_scale(in_yscale), z_scale(in_zscale),
        coll_checker(in_coll_checker), conf_mapper(in_conf_mapper), gv(in_gv),
        extra_factor(in_extra_factor)

        
{
    x_min = in_xmin*x_scale;
    x_max = in_xmax*x_scale;
    y_min = in_ymin*y_scale;
    y_max = in_ymax*y_scale;
    z_min = in_zmin*z_scale;
    z_max = in_zmax*z_scale;
    underlying_prover = NULL;
}		/* -----  end of function Triangulation_Planner::Triangulation_Planner  ----- */





/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Triangulation_Planner::Clear_Stored_Data
 *  Description:  This class clears all of the data associated with the Planner.  After
 *                this function is called, the planner will be reset to initialized state.
 * =====================================================================================
 */
    void
Triangulation_Planner::Clear_Stored_Data ()
{
    underlying_triangulation.clear();
    cell_classification.clear();
    point_in_collision.clear();
    delete underlying_prover;
}		/* -----  end of function Triangulation_Planner::Clear_Stored_Data  ----- */




/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Triangulation_Planner::Sample_Random_Configuration
 *  Description:  This function chooses a random configuration and checks whether or not
 *                it places the robot in collision with the obstacles.  It then calculates
 *                either penetration depth or separation depth and stores this data in a map.
 * =====================================================================================
 */
    void
Triangulation_Planner::Sample_Random_Configuration ()
{
    double x, y, z;
    Configuration conf;
    x = ((double)rand()/(double)RAND_MAX)*(extra_factor*x_max-extra_factor*x_min)+extra_factor*x_min;
    y = ((double)rand()/(double)RAND_MAX)*(extra_factor*y_max-extra_factor*y_min)+extra_factor*y_min;
    z = ((double)rand()/(double)RAND_MAX)*(extra_factor*z_max-extra_factor*z_min)+extra_factor*z_min;
    conf.push_back(x);
    conf.push_back(y);
    conf.push_back(z);
    Sample_Configuration(conf);
}		/* -----  end of function Triangulation_Planner::Sample_Random_Configuration  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Triangulation_Planner::Sample_Tetrahedron_Configurations
 *  Description:  This function takes in a Cell_path, assumed to be composed of all mixed
 *                cells, and then subdivides in each one by sampling a random point within
 *                each cell.
 * =====================================================================================
 */
    void
Triangulation_Planner::Sample_Tetrahedron_Configurations (Cell_path& in_Path)
{
    std::list<Configuration> sampled_points;
    //sample all the points first so that the cells are still in the triangulation
    for(Cell_path::iterator iter = in_Path.begin(); iter != in_Path.end(); iter++)
    {
        sampled_points.push_back(Random_Point_In_Tetrahedron(*iter));
    }
    //then insert all of them into the triangulation
    for(std::list<Configuration>::iterator iter = sampled_points.begin(); iter != sampled_points.end(); iter++)
    {
        Sample_Configuration(*iter);
    }
}		/* -----  end of function Triangulation_Planner::Sample_Tetrahedron_Configurations  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Triangulation_Planner::Random_Point_In_Tetrahedron
 *  Description:  Generates a uniform point in a tetrahedron using barycentric coordinates.
 *                Method from Generating Random Points in a Tetrahedron by C. Rocchini,
 *                P. Cignoni
 * =====================================================================================
 */
    Configuration
Triangulation_Planner::Random_Point_In_Tetrahedron ( Cell_handle in_cell )
{
    double s, t, u;
    Configuration conf;
    //generate 3 random values between 0 and 1
    s = ((double)rand()/(double)RAND_MAX);
    t = ((double)rand()/(double)RAND_MAX);
    u = ((double)rand()/(double)RAND_MAX);


    double new_s, new_t, new_u;
    //fold the sampled point into a tetrahedron from a cube
    //resulting doubles are barycentric coordinates for a tetrahedron
    if(s+t+u > 1)
    {
        if(t+u > 1)
        {
            new_s = s;
            new_t = 1-u;
            new_u = 1-s-t;
        }
        else
        {
            new_s = 1-t-u;
            new_t = t;
            new_u = s+t+u-1;
        }
    }
    else
    {
        new_s = s;
        new_t = t;
        new_u = u;
    }


    Bare_point base(in_cell->vertex(0)->point());
    Bare_point a(in_cell->vertex(1)->point());
    Bare_point b(in_cell->vertex(2)->point());
    Bare_point c(in_cell->vertex(3)->point());

    Bare_point end = base + (new_s*(a-base) + new_t*(b-base) + new_u*(c-base));
    conf.push_back(end[0]);
    conf.push_back(end[1]);
    conf.push_back(end[2]);
    return conf;
}		/* -----  end of function Triangulation_Planner::Random_Point_In_Tetrahedron  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Triangulation_Planner::Sample_Configuration
 *  Description:  This function takes in a sampled configuration and inserts it into the
 *                triangulation structure after determining appropriate distance.
 *  * =====================================================================================
 */
    void
Triangulation_Planner::Sample_Configuration (Configuration& c)
{
    
    Bare_point point = Bare_point(c[0], c[1], c[2]);
    //Cell_handle cell = underlying_triangulation.locate(point);
    //Alpha_Classification classification = cell_classification[cell];
    //if(classification == EMPTY or classification == FULL)
    //{
    //    return;
    //}
    double distance_squared_to_valid_region = 0.0;

    if(c[0] < x_min)
    {
        distance_squared_to_valid_region += (x_min-c[0])*(x_min-c[0]);
    } 
    else if(c[0] > x_max)
    {
        distance_squared_to_valid_region += (c[0]-x_max)*(c[0]-x_max);
    } 
    if(c[1] < y_min)
    {
        distance_squared_to_valid_region += (y_min-c[1])*(y_min-c[1]);
    } 
    else if(c[1] > y_max)
    {
        distance_squared_to_valid_region += (c[1]-y_max)*(c[1]-y_max);
    } 
    if(c[2] < z_min)
    {
        distance_squared_to_valid_region += (z_min-c[2])*(z_min-c[2]);
    } 
    else if(c[2] > z_max)
    {
        distance_squared_to_valid_region += (c[2]-z_max)*(c[2]-z_max);
    } 


    CGAL::Triple<bool, DT_Vector3, DT_Vector3> collision_data;
    collision_data = Check_Collision(c);
    bool collision = collision_data.first;
    double distance_squared = 0.0;
    //std::cout << "Collision?: " << (collision ? 1 : 0 ) << std::endl;
    if(collision) //reject non-collision points
    {
        distance_squared = Distance_Squared(collision_data.second, collision_data.third);
    }

    //take max of joint limit collision and obstacle collision
    if(distance_squared_to_valid_region > distance_squared)
    {
        distance_squared = distance_squared_to_valid_region;
    }

    //collision exists in some form
    if(distance_squared > 0.0)
    {
        Weighted_point new_point(Bare_point(c[0],c[1],c[2]), distance_squared);
        underlying_triangulation.insert(new_point);

#ifdef DISPLAY_SPHERES

        gv << CGAL::PURPLE;
        gv << K::Sphere_3 (Bare_point(c[0],c[1],c[2]), distance_squared);

#endif
        
    }


}		/* -----  end of function Triangulation_Planner::Sample_Configuration  ----- */




        /* 
         * ===  FUNCTION  ======================================================================
         *         Name:  Triangulation_Planner::Find_Path
         *  Description:  This function computes a path and returns that no path exists if it
         *                detects that none exists.
         * =====================================================================================
         */
    Path_with_exist
Triangulation_Planner::Find_Path (Configuration& c1, Configuration& c2)
{

    gv.clear();
    //get initial random samples going if they aren't
    int size = underlying_triangulation.number_of_vertices();
    if (size == 0)
    {
        //sample 8 corners of rectangloid, this is so that any points in infinite
        //cells are not valid samples
        //Weighted_point new_point(Bare_point(x_min,y_min,z_min), 0);
        //underlying_triangulation.insert(new_point);

        //new_point = Weighted_point(Bare_point(x_max,y_min,z_min), 0);
        //underlying_triangulation.insert(new_point);
        //
        //new_point = Weighted_point(Bare_point(x_min,y_max,z_min), 0);
        //underlying_triangulation.insert(new_point);

        //new_point = Weighted_point(Bare_point(x_min,y_min,z_max), 0);
        //underlying_triangulation.insert(new_point);

        //new_point = Weighted_point(Bare_point(x_max,y_max,z_min), 0);
        //underlying_triangulation.insert(new_point);

        //new_point = Weighted_point(Bare_point(x_max,y_min,z_max), 0);
        //underlying_triangulation.insert(new_point);

        //new_point = Weighted_point(Bare_point(x_min,y_max,z_max), 0);
        //underlying_triangulation.insert(new_point);

        //new_point = Weighted_point(Bare_point(x_max,y_max,z_max), 0);
        //underlying_triangulation.insert(new_point);


        //sample a bunch of random points inside the rectangloid to begin the search
        timer.start();
        for(int i = 0; i < 6000; i++)
        {
            Sample_Random_Configuration();
        }
        timer.stop();
        std::cout << "6000 points sampled in: " << timer.time() << " seconds." << std::endl;
        timer.reset();
    }

    //initialize new variables
    Cell_handle start, end;

    //c*_PLAN is in PLAN_CONF units
    Configuration c1_PLAN, c2_PLAN;
    c1_PLAN = Convert_Config_From_PRM(c1);
    c2_PLAN = Convert_Config_From_PRM(c2);
    Bare_point point1 = Bare_point(c1_PLAN[0], c1_PLAN[1], c1_PLAN[2]);
    Bare_point point2 = Bare_point(c2_PLAN[0], c2_PLAN[1], c2_PLAN[2]);
    
    bool adaptive_subdivision = false;

    bool Path_exists = true;
    Path_with_exist possible_path;
    //main loop of iterate, subdivide until disconnection is proven or a path is found
    //GO TO OLD CODE FOLDER FOR AN ACTUAL WORKING PATH PLANNER.  THIS IS ONLY DISCONNECTION PROVING
    for(int i = 0; i < 100; i++)
    {
#ifdef DISPLAY_TETRA
        //display debug codes
        //gv.clear();
        //gv << CGAL::GREEN;
        //gv << Convert_Cell_To_Tetrahedron(start);
        //gv << CGAL::RED;
        //gv << Convert_Cell_To_Tetrahedron(end);
#endif

        start = underlying_triangulation.locate(point1); //find start cell //TODO: possible bug here to due with the locate query returning something 
        end = underlying_triangulation.locate(point2); //find end cell   //other than a cell.  also, likely has a bug when the cell returned is a 
        std::cout << start->get_alpha();
        //check if path exists using alpha shapes disconnection.
        delete underlying_prover;
        timer.start();
        underlying_prover = new Alpha_shapes_disconnection(underlying_triangulation);
        timer.stop();
        std::cout << "Alpha Shape constructed in: " << timer.time() << " seconds." << std::endl;
        timer.reset();

        Path_exists = underlying_prover->Query_disconnection(point1, point2);
        std::cout << "Path exists?: " << (Path_exists ? "yes" : "no" )<< std::endl;

        Cell_vector cell_path;
        timer.start();
        cell_path = Find_Path_Between_And_Classify(point1, point2);
        timer.stop();
        std::cout << "Alpha structure on regular triangulation constructed in: " << timer.time() << " seconds." << std::endl;
        timer.reset();
        if(cell_path.empty())
        {
            if(!Path_exists)
            {
                std::cout << "New method does not agree with old method" << std::endl;
            }
            else
            {
                std::cout << "New method does agree with old method" << std::endl;
            }
        }
        else
        {
            if(Path_exists)
            {
                std::cout << "New method does not agree with old method" << std::endl;
            }
            else
            {
                std::cout << "New method does agree with old method" << std::endl;
            }

        }


        if( i % 1 == 0)
        {
            std::cout << underlying_triangulation.number_of_cells() << " cells in the triangulation." << std::endl;
            std::cout << underlying_prover->Num_sets() << " sets in the disconnection prover." << std::endl;
            underlying_prover->Num_sets();

        }

#ifdef DISPLAY_AS

        gv.clear();
        gv << CGAL::DEEPBLUE;
        gv << underlying_prover->as;
        gv << CGAL::PURPLE;

        //if(Path_exists)
        //{
        //    Cell_vector cell_path(underlying_prover->Find_Path_Between(point1, point2));
        //    while(!cell_path.empty())
        //    {
        //        gv << Convert_Cell_To_Tetrahedron(cell_path.back());
        //        cell_path.pop_back();
        //    }

        //}
        gv << CGAL::GREEN;
        gv << Convert_Cell_To_Tetrahedron(start);
        gv << CGAL::RED;
        gv << Convert_Cell_To_Tetrahedron(end);
        
        //as.get_alpha_shape_facets(std::back_inserter(external_cells), 
	    //      Alpha_shape_3::EXTERIOR);

#endif
        
        if(!Path_exists)
        {
            possible_path.second = false;
            return possible_path;
        }
        else
        {
            if(adaptive_subdivision)//mixed, subdivide
            {
                //do nothing currently
            }
            else
            {
                //else, sample a bunch more points
                timer.start();
                for(int j = 0; j < 6000; j++)
                {
                    Sample_Random_Configuration();
                }
                timer.stop();
                std::cout << "6000 points sampled in: " << timer.time() << " seconds." << std::endl;
                timer.reset();
            }
        }
    }
    possible_path.second = true;
    std::cout << "Computation timed out.  Need more iterations or something is wrong." << std::endl;
    std::cout << underlying_triangulation.number_of_finite_cells() << " cells in the triangulation." << std::endl;
    return possible_path;
}		/* -----  end of function Triangulation_Planner::Find_Path  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Triangulation_Planner::Check_Collision
 *  Description:  Links the CollisionChecker and ConfigurationMapper
 * =====================================================================================
 */
CGAL::Triple<bool, DT_Vector3, DT_Vector3>
Triangulation_Planner::Check_Collision (Configuration& configuration)
{
    Configuration PRM_CONF = Convert_Config_To_PRM(configuration);
    return coll_checker.Check_Intersection(conf_mapper->Collision_Transformations(PRM_CONF));
}		/* -----  end of function Triangulation_Planner::Check_Collision  ----- */


Cell_vector 
Triangulation_Planner::Find_Path_Between_And_Classify(Bare_point qs, Bare_point qg)
{
  Cell_handle cs, cg; //start and goal cells

  cs = underlying_triangulation.locate(qs); //find start cell //TODO: possible bug here to due with the locate query returning something 
  cg = underlying_triangulation.locate(qg); //find end cell   //other than a cell.  also, likely has a bug when the cell returned is a 

  Cell_vector  cell_path;

  //if either cell is infinite, query the ds with the infinite cell set
  if(underlying_triangulation.is_infinite(cs))
  {
      std::cout << "Start is infinite" << std::endl;
      return cell_path;
  }


  if(underlying_triangulation.is_infinite(cg))
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

    if(underlying_triangulation.is_infinite(current_cell))
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
          if(underlying_triangulation.is_infinite(current_neighbor))  //if infinite
          {
            //std::cout << "Infinite cell skipped in disjoint set initialization" << std::endl;
            //do nothing if infinite
          }
          else if(classify_cell(current_neighbor) == EXTERIOR)  // check if neighbor is an exterior cell
          {
            if(Compute_squared_radius_facet(underlying_triangulation, current_cell, neighbor) > 0 ) //check if the joining face is exterior
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


#endif //TRIANGULATION_PLANNER_H
