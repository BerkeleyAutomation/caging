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

#include "Planner.h"
#include "CollisionChecker.h"
#include "ConfigurationMapper.h"


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
#include <math.h>

typedef std::vector<double>                                 Configuration;
typedef std::list< Configuration >                          Path;
typedef std::pair< Path, bool >                             Path_with_exist;




typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Regular_triangulation_euclidean_traits_3<K>   Gt;

typedef CGAL::Timer                                         Timer;

typedef CGAL::Regular_triangulation_3<Gt>                   Triangulation_3;

typedef Triangulation_3::Cell_handle                        Cell_handle;
typedef Triangulation_3::Vertex_handle                      Vertex_handle;
typedef Triangulation_3::Facet                              Facet;
typedef Triangulation_3::Edge                               Edge;
typedef Gt::Weighted_point                                  Weighted_point;
typedef Gt::Bare_point                                      Bare_point;

typedef K::Tetrahedron_3                                    Tetrahedron;


enum Alpha_Classification {EMPTY, FULL, MIXED};


typedef std::list< Cell_handle >                            Cell_path;
typedef std::pair< Cell_path, Alpha_Classification >        Cell_path_with_exist;

typedef CGAL::Unique_hash_map<Vertex_handle, bool>         Point_map;
typedef CGAL::Unique_hash_map<Cell_handle, Alpha_Classification>  Cell_map;





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
        CollisionChecker& in_coll_checker, ConfigurationMapper* in_conf_mapper, CGAL::Geomview_stream& in_gv);          /* constructor */
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



        //currently only allows for linear scaling
        //scaling is: PLAN_CONF = _scale*PRM_CONF
        double x_scale, y_scale, z_scale;

        // min/max are in PLAN_CONF units
        double x_min, x_max, y_min, y_max, z_min, z_max;

        Cell_map cell_classification; 

        CollisionChecker&            coll_checker;

        ConfigurationMapper*        conf_mapper;
        
        /* ====================  MEMBER FUNCTIONS  =================================== */
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
                //implements A* with euclidean distance between centroids of cells
                Cell_path_with_exist Search_For_Path()
                {
                    Path_Tree_Node* current_node;
                    Path_Tree_Node* end_node = NULL;
                    Cell_path_with_exist return_path;

                    //classify start
                    root = Classify(start, NULL);
                    if(planner.cell_classification[root->Give_Cell()] == FULL)
                    {
                        std::cout << "Error.  Start Cell is full." << std::endl;
                    }

                    //check for a wholly empty path first
                    while((!unchecked_empty.empty()) && (end_node == NULL))
                    {
                        current_node = unchecked_empty.back();
                        unchecked_empty.pop_back();
                        end_node = Add_Children_Nodes(current_node);
                    }
                    while((!unchecked_mixed.empty()) && (end_node == NULL))
                    {
                        while((!unchecked_empty.empty()) && (end_node == NULL))
                        {
                            current_node = unchecked_empty.back();
                            unchecked_empty.pop_back();
                            end_node = Add_Children_Nodes(current_node);
                        }
                        if(end_node == NULL)
                        {
                            current_node = unchecked_mixed.back();
                            unchecked_mixed.pop_back();
                            end_node = Add_Children_Nodes(current_node);
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

                    private:
                        Cell_handle cell;
                        Path_Tree_Node* parent;
                        std::list<Path_Tree_Node*> children;
                };

                Path_Tree_Node* root;
                typedef CGAL::Unique_hash_map<Cell_handle, Path_Tree_Node*>         Node_map;
                Node_map                                                            closed_set, open_set;

                typedef CGAL::Unique_hash_map<Cell_handle, double>                  Cost_map;
                Cost_map                                                            g_score, h_score, f_score;

                typedef std::vector<Path_Tree_Node*> Node_Vector;
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
                            //check if its been classified
                            if(!planner.cell_classification.is_defined(new_cell)) //not been classified
                            {
                                new_node = Classify(new_cell, cur_root);
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
                                if(planner.cell_classification[new_cell] == MIXED)
                                {
                                    new_node = new Path_Tree_Node(new_cell, cur_root);
                                    unchecked_mixed.push_back(new_node);
                                    closed_set[new_cell] = new_node;
                                    if(new_cell == end)
                                    {
                                        return new_node;
                                    }
                                }
                                else if(planner.cell_classification[new_cell] == EMPTY)
                                {
                                    new_node = new Path_Tree_Node(new_cell, cur_root);
                                    unchecked_empty.push_back(new_node);
                                    closed_set[new_cell] = new_node;
#ifdef DEBUG
                                    planner.gv << CGAL::PURPLE;
                                    planner.gv << Convert_Cell_To_Tetrahedron(new_cell);
#endif
                                    if(new_cell == end)
                                    {
                                        return new_node;
                                    }
                                }
                                else //full
                                {
                                    closed_set[new_cell] = NULL;
#ifdef DEBUG
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
                            unchecked_mixed.push_back(new_node);
                            closed_set[new_cell] = new_node;
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
                                closed_set[new_cell] = NULL;
                                //do not add into tree structure
                                //std::cout << "Full cell found." << std::endl;
#ifdef DEBUG
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
                                unchecked_empty.push_back(new_node);
                                closed_set[new_cell] = new_node;
#ifdef DEBUG
                                planner.gv << CGAL::PURPLE;
                                planner.gv << Convert_Cell_To_Tetrahedron(new_cell);
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
        CollisionChecker& in_coll_checker, ConfigurationMapper* in_conf_mapper, CGAL::Geomview_stream& in_gv)
        : x_scale(in_xscale), y_scale(in_yscale), z_scale(in_zscale),
        coll_checker(in_coll_checker), conf_mapper(in_conf_mapper), gv(in_gv)

        
{
    x_min = in_xmin*x_scale;
    x_max = in_xmax*x_scale;
    y_min = in_ymin*y_scale;
    y_max = in_ymax*y_scale;
    z_min = in_zmin*z_scale;
    z_max = in_zmax*z_scale;
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
    x = ((double)rand()/(double)RAND_MAX)*(x_max-x_min)+x_min;
    y = ((double)rand()/(double)RAND_MAX)*(y_max-y_min)+y_min;
    z = ((double)rand()/(double)RAND_MAX)*(z_max-z_min)+z_min;
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
    if(c[0] < x_min || c[0] > x_max)
    {
        //std::cout << "Sample violates x bounds" << std::endl;
        return;
    }
    if(c[1] < y_min || c[1] > y_max)
    {
        //std::cout << "Sample violates y bounds" << std::endl;
        return;
    }
    if(c[2] < z_min || c[2] > z_max)
    {
        //std::cout << "Sample violates z bounds" << std::endl;
        return;
    }

    CGAL::Triple<bool, DT_Vector3, DT_Vector3> collision_data;
    collision_data = Check_Collision(c);
    bool collision = collision_data.first;

    double distance = Distance_Squared(collision_data.second, collision_data.third);

    Weighted_point new_point(Bare_point(c[0],c[1],c[2]), distance);
    
    //map point to whether or not it is in collision or not.

    Vertex_handle new_vertex_handle;
    new_vertex_handle = underlying_triangulation.insert(new_point);

    point_in_collision[new_vertex_handle] = collision;


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

    //get initial random samples going if they aren't
    int size = underlying_triangulation.number_of_vertices();
    if (size == 0)
    {
        //sample 8 corners of rectangloid, this is so that any points in infinite
        //cells are not valid samples
        Configuration conf;

        conf.push_back(x_min);
        conf.push_back(y_min);
        conf.push_back(z_min);
        Sample_Configuration(conf);
        conf.clear();

        conf.push_back(x_max);
        conf.push_back(y_min);
        conf.push_back(z_min);
        Sample_Configuration(conf);
        conf.clear();
        

        conf.push_back(x_min);
        conf.push_back(y_max);
        conf.push_back(z_min);
        Sample_Configuration(conf);
        conf.clear();


        conf.push_back(x_min);
        conf.push_back(y_min);
        conf.push_back(z_max);
        Sample_Configuration(conf);
        conf.clear();

        conf.push_back(x_max);
        conf.push_back(y_max);
        conf.push_back(z_min);
        Sample_Configuration(conf);
        conf.clear();


        conf.push_back(x_max);
        conf.push_back(y_min);
        conf.push_back(z_max);
        Sample_Configuration(conf);
        conf.clear();


        conf.push_back(x_min);
        conf.push_back(y_max);
        conf.push_back(z_max);
        Sample_Configuration(conf);
        conf.clear();


        conf.push_back(x_max);
        conf.push_back(y_max);
        conf.push_back(z_max);
        Sample_Configuration(conf);
        conf.clear();


        //sample a bunch of random points inside the rectangloid to begin the search
        for(int i = 0; i < 100; i++)
        {
            Sample_Random_Configuration();
        }
    }

    //initialize new variables
    Cell_handle start, end;
    Path_Tree* search_tree = NULL;
    Cell_path_with_exist possible_cell_path;
    Path_with_exist possible_path;

    //c*_PLAN is in PLAN_CONF units
    Configuration c1_PLAN, c2_PLAN;
    c1_PLAN = Convert_Config_From_PRM(c1);
    c2_PLAN = Convert_Config_From_PRM(c2);
    Bare_point point1 = Bare_point(c1_PLAN[0], c1_PLAN[1], c1_PLAN[2]);
    Bare_point point2 = Bare_point(c2_PLAN[0], c2_PLAN[1], c2_PLAN[2]);

    //main loop of iterate, subdivide until disconnection is proven or a path is found
    for(int i = 0; i < 350; i++)
    {
        start = underlying_triangulation.locate(point1);
        end = underlying_triangulation.locate(point2);
#ifdef DEBUG
        //display debug codes
        gv.clear();
        gv << CGAL::GREEN;
        gv << Convert_Cell_To_Tetrahedron(start);
        gv << CGAL::RED;
        gv << Convert_Cell_To_Tetrahedron(end);
#endif
        search_tree = new Path_Tree(start, end, this);
        possible_cell_path = search_tree->Search_For_Path();
        if(possible_cell_path.second == FULL) //proven that no path exists
        {
            possible_path.first.clear();
            possible_path.second = false;
            delete search_tree;
            std::cout << underlying_triangulation.number_of_finite_cells() << " cells in the triangulation." << std::endl;
            return possible_path;
        }
        else if(possible_cell_path.second == EMPTY) // path found
        {  //connect a path from the start to finish via the centroids of empty cells
            possible_path.first.clear();
            possible_path.second = true;
            possible_path.first.push_back(c1);
            for(Cell_path::iterator iter = possible_cell_path.first.begin(); iter != possible_cell_path.first.end(); iter++)
            {
                Configuration temp_config = Centroid(*iter);
                possible_path.first.push_back(Convert_Config_To_PRM(temp_config));
            }
            possible_path.first.push_back(c2);
            delete search_tree;
            std::cout << underlying_triangulation.number_of_finite_cells() << " cells in the triangulation." << std::endl;
            return possible_path;
        }
        else //mixed, subdivide
        {
            Sample_Tetrahedron_Configurations(possible_cell_path.first);
            delete search_tree;
        }
    }
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




#endif //TRIANGULATION_PLANNER_H
