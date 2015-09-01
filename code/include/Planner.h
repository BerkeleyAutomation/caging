/*
 * =====================================================================================
 *
 *       Filename:  Planner.h
 *
 *    Description:  This file contains the code for the Planner class, which is an
 *                  abstract class for PRM path planners.
 *
 *        Version:  1.0
 *        Created:  07/26/2011 09:17:10 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Zoe McCarthy (zm), ZoeMcCarthy12@gmail.com
 *        Company:  University of Illinois at Urbana-Champaign
 *
 * =====================================================================================
 */


#ifndef PLANNER_H
#define PLANNER_H

#include <list>
#include <vector>

typedef std::vector<double>                                 Configuration;
typedef std::list< Configuration >                          Path;
typedef std::pair< Path, bool >                             Path_with_exist;


/*
 * =====================================================================================
 *        Class:  Planner
 *  Description:  This class contains virtual member functions for sampling points randomly.
 *                It also maintains a structure (graph) for storing sampled points and their
 *                topology.  
 * =====================================================================================
 */
class Planner
{
    public:
        /* ====================  LIFECYCLE     ======================================= */
        Planner(){}                             /* constructor */

        /* ====================  ACCESSORS     ======================================= */

        /* ====================  MUTATORS      ======================================= */
        virtual void Clear_Stored_Data() = 0;

        /* ====================  OPERATORS     ======================================= */
        virtual Path_with_exist Find_Path(Configuration& c1, Configuration& c2) = 0; //find a path between configurations

    protected:
        /* ====================  DATA MEMBERS  ======================================= */

    private:
        /* ====================  DATA MEMBERS  ======================================= */

}; /* -----  end of class Planner  ----- */

#endif //PLANNER_H
