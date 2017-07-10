/*
 * =====================================================================================
 *
 *       Filename:  ConfigurationMapper.h
 *
 *    Description:  This class is a class that is contained in a PRM instantiation.
 *                  It maps configurations to lists of affine transformations for
 *                  each section of the robot.
 *
 *        Version:  1.0
 *        Created:  07/22/2011 01:21:21 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Zoe McCarthy (zm), ZoeMcCarthy12@gmail.com
 *        Company:  University of Illinois at Urbana-Champaign
 *
 * =====================================================================================
 */


#ifndef CONF_MAPPER_H
#define CONF_MAPPER_H

#include "Typedef.h"

/*
 * =====================================================================================
 *        Class:  ConfigurationMapper
 *  Description:  This class is a member in the PRM class.  It provides an interface
 *                for the PRM to convert configurations into a list of affine transformations
 *                for the robot's different polyhedra.
 * =====================================================================================
 */
class ConfigurationMapper
{
    public:
        /* ====================  LIFECYCLE     ======================================= */

        /* ====================  ACCESSORS     ======================================= */

        /* ====================  MUTATORS      ======================================= */

        /* ====================  OPERATORS     ======================================= */
        virtual std::list<CGAL_Aff_Transform> Display_Transformations(Configuration& configuration) = 0;

        virtual double** Collision_Transformations(Configuration& configuration) = 0;


    protected:
        /* ====================  DATA MEMBERS  ======================================= */

    private:
        /* ====================  DATA MEMBERS  ======================================= */

}; /* -----  end of class ConfigurationMapper  ----- */



#endif // CONF_MAPPER_H
