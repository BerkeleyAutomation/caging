// Utility functions for the caging project

#ifndef UTIL_H
#define UTIL_H

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/number_utils.h>
#include <CGAL/utility.h>
#include <SOLID/SOLID.h>
#include <SOLID/MT_Scalar.h>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

#include <math.h>

#include "fcl/traversal/traversal_node_bvhs.h"
#include "fcl/traversal/traversal_node_setup.h"
#include "fcl/collision_node.h"
#include "fcl/collision.h"
#include "fcl/BV/BV.h"
#include "fcl/shape/geometric_shapes.h"
#include "fcl/narrowphase/narrowphase.h"

#include "CageConfigurationParser.h"

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
typedef CGAL::Aff_transformation_3<Kernel> CGAL_Aff_Transform;

double Distance_Squared(DT_Vector3 p1, DT_Vector3 p2);

double Magnitude_Squared(const DT_Vector3 p);
double Magnitude_Squared(DT_Vector3 p);

Eigen::Matrix4f CreatePose(float tx, float ty, float theta);
Eigen::Matrix3f CreatePose2D(float tx, float ty, float theta);

void Pose2DToParams2D(Eigen::Matrix3f pose, float& tx, float& ty, float& theta);
void Pose3DToParams2D(Eigen::Matrix4f pose, float& tx, float& ty, float& theta);


void UniformRandomPose(float& tx, float& ty, float& theta,
                       float min_tx = 0, float max_tx = 1,
                       float min_ty = 0, float max_ty = 1,
                       float min_theta = 0, float max_theta = M_PI);

CGAL_Aff_Transform EigenPoseToCGAL(Eigen::Matrix4f pose);

// Modified from FCL test utility
bool LoadOBJFile(const char* filename, std::vector<fcl::Vec3f>& points, std::vector<fcl::Triangle>& triangles);

#endif // UTIL_H
