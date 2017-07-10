// Utility functions for the caging project
#ifndef UTIL_H
#define UTIL_H

#include <fcl/traversal/traversal_node_bvhs.h>
#include <fcl/traversal/traversal_node_setup.h>
#include <fcl/collision_node.h>
#include <fcl/collision.h>
#include <fcl/BV/BV.h>
#include <fcl/shape/geometric_shapes.h>
#include <fcl/narrowphase/narrowphase.h>

#include "CageConfigurationParser.h"
#include "Typedef.h"
#include "Simplex.h"
#include <boost/filesystem.hpp>

void wait ( int seconds );
bool file_exists(const char* file_name); 
std::string get_filename(std::string full_path);
Eigen::Vector2d rotate_vector(Eigen::Vector2d input_vector, float theta_offset);
//Given vectors of points in a 3D triangle, return its area
float triangle_area(std::vector< Eigen::Vector3d > triangle_points);

double Distance_Squared(Eigen::Vector3d p1, Eigen::Vector3d p2);
double Magnitude_Squared(Eigen::Vector3d p);

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

Tetrahedron_3 Convert_Cell_To_Tetrahedron ( Cell_handle in_cell );


// Variants for simplices!
Tetrahedron_3 Convert_Simplex_To_Tetrahedron ( ASE::Simplex in_simplex );

Tri_3 Convert_Simplex_To_Triangle ( ASE::Simplex in_simplex );

struct DifferentialTriangle
{
  Bare_point center;
  Dir_3 dir;
  Tri_3 tri;
  float area;
};
#endif // UTIL_H
