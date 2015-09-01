/*
 * Mesh class for collision checking
 *
 */
#ifndef MESH_H
#define MESH_H

#include "CageConfigurationParser.h"
#include "ShapeFactory.hpp"
#include "SDF.h"
#include "Util.h"

#include <fstream>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <vector>

// Box2d includes
#include <Box2D/Box2D.h>

#include <opencv2/opencv.hpp>

// CGAL includes
#include <CGAL/Cartesian.h>
#include <CGAL/Aff_transformation_3.h>

#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/Triangulation_geomview_ostream_3.h>
#include <CGAL/IO/Polyhedron_geomview_ostream.h>
#include <CGAL/utility.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Timer.h>

// SOLID includes
#include <SOLID/SOLID.h>
#include <SOLID/MT_Point3.h>
#include <SOLID/MT_Vector3.h>

// PolyDepth includes
#include <PQP.h>
#include <C2A/LinearMath.h>
#include <C2A/C2A.h>
#include <PolyDepth/PolyDepth.h>
#include <PolyDepth/MeshObject.h>

// FCL includes
#include "fcl/traversal/traversal_node_bvhs.h"
#include "fcl/traversal/traversal_node_setup.h"
#include "fcl/collision_node.h"
#include "fcl/collision.h"
#include "fcl/BV/BV.h"
#include "fcl/shape/geometric_shapes.h"
#include "fcl/narrowphase/narrowphase.h"

// Eigen includes
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

typedef CGAL::Timer Timer;

struct Poly {
  std::vector<Eigen::Vector3f> vertices;
  Eigen::Matrix4f pose;
};

namespace ASE {
  typedef Eigen::Vector2f Vertex2d;
  typedef Eigen::Vector2f Direction2d;
  typedef Eigen::Vector2i Edge2d;  
  typedef Eigen::Vector3f Vertex3d;
  typedef Eigen::Vector3f Direction3d;
  typedef Eigen::Vector3i Edge3d;  
};


struct MeshCollisionResult
{
  bool collision;
  float distance;
  float time;
  int coll_i; // index of component of 1st mesh in collision 
  int coll_j; // index of component of 2nd mesh in collision
};

// holds refs to two vertices in a list
struct MeshEdge
{
  unsigned int i;
  unsigned int j;
};

class Mesh
{
 public:
  Mesh();
  Mesh(ComponentConfig config, Eigen::Matrix4f root_pose = Eigen::Matrix4f::Identity());
  Mesh(std::vector<ComponentConfig> configs, Eigen::Matrix4f root_pose = Eigen::Matrix4f::Identity());
  ~Mesh();

 public:
  void RenderToGeomview(CGAL::Geomview_stream& gv);
  void RenderToImage(cv::Mat& image, char color = 'r');
  void RenderTrisToImage(cv::Mat& image, float scale = 1.0f, char color = 'r', bool draw_com = false);

  // returns a lower bound on the penetration depth if a collision occurs
  static MeshCollisionResult LowerBoundCollision(Mesh* mesh1, Mesh* mesh2, bool print=false);

  // returns a upper bound (tighter than LB) on the penetration depth if a collision occurs
  static MeshCollisionResult UpperBoundCollision(Mesh* mesh1, Mesh* mesh2, bool compute_separation = false);

  // returns the minimum distance between the meshes
  static MeshCollisionResult Distance(Mesh* mesh1, Mesh* mesh2);

 public:
  std::vector<std::vector<fcl::Vec3f> > BoundaryVertices();
  void ConcaveVertices(std::vector<ASE::Vertex2d>& verts, std::vector<std::pair<ASE::Direction2d, ASE::Direction2d> >& directions);
  void BoundaryEdges(std::vector<ASE::Vertex2d>& boundary_edges, std::vector<ASE::Direction2d>& directions);
  
  float MaxMomentArm();

  void SetPose(Eigen::Matrix4f pose, bool update_tris = false, bool print = false);
  bool SetRelativePoses(std::vector<ComponentConfig> configs);
  //  void SetComponents(std::vector<ComponentConfig> configs);

  SDF* GetSDF();
  std::vector<SDF*> SDFs();
  std::vector<SDFCollisionCache*> SDFCollCaches();
  std::vector<MeshObject*> PolyDepthObjects();
  std::vector<Transform> PolyDepthPoses();

  std::vector<std::vector<fcl::Vec3f> > FCLVertices();
  std::vector<std::vector<fcl::Triangle> > FCLTriangles();

  Eigen::Matrix4f Pose();
  std::vector<Eigen::Matrix4f> PosesWorldFrame();
  std::vector<Eigen::Matrix4f> PosesRootFrame();
  unsigned int NumComponents();

  float Mass();
  float MomentOfInertia();
  std::vector< std::vector<Poly> > ConvexPieces2D();
  bool Initialized() { return initialized_; }

 protected:
  void Initialize(ComponentConfig config, Eigen::Matrix4f root_pose = Eigen::Matrix4f::Identity());
  void Initialize(std::vector<ComponentConfig> configs, Eigen::Matrix4f root_pose = Eigen::Matrix4f::Identity());

 private:
  void InitializeSDF(std::vector<ComponentConfig> configs);
  void InitializeMeshObject(std::vector<ComponentConfig> configs);
  void InitializeFCLObject(std::vector<ComponentConfig> configs);
  void InitializeRendering(bool update_tris = true);
  void InitializeMass();
  void CenterFCLVertices();
  void ConvexDecomposition2D(float extrusion = 50.0f);

  // 2D and 3D pose conversions
  Eigen::Matrix3f RootFrameToWorldFrame(Eigen::Matrix3f pose_root_frame);
  Eigen::Matrix3f WorldFrameToRootFrame(Eigen::Matrix3f pose_world_frame);
  Eigen::Matrix4f RootFrameToWorldFrame(Eigen::Matrix4f pose_root_frame);
  Eigen::Matrix4f WorldFrameToRootFrame(Eigen::Matrix4f pose_world_frame);
  Eigen::Matrix3f Pose3DTo2D(Eigen::Matrix4f pose_3d);
  Transform Pose3DToPolyDepth(Eigen::Matrix4f pose);

 private:
  //  std::vector<fcl::Vec3f> vertices_; // might not need this
  std::vector<float> scales_;
  std::vector<SDF*> sdfs_;                // signed distance field for lower bound on PD
  std::vector<SDFCollisionCache*> sdf_coll_caches_;
  std::vector<MeshObject*> mesh_objects_; // polydepth mesh for upper bound on PD (more accurate but can't prove caging
  Eigen::Matrix4f pose_root_world_; // the pose of the root in the world frame
  std::vector<Eigen::Matrix4f> poses_comp_root_; // poses of components relative to the "root" world frame

  std::vector<Transform> pd_poses_; // polydepth version of poses in world frame
  std::vector<std::vector<Triangle_3> > render_triangles_; // triangle rep for rendering (in world frame)
  bool tris_need_update_;
  bool initialized_;

  float mass_;
  float moment_of_inertia_;

  std::vector<fcl::Vec3f> centroids_;
  std::vector<std::vector<fcl::Vec3f> > fcl_vertices_;
  std::vector<std::vector<fcl::Triangle> > fcl_triangles_;

  std::vector<std::vector<Poly> > convex_pieces_2d_;
};

#endif // MESH_H
