#include "Mesh.hpp"

#include "fcl/traversal/traversal_node_bvhs.h"
#include "fcl/traversal/traversal_node_setup.h"

#include "Simulation/Box2DShapeFactory.hpp"
#include "Simulation/ConvergenceTest.h"

#include <ccd/ccd.h>
#include <SOLID/MT_Quaternion.h>

#include <glog/logging.h>

bool g_p;

struct Fake {
  float fake;
};

Fake pt_at_me;

// Support function for convex polygon
void support(const void* obj, const ccd_vec3_t *dir, ccd_vec3_t* vec)
{
  Poly* poly = (Poly*)(obj);

  // get direction
  Eigen::Vector3f d;
  d << dir->v[0], dir->v[1], dir->v[2];
  if (d.norm() == 0.0f) {
    vec->v[0] = 0.0f;
    vec->v[1] = 0.0f;
    vec->v[2] = 0.0f;
    return;
  }
  d = d / d.norm();
 
  float max_dot = -FLT_MAX;
  int max_ind = -1;
  Eigen::Vector4f max_v;

  // get vertex with max x coordinate
  // std::cout << "Dir " << d.transpose() << std::endl;

  if (g_p) {
    std::cout << "Dir " << d.transpose() << std::endl;
    std::cout << "Pose " << poly->pose << std::endl;
  }

  for (unsigned int i = 0; i < poly->vertices.size(); i++) {
    //    std::cout << "V " << i << " " << poly->vertices[i].transpose() << std::endl;

    Eigen::Vector4f v;
    v << poly->vertices[i], 1.0f;
    Eigen::Vector4f p = poly->pose * v;

    float dot = d.dot(p.head<3>());
    // std::cout << "P " << i << " " << p.transpose() << " " << dot << std::endl;
    if (dot > max_dot) {
      max_dot = dot;//p(0);
      max_ind = i;
      max_v = p;
    }
  }

  if (g_p) {
    std::cout << "Chose index " << max_ind << " " << max_v.transpose() << std::endl;
    std::cout << std::endl;
  }  
  // assign support vertex
  ccdVec3Set(vec, max_v(0), max_v(1), max_v(2));
}

// empty constructor
Mesh::Mesh()
  : initialized_(false)
{
}

Mesh::Mesh(ComponentConfig config, Eigen::Matrix4f root_pose)
  : pose_root_world_(root_pose),
    tris_need_update_(false),
    initialized_(true)
{
  std::vector<ComponentConfig> configs;
  configs.push_back(config);
  InitializeSDF(configs);
  InitializeMeshObject(configs);
  InitializeFCLObject(configs);
  InitializeRendering(true);
  InitializeMass();
  ConvexDecomposition2D();
}

Mesh::Mesh(std::vector<ComponentConfig> configs, Eigen::Matrix4f root_pose)
  : pose_root_world_(root_pose),
    tris_need_update_(false),
    initialized_(true)
{
  InitializeSDF(configs);
  InitializeMeshObject(configs);
  InitializeFCLObject(configs);
  InitializeRendering(true);
  InitializeMass();
  ConvexDecomposition2D();
}

Mesh::~Mesh()
{
  for (unsigned int i = 0; i < sdfs_.size(); i++) {
    delete sdfs_[i];
    delete sdf_coll_caches_[i];
    delete mesh_objects_[i];
  }
}

void Mesh::Initialize(ComponentConfig config, Eigen::Matrix4f root_pose)
{
  std::vector<ComponentConfig> configs;
  configs.push_back(config);
  Initialize(configs, root_pose);
}

void Mesh::Initialize(std::vector<ComponentConfig> configs, Eigen::Matrix4f root_pose)
{
  initialized_ = true;
  pose_root_world_ = root_pose;
  tris_need_update_ = false;
  
  scales_.clear();
  sdfs_.clear();
  sdf_coll_caches_.clear();
  mesh_objects_.clear();
  poses_comp_root_.clear();
  pd_poses_.clear();
  render_triangles_.clear();
  fcl_vertices_.clear();
  fcl_triangles_.clear();

  InitializeSDF(configs);
  InitializeMeshObject(configs);
  InitializeFCLObject(configs);
  InitializeRendering(true);
  InitializeMass();
  ConvexDecomposition2D();
}

void Mesh::InitializeSDF(std::vector<ComponentConfig> configs)
{
  ComponentConfig config;
  float theta;
  float tx;
  float ty;
  Eigen::Matrix3f pose_comp_root; // pose of component wrt root
  Eigen::Matrix3f pose_comp_world; // pose of component wrt world

  for (unsigned int i = 0; i < configs.size(); i++) {
    config = configs[i];
    scales_.push_back(1.0f);//config.scale);
    theta = config.theta;
    tx = config.cx;
    ty = config.cy;
    pose_comp_root = CreatePose2D(tx, ty, theta);
    pose_comp_world = RootFrameToWorldFrame(pose_comp_root);

    SDF* s = new SDF(config.sdf_filename, pose_comp_world, scales_[i]);
    sdfs_.push_back(s);
    if (!s->Initialized()) {
      initialized_ = false;
    }
    sdf_coll_caches_.push_back(new SDFCollisionCache());
  }
}

void Mesh::InitializeMeshObject(std::vector<ComponentConfig> configs)
{
  ComponentConfig config;
  float theta = config.theta;
  float tx = config.cx;
  float ty = config.cy;
  Eigen::Matrix4f relative_pose;

  for (unsigned int i = 0; i < configs.size(); i++) {
    config = configs[i];
    theta = config.theta;
    tx = config.cx;
    ty = config.cy;
    
    mesh_objects_.push_back(new MeshObject(config.obj_filename.c_str(), scales_[i]));
    relative_pose = CreatePose(tx, ty, theta);
    poses_comp_root_.push_back(relative_pose);
  }
}

void Mesh::InitializeFCLObject(std::vector<ComponentConfig> configs)
{
  ComponentConfig config;
  fcl_vertices_.resize(configs.size());
  fcl_triangles_.resize(configs.size());
  
  for (unsigned int i = 0; i < configs.size(); i++) {
    config = configs[i];
    if (!LoadOBJFile(config.obj_filename.c_str(), fcl_vertices_[i], fcl_triangles_[i])) {
      initialized_ = false;
    }
  }
  CenterFCLVertices();
}

void Mesh::InitializeRendering(bool update_tris)
{
  if (update_tris) {
    render_triangles_.clear();
    render_triangles_.resize(NumComponents());
  }
  if (pd_poses_.size() == 0)
    pd_poses_.resize(NumComponents());

  for (unsigned int k = 0; k < NumComponents(); k++) {    
    // create pd pose
    Eigen::Matrix4f pose = RootFrameToWorldFrame(poses_comp_root_[k]);

    Transform pd_pose = Pose3DToPolyDepth(pose);
    pd_poses_[k] = pd_pose;

    // use pd pose to create rendering triangles
    if (update_tris) {
      render_triangles_[k].clear();
      int num_tris = mesh_objects_[k]->triangles().size();
      for (unsigned int i = 0; i < num_tris; i += 3) {
        Point_3 points[3];
        for (unsigned int j = 0; j < 3; j++) {
          int ind = mesh_objects_[k]->triangles()[i+j];
          Coord3D v = mesh_objects_[k]->vertices()[ind];
          v = pd_pose * v;
          points[j] = Point_3(v.X(), v.Y(), v.Z()); 
        }
        Triangle_3 local_t(points[0], points[1], points[2]);
        render_triangles_[k].push_back(local_t);
      }
    }
  }

  tris_need_update_ = (!update_tris);
}

void Mesh::InitializeMass()
{
  // compute mass using box2D to make things consistent with simulation 
  b2Vec2 position(0.0f, 0.0f);
  float angle = 0.0f;
  b2World* world = new b2World(position);
  b2Body* object_body;
  b2FixtureDef object_fixture_def;

  // create dynamic circle
  b2Separator object_sep;
  b2BodyDef body_def;
  body_def.type = b2_dynamicBody;
  body_def.position.Set(0.0f, 0.0f);
  body_def.angle = 0.0f;
  object_body = world->CreateBody(&body_def);

  Box2DShapeFactory::CreateNonconvexBody(this, Box2DShapeFactory::default_density, Box2DShapeFactory::default_friction, world,
                                         object_body, object_fixture_def, object_sep, position, angle);
  mass_ = object_body->GetMass();
  moment_of_inertia_ = object_body->GetInertia();
}

void Mesh::CenterFCLVertices()
{
  centroids_.resize(fcl_vertices_.size());
  for (unsigned int i = 0; i < fcl_vertices_.size(); i++) {
    // compute centroid
    centroids_[i] = fcl::Vec3f(0, 0, 0);
    for (unsigned int j = 0; j < fcl_vertices_[i].size(); j++) {
      centroids_[i] += fcl_vertices_[i][j];
    }
    centroids_[i] = centroids_[i] / fcl_vertices_[i].size();

    // subtract out centroid
    for (unsigned int j = 0; j < fcl_vertices_[i].size(); j++) {
      fcl_vertices_[i][j] -= centroids_[i];
    }
  }
}

void Mesh::ConvexDecomposition2D(float extrusion)
{
  convex_pieces_2d_.clear();
  convex_pieces_2d_.resize(NumComponents());

  // first 2D convex decomposition of the vertices  
  b2Separator object_sep;
  std::vector< std::vector< std::vector<b2Vec2> > > all_convex_pieces_2d;
  Box2DShapeFactory::ConvexDecomposition(this, object_sep, all_convex_pieces_2d);

  // add convex pieces for each component
  for (unsigned int k = 0; k < NumComponents(); k++) {
    LOG(INFO) << "Adding convex pieces for component " << k;
    std::vector< std::vector<b2Vec2> > convex_pieces_2d = all_convex_pieces_2d[k];

    // set up the goddamn polygons
    convex_pieces_2d_[k].resize(convex_pieces_2d.size());
    for (unsigned int i = 0; i < convex_pieces_2d.size(); i++) {

      // add vertices
      for (unsigned int j = 0; j < convex_pieces_2d[i].size(); j++) {
        // front face
        Eigen::Vector3f v;
        v << convex_pieces_2d[i][j].x, convex_pieces_2d[i][j].y, extrusion;
        convex_pieces_2d_[k][i].vertices.push_back(v);

        // back face
        Eigen::Vector3f w;
        w << convex_pieces_2d[i][j].x, convex_pieces_2d[i][j].y, -extrusion;
        convex_pieces_2d_[k][i].vertices.push_back(w);
      }

      // add pose
      convex_pieces_2d_[k][i].pose = RootFrameToWorldFrame(poses_comp_root_[k]);
    }
  }
}

Eigen::Matrix3f Mesh::RootFrameToWorldFrame(Eigen::Matrix3f pose_root_frame)
{  
  return Pose3DTo2D(pose_root_world_) * pose_root_frame;
}

Eigen::Matrix3f Mesh::WorldFrameToRootFrame(Eigen::Matrix3f pose_world_frame)
{
  return Pose3DTo2D(pose_root_world_.inverse()) * pose_world_frame;
}

Eigen::Matrix4f Mesh::RootFrameToWorldFrame(Eigen::Matrix4f pose_root_frame)
{
  return pose_root_world_ * pose_root_frame;
}

Eigen::Matrix4f Mesh::WorldFrameToRootFrame(Eigen::Matrix4f pose_world_frame)
{
  return pose_root_world_.inverse() * pose_world_frame;
}

Eigen::Matrix3f Mesh::Pose3DTo2D(Eigen::Matrix4f pose)
{
  // negatives for opposite y axis in image basis
  float tx, ty, theta;
  Pose3DToParams2D(pose, tx, ty, theta);
  Eigen::Matrix3f pose_2d = CreatePose2D(tx, ty, theta);
  return pose_2d;
}

Transform Mesh::Pose3DToPolyDepth(Eigen::Matrix4f pose)
{
  Transform pd_pose;
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      pd_pose.Rotation()[i][j] = pose(i,j);
    }
  }
  pd_pose.Translation()[0] = pose(0,3);
  pd_pose.Translation()[1] = pose(1,3);
  pd_pose.Translation()[2] = pose(2,3);
  return pd_pose;
}

Eigen::Matrix4f Mesh::Pose()
{
  return pose_root_world_;
}

// Adds the mesh boundary vertices in clockwise order
std::vector<std::vector<fcl::Vec3f> > Mesh::BoundaryVertices()
{
  std::vector<std::vector<fcl::Vec3f> > boundary_verts;
  boundary_verts.resize(NumComponents());

  // iterate through the triangles to find boundary
  for (unsigned int i = 0; i < fcl_triangles_.size(); i++) {

    std::map<unsigned int, std::vector<unsigned int> > boundary_edge_map;
    std::vector<MeshEdge> boundary_edges;

    // add all vertices on the front boundary
    for (unsigned int j = 0; j < fcl_triangles_[i].size(); j++) {
      fcl::Triangle t = fcl_triangles_[i][j];
      std::vector<int> front_indices;

      // NOTE: assumes 2D extruded polygon, will need modification in the future 
      for (unsigned int k = 0; k < 3; k++) {
        fcl::Vec3f vert = fcl_vertices_[i][t[k]];
        if (vert[2] > 0) {
          front_indices.push_back(t[k]);
        }
      }      

      // add vertices to pool if exactly two frontal vertices
      if (front_indices.size() == 2) {
        MeshEdge e;
        e.i = front_indices[0];
        e.j = front_indices[1];
        boundary_edges.push_back(e);

        unsigned int edge_ind = boundary_edges.size()-1;
        boundary_edge_map[e.i].push_back(edge_ind);
        boundary_edge_map[e.j].push_back(edge_ind);
      }
    }

    // go around the boundary, adding vertices in clockwise order
    std::map<unsigned int, std::vector<unsigned int> >::iterator it = boundary_edge_map.begin();
    std::vector<bool> edge_visited(boundary_edges.size(), false);
    unsigned int first_vert_ind = it->first;
    unsigned int cur_vert_ind = first_vert_ind;
    unsigned int prev_vert_ind = first_vert_ind;
    std::vector<unsigned int> cur_edges = it->second;
    boundary_verts[i].push_back(fcl_vertices_[i][cur_vert_ind]);                                   
    
    // assume that it forms a loop
    do {
      // get next edge
      bool found_next = false;
      for (unsigned int k = cur_edges.size() - 1; k >= 0 && !found_next; k--) { // index backwards for CW
        unsigned int edge_ind = cur_edges[k];
        unsigned int next_vert_ind;

        // next edge is next unvisited
        if (!edge_visited[edge_ind]) {
          found_next = true;

          // update cur vertex (assign current to opposite from edge)
          if (cur_vert_ind == boundary_edges[edge_ind].i) {
            next_vert_ind = boundary_edges[edge_ind].j;
          }
          else if (cur_vert_ind == boundary_edges[edge_ind].j) {
            next_vert_ind = boundary_edges[edge_ind].i;
          }

          // update edges
          if (next_vert_ind != first_vert_ind) {
            cur_edges = boundary_edge_map[next_vert_ind];
            edge_visited[edge_ind] = true;

            // check if on line segment
            fcl::Vec3f cur_v = fcl_vertices_[i][next_vert_ind]; 
            fcl::Vec3f prev_v = fcl_vertices_[i][cur_vert_ind]; 
            fcl::Vec3f prev_prev_v = fcl_vertices_[i][prev_vert_ind]; 
            Eigen::Matrix3f points;
            points << prev_prev_v[0], prev_prev_v[1], 1,
              prev_v[0], prev_v[1], 1,
              cur_v[0], cur_v[1], 1;

            // add if not on same line segment
            if (fabs(points.determinant()) > 1e-2) {
              prev_vert_ind = cur_vert_ind;
              boundary_verts[i].push_back(fcl_vertices_[i][cur_vert_ind]);                                   
            }
          }
          cur_vert_ind = next_vert_ind;
        }
      }
    } while (cur_vert_ind != first_vert_ind);

    // check CW, CCW
    float edge_sum = 0;
    for (unsigned int k = 0; k < boundary_verts[i].size(); k++) {
      if (k < boundary_verts[i].size()- 1) {
        edge_sum += (boundary_verts[i][k+1][0] - boundary_verts[i][k][0]) * (boundary_verts[i][k+1][1] + boundary_verts[i][k][1]);
      }
      else {
        edge_sum += (boundary_verts[i][0][0] - boundary_verts[i][k][0]) * (boundary_verts[i][0][1] + boundary_verts[i][k][1]);
      }
    }

    // if CCW, reverse order
    if (edge_sum > 0) {
      std::reverse(boundary_verts[i].begin(), boundary_verts[i].end());
    }
  }
  return boundary_verts;
}

void Mesh::ConcaveVertices(std::vector<ASE::Vertex2d>& verts, std::vector<std::pair<ASE::Direction2d, ASE::Direction2d> >& directions)
{
  // boundary vertices
  std::vector<std::vector<fcl::Vec3f> > boundary_vertices = BoundaryVertices();
  verts.clear();
  directions.clear();

  // previous indices
  int prev_i;
  int next_i;

  // edges and vertices
  ASE::Vertex2d v;
  ASE::Vertex2d prev_v;
  ASE::Vertex2d next_v;
  fcl::Vec3f fv;
  fcl::Vec3f prev_fv;
  fcl::Vec3f next_fv;
  ASE::Direction3d edge1;
  ASE::Direction3d edge2;
  ASE::Direction3d dir1;
  ASE::Direction3d dir2;
  ASE::Direction2d d1;
  ASE::Direction2d d2;

  // z down
  ASE::Direction3d n;
  ASE::Direction3d z;
  z << 0, 0, -1.0f;

  // loop through components
  for (unsigned int j = 0; j < boundary_vertices.size(); j++) {

    int num_verts = boundary_vertices[j].size(); 
    if (num_verts == 3)
      continue;

    for (int i = 0; i < num_verts; i++) {
      // get vertices in ordere going clockwise
      fv = boundary_vertices[j][i];
      v << fv[0], fv[1];

      prev_i = (i - 1) % num_verts;
      prev_fv = boundary_vertices[j][prev_i];
      prev_v << prev_fv[0], prev_fv[1];

      next_i = (i + 1) % num_verts;
      next_fv = boundary_vertices[j][next_i];
      next_v << next_fv[0], next_fv[1];    

      // check concavity via cross product trick
      edge1 << v(0) - prev_v(0), v(1) - prev_v(1), 0;
      edge2 << next_v(0) - v(0), next_v(1) - v(1), 0;
      n = edge1.cross(edge2);

      // concave, so get directions
      if (n(2) < 0) {
        dir1 = edge1.cross(z);
        dir2 = edge2.cross(z);
        d1 << dir1(0), dir1(1);
        d2 << dir2(0), dir2(1);
        d1 = d1 / d1.norm();
        d2 = d2 / d2.norm();
        
        // add back in the centroids
        verts.push_back(v);
        directions.push_back(std::make_pair<ASE::Direction2d, ASE::Direction2d>(d1, d2));
      } 
    }
  }
}

void Mesh::BoundaryEdges(std::vector<ASE::Vertex2d>& boundary_edges, std::vector<ASE::Direction2d>& directions)
{
  std::vector<std::vector<fcl::Vec3f> > boundary_vertices = BoundaryVertices();
  boundary_edges.clear();
  directions.clear();

  // previous indices
  int prev_i;
  int next_i;

  // edges and vertices
  ASE::Vertex2d v;
  ASE::Vertex2d next_v;
  fcl::Vec3f fv;
  fcl::Vec3f next_fv;
  ASE::Direction3d edge;
  ASE::Direction3d dir;
  ASE::Direction2d d;

  // z down
  ASE::Direction3d z;
  z << 0, 0, -1.0f;

  // loop through components
  for (unsigned int j = 0; j < boundary_vertices.size(); j++) {

    int num_verts = boundary_vertices[j].size(); 
    if (num_verts == 3)
      continue;

    for (int i = 0; i < num_verts; i++) {
      // get vertices in ordere going clockwise
      fv = boundary_vertices[j][i];
      v << fv[0], fv[1];

      next_i = (i + 1) % num_verts;
      next_fv = boundary_vertices[j][next_i];
      next_v << next_fv[0], next_fv[1];    

      // check concavity via cross product trick
      edge << next_v(0) - v(0), next_v(1) - v(1), 0;

      // concave, so get directions
      dir = edge.cross(z);
      d << dir(0), dir(1);
      d = d / d.norm();
      v = (v + next_v) / 2.0f;

      // add back in the centroids
      boundary_edges.push_back(v);
      directions.push_back(d);
    }
  }
}

float Mesh::MaxMomentArm()
{
  // take max over the individual comopnents
  // TODO: technically should be max moment from all points combined but for now its ok bc object is single components
  float max_moment_arm = 0.0f;
  for (unsigned int i = 0; i < NumComponents(); i++) {
    for (unsigned int j = 0; j < fcl_vertices_[i].size(); j++) {
      fcl::Vec3f v = fcl_vertices_[i][j];
      // assume centered vertices and don't use the z coordinate
      float m_i = sqrt(v[0]*v[0] + v[1]*v[1]);
      if (m_i > max_moment_arm) {
        max_moment_arm = m_i;
      }
    }
  }
  return max_moment_arm;
}

// pose of the object in world frame
void Mesh::SetPose(Eigen::Matrix4f pose, bool update_tris, bool print)
{
  g_p = false;
  pose_root_world_ = pose;
  float tx, ty, theta;
  // std::cout << "Mesh" << std::endl;
  //  std::cout << pose_root_world_ << std::endl;

  // create pd pose
  for (unsigned int k = 0; k < NumComponents(); k++) {
    Eigen::Matrix4f pose = RootFrameToWorldFrame(poses_comp_root_[k]);
    double transform[16];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 4; j++) {
        transform[i+4*j] = pose(i,j);
      } 
    } 
    for (int j = 0; j < 3; j++) {
      transform[3+4*j] = 0.0f;
    }
    transform[15] = 1.0f;

    if (print) {
      std::cout << "Pose world" << std::endl;
      std::cout << pose_root_world_ << std::endl;
      g_p = true;
    }
    // std::cout << "Tf" << std::endl;
    // for (unsigned int a = 0; a < 4; a++) {
    //   for (unsigned int b = 0; b < 4; b++)
    //     std::cout << transform[a+4*b] << " ";
    //   std::cout << std::endl;
    // }

    // negative because x and y axes are reverse for SDF coords
    Eigen::Matrix3f pose_2d = Pose3DTo2D(pose);
    sdfs_[k]->SetPose(pose_2d);

    // update collision cache
    sdfs_[k]->RenderToWorldFrame(sdf_coll_caches_[k]->prerendered_sdf, false);

    // update convex pieces
    for (unsigned int i = 0; i < convex_pieces_2d_[k].size(); i++) {
      convex_pieces_2d_[k][i].pose = pose;
    }
  }

  InitializeRendering(update_tris);
}

bool Mesh::SetRelativePoses(std::vector<ComponentConfig> configs)
{
  if (configs.size() != NumComponents()) {
    return false;
  }

  ComponentConfig config;
  float theta;
  float tx;
  float ty;
  Eigen::Matrix3f pose_comp_root_2d; // pose of component wrt root

  for (unsigned int k = 0; k < NumComponents(); k++) {
    config = configs[k];
    theta = config.theta;
    tx = config.cx;
    ty = config.cy;

    pose_comp_root_2d = CreatePose2D(tx, ty, theta);
    poses_comp_root_[k] = CreatePose(tx, ty, theta);
    sdfs_[k]->SetPose(pose_comp_root_2d);
  //   std::cout << "Pose " << k << "\n" << pose_comp_root_2d << std::endl;
  }

  InitializeRendering(true);
  return true;
}

std::vector<SDF*> Mesh::SDFs()
{
  return sdfs_;
}

std::vector<SDFCollisionCache*> Mesh::SDFCollCaches()
{
  return sdf_coll_caches_;
}

std::vector<MeshObject*> Mesh::PolyDepthObjects()
{
  return mesh_objects_;
}

std::vector<Transform> Mesh::PolyDepthPoses()
{
  return pd_poses_;
}

std::vector<Eigen::Matrix4f> Mesh::PosesWorldFrame()
{
  std::vector<Eigen::Matrix4f> poses;
  for (unsigned int i = 0; i < NumComponents(); i++) {
    poses.push_back(RootFrameToWorldFrame(poses_comp_root_[i]));
  }
  return poses;
}

std::vector<Eigen::Matrix4f> Mesh::PosesRootFrame()
{
  return poses_comp_root_;
}

unsigned int Mesh::NumComponents()
{
  return sdfs_.size();
}

std::vector<std::vector<fcl::Vec3f> > Mesh::FCLVertices()
{
  return fcl_vertices_;
}

std::vector<std::vector<fcl::Triangle> > Mesh::FCLTriangles()
{
  return fcl_triangles_;
}

float Mesh::Mass()
{
  return mass_;
}

float Mesh::MomentOfInertia()
{
  return moment_of_inertia_;
}

std::vector< std::vector<Poly> > Mesh::ConvexPieces2D()
{
  return convex_pieces_2d_;
}

void Mesh::RenderToGeomview(CGAL::Geomview_stream& gv)
{
  if (tris_need_update_) {
    InitializeRendering(true);
  }

  for (unsigned int i = 0; i < NumComponents(); i++) {
    gv.draw_triangles(render_triangles_[i].begin(), render_triangles_[i].end());
  }
}

void Mesh::RenderToImage(cv::Mat& image, char color)
{
  for (unsigned int i = 0; i < NumComponents(); i++) {
    sdfs_[i]->RenderToRGBImage(image, color);
  }
}

void Mesh::RenderTrisToImage(cv::Mat& image, float scale, char color, bool draw_com)
{
  std::vector<std::vector<cv::Point2i> > contours;
  std::vector<Eigen::Matrix4f> poses = PosesWorldFrame();

  int num_points = 3;
  cv::Size image_size = image.size();
  Eigen::Matrix4f T_image_to_image_center;
  T_image_to_image_center << 1, 0, 0, image_size.width / 2.0f,
    0, 1, 0, (float)image_size.height / 2.0f,
    0, 0, 1, 0,
    0, 0, 0, 1;
  Eigen::Matrix4f T_object_rot_to_object;
  T_object_rot_to_object << 1, 0, 0, 0,
    0, -1, 0, 0,
    0, 0, -1, 0,
    0, 0, 0, 1;

  std::vector<Eigen::Matrix4f> centers;

  for (unsigned int i = 0; i < NumComponents(); i++) {
    Eigen::Matrix4f pose = poses[i];
    MeshObject* mesh_object = mesh_objects_[i];
    Eigen::Matrix4f T_image_to_object;
    T_image_to_object << 1, 0, 0, 0, //-mesh_object->centroid()[0],
      0, -1, 0, 0, //-mesh_object->centroid()[1],
      0, 0, -1, 0,
      0, 0, 0, 1;
    T_image_to_object = T_object_rot_to_object * T_image_to_object;

    // get color
    cv::Scalar act_color(255, 0, 0);
    if (color == 'g') {
      act_color = cv::Scalar(0, 255, 0);
    }
    else if (color == 'r') {
      act_color = cv::Scalar(0, 0, 255);
    }


    // draw center of mass
    if (draw_com) {
      Eigen::Matrix4f v;
      v <<
        1, -1, 0, 0,
        0, 0, 1, -1,
        0, 0, 0, 0,
        1, 1, 1, 1; // assume vertices are locally centered
      v = T_image_to_object * v; // convert from image to object
      v = pose * v;
      v.topLeftCorner(3,4) = scale * v.topLeftCorner(3,4);
      v = T_image_to_image_center * T_object_rot_to_object.inverse() * v; // convert from object to image
      centers.push_back(v);
    }

    cv::Point2i* tris = new cv::Point2i[num_points];

    for (unsigned int j = 0; j < fcl_triangles_[i].size(); j++) {
      fcl::Triangle t = fcl_triangles_[i][j];

      contours.push_back(std::vector<cv::Point2i>());
      unsigned int a = contours.size() - 1;

      for (unsigned int k = 0; k < 3; k++) {
        fcl::Vec3f vert = fcl_vertices_[i][t[k]];
        Eigen::Vector4f v;
        v << vert[0], vert[1], 0, 1;
        v = T_image_to_object * v; // convert from image to object
        v = pose * v;
        v(0) = scale * v(0); // scale verts
        v(1) = scale * v(1);
        v = T_image_to_image_center * T_object_rot_to_object.inverse() * v; // convert from object to image
        contours[a].push_back(cv::Point2i((int)v(0), (int)v(1)));
        tris[k] = cv::Point2i((int)v(0), (int)v(1));
      } 

      cv::fillConvexPoly(image, (const cv::Point2i*)tris, num_points, act_color, CV_AA);
    }

    delete tris;
  }

  // draw center of mass circle
  int radius = 20;
  int thickness = 10;
  cv::Scalar center_color(0, 0, 255);
  if (draw_com) {
    for (unsigned int i = 0; i < NumComponents(); i++) {
      cv::Point2i pt_up(centers[i](0,0), centers[i](1,0));
      cv::Point2i pt_down(centers[i](0,1), centers[i](1,1));
      cv::Point2i pt_left(centers[i](0,2), centers[i](1,2));
      cv::Point2i pt_right(centers[i](0,3), centers[i](1,3));

      cv::line(image, pt_up, pt_down, center_color, thickness);
      cv::line(image, pt_left, pt_right, center_color, thickness);
      //      cv::circle(image, centers[i], radius, center_color, thickness);
    }
  }


  // draw contours with color
  // int idx = -1;
  // for (int i = 0; i < contours.size(); i++) {
  //   if (contours[i].size() != 3) {
  //     std::cout << "Contours fkd at " << i << ": " << contours[i].size() << std::endl;
  //   }
  // }
  // cv::drawContours(image, contours, idx, act_color, CV_FILLED);
}

MeshCollisionResult Mesh::LowerBoundCollision(Mesh* mesh1, Mesh* mesh2, bool print)
{
  bool collision = true;
  Timer timer;
  timer.start();

  float pen_depth = 0.0f;
  float pd;
  int coll_i = -1;
  int coll_j = -1;  
  int coll_k = -1;
  int coll_l = -1;  
  std::vector<std::vector<Poly> > convex_pieces1 = mesh1->ConvexPieces2D();
  std::vector<std::vector<Poly> > convex_pieces2 = mesh2->ConvexPieces2D();

  // loop through components to find max penetration
  for (unsigned int i = 0; i < mesh1->NumComponents(); i++) {
    for (unsigned int j = 0; j < mesh2->NumComponents(); j++) {

      // loop through all convex piece combinations 
      for (unsigned int k = 0; k < convex_pieces1[i].size(); k++) {
        for (unsigned int l = 0; l < convex_pieces2[j].size(); l++) {

          if (print) {
            std::cout << "Indices " << i << " " << j << " " << k << " " << l << std::endl;
          }
          
          ccd_t ccd;
          CCD_INIT(&ccd);
          ccd.support1 = support;
          ccd.support2 = support;
          ccd.max_iterations = 100;
          ccd.epa_tolerance  = 1e-6;
          ccd.dist_tolerance  = 0.0f;
          ccd.mpr_tolerance  = 1e-6;
          ccd_real_t depth = 0.0;
          ccd_vec3_t dir;
          ccd_vec3_t pos;
          int intersect = ccdGJKPenetration((const void*)(&convex_pieces1[i][k]), 
                                            (const void*)(&convex_pieces2[j][l]),
                                            &ccd, &depth, &dir, &pos);
          bool coll = (intersect == 0);
          pd = (float)depth;

          if (print) {
            std::cout << "Collision " << coll << std::endl;
            std::cout << "Pen Depth " << pd << std::endl;
            // std::cout << "Dir " << dir.v[0] << " " << dir.v[1] << " " << dir.v[2] << std::endl;
            // std::cout << "Pose k " << convex_pieces1[i][k].pose << std::endl;
            // std::cout << "Pose l " << convex_pieces2[j][l].pose << std::endl;
          }

          // update penetration depth if larger
          if (coll && pd > pen_depth) {
            pen_depth = pd;
            coll_i = i;
            coll_j = j;
            coll_k = k;
            coll_l = l;
          }
        }
      }
    }
  }
  //  while(true);
  if (print) {
    std::cout << "Collision " << collision << std::endl;
    std::cout << "Distance " << pen_depth << std::endl;
    std::cout << "Coll I " << coll_i << std::endl;
    std::cout << "Coll J " << coll_j << std::endl;
    std::cout << "Coll K " << coll_k << std::endl;
    std::cout << "Coll L " << coll_l << std::endl;
  }

  // if pen depth is 0 then we calculate the closest distance
  if (pen_depth == 0.0f) {
    collision = false;
    pen_depth = 0.0f;
    //    return Distance(mesh1, mesh2);
  }
  timer.stop();
  float time_elapsed = timer.time();
  
  MeshCollisionResult result;
  result.collision = collision;
  result.distance = pen_depth;
  result.time = time_elapsed;
  result.coll_i = coll_i;
  result.coll_j = coll_j;
  return result;
}

// returns a upper bound (tighter than LB) on the penetration depth if a collision occurs
MeshCollisionResult Mesh::UpperBoundCollision(Mesh* mesh1, Mesh* mesh2, bool compute_separation)
{
  bool use_motion_coherance = true;
  bool use_maximal_clear_conf = false;
  Tri *t1;
  Tri *t2;

  std::vector<Transform> trans1 = mesh1->PolyDepthPoses();
  std::vector<Transform> trans2 = mesh2->PolyDepthPoses();	
  std::vector<MeshObject*> objects1 = mesh1->PolyDepthObjects();
  std::vector<MeshObject*> objects2 = mesh2->PolyDepthObjects();
  float pen_depth = 0.0f;
  float local_pd;
  bool collision = false;
  int coll_i = -1;
  int coll_j = -1;

  Timer timer;
  timer.start();

  for (unsigned int i = 0; i < mesh1->NumComponents(); i++) {
    for (unsigned int j = 0; j < mesh2->NumComponents(); j++) {
      t1 = objects1[i]->c2a_model()->last_tri;
      t2 = objects2[i]->c2a_model()->last_tri;

      std::vector<Coord3D> local_penetration_depth;
      std::vector<Coord3D> local_penetration_features1;
      std::vector<Coord3D> local_penetration_features2;
      Coord3D global_penetration_depth;
      Transform trans1_after;
      Transform trans2_after;
      int n_itr = 0;
      PolyDepthReturnValue result = PolyDepth(&trans1[i], objects1[i]->c2a_model(), 
                                              &trans2[j], objects2[j]->c2a_model(), 
                                              t1,t2,
                                              objects1[i]->vertices(), objects2[j]->vertices(),
                                              trans1_after, trans2_after, 
                                              local_penetration_depth,
                                              local_penetration_features1,
                                              local_penetration_features2,
                                              global_penetration_depth,
                                              n_itr,
                                              objects1[i]->clear_conf(),
                                              use_motion_coherance,
                                              use_maximal_clear_conf);

      bool local_collision = (result.result_ == kPenetration);
      local_pd = sqrt(global_penetration_depth.X() * global_penetration_depth.X() +
                global_penetration_depth.Y() * global_penetration_depth.Y() + 
                global_penetration_depth.Z() * global_penetration_depth.Z());
      if (local_collision && local_pd > pen_depth) {
        collision = true;
        pen_depth = local_pd;
        coll_i = i;
        coll_j = j;
      }
    }
  }
  timer.stop();
  float time_elapsed = timer.time();

  if (!collision && compute_separation) {
    return Distance(mesh1, mesh2);
  }

  MeshCollisionResult collision_result;
  collision_result.collision = collision;
  collision_result.distance = pen_depth;
  collision_result.time = time_elapsed;
  collision_result.coll_i = coll_i;
  collision_result.coll_j = coll_j;
  return collision_result;
}

// returns the minimum distance between two meshes
MeshCollisionResult Mesh::Distance(Mesh* mesh1, Mesh* mesh2)
{
  float min_distance = FLT_MAX;
  int coll_i = -1;
  int coll_j = -1;

  Timer timer;
  timer.start();

  std::vector<std::vector<fcl::Vec3f> > vertices1 = mesh1->FCLVertices();
  std::vector<std::vector<fcl::Triangle> > triangles1 = mesh1->FCLTriangles();
  std::vector<std::vector<fcl::Vec3f> > vertices2 = mesh2->FCLVertices();
  std::vector<std::vector<fcl::Triangle> > triangles2 = mesh2->FCLTriangles();
  std::vector<Eigen::Matrix4f> poses1 = mesh1->PosesWorldFrame();
  std::vector<Eigen::Matrix4f> poses2 = mesh2->PosesWorldFrame();

  for (unsigned int i = 0; i < mesh1->NumComponents(); i++) {
    for (unsigned int j = 0; j < mesh2->NumComponents(); j++) {

      // set up the fcl geometries
      fcl::BVHModel<fcl::RSS> m1;
      m1.bv_splitter.reset(new fcl::BVSplitter<fcl::RSS>(fcl::SPLIT_METHOD_MEAN));
      fcl::BVHModel<fcl::RSS> m2;
      m2.bv_splitter.reset(new fcl::BVSplitter<fcl::RSS>(fcl::SPLIT_METHOD_MEAN));

      m1.beginModel();
      m1.addSubModel(vertices1[i], triangles1[i]);
      m1.endModel();

      m2.beginModel();
      m2.addSubModel(vertices2[j], triangles2[j]);
      m2.endModel();

      // convert poses to FCL
      Eigen::Matrix4f pose1 = poses1[i];
      fcl::Matrix3f R1(pose1(0,0), pose1(0,1), pose1(0,2),
                       pose1(1,0), pose1(1,1), pose1(1,2),
                       pose1(2,0), pose1(2,1), pose1(2,2));
      fcl::Vec3f t1(pose1(0,3), pose1(1,3), pose1(2,3));
      const fcl::Transform3f fcl_pose1(R1, t1);

      Eigen::Matrix4f pose2 = poses2[j];
      fcl::Matrix3f R2(pose2(0,0), pose2(0,1), pose2(0,2),
                       pose2(1,0), pose2(1,1), pose2(1,2),
                       pose2(2,0), pose2(2,1), pose2(2,2));
      fcl::Vec3f t2(pose2(0,3), pose2(1,3), pose2(2,3));
      const fcl::Transform3f fcl_pose2(R2, t2);

      // request distance
      int qsize = 2;
      fcl::DistanceResult result;
      fcl::MeshDistanceTraversalNodeRSS node;
      if(!fcl::initialize(node, (const fcl::BVHModel<fcl::RSS>&)m1, fcl_pose1, (const fcl::BVHModel<fcl::RSS>&)m2, fcl_pose2,
                          fcl::DistanceRequest(true), result)) {
        std::cout << "Error initializing FCL distance request" << std::endl;
      }
      fcl::distance(&node, NULL, qsize);

      // update min distance
      if (result.min_distance < min_distance) {
        min_distance = result.min_distance; 
        coll_i = i;
        coll_j = j;
      }
    }
  }
  timer.stop();
  // std::cout << "Min " << min_distance << std::endl;
  // std::cout << "Comp1 " << mesh1->NumComponents() << std::endl;
  // std::cout << "Comp2 " << mesh2->NumComponents() << std::endl;
  float time_elapsed = timer.time();
  bool collision = (min_distance <= 0);
  
  MeshCollisionResult collision_result;
  collision_result.collision = collision;
  collision_result.distance = min_distance;
  collision_result.time = time_elapsed;
  collision_result.coll_i = coll_i;
  collision_result.coll_j = coll_j;
  return collision_result;
}
