#include "GraspGenerator.hpp"

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <glog/logging.h>

unsigned int random_index(unsigned int max_ind)
{
  return (unsigned int)rand() % max_ind;
}

float random_angle(float theta, float sigma = 1.0f)
{
  boost::mt19937 rng;
  boost::normal_distribution<> nd(theta, sigma);
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_nor(rng, nd);
  return var_nor();  
}

ASE::Vertex2d perturb_point(ASE::Vertex2d point, float sigma = 1.0f)
{
  boost::mt19937 rng;
  boost::normal_distribution<> ndx(point(0), sigma);
  boost::normal_distribution<> ndy(point(1), sigma);
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_nor_x(rng, ndx);
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_nor_y(rng, ndy);

  ASE::Vertex2d perturbed_point;
  perturbed_point(0) = var_nor_x();
  perturbed_point(1) = var_nor_y();
  return perturbed_point;
}

GraspGenerator::GraspGenerator(Mesh* object)
  : object_(object)
{
}

ParallelJawGraspGenerator::ParallelJawGraspGenerator(Mesh* object)
  : GraspGenerator(object)
{
}

bool ParallelJawGraspGenerator::RandomGrasp(ParallelJawGrasp2D* grasp, GripperConfig& config, float max_grasp_width, float delta_jaw, float p_concave)
{
  float r = (float)(rand()) / (float)(RAND_MAX);
  if (r < p_concave) {
    LOG(INFO) << "Returning concave grasp";
    return ConcaveGrasp(grasp, config, max_grasp_width, delta_jaw);
  }
  LOG(INFO) << "Returning antipodal grasp";
  return AntipodalGrasp(grasp, config, max_grasp_width, delta_jaw);
}
 
bool ParallelJawGraspGenerator::ConcaveGrasp(ParallelJawGrasp2D* grasp, GripperConfig& config, float max_grasp_width, float delta_jaw)
{
  // set up variables, get concave vertices
  std::vector<ASE::Vertex2d> concave_vertices;
  std::vector<std::pair<ASE::Direction2d, ASE::Direction2d> > vertex_dirs;

  object_->ConcaveVertices(concave_vertices, vertex_dirs);
  unsigned int num_vertices = concave_vertices.size();

  if (concave_vertices.size() == 0) {
    return false;//AntipodalGrasp(grasp, config);
  }

  // get random vertex to grasp at
  unsigned int index = random_index(num_vertices);
  ASE::Vertex2d v1 = perturb_point(concave_vertices[index]);
  ASE::Direction2d d1 = vertex_dirs[index].first;
  ASE::Direction2d d2 = vertex_dirs[index].second;
  
  // sample the directions to step along
  float w = (float)rand() / (float)RAND_MAX;
  ASE::Direction2d dir = w * d1 + (1.0 - w) * d2;
  dir = dir / dir.norm();

  float width = max_grasp_width;
  ASE::Vertex2d g1 = v1 - delta_jaw * dir;
  ASE::Vertex2d g2 = g1 + width * dir;
  float dist_from_start = (g1 - v1).norm();

  // check collisions on initial jaw placement
  ParallelJawGrasp2D::EndpointsToGrasp(g1, g2, config);
  grasp->SetState(config);
  MeshCollisionResult c = Mesh::UpperBoundCollision(object_, grasp);

  // step first jaw backwards until out of collision
  while (c.collision && dist_from_start < max_grasp_width) {
    g1 = g1 - delta_jaw * dir;
    g2 = g2 + width * dir;
    ParallelJawGrasp2D::EndpointsToGrasp(g1, g2, config);
    grasp->SetState(config);

    // check collisions
    c = Mesh::UpperBoundCollision(object_, grasp);
  }

  // if still in collision, then we failed
  if (c.collision) {
    return false;
  }

  // close until collision found
  while (!c.collision && width > delta_jaw) {
    // shorten width
    width = width - delta_jaw;
    g2 = g1 + width * dir;
    ParallelJawGrasp2D::EndpointsToGrasp(g1, g2, config);
    grasp->SetState(config);

    // check collisions
    c = Mesh::UpperBoundCollision(object_, grasp);
  }

  // set the grasp to the previous state
  width = width + delta_jaw;
  g2 = g1 + width * dir;
  ParallelJawGrasp2D::EndpointsToGrasp(g1, g2, config);
  grasp->SetState(config);
  c = Mesh::UpperBoundCollision(object_, grasp);

  // report failures
  if (c.collision) {
    return false;
  }
  return true;
}

bool ParallelJawGraspGenerator::AntipodalGrasp(ParallelJawGrasp2D* grasp, GripperConfig& config, float max_grasp_width, float delta_jaw)
{
  std::vector<ASE::Vertex2d> boundary_edge_centers;
  std::vector<ASE::Direction2d> edge_dirs;

  object_->BoundaryEdges(boundary_edge_centers, edge_dirs);
  unsigned int num_vertices = boundary_edge_centers.size();

  if (boundary_edge_centers.size() == 0) {
    return false;
  }

  // get grasps
  unsigned int index = random_index(num_vertices);
  ASE::Vertex2d v1 = perturb_point(boundary_edge_centers[index]);
  ASE::Direction2d dir = edge_dirs[index];

  // random angle perturbation
  float theta = acos(dir(0));
  if (dir(1) < 0) {
    theta = 2 * M_PI - theta;  
  }
  theta = random_angle(theta);
  dir(0) = cos(theta);
  dir(1) = sin(theta);

  // initial grasp
  float width = max_grasp_width;
  ASE::Vertex2d g1 = v1 - delta_jaw * dir;
  ASE::Vertex2d g2 = g1 + width * dir;
  float dist_from_start = (g1 - v1).norm();

  // check collisions on initial jaw placement
  ParallelJawGrasp2D::EndpointsToGrasp(g1, g2, config);
  grasp->SetState(config);
  MeshCollisionResult c = Mesh::UpperBoundCollision(object_, grasp);

  // step first jaw backwards until out of collision
  while (c.collision && dist_from_start < max_grasp_width) {
    g1 = g1 - delta_jaw * dir;
    g2 = g2 + width * dir;
    ParallelJawGrasp2D::EndpointsToGrasp(g1, g2, config);
    grasp->SetState(config);

    // check collisions
    c = Mesh::UpperBoundCollision(object_, grasp);
  }

  // if still in collision, then we failed
  if (c.collision) {
    return false;
  }

  // close until collision found
  while (!c.collision && width > delta_jaw) {
    // shorten width
    width = width - delta_jaw;
    g2 = g1 + width * dir;
    ParallelJawGrasp2D::EndpointsToGrasp(g1, g2, config);
    grasp->SetState(config);

    // check collisions
    c = Mesh::UpperBoundCollision(object_, grasp);
  }

  // set the grasp to the previous state
  width = width + delta_jaw;
  g2 = g1 + width * dir;
  ParallelJawGrasp2D::EndpointsToGrasp(g1, g2, config);
  grasp->SetState(config);
  c = Mesh::UpperBoundCollision(object_, grasp);

  // report failures
  if (c.collision) {
    return false;
  }
  return true;
}
