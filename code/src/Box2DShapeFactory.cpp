#include "Simulation/Box2DShapeFactory.hpp"

#include <iostream>
#include "Util.h"

float Box2DShapeFactory::default_density = 0.01f;
float Box2DShapeFactory::default_friction = 0.0f;

bool Box2DShapeFactory::CreateNonconvexBody(Mesh* object, float density, float friction, b2World* world, b2Body* object_body,
                                            b2FixtureDef& object_fixture_def, b2Separator& object_sep,
                                            b2Vec2& object_start_pos, float& object_start_angle)
{
  std::vector<Eigen::Matrix4f> object_poses_3d = object->PosesWorldFrame();    
  float object_tx, object_ty, object_theta;
  Pose3DToParams2D(object_poses_3d[0], object_tx, object_ty, object_theta);

  object_start_pos.Set(object_ty, object_tx);
  object_start_angle = object_theta;

  // add vertices from object, obstacles
  std::vector<std::vector<fcl::Vec3f> > object_vertices_3d = object->BoundaryVertices();
  std::vector<b2Vec2> object_vertices_2d;
  if (object_vertices_3d.size() == 0) {
    std::cout << "Error: No object vertices" << std::endl;
    return false;
  }

  b2Vec2 object_com(0.0f, 0.0f);
  int num_v = 0;
  for (unsigned int i = 0; i < object_vertices_3d[0].size(); i++) {
    fcl::Vec3f v = object_vertices_3d[0][i];;

    // only add extruded fronts
    if (v[2] > 0) {
      object_vertices_2d.push_back(b2Vec2(v[0], v[1]));
      object_com += b2Vec2(v[0], v[1]);
      num_v++;
    }
  }
  object_com = (1.0f / num_v) * object_com;

  // convex decomposition of vertices for simulation
  float scale = 1.0f; // some random param for separator that seems to work
  object_sep.Separate(object_body, &object_fixture_def, &object_vertices_2d, scale);

  // set circle density and friction
  object_fixture_def.density = density;
  object_fixture_def.friction = friction;

  // add the shape to the body.
  object_body->CreateFixture(&object_fixture_def);

  // set object mass at true center (b2 separator messes it up badly)
  b2MassData mass_data;
  object_body->GetMassData(&mass_data);
  mass_data.center = object_com;
  object_body->SetMassData(&mass_data);  

  return true;
}

bool Box2DShapeFactory::ConvexDecomposition(Mesh* object, b2Separator& object_sep,
                                            std::vector<std::vector<std::vector<b2Vec2> > >& convex_pieces)
{
  std::vector<Eigen::Matrix4f> object_poses_3d = object->PosesWorldFrame();    

  // add vertices from object, obstacles
  std::vector<std::vector<fcl::Vec3f> > object_vertices_3d = object->BoundaryVertices();
  if (object_vertices_3d.size() == 0) {
    std::cout << "Error: No object vertices" << std::endl;
    return false;
  }

  convex_pieces.clear();
  convex_pieces.resize(object_vertices_3d.size());

  // convert to 2D
  for (unsigned int k = 0; k < object_vertices_3d.size(); k++) {
    std::vector<b2Vec2> object_vertices_2d;
    int num_v = 0;
    for (unsigned int i = 0; i < object_vertices_3d[k].size(); i++) {
      // get object vertex
      fcl::Vec3f v = object_vertices_3d[k][i];

      // only add extruded fronts
      if (v[2] > 0) {
        object_vertices_2d.push_back(b2Vec2(v[0], v[1]));
        num_v++;
      }
    }

    // convex decomposition of vertices for simulation
    object_sep.ConvexDecomp(object_vertices_2d, convex_pieces[k]);
  }

  return true;
}
