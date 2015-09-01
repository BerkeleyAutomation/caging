#pragma once

// Small class to create box2d bodies
#include "Mesh.hpp"
#include "Simulation/Test.h"
#include "Simulation/Render.h"
#include <Box2D/Dynamics/b2Separator.h>

class Box2DShapeFactory
{
 public:
  static bool CreateNonconvexBody(Mesh* object, float density, float friction, b2World* world, b2Body* object_body,
                                  b2FixtureDef& object_fixture_def, b2Separator& object_sep, b2Vec2& object_start_pos, float& object_start_angle);

  static bool ConvexDecomposition(Mesh* object, b2Separator& object_sep,
                                  std::vector<std::vector<std::vector<b2Vec2> > >& convex_pieces);

public:
  static float default_density;
  static float default_friction;
};
