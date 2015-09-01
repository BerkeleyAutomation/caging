#pragma once

#include "Grasp.hpp"
#include "Mesh.hpp"
#include "ShapeFactory.hpp"

class GraspGenerator
{
 public:
  GraspGenerator(Mesh* object);

 protected:
  Mesh* object_;
};

class ParallelJawGraspGenerator : public GraspGenerator
{
 public:
  ParallelJawGraspGenerator(Mesh* object);
  
 public:
  bool RandomGrasp(ParallelJawGrasp2D* grasp, GripperConfig& config, float max_grasp_width = 10.0f, float delta_jaw = 1.0f, float p_concave = 0.5f);  
  bool ConcaveGrasp(ParallelJawGrasp2D* grasp, GripperConfig& config, float max_grasp_width = 10.0f, float delta_jaw = 1.0f);
  bool AntipodalGrasp(ParallelJawGrasp2D* grasp, GripperConfig& config, float max_grasp_width = 10.0f, float delta_jaw = 1.0f);
};
