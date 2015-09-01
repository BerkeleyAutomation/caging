#pragma once

#include "Mesh.hpp"

struct GripperConfig : public ComponentConfig
{
  float width;
  float angle; // rotations of the space
};

void LoadGrasps(const std::string& filename, GripperConfig config_template, std::vector<GripperConfig>& configs);
void LoadGraspsAndFilenames(const std::string& filename, GripperConfig config_template, std::vector<GripperConfig>& configs);

class ParallelJawGrasp2D : public Mesh
{
 public:
  // object containing pose, width, and filename for obj and sdfs
  ParallelJawGrasp2D(GripperConfig config); 
  ~ParallelJawGrasp2D() {}

public:
  static void EndpointsToGrasp(ASE::Vertex2d g1, ASE::Vertex2d g2, GripperConfig& config);
  void SetState(GripperConfig config);

private:
  GripperConfig config_;
};

class ThreeFingerGrasp2D : public Mesh
{
 public:
  // object containing pose, width, and filename for obj and sdfs
  ThreeFingerGrasp2D(GripperConfig config); 
  ~ThreeFingerGrasp2D() {}

public:
  //  static void EndpointsToGrasp(ASE::Vertex2d g1, ASE::Vertex2d g2, GripperConfig& config);
  void SetState(GripperConfig config);

private:
  GripperConfig config_;
};

class OneFingerGrasp2D : public Mesh
{
 public:
  // object containing pose, width, and filename for obj and sdfs
  OneFingerGrasp2D(GripperConfig config); 
  ~OneFingerGrasp2D() {}

public:
  //  static void EndpointsToGrasp(ASE::Vertex2d g1, ASE::Vertex2d g2, GripperConfig& config);
  void SetState(GripperConfig config);

private:
  GripperConfig config_;
};
