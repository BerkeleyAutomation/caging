/*
 * Parses YAML configuration files for the caging experiments
 *
 */
#ifndef CAGE_CONFIG_H
#define CAGE_CONFIG_H

#include <yaml-cpp/yaml.h>
#include "Typedef.h"

#include "ShapeFactory.hpp"

enum GripperType {
  PARALLEL_JAW = 0,
  RIGID = 1,
  QUAD = 2,
  TRI
};

class CageConfigParser
{
 public:
  CageConfigParser(std::string filename);

 public:
  bool Parse();
  void Print();

 public:
  float Extrusion();
  float Bounding();
  bool Debug();
  bool RandomGrasps();
  bool CloseGrippers();
  unsigned int NumSamples();
  float SamplingScale();
  GripperType Gripper();
  float GripperRadius();
  std::string GripperPath(); // for composite built in gripper models
  std::vector<std::string> ObjectList();

 private:
  void ParseObjectFilenames(std::string list_file);

 private:
  std::string config_filename_;
  float extrusion_;
  float bounding_;
  bool debug_;
  bool random_grasps_;
  bool close_grippers_;
  unsigned int num_samples_;
  float sampling_scale_;
  GripperType gripper_type_;
  float gripper_radius_;
  std::string gripper_path_;
  std::vector<std::string> object_list_;
};

#endif // CAGE_CONFIG_H
