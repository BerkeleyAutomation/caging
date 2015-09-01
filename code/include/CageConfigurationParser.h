/*
 * Parses YAML configuration files for the caging experiments
 *
 */
#ifndef CAGE_CONFIG_H
#define CAGE_CONFIG_H

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>
#include <sstream>
#include <vector>

#include <yaml-cpp/yaml.h>

#include "ShapeFactory.hpp"

enum GripperType {
  PARALLEL_JAW = 0,
  BARRETT_HAND
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
  MultiObjectConfig GripperConfig();
  MultiObjectConfig ObjectConfig();

 private:
  void ParseCompositeConfig(const YAML::Node& node, MultiObjectConfig& config);
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
  MultiObjectConfig gripper_config_;
  MultiObjectConfig object_config_;  
};

#endif // CAGE_CONFIG_H
