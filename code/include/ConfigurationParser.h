/*
 * Parses YAML configuration files for the caging experiments
 *
 * Now deprecated! This works only for SOLID library collision detection and whatnot
 */
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>
#include <sstream>
#include <vector>

#include <yaml-cpp/yaml.h>

#include "ShapeFactory.hpp"

class ConfigParser
{
 public:
  ConfigParser(std::string filename);

 public:
  bool Parse();
  void Print();

 public:
  float Extrusion();
  float Bounding();
  bool Debug();
  unsigned int NumSamples();
  CompositeObjectConfig GripperConfig();
  CompositeObjectConfig ObjectConfig();

 private:
  void ParseCompositeConfig(const YAML::Node& node, CompositeObjectConfig& config);

 private:
  std::string config_filename_;
  float extrusion_;
  float bounding_;
  bool debug_;
  unsigned int num_samples_;
  CompositeObjectConfig gripper_config_;
  CompositeObjectConfig object_config_;  
};
