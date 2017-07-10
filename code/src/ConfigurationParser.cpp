#include "ConfigurationParser.h"

ConfigParser::ConfigParser(std::string filename) : config_filename_(filename)
{
  Parse();
}                                                         

// Parses the YAML config file
bool ConfigParser::Parse()
{
  YAML::Node node = YAML::LoadFile(config_filename_.c_str());

  // attempt to parse config
  // NOTE: assumes given keys are there
  try {
    // get constants
    extrusion_ = node["extrusion"].as<float>();
    bounding_ = node["bounding"].as<float>();
    debug_ = node["debug"].as<bool>();
    num_samples_ = node["samples"].as<unsigned int>();

    // get gripper subnode
    const YAML::Node& grip_subnode = node["gripper"];
    ParseCompositeConfig(grip_subnode, gripper_config_);

    // get object subnode
    const YAML::Node& obj_subnode = node["object"];
    ParseCompositeConfig(obj_subnode, object_config_);    
  }
  catch (YAML::ParserException& e) {
    std::cout << e.what() << std::endl;
    return false;
  }

  return true;
}

void ConfigParser::Print()
{
  std::cout << "CONFIGURATION VALUES" << std::endl;
  std::cout << "Extrusion: " << extrusion_ << std::endl;
  std::cout << "Bounding: " << bounding_ << std::endl;
  std::cout << "Debug: " << debug_ << std::endl;
  std::cout << "Num Samples: " << num_samples_ << std::endl;
  std::cout << "Num Gripper Boxes: " << gripper_config_.boxes.size() << std::endl;
  std::cout << "Num Gripper Tris: " << gripper_config_.tris.size() << std::endl;
  std::cout << "Num Object Boxes: " << object_config_.boxes.size() << std::endl;
  std::cout << "Num Object Tris: " << object_config_.tris.size() << std::endl;
  std::cout << std::endl;
}

void ConfigParser::ParseCompositeConfig(const YAML::Node& node, CompositeObjectConfig& comp_config)
{
  // parse boxes
  for (unsigned int i = 0; i < node["boxes"].size(); i++) {
    BoxConfig box;
    box.width = node["boxes"][i]["width"].as<float>();
    box.height = node["boxes"][i]["height"].as<float>();
    box.cx = node["boxes"][i]["cx"].as<float>();
    box.cy = node["boxes"][i]["cy"].as<float>();
    box.theta = node["boxes"][i]["theta"].as<float>();
    comp_config.boxes.push_back(box);
  }

  // parse triangles
  std::vector<TriangleConfig> tri_configs;
  for (unsigned int i = 0; i < node["triangles"].size(); i++) {
    TriangleConfig tri;
    tri.scale = node["triangles"][i]["scale"].as<float>();
    tri.cx = node["triangles"][i]["cx"].as<float>();
    tri.cy = node["triangles"][i]["cy"].as<float>();
    tri.theta = node["triangles"][i]["theta"].as<float>();
    comp_config.tris.push_back(tri);
  }
}

float ConfigParser::Extrusion()
{
  return extrusion_;
}

float ConfigParser::Bounding()
{
  return bounding_;
}

bool ConfigParser::Debug()
{
  return debug_;
}

unsigned int ConfigParser::NumSamples()
{
  return num_samples_;
}

CompositeObjectConfig ConfigParser::GripperConfig()
{
  return gripper_config_;
}

CompositeObjectConfig ConfigParser::ObjectConfig()
{
  return object_config_;
}


