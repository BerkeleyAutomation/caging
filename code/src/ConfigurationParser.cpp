#include "ConfigurationParser.h"

ConfigParser::ConfigParser(std::string filename) : config_filename_(filename)
{
  Parse();
}                                                         

// Parses the YAML config file
bool ConfigParser::Parse()
{
  std::ifstream yaml_in(config_filename_.c_str());
  YAML::Parser parser(yaml_in);
  YAML::Node node;

  // attempt to parse config
  // NOTE: assumes given keys are there
  try {
    parser.GetNextDocument(node);

    // get constants
    node["extrusion"] >> extrusion_;
    node["bounding"] >> bounding_;
    node["debug"] >> debug_;
    node["samples"] >> num_samples_;

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
    node["boxes"][i]["width"] >> box.width;
    node["boxes"][i]["height"] >> box.height;
    node["boxes"][i]["cx"] >> box.cx;
    node["boxes"][i]["cy"] >> box.cy;
    node["boxes"][i]["theta"] >> box.theta;
    comp_config.boxes.push_back(box);
  }

  // parse triangles
  std::vector<TriangleConfig> tri_configs;
  for (unsigned int i = 0; i < node["triangles"].size(); i++) {
    TriangleConfig tri;
    node["triangles"][i]["scale"] >> tri.scale;
    node["triangles"][i]["cx"] >> tri.cx;
    node["triangles"][i]["cy"] >> tri.cy;
    node["triangles"][i]["theta"] >> tri.theta;
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


