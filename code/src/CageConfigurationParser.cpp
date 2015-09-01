#include "CageConfigurationParser.h"

CageConfigParser::CageConfigParser(std::string filename) : config_filename_(filename)
{
  Parse();
}                                                         

// Parses the YAML config file
bool CageConfigParser::Parse()
{
  std::ifstream yaml_in(config_filename_.c_str());
  YAML::Parser parser(yaml_in);
  YAML::Node node;

  // attempt to parse config
  // NOTE: assumes given keys are there
  try {
    parser.GetNextDocument(node);

    // get constants
    int gt;
    std::string list_file;
    node["extrusion"] >> extrusion_;
    node["bounding"] >> bounding_;
    node["debug"] >> debug_;
    node["random_grasps"] >> random_grasps_;
    node["close_grippers"] >> close_grippers_;
    node["samples"] >> num_samples_;
    node["sampling_scale"] >> sampling_scale_;
    node["gripper_type"] >> gt;
    node["gripper_radius"] >> gripper_radius_;
    node["gripper_path"] >> gripper_path_;
    node["object_list"] >> list_file;

    if (gt == 0) {
      gripper_type_ = PARALLEL_JAW;
    }
    else {
      gripper_type_ = BARRETT_HAND;
    }
    ParseObjectFilenames(list_file);

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

void CageConfigParser::Print()
{
  std::cout << "CONFIGURATION VALUES" << std::endl;
  std::cout << "Extrusion: " << extrusion_ << std::endl;
  std::cout << "Bounding: " << bounding_ << std::endl;
  std::cout << "Debug: " << debug_ << std::endl;
  std::cout << "Num Samples: " << num_samples_ << std::endl;
  std::cout << "Gripper Path: " << gripper_path_ << std::endl;
  std::cout << "Num Gripper Components: " << gripper_config_.size() << std::endl;
  std::cout << "Num Object Components: " << object_config_.size() << std::endl;
  std::cout << std::endl;
}

void CageConfigParser::ParseCompositeConfig(const YAML::Node& node, MultiObjectConfig& comp_config)
{
  comp_config.clear();

  // parse boxes
  for (unsigned int i = 0; i < node.size(); i++) {
    ComponentConfig c;
    node[i]["scale"] >> c.scale;
    node[i]["cx"] >> c.cx;
    node[i]["cy"] >> c.cy;
    node[i]["theta"] >> c.theta;
    node[i]["sdf_filename"] >> c.sdf_filename;
    node[i]["obj_filename"] >> c.obj_filename;
    comp_config.push_back(c);
  }

}

void CageConfigParser::ParseObjectFilenames(std::string list_file)
{
  std::ifstream ifile(list_file.c_str());
  std::string filename;
  object_list_.clear();
  while (!ifile.eof()) {
    ifile >> filename;
    object_list_.push_back(filename);
  }
  object_list_.pop_back();
}

float CageConfigParser::Extrusion()
{
  return extrusion_;
}

float CageConfigParser::Bounding()
{
  return bounding_;
}

bool CageConfigParser::Debug()
{
  return debug_;
}

bool CageConfigParser::RandomGrasps()
{
  return random_grasps_;
}

bool CageConfigParser::CloseGrippers()
{
  return close_grippers_;
}

unsigned int CageConfigParser::NumSamples()
{
  return num_samples_;
}

float CageConfigParser::SamplingScale()
{
  return sampling_scale_;
}

GripperType CageConfigParser::Gripper()
{
  return gripper_type_;
}

float CageConfigParser::GripperRadius()
{
  return gripper_radius_;
}

std::string CageConfigParser::GripperPath()
{
  return gripper_path_;
}

std::vector<std::string> CageConfigParser::ObjectList()
{
  return object_list_;
}

MultiObjectConfig CageConfigParser::GripperConfig()
{
  return gripper_config_;
}

MultiObjectConfig CageConfigParser::ObjectConfig()
{
  return object_config_;
}


