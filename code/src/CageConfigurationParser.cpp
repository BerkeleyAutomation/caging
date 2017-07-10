#include "CageConfigurationParser.h"

CageConfigParser::CageConfigParser(std::string filename) : config_filename_(filename)
{
  Parse();
}                                                         

// Parses the YAML config file
bool CageConfigParser::Parse()
{
  YAML::Node node = YAML::LoadFile(config_filename_.c_str());

  // attempt to parse config
  // NOTE: assumes given keys are there
  try {
    // get constants
    int gt;
    std::string list_file;
    extrusion_ = node["extrusion"].as<float>();
    bounding_ = node["bounding"].as<float>();
    debug_ = node["debug"].as<bool>();
    random_grasps_ = node["random_grasps"].as<bool>();
    close_grippers_ = node["close_grippers"].as<bool>();
    num_samples_ = node["samples"].as<unsigned int>();
    sampling_scale_ = node["sampling_scale"].as<float>();
    gt = node["gripper_type"].as<int>();
    gripper_radius_ = node["gripper_radius"].as<float>();
    gripper_path_ = node["gripper_path"].as<std::string>();
    list_file = node["object_list"].as<std::string>();

    if (gt == 0) {
      gripper_type_ = PARALLEL_JAW;
    }
    else if (gt == 1) {
      gripper_type_ = RIGID;
    } 
    else if (gt == 2) {
      gripper_type_ = QUAD;
    } 
    else {
      gripper_type_ = TRI;
    }
    ParseObjectFilenames(list_file);
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
  std::cout << std::endl;
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
