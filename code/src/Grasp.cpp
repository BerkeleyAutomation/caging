#include "Grasp.hpp"

ParallelJawGrasp2D::ParallelJawGrasp2D(GripperConfig config)
  : config_(config)
{
  SetState(config);
}

void ParallelJawGrasp2D::SetState(GripperConfig config)
{
  const unsigned int num_gripper_fingers = 2;
  std::vector<ComponentConfig> configs;
  Eigen::Matrix3f root_pose = CreatePose2D(config.cx, config.cy, config.theta);
  Eigen::Matrix4f root_pose_3d = CreatePose(config.cx, config.cy, config.theta);
  Eigen::Matrix3f grip_pose;
  float tx, ty, theta;
  configs.resize(num_gripper_fingers);

  // gripper 1
  configs[0].scale = config.scale;
  configs[0].cx = 0;
  configs[0].cy = -config.width / 2;
  configs[0].theta = 0;//M_PI/ 2;
  configs[0].sdf_filename = config.sdf_filename;
  configs[0].obj_filename = config.obj_filename;

  // gripper 2
  configs[1].scale = config.scale;
  configs[1].cx = 0;
  configs[1].cy = config.width / 2;
  configs[1].theta = M_PI;
  configs[1].sdf_filename = config.sdf_filename;
  configs[1].obj_filename = config.obj_filename;

  // init root object
  if (!Initialized()) {
    Mesh::Initialize(configs, root_pose_3d);
  }
  else {
    // std::cout << "Root pose " << config.theta << " " << config.width << std::endl;
    // std::cout << root_pose_3d << std::endl;
    Mesh::SetPose(root_pose_3d);
    Mesh::SetRelativePoses(configs);
  }
}

void ParallelJawGrasp2D::EndpointsToGrasp(ASE::Vertex2d g1, ASE::Vertex2d g2, GripperConfig& config)
{
  // directions
  ASE::Direction2d dir = g2 - g1;
  config.width = dir.norm();
  if (config.width == 0) {
    config.width = 1;
  }
  dir = dir / config.width;

  // compute the grasp center
  ASE::Vertex2d center = (g1 + g2) / 2.0f;
  config.cx = center(1);
  config.cy = center(0);
  //  std::cout << "center " << config.cx << " " << config.cy << std::endl;
  config.theta = acos(dir(0));
  if (dir(1) < 0) {
    config.theta = 2 * M_PI - config.theta;  
  }
}

ThreeFingerGrasp2D::ThreeFingerGrasp2D(GripperConfig config)
  : config_(config)
{
  SetState(config);
}

void ThreeFingerGrasp2D::SetState(GripperConfig config)
{
  const unsigned int num_gripper_fingers = 3;
  std::vector<ComponentConfig> configs;
  Eigen::Matrix3f root_pose = CreatePose2D(config.cx, config.cy, config.theta);
  Eigen::Matrix4f root_pose_3d = CreatePose(config.cx, config.cy, config.theta);
  Eigen::Matrix3f grip_pose;
  configs.resize(num_gripper_fingers);

  float arc = 1.2f;
  Eigen::Vector3f v;
  v << 0, -config.width, 1.0f;
  float theta0 = 0.0f * 2.0f * M_PI / 3.0f;
  float theta1 = M_PI - arc * M_PI / 3.0f;
  float theta2 = M_PI + arc * M_PI / 3.0f;
  Eigen::Matrix3f T0 = CreatePose2D(0.0f, 0.0f, -theta0);
  Eigen::Matrix3f T1 = CreatePose2D(0.0f, 0.0f, -theta1);
  Eigen::Matrix3f T2 = CreatePose2D(0.0f, 0.0f, -theta2);
  Eigen::Vector3f w;

  // gripper 1
  w = T0 * v;
  configs[0].scale = config.scale;
  configs[0].cx = w(0);
  configs[0].cy = w(1);
  configs[0].theta = theta0;//M_PI/ 2;
  configs[0].sdf_filename = config.sdf_filename;
  configs[0].obj_filename = config.obj_filename;

  // gripper 2
  w = T1 * v;
  configs[1].scale = config.scale;
  configs[1].cx = w(0);
  configs[1].cy = w(1);
  configs[1].theta = theta1;
  configs[1].sdf_filename = config.sdf_filename;
  configs[1].obj_filename = config.obj_filename;

  // gripper 2
  w = T2 * v;
  configs[2].scale = config.scale;
  configs[2].cx = w(0);
  configs[2].cy = w(1);
  configs[2].theta = theta2;
  configs[2].sdf_filename = config.sdf_filename;
  configs[2].obj_filename = config.obj_filename;

  // init root object
  if (!Initialized()) {
    Mesh::Initialize(configs, root_pose_3d);
  }
  else {
    // std::cout << "Root pose " << config.theta << " " << config.width << std::endl;
    // std::cout << root_pose_3d << std::endl;
    Mesh::SetPose(root_pose_3d);
    Mesh::SetRelativePoses(configs);
  }
}

OneFingerGrasp2D::OneFingerGrasp2D(GripperConfig config)
  : config_(config)
{
  SetState(config);
}

void OneFingerGrasp2D::SetState(GripperConfig config)
{
  const unsigned int num_gripper_fingers = 1;
  std::vector<ComponentConfig> configs;
  Eigen::Matrix3f root_pose = CreatePose2D(config.cx, config.cy, config.theta);
  Eigen::Matrix4f root_pose_3d = CreatePose(config.cx, config.cy, config.theta);
  Eigen::Matrix3f grip_pose;
  configs.resize(num_gripper_fingers);

  // gripper 1
  configs[0].scale = config.scale;
  configs[0].cx = 0.0f;
  configs[0].cy = 0.0f;
  configs[0].theta = 0.0f;
  configs[0].sdf_filename = config.sdf_filename;
  configs[0].obj_filename = config.obj_filename;

  // init root object
  if (!Initialized()) {
    Mesh::Initialize(configs, root_pose_3d);
  }
  else {
    Mesh::SetPose(root_pose_3d);
    Mesh::SetRelativePoses(configs);
  }
}

void LoadGrasps(const std::string& filename, GripperConfig config_template, std::vector<GripperConfig>& configs)
{
  std::ifstream ifile(filename.c_str());

  char line_buffer[2000];
  while(!ifile.eof()) { 
    ifile.getline(line_buffer, 2000);
    char first_token = line_buffer[0];
    if (first_token == '#' || first_token == 0)
      continue;

    // add to config list
    GripperConfig config;
    float width;
    std::stringstream line_parser(line_buffer);
    line_parser >> config.cx;
    line_parser >> config.cy >> config.theta >> config.width >> config.angle;
    config.scale = config_template.scale;
    config.sdf_filename = config_template.sdf_filename;
    config.obj_filename = config_template.obj_filename;
    config.theta = config.theta;// + M_PI / 2;
    configs.push_back(config);
  }  
}

void LoadGraspsAndFilenames(const std::string& filename, GripperConfig config_template, std::vector<GripperConfig>& configs)
{
  std::ifstream ifile(filename.c_str());

  char line_buffer[2000];
  while(!ifile.eof()) { 
    ifile.getline(line_buffer, 2000);
    char first_token = line_buffer[0];
    if (first_token == '#' || first_token == 0)
      continue;

    // add to config list
    GripperConfig config;
    float width;
    std::stringstream line_parser(line_buffer);
    line_parser >> config.cx;
    line_parser >> config.cy >> config.theta >> config.width >> config.angle;
    config.scale = config_template.scale;
    
    std::string gripper_path;
    line_parser >> gripper_path;
    std::stringstream grip_sdf_filename;
    std::stringstream grip_obj_filename;
    grip_sdf_filename << gripper_path << ".csv";
    grip_obj_filename << gripper_path << ".obj";

    config.sdf_filename = grip_sdf_filename.str();
    config.obj_filename = grip_obj_filename.str();
    config.theta = config.theta;// + M_PI / 2;
    configs.push_back(config);
  }  
}
