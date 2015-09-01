/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
 *
 *    Description:  this file is the main file for my hybrid PRM path planner
 *
 *        Version:  1.0
 *        Created:  07/14/2011 06:54:39 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Zoe McCarthy (zm), ZoeMcCarthy12@gmail.com
 *        Company:  University of Illinois at Urbana-Champaign
 *
 * =====================================================================================
 */


#include    <stdlib.h>
#include    <string>
#include    <iostream>
#include    <fstream>
#include    <boost/thread/thread.hpp>
#include    "CageConfigurationParser.h"
#include    "CageChecker.h"
#include    "ShapeFactory.hpp"
#include    "Mesh.hpp"
#include    "Grasp.hpp"
#include    "GraspGenerator.hpp"

#include "glui/glui.h"

#include <glog/logging.h>
#include <opencv2/opencv.hpp>

// FCL includes
#include "fcl/traversal/traversal_node_bvhs.h"
#include "fcl/traversal/traversal_node_setup.h"
#include "fcl/collision_node.h"
#include "fcl/collision.h"
#include "fcl/BV/BV.h"
#include "fcl/shape/geometric_shapes.h"
#include "fcl/narrowphase/narrowphase.h"

// PolyDepth includes
#include <PQP.h>
#include <C2A/LinearMath.h>
#include <C2A/C2A.h>
#include <PolyDepth/PolyDepth.h>
#include <PolyDepth/MeshObject.h>

#include    <math.h>
#include    <time.h>

//#define VIS
using namespace std;

//helper wait function
void wait ( int seconds )
{
  clock_t endwait;
  endwait = clock () + seconds * CLOCKS_PER_SEC ;
  while (clock() < endwait) {}
}

// prints out usage statement
void print_help()
{
  std::cout << "Usage for alpha shapes caging program" << std::endl;
  std::cout << "alpha <config_filename>" << std::endl;
}

// save image to the mat
cv::Mat save_image(CageChecker& cage_checker, std::string image_filename, unsigned int width = 1500, unsigned int height = 1500)
{
  cv::Mat image(width, height, CV_8UC3);
  image.setTo(cv::Scalar(255, 255, 255));
  cage_checker.Render_Object_Tris(image);
  cage_checker.Render_Gripper_Tris(image);
  cv::imwrite(image_filename, image);
  return image;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  This function is the main function for my PRM planner.  It will be the
 *                command line interface to the geomview stream.
 * =====================================================================================
 */
int main ( int argc, char *argv[] )
{
  char** fake_argv;
  int fake_argc = 0;
  glutInit(&argc, argv);

  // get args
  if (argc < 2) {
    print_help();
  }

  std::string config_filename(argv[1]);
  //  srand(time(NULL));
  srand(110);

  float max_scale = -1.0f;
  if (argc > 2)
    max_scale = std::atof(argv[2]);
    
  CageConfigParser config(config_filename.c_str());
  if (config.Debug()) {
    config.Print();
  }
  bool random_grasps = config.RandomGrasps();
  bool close_grippers = config.CloseGrippers();
  float bounding = config.Bounding();
  float extrusion = config.Extrusion();
  GripperType grip_type = config.Gripper();
  std::string gripper_path = config.GripperPath();
  std::vector<std::string> object_list = config.ObjectList();

  MultiObjectConfig gc = config.GripperConfig();
  MultiObjectConfig oc = config.ObjectConfig();

  // render the sdfs for visualization
#ifdef VIS
  CGAL::Geomview_stream gv(CGAL::Bbox_3(-extrusion, -extrusion, -extrusion, extrusion, extrusion, extrusion));
  gv.set_line_width(4);
  gv.set_bg_color(CGAL::Color(0, 200, 200));
  gv.set_wired(true);
#endif

  // set up gripper config
  std::stringstream grip_sdf_filename;
  std::stringstream grip_obj_filename;
  grip_sdf_filename << gripper_path << ".csv";
  grip_obj_filename << gripper_path << ".obj";

  float grip_diff = config.GripperRadius();
  float sampling_scale = config.SamplingScale();
  float max_grip_width = 2*bounding;
  float sample_grip_width = 2*bounding;
  float delta_jaw = 1.0f;
  float grip_inc = 2.5f;
  GripperConfig gripper_config;
  gripper_config.scale = 1.0f;
  gripper_config.cx = 0;
  gripper_config.cy = 0;
  gripper_config.theta = 0;
  gripper_config.width = max_grip_width;
  gripper_config.sdf_filename = grip_sdf_filename.str();
  gripper_config.obj_filename = grip_obj_filename.str();
  std::vector<GripperConfig> gripper_configs;

  ComponentConfig object_config;
  object_config.scale = 1.0f;
  object_config.cx = 0;
  object_config.cy = 0;
  object_config.theta = 0;

  // write file header
  std::string results_filename = "caging_output_compound.csv";
  std::string tmp_results_filename = "caging_output_compound.csv";
  std::ofstream of_header(tmp_results_filename.c_str());  
  of_header << "# CSV of caging experiments run with " << config_filename << "\n";
  of_header << "# ObjectNum Samples Alpha CollRatio SampleTime AlphaTime IterTime TotalTime Volume PPExists Tx Ty Theta\n";
  of_header.close();

  // loop through objects and number of samples, recording the values for each
  int snapshot_interval = 1;
  int num_powers = 1;
  int num_trials = 3;
  int num_grasps = 10;
  int num_objects = object_list.size();
  int num_sim_trials = 10;
  int max_num_grasp_attempts = 10;

  float max_samples = config.NumSamples();
  float theta_inc = 4 * M_PI;// / 4.0f;
  float p_concave = 0.75f;
  float timeout = 120.0f;
  unsigned int grasp_id;
  Timer timer;
  if (!random_grasps)
    theta_inc = 2 * M_PI;

  float grasp_x, grasp_y, grasp_a;
  float grip_width;
      
  std::vector<float> samples;
  std::vector<CageResult> cage_results;
  std::vector<float> txs;
  std::vector<float> tys;
  std::vector<float> thetas;
  std::vector<float> d_thetas;// the rotation of the world that we use nameen
  std::vector<float> times;
  std::vector<PathPlanningResult> planned_path_results;
  std::vector<SimulationResult> simulation_results;
  std::vector<int> objects;
  std::vector<unsigned int> grasp_ids;
  std::vector<float> all_grip_widths;

  CfApproxConfig cf_config;
  cf_config.x_scale = 1.0f;
  cf_config.y_scale = 1.0f;
  cf_config.theta_scale = 1.0f;
  cf_config.x_min = -config.Bounding();
  cf_config.x_max =  config.Bounding();
  cf_config.y_min = -config.Bounding();
  cf_config.y_max =  config.Bounding();
  cf_config.theta_min = 0;
  cf_config.theta_max = 2 * M_PI;
  cf_config.num_rots = 1;

  EscapeEnergyConfig energy_config;
  energy_config.energy_min = -FLT_MAX;
  energy_config.energy_max = 0.0f;
  energy_config.energy_res = 0.001f;
  energy_config.prove_exists = true;

  Pose2D object_pose;
  object_pose.x = object_config.cx;
  object_pose.y = object_config.cy;
  object_pose.theta = object_config.theta;

  cv::Mat im1;  
  cv::Mat im2;  

  for (int k = 0; k < num_objects; k++) {
    LOG(INFO) << "Testing object " << object_list[k];
    gripper_configs.clear();

    // set up object config
    std::stringstream object_sdf_filename;
    std::stringstream object_obj_filename;
    std::stringstream object_grasps_filename;
    object_sdf_filename << object_list[k] << ".csv";
    object_obj_filename << object_list[k] << ".obj";
    object_grasps_filename << object_list[k] << "_grasps.txt";
    object_config.sdf_filename = object_sdf_filename.str();
    object_config.obj_filename = object_obj_filename.str();

    // create meshes
    Mesh* object = new Mesh(object_config);
    if (!object->Initialized()) {
      LOG(ERROR) << "Failed to initialize object";
      delete object;
      continue;
      // TODO: handle
    }
    cf_config.theta_scale = object->MaxMomentArm();

    // load grasp configurations
    int grasp_index = 0;
    LoadGraspsAndFilenames(object_grasps_filename.str(), gripper_config, gripper_configs);
    num_grasps = gripper_configs.size();

    // get the number of samples
    float cur_samples = max_samples / std::pow(2, num_powers - 1);

    // get the result for multiple numbers of samples
    while (cur_samples <= max_samples) {
      LOG(INFO) << "ATTEMPTING WITH " << cur_samples << " SAMPLES";

      // loop through angles
      for (float d_theta = 0.0f; d_theta <= 0.0f; d_theta += theta_inc) { 
        d_theta = gripper_configs[0].angle;
        LOG(INFO) << "USING ANGLE = " << d_theta;

        // iteratively try the same grasp to average out random effects
        for (int j = 0; j < num_trials; j++) {
          LOG(INFO) << "GRASP TRIAL " << j << " of " << num_trials;

          // create gripper
          std::vector<Mesh*> gripper_fingers;
          for (int g = 0; g < num_grasps; g++) {
            OneFingerGrasp2D* gripper = new OneFingerGrasp2D(gripper_configs[g]);
            if (gripper->Initialized()) {
              gripper_fingers.push_back(gripper);
            }
            else {
              LOG(ERROR) << "Failed to initialize gripper";
              delete gripper;
              continue;
            }
          }
          
          // initialize new cage checker with optional vis
#ifdef VIS
          CageChecker cage_checker(&gv, cf_config, object, gripper_fingers);  
#else
          CageChecker cage_checker(cf_config, object, gripper_fingers);  
#endif
          // reset poses
          cage_checker.Set_Object_Pose(object_config.cx, object_config.cy, object_config.theta);
          
          // rotate world by theta for different gravity vectors
          cage_checker.Rotate_World(d_theta);

          // draw the grippers
#ifdef VIS            
          gv.clear();
          cage_checker.Draw_Gripper();
          cage_checker.Draw_Object();
#endif

          // draw and save an image
          if (j >= 0){ 
            grasp_id = rand();
            std::stringstream image_filename;
            image_filename << "../results/icra_compound_examples/object_" << k << "_cage_width_" << grip_width << "_" << grasp_id << "_" << d_theta << ".jpg"; 
            im1 = save_image(cage_checker, image_filename.str());
            LOG(INFO) << "Saved image to " << image_filename.str();            
          }

          // check for initial collision
          MeshCollisionResult r = Mesh::LowerBoundCollision(object, gripper_fingers[0], false);
          if (r.collision) {
            LOG(ERROR) << "Invalid Configuration. Initial object pose is in collision.";
            LOG(INFO) << "Dist " << r.distance;
            continue;
          }
          else {
            LOG(INFO) << "No collision";
          }
                
          // init results
          CageResult cage_result;
          SimulationResult sr;
          PathPlanningResult ppr;              

          LOG(INFO) << "Object mass: " << object->Mass();
          LOG(INFO) << "Object moment arm: " << object->MaxMomentArm();

          // check for a cage
          LOG(INFO) << "Checking cage with " << cur_samples << " samples on trial " << j;
          timer.reset();
          timer.start();
          energy_config.num_samples = cur_samples;
          cage_result = cage_checker.Check_Cage(energy_config);
          timer.stop();
          LOG(INFO) << "Object caged?: " << !(cage_result.energy_result.path_exists);
          LOG(INFO) << "Cage energy: " << cage_result.energy_result.energy;
          LOG(INFO) << "Normalized energy: " << cage_result.energy_result.normalized_energy;

          // simulate to get empirical cage ratio (have to reset cage checker - silly but need to just get it running)
          Eigen::Matrix4f rot_pose = CreatePose(0, 0, d_theta);
          Eigen::Matrix4f obj_pose2 = CreatePose(object_config.cx, object_config.cy, object_config.theta);
          object->SetPose(rot_pose * obj_pose2);
          timer.reset();
          timer.start();
          sr = cage_checker.Ratio_Escape_Paths(object_config.cx, object_config.cy, object_config.theta,
                                               num_sim_trials, -cage_result.energy_result.energy);
          timer.stop();
          LOG(INFO) << "Thresholded cage ratio: " << sr.cage_rate;
          LOG(INFO) << "Raw cage ratio: " << sr.raw_cage_rate;
          LOG(INFO) << "Simulation took (sec): " << timer.time();

          // check for an escape path with rrt
          object->SetPose(rot_pose * obj_pose2);
          Eigen::Matrix4f obj_pose = object->Pose();
          float gx, gy, gtheta;
          Pose3DToParams2D(obj_pose, gx, gy, gtheta);
                
          LOG(INFO) << "Planning escape path";
          timer.reset();
          timer.start();
          float range = 0.025f;
          ppr = cage_checker.Find_Escape_Path(gx, gy, gtheta,
                                              std::max<float>(0.0f, -cage_result.energy_result.energy),
                                              timeout, range);
          timer.stop();
          LOG(INFO) << "Object able to escape?: " << ppr.path_exists;
          LOG(INFO) << "Escape Energy: " << ppr.max_energy << " at state " << ppr.max_state;
          LOG(INFO) << "Planning took (sec): " << timer.time();

          // visualization
          if (ppr.path_exists) {
            for (unsigned int a = 0; a < ppr.states.size(); a++) {
              cage_checker.Set_Object_Pose(ppr.states[a](0), ppr.states[a](1), ppr.states[a](2));
#ifdef VIS
              if (a % 5 == 0)
                LOG(INFO) << "State " << a << ": " << ppr.states[a](0) << " " << ppr.states[a](1) << " " << ppr.states[a](2);

              gv.clear();
              cage_checker.Draw_Gripper();
              cage_checker.Draw_Object();
              boost::this_thread::sleep(boost::posix_time::seconds(0.25f));
#endif

              if (a == ppr.max_state || (abs((int)a - (int)ppr.max_state) < 10 && a % 2 == 0)) {
                std::stringstream image_filename;
                image_filename << "../results/icra_compound_examples/object_" << k << "_cage_width_" << grip_width << "_" << grasp_id << "_path_state_" << a << ".jpg"; 
                im2 = save_image(cage_checker, image_filename.str());

                if (a == ppr.max_state) {
                  float alpha = 0.7f;
                  cv::addWeighted(im1, alpha, im2, 1.0 - alpha, 0.0, im2);
                  std::stringstream image_filename2;
                  image_filename2 << "../results/icra_compound_examples/object_" << k << "_cage_width_" << grip_width << "_" << grasp_id << "_path_state_translucent_" << a << ".jpg"; 
                  cv::imwrite(image_filename2.str(), im2);
                }                      
                LOG(INFO) << "Saved image to " << image_filename.str();            
              }
            }
          }

          // save results of the run
          txs.push_back(gripper_config.cx);
          tys.push_back(gripper_config.cy);
          thetas.push_back(gripper_config.theta);
          d_thetas.push_back(d_theta);
          all_grip_widths.push_back(grip_width);
          samples.push_back(cur_samples);
          cage_results.push_back(cage_result);
          planned_path_results.push_back(ppr);
          simulation_results.push_back(sr);
          times.push_back(timer.time());
          objects.push_back(k);
          grasp_ids.push_back(grasp_id);
          grasp_index++;

          if (grasp_index % snapshot_interval == 0) {
            // write the result to file
            std::ofstream of(tmp_results_filename.c_str(), std::ios::app);  
            for (unsigned int j = txs.size() - snapshot_interval; j < samples.size(); j++) {
              of << objects[j] << ", " << grasp_ids[j] << ", " << samples[j] << ", " << cage_results[j].energy_result.energy << ", " << cage_results[j].energy_result.normalized_energy << \
                ", " << cage_results[j].total_time << ", " << cage_results[j].energy_result.sample_time << ", " << cage_results[j].energy_result.triangulation_time << \
                ", " << cage_results[j].energy_result.iter_time << ", " << planned_path_results[j].path_exists << ", " << planned_path_results[j].max_energy << ", " << simulation_results[j].cage_rate << ", " << simulation_results[j].raw_cage_rate << \
                ", " << txs[j] << ", " << tys[j] << ", " << thetas[j] << ", " << d_thetas[j] << ", " << all_grip_widths[j] << "\n";
            }
            of.close();
          }

          // clean up grippers
          for (int g = 0; g < num_grasps; g++) {
            delete gripper_fingers[g];
          }
        }
      }
      cur_samples = cur_samples * 2;
    }
    // clean up
    delete object;
  }


  // write all results to file
  std::ofstream of(results_filename.c_str());  
  of << "# CSV of caging experiments run with " << config_filename << "\n";
  of << "# ObjectNum Samples Alpha CollRatio SampleTime AlphaTime IterTime TotalTime Volume Tx Ty Theta\n";
  for (unsigned int j = 0; j < samples.size(); j++) {
    of << objects[j] << ", " << grasp_ids[j] << ", " << samples[j] << ", " << cage_results[j].energy_result.energy << ", " << cage_results[j].energy_result.normalized_energy << \
      ", " << cage_results[j].total_time << ", " << cage_results[j].energy_result.sample_time << ", " << cage_results[j].energy_result.triangulation_time << \
      ", " << cage_results[j].energy_result.iter_time << ", " << planned_path_results[j].path_exists << ", " << planned_path_results[j].max_energy << ", " << simulation_results[j].cage_rate << ", " << simulation_results[j].raw_cage_rate << \
      ", " << txs[j] << ", " << tys[j] << ", " << thetas[j] << ", " << d_thetas[j] << ", " << all_grip_widths[j] << "\n";
  }
  of.close();

  return EXIT_SUCCESS;
}

