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
#include    "Util.h"

#include "glui/glui.h"

#include <glog/logging.h>
#include <opencv2/opencv.hpp>

// FCL includes
#include <fcl/traversal/traversal_node_bvhs.h>
#include <fcl/traversal/traversal_node_setup.h>
#include <fcl/collision_node.h>
#include <fcl/collision.h>
#include <fcl/BV/BV.h>
#include <fcl/shape/geometric_shapes.h>
#include <fcl/narrowphase/narrowphase.h>

// PolyDepth includes
#include <PQP.h>
#include <C2A/LinearMath.h>
#include <C2A/C2A.h>
#include <PolyDepth/PolyDepth.h>
#include <PolyDepth/MeshObject.h>

#include    <math.h>
#include    <time.h>

#define VIS
using namespace std;

// prints out usage statement
void print_help()
{
  std::cout << "Usage for alpha shapes caging program" << std::endl;
  std::cout << "alpha <config_filename>" << std::endl;
}

// save image to the mat
cv::Mat save_image(CageChecker& cage_checker, std::string image_filename, float theta_offset, unsigned int width = 1500, unsigned int height = 1500)
{
  cv::Mat image(width, height, CV_8UC3);
  image.setTo(cv::Scalar(255, 255, 255));
  cage_checker.Render_Object_Tris(image);
  cage_checker.Render_Gripper_Tris(image);
  //But we also want to draw the direction of force 
  Eigen::Vector2d rotation_vector(0, -((float) height)/8);
  //Use negative theta offset since our fram of refernce for gemview flips rotations
  rotation_vector = rotate_vector(rotation_vector, -theta_offset); 
  cv::Point2i center_point(height/8, width/8);
  //Note that for some odd reason, "up" is negative in OpenCV
  cv::Point2i end_point(height/8 + rotation_vector(0), width/8 - rotation_vector(1));
  cv::circle(image, center_point, 10, cv::Scalar(0, 255, 0), 4);
  cv::line(image, center_point, end_point, cv::Scalar(0, 255, 0), 4);
  cv::imwrite(image_filename, image);
  return image;
}

//Version of image saving for only the object (mainly for physical experiment templates).
cv::Mat save_image_object_only(CageChecker& cage_checker, std::string image_filename, unsigned int width = 1500, unsigned int height = 1500)
{
  cage_checker.Set_Object_Pose(0, 0, 0);
  cv::Mat image(width, height, CV_8UC3);
  image.setTo(cv::Scalar(255, 255, 255));
  cage_checker.Render_Object_Tris(image);
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

  bool debug = false;

  // get args
  if (argc < 2) {
    print_help();
    exit(0);
  }
  std::string config_filename(argv[1]);
  std::string result_dir(argv[2]);
  //The results directory must end in "/"
  if (result_dir.at(result_dir.length() - 1) != '/') {
    result_dir = result_dir + "/";
  }

  //srand(time(NULL));
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
  //float angle_sweep = M_PI / 4.0f;
  float angle_sweep = 0.0f / 4.0f;
  float angle_disc = M_PI / 8.0f;
  float max_push_force = 1.0f;
  int num_searches = 25;
  float rrt_range = 0.25f;
  float rrt_vis_range = 0.25f;
  float rrt_timeout = 120.0f;
  float rrt_short_timeout = 30.0f;
  float alpha = 0.7f;
  GripperConfig gripper_config;
  gripper_config.scale = 1.0f;
  gripper_config.cx = 0;
  gripper_config.cy = 0;
  gripper_config.theta = 0;
  gripper_config.angle = 0;
  gripper_config.width = max_grip_width;
  gripper_config.sdf_filename = grip_sdf_filename.str();
  gripper_config.obj_filename = grip_obj_filename.str();
  std::vector<GripperConfig> gripper_configs;
  gripper_configs.push_back(gripper_config);

  ComponentConfig object_config;
  object_config.scale = 1.0f;
  object_config.cx = 0;
  object_config.cy = 0;
  object_config.theta = 0;

  // loop through objects and number of samples, recording the values for each
  int snapshot_interval = 1;
  int num_powers = 1;
  int num_grasps = 1;
  int num_grasp_trials = 1;
  int num_objects = object_list.size();
  int num_sim_trials = 10;
  int max_num_grasp_attempts = 10;

  float max_samples = config.NumSamples();
  float theta_inc = 4 * M_PI;// / 4.0f;
  float p_concave = 0.75f;
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

  // create gripper
  Mesh* gripper;
  if (grip_type == PARALLEL_JAW) {
    gripper = new ParallelJawGrasp2D(gripper_config);
  }
  else if (grip_type == RIGID) {
    gripper = new OneFingerGrasp2D(gripper_config);
  } 
  else if (grip_type == QUAD){
    gripper = new QuadGrasp2D(gripper_config);
  }
  else {
    gripper = new TriGrasp2D(gripper_config);
  }
  if (!gripper->Initialized()) {
    LOG(INFO) << "Failed to initialize object";
    delete gripper;
    // TODO: handle
  }
  std::vector<Mesh*> gripper_fingers;
  gripper_fingers.push_back(gripper);

  cv::Mat im1;  
  cv::Mat im2;  

  for (int k = 0; k < num_objects; k++) {
    LOG(INFO) << "Using object " << object_list[k];

    // open file
    std::stringstream results;
    results << object_list[k]  << "_synthesis_output.csv";
    
    std::string results_filename = results.str();
    results_filename = get_filename(results_filename);
    results_filename = result_dir +  results_filename;
    std::ofstream of_header(results_filename.c_str());  
    LOG(INFO) << "Saving results to " << results_filename;
    of_header << "trial" << ", " << "gripper_width" << ", " << "pose_x"  << ", " << "pose_y" << ", " << "pose_theta" << ", "\ 
              << "energy"  << ", "  << "norm_energy" << ", " << "void_volume" <<  ", " "void_area" << ", " << "sample_time" << ", "\ 
              << "triangulation_time" <<  ", " <<  "filtration_time" << ", " << "persistence_time" <<  ", " \
              << "void_time" <<  ", " << "push_direction" << ", " << "image_filename" << ", " <<  "rrt* found path" \
              << ", " << "rrt* max energy" << ", " << "num_samples" << "\n";
    of_header.close();

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
    LOG(INFO) << "Initializing object";
    Mesh* object = new Mesh(object_config);
    if (!object->Initialized()) {
      LOG(INFO) << "Failed to initialize object";
      delete object;
      continue;
      // TODO: handle
    }
    LOG(INFO) << "Object initialized";
    cf_config.theta_scale = object->MaxMomentArm()*1.5f; // for numeric tolerancing
    LOG(INFO) << "Max moment: " << object->MaxMomentArm();
    LOG(INFO) << "Mass: " << object->Mass();    

    num_grasps = gripper_configs.size();

    // iterate for a number of individual grasps
    for (int g = 0; g < num_grasps; g++) {
        LOG(INFO) << "GRASP " << g+1 << " of " << num_grasps;

        // get the random grasps
        float cur_samples = max_samples / std::pow(2, num_powers - 1);

        // set up next grasp
        gripper_config = gripper_configs[g];
        if (grip_type == PARALLEL_JAW) {
          static_cast<ParallelJawGrasp2D*>(gripper_fingers[0])->SetState(gripper_config);
        }
        else {
          static_cast<OneFingerGrasp2D*>(gripper_fingers[0])->SetState(gripper_config);
        }

        // set grasp width based on object shape
        gripper_config.width = grip_diff; //object->MaxMomentArm() / 2.0f;
        float orig_grip_width = gripper_config.width;
        LOG(INFO) << "CX " << gripper_config.cx << " CY " << gripper_config.cy << " THETA " << gripper_config.theta << " WIDTH " << gripper_config.width;

        // get the result for multiple numbers of samples
        while (cur_samples <= max_samples) {
          LOG(INFO) << "ATTEMPTING WITH " << cur_samples << " SAMPLES";

          // loop through angles
          for (int trial = 0; trial < num_grasp_trials; trial++) {
            float d_theta = gripper_config.angle;
            LOG(INFO) << "TRIAL " << trial;
            LOG(INFO) << "USING ANGLE = " << d_theta;
		    
            grip_inc = orig_grip_width;
	    grip_width = orig_grip_width;
	    if (grip_type == PARALLEL_JAW) {
	      //max_grip_width = 2 * object->MaxMomentArm();
	      max_grip_width = grip_width + grip_inc; 
	    } else {
	      //If we're using the Zeke gripper, width iteration makes no sense!
	      max_grip_width = grip_width + grip_inc; 
	    }
            int grasp_num = 0;
	    while (grip_width < max_grip_width) {
              LOG(INFO) << "GRASP WIDTH " << grip_width;

              //Update width of the grippers that we use...
	      gripper_config.width = grip_width;
              if (grip_type == PARALLEL_JAW) {
                static_cast<ParallelJawGrasp2D*>(gripper_fingers[0])->SetState(gripper_config);
              }
              else {
                static_cast<OneFingerGrasp2D*>(gripper_fingers[0])->SetState(gripper_config);
              }
              
              // initialize new cage checker with optional vis
#ifdef VIS
              CageChecker cage_checker(&gv, cf_config, object, gripper_fingers);  
#else
              CageChecker cage_checker(cf_config, object, gripper_fingers);  
#endif
              // reset poses
              if (!debug) {
                cage_checker.Set_Object_Pose(object_config.cx, object_config.cy, object_config.theta);
              }
              cage_checker.Set_Gripper_Pose(gripper_config.cx, gripper_config.cy, gripper_config.theta);        

              // rotate world by theta for different gravity vectors
              cage_checker.Rotate_World(d_theta);
              
              // draw the grippers
#ifdef VIS
              gv.clear();
              cage_checker.Draw_Gripper();
              cage_checker.Draw_Object();
#endif        
	      ////OUR CORE ALGORITHM LIVES HERE 
	      timer.reset();
              timer.start();
              energy_config.num_samples = cur_samples;
              //Checking this doesn't make sense with the quad gripper because of complete cages.
	      bool check_reachability = grip_type != QUAD; 
	      std::vector< std::vector<synthesis_info> > all_poses =  \
                cage_checker.synthesize_grasps(energy_config, num_searches, 
					       angle_sweep, angle_disc, check_reachability, max_push_force);
              timer.stop();
              
	      //LET's SAVE SOME RESULTS	
              // write file header
 
	      std::ofstream of(results_filename.c_str(), std::ios::app);  	
	      
	      for (int r = 0; r < all_poses.size(); r++) { 
	        std::vector<synthesis_info> extracted_poses = all_poses[r];

	        for (unsigned int j = 0 ; j < extracted_poses.size(); j++) { 
		  float seen_rot = roundf(extracted_poses[j].get_rotation() * 100)/100;
	          std::ostringstream ss;
	          ss << object_list[k] << "_rot_" << seen_rot <<  "_grasp_" << grasp_num << "_" << j\
                     << "_width_" << grip_width << "_.png";
	          std::string pic_string = ss.str();
	          pic_string = get_filename(pic_string);
	          pic_string = result_dir + pic_string; 
	          Pose2D cur_pose;
	          cur_pose = extracted_poses[j].get_pose();
                  float delta_lower_bound = extracted_poses[j].get_delta();
                  float energy_dir = extracted_poses[j].get_rotation();

                  cage_checker.Set_Object_Pose(cur_pose.x, cur_pose.y, cur_pose.theta);
                  gv.clear();
                  cage_checker.Draw_Gripper();
                  cage_checker.Draw_Object();
	          
                  //DRAW THE BIRTH POSE 
	          ss << object_list[k] << "_BIRTH" << "_rot_" << seen_rot <<  "_grasp_" << grasp_num << "_" << j\
                     << "_width_" << grip_width << "_.png";
	          std::string birth_string = ss.str();
	          birth_string = get_filename(birth_string);
	          birth_string = result_dir + birth_string; 
	          Pose2D birth_pose;
	          birth_pose = extracted_poses[j].get_birth_pose();
		  cage_checker.Set_Object_Pose(birth_pose.x, birth_pose.y, birth_pose.theta);
		  save_image(cage_checker, birth_string, extracted_poses[j].get_rotation());
		  
		  //Check in case the algorithm gives ys a bad bound or bad cage...
		  cage_checker.Set_Object_Pose(cur_pose.x, cur_pose.y, cur_pose.theta);
                  PathPlanningResult ppr;  //= cage_checker.Upper_Bound_Escape_Energy(cur_pose.x, cur_pose.y, cur_pose.theta,
                  //                                                                 std::max<float>(0.0f,
                  //                                                                                delta_lower_bound),
                  //                                                                energy_dir,
                  //                                                                rrt_timeout, rrt_range, max_push_force);
                  //
                  //LOG(INFO) << "Object able to escape?: " << ppr.path_exists;
                  //LOG(INFO) << "Escape Energy: " << ppr.max_energy << " at state " << ppr.max_state;
                  //LOG(INFO) << "Normalized Escape Energy: " << ppr.normalized_max_energy << " at state " << ppr.max_state;

                  //// visualization
                  //if (ppr.path_exists) {
                  //  for (unsigned int a = 0; a < ppr.states.size(); a++) {
                  //    cage_checker.Set_Object_Pose(ppr.states[a](0), ppr.states[a](1), ppr.states[a](2));
                  //    MeshCollisionResult r = Mesh::LowerBoundCollision(object, gripper, false);
                  //    if (r.collision) {
                  //      LOG(INFO) << "Collision at state " << a;
                  //    }

                  //    if (a % 1 == 0)
                  //      LOG(INFO) << "State " << a << ": " << ppr.states[a](0) << " " << ppr.states[a](1) << " " << ppr.states[a](2);

                  //    gv.clear();
                  //    cage_checker.Draw_Gripper();
                  //    cage_checker.Draw_Object();
                  //    boost::this_thread::sleep(boost::posix_time::seconds(1.0f));

                  //    std::stringstream image_filename;
                  //    image_filename << "state_" << a << ".png";
		  //      
                  //    std::string pic_string2 = image_filename.str();
                  //    pic_string2 = get_filename(pic_string2);
                  //    pic_string2 = result_dir + pic_string2; 
                  //    im2 = save_image(cage_checker, pic_string2, energy_dir);                  
                  //  }
                  //}

		  //Write to reults .csv
		  of << trial << ", " << grip_width << ", " << cur_pose.x << ", " << cur_pose.y << ", " << cur_pose.theta << ", " << extracted_poses[j].get_delta() \
                     << ", " << extracted_poses[j].get_norm_delta() << ", " << extracted_poses[j].get_void_volume() << ", " << extracted_poses[j].v_area_ << ", " << extracted_poses[j].get_sample_time() \
                     << ", " << extracted_poses[j].get_triangulation_time() <<  ", "\
                     << extracted_poses[j].get_filtration_time() << ", " << extracted_poses[j].get_persistence_time()\
                     <<  ", " << extracted_poses[j].get_void_time() << ", " << extracted_poses[j].get_rotation()\
                     << ", " << pic_string << ", " << ppr.path_exists << ", " << ppr.max_energy << ", " << cur_samples << "\n";

	        
	          //Also write out pictures for each grasp...
		  cage_checker.Set_Object_Pose(cur_pose.x, cur_pose.y, cur_pose.theta);
		  im1 = save_image(cage_checker, pic_string, extracted_poses[j].get_rotation());

                }
	      }              
	      of.close();

              grip_width = grip_width + grip_inc;
	      grasp_num++;
            }
          }
          cur_samples = cur_samples * 2;
        }
    }
    delete object;
  }
  delete gripper;

  return EXIT_SUCCESS;
}

