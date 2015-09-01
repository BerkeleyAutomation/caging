/*
 * SDF class for collision checking
 *
 */
#ifndef SDF_H
#define SDF_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include <opencv2/opencv.hpp>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

#include <CGAL/Timer.h>
typedef CGAL::Timer Timer;

class SDFCollisionCache;

struct Point2ui
{
  unsigned int i;
  unsigned int j;

  Point2ui(unsigned int i_set, unsigned int j_set) {
    i = i_set;
    j = j_set;
  }
};

class SDF
{
 public:
  // create an empty sdf of the given size
  SDF(unsigned int height, unsigned int width);
  // create an sdf from an image file
  SDF(const std::string& filename, Eigen::Matrix3f pose = Eigen::Matrix3f::Identity(), float scale = 1.0f);  
  ~SDF();

 public:
  // renders the SDF to a new one in identity world frame (with interior-only option) 
  bool RenderToWorldFrame(SDF* render_sdf, bool only_interior = false);

  void RenderToRGBImage(cv::Mat& image, char color = 'r');
  
  // check collision between two SDFs
  static float PenetrationDepth(SDF* sdf1, SDF* sdf2, SDFCollisionCache* cache1, SDFCollisionCache* cache2);
  static float Distance(SDF* sdf1, SDF* sdf2, SDFCollisionCache* cache1, SDFCollisionCache* cache2);

 public:
  // return the signed distance at a location
  float At(unsigned int i, unsigned int j);

  // return the signed distance at a non-integer location (using square interpolation)
  float At(float y, float x);

  // set the signed distance at a location
  void Set(float val, unsigned int i, unsigned int j);
  void ComputeMoments();
  void PrerenderSDFToWorldFrame();

  float MaxMomentArm();
  unsigned int Height();
  unsigned int Width();
  Eigen::Vector3f Center();
  Eigen::Matrix3f Pose();
  float Scale();
  bool Initialized() { return initialized_; }

  std::vector<Point2ui> InterestPoints();
  void SetPose(Eigen::Matrix3f pose);
  void ClearInterestPoints();

 private:
  void Initialize();

 private:
  float** sdf_vals_;
  //  bool** medial_axis_;
  float max_moment_arm_;
  unsigned int height_;
  unsigned int width_;
  std::vector<Point2ui> interest_points_;
  Eigen::Vector2f object_center_; // object center in grid
  Eigen::Vector3f center_; // grid center
  Eigen::Matrix3f pose_; // in 'world' frame
  float scale_;
  bool initialized_;
};

struct SDFCollisionCache
{
  bool cache_init;
  SDF* prerendered_sdf;  
  
  SDFCollisionCache() { cache_init = false;  prerendered_sdf = NULL; }
  ~SDFCollisionCache() { delete prerendered_sdf; }
};

#endif // SDF_H
