#include "SDF.h"

#include <fstream>
#include <sstream>

SDF::SDF(unsigned int height, unsigned int width)
  : height_(height),
    width_(width),
    pose_(Eigen::Matrix3f::Identity()),
    scale_(1.0f),
    initialized_(false)
{
  Initialize();
  ComputeMoments();
}

SDF::SDF(const std::string& filename, Eigen::Matrix3f pose, float scale)
  : height_(0),
    width_(0),
    pose_(pose),
    scale_(scale),
    initialized_(false)
{
  std::ifstream ifile(filename.c_str());
  if (!ifile.is_open()) {
    std::cout << "Error: failed to open file " << filename << std::endl; 
    return;
  }

  int max_len = 100000;
  char dummy;
  char buffer[max_len];
  std::string buffer2;
  interest_points_.clear();

  ifile >> dummy >> height_ >> width_;
  ifile.getline(buffer, max_len); // remove next line

  Initialize();

  // read csv to file
  unsigned int i = 0;
  float sdf_val;
  while (!ifile.eof() && i < height_) {
    ifile.getline(buffer, max_len);
    std::stringstream stream(buffer);

    for (unsigned int j = 0; j < width_; j++) {
      std::getline(stream, buffer2, ',');
      sdf_vals_[i][j] = atof(buffer2.c_str());
      
      if (sdf_vals_[i][j] < 0.0f) {
        interest_points_.push_back(Point2ui(i,j));
      }
    }
    i++;
  }
  ifile.close();

  ComputeMoments();
  initialized_ = true;
}

SDF::~SDF()
{
  for (unsigned int i = 0; i < height_; i++) {
    delete [] sdf_vals_[i];
  }
  delete [] sdf_vals_;
}

void SDF::Initialize()
{
  // center of grid
  center_ << (float)(width_ - 1) / 2.0f, (float)(height_ - 1) / 2.0f, 0;

  // index rows, then column
  sdf_vals_ = new float*[height_];

  for (unsigned int i = 0; i < height_; i++) {
    sdf_vals_[i] = new float[width_];
  }
}

float SDF::At(unsigned int i, unsigned int j)
{
  // snap to boundary
  if (i >= height_) {
    i = height_ - 1;
  }
  if (j >= width_) {
    j = width_ - 1;
  }

  return sdf_vals_[i][j];
}

float SDF::At(float y, float x)
{
  int x_min = std::min((float)width_-1.0f, std::max(0.0f, std::floor(x)));
  int y_min = std::min((float)height_-1.0f, std::max(0.0f, std::floor(y)));
  int x_max = std::max(0.0f, std::min((float)width_-1.0f, std::ceil(x)));
  int y_max = std::max(0.0f, std::min((float)height_-1.0f, std::ceil(y)));
  float value = 0.0f;
  float weight_sum = 0.0f;

  int x_coords[4] = {x_min, x_min, x_max, x_max};
  int y_coords[4] = {y_min, y_max, y_min, y_max};
  float weights[4];
  float eps = 1e-2;

  // weighted sum of four corner points
  for (int i = 0; i < 4; i++) {
    float dist = (x_coords[i] - x) * (x_coords[i] - x) + (y_coords[i] - y) * (y_coords[i] - y);
    float weight = 1.0f / eps; // max_weight
    if (dist > eps) {
      weight = 1.0f / dist;
    }
    //    std::cout << y_coords[i] << " " << x_coords[i] << std::endl;
    value += weight * sdf_vals_[y_coords[i]][x_coords[i]];
    weight_sum += weight;
  }

  value = value / weight_sum;
  return value;
}

void SDF::Set(float val, unsigned int i, unsigned int j)
{
  if (i >= height_ || j >= width_) {
    return;
  }
  sdf_vals_[i][j] = val;

  if (sdf_vals_[i][j] < 0.0f) {
    interest_points_.push_back(Point2ui(i,j));
  }
}

void SDF::ComputeMoments()
{
  int num_interior = 0;
  Eigen::Vector2f p;
  object_center_ << 0.0f, 0.0f;

  // compute the object center
  for (unsigned int i = 0; i < height_; i++) {
    for (unsigned int j = 0; j < width_; j++) {
      if (At(i,j) <= 0.0f) {
        p << (float)j - center_(0), (float)i - center_(1);
        object_center_ = object_center_ + p;
        num_interior++;
      }
    }
  }
  object_center_ = object_center_ / num_interior;

  max_moment_arm_ = 0.0f;
  for (unsigned int i = 0; i < height_; i++) {
    for (unsigned int j = 0; j < width_; j++) {
      if (At(i,j) <= 0.0f) {
        p << (float)j - center_(0), (float)i - center_(1);
        p = p - object_center_;
        if (p.norm() > max_moment_arm_) {
          max_moment_arm_ = p.norm();
        }
      }
    }
  }
}

float SDF::MaxMomentArm()
{
  return max_moment_arm_;
}

unsigned int SDF::Height()
{
  return height_;
}

unsigned int SDF::Width()
{
  return width_;
}

Eigen::Vector3f SDF::Center()
{
  return center_;
}

Eigen::Matrix3f SDF::Pose()
{
  return pose_;
}

float SDF::Scale()
{
  return scale_;
}

std::vector<Point2ui> SDF::InterestPoints()
{
  return interest_points_;
}

void SDF::SetPose(Eigen::Matrix3f pose)
{
  pose_ = pose;
}

void SDF::ClearInterestPoints() {
  interest_points_.clear();
}

bool SDF::RenderToWorldFrame(SDF* render_sdf, bool only_interior)
{
  Timer timer;
  timer.start();
  if (render_sdf == NULL || render_sdf->Height() != height_ || render_sdf->Width() != width_) {
    return false;
  }

  float val;
  Eigen::Vector3f v;
  Eigen::Vector3f v_trans;
  render_sdf->ClearInterestPoints();

  if (only_interior) {
    // only transform the interior points
    // useful for fast coll checking but less accurate due to rounding
    for (unsigned int i = 0; i < interest_points_.size(); i++) {
      Point2ui p = interest_points_[i];
      v << p.j, p.i, 1;
      v = v - center_; // normalize to center of grid
      v_trans = scale_ *  pose_ * v;
      v_trans  = v_trans + render_sdf->Center();      

      val = At(p.i, p.j);
      render_sdf->Set(val, v_trans(1), v_trans(0));
    }
  }
  else {
    // transform the entire grid
    for (unsigned int i = 0; i < height_; i++) {
      for (unsigned int j = 0; j < width_; j++) {
        // location in render sdf
        v << j, i, 1;
        v = v - render_sdf->Center(); // normalize to center of grid
        v_trans = (1.0f / scale_) *  pose_.inverse() * v;
        v_trans  = v_trans + center_;      

        val = At(v_trans(1), v_trans(0));
        render_sdf->Set(val, i, j);
      }    
    }
  }
  timer.stop();
  //  std::cout << "Rendering world frame took " << timer.time() << " sec" << std::endl;
  return true;
}

void SDF::RenderToRGBImage(cv::Mat& image, char color)
{
  SDF* sdf_render = new SDF(height_, width_);
  RenderToWorldFrame(sdf_render);
  int image_height = image.size().height;
  int image_width = image.size().width;
  Eigen::Vector3f image_center;
  Eigen::Vector3f v;
  image_center << ((float)image_width - 1.0f) / 2.0f, ((float)image_height - 1.0f) / 2.0f, 0;

  for (unsigned int i = 0; i < image_height; i++) {
    for (unsigned int j = 0; j < image_width; j++) {
      
      // render color where inside shape (TODO: change to negative at some point)
      v << j, i, 1;
      v = v - image_center; // normalize to center of grid
      v = v + center_;
      if (sdf_render->At(v(1), v(0)) < 0.0f) {
        if (color == 'r') {
          image.at<cv::Vec3b>(cv::Point(j,i)) = cv::Vec3b(0, 0, 255);
        }
        else if (color == 'g') {
          image.at<cv::Vec3b>(cv::Point(j,i)) = cv::Vec3b(0, 255, 0);
        }
        else {
          image.at<cv::Vec3b>(cv::Point(j,i)) = cv::Vec3b(255, 0, 0);
        }
      }
    }
  }
  delete sdf_render;
}

float SDF::PenetrationDepth(SDF* sdf1, SDF* sdf2, SDFCollisionCache* cache1, SDFCollisionCache* cache2)
{
  Timer timer;
  timer.start();
  if (!cache1->cache_init) {
    cache1->prerendered_sdf = new SDF(sdf1->Height(), sdf1->Width());
    if (sdf1->RenderToWorldFrame(cache1->prerendered_sdf, false)) {
      cache1->cache_init = true;
    }
  }
  if (!cache2->cache_init) {
    cache2->prerendered_sdf = new SDF(sdf2->Height(), sdf2->Width());
    if (sdf2->RenderToWorldFrame(cache2->prerendered_sdf, false)) {
      cache2->cache_init = true;
    }
  }
  float pen_depth = 0.0f;
  std::vector<Point2ui> interest_points = cache1->prerendered_sdf->InterestPoints();

  for (unsigned int i = 0; i < interest_points.size(); i++) {
    Point2ui p = interest_points[i];
    float sd2 = cache2->prerendered_sdf->At(p.i,p.j);
  
    // check both inside of shape (TODO: change to negative at some point)
    if (-sd2 > pen_depth) {
        pen_depth = -sd2; 
    }
  }
  timer.stop();
  //  std::cout << "Pen depth took " << timer.time() << " sec" << std::endl;

  return pen_depth;
}

float SDF::Distance(SDF* sdf1, SDF* sdf2, SDFCollisionCache* cache1, SDFCollisionCache* cache2)
{
  if (!cache1->cache_init) {
    cache1->prerendered_sdf = new SDF(sdf1->Height(), sdf1->Width());
    if (sdf1->RenderToWorldFrame(cache1->prerendered_sdf, false)) {
      cache1->cache_init = true;
    }
  }
  if (!cache2->cache_init) {
    cache2->prerendered_sdf = new SDF(sdf2->Height(), sdf2->Width());
    if (sdf2->RenderToWorldFrame(cache2->prerendered_sdf, false)) {
      cache2->cache_init = true;
    }
  }
  float dist = FLT_MAX;
  std::vector<Point2ui> interest_points = cache1->prerendered_sdf->InterestPoints();

  for (unsigned int i = 0; i < interest_points.size(); i++) {
    Point2ui p = interest_points[i];
    float sd2 = cache2->prerendered_sdf->At(p.i,p.j);

    // check both inside of shape (TODO: change to negative at some point)
    if (sd2 > 0.0f && sd2 < dist) {
      dist = sd2; 
    }
  }
  return dist;
}


