#include "Util.h"

// helper functions
double Distance_Squared(DT_Vector3 p1, DT_Vector3 p2)
{
  double ret_dist = 0.0;
  for(int i = 0; i < 3; i++) {
    ret_dist += (p1[i]-p2[i])*(p1[i]-p2[i]);
  }
  return ret_dist;
}

double Magnitude_Squared(const DT_Vector3 p)
{
  double ret_dist = 0.0;
  for(int i = 0; i < 3; i++) {
    ret_dist += p[i]*p[i];
  }
  return ret_dist;
}

double Magnitude_Squared(DT_Vector3 p)
{
  double ret_dist = 0.0;
  for(int i = 0; i < 3; i++) {
    ret_dist += p[i]*p[i];
  }
  return ret_dist;
}

Eigen::Matrix4f CreatePose(float tx, float ty, float theta)
{
  Eigen::Matrix4f pose;
  pose << cos(theta), -sin(theta), 0, ty,
    sin(theta), cos(theta), 0 , tx,
    0, 0, 1, 0,
    0, 0, 0, 1;
  return pose;
}

Eigen::Matrix3f CreatePose2D(float tx, float ty, float theta)
{
  Eigen::Matrix3f pose;
  pose << cos(theta), -sin(theta), tx,
    sin(theta), cos(theta), ty,
    0, 0, 1;
  return pose;
}

void Pose2DToParams2D(Eigen::Matrix3f pose, float& tx, float& ty, float& theta)
{
  tx = pose(0,2);
  ty = pose(1,2);

  Eigen::Vector2f v, w;
  v << pose(0,0), pose(1,0);
  w << 1, 0;
  theta = acos(std::max<float>(std::min<float>(v.dot(w), 1.0f), -1.0f));
  if (v(1) < 0) {
    theta = -theta;
  }
}

void Pose3DToParams2D(Eigen::Matrix4f pose, float& tx, float& ty, float& theta)
{
  tx = pose(1,3);
  ty = pose(0,3);

  Eigen::Vector2f v, w;
  v << pose(0,0), pose(1,0);
  w << 1, 0;
  theta = acos(std::max<float>(std::min<float>(v.dot(w), 1.0f), -1.0f));
  if (v(1) < 0) {
    theta = -theta;
  }
}

void UniformRandomPose(float& tx, float& ty, float& theta,
                       float min_tx, float max_tx,
                       float min_ty, float max_ty,
                       float min_theta, float max_theta)
{
  tx = (max_tx - min_tx) * ((float)rand() / (float)RAND_MAX) + min_tx;
  ty = (max_ty - min_ty) * ((float)rand() / (float)RAND_MAX) + min_ty;
  theta = (max_theta - min_theta) * ((float)rand() / (float)RAND_MAX) + min_theta;
}


CGAL_Aff_Transform EigenPoseToCGAL(Eigen::Matrix4f pose)
{
  CGAL_Aff_Transform ret_trans(pose(0,0), pose(0,1), pose(0,2), pose(0,3),
                               pose(1,0), pose(1,1), pose(1,2), pose(1,3),
                               pose(2,0), pose(2,1), pose(2,2), pose(2,3));
  return ret_trans;
}

// Modified from FCL test utility
bool LoadOBJFile(const char* filename, std::vector<fcl::Vec3f>& points, std::vector<fcl::Triangle>& triangles)
{
  FILE* file = fopen(filename, "rb");
  if(!file) {
    std::cerr << "File " << filename << " does not exist" << std::endl;
    return false;
  }

  bool has_normal = false;
  bool has_texture = false;
  char line_buffer[2000];
  while(fgets(line_buffer, 2000, file)) {
    char* first_token = strtok(line_buffer, "\r\n\t ");
    if(!first_token || first_token[0] == '#' || first_token[0] == 0)
      continue;

    switch(first_token[0])
    {
    case 'v':
      {
        if(first_token[1] == 'n') {
          strtok(NULL, "\t ");
          strtok(NULL, "\t ");
          strtok(NULL, "\t ");
          has_normal = true;
        }
        else if(first_token[1] == 't') {
          strtok(NULL, "\t ");
          strtok(NULL, "\t ");
          has_texture = true;
        }
        else {
          float x = (float)atof(strtok(NULL, "\t "));
          float y = (float)atof(strtok(NULL, "\t "));
          float z = (float)atof(strtok(NULL, "\t "));
          fcl::Vec3f p(x, y, z);
          points.push_back(p);
        }
      }
      break;
    case 'f':
      {
        fcl::Triangle tri;
        char* data[30];
        int n = 0;
        while((data[n] = strtok(NULL, "\t \r\n")) != NULL) {
          if(strlen(data[n]))
            n++;
        }

        for(int t = 0; t < (n - 2); ++t) {
          if((!has_texture) && (!has_normal)) {
            tri[0] = atoi(data[0]) - 1;
            tri[1] = atoi(data[1]) - 1;
            tri[2] = atoi(data[2]) - 1;
          }
          else {
            const char *v1;
            for(int i = 0; i < 3; i++) {
              // vertex ID
              if(i == 0)
                v1 = data[0];
              else
                v1 = data[t + i];
              
              tri[i] = atoi(v1) - 1;
            }
          }
          triangles.push_back(tri);
        }
      }
    }
  }

  int closed = fclose(file);
  if (closed != 0) {
    std::cerr << "Failed to close " << filename << std::endl;
    return false;
  }
  return true;
}
