// Class for level set probing

#ifndef __LEVEL_SET_H__
#define __LEVEL_SET_H__

#define _USE_MATH_DEFINES
#include <cmath>          
#include <algorithm>    

#include <boost/function.hpp>
#include <bayesopt/bayesopt.hpp>

class LevelSetProber : public bayesopt::ContinuousModel
{
 public:
  LevelSetProber(bopt_params par, float h, boost::function<double (const vectord& x)> f, int d,
                 boost::function<void (const vectord& x)> cb);
 
 public:
  double evaluateSample(const vectord& xin);
  bool checkReachability(const vectord& query);
  void findOptimal(vectord& xOpt);

  void saveSequenceToCsv(std::string filename);

private:
  float level_set_; // level set to probe
  boost::function<float (const vectord& x)> eval_f_; // function to evaluate a point
  boost::function<void (const vectord& x)> cb_f_;   // function to call back when a poitn is evaluated
  int dim_; // dimension of acceptable vectors
  std::vector<vectord> x_seq_;
};

// Testing
bool RunLevelSetTest();

#endif // LEVEL_SET_H
