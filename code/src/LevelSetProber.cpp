#include "LevelSetProber.hpp"

#include <boost/bind.hpp>
#include <boost/numeric/ublas/assignment.hpp>
#include <fstream>
#include <iostream>

LevelSetProber::LevelSetProber(bopt_params par, float h, boost::function<double (const vectord& x)> f, int d,
                               boost::function<void (const vectord& x)> opt_cb)
  : ContinuousModel(d, par),
    level_set_(h),
    eval_f_(f),
    cb_f_(opt_cb),
    dim_(d)
{
}

double LevelSetProber::evaluateSample(const vectord& xin)
{
  if (xin.size() != dim_) {
    std::cout << "Warning: Dim " << dim_
              << " specified but queried with " << xin.size()
              << "-dimensional vector" << std::endl;
  }

  double f_x = eval_f_(xin);
  return std::abs(f_x - level_set_);
}

bool LevelSetProber::checkReachability(const vectord& query)
{
  return true; // for now assume all reachable
}

void LevelSetProber::findOptimal(vectord& xOpt)
{
  ContinuousModel::findOptimal(xOpt);
  cb_f_(xOpt); // callback function to update model with optimal probes
  x_seq_.push_back(xOpt);
}

void LevelSetProber::saveSequenceToCsv(std::string filename)
{
  std::ofstream of(filename.c_str());
  for (unsigned int i = 0; i < x_seq_.size(); i++) {
    for (unsigned int j = 0; j < dim_; j++) {
      if (j < dim_ - 1)
        of << x_seq_[i](j) << ",";
      else
        of << x_seq_[i](j) << "\n";          
    }
  }
  of.close();
}

double quad_func(const vectord& x)
{
  float a = 1.0f;
  float b = 1.0f;
  float c = 1.0f;
  return a*x(0)*x(0) + b*x(1)*x(1);// + c*x(2)*x(2);
}

void print_probe(const vectord& x)
{
  std::cout << "Probed " << x(0) << " " << x(1) << " " << x(2) << std::endl;
}

bool RunLevelSetTest()
{
  bopt_params par = initialize_parameters_to_default();
  par.n_iterations = 190;
  par.random_seed = 0;
  par.verbose_level = 1;
  par.noise = 1e-10;
  par.n_init_samples = 20;

  // set gp kernel priors
  set_kernel(&par,"kSEISO");
  par.kernel.hp_mean[0] = 0.0f;
  par.kernel.hp_std[0] = 10.0f;
  par.kernel.n_hp = 1;

  int d = 2;
  float h = 1.0f;
  boost::function<double (const vectord& x)> f = boost::bind(&quad_func, _1);
  boost::function<void (const vectord& x)> cb = boost::bind(&print_probe, _1);
  LevelSetProber ols(par, h, f, d, cb);

  vectord lower_bound(d);
  vectord upper_bound(d);
  lower_bound(0) = -2;
  lower_bound(1) = -2;
  upper_bound(0) = 2;
  upper_bound(1) = 2;
  ols.setBoundingBox(lower_bound, upper_bound);

  vectord result(d);
  ols.optimize(result);
  std::cout << "Result: " << result << "->" 
	    << ols.evaluateSample(result) << std::endl;

  std::string filename = "test.csv";
  ols.saveSequenceToCsv(filename);
  return true;
}
