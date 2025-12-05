#include <vector>
#include <math.h>
#include <Eigen/Dense>
#include <functional>
#include "functions.h"

using Eigen::Vector3d;

double ANGLE (Vector3d vec1, Vector3d vec2) {
  return 2 * std::atan2((vec2.normalized() - vec1.normalized()).norm(), (vec2.normalized() + vec1.normalized()).norm());
} 

double Arcsinh(double x) {
  return std::log(x + std::sqrt(std::pow(x, 2.0) + 1.0));
}

double sgn (double value) {
  if (value >= 0.0) {
    return 1.0;
  } else {
    return -1.0;
  }
}

double get_derivative_along_the_ray(std::function<double(double)> func, double l, double dl)
{
  return (func(l + dl) - func(l)) / dl; // Very approximate, but high accuracy is not needed (probably better still to make at least 2nd order)
}

double DX(std::function<double(double, double)> func, double x, double y) {
    double h = 0.00001;
    double fm2 = func(x - 2 * h, y);
    double fp2 = func(x + 2 * h, y);
    double fm1 = func(x - h, y);
    double fp1 = func(x + h, y);
    return (fm2 - 8 * fm1 + 8 * fp1 - fp2) / (12 * h);
}
double DY (std::function<double(double, double)> func, double x, double y) {
    double h = 0.00001;
    double fm2 = func(x, y - 2 * h);
    double fp2 = func(x, y + 2 * h);
    double fm1 = func(x, y - h);
    double fp1 = func(x, y + h);
    return (fm2 - 8 * fm1 + 8 * fp1 - fp2) / (12 * h);
}

pair<double, double> split_phases(double phi_start_global, double phi_end_global, double phi_step, int size, int rank){
  /*
  function to find whitch phases should be processed by this process 
  returns phi_initial and phi_final for process
  */
  pair <double, double> phases;
  int Nsteps_global = int((phi_end_global - phi_start_global) / phi_step);
  int Nsteps = Nsteps_global / size;
  int residual = Nsteps_global % size;
  if(rank < residual){
    phases.first = phi_start_global + rank * (Nsteps + 1) * phi_step;
    phases.second = phi_start_global + (rank + 1) * (Nsteps + 1) * phi_step;
  }
  else{
    phases.first = phi_start_global + ((Nsteps + 1) * residual + Nsteps * (rank - residual)) * phi_step;
    phases.second = phi_start_global + ((Nsteps + 1) * residual + Nsteps * (rank - residual + 1)) * phi_step;
  }   
  return phases;
}