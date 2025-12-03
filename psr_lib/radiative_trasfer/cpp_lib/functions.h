#pragma once

#include <functional>
#include <Eigen/Dense>
using namespace std;
using Eigen::Vector3d;

/// @brief returns angle between two vecors using numerically stable routine for all possible vector directions
/// @param vec1 
/// @param vec2 
/// @return angle between two 3D vectors
double ANGLE (Vector3d vec1, Vector3d vec2);

/// @brief just simple hyperbolic acrsin function
/// @param x 
/// @return arcsinh(x)
double Arcsinh(double x);

/// @brief just simple signum function
/// @param value
/// @return -1 if x < 0, 1 if x >= 0, 
double sgn (double value);

double get_derivative_along_the_ray(std::function<double(double)> func, double l, double dl=1);
double DX(std::function<double(double, double)> func, double x, double y);
double DY(std::function<double(double, double)> func, double x, double y);


/// @brief split array of pahses into [size] number of chunks for multiprocessing and return 
/// initial and final phases for chunk corrsponding to [rank]
/// @param phi_start_global inital phase for whole computation
/// @param phi_end_global final phase for whole computation
/// @param phi_step phase step
/// @param size number of processes
/// @param rank rank of proces of interest
/// @return initial and final phases. Final phase should not be included into computation
pair<double, double> split_phases(double phi_start_global, double phi_end_global, double phi_step, int size, int rank);