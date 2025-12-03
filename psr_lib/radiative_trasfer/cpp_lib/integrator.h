#pragma once

#include <functional>

// double integrate (double (*func)(double), double a, double b);
double integrate (std::function<double(double)> func, double a, double b);
