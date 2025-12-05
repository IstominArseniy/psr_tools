#include <Eigen/Dense>
#include <map>
#include <string>

#pragma once

using Eigen::Vector3d;

class ObservedRadioPulsar{
    public:
    ObservedRadioPulsar();
    ObservedRadioPulsar(std::map<std::string, double> psr_dict); // constructor from dictionary
    void init_from_file(std::string filename);
    double B12;
    double freqGHz;
    double omega_obs;
    double Period;
    double Omega;
    double chi_deg;
    double beta_deg;
    double chi;
    double beta;
    double dzeta;
    double Rs;
    double RLC;
    double Rpc;
    Vector3d observer_vec;
    Vector3d Omega_vec;
    private:
    void initialize_derived_params();
};
