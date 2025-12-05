#include <map>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

#include "constants.h"
#include "functions.h"
#include "CppInterfaceClass.h"
#include "SolverClass.h"
#include "read_write.h"


using Eigen::Vector3d;

CppInterface::CppInterface(std::string log_path){
    log_path_ = log_path;
}

CppInterface::CppInterface(std::map<std::string, double> psr_dict, std::map<std::string, double> param_dict, std::string log_path)
{ // Constructor form dict
    PSR_ = ObservedRadioPulsar(psr_dict);
    model_ = FixedHeightModel(param_dict, PSR_);
    log_path_ = log_path;
}

CppInterface::CppInterface(ObservedRadioPulsar PSR, std::map<std::string, double> param_dict, std::string log_path)
{
    PSR_ = PSR;
    model_ = FixedHeightModel(param_dict, PSR_);
    log_path_ = log_path;
}

void CppInterface::init_from_file(std::string filename)
{    
    std::map<std::string, double> psr_dict;
    psr_dict["B12"] = read_from_file(filename, "B12"); // Surface B-field in 10^12 Gs
    psr_dict["Period"] = read_from_file(filename, "Period"); // Rotation period in sec
    psr_dict["freqGHz"] = read_from_file(filename, "freqGHz"); // Radiation frequency in GHz
    psr_dict["chi_deg"] = read_from_file(filename, "alpha_deg"); // alpha_deg name is for backward compatability ->change to chi_deg
    psr_dict["beta_deg"] = read_from_file(filename, "beta_deg");
    psr_dict["Rs"] = 1.2e6; // Radius of the star. Almost never changed
    std::map<std::string, double> model_dict;
    model_dict["fr"] = read_from_file(filename, "fr");
    model_dict["fphi"] = read_from_file(filename, "fphi");
    model_dict["Rem"] = read_from_file(filename, "R_em");
    model_dict["lambda"] = read_from_file(filename, "lambda");
    model_dict["gamma0"] = read_from_file(filename, "gamma0");
    PSR_ = ObservedRadioPulsar(psr_dict);
    model_ = FixedHeightModel(model_dict, PSR_);
}

std::map<std::string, double> CppInterface::find_ILVPA(double phi, int mode, bool with_absorption)
{
    FixedHeightSolver solver(phi, mode, PSR_, model_);
    double l1 = solver.find_initial_point(false);
    double l2 = std::min(2.5 * get_R_escape(), 2*PSR_.RLC);
    solver.write_params_on_ray(log_path_);
    std::vector<double> theta_init = solver.find_approximate_KO_solution(l1);
    std::vector<double> theta_final = solver.solve_KO_equations(theta_init, l1, l2, log_path_);
    double I = solver.find_intensity();
    if(with_absorption==true){
        double tau = solver.get_tau();
        I *= std::exp(-tau);
    }
    double V = I * std::tanh(2.0 * theta_final[1]);
    double PA = theta_final[0] * 180.0 / constants::PI;
    std::map<std::string, double> result;
    result["phase"] = phi;
    result["I"] = I;
    result["L"] = std::sqrt(I*I - V*V);
    result["V"] = -V;
    result["PA"] = PA;
    return result;
}

double CppInterface::find_I(double phi, int mode, bool with_absorption)
{
    FixedHeightSolver solver(phi, mode, PSR_, model_);
    double I = solver.find_intensity();
    if(with_absorption){
        double tau = solver.get_tau();
        I *= std::exp(-tau);
    }
    return I;
}

std::vector<std::map<std::string, double> > CppInterface::calculate_profile(double phi1, double phi2, double phi_step, int mode, bool with_absorption){
    std::vector<std::map<std::string, double> > profile;
    double phi_tmp = phi1;
    while(phi_tmp < phi2){
        auto result = find_ILVPA(phi_tmp, mode, with_absorption);
        profile.push_back(result);
        phi_tmp += phi_step;
    }
    return profile;
}



std::vector<double> CppInterface::vectorize_data(std::vector<std::map<std::string, double> > data)
{
    int N_counts = data.size();
    std::vector<double> vectorized_data(N_counts * 5);
    for(int i=0; i<N_counts; i++){
        vectorized_data[i*5] = data[i]["phase"];
        vectorized_data[i*5+1] = data[i]["I"];
        vectorized_data[i*5+2] = data[i]["L"];
        vectorized_data[i*5+3] = data[i]["V"];
        vectorized_data[i*5+4] = data[i]["PA"];
    }
    return vectorized_data;
}

std::vector<std::map<std::string, double> > CppInterface::restore_data(std::vector<double> vectorized_data)
{   
    int N_counts = vectorized_data.size() / 5;
    std::vector<std::map<std::string, double> > restored_data(N_counts);
    for(int i=0; i<N_counts; i++){
        restored_data[i]["phase"] = vectorized_data[i*5];
        restored_data[i]["I"] = vectorized_data[i*5+1];
        restored_data[i]["L"] = vectorized_data[i*5+2];
        restored_data[i]["V"] = vectorized_data[i*5+3];
        restored_data[i]["PA"] = vectorized_data[i*5+4];
    }
    return restored_data;
}

double CppInterface::get_R_escape()
{
    return 1.0e3 * std::pow(model_.lambda / 1.0e4, 1.0/3.0) * std::pow(model_.gamma0 / 100.0, -6.0/5.0) * std::pow(PSR_.B12, 2.0/5.0) * std::pow(PSR_.freqGHz, -2.0/5.0) * std::pow(PSR_.Period, -1.0/5.0);
}

double CppInterface::get_RLC()
{
    return PSR_.RLC;
}

double CppInterface::get_rho(){
    return 3/2.0*std::sqrt(model_.Rem)*PSR_.Rpc * 180 / constants::PI;
}
