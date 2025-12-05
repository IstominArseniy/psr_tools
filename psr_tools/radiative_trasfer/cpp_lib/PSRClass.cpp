#include "PSRClass.h"
#include "constants.h"
#include "read_write.h"

using Eigen::Vector3d;

ObservedRadioPulsar::ObservedRadioPulsar(){}

ObservedRadioPulsar::ObservedRadioPulsar(std::map<std::string, double> psr_dict)
{
    B12 = psr_dict["B12"];
    freqGHz = psr_dict["freqGHz"];
    Period = psr_dict["Period"];
    chi_deg = psr_dict["chi_deg"];
    beta_deg = psr_dict["beta_deg"];
    Rs = psr_dict["Rs"];
    initialize_derived_params();
}

void ObservedRadioPulsar::init_from_file(std::string filename)
{
    std::map<std::string, double> psr_dict;
    B12 = read_from_file(filename, "B12"); // Surface B-field in 10^12 Gs
    Period = read_from_file(filename, "Period"); // Rotation period in sec
    freqGHz = read_from_file(filename, "freqGHz"); // Radiation frequency in GHz
    chi_deg = read_from_file(filename, "alpha_deg"); // alpha_deg name is for backward compatability ->change to chi_deg
    beta_deg = read_from_file(filename, "beta_deg");
    Rs = 1.2e6; // Radius of the star. Almost never changed
    initialize_derived_params();
}

void ObservedRadioPulsar::initialize_derived_params()
{
    omega_obs = 2 * constants::PI * freqGHz * 1e9;
    Omega = 2 * constants::PI / Period;
    chi = chi_deg * constants::PI / 180;
    beta = beta_deg * constants::PI / 180;
    dzeta = chi - beta;
    RLC = (constants::c / Omega) / Rs;
    Rpc = std::sqrt(1/RLC);
    Omega_vec << 0, 0, Omega;
    observer_vec << std::sin(dzeta), 0, std::cos(dzeta);
}


