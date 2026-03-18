#include "FixedHeightModelClass.h"
#include "constants.h"
FixedHeightModel::FixedHeightModel(){}

FixedHeightModel::FixedHeightModel(std::map<std::string, double> model_dict, ObservedRadioPulsar PSR)
{
    fr = model_dict["fr"];
    fphi = model_dict["fphi"];
    Rem = model_dict["Rem"];
    lambda = model_dict["lambda"]; // questionable
    gamma0 = model_dict["gamma0"]; // questionable
    PSR_ = PSR;

}

double FixedHeightModel::density_profile(double x_pc, double phi_pc)
{
    // double f = std::pow(x_pc, 2);
    // return (pow(f, 2.5) * exp(-f * f) / (pow(f, 2.5) + pow(0.5, 2.5))) * 1;
    return 100 * std::pow(x_pc, 3) * std::exp(-7 * std::pow(x_pc, 2));
}

Vector3d FixedHeightModel::Bfield(Vector3d vR, Vector3d m){
    double Rdist = vR.norm();
    Vector3d n = vR.normalized();
    Vector3d Bdipole = 3 * m.dot(n) * n - m; // Dipole component
    Vector3d Bwind; //Wind component
    Bwind(0) = Rdist/PSR_.RLC * fr * n(0) - std::pow(Rdist/PSR_.RLC, 2)* 
    fphi * fr * (-n(1));
    Bwind(1) = Rdist/PSR_.RLC * fr * n(1) - std::pow(Rdist/PSR_.RLC, 2)* 
    fphi * fr * n(0);
    Bwind(2) = Rdist/PSR_.RLC * fr * n(2);
    return Bdipole + Bwind;
}

double FixedHeightModel::fDist (double gamma) {
  return ((6.0 * gamma0) / (std::pow(2.0, 1.0/6.0) * constants::PI)) * (std::pow(gamma, 4) / (2.0 * std::pow(gamma, 6) + std::pow(gamma0, 6)));
}

double FixedHeightModel::INTEGRAL (double gamma, double A) { // THIS IS TRASH. REDO AS FAST AS POSSIBLE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  return -(std::pow(2, 2.0 / 3.0)*(-2*std::sqrt(3)*std::atan((2*std::pow(2, 1.0 / 3.0)*std::pow(gamma,2) - std::pow(gamma0,2))/
          (std::sqrt(3)*std::pow(gamma0,2))) - 2*std::log(std::pow(2, 1.0 / 3.0)*std::pow(gamma,2) + std::pow(gamma0,2)) +
          std::log(std::pow(2, 2.0 / 3.0)*std::pow(gamma,4) - std::pow(2, 1.0 / 3.0)*std::pow(gamma,2)*std::pow(gamma0,2) +
          std::pow(gamma0,4))) - std::pow(2, 1.0 / 3.0)*std::pow(gamma0,2)*A*
          (2*std::sqrt(3)*std::atan((2*std::pow(2, 1.0 / 3.0)*std::pow(gamma,2) - std::pow(gamma0,2))/(std::sqrt(3)*std::pow(gamma0,2))) -
          2*std::log(std::pow(2, 1.0 / 3.0)*std::pow(gamma,2) + std::pow(gamma0,2)) +
          std::log(std::pow(2, 2.0 / 3.0)*std::pow(gamma,4) - std::pow(2, 1.0 / 3.0)*std::pow(gamma,2)*std::pow(gamma0,2) +
          std::pow(gamma0,4))) - 2*std::pow(gamma0,4)*std::pow(A,2)*
          (std::log(2*std::pow(gamma,6) + std::pow(gamma0,6)) - 3*std::log(fabs((1 - gamma*std::sqrt(A))*(1 + gamma*std::sqrt(A))))))/
          (2.*std::pow(2, 1.0 / 6.0)*constants::PI*std::pow(gamma0,3)*(2 + std::pow(gamma0,6)*std::pow(A,3)));
}


double FixedHeightModel::Lambda_type_avrg(double A)
{
    return INTEGRAL(1000000, A) - INTEGRAL(0, A); // TRASH CONTINUES. REDO AS FAST AS POSSIBLE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
}

double FixedHeightModel::ImEps_type_avrg(double A)
{
    return gamma0;
}

double FixedHeightModel::Q_type_avrg(double A)
{
    // should be Lambda_type_avrg / ImEps_type_avrg, but now for test it ist just 1/gamma^3 (cold plasma)
    return std::pow(1/gamma0, 3);
}
