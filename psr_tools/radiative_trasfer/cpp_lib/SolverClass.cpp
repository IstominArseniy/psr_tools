#include <boost/numeric/odeint.hpp>

#include "SolverClass.h"
#include "functions.h"
#include "integrator.h"
#include "constants.h"

#include <fstream>

using Eigen::Vector3d;
using namespace boost::numeric::odeint;

FixedHeightSolver::FixedHeightSolver(double phi_deg, int mode, ObservedRadioPulsar &PSR, FixedHeightModel &model)
{
    model_ = model;
    PSR_ = PSR;
    phi_ = phi_deg * constants::PI / 180;
    mode_ = mode;
    std::pair<double, double> emission_coords = find_emission_point();
    theta_em_ = emission_coords.first;
    phi_em_ = emission_coords.second;
}

std::vector<double> FixedHeightSolver::solve_KO_equations(std::vector<double> theta_initial, double l1, double l2, std::string log_path)
{
  std::vector<double> dep_vars = theta_initial;
  double eps_abs = 1e-6, eps_rel = 1e-6, h_init = 1.0e-8;
  auto addaptive_stepper = make_controlled(eps_abs, eps_rel, runge_kutta_dopri5 < std::vector<double> >()); // make stepper for ODE integration
  auto ode_range = make_adaptive_time_range(addaptive_stepper, [this](const std::vector<double>& f, std::vector<double> &dydx, double l){this->RHS_for_boost(f, dydx, l);}, dep_vars, l1, l2, h_init); // make range for integration
  auto it=ode_range.first;
  if (log_path != ""){ // with logs
    ofstream log_stream(log_path + "/log_" + std::to_string(std::round(100 * phi_*180/constants::PI)/100));
    while(it != ode_range.second){ // integration
    if(stop_condition(it->second)){ // stop integration if Udr > c(?!)
      break;
    }
    log_stream << it->second <<  ", " << it->first[0] << ", " << it->first[1] << ", " << PSR_.Rs * PSR_.omega_obs / (2.0 * constants::c) * Lambda(it->second) << ", " << BetaB (it->second) + delta (it->second) << ", " <<BetaB (it->second) << ", " << delta (it->second) << ", " << vUdr(it->second)(0) << ", " << vUdr(it->second)(1) << std::endl;
    // log_stream << -2*vUdr(it->second)(1)*std::cos(theta_kb(it->second))*(std::sin(theta_kb(it->second))-vUdr(it->second)(0)) / ((std::sin(theta_kb(it->second))-vUdr(it->second)(0))*(std::sin(theta_kb(it->second))-vUdr(it->second)(0)) - std::cos(theta_kb(it->second))*std::cos(theta_kb(it->second))*vUdr(it->second)(1)*vUdr(it->second)(1)) << std::endl;
    it++; // make integration step
    }
    log_stream.close();
  }
  else{ //no logs
    while(it != ode_range.second){ // integration 
      if(stop_condition(it->second)){ // stop integration if Udr > c(?!)
        break;
      }
      it++; // make integration step
    }
  }
 return dep_vars;
}


bool FixedHeightSolver::stop_condition(double l){
  return false;
}

double FixedHeightSolver::Q_denominator(double l)
{
  double vx = vUdr(l)(0);
  double vy = vUdr(l)(1);
  double vz = vUdr(l)(2);
  double sinth = std::sin(theta_kb(l));
  double costh = std::cos(theta_kb(l));
  return (costh * (1 - vx * vx - vy * vy) - vz * (1.0 - sinth * vx));
}

void FixedHeightSolver::RHS_for_boost(const std::vector<double>& f, std::vector<double> &dydx, double l) {
	double coeff = PSR_.Rs * PSR_.omega_obs / (2.0 * constants::c);

	double LL = Lambda (l);
	double QQ = Q (l);
	double BB = BetaB (l);
	double DD = delta (l);

	dydx[0] = coeff * (-LL / QQ - LL * std::cos(2 * f[0] - 2 * BB - 2 * DD) * std::sinh(2 * f[1]));
	dydx[1] = coeff * LL * std::sin(2 * f[0] - 2 * BB - 2 * DD) * std::cosh(2 * f[1]);
}

std::vector<double> FixedHeightSolver::find_approximate_KO_solution(std::vector<double> theta_initial)
{
  // SEEMS TO BE IMPOSSIBLE
  return std::vector<double>();
}

std::vector<double> FixedHeightSolver::find_approximate_KO_solution(double l)
{
  std::vector<double> theta_final(2);
  double Lambda_integral;
  double starting_point=0;
  Lambda_integral = integrate([this](double l){return this->Lambda(l);}, starting_point, l);
  double coefficient = constants::c / PSR_.Rs / PSR_.omega_obs;
  double delta_theta = - coefficient / Lambda(starting_point) * (BetaB_derivative(starting_point) + delta_derivative(starting_point)) * std::sin(Lambda_integral / coefficient);
  double theta2 = -coefficient / Lambda(l) * (BetaB_derivative(l) + delta_derivative(l)) - 1 / 2 / Q(l) + 
  coefficient / Lambda(starting_point) * (BetaB_derivative(0) + delta_derivative(starting_point)) * std::cos(Lambda_integral / coefficient);
  if (mode_ == 0){
    theta_final[0] = constants::PI / 2 + BetaB(l) + delta(l) + delta_theta;
    theta_final[1] = -theta2;
  }
  else{
    theta_final[0] = BetaB(l) + delta(l) + delta_theta;
    theta_final[1] = theta2;
  }
  return theta_final;
}


double FixedHeightSolver::find_initial_point(bool use_binary_search) {
  double freq0 = 0.1;
  if(use_binary_search){ // binary search
    double n_iter=0;

    if(std::abs(Lambda_derivative(0) / std::pow(Lambda(0), 2) * 2 * constants::c / PSR_.Rs / PSR_.omega_obs) > freq0)
      return 0; 
    
    double l_left = 0, l_right = PSR_.RLC / 10, l_cur;   
    l_cur = (l_left + l_right) / 2;

    while(std::abs(std::abs(Lambda_derivative(l_cur) / std::pow(Lambda(l_cur), 2) * 2 * constants::c  / PSR_.Rs/ PSR_.omega_obs)  - freq0) > 0.01 && n_iter < 30){
      l_cur = (l_left + l_right) / 2;
      n_iter++;
      if(std::abs(Lambda_derivative(l_cur) / std::pow(Lambda(l_cur), 2) * 2 * constants::c / PSR_.Rs / PSR_.omega_obs) > freq0){
        l_right = l_cur;
      }
      else{
        l_left = l_cur; 
      }
    }
    return l_cur;
  }
  else{ // linear search
    double cr_l = 0, step = 10;
    while(std::abs(Lambda_derivative(cr_l + step) / std::pow(Lambda(cr_l+step), 2) * 2 * constants::c / PSR_.Rs / PSR_.omega_obs) < freq0){
      cr_l += step;
    }
    return cr_l;
  }
}

double FixedHeightSolver::find_intensity()
{
    return gFunc(0);
}

void FixedHeightSolver::write_params_on_ray(std::string log_path)
{
    if (log_path == "")
      return;
    std::ofstream param_stream(log_path + "/prop_params_" + std::to_string(std::round(100 * phi_*180/constants::PI)/100));
    double l = 0;
    double dl = 2;
    while (l <= 1.5 * PSR_.RLC){
        param_stream << l << ", " << Lambda(l) << ", " << delta(l) << ", " << BetaB(l) << ", " << vUdr(l)(0) << ", " << vUdr(l)(1) << ", " << theta_kb(l) << ", " << x_pc(l) << ", " << Q(l) << std::endl;
        l += dl;
    }
}


// ------------------------------------------------------------"on ray" functions ----------------------------------------------------------
Vector3d FixedHeightSolver::vMoment (double l) {
  /*
  Magnetic momentum unit vector
  */
  Vector3d mvec;
  mvec(0) = std::sin(PSR_.chi) * std::cos(phi_ + l / PSR_.RLC);
  mvec(1) = std::sin(PSR_.chi) * std::sin(phi_ + l / PSR_.RLC);
  mvec(2) = std::cos(PSR_.chi);
  return mvec;
}

std::pair<double, double> FixedHeightSolver::find_emission_point()
{ 
  double theta_em = PSR_.chi;
  double phi_em = phi_;
  Vector3d vM = vMoment(0);
  // Vector3d rand_vector = (Vector3d::Random().cross(PSR_.observer_vec)).normalized(); // regularisation of theta_kb
  Vector3d target_vector = PSR_.observer_vec;
  auto func1 = [&](double theta, double phi){
    Vector3d vPoint;
    vPoint << model_.Rem * std::sin(theta) * std::cos(phi), model_.Rem * std::sin(theta) * std::sin(phi), model_.Rem * std::cos(theta);
    return (model_.Bfield(vPoint, vM)).cross(target_vector)(0);
  };
  auto func2 = [&](double theta, double phi){
    Vector3d vPoint;
    vPoint << model_.Rem * std::sin(theta) * std::cos(phi), model_.Rem * std::sin(theta) * std::sin(phi), model_.Rem * std::cos(theta);
    return (model_.Bfield(vPoint, vM)).cross(target_vector)(1);
  };
  for(int i = 0; i < 15; i ++) {
      double f1x = DX(func1, theta_em, phi_em);
      double f2x = DX(func2, theta_em, phi_em);
      double f1y = DY(func1, theta_em, phi_em);
      double f2y = DY(func2, theta_em, phi_em);
      double f1 = func1(theta_em, phi_em);
      double f2 = func2(theta_em, phi_em);
      double dX = (f1y * f2 - f1 * f2y) / (f1x * f2y - f1y * f2x);
      double dY = (f1x * f2 - f1 * f2x) / (f1y * f2x - f2y * f1x);
      theta_em += dX;
      phi_em += dY;
  }
  std::pair<double, double> em_point = {theta_em, phi_em};
  return em_point;
}

Vector3d FixedHeightSolver::vR (double l) {
  /*
  Propagation radius vector
  Here strightforward ray propagation is implemented (but more complex cases of refraction can be considere here too) 
  */
  Vector3d n0(3); // unit vector along the ray
  n0(0) = std::sin(theta_em_) * std::cos(phi_em_); 
  n0(1) = std::sin(theta_em_) * std::sin(phi_em_);
  n0(2) = std::cos(theta_em_);
  return model_.Rem * n0 + l * PSR_.observer_vec;
} 

double FixedHeightSolver::psi_m (double l) {
  return ANGLE(vR(l), vMoment(l));
}

Vector3d FixedHeightSolver::vB (double l) {
  return model_.Bfield(vR(l), vMoment(l));
}

Vector3d FixedHeightSolver::vb (double l) {
  return vB(l).normalized();
}

double FixedHeightSolver::theta_kb (double l) {
  return ANGLE(vB(l), PSR_.observer_vec);
}

Vector3d FixedHeightSolver::vBetaR (double l) {
  return PSR_.Rs / constants::c * PSR_.Omega_vec.cross(vR(l)); 
}

Vector3d FixedHeightSolver::vUdr (double l) {  
  Vector3d vn;
  Vector3d vm;
  if(l < 1){
    vn = (PSR_.observer_vec - PSR_.observer_vec.dot(vb(1)) * vb(1)).normalized();
    vm = (vb(1).cross(vn)).normalized();    
  }
  else{
    vn = (PSR_.observer_vec - PSR_.observer_vec.dot(vb(l)) * vb(l)).normalized();
    vm = (vb(l).cross(vn)).normalized();
  }

  Vector3d temp;
  temp(0) = vBetaR(l).dot(vn);
  temp(1) = vBetaR(l).dot(vm);
  if (temp(0)*temp(0) + temp(1) * temp(1) >= 1.0) {
    double coef = 1/std::sqrt(temp(0)*temp(0) + temp(1) * temp(1))*(1e-3);
    temp(0) *= coef;
    temp(1) *= coef;
    temp(2) = 0;
    // throw_error("ERROR: vUdr > 1.");
    std::cout << "Udr > 1" << std::endl;
    return temp;
  }
  temp(2) = std::sqrt(1 - std::pow(temp(0), 2) - std::pow(temp(1), 2));
  return temp;
}

double FixedHeightSolver::gammaU (double l) {
  double vx = vUdr(l)(0);
  double vy = vUdr(l)(1);
  return std::pow(1 - vx * vx - vy * vy, -0.5);
}

double FixedHeightSolver::x_pc(double l){
  return std::abs(std::sin(psi_m(l))) * std::sqrt(PSR_.RLC / vR(l).norm());
}

double FixedHeightSolver::phi_pc(double l){
  Vector3d m_perp; // Vector, perpendicular to the magnetic axis and e_phi basis vector
  m_perp(0) = -vMoment(l)(2) *  vMoment(l)(0) / std::sqrt(std::pow(vMoment(l)(0), 2) + std::pow(vMoment(l)(1), 2));
  m_perp(1) = -vMoment(l)(2) *  vMoment(l)[1] / std::sqrt(pow(vMoment(l)(0), 2) + std::pow(vMoment(l)(1), 2));
  m_perp(2) = std::sqrt(std::pow(vMoment(l)[0], 2) + std::pow(vMoment(l)[1], 2));
  Vector3d v_perp = vR(l) - vR(l).dot(vMoment(l)) * vMoment(l); // projection of vR, perpendicular to the magntic axis
  if((m_perp.cross(v_perp)).norm() >= 0)
    return constants::PI / 2 + ANGLE(v_perp, m_perp); //REDO ANGLE (?!)
  else
    return constants::PI / 2 - ANGLE(v_perp, m_perp);
}

double FixedHeightSolver::gFunc (double l) {
    return model_.density_profile(x_pc(l), phi_pc(l)); // TODO make 1D/2D differentiation (performace issue)!!!
}

double FixedHeightSolver::Ne(double l) {
  double nGJ = PSR_.Omega_vec.dot(vB(l)) * PSR_.B12 * 1e12 / std::pow(vR(l).norm(), 3) / 
  (2 * constants::PI * constants::c * constants::e);
  return model_.lambda * gFunc(l) * nGJ;
}

double FixedHeightSolver::omegaB (double l) {
  return -constants::e * vB(l).norm() * (PSR_.B12*1e12 / std::pow(vR(l).norm(), 3)) / (constants::me * constants::c);
}

double FixedHeightSolver::omegaW (double l) {
  double vx = vUdr(l)(0);
  double vz = vUdr(l)(2);
  double sinth = std::sin(theta_kb(l));
  double costh = std::cos(theta_kb(l));
  return PSR_.omega_obs * (1 - sinth * vx - costh * vz);
}

double FixedHeightSolver::omegaP (double l) {
  return std::sqrt(4 * constants::PI * constants::e * constants::e * std::abs(Ne(l)) / constants::me);
}

// --------------------------Main Kravtsov-Orlov equation functions------------------------------------

double FixedHeightSolver::delta (double l) {
  double vx = vUdr(l) (0);
  double vy = vUdr(l) (1);
  double sinth = std::sin(theta_kb(l));
  double costh = std::cos(theta_kb(l));
  // std::cout << phi_ * 180 /constants::PI << " " << -vUdr(0) (1) << " " << -vUdr(0) (0) << " " << std::atan2(-std::cos(theta_kb(0))*vUdr(0) (1), std::sin(theta_kb(0))- vUdr(0) (0)) << std::endl;
  // return std::atan2(-costh * vy, sinth-vx);
  // return 0.5 * std::atan(-2*vy*costh*(sinth-vx) / ((sinth-vx)*(sinth-vx) - costh*costh*vy*vy));
  return 0.5 * std::atan2(-2*vy*costh*(sinth-vx), ((sinth-vx)*(sinth-vx) - costh*costh*vy*vy));
}

double FixedHeightSolver::BetaB (double l) {
  Vector3d XX;
  Vector3d YY;
  XX = (PSR_.Omega_vec - PSR_.observer_vec.dot(PSR_.Omega_vec) * PSR_.observer_vec).normalized();
  YY = PSR_.observer_vec.cross(XX);
  if (l < 1e-1){
    l = 1e-1;
  }
  double bx = XX.dot(vB(l));
  double by = YY.dot(vB(l));
  return std::atan2(by, bx);  // ??? was atan (by/bx)
}

double FixedHeightSolver::Q (double l) {
  double vx = vUdr(l)(0);
  double vy = vUdr(l)(1);
  double vz = vUdr(l)(2);
  double sinth = std::sin(theta_kb(l));
  double costh = std::cos(theta_kb(l));
  return  model_.lambda * omegaB(l) * PSR_.omega_obs * (std::pow(sinth - vx, 2) + std::pow(vy * costh, 2)) * model_.Q_type_avrg(std::pow(gammaU(l) * omegaW(l) / omegaB(l), 2))
   / (2 * std::pow(omegaW(l), 2) * (costh * (1 - vx * vx - vy * vy) - vz * (1.0 - sinth * vx))); // ONLY FOR ZERO TEMPERATURE !!! 
}


double FixedHeightSolver::Lambda (double l) {
  double vx = vUdr(l)(0);
  double vy = vUdr(l)(1);
  double sinth = std::sin(theta_kb(l));
  double costh = std::cos(theta_kb(l));
  double avrg = model_.Lambda_type_avrg(std::pow(gammaU(l) * omegaW(l) / omegaB(l), 2));
  return sgn(avrg) * (-1.0 / 2.0) * std::pow(omegaP(l) * gammaU(l) / omegaW(l), 2) * avrg *
   (std::pow(sinth - vx, 2) + std::pow(vy * costh, 2));
}

double FixedHeightSolver::BetaB_derivative(double l)
{
  double dl = 1;
  return (BetaB(l + dl) - BetaB(l)) / dl;
}

double FixedHeightSolver::Q_derivative(double l)
{
  double dl = 1;
  return (Q(l + dl) - Q(l)) / dl;
}

double FixedHeightSolver::Lambda_derivative(double l)
{
  double dl = 1;
  return (Lambda(l + dl) - Lambda(l)) / dl;
}

double FixedHeightSolver::delta_derivative(double l)
{
  double dl = 1;
  return (delta(l + dl) - delta(l)) / dl;
}

// -----------------------Cyclotron Absorption------------------------

double FixedHeightSolver::dtau (double l) {
  double coef = constants::PI * PSR_.Rs / (constants::c * PSR_.omega_obs);
  return coef * std::pow(omegaP(l), 2) * model_.fDist(std::abs(omegaB(l)) / (omegaW(l) * gammaU(l)));
}

double FixedHeightSolver::get_tau(){
  return integrate([this](double l){return this->dtau(l);}, 1, PSR_.RLC);
}
