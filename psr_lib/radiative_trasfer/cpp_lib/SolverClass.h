#include "PSRClass.h"
#include "FixedHeightModelClass.h"
#include <Eigen/Dense>


#pragma once

using Eigen::Vector3d;

class FixedHeightSolver{
    private:
    double phi_;
    double mode_;
    ObservedRadioPulsar PSR_;
    double theta_em_;
    double phi_em_;
    FixedHeightModel model_;

    public:
    FixedHeightSolver(double phi_deg, int mode, ObservedRadioPulsar &PSR, FixedHeightModel &model);

    /// @brief solve Kravtsov-Orlov equations 
    /// @param theta_inital inital values for comlex polarisation angle (as length 2 vector)
    /// @param l1 - initial point
    /// @param l2 - final point
    /// @param mode emission mode (1 - Omode, 0 - Xmode)
    /// @return theta_final -  final complex polarisation angle (as length 2 vector)
    std::vector<double> solve_KO_equations(std::vector<double> theta_initial, double l1, double l2, std::string log_path="");

    /// @brief RHS for solving Kravtsov-Orlov equations using <boost/numeric/odeint.hpp> integration routines
    /// @param f RHS function
    /// @param dydx derivative by reference
    /// @param l independent variable by reference
    void RHS_for_boost(const std::vector<double>& f, std::vector<double> &dydx, double l);

    /// @brief find approximate solution for Kravtsov-Orlov equations in dense plasma regions 
    /// @param theta_inital inital values for comlex polarisation angle (as length 2 vector)
    /// @param mode emission mode (1 - Omode, 0 - Xmode)
    /// @return theta_final -  final complex polarisation angle (as length 2 vector)
    std::vector<double> find_approximate_KO_solution(std::vector<double> theta_initial);

    /// @brief find approximate solution for Kravtsov-Orlov equations in dense plasma regions 
    /// @param l - distance from the star. Theta_inintal is determined by emission mode
    /// @param mode emission mode (1 - Omode, 0 - Xmode)
    /// @return theta_final -  final complex polarisation angle (as length 2 vector)
    std::vector<double> find_approximate_KO_solution(double l);

    /// @brief This function finds a distance from emission point where oscillations fade out but p.a. is still strictly
    /// following beta + delta. 
    /// This point is determined from the condition |Lambda_derivative / Lambda^2 * 2 * c / omega| ~ 1 
    /// @param use_binary_search - flag, which determines wither to use binary search or not. Binary search can fail when 
    /// Lambda is not monotonous
    /// @return point, where integration of Kravtsov-Orlov equations can be started
    double find_initial_point(bool use_binary_search);

    /// @brief find emission intensity
    /// @return emission intensity in arbitrary uits
    double find_intensity();

    /// @return tau - optical thickness (I = I0 * e^(-tau))
    double get_tau();

    void write_params_on_ray(std::string log_path);



    // --------------------------Main Kravtsov-Orlov equation functions------------------------------------

    /// @brief additional phase parameter related to ExB drift 
    /// @param l distance alnog the ray 
    /// @return delta phase 
    double delta (double l);

    /// @brief polar angle of magnetic filed projected to the plane perpendicular to the ray
    /// @param l distnace along the ray
    /// @return BetaB angle in radians
    double BetaB (double l);

    /// @param l distance alnog the ray
    /// @return Q parameter from Kravtsov-Orlov equations
    double Q (double l);

    /// @param l distance along the ray
    /// @return Lambda parameter form Kravtsov-Orlov equations
    double Lambda (double l);
    // -----------------------------------------------------------------"On ray" functions-------------------------------------------------
    private:

    /// @param l distance along the ray
    /// @return magnetic moment unit vector (Egien 3d Vector)
    Vector3d vMoment (double l);

    /// @brief finds spherical angular coordinates of the emission point (on the distance Rem) using ...
    /// @return pair theta_em, phi_em
    std::pair<double, double> find_emission_point();

    /// @param l distance along the ray
    /// @return vector to the point on the ray as Eigen 3d Vector
    Vector3d vR (double l);

    /// @param l distance along the ray
    /// @return angle between magnetic axis and r vector
    double psi_m (double l);

    /// @param l distance along the ray
    /// @return Magnetic field vector as Eigen 3d Vector
    Vector3d vB (double l);

    /// @param l distance along the ray
    /// @return normalized magnetic field vector
    Vector3d vb (double l);

    /// @param l distance along the ray
    /// @return angle between wave vecotr and magnetic field (in radians)
    double theta_kb (double l);

    /// @param l distance along the ray
    /// @return Omega_vec corss r_vec 
    Vector3d vBetaR (double l);

    /// @brief ExB drift particle velocity component (V = V_|| * b_vec + Udr)
    /// @param l distance along the ray
    /// @return drift velocity Eigen 3d Vector  
    Vector3d vUdr (double l);

    /// @param l distance along the ray
    /// @return drift velocity gamma factor
    double gammaU (double l);

    /// @param l distance along the ray
    /// @return distance from the field line footpoint to the polar cap center
    /// @note Only for dipolar magnetic field!!!
    double x_pc(double l);

    /// @param l distance along the ray
    /// @return polar cap angle cooridnate of the magnetic field line footpoint.
    /// Angle is counted from the East-West line on the polar cap
    /// @note Applicable only for magnetic fields with zero torsion
    double phi_pc(double l);

    /// @brief Plasma density transvers profile
    /// @param l distance along the ray
    /// @return normalized plasma density
    double gFunc (double l);

    /// @param l  distance along the ray
    /// @return Plasma density in physical units (g/cm^3)
    double Ne(double l);

    /// @param l distance along the ray
    /// @return Local cyclotron frequency (s^-1)
    double omegaB (double l);

    /// @param l distance along the ray
    /// @return wave frequency in plasma rest frame
    double omegaW (double l);

    /// @param l distance along the ray
    /// @return Local plasma frequency (s^-1)
    double omegaP (double l);

    /// @brief ...
    /// @param l 
    /// @return ...
    bool stop_condition(double l);


    double Q_denominator(double l);



    /// @brief polar angle of magnetic filed projected to the plane perpendicular to the ray
    /// @param l distnace along the ray
    /// @return Derivative of BetaB angle in radians along the ray
    double BetaB_derivative (double l);

    /// @param l distance alnog the ray
    /// @return Derivative of Q parameter from Kravtsov-Orlov equations along the ray
    double Q_derivative (double l);

    /// @param l distance along the ray
    /// @return Derivative of Lambda parameter form Kravtsov-Orlov equations along the ray
    double Lambda_derivative (double l);

    /// @param l distance along the ray
    /// @return Derivative of delta parameter form Kravtsov-Orlov equations along the ray
    double delta_derivative(double l);

    // -----------------------Cyclotron Absorption------------------------

    /// @param l - distance along the ray
    /// @return dtau - differential optical thickness (I = I0 * e^(-tau), tau = integral dtau)
    double dtau (double l);

};

