#include "interpolator.h"
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <functional> 
#include <fstream>

double bilinear_interp(std::vector<double>& xs, std::vector<double>& ys, std::vector<std::vector<double> >& data, double x, double y){
	double x1, y1, x2, y2;
	double x_start, x_end, x_step;
	double y_start, y_end, y_step;
	x_start = xs[0];
	x_end = xs[xs.size() - 1];
	x_step = xs[1] - xs[0];
	y_start = ys[0];
	y_end = ys[ys.size() - 1];
	y_step = ys[1] - ys[0];
	//std::cout << x_start << " " << x_end << " " << x_step << " " << y_start << " " << y_end << " " << y_step << std::endl;
	int x_ind1, y_ind1, x_ind2, y_ind2;
	x_ind1 = (x - x_start) / x_step;
	x_ind2 = x_ind1 + 1;
	y_ind1 = (y - y_start) / y_step;
	y_ind2 = y_ind1 + 1;
	if(x_ind2 >= xs.size() || x_ind1 < 0 || y_ind2 >= ys.size() || y_ind1 < 0){
		return 0.01;
	}
	x1 = xs[x_ind1];
	x2 = xs[x_ind2];
	y1 = ys[y_ind1];
	y2 = ys[y_ind2];
	double w11, w12, w21, w22;
	double denominator = x_step * y_step;
	w11 = (x2 - x) * (y2 - y) / denominator;
	w12 = (x2 - x) * (y - y1) / denominator;
	w21 = (x - x1) * (y2 - y) / denominator;
	w22 = (x - x1) * (y - y1) / denominator;
	return w11 * data[x_ind1][y_ind1] + w12 * data[x_ind1][y_ind2] + w21 * data[x_ind2][y_ind1] + w22 * data[x_ind2][y_ind2];
}


interpolator2D::interpolator2D(){
	Nrs = 0;
	Nphis = 0;
}

void interpolator2D::init(std::string fname){
	std::ifstream in (fname);
	in >> Nrs >> Nphis;
	rs.resize(Nrs);
	phis.resize(Nphis);
	fs.resize(Nrs, std::vector<double>(Nphis));
	for(int i=0; i<Nrs; i++){
		in >> rs[i];
	}
	for(int i=0; i<Nphis; i++){
		in >> phis[i];
	}
	for(int i=0; i<Nrs; i++){
		for(int j=0; j<Nphis; j++){
			in>>fs[i][j];
		}
	}
	in.close();
}

void interpolator2D::init(std::string fname, double r0, double r1, double phi0, double phi1){
	std::ifstream in (fname);
	in >> Nrs >> Nphis;
	rs.resize(Nrs);
	phis.resize(Nphis);
	fs.resize(Nrs, std::vector<double>(Nphis));
	for(int i=0; i<Nrs; i++){
		rs[i] = r0 + (r1 - r0) / Nrs * i;
	}
	for(int i=0; i<Nphis; i++){
		phis[i] = phi0 + (phi1 - phi0) / Nphis * i;
	}
	for(int i=0; i<Nrs; i++){
		for(int j=0; j<Nphis; j++){
			in>>fs[i][j];
		}
	}
	in.close();
}

double interpolator2D::get_f(double r, double phi){
	return bilinear_interp(rs, phis, fs, r, phi);
}