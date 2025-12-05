#include <functional>
#include <vector>
#include <string>

double bilinear_interp(std::vector<double>& xs, std::vector<double>& ys, std::vector<std::vector<double> >& data, double x, double y);


class interpolator2D{
private:
        int Nrs, Nphis;
        std::vector<double> rs;
        std::vector<double> phis;
        std::vector<std::vector<double> > fs;
public:
        interpolator2D ();
        void init(std::string fname);
        void init(std::string fname, double r0, double r1, double phi0, double phi1);
        double get_f(double r, double phi);
};
