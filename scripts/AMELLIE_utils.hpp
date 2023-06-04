// header guard:
#ifndef AMELLIE_utils
#define AMELLIE_utils

// include
#include <iostream>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>


/**
 * @brief Class to encode the direct and reflected regions
 * IDEA: make the coords in bin number + save hist binning info so that it can take normal coords
 * as arguments and compute bin coords. That way bin coords are only computed once for all hists.
 * 
 */
class rectangle {
    private:
        double x_max;
        double x_min;
        double t_max;
        double t_min;

    public:
        //constructor
        rectangle(double X_maxi, double X_mini, double T_maxi, double T_mini);

        // Memeber function
        const double& X_max() const {return x_max;}
        const double& X_min() const {return x_min;}
        const double& T_max() const {return t_max;}
        const double& T_min() const {return t_min;}

        double& X_max() {return x_max;}
        double& X_min() {return x_min;}
        double& T_max() {return t_max;}
        double& T_min() {return t_min;}

        bool check_point_inside_rectangle(const double point_x, const double point_t) const;
};

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

unsigned int Xbin_max(const rectangle& rect, const double Xbin_size, const double Xbin_max, const double Xbin_min);
unsigned int Xbin_min(const rectangle& rect, const double Xbin_size, const double Xbin_max, const double Xbin_min);
unsigned int Tbin_max(const rectangle& rect, const double Tbin_size, const double Tbin_max, const double Tbin_min);
unsigned int Tbin_min(const rectangle& rect, const double Tbin_size, const double Tbin_max, const double Tbin_min);


TH2F* GetHist(std::string root_filename, std::string hist_name = "hPmtResTimeVsCosTheta");

//end header guard
#endif