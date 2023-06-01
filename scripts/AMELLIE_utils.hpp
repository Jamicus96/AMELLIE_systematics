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


TH2F* GetHist(std::string root_filename, std::string hist_name = "hPmtResTimeVsCosTheta");

//end header guard
#endif