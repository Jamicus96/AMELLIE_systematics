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
        double y_max;
        double y_min;

    public:
        //constructor
        rectangle(double X_maxi, double X_mini, double Y_maxi, double Y_mini);

        // Memeber function
        double& X_max();
        double& X_min();
        double& Y_max();
        double& Y_min();

        bool check_point_inside_rectangle(const double point_x, const double point_y);
};

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


TH2F* GetHist(std::string root_filename, std::string hist_name = "hPmtResTimeVsCosTheta");

//end header guard
#endif