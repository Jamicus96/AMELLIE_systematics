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
        // Rectangle coordinates (in bin index, recall: bins 0 and N+1 are underflow and overflow)
        unsigned int x_min_bin;
        unsigned int x_max_bin;
        unsigned int t_min_bin;
        unsigned int t_max_bin;

        // Histogram info (min and max are from lower lim of bin 1, to upper lim of bin N)
        double hist_x_min;
        double hist_x_max;
        double hist_t_min;
        double hist_t_max;

        double hist_x_bin_size;
        double hist_t_bin_size;

    public:
        //constructors
        rectangle(TH2F* hist);
        rectangle(TH2F* hist, double x_max, double x_min, double t_max, double t_min);

        // Memeber function
        const unsigned int& X_max_bin() const {return x_max_bin;}
        const unsigned int& X_min_bin() const {return x_min_bin;}
        const unsigned int& T_max_bin() const {return t_max_bin;}
        const unsigned int& T_min_bin() const {return t_min_bin;}

        unsigned int& X_max_bin() {return x_max_bin;}
        unsigned int& X_min_bin() {return x_min_bin;}
        unsigned int& T_max_bin() {return t_max_bin;}
        unsigned int& T_min_bin() {return t_min_bin;}

        void Set_x_max(const double& x_max);
        void Set_x_min(const double& x_min);
        void Set_t_max(const double& t_max);
        void Set_t_min(const double& t_min);
};

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


TH2F* GetHist(std::string root_filename, std::string hist_name = "hPmtResTimeVsCosTheta");

//end header guard
#endif