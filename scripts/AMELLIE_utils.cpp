#include "AMELLIE_utils.hpp"



/**
 * @brief Construct a new rectangle::rectangle object.
 * 
 * @param X_max 
 * @param X_min 
 * @param Y_max 
 * @param Y_min 
 */
rectangle::rectangle(double X_max, double X_min, double Y_max, double Y_min) {
    if (X_max < X_min or Y_max < Y_min) {
        std::cout << "Wrong limits: max less than min." << std::endl;
        exit(1);
    }
    x_max = X_max; x_min = X_min; y_max = Y_max; y_min = Y_min;
}

const double& rectangle::X_max() const {return x_max;}
const double& rectangle::X_min() const {return x_min;}
const double& rectangle::Y_max() const {return y_max;}
const double& rectangle::Y_min() const {return y_min;}

double& rectangle::X_max() {return x_max;}
double& rectangle::X_min() {return x_min;}
double& rectangle::Y_max() {return y_max;}
double& rectangle::Y_min() {return y_min;}

/**
 * @brief Checks if provided point is inside the rectangular region.
 * 
 * @param point_x 
 * @param point_x 
 * @return true 
 * @return false 
 */
bool rectangle::check_point_inside_rectangle(const double point_x, const double point_y) const {
    if (point_x > x_max or point_x < x_min or point_y > y_max or point_y < y_min) {
        return false;
    } else {
        return true;
    }
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


TH2F* GetHist(std::string root_filename, std::string hist_name) {
    // Read in root file
    TFile* fin = new TFile(root_filename.c_str());
    if (!fin->IsOpen()) {
        std::cout << "Cannot open input file " << root_filename << std::endl;
        exit(1);
    }

    // Get histogram
    TH2F* h = (TH2F*)fin->Get(hist_name.c_str());

    return h;
}