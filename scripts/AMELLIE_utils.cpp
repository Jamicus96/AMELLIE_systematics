#include "AMELLIE_utils.hpp"



/**
 * @brief Construct a new rectangle::rectangle object.
 * 
 * @param X_max 
 * @param X_min 
 * @param T_max 
 * @param T_min 
 */
rectangle::rectangle(double X_max, double X_min, double T_max, double T_min) {
    if (X_max < X_min or T_max < T_min) {
        std::cout << "Wrong limits: max less than min." << std::endl;
        exit(1);
    }
    x_max = X_max; x_min = X_min; t_max = T_max; t_min = T_min;
}

unsigned int Xbin_max(const rectangle& rect, const double Xbin_size, const double Xbin_max, const double Xbin_min) {
    if (rect.X_max() >= Xbin_max) {
        return (unsigned int) ((Xbin_max - Xbin_min) / Xbin_size);
    } else {
        return (unsigned int) ((rect.X_max() - Xbin_min) / Xbin_size);
    }
}
unsigned int Xbin_min(const rectangle& rect, const double Xbin_size, const double Xbin_max, const double Xbin_min) {
    if (rect.X_min() <= Xbin_min) {
        return 1;
    } else {
        return (unsigned int) ((rect.X_min() - Xbin_min) / Xbin_size);
    }
}
unsigned int Tbin_max(const rectangle& rect, const double Tbin_size, const double Tbin_max, const double Tbin_min) {
    if (rect.T_max() >= Tbin_max) {
        return (unsigned int) ((Tbin_max - Tbin_min) / Tbin_size);
    } else {
        return (unsigned int) ((rect.T_max() - Tbin_min) / Tbin_size);
    }
}
unsigned int Tbin_min(const rectangle& rect, const double Tbin_size, const double Tbin_max, const double Tbin_min) {
    if (rect.T_min() <= Tbin_min) {
        return 1;
    } else {
        return (unsigned int) ((rect.T_min() - Tbin_min) / Tbin_size);
    }
}

/**
 * @brief Checks if provided point is inside the rectangular region.
 * 
 * @param point_x 
 * @param point_t 
 * @return true 
 * @return false 
 */
bool rectangle::check_point_inside_rectangle(const double point_x, const double point_t) const {
    if (point_x > x_max or point_x < x_min or point_t > t_max or point_t < t_min) {
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