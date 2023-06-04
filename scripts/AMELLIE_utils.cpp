#include "AMELLIE_utils.hpp"



rectangle::rectangle(TH2F* hist) {
    // Set histogram binning info
    hist_x_min = hist->GetXaxis()->GetXmin();
    hist_x_max = hist->GetXaxis()->GetXmax();
    hist_t_min = hist->GetYaxis()->GetXmin();
    hist_t_max = hist->GetYaxis()->GetXmax();
    hist_x_bin_size = (hist_x_max - hist_x_min) / (double)(hist->GetNbinsX());
    hist_t_bin_size = (hist_t_max - hist_t_min) / (double)(hist->GetNbinsY());

    // Set default rectangle bin coordinates (max histogram coordinates)
    x_min_bin = 1;
    x_max_bin = round(0.5 + ((hist_x_max - hist_x_min) / hist_x_bin_size));
    t_min_bin = 1;
    t_max_bin = round(0.5 + ((hist_t_max - hist_t_min) / hist_t_bin_size));
}

/**
 * @brief Construct a new rectangle::rectangle object.
 * 
 * @param x_max 
 * @param x_min 
 * @param t_max 
 * @param t_min 
 */
rectangle::rectangle(TH2F* hist, double x_max, double x_min, double t_max, double t_min) {
    if (x_max < x_min or t_max < t_min) {
        std::cout << "Wrong limits: max less than min." << std::endl;
        exit(1);
    }

    // Set histogram binning info
    hist_x_min = hist->GetXaxis()->GetXmin();
    hist_x_max = hist->GetXaxis()->GetXmax();
    hist_t_min = hist->GetYaxis()->GetXmin();
    hist_t_max = hist->GetYaxis()->GetXmax();
    hist_x_bin_size = (hist_x_max - hist_x_min) / (double)(hist->GetNbinsX());
    hist_t_bin_size = (hist_t_max - hist_t_min) / (double)(hist->GetNbinsY());

    // Set rectangle bin coordinates
    this->Set_x_min(x_min);
    this->Set_x_max(x_max);
    this->Set_t_min(t_min);
    this->Set_t_max(t_max);
}

void rectangle::Set_x_max(const double& x_max) {
    if (hist_x_max <= x_max) {
        x_max_bin = round(0.5 + ((hist_x_max - hist_x_min) / hist_x_bin_size));
    } else {
        x_max_bin = round(0.5 + ((x_max - hist_x_min) / hist_x_bin_size));
    }
}
void rectangle::Set_x_min(const double& x_min) {
    if (hist_x_min <= x_min) {
        x_min_bin = 1;
    } else {
        x_min_bin = round(0.5 + ((x_min - hist_x_min) / hist_x_bin_size));
    }
}
void rectangle::Set_t_max(const double& t_max) {
    if (hist_t_max <= t_max) {
        t_max_bin = round(0.5 + ((hist_t_max - hist_t_min) / hist_t_bin_size));
    } else {
        t_max_bin = round(0.5 + ((t_max - hist_t_min) / hist_t_bin_size));
    }
}
void rectangle::Set_t_min(const double& t_min) {
    if (hist_t_min <= t_min) {
        t_min_bin = 1;
    } else {
        t_min_bin = round(0.5 + ((t_min - hist_t_min) / hist_t_bin_size));
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