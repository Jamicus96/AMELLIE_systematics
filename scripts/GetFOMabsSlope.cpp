#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLine.h>
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <vector>
#include "AMELLIE_utils.hpp"


std::vector<double> create_lims(const double min, const double max, const unsigned int N);
double GetFOMabsSlope(const rectangle& direct_region, const rectangle& reflected_region, const std::vector<TH2F*>& hists, const std::vector<double>& abs_scalings, const std::vector<double>& bin_info, const int abs1_idx);
std::vector<double> getRatio(const rectangle& direct_region, const rectangle& reflected_region, TH2F* allPathsHist, const std::vector<double>& bin_info);
std::vector<std::vector<std::string> > readInfoFile(std::string tracking_hist_repo, std::string info_file, const unsigned int first_line, const unsigned int last_line, const bool verbose);



int main(int argc, char** argv) {
    // Read in text file with list of stats that were simulated (different abs, the rest the same for now)
    std::string info_file = argv[1];
    unsigned int first_line = std::stoi(argv[2]);
    unsigned int last_line = std::stoi(argv[3]);
    std::string tracking_hist_repo = argv[4];
    // Read in output filename (txt)
    std::string output_file = argv[5];
    // record extra info?
    bool verbose = std::stoi(argv[6]);

    // Create list of tracking hist file names from info file, and list abs_scalings
    std::vector<std::vector<std::string> > file_info = readInfoFile(tracking_hist_repo, info_file, first_line, last_line, verbose);

    // Unpack info
    unsigned int abs1_idx = std::stoi(file_info.at(2).at(0));
    std::vector<double> abs_scalings;
    std::vector<TH2F*> hists;
    for (unsigned int i = 0; i < file_info.at(0).size(); ++i) {
        abs_scalings.push_back(std::stod(file_info.at(1).at(i)));
        hists.push_back(GetHist(file_info.at(0).at(i), "hPmtResTimeVsCosTheta"));
    }

    // Create reagion lims to iterate over (6 values):
    double dir_x_max_min = -0.99; double dir_x_max_max = -0.7;
    double ref_x_min_min = 0.7; double ref_x_min_max = 0.99;
    double dir_t_cen_min = -10.0; double dir_t_cen_max = 10.0;
    double ref_t_cen_min = -10.0; double ref_t_cen_max = 10.0;
    double dir_t_wid_min = 1.0; double dir_t_wid_max = 30.0;
    double ref_t_wid_min = 1.0; double ref_t_wid_max = 30.0;

    unsigned int N = 10;
    std::vector<double> dir_x_max_lims = create_lims(dir_x_max_min, dir_x_max_max, N);
    std::vector<double> ref_x_min_lims = create_lims(ref_x_min_min, ref_x_min_max, N);
    std::vector<double> dir_t_cen_lims = create_lims(dir_t_cen_min, dir_t_cen_max, N);
    std::vector<double> ref_t_cen_lims = create_lims(ref_t_cen_min, ref_t_cen_max, N);
    std::vector<double> dir_t_wid_lims = create_lims(dir_t_wid_min, dir_t_wid_max, N);
    std::vector<double> ref_t_wid_lims = create_lims(ref_t_wid_min, ref_t_wid_max, N);

    // Set up regions
    rectangle direct_region = rectangle(dir_x_max_lims.at(0), -1.0, dir_t_cen_lims.at(0) + 0.5 * dir_t_wid_lims.at(0), dir_t_cen_lims.at(0) - 0.5 * dir_t_wid_lims.at(0));
    rectangle reflected_region = rectangle(1.0, ref_x_min_lims.at(0), ref_t_cen_lims.at(0) + 0.5 * ref_t_wid_lims.at(0), ref_t_cen_lims.at(0) - 0.5 * ref_t_wid_lims.at(0));

    // Get hist binning info (Assuming they are all the same, which they should be)
    double X_min = hists.at(0)->GetXaxis()->GetXmin();
    double X_max = hists.at(0)->GetXaxis()->GetXmax();
    double T_min = hists.at(0)->GetYaxis()->GetXmin();
    double T_max = hists.at(0)->GetYaxis()->GetXmax();

    double Xbin_size = (X_max - X_min) / (double)(hists.at(0)->GetNbinsX());
    double Tbin_size = (T_max - T_min) / (double)(hists.at(0)->GetNbinsY());
    std::vector<double> bin_info = {X_min, X_max, Xbin_size, T_min, T_max, Tbin_size}; // Package up to pass to functions

    // Loop through all region lims, and compute FOM, and print them to file
    std::ofstream datafile;
    datafile.open(output_file.c_str(), std::ios::trunc);
    double slope;
    for (unsigned int i = 0; i < N; ++i) {
        direct_region.X_max() = dir_x_max_lims.at(i);
        for (unsigned int j = 0; j < N; ++j) {
            for (unsigned int k = 0; k < N; ++k) {
                direct_region.T_max() = dir_t_cen_lims.at(j) + 0.5 * dir_t_wid_lims.at(k);
                direct_region.T_min() = dir_t_cen_lims.at(j) - 0.5 * dir_t_wid_lims.at(k);
                for (unsigned int l = 0; l < N; ++l) {
                    reflected_region.X_min() = ref_x_min_lims.at(l);
                    for (unsigned int m = 0; m < N; ++m) {
                        for (unsigned int n = 0; n < N; ++n) {
                            direct_region.T_max() = ref_t_cen_lims.at(m) + 0.5 * ref_t_wid_lims.at(n);
                            direct_region.T_min() = ref_t_cen_lims.at(m) - 0.5 * ref_t_wid_lims.at(n);

                            // Compute slope
                            slope = GetFOMabsSlope(direct_region, reflected_region, hists, abs_scalings, bin_info, abs1_idx);

                            // print region limits first, and then slope. All space-separated
                            datafile << dir_x_max_lims.at(i) << " " << dir_t_cen_lims.at(j) << " " << dir_t_wid_lims.at(k) << " " << ref_x_min_lims.at(l)
                                     << " " << ref_t_cen_lims.at(m) << " " << ref_t_wid_lims.at(n) << " " << slope << std::endl;
                        }
                    }
                }
            }
        }
    }

    return 0;
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

std::vector<double> create_lims(const double min, const double max, const unsigned int N) {
    const double step = (max - min) / (double)(N - 1);

    std::vector<double> lims = {min};
    for (unsigned int n = 0; n < N; ++n) {
        lims.push_back(min + n * step);
    }

    return lims;
}


/**
 * @brief Returns FOM/abs slope (called y)
 * 
 * @param direct_region 
 * @param reflected_region 
 * @param hists 
 * @param abs_scalings 
 * @param bin_info 
 * @param abs1_idx 
 * @return double 
 */
double GetFOMabsSlope(const rectangle& direct_region, const rectangle& reflected_region, const std::vector<TH2F*>& hists, const std::vector<double>& abs_scalings, const std::vector<double>& bin_info, const int abs1_idx) {

    // Get ratio of hits in both triangular regions, and associated error, then normalise all to the ratio
    // found at abs = 1. Then compute best fit slope of normalised ratio vs abs:

    // abs = 1 ratio
    std::vector<double> ratio_abs1 = getRatio(direct_region, reflected_region, hists.at(abs1_idx), bin_info);

    double ratio;
    double ratio_err;
    double S_xy;
    double S_xx;
    std::vector<double> ratio_res = {0.0, 0.0};
    for (unsigned int n = 0; n < abs_scalings.size(); ++n) {
        if (n == abs1_idx) continue;  // values are just zero

        // ratio for other abs
        ratio_res = getRatio(direct_region, reflected_region, hists.at(n), bin_info);
        if (ratio_res.at(0) == -1.0) {
            std::cout << "ERROR: Region with no hits :/" << std::endl;
            exit(1);
        }
        ratio = ratio_res.at(0) / ratio_abs1.at(0); // normalise ratio to value at abs=1.0
        ratio_err = sqrt(((ratio_res.at(1)*ratio_res.at(1)) / (ratio_abs1.at(0)*ratio_abs1.at(0)))
                        + ((ratio_res.at(0)*ratio_res.at(0) * ratio_abs1.at(1)*ratio_abs1.at(1))
                        / (ratio_abs1.at(0)*ratio_abs1.at(0)*ratio_abs1.at(0)*ratio_abs1.at(0)))); // update errors accordingly

        // Add to quantities used to find slope (FOM) of ratio vs abs
        S_xy += ((abs_scalings.at(n) - 1.0) * (ratio - 1.0)) / (ratio_err*ratio_err);
        S_xx += ((abs_scalings.at(n) - 1.0) * (abs_scalings.at(n) - 1.0)) / (ratio_err*ratio_err);
    }

    return S_xy / S_xx;  // return best fit slope (FOM)
}


/**
 * @brief Get the Ratio [#hits in reflected region / #hits in direct region], along with the associated statistical error.
 * 
 * @param direct_region 
 * @param reflected_region 
 * @param hist 
 * @param bin_info = {X_min, X_max, Xbin_size, T_min, T_max, Tbin_size} Histogram binning info
 * @return std::vector<double> 
 */
std::vector<double> getRatio(const rectangle& direct_region, const rectangle& reflected_region, TH2F* hist, const std::vector<double>& bin_info) {

    // Count number of hits in direct and reflected regions
    double direct_count = 0.0;
    double reflected_count = 0.0;
    double xBinCenter;
    double tBinCenter;

    // Loop over bins in direct region
    for (unsigned int x = Xbin_min(direct_region, bin_info.at(2), bin_info.at(1), bin_info.at(0)); x < Xbin_max(direct_region, bin_info.at(2), bin_info.at(1), bin_info.at(0)); ++x) {  // loop over histogram bins
        xBinCenter = hist->GetXaxis()->GetBinCenter(x);
        for (unsigned int t = Tbin_min(direct_region, bin_info.at(5), bin_info.at(4), bin_info.at(3)); t < Tbin_max(direct_region, bin_info.at(5), bin_info.at(4), bin_info.at(3)); ++t) {
            direct_count += hist->GetBinContent(x, hist->GetYaxis()->GetBinCenter(t));
        }
    }

    // Loop over bins in reflected region
    for (unsigned int x = Xbin_min(reflected_region, bin_info.at(2), bin_info.at(1), bin_info.at(0)); x < Xbin_max(reflected_region, bin_info.at(2), bin_info.at(1), bin_info.at(0)); ++x) {  // loop over histogram bins
        xBinCenter = hist->GetXaxis()->GetBinCenter(x);
        for (unsigned int t = Tbin_min(reflected_region, bin_info.at(5), bin_info.at(4), bin_info.at(3)); t < Tbin_max(reflected_region, bin_info.at(5), bin_info.at(4), bin_info.at(3)); ++t) {
            reflected_count += hist->GetBinContent(x, hist->GetYaxis()->GetBinCenter(t));
        }
    }

    // Compute ratio of both, and associated stat error
    std::vector<double> ratio;
    if (direct_count == 0 || reflected_count == 0) {
        ratio.push_back(-1.0);
        ratio.push_back(-1.0);
    } else {
        ratio.push_back(reflected_count / direct_count);  // ratio value
        ratio.push_back((sqrt(reflected_count) / direct_count) * sqrt(1.0 + (reflected_count / direct_count)));  // ratio error (assuming stat errors for regions of sqrt(N))
    }
    return ratio;
}


/**
 * @brief Reads information from text file with simulation info into.
 * 
 * @param tracking_hist_repo folder containing root file with histograms to analyse
 * @param info_file File with analysis info this will read
 * @param first_line First line to read info from
 * @param last_line Last line to read info from
 * @param verbose
 * @return std::vector<std::string> {traking_files, abs_scalings, abs1_idx}
 */
std::vector<std::vector<std::string> > readInfoFile(std::string tracking_hist_repo, std::string info_file, const unsigned int first_line, const unsigned int last_line, const bool verbose) {

    std::ifstream file(info_file);
    std::string str;
    std::string geo_file;  // sans the ".geo" part
    std::string inner_av_material;
    std::string wavelength;  // nm
    std::string fibre;
    std::string reemission;  // reemission fraction of absorbed light
    std::vector<std::string> abs_scalings;  // absorption coefficient normalised so that 1.0 = current value
    std::vector<std::string> traking_files;
    std::vector<std::string> abs1_idx = {"-1"};  // Index in abs_scalings where abs=1.0

    unsigned int i = 0;
    std::string delimiter = ", ";
    std::size_t pos; std::size_t pos_1;
    std::string substr;
    std::string full_geo;
    while (std::getline(file, str)) {
        // format in info file: geo_file.geo, LEDnum, fibre, reemis, abs
        // format in tracking hist file name: tot_AMELLIE_geoFile_LEDnum_fibre_reemis_abs.root

        if (i < first_line || i > last_line) {continue;} // only read between selected lines
        
        pos = str.find(delimiter);
        full_geo = str.substr(0, pos);
        pos_1 = full_geo.find(".geo");
        geo_file = full_geo.substr(0, pos_1);
        str = str.substr(pos + delimiter.size());

        pos = str.find(delimiter);
        substr = str.substr(0, pos);
        inner_av_material = str.substr(0, pos);
        str = str.substr(pos + delimiter.size());

        pos = str.find(delimiter);
        substr = str.substr(0, pos);
        wavelength = str.substr(0, pos);
        str = str.substr(pos + delimiter.size());

        pos = str.find(delimiter);
        substr = str.substr(0, pos);
        fibre = str.substr(0, pos);
        str = str.substr(pos + delimiter.size());

        pos = str.find(delimiter);
        substr = str.substr(0, pos);
        reemission = str.substr(0, pos);
        str = str.substr(pos + delimiter.size());

        abs_scalings.push_back(str);
        traking_files.push_back(tracking_hist_repo + "tot_hists_AMELLIE_" + geo_file + "_" + inner_av_material + "_" + wavelength + "_" + fibre + "_reemis" + reemission + "_abs" + str + ".root");

        if (std::stof(str) == 1.0) {abs1_idx.at(0) = std::to_string(i);}  // found absorption=1 case
        if (verbose) {
            std::cout << "abs = " << abs_scalings.at(i) << ", tracking file = " << traking_files.at(i) << ", abs1_idx = " << abs1_idx.at(0) << std::endl;
        }
        ++i;
    }

    if (abs1_idx.at(0) == "-1") {
        std::cout << "ERROR: Absorption 1.0 case not found." << std::endl;
        exit(1);
    }

    std::vector<std::vector<std::string> > results = {traking_files, abs_scalings, abs1_idx};
    return results;
}
