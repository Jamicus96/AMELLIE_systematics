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
void GetFOMabsSlope(std::vector<double>& slope_info, const rectangle& direct_region, const rectangle& reflected_region, const std::vector<TH2F*>& hists, const std::vector<double>& abs_scalings_min1, const bool verbose);
void getRatio(std::vector<double>& ratio_res, const rectangle& direct_region, const rectangle& reflected_region, TH2F* allPathsHist, const bool verbose);
// std::vector<std::vector<std::string>> readInfoFile(std::string tracking_hist_repo, std::string info_file, const unsigned int first_line, const unsigned int last_line, const bool verbose);



int main(int argc, char** argv) {
    // Read in text file with list of stats that were simulated (different abs, the rest the same for now)
    std::string tracking_hist_repo = argv[1];
    std::string output_file = argv[2];
    bool verbose = std::stoi(argv[3]);

    // Read in region limit ranges
    double dir_x_max_min = std::stod(argv[4]);
    double dir_x_max_max = std::stod(argv[5]);
    double ref_x_min_min = std::stod(argv[6]);
    double ref_x_min_max = std::stod(argv[7]);
    double dir_t_cen_min = std::stod(argv[8]);
    double dir_t_cen_max = std::stod(argv[9]);
    double ref_t_cen_min = std::stod(argv[10]);
    double ref_t_cen_max = std::stod(argv[11]);
    double dir_t_wid_min = std::stod(argv[12]);
    double dir_t_wid_max = std::stod(argv[13]);
    double ref_t_wid_min = std::stod(argv[14]);
    double ref_t_wid_max = std::stod(argv[15]);
    unsigned int N = std::stoi(argv[16]);

    // Read in hist files and corresponsing abs
    std::vector<double> abs_scalings_min1;
    std::vector<std::string> traking_files;
    bool write_abs = true;
    for (unsigned int i = 17; i < argc; ++i) {
        if (write_abs) {
            abs_scalings_min1.push_back(std::stod(argv[i]) - 1.0);
            write_abs = false;
        } else {
            traking_files.push_back(argv[i]);
            write_abs = true;
        }
    }

    // Unpack info
    std::vector<TH2F*> hists;
    for (unsigned int i = 0; i < traking_files.size(); ++i) {
        hists.push_back(GetHist(traking_files.at(i), "hPmtResTimeVsCosTheta"));
    }

    // Create reagion lims to iterate over (6 values):
    std::vector<double> dir_x_max_lims = create_lims(dir_x_max_min, dir_x_max_max, N);
    std::vector<double> ref_x_min_lims = create_lims(ref_x_min_min, ref_x_min_max, N);
    std::vector<double> dir_t_cen_lims = create_lims(dir_t_cen_min, dir_t_cen_max, N);
    std::vector<double> ref_t_cen_lims = create_lims(ref_t_cen_min, ref_t_cen_max, N);
    std::vector<double> dir_t_wid_lims = create_lims(dir_t_wid_min, dir_t_wid_max, N);
    std::vector<double> ref_t_wid_lims = create_lims(ref_t_wid_min, ref_t_wid_max, N);

    // Set up regions
    rectangle direct_region = rectangle(hists.at(0));
    direct_region.Set_x_min(-1.0);  // harcoded at histgram left edge
    rectangle reflected_region = rectangle(hists.at(0));
    direct_region.Set_x_max(1.0);  // harcoded at histgram right edge

    // Loop through all region lims, and compute FOM, and print them to file
    std::ofstream datafile;
    datafile.open(output_file.c_str(), std::ios::trunc);
    std::vector<double> slope_info = {0.0, 0.0};
    for (unsigned int i = 0; i < N; ++i) {
        direct_region.Set_x_max(dir_x_max_lims.at(i));
        for (unsigned int j = 0; j < N; ++j) {
            for (unsigned int k = 0; k < N; ++k) {
                direct_region.Set_t_max(dir_t_cen_lims.at(j) + 0.5 * dir_t_wid_lims.at(k));
                direct_region.Set_t_min(dir_t_cen_lims.at(j) - 0.5 * dir_t_wid_lims.at(k));
                if (verbose) std::cout << "direct region: x bins € [" << direct_region.X_min_bin() << ", " << direct_region.X_max_bin() << "], t bins € [" << direct_region.T_min_bin() << ", " << direct_region.T_max_bin() << "]" << std::endl;

                for (unsigned int l = 0; l < N; ++l) {
                    reflected_region.Set_x_min(ref_x_min_lims.at(l));
                    for (unsigned int m = 0; m < N; ++m) {
                        for (unsigned int n = 0; n < N; ++n) {
                            reflected_region.Set_t_max(ref_t_cen_lims.at(m) + 0.5 * ref_t_wid_lims.at(n));
                            reflected_region.Set_t_min(ref_t_cen_lims.at(m) - 0.5 * ref_t_wid_lims.at(n));
                            if (verbose) std::cout << "reflected region: x bins € [" << reflected_region.X_min_bin() << ", " << reflected_region.X_max_bin() << "], t bins € [" << reflected_region.T_min_bin() << ", " << reflected_region.T_max_bin() << "]" << std::endl;

                            // Compute slope
                            GetFOMabsSlope(slope_info, direct_region, reflected_region, hists, abs_scalings_min1, verbose);
                            if (verbose) std::cout << "Slope = " << slope_info.at(0) << std::endl;

                            // print region limits first, and then slope. All space-separated
                            datafile << dir_x_max_lims.at(i) << " " << dir_t_cen_lims.at(j) << " " << dir_t_wid_lims.at(k) << " " << ref_x_min_lims.at(l)
                                     << " " << ref_t_cen_lims.at(m) << " " << ref_t_wid_lims.at(n) << " " << slope_info.at(0) << " " << slope_info.at(1) << std::endl;
                        }
                    }
                }

            }
        }
    }

    return 0;
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/**
 * @brief Create a list of region limits to use between min and max
 * 
 * @param min 
 * @param max 
 * @param N 
 * @return std::vector<double> 
 */
std::vector<double> create_lims(const double min, const double max, const unsigned int N) {
    const double step = (max - min) / (double)(N - 1);

    std::vector<double> lims = {min};
    for (unsigned int n = 1; n < N; ++n) {
        lims.push_back(min + n * step);
    }

    return lims;
}


/**
 * @brief Returns FOM/abs slope (called y), and reduced Chi squared
 * 
 * @param slope_info = {fit slope (FOM), reduced Chi squared}, what gets modified to return value
 * @param direct_region 
 * @param reflected_region 
 * @param hists 
 * @param abs_scalings_min1 
 * @param bin_info 
 * @return std::vector<double> 
 */
void GetFOMabsSlope(std::vector<double>& slope_info, const rectangle& direct_region, const rectangle& reflected_region, const std::vector<TH2F*>& hists, const std::vector<double>& abs_scalings_min1, const bool verbose) {

    // Get ratio of hits in both triangular regions, and associated error, then normalise all to the ratio
    // found at abs = 1. Then compute best fit slope of normalised ratio vs abs:

    // abs = 1 ratio
    // if (verbose) std::cout << "Getting ratio for idx = " << 0 << " (abs = 1)" << std::endl;
    std::vector<double> ratio_abs1 = {0.0, 0.0};
    getRatio(ratio_abs1, direct_region, reflected_region, hists.at(0), verbose);

    double ratio_min1;
    double ratio_err2;
    double S_xy = 0.0;
    double S_xx = 0.0;
    double S_yy = 0.0;
    std::vector<double> ratio_res = {0.0, 0.0};
    for (unsigned int n = 1; n < abs_scalings_min1.size(); ++n) {
        // ratio for other abs
        // if (verbose) std::cout << "Getting ratio for idx = " << n << std::endl;
        getRatio(ratio_res, direct_region, reflected_region, hists.at(n), verbose);
        ratio_min1 = (ratio_res.at(0) / ratio_abs1.at(0)) - 1.0; // normalise ratio to value at abs=1.0, then subtract 1
        ratio_err2 = (ratio_res.at(1)*ratio_res.at(1)+ (ratio_res.at(0)*ratio_res.at(0)
                     / (ratio_abs1.at(0)*ratio_abs1.at(0))) * ratio_abs1.at(1)*ratio_abs1.at(1))
                     / (ratio_abs1.at(0)*ratio_abs1.at(0)); // update errors (squared) accordingly

        // Add to quantities used to find slope (FOM) of ratio vs abs
        S_xy += (abs_scalings_min1.at(n) * ratio_min1) / ratio_err2;
        S_xx += (abs_scalings_min1.at(n) * abs_scalings_min1.at(n)) / ratio_err2;
        S_yy += (ratio_min1 * ratio_min1) / ratio_err2;
    }

    // Return results
    slope_info.at(0) = S_xy / S_xx;  // fit slope (FOM)
    slope_info.at(1) = (S_yy - (S_xy*S_xy / S_xx)) / ((double)abs_scalings_min1.size() - 1.0);  // reduced Chi squared
}


/**
 * @brief Get the Ratio [#hits in reflected region / #hits in direct region], along with the associated statistical error.
 * 
 * @param ratio_res = {ratio, stats error}, what gets modified to return value
 * @param direct_region 
 * @param reflected_region 
 * @param hist 
 * @param bin_info = {X_min, X_max, Xbin_size, T_min, T_max, Tbin_size} Histogram binning info
 * @return std::vector<double> 
 */
void getRatio(std::vector<double>& ratio_res, const rectangle& direct_region, const rectangle& reflected_region, TH2F* hist, const bool verbose) {

    // Count number of hits in direct and reflected regions
    double direct_count = 0.0;
    double reflected_count = 0.0;
    double xBinCenter;

    // Loop over bins in direct region
    for (unsigned int x = direct_region.X_min_bin(); x <= direct_region.X_max_bin(); ++x) {  // loop over histogram bins
        xBinCenter = hist->GetXaxis()->GetBinCenter(x);
        for (unsigned int t = direct_region.T_min_bin(); t < direct_region.T_max_bin(); ++t) {
            direct_count += hist->GetBinContent(x, hist->GetYaxis()->GetBinCenter(t));
        }
    }

    // Loop over bins in reflected region
    for (unsigned int x = reflected_region.X_min_bin(); x <= reflected_region.X_max_bin(); ++x) {  // loop over histogram bins
        xBinCenter = hist->GetXaxis()->GetBinCenter(x);
        for (unsigned int t = reflected_region.T_min_bin(); t < reflected_region.T_max_bin(); ++t) {
            reflected_count += hist->GetBinContent(x, hist->GetYaxis()->GetBinCenter(t));
        }
    }

    // Error handling
    if (direct_count == 0) {
        std::cout << "ERROR: direct region with no hits :/" << std::endl;
        exit(1);
    }
    if (reflected_count == 0) {
        std::cout << "ERROR: reflected region with no hits :/" << std::endl;
        exit(1);
    }

    // Compute ratio of both, and associated stat error (assuming stat errors for regions of sqrt(N))
    ratio_res.at(0) = reflected_count / direct_count;
    ratio_res.at(1) = (sqrt(reflected_count) / direct_count) * sqrt(1.0 + (reflected_count / direct_count));
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
