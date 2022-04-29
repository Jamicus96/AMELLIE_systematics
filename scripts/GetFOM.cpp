#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include "AMELLIE_utils.hpp"


double GetFOM(std::vector<std::string> traking_files, std::vector<double> abs_scalings, int abs1_idx, double points[8]);
std::vector<double> getRatio(rectangle direct_region, rectangle reflected_region, TH2F *allPathsHist);



int main(int argc, char** argv){
    // Read in text file with list of stats that were simulated (different abs, the rest the same for now)
    std::string info_file = argv[1];
    std::string tracking_hist_repo = argv[2];
    // Read in output filename (txt)
    std::string output_file = argv[3];
    // Read in region limits
    double direct_x_max = std::stod(argv[4]);
    double direct_x_min = std::stod(argv[5]);
    double direct_y_max = std::stod(argv[6]);
    double direct_y_min = std::stod(argv[7]);
    double reflected_x_max = std::stod(argv[8]);
    double reflected_x_min = std::stod(argv[9]);
    double reflected_y_max = std::stod(argv[10]);
    double reflected_y_min = std::stod(argv[11]);
    double points[8] = {direct_x_max, direct_x_min, direct_y_max, direct_y_min,
                            reflected_x_max, reflected_x_min, reflected_y_max, reflected_y_min};

    /* ~~~~~~~~~ Create list of tracking hist file name from info file, and list abs_scalings ~~~~~~~~~ */

    std::ifstream file(info_file);
    std::string str;
    std::vector<std::string> traking_files;
    std::vector<double> abs_scalings;
    int abs1_idx = -1;  // Index in abs_scalings where abs=1.0
    int i = 0;
    while (std::getline(file, str)) {
        // format in info file: geo_file.geo, LEDnum, fibre, reemis, abs
        // format in tracking hist file name: tot_AMELLIE_geoFile_LEDnum_fibre_reemis_abs.root
        std::vector<std::string> str_lst;
        std::size_t pos;
        std::string delimiter = ", ";
        bool splitting = true;
        while (splitting) {
            pos = str.find(delimiter);
            if ((unsigned int)pos < str.size()) {
                str_lst.push_back(str.substr(0, pos));
                str = str.substr(pos + delimiter.size());
            } else {
                splitting = false;
            }
        }
        str_lst.push_back(str);

        pos = str_lst.at(0).find(".geo");
        std::string geo_name = str_lst.at(0).substr(0, pos);
        std::string file_address = tracking_hist_repo + "tot_AMELLIE_" + geo_name + "_" + str_lst.at(1) + "_" + str_lst.at(2) + "_reemis" + str_lst.at(3) + "_abs" + str_lst.at(4) + ".root";
        traking_files.push_back(file_address);
        std::cout << file_address << std::endl;
        abs_scalings.push_back(std::stod(str_lst.at(4)));
        std::cout << std::stod(str_lst.at(4)) << std::endl;
        if (str_lst.at(4) == "1.0" || abs_scalings.at(i) == 1.0) {
            abs1_idx = i;
        }
        str_lst.clear();
        ++i;
    }

    if (abs1_idx == -1){
        std::cout << "ERROR: Absorption 1.0 case not found." << std::endl;
        exit(1);
    }

    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

    // Get FOM
    double FOM = GetFOM(traking_files, abs_scalings, abs1_idx, points);

    // Print it to file, after region info
    std::ofstream datafile;
    datafile.open(output_file.c_str(), std::ios::app);
    for (int i = 0; i < 8; ++i) {
        datafile << points[i] << " ";
    }
    datafile << FOM << std::endl;
}



/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


double GetFOM(std::vector<std::string> traking_files, std::vector<double> abs_scalings, int abs1_idx, double points[8]) {

    // Read in histograms
    std::vector<HistList> hist_lists;
    for (unsigned int i = 0; i < traking_files.size(); ++i) {
        hist_lists.push_back(HistList(traking_files.at(i)));
    }

    // Set up rectangle regions
    rectangle direct_region = rectangle(points[0], points[1], points[2], points[3]);
    rectangle reflected_region = rectangle(points[4], points[5], points[6], points[7]);

    // Get ratio of hits in both triangular regions, and associated error, then normalise all to the ratio
    // found at abs = 1. Then compute best fit slope of normalised ratio vs abs:

    // abs = 1 ratio
    std::vector<double> ratio_abs1 = getRatio(direct_region, reflected_region, hist_lists.at(abs1_idx).Tracking_Hists().at(1));

    double ratio;
    double ratio_err;
    double S_xy;
    double S_xx;
    std::vector<double> ratio_res = {0.0, 0.0};
    for (unsigned int n = 0; n < abs_scalings.size(); ++n) {
        // ratio for other abs
        ratio_res = getRatio(direct_region, reflected_region, hist_lists.at(n).Tracking_Hists().at(1));
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
 * @param allPathsHist 
 * @return std::vector<double> 
 */
std::vector<double> getRatio(rectangle direct_region, rectangle reflected_region, TH2F *allPathsHist) {

    // Count number of hits in direct and reflected regions
    double direct_count = 0.0;
    double reflected_count = 0.0;
    for(int x=1; x<allPathsHist->GetNbinsX()+1; x++){ //loop over histogram bins
        double xBinCenter = allPathsHist->GetXaxis()->GetBinCenter(x);
        for(int y=1; y<allPathsHist->GetNbinsY()+1; y++){
            double yBinCenter = allPathsHist->GetYaxis()->GetBinCenter(y);
            if(direct_region.check_point_inside_rectangle(xBinCenter, yBinCenter)){
                direct_count += allPathsHist->GetBinContent(x,y);
                // std::cout << "+d" << std::endl;
            }
            if(reflected_region.check_point_inside_rectangle(xBinCenter, yBinCenter)){
                reflected_count += allPathsHist->GetBinContent(x,y);
                // std::cout << "+r" << std::endl;
            }
        }
    }

    // std::cout << "direct_count = " << direct_count << std::endl;
    // std::cout << "reflected_count = " << reflected_count << std::endl;

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