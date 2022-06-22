#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLine.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include "AMELLIE_utils.hpp"


std::vector<std::vector<double> > GetFOMabsSlope(rectangle direct_region, rectangle reflected_region, std::vector<std::string> traking_files, std::vector<double> abs_scalings, int abs1_idx, bool record_info = true,  std::string outRoot_filename = "Regions.root");
std::vector<double> getRatio(rectangle direct_region, rectangle reflected_region, TH2F *allPathsHist);
void DrawRegionLims(rectangle direct_region, rectangle reflected_region, TH2F *allPathsHist, double abs);
std::vector<std::vector<std::string> > readInfoFile(std::string tracking_hist_repo, std::string info_file);



int main(int argc, char** argv){
    // Read in text file with list of stats that were simulated (different abs, the rest the same for now)
    std::string info_file = argv[1];
    std::string tracking_hist_repo = argv[2];
    // Read in output filename (txt)
    std::string output_file = argv[3];
    // record extra info?
    bool record_info = std::stoi(argv[4]);
    std::string output_root_filename = argv[5];
    // Read in region limits
    double direct_x_max = std::stod(argv[6]);
    double direct_x_min = std::stod(argv[7]);
    double direct_y_max = std::stod(argv[8]);
    double direct_y_min = std::stod(argv[9]);
    double reflected_x_max = std::stod(argv[10]);
    double reflected_x_min = std::stod(argv[11]);
    double reflected_y_max = std::stod(argv[12]);
    double reflected_y_min = std::stod(argv[13]);

    // Set up regions
    rectangle direct_region = rectangle(direct_x_max, direct_x_min, direct_y_max, direct_y_min);
    rectangle reflected_region = rectangle(reflected_x_max, reflected_x_min, reflected_y_max, reflected_y_min);

    //Create list of tracking hist file name from info file, and list abs_scalings

    std::vector<std::vector<std::string> > file_info = readInfoFile(tracking_hist_repo, info_file);

    std::vector<std::string> traking_files = file_info.at(0);
    unsigned int abs1_idx = std::stoi(file_info.at(2).at(0));
    std::vector<double> abs_scalings;
    for (unsigned int i = 0; i < file_info.size(); ++i) {
        abs_scalings.push_back(std::stod(file_info.at(1).at(i)));
    }

    // Get FOM/abs slope + FOM points
    std::vector<std::vector<double> > results = GetFOMabsSlope(direct_region, reflected_region, traking_files, abs_scalings, abs1_idx, record_info, output_root_filename);
    double slope = results.at(0).at(0);
    std::vector<double> FOM = results.at(1);  // normalised FOM (reflected / direct)
    std::vector<double> FOM_err = results.at(2);  // normalised FOM stat errors

    // Print it to file, after region info
    std::ofstream datafile;
    datafile.open(output_file.c_str(), std::ios::app);
    // print region limits first
    datafile << direct_x_max << " " << direct_x_min << " " << direct_y_max << " " << direct_y_min << " " << reflected_x_max
             << " " << reflected_x_min << " " << reflected_y_max << " " << reflected_y_min << std::endl;
    // then print slope (final results)
    datafile << slope << std::endl;
    if (record_info) {
        // print normalised FOM for each absorption
        for (unsigned int i = 0; i < abs_scalings.size(); ++i) {
            datafile << abs_scalings.at(i) << " " << FOM.at(i) << " " << FOM_err.at(i) << std::endl;
        }
        std::cout << "num_abs = " << abs_scalings.size() << ", num_FOM" << FOM.size() << std::endl;
    }
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


/**
 * @brief Returns FOM/abs slope (called y)
 * 
 * @param traking_files 
 * @param abs_scalings 
 * @param abs1_idx 
 * @param points 
 * @return double 
 */
std::vector<std::vector<double> > GetFOMabsSlope(rectangle direct_region, rectangle reflected_region, std::vector<std::string> traking_files, std::vector<double> abs_scalings, int abs1_idx, bool record_info, std::string outRoot_filename) {

    // Read in histograms
    std::vector<TH2F*> hists;
    for (unsigned int i = 0; i < traking_files.size(); ++i) {
        hists.push_back(GetHist(traking_files.at(i), "hPmtResTimeVsCosTheta"));
    }

    // Set up output histogram
    TFile *output_file = new TFile(outRoot_filename.c_str(), "RECREATE");
    output_file->cd();

    // Get ratio of hits in both triangular regions, and associated error, then normalise all to the ratio
    // found at abs = 1. Then compute best fit slope of normalised ratio vs abs:

    // abs = 1 ratio
    std::vector<double> ratio_abs1 = getRatio(direct_region, reflected_region, hists.at(abs1_idx));

    double ratio;
    double ratio_err;
    double S_xy;
    double S_xx;
    std::vector<double> ratio_res = {0.0, 0.0};
    std::vector<double> ratios;
    std::vector<double> ratio_errs;
    for (unsigned int n = 0; n < abs_scalings.size(); ++n) {
        // ratio for other abs
        ratio_res = getRatio(direct_region, reflected_region, hists.at(n));
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

        if (record_info) {
            // Create histogram to show regions on 2D t_res vs cos(theta) histogram
            DrawRegionLims(direct_region, reflected_region, hists.at(n), abs_scalings.at(n));
            // Add of vectors
            ratios.push_back(ratio);
            ratio_errs.push_back(ratio_err);
        }
    }

    // write to output file and close
    output_file->Write();
    output_file->Close();
    if (!record_info) {
        // If root file wasn't used, delete it
        int status = remove(output_rootFilename.c_str());
    }

    // package results (main result is slope, but record other things too)
    std::vector<double> slope = {S_xy / S_xx};
    std::vector<std::vector<double> > results = {slope, ratios, ratio_errs};

    return results;  // return best fit slope (FOM)
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


/**
 * @brief Draw 2-D histogram from tracking info with regions overlain.
 * 
 * @param direct_region 
 * @param reflected_region 
 * @param allPathsHist 
 */
void DrawRegionLims(rectangle direct_region, rectangle reflected_region, TH2F *allPathsHist, double abs) {

    // Draw box cuts on resthit vs costheta hist
    std::string canvas_name = "hHitTimeResiduals_regions_" + std::to_string(abs);
    std::string canvas_title = "Hit time residuals vs cos(theta) with Regions Overlain, for abs = " + std::to_string(abs);
    TCanvas* c1 = new TCanvas(canvas_name.c_str(), canvas_title.c_str());  // Create output canvas to be saved in output file
    TH2F* h = (TH2F*)allPathsHist->Clone();  // Make a copy to use in canvas
    h->Draw("colz");  // Draw histogram

    // create lines
    std::vector<TLine> lines;
    // direct box
    lines.push_back(TLine(direct_region.X_min(), direct_region.Y_min(), direct_region.X_min(), direct_region.Y_max()));
    lines.push_back(TLine(direct_region.X_max(), direct_region.Y_min(), direct_region.X_max(), direct_region.Y_max()));
    lines.push_back(TLine(direct_region.X_min(), direct_region.Y_min(), direct_region.X_max(), direct_region.Y_min()));
    lines.push_back(TLine(direct_region.X_min(), direct_region.Y_max(), direct_region.X_max(), direct_region.Y_max()));
    // reflected box
    lines.push_back(TLine(reflected_region.X_min(), reflected_region.Y_min(), reflected_region.X_min(), reflected_region.Y_max()));
    lines.push_back(TLine(reflected_region.X_max(), reflected_region.Y_min(), reflected_region.X_max(), reflected_region.Y_max()));
    lines.push_back(TLine(reflected_region.X_min(), reflected_region.Y_min(), reflected_region.X_max(), reflected_region.Y_min()));
    lines.push_back(TLine(reflected_region.X_min(), reflected_region.Y_max(), reflected_region.X_max(), reflected_region.Y_max()));

    // draw lines
    for(unsigned int i = 0; i < lines.size(); ++i){
        lines[i].SetLineColor(kBlack);
        lines[i].Draw("SAME");
    }

    c1->Write();
    delete c1;
}


/**
 * @brief Reads information from text file with simulation info into.
 * 
 * @param info_file 
 * @return std::vector<std::string> {traking_files, abs_scalings, abs1_idx}
 */
std::vector<std::vector<std::string> > readInfoFile(std::string tracking_hist_repo, std::string info_file) {

    std::ifstream file(info_file);
    std::string str;
    std::string geo_file;  // sans the ".geo" part
    std::string wavelength;  // nm
    std::string fibre;
    std::string reemission;  // reemission fraction of absorbed light
    std::vector<std::string> abs_scalings;  // absorption coefficient normalised so that 1.0 = current value
    std::vector<std::string> traking_files;
    std::vector<std::string> abs1_idx = {"-1"};  // Index in abs_scalings where abs=1.0

    int i = 0;
    std::string delimiter = ", ";
    std::size_t pos;
    std::string substr;
    std::string full_geo;
    while (std::getline(file, str)) {
        // format in info file: geo_file.geo, LEDnum, fibre, reemis, abs
        // format in tracking hist file name: tot_AMELLIE_geoFile_LEDnum_fibre_reemis_abs.root
        
        pos = str.find(delimiter);
        substr = str.substr(0, pos);
        full_geo = str.substr(0, pos);
        pos = full_geo.find(".geo");
        geo_file = full_geo.substr(0, pos);
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
        traking_files.push_back(tracking_hist_repo + "tot_AMELLIE_" + geo_file + "_" + wavelength + "_" + fibre + "_reemis" + reemission + "_abs" + str + ".root");

        if (std::stof(str) == 1.0) {abs1_idx.at(0) = std::to_string(i);}  // found absorption=1 case
        ++i;
    }

    if (abs1_idx.at(0) == "-1") {
        std::cout << "ERROR: Absorption 1.0 case not found." << std::endl;
        exit(1);
    }

    std::vector<std::vector<std::string> > results = {traking_files, abs_scalings, abs1_idx};
    return results;
}
