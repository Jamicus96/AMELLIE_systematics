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


std::vector<std::vector<unsigned int> > GetStats(const rectangle& direct_region, const rectangle& reflected_region, const std::vector<std::string>& hist_file_addresses, bool verbose);
void DrawRegionLims(const rectangle& direct_region, const rectangle& reflected_region, const std::vector<std::string>& hist_file_addresses, const std::string& output_file);


int main(int argc, char** argv) {
    // Read in output file address (txt)
    std::string output_file = argv[1];
    // Read in region limits
    double direct_x_max = std::stod(argv[2]);
    double direct_x_min = std::stod(argv[3]);
    double direct_y_max = std::stod(argv[4]);
    double direct_y_min = std::stod(argv[5]);
    double reflected_x_max = std::stod(argv[6]);
    double reflected_x_min = std::stod(argv[7]);
    double reflected_y_max = std::stod(argv[8]);
    double reflected_y_min = std::stod(argv[9]);
    // record extra info?
    bool verbose = std::stoi(argv[10]);
    // Addresses of hist files to be analysed
    std::vector<std::string> hist_file_addresses;
    for (unsigned int i = 11; i < argc; ++i) {
        hist_file_addresses.push_back(argv[i]);
    }
    
    if (verbose) {std::cout << "Number of files being analysed: " << hist_file_addresses.size() << std::endl;}

    // Set up regions
    rectangle direct_region = rectangle(direct_x_max, direct_x_min, direct_y_max, direct_y_min);
    rectangle reflected_region = rectangle(reflected_x_max, reflected_x_min, reflected_y_max, reflected_y_min);

    // Get number of hits in different regions of each hist
    std::vector<std::vector<unsigned int> > results = GetStats(direct_region, reflected_region, hist_file_addresses, verbose);
    std::vector<unsigned int> tot_hits = results.at(0);
    std::vector<unsigned int> direct_hits = results.at(1);
    std::vector<unsigned int> reflected_hits = results.at(2);

    // Save hists with regions overlain
    if (verbose) {DrawRegionLims(direct_region, reflected_region, hist_file_addresses, output_file);}

    // Print it to file, after region info
    std::ofstream datafile;
    datafile.open(output_file.c_str(), std::ios::trunc);
    // then print stats
    for (unsigned int i = 0; i < tot_hits.size(); ++i) {
        datafile << tot_hits.at(i) << " " << direct_hits.at(i) << " " << reflected_hits.at(i) << std::endl;
    }
    datafile.close();
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


/**
 * @brief Returns number of hits in direct and reflected regions, and in overall histogram
 * 
 * @param direct_region 
 * @param reflected_region 
 * @param hist_file_addresses 
 * @param verbose 
 * @return std::vector<std::vector<unsigned int> > 
 */
std::vector<std::vector<unsigned int> > GetStats(const rectangle& direct_region, const rectangle& reflected_region, const std::vector<std::string>& hist_file_addresses, bool verbose) {

    std::vector<unsigned int> tot_hits, direct_hits, reflected_hits;
    unsigned int NbinsX, NbinsY;
    TH2F* hist;
    for (unsigned int n = 0; n < hist_file_addresses.size(); ++n) {
        // Read in histogram
        hist = GetHist(hist_file_addresses.at(n), "hPmtResTimeVsCosTheta");

        // Count number of hits in direct and reflected regions
        tot_hits.push_back(0);
        direct_hits.push_back(0);
        reflected_hits.push_back(0);
        NbinsX = hist->GetNbinsX();
        NbinsY = hist->GetNbinsY();
        if (verbose) {std::cout << "NbinsX = " << NbinsX << ", NbinsY = " << NbinsY << std::endl;}
        for (int x = 1; x < NbinsX + 1; x++) {  // loop over histogram bins
            double xBinCenter = hist->GetXaxis()->GetBinCenter(x);
            for (int y = 1; y < NbinsY + 1; y++) {
                double yBinCenter = hist->GetYaxis()->GetBinCenter(y);
                tot_hits.at(n) += hist->GetBinContent(x, y);
                if (direct_region.check_point_inside_rectangle(xBinCenter, yBinCenter)) {
                    direct_hits.at(n) += hist->GetBinContent(x, y);
                }
                if (reflected_region.check_point_inside_rectangle(xBinCenter, yBinCenter)) {
                    reflected_hits.at(n) += hist->GetBinContent(x, y);
                }
            }
        }
    }

    // package results
    std::vector<std::vector<unsigned int> > results = {tot_hits, direct_hits, reflected_hits};
    return results;
}


/**
 * @brief Draw 2-D histogram from tracking info with regions overlain.
 * 
 * @param direct_region 
 * @param reflected_region 
 * @param hist_file_addresses 
 * @param output_file 
 */
void DrawRegionLims(const rectangle& direct_region, const rectangle& reflected_region, const std::vector<std::string>& hist_file_addresses, const std::string& output_file) {

    // Create output root file for drawn regions
    std::string root_fileName = output_file.substr(0, output_file.find_last_of(".")) + ".root";
    TFile *root_file = new TFile(root_fileName.c_str(), "RECREATE");
    root_file->cd();

    TH2F* hist;
    std::size_t slash_idx, dot_idx;
    std::string hist_filename_str;
    for (unsigned int n = 0; n < hist_file_addresses.size(); ++n) {
        // Read in histogram
        hist = GetHist(hist_file_addresses.at(n), "hPmtResTimeVsCosTheta");

        // Get file name for texting (without extension)
        slash_idx = hist_file_addresses.at(n).find_last_of("/");
        dot_idx = hist_file_addresses.at(n).find_last_of(".");
        hist_filename_str = hist_file_addresses.at(n).substr(slash_idx+1, dot_idx);

        // Draw box cuts on resthit vs costheta hist
        TCanvas* c1 = new TCanvas(hist_filename_str.c_str(), hist_filename_str.c_str());  // Create output canvas to be saved in output file
        TH2F* h = (TH2F*)hist->Clone();  // Make a copy to use in canvas
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

        // Write to file
        c1->Write();
        root_file->Write();
        delete c1;
    }

    // Close
    root_file->Close();
}
