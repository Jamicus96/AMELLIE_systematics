#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include "AMELLIE_utils.hpp"


std::vector<double> OptimiseDivideAndConquer(std::vector<std::string> traking_files, std::vector<double> abs_scalings, int abs1_idx, bool verbose);
std::vector<double> GetThreePoints(double bestPoint, double worstPoint, std::vector<double> originalPoints);
std::vector<double> GetFOMs(std::vector<double> points, std::vector<double> fixedPoints, int numVar, std::vector<HistList> hist_lists, std::vector<double> abs_scalings, int abs1_idx);
std::vector<double> getRatio(rectangle direct_region, rectangle reflected_region, TH2F *allPathsHist);
std::vector<double> GetBestFOM(std::vector<double> FOMs, std::vector<double> points);

int main(int argc, char** argv){
    // Read in text file with list of stats that were simulated (different abs, the rest the same for now)
    std::string info_file = argv[1];
    std::string info_file_repo = argv[2];
    bool verbose = std::stoi(argv[3]);

    // Create list of tracking hist file name from info file, and list abs_scalings
    std::ifstream file(info_file);
    std::string str;
    std::vector<std::string> traking_files;
    std::vector<double> abs_scalings;
    int abs1_idx = -1;
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
        std::string file_address = info_file_repo + "tot_AMELLIE_" + geo_name + "_" + str_lst.at(1) + "_" + str_lst.at(2) + "_reemis" + str_lst.at(3) + "_abs" + str_lst.at(4) + ".root";
        traking_files.push_back(file_address);
        abs_scalings.push_back(std::stod(str_lst.at(4)));
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

    auto t1 = std::chrono::high_resolution_clock::now();
    std::vector<double> lims_FOM = OptimiseDivideAndConquer(traking_files, abs_scalings, abs1_idx, verbose); 
    auto t2 = std::chrono::high_resolution_clock::now();

    std::cout << "Script took " << std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count() / 1E6 << " s to execute" << std::endl;

    std::cout << "The best region limits , followed by the FOM (slope) are:" << std::endl;    
    for (unsigned int i = 0; i < lims_FOM.size(); ++i) {
        std::cout << lims_FOM.at(i) << std::endl;
    }
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

    return 0;
}


std::vector<double> OptimiseDivideAndConquer(std::vector<std::string> traking_files, std::vector<double> abs_scalings, int abs1_idx, bool verbose){

    // Read in histograms
    std::vector<HistList> hist_lists;
    for (unsigned int i = 0; i < traking_files.size(); ++i) {
        hist_lists.push_back(HistList(traking_files.at(i)));
    }

    //set up constants (note dx and dy denote HALF the width and height of the rectangular regions, respectively)
    std::string point_names[8] = {"direct_x_c", "direct_y_c", "direct_dx", "direct_dy", "reflected_x_c", "reflected_y_c", "reflected_dx", "reflected_dy"};
    double y_min = hist_lists.at(0).Tracking_Hists().at(0)->GetYaxis()->GetXmin();
    double y_max = hist_lists.at(0).Tracking_Hists().at(0)->GetYaxis()->GetXmax();
    double point_mins[8] = {-1.0, y_min, 0.1, 0.1, 0.0, y_min, 0.1, 0.1};
    double point_maxs[8] = {0.0, y_max, 1.0, 0.25 * (y_max - y_min), 0.5, y_max, 0.5, 0.25 * (y_max - y_min)};

    std::vector<double> point_tolerances = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double xBinWidth = hist_lists.at(0).Tracking_Hists().at(0)->GetXaxis()->GetBinCenter(2) - hist_lists.at(0).Tracking_Hists().at(0)->GetXaxis()->GetBinCenter(1);
    double yBinWidth = hist_lists.at(0).Tracking_Hists().at(0)->GetYaxis()->GetBinCenter(2) - hist_lists.at(0).Tracking_Hists().at(0)->GetYaxis()->GetBinCenter(1);
    for (unsigned int i = 0; i < 4; ++i) {
        point_tolerances.at(2*i) = xBinWidth;
        point_tolerances.at(2*i + 1) = yBinWidth;
    }

    std::vector<double> points_diffs = {9999999999, 9999999999, 9999999999, 9999999999, 9999999999, 9999999999, 9999999999, 9999999999};
    std::vector<double> points_temp_diffs = {9999999999, 9999999999, 9999999999, 9999999999, 9999999999, 9999999999, 9999999999, 9999999999};

    // Initialise variables:
    // fixedPoints = {(direct_x_c_max - direct_x_c_min) / 2, (direct_y_c_max - direct_y_c_min) / 2, direct_dx_max, direct_dy_max,
    //          (reflected_x_c_max - reflected_x_c_min) / 2, (reflected_y_c_max - reflected_y_c_min) / 2, reflected_dx_max, reflected_dy_max};
    std::vector<double> fixedPoints = {0.5 * (point_maxs[0] - point_mins[0]), 0.5 * (point_maxs[1] - point_mins[1]), point_maxs[2], point_mins[3],
                                    0.5 * (point_maxs[4] - point_mins[4]), 0.5 * (point_maxs[5] - point_mins[5]), point_maxs[6], point_mins[7]};
    bool firstRun = true;

    
    int numMainLoopIterations = 0;
    double prev_best_points_main[8];
    bool first_runs[8];
    double prevBestFOMPoints[8];
    double temp_diffs[8];
    std::vector<double> FOMs[8] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    std::vector<double> bestworstFOMPoints[8] = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};

    std::vector<double> points[8];
    for (unsigned int i = 0; i < 8; ++i) {
        points[i] = {point_mins[i], (point_maxs[i] + point_mins[i])/2., point_maxs[i]};
    }

    bool loop_condition = true; // If any point difference is still larger than its associated tolerance, keep looping. See end of while loop.
    unsigned int N = 8;
    while(loop_condition or firstRun){
        // Set/Reset variables
        for (unsigned int i = 0; i < N; ++i) {
            first_runs[i] = true;
            prevBestFOMPoints[i] = 9999999999;
            temp_diffs[i] = 9999999999;
        }

        if(verbose) std::cout << "In main while loop with direct_x_c_temp_diff: " << std::abs(temp_diffs[0]) << std::endl;

        for (unsigned int i = 0; i < N; ++i) {
            while(std::abs(temp_diffs[i]) > point_tolerances.at(i)){
                if(first_runs[i]){
                    points[i] = {point_mins[i], (point_maxs[i] + point_mins[i])/2., point_maxs[i]}; //l, c, r
                    first_runs[i] = false;
                    if(verbose) std::cout << "Set up first " << point_names[i] << "_run points: " << points[i].at(0)
                                        << ", " << points[i].at(1) << ", " << points[i].at(2) << std::endl;
                }
                FOMs[i] = GetFOMs(points[i], fixedPoints, i, hist_lists, abs_scalings, abs1_idx);  //hAllPaths, hReEmittedPaths, hSingleScatterPaths
                if(verbose) std::cout << "Got FOMs " << FOMs[i].at(0) << ", " << FOMs[i].at(1) << ", " << FOMs[i].at(2) << std::endl;
                
                bestworstFOMPoints[i] = GetBestFOM(FOMs[i], points[i]);
                if(verbose) std::cout << "Got best FOMs" << std::endl;
                if(verbose) std::cout << "Best slope: " << bestworstFOMPoints[0] << ", best FOM: " << bestworstFOMPoints[2] << std::endl;
                
                points[i] = GetThreePoints(bestworstFOMPoints[i].at(0), bestworstFOMPoints[i].at(1), points[i]);
                if(verbose) std::cout << "Got new points: " << points[i].at(0) << ", " << points[i].at(1) << ", " << points[i].at(2) << std::endl;
                
                temp_diffs[i] = std::abs(points[i].at(1) - prevBestFOMPoints[i]);
                if(verbose) std::cout << "Got difference: " << temp_diffs[i] << std::endl;
                
                prevBestFOMPoints[i] = points[i].at(1);
            }
            if(verbose) std::cout << "Setting fixed point: " << bestworstFOMPoints[i].at(0) << std::endl;
            fixedPoints.at(i) = bestworstFOMPoints[i].at(0);
            if(verbose) std::cout << "" << std::endl;
        }

        if(verbose) std::cout << "Getting main differences" << std::endl;

        for (unsigned int i = 0; i < N; ++i) {
            if(points_diffs.at(i) == 9999999999){
                points_diffs.at(i) = std::abs(bestworstFOMPoints[i].at(0));
                prev_best_points_main[i] = bestworstFOMPoints[i].at(0);
            }
            else{
                points_diffs.at(i) = std::abs(prev_best_points_main[i] - bestworstFOMPoints[i].at(0));
                prev_best_points_main[i] = bestworstFOMPoints[i].at(0);
                firstRun = false;
            }
        }
        
        if(verbose) std::cout << "Got main differences" << std::endl;
        // If any point difference is still larger than its associated tolerance, keep looping.
        loop_condition = false;
        for (unsigned int i = 0; i < N; ++i) {
            if (points_diffs.at(i) > point_tolerances.at(i)) {
                loop_condition = true;
                break;
            }
        }
        numMainLoopIterations++;
    }

    // Compute final region limits
    double direct_x_max = fixedPoints[0] + fixedPoints[2];
    double direct_x_min = fixedPoints[0] - fixedPoints[2];
    double direct_y_max = fixedPoints[1] + fixedPoints[3];
    double direct_y_min = fixedPoints[1] - fixedPoints[3];
    double reflected_x_max = fixedPoints[4] + fixedPoints[6];
    double reflected_x_min = fixedPoints[4] - fixedPoints[6];
    double reflected_y_max = fixedPoints[5] + fixedPoints[7];
    double reflected_y_min = fixedPoints[5] - fixedPoints[7];
    std::vector<double> lims_FOM = {direct_x_max, direct_x_min, direct_y_max, direct_y_min,
                                reflected_x_max, reflected_x_min, reflected_y_max, reflected_y_min, bestworstFOMPoints[N-1].at(2)};

    return lims_FOM;
}

/**
 * @brief Remove the point in points with the worst FOM, and recalculate the mid point from the two remaining.
 * If the middle point is the worst, remove the second worst instead.
 * 
 * @param bestPoint Point with best FOM
 * @param worstPoint Point with worst FOM
 * @param originalPoints (z_i_min, z_i_mid, z_i_max), for z=x,y, and i in {a, b, c}
 * @return std::vector<double> 
 */
std::vector<double> GetThreePoints(double bestPoint, double worstPoint, std::vector<double> originalPoints){
    std::vector<double> newPoints;
    if(worstPoint == originalPoints.at(0)){ //left point worst
        newPoints = {originalPoints.at(1), (originalPoints.at(1) + originalPoints.at(2)) / 2., originalPoints.at(2)};
    }
    else if(worstPoint == originalPoints.at(2)){ //right point worst
        newPoints = {originalPoints.at(0), (originalPoints.at(0) + originalPoints.at(1)) / 2., originalPoints.at(1)};
    }
    else if(worstPoint == originalPoints.at(1)){ //middle point worst, remove next worse point
        if(originalPoints.at(0) == bestPoint){ //right point worst
            newPoints = {originalPoints.at(0), (originalPoints.at(0) + originalPoints.at(1)) / 2., originalPoints.at(1)};
        }
        else if(originalPoints.at(2) == bestPoint){ //left point worst
            newPoints = {originalPoints.at(1), (originalPoints.at(1) + originalPoints.at(2)) / 2., originalPoints.at(2)};
        }
    }
    else{
        std::cout << "Error! No worst point" <<std::endl;
    }

    return newPoints;
}

/**
 * @brief Get Figures of Merit (FOMs) from replacing one fixed point of a rectangle with three different options from points. FOM is the slope
 * found from plotting the ratio [#hits in reflected region / #hits in direct region] (normalised to value at abs=1.0) vs absorptions.
 * 
 * @param points 
 * @param fixedPoints 
 * @param numVar 
 * @param hist_lists 
 * @param abs_scalings 
 * @param abs1_idx 
 * @return std::vector<double> 
 */
std::vector<double> GetFOMs(std::vector<double> points, std::vector<double> fixedPoints, int numVar, std::vector<HistList> hist_lists, std::vector<double> abs_scalings, int abs1_idx){

    // Create rectangles with fixed points (x_max = x_c + dx, x_min = x_c - dx, and likewise for y)
    rectangle direct_region = rectangle(fixedPoints[0] + fixedPoints[2], fixedPoints[0] - fixedPoints[2],
                                        fixedPoints[1] + fixedPoints[3], fixedPoints[1] - fixedPoints[3]);
    rectangle reflected_region = rectangle(fixedPoints[4] + fixedPoints[6], fixedPoints[4] - fixedPoints[6],
                                        fixedPoints[5] + fixedPoints[7], fixedPoints[5] - fixedPoints[7]);

    // Replace the appropriate point in the rectangles with each point in points and get the ratio of reflected/direct number of hits.
    unsigned int N = abs_scalings.size();
    std::vector<double> temp_points = fixedPoints;
    std::vector<double> ratio_res = {0.0, 0.0};
    double ratio;
    double ratio_err;
    std::vector<double> ratio_abs1;
    double S_xy;
    double S_xx;
    std::vector<double> FOMs = {0.0, 0.0, 0.0}; // slopes
    bool zero = false;
    for (unsigned int i = 0; i < 3; ++i){
        temp_points = fixedPoints;
        temp_points.at(numVar) = points.at(i);
        direct_region = rectangle(temp_points[0] + temp_points[2], temp_points[0] - temp_points[2],
                                            temp_points[1] + temp_points[3], temp_points[1] - temp_points[3]);
        reflected_region = rectangle(temp_points[4] + temp_points[6], temp_points[4] - temp_points[6],
                                            temp_points[5] + temp_points[7], temp_points[5] - temp_points[7]);

        ratio_abs1 = getRatio(direct_region, reflected_region, hist_lists.at(abs1_idx).Tracking_Hists().at(0));
        S_xy = 0.0;
        S_xx = 0.0;

        for (unsigned int n = 0; n < N; ++n) {
            ratio_res = getRatio(direct_region, reflected_region, hist_lists.at(n).Tracking_Hists().at(0));
            if (ratio_res.at(0) == -1.0) {
                zero = true;
                break;
            } else {
                ratio = ratio_res.at(0) / ratio_abs1.at(0); // normalise ratio to value at abs=1.0
                ratio_err = sqrt(((ratio_res.at(1)*ratio_res.at(1)) / (ratio_abs1.at(0)*ratio_abs1.at(0)))
                                        + ((ratio_res.at(0)*ratio_res.at(0) * ratio_abs1.at(1)*ratio_abs1.at(1))
                                        / (ratio_abs1.at(0)*ratio_abs1.at(0)*ratio_abs1.at(0)*ratio_abs1.at(0)))); // update errors accordingly

                // Add to quantities used to find slope of ratio vs abs
                S_xy += ((abs_scalings.at(n) - 1.0) * (ratio - 1.0)) / (ratio_err*ratio_err);
                S_xx += ((abs_scalings.at(n) - 1.0) * (abs_scalings.at(n) - 1.0)) / (ratio_err*ratio_err);
            }
        }
        if (zero) {
            FOMs.at(i) = -999999.0; // Bad region: one of them had zero hits.
            zero = false;
        } else {
            FOMs.at(i) = S_xy / S_xx;
        }
    }

    return FOMs;
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
            }
            if(reflected_region.check_point_inside_rectangle(xBinCenter, yBinCenter)){
                reflected_count += allPathsHist->GetBinContent(x,y);
            }
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
 * @brief Compares the FOM of each point (of 3) in points, and returns a vector with:
 * (point with best (highest) FOM, point with worst (lowest) FOM, FOM of best point, FOM of worst point).
 * 
 * @param FOMs Figures of merit, See GetFOM() function.
 * @param points (z_i_min, z_i_mid, z_i_max), for z=x,y, and i in {a, b, c}
 * @return std::vector<double> 
 */
std::vector<double> GetBestFOM(std::vector<double> FOMs, std::vector<double> points){
    std::vector<double> output;
    if(FOMs.at(0) >= FOMs.at(1) and FOMs.at(0) >= FOMs.at(2)){ //left side best point
        output.push_back(points.at(0));
        if(FOMs.at(1) >= FOMs.at(2)){ //right point worst
            output.push_back(points.at(2));
            output.push_back(FOMs.at(0));
            output.push_back(FOMs.at(2));
        }
        else{ // middle point worst
            output.push_back(points.at(1));
            output.push_back(FOMs.at(0));
            output.push_back(FOMs.at(1));
        }
    }
    else if(FOMs.at(1) >= FOMs.at(0) and FOMs.at(1) >= FOMs.at(2)){ //middle best point
        output.push_back(points.at(1));
        if(FOMs.at(0) >= FOMs.at(2)){ //right point worst
            output.push_back(points.at(2));
            output.push_back(FOMs.at(1));
            output.push_back(FOMs.at(2));
        }
        else{ // left point worst
            output.push_back(points.at(0));
            output.push_back(FOMs.at(1));
            output.push_back(FOMs.at(0));
        }
    }
    else if(FOMs.at(2) >= FOMs.at(1) and FOMs.at(2) >= FOMs.at(0)){ //right side best point
        output.push_back(points.at(2));
        if(FOMs.at(1) >= FOMs.at(0)){ //left point worst
            output.push_back(points.at(0));
            output.push_back(FOMs.at(2));
            output.push_back(FOMs.at(0));
        }
        else{ // middle point worst
            output.push_back(points.at(1));
            output.push_back(FOMs.at(2));
            output.push_back(FOMs.at(1));
        }
    }
    else{
        std::cout << "ERROR! No best fit point, points: ";
        for(unsigned int i = 0; i< points.size();i++){
            std::cout << points.at(i) << ", ";
        }
        std::cout << std::endl << "FOMs: ";
        for(unsigned int i = 0; i< FOMs.size();i++){
            std::cout << FOMs.at(i) << ", ";
        }
    }
    return output;
}