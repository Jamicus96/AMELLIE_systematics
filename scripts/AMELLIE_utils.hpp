// header guard:
#ifndef AMELLIE_utils
#define AMELLIE_utils

// include
#include <iostream>
#include <fstream>
#include <assert.h>
#include <algorithm>
#include <cmath> 
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TKey.h>
#include <TObject.h>
#include <TList.h>
#include <TVectorD.h>

/**
 * @brief Class to keep track of vertices of a defined triangle on a 2D plane.
 * A member function can thus be used to automatically check whether a point lies
 * inside the triangle.
 * Vertices can be easily accessed and modified.
 * 
 */
class triangle {
    private:
        double points[6]; //x_a, x_b, x_c, y_a, y_b, y_c
        //double Area2;

    public:
        // constructors
        triangle(double x_a, double x_b, double x_c, double y_a, double y_b, double y_c);

        // member functions
        double& X_a();
        double& X_b();
        double& X_c();
        double& Y_a();
        double& Y_b();
        double& Y_c();

        bool check_point_inside_triangle(const double point_x, const double point_y, const bool print=false);

    /* ~~~~~~~~~ operator overload ~~~~~~~~~ */

        // access elements easier
        double& operator [] (int i);
};

/**
 * @brief Class to encode the direct and reflected regions
 * 
 */
class rectangle {
    private:
        double x_max;
        double x_min;
        double y_max;
        double y_min;

    public:
        //constructor
        rectangle(const std::string region_type, const std::string fibre, double max_y_bin);
        rectangle(double X_maxi, double X_mini, double Y_maxi, double Y_mini);

        // Memeber function
        double& X_max();
        double& X_min();
        double& Y_max();
        double& Y_min();

        bool check_point_inside_rectangle(const double point_x, const double point_y);
};

/**
 * @brief Class to read in desired 2D hists of tracking root file and clone 3 more
 * copies of these for the region, direct and reflected histograms. All fours lists
 * are in the HistList object and can be easily accessed and modified.
 * 
 */
class HistList {
    private:
        unsigned int length;
        std::vector<TH2F*> tracking_hists;
        std::vector<TH2F*> region_hists;
        std::vector<TH2F*> direct_hists;
        std::vector<TH2F*> reflected_hists;

    public:
    // constructors
    HistList();
    HistList(std::string tracking_file, std::vector<std::string> name_list
                = {"hReemissionResTimeVsCosTheta", "hPmtResTimeVsCosTheta", "hNoiseResTimeVsCosTheta", 
                "hSingleScatterResTimeVsCosTheta", "hOtherEffectResTimeVsCosTheta", "hNoEffectResTimeVsCosTheta", 
                "hNearReflectResTimeVsCosTheta", "hRopesResTimeVsCosTheta", "hPMTReflectionResTimeVsCosTheta", 
                "hExtWaterScatterResTimeVsCosTheta", "hInnerAvReflectionResTimeVsCosTheta", "hMultipleEffectResTimeVsCosTheta",
                "hAVPipesResTimeVsCosTheta", "hAcrylicScatterResTimeVsCosTheta", "OtherScatterResTimeVsCosTheta"});

    // read_file (does the same as previous constructor)
    void Read_File(std::string tracking_file, std::vector<std::string> name_list
                = {"hReemissionResTimeVsCosTheta", "hPmtResTimeVsCosTheta", "hNoiseResTimeVsCosTheta", 
                "hSingleScatterResTimeVsCosTheta", "hOtherEffectResTimeVsCosTheta", "hNoEffectResTimeVsCosTheta", 
                "hNearReflectResTimeVsCosTheta", "hRopesResTimeVsCosTheta", "hPMTReflectionResTimeVsCosTheta", 
                "hExtWaterScatterResTimeVsCosTheta", "hInnerAvReflectionResTimeVsCosTheta", "hMultipleEffectResTimeVsCosTheta",
                "hAVPipesResTimeVsCosTheta", "hAcrylicScatterResTimeVsCosTheta", "OtherScatterResTimeVsCosTheta"});

    // return suitable histograms
    unsigned int len();
    std::vector<TH2F*>& Tracking_Hists();
    std::vector<TH2F*>& Region_Hists();
    std::vector<TH2F*>& Direct_Hists();
    std::vector<TH2F*>& Reflected_Hists();

    // Write all histograms to fiel
    void Write();
};

//end header guard
#endif