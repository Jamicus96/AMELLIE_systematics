#include "AMELLIE_utils.hpp"


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TRIANGLE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/**
 * @brief Construct a new triangle::triangle object
 * 
 * @param x_a 
 * @param x_b 
 * @param x_c 
 * @param y_a 
 * @param y_b 
 * @param y_c 
 */
triangle::triangle(double x_a, double x_b, double x_c, double y_a, double y_b, double y_c) {
    points[0] = x_a; points[1] = x_b; points[2] = x_c;
    points[3] = y_a; points[4] = y_b; points[5] = y_c;
    // // Twice the area of the triangle (twice to save on computation later)
    // Area2 = fabs(points[0] * (points[4] - points[5]) + points[1] * (points[5] - points[3])
    //             + points[2] * (points[3] - points[4]));
}

// member functions

// Get/Set triangle vertex coordinates
double& triangle::X_a() {return points[0];}
double& triangle::X_b() {return points[1];}
double& triangle::X_c() {return points[2];}
double& triangle::Y_a() {return points[3];}
double& triangle::Y_b() {return points[4];}
double& triangle::Y_c() {return points[5];}

/**
 * @brief Checks if point is inside triangle by using the sum of triangles method:
 * if the sum of the area of the three triangles formed by the point and two vertices of the triangle add up to
 * the area of the original triangle, the point is inside. It is outside if the sum is larger. 
 * 
 * @param point_x 
 * @param point_y 
 * @param lim_count true (default) = point on triangle edge counts as inside, false = point on edge counts as outside.
 * @return true 
 * @return false 
 */
bool triangle::check_point_inside_triangle(const double point_x, const double point_y, const bool print) {
    // Twice the area of the triangle (twice to save on computation later)
    double Area2 = fabs(points[0] * (points[4] - points[5]) + points[1] * (points[5] - points[3])
                    + points[2] * (points[3] - points[4]));
    if (Area2 == 0.0) {
        // std::cout << "Triangle area is zero." << std::endl;
        return false;
    }

    // Calculate (twice) the area of each triangle formed by the point, and check if it is zero successively
    double A_a = fabs(point_x * (points[4] - points[5]) + points[1] * (points[5] - point_y)
                    + points[2] * (point_y - points[4]));
    if (A_a == 0.0) {return false;} // if area is zero, point is on edge of triangle.

    double A_b = fabs(points[0] * (point_y - points[5]) + point_x * (points[5] - points[3])
                    + points[2] * (points[3] - point_y));
    if (A_b == 0.0) {return false;} // if area is zero, point is on edge of triangle.

    double A_c = fabs(points[0] * (points[4] - point_y) + points[1] * (point_y - points[3])
                    + point_x * (points[3] - points[4]));
    if (A_c == 0.0) {return false;} // if area is zero, point is on edge of triangle.

    // Check sum of areas of point triangles adds up to original triangle aread (tolerance of 0.1%).
    // If they are equal, the point is inside the triangle. If the sum is larger, the point is outside.
    double frac_diff = ((A_a + A_b + A_c) / Area2) - 1.0;

    // print info if desired
    if (print) {
        std::cout << "2*Area = " << Area2 << std::endl;
        std::cout << "2*Area_a = " << A_a << std::endl;
        std::cout << "2*Area_b = " << A_b << std::endl;
        std::cout << "2*Area_c = " << A_c << std::endl;
        std::cout << "frac_diff = " << frac_diff << std::endl;
    }
    
    if (frac_diff <= 0.001) {
        return true;
    } else {
        return false;
    }
}


/* ~~~~~~~~~ operator overload ~~~~~~~~~ */
double& triangle::operator [] (int i) {
    assert((i < 6 && i >= 0) && "Index out of range");
    return points[i];
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ RECTANGLE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/**
 * @brief Construct a new rectangle::rectangle object.
 * 
 * @param region_type = "direct" or "reflected".
 * @param fibre to adjust direct region x-limit.
 * @param reflected_max_time to set y-limits.
 */
rectangle::rectangle(const std::string region_type, const std::string fibre, double max_y_bin) {
    y_max = max_y_bin + 10;
    y_min = max_y_bin - 10;
    if (region_type == "direct") {
        x_min = -1.0;
        double dx;
        if (fibre == "FA089") {  // 10deg off-axis
            dx = 0.15;
        } else if (fibre == "FA173" or fibre == "FA150" or fibre == "FA093") {  // 20deg off-axis
            dx = 0.4;
        } else {  // on-axis
            dx = 0.1;
        }
        x_max = x_min + dx;
    } else if (region_type == "reflected") {
        x_max = 1.0;
        x_min = 0.95;
    } else {
        std::cout << "Wrong type, should be 'direct' or 'reflected', not '" << region_type << "'." << std::endl;
        exit(1);
    }
}


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
bool rectangle::check_point_inside_rectangle(const double point_x, const double point_y) {
    if (point_x > x_max or point_x < x_min or point_y > y_max or point_y < y_min) {
        return false;
    } else {
        return true;
    }
}



/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ HISTLIST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/**
 * @brief Construct a new Hist List:: Hist List object.
 * Default constructor
 * 
 */
HistList::HistList() {length = 0;}

/**
 * @brief Construct a new Hist List:: Hist List object.
 * Creates vector of 2D Tracking Hists in root file, and clones them to
 * 3 new vectors of 2D hists with modified names.
 * 
 * @param tracking_file Root file containing hists.
 * @param name_list List of histogram names the user wishes to read in. Default list defined
 * in header file. (FIXME: can add option to just read in all 2D hists)
 */
HistList::HistList(std::string tracking_file, std::vector<std::string> name_list) {
    Read_File(tracking_file, name_list);
}

/**
 * @brief Basically constructs a new HistList object.
 * Creates vector of 2D Tracking Hists in root file, and clones them to
 * 3 new vectors of 2D hists with modified names.
 * 
 * @param tracking_file Root file containing hists.
 * @param name_list List of histogram names the user wishes to read in. Default list defined
 * in header file. (FIXME: can add option to just read in all 2D hists)
 */
void HistList::Read_File(std::string tracking_file, std::vector<std::string> name_list) {
    // Read in root file
    TFile* fin = new TFile(tracking_file.c_str());
    if (!fin->IsOpen()) {
        std::cout << "Cannot open input file " << tracking_file << std::endl;
        exit(1);
    }
    length = name_list.size();
    
    // Iterate through list of objects in root file
    TList* list = fin->GetListOfKeys() ;
    if (!list) {std::cout << "No keys found in file\n" << std::endl; exit(1);}
    TIter next(list);
    TObject* obj;
    TKey* key;
    int j = 0;
    std::string name;
    int order[name_list.size()];
    for (unsigned int i = 0; i < length; ++i) {
        order[i] = -1;
    }
    std::vector<TH2F*> temp_HistList;

    // Go through list of histograms and add them to temp list if they are included in name list
    while((key = (TKey*)next())){
        obj = key->ReadObj() ;
        if(obj->InheritsFrom(TH2::Class())){
            // Check which histogram in file matches name
            name = (std::string)(((TH2F*)obj)->GetName());
            for (unsigned int i = 0; i < length; ++i) {
                if (name == name_list.at(i)) {
                    // Add tracking histograms to temporary list
                    temp_HistList.push_back((TH2F*)obj);
                    order[i] = j;
                    ++j;
                    break;
                }
            }
        }
    }

    // Reorder histograms, and clone new ones in that order
    for (unsigned int i = 0; i < length; ++i) {
        if (order[i] == -1) {
            std::cout << "Histogram '" << name_list.at(i) << "' not found." << std::endl;
            exit(1);
        } else {
            // Add tracking histograms to list
            tracking_hists.push_back(temp_HistList.at(order[i]));

            // Clone tracking hists to other lists
            name = name_list.at(i).erase(0,1); // removes first character ("h" in this case)
            region_hists.push_back((TH2F*)(tracking_hists.at(i))->Clone());
            region_hists.at(i)->SetName(((std::string)"hRegion" + name).c_str());
            direct_hists.push_back((TH2F*)(tracking_hists.at(i))->Clone());
            direct_hists.at(i)->SetName(((std::string)"hDirect" + name).c_str());
            reflected_hists.push_back((TH2F*)(tracking_hists.at(i))->Clone());
            reflected_hists.at(i)->SetName(((std::string)"hReflected" + name).c_str());
        }
    }
}

// return suitable histograms / histogram names
unsigned int HistList::len() {return length;}
std::vector<TH2F*>& HistList::Tracking_Hists() {return tracking_hists;}
std::vector<TH2F*>& HistList::Region_Hists() {return region_hists;}
std::vector<TH2F*>& HistList::Direct_Hists() {return direct_hists;}
std::vector<TH2F*>& HistList::Reflected_Hists() {return reflected_hists;}

/**
 * @brief Write all  new (Region, Direct, Reflected) hists to file
 * 
 */
void HistList::Write() {
    for (unsigned int i = 0; i < length; ++i) {
        region_hists.at(i)->Write();
        direct_hists.at(i)->Write();
        reflected_hists.at(i)->Write();
    }
}