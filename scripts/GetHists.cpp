#include <iostream>
#include <string>
#include <TCanvas.h>
#include "RAT/DU/DSReader.hh"
#include "RAT/DU/PMTInfo.hh"
#include "RAT/DU/Utility.hh"
#include "RAT/DS/Run.hh"
#include "RAT/DS/Entry.hh"
#include <RAT/PhysicsUtil.hh>
#include <TH2F.h>
#include <TH1F.h>
#include <TF1.h>
#include <TMath.h>
#include <TFile.h>
#include <TVector.h>
#include <TGraph2D.h>


void GetLightPaths(std::string input_file, std::string output_file, std::string fibre, double wavelength, bool verbose);


int main(int argc, char** argv){
    std::string input_file = argv[1];
    std::string output_file = argv[2];
    std::string fibre = argv[3];
    double wavelength = std::stod(argv[4]) * 1E-6; // Conver wavelength from nm to mm
    bool verbose = std::stoi(argv[5]);
    GetLightPaths(input_file, output_file, fibre, wavelength, verbose);
    return 0;
}


/**
 * @brief Get the Light Paths object
 * 
 * @param input_file Input simulation root file
 * @param output_file Output root file
 * @param fibre AMELLIE fibre that was fired in the simulation
 * @param wavelength in mm I think (403E-6)
 * @param verbose
 */
void GetLightPaths(std::string input_file, std::string output_file, std::string fibre, double wavelength, bool verbose){

    // Initialise variables and histograms
    int pmtcount = 0;
    double locality = 10.0; //lpc sensitivity
    double energy = RAT::util::WavelengthToEnergy(wavelength); //FIXME: could be input argument as easy to forget

    TFile *rootfile = new TFile(output_file.c_str(),"RECREATE");

    if (verbose) {std::cout << "Initialising RAT" << std::endl;}

    // Initialise RAT
    RAT::DU::DSReader dsreader(input_file);
    RAT::DU::GroupVelocity groupVelocity = RAT::DU::Utility::Get()->GetGroupVelocity();
    RAT::DU::LightPathCalculator lightPath = RAT::DU::Utility::Get()->GetLightPathCalculator();
    const RAT::DU::PMTCalStatus& PMTCalStatus = RAT::DU::Utility::Get()->GetPMTCalStatus();  // Needed for GetHitStatus() later
    //RAT::DU::TimeResidualCalculator fTRCalc = RAT::DU::Utility::Get()->GetTimeResidualCalculator();

    lightPath.SetELLIEEvent(true); // event originates outside AV (see PR #2621)
    const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
    const int NPMTS = pmtinfo.GetCount();
    if (verbose) {std::cout << NPMTS << " PMTs found" << std::endl;}
    double transitTime[NPMTS];
    double bucketTime[NPMTS];
    double cosTheta[NPMTS];

    RAT::DB *db = RAT::DB::Get();
    RAT::DBLinkPtr entry = db->GetLink("FIBRE", fibre);
    TVector3 fibrePos(entry->GetD("x"), entry->GetD("y"), entry->GetD("z")); // position of fibre [mm]
    TVector3 fibreDir(entry->GetD("u"), entry->GetD("v"), entry->GetD("w")); // direction of fibre
    if (verbose) {std::cout << "RATDB: fibre " << fibre << ", pos: (" << fibrePos.X() << "," << fibrePos.Y() << "," << fibrePos.Z() << "), dir: (" << fibreDir.X() << "," << fibreDir.Y() << "," << fibreDir.Z() << ")" << std::endl;}

    for (int it = 0; it < NPMTS; it++) { //only need to calculate transit time for each PMT once, since fibre and pmt locations are fixed
        // Use normal or HQE PMTs only
        if (pmtinfo.GetType(it) != 1 && pmtinfo.GetType(it) != 7) {continue;}
        pmtcount++;

        // Get PMT information
        TVector3 pmtPos = pmtinfo.GetPosition(it);       // position [mm]
        TVector3 pmtDir = pmtinfo.GetDirection(it);      // direction

        // Calculate light path
        lightPath.CalcByPosition(fibrePos,pmtPos,energy,locality);

        // Get light travel time (as done in BiPoLikelihoodDiff.cc)
        double distInInnerAV = lightPath.GetDistInInnerAV();
        double distInAV = lightPath.GetDistInAV();
        double distInWater = lightPath.GetDistInWater();
        double tof = groupVelocity.CalcByDistance(distInInnerAV, distInAV, distInWater, energy);

        // Get light bucket time (as done in DQLaserBallProc.cc)
        TVector3 endDir = lightPath.GetIncidentVecOnPMT();        // end direction at PMT
        double thetaAtPMT = endDir.Angle(pmtDir)*180./TMath::Pi();   // incident angle with bucket face
        double timeInBucket = groupVelocity.PMTBucketTime(thetaAtPMT);   // DocDB 3138

        double cosTheta_calc = (fibrePos * pmtPos)/(fibrePos.Mag() * pmtPos.Mag());

        //fill arrays
        transitTime[it] = tof;
        bucketTime[it] = timeInBucket;
        cosTheta[it] = cosTheta_calc;
    }

    // Create histograms
    TH1D *h1DResTimeAll_raw = new TH1D("g", "Residual Hit Time, not reajusted", 1000, -50., 500.);
    TH1D *h1DResTimeAll = new TH1D("h1DResTimeAll", "Residual Hit Time", 1000, -50., 250.);
    TH2F *hPMTResTimeCosTheta = new TH2F("hPmtResTimeVsCosTheta", "title",1000, -1., 1., 1000, -50., 250.);

    size_t entryCount = dsreader.GetEntryCount(); //number of entries, want to loop over each one
    std::vector<Double_t> evPMTTimes;
    std::vector<UInt_t> pmtID;
    int index = 0;
    if (verbose) {std::cout << "No of entries in run: " << entryCount << " events" << std::endl;}
    for (size_t iEntry = 0; iEntry < entryCount; ++iEntry) {
        if (iEntry %100 == 0 and verbose) {std::cout << "Entry no " << iEntry << std::endl;}
        const RAT::DS::Entry &rDS = dsreader.GetEntry(iEntry);
        for (size_t i_ev = 0; i_ev< rDS.GetEVCount(); ++i_ev){
            if (i_ev %100 == 0 and verbose) {std::cout << "Event no " << i_ev << std::endl;}
            const RAT::DS::EV &rEV = rDS.GetEV(i_ev);
            Int_t triggerWord = rEV.GetTrigType();
            if (!(triggerWord & (1 << 15))) {continue;} // EXTA cut
            std::cout << "passed" << std::endl;

            const RAT::DS::CalPMTs &calPMTs = rEV.GetCalPMTs(); 
            size_t calPMT_count = calPMTs.GetCount();
            
            // calculate time residuals (not ajusted for peak hit time yet), to fit gaussian
            for (size_t i_evpmt = 0; i_evpmt < calPMT_count; ++i_evpmt) {
                const RAT::DS::PMTCal& pmtCal = calPMTs.GetPMT(i_evpmt);
                // Check if PMT passes data cleaning
                if (PMTCalStatus.GetHitStatus(pmtCal) != 0) {continue;}

                pmtID.push_back(pmtCal.GetID());
                //evPMTTimes.push_back(fTRCalc.CalcTimeResidual(pmtCal, fibrePos, 0.0, true, energy, true, locality));  // Can't set ELLIE event in CalcTimeResidual.
                evPMTTimes.push_back(pmtCal.GetTime() - transitTime[pmtID[index]] - bucketTime[pmtID[index]]);
                h1DResTimeAll_raw->Fill(evPMTTimes[index]);
                ++index;
            }
        }
    }

    Int_t binmax = h1DResTimeAll_raw->GetMaximumBin();
    Double_t peak_time = h1DResTimeAll_raw->GetXaxis()->GetBinCenter(binmax);

    // calculate time residuals (ajusted)
    //std::cout << "peak_time = " << peak_time << std::endl;
    for (size_t i_evpmt = 0; i_evpmt < pmtID.size(); ++i_evpmt) {
        h1DResTimeAll->Fill(evPMTTimes[i_evpmt] - peak_time);
        // cos(theta) hist
        hPMTResTimeCosTheta->Fill(cosTheta[pmtID[i_evpmt]], evPMTTimes[i_evpmt] - peak_time);
    }

    //now write everything
    rootfile->cd();
    hPMTResTimeCosTheta->Write();
    h1DResTimeAll->Write();
    h1DResTimeAll_raw->Write();
    rootfile->Write();
    rootfile->Close();
}


// /**
//  * @brief Get the Light Paths object
//  * 
//  * @param input_file Input simulation root file
//  * @param output_file Output root file
//  * @param fibre AMELLIE fibre that was fired in the simulation
//  * @param wavelength in mm I think (403E-6)
//  * @param verbose
//  */
// void GetLightPaths(std::string input_file, std::string output_file, std::string fibre, double wavelength, bool verbose){

//     // Initialise variables and histograms
//     int pmtcount = 0;
//     double locality = 10.0; // lpc sensitivity
//     //double locality = 0.0; // lpc sensitivity
//     double energy = RAT::util::WavelengthToEnergy(wavelength);

//     TFile *rootfile = new TFile(output_file.c_str(),"RECREATE");

//     if(verbose) {std::cout << "Initialising RAT" << std::endl;}
//     // Initialise RAT
//     RAT::DU::DSReader dsreader(input_file);
//     RAT::DU::GroupVelocity groupVelocity = RAT::DU::Utility::Get()->GetGroupVelocity();
//     RAT::DU::LightPathCalculator lightPath = RAT::DU::Utility::Get()->GetLightPathCalculator();
//     lightPath.SetELLIEEvent(true); // event originates outside AV (see PR #2621)
//     const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
//     const int NPMTS = pmtinfo.GetCount();
//     if(verbose) {std::cout << NPMTS << " PMTs found" << std::endl;}

//     RAT::DB *db = RAT::DB::Get();
//     RAT::DBLinkPtr entry = db->GetLink("FIBRE", fibre);
//     TVector3 fibrePos(entry->GetD("x"), entry->GetD("y"), entry->GetD("z")); // position of fibre [mm]
//     TVector3 fibreDir(entry->GetD("u"), entry->GetD("v"), entry->GetD("w")); // direction of fibre
//     if(verbose) {std::cout << "RATDB: fibre " << fibre << ", pos: (" << fibrePos.X() << "," << fibrePos.Y() << "," << fibrePos.Z() << "), dir: (" << fibreDir.X() << "," << fibreDir.Y() << "," << fibreDir.Z() << ")" << std::endl;}

//     // Create histograms
//     TH2F *hPMTResTimeCosTheta = new TH2F("hPmtResTimeVsCosTheta", "title", 1000, -1., 1., 1000, -50., 250.);
//     TH1D *h1DResTimeAll = new TH1D("h1DResTimeAll", "Hit Time (not residual, would need to subtract from peak time)", 1000, -50., 500.);

//     const RAT::DU::PMTCalStatus& PMTCalStatus = RAT::DU::Utility::Get()->GetPMTCalStatus();  // Needed for GetHitStatus() later
//     RAT::DU::TimeResidualCalculator fTRCalc = RAT::DU::Utility::Get()->GetTimeResidualCalculator();

//     size_t entryCount = dsreader.GetEntryCount(); //number of entries, want to loop over each one
//     std::vector<Double_t> evPMTTimes;
//     std::vector<UInt_t> pmtID;
//     std::vector<Double_t> cosTheta;
//     int index = 0;
//     if(verbose) {std::cout << "No of entries in run: " << entryCount << " events" << std::endl;}
//     for (size_t iEntry = 0; iEntry < entryCount; ++iEntry) {
//         if (iEntry %100 == 0 and verbose) {std::cout << "Entry no " << iEntry << std::endl;}
//         const RAT::DS::Entry &rDS = dsreader.GetEntry(iEntry);
//         for(size_t i_ev = 0; i_ev< rDS.GetEVCount(); ++i_ev) {
//             if (i_ev %100 == 0 and verbose) {std::cout << "Event no " << i_ev << std::endl;}
//             const RAT::DS::EV &rEV = rDS.GetEV(i_ev);
//             Int_t triggerWord = rEV.GetTrigType();
//             if(!(triggerWord & (1 << 15))) {continue;} // EXTA cut
//             std::cout << "passed" << std::endl;

//             const RAT::DS::CalPMTs &calPMTs = rEV.GetCalPMTs(); 
//             size_t calPMT_count = calPMTs.GetCount();
            
            
//             for (size_t i_evpmt = 0; i_evpmt < calPMT_count; ++i_evpmt) {
//                 // calculate time residuals (not ajusted for peak hit time yet), to fit gaussian
//                 const RAT::DS::PMTCal& pmtCal = calPMTs.GetPMT(i_evpmt);

//                 // Check if PMT passes data cleaning
//                 if (PMTCalStatus.GetHitStatus(pmtCal) != 0) {continue;}

//                 // Get PMT ID, angle (cos(theta)) and residual hit time
//                 pmtID.push_back(pmtCal.GetID());
//                 TVector3 pmtPos = pmtinfo.GetPosition(pmtID.at(index));       // position [mm]
//                 cosTheta.push_back((fibrePos * pmtPos)/(fibrePos.Mag() * pmtPos.Mag()));
//                 evPMTTimes.push_back(fTRCalc.CalcTimeResidual(pmtCal, fibrePos, 0.0, true, energy, false, locality));  // set event time to zero for now. Will subtract in hist later.
//                 h1DResTimeAll->Fill(evPMTTimes[index]);
//                 ++index;
//             }            
//         }
//     }

//     // Find bin with highest count -> peak
//     Int_t binmax = h1DResTimeAll->GetMaximumBin();
//     Double_t peak_time = h1DResTimeAll->GetXaxis()->GetBinCenter(binmax);

//     // calculate time residuals (ajusted)
//     //std::cout << "peak_time = " << peak_time << std::endl;
//     for (size_t i_evpmt = 0; i_evpmt < pmtID.size(); ++i_evpmt) {
//         hPMTResTimeCosTheta->Fill(cosTheta[pmtID[i_evpmt]], evPMTTimes[i_evpmt] - peak_time);
//     }

//     //now write everything and close
//     rootfile->cd();
//     h1DResTimeAll->Write();
//     hPMTResTimeCosTheta->Write();
//     rootfile->Write();
//     rootfile->Close();
// }