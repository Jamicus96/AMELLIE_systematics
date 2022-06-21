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
    double wavelength = std::stod(argv[4]);
    bool verbose = std::stoi(argv[5]);
    GetLightPaths(input_file, output_file, fibre, wavelength, verbose);
    return 0;
}


/**
 * @brief Get the Light Paths object
 * 
 * @param file Input simulation root file
 * @param fibre AMELLIE fibre that was fired in the simulation
 * @param wavelength in mm I think (403E-6)
 */
void GetLightPaths(std::string input_file, std::string output_file, std::string fibre, double wavelength, bool verbose){

    // Initialise variables and histograms
    int pmtcount = 0;
    //double locality = 10.0; // lpc sensitivity
    double locality = 0.0; // lpc sensitivity
    double energy = RAT::util::WavelengthToEnergy(wavelength);

    TFile *rootfile = new TFile(output_file.c_str(),"RECREATE");

    if(verbose) {std::cout << "Initialising RAT" << std::endl;}
    // Initialise RAT
    RAT::DU::DSReader dsreader(input_file);
    RAT::DU::GroupVelocity groupVelocity = RAT::DU::Utility::Get()->GetGroupVelocity();
    RAT::DU::LightPathCalculator lightPath = RAT::DU::Utility::Get()->GetLightPathCalculator();
    lightPath.SetELLIEEvent(true); // event originates outside AV (see PR #2621)
    const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
    const int NPMTS = pmtinfo.GetCount();
    if(verbose) {std::cout << NPMTS << " PMTs found" << std::endl;}

    RAT::DB *db = RAT::DB::Get();
    RAT::DBLinkPtr entry = db->GetLink("FIBRE", fibre);
    TVector3 fibrePos(entry->GetD("x"), entry->GetD("y"), entry->GetD("z")); // position of fibre [mm]
    TVector3 fibreDir(entry->GetD("u"), entry->GetD("v"), entry->GetD("w")); // direction of fibre
    if(verbose) {std::cout << "RATDB: fibre " << fibre << ", pos: (" << fibrePos.X() << "," << fibrePos.Y() << "," << fibrePos.Z() << "), dir: (" << fibreDir.X() << "," << fibreDir.Y() << "," << fibreDir.Z() << ")" << std::endl;}

    // Create histograms
    TH2F *hPMTResTimeCosTheta = new TH2F("hPmtResTimeVsCosTheta", "title", 1000, -1., 1., 1000, -50., 250.);
    TH1D *h1DResTimeAll = new TH1D("h1DResTimeAll", "Residual Hit Time", 1000, -50., 250.);

    const RAT::DU::PMTCalStatus& PMTCalStatus = RAT::DU::Utility::Get()->GetPMTCalStatus();  // Needed for GetHitStatus() later
    RAT::DU::TimeResidualCalculator fTRCalc = RAT::DU::Utility::Get()->GetTimeResidualCalculator();

    size_t entryCount = dsreader.GetEntryCount(); //number of entries, want to loop over each one
    std::vector<Double_t> evPMTTimes;
    std::vector<UInt_t> pmtID;
    std::vector<Double_t> cosTheta;
    int index = 0;
    if(verbose) {std::cout << "No of entries in run: " << entryCount << " events" << std::endl;}
    for (size_t iEntry = 0; iEntry < entryCount; ++iEntry) {
        if (iEntry %100 == 0 and verbose) {std::cout << "Entry no " << iEntry << std::endl;}
        const RAT::DS::Entry &rDS = dsreader.GetEntry(iEntry);
        for(size_t i_ev = 0; i_ev< rDS.GetEVCount(); ++i_ev) {
            if (i_ev %100 == 0 and verbose) {std::cout << "Event no " << i_ev << std::endl;}
            const RAT::DS::EV &rEV = rDS.GetEV(i_ev);
            Int_t triggerWord = rEV.GetTrigType();
            if(!(triggerWord & (1 << 15))) {continue;} // EXTA cut
            std::cout << "passed" << std::endl;

            const RAT::DS::CalPMTs &calPMTs = rEV.GetCalPMTs(); 
            size_t calPMT_count = calPMTs.GetCount();
            
            
            for (size_t i_evpmt = 0; i_evpmt < calPMT_count; ++i_evpmt) {
                // calculate time residuals (not ajusted for peak hit time yet), to fit gaussian
                const RAT::DS::PMTCal& pmtCal = calPMTs.GetPMT(i_evpmt);

                // Check if PMT passes data cleaning
                if (PMTCalStatus.GetHitStatus(pmtCal) != 0) {continue;}

                // Get PMT ID, angle (cos(theta)) and residual hit time
                pmtID.push_back(pmtCal.GetID());
                TVector3 pmtPos = pmtinfo.GetPosition(pmtID.at(index));       // position [mm]
                cosTheta.push_back((fibrePos * pmtPos)/(fibrePos.Mag() * pmtPos.Mag()));
                evPMTTimes.push_back(fTRCalc.CalcTimeResidual(pmtCal, fibrePos, 0.0, true, energy, false, locality));  // set event time to zero for now. Will subtract in hist later.
                h1DResTimeAll->Fill(evPMTTimes[index]);
                ++index;
            }            
        }
    }

    // Find bin with highest count -> peak
    Int_t binmax = h1DResTimeAll->GetMaximumBin();
    Double_t peak_time = h1DResTimeAll->GetXaxis()->GetBinCenter(binmax);

    // calculate time residuals (ajusted)
    //std::cout << "peak_time = " << peak_time << std::endl;
    for (size_t i_evpmt = 0; i_evpmt < pmtID.size(); ++i_evpmt) {
        hPMTResTimeCosTheta->Fill(cosTheta[pmtID[i_evpmt]], evPMTTimes[i_evpmt] - peak_time);
    }

    //now write everything and close
    rootfile->cd();
    h1DResTimeAll->Write();
    hPMTResTimeCosTheta->Write();
    rootfile->Write();
    rootfile->Close();
}