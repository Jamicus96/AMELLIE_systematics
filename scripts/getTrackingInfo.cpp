#include <iostream>
#include <string>
#include <TCanvas.h>
#include "RAT/DU/DSReader.hh"
#include "RAT/DU/PMTInfo.hh"
#include "RAT/DU/Utility.hh"
#include "RAT/DS/Run.hh"
#include "RAT/DS/Entry.hh"
#include "RAT/DS/MC.hh"
#include <RAT/PhysicsUtil.hh>
#include <TH2F.h>
#include <TH1F.h>
#include <TF1.h>
#include <TMath.h>
#include <TFile.h>
#include <RAT/TrackNav.hh>
#include <RAT/TrackCursor.hh>
#include <RAT/TrackNode.hh>
#include <RAT/DS/MCPhoton.hh>
#include <TVector.h>
#include <TGraph2D.h>
#include <regex>

/*
Photon processes:
0 - Single scatter in scint
1 - Single absorption and re-emission in scint
2 - No effects
3 - Multiple effect
4 - Unknown effect
5 - Near reflection
6 - Ropes
7 - PMT reflection
8 - single scatter in external water
9 - inner av reflection
10 - AV pipes
11 - acrylic scatter
12 - other scattering
*/

std::vector<std::string> ESummaryFlag = {"OpScintillation", "OpReemission", "OpCerenkov", "OpAbsorption", "OpRayleigh", "OpMie", "HitConc", "HitPMT", "Compton", "OpRayleighH2O", "OpRayleighAV", "OpRayleighInnerAV", "OpMieH2O", "OpMieAV", "OpMieInnerAV", "CoulombScatter"};

std::vector<RAT::DS::MCTrack::ESummaryFlag> flags = {RAT::DS::MCTrack::OpScintillation ,RAT::DS::MCTrack::OpReemission, RAT::DS::MCTrack::OpCerenkov, RAT::DS::MCTrack::OpAbsorption, RAT::DS::MCTrack::OpRayleigh, RAT::DS::MCTrack::OpMie, RAT::DS::MCTrack::HitConc, RAT::DS::MCTrack::HitPMT, RAT::DS::MCTrack::Compton, RAT::DS::MCTrack::OpRayleighH2O, RAT::DS::MCTrack::OpRayleighAV, RAT::DS::MCTrack::OpRayleighInnerAV, RAT::DS::MCTrack::OpMieH2O, RAT::DS::MCTrack::OpMieAV, RAT::DS::MCTrack::OpMieInnerAV, RAT::DS::MCTrack::CoulombScatter};

int GetPhotonProcess(const RAT::DS::MCTrack &inputPhotonTrack, const RAT::DS::MC &inputMC, int event_num, bool debug);

int GetLightPaths(std::string file, std::string fibre, std::string data_type);

void WriteProcess(const RAT::DS::MCTrack &inputPhotonTrack);

double ReflectionAngle(const TVector3& eventPos);

TVector2 IcosProject(TVector3 pmtPos);

TVector2 TransformCoord( const TVector3& V1, const TVector3& V2, const TVector3& V3, const TVector2& A1, const TVector2& A2, const TVector2& A3,const TVector3& P);

int GetPMTID(std::string input);

int main(int argc, char** argv){
    std::string file = argv[1];
    std::string fibre = argv[2];
    std::string data_type = argv[3];  // MC or raw
    int returnCode = GetLightPaths(file, fibre, data_type);
    return 0;
}


double ReflectionAngle(const TVector3& eventPos){
    double angle = 2.0 * ( TMath::ACos(6060.418 / eventPos.Mag()));
    return angle;
}


TVector2 TransformCoord( const TVector3& V1, const TVector3& V2, const TVector3& V3, const TVector2& A1, const TVector2& A2, const TVector2& A3,const TVector3& P ) {
    TVector3 xV = V2 - V1;
    TVector3 yV = ( ( V3 - V1 ) + ( V3 - V2 ) ) * 0.5;
    TVector3 zV = xV.Cross( yV ).Unit();

    double planeD = V1.Dot( zV );

    double t = planeD / P.Dot( zV );

    TVector3 localP = t*P - V1;

    TVector2 xA = A2 - A1;
    TVector2 yA = ( ( A3 - A1 ) +( A3 - A2 ) ) * 0.5;

    double convUnits = xA.Mod() / xV.Mag();

    TVector2 result;
    result = localP.Dot( xV.Unit() ) * xA.Unit() * convUnits;
    result += localP.Dot( yV.Unit() ) * yA.Unit() * convUnits + A1;
    return result;
}


TVector2 IcosProject(TVector3 pmtPos) {

    double fa = 1.0 / 5.5;
    double fb = fa * sqrt( 3.0 ) / 2.0;
    TVector2 *fA12a = new TVector2( fa / 2.0, 0.0 );
    TVector2 *fA12b = new TVector2( 3.0 * fa / 2.0, 0.0 );
    TVector2 *fA12c = new TVector2( 5.0 * fa / 2.0, 0.0 );
    TVector2 *fA12d = new TVector2( 7.0 *fa / 2.0, 0.0 );
    TVector2 *fA12e = new TVector2( 9.0 * fa / 2.0, 0.0 );
    TVector2 *fA2a = new TVector2( 0.0, fb );
    TVector2 *fA2b = new TVector2( 5.0 * fa, fb );
    TVector2 *fA17a = new TVector2( fa / 2.0 , 2.0 * fb );
    TVector2 *fA17b = new TVector2( 11.0 * fa / 2.0 , 2.0 * fb );
    TVector2 *fA51a = new TVector2( fa, 3.0 * fb );
    TVector2 *fA51b = new TVector2( 2.0 * fa, 3.0 * fb );
    TVector2 *fA51c = new TVector2( 3.0 * fa, 3.0 * fb );
    TVector2 *fA51d = new TVector2( 4.0 * fa, 3.0 * fb );
    TVector2 *fA51e = new TVector2( 5.0 * fa, 3.0 * fb );
    TVector2 *fA27 = new TVector2( 4.0 * fa, fb );
    TVector2 *fA46 = new TVector2( 3.0 * fa, fb );
    TVector2 *fA31 = new TVector2( 2.0 * fa, fb );
    TVector2 *fA6 = new TVector2( fa, fb );
    TVector2 *fA37 = new TVector2( 9.0 * fa / 2.0 , 2.0 * fb );
    TVector2 *fA33 = new TVector2( 3.0 * fa / 2.0 , 2.0 * fb );
    TVector2 *fA58 = new TVector2( 5.0 * fa / 2.0 , 2.0 * fb );
    TVector2 *fA54 = new TVector2( 7.0 * fa / 2.0 , 2.0 * fb );

    TVector3 pointOnSphere( pmtPos.X(), pmtPos.Y(), pmtPos.Z() );
    pointOnSphere = pointOnSphere.Unit();
    pointOnSphere.RotateX( -45.0 );
    // From http://www.rwgrayprojects.com/rbfnotes/polyhed/PolyhedraData/Icosahedralsahedron/Icosahedralsahedron.pdf
    const double t = ( 1.0 + sqrt( 5.0 ) ) / 2.0;
    const TVector3 V2 = TVector3( t * t, 0.0, t * t * t ).Unit();
    const TVector3 V6 = TVector3( -t * t, 0.0, t * t * t ).Unit();
    const TVector3 V12 = TVector3( 0.0, t * t * t, t * t ).Unit();
    const TVector3 V17 = TVector3( 0.0, -t * t * t, t * t ).Unit();
    const TVector3 V27 = TVector3( t * t * t, t * t, 0.0 ).Unit();
    const TVector3 V31 = TVector3( -t * t * t, t * t, 0.0 ).Unit();
    const TVector3 V33 = TVector3( -t * t * t, -t * t, 0.0 ).Unit();
    const TVector3 V37 = TVector3( t * t * t, -t * t, 0.0 ).Unit();
    const TVector3 V46 = TVector3( 0.0, t * t * t, -t * t ).Unit();
    const TVector3 V51 = TVector3( 0.0, -t * t * t, -t * t ).Unit();
    const TVector3 V54 = TVector3( t * t, 0.0, -t * t * t ).Unit();
    const TVector3 V58 = TVector3( -t * t, 0.0, -t * t * t ).Unit();

    // Faces {{ 2, 6, 17}, { 2, 12, 6}, { 2, 17, 37}, { 2, 37, 27}, { 2, 27, 12}, {37, 54, 27},
    // {27, 54, 46}, {27, 46, 12}, {12, 46, 31}, {12, 31, 6}, { 6, 31, 33}, { 6, 33, 17},
    // {17, 33, 51}, {17, 51, 37}, {37, 51, 54}, {58, 54, 51}, {58, 46, 54}, {58, 31, 46},
    // {58, 33, 31}, {58, 51, 33}}

    std::vector<TVector3> IcosahedralCentres;
    IcosahedralCentres.push_back( ( V2 + V6 + V17 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V2 + V12 + V6 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V2 + V17 + V37 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V2 + V37 + V27 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V2 + V27 + V12 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V37 + V54 + V27 ) * ( 1.0 / 3.0 ) );

    IcosahedralCentres.push_back( ( V27 + V54 + V46 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V27 + V46 + V12 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V12 + V46 + V31 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V12 + V31 + V6 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V6 + V31 + V33 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V6 + V33 + V17 ) * ( 1.0 / 3.0 ) );

    IcosahedralCentres.push_back( ( V17 + V33 + V51 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V17 + V51 + V37 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V37 + V51 + V54 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V58 + V54 + V51 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V58 + V46 + V54 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V58 + V31 + V46 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V58 + V33 + V31 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V58 + V51 + V33 ) * ( 1.0 / 3.0 ) );

    std::vector<double> distFromCentre;
    unsigned int uLoop;
    for( uLoop = 0; uLoop < IcosahedralCentres.size(); uLoop++ ){
      distFromCentre.push_back( ( IcosahedralCentres[uLoop] - pointOnSphere ).Mag() );
    }
    const int face = min_element( distFromCentre.begin(), distFromCentre.end() ) - distFromCentre.begin() + 1;
    TVector2 resultPosition;
    switch(face){
    case 1://{ 2, 6, 17}
      resultPosition = TransformCoord( V2, V6, V17, *fA2a, *fA6, *fA17a, pointOnSphere );
      break;
    case 2://{ 2, 12, 6}
      resultPosition = TransformCoord( V2, V12, V6, *fA2a, *fA12a, *fA6, pointOnSphere );
      break;
    case 3://{ 2, 17, 37}
      resultPosition = TransformCoord( V2, V17, V37, *fA2b, *fA17b, *fA37, pointOnSphere );
      break;
    case 4://{ 2, 37, 27}
      resultPosition = TransformCoord( V2, V37, V27, *fA2b, *fA37, *fA27, pointOnSphere );
      break;
    case 5://{ 2, 27, 12}
      resultPosition = TransformCoord( V2, V27, V12, *fA2b, *fA27, *fA12e, pointOnSphere );
      break;
    case 6://{37, 54, 27}
      resultPosition = TransformCoord( V37, V54, V27, *fA37, *fA54, *fA27, pointOnSphere );
      break;
    case 7://{27, 54, 46}
      resultPosition = TransformCoord( V27, V54, V46, *fA27, *fA54, *fA46, pointOnSphere );
      break;
    case 8://{27, 46, 12}
      resultPosition = TransformCoord( V27, V46, V12, *fA27, *fA46, *fA12d, pointOnSphere );
      break;
    case 9://{12, 46, 31}
      resultPosition = TransformCoord( V12, V46, V31, *fA12c, *fA46, *fA31, pointOnSphere );
      break;
    case 10://{12, 31, 6}
      resultPosition = TransformCoord( V12, V31, V6, *fA12b, *fA31, *fA6, pointOnSphere );
      break;
    case 11://{ 6, 31, 33}
      resultPosition = TransformCoord( V6, V31, V33, *fA6, *fA31, *fA33, pointOnSphere );
      break;
    case 12://{ 6, 33, 17}
      resultPosition = TransformCoord( V6, V33, V17, *fA6, *fA33, *fA17a, pointOnSphere );
      break;
    case 13://{17, 33, 51}
      resultPosition = TransformCoord( V17, V33, V51, *fA17a, *fA33, *fA51a, pointOnSphere );
      break;
    case 14://{17, 51, 37}
      resultPosition = TransformCoord( V17, V51, V37, *fA17b, *fA51e, *fA37, pointOnSphere );
      break;
    case 15://{37, 51, 54}
      resultPosition = TransformCoord( V37, V51, V54, *fA37, *fA51d, *fA54, pointOnSphere );
      break;
    case 16://{58, 54, 51}
      resultPosition = TransformCoord( V58, V54, V51, *fA58, *fA54, *fA51c, pointOnSphere );
      break;
    case 17://{58, 46, 54}
      resultPosition = TransformCoord( V58, V46, V54, *fA58, *fA46, *fA54, pointOnSphere );
      break;
    case 18://{58, 31, 46}
      resultPosition = TransformCoord( V58, V31, V46, *fA58, *fA31, *fA46, pointOnSphere );
      break;
    case 19://{58, 33, 31}
      resultPosition = TransformCoord( V58, V33, V31, *fA58, *fA33, *fA31, pointOnSphere );
      break;
    case 20://{58, 51, 33}
      resultPosition = TransformCoord( V58, V51, V33, *fA58, *fA51b, *fA33, pointOnSphere );
      break;
    }
    // output face, if needed
    //segment = face;
    // 1 - x,y pos to project the same as node map
    return TVector2(1.0 - resultPosition.X(), 1.0 - 2.0 * resultPosition.Y() );
}


int GetPhotonProcess(const RAT::DS::MCTrack &inputPhotonTrack, const RAT::DS::MC &inputMC, int event_num, bool debug){

    //make flags

    int num_abs = 0;
    int num_av_scat = 0;
    int num_other = 0;
    int num_psup_scat = 0;
    int num_near_reflect = 0;
    int num_inner_av_reflect = 0;
    int num_PMT_reflect = 0;
    int num_abs_other = 0;
    int num_scat_other = 0;
    int num_ropes = 0;
    int num_extwater_scat = 0;
    int num_av_pipes = 0;
    int num_acrylic_scatter = 0;
    int num_outer_av_reflect = 0;
    int sum_effects = 0;

    bool abs_phot = false;

    TVector3 prev_pos;
    TVector3 first_pos;
    std::string prev_end_vol = "";
    std::string prev_start_vol = "";

    RAT::DS::MCTrack::ESummaryFlag flag = RAT::DS::MCTrack::OpReemission;

    if(debug){ //probs better way of doing this
        for(int i=0; i<16; i++){
            RAT::DS::MCTrack::ESummaryFlag flag = flags.at(i);
            if(inputPhotonTrack.GetSummaryFlag(flag)){
                std::cout << "Track underwent: " << ESummaryFlag.at(i) << std::endl;
            }
        }
    }

    //loop through

    for(int i_step = 0; i_step < inputPhotonTrack.GetMCTrackStepCount();i_step++){
        const RAT::DS::MCTrackStep &rTrackStep = inputPhotonTrack.GetMCTrackStep(i_step);

        if(debug) std::cout << "Step: " << i_step << ", Start vol: " << rTrackStep.GetStartVolume() << ", End Vol: " << rTrackStep.GetEndVolume() << ", process: " << rTrackStep.GetProcess() << ", Start Pos: " << rTrackStep.GetPosition().X() << ", " << rTrackStep.GetPosition().Y() << ", " << rTrackStep.GetPosition().Z() << ", time: " << rTrackStep.GetGlobalTime() << std::endl;

        if(i_step == 0 and rTrackStep.GetStartVolume() != "cavity"){ //reemission - I don't think there's another case where the photon is not in cavity at step 0
            num_abs++;
            sum_effects++;
            abs_phot = true;
        }

        if(rTrackStep.GetEndVolume() == "snorope" or rTrackStep.GetEndVolume() == "internal_ropes_3"){ //photon hits rope
            num_ropes++;
            sum_effects++;
        }
        if(rTrackStep.GetEndVolume() == "world"){ //probably should never happen
            num_other++;
            sum_effects++;
        }

        if(rTrackStep.GetEndVolume().find("av_pipes") != std::string::npos){ //photon hits AV pipe
            num_av_pipes++;
            sum_effects++;
        }

        if(rTrackStep.GetStartVolume() == "av" and rTrackStep.GetEndVolume() == "cavity" and prev_start_vol == "cavity" and prev_end_vol == "av"){
            if(rTrackStep.GetPosition().Angle(first_pos) * (180/TMath::Pi()) > 90.){ //far AV - since angle is defined as seen from center of detector
                num_outer_av_reflect++;
                sum_effects++;
            }
            else{ //near AV
                num_near_reflect++;
                sum_effects++;
            }
        }

        if((rTrackStep.GetStartVolume() == "av" and rTrackStep.GetEndVolume() == "inner_av") and ((sum_effects == 0 and i_step != 2) or (sum_effects != 0))){ //FIXME: This might not be robust, e.g if external water scatter before entering acrylic
            if(prev_start_vol != "cavity"){
                num_inner_av_reflect++;
                sum_effects++;
            }
            else{
                const RAT::DS::MCTrackStep &prevPrevTrackStep = inputPhotonTrack.GetMCTrackStep(i_step-2); //don't want to include light that returns into the AV
                if(prevPrevTrackStep.GetStartVolume() == "av"){
                    num_inner_av_reflect++;
                    sum_effects++;
                }
            }
        }

        if(rTrackStep.GetProcess() == "OpRayleigh"){ //scattering
            if(rTrackStep.GetStartVolume() == "inner_av" and rTrackStep.GetEndVolume() == "inner_av"){
                num_av_scat++;
                sum_effects++;
            }
            else if(rTrackStep.GetStartVolume() == "cavity" and rTrackStep.GetEndVolume() == "cavity"){
                num_extwater_scat++;
                sum_effects++;
            }
            else if(rTrackStep.GetStartVolume() == "av" and rTrackStep.GetEndVolume() == "av"){
                num_acrylic_scatter++;
                sum_effects++;
            }
            else{
                num_scat_other++;
                sum_effects++;
            }
        }
        if(rTrackStep.GetProcess() == "OpAbsorption"){ //absorption - this might be redundant since we are looking backwards and so we see the reemitted photon first
            if(rTrackStep.GetStartVolume() == "inner_av" and rTrackStep.GetEndVolume() == "inner_av"){
                num_abs++;
                sum_effects++;
            }
            else{ //shouldn't have absorption that isn't in the scintillator
                num_abs_other++;
                sum_effects++;
            }
        }

        if((rTrackStep.GetEndVolume() == "cavity" and rTrackStep.GetStartVolume().find("innerPMT") != std::string::npos) or (rTrackStep.GetStartVolume() == "cavity" and prev_end_vol.find("innerPMT") != std::string::npos)){ //3d pmt or grey disk
            num_PMT_reflect++;
            sum_effects++;
        }

        if(i_step != 0){
            if(debug) std::cout << "Angle between this direction and last direction is " << rTrackStep.GetPosition().Angle(prev_pos) * (180/TMath::Pi()) << std::endl;
        }
        else{
            first_pos = rTrackStep.GetPosition();
        }
        prev_pos = rTrackStep.GetPosition(); //keep record of current step to use in next step logic
        prev_end_vol = rTrackStep.GetEndVolume();
        prev_start_vol = rTrackStep.GetStartVolume();
    }

    if(abs_phot and num_other == 0 and num_av_scat == 0 and num_abs == 1 and num_near_reflect == 0 and num_ropes == 0 and num_av_pipes == 0 and num_extwater_scat == 0){ //might need to check parent track
        //find parent - FIXME: Multiple absorptions?
        UInt_t rparentTrackID = inputPhotonTrack.GetParentID();
        const RAT::DS::MCTrack &parentPhotonTrack = inputMC.GetMCTrack(rparentTrackID);

        for(int i_step = 0; i_step < parentPhotonTrack.GetMCTrackStepCount();i_step++){
            const RAT::DS::MCTrackStep &rParentTrackStep = parentPhotonTrack.GetMCTrackStep(i_step);
            if(debug) std::cout << "Step: " << i_step << ", Start vol: " << rParentTrackStep.GetStartVolume() << ", End Vol: " << rParentTrackStep.GetEndVolume() << ", process: " << rParentTrackStep.GetProcess() << std::endl;

            if(i_step == 0 and rParentTrackStep.GetStartVolume() != "cavity"){ //absorption
                num_abs++;
                sum_effects++;
            }

            if(rParentTrackStep.GetEndVolume() == "snorope"){
                num_ropes++;
                sum_effects++;
            }

            if(rParentTrackStep.GetEndVolume().find("innerPMT") != std::string::npos){
                num_av_pipes++;
                sum_effects++;
            }


            if(rParentTrackStep.GetProcess() == "OpRayleigh"){ //scattering
                if(rParentTrackStep.GetStartVolume() == "inner_av" and rParentTrackStep.GetEndVolume() == "inner_av"){
                    num_av_scat++;
                    sum_effects++;
                }
                else if(rParentTrackStep.GetStartVolume() == "cavity" and rParentTrackStep.GetEndVolume() == "cavity"){
                    num_extwater_scat++;
                    sum_effects++;
                }
                else if(rParentTrackStep.GetStartVolume() == "av" and rParentTrackStep.GetEndVolume() == "av"){
                    num_acrylic_scatter++;
                    sum_effects++;
                }
                else{
                    num_scat_other++;
                    sum_effects++;
                }
            }
            if(rParentTrackStep.GetProcess() == "OpAbsorption"){ //scattering
                if(i_step != parentPhotonTrack.GetMCTrackStepCount() - 1){
                    if(rParentTrackStep.GetStartVolume() == "inner_av" and rParentTrackStep.GetEndVolume() == "inner_av"){
                        num_abs++;
                        sum_effects++;
                    }
                    else{
                        num_abs_other++;
                        sum_effects++;
                    }
                }
            }
        }

    }

    if(debug) std::cout << "Angle between first pos and end pos is " << prev_pos.Angle(first_pos) * (180/TMath::Pi()) << std::endl;
    //now decide what the process was and what to return


    if(num_abs == 0 and num_av_scat == 0 and num_other == 0 and num_near_reflect == 0 and num_ropes == 0 and num_PMT_reflect == 0 and num_scat_other == 0 and num_inner_av_reflect == 0 and num_av_pipes == 0 and num_acrylic_scatter == 0 and num_extwater_scat == 0){
        if(debug){
            std:: cout << "Straight through transmission" << std::endl;
            std::cout << "" << std::endl;
        }
        return 2;
    }
    else if(num_abs == 1 and num_av_scat == 0 and num_other == 0 and num_near_reflect == 0 and num_PMT_reflect == 0 and num_inner_av_reflect == 0 and num_av_pipes == 0 and num_acrylic_scatter == 0 and num_extwater_scat == 0){
        if(debug){
            std:: cout << "Absorbed and reemitted" << std::endl;
            std::cout << "" << std::endl;
        }
        return 1;
    }
    else if(num_abs == 0 and num_av_scat == 1 and num_other == 0 and num_near_reflect == 0 and num_PMT_reflect == 0 and num_inner_av_reflect == 0 and num_av_pipes == 0 and num_acrylic_scatter == 0 and num_extwater_scat == 0){
        if(debug){
            std:: cout << "Single scatter in AV" << std::endl;
            std::cout << "" << std::endl;
        }
        return 0;
    }
    else if(num_abs == 0 and num_av_scat == 0 and num_other == 0 and num_near_reflect == 1 and num_PMT_reflect == 0 and num_inner_av_reflect == 0 and num_av_pipes == 0 and num_acrylic_scatter == 0 and num_extwater_scat == 0){
        if(debug){
            std:: cout << "near reflect" << std::endl;
            std::cout << "" << std::endl;
        }
        return 5;
    }
    else if(num_abs > 0 and num_av_scat > 0){
        if(debug){
            std::cout << "Absorbed and scattered" << std::endl;
            std::cout << "" << std::endl;
        }
        return 3;
    }
    else if(num_av_scat > 0 and num_near_reflect > 0 and num_PMT_reflect == 0){
        if(debug){
            std::cout << "Near reflection, scattering" << std::endl;
            std::cout << "" << std::endl;
        }
        return 3;
    }
    else if(num_abs == 0 and num_av_scat == 0 and num_other == 0 and num_near_reflect == 0 and num_PMT_reflect == 1 and num_inner_av_reflect == 0 and num_av_pipes == 0 and num_acrylic_scatter == 0 and num_extwater_scat == 0 and num_ropes == 0){
        if(num_outer_av_reflect == 1){
            if(debug){
                std::cout << "PMT/AV reflection" << std::endl;
                std::cout << "" << std::endl;
            }
            return 13;
        }
        else{
            if(debug){
                std::cout << "PMT reflection" << std::endl;
                std::cout << "" << std::endl;
            }
            return 7;
        }
    }
    else if(num_abs == 0 and num_av_scat == 0 and num_other == 0 and num_near_reflect == 0 and num_PMT_reflect == 0 and num_inner_av_reflect == 1 and num_av_pipes == 0 and num_acrylic_scatter == 0 and num_extwater_scat == 0){
        if(debug){
            std::cout << "inner av reflection" << std::endl;
            std::cout << "" << std::endl;
        }
        return 9;
    }
    else if((num_abs > 1 or num_av_scat > 1 or num_near_reflect > 1) or sum_effects > 1){
        if(num_other == 0){
            if(debug){
                std:: cout << "Multiple effect" << std::endl;
                std::cout << "" << std::endl;
            }
            return 3;
        }
        else{
            // if(debug){
                std:: cout << "Unknown - num_abs: " << num_abs << ", num_av_scat: " << num_av_scat << ", num_other: " << num_other << ", num_near_reflect: " << num_near_reflect << ", num pmt reflect: " << num_PMT_reflect << ", sum effects: " << sum_effects << ", event number: " << event_num << std::endl;
                std::cout << "" << std::endl;
            // }
            return 4;
        }
    }

    else if(num_abs == 0 and num_av_scat == 0 and num_other == 0 and num_near_reflect == 0 and num_ropes > 0  and num_av_pipes == 0  and num_acrylic_scatter == 0){
        if(debug){
            std::cout << "ropes" << std::endl;
            std::cout << "" << std::endl;
        }
        return 6;
    }
    else if(num_abs == 0 and num_av_scat == 0 and num_other == 0 and num_near_reflect == 0 and num_ropes == 0 and num_extwater_scat == 1 and num_av_pipes == 0 and num_acrylic_scatter == 0){
        if(debug){
            std::cout << "external water scatter" << std::endl;
            std::cout << "" << std::endl;
        }
        return 8;
    }  
    else if(num_abs == 0 and num_av_scat == 0 and num_other == 0 and num_near_reflect == 0 and num_ropes == 0 and num_extwater_scat == 0 and num_av_pipes == 1 and num_acrylic_scatter == 0){
        if(debug){
            std::cout << "av pipes" << std::endl;
            std::cout << "" << std::endl;
        }
        return 10;
    }
    else if(num_abs == 0 and num_av_scat == 0 and num_other == 0 and num_near_reflect == 0 and num_ropes == 0 and num_extwater_scat == 0 and num_av_pipes == 0 and num_acrylic_scatter == 1){
        if(debug){
            std::cout << "acrylic scatter" << std::endl;
            std::cout << "" << std::endl;
        }
        return 11;
    }
    else if(num_abs == 0 and num_av_scat == 0 and num_other == 0 and num_near_reflect == 0 and num_ropes == 0 and num_extwater_scat == 0 and num_av_pipes == 0 and num_acrylic_scatter == 0 and num_scat_other == 1){
        if(debug){
            std::cout << "other scatter" << std::endl;
            std::cout << "" << std::endl;
        }
        return 12;
    }              
    else{
        // if(debug){
            std:: cout << "Unknown - num_abs: " << num_abs << ", num_av_scat: " << num_av_scat << ", num_other: " << num_other << ", num_near_reflect: " << num_near_reflect << ", num pmt reflect: " << num_PMT_reflect << ", num acrylic scatter: " << num_acrylic_scatter << ", av pipes: " << num_av_pipes << ", inner av reflect: " << num_inner_av_reflect << ", other absorption: " << num_abs_other << ", other scattering: " << num_scat_other << ", sum effects: " << sum_effects <<  ", event number: " << event_num << std::endl;
            std::cout << "" << std::endl;
        // }
        return 4;
    }
}

int GetMultipleEffect(const RAT::DS::MCTrack &inputPhotonTrack, const RAT::DS::MC &inputMC){

    //make flags
    bool debug = false;

    int num_abs = 0;
    int num_av_scat = 0;
    int num_other = 0;
    int num_psup_scat = 0;
    int num_near_reflect = 0;
    int num_inner_av_reflect = 0;
    int num_acrylic_bounce = 0;
    int num_PMT_reflect = 0;
    int num_abs_other = 0;
    int num_scat_other = 0;
    int num_ropes = 0;
    int num_extwater_scat = 0;
    int num_av_pipes = 0;
    int num_acrylic_scatter = 0;
    int num_outer_av_reflect = 0;
    int sum_effects = 0;

    bool abs_phot = false;

    TVector3 prev_pos;
    TVector3 first_pos;
    std::string prev_end_vol = "";
    std::string prev_start_vol = "";

    RAT::DS::MCTrack::ESummaryFlag flag = RAT::DS::MCTrack::OpReemission;

    for(int i_step = 0; i_step < inputPhotonTrack.GetMCTrackStepCount();i_step++){
        const RAT::DS::MCTrackStep &rTrackStep = inputPhotonTrack.GetMCTrackStep(i_step);

        if(i_step == 0 and rTrackStep.GetStartVolume() != "cavity"){ //reemission - don't double count
            num_abs++;
            sum_effects++;
            abs_phot = true;
        }

        if(rTrackStep.GetEndVolume() == "snorope" or rTrackStep.GetEndVolume() == "internal_ropes_3"){ //scatter off rope
            num_ropes++;
            sum_effects++;
        }
        if(rTrackStep.GetEndVolume() == "world"){ 
            num_other++;
            sum_effects++;
        }
        if(rTrackStep.GetEndVolume().find("av_pipes") != std::string::npos){
            num_av_pipes++;
            sum_effects++;
        }

        if(rTrackStep.GetStartVolume() == "av" and rTrackStep.GetEndVolume() == "cavity" and prev_start_vol == "cavity" and prev_end_vol == "av"){
            const RAT::DS::MCTrackStep &prevPrevTrackStep = inputPhotonTrack.GetMCTrackStep(i_step-2);
            if(prevPrevTrackStep.GetStartVolume() != "av"){    
                if(rTrackStep.GetPosition().Angle(first_pos) * (180/TMath::Pi()) > 90.){ //far AV
                    num_outer_av_reflect++;
                    sum_effects++;
                }
                else{ //near AV
                    num_near_reflect++;
                    sum_effects++;
                }
            }
            else{
                num_acrylic_bounce++;
                sum_effects++;
            }
        }

        if((rTrackStep.GetStartVolume() == "av" and rTrackStep.GetEndVolume() == "inner_av") and ((sum_effects == 0 and i_step != 2) or (sum_effects != 0))){
            if(prev_start_vol != "cavity"){
                num_inner_av_reflect++;
                sum_effects++;
            }
            else{
                const RAT::DS::MCTrackStep &prevPrevTrackStep = inputPhotonTrack.GetMCTrackStep(i_step-2);
                if(prevPrevTrackStep.GetStartVolume() == "av"){
                    num_inner_av_reflect++;
                    sum_effects++;
                }
            }
        }

        if(rTrackStep.GetProcess() == "OpRayleigh"){ //scattering
            if(rTrackStep.GetStartVolume() == "inner_av" and rTrackStep.GetEndVolume() == "inner_av"){
                num_av_scat++;
                sum_effects++;
            }
            else if(rTrackStep.GetStartVolume() == "cavity" and rTrackStep.GetEndVolume() == "cavity"){
                num_extwater_scat++;
                sum_effects++;
            }
            else if(rTrackStep.GetStartVolume() == "av" and rTrackStep.GetEndVolume() == "av"){
                num_acrylic_scatter++;
                sum_effects++;
            }
            else{
                num_scat_other++;
                sum_effects++;
            }
        }
        if(rTrackStep.GetProcess() == "OpAbsorption"){ //absorption
            if(rTrackStep.GetStartVolume() == "inner_av" and rTrackStep.GetEndVolume() == "inner_av"){
                num_abs++;
                sum_effects++;
            }
            else{
                num_abs_other++;
                sum_effects++;
            }
        }

        if((rTrackStep.GetEndVolume() == "cavity" and rTrackStep.GetStartVolume().find("innerPMT") != std::string::npos) or (rTrackStep.GetStartVolume() == "cavity" and prev_end_vol.find("innerPMT") != std::string::npos)){ //3d pmt or grey disk
            num_PMT_reflect++;
            sum_effects++;
        }

        if(i_step == 0){
            first_pos = rTrackStep.GetPosition();
        }
        prev_pos = rTrackStep.GetPosition();
        prev_end_vol = rTrackStep.GetEndVolume();
        prev_start_vol = rTrackStep.GetStartVolume();
    }

    bool parentAbsPhot = false;

    if(abs_phot and sum_effects < 3){ //might need to check parent track
        //find parent
        UInt_t rparentTrackID = inputPhotonTrack.GetParentID();
        const RAT::DS::MCTrack &parentPhotonTrack = inputMC.GetMCTrack(rparentTrackID);

        for(int i_step = 0; i_step < parentPhotonTrack.GetMCTrackStepCount();i_step++){
            const RAT::DS::MCTrackStep &rParentTrackStep = parentPhotonTrack.GetMCTrackStep(i_step);

            if(i_step == 0 and rParentTrackStep.GetStartVolume() != "cavity"){ //absorption
                num_abs++;
                sum_effects++;
                parentAbsPhot = true;
            }

            if(rParentTrackStep.GetEndVolume() == "snorope"){
                num_ropes++;
                sum_effects++;
            }

            if((rParentTrackStep.GetEndVolume() == "cavity" and rParentTrackStep.GetStartVolume().find("innerPMT") != std::string::npos) or (rParentTrackStep.GetStartVolume() == "cavity" and prev_end_vol.find("innerPMT") != std::string::npos)){
                num_PMT_reflect++;
                sum_effects++;
            }

            if(rParentTrackStep.GetEndVolume().find("av_pipes") != std::string::npos){
                num_av_pipes++;
                sum_effects++;
            }

            if(rParentTrackStep.GetProcess() == "OpRayleigh"){ //scattering
                if(rParentTrackStep.GetStartVolume() == "inner_av" and rParentTrackStep.GetEndVolume() == "inner_av"){
                    num_av_scat++;
                    sum_effects++;
                }
                else if(rParentTrackStep.GetStartVolume() == "cavity" and rParentTrackStep.GetEndVolume() == "cavity"){
                    num_extwater_scat++;
                    sum_effects++;
                }
                else if(rParentTrackStep.GetStartVolume() == "av" and rParentTrackStep.GetEndVolume() == "av"){
                    num_acrylic_scatter++;
                    sum_effects++;
                }
                else{
                    num_scat_other++;
                    sum_effects++;
                }
            }
        }
    }
    if(parentAbsPhot and sum_effects < 3){ //FIXME: Could be more than two absorptions, here I'm simply ignoring due to multiple effects
        //find parent of parent
        UInt_t rparentTrackID = inputPhotonTrack.GetParentID();
        const RAT::DS::MCTrack &parentPhotonTrack = inputMC.GetMCTrack((inputMC.GetMCTrack(rparentTrackID).GetParentID()));

        for(int i_step = 0; i_step < parentPhotonTrack.GetMCTrackStepCount();i_step++){
            const RAT::DS::MCTrackStep &rParentTrackStep = parentPhotonTrack.GetMCTrackStep(i_step);

            if(i_step == 0 and rParentTrackStep.GetStartVolume() != "cavity"){ //absorption
                num_abs++;
                sum_effects++;
            }

            if(rParentTrackStep.GetEndVolume() == "snorope"){
                num_ropes++;
                sum_effects++;
            }

            if((rParentTrackStep.GetEndVolume() == "cavity" and rParentTrackStep.GetStartVolume().find("innerPMT") != std::string::npos) or (rParentTrackStep.GetStartVolume() == "cavity" and prev_end_vol.find("innerPMT") != std::string::npos)){
                num_PMT_reflect++;
                sum_effects++;
            }

            if(rParentTrackStep.GetEndVolume().find("av_pipes") != std::string::npos){
                num_av_pipes++;
                sum_effects++;
            }

            if(rParentTrackStep.GetProcess() == "OpRayleigh"){ //scattering
                if(rParentTrackStep.GetStartVolume() == "inner_av" and rParentTrackStep.GetEndVolume() == "inner_av"){
                    num_av_scat++;
                    sum_effects++;
                }
                else if(rParentTrackStep.GetStartVolume() == "cavity" and rParentTrackStep.GetEndVolume() == "cavity"){
                    num_extwater_scat++;
                    sum_effects++;
                }
                else if(rParentTrackStep.GetStartVolume() == "av" and rParentTrackStep.GetEndVolume() == "av"){
                    num_acrylic_scatter++;
                    sum_effects++;
                }
                else{
                    num_scat_other++;
                    sum_effects++;
                }
            }
        }
    }
    
    /*
    0 - three effects or more
    1 - double absorption
    2 - double scattering
    3 - absorption and scatter
    4 - pmt reflect and absorption
    5 - pmt reflect and scattering
    6 - double pmt reflect
    7 - double inner av reflect
    8 - inner av reflect and abs
    9 - inner av reflect and scatter
    10 - double external scatter
    11 - external scatter and abs
    12 - external scatter and internal scatter
    13 - external scatter and inner av reflect
    14 - external scatter and pmt reflect
    15 - external scatter and near reflect
    16 - Near reflect and PMT reflect
    17 - PMT reflect and ropes
    18 - External scatter and ropes
    19 - Double ropes
    20 - Scattering and pipes
    21 - Double Pipes
    22 - Inner AV reflect and acylic scatter
    23 - Inner AV reflect and PMT reflect
    24 - External scatter and outer AV reflect
    25 - PMT reflect and pipes
    26 - Ropes and acrylic scatter
    27 - PMT reflect and acrylic scatter
    28 - Scatter and acrylic scatter
    29 - ropes and other scatter
    30 - scatter and acrylic bounce
    31 - external scatter and acrylic bounce
    32 - PMT reflect and acrylic bounce
    33 - other
    */

    if(sum_effects < 3){ //FIXME: Could be nicer instead of brute force checking everything
        if(num_abs == 2 and num_av_scat == 0 and num_other == 0 and num_psup_scat == 0 and num_near_reflect == 0 and num_inner_av_reflect == 0 and num_PMT_reflect == 0 and num_abs_other == 0 and num_scat_other == 0 and num_ropes == 0 and num_extwater_scat == 0 and num_av_pipes == 0 and num_acrylic_scatter == 0 and num_outer_av_reflect == 0 and num_acrylic_bounce == 0){
            return 1;
        }
        else if(num_abs == 0 and num_av_scat == 2 and num_other == 0 and num_psup_scat == 0 and num_near_reflect == 0 and num_inner_av_reflect == 0 and num_PMT_reflect == 0 and num_abs_other == 0 and num_scat_other == 0 and num_ropes == 0 and num_extwater_scat == 0 and num_av_pipes == 0 and num_acrylic_scatter == 0 and num_outer_av_reflect == 0 and num_acrylic_bounce == 0){
            return 2;
        }
        else if(num_abs == 1 and num_av_scat == 1 and num_other == 0 and num_psup_scat == 0 and num_near_reflect == 0 and num_inner_av_reflect == 0 and num_PMT_reflect == 0 and num_abs_other == 0 and num_scat_other == 0 and num_ropes == 0 and num_extwater_scat == 0 and num_av_pipes == 0 and num_acrylic_scatter == 0 and num_outer_av_reflect == 0 and num_acrylic_bounce == 0){
            return 3;
        }
        else if(num_abs == 1 and num_av_scat == 0 and num_other == 0 and num_psup_scat == 0 and num_near_reflect == 0 and num_inner_av_reflect == 0 and num_PMT_reflect == 1 and num_abs_other == 0 and num_scat_other == 0 and num_ropes == 0 and num_extwater_scat == 0 and num_av_pipes == 0 and num_acrylic_scatter == 0 and num_outer_av_reflect == 0 and num_acrylic_bounce == 0){
            return 4;
        }
        else if(num_abs == 0 and num_av_scat == 1 and num_other == 0 and num_psup_scat == 0 and num_near_reflect == 0 and num_inner_av_reflect == 0 and num_PMT_reflect == 1 and num_abs_other == 0 and num_scat_other == 0 and num_ropes == 0 and num_extwater_scat == 0 and num_av_pipes == 0 and num_acrylic_scatter == 0 and num_outer_av_reflect == 0 and num_acrylic_bounce == 0){
            return 5;
        }
        else if(num_abs == 0 and num_av_scat == 0 and num_other == 0 and num_psup_scat == 0 and num_near_reflect == 0 and num_inner_av_reflect == 0 and num_PMT_reflect == 2 and num_abs_other == 0 and num_scat_other == 0 and num_ropes == 0 and num_extwater_scat == 0 and num_av_pipes == 0 and num_acrylic_scatter == 0 and num_outer_av_reflect == 0 and num_acrylic_bounce == 0){
            return 6;
        }
        else if(num_abs == 0 and num_av_scat == 0 and num_other == 0 and num_psup_scat == 0 and num_near_reflect == 0 and num_inner_av_reflect == 2 and num_PMT_reflect == 0 and num_abs_other == 0 and num_scat_other == 0 and num_ropes == 0 and num_extwater_scat == 0 and num_av_pipes == 0 and num_acrylic_scatter == 0 and num_outer_av_reflect == 0 and num_acrylic_bounce == 0){
            return 7;
        }
        else if(num_abs == 1 and num_av_scat == 0 and num_other == 0 and num_psup_scat == 0 and num_near_reflect == 0 and num_inner_av_reflect == 1 and num_PMT_reflect == 0 and num_abs_other == 0 and num_scat_other == 0 and num_ropes == 0 and num_extwater_scat == 0 and num_av_pipes == 0 and num_acrylic_scatter == 0 and num_outer_av_reflect == 0 and num_acrylic_bounce == 0){
            return 8;
        }
        else if(num_abs == 0 and num_av_scat == 1 and num_other == 0 and num_psup_scat == 0 and num_near_reflect == 0 and num_inner_av_reflect == 1 and num_PMT_reflect == 0 and num_abs_other == 0 and num_scat_other == 0 and num_ropes == 0 and num_extwater_scat == 0 and num_av_pipes == 0 and num_acrylic_scatter == 0 and num_outer_av_reflect == 0 and num_acrylic_bounce == 0){
            return 9;
        }
        else if(num_abs == 0 and num_av_scat == 0 and num_other == 0 and num_psup_scat == 0 and num_near_reflect == 0 and num_inner_av_reflect == 0 and num_PMT_reflect == 0 and num_abs_other == 0 and num_scat_other == 0 and num_ropes == 0 and num_extwater_scat == 2 and num_av_pipes == 0 and num_acrylic_scatter == 0 and num_outer_av_reflect == 0 and num_acrylic_bounce == 0){
            return 10;
        }
        else if(num_abs == 1 and num_av_scat == 0 and num_other == 0 and num_psup_scat == 0 and num_near_reflect == 0 and num_inner_av_reflect == 0 and num_PMT_reflect == 0 and num_abs_other == 0 and num_scat_other == 0 and num_ropes == 0 and num_extwater_scat == 1 and num_av_pipes == 0 and num_acrylic_scatter == 0 and num_outer_av_reflect == 0 and num_acrylic_bounce == 0){
            return 11;
        }
        else if(num_abs == 0 and num_av_scat == 1 and num_other == 0 and num_psup_scat == 0 and num_near_reflect == 0 and num_inner_av_reflect == 0 and num_PMT_reflect == 0 and num_abs_other == 0 and num_scat_other == 0 and num_ropes == 0 and num_extwater_scat == 1 and num_av_pipes == 0 and num_acrylic_scatter == 0 and num_outer_av_reflect == 0 and num_acrylic_bounce == 0){
            return 12;
        }
        else if(num_abs == 0 and num_av_scat == 0 and num_other == 0 and num_psup_scat == 0 and num_near_reflect == 0 and num_inner_av_reflect == 1 and num_PMT_reflect == 0 and num_abs_other == 0 and num_scat_other == 0 and num_ropes == 0 and num_extwater_scat == 1 and num_av_pipes == 0 and num_acrylic_scatter == 0 and num_outer_av_reflect == 0 and num_acrylic_bounce == 0){
            return 13;
        }
        else if(num_abs == 0 and num_av_scat == 0 and num_other == 0 and num_psup_scat == 0 and num_near_reflect == 0 and num_inner_av_reflect == 0 and num_PMT_reflect == 1 and num_abs_other == 0 and num_scat_other == 0 and num_ropes == 0 and num_extwater_scat == 1 and num_av_pipes == 0 and num_acrylic_scatter == 0 and num_outer_av_reflect == 0 and num_acrylic_bounce == 0){
            return 14;
        }
        else if(num_abs == 0 and num_av_scat == 0 and num_other == 0 and num_psup_scat == 0 and num_near_reflect == 1 and num_inner_av_reflect == 0 and num_PMT_reflect == 0 and num_abs_other == 0 and num_scat_other == 0 and num_ropes == 0 and num_extwater_scat == 1 and num_av_pipes == 0 and num_acrylic_scatter == 0 and num_outer_av_reflect == 0 and num_acrylic_bounce == 0){
            return 15;
        }
        else if(num_abs == 0 and num_av_scat == 0 and num_other == 0 and num_psup_scat == 0 and num_near_reflect == 1 and num_inner_av_reflect == 0 and num_PMT_reflect == 1 and num_abs_other == 0 and num_scat_other == 0 and num_ropes == 0 and num_extwater_scat == 0 and num_av_pipes == 0 and num_acrylic_scatter == 0 and num_outer_av_reflect == 0 and num_acrylic_bounce == 0){
            return 16;
        }
        else if(num_abs == 0 and num_av_scat == 0 and num_other == 0 and num_psup_scat == 0 and num_near_reflect == 0 and num_inner_av_reflect == 0 and num_PMT_reflect == 1 and num_abs_other == 0 and num_scat_other == 0 and num_ropes == 1 and num_extwater_scat == 0 and num_av_pipes == 0 and num_acrylic_scatter == 0 and num_outer_av_reflect == 0 and num_acrylic_bounce == 0){
            return 17;
        }
        else if(num_abs == 0 and num_av_scat == 0 and num_other == 0 and num_psup_scat == 0 and num_near_reflect == 0 and num_inner_av_reflect == 0 and num_PMT_reflect == 0 and num_abs_other == 0 and num_scat_other == 0 and num_ropes == 1 and num_extwater_scat == 1 and num_av_pipes == 0 and num_acrylic_scatter == 0 and num_outer_av_reflect == 0 and num_acrylic_bounce == 0){
            return 18;
        }
        else if(num_abs == 0 and num_av_scat == 0 and num_other == 0 and num_psup_scat == 0 and num_near_reflect == 0 and num_inner_av_reflect == 0 and num_PMT_reflect == 0 and num_abs_other == 0 and num_scat_other == 0 and num_ropes == 2 and num_extwater_scat == 0 and num_av_pipes == 0 and num_acrylic_scatter == 0 and num_outer_av_reflect == 0 and num_acrylic_bounce == 0){
            return 19;
        }
        else if(num_abs == 0 and num_av_scat == 1 and num_other == 0 and num_psup_scat == 0 and num_near_reflect == 0 and num_inner_av_reflect == 0 and num_PMT_reflect == 0 and num_abs_other == 0 and num_scat_other == 0 and num_ropes == 0 and num_extwater_scat == 0 and num_av_pipes == 1 and num_acrylic_scatter == 0 and num_outer_av_reflect == 0 and num_acrylic_bounce == 0){
            return 20;
        }
        else if(num_abs == 0 and num_av_scat == 0 and num_other == 0 and num_psup_scat == 0 and num_near_reflect == 0 and num_inner_av_reflect == 0 and num_PMT_reflect == 0 and num_abs_other == 0 and num_scat_other == 0 and num_ropes == 0 and num_extwater_scat == 0 and num_av_pipes == 2 and num_acrylic_scatter == 0 and num_outer_av_reflect == 0 and num_acrylic_bounce == 0){
            return 21;
        }
        else if(num_abs == 0 and num_av_scat == 0 and num_other == 0 and num_psup_scat == 0 and num_near_reflect == 0 and num_inner_av_reflect == 1 and num_PMT_reflect == 0 and num_abs_other == 0 and num_scat_other == 0 and num_ropes == 0 and num_extwater_scat == 0 and num_av_pipes == 0 and num_acrylic_scatter == 1 and num_outer_av_reflect == 0 and num_acrylic_bounce == 0){
            return 22;
        }
        else if(num_abs == 0 and num_av_scat == 0 and num_other == 0 and num_psup_scat == 0 and num_near_reflect == 0 and num_inner_av_reflect == 1 and num_PMT_reflect == 1 and num_abs_other == 0 and num_scat_other == 0 and num_ropes == 0 and num_extwater_scat == 0 and num_av_pipes == 0 and num_acrylic_scatter == 0 and num_outer_av_reflect == 0 and num_acrylic_bounce == 0){
            return 23;
        }
        else if(num_abs == 0 and num_av_scat == 0 and num_other == 0 and num_psup_scat == 0 and num_near_reflect == 0 and num_inner_av_reflect == 0 and num_PMT_reflect == 0 and num_abs_other == 0 and num_scat_other == 0 and num_ropes == 0 and num_extwater_scat == 1 and num_av_pipes == 0 and num_acrylic_scatter == 0 and num_outer_av_reflect == 1 and num_acrylic_bounce == 0){
            return 24;
        }
        else if(num_abs == 0 and num_av_scat == 0 and num_other == 0 and num_psup_scat == 0 and num_near_reflect == 0 and num_inner_av_reflect == 0 and num_PMT_reflect == 1 and num_abs_other == 0 and num_scat_other == 0 and num_ropes == 0 and num_extwater_scat == 0 and num_av_pipes == 1 and num_acrylic_scatter == 0 and num_outer_av_reflect == 0 and num_acrylic_bounce == 0){
            return 25;
        }
        else if(num_abs == 0 and num_av_scat == 0 and num_other == 0 and num_psup_scat == 0 and num_near_reflect == 0 and num_inner_av_reflect == 0 and num_PMT_reflect == 0 and num_abs_other == 0 and num_scat_other == 0 and num_ropes == 1 and num_extwater_scat == 0 and num_av_pipes == 0 and num_acrylic_scatter == 1 and num_outer_av_reflect == 0 and num_acrylic_bounce == 0){
            return 26;
        }
        else if(num_abs == 0 and num_av_scat == 0 and num_other == 0 and num_psup_scat == 0 and num_near_reflect == 0 and num_inner_av_reflect == 0 and num_PMT_reflect == 1 and num_abs_other == 0 and num_scat_other == 0 and num_ropes == 0 and num_extwater_scat == 0 and num_av_pipes == 0 and num_acrylic_scatter == 1 and num_outer_av_reflect == 0 and num_acrylic_bounce == 0){
            return 27;
        }
        else if(num_abs == 0 and num_av_scat == 1 and num_other == 0 and num_psup_scat == 0 and num_near_reflect == 0 and num_inner_av_reflect == 0 and num_PMT_reflect == 0 and num_abs_other == 0 and num_scat_other == 0 and num_ropes == 0 and num_extwater_scat == 0 and num_av_pipes == 0 and num_acrylic_scatter == 1 and num_outer_av_reflect == 0 and num_acrylic_bounce == 0){
            return 28;
        }
        else if(num_abs == 0 and num_av_scat == 0 and num_other == 0 and num_psup_scat == 0 and num_near_reflect == 0 and num_inner_av_reflect == 0 and num_PMT_reflect == 0 and num_abs_other == 0 and num_scat_other == 1 and num_ropes == 1 and num_extwater_scat == 0 and num_av_pipes == 0 and num_acrylic_scatter == 0 and num_outer_av_reflect == 0 and num_acrylic_bounce == 0){
            return 29;
        }
        else if(num_abs == 0 and num_av_scat == 1 and num_other == 0 and num_psup_scat == 0 and num_near_reflect == 0 and num_inner_av_reflect == 0 and num_PMT_reflect == 0 and num_abs_other == 0 and num_scat_other == 0 and num_ropes == 0 and num_extwater_scat == 0 and num_av_pipes == 0 and num_acrylic_scatter == 0 and num_outer_av_reflect == 0 and num_acrylic_bounce == 1){
            return 30;
        }
        else if(num_abs == 0 and num_av_scat == 0 and num_other == 0 and num_psup_scat == 0 and num_near_reflect == 0 and num_inner_av_reflect == 0 and num_PMT_reflect == 0 and num_abs_other == 0 and num_scat_other == 0 and num_ropes == 0 and num_extwater_scat == 1 and num_av_pipes == 0 and num_acrylic_scatter == 0 and num_outer_av_reflect == 0 and num_acrylic_bounce == 1){
            return 31;
        }
        else if(num_abs == 0 and num_av_scat == 0 and num_other == 0 and num_psup_scat == 0 and num_near_reflect == 0 and num_inner_av_reflect == 0 and num_PMT_reflect == 1 and num_abs_other == 0 and num_scat_other == 0 and num_ropes == 0 and num_extwater_scat == 0 and num_av_pipes == 0 and num_acrylic_scatter == 0 and num_outer_av_reflect == 0 and num_acrylic_bounce == 1){
            return 32;
        }
        else{ //should not get here
            std::string output_string = "";
            if(num_abs != 0){
                output_string = output_string + "Abs: " + std::to_string(num_abs) + ", ";
            }
            if(num_av_scat != 0){
                output_string = output_string + "Scat: " + std::to_string(num_av_scat) + ", ";
            }
            if(num_other != 0){
                output_string = output_string + "Other: " + std::to_string(num_other) + ", ";
            }
            if(num_psup_scat != 0){
                output_string = output_string + "PSUP scat: " + std::to_string(num_psup_scat) + ", ";
            }
            if(num_near_reflect != 0){
                output_string = output_string + "Near reflect: " + std::to_string(num_near_reflect) + ", ";
            }
            if(num_inner_av_reflect != 0){
                output_string = output_string + "Inner av reflect: " + std::to_string(num_inner_av_reflect) + ", ";
            }
            if(num_PMT_reflect != 0){
                output_string = output_string + "PMT reflect: " + std::to_string(num_PMT_reflect) + ", ";
            }
            if(num_abs_other != 0){
                output_string = output_string + "Abs other: " + std::to_string(num_abs_other) + ", ";
            }
            if(num_scat_other != 0){
                output_string = output_string + "Other scat: " + std::to_string(num_scat_other) + ", ";
            }
            if(num_ropes != 0){
                output_string = output_string + "Ropes: " + std::to_string(num_ropes) + ", ";
            }
            if(num_extwater_scat != 0){
                output_string = output_string + "External scatter: " + std::to_string(num_extwater_scat) + ", ";
            }
            if(num_av_pipes != 0){
                output_string = output_string + "Pipes: " + std::to_string(num_av_pipes) + ", ";
            }
            if(num_acrylic_scatter != 0){
                output_string = output_string + "Acrylic scatter: " + std::to_string(num_acrylic_scatter) + ", ";
            }            
            if(num_outer_av_reflect != 0){
                output_string = output_string + "Outer AV reflect: " + std::to_string(num_outer_av_reflect) + ", ";
            }
            if(num_acrylic_bounce != 0){
                output_string = output_string + "Acrylic bounce: " + std::to_string(num_acrylic_bounce) + ", ";
            }

            std::cout << output_string << std::endl;

            return 33;
        }
    }
    else{
        return 0;
    }
}

int GetPMTID(std::string input){
    if(std::isdigit(input.substr(input.length()-4, 4)[0])) return std::stoi(input.substr(input.length()-4, 4));
    else if(std::isdigit(input.substr(input.length()-3, 3)[0])) return std::stoi(input.substr(input.length()-3, 3));
    else if(std::isdigit(input.substr(input.length()-2, 2)[0])) return std::stoi(input.substr(input.length()-2, 2));
    else if(std::isdigit(input.substr(input.length()-1, 1)[0])) return std::stoi(input.substr(input.length()-1, 1));
}
    
void WriteProcess(const RAT::DS::MCTrack &inputPhotonTrack){
    TVector3 prev_pos;
    for(int i_step = 0; i_step < inputPhotonTrack.GetMCTrackStepCount();i_step++){
        const RAT::DS::MCTrackStep &rTrackStep = inputPhotonTrack.GetMCTrackStep(i_step);
        std::cout << "Step: " << i_step << ", Start vol: " << rTrackStep.GetStartVolume() << ", End Vol: " << rTrackStep.GetEndVolume() << ", process: " << rTrackStep.GetProcess() << ", Start Pos: " << rTrackStep.GetPosition().X() << ", " << rTrackStep.GetPosition().Y() << ", " << rTrackStep.GetPosition().Z() << ", time: " << rTrackStep.GetGlobalTime() << std::endl;
        if(i_step != 0) std::cout << "Angle between this direction and last direction is " << rTrackStep.GetPosition().Angle(prev_pos) * (180/TMath::Pi()) << std::endl;
        prev_pos = rTrackStep.GetPosition();
    }
    std::cout << " " << std::endl;
}

int GetLightPaths(std::string file, std::string fibre, std::string data_type){

    bool verbose = true;
    bool debug = false;

    // Initialise variables and histograms
    int pmtcount = 0;
    double locality = 10.0; //lpc sensitivity
    double energy = RAT::util::WavelengthToEnergy(403E-6); //FIXME: could be input argument as easy to forget

    // get file name from path+filename string
    std::size_t botDirPos = file.find_last_of("/");
    std::string filename = file.substr(botDirPos+1, file.length());
    std::string saveroot = "Tracking_ResHitCosTheta_" + filename;
    TFile *rootfile = new TFile(saveroot.c_str(),"RECREATE");

    TH2F *h2DLPCTIR = new TH2F("h2DLPCTIR", "title", 1000, -1., 1., 20, 0, 2);
    TH2F *h2DLPCStraightLinePath = new TH2F("h2DLPCStraightLinePath", "title", 1000, -1., 1., 20, 0, 2);
    TH2F *h2DLPCTransitTime = new TH2F("h2DLPCTransitTime", "title", 1000, -1., 1., 900, 0, 90);

    TH1D *h1DWreflPaths = new TH1D("h1DWreflPaths", "AV reflect path", 1000, -1., 1.);
    TH1D *h1DWASAWPaths = new TH1D("h1DWASAWPaths", "Straight through path", 1000, -1., 1.);
    TH1D *h1DWAWPaths = new TH1D("h1DWAWPaths", "Water AV water path", 1000, -1., 1.);
    TH1D *h1DWaterPaths = new TH1D("h1DWaterPaths", "Water AV water path", 1000, -1., 1.);
    TH1D *h1DBelowMaxReflectionAngle = new TH1D("h1DBelowMaxReflectionAngle", "Below max angle", 1000, -1., 1.);

    if(verbose) std::cout << "Initialising RAT" << std::endl;

    // Initialise RAT
    RAT::DU::DSReader dsreader(file);
    RAT::DU::GroupVelocity groupVelocity = RAT::DU::Utility::Get()->GetGroupVelocity();
    RAT::DU::LightPathCalculator lightPath = RAT::DU::Utility::Get()->GetLightPathCalculator();
    lightPath.SetELLIEEvent(true); // event originates outside AV (see PR #2621)
    const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
    const int NPMTS = pmtinfo.GetCount();
    if(verbose) std::cout << NPMTS << " PMTs found" << std::endl;
    double transitTime[NPMTS];
    double bucketTime[NPMTS];
    double cosTheta[NPMTS];
    TVector2 pmtPosFlat[NPMTS];

    RAT::DB *db = RAT::DB::Get();
    RAT::DBLinkPtr entry = db->GetLink("FIBRE", fibre);
    TVector3 fibrePos(entry->GetD("x"), entry->GetD("y"), entry->GetD("z")); // position of fibre [mm]
    TVector3 fibreDir(entry->GetD("u"), entry->GetD("v"), entry->GetD("w")); // direction of fibre
    if(verbose) std::cout << "RATDB: fibre " << fibre << ", pos: (" << fibrePos.X() << "," << fibrePos.Y() << "," << fibrePos.Z() << "), dir: (" << fibreDir.X() << "," << fibreDir.Y() << "," << fibreDir.Z() << ")" << std::endl; 

    double maxReflectionAngle = ReflectionAngle(fibrePos);

    for (int it=0; it<NPMTS; it++){ //only need to calculate transit time for each PMT once, since fibre and pmt locations are fixed
        // Use normal or HQE PMTs only
        if (pmtinfo.GetType(it) != 1 && pmtinfo.GetType(it) != 7) continue;
        pmtcount++;

        // Get PMT information
        TVector3 pmtPos = pmtinfo.GetPosition(it);       // position [mm]
        TVector3 pmtDir = pmtinfo.GetDirection(it);      // direction
        TVector2 PMTPosFlat = IcosProject(pmtPos);

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
        pmtPosFlat[it] = PMTPosFlat;

        //fill histogram
        h2DLPCTIR->Fill(cosTheta_calc, lightPath.GetTIR()); //check total internal reflection flag
        h2DLPCStraightLinePath->Fill(cosTheta_calc, lightPath.GetStraightLine()); //check if straight line calculation was used
        h2DLPCTransitTime->Fill(cosTheta_calc, tof);
        if(lightPath.GetLightPathType() == "Water->AV->Scint/InnerAV->AV->Water"){
            h1DWASAWPaths->Fill(cosTheta_calc);
        }
        else if(lightPath.GetLightPathType() == "Water->Reflection->Water"){
            h1DWreflPaths->Fill(cosTheta_calc);
        }
        else if(lightPath.GetLightPathType() == "Water->AV->Water"){
            h1DWAWPaths->Fill(cosTheta_calc);
        }
        else if(lightPath.GetLightPathType() == "Water"){
            h1DWaterPaths->Fill(cosTheta_calc);
        }
        else{
            std::cout << "Not accounted for path: " << lightPath.GetLightPathType() << std::endl;
        }
        if(fibrePos.Angle(pmtPos) < maxReflectionAngle){
            h1DBelowMaxReflectionAngle->Fill(cosTheta_calc);
        }
    }

    if (data_type == "MC") {
        // Create histograms
        TH1D* hNhits = new TH1D("hNhits", "nhits", 101, 0, 100);

        TH2F *hPMTResTimeCosTheta = new TH2F("hPmtResTimeVsCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hNoiseResTimeCosTheta = new TH2F("hNoiseResTimeVsCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hSingleScatterResTimeCosTheta = new TH2F("hSingleScatterResTimeVsCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hReemissionResTimeCosTheta = new TH2F("hReemissionResTimeVsCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hOtherEffectResTimeCosTheta = new TH2F("hOtherEffectResTimeVsCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hNoEffectResTimeCosTheta = new TH2F("hNoEffectResTimeVsCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hNearReflectResTimeCosTheta = new TH2F("hNearReflectResTimeVsCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hRopesResTimeCosTheta = new TH2F("hRopesResTimeVsCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hPMTReflectionResTimeCosTheta = new TH2F("hPMTReflectionResTimeVsCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hPMTReflectionAVReflectionResTimeCosTheta = new TH2F("hPMTReflectionAVReflectionResTimeCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hExtWaterScatterResTimeCosTheta = new TH2F("hExtWaterScatterResTimeVsCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hInnerAvReflectionResTimeCosTheta = new TH2F("hInnerAvReflectionResTimeVsCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hMultipleEffectResTimeCosTheta = new TH2F("hMultipleEffectResTimeVsCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hAVPipesResTimeCosTheta = new TH2F("hAVPipesResTimeVsCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hAcrylicScatterResTimeCosTheta = new TH2F("hAcrylicScatterResTimeVsCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hOtherScatterResTimeCosTheta = new TH2F("OtherScatterResTimeVsCosTheta", "title",1000, -1., 1., 1000, -50., 250.);

        TH1D *h1DResTimeAll = new TH1D("h1DResTimeAll", "Residual Hit Time", 1000, -50., 250.);
        TH1D *h1DNoEffectResTime = new TH1D("h1DNoEffectResTime", "MC Residual Front End Time", 1000, -50., 250.);
        TH1D *h1DNoEffectTime = new TH1D("h1DNoEffectTime", "MC Front End Time", 1000, -50., 250.);
        TH1D *h1DNoEffectMCTransitTime = new TH1D("h1DNoEffectMCTransitTime", "MC Transit Time", 1000, -50., 250.);
        TH1D *h1DNoEffectElectronicsTime = new TH1D("h1DNoEffectElectronicsTime", "Front End Time - Transit Time", 1000, -25., 200.);

        TH2D *h1DDiscontinuityPMTs = new TH2D("h1DDiscontinuityPMTs", "title",9729, 0, 9728, 1000, -50., 250.);

        TH2F *hMultipleDoubleAbsResTimeCosTheta = new TH2F("hMultipleDoubleAbsResTimeCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hMultipleDoubleScatResTimeCosTheta = new TH2F("hMultipleDoubleScatResTimeCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hMultipleAbsScatScatterResTimeCosTheta = new TH2F("hMultipleAbsScatScatterResTimeCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hMultiplePMTReflectAbsAbsResTimeCosTheta = new TH2F("hMultiplePMTReflectAbsAbsResTimeCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hMultiplePMTReflectScatResTimeCosTheta = new TH2F("hMultiplePMTReflectScatResTimeCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hMultipleDoublePMTReflectResTimeCosTheta = new TH2F("hMultipleDoublePMTReflectResTimeCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hMultipleDoubleInnerAVReflectResTimeCosTheta = new TH2F("hMultipleDoubleInnerAVReflectResTimeCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hMultipleInnerAVReflectAbsResTimeCosTheta = new TH2F("hMultipleInnerAVReflectAbsResTimeCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hMultipleInnerAVReflectScatResTimeCosTheta = new TH2F("hMultipleInnerAVReflectScatResTimeCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hMultipleDoubleExternalScatResTimeCosTheta = new TH2F("hMultipleDoubleExternalScatResTimeCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hMultipleExtenalScatAbsResTimeCosTheta = new TH2F("hMultipleExtenalScatAbsResTimeCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hMultipleExternalScatInternalScatResTimeCosTheta = new TH2F("hMultipleExternalScatInternalScatResTimeCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hMultipleExternalScatInnerAvReflectResTimeCosTheta = new TH2F("hMultipleExternalScatInnerAvReflectResTimeCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hMultipleExternalScatPMTReflectResTimeCosTheta = new TH2F("hMultipleExternalScatPMTReflectResTimeCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hMultipleExternalScatNearReflectResTimeCosTheta = new TH2F("hMultipleExternalScatNearReflectResTimeCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hMultipleNearReflectPMTReflectResTimeCosTheta = new TH2F("hMultipleNearReflectPMTReflectResTimeCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hMultiplePMTReflectRopesResTimeCosTheta = new TH2F("hMultiplePMTReflectRopesResTimeCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hMultipleExternalScatRopesResTimeCosTheta = new TH2F("hMultipleExternalScatRopesResTimeCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hMultipleDoubleRopesResTimeCosTheta = new TH2F("hMultipleDoubleRopesResTimeCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hMultipleScatterPipesResTimeCosTheta = new TH2F("hMultipleScatterPipesResTimeCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hMultipleDoublePipesResTimeCosTheta = new TH2F("hMultipleDoublePipesResTimeCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hMultipleInnerAVReflectAcrylicScatResTimeCosTheta = new TH2F("hMultipleInnerAVReflectAcrylicScatResTimeCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hMultipleInnerAVReflectPMTReflectResTimeCosTheta = new TH2F("hMultipleInnerAVReflectPMTReflectResTimeCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hMultipleOuterAVReflectExternalScatResTimeCosTheta = new TH2F("hMultipleOuterAVReflectExternalScatResTimeCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hMultiplePMTReflectPipesResTimeCosTheta = new TH2F("hMultiplePMTReflectPipesResTimeCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hMultipleRopesAcrylicScatResTimeCosTheta = new TH2F("hMultipleRopesAcrylicScatResTimeCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hMultiplePMTReflectAcrylicScatResTimeCosTheta = new TH2F("hMultiplePMTReflectAcrylicScatResTimeCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hMultipleScatterAcrylicScatResTimeCosTheta = new TH2F("hMultipleScatterAcrylicScatResTimeCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hMultipleOtherScatRopesResTimeCosTheta = new TH2F("hMultipleOtherScatRopesResTimeCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hMultipleScatterAcrylicBounceResTimeCosTheta = new TH2F("hMultipleScatterAcrylicBounceResTimeCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hMultipleExternalScatterAcrylicBounceResTimeCosTheta = new TH2F("hMultipleExternalScatterAcrylicBounceResTimeCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hMultiplePMTReflectAcrylicBounceResTimeCosTheta = new TH2F("hMultiplePMTReflectAcrylicBounceResTimeCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hMultipleOtherTimeCosTheta = new TH2F("hMultipleOtherTimeCosTheta", "title",1000, -1., 1., 1000, -50., 250.);
        TH2F *hMultipleMoreThan2EffectResTimeCosTheta = new TH2F("hMultipleMoreThan2EffectResTimeCosTheta", "title",1000, -1., 1., 1000, -50., 250.);


        std::vector<Double_t> mctimes;
        std::vector<Double_t> evtimes;
        std::vector<int> mcids;
        std::vector<int> evids;
        Double_t gttime;
        bool print_status = false;

        size_t nEV = 0;
        size_t entryCount = dsreader.GetEntryCount(); //number of entries, want to loop over each one
        if(verbose) std::cout << "No of entries in run: " << entryCount << " events" << std::endl;
        for (size_t iEntry = 0; iEntry < entryCount; ++iEntry){
            if (iEntry %100 == 0 and verbose) std::cout << "Entry no " << iEntry << std::endl;
            const RAT::DS::Entry &rDS = dsreader.GetEntry(iEntry);
            const RAT::DS::MC &rMC = rDS.GetMC();
            for(size_t i_ev = 0; i_ev< rDS.GetMCEVCount(); ++i_ev){
                const RAT::DS::MCEV &rMCEV = rDS.GetMCEV(i_ev);
                const RAT::DS::MCHits &rMCHits = rMCEV.GetMCHits();
                const RAT::DS::EV &rEV = rDS.GetEV(i_ev);
                const RAT::DS::CalPMTs &calPMTs = rEV.GetCalPMTs();
                size_t n_hits_MCHits = rMCHits.GetCount();
                size_t calPMT_count = calPMTs.GetCount();
                hNhits->Fill(n_hits_MCHits);

                //get ev PMT ids

                std::vector<int> evPMTIDs;
                std::vector<double> evPMTTimes;

                for(size_t i_evpmt = 0; i_evpmt < calPMT_count; ++i_evpmt){
                    evPMTIDs.push_back(calPMTs.GetPMT(i_evpmt).GetID());
                    evPMTTimes.push_back(calPMTs.GetPMT(i_evpmt).GetTime() - transitTime[i_evpmt] - 390 + rMCEV.GetGTTime());
                    h1DResTimeAll->Fill(calPMTs.GetPMT(i_evpmt).GetTime() - transitTime[i_evpmt] - bucketTime[i_evpmt] - 390 + rMCEV.GetGTTime());
                }

                //now do main tracking
                std::vector<int> MCPMTIDs; 

                for(size_t i_mcpmt = 0; i_mcpmt < rMC.GetMCPMTCount(); ++i_mcpmt){
                    const RAT::DS::MCPMT &rMCPMT = rMC.GetMCPMT(i_mcpmt);
                    if(rMCPMT.GetMCPECount() == 0) continue; //should be a photo electron produced at pmt

                    UInt_t pmtID = rMCPMT.GetID();

                    if(std::find(evPMTIDs.begin(), evPMTIDs.end(), pmtID) == evPMTIDs.end()) continue; //PMT not in EV branch

                    MCPMTIDs.push_back(pmtID);

                    for (int z=0;z<rMCPMT.GetMCPECount();z++) { //FIXME: multiple PEs?
                        Double_t t = rMCPMT.GetMCPE(z).GetFrontEndTime();
                        Double_t t_res = t - transitTime[pmtID] - bucketTime[pmtID];

                        hPMTResTimeCosTheta->Fill(cosTheta[pmtID], t_res);

                        //now, find photon process and fill relevant histogram
                        //first check is for pmts with no photons, these are noise
                        if (rMCPMT.GetMCPhotonCount() == 0) {
                            hNoiseResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                        } else {
                            //get the photon track associated with the PE
                            const RAT::DS::MCPE &rPE = rMCPMT.GetMCPE(z);
                            if(rPE.GetNoise() != 1 and rPE.GetAfterPulse() != 1){ //not noise or afterpulse
                                UInt_t trackID = rPE.GetPhotonTrackID();
                                const RAT::DS::MCTrack &rPhotonTrack = rMC.GetMCTrack(trackID);
                                int photonProcess = GetPhotonProcess(rPhotonTrack, rMC, iEntry, debug);

                                h1DNoEffectElectronicsTime->Fill(t - (rPhotonTrack.GetMCTrackStep(rPhotonTrack.GetMCTrackStepCount()-1).GetGlobalTime() - rPhotonTrack.GetMCTrackStep(0).GetGlobalTime())); //time from PE to front end readout

                                if(photonProcess == 0){
                                    hSingleScatterResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                }
                                else if(photonProcess == 1){
                                    hReemissionResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                }
                                else if(photonProcess == 2){
                                    hNoEffectResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                    h1DNoEffectResTime->Fill(t_res);
                                    h1DNoEffectTime->Fill(t);
                                    h1DNoEffectMCTransitTime->Fill(rPhotonTrack.GetMCTrackStep(rPhotonTrack.GetMCTrackStepCount()-1).GetGlobalTime() - rPhotonTrack.GetMCTrackStep(0).GetGlobalTime());
                                }
                                else if(photonProcess == 3){
                                    hMultipleEffectResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                    int multipleStatus = GetMultipleEffect(rPhotonTrack, rMC);
                                    if(multipleStatus == 0){
                                        hMultipleMoreThan2EffectResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                    }
                                    if(multipleStatus == 1){
                                        hMultipleDoubleAbsResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                    }
                                    if(multipleStatus == 2){
                                        hMultipleDoubleScatResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                    }
                                    if(multipleStatus == 3){
                                        hMultipleAbsScatScatterResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                    }
                                    if(multipleStatus == 4){
                                        hMultiplePMTReflectAbsAbsResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                    }
                                    if(multipleStatus == 5){
                                        hMultiplePMTReflectScatResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                    }
                                    if(multipleStatus == 6){
                                        hMultipleDoublePMTReflectResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                    }
                                    if(multipleStatus == 7){
                                        hMultipleDoubleInnerAVReflectResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                    }
                                    if(multipleStatus == 8){
                                        hMultipleInnerAVReflectAbsResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                    }
                                    if(multipleStatus == 9){
                                        hMultipleInnerAVReflectScatResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                    }
                                    if(multipleStatus == 10){
                                        hMultipleDoubleExternalScatResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                    }
                                    if(multipleStatus == 11){
                                        hMultipleExtenalScatAbsResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                    }
                                    if(multipleStatus == 12){
                                        hMultipleExternalScatInternalScatResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                    }
                                    if(multipleStatus == 13){
                                        hMultipleExternalScatInnerAvReflectResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                    }
                                    if(multipleStatus == 14){
                                        hMultipleExternalScatPMTReflectResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                    }
                                    if(multipleStatus == 15){
                                        hMultipleExternalScatNearReflectResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                    }
                                    if(multipleStatus == 16){
                                        hMultipleNearReflectPMTReflectResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                    }
                                    if(multipleStatus == 17){
                                        hMultiplePMTReflectRopesResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                    }
                                    if(multipleStatus == 18){
                                        hMultipleExternalScatRopesResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                    }
                                    if(multipleStatus == 19){
                                        hMultipleDoubleRopesResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                    }
                                    if(multipleStatus == 20){
                                        hMultipleScatterPipesResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                    }
                                    if(multipleStatus == 21){
                                        hMultipleDoublePipesResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                    }
                                    if(multipleStatus == 22){
                                        hMultipleInnerAVReflectAcrylicScatResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                    }
                                    if(multipleStatus == 23){
                                        hMultipleInnerAVReflectPMTReflectResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                    }
                                    if(multipleStatus == 24){
                                        hMultipleOuterAVReflectExternalScatResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                    }
                                    if(multipleStatus == 25){
                                        hMultiplePMTReflectPipesResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                    }
                                    if(multipleStatus == 26){
                                        hMultipleRopesAcrylicScatResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                    }
                                    if(multipleStatus == 27){
                                        hMultiplePMTReflectAcrylicScatResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                    }
                                    if(multipleStatus == 28){
                                        hMultipleScatterAcrylicScatResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                    }
                                    if(multipleStatus == 29){
                                        hMultipleOtherScatRopesResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                    }
                                    if(multipleStatus == 30){
                                        hMultipleScatterAcrylicBounceResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                    }
                                    if(multipleStatus == 31){
                                        hMultipleExternalScatterAcrylicBounceResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                    }
                                    if(multipleStatus == 32){
                                        hMultiplePMTReflectAcrylicBounceResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                    }
                                    if(multipleStatus == 33){
                                        hMultipleOtherTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                    }
                                }
                                else if(photonProcess == 4){
                                    hOtherEffectResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                }
                                else if(photonProcess == 5){
                                    hNearReflectResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                }
                                else if (photonProcess == 6){
                                    hRopesResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                }
                                else if(photonProcess == 7){
                                    hPMTReflectionResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                }
                                else if(photonProcess == 8){
                                    hExtWaterScatterResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                }
                                else if(photonProcess == 9){
                                    hInnerAvReflectionResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                }
                                else if(photonProcess == 10){
                                    hAVPipesResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                }
                                else if(photonProcess == 11){
                                    hAcrylicScatterResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                }
                                else if(photonProcess == 12){
                                    hOtherScatterResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                }
                                else if(photonProcess == 13){
                                    hPMTReflectionAVReflectionResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                                }
                            }
                            else{
                                hNoiseResTimeCosTheta->Fill(cosTheta[pmtID], t_res);
                            }
                        }
                    }
                }

                //compare and add others to noise
                for(int a=0;a<evPMTIDs.size();a++){
                    if(std::find(MCPMTIDs.begin(), MCPMTIDs.end(), evPMTIDs.at(a)) == MCPMTIDs.end()){
                        hNoiseResTimeCosTheta->Fill(cosTheta[evPMTIDs.at(a)], evPMTTimes.at(a));
                    }
                }
            }
        }

        //now write everything
        rootfile->cd();
        hNhits->Write();
        hPMTResTimeCosTheta->Write();
        hNoiseResTimeCosTheta->Write();
        hSingleScatterResTimeCosTheta->Write();
        hReemissionResTimeCosTheta->Write();
        hOtherEffectResTimeCosTheta->Write();
        hNoEffectResTimeCosTheta->Write();
        hMultipleEffectResTimeCosTheta->Write();
        hNearReflectResTimeCosTheta->Write();
        hRopesResTimeCosTheta->Write();
        hPMTReflectionResTimeCosTheta->Write();
        hPMTReflectionAVReflectionResTimeCosTheta->Write();
        hExtWaterScatterResTimeCosTheta->Write();
        hInnerAvReflectionResTimeCosTheta->Write();
        hAVPipesResTimeCosTheta->Write();
        hAcrylicScatterResTimeCosTheta->Write();
        hOtherScatterResTimeCosTheta->Write();
        h1DNoEffectResTime->Write();
        h1DNoEffectTime->Write();
        h1DNoEffectMCTransitTime->Write();
        h1DNoEffectElectronicsTime->Write();
        h1DResTimeAll->Write();
        h2DLPCTIR->Write();
        h2DLPCStraightLinePath->Write();
        h2DLPCTransitTime->Write();
        h1DWASAWPaths->Write();
        h1DWreflPaths->Write();
        h1DWAWPaths->Write();
        h1DWaterPaths->Write();
        h1DBelowMaxReflectionAngle->Write();
        h1DDiscontinuityPMTs->Write();
        hMultipleMoreThan2EffectResTimeCosTheta->Write();
        hMultipleDoubleAbsResTimeCosTheta->Write();
        hMultipleDoubleScatResTimeCosTheta->Write();
        hMultipleAbsScatScatterResTimeCosTheta->Write();
        hMultiplePMTReflectAbsAbsResTimeCosTheta->Write();
        hMultiplePMTReflectScatResTimeCosTheta->Write();
        hMultipleDoublePMTReflectResTimeCosTheta->Write();
        hMultipleDoubleInnerAVReflectResTimeCosTheta->Write();
        hMultipleInnerAVReflectAbsResTimeCosTheta->Write();
        hMultipleInnerAVReflectScatResTimeCosTheta->Write();
        hMultipleDoubleExternalScatResTimeCosTheta->Write();  
        hMultipleExtenalScatAbsResTimeCosTheta->Write();  
        hMultipleExternalScatInternalScatResTimeCosTheta->Write();
        hMultipleExternalScatInnerAvReflectResTimeCosTheta->Write();  
        hMultipleExternalScatPMTReflectResTimeCosTheta->Write();
        hMultipleExternalScatNearReflectResTimeCosTheta->Write();
        hMultipleNearReflectPMTReflectResTimeCosTheta->Write();
        hMultiplePMTReflectRopesResTimeCosTheta->Write();
        hMultipleExternalScatRopesResTimeCosTheta->Write();
        hMultipleDoubleRopesResTimeCosTheta->Write();
        hMultipleScatterPipesResTimeCosTheta->Write();
        hMultipleDoublePipesResTimeCosTheta->Write();
        hMultipleInnerAVReflectAcrylicScatResTimeCosTheta->Write();
        hMultipleInnerAVReflectPMTReflectResTimeCosTheta->Write();
        hMultipleOuterAVReflectExternalScatResTimeCosTheta->Write();
        hMultiplePMTReflectPipesResTimeCosTheta->Write();
        hMultipleRopesAcrylicScatResTimeCosTheta->Write();
        hMultiplePMTReflectAcrylicScatResTimeCosTheta->Write();
        hMultipleScatterAcrylicScatResTimeCosTheta->Write();
        hMultipleOtherScatRopesResTimeCosTheta->Write();
        hMultipleScatterAcrylicBounceResTimeCosTheta->Write();
        hMultipleExternalScatterAcrylicBounceResTimeCosTheta->Write();
        hMultiplePMTReflectAcrylicBounceResTimeCosTheta->Write();
        hMultipleOtherTimeCosTheta->Write();
        rootfile->Write();

    } else if (data_type == "raw") {
        // Create histograms
        TH1D* hNhits = new TH1D("hNhits", "nhits", 101, 0, 100);

        TH1D *h1DResTimeAll_raw = new TH1D("g", "Residual Hit Time, not reajusted", 1000, -50., 500.);
        TH1D *h1DResTimeAll = new TH1D("h1DResTimeAll", "Residual Hit Time", 1000, -50., 250.);
        TH2F *hPMTResTimeCosTheta = new TH2F("hPmtResTimeVsCosTheta", "title",1000, -1., 1., 1000, -50., 250.);

        size_t entryCount = dsreader.GetEntryCount(); //number of entries, want to loop over each one
        std::vector<Double_t> evPMTTimes;
        std::vector<UInt_t> pmtID;
        int index = 0;
        if(verbose) std::cout << "No of entries in run: " << entryCount << " events" << std::endl;
        for (size_t iEntry = 0; iEntry < entryCount; ++iEntry){
            if (iEntry %100 == 0 and verbose) std::cout << "Entry no " << iEntry << std::endl;
            const RAT::DS::Entry &rDS = dsreader.GetEntry(iEntry);
            for(size_t i_ev = 0; i_ev< rDS.GetEVCount(); ++i_ev){
                if (i_ev %100 == 0 and verbose) std::cout << "Event no " << i_ev << std::endl;
                const RAT::DS::EV &rEV = rDS.GetEV(i_ev);
                Int_t triggerWord = rEV.GetTrigType();
                if(!(triggerWord & (1 << 15))) continue; // EXTA cut
                std::cout << "passed" << std::endl;

                const RAT::DS::CalPMTs &calPMTs = rEV.GetCalPMTs(); 
                size_t calPMT_count = calPMTs.GetCount();
                hNhits->Fill(calPMT_count);
                
                // calculate time residuals (not ajusted for peak hit time yet), to fit gaussian
                for (size_t i_evpmt = 0; i_evpmt < calPMT_count; ++i_evpmt) {
                    const RAT::DS::PMTCal& pmtCal = calPMTs.GetPMT(i_evpmt);
                    pmtID.push_back(pmtCal.GetID());
                    evPMTTimes.push_back(pmtCal.GetTime() - transitTime[pmtID[index]] - bucketTime[pmtID[index]]);
                    h1DResTimeAll_raw->Fill(evPMTTimes[index]);
                    ++index;
                }

                
            }
        }

        // // Fit hist to gaussian
        // h1DResTimeAll_raw->Fit("gaus");
        // TF1 *fit = h1DResTimeAll_raw->GetFunction("gaus");
        // //Double_t chi2 = fit->GetChisquare();

        // // Value of the first parameter:
        // Double_t peak_time = fit->GetParameter(1);
        // //Double_t e1 = fit->GetParError(0);

        // Find bin with highest count -> peak
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
        hNhits->Write();
        hPMTResTimeCosTheta->Write();
        h1DResTimeAll->Write();
        h1DResTimeAll_raw->Write();
        rootfile->Write();
        rootfile->Close();

    } else {
        std::cout << "Wrong data type. Should be MC or raw" << std::endl;
        throw;
    }

    // Close everything
    rootfile->Close();
}