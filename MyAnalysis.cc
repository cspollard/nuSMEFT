// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PromptFinalState.hh" 
#include "Rivet/Projections/InvisibleFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/NeutralFinalState.hh"
#include <iostream>
#include <random>
#include <cmath>
#include <string>
#include <fstream>

// //for CDF do
// const bool CDF = true;
// const bool lhcrun3 = false;

// //for Atlas do:
// const bool CDF = false;
// const bool lhcrun3 = false;

//for Atlas at 13.6TeV do:
const bool CDF = false;
const bool lhcrun3 = true;

namespace Rivet {

 
  class MyAnalysis : public Analysis {
  public:
 
    /// Constructor
    MyAnalysis()
      : Analysis("MyAnalysis")
    { }
 

    /// @name Analysis methods
    ///@{
 
    /// Book histograms and initialize projections before the run
    void init() {
 
      // projection to find the electrons
      IdentifiedFinalState elecs(Cuts::abseta < 2.47 && Cuts::pT > 10*GeV);
      elecs.acceptIdPair(PID::ELECTRON);
      declare(elecs, "elecs");

      // projection to find the muons
      IdentifiedFinalState muons(Cuts::abseta < 2.4 && Cuts::pT > 10*GeV);
      muons.acceptIdPair(PID::MUON);
      declare(muons, "muons");
      
      // projection to find the electron neutrino
      IdentifiedFinalState neut(Cuts::pT > 0*GeV);
      neut.acceptIdPair(PID::NU_E);
      declare(neut, "neut");

      // projection to find the sterile neutrino
      IdentifiedFinalState sneut(Cuts::pT > 0*GeV);
      sneut.acceptIdPair(9900012);
      declare(sneut, "sneut");

      // // Jet Finder
      // const FinalState fs;
      // FastJets jetfs(fs, JetAlg::ANTIKT, 1, JetMuons::NONE, JetInvisibles::NONE);
      // declare(jetfs, "jets");

      // all tracks (to do deltaR with leptons)
      declare(ChargedFinalState(Cuts::abseta < 3 && Cuts::pT > 0.5*GeV), "cfs");

      // for pTmiss
      declare(VisibleFinalState(Cuts::abseta < 4.9 && Cuts::abspid!=9900012 && Cuts::abspid!=12), "vfs");
      declare(MissingMomentum(), "MissingMomentum");
      declare(NeutralFinalState(Cuts::abseta < 4.9), "nfs");

      double ptbinlow = 30.0;
      double ptbinhigh = 55.0;
      int ptbincount = 50;

      double mtbinlow = 60.0;
      double mtbinhigh = 100.0;
      int mtbincount = 80;

      //different binning for CDF vs atlas
      if (CDF == false) {
        ptbinlow = 30.;
        ptbinhigh = 60.;
        ptbincount = 60;

        mtbinlow = 60.;
        mtbinhigh = 120.;
        mtbincount = 120;
      }

      // // Each histogram is commented to describe what it does.
      // // Note that where I say 'no kinematic cuts' I mean none further than the generator level cuts and the global 10GeV pt
      // // restriction on electrons in this Rivet analysis, as well as the requirement each event have one electron and
      // // one neutrino.

      // reconstructed quantities, after kinematic cuts corresponding to combination paper (arXiv:2308.09417)
      book(_h["e pt"], "hist_e_pT_reconstructed", ptbincount, ptbinlow, ptbinhigh);      //reconstructed pT of electron
      book(_h["mt"], "hist_mT_reconstructed",mtbincount, mtbinlow, mtbinhigh);           //mT of electron and N (or n) using reconstructed missing energy and reconstructed electron data
      book(_h["pt miss"], "hist_pT_miss_reconstructed", ptbincount, ptbinlow, ptbinhigh);//reconstructed missing energy

      // reconstructed quantities, weighting to CDF PDF reweighting, after kinematic cuts corresponding to combination paper (arXiv:2308.09417)
      book(_hw2["e pt"], "hist_e_pT_reconstructed_weight2", ptbincount, ptbinlow, ptbinhigh);      //reconstructed pT of electron
      book(_hw2["mt"], "hist_mT_reconstructed_weight2",mtbincount, mtbinlow, mtbinhigh);           //mT of electron and N (or n) using reconstructed missing energy and reconstructed electron data
      book(_hw2["pt miss"], "hist_pT_miss_reconstructed_weight2", ptbincount, ptbinlow, ptbinhigh);//reconstructed missing energy
      

      // truth level histograms, after kinematic cuts corresponding to combination paper (arXiv:2308.09417)
      // the two mt histograms are either calculated with truth neutrino pt or truth missing pt
      book(_h["e pt truth"], "hist_e_pT_true", 100, 0.0, 100.0);                //true pT of electron
      book(_h["n pt truth"], "hist_N_pT_true", 100, 0.0, 100.0);                //true pT of N (or n)
      book(_h["pt miss truth"], "hist_pT_miss_true", 100, 0.0, 100.0);          //true missing energy
      book(_h["mt e pt miss truth"], "hist_mT_epT_miss_true", 100, 0.0, 100.0);  //mT of electron and N (or n) using true missing energy and true electron data
      book(_h["mt truth"], "hist_mT_eN_true", 100, 0.0, 100.0);              //mT of electron and N (or n) using only true data
      book(_h["m inv"], "hist_m_inv_true", 100, 0.0, 100.0);              //invariant mass using all only true data

      // reconstructed quantities, no cuts besides requiring one electron and one neutrino
      book(_h["e pt before"], "hist_e_pT_before_reconstructed", 100, 0.0, 100.0);          //reconstructed pT of electron
      book(_h["mt before"], "hist_mT_before_reconstructed",100, 0.0, 100.0);            //mT of electron and N (or n) using reconstructed missing energy and reconstructed electron data
      book(_h["pt miss before"], "hist_pT_miss_before_reconstructed", 100, 0.0, 100.0); //reconstructed missing energy
      book(_h["e abseta before"], "hist_e_abseta_before", 50, 0, 2.5);
      book(_h["ETU"], "hist_E_T_U", 100, 0.0, 100.0);

      double mw_count = 150;
      double mw_low = 79.750;
      double mw_high = 79.900;

      // 2d histos around mw
      // x variable is given in name, y variable is Breit-Wigner reweighting 
      // after kinematic cuts corresponding to combination paper (arXiv:2308.09417)
      book(_h2["e pt"], "hist_e_pT_2d", mw_count, mw_low, mw_high, ptbincount, ptbinlow, ptbinhigh);
      book(_h2["mt"], "hist_mT_2d", mw_count, mw_low, mw_high, mtbincount, mtbinlow, mtbinhigh);
      book(_h2["pt miss"], "hist_pT_miss_2d", mw_count, mw_low, mw_high, ptbincount, ptbinlow, ptbinhigh);
      // book(_h2w1["e pt pdf"], "hist_e_pT_2d_weighted", 150, 79.75, 79.90, ptbincount, ptbinlow, ptbinhigh);
      // book(_h2w1["mt pdf"], "hist_mT_2d_weighted", 150, 79.75, 79.90, mtbincount, mtbinlow, mtbinhigh);
      // book(_h2w1["pt miss pdf"], "hist_pT_miss_2d_weighted", 150, 79.75, 79.90, ptbincount, ptbinlow, ptbinhigh);
      book(_h2w2["e pt pdf 2"], "hist_e_pT_2d_weight2",mw_count, mw_low, mw_high, ptbincount, ptbinlow, ptbinhigh);
      book(_h2w2["mt pdf 2"], "hist_mT_2d_weight2", mw_count, mw_low, mw_high, mtbincount, mtbinlow, mtbinhigh);
      book(_h2w2["pt miss pdf 2"], "hist_pT_miss_2d_weight2", mw_count, mw_low, mw_high, ptbincount, ptbinlow, ptbinhigh);
      // book(_h2w3["e pt pdf 3"], "hist_e_pT_2d_weight3", 150, 79.75, 79.90, ptbincount, ptbinlow, ptbinhigh);
      // book(_h2w3["mt pdf 3"], "hist_mT_2d_weight3", 150, 79.75, 79.90, mtbincount, mtbinlow, mtbinhigh);
      // book(_h2w3["pt miss pdf 3"], "hist_pT_miss_2d_weight3", 150, 79.75, 79.90, ptbincount, ptbinlow, ptbinhigh);

      // mass fit in four lepton rapidity bins for Atlas
      // 2d histograms of ept or mt
      // after kinematic cuts corresponding to combination paper (arXiv:2308.09417)
      book(_h2["e pt 0-0.8"], "e_pT_2d_0_08", mw_count, mw_low, mw_high, ptbincount, ptbinlow, ptbinhigh);
      book(_h2["e pt 0.8-1.4"], "e_pT_2d_08_14", mw_count, mw_low, mw_high, ptbincount, ptbinlow, ptbinhigh);
      book(_h2["e pt 1.4-2.0"], "e_pT_2d_14_20", mw_count, mw_low, mw_high, ptbincount, ptbinlow, ptbinhigh);
      book(_h2["e pt 2.0-2.4"], "e_pT_2d_20_24", mw_count, mw_low, mw_high, ptbincount, ptbinlow, ptbinhigh);
      book(_h2["mt 0-0.8"], "hist_mT_2d_0_08", mw_count, mw_low, mw_high, mtbincount, mtbinlow, mtbinhigh);
      book(_h2["mt 0.8-1.4"], "hist_mT_2d_08_14", mw_count, mw_low, mw_high, mtbincount, mtbinlow, mtbinhigh);
      book(_h2["mt 1.4-2.0"], "hist_mT_2d_14_20", mw_count, mw_low, mw_high, mtbincount, mtbinlow, mtbinhigh);
      book(_h2["mt 2.0-2.4"], "hist_mT_2d_20_24", mw_count, mw_low, mw_high, mtbincount, mtbinlow, mtbinhigh);

      //quark histos
      // u quark pz, no kinematic cuts
      book(_h["upzbefore"], "abs_up_pz_before", 100, 0.0, 100.0); 
      // d quark pz, no kinematic cuts
      book(_h["dpzbefore"], "abs_down_pz_before", 100, 0.0, 100.0);
      // log_10 (u quark pz/d quark pz), no kinematic cuts
      book(_h["udpzbefore"], "abs_u_pz_over_abs_d_pz_before", 100,-2.5,2.5);            //log base 10
      // log_10 (u quark pz/d quark pz), no kinematic cuts, PDF reweighted for CDF
      book(_hw2["udpzbefore2"], "abs_u_pz_over_abs_d_pz_weight2_before", 100,-2.5,2.5);  //log base 10

      // book(_hw["udpzbefore"], "abs_u_pz_over_abs_d_pz_weighted_before", 100,-2.5,2.5);  //log base 10
      // book(_hw3["udpzbefore3"], "abs_u_pz_over_abs_d_pz_weight3_before", 100,-2.5,2.5);  //log base 10

      //same histos as above but after kinematic cuts corresponding to combination paper
      book(_h["upz"], "abs_up_pz", 100, 0.0, 100.0);
      book(_h["dpz"], "abs_down_pz", 100, 0.0, 100.0);
      book(_h["udpz"], "abs_u_pz_over_abs_d_pz", 100,-2.5,2.5);
      book(_hw2["udpz"], "abs_u_pz_over_abs_d_pz_weight2", 100,-2.5,2.5);

      // histograms of x (quark pz/ total proton energy)
      book(_h["xubefore"], "log_x_u_before",100,-3,0.0);
      book(_h["xdbefore"], "log_x_d_before",100,-3,0.0);
      // book(_hw["xubefore"], "log_x_u_weighted_before",100,-3,0.0);
      // book(_hw["xdbefore"], "log_x_d_weighted_before",100,-3,0.0);

      double abseta_bin = 2.4;

      if (CDF) {
        abseta_bin = 1.0;
      }

      // lepton rapidity for quark events, before histograms have no cuts, after have kinematic cuts corresponding to combination paper
      book(_h["eetabefore"],"e_eta_before",50,-2.5,2.5);
      book(_h["eeta"],"e_eta",50,-1*abseta_bin,abseta_bin);
      // book(_hw["eetabefore"],"e_eta_before_weighted",50,-2.5,2.5);
      // book(_hw["eeta"],"e_eta_weighted",50,-abseta_bin,abseta_bin);

      // same as above but reweighted
      book(_hw2["eeta"],"e_eta_weight2",50,-1*abseta_bin,abseta_bin);
      book(_hw2["eetabefore2"],"e_eta_weight2_before",50,-2.5,2.5);
      // book(_hw3["eeta"],"e_eta_weight3",50,-abseta_bin,abseta_bin);
      // book(_hw3["eetabefore3"],"e_eta_weight3_before",50,-2.5,2.5);

      //Ai histograms, truth level, no kinematic cuts
      book(_g["hA0"], "A0_Histogram", 25, 0.0, 50);
      book(_g["hA1"], "A1_Histogram", 25, 0.0, 50);
      book(_g["hA2"], "A2_Histogram", 25, 0.0, 50.0);
      book(_g["hA3"], "A3_Histogram", 25, 0.0, 50.0);
      book(_g["hA4"], "A4_Histogram", 25, 0.0, 50.0);
      // book(_h["hA5"], "A5_Histogram", 25, 0.0, 50.0);
      // book(_h["hA6"], "A6_Histogram", 25, 0.0, 50.0);
      // book(_h["hA7"], "A7_Histogram", 25, 0.0, 50.0);
      book(_h["hNormal"],"Normal_Histogram",25,0.0,50.0);       
            
      //angular histograms, truth level, no kinematic cuts
      book(_h["phi0_15"],"e_phi_0_15",25,0.,6.28);
      book(_h["phi15_25"],"e_phi_15_25",25,0.,6.28);
      book(_h["phi25_35"],"e_phi_25_35",25,0.,6.28);
      book(_h["phi35_65"],"e_phi_35_65",25,0., 6.28);
      book(_h["phi65_105"],"e_phi_65_105",25,0.,6.28);

      book(_h["theta0_15"],"e_costheta_0_15",25,-1.,1.);
      book(_h["theta15_25"],"e_costheta_15_25",25,-1.,1.);
      book(_h["theta25_35"],"e_costheta_25_35",25,-1.,1.);
      book(_h["theta35_65"],"e_costheta_35_65",25,-1.,1.);
      book(_h["theta65_105"],"e_costheta_65_105",25,-1.,1.);

      //reconstructed Ais, no kinematic cuts
      // filled twice - once with each branch of the equation required to reconstructed neutrino pz
      book(_g["hA0reco"], "A0_Histogram_reconstructed", 25, 0.0, 50);
      book(_g["hA1reco"], "A1_Histogram_reconstructed", 25, 0.0, 50);
      book(_g["hA2reco"], "A2_Histogram_reconstructed", 25, 0.0, 50.0);
      book(_g["hA3reco"], "A3_Histogram_reconstructed", 25, 0.0, 50.0);
      book(_g["hA4reco"], "A4_Histogram_reconstructed", 25, 0.0, 50.0);
      book(_h["hNormalreco"],"Normal_Histogram_reconstructed",25,0.0,50.0);

      //reconstructed Ais, no kinematic cuts
      // filled only with high neutrino pz solution
      book(_g["hA0recohigh"], "A0_Histogram_reconstructed_high_solution", 25, 0.0, 50);
      book(_g["hA1recohigh"], "A1_Histogram_reconstructed_high_solution", 25, 0.0, 50);
      book(_g["hA2recohigh"], "A2_Histogram_reconstructed_high_solution", 25, 0.0, 50.0);
      book(_g["hA3recohigh"], "A3_Histogram_reconstructed_high_solution", 25, 0.0, 50.0);
      book(_g["hA4recohigh"], "A4_Histogram_reconstructed_high_solution", 25, 0.0, 50.0);
      book(_h["hNormalrecohigh"],"Normal_Histogram_reconstructed_high_solution",25,0.0,50.0);

      //reconstructed phi, costheta, no cuts
      book(_h["phireco0_15"],"e_phi_0_15_reco",25,0.,6.28);
      book(_h["phireco15_25"],"e_phi_15_25_reco",25,0.,6.28);
      book(_h["phireco25_35"],"e_phi_25_35_reco",25,0.,6.28);
      book(_h["phireco35_65"],"e_phi_35_65_reco",25,0., 6.28);
      book(_h["phireco65_105"],"e_phi_65_105_reco",25,0.,6.28);

      book(_h["thetareco0_15"],"e_abscostheta_0_15_reco",25,0.,1.);
      book(_h["thetareco15_25"],"e_abscostheta_15_25_reco",25,0.,1.);
      book(_h["thetareco25_35"],"e_abscostheta_25_35_reco",25,0.,1.);
      book(_h["thetareco35_65"],"e_abscostheta_35_65_reco",25,0.,1.);
      book(_h["thetareco65_105"],"e_abscostheta_65_105_reco",25,0.,1.);

      //reconstructed Ais, with cuts corresponding to Errede and Strogolas (arXiv:hep-ph/0503291)
      book(_g["hA0reco_2"], "A0_Histogram_reconstructed_25", 25, 0.0, 50);
      book(_g["hA1reco_2"], "A1_Histogram_reconstructed_25", 25, 0.0, 50);
      book(_g["hA2reco_2"], "A2_Histogram_reconstructed_25", 25, 0.0, 50.0);
      book(_g["hA3reco_2"], "A3_Histogram_reconstructed_25", 25, 0.0, 50.0);
      book(_g["hA4reco_2"], "A4_Histogram_reconstructed_25", 25, 0.0, 50.0);
      book(_h["hNormalreco_2"],"Normal_Histogram_reconstructed_25",25,0.0,50.0);

      //reconstructed phi, costheta, with cuts corresponding to Errede and Strogolas (arXiv:hep-ph/0503291)
      book(_h["phireco0_15_2"],"e_phi_0_15_2_reco",25,0.,6.28);
      book(_h["phireco15_25_2"],"e_phi_15_25_2_reco",25,0.,6.28);
      book(_h["phireco25_35_2"],"e_phi_25_35_2_reco",25,0.,6.28);
      book(_h["phireco35_65_2"],"e_phi_35_65_2_reco",25,0., 6.28);
      book(_h["phireco65_105_2"],"e_phi_65_105_2_reco",25,0.,6.28);

      book(_h["thetareco0_15_2"],"e_abscostheta_0_15_2_reco",25,0.,1.);
      book(_h["thetareco15_25_2"],"e_abscostheta_15_25_2_reco",25,0.,1.);
      book(_h["thetareco25_35_2"],"e_abscostheta_25_35_2_reco",25,0.,1.);
      book(_h["thetareco35_65_2"],"e_abscostheta_35_65_2_reco",25,0.,1.);
      book(_h["thetareco65_105_2"],"e_abscostheta_65_105_2_reco",25,0.,1.);

      // same costheta as above but PDF reweighted for CDF
      book(_hw2["thetareco0_15_2"],"e_abscostheta_0_15_2_reco_weight2",25,0.,1.);
      book(_hw2["thetareco15_25_2"],"e_abscostheta_15_25_2_reco_weight2",25,0.,1.);
      book(_hw2["thetareco25_35_2"],"e_abscostheta_25_35_2_reco_weight2",25,0.,1.);
      book(_hw2["thetareco35_65_2"],"e_abscostheta_35_65_2_reco_weight2",25,0.,1.);
      book(_hw2["thetareco65_105_2"],"e_abscostheta_65_105_2_reco_weight2",25,0.,1.);

      //reconstructed Ais, with cuts corresponding to Errede and Strogolas (arXiv:hep-ph/0503291)
      // filled only with high neutrino pz solution
      book(_g["hA0recohigh_2"], "A0_Histogram_reconstructed_high_solution_25", 25, 0.0, 50);
      book(_g["hA1recohigh_2"], "A1_Histogram_reconstructed_high_solution_25", 25, 0.0, 50);
      book(_g["hA2recohigh_2"], "A2_Histogram_reconstructed_high_solution_25", 25, 0.0, 50.0);
      book(_g["hA3recohigh_2"], "A3_Histogram_reconstructed_high_solution_25", 25, 0.0, 50.0);
      book(_g["hA4recohigh_2"], "A4_Histogram_reconstructed_high_solution_25", 25, 0.0, 50.0);
      book(_h["hNormalrecohigh_2"],"Normal_Histogram_reconstructed_high_solution_25",25,0.0,50.0);

      //Ai histograms, truth level,  with cuts corresponding to Errede and Strogolas (arXiv:hep-ph/0503291)
      book(_g["hA0_2"], "A0_Histogram_2", 25, 0.0, 50);
      book(_g["hA1_2"], "A1_Histogram_2", 25, 0.0, 50);
      book(_g["hA2_2"], "A2_Histogram_2", 25, 0.0, 50.0);
      book(_g["hA3_2"], "A3_Histogram_2", 25, 0.0, 50.0);
      book(_g["hA4_2"], "A4_Histogram_2", 25, 0.0, 50.0);
      // book(_h["hA5_2"], "A5_Histogram_2", 25, 0.0, 50.0);
      // book(_h["hA6_2"], "A6_Histogram_2", 25, 0.0, 50.0);
      // book(_h["hA7_2"], "A7_Histogram_2", 25, 0.0, 50.0);
      book(_h["hNormal_2"],"Normal_Histogram_2",25,0.0,50.0);

      //angular histograms, truth level, with cuts corresponding to Errede and Strogolas (arXiv:hep-ph/0503291)
      book(_h["phi15_25_2"],"e_phi_15_25_2",25,0.,6.28);
      book(_h["phi25_35_2"],"e_phi_25_35_2",25,0.,6.28);
      book(_h["phi35_65_2"],"e_phi_35_65_2",25,0., 6.28);
      book(_h["phi65_105_2"],"e_phi_65_105_2",25,0.,6.28);

      book(_h["theta15_25_2"],"e_costheta_15_25_2",25,-1.,1.);
      book(_h["theta25_35_2"],"e_costheta_25_35_2",25,-1.,1.);
      book(_h["theta35_65_2"],"e_costheta_35_65_2",25,-1.,1.);
      book(_h["theta65_105_2"],"e_costheta_65_105_2",25,-1.,1.);

      //Ais in lepton rapidity bins, truth level, no kinematic cuts
      book(_g["hA0_12_less"], "A0_Histogram_less_than_-1.2", 25, 0.0, 50);
      book(_g["hA1_12_less"], "A1_Histogram_less_than_-1.2", 25, 0.0, 50);
      book(_g["hA2_12_less"], "A2_Histogram_less_than_-1.2", 25, 0.0, 50.0);
      book(_g["hA3_12_less"], "A3_Histogram_less_than_-1.2", 25, 0.0, 50.0);
      book(_g["hA4_12_less"], "A4_Histogram_less_than_-1.2", 25, 0.0, 50.0);

      book(_g["hA0_12_0"], "A0_Histogram_-1.2_0", 25, 0.0, 50);
      book(_g["hA1_12_0"], "A1_Histogram_-1.2_0", 25, 0.0, 50);
      book(_g["hA2_12_0"], "A2_Histogram_-1.2_0", 25, 0.0, 50.0);
      book(_g["hA3_12_0"], "A3_Histogram_-1.2_0", 25, 0.0, 50.0);
      book(_g["hA4_12_0"], "A4_Histogram_-1.2_0", 25, 0.0, 50.0);

      book(_g["hA0_0_12"], "A0_Histogram_0_1.2", 25, 0.0, 50);
      book(_g["hA1_0_12"], "A1_Histogram_0_1.2", 25, 0.0, 50);
      book(_g["hA2_0_12"], "A2_Histogram_0_1.2", 25, 0.0, 50.0);
      book(_g["hA3_0_12"], "A3_Histogram_0_1.2", 25, 0.0, 50.0);
      book(_g["hA4_0_12"], "A4_Histogram_0_1.2", 25, 0.0, 50.0);

      book(_g["hA0_12_gtr"], "A0_Histogram_greater_than_1.2", 25, 0.0, 50);
      book(_g["hA1_12_gtr"], "A1_Histogram_greater_than_1.2", 25, 0.0, 50);
      book(_g["hA2_12_gtr"], "A2_Histogram_greater_than_1.2", 25, 0.0, 50.0);
      book(_g["hA3_12_gtr"], "A3_Histogram_greater_than_1.2", 25, 0.0, 50.0);
      book(_g["hA4_12_gtr"], "A4_Histogram_greater_than_1.2", 25, 0.0, 50.0);

      // // histos for 13.6TeV LHC data, binned by PTW
      // // histograms with names ending in 'cuts' have cuts corresponding to A Sparrow's PhD thesis.

      // 2d histograms to show correlation of 2Lp-1 with truth costheta in the collins-soper frame
      book(_h2["Lp correlation 50GeV truth"],"Lp_correlation_50_truth", 28, -1.4,1.4,20,-1.,1.);
      book(_h2["Lp correlation 50GeV reconstructed"],"Lp_correlation_50_reconstructed", 28, -1.4,1.4,20,-1.,1.);
      book(_h2["Lp correlation 50GeV reconstructed cuts"],"Lp_correlation_50_reconstructed_cuts", 28, -1.4,1.4,20,-1.,1.);
      book(_h2["Lp correlation 100GeV truth"],"Lp_correlation_100_truth", 28, -1.4,1.4,20,-1.,1.);
      book(_h2["Lp correlation 100GeV reconstructed"],"Lp_correlation_100_reconstructed", 28, -1.4,1.4,20,-1.,1.);
      book(_h2["Lp correlation 100GeV reconstructed cuts"],"Lp_correlation_100_reconstructed_cuts", 28, -1.4,1.4,20,-1.,1.);
      book(_h2["Lp correlation 200GeV truth"],"Lp_correlation_200_truth", 28, -1.4,1.4,20,-1.,1.);
      book(_h2["Lp correlation 200GeV reconstructed"],"Lp_correlation_200_reconstructed", 28, -1.4,1.4,20,-1.,1.);
      book(_h2["Lp correlation 200GeV reconstructed cuts"],"Lp_correlation_200_reconstructed_cuts", 28, -1.4,1.4,20,-1.,1.);

      //histograms of 2Lp-1
      book(_h["Lp histogram 50GeV truth"],"Lp_histogram_50_truth", 25, -1.,1.);
      book(_h["Lp histogram 50GeV reconstructed"],"Lp_histogram_50_reconstructed", 25, -1.,1.);
      book(_h["Lp histogram 50GeV reconstructed cuts"],"Lp_histogram_50_reconstructed_cuts", 25, -1.,1.);
      book(_h["Lp histogram 100GeV truth"],"Lp_histogram_100_truth", 25, -1.,1.);
      book(_h["Lp histogram 100GeV reconstructed"],"Lp_histogram_100_reconstructed", 25, -1.,1.);
      book(_h["Lp histogram 100GeV reconstructed cuts"],"Lp_histogram_100_reconstructed_cuts", 25, -1.,1.);
      book(_h["Lp histogram 200GeV truth"],"Lp_histogram_200_truth", 25, -1.,1.);
      book(_h["Lp histogram 200GeV reconstructed"],"Lp_histogram_200_reconstructed", 25, -1.,1.);
      book(_h["Lp histogram 200GeV reconstructed cuts"],"Lp_histogram_200_reconstructed_cuts", 25, -1.,1.);
      
      // 2d histograms for correlation with truth costheta in the helicity frame
      book(_h2["Lp helicity correlation 50GeV truth"],"Lp_helicity_correlation_50_truth", 28, -1.4,1.4,20,-1.,1.);
      book(_h2["Lp helicity correlation 50GeV reconstructed"],"Lp_helicity_correlation_50_reconstructed", 28, -1.4,1.4,20,-1.,1.);
      book(_h2["Lp helicity correlation 50GeV reconstructed cuts"],"Lp_helicity_correlation_50_reconstructed_cuts", 28, -1.4,1.4,20,-1.,1.);
      book(_h2["Lp helicity correlation 100GeV truth"],"Lp_helicity_correlation_100_truth", 28, -1.4,1.4,20,-1.,1.);
      book(_h2["Lp helicity correlation 100GeV reconstructed"],"Lp_helicity_correlation_100_reconstructed", 28, -1.4,1.4,20,-1.,1.);
      book(_h2["Lp helicity correlation 100GeV reconstructed cuts"],"Lp_helicity_correlation_100_reconstructed_cuts", 28, -1.4,1.4,20,-1.,1.);
      book(_h2["Lp helicity correlation 200GeV truth"],"Lp_helicity_correlation_200_truth", 28, -1.4,1.4,20,-1.,1.);
      book(_h2["Lp helicity correlation 200GeV reconstructed"],"Lp_helicity_correlation_200_reconstructed", 28, -1.4,1.4,20,-1.,1.);
      book(_h2["Lp helicity correlation 200GeV reconstructed cuts"],"Lp_helicity_correlation_200_reconstructed_cuts", 28, -1.4,1.4,20,-1.,1.);

      //weight counters
      book(event_counter,"_event_counter");
      // book(counter_1,"_weight_1_counter");
      book(counter_2,"_weight_2_counter");
      // book(counter_3,"_weight_3_counter");

    }

    //Breit-Wigner
    double bw(double m, double mp, double ww) {
      double ww2 = ww*ww;
      double mp2 = mp*mp;
      double m2 = m*m;
      double d = m2 - mp2;

      return 1/(d*d + m2*m2*ww2/mp2);
    }

    // //Calculate event weight
    // double event_weighting(double logratio, int charge) {
    //   double weight = 0;
    //   if (charge == 1) {
    //     weight = exp( -1 * (3.52 * logratio + 5.92) ) - 0.585*(1 + tanh(2.15 * logratio - 3.65)) + 1.25;
    //   } else {
    //     weight = exp( -1 * (3.28 * logratio + 5.46) ) - 0.558*(1 + tanh(2.35 * logratio - 4.01)) + 1.22;
    //   }

    //   return 1 / weight;
    // }

    //second candidate event weight
    double weighting_2(double logratio, int charge) {
      double weight = 0;
      if (charge == 1) {
        weight = exp( -1 * (3.52 * logratio + 5.92) ) + 1.25;
      } else {
        weight = exp( -1 * (3.28 * logratio + 5.46) ) + 1.22;
      }

      return 1 / weight;
    }

    // //third candidate event weight
    // double weighting_3(double logratio, int charge) {
      
    //   double weight = 0.8 * (logratio - 0.75) * (logratio - 0.75) * (logratio + 1.2) * (logratio + 1.2) + 1.25;

    //   return 1 / weight;
    // }

    #include "AngularFunction.cc"
      
    // Per event analysis
    void analyze(const Event& event) {
      MSG_INFO("(new event)");

      MSG_INFO("evt weight = " + std::to_string(event.weights()[0]));

      const double ww = 2.085*GeV;
      const double mw_ref = 79.82436*GeV; //value used in param cards for MC

      int ncount = 0;
      int ecount = 0;
      int charge = 0;

      if (CDF) {
        MSG_INFO("CDF");
      } else {
        MSG_INFO("ATLAS");
      }

      int qcount = 0;
      int ucount = 0;
      int dcount = 0;

      double u_pz = 0.;
      double d_pz = 0.;

      int q1pid = 0;
      double q1pz = 0;

      FourMomentum nFourMomentum;
      FourMomentum eFourMomentum;
      FourMomentum pTmissA;

      //find the neutrinos (either n or N)
      for( const Particle &  n : apply<IdentifiedFinalState>(event, "neut").particlesByPt()) {
        ncount += 1 ;
        nFourMomentum = n.momentum();
      }

      for( const Particle & N : apply<IdentifiedFinalState>(event, "sneut").particlesByPt()) {
        ncount += 1;
        nFourMomentum = N.momentum();
      }

        
      //determine 'true' missing energy
      for ( const Particle & p : apply<VisibleFinalState>(event, "vfs").particles() ) {
        pTmissA -= p.momentum();
      } 

        
      // find the electrons
      std::random_device rd;
      std::mt19937 gen(rd());           
      for( const Particle & e : apply<IdentifiedFinalState>(event, "elecs").particlesByPt()) {
        eFourMomentum  = e.momentum();
        ecount += 1;
        charge = e.charge();


        //loop through all ancestors which are quarks (quarks have abs(charge) 1/3 or 2/3)
        for( Particle q : e.ancestors( Cuts::abscharge3 == 2 || Cuts::abscharge3 == 1 , false) ) {

          //choose only ancestors which have no quark children
          if ( q.hasChildWith( Cuts::abscharge3 == 2 || Cuts::abscharge3 == 1 ) == false ) {

            //count total number of quarks found just in case (should be two)
            qcount += 1;

            //check interaction is ud -> w
            if (abs(q.pid()) == 2) {
              ucount++;
              u_pz = q.momentum().pz();
            } else if (abs(q.pid()) == 1) {
              dcount++;
              d_pz = q.momentum().pz();
            }

            if (qcount == 1) {
              q1pid = q.pid();
              q1pz = q.momentum().pz();
            }
            
            
          }
        }
      }

      // // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      // Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 50*GeV);

      // MSG_INFO("found " + std::to_string(jets.size()) + " jets");
      // for (const Jet& j : jets) {
      //   MSG_INFO("jet pt " + std::to_string(j.pt()));
      // }

      if (ecount == 1 && qcount !=2) {
        MSG_INFO("found electron but not two quarks!");
      }



      MSG_INFO("ecount = " + std::to_string(ecount));
      MSG_INFO("ncount = " + std::to_string(ncount));

      //only analyses events with one electron and one neutrino
      if ( ecount == 1 && ncount == 1) {
        double e_pT = eFourMomentum.pT();
        double e_abseta = eFourMomentum.abseta();
        double e_eta = eFourMomentum.eta();


        //smearing parameters
        double recoil_scale = 0.;
        float a = 0.;
        float b = 0.;
        float c = 0.;
        float underlying_energy = 0.;
        double C = 0.;

        if (CDF){
          recoil_scale = 0.79;
          a = 0.632;
          b = 17;
          c = 3;
          underlying_energy = 5.17;  
          C = 0.009;
        } else {
          //Atlas
          recoil_scale = 1.;
          a = 0.75;
          b = 0.002;
          c = 2;
          underlying_energy = 19;  
          C = 0.011;
        }
        double sigmaP = sqrt(0.126*0.126*e_pT + C*C*e_pT*e_pT);

        std::normal_distribution<double> dist(0, sigmaP);

        double dE = dist(gen);

        //smears tranvserse momentum
        double e_pTD = e_pT + dE;

        //create four vector for reconstructed neutrino mometum
        //looks strange but first entry in array is high-pz solution, second is low-pz solution
        std::array<FourMomentum, 2> missReconstructed;

        double e_pTD_x = eFourMomentum.px() + dE * cos(eFourMomentum.phi());
        double e_pTD_y = eFourMomentum.py() + dE * sin(eFourMomentum.phi());
        double e_pTD_z = eFourMomentum.pz() + dE / tan(eFourMomentum.theta());

        //create 4-vector for reconstructed charged lepton momentum
        FourMomentum eReconstructed;
        eReconstructed.setPM(e_pTD_x, e_pTD_y, e_pTD_z, 0.);
                    
        double recoil_x =  - eFourMomentum.px() - pTmissA.px();
        double recoil_y = - eFourMomentum.py() - pTmissA.py();
        double E_T_U = sqrt(recoil_x*recoil_x + recoil_y*recoil_y);


        double sigmaU = recoil_scale*sqrt(E_T_U);

        double R = a*log(b*E_T_U + c)/log(b*15 + c);
        std::normal_distribution<double> dist2(0, sigmaU);
        double E_T_UReco1 = E_T_U*R + dist2(gen);
        double E_T_U_XReco1 = 0;
        double E_T_U_YReco1 = 0;

        if (E_T_U != 0) {
          E_T_U_XReco1 = E_T_UReco1*recoil_x/(sqrt(recoil_x*recoil_x + recoil_y*recoil_y));
          E_T_U_YReco1 = E_T_UReco1*recoil_y/(sqrt(recoil_x*recoil_x + recoil_y*recoil_y));
        }

        std::normal_distribution<double> dist3(0, underlying_energy);
        double E_T_U_ue = dist3(gen);

        //Random angle for underlying energy
        std::uniform_real_distribution<> dist4(0, 2*M_PI);
        double angle = dist4(gen);
        double E_T_U_X_ue = E_T_U_ue*cos(angle);
        double E_T_U_Y_ue = E_T_U_ue*sin(angle);
          
        //double E_T_UReco = sqrt((E_T_U_X_ue + E_T_U_XReco1)*(E_T_U_X_ue + E_T_U_XReco1) + (E_T_U_Y_ue + E_T_U_YReco1)*(E_T_U_Y_ue + E_T_U_YReco1));    
        double E_T_U_XReco = E_T_U_X_ue + E_T_U_XReco1;
        double E_T_U_YReco = E_T_U_Y_ue + E_T_U_YReco1;

        //reconstruct neutrino pT (MET)
        double pxmiss_reco = - E_T_U_XReco - e_pTD_x;
        double pymiss_reco = - E_T_U_YReco - e_pTD_y;
        double MET = sqrt( pxmiss_reco * pxmiss_reco + pymiss_reco * pymiss_reco );

        //reconstruct neutrino pz (there are two solutions)
        //returns nan quite often, be careful
        //equation is in arXiv:hep-ph/0503291 but is a factor of two out
        double xi2 = 0.5 * mw_ref * mw_ref + e_pTD_x * pxmiss_reco + e_pTD_y * pymiss_reco; //intermediate quantity in calculating reconstructed pz_miss
        double Ee_reco = sqrt(e_pTD * e_pTD + e_pTD_z * e_pTD_z);
        double pzmiss_reco_high = ( xi2*e_pTD_z + Ee_reco*sqrt( xi2*xi2 - e_pTD*e_pTD*MET*MET ) ) / ( e_pTD*e_pTD );
        double pzmiss_reco_low = ( xi2*e_pTD_z - Ee_reco*sqrt( xi2*xi2 - e_pTD*e_pTD*MET*MET ) ) / ( e_pTD*e_pTD );

        // MSG_INFO("pzmiss low = " + std::to_string(pzmiss_reco_low));
        // MSG_INFO("pzmiss high = " + std::to_string(pzmiss_reco_high));

        //fill reco neutrino four momenta
        missReconstructed[0].setPM(pxmiss_reco,pymiss_reco,pzmiss_reco_high, 0.);
        missReconstructed[1].setPM(pxmiss_reco,pymiss_reco,pzmiss_reco_low, 0.);

        //Reconstruct mT
        double mT = e_pTD*MET - e_pTD_x*pxmiss_reco - e_pTD_y*pymiss_reco ;

        mT = sqrt(2 * mT);

        // MSG_INFO("e_pTD = " + std::to_string(e_pTD/GeV));
        // MSG_INFO("e_abseta = " + std::to_string(e_abseta/GeV));
        // MSG_INFO("MET = " + std::to_string(MET/GeV));
        // MSG_INFO("mT = " + std::to_string(mT/GeV));
        // MSG_INFO("E_T_U = " + std::to_string(E_T_U/GeV));
        // MSG_INFO("---");

        //plot cut criteria values
        _h["e pt before"]->fill(e_pTD);
        _h["mt before"]->fill(mT);
        _h["pt miss before"]->fill(MET);
        _h["e abseta before"]->fill(e_abseta);
        _h["ETU"]->fill(E_T_U);

        double logratio = 0.;
        // double weight = 1.;
        double weight_2 = 1.;
        // double weight_3 = 1.;

        // analyse quark stuff for only ud events, and only event where both quarks are valence
        // for CDF, valence strongly indicated when quark pid sign matches pz sign
        // for Atlas, valence strongly indicated when quark with positive pid has charge matching w
        bool quarkevent = (qcount ==2 && ucount ==1 && dcount ==1 && q1pid * q1pz > 0);

        if (!CDF) {
          // Atlas quark weighting event criteria
          quarkevent = (qcount==2 && ucount==1 && dcount ==1);
        }

        if ( quarkevent ) {
          logratio = log10( std::abs( u_pz / d_pz ) );
          // weight = event_weighting(logratio,charge);
          weight_2 = weighting_2(logratio,charge);
          // weight_3 = weighting_3(logratio,charge);

          _h["upzbefore"]->fill( std::abs(u_pz) );
          _h["dpzbefore"]->fill( std::abs(d_pz) );
          _h["udpzbefore"]->fill( logratio );
          // _hw["udpzbefore"]->fill( logratio , weight );
          _hw2["udpzbefore2"]->fill( logratio , weight_2 );
          // _hw3["udpzbefore3"]->fill( logratio , weight_3 );

          //CM energy
          double Ebeam = beamEnergies().first;

          //fill pdf histos
          double u_x = log10( std::abs(u_pz) / Ebeam );
          double d_x = log10( std::abs(d_pz) / Ebeam );
          _h["xubefore"]->fill( u_x );
          _h["xdbefore"]->fill( d_x );
          // _hw["xubefore"]->fill( u_x , weight );
          // _hw["xdbefore"]->fill( d_x , weight );
        }
        

        _h["eetabefore"]->fill(e_eta);
        // _hw["eetabefore"]->fill(e_eta,weight);
        _hw2["eetabefore2"]->fill(e_eta,weight_2);
        // _hw3["eetabefore3"]->fill(e_eta,weight_3);




        // MSG_INFO("* " + std::to_string(std::abs(u_pz / d_pz)));
        // MSG_INFO("  quark 1, id: " + std::to_string(q1pid));
        // MSG_INFO("  quark 1, pz: " + std::to_string(q1pz));
        // MSG_INFO("  quark 2, pid: " + std::to_string(q2pid));
        // MSG_INFO("  quark 2, pz: " + std::to_string(q2pz));


        double m_inv = eFourMomentum.p() * nFourMomentum.p() - eFourMomentum.px() * nFourMomentum.px() - eFourMomentum.py() * nFourMomentum.py() - eFourMomentum.pz() * nFourMomentum.pz() ;
        m_inv = sqrt(2 * m_inv);

        double m_inv_reco = 0.;

        //find angular coefficients before selection

        double costheta = 0;
        double sintheta = 0;
        double cosphi = 0;
        double sinphi = 0;

        // MSG_INFO("(calculating truth level angular coefficients)");
        //all angles in Collins-Soper frame
        CalculateAngles(eFourMomentum, nFourMomentum, m_inv, charge, costheta, sintheta, cosphi, sinphi);
        double sin2theta = 2. * sintheta * costheta;
        double sin2phi = 2. * sinphi * cosphi;
        double cos2phi = cosphi * cosphi - sinphi * sinphi;

        double a0 = 10.*(1. - 3. * costheta * costheta)/3. + 2. /3. ;
        double a1 = 5. * sin2theta * cosphi;
        double a2 = 10. * sintheta * sintheta * cos2phi;
        double a3 = 4. * sintheta * cosphi;
        double a4 = 4. * costheta;
        // double a5 = 5. * sintheta * sintheta * sin2phi;
        // double a6 = 5. * sin2theta * sinphi;
        // double a7 = 4. * sintheta * sinphi;

        double w_px = eFourMomentum.px() + nFourMomentum.px();
        double w_py = eFourMomentum.py() + nFourMomentum.py();
        double w_pz = eFourMomentum.pz() + nFourMomentum.pz();
        double P_T_W_truth =  sqrt(w_px*w_px + w_py*w_py);
        MSG_INFO("w pt is " + std::to_string(P_T_W_truth));

        //calculate w transvserse momentum
        double w_px_reco = missReconstructed[0].px() + eReconstructed.px();
        double w_py_reco = missReconstructed[0].py() + eReconstructed.py();
        double P_T_W_reco = sqrt( w_px_reco*w_px_reco + w_py_reco*w_py_reco );

        _g["hA0"]->fill(P_T_W_truth,a0);
        _g["hA1"]->fill(P_T_W_truth,a1);
        _g["hA2"]->fill(P_T_W_truth,a2);
        _g["hA3"]->fill(P_T_W_truth,a3);
        _g["hA4"]->fill(P_T_W_truth,a4);
        // _h["hA5"]->fill(P_T_W_truth,a5);
        // _h["hA6"]->fill(P_T_W_truth,a6);
        // _h["hA7"]->fill(P_T_W_truth,a7);
        _h["hNormal"]->fill(P_T_W_truth);

        //find phi to plot
        double phi = 0;

        if (sinphi > 0) {
          phi = acos(cosphi);
        } else {
          phi = 2 * pi - acos(cosphi);
        }
        
        // if (tan(phi) < 0.11) {
        //   MSG_INFO("phi is less than 0.11!");
        // }

        // if (phi > 3.015 && phi < 3.268) {
        //   MSG_INFO("phi is in pi bin!");
        // }


        if (P_T_W_truth < 15) {
          _h["phi0_15"]->fill(phi);
          _h["theta0_15"]->fill(costheta);
        } else if (P_T_W_truth < 25) {
          _h["phi15_25"]->fill(phi);
          _h["theta15_25"]->fill(costheta);
        } else if (P_T_W_truth < 35) {
          _h["phi25_35"]->fill(phi);
          _h["theta25_35"]->fill(costheta);
        } else if (P_T_W_truth < 65) {
          _h["phi35_65"]->fill(phi);
          _h["theta35_65"]->fill(costheta);
        } else if (P_T_W_truth < 105) {
          _h["phi65_105"]->fill(phi);
          _h["theta65_105"]->fill(costheta);
        }

        //Ais in lepton rapidity bins
        if (e_eta < -1.2) {
          _g["hA0_12_less"]->fill(P_T_W_truth,a0);
          _g["hA1_12_less"]->fill(P_T_W_truth,a1);
          _g["hA2_12_less"]->fill(P_T_W_truth,a2);
          _g["hA3_12_less"]->fill(P_T_W_truth,a3);
          _g["hA4_12_less"]->fill(P_T_W_truth,a4);
        } else if (e_eta < 0) {
          _g["hA0_12_0"]->fill(P_T_W_truth,a0);
          _g["hA1_12_0"]->fill(P_T_W_truth,a1);
          _g["hA2_12_0"]->fill(P_T_W_truth,a2);
          _g["hA3_12_0"]->fill(P_T_W_truth,a3);
          _g["hA4_12_0"]->fill(P_T_W_truth,a4);
        } else if (e_eta < 1.2) {
          _g["hA0_0_12"]->fill(P_T_W_truth,a0);
          _g["hA1_0_12"]->fill(P_T_W_truth,a1);
          _g["hA2_0_12"]->fill(P_T_W_truth,a2);
          _g["hA3_0_12"]->fill(P_T_W_truth,a3);
          _g["hA4_0_12"]->fill(P_T_W_truth,a4);
        } else {
          _g["hA0_12_gtr"]->fill(P_T_W_truth,a0);
          _g["hA1_12_gtr"]->fill(P_T_W_truth,a1);
          _g["hA2_12_gtr"]->fill(P_T_W_truth,a2);
          _g["hA3_12_gtr"]->fill(P_T_W_truth,a3);
          _g["hA4_12_gtr"]->fill(P_T_W_truth,a4);
        }

        if (lhcrun3) {

          //lepton projection variable (see Measurement of the Polarisation of the W Boson and Application to Supersymmetry Searches at the Large Hadron Collider, Alexander Sparrow PhD thesis)
          //note the two reconstructed correlation histograms are compared to truth costheta, as only |costheta| can be reconstructed
          double Lp_truth = (eFourMomentum.px()*w_px + eFourMomentum.py()*w_py) / (P_T_W_truth*P_T_W_truth);
          double Lp_reco = (eReconstructed.px()*w_px_reco + eReconstructed.py()*w_py_reco) / (P_T_W_reco*P_T_W_reco);
          double w_p = sqrt(w_px*w_px + w_py+w_py + w_pz*w_pz);
          double w_e = eFourMomentum.E() + nFourMomentum.E();
          double costheta_helicityframe = (2*eFourMomentum.p() - w_e)/w_p;
          
          if (P_T_W_truth > 50) {


            _h2["Lp correlation 50GeV truth"]->fill( 2*Lp_truth - 1,costheta );
            _h2["Lp helicity correlation 50GeV truth"]->fill( 2*Lp_truth - 1,costheta_helicityframe );
            _h["Lp histogram 50GeV truth"]->fill( 2*Lp_truth - 1 );

            if (P_T_W_truth > 100) {
              _h2["Lp correlation 100GeV truth"]->fill( 2*Lp_truth - 1,costheta );
              _h2["Lp helicity correlation 100GeV truth"]->fill( 2*Lp_truth - 1,costheta_helicityframe );
              _h["Lp histogram 100GeV truth"]->fill( 2*Lp_truth - 1 );

              if (P_T_W_truth > 200){
                _h2["Lp correlation 200GeV truth"]->fill( 2*Lp_truth - 1,costheta );
                _h2["Lp helicity correlation 200GeV truth"]->fill( 2*Lp_truth - 1,costheta_helicityframe );
                _h["Lp histogram 200GeV truth"]->fill( 2*Lp_truth - 1 );
              }
            }
          }

        

          if (P_T_W_reco > 50) {
            _h2["Lp correlation 50GeV reconstructed"]->fill( 2*Lp_reco - 1,costheta );
            _h2["Lp helicity correlation 50GeV reconstructed"]->fill( 2*Lp_reco - 1,costheta_helicityframe );
            _h["Lp histogram 50GeV reconstructed"]->fill( 2*Lp_reco - 1 );

            if (P_T_W_reco > 100) {
              _h2["Lp correlation 100GeV reconstructed"]->fill( 2*Lp_reco - 1,costheta );
              _h2["Lp helicity correlation 100GeV reconstructed"]->fill( 2*Lp_reco - 1,costheta_helicityframe );
              _h["Lp histogram 100GeV reconstructed"]->fill( 2*Lp_reco - 1 );

              if (P_T_W_reco > 200){
                _h2["Lp correlation 200GeV reconstructed"]->fill( 2*Lp_reco - 1,costheta );
                _h2["Lp helicity correlation 200GeV reconstructed"]->fill( 2*Lp_reco - 1,costheta_helicityframe );
                _h["Lp histogram 200GeV reconstructed"]->fill( 2*Lp_reco - 1 );
              }
            }
          }

          if ( e_pTD > 25 && e_abseta < 2.4 && mT > 50) {
            if (P_T_W_reco > 50) {
              _h2["Lp correlation 50GeV reconstructed cuts"]->fill( 2*Lp_reco - 1,costheta );
              _h2["Lp helicity correlation 50GeV reconstructed cuts"]->fill( 2*Lp_reco - 1,costheta_helicityframe );
              _h["Lp histogram 50GeV reconstructed cuts"]->fill( 2*Lp_reco - 1 );

              if (P_T_W_reco > 100) {
                _h2["Lp correlation 100GeV reconstructed cuts"]->fill( 2*Lp_reco - 1,costheta );
                _h2["Lp helicity correlation 100GeV reconstructed cuts"]->fill( 2*Lp_reco - 1,costheta_helicityframe );
                _h["Lp histogram 100GeV reconstructed cuts"]->fill( 2*Lp_reco - 1 );

                if (P_T_W_reco > 200){
                  _h2["Lp correlation 200GeV reconstructed cuts"]->fill( 2*Lp_reco - 1,costheta );
                  _h2["Lp helicity correlation 200GeV reconstructed cuts"]->fill( 2*Lp_reco - 1,costheta_helicityframe );
                  _h["Lp histogram 200GeV reconstructed cuts"]->fill( 2*Lp_reco - 1 );
                }
              }
            }
          }
        
        }


        //find and plot reconstructed angular coefficients (fill twice, once for each solution to pzmiss quadratic)
        //thus for loop as each entry in vector contains a different pz solution.
        if (!isnan(missReconstructed[0].pz())) {

          for (int i=0; i<2;i++) {

            //only calculate reconstructed invariant mass if used in further calculations
            m_inv_reco = eReconstructed.p() * missReconstructed[i].p() - eReconstructed.px() * missReconstructed[i].px() - eReconstructed.py() * missReconstructed[i].py() - eReconstructed.pz() * missReconstructed[i].pz() ;
            m_inv_reco = sqrt( 2*m_inv_reco );


            double costheta_reco = 0.;
            double sintheta_reco = 0.;
            double cosphi_reco = 0.;
            double sinphi_reco = 0.;
            CalculateAngles(eReconstructed, missReconstructed[i], m_inv_reco, charge, costheta_reco, sintheta_reco, cosphi_reco, sinphi_reco);
            double cos2phi_reco = cosphi_reco * cosphi_reco - sinphi_reco * sinphi_reco;
            double sin2theta_reco = 2. * sintheta_reco * costheta_reco;

            //find coefficients
            double a0reco = 10.*(1. - 3. * costheta_reco * costheta_reco)/3. + 2. /3. ;
            double a1reco = 5. * sin2theta_reco * cosphi_reco;
            double a2reco = 10. * sintheta_reco * sintheta_reco * cos2phi_reco;
            double a3reco = 4. * sintheta_reco * cosphi_reco;
            double a4reco = 4. * costheta_reco;

            double phi_reco = 0.;

            if (sinphi_reco > 0) {
              phi_reco = acos(cosphi_reco);
            } else {
              phi_reco = 2 * pi - acos(cosphi_reco);
            }

            // MSG_INFO("costheta: " + std::to_string(costheta_reco));
            // MSG_INFO("sintheta: " + std::to_string(sintheta_reco));
            // MSG_INFO("cosphi: " + std::to_string(cosphi_reco));
            // MSG_INFO("sinphi: " + std::to_string(sinphi_reco));

            // reco collins-soper angles

            _g["hA0reco"]->fill(P_T_W_reco,a0reco);
            _g["hA1reco"]->fill(P_T_W_reco,a1reco);
            _g["hA2reco"]->fill(P_T_W_reco,a2reco);
            _g["hA3reco"]->fill(P_T_W_reco,a3reco);
            _g["hA4reco"]->fill(P_T_W_reco,a4reco);
            _h["hNormalreco"]->fill(P_T_W_reco);


            if (i==0) {

              _g["hA0recohigh"]->fill(P_T_W_reco,a0reco);
              _g["hA1recohigh"]->fill(P_T_W_reco,a1reco);
              _g["hA2recohigh"]->fill(P_T_W_reco,a2reco);
              _g["hA3recohigh"]->fill(P_T_W_reco,a3reco);
              _g["hA4recohigh"]->fill(P_T_W_reco,a4reco);
              _h["hNormalrecohigh"]->fill(P_T_W_reco);

              if (P_T_W_reco < 15) {
                _h["phireco0_15"]->fill(phi_reco);
                _h["thetareco0_15"]->fill(abs(costheta_reco));
              } else if (P_T_W_reco < 25) {
                _h["phireco15_25"]->fill(phi_reco);
                _h["thetareco15_25"]->fill(abs(costheta_reco));
              } else if (P_T_W_reco < 35) {
                _h["phireco25_35"]->fill(phi_reco);
                _h["thetareco25_35"]->fill(abs(costheta_reco));
              } else if (P_T_W_reco < 65) {
                _h["phireco35_65"]->fill(phi_reco);
                _h["thetareco35_65"]->fill(abs(costheta_reco));
              } else if (P_T_W_reco < 105) {
                _h["phireco65_105"]->fill(phi_reco);
                _h["thetareco65_105"]->fill(abs(costheta_reco));
              }
            }
        

            //mirrors experimental distribution cuts in Errede and Strogolas
            if ( e_pTD > 20 && MET > 20 && e_abseta < 1 ) {

              _g["hA0reco_2"]->fill(P_T_W_reco,a0reco);
              _g["hA1reco_2"]->fill(P_T_W_reco,a1reco);
              _g["hA2reco_2"]->fill(P_T_W_reco,a2reco);
              _g["hA3reco_2"]->fill(P_T_W_reco,a3reco);
              _g["hA4reco_2"]->fill(P_T_W_reco,a4reco);
              _h["hNormalreco_2"]->fill(P_T_W_reco);

              //fill high-solution only histos

              if (i==0) {

                _g["hA0recohigh_2"]->fill(P_T_W_reco,a0reco);
                _g["hA1recohigh_2"]->fill(P_T_W_reco,a1reco);
                _g["hA2recohigh_2"]->fill(P_T_W_reco,a2reco);
                _g["hA3recohigh_2"]->fill(P_T_W_reco,a3reco);
                _g["hA4recohigh_2"]->fill(P_T_W_reco,a4reco);
                _h["hNormalrecohigh_2"]->fill(P_T_W_reco);

                if (P_T_W_reco < 15) {
                  _h["phireco0_15_2"]->fill(phi_reco);
                  _h["thetareco0_15_2"]->fill(abs(costheta_reco));
                  _hw2["thetareco0_15_2"]->fill(abs(costheta_reco),weight_2);
                } else if (P_T_W_reco < 25) {
                  _h["phireco15_25_2"]->fill(phi_reco);
                  _h["thetareco15_25_2"]->fill(abs(costheta_reco));
                  _hw2["thetareco15_25_2"]->fill(abs(costheta_reco),weight_2);
                } else if (P_T_W_reco < 35) {
                  _h["phireco25_35_2"]->fill(phi_reco);
                  _h["thetareco25_35_2"]->fill(abs(costheta_reco));
                  _hw2["thetareco25_35_2"]->fill(abs(costheta_reco),weight_2);
                } else if (P_T_W_reco < 65) {
                  _h["phireco35_65_2"]->fill(phi_reco);
                  _h["thetareco35_65_2"]->fill(abs(costheta_reco));
                  _hw2["thetareco35_65_2"]->fill(abs(costheta_reco),weight_2);
                } else if (P_T_W_reco < 105) {
                  _h["phireco65_105_2"]->fill(phi_reco);
                  _h["thetareco65_105_2"]->fill(abs(costheta_reco));
                  _hw2["thetareco65_105_2"]->fill(abs(costheta_reco),weight_2);
                }
                
              }
            }
          
          }
        
        }


        
        //only add to histos if event meets thresholds
        if ( ( CDF && ( e_pTD < 55*GeV && e_pTD > 30*GeV && e_abseta < 1 && MET > 30*GeV && MET < 55*GeV && mT > 60*GeV && mT < 100*GeV && E_T_U < 15*GeV ) ) || 
            ( (CDF == false) && ( e_pTD > 30*GeV && e_abseta < 2.4 && MET > 30*GeV && mT > 60*GeV &&  E_T_U < 30*GeV ) ) ) {



          double n_pT = nFourMomentum.pT();
          double n_px = nFourMomentum.px();
          double n_py = nFourMomentum.py();

          double e_px = eFourMomentum.px();
          double e_py = eFourMomentum.py();

          double miss_pT = pTmissA.pT();
          double miss_px = pTmissA.px();
          double miss_py = pTmissA.py();

          _h["e pt"]->fill(e_pTD);
          _h["mt"]->fill(mT);
          _h["pt miss"]->fill(MET);
          _h["e pt truth"]->fill(e_pT);
          _h["pt miss truth"]->fill(pTmissA.pT());
          _h["n pt truth"]->fill(n_pT);
          
          double mT_true =  n_pT * e_pT - n_px * e_px - n_py * e_py;
          mT_true = sqrt(2 * mT_true);

          double mT_miss = miss_pT * e_pT - miss_px * e_px - miss_py * e_py ;
          mT_miss = sqrt(2 * mT_miss);

          _h["mt truth"]->fill(mT_true);
          _h["mt e pt miss truth"]->fill(mT_miss);
          _h["m inv"]->fill(m_inv);

          _h["eeta"]->fill(e_eta);
          // _hw["eeta"]->fill(e_eta, weight);
          _hw2["eeta"]->fill(e_eta, weight_2);
          // _hw3["eeta"]->fill(e_eta, weight_3);
          
          if (quarkevent) {
            _h["upz"]->fill(u_pz);
            _h["dpz"]->fill(d_pz);
            _h["udpz"]->fill(logratio);
            _hw2["udpz"]->fill(logratio,weight_2);
          }

          _hw2["e pt"]->fill(e_pTD,weight_2);
          _hw2["mt"]->fill(mT,weight_2);
          _hw2["pt miss"]->fill(MET,weight_2);

          //2d histos w/ breit wigner reweighting

          for (double mw = 79.750; mw <= 79.900; mw += 0.001) {

            double bwreweight = bw(m_inv,mw,ww)/bw(m_inv,mw_ref,ww);

            if (bwreweight > 100) {
              MSG_INFO("bwreweight = " + std::to_string(bwreweight));
            }

            _h2["e pt"]->fill(mw,e_pTD,bwreweight);
            _h2["mt"]->fill(mw,mT,bwreweight);
            _h2["pt miss"]->fill(mw,MET,bwreweight);

            if (CDF) {
              // _h2w1["e pt pdf"]->fill(mw,e_pTD,bwreweight*weight);
              // _h2w1["mt pdf"]->fill(mw,mT,bwreweight*weight);
              // _h2w1["pt miss pdf"]->fill(mw,MET,bwreweight*weight);

              _h2w2["e pt pdf 2"]->fill(mw,e_pTD,bwreweight*weight_2);
              _h2w2["mt pdf 2"]->fill(mw,mT,bwreweight*weight_2);
              _h2w2["pt miss pdf 2"]->fill(mw,MET,bwreweight*weight_2);

              // _h2w3["e pt pdf 3"]->fill(mw,e_pTD,bwreweight*weight_3);
              // _h2w3["mt pdf 3"]->fill(mw,mT,bwreweight*weight_3);
              // _h2w3["pt miss pdf 3"]->fill(mw,MET,bwreweight*weight_3);
            }

            // mass fit in four rapidity bins for atlas
            if (CDF == false) {
              if (e_abseta < 0.8) {
                _h2["e pt 0-0.8"]->fill(mw,e_pTD,bwreweight);
                _h2["mt 0-0.8"]->fill(mw,mT,bwreweight);
              } else if (e_abseta < 1.4) {
                _h2["e pt 0.8-1.4"]->fill(mw,e_pTD,bwreweight);
                _h2["mt 0.8-1.4"]->fill(mw,mT,bwreweight);
              } else if (e_abseta < 2.0) {
                _h2["e pt 1.4-2.0"]->fill(mw,e_pTD,bwreweight);
                _h2["mt 1.4-2.0"]->fill(mw,mT,bwreweight);
              } else {
                // (i.e. <2.4 as max selection cut)
                _h2["e pt 2.0-2.4"]->fill(mw,e_pTD,bwreweight);
                _h2["mt 2.0-2.4"]->fill(mw,mT,bwreweight);
              }
            }

          }
          
        }

        //mirrors experimental distribution cuts in Errede and Strogolas
        if ( e_pTD > 20 && MET > 20 && e_abseta < 1 ) {

          _g["hA0_2"]->fill(P_T_W_truth,a0);
          _g["hA1_2"]->fill(P_T_W_truth,a1);
          _g["hA2_2"]->fill(P_T_W_truth,a2);
          _g["hA3_2"]->fill(P_T_W_truth,a3);
          _g["hA4_2"]->fill(P_T_W_truth,a4);
          // _h["hA5_2"]->fill(P_T_W_truth,a5);
          // _h["hA6_2"]->fill(P_T_W_truth,a6);
          // _h["hA7_2"]->fill(P_T_W_truth,a7);
          _h["hNormal_2"]->fill(P_T_W_truth);


          if (P_T_W_truth > 15) {
            if (P_T_W_truth < 25) {
              _h["phi15_25_2"]->fill(phi);
              _h["theta15_25_2"]->fill(costheta);
            } else if (P_T_W_truth < 35) {
              _h["phi25_35_2"]->fill(phi);
              _h["theta25_35_2"]->fill(costheta);
            } else if (P_T_W_truth < 65) {
              _h["phi35_65_2"]->fill(phi);
              _h["theta35_65_2"]->fill(costheta);
            } else if (P_T_W_truth < 105) {
              _h["phi65_105_2"]->fill(phi);
              _h["theta65_105_2"]->fill(costheta);
            }
            
          }
        }

        //fill event weight counters for all one electron one neutrino events
        // counter_1->fill(weight);
        counter_2->fill(weight_2);
        // counter_3->fill(weight_3);
        event_counter->fill(1);
      }

    }
   


    void finalize() {
      
      double luminosity = 0;

      //luminosity in inverse picobarns to match xsection dimensions
      if (CDF == true){ 
        luminosity = 8800;
      } else {
        luminosity = 4700; //atlas
      }


      // //to use crossSection() pass the cross section to rivet using rivet -x <cross section in pb>
      // double normVal = crossSection() * luminosity / numEvents();

      // for (auto& [name, histo] : _h) {
      //   scale(histo, histo->numEntries() * normVal);
      // }


      scale(_h, luminosity);
      // scale(_hw, luminosity * event_counter->val()/counter_1->val());
      scale(_hw2, luminosity * event_counter->val()/counter_2->val());
      // scale(_hw3, luminosity * event_counter->val()/counter_3->val());
      scale(_h2, luminosity);
      // scale(_h2w1, luminosity * event_counter->val()/counter_1->val());
      scale(_h2w2, luminosity * event_counter->val()/counter_2->val());
      // scale(_h2w3, luminosity * event_counter->val()/counter_3->val());

      

      MSG_INFO("e pt : " + std::to_string(_h["e pt"]->numEntries() ));
      MSG_INFO("e pt before : " + std::to_string(_h["e pt before"]->numEntries() ));
      MSG_INFO("e pt truth : " + std::to_string(_h["e pt truth"]->numEntries() ));


      // for (int i = 0; i < _h.size(); i++) {
      //   normalize(_h[i].second, _h[i].second->numEntries() * NormVal);
      // }



      // double norm = luminosity * crossSection() / ( sumOfWeights() ) / femtobarn;
      // double norm1 = norm  * quark_event_counter / counter_1->val();
      // double norm2 = norm  * quark_event_counter / counter_2->val();
      // double norm3 = norm  * quark_event_counter / counter_3->val();

      // scale(_h, norm);
      // scale(_g, norm);
      // scale(_h2, norm);
      // scale(_h2w1, norm1);
      // scale(_h2w2, norm2);
      // scale(_h2w3, norm3);
      // scale(_hw, norm1);
      // scale(_hw2, norm2);
      // scale(_hw3, norm3);

      if (CDF) {
        MSG_WARNING("analysis performed with CDF parameters");
      } else {
        //Atlas
        MSG_WARNING("analysis performed with Atlas parameters");
      }

    }

    ///@}

    /// @name Histograms
    ///@{

    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _g;
    // map<string, Histo1DPtr> _hw;
    map<string, Histo1DPtr> _hw2;
    // map<string, Histo1DPtr> _hw3;
    map<string, Histo2DPtr> _h2;
    // map<string, Histo2DPtr> _h2w1;
    map<string, Histo2DPtr> _h2w2;
    // map<string, Histo2DPtr> _h2w3;
    // CounterPtr counter_1;
    CounterPtr counter_2;
    // CounterPtr counter_3;
    CounterPtr event_counter;


    ///@}


  };

  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(MyAnalysis);

} 
