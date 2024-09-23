#include<iostream>
#include<fstream>
#include<vector>
#include"TH1.h"
#include"TH2.h"
#include"TF1.h"
#include"TMath.h"
#include"TCanvas.h"
#include"TFile.h"
#include"TLatex.h"
#include"TStyle.h"
#include"TPad.h"
#include"TLegend.h"
#include"TPaveText.h"
#include"TAxis.h"
#include"TTree.h"
#include"TFitResultPtr.h"
#include"TFitResult.h"
#include"TString.h"

using namespace std;

//const int energy = 510 is energy in GeV - for fig. labels
void Closure_test_new_US_LS(const int energy = 510)
{
  const float L0_alpha = 0.732; //decay parameter of L0
  const float L0_alpha_relat_err = 0.014/L0_alpha; //relative error of decay parameter

  const float L0bar_alpha = -0.758; //decay paramteter of L0bar
  const float L0bar_alpha_relat_err = 0.012/fabs(L0bar_alpha); //relative error of decay paramteter


  const int nPtBins = 8;
  float const pT_bins[nPtBins+1] = { 0., 0.5, 1.,1.5, 2., 2.5, 3., 4., 5.};

  const int nPtBins_corr = 2;
  float const pT_bins_corr[nPtBins_corr+1] = { 0.5, 1.5, 5.};

  const int nEtaBins = 3;
  //float const eta_bins[nEtaBins+1] = { -1, -0.2, 0.2, 1 };
  float const eta_bins[nEtaBins+1] = { -1, -0.4, 0.4, 1 };

  //TFile *outFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/output/L_cosThetaStar_closure_test_work.root", "recreate");

  //load all files
  TFile *inFile;


  if(energy == 510)
  {
    inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run17/output_Lambda_pp_510_MB_1B_events_new_hists_Minv.root", "READ");
  }
  else if(energy == 200)
  {
    inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/old_2/output_Lambda_pp_200_MB_1B_events_hists_work.root", "READ");

    //ME tests
    //inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/ME_tests/output_Lambda_pp_200_MB_1B_events_hists_standard_ME.root", "READ");
    //inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/ME_tests/output_Lambda_pp_200_MB_1B_events_hists_single_L_ME.root", "READ");
    //inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/ME_tests/output_Lambda_pp_200_MB_1B_events_hists_delta_phi_cut_ME.root", "READ");
    //inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/ME_tests/output_Lambda_pp_200_MB_1B_events_hists_smear_phi.root", "READ");
    //inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/ME_tests/output_Lambda_pp_200_MB_1B_events_hists_smear_phi_2.root", "READ");
    //inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/ME_tests/output_Lambda_pp_200_MB_1B_events_hists_new_ME_with_SE_kine.root", "READ");

    //no ME weight test
    //inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/old_ME_tests/output_Lambda_pp_200_MB_1B_events_hists_old_ME_work.root", "READ");

    //inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/mom_smear_new/output_Lambda_pp_200_MB_1B_events_hists_pair_weight_mom_smear_full.root", "READ");
    //inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/mom_smear_pT/output_Lambda_pp_200_MB_1B_events_hists_pT_times_20.root", "READ");

  }
  else
  {
    cout<<"Not a valid collision energy! Abborting!"<<endl;
    return;
  }


  TFile *sysErrFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/output/sys_err/SysErrSlope.root", "recreate");

  //histogram to store systematic error from slope difference in closure test
  //each bin is for one L charge combination
  TH1F *SysErrSlope_hist = new TH1F("SysErrSlope_hist", "SysErrSlope_hist", 3, 0 , 3);

  //histograms

  //p pi pairs - to simulate combinatorial background
  //invariant mass before cuts
  TH2F *L0_inv_mass_vs_L0bar_inv_mass_US[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0_inv_mass_vs_L0bar_inv_mass_US_LS[nPtBins_corr][nPtBins_corr]; //for US-LS Lambda pairs

  TH2F *L0_inv_mass_vs_L0_inv_mass_US[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0_inv_mass_vs_L0_inv_mass_US_LS[nPtBins_corr][nPtBins_corr]; //for US-LS Lambda pairs

  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_US[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[nPtBins_corr][nPtBins_corr]; //for US-LS Lambda pairs

  //invariant mass after cuts
  TH2F *L0_inv_mass_vs_L0bar_inv_mass_US_cuts[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0_inv_mass_vs_L0bar_inv_mass_US_LS_cuts[nPtBins_corr][nPtBins_corr]; //for US-LS Lambda pairs

  TH2F *L0_inv_mass_vs_L0_inv_mass_US_cuts[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0_inv_mass_vs_L0_inv_mass_US_LS_cuts[nPtBins_corr][nPtBins_corr]; //for US-LS Lambda pairs

  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_US_cuts[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_cuts[nPtBins_corr][nPtBins_corr]; //for US-LS Lambda pairs

  //before cuts
  //True MC
  TH1F *L0_L0bar_cosThetaProdPlane = (TH1F*)inFile->Get("L0_L0bar_cosThetaProdPlane");
  TH1F *L0_L0bar_cosThetaProdPlane_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0bar_cosThetaProdPlane_eta_hist[nEtaBins][nEtaBins];


  TH1F *L0_L0_cosThetaProdPlane = (TH1F*)inFile->Get("L0_L0_cosThetaProdPlane");
  TH1F *L0_L0_cosThetaProdPlane_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0_cosThetaProdPlane_eta_hist[nEtaBins][nEtaBins];


  TH1F *L0bar_L0bar_cosThetaProdPlane = (TH1F*)inFile->Get("L0bar_L0bar_cosThetaProdPlane");
  TH1F *L0bar_L0bar_cosThetaProdPlane_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar_L0bar_cosThetaProdPlane_eta_hist[nEtaBins][nEtaBins];

  //pi kinematics
  TH2F *L0_L0bar_pi_pT1_vs_pi_pT2_hist = (TH2F*)inFile->Get("L0_L0bar_pi_pT1_vs_pi_pT2_hist");

  //----------------------------------------------------------------------------------------------

  //unlike-sign (signal + background)
  TH1F *L0_L0bar_cosThetaProdPlane_US = (TH1F*)inFile->Get("L0_L0bar_cosThetaProdPlane_US");
  TH1F *L0_L0bar_cosThetaProdPlane_US_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0bar_cosThetaProdPlane_US_eta_hist[nEtaBins][nEtaBins];


  TH1F *L0_L0_cosThetaProdPlane_US = (TH1F*)inFile->Get("L0_L0_cosThetaProdPlane_US");
  TH1F *L0_L0_cosThetaProdPlane_US_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0_cosThetaProdPlane_US_eta_hist[nEtaBins][nEtaBins];


  TH1F *L0bar_L0bar_cosThetaProdPlane_US = (TH1F*)inFile->Get("L0bar_L0bar_cosThetaProdPlane_US");
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_eta_hist[nEtaBins][nEtaBins];


  //US matched to LS (background)
  TH1F *L0_L0bar_cosThetaProdPlane_US_LS = (TH1F*)inFile->Get("L0_L0bar_cosThetaProdPlane_US_LS");
  TH1F *L0_L0bar_cosThetaProdPlane_US_LS_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0bar_cosThetaProdPlane_US_LS_eta_hist[nEtaBins][nEtaBins];

  TH1F *L0_L0_cosThetaProdPlane_US_LS = (TH1F*)inFile->Get("L0_L0_cosThetaProdPlane_US_LS");
  TH1F *L0_L0_cosThetaProdPlane_US_LS_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0_cosThetaProdPlane_US_LS_eta_hist[nEtaBins][nEtaBins];


  TH1F *L0bar_L0bar_cosThetaProdPlane_US_LS = (TH1F*)inFile->Get("L0bar_L0bar_cosThetaProdPlane_US_LS");
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_LS_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_LS_eta_hist[nEtaBins][nEtaBins];

  //______________________________________________________________________________________________________

  //mixed event
  //True MC
  TH1F *L0_L0bar_cosThetaProdPlane_ME_weight = (TH1F*)inFile->Get("L0_L0bar_cosThetaProdPlane_ME_weight");
  TH1F *L0_L0bar_cosThetaProdPlane_ME_weight_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0bar_cosThetaProdPlane_ME_weight_eta_hist[nEtaBins][nEtaBins];


  TH1F *L0_L0_cosThetaProdPlane_ME_weight = (TH1F*)inFile->Get("L0_L0_cosThetaProdPlane_ME_weight");
  TH1F *L0_L0_cosThetaProdPlane_ME_weight_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0_cosThetaProdPlane_ME_weight_eta_hist[nEtaBins][nEtaBins];


  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_weight = (TH1F*)inFile->Get("L0bar_L0bar_cosThetaProdPlane_ME_weight");
  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_weight_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_weight_eta_hist[nEtaBins][nEtaBins];

  //pi kinematics
  TH2F *L0_L0bar_pi_pT1_vs_pi_pT2_ME_hist = (TH2F*)inFile->Get("L0_L0bar_pi_pT1_vs_pi_pT2_ME_hist");

  //---------------------------------------------------------------------------------------------------------

  //US
  TH1F *L0_L0bar_cosThetaProdPlane_US_ME_weight = (TH1F*)inFile->Get("L0_L0bar_cosThetaProdPlane_US_ME_weight");
  TH1F *L0_L0bar_cosThetaProdPlane_US_ME_weight_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0bar_cosThetaProdPlane_US_ME_weight_eta_hist[nEtaBins][nEtaBins];


  TH1F *L0_L0_cosThetaProdPlane_US_ME_weight = (TH1F*)inFile->Get("L0_L0_cosThetaProdPlane_US_ME_weight");
  TH1F *L0_L0_cosThetaProdPlane_US_ME_weight_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0_cosThetaProdPlane_US_ME_weight_eta_hist[nEtaBins][nEtaBins];


  TH1F *L0bar_L0bar_cosThetaProdPlane_US_ME_weight = (TH1F*)inFile->Get("L0bar_L0bar_cosThetaProdPlane_US_ME_weight");
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_ME_weight_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_ME_weight_eta_hist[nEtaBins][nEtaBins];


  //US matched to LS (background)
  TH1F *L0_L0bar_cosThetaProdPlane_US_LS_ME_weight = (TH1F*)inFile->Get("L0_L0bar_cosThetaProdPlane_US_LS_ME_weight");
  TH1F *L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_hist[nEtaBins][nEtaBins];


  TH1F *L0_L0_cosThetaProdPlane_US_LS_ME_weight = (TH1F*)inFile->Get("L0_L0_cosThetaProdPlane_US_LS_ME_weight");
  TH1F *L0_L0_cosThetaProdPlane_US_LS_ME_weight_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0_cosThetaProdPlane_US_LS_ME_weight_eta_hist[nEtaBins][nEtaBins];


  TH1F *L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight = (TH1F*)inFile->Get("L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight");
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_hist[nEtaBins][nEtaBins];

  //_____________________________________________________________________________________________________________________________

  //after cuts

  //True MC
  TH1F *L0_L0bar_cosThetaProdPlane_cuts = (TH1F*)inFile->Get("L0_L0bar_cosThetaProdPlane_cuts");
  TH1F *L0_L0bar_cosThetaProdPlane_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0bar_cosThetaProdPlane_eta_cuts_hist[nEtaBins][nEtaBins];


  TH1F *L0_L0_cosThetaProdPlane_cuts = (TH1F*)inFile->Get("L0_L0_cosThetaProdPlane_cuts");
  TH1F *L0_L0_cosThetaProdPlane_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0_cosThetaProdPlane_eta_cuts_hist[nEtaBins][nEtaBins];


  TH1F *L0bar_L0bar_cosThetaProdPlane_cuts = (TH1F*)inFile->Get("L0bar_L0bar_cosThetaProdPlane_cuts");
  TH1F *L0bar_L0bar_cosThetaProdPlane_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar_L0bar_cosThetaProdPlane_eta_cuts_hist[nEtaBins][nEtaBins];

  //pi kinematics
  TH2F *L0_L0bar_pi_pT1_vs_pi_pT2_cuts_hist = (TH2F*)inFile->Get("L0_L0bar_pi_pT1_vs_pi_pT2_cuts_hist");

  //---------------------------------------------------------------------------------------------------------

  //unlike-sign (signal + background)
  TH1F *L0_L0bar_cosThetaProdPlane_US_cuts = (TH1F*)inFile->Get("L0_L0bar_cosThetaProdPlane_US_cuts");
  TH1F *L0_L0bar_cosThetaProdPlane_US_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0bar_cosThetaProdPlane_US_eta_cuts_hist[nEtaBins][nEtaBins];


  TH1F *L0_L0_cosThetaProdPlane_US_cuts = (TH1F*)inFile->Get("L0_L0_cosThetaProdPlane_US_cuts");
  TH1F *L0_L0_cosThetaProdPlane_US_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0_cosThetaProdPlane_US_eta_cuts_hist[nEtaBins][nEtaBins];


  TH1F *L0bar_L0bar_cosThetaProdPlane_US_cuts = (TH1F*)inFile->Get("L0bar_L0bar_cosThetaProdPlane_US_cuts");
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_eta_cuts_hist[nEtaBins][nEtaBins];


  //US matched to LS (background)
  TH1F *L0_L0bar_cosThetaProdPlane_US_LS_cuts = (TH1F*)inFile->Get("L0_L0bar_cosThetaProdPlane_US_LS_cuts");
  TH1F *L0_L0bar_cosThetaProdPlane_US_LS_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0bar_cosThetaProdPlane_US_LS_eta_cuts_hist[nEtaBins][nEtaBins];


  TH1F *L0_L0_cosThetaProdPlane_US_LS_cuts = (TH1F*)inFile->Get("L0_L0_cosThetaProdPlane_US_LS_cuts");
  TH1F *L0_L0_cosThetaProdPlane_US_LS_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0_cosThetaProdPlane_US_LS_eta_cuts_hist[nEtaBins][nEtaBins];


  TH1F *L0bar_L0bar_cosThetaProdPlane_US_LS_cuts = (TH1F*)inFile->Get("L0bar_L0bar_cosThetaProdPlane_US_LS_cuts");
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_LS_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_LS_eta_cuts_hist[nEtaBins][nEtaBins];

  //_____________________________________________________________________________________________________________________________

  //mixed event
  //true MC
  //weight from distributions after cuts
  TH1F *L0_L0bar_cosThetaProdPlane_ME_weight_cuts = (TH1F*)inFile->Get("L0_L0bar_cosThetaProdPlane_ME_weight_cuts");
  TH1F *L0_L0bar_cosThetaProdPlane_ME_weight_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0bar_cosThetaProdPlane_ME_weight_eta_cuts_hist[nEtaBins][nEtaBins];


  TH1F *L0_L0_cosThetaProdPlane_ME_weight_cuts = (TH1F*)inFile->Get("L0_L0_cosThetaProdPlane_ME_weight_cuts");
  TH1F *L0_L0_cosThetaProdPlane_ME_weight_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0_cosThetaProdPlane_ME_weight_eta_cuts_hist[nEtaBins][nEtaBins];


  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_weight_cuts = (TH1F*)inFile->Get("L0bar_L0bar_cosThetaProdPlane_ME_weight_cuts");
  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_weight_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_weight_eta_cuts_hist[nEtaBins][nEtaBins];

  //pi kinematics
  TH2F *L0_L0bar_pi_pT1_vs_pi_pT2_ME_cuts_hist = (TH2F*)inFile->Get("L0_L0bar_pi_pT1_vs_pi_pT2_ME_cuts_hist");


  //weight from distributions befire cuts
  TH1F *L0_L0bar_cosThetaProdPlane_ME_weight_2_cuts = (TH1F*)inFile->Get("L0_L0bar_cosThetaProdPlane_ME_weight_2_cuts");
  TH1F *L0_L0bar_cosThetaProdPlane_ME_weight_2_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0bar_cosThetaProdPlane_ME_weight_2_eta_cuts_hist[nEtaBins][nEtaBins];


  TH1F *L0_L0_cosThetaProdPlane_ME_weight_2_cuts = (TH1F*)inFile->Get("L0_L0_cosThetaProdPlane_ME_weight_2_cuts");
  TH1F *L0_L0_cosThetaProdPlane_ME_weight_2_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0_cosThetaProdPlane_ME_weight_2_eta_cuts_hist[nEtaBins][nEtaBins];


  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_weight_2_cuts = (TH1F*)inFile->Get("L0bar_L0bar_cosThetaProdPlane_ME_weight_2_cuts");
  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_weight_2_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_weight_2_eta_cuts_hist[nEtaBins][nEtaBins];

  //pi kinematics
  TH2F *L0_L0bar_pi_pT1_vs_pi_pT2_ME_cuts_2_hist = (TH2F*)inFile->Get("L0_L0bar_pi_pT1_vs_pi_pT2_ME_cuts_2_hist");


  //--------------------------------------------------------------------------------------------------------------------------

  //US
  TH1F *L0_L0bar_cosThetaProdPlane_US_ME_weight_cuts = (TH1F*)inFile->Get("L0_L0bar_cosThetaProdPlane_US_ME_weight_cuts");
  TH1F *L0_L0bar_cosThetaProdPlane_US_ME_weight_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0bar_cosThetaProdPlane_US_ME_weight_eta_cuts_hist[nEtaBins][nEtaBins];


  TH1F *L0_L0_cosThetaProdPlane_US_ME_weight_cuts = (TH1F*)inFile->Get("L0_L0_cosThetaProdPlane_US_ME_weight_cuts");
  TH1F *L0_L0_cosThetaProdPlane_US_ME_weight_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0_cosThetaProdPlane_US_ME_weight_eta_cuts_hist[nEtaBins][nEtaBins];


  TH1F *L0bar_L0bar_cosThetaProdPlane_US_ME_weight_cuts = (TH1F*)inFile->Get("L0bar_L0bar_cosThetaProdPlane_US_ME_weight_cuts");
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_ME_weight_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_ME_weight_eta_cuts_hist[nEtaBins][nEtaBins];


  //US matched to LS (background)
  TH1F *L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts = (TH1F*)inFile->Get("L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts");
  TH1F *L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_cuts_hist[nEtaBins][nEtaBins];


  TH1F *L0_L0_cosThetaProdPlane_US_LS_ME_weight_cuts = (TH1F*)inFile->Get("L0_L0_cosThetaProdPlane_US_LS_ME_weight_cuts");
  TH1F *L0_L0_cosThetaProdPlane_US_LS_ME_weight_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0_cosThetaProdPlane_US_LS_ME_weight_eta_cuts_hist[nEtaBins][nEtaBins];


  TH1F *L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts = (TH1F*)inFile->Get("L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts");
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_cuts_hist[nEtaBins][nEtaBins];
  //_________________________________________________________________________________________________________________________________________________________________________________

  //True MC before cuts
  TCanvas *L0_L0bar_cosThetaProdPlane_can = new TCanvas("L0_L0bar_cosThetaProdPlane_can", "L0_L0bar_cosThetaProdPlane_can", 1200, 1000);
  TCanvas *L0_L0_cosThetaProdPlane_can = new TCanvas("L0_L0_cosThetaProdPlane_can", "L0_L0_cosThetaProdPlane_can", 1200, 1000);
  TCanvas *L0bar_L0bar_cosThetaProdPlane_can = new TCanvas("L0bar_L0bar_cosThetaProdPlane_can", "L0bar_L0bar_cosThetaProdPlane_can", 1200, 1000);

  //True MC after cuts
  TCanvas *L0_L0bar_cosThetaProdPlane_cuts_can = new TCanvas("L0_L0bar_cosThetaProdPlane_cuts_can", "L0_L0bar_cosThetaProdPlane_cuts_can", 1200, 1000);
  TCanvas *L0_L0bar_cosThetaProdPlane_cuts_2_can = new TCanvas("L0_L0bar_cosThetaProdPlane_cuts_2_can", "L0_L0bar_cosThetaProdPlane_cuts_2_can", 1200, 1000);

  TCanvas *L0_L0_cosThetaProdPlane_cuts_can = new TCanvas("L0_L0_cosThetaProdPlane_cuts_can", "L0_L0_cosThetaProdPlane_cuts_can", 1200, 1000);
  TCanvas *L0_L0_cosThetaProdPlane_cuts_2_can = new TCanvas("L0_L0_cosThetaProdPlane_cuts_2_can", "L0_L0_cosThetaProdPlane_cuts_2_can", 1200, 1000);

  TCanvas *L0bar_L0bar_cosThetaProdPlane_cuts_can = new TCanvas("L0bar_L0bar_cosThetaProdPlane_cuts_can", "L0bar_L0bar_cosThetaProdPlane_cuts_can", 1200, 1000);
  TCanvas *L0bar_L0bar_cosThetaProdPlane_cuts_2_can = new TCanvas("L0bar_L0bar_cosThetaProdPlane_cuts_2_can", "L0bar_L0bar_cosThetaProdPlane_cuts_2_can", 1200, 1000);


  //US pairs before cuts
  TCanvas *L0_L0bar_cosThetaProdPlane_US_no_corr_can = new TCanvas("L0_L0bar_cosThetaProdPlane_US_no_corr_can", "L0_L0bar_cosThetaProdPlane_US_no_corr_can", 1200, 1000);
  TCanvas *L0_L0bar_cosThetaProdPlane_US_separate_can = new TCanvas("L0_L0bar_cosThetaProdPlane_US_separate_can", "L0_L0bar_cosThetaProdPlane_US_separate_can", 1200, 1000);
  TCanvas *L0_L0bar_cosThetaProdPlane_US_can = new TCanvas("L0_L0bar_cosThetaProdPlane_US_can", "L0_L0bar_cosThetaProdPlane_US_can", 1200, 1000);

  TCanvas *L0_L0_cosThetaProdPlane_US_can = new TCanvas("L0_L0_cosThetaProdPlane_US_can", "L0_L0_cosThetaProdPlane_US_can", 1200, 1000);
  TCanvas *L0bar_L0bar_cosThetaProdPlane_US_can = new TCanvas("L0bar_L0bar_cosThetaProdPlane_US_can", "L0bar_L0bar_cosThetaProdPlane_US_can", 1200, 1000);

  //US pairs after cuts
  TCanvas *L0_L0bar_cosThetaProdPlane_US_no_corr_cuts_can = new TCanvas("L0_L0bar_cosThetaProdPlane_US_no_corr_cuts_can", "L0_L0bar_cosThetaProdPlane_US_no_corr_cuts_can", 1200, 1000);
  TCanvas *L0_L0bar_cosThetaProdPlane_US_separate_cuts_can = new TCanvas("L0_L0bar_cosThetaProdPlane_US_separate_cuts_can", "L0_L0bar_cosThetaProdPlane_US_separate_cuts_can", 1200, 1000);
  TCanvas *L0_L0bar_cosThetaProdPlane_US_cuts_can = new TCanvas("L0_L0bar_cosThetaProdPlane_US_cuts_can", "L0_L0bar_cosThetaProdPlane_US_cuts_can", 1200, 1000);

  TCanvas *L0_L0_cosThetaProdPlane_US_cuts_can = new TCanvas("L0_L0_cosThetaProdPlane_US_cuts_can", "L0_L0_cosThetaProdPlane_US_cuts_can", 1200, 1000);
  TCanvas *L0bar_L0bar_cosThetaProdPlane_US_cuts_can = new TCanvas("L0bar_L0bar_cosThetaProdPlane_US_cuts_can", "L0bar_L0bar_cosThetaProdPlane_US_cuts_can", 1200, 1000);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  //_______________________________________________________________________________________________________________________________________________________________________________

  //outFile->cd();

  //closure test

  //true MC corrected with ME
  //L-Lbar
  //before cuts
  TF1 *fit_L0_L0bar_before_cuts = new TF1("fit_L0_L0bar_before_cuts", "[0]*(1+[1]*x)", -1, 1);
  fit_L0_L0bar_before_cuts->SetParameters(1000, 0.10);

  L0_L0bar_cosThetaProdPlane_can->cd();

  float nL0L0bar = L0_L0bar_cosThetaProdPlane->Integral();

  L0_L0bar_cosThetaProdPlane_ME_weight->Sumw2();
  //L0_L0bar_cosThetaProdPlane_ME_weight->Add(L0_L0bar_cosThetaProdPlane_LS_ME_weight);
  L0_L0bar_cosThetaProdPlane_ME_weight->Scale(nL0L0bar/L0_L0bar_cosThetaProdPlane_ME_weight->Integral());
  //L0_L0bar_cosThetaProdPlane_ME_weight->Draw("p e");

  L0_L0bar_cosThetaProdPlane->Sumw2();
  L0_L0bar_cosThetaProdPlane->SetMarkerStyle(20);
  L0_L0bar_cosThetaProdPlane->SetMarkerColor(kRed);
  L0_L0bar_cosThetaProdPlane->SetLineColor(kRed);
  L0_L0bar_cosThetaProdPlane->GetXaxis()->SetTitle("cos(#theta*)");
  L0_L0bar_cosThetaProdPlane->GetXaxis()->CenterTitle();
  L0_L0bar_cosThetaProdPlane->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0_L0bar_cosThetaProdPlane->GetYaxis()->CenterTitle();
  L0_L0bar_cosThetaProdPlane->Divide(L0_L0bar_cosThetaProdPlane_ME_weight); //correct using ME
  L0_L0bar_cosThetaProdPlane->Scale(nL0L0bar/L0_L0bar_cosThetaProdPlane->Integral()); //scale back
  L0_L0bar_cosThetaProdPlane->Scale(1./L0_L0bar_cosThetaProdPlane->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
  L0_L0bar_cosThetaProdPlane->Fit(fit_L0_L0bar_before_cuts, "i 0 r");
  L0_L0bar_cosThetaProdPlane->SetMinimum(0);
  //L0_L0bar_cosThetaProdPlane->GetYaxis()->SetRangeUser(-20e9, 0);
  L0_L0bar_cosThetaProdPlane->Draw("p e");

  //L0_L0bar_cosThetaProdPlane_ME_weight->Draw("p e same");

  fit_L0_L0bar_before_cuts->SetLineColor(1);
  fit_L0_L0bar_before_cuts->Draw("same");


  float L0_L0bar_slope = fit_L0_L0bar_before_cuts->GetParameter(1)/(L0_alpha*L0bar_alpha);
  float L0_L0bar_slope_err = fit_L0_L0bar_before_cuts->GetParError(1)/(L0_alpha*L0bar_alpha);


  TPaveText *L0_L0bar_text_MC = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
  L0_L0bar_text_MC->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  L0_L0bar_text_MC->AddText("Minimum bias");
  L0_L0bar_text_MC->AddText("#Lambda-#bar{#Lambda}");
  L0_L0bar_text_MC->AddText("True MC");
  L0_L0bar_text_MC->AddText("No cuts");
  L0_L0bar_text_MC->AddText(Form("Fit slope: %.3f #pm %.3f", L0_L0bar_slope, L0_L0bar_slope_err));
  L0_L0bar_text_MC->SetFillColorAlpha(0, 0.01);
  L0_L0bar_text_MC->Draw("same");

  L0_L0bar_cosThetaProdPlane_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/L0_L0bar_cosThetaProdPlane.png");


  //pi kinematics, before cuts

  //int firstBin = L0_L0bar_pi_pT1_vs_pi_pT2_hist->GetXaxis()->FindBin(0.151);
  //int lastBin = L0_L0bar_pi_pT1_vs_pi_pT2_hist->GetXaxis()->FindBin(0.99);

  TCanvas *pi_plus_SE_vs_ME_MC = new TCanvas("pi_plus_SE_vs_ME_MC", "pi_plus_SE_vs_ME_MC", 1200, 1000);

  pi_plus_SE_vs_ME_MC->cd();

  TH1D *pi_plus_SE_MC_hist = L0_L0bar_pi_pT1_vs_pi_pT2_hist->ProjectionY("pi_plus_SE_MC_hist"); //projection to Lbar - i.e. pi+
  pi_plus_SE_MC_hist->SetLineColor(1);
  //pi_plus_SE_MC_hist->Scale(1./pi_plus_SE_MC_hist->Integral(firstBin, lastBin));
  pi_plus_SE_MC_hist->Scale(1./pi_plus_SE_MC_hist->Integral());
  pi_plus_SE_MC_hist->GetXaxis()->SetTitle("p_{T} (GeV/#it{c})");
  pi_plus_SE_MC_hist->GetXaxis()->CenterTitle();
  pi_plus_SE_MC_hist->GetXaxis()->SetRangeUser(0,1);
  pi_plus_SE_MC_hist->GetYaxis()->SetRangeUser(0,0.11);
  pi_plus_SE_MC_hist->Draw("hist");

  TH1D *pi_plus_ME_MC_hist = L0_L0bar_pi_pT1_vs_pi_pT2_ME_hist->ProjectionY("pi_plus_ME_MC_hist");
  pi_plus_ME_MC_hist->SetLineColor(kRed);
  pi_plus_ME_MC_hist->Scale(1./pi_plus_ME_MC_hist->Integral());
  pi_plus_ME_MC_hist->Draw("hist same");

  TLegend *SE_vs_ME_leg = new TLegend(0.65, 0.5, 0.89, 0.69);
  SE_vs_ME_leg->AddEntry(pi_plus_SE_MC_hist, "Same event");
  SE_vs_ME_leg->AddEntry(pi_plus_ME_MC_hist, "Mixed event");
  SE_vs_ME_leg->SetBorderSize(0);
  SE_vs_ME_leg->SetFillColorAlpha(0, 0.01);
  SE_vs_ME_leg->Draw("same");


  TPaveText *text_pi_plus_MC = new TPaveText(0.65, 0.69, 0.89, 0.89, "NDC");
  text_pi_plus_MC->SetTextFont(42);
  text_pi_plus_MC->AddText("#pi^{+}");
  text_pi_plus_MC->AddText("PYTHIA before cuts");
  text_pi_plus_MC->AddText(Form("SE mean = %.3f", pi_plus_SE_MC_hist->GetMean()));
  text_pi_plus_MC->AddText(Form("ME mean = %.3f", pi_plus_ME_MC_hist->GetMean()));
  text_pi_plus_MC->SetFillColorAlpha(0, 0.01);
  text_pi_plus_MC->Draw("same");

  pi_plus_SE_vs_ME_MC->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/pi_pT/pi_plus_SE_vs_ME_MC.png");


  TCanvas *pi_minus_SE_vs_ME_MC = new TCanvas("pi_minus_SE_vs_ME_MC", "pi_minus_SE_vs_ME_MC", 1200, 1000);

  pi_minus_SE_vs_ME_MC->cd();

  TH1D *pi_minus_SE_MC_hist = L0_L0bar_pi_pT1_vs_pi_pT2_hist->ProjectionX("pi_minus_SE_MC_hist"); //projection to Lbar - i.e. pi+
  pi_minus_SE_MC_hist->SetLineColor(1);
  //pi_minus_SE_MC_hist->Scale(1./pi_minus_SE_MC_hist->Integral(firstBin, lastBin));
  pi_minus_SE_MC_hist->Scale(1./pi_minus_SE_MC_hist->Integral());
  pi_minus_SE_MC_hist->GetXaxis()->SetTitle("p_{T} (GeV/#it{c})");
  pi_minus_SE_MC_hist->GetXaxis()->CenterTitle();
  pi_minus_SE_MC_hist->GetXaxis()->SetRangeUser(0,1);
  pi_minus_SE_MC_hist->GetYaxis()->SetRangeUser(0,0.11);
  pi_minus_SE_MC_hist->Draw("hist");

  TH1D *pi_minus_ME_MC_hist = L0_L0bar_pi_pT1_vs_pi_pT2_ME_hist->ProjectionX("pi_minus_ME_MC_hist");
  pi_minus_ME_MC_hist->SetLineColor(kRed);
  pi_minus_ME_MC_hist->Scale(1./pi_minus_ME_MC_hist->Integral());
  pi_minus_ME_MC_hist->Draw("hist same");


  SE_vs_ME_leg->Draw("same");


  TPaveText *text_pi_minus_MC = new TPaveText(0.65, 0.69, 0.89, 0.89, "NDC");
  text_pi_minus_MC->SetTextFont(42);
  text_pi_minus_MC->AddText("#pi^{-}");
  text_pi_minus_MC->AddText("PYTHIA before cuts");
  text_pi_minus_MC->AddText(Form("SE mean = %.3f", pi_minus_SE_MC_hist->GetMean()));
  text_pi_minus_MC->AddText(Form("ME mean = %.3f", pi_minus_ME_MC_hist->GetMean()));
  text_pi_minus_MC->SetFillColorAlpha(0, 0.01);
  text_pi_minus_MC->Draw("same");

  pi_minus_SE_vs_ME_MC->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/pi_pT/pi_minus_SE_vs_ME_MC.png");


  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  //after cuts
  TF1 *fit_L0_L0bar_after_cuts = new TF1("fit_L0_L0bar_after_cuts", "[0]*(1+[1]*x)", -1, 1);
  fit_L0_L0bar_after_cuts->SetParameters(1000, 0.10);

  TF1 *fit_L0_L0bar_after_cuts_ME = new TF1("fit_L0_L0bar_after_cuts_ME", "[0]*(1+[1]*x)", -1, 1);
  fit_L0_L0bar_after_cuts_ME->SetParameters(1000, 0.10);

  L0_L0bar_cosThetaProdPlane_cuts_can->cd();

  float nL0L0bar_cuts = L0_L0bar_cosThetaProdPlane_cuts->Integral();

  L0_L0bar_cosThetaProdPlane_ME_weight_cuts->Sumw2();
  L0_L0bar_cosThetaProdPlane_ME_weight_cuts->Scale(nL0L0bar_cuts/L0_L0bar_cosThetaProdPlane_ME_weight_cuts->Integral());
  L0_L0bar_cosThetaProdPlane_ME_weight_cuts->Scale(1./L0_L0bar_cosThetaProdPlane_ME_weight_cuts->GetXaxis()->GetBinWidth(1));
  L0_L0bar_cosThetaProdPlane_ME_weight_cuts->Fit(fit_L0_L0bar_after_cuts_ME, "i 0 r");
  //L0_L0bar_cosThetaProdPlane_ME_weight->Draw("p e");

  L0_L0bar_cosThetaProdPlane_cuts->Sumw2();
  L0_L0bar_cosThetaProdPlane_cuts->SetMarkerStyle(20);
  L0_L0bar_cosThetaProdPlane_cuts->SetMarkerColor(kRed);
  L0_L0bar_cosThetaProdPlane_cuts->SetLineColor(kRed);
  L0_L0bar_cosThetaProdPlane_cuts->GetXaxis()->SetTitle("cos(#theta*)");
  L0_L0bar_cosThetaProdPlane_cuts->GetXaxis()->CenterTitle();
  L0_L0bar_cosThetaProdPlane_cuts->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0_L0bar_cosThetaProdPlane_cuts->GetYaxis()->CenterTitle();
  //L0_L0bar_cosThetaProdPlane_cuts->Divide(L0_L0bar_cosThetaProdPlane_ME_weight_cuts); //correct using ME
  //L0_L0bar_cosThetaProdPlane_cuts->Scale(nL0L0bar_cuts/L0_L0bar_cosThetaProdPlane_cuts->Integral()); //scale back
  L0_L0bar_cosThetaProdPlane_cuts->Scale(1./L0_L0bar_cosThetaProdPlane_cuts->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
  L0_L0bar_cosThetaProdPlane_cuts->Fit(fit_L0_L0bar_after_cuts, "i 0 r");
  L0_L0bar_cosThetaProdPlane_cuts->SetMinimum(0);
  //L0_L0bar_cosThetaProdPlane_cuts->GetYaxis()->SetRangeUser(-20e9, 0);
  L0_L0bar_cosThetaProdPlane_cuts->Draw("p e");

  L0_L0bar_cosThetaProdPlane_ME_weight_cuts->Draw("p e same");

  fit_L0_L0bar_after_cuts->SetLineColor(1);
  fit_L0_L0bar_after_cuts->Draw("same");

  fit_L0_L0bar_after_cuts_ME->SetLineColor(kBlue);
  fit_L0_L0bar_after_cuts_ME->Draw("same");


  float L0_L0bar_slope_cuts = fit_L0_L0bar_after_cuts->GetParameter(1)/(L0_alpha*L0bar_alpha);
  float L0_L0bar_slope_cuts_err = fit_L0_L0bar_after_cuts->GetParError(1)/(L0_alpha*L0bar_alpha);


  float L0_L0bar_slope_cuts_sys_err = fabs( fit_L0_L0bar_after_cuts->GetParameter(1) - fit_L0_L0bar_after_cuts_ME->GetParameter(1) )/fit_L0_L0bar_after_cuts->GetParameter(1);

  SysErrSlope_hist->SetBinContent(1, L0_L0bar_slope_cuts_sys_err);

  TPaveText *L0_L0bar_text_MC_cuts = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
  L0_L0bar_text_MC_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  L0_L0bar_text_MC_cuts->AddText("Minimum bias");
  L0_L0bar_text_MC_cuts->AddText("#Lambda-#bar{#Lambda}");
  L0_L0bar_text_MC_cuts->AddText("True MC");
  L0_L0bar_text_MC_cuts->AddText("Analysis cuts");
  L0_L0bar_text_MC_cuts->AddText(Form("P = %.3f #pm %.3f", L0_L0bar_slope_cuts, L0_L0bar_slope_cuts_err ));
  L0_L0bar_text_MC_cuts->AddText(Form("P_{ME} = %.3f", fit_L0_L0bar_after_cuts_ME->GetParameter(1)/(L0_alpha*L0bar_alpha) ));
  L0_L0bar_text_MC_cuts->AddText(Form("Slope sys. err: %.1f %%", L0_L0bar_slope_cuts_sys_err*100 ));
  L0_L0bar_text_MC_cuts->SetFillColorAlpha(0, 0.01);
  L0_L0bar_text_MC_cuts->Draw("same");

  L0_L0bar_cosThetaProdPlane_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/L0_L0bar_cosThetaProdPlane_cuts.png");

  //pi kinematics
  TCanvas *pi_plus_SE_vs_ME_MC_cuts = new TCanvas("pi_plus_SE_vs_ME_MC_cuts", "pi_plus_SE_vs_ME_MC_cuts", 1200, 1000);

  pi_plus_SE_vs_ME_MC_cuts->cd();

  TH1D *pi_plus_SE_MC_cuts_hist = L0_L0bar_pi_pT1_vs_pi_pT2_cuts_hist->ProjectionY("pi_plus_SE_MC_cuts_hist"); //projection to Lbar - i.e. pi+
  pi_plus_SE_MC_cuts_hist->SetLineColor(1);
  //pi_plus_SE_MC_cuts_hist->Scale(1./pi_plus_SE_MC_cuts_hist->Integral(firstBin, lastBin));
  pi_plus_SE_MC_cuts_hist->Scale(1./pi_plus_SE_MC_cuts_hist->Integral());
  pi_plus_SE_MC_cuts_hist->GetXaxis()->SetTitle("p_{T} (GeV/#it{c})");
  pi_plus_SE_MC_cuts_hist->GetXaxis()->CenterTitle();
  pi_plus_SE_MC_cuts_hist->GetXaxis()->SetRangeUser(0,1);
  pi_plus_SE_MC_cuts_hist->GetYaxis()->SetRangeUser(0,0.16);
  pi_plus_SE_MC_cuts_hist->Draw("hist");

  TH1D *pi_plus_ME_MC_cuts_hist = L0_L0bar_pi_pT1_vs_pi_pT2_ME_cuts_hist->ProjectionY("pi_plus_ME_MC_cuts_hist");
  pi_plus_ME_MC_cuts_hist->SetLineColor(kRed);
  pi_plus_ME_MC_cuts_hist->Scale(1./pi_plus_ME_MC_cuts_hist->Integral());
  pi_plus_ME_MC_cuts_hist->Draw("hist same");

  SE_vs_ME_leg->Draw("same");

  TPaveText *text_pi_plus_MC_cuts = new TPaveText(0.65, 0.69, 0.89, 0.89, "NDC");
  text_pi_plus_MC_cuts->SetTextFont(42);
  text_pi_plus_MC_cuts->AddText("#pi^{+}");
  text_pi_plus_MC_cuts->AddText("PYTHIA after cuts");
  text_pi_plus_MC_cuts->AddText("Weight after cuts");
  text_pi_plus_MC_cuts->AddText(Form("SE mean = %.3f", pi_plus_SE_MC_cuts_hist->GetMean()));
  text_pi_plus_MC_cuts->AddText(Form("ME mean = %.3f", pi_plus_ME_MC_cuts_hist->GetMean()));
  text_pi_plus_MC_cuts->SetFillColorAlpha(0, 0.01);
  text_pi_plus_MC_cuts->Draw("same");

  pi_plus_SE_vs_ME_MC_cuts->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/pi_pT/pi_plus_SE_vs_ME_MC_cuts.png");


  TCanvas *pi_minus_SE_vs_ME_MC_cuts = new TCanvas("pi_minus_SE_vs_ME_MC_cuts", "pi_minus_SE_vs_ME_MC_cuts", 1200, 1000);

  pi_minus_SE_vs_ME_MC_cuts->cd();

  TH1D *pi_minus_SE_MC_cuts_hist = L0_L0bar_pi_pT1_vs_pi_pT2_cuts_hist->ProjectionX("pi_minus_SE_MC_cuts_hist"); //projection to Lbar - i.e. pi+
  pi_minus_SE_MC_cuts_hist->SetLineColor(1);
  //pi_minus_SE_MC_cuts_hist->Scale(1./pi_minus_SE_MC_cuts_hist->Integral(firstBin, lastBin));
  pi_minus_SE_MC_cuts_hist->Scale(1./pi_minus_SE_MC_cuts_hist->Integral());
  pi_minus_SE_MC_cuts_hist->GetXaxis()->SetTitle("p_{T} (GeV/#it{c})");
  pi_minus_SE_MC_cuts_hist->GetXaxis()->CenterTitle();
  pi_minus_SE_MC_cuts_hist->GetXaxis()->SetRangeUser(0,1);
  pi_minus_SE_MC_cuts_hist->GetYaxis()->SetRangeUser(0,0.16);
  pi_minus_SE_MC_cuts_hist->Draw("hist");

  TH1D *pi_minus_ME_MC_cuts_hist = L0_L0bar_pi_pT1_vs_pi_pT2_ME_cuts_hist->ProjectionX("pi_minus_ME_MC_cuts_hist");
  pi_minus_ME_MC_cuts_hist->SetLineColor(kRed);
  pi_minus_ME_MC_cuts_hist->Scale(1./pi_minus_ME_MC_cuts_hist->Integral());
  pi_minus_ME_MC_cuts_hist->Draw("hist same");

  SE_vs_ME_leg->Draw("same");

  TPaveText *text_pi_minus_MC_cuts = new TPaveText(0.65, 0.69, 0.89, 0.89, "NDC");
  text_pi_minus_MC_cuts->SetTextFont(42);
  text_pi_minus_MC_cuts->AddText("#pi^{-}");
  text_pi_minus_MC_cuts->AddText("PYTHIA after cuts");
  text_pi_minus_MC_cuts->AddText("Weight after cuts");
  text_pi_minus_MC_cuts->AddText(Form("SE mean = %.3f", pi_minus_SE_MC_cuts_hist->GetMean()));
  text_pi_minus_MC_cuts->AddText(Form("ME mean = %.3f", pi_minus_ME_MC_cuts_hist->GetMean()));
  text_pi_minus_MC_cuts->SetFillColorAlpha(0, 0.01);
  text_pi_minus_MC_cuts->Draw("same");

  pi_minus_SE_vs_ME_MC_cuts->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/pi_pT/pi_minus_SE_vs_ME_MC_cuts.png");

  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  //after cuts, ME with weight from distributions before cuts

  TF1 *fit_L0_L0bar_after_cuts_ME_2 = new TF1("fit_L0_L0bar_after_cuts_ME_2", "[0]*(1+[1]*x)", -1, 1);
  fit_L0_L0bar_after_cuts_ME_2->SetParameters(1000, 0.10);

  L0_L0bar_cosThetaProdPlane_cuts_2_can->cd();

  //float nL0L0bar_cuts = L0_L0bar_cosThetaProdPlane_cuts->Integral();

  L0_L0bar_cosThetaProdPlane_ME_weight_2_cuts->Sumw2();
  L0_L0bar_cosThetaProdPlane_ME_weight_2_cuts->Scale(nL0L0bar_cuts/L0_L0bar_cosThetaProdPlane_ME_weight_2_cuts->Integral());
  L0_L0bar_cosThetaProdPlane_ME_weight_2_cuts->Scale(1./L0_L0bar_cosThetaProdPlane_ME_weight_2_cuts->GetXaxis()->GetBinWidth(1));
  L0_L0bar_cosThetaProdPlane_ME_weight_2_cuts->Fit(fit_L0_L0bar_after_cuts_ME_2, "i 0 r");
  //L0_L0bar_cosThetaProdPlane_ME_weight_2->Draw("p e");

  L0_L0bar_cosThetaProdPlane_cuts->Draw("p e");

  L0_L0bar_cosThetaProdPlane_ME_weight_2_cuts->Draw("p e same");

  //fit_L0_L0bar_after_cuts->SetLineColor(1);
  fit_L0_L0bar_after_cuts->Draw("same");

  fit_L0_L0bar_after_cuts_ME_2->SetLineColor(kBlue);
  fit_L0_L0bar_after_cuts_ME_2->Draw("same");


  float L0_L0bar_slope_cuts_sys_err_2 = fabs( fit_L0_L0bar_after_cuts->GetParameter(1) - fit_L0_L0bar_after_cuts_ME_2->GetParameter(1) )/fit_L0_L0bar_after_cuts->GetParameter(1);

  //SysErrSlope_hist->SetBinContent(1, L0_L0bar_slope_cuts_sys_err);

  TPaveText *L0_L0bar_text_MC_cuts_2 = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
  L0_L0bar_text_MC_cuts_2->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  L0_L0bar_text_MC_cuts_2->AddText("Minimum bias");
  L0_L0bar_text_MC_cuts_2->AddText("#Lambda-#bar{#Lambda}");
  L0_L0bar_text_MC_cuts_2->AddText("True MC");
  L0_L0bar_text_MC_cuts_2->AddText("Analysis cuts");
  L0_L0bar_text_MC_cuts_2->AddText(Form("P = %.3f #pm %.3f", L0_L0bar_slope_cuts, L0_L0bar_slope_cuts_err ));
  L0_L0bar_text_MC_cuts_2->AddText(Form("P_{ME} = %.3f", fit_L0_L0bar_after_cuts_ME_2->GetParameter(1)/(L0_alpha*L0bar_alpha) ));
  L0_L0bar_text_MC_cuts_2->AddText(Form("Slope sys. err: %.1f %%", L0_L0bar_slope_cuts_sys_err_2*100 ));
  L0_L0bar_text_MC_cuts_2->SetFillColorAlpha(0, 0.01);
  L0_L0bar_text_MC_cuts_2->Draw("same");

  L0_L0bar_cosThetaProdPlane_cuts_2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/L0_L0bar_cosThetaProdPlane_cuts_weight_no_cuts.png");

  //pi kinematics
  TCanvas *pi_plus_SE_vs_ME_MC_cuts_2 = new TCanvas("pi_plus_SE_vs_ME_MC_cuts_2", "pi_plus_SE_vs_ME_MC_cuts_2", 1200, 1000);

  pi_plus_SE_vs_ME_MC_cuts_2->cd();

  pi_plus_SE_MC_cuts_hist->Draw("hist");

  TH1D *pi_plus_ME_MC_cuts_2_hist = L0_L0bar_pi_pT1_vs_pi_pT2_ME_cuts_2_hist->ProjectionY("pi_plus_ME_MC_cuts_2_hist");
  pi_plus_ME_MC_cuts_2_hist->SetLineColor(kRed);
  pi_plus_ME_MC_cuts_2_hist->Scale(1./pi_plus_ME_MC_cuts_2_hist->Integral());
  pi_plus_ME_MC_cuts_2_hist->Draw("hist same");

  SE_vs_ME_leg->Draw("same");

  TPaveText *text_pi_plus_MC_cuts_2 = new TPaveText(0.65, 0.69, 0.89, 0.89, "NDC");
  text_pi_plus_MC_cuts_2->SetTextFont(42);
  text_pi_plus_MC_cuts_2->AddText("#pi^{+}");
  text_pi_plus_MC_cuts_2->AddText("PYTHIA after cuts");
  text_pi_plus_MC_cuts_2->AddText("Weight before cuts");
  text_pi_plus_MC_cuts_2->AddText(Form("SE mean = %.3f", pi_plus_SE_MC_cuts_hist->GetMean()));
  text_pi_plus_MC_cuts_2->AddText(Form("ME mean = %.3f", pi_plus_ME_MC_cuts_2_hist->GetMean()));
  text_pi_plus_MC_cuts_2->SetFillColorAlpha(0, 0.01);
  text_pi_plus_MC_cuts_2->Draw("same");

  pi_plus_SE_vs_ME_MC_cuts_2->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/pi_pT/pi_plus_SE_vs_ME_MC_cuts_before_cut_weight.png");


  TCanvas *pi_minus_SE_vs_ME_MC_cuts_2 = new TCanvas("pi_minus_SE_vs_ME_MC_cuts_2", "pi_minus_SE_vs_ME_MC_cuts_2", 1200, 1000);

  pi_minus_SE_vs_ME_MC_cuts_2->cd();

  pi_minus_SE_MC_cuts_hist->Draw("hist");

  TH1D *pi_minus_ME_MC_cuts_2_hist = L0_L0bar_pi_pT1_vs_pi_pT2_ME_cuts_2_hist->ProjectionX("pi_minus_ME_MC_cuts_2_hist");
  pi_minus_ME_MC_cuts_2_hist->SetLineColor(kRed);
  pi_minus_ME_MC_cuts_2_hist->Scale(1./pi_minus_ME_MC_cuts_2_hist->Integral());
  pi_minus_ME_MC_cuts_2_hist->Draw("hist same");

  SE_vs_ME_leg->Draw("same");

  TPaveText *text_pi_minus_MC_cuts_2 = new TPaveText(0.65, 0.69, 0.89, 0.89, "NDC");
  text_pi_minus_MC_cuts_2->SetTextFont(42);
  text_pi_minus_MC_cuts_2->AddText("#pi^{-}");
  text_pi_minus_MC_cuts_2->AddText("PYTHIA after cuts");
  text_pi_minus_MC_cuts_2->AddText("Weight before cuts");
  text_pi_minus_MC_cuts_2->AddText(Form("SE mean = %.3f", pi_minus_SE_MC_cuts_hist->GetMean()));
  text_pi_minus_MC_cuts_2->AddText(Form("ME mean = %.3f", pi_minus_ME_MC_cuts_2_hist->GetMean()));
  text_pi_minus_MC_cuts_2->SetFillColorAlpha(0, 0.01);
  text_pi_minus_MC_cuts_2->Draw("same");

  pi_minus_SE_vs_ME_MC_cuts_2->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/pi_pT/pi_minus_SE_vs_ME_MC_cuts_before_cut_weight.png");

  //___________________________________________________________________________________________________________________________________________________________________________________________________________________________

  //L-L
  //before cuts
  TF1 *fit_L0_L0_before_cuts = new TF1("fit_L0_L0_before_cuts", "[0]*(1+[1]*x)", -1, 1);
  fit_L0_L0_before_cuts->SetParameters(1000, 0.10);

  L0_L0_cosThetaProdPlane_can->cd();

  float nL0L0 = L0_L0_cosThetaProdPlane->Integral();

  L0_L0_cosThetaProdPlane_ME_weight->Sumw2();
  //L0_L0_cosThetaProdPlane_ME_weight->Add(L0_L0_cosThetaProdPlane_LS_ME_weight);
  L0_L0_cosThetaProdPlane_ME_weight->Scale(nL0L0/L0_L0_cosThetaProdPlane_ME_weight->Integral());
  //L0_L0_cosThetaProdPlane_ME_weight->Draw("p e");

  L0_L0_cosThetaProdPlane->Sumw2();
  L0_L0_cosThetaProdPlane->SetMarkerStyle(20);
  L0_L0_cosThetaProdPlane->SetMarkerColor(kRed);
  L0_L0_cosThetaProdPlane->SetLineColor(kRed);
  L0_L0_cosThetaProdPlane->GetXaxis()->SetTitle("cos(#theta*)");
  L0_L0_cosThetaProdPlane->GetXaxis()->CenterTitle();
  L0_L0_cosThetaProdPlane->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0_L0_cosThetaProdPlane->GetYaxis()->CenterTitle();
  L0_L0_cosThetaProdPlane->Divide(L0_L0_cosThetaProdPlane_ME_weight); //correct using ME
  L0_L0_cosThetaProdPlane->Scale(nL0L0/L0_L0_cosThetaProdPlane->Integral()); //scale back
  //L0_L0_cosThetaProdPlane->Scale(1./L0_L0_cosThetaProdPlane->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
  L0_L0_cosThetaProdPlane->Fit(fit_L0_L0_before_cuts, "i 0 r");
  L0_L0_cosThetaProdPlane->SetMinimum(0);
  //L0_L0_cosThetaProdPlane->GetYaxis()->SetRangeUser(-20e9, 0);
  L0_L0_cosThetaProdPlane->Draw("p e");

  //L0_L0_cosThetaProdPlane_ME_weight->Draw("p e same");

  fit_L0_L0_before_cuts->SetLineColor(1);
  fit_L0_L0_before_cuts->Draw("same");

  float L0_L0_slope = fit_L0_L0_before_cuts->GetParameter(1)/(L0_alpha*L0_alpha);
  float L0_L0_slope_err = fit_L0_L0_before_cuts->GetParError(1)/(L0_alpha*L0_alpha);

  TPaveText *L0_L0_text_MC = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
  L0_L0_text_MC->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  L0_L0_text_MC->AddText("Minimum bias");
  L0_L0_text_MC->AddText("#Lambda-#Lambda");
  L0_L0_text_MC->AddText("True MC");
  L0_L0_text_MC->AddText("No cuts");
  L0_L0_text_MC->AddText(Form("P = %.3f #pm %.3f", L0_L0_slope, L0_L0_slope_err));
  L0_L0_text_MC->SetFillColorAlpha(0, 0.01);
  L0_L0_text_MC->Draw("same");

  L0_L0_cosThetaProdPlane_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_L/L0_L0_cosThetaProdPlane.png");

  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  //after cuts
  TF1 *fit_L0_L0_after_cuts = new TF1("fit_L0_L0_after_cuts", "[0]*(1+[1]*x)", -1, 1);
  fit_L0_L0_after_cuts->SetParameters(1000, 0.10);

  TF1 *fit_L0_L0_after_cuts_ME = new TF1("fit_L0_L0_after_cuts_ME", "[0]*(1+[1]*x)", -1, 1);
  fit_L0_L0_after_cuts_ME->SetParameters(1000, 0.10);

  L0_L0_cosThetaProdPlane_cuts_can->cd();

  float nL0L0_cuts = L0_L0_cosThetaProdPlane_cuts->Integral();

  L0_L0_cosThetaProdPlane_ME_weight_cuts->Sumw2();
  L0_L0_cosThetaProdPlane_ME_weight_cuts->Scale(nL0L0_cuts/L0_L0_cosThetaProdPlane_ME_weight_cuts->Integral());
  L0_L0_cosThetaProdPlane_ME_weight_cuts->Fit(fit_L0_L0_after_cuts_ME, "i 0 r");
  //L0_L0_cosThetaProdPlane_ME_weight->Draw("p e");

  L0_L0_cosThetaProdPlane_cuts->Sumw2();
  L0_L0_cosThetaProdPlane_cuts->SetMarkerStyle(20);
  L0_L0_cosThetaProdPlane_cuts->SetMarkerColor(kRed);
  L0_L0_cosThetaProdPlane_cuts->SetLineColor(kRed);
  L0_L0_cosThetaProdPlane_cuts->GetXaxis()->SetTitle("cos(#theta*)");
  L0_L0_cosThetaProdPlane_cuts->GetXaxis()->CenterTitle();
  L0_L0_cosThetaProdPlane_cuts->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0_L0_cosThetaProdPlane_cuts->GetYaxis()->CenterTitle();
  //L0_L0_cosThetaProdPlane_cuts->Divide(L0_L0_cosThetaProdPlane_ME_weight_cuts); //correct using ME
  //L0_L0_cosThetaProdPlane_cuts->Scale(nL0L0_cuts/L0_L0_cosThetaProdPlane_cuts->Integral()); //scale back
  //L0_L0_cosThetaProdPlane_cuts->Scale(1./L0_L0_cosThetaProdPlane_cuts->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
  L0_L0_cosThetaProdPlane_cuts->Fit(fit_L0_L0_after_cuts, "i 0 r");
  L0_L0_cosThetaProdPlane_cuts->SetMinimum(0);
  //L0_L0_cosThetaProdPlane_cuts->GetYaxis()->SetRangeUser(-20e9, 0);
  L0_L0_cosThetaProdPlane_cuts->Draw("p e");

  L0_L0_cosThetaProdPlane_ME_weight_cuts->Draw("p e same");

  fit_L0_L0_after_cuts->SetLineColor(1);
  fit_L0_L0_after_cuts->Draw("same");

  float L0_L0_slope_cuts = fit_L0_L0_after_cuts->GetParameter(1)/(L0_alpha*L0_alpha);
  float L0_L0_slope_cuts_err = fit_L0_L0_after_cuts->GetParError(1)/(L0_alpha*L0_alpha);

  float L0_L0_slope_cuts_sys_err = fabs( fit_L0_L0_after_cuts->GetParameter(1) - fit_L0_L0_after_cuts_ME->GetParameter(1) )/fabs(fit_L0_L0_after_cuts->GetParameter(1));

  SysErrSlope_hist->SetBinContent(2, L0_L0_slope_cuts_sys_err);

  TPaveText *L0_L0_text_MC_cuts = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
  L0_L0_text_MC_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  L0_L0_text_MC_cuts->AddText("Minimum bias");
  L0_L0_text_MC_cuts->AddText("#Lambda-#Lambda");
  L0_L0_text_MC_cuts->AddText("True MC");
  L0_L0_text_MC_cuts->AddText("Analysis cuts");
  L0_L0_text_MC_cuts->AddText(Form("P = %.3f #pm %.3f", L0_L0_slope_cuts, L0_L0_slope_cuts_err));
  L0_L0_text_MC_cuts->AddText(Form("Sys. err: %.1f %%", L0_L0_slope_cuts_sys_err*100 ));
  L0_L0_text_MC_cuts->SetFillColorAlpha(0, 0.01);
  L0_L0_text_MC_cuts->Draw("same");

  L0_L0_cosThetaProdPlane_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_L/L0_L0_cosThetaProdPlane_cuts.png");

  //--------------------------------------------------------------------

  //after cuts, ME with weight from distributions before cuts

  TF1 *fit_L0_L0_after_cuts_ME_2 = new TF1("fit_L0_L0_after_cuts_ME_2", "[0]*(1+[1]*x)", -1, 1);
  fit_L0_L0_after_cuts_ME_2->SetParameters(1000, 0.10);

  L0_L0_cosThetaProdPlane_cuts_2_can->cd();

  //float nL0L0_cuts = L0_L0_cosThetaProdPlane_cuts->Integral();

  L0_L0_cosThetaProdPlane_ME_weight_2_cuts->Sumw2();
  L0_L0_cosThetaProdPlane_ME_weight_2_cuts->Scale(nL0L0_cuts/L0_L0_cosThetaProdPlane_ME_weight_2_cuts->Integral());
  //L0_L0_cosThetaProdPlane_ME_weight_2_cuts->Scale(1./L0_L0_cosThetaProdPlane_ME_weight_2_cuts->GetXaxis()->GetBinWidth(1));
  L0_L0_cosThetaProdPlane_ME_weight_2_cuts->Fit(fit_L0_L0_after_cuts_ME_2, "i 0 r");
  //L0_L0_cosThetaProdPlane_ME_weight_2->Draw("p e");

  L0_L0_cosThetaProdPlane_cuts->Draw("p e");

  L0_L0_cosThetaProdPlane_ME_weight_2_cuts->Draw("p e same");

  //fit_L0_L0_after_cuts->SetLineColor(1);
  fit_L0_L0_after_cuts->Draw("same");

  fit_L0_L0_after_cuts_ME_2->SetLineColor(kBlue);
  fit_L0_L0_after_cuts_ME_2->Draw("same");


  float L0_L0_slope_cuts_sys_err_2 = fabs( fit_L0_L0_after_cuts->GetParameter(1) - fit_L0_L0_after_cuts_ME_2->GetParameter(1) )/fabs(fit_L0_L0_after_cuts->GetParameter(1));

  //SysErrSlope_hist->SetBinContent(1, L0_L0_slope_cuts_sys_err);

  TPaveText *L0_L0_text_MC_cuts_2 = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
  L0_L0_text_MC_cuts_2->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  L0_L0_text_MC_cuts_2->AddText("Minimum bias");
  L0_L0_text_MC_cuts_2->AddText("#Lambda-#bar{#Lambda}");
  L0_L0_text_MC_cuts_2->AddText("True MC");
  L0_L0_text_MC_cuts_2->AddText("Analysis cuts");
  L0_L0_text_MC_cuts_2->AddText(Form("P = %.3f #pm %.3f", L0_L0_slope_cuts, L0_L0_slope_cuts_err ));
  L0_L0_text_MC_cuts_2->AddText(Form("P_{ME} = %.3f", fit_L0_L0_after_cuts_ME_2->GetParameter(1)/(L0_alpha*L0_alpha) ));
  L0_L0_text_MC_cuts_2->AddText(Form("Slope sys. err: %.1f %%", L0_L0_slope_cuts_sys_err_2*100 ));
  L0_L0_text_MC_cuts_2->SetFillColorAlpha(0, 0.01);
  L0_L0_text_MC_cuts_2->Draw("same");

  L0_L0_cosThetaProdPlane_cuts_2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_L/L0_L0_cosThetaProdPlane_cuts_weight_no_cuts.png");

  //_________________________________________________________________________

  //Lbar-Lbar
  //before cuts
  TF1 *fit_L0bar_L0bar_before_cuts = new TF1("fit_L0bar_L0bar_before_cuts", "[0]*(1+[1]*x)", -1, 1);
  fit_L0bar_L0bar_before_cuts->SetParameters(1000, 0.10);

  L0bar_L0bar_cosThetaProdPlane_can->cd();

  float nL0barL0bar = L0bar_L0bar_cosThetaProdPlane->Integral();

  L0bar_L0bar_cosThetaProdPlane_ME_weight->Sumw2();
  //L0bar_L0bar_cosThetaProdPlane_ME_weight->Add(L0bar_L0bar_cosThetaProdPlane_LS_ME_weight);
  L0bar_L0bar_cosThetaProdPlane_ME_weight->Scale(nL0barL0bar/L0bar_L0bar_cosThetaProdPlane_ME_weight->Integral());
  //L0bar_L0bar_cosThetaProdPlane_ME_weight->Draw("p e");

  L0bar_L0bar_cosThetaProdPlane->Sumw2();
  L0bar_L0bar_cosThetaProdPlane->SetMarkerStyle(20);
  L0bar_L0bar_cosThetaProdPlane->SetMarkerColor(kRed);
  L0bar_L0bar_cosThetaProdPlane->SetLineColor(kRed);
  L0bar_L0bar_cosThetaProdPlane->GetXaxis()->SetTitle("cos(#theta*)");
  L0bar_L0bar_cosThetaProdPlane->GetXaxis()->CenterTitle();
  L0bar_L0bar_cosThetaProdPlane->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0bar_L0bar_cosThetaProdPlane->GetYaxis()->CenterTitle();
  L0bar_L0bar_cosThetaProdPlane->Divide(L0bar_L0bar_cosThetaProdPlane_ME_weight); //correct using ME
  L0bar_L0bar_cosThetaProdPlane->Scale(nL0barL0bar/L0bar_L0bar_cosThetaProdPlane->Integral()); //scale back
  //L0bar_L0bar_cosThetaProdPlane->Scale(1./L0bar_L0bar_cosThetaProdPlane->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
  L0bar_L0bar_cosThetaProdPlane->Fit(fit_L0bar_L0bar_before_cuts, "i 0 r");
  L0bar_L0bar_cosThetaProdPlane->SetMinimum(0);
  //L0bar_L0bar_cosThetaProdPlane->GetYaxis()->SetRangeUser(-20e9, 0);
  L0bar_L0bar_cosThetaProdPlane->Draw("p e");

  //L0bar_L0bar_cosThetaProdPlane_ME_weight->Draw("p e same");

  fit_L0bar_L0bar_before_cuts->SetLineColor(1);
  fit_L0bar_L0bar_before_cuts->Draw("same");


  float L0bar_L0bar_slope = fit_L0bar_L0bar_before_cuts->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
  float L0bar_L0bar_slope_err = fit_L0bar_L0bar_before_cuts->GetParError(1)/(L0bar_alpha*L0bar_alpha);

  TPaveText *L0bar_L0bar_text_MC = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
  L0bar_L0bar_text_MC->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  L0bar_L0bar_text_MC->AddText("Minimum bias");
  L0bar_L0bar_text_MC->AddText("#bar{#Lambda}-#bar{#Lambda}");
  L0bar_L0bar_text_MC->AddText("True MC");
  L0bar_L0bar_text_MC->AddText("No cuts");
  L0bar_L0bar_text_MC->AddText(Form("P = %.3f #pm %.3f", L0bar_L0bar_slope, L0bar_L0bar_slope_err));
  L0bar_L0bar_text_MC->SetFillColorAlpha(0, 0.01);
  L0bar_L0bar_text_MC->Draw("same");

  L0bar_L0bar_cosThetaProdPlane_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/Lbar_Lbar/L0bar_L0bar_cosThetaProdPlane.png");

  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  //after cuts
  TF1 *fit_L0bar_L0bar_after_cuts = new TF1("fit_L0bar_L0bar_after_cuts", "[0]*(1+[1]*x)", -1, 1);
  fit_L0bar_L0bar_after_cuts->SetParameters(1000, 0.10);

  TF1 *fit_L0bar_L0bar_after_cuts_ME = new TF1("fit_L0bar_L0bar_after_cuts_ME", "[0]*(1+[1]*x)", -1, 1);
  fit_L0bar_L0bar_after_cuts_ME->SetParameters(1000, 0.10);

  L0bar_L0bar_cosThetaProdPlane_cuts_can->cd();

  float nL0barL0bar_cuts = L0bar_L0bar_cosThetaProdPlane_cuts->Integral();

  L0bar_L0bar_cosThetaProdPlane_ME_weight_cuts->Sumw2();
  L0bar_L0bar_cosThetaProdPlane_ME_weight_cuts->Scale(nL0barL0bar_cuts/L0bar_L0bar_cosThetaProdPlane_ME_weight_cuts->Integral());
  L0bar_L0bar_cosThetaProdPlane_ME_weight_cuts->Fit(fit_L0bar_L0bar_after_cuts_ME, "i 0 r");
  //L0bar_L0bar_cosThetaProdPlane_ME_weight->Draw("p e");

  L0bar_L0bar_cosThetaProdPlane_cuts->Sumw2();
  L0bar_L0bar_cosThetaProdPlane_cuts->SetMarkerStyle(20);
  L0bar_L0bar_cosThetaProdPlane_cuts->SetMarkerColor(kRed);
  L0bar_L0bar_cosThetaProdPlane_cuts->SetLineColor(kRed);
  L0bar_L0bar_cosThetaProdPlane_cuts->GetXaxis()->SetTitle("cos(#theta*)");
  L0bar_L0bar_cosThetaProdPlane_cuts->GetXaxis()->CenterTitle();
  L0bar_L0bar_cosThetaProdPlane_cuts->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0bar_L0bar_cosThetaProdPlane_cuts->GetYaxis()->CenterTitle();
  //L0bar_L0bar_cosThetaProdPlane_cuts->Divide(L0bar_L0bar_cosThetaProdPlane_ME_weight_cuts); //correct using ME
  //L0bar_L0bar_cosThetaProdPlane_cuts->Scale(nL0barL0bar_cuts/L0bar_L0bar_cosThetaProdPlane_cuts->Integral()); //scale back
  //L0bar_L0bar_cosThetaProdPlane_cuts->Scale(1./L0bar_L0bar_cosThetaProdPlane_cuts->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
  L0bar_L0bar_cosThetaProdPlane_cuts->Fit(fit_L0bar_L0bar_after_cuts, "i 0 r");
  L0bar_L0bar_cosThetaProdPlane_cuts->SetMinimum(0);
  //L0bar_L0bar_cosThetaProdPlane_cuts->GetYaxis()->SetRangeUser(-20e9, 0);
  L0bar_L0bar_cosThetaProdPlane_cuts->Draw("p e");

  L0bar_L0bar_cosThetaProdPlane_ME_weight_cuts->Draw("p e same");

  fit_L0bar_L0bar_after_cuts->SetLineColor(1);
  fit_L0bar_L0bar_after_cuts->Draw("same");

  float L0bar_L0bar_slope_cuts = fit_L0bar_L0bar_after_cuts->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
  float L0bar_L0bar_slope_cuts_err = fit_L0bar_L0bar_after_cuts->GetParError(1)/(L0bar_alpha*L0bar_alpha);

  float L0bar_L0bar_slope_cuts_sys_err = fabs( fit_L0bar_L0bar_after_cuts->GetParameter(1) - fit_L0bar_L0bar_after_cuts_ME->GetParameter(1) )/fabs(fit_L0bar_L0bar_after_cuts->GetParameter(1));

  SysErrSlope_hist->SetBinContent(3, L0bar_L0bar_slope_cuts_sys_err);

  TPaveText *L0bar_L0bar_text_MC_cuts = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
  L0bar_L0bar_text_MC_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  L0bar_L0bar_text_MC_cuts->AddText("Minimum bias");
  L0bar_L0bar_text_MC_cuts->AddText("#bar{#Lambda}-#bar{#Lambda}");
  L0bar_L0bar_text_MC_cuts->AddText("True MC");
  L0bar_L0bar_text_MC_cuts->AddText("Analysis cuts");
  L0bar_L0bar_text_MC_cuts->AddText(Form("P = %.3f #pm %.3f", L0bar_L0bar_slope_cuts, L0bar_L0bar_slope_cuts_err));
  L0bar_L0bar_text_MC_cuts->AddText(Form("Slope sys. err: %.1f %%", L0bar_L0bar_slope_cuts_sys_err*100 ));
  L0bar_L0bar_text_MC_cuts->SetFillColorAlpha(0, 0.01);
  L0bar_L0bar_text_MC_cuts->Draw("same");

  L0bar_L0bar_cosThetaProdPlane_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/Lbar_Lbar/L0bar_L0bar_cosThetaProdPlane_cuts.png");

  //--------------------------------------------------------------------

  //after cuts, ME with weight from distributions before cuts

  TF1 *fit_L0bar_L0bar_after_cuts_ME_2 = new TF1("fit_L0bar_L0bar_after_cuts_ME_2", "[0]*(1+[1]*x)", -1, 1);
  fit_L0bar_L0bar_after_cuts_ME_2->SetParameters(1000, 0.10);

  L0bar_L0bar_cosThetaProdPlane_cuts_2_can->cd();

  //float nL0barL0bar_cuts = L0bar_L0bar_cosThetaProdPlane_cuts->Integral();

  L0bar_L0bar_cosThetaProdPlane_ME_weight_2_cuts->Sumw2();
  L0bar_L0bar_cosThetaProdPlane_ME_weight_2_cuts->Scale(nL0barL0bar_cuts/L0bar_L0bar_cosThetaProdPlane_ME_weight_2_cuts->Integral());
  //L0bar_L0bar_cosThetaProdPlane_ME_weight_2_cuts->Scale(1./L0bar_L0bar_cosThetaProdPlane_ME_weight_2_cuts->GetXaxis()->GetBinWidth(1));
  L0bar_L0bar_cosThetaProdPlane_ME_weight_2_cuts->Fit(fit_L0bar_L0bar_after_cuts_ME_2, "i 0 r");
  //L0bar_L0bar_cosThetaProdPlane_ME_weight_2->Draw("p e");

  L0bar_L0bar_cosThetaProdPlane_cuts->Draw("p e");

  L0bar_L0bar_cosThetaProdPlane_ME_weight_2_cuts->Draw("p e same");

  //fit_L0bar_L0bar_after_cuts->SetLineColor(1);
  fit_L0bar_L0bar_after_cuts->Draw("same");

  fit_L0bar_L0bar_after_cuts_ME_2->SetLineColor(kBlue);
  fit_L0bar_L0bar_after_cuts_ME_2->Draw("same");


  float L0bar_L0bar_slope_cuts_sys_err_2 = fabs( fit_L0bar_L0bar_after_cuts->GetParameter(1) - fit_L0bar_L0bar_after_cuts_ME_2->GetParameter(1) )/fabs(fit_L0bar_L0bar_after_cuts->GetParameter(1));

  //SysErrSlope_hist->SetBinContent(1, L0bar_L0bar_slope_cuts_sys_err);

  TPaveText *L0bar_L0bar_text_MC_cuts_2 = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
  L0bar_L0bar_text_MC_cuts_2->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  L0bar_L0bar_text_MC_cuts_2->AddText("Minimum bias");
  L0bar_L0bar_text_MC_cuts_2->AddText("#Lambda-#bar{#Lambda}");
  L0bar_L0bar_text_MC_cuts_2->AddText("True MC");
  L0bar_L0bar_text_MC_cuts_2->AddText("Analysis cuts");
  L0bar_L0bar_text_MC_cuts_2->AddText(Form("P = %.3f #pm %.3f", L0bar_L0bar_slope_cuts, L0bar_L0bar_slope_cuts_err ));
  L0bar_L0bar_text_MC_cuts_2->AddText(Form("P_{ME} = %.3f", fit_L0bar_L0bar_after_cuts_ME_2->GetParameter(1)/(L0bar_alpha*L0bar_alpha) ));
  L0bar_L0bar_text_MC_cuts_2->AddText(Form("Slope sys. err: %.1f %%", L0bar_L0bar_slope_cuts_sys_err_2*100 ));
  L0bar_L0bar_text_MC_cuts_2->SetFillColorAlpha(0, 0.01);
  L0bar_L0bar_text_MC_cuts_2->Draw("same");

  L0bar_L0bar_cosThetaProdPlane_cuts_2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/Lbar_Lbar/L0bar_L0bar_cosThetaProdPlane_cuts_weight_no_cuts.png");

  //_________________________________________________________________________________________________________________________________________________________________________

  //US matched with US
  //L-Lbar
  L0_L0bar_cosThetaProdPlane_US_no_corr_can->cd();

  L0_L0bar_cosThetaProdPlane_US->Sumw2();
  L0_L0bar_cosThetaProdPlane_US->SetMarkerSize(1.5);
  L0_L0bar_cosThetaProdPlane_US->SetMarkerStyle(20);
  L0_L0bar_cosThetaProdPlane_US->SetMarkerColor(kRed);
  L0_L0bar_cosThetaProdPlane_US->SetLineColor(kRed);
  L0_L0bar_cosThetaProdPlane_US->GetXaxis()->SetTitle("cos(#theta*)");
  L0_L0bar_cosThetaProdPlane_US->GetXaxis()->CenterTitle();
  L0_L0bar_cosThetaProdPlane_US->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0_L0bar_cosThetaProdPlane_US->GetYaxis()->CenterTitle();
  float nL0L0bar_US = L0_L0bar_cosThetaProdPlane_US->Integral();
  L0_L0bar_cosThetaProdPlane_US->Scale(1./L0_L0bar_cosThetaProdPlane_US->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
  L0_L0bar_cosThetaProdPlane_US->SetMinimum(0);
  L0_L0bar_cosThetaProdPlane_US->Draw("p e");

  L0_L0bar_cosThetaProdPlane_US_ME_weight->Sumw2();
  L0_L0bar_cosThetaProdPlane_US_ME_weight->SetMarkerSize(1.5);
  L0_L0bar_cosThetaProdPlane_US_ME_weight->SetMarkerStyle(24);
  L0_L0bar_cosThetaProdPlane_US_ME_weight->SetMarkerColor(1);
  L0_L0bar_cosThetaProdPlane_US_ME_weight->SetLineColor(1);
  //L0_L0bar_cosThetaProdPlane_US_ME_weight->Add(L0_L0bar_cosThetaProdPlane_US_LS_ME_weight);
  L0_L0bar_cosThetaProdPlane_US_ME_weight->Scale(nL0L0bar_US/L0_L0bar_cosThetaProdPlane_US_ME_weight->Integral());
  L0_L0bar_cosThetaProdPlane_US_ME_weight->Scale(1./L0_L0bar_cosThetaProdPlane_US_ME_weight->GetXaxis()->GetBinWidth(1));
  L0_L0bar_cosThetaProdPlane_US_ME_weight->Draw("p e same");


  L0_L0bar_cosThetaProdPlane_US_LS->Sumw2();
  L0_L0bar_cosThetaProdPlane_US_LS->SetMarkerSize(1.5);
  L0_L0bar_cosThetaProdPlane_US_LS->SetMarkerStyle(21);
  L0_L0bar_cosThetaProdPlane_US_LS->SetMarkerColor(kBlue);
  L0_L0bar_cosThetaProdPlane_US_LS->SetLineColor(kBlue);
  float nL0L0bar_US_LS = L0_L0bar_cosThetaProdPlane_US_LS->Integral();
  L0_L0bar_cosThetaProdPlane_US_LS->Scale(1./L0_L0bar_cosThetaProdPlane_US_LS->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
  L0_L0bar_cosThetaProdPlane_US_LS->Draw("p e same");

  L0_L0bar_cosThetaProdPlane_US_LS_ME_weight->Sumw2();
  L0_L0bar_cosThetaProdPlane_US_LS_ME_weight->SetMarkerSize(1.5);
  L0_L0bar_cosThetaProdPlane_US_LS_ME_weight->SetMarkerStyle(25);
  L0_L0bar_cosThetaProdPlane_US_LS_ME_weight->SetMarkerColor(kMagenta+1);
  L0_L0bar_cosThetaProdPlane_US_LS_ME_weight->SetLineColor(kMagenta+1);
  L0_L0bar_cosThetaProdPlane_US_LS_ME_weight->Scale(nL0L0bar_US_LS/L0_L0bar_cosThetaProdPlane_US_LS_ME_weight->Integral()); //scale ME
  L0_L0bar_cosThetaProdPlane_US_LS_ME_weight->Scale(1./L0_L0bar_cosThetaProdPlane_US_LS_ME_weight->GetXaxis()->GetBinWidth(1));
  L0_L0bar_cosThetaProdPlane_US_LS_ME_weight->Draw("p e same");


  TPaveText *L0_L0bar_text_no_corr = new TPaveText(0.6, 0.35, 0.8, 0.55, "NDC");
  L0_L0bar_text_no_corr->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  L0_L0bar_text_no_corr->AddText("Minimum bias");
  L0_L0bar_text_no_corr->AddText("#Lambda-#bar{#Lambda}");
  L0_L0bar_text_no_corr->AddText("No correction");
  L0_L0bar_text_no_corr->AddText("No cuts");
  L0_L0bar_text_no_corr->SetFillColorAlpha(0, 0.01);
  L0_L0bar_text_no_corr->Draw("same");


  L0_L0bar_cosThetaProdPlane_US_no_corr_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/L0_L0bar_cosThetaProdPlane_US_no_corr.png");

  //---------------------------------------

  L0_L0bar_cosThetaProdPlane_US_separate_can->cd();

  L0_L0bar_cosThetaProdPlane_US->Divide(L0_L0bar_cosThetaProdPlane_US_ME_weight); //correct using ME
  L0_L0bar_cosThetaProdPlane_US->Scale(nL0L0bar_US/L0_L0bar_cosThetaProdPlane_US->Integral()); //scale back
  L0_L0bar_cosThetaProdPlane_US->Scale(1./L0_L0bar_cosThetaProdPlane_US->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
  L0_L0bar_cosThetaProdPlane_US->SetMinimum(0);
  L0_L0bar_cosThetaProdPlane_US->Draw("p e");


  L0_L0bar_cosThetaProdPlane_US_LS->Divide(L0_L0bar_cosThetaProdPlane_US_LS_ME_weight);
  L0_L0bar_cosThetaProdPlane_US_LS->Scale(nL0L0bar_US_LS/L0_L0bar_cosThetaProdPlane_US_LS->Integral()); //scale back
  L0_L0bar_cosThetaProdPlane_US_LS->Scale(1./L0_L0bar_cosThetaProdPlane_US_LS->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
  //L0_L0bar_cosThetaProdPlane_US_LS->Draw("p e same");

  TPaveText *L0_L0bar_text_separate = new TPaveText(0.6, 0.35, 0.8, 0.55, "NDC");
  L0_L0bar_text_separate->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  L0_L0bar_text_separate->AddText("Minimum bias");
  L0_L0bar_text_separate->AddText("#Lambda-#bar{#Lambda}");
  L0_L0bar_text_separate->AddText("ME corrected");
  L0_L0bar_text_separate->AddText("No cuts");
  L0_L0bar_text_separate->SetFillColorAlpha(0, 0.01);
  L0_L0bar_text_separate->Draw("same");


  L0_L0bar_cosThetaProdPlane_US_separate_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/L0_L0bar_cosThetaProdPlane_US_separate.png");

  //---------------------------------------

  L0_L0bar_cosThetaProdPlane_US_can->cd();

  //correct signal+background

  L0_L0bar_cosThetaProdPlane_US->Add(L0_L0bar_cosThetaProdPlane_US_LS, -1); //subtract combinatorial background
  //L0_L0bar_cosThetaProdPlane_US->Scale(1./L0_L0bar_cosThetaProdPlane_US->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
  L0_L0bar_cosThetaProdPlane_US->SetMinimum(0);
  L0_L0bar_cosThetaProdPlane_US->Draw("p e");


  TPaveText *L0_L0bar_text = new TPaveText(0.6, 0.35, 0.8, 0.55, "NDC");
  L0_L0bar_text->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  L0_L0bar_text->AddText("Minimum bias");
  L0_L0bar_text->AddText("#Lambda-#bar{#Lambda}");
  L0_L0bar_text->AddText("US - Background");
  L0_L0bar_text->AddText("No cuts");
  L0_L0bar_text->SetFillColorAlpha(0, 0.01);
  L0_L0bar_text->Draw("same");

  L0_L0bar_cosThetaProdPlane_US_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/L0_L0bar_cosThetaProdPlane_US.png");

  //-------------------------------------------------------------------------



  L0_L0bar_cosThetaProdPlane_US_no_corr_cuts_can->cd();

  L0_L0bar_cosThetaProdPlane_US_cuts->Sumw2();
  L0_L0bar_cosThetaProdPlane_US_cuts->SetMarkerSize(1.5);
  L0_L0bar_cosThetaProdPlane_US_cuts->SetMarkerStyle(20);
  L0_L0bar_cosThetaProdPlane_US_cuts->SetMarkerColor(kRed);
  L0_L0bar_cosThetaProdPlane_US_cuts->SetLineColor(kRed);
  L0_L0bar_cosThetaProdPlane_US_cuts->GetXaxis()->SetTitle("cos(#theta*)");
  L0_L0bar_cosThetaProdPlane_US_cuts->GetXaxis()->CenterTitle();
  L0_L0bar_cosThetaProdPlane_US_cuts->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0_L0bar_cosThetaProdPlane_US_cuts->GetYaxis()->CenterTitle();
  float nL0L0bar_US_cuts = L0_L0bar_cosThetaProdPlane_US_cuts->Integral();
  L0_L0bar_cosThetaProdPlane_US_cuts->Scale(1./L0_L0bar_cosThetaProdPlane_US_cuts->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
  L0_L0bar_cosThetaProdPlane_US_cuts->SetMinimum(0);
  L0_L0bar_cosThetaProdPlane_US_cuts->Draw("p e");

  L0_L0bar_cosThetaProdPlane_US_ME_weight_cuts->Sumw2();
  L0_L0bar_cosThetaProdPlane_US_ME_weight_cuts->SetMarkerSize(1.5);
  L0_L0bar_cosThetaProdPlane_US_ME_weight_cuts->SetMarkerStyle(24);
  L0_L0bar_cosThetaProdPlane_US_ME_weight_cuts->SetMarkerColor(1);
  L0_L0bar_cosThetaProdPlane_US_ME_weight_cuts->SetLineColor(1);
  //L0_L0bar_cosThetaProdPlane_US_ME_weight_cuts->Add(L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts);
  L0_L0bar_cosThetaProdPlane_US_ME_weight_cuts->Scale(nL0L0bar_US_cuts/L0_L0bar_cosThetaProdPlane_US_ME_weight_cuts->Integral());
  L0_L0bar_cosThetaProdPlane_US_ME_weight_cuts->Scale(1./L0_L0bar_cosThetaProdPlane_US_ME_weight_cuts->GetXaxis()->GetBinWidth(1));
  L0_L0bar_cosThetaProdPlane_US_ME_weight_cuts->Draw("p e same");


  L0_L0bar_cosThetaProdPlane_US_LS_cuts->Sumw2();
  L0_L0bar_cosThetaProdPlane_US_LS_cuts->SetMarkerSize(1.5);
  L0_L0bar_cosThetaProdPlane_US_LS_cuts->SetMarkerStyle(21);
  L0_L0bar_cosThetaProdPlane_US_LS_cuts->SetMarkerColor(kBlue);
  L0_L0bar_cosThetaProdPlane_US_LS_cuts->SetLineColor(kBlue);
  float nL0L0bar_US_LS_cuts = L0_L0bar_cosThetaProdPlane_US_LS_cuts->Integral();
  L0_L0bar_cosThetaProdPlane_US_LS_cuts->Scale(1./L0_L0bar_cosThetaProdPlane_US_LS_cuts->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
  L0_L0bar_cosThetaProdPlane_US_LS_cuts->Draw("p e same");

  L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts->Sumw2();
  L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts->SetMarkerSize(1.5);
  L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts->SetMarkerStyle(25);
  L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts->SetMarkerColor(kMagenta+1);
  L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts->SetLineColor(kMagenta+1);
  L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts->Scale(nL0L0bar_US_LS_cuts/L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts->Integral()); //scale ME
  L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts->Scale(1./L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts->GetXaxis()->GetBinWidth(1));
  L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts->Draw("p e same");

  TPaveText *L0_L0bar_text_no_corr_cuts = new TPaveText(0.6, 0.35, 0.8, 0.55, "NDC");
  L0_L0bar_text_no_corr_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  L0_L0bar_text_no_corr_cuts->AddText("Minimum bias");
  L0_L0bar_text_no_corr_cuts->AddText("#Lambda-#bar{#Lambda}");
  L0_L0bar_text_no_corr_cuts->AddText("No correction");
  L0_L0bar_text_no_corr_cuts->AddText("Analysis cuts");
  L0_L0bar_text_no_corr_cuts->SetFillColorAlpha(0, 0.01);
  L0_L0bar_text_no_corr_cuts->Draw("same");

  L0_L0bar_cosThetaProdPlane_US_no_corr_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/L0_L0bar_cosThetaProdPlane_US_no_corr_cuts.png");

  //---------------------------------------

  L0_L0bar_cosThetaProdPlane_US_separate_cuts_can->cd();

  L0_L0bar_cosThetaProdPlane_US_cuts->Divide(L0_L0bar_cosThetaProdPlane_US_ME_weight_cuts); //correct using ME
  L0_L0bar_cosThetaProdPlane_US_cuts->Scale(nL0L0bar_US_cuts/L0_L0bar_cosThetaProdPlane_US_cuts->Integral()); //scale back
  L0_L0bar_cosThetaProdPlane_US_cuts->Scale(1./L0_L0bar_cosThetaProdPlane_US_cuts->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
  L0_L0bar_cosThetaProdPlane_US_cuts->SetMinimum(0);
  L0_L0bar_cosThetaProdPlane_US_cuts->Draw("p e");


  L0_L0bar_cosThetaProdPlane_US_LS_cuts->Divide(L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts);
  L0_L0bar_cosThetaProdPlane_US_LS_cuts->Scale(nL0L0bar_US_LS_cuts/L0_L0bar_cosThetaProdPlane_US_LS_cuts->Integral()); //scale back
  L0_L0bar_cosThetaProdPlane_US_LS_cuts->Scale(1./L0_L0bar_cosThetaProdPlane_US_LS_cuts->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
  L0_L0bar_cosThetaProdPlane_US_LS_cuts->Draw("p e same");

  TPaveText *L0_L0bar_text_separate_cuts = new TPaveText(0.6, 0.35, 0.8, 0.55, "NDC");
  L0_L0bar_text_separate_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  L0_L0bar_text_separate_cuts->AddText("Minimum bias");
  L0_L0bar_text_separate_cuts->AddText("#Lambda-#bar{#Lambda}");
  L0_L0bar_text_separate_cuts->AddText("ME corrected");
  L0_L0bar_text_separate_cuts->AddText("Analysis cuts");
  L0_L0bar_text_separate_cuts->SetFillColorAlpha(0, 0.01);
  L0_L0bar_text_separate_cuts->Draw("same");


  L0_L0bar_cosThetaProdPlane_US_separate_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/L0_L0bar_cosThetaProdPlane_US_separate_cuts.png");

  //---------------------------------------

  L0_L0bar_cosThetaProdPlane_US_cuts_can->cd();

  //correct signal+background

  L0_L0bar_cosThetaProdPlane_US_cuts->Sumw2();
  L0_L0bar_cosThetaProdPlane_US_cuts->SetMarkerStyle(20);
  L0_L0bar_cosThetaProdPlane_US_cuts->SetMarkerColor(kRed);
  L0_L0bar_cosThetaProdPlane_US_cuts->SetLineColor(kRed);
  L0_L0bar_cosThetaProdPlane_US_cuts->GetXaxis()->SetTitle("cos(#theta*)");
  L0_L0bar_cosThetaProdPlane_US_cuts->GetXaxis()->CenterTitle();
  L0_L0bar_cosThetaProdPlane_US_cuts->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0_L0bar_cosThetaProdPlane_US_cuts->GetYaxis()->CenterTitle();
  L0_L0bar_cosThetaProdPlane_US_cuts->Add(L0_L0bar_cosThetaProdPlane_US_LS_cuts, -1); //subtract combinatorial background
  //L0_L0bar_cosThetaProdPlane_US_cuts->Scale(1./L0_L0bar_cosThetaProdPlane_US_cuts->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
  L0_L0bar_cosThetaProdPlane_US_cuts->SetMinimum(0);
  L0_L0bar_cosThetaProdPlane_US_cuts->Draw("p e");



  TPaveText *L0_L0bar_text_cuts = new TPaveText(0.6, 0.35, 0.8, 0.55, "NDC");
  L0_L0bar_text_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  L0_L0bar_text_cuts->AddText("Minimum bias");
  L0_L0bar_text_cuts->AddText("#Lambda-#bar{#Lambda}");
  L0_L0bar_text_cuts->AddText("US - Background");
  L0_L0bar_text_cuts->AddText("Analysis cuts");
  L0_L0bar_text_cuts->SetFillColorAlpha(0, 0.01);
  L0_L0bar_text_cuts->Draw("same");

  L0_L0bar_cosThetaProdPlane_US_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/L0_L0bar_cosThetaProdPlane_US_cuts.png");


  //_________________________________________________________________________

  //L-L
  L0_L0_cosThetaProdPlane_US_can->cd();

  //correct background first
  float nL0L0_US_LS = L0_L0_cosThetaProdPlane_US_LS->Integral();

  L0_L0_cosThetaProdPlane_US_LS_ME_weight->Sumw2();
  L0_L0_cosThetaProdPlane_US_LS_ME_weight->Scale(nL0L0_US_LS/L0_L0_cosThetaProdPlane_US_LS_ME_weight->Integral());

  L0_L0_cosThetaProdPlane_US_LS->Sumw2();
  L0_L0_cosThetaProdPlane_US_LS->Divide(L0_L0_cosThetaProdPlane_US_LS_ME_weight); //correct background using ME
  L0_L0_cosThetaProdPlane_US_LS->Scale(nL0L0_US_LS/L0_L0_cosThetaProdPlane_US_LS->Integral()); //correct background using ME
  //L0_L0_cosThetaProdPlane_US_LS->Scale(1./L0_L0_cosThetaProdPlane_US_LS->GetXaxis()->GetBinWidth(1)); //scale by bin width
  L0_L0_cosThetaProdPlane_US_LS->SetMinimum(0);

  //correct signal+background
  float nL0L0_US = L0_L0_cosThetaProdPlane_US->Integral();

  L0_L0_cosThetaProdPlane_US_ME_weight->Sumw2();
  L0_L0_cosThetaProdPlane_US_ME_weight->Scale(nL0L0_US/L0_L0_cosThetaProdPlane_US_ME_weight->Integral());

  L0_L0_cosThetaProdPlane_US->Sumw2();
  L0_L0_cosThetaProdPlane_US->SetMarkerStyle(20);
  L0_L0_cosThetaProdPlane_US->SetMarkerColor(kRed);
  L0_L0_cosThetaProdPlane_US->SetLineColor(kRed);
  L0_L0_cosThetaProdPlane_US->GetXaxis()->SetTitle("cos(#theta*)");
  L0_L0_cosThetaProdPlane_US->GetXaxis()->CenterTitle();
  L0_L0_cosThetaProdPlane_US->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0_L0_cosThetaProdPlane_US->GetYaxis()->CenterTitle();
  L0_L0_cosThetaProdPlane_US->Divide(L0_L0_cosThetaProdPlane_US_ME_weight); //correct using ME
  L0_L0_cosThetaProdPlane_US->Scale(nL0L0_US/L0_L0_cosThetaProdPlane_US->Integral()); //scale back to correct integral
  L0_L0_cosThetaProdPlane_US->Add(L0_L0_cosThetaProdPlane_US_LS, -1); //subtract combinatorial background
  L0_L0_cosThetaProdPlane_US->Scale(1./L0_L0_cosThetaProdPlane_US->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
  L0_L0_cosThetaProdPlane_US->SetMinimum(0);
  L0_L0_cosThetaProdPlane_US->Draw("p e");

  TPaveText *L0_L0_text = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
  L0_L0_text->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  L0_L0_text->AddText("Minimum bias");
  L0_L0_text->AddText("#Lambda-#Lambda");
  L0_L0_text->AddText("US - Background");
  L0_L0_text->AddText("No cuts");
  L0_L0_text->SetFillColorAlpha(0, 0.01);
  L0_L0_text->Draw("same");

  L0_L0_cosThetaProdPlane_US_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_L/L0_L0_cosThetaProdPlane_US.png");

  //-------------------------------------------------------------------------

  L0_L0_cosThetaProdPlane_US_cuts_can->cd();

  //correct background first
  float nL0L0_US_LS_cuts = L0_L0_cosThetaProdPlane_US_LS_cuts->Integral();

  L0_L0_cosThetaProdPlane_US_LS_ME_weight_cuts->Sumw2();
  L0_L0_cosThetaProdPlane_US_LS_ME_weight_cuts->Scale(nL0L0_US_LS_cuts/L0_L0_cosThetaProdPlane_US_LS_ME_weight_cuts->Integral());

  L0_L0_cosThetaProdPlane_US_LS_cuts->Sumw2();
  L0_L0_cosThetaProdPlane_US_LS_cuts->Divide(L0_L0_cosThetaProdPlane_US_LS_ME_weight_cuts); //correct background using ME
  L0_L0_cosThetaProdPlane_US_LS_cuts->Scale(nL0L0_US_LS_cuts/L0_L0_cosThetaProdPlane_US_LS_cuts->Integral()); //scale back
  //L0_L0_cosThetaProdPlane_US_LS_cuts->Scale(1./L0_L0_cosThetaProdPlane_US_LS_cuts->GetXaxis()->GetBinWidth(1)); //scale by bin width
  L0_L0_cosThetaProdPlane_US_LS_cuts->SetMinimum(0);

  //correct signal+background
  float nL0L0_US_cuts = L0_L0_cosThetaProdPlane_US_cuts->Integral();

  L0_L0_cosThetaProdPlane_US_ME_weight_cuts->Sumw2();
  L0_L0_cosThetaProdPlane_US_ME_weight_cuts->Scale(nL0L0_US_cuts/L0_L0_cosThetaProdPlane_US_ME_weight_cuts->Integral());

  L0_L0_cosThetaProdPlane_US_cuts->Sumw2();
  L0_L0_cosThetaProdPlane_US_cuts->SetMarkerStyle(20);
  L0_L0_cosThetaProdPlane_US_cuts->SetMarkerColor(kRed);
  L0_L0_cosThetaProdPlane_US_cuts->SetLineColor(kRed);
  L0_L0_cosThetaProdPlane_US_cuts->GetXaxis()->SetTitle("cos(#theta*)");
  L0_L0_cosThetaProdPlane_US_cuts->GetXaxis()->CenterTitle();
  L0_L0_cosThetaProdPlane_US_cuts->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0_L0_cosThetaProdPlane_US_cuts->GetYaxis()->CenterTitle();
  L0_L0_cosThetaProdPlane_US_cuts->Divide(L0_L0_cosThetaProdPlane_US_ME_weight_cuts); //correct using ME
  L0_L0_cosThetaProdPlane_US_cuts->Scale(nL0L0_US_cuts/L0_L0_cosThetaProdPlane_US_cuts->Integral()); //scale back
  L0_L0_cosThetaProdPlane_US_cuts->Add(L0_L0_cosThetaProdPlane_US_LS_cuts, -1); //subtract combinatorial background
  L0_L0_cosThetaProdPlane_US_cuts->Scale(1./L0_L0_cosThetaProdPlane_US_cuts->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
  L0_L0_cosThetaProdPlane_US_cuts->SetMinimum(0);
  L0_L0_cosThetaProdPlane_US_cuts->Draw("p e");

  //L0_L0_cosThetaProdPlane_US_ME_weight_cuts->Draw("p e same");

  L0_L0_cosThetaProdPlane_US_LS_cuts->SetMarkerStyle(20);
  L0_L0_cosThetaProdPlane_US_LS_cuts->SetMarkerColor(kBlue);
  L0_L0_cosThetaProdPlane_US_LS_cuts->SetLineColor(kBlue);
  //L0_L0_cosThetaProdPlane_US_LS_cuts->Draw("p e");
  //L0_L0_cosThetaProdPlane_US_LS_cuts->Draw("p e same");

  //L0_L0_cosThetaProdPlane_US_LS_ME_weight_cuts->Draw("p e same");

  TPaveText *L0_L0_text_cuts = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
  L0_L0_text_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  L0_L0_text_cuts->AddText("Minimum bias");
  L0_L0_text_cuts->AddText("#Lambda-#Lambda");
  L0_L0_text_cuts->AddText("US - Background");
  L0_L0_text_cuts->AddText("Analysis cuts");
  L0_L0_text_cuts->SetFillColorAlpha(0, 0.01);
  L0_L0_text_cuts->Draw("same");

  L0_L0_cosThetaProdPlane_US_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_L/L0_L0_cosThetaProdPlane_US_cuts.png");

  //_________________________________________________________________________

  //Lbar-Lbar
  L0bar_L0bar_cosThetaProdPlane_US_can->cd();

  //correct background first
  float nL0barL0bar_US_LS = L0bar_L0bar_cosThetaProdPlane_US_LS->Integral();

  L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight->Sumw2();
  L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight->Scale(nL0barL0bar_US_LS/L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight->Integral());

  L0bar_L0bar_cosThetaProdPlane_US_LS->Sumw2();
  L0bar_L0bar_cosThetaProdPlane_US_LS->Divide(L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight); //correct background using ME
  L0bar_L0bar_cosThetaProdPlane_US_LS->Scale(nL0barL0bar_US_LS/L0bar_L0bar_cosThetaProdPlane_US_LS->Integral()); //scale by bin width
  //L0bar_L0bar_cosThetaProdPlane_US_LS->Scale(1./L0bar_L0bar_cosThetaProdPlane_US_LS->GetXaxis()->GetBinWidth(1)); //scale by bin width
  L0bar_L0bar_cosThetaProdPlane_US_LS->SetMinimum(0);

  //correct signal+background
  float nL0barL0bar_US = L0bar_L0bar_cosThetaProdPlane_US->Integral();

  L0bar_L0bar_cosThetaProdPlane_US_ME_weight->Sumw2();
  L0bar_L0bar_cosThetaProdPlane_US_ME_weight->Scale(nL0barL0bar_US/L0bar_L0bar_cosThetaProdPlane_US_ME_weight->Integral());

  L0bar_L0bar_cosThetaProdPlane_US->Sumw2();
  L0bar_L0bar_cosThetaProdPlane_US->SetMarkerStyle(20);
  L0bar_L0bar_cosThetaProdPlane_US->SetMarkerColor(kRed);
  L0bar_L0bar_cosThetaProdPlane_US->SetLineColor(kRed);
  L0bar_L0bar_cosThetaProdPlane_US->GetXaxis()->SetTitle("cos(#theta*)");
  L0bar_L0bar_cosThetaProdPlane_US->GetXaxis()->CenterTitle();
  L0bar_L0bar_cosThetaProdPlane_US->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0bar_L0bar_cosThetaProdPlane_US->GetYaxis()->CenterTitle();
  L0bar_L0bar_cosThetaProdPlane_US->Divide(L0bar_L0bar_cosThetaProdPlane_US_ME_weight); //correct using ME
  L0bar_L0bar_cosThetaProdPlane_US->Scale(nL0barL0bar_US/L0bar_L0bar_cosThetaProdPlane_US->Integral()); //scale back
  L0bar_L0bar_cosThetaProdPlane_US->Add(L0bar_L0bar_cosThetaProdPlane_US_LS, -1); //subtract combinatorial background
  L0bar_L0bar_cosThetaProdPlane_US->Scale(1./L0bar_L0bar_cosThetaProdPlane_US->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
  L0bar_L0bar_cosThetaProdPlane_US->SetMinimum(0);
  L0bar_L0bar_cosThetaProdPlane_US->Draw("p e");


  TPaveText *L0bar_L0bar_text = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
  L0bar_L0bar_text->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  L0bar_L0bar_text->AddText("Minimum bias");
  L0bar_L0bar_text->AddText("#bar{#Lambda}-#bar{#Lambda}");
  L0bar_L0bar_text->AddText("US - Background");
  L0bar_L0bar_text->AddText("No cuts");
  L0bar_L0bar_text->SetFillColorAlpha(0, 0.01);
  L0bar_L0bar_text->Draw("same");

  L0bar_L0bar_cosThetaProdPlane_US_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/Lbar_Lbar/L0bar_L0bar_cosThetaProdPlane_US.png");

  //-------------------------------------------------------------------------

  L0bar_L0bar_cosThetaProdPlane_US_cuts_can->cd();

  //correct background first
  float nL0barL0bar_US_LS_cuts = L0bar_L0bar_cosThetaProdPlane_US_LS_cuts->Integral();

  L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts->Sumw2();
  L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts->Scale(nL0barL0bar_US_LS_cuts/L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts->Integral());

  L0bar_L0bar_cosThetaProdPlane_US_LS_cuts->Sumw2();
  L0bar_L0bar_cosThetaProdPlane_US_LS_cuts->Divide(L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts); //correct background using ME
  L0bar_L0bar_cosThetaProdPlane_US_LS_cuts->Scale(nL0barL0bar_US_LS_cuts/L0bar_L0bar_cosThetaProdPlane_US_LS_cuts->Integral()); //scale by bin width
  //L0bar_L0bar_cosThetaProdPlane_US_LS_cuts->Scale(1./L0bar_L0bar_cosThetaProdPlane_US_LS_cuts->GetXaxis()->GetBinWidth(1)); //scale by bin width
  L0bar_L0bar_cosThetaProdPlane_US_LS_cuts->SetMinimum(0);

  //correct signal+background
  float nL0barL0bar_US_cuts = L0bar_L0bar_cosThetaProdPlane_US_cuts->Integral();

  L0bar_L0bar_cosThetaProdPlane_US_ME_weight_cuts->Sumw2();
  L0bar_L0bar_cosThetaProdPlane_US_ME_weight_cuts->Scale(nL0barL0bar_US_cuts/L0bar_L0bar_cosThetaProdPlane_US_ME_weight_cuts->Integral());

  L0bar_L0bar_cosThetaProdPlane_US_cuts->Sumw2();
  L0bar_L0bar_cosThetaProdPlane_US_cuts->SetMarkerStyle(20);
  L0bar_L0bar_cosThetaProdPlane_US_cuts->SetMarkerColor(kRed);
  L0bar_L0bar_cosThetaProdPlane_US_cuts->SetLineColor(kRed);
  L0bar_L0bar_cosThetaProdPlane_US_cuts->GetXaxis()->SetTitle("cos(#theta*)");
  L0bar_L0bar_cosThetaProdPlane_US_cuts->GetXaxis()->CenterTitle();
  L0bar_L0bar_cosThetaProdPlane_US_cuts->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0bar_L0bar_cosThetaProdPlane_US_cuts->GetYaxis()->CenterTitle();
  L0bar_L0bar_cosThetaProdPlane_US_cuts->Divide(L0bar_L0bar_cosThetaProdPlane_US_ME_weight_cuts); //correct using ME
  L0bar_L0bar_cosThetaProdPlane_US_cuts->Scale(nL0barL0bar_US_cuts/L0bar_L0bar_cosThetaProdPlane_US_cuts->Integral()); //scale back
  L0bar_L0bar_cosThetaProdPlane_US_cuts->Add(L0bar_L0bar_cosThetaProdPlane_US_LS_cuts, -1); //subtract combinatorial background
  L0bar_L0bar_cosThetaProdPlane_US_cuts->Scale(1./L0bar_L0bar_cosThetaProdPlane_US_cuts->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
  L0bar_L0bar_cosThetaProdPlane_US_cuts->SetMinimum(0);
  //L0bar_L0bar_cosThetaProdPlane_US_cuts->GetYaxis()->SetRangeUser(0, 12000);
  //L0bar_L0bar_cosThetaProdPlane_US_cuts->GetYaxis()->SetRangeUser(-8000, 0);
  L0bar_L0bar_cosThetaProdPlane_US_cuts->Draw("p e");

  //L0bar_L0bar_cosThetaProdPlane_US_LS_cuts->Draw("p e same");

  TPaveText *L0bar_L0bar_text_cuts = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
  L0bar_L0bar_text_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  L0bar_L0bar_text_cuts->AddText("Minimum bias");
  L0bar_L0bar_text_cuts->AddText("#bar{#Lambda}-#bar{#Lambda}");
  L0bar_L0bar_text_cuts->AddText("US - Background");
  L0bar_L0bar_text_cuts->AddText("Analysis cuts");
  L0bar_L0bar_text_cuts->SetFillColorAlpha(0, 0.01);
  L0bar_L0bar_text_cuts->Draw("same");

  L0bar_L0bar_cosThetaProdPlane_US_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/Lbar_Lbar/L0bar_L0bar_cosThetaProdPlane_US_cuts.png");

  //_________________________________________________________________________________________________________________________

  float nLLbar_pT_US[nPtBins_corr][nPtBins_corr];
  float nLL_pT_US[nPtBins_corr][nPtBins_corr];
  float nLbarLbar_pT_US[nPtBins_corr][nPtBins_corr];

  float nLLbar_pT_US_LS[nPtBins_corr][nPtBins_corr];
  float nLL_pT_US_LS[nPtBins_corr][nPtBins_corr];
  float nLbarLbar_pT_US_LS[nPtBins_corr][nPtBins_corr];


  float nLLbar_pT_US_cuts[nPtBins_corr][nPtBins_corr];
  float nLL_pT_US_cuts[nPtBins_corr][nPtBins_corr];
  float nLbarLbar_pT_US_cuts[nPtBins_corr][nPtBins_corr];

  float nLLbar_pT_US_LS_cuts[nPtBins_corr][nPtBins_corr];
  float nLL_pT_US_LS_cuts[nPtBins_corr][nPtBins_corr];
  float nLbarLbar_pT_US_LS_cuts[nPtBins_corr][nPtBins_corr];

  //L-L correlation histograms in bins
  for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
  {
    for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
    {
      //L-Lbar
      L0_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2] = (TH1F*)inFile->Get(Form("L0_L0bar_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
      L0_L0bar_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2] = (TH1F*)inFile->Get(Form("L0_L0bar_cosThetaProdPlane_US_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));

      L0_L0bar_cosThetaProdPlane_US_ME_weight_pT_hist[pTbin1][pTbin2] = (TH1F*)inFile->Get(Form("L0_L0bar_cosThetaProdPlane_US_ME_weight_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
      L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_hist[pTbin1][pTbin2] = (TH1F*)inFile->Get(Form("L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));


      TCanvas *L0_L0bar_cosThetaProdPlane_US_pT_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_US_pT_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_US_pT_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);

      //correct background first
      nLLbar_pT_US_LS[pTbin1][pTbin2] = L0_L0bar_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2]->Integral();

      L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_hist[pTbin1][pTbin2]->Sumw2();
      L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_hist[pTbin1][pTbin2]->Scale(nLLbar_pT_US_LS[pTbin1][pTbin2]/L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_hist[pTbin1][pTbin2]->Integral());


      L0_L0bar_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2]->Sumw2();
      L0_L0bar_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2]->Divide(L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_hist[pTbin1][pTbin2]);
      L0_L0bar_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2]->Scale(nLLbar_pT_US_LS[pTbin1][pTbin2]/L0_L0bar_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2]->Integral()); //scale back
      //L0_L0bar_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2]->Scale(1./L0_L0bar_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2]->GetXaxis()->GetBinWidth(1));


      //correct signal+background
      nLLbar_pT_US[pTbin1][pTbin2] = L0_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->Integral();

      L0_L0bar_cosThetaProdPlane_US_ME_weight_pT_hist[pTbin1][pTbin2]->Sumw2();
      L0_L0bar_cosThetaProdPlane_US_ME_weight_pT_hist[pTbin1][pTbin2]->Scale(nLLbar_pT_US[pTbin1][pTbin2]/L0_L0bar_cosThetaProdPlane_US_ME_weight_pT_hist[pTbin1][pTbin2]->Integral());

      L0_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->Sumw2();
      L0_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->SetLineColor(kRed);
      L0_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      L0_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->Divide(L0_L0bar_cosThetaProdPlane_US_ME_weight_pT_hist[pTbin1][pTbin2]);
      L0_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->Scale(nLLbar_pT_US[pTbin1][pTbin2]/L0_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->Integral());
      L0_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->Add(L0_L0bar_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2], -1);
      L0_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->Scale(1./L0_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->GetXaxis()->GetBinWidth(1));
      L0_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->SetMinimum(0);
      L0_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->Draw("p e");


      TPaveText *L0_L0bar_text_pT = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
      L0_L0bar_text_pT->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
      L0_L0bar_text_pT->AddText("Minimum bias");
      L0_L0bar_text_pT->AddText("#Lambda-#bar{#Lambda}");
      L0_L0bar_text_pT->AddText("US - Background");
      L0_L0bar_text_pT->AddText("No cuts");
      L0_L0bar_text_pT->AddText("|#eta| < 1");
      L0_L0bar_text_pT->AddText(Form("%0.1f < p_{T}^{1} < %0.1f", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0_L0bar_text_pT->AddText(Form("%0.1f < p_{T}^{2} < %0.1f", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0_L0bar_text_pT->SetFillColorAlpha(0, 0.01);
      L0_L0bar_text_pT->Draw("same");

      L0_L0bar_cosThetaProdPlane_US_pT_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/L0_L0bar_cosThetaProdPlane_US_pT1_%i_pT2_%i.png", pTbin1, pTbin2));

      //--------------------------------------------------------------------------------

      L0_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2] = (TH1F*)inFile->Get(Form("L0_L0bar_cosThetaProdPlane_US_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
      L0_L0bar_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2] = (TH1F*)inFile->Get(Form("L0_L0bar_cosThetaProdPlane_US_LS_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));

      L0_L0bar_cosThetaProdPlane_US_ME_weight_pT_cuts_hist[pTbin1][pTbin2] = (TH1F*)inFile->Get(Form("L0_L0bar_cosThetaProdPlane_US_ME_weight_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
      L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_cuts_hist[pTbin1][pTbin2] = (TH1F*)inFile->Get(Form("L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));


      TCanvas *L0_L0bar_cosThetaProdPlane_US_pT_cuts_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_US_cuts_pT_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_US_cuts_pT_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);

      //correct background first
      nLLbar_pT_US_LS_cuts[pTbin1][pTbin2] = L0_L0bar_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2]->Integral();

      L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_cuts_hist[pTbin1][pTbin2]->Sumw2();
      L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_cuts_hist[pTbin1][pTbin2]->Scale(nLLbar_pT_US_LS_cuts[pTbin1][pTbin2]/L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_cuts_hist[pTbin1][pTbin2]->Integral());


      L0_L0bar_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2]->Sumw2();
      L0_L0bar_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2]->Divide(L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_cuts_hist[pTbin1][pTbin2]);
      L0_L0bar_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2]->Scale(nLLbar_pT_US_LS_cuts[pTbin1][pTbin2]/L0_L0bar_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2]->Integral());
      //L0_L0bar_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2]->Scale(1./L0_L0bar_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2]->GetXaxis()->GetBinWidth(1));


      //correct signal+background
      nLLbar_pT_US_cuts[pTbin1][pTbin2] = L0_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->Integral();

      L0_L0bar_cosThetaProdPlane_US_ME_weight_pT_cuts_hist[pTbin1][pTbin2]->Sumw2();
      L0_L0bar_cosThetaProdPlane_US_ME_weight_pT_cuts_hist[pTbin1][pTbin2]->Scale(nLLbar_pT_US_cuts[pTbin1][pTbin2]/L0_L0bar_cosThetaProdPlane_US_ME_weight_pT_cuts_hist[pTbin1][pTbin2]->Integral());

      L0_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->Sumw2();
      L0_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->SetLineColor(kRed);
      L0_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      L0_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->Divide(L0_L0bar_cosThetaProdPlane_US_ME_weight_pT_cuts_hist[pTbin1][pTbin2]);
      L0_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->Scale(nLLbar_pT_US_cuts[pTbin1][pTbin2]/L0_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->Integral());
      L0_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->Add(L0_L0bar_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2], -1);
      L0_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->Scale(1./L0_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->GetXaxis()->GetBinWidth(1));
      L0_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->SetMinimum(0);
      L0_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->Draw("p e");


      TPaveText *L0_L0bar_text_pT_cuts = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
      L0_L0bar_text_pT_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
      L0_L0bar_text_pT_cuts->AddText("Minimum bias");
      L0_L0bar_text_pT_cuts->AddText("#Lambda-#bar{#Lambda}");
      L0_L0bar_text_pT_cuts->AddText("US - Background");
      L0_L0bar_text_pT_cuts->AddText("Analysis cuts");
      L0_L0bar_text_pT_cuts->AddText("|#eta| < 1");
      L0_L0bar_text_pT_cuts->AddText(Form("%0.1f < p_{T}^{1} < %0.1f", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0_L0bar_text_pT_cuts->AddText(Form("%0.1f < p_{T}^{2} < %0.1f", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0_L0bar_text_pT_cuts->SetFillColorAlpha(0, 0.01);
      L0_L0bar_text_pT_cuts->Draw("same");

      L0_L0bar_cosThetaProdPlane_US_pT_cuts_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/L0_L0bar_cosThetaProdPlane_US_cuts_pT1_%i_pT2_%i.png", pTbin1, pTbin2));

      //_________________________________________________________________________________________________________________________________________________________________________________________________________________________

      //L-L
      L0_L0_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2] = (TH1F*)inFile->Get(Form("L0_L0_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
      L0_L0_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2] = (TH1F*)inFile->Get(Form("L0_L0_cosThetaProdPlane_US_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));

      L0_L0_cosThetaProdPlane_US_ME_weight_pT_hist[pTbin1][pTbin2] = (TH1F*)inFile->Get(Form("L0_L0_cosThetaProdPlane_US_ME_weight_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
      L0_L0_cosThetaProdPlane_US_LS_ME_weight_pT_hist[pTbin1][pTbin2] = (TH1F*)inFile->Get(Form("L0_L0_cosThetaProdPlane_US_LS_ME_weight_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));



      TCanvas *L0_L0_cosThetaProdPlane_US_pT_can = new TCanvas(Form("L0_L0_cosThetaProdPlane_US_pT_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_US_pT_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);

      //correct background first
      nLL_pT_US_LS[pTbin1][pTbin2] = L0_L0_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2]->Integral();

      L0_L0_cosThetaProdPlane_US_LS_ME_weight_pT_hist[pTbin1][pTbin2]->Sumw2();
      L0_L0_cosThetaProdPlane_US_LS_ME_weight_pT_hist[pTbin1][pTbin2]->Scale(nLL_pT_US_LS[pTbin1][pTbin2]/L0_L0_cosThetaProdPlane_US_LS_ME_weight_pT_hist[pTbin1][pTbin2]->Integral());


      L0_L0_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2]->Sumw2();
      L0_L0_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2]->Divide(L0_L0_cosThetaProdPlane_US_LS_ME_weight_pT_hist[pTbin1][pTbin2]);
      L0_L0_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2]->Scale(nLL_pT_US_LS[pTbin1][pTbin2]/L0_L0_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2]->Integral());
      //L0_L0_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2]->Scale(1./L0_L0_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2]->GetXaxis()->GetBinWidth(1));


      //correct signal+background
      nLL_pT_US[pTbin1][pTbin2] = L0_L0_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->Integral();

      L0_L0_cosThetaProdPlane_US_ME_weight_pT_hist[pTbin1][pTbin2]->Sumw2();
      L0_L0_cosThetaProdPlane_US_ME_weight_pT_hist[pTbin1][pTbin2]->Scale(nLL_pT_US[pTbin1][pTbin2]/L0_L0_cosThetaProdPlane_US_ME_weight_pT_hist[pTbin1][pTbin2]->Integral());

      L0_L0_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->Sumw2();
      L0_L0_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0_L0_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0_L0_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->SetLineColor(kRed);
      L0_L0_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0_L0_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_L0_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      L0_L0_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_L0_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->Divide(L0_L0_cosThetaProdPlane_US_ME_weight_pT_hist[pTbin1][pTbin2]);
      L0_L0_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->Scale(nLL_pT_US[pTbin1][pTbin2]/L0_L0_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->Integral());
      L0_L0_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->Add(L0_L0_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2], -1);
      L0_L0_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->Scale(1./L0_L0_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->GetXaxis()->GetBinWidth(1));
      L0_L0_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->SetMinimum(0);
      L0_L0_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->Draw("p e");

      TPaveText *L0_L0_text_pT = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
      L0_L0_text_pT->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
      L0_L0_text_pT->AddText("Minimum bias");
      L0_L0_text_pT->AddText("#Lambda-#Lambda");
      L0_L0_text_pT->AddText("US - Background");
      L0_L0_text_pT->AddText("No cuts");
      L0_L0_text_pT->AddText("|#eta| < 1");
      L0_L0_text_pT->AddText(Form("%0.1f < p_{T}^{1} < %0.1f", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0_L0_text_pT->AddText(Form("%0.1f < p_{T}^{2} < %0.1f", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0_L0_text_pT->SetFillColorAlpha(0, 0.01);
      L0_L0_text_pT->Draw("same");

      L0_L0_cosThetaProdPlane_US_pT_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_L/L0_L0_cosThetaProdPlane_US_pT1_%i_pT2_%i.png", pTbin1, pTbin2));

      //--------------------------------------------------------------------------------

      L0_L0_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2] = (TH1F*)inFile->Get(Form("L0_L0_cosThetaProdPlane_US_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
      L0_L0_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2] = (TH1F*)inFile->Get(Form("L0_L0_cosThetaProdPlane_US_LS_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));

      L0_L0_cosThetaProdPlane_US_ME_weight_pT_cuts_hist[pTbin1][pTbin2] = (TH1F*)inFile->Get(Form("L0_L0_cosThetaProdPlane_US_ME_weight_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
      L0_L0_cosThetaProdPlane_US_LS_ME_weight_pT_cuts_hist[pTbin1][pTbin2] = (TH1F*)inFile->Get(Form("L0_L0_cosThetaProdPlane_US_LS_ME_weight_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));

      TCanvas *L0_L0_cosThetaProdPlane_US_pT_cuts_can = new TCanvas(Form("L0_L0_cosThetaProdPlane_US_pT_cuts_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_US_pT_cuts_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);

      //correct background first
      nLL_pT_US_LS_cuts[pTbin1][pTbin2] = L0_L0_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2]->Integral();

      L0_L0_cosThetaProdPlane_US_LS_ME_weight_pT_cuts_hist[pTbin1][pTbin2]->Sumw2();
      L0_L0_cosThetaProdPlane_US_LS_ME_weight_pT_cuts_hist[pTbin1][pTbin2]->Scale(nLL_pT_US_LS_cuts[pTbin1][pTbin2]/L0_L0_cosThetaProdPlane_US_LS_ME_weight_pT_cuts_hist[pTbin1][pTbin2]->Integral());


      L0_L0_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2]->Sumw2();
      L0_L0_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2]->Divide(L0_L0_cosThetaProdPlane_US_LS_ME_weight_pT_cuts_hist[pTbin1][pTbin2]);
      L0_L0_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2]->Scale(nLL_pT_US_LS_cuts[pTbin1][pTbin2]/L0_L0_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2]->GetXaxis()->GetBinWidth(1));
      //L0_L0_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2]->Scale(1./L0_L0_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2]->GetXaxis()->GetBinWidth(1));

      //correct signal+background
      nLL_pT_US_cuts[pTbin1][pTbin2] = L0_L0_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->Integral();

      L0_L0_cosThetaProdPlane_US_ME_weight_pT_cuts_hist[pTbin1][pTbin2]->Sumw2();
      L0_L0_cosThetaProdPlane_US_ME_weight_pT_cuts_hist[pTbin1][pTbin2]->Scale(nLL_pT_US_cuts[pTbin1][pTbin2]/L0_L0_cosThetaProdPlane_US_ME_weight_pT_cuts_hist[pTbin1][pTbin2]->Integral());

      L0_L0_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->Sumw2();
      L0_L0_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0_L0_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0_L0_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->SetLineColor(kRed);
      L0_L0_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0_L0_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_L0_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      L0_L0_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_L0_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->Divide(L0_L0_cosThetaProdPlane_US_ME_weight_pT_cuts_hist[pTbin1][pTbin2]);
      L0_L0_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->Scale(nLL_pT_US_cuts[pTbin1][pTbin2]/L0_L0_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->Integral());
      L0_L0_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->Add(L0_L0_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2], -1);
      L0_L0_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->Scale(1./L0_L0_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->GetXaxis()->GetBinWidth(1));
      L0_L0_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->SetMinimum(0);
      L0_L0_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->Draw("p e");

      TPaveText *L0_L0_text_pT_cuts = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
      L0_L0_text_pT_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
      L0_L0_text_pT_cuts->AddText("Minimum bias");
      L0_L0_text_pT_cuts->AddText("#Lambda-#Lambda");
      L0_L0_text_pT_cuts->AddText("US - Background");
      L0_L0_text_pT_cuts->AddText("Analysis cuts");
      L0_L0_text_pT_cuts->AddText("|#eta| < 1");
      L0_L0_text_pT_cuts->AddText(Form("%0.1f < p_{T}^{1} < %0.1f", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0_L0_text_pT_cuts->AddText(Form("%0.1f < p_{T}^{2} < %0.1f", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0_L0_text_pT_cuts->SetFillColorAlpha(0, 0.01);
      L0_L0_text_pT_cuts->Draw("same");

      L0_L0_cosThetaProdPlane_US_pT_cuts_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_L/L0_L0_cosThetaProdPlane_US_cuts_pT1_%i_pT2_%i.png", pTbin1, pTbin2));


      //_________________________________________________________________________________________________

      //Lbar-Lbar
      L0bar_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2] = (TH1F*)inFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
      L0bar_L0bar_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2] = (TH1F*)inFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_US_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));

      L0bar_L0bar_cosThetaProdPlane_US_ME_weight_pT_hist[pTbin1][pTbin2] = (TH1F*)inFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_US_ME_weight_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
      L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_hist[pTbin1][pTbin2] = (TH1F*)inFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));


      TCanvas *L0bar_L0bar_cosThetaProdPlane_US_pT_can = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_US_pT_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_US_pT_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);

      //correct background first
      nLbarLbar_pT_US_LS[pTbin1][pTbin2] = L0bar_L0bar_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2]->Integral();

      L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_hist[pTbin1][pTbin2]->Sumw2();
      L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_hist[pTbin1][pTbin2]->Scale(nLbarLbar_pT_US_LS[pTbin1][pTbin2]/L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_hist[pTbin1][pTbin2]->Integral());


      L0bar_L0bar_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2]->Sumw2();
      L0bar_L0bar_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2]->Divide(L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_hist[pTbin1][pTbin2]);
      L0bar_L0bar_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2]->Scale(nLbarLbar_pT_US_LS[pTbin1][pTbin2]/L0bar_L0bar_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2]->Integral());
      //L0bar_L0bar_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2]->Scale(1./L0bar_L0bar_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2]->GetXaxis()->GetBinWidth(1));


      //correct signal+background
      nLbarLbar_pT_US[pTbin1][pTbin2] = L0bar_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->Integral();

      L0bar_L0bar_cosThetaProdPlane_US_ME_weight_pT_hist[pTbin1][pTbin2]->Sumw2();
      L0bar_L0bar_cosThetaProdPlane_US_ME_weight_pT_hist[pTbin1][pTbin2]->Scale(nLbarLbar_pT_US[pTbin1][pTbin2]/L0bar_L0bar_cosThetaProdPlane_US_ME_weight_pT_hist[pTbin1][pTbin2]->Integral());

      L0bar_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->Sumw2();
      L0bar_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0bar_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0bar_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->SetLineColor(kRed);
      L0bar_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0bar_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0bar_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      L0bar_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0bar_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->Divide(L0bar_L0bar_cosThetaProdPlane_US_ME_weight_pT_hist[pTbin1][pTbin2]);
      L0bar_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->Scale(nLbarLbar_pT_US[pTbin1][pTbin2]/L0bar_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->Integral());
      L0bar_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->Add(L0bar_L0bar_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2], -1);
      L0bar_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->Scale(1./L0bar_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->GetXaxis()->GetBinWidth(1));
      L0bar_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->SetMinimum(0);
      L0bar_L0bar_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->Draw("p e");

      TPaveText *L0bar_L0bar_text_pT = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
      L0bar_L0bar_text_pT->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
      L0bar_L0bar_text_pT->AddText("Minimum bias");
      L0bar_L0bar_text_pT->AddText("#bar{#Lambda}-#bar{#Lambda}");
      L0bar_L0bar_text_pT->AddText("US - Background");
      L0bar_L0bar_text_pT->AddText("No cuts");
      L0bar_L0bar_text_pT->AddText("|#eta| < 1");
      L0bar_L0bar_text_pT->AddText(Form("%0.1f < p_{T}^{1} < %0.1f", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0bar_L0bar_text_pT->AddText(Form("%0.1f < p_{T}^{2} < %0.1f", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0bar_L0bar_text_pT->SetFillColorAlpha(0, 0.01);
      L0bar_L0bar_text_pT->Draw("same");



      L0bar_L0bar_cosThetaProdPlane_US_pT_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/Lbar_Lbar/L0bar_L0bar_cosThetaProdPlane_US_pT1_%i_pT2_%i.png", pTbin1, pTbin2));

      //--------------------------------------------------------------------------------

      L0bar_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2] = (TH1F*)inFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_US_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
      L0bar_L0bar_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2] = (TH1F*)inFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_US_LS_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));

      L0bar_L0bar_cosThetaProdPlane_US_ME_weight_pT_cuts_hist[pTbin1][pTbin2] = (TH1F*)inFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_US_ME_weight_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
      L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_cuts_hist[pTbin1][pTbin2] = (TH1F*)inFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));

      TCanvas *L0bar_L0bar_cosThetaProdPlane_US_pT_cuts_can = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_US_pT_cuts_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_US_pT_cuts_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);

      //correct background first
      nLbarLbar_pT_US_LS_cuts[pTbin1][pTbin2] = L0bar_L0bar_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2]->Integral();

      L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_cuts_hist[pTbin1][pTbin2]->Sumw2();
      L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_cuts_hist[pTbin1][pTbin2]->Scale(nLbarLbar_pT_US_LS_cuts[pTbin1][pTbin2]/L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_cuts_hist[pTbin1][pTbin2]->Integral());


      L0bar_L0bar_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2]->Sumw2();
      L0bar_L0bar_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2]->Divide(L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_pT_cuts_hist[pTbin1][pTbin2]);
      L0bar_L0bar_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2]->Scale(nLbarLbar_pT_US_LS_cuts[pTbin1][pTbin2]/L0bar_L0bar_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2]->Integral());
      //L0bar_L0bar_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2]->Scale(1./L0bar_L0bar_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2]->GetXaxis()->GetBinWidth(1));

      //correct signal+background
      nLbarLbar_pT_US_cuts[pTbin1][pTbin2] = L0bar_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->Integral();

      L0bar_L0bar_cosThetaProdPlane_US_ME_weight_pT_cuts_hist[pTbin1][pTbin2]->Sumw2();
      L0bar_L0bar_cosThetaProdPlane_US_ME_weight_pT_cuts_hist[pTbin1][pTbin2]->Scale(nLbarLbar_pT_US_cuts[pTbin1][pTbin2]/L0bar_L0bar_cosThetaProdPlane_US_ME_weight_pT_cuts_hist[pTbin1][pTbin2]->Integral());

      L0bar_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->Sumw2();
      L0bar_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0bar_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0bar_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->SetLineColor(kRed);
      L0bar_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0bar_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0bar_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      L0bar_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0bar_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->Divide(L0bar_L0bar_cosThetaProdPlane_US_ME_weight_pT_cuts_hist[pTbin1][pTbin2]);
      L0bar_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->Scale(nLbarLbar_pT_US_cuts[pTbin1][pTbin2]/L0bar_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->Integral());
      L0bar_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->Add(L0bar_L0bar_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2], -1);
      L0bar_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->Scale(1./L0bar_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->GetXaxis()->GetBinWidth(1));
      L0bar_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->SetMinimum(0);
      L0bar_L0bar_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->Draw("p e");

      TPaveText *L0bar_L0bar_text_pT_cuts = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
      L0bar_L0bar_text_pT_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
      L0bar_L0bar_text_pT_cuts->AddText("Minimum bias");
      L0bar_L0bar_text_pT_cuts->AddText("#bar{#Lambda}-#bar{#Lambda}");
      L0bar_L0bar_text_pT_cuts->AddText("US - Background");
      L0bar_L0bar_text_pT_cuts->AddText("Analysis cuts");
      L0bar_L0bar_text_pT_cuts->AddText("|#eta| < 1");
      L0bar_L0bar_text_pT_cuts->AddText(Form("%0.1f < p_{T}^{1} < %0.1f", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0bar_L0bar_text_pT_cuts->AddText(Form("%0.1f < p_{T}^{2} < %0.1f", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0bar_L0bar_text_pT_cuts->SetFillColorAlpha(0, 0.01);
      L0bar_L0bar_text_pT_cuts->Draw("same");

      L0bar_L0bar_cosThetaProdPlane_US_pT_cuts_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/Lbar_Lbar/L0bar_L0bar_cosThetaProdPlane_US_cuts_pT1_%i_pT2_%i.png", pTbin1, pTbin2));

      //_________________________________________________________________________________________________________________________________________________________________________________________________________________________________


      //invariant mass histograms after cuts
      L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2] = (TH2F*)inFile->Get(Form("L0_inv_mass_vs_L0bar_inv_mass_US_pT1_%i_pT2_%i", pTbin1, pTbin2)); //for US-US Lambda pairs
      L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2] = (TH2F*)inFile->Get(Form("L0_inv_mass_vs_L0bar_inv_mass_US_LS_pT1_%i_pT2_%i", pTbin1, pTbin2)); //for US-LS Lambda pairs

      TCanvas *L0_inv_mass_vs_L0bar_inv_mass_US_can = new TCanvas(Form("L0_inv_mass_vs_L0bar_inv_mass_US_can_%i_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0bar_inv_mass_US_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0_inv_mass_vs_L0bar_inv_mass_US_can->cd();

      L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(2);
      L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->SetTitleOffset(2.1);
      L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      //L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetZaxis()->SetRangeUser(0,1000);
      L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->Draw("surf1");

      TPaveText *L0_L0bar_text_Minv_2D_US = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
      L0_L0bar_text_Minv_2D_US->SetTextFont(42);
      //L0_L0bar_text_no_corr->AddText("STAR Internal");
      //L0_L0bar_text_no_corr->AddText("STAR preliminary");
      //((TText*)L0_L0bar_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
      //L0_L0bar_text_Minv_2D_US->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0_L0bar_text_Minv_2D_US->AddText("Minimum bias");
      L0_L0bar_text_Minv_2D_US->AddText("#Lambda^{0}-#bar{#Lambda^{0}}");
      L0_L0bar_text_Minv_2D_US->SetFillColorAlpha(0, 0.01);
      L0_L0bar_text_Minv_2D_US->Draw("same");


      TPaveText *L0_L0bar_text_Minv_2D_kine = new TPaveText(0.75, 0.83, 0.95, 0.98, "NDC");
      L0_L0bar_text_Minv_2D_kine->SetTextFont(42);
      L0_L0bar_text_Minv_2D_kine->AddText("|#it{y}| < 1");
      L0_L0bar_text_Minv_2D_kine->AddText(Form("%.2f < p_{T}^{1} < %.2f", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0_L0bar_text_Minv_2D_kine->AddText(Form("%.2f < p_{T}^{2} < %.2f", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0_L0bar_text_Minv_2D_kine->SetFillColorAlpha(0, 0.01);
      L0_L0bar_text_Minv_2D_kine->Draw("same");

      L0_inv_mass_vs_L0bar_inv_mass_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/Minv/L0_inv_mass_vs_L0bar_inv_mass_US_pT1_%i_pT2_%i.png", pTbin1, pTbin2));



      TCanvas *L0_inv_mass_vs_L0bar_inv_mass_US_LS_can = new TCanvas(Form("L0_inv_mass_vs_L0bar_inv_mass_US_LS_can_%i_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0bar_inv_mass_US_LS_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0_inv_mass_vs_L0bar_inv_mass_US_LS_can->cd();

      L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->Sumw2();
      L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(2);
      L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetYaxis()->SetTitleOffset(2);
      L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->Draw("surf1");

      TPaveText *L0_L0bar_text_Minv_2D_US_LS = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
      L0_L0bar_text_Minv_2D_US_LS->SetTextFont(42);
      //L0_L0bar_text_no_corr->AddText("STAR Internal");
      //L0_L0bar_text_no_corr->AddText("STAR preliminary");
      //((TText*)L0_L0bar_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
      //L0_L0bar_text_Minv_2D_US_LS->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0_L0bar_text_Minv_2D_US_LS->AddText("Minimum bias");
      L0_L0bar_text_Minv_2D_US_LS->AddText("#Lambda^{0}-#bar{#Lambda^{0}}");
      L0_L0bar_text_Minv_2D_US_LS->SetFillColorAlpha(0, 0.01);
      L0_L0bar_text_Minv_2D_US_LS->Draw("same");

      L0_L0bar_text_Minv_2D_kine->Draw("same");

      L0_inv_mass_vs_L0bar_inv_mass_US_LS_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/Minv/L0_inv_mass_vs_L0bar_inv_mass_US_LS_pT1_%i_pT2_%i.png", pTbin1, pTbin2));

      //-------------------------------------------------------------------------

      L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2] = (TH2F*)inFile->Get(Form("L0_inv_mass_vs_L0_inv_mass_US_pT1_%i_pT2_%i", pTbin1, pTbin2)); //for US-US Lambda pairs
      L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2] = (TH2F*)inFile->Get(Form("L0_inv_mass_vs_L0_inv_mass_US_LS_pT1_%i_pT2_%i", pTbin1, pTbin2)); //for US-LS Lambda pairs

      TCanvas *L0_inv_mass_vs_L0_inv_mass_US_can = new TCanvas(Form("L0_inv_mass_vs_L0_inv_mass_US_can_%i_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0_inv_mass_US_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0_inv_mass_vs_L0_inv_mass_US_can->cd();

      L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(2);
      L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->SetTitleOffset(2.1);
      L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      //L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->GetZaxis()->SetRangeUser(0,1000);
      L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->Draw("surf1");

      TPaveText *L0_L0_text_Minv_2D_US = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
      L0_L0_text_Minv_2D_US->SetTextFont(42);
      //L0_L0_text_no_corr->AddText("STAR Internal");
      //L0_L0_text_no_corr->AddText("STAR preliminary");
      //((TText*)L0_L0_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
      //L0_L0_text_Minv_2D_US->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0_L0_text_Minv_2D_US->AddText("Minimum bias");
      L0_L0_text_Minv_2D_US->AddText("#Lambda^{0}-#Lambda^{0}");
      L0_L0_text_Minv_2D_US->SetFillColorAlpha(0, 0.01);
      L0_L0_text_Minv_2D_US->Draw("same");


      TPaveText *L0_L0_text_Minv_2D_kine = new TPaveText(0.75, 0.83, 0.95, 0.98, "NDC");
      L0_L0_text_Minv_2D_kine->SetTextFont(42);
      L0_L0_text_Minv_2D_kine->AddText("|#it{y}| < 1");
      L0_L0_text_Minv_2D_kine->AddText(Form("%.2f < p_{T}^{1} < %.2f", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0_L0_text_Minv_2D_kine->AddText(Form("%.2f < p_{T}^{2} < %.2f", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0_L0_text_Minv_2D_kine->SetFillColorAlpha(0, 0.01);
      L0_L0_text_Minv_2D_kine->Draw("same");

      L0_inv_mass_vs_L0_inv_mass_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_L/Minv/L0_inv_mass_vs_L0_inv_mass_US_pT1_%i_pT2_%i.png", pTbin1, pTbin2));



      TCanvas *L0_inv_mass_vs_L0_inv_mass_US_LS_can = new TCanvas(Form("L0_inv_mass_vs_L0_inv_mass_US_LS_can_%i_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0_inv_mass_US_LS_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0_inv_mass_vs_L0_inv_mass_US_LS_can->cd();

      L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2]->Sumw2();
      L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(2);
      L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2]->GetYaxis()->SetTitleOffset(2);
      L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2]->Draw("surf1");

      TPaveText *L0_L0_text_Minv_2D_US_LS = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
      L0_L0_text_Minv_2D_US_LS->SetTextFont(42);
      //L0_L0_text_no_corr->AddText("STAR Internal");
      //L0_L0_text_no_corr->AddText("STAR preliminary");
      //((TText*)L0_L0_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
      //L0_L0_text_Minv_2D_US_LS->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0_L0_text_Minv_2D_US_LS->AddText("Minimum bias");
      L0_L0_text_Minv_2D_US_LS->AddText("#Lambda^{0}-#Lambda^{0}");
      L0_L0_text_Minv_2D_US_LS->SetFillColorAlpha(0, 0.01);
      L0_L0_text_Minv_2D_US_LS->Draw("same");

      L0_L0_text_Minv_2D_kine->Draw("same");

      L0_inv_mass_vs_L0_inv_mass_US_LS_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_L/Minv/L0_inv_mass_vs_L0_inv_mass_US_LS_pT1_%i_pT2_%i.png", pTbin1, pTbin2));

      //-------------------------------------------------------------------------

      L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2] = (TH2F*)inFile->Get(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_pT1_%i_pT2_%i", pTbin1, pTbin2)); //for US-US Lambda pairs
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2] = (TH2F*)inFile->Get(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_pT1_%i_pT2_%i", pTbin1, pTbin2)); //for US-LS Lambda pairs

      TCanvas *L0bar_inv_mass_vs_L0bar_inv_mass_US_can = new TCanvas(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_can_%i_%i", pTbin1, pTbin2), Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_can->cd();

      L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(2);
      L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->SetTitleOffset(2.1);
      L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      //L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetZaxis()->SetRangeUser(0,1000);
      L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->Draw("surf1");

      TPaveText *L0bar_L0bar_text_Minv_2D_US = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
      L0bar_L0bar_text_Minv_2D_US->SetTextFont(42);
      //L0bar_L0bar_text_no_corr->AddText("STAR Internal");
      //L0bar_L0bar_text_no_corr->AddText("STAR preliminary");
      //((TText*)L0bar_L0bar_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
      //L0bar_L0bar_text_Minv_2D_US->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0bar_L0bar_text_Minv_2D_US->AddText("Minimum bias");
      L0bar_L0bar_text_Minv_2D_US->AddText("#bar{#Lambda^{0}}-#bar{#Lambda^{0}}");
      L0bar_L0bar_text_Minv_2D_US->SetFillColorAlpha(0, 0.01);
      L0bar_L0bar_text_Minv_2D_US->Draw("same");


      TPaveText *L0bar_L0bar_text_Minv_2D_kine = new TPaveText(0.75, 0.83, 0.95, 0.98, "NDC");
      L0bar_L0bar_text_Minv_2D_kine->SetTextFont(42);
      L0bar_L0bar_text_Minv_2D_kine->AddText("|#it{y}| < 1");
      L0bar_L0bar_text_Minv_2D_kine->AddText(Form("%.2f < p_{T}^{1} < %.2f", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0bar_L0bar_text_Minv_2D_kine->AddText(Form("%.2f < p_{T}^{2} < %.2f", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0bar_L0bar_text_Minv_2D_kine->SetFillColorAlpha(0, 0.01);
      L0bar_L0bar_text_Minv_2D_kine->Draw("same");

      L0bar_inv_mass_vs_L0bar_inv_mass_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/Lbar_Lbar/Minv/L0bar_inv_mass_vs_L0bar_inv_mass_US_pT1_%i_pT2_%i.png", pTbin1, pTbin2));



      TCanvas *L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_can = new TCanvas(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_can_%i_%i", pTbin1, pTbin2), Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_can->cd();

      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->Sumw2();
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(2);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetYaxis()->SetTitleOffset(2);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->Draw("surf1");

      TPaveText *L0bar_L0bar_text_Minv_2D_US_LS = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
      L0bar_L0bar_text_Minv_2D_US_LS->SetTextFont(42);
      //L0bar_L0bar_text_no_corr->AddText("STAR Internal");
      //L0bar_L0bar_text_no_corr->AddText("STAR preliminary");
      //((TText*)L0bar_L0bar_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
      //L0bar_L0bar_text_Minv_2D_US_LS->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0bar_L0bar_text_Minv_2D_US_LS->AddText("Minimum bias");
      L0bar_L0bar_text_Minv_2D_US_LS->AddText("#bar{#Lambda^{0}}-#bar{#Lambda^{0}}");
      L0bar_L0bar_text_Minv_2D_US_LS->SetFillColorAlpha(0, 0.01);
      L0bar_L0bar_text_Minv_2D_US_LS->Draw("same");

      L0bar_L0bar_text_Minv_2D_kine->Draw("same");

      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/Lbar_Lbar/Minv/L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_pT1_%i_pT2_%i.png", pTbin1, pTbin2));

      //________________________________________________________________________________________________________________________________________________

      //invariant mass histograms after cuts
      L0_inv_mass_vs_L0bar_inv_mass_US_cuts[pTbin1][pTbin2] = (TH2F*)inFile->Get(Form("L0_inv_mass_vs_L0bar_inv_mass_US_cuts_pT1_%i_pT2_%i", pTbin1, pTbin2)); //for US-US Lambda pairs
      L0_inv_mass_vs_L0bar_inv_mass_US_LS_cuts[pTbin1][pTbin2] = (TH2F*)inFile->Get(Form("L0_inv_mass_vs_L0bar_inv_mass_US_LS_cuts_pT1_%i_pT2_%i", pTbin1, pTbin2)); //for US-LS Lambda pairs

      TCanvas *L0_inv_mass_vs_L0bar_inv_mass_US_cuts_can = new TCanvas(Form("L0_inv_mass_vs_L0bar_inv_mass_US_cuts_can_%i_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0bar_inv_mass_US_cuts_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0_inv_mass_vs_L0bar_inv_mass_US_cuts_can->cd();

      L0_inv_mass_vs_L0bar_inv_mass_US_cuts[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0_inv_mass_vs_L0bar_inv_mass_US_cuts[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(2);
      L0_inv_mass_vs_L0bar_inv_mass_US_cuts[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_inv_mass_vs_L0bar_inv_mass_US_cuts[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0_inv_mass_vs_L0bar_inv_mass_US_cuts[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0_inv_mass_vs_L0bar_inv_mass_US_cuts[pTbin1][pTbin2]->GetYaxis()->SetTitleOffset(2.1);
      L0_inv_mass_vs_L0bar_inv_mass_US_cuts[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_inv_mass_vs_L0bar_inv_mass_US_cuts[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      //L0_inv_mass_vs_L0bar_inv_mass_US_cuts[pTbin1][pTbin2]->GetZaxis()->SetRangeUser(0,1000);
      L0_inv_mass_vs_L0bar_inv_mass_US_cuts[pTbin1][pTbin2]->Draw("surf1");

      TPaveText *L0_L0bar_text_Minv_2D_US_cuts = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
      L0_L0bar_text_Minv_2D_US_cuts->SetTextFont(42);
      //L0_L0bar_text_no_corr->AddText("STAR Internal");
      //L0_L0bar_text_no_corr->AddText("STAR preliminary");
      //((TText*)L0_L0bar_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
      //L0_L0bar_text_Minv_2D_US_cuts->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0_L0bar_text_Minv_2D_US_cuts->AddText("Minimum bias");
      L0_L0bar_text_Minv_2D_US_cuts->AddText("#Lambda^{0}-#bar{#Lambda^{0}}");
      L0_L0bar_text_Minv_2D_US_cuts->SetFillColorAlpha(0, 0.01);
      L0_L0bar_text_Minv_2D_US_cuts->Draw("same");

      L0_L0bar_text_Minv_2D_kine->Draw("same");

      L0_inv_mass_vs_L0bar_inv_mass_US_cuts_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/Minv/L0_inv_mass_vs_L0bar_inv_mass_US_cuts_pT1_%i_pT2_%i.png", pTbin1, pTbin2));



      TCanvas *L0_inv_mass_vs_L0bar_inv_mass_US_LS_cuts_can = new TCanvas(Form("L0_inv_mass_vs_L0bar_inv_mass_US_LS_cuts_can_%i_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0bar_inv_mass_US_LS_cuts_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0_inv_mass_vs_L0bar_inv_mass_US_LS_cuts_can->cd();

      L0_inv_mass_vs_L0bar_inv_mass_US_LS_cuts[pTbin1][pTbin2]->Sumw2();
      L0_inv_mass_vs_L0bar_inv_mass_US_LS_cuts[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0_inv_mass_vs_L0bar_inv_mass_US_LS_cuts[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(2);
      L0_inv_mass_vs_L0bar_inv_mass_US_LS_cuts[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_inv_mass_vs_L0bar_inv_mass_US_LS_cuts[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0_inv_mass_vs_L0bar_inv_mass_US_LS_cuts[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0_inv_mass_vs_L0bar_inv_mass_US_LS_cuts[pTbin1][pTbin2]->GetYaxis()->SetTitleOffset(2);
      L0_inv_mass_vs_L0bar_inv_mass_US_LS_cuts[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_inv_mass_vs_L0bar_inv_mass_US_LS_cuts[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      L0_inv_mass_vs_L0bar_inv_mass_US_LS_cuts[pTbin1][pTbin2]->Draw("surf1");

      TPaveText *L0_L0bar_text_Minv_2D_US_LS_cuts = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
      L0_L0bar_text_Minv_2D_US_LS_cuts->SetTextFont(42);
      //L0_L0bar_text_no_corr->AddText("STAR Internal");
      //L0_L0bar_text_no_corr->AddText("STAR preliminary");
      //((TText*)L0_L0bar_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
      //L0_L0bar_text_Minv_2D_US_LS_cuts->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0_L0bar_text_Minv_2D_US_LS_cuts->AddText("Minimum bias");
      L0_L0bar_text_Minv_2D_US_LS_cuts->AddText("#Lambda^{0}-#bar{#Lambda^{0}}");
      L0_L0bar_text_Minv_2D_US_LS_cuts->SetFillColorAlpha(0, 0.01);
      L0_L0bar_text_Minv_2D_US_LS_cuts->Draw("same");

      L0_L0bar_text_Minv_2D_kine->Draw("same");

      L0_inv_mass_vs_L0bar_inv_mass_US_LS_cuts_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/Minv/L0_inv_mass_vs_L0bar_inv_mass_US_LS_cuts_pT1_%i_pT2_%i.png", pTbin1, pTbin2));

      //-------------------------------------------------------------------------

      L0_inv_mass_vs_L0_inv_mass_US_cuts[pTbin1][pTbin2] = (TH2F*)inFile->Get(Form("L0_inv_mass_vs_L0_inv_mass_US_cuts_pT1_%i_pT2_%i", pTbin1, pTbin2)); //for US-US Lambda pairs
      L0_inv_mass_vs_L0_inv_mass_US_LS_cuts[pTbin1][pTbin2] = (TH2F*)inFile->Get(Form("L0_inv_mass_vs_L0_inv_mass_US_LS_cuts_pT1_%i_pT2_%i", pTbin1, pTbin2)); //for US-LS Lambda pairs

      TCanvas *L0_inv_mass_vs_L0_inv_mass_US_cuts_can = new TCanvas(Form("L0_inv_mass_vs_L0_inv_mass_US_cuts_can_%i_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0_inv_mass_US_cuts_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0_inv_mass_vs_L0_inv_mass_US_cuts_can->cd();

      L0_inv_mass_vs_L0_inv_mass_US_cuts[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0_inv_mass_vs_L0_inv_mass_US_cuts[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(2);
      L0_inv_mass_vs_L0_inv_mass_US_cuts[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_inv_mass_vs_L0_inv_mass_US_cuts[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0_inv_mass_vs_L0_inv_mass_US_cuts[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0_inv_mass_vs_L0_inv_mass_US_cuts[pTbin1][pTbin2]->GetYaxis()->SetTitleOffset(2.1);
      L0_inv_mass_vs_L0_inv_mass_US_cuts[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_inv_mass_vs_L0_inv_mass_US_cuts[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      //L0_inv_mass_vs_L0_inv_mass_US_cuts[pTbin1][pTbin2]->GetZaxis()->SetRangeUser(0,1000);
      L0_inv_mass_vs_L0_inv_mass_US_cuts[pTbin1][pTbin2]->Draw("surf1");

      TPaveText *L0_L0_text_Minv_2D_US_cuts = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
      L0_L0_text_Minv_2D_US_cuts->SetTextFont(42);
      //L0_L0_text_no_corr->AddText("STAR Internal");
      //L0_L0_text_no_corr->AddText("STAR preliminary");
      //((TText*)L0_L0_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
      //L0_L0_text_Minv_2D_US_cuts->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0_L0_text_Minv_2D_US_cuts->AddText("Minimum bias");
      L0_L0_text_Minv_2D_US_cuts->AddText("#Lambda^{0}-#Lambda^{0}");
      L0_L0_text_Minv_2D_US_cuts->SetFillColorAlpha(0, 0.01);
      L0_L0_text_Minv_2D_US_cuts->Draw("same");

      L0_L0_text_Minv_2D_kine->Draw("same");

      L0_inv_mass_vs_L0_inv_mass_US_cuts_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_L/Minv/L0_inv_mass_vs_L0_inv_mass_US_cuts_pT1_%i_pT2_%i.png", pTbin1, pTbin2));



      TCanvas *L0_inv_mass_vs_L0_inv_mass_US_LS_cuts_can = new TCanvas(Form("L0_inv_mass_vs_L0_inv_mass_US_LS_cuts_can_%i_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0_inv_mass_US_LS_cuts_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0_inv_mass_vs_L0_inv_mass_US_LS_cuts_can->cd();

      L0_inv_mass_vs_L0_inv_mass_US_LS_cuts[pTbin1][pTbin2]->Sumw2();
      L0_inv_mass_vs_L0_inv_mass_US_LS_cuts[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0_inv_mass_vs_L0_inv_mass_US_LS_cuts[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(2);
      L0_inv_mass_vs_L0_inv_mass_US_LS_cuts[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_inv_mass_vs_L0_inv_mass_US_LS_cuts[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0_inv_mass_vs_L0_inv_mass_US_LS_cuts[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0_inv_mass_vs_L0_inv_mass_US_LS_cuts[pTbin1][pTbin2]->GetYaxis()->SetTitleOffset(2);
      L0_inv_mass_vs_L0_inv_mass_US_LS_cuts[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_inv_mass_vs_L0_inv_mass_US_LS_cuts[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      L0_inv_mass_vs_L0_inv_mass_US_LS_cuts[pTbin1][pTbin2]->Draw("surf1");

      TPaveText *L0_L0_text_Minv_2D_US_LS_cuts = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
      L0_L0_text_Minv_2D_US_LS_cuts->SetTextFont(42);
      //L0_L0_text_no_corr->AddText("STAR Internal");
      //L0_L0_text_no_corr->AddText("STAR preliminary");
      //((TText*)L0_L0_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
      //L0_L0_text_Minv_2D_US_LS_cuts->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0_L0_text_Minv_2D_US_LS_cuts->AddText("Minimum bias");
      L0_L0_text_Minv_2D_US_LS_cuts->AddText("#Lambda^{0}-#Lambda^{0}");
      L0_L0_text_Minv_2D_US_LS_cuts->SetFillColorAlpha(0, 0.01);
      L0_L0_text_Minv_2D_US_LS_cuts->Draw("same");

      L0_L0_text_Minv_2D_kine->Draw("same");

      L0_inv_mass_vs_L0_inv_mass_US_LS_cuts_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_L/Minv/L0_inv_mass_vs_L0_inv_mass_US_LS_cuts_pT1_%i_pT2_%i.png", pTbin1, pTbin2));

      //----------------------------------------------------------------------------------------------------

      L0bar_inv_mass_vs_L0bar_inv_mass_US_cuts[pTbin1][pTbin2] = (TH2F*)inFile->Get(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_cuts_pT1_%i_pT2_%i", pTbin1, pTbin2)); //for US-US Lambda pairs
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_cuts[pTbin1][pTbin2] = (TH2F*)inFile->Get(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_cuts_pT1_%i_pT2_%i", pTbin1, pTbin2)); //for US-LS Lambda pairs

      TCanvas *L0bar_inv_mass_vs_L0bar_inv_mass_US_cuts_can = new TCanvas(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_cuts_can_%i_%i", pTbin1, pTbin2), Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_cuts_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_cuts_can->cd();

      L0bar_inv_mass_vs_L0bar_inv_mass_US_cuts[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0bar_inv_mass_vs_L0bar_inv_mass_US_cuts[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(2);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_cuts[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0bar_inv_mass_vs_L0bar_inv_mass_US_cuts[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_cuts[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0bar_inv_mass_vs_L0bar_inv_mass_US_cuts[pTbin1][pTbin2]->GetYaxis()->SetTitleOffset(2.1);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_cuts[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0bar_inv_mass_vs_L0bar_inv_mass_US_cuts[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      //L0bar_inv_mass_vs_L0bar_inv_mass_US_cuts[pTbin1][pTbin2]->GetZaxis()->SetRangeUser(0,1000);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_cuts[pTbin1][pTbin2]->Draw("surf1");

      TPaveText *L0bar_L0bar_text_Minv_2D_US_cuts = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
      L0bar_L0bar_text_Minv_2D_US_cuts->SetTextFont(42);
      //L0bar_L0bar_text_no_corr->AddText("STAR Internal");
      //L0bar_L0bar_text_no_corr->AddText("STAR preliminary");
      //((TText*)L0bar_L0bar_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
      //L0bar_L0bar_text_Minv_2D_US_cuts->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0bar_L0bar_text_Minv_2D_US_cuts->AddText("Minimum bias");
      L0bar_L0bar_text_Minv_2D_US_cuts->AddText("#bar{#Lambda^{0}}-#bar{#Lambda^{0}}");
      L0bar_L0bar_text_Minv_2D_US_cuts->SetFillColorAlpha(0, 0.01);
      L0bar_L0bar_text_Minv_2D_US_cuts->Draw("same");

      L0bar_L0bar_text_Minv_2D_kine->Draw("same");

      L0bar_inv_mass_vs_L0bar_inv_mass_US_cuts_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/Lbar_Lbar/Minv/L0bar_inv_mass_vs_L0bar_inv_mass_US_cuts_pT1_%i_pT2_%i.png", pTbin1, pTbin2));



      TCanvas *L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_cuts_can = new TCanvas(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_cuts_can_%i_%i", pTbin1, pTbin2), Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_cuts_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_cuts_can->cd();

      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_cuts[pTbin1][pTbin2]->Sumw2();
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_cuts[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_cuts[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(2);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_cuts[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_cuts[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_cuts[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_cuts[pTbin1][pTbin2]->GetYaxis()->SetTitleOffset(2);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_cuts[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_cuts[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_cuts[pTbin1][pTbin2]->Draw("surf1");

      TPaveText *L0bar_L0bar_text_Minv_2D_US_LS_cuts = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
      L0bar_L0bar_text_Minv_2D_US_LS_cuts->SetTextFont(42);
      //L0bar_L0bar_text_no_corr->AddText("STAR Internal");
      //L0bar_L0bar_text_no_corr->AddText("STAR preliminary");
      //((TText*)L0bar_L0bar_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
      //L0bar_L0bar_text_Minv_2D_US_LS_cuts->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0bar_L0bar_text_Minv_2D_US_LS_cuts->AddText("Minimum bias");
      L0bar_L0bar_text_Minv_2D_US_LS_cuts->AddText("#bar{#Lambda^{0}}-#bar{#Lambda^{0}}");
      L0bar_L0bar_text_Minv_2D_US_LS_cuts->SetFillColorAlpha(0, 0.01);
      L0bar_L0bar_text_Minv_2D_US_LS_cuts->Draw("same");

      L0bar_L0bar_text_Minv_2D_kine->Draw("same");

      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_cuts_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/Lbar_Lbar/Minv/L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_cuts_pT1_%i_pT2_%i.png", pTbin1, pTbin2));

      //__________________________________________________________________________________________________________________________________________

    }
  }

  float nLLbar_eta_US[nEtaBins][nEtaBins];
  float nLL_eta_US[nEtaBins][nEtaBins];
  float nLbarLbar_eta_US[nEtaBins][nEtaBins];

  float nLLbar_eta_US_LS[nEtaBins][nEtaBins];
  float nLL_eta_US_LS[nEtaBins][nEtaBins];
  float nLbarLbar_eta_US_LS[nEtaBins][nEtaBins];


  float nLLbar_eta_US_cuts[nEtaBins][nEtaBins];
  float nLL_eta_US_cuts[nEtaBins][nEtaBins];
  float nLbarLbar_eta_US_cuts[nEtaBins][nEtaBins];

  float nLLbar_eta_US_LS_cuts[nEtaBins][nEtaBins];
  float nLL_eta_US_LS_cuts[nEtaBins][nEtaBins];
  float nLbarLbar_eta_US_LS_cuts[nEtaBins][nEtaBins];

  for(unsigned int etaBin1 = 0; etaBin1 < nEtaBins; etaBin1++)
  {
    for(unsigned int etaBin2 = 0; etaBin2 < nEtaBins; etaBin2++)
    {
      //L-Lbar
      L0_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2] = (TH1F*)inFile->Get(Form("L0_L0bar_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));
      L0_L0bar_cosThetaProdPlane_US_LS_eta_hist[etaBin1][etaBin2] = (TH1F*)inFile->Get(Form("L0_L0bar_cosThetaProdPlane_US_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));

      L0_L0bar_cosThetaProdPlane_US_ME_weight_eta_hist[etaBin1][etaBin2] = (TH1F*)inFile->Get(Form("L0_L0bar_cosThetaProdPlane_US_ME_weight_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));
      L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_hist[etaBin1][etaBin2] = (TH1F*)inFile->Get(Form("L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));


      TCanvas *L0_L0bar_cosThetaProdPlane_US_eta_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_US_eta_can_eta1_%i_eta2_%i", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_US_eta_can_eta1_%i_eta2_%i", etaBin1, etaBin2), 1200, 1000);

      //correct background first
      nLLbar_eta_US_LS[etaBin1][etaBin2] = L0_L0bar_cosThetaProdPlane_US_LS_eta_hist[etaBin1][etaBin2]->Integral();

      L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_hist[etaBin1][etaBin2]->Sumw2();
      L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_hist[etaBin1][etaBin2]->Scale(nLLbar_eta_US_LS[etaBin1][etaBin2]/L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_hist[etaBin1][etaBin2]->Integral());

      L0_L0bar_cosThetaProdPlane_US_LS_eta_hist[etaBin1][etaBin2]->Sumw2();
      L0_L0bar_cosThetaProdPlane_US_LS_eta_hist[etaBin1][etaBin2]->Divide(L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_hist[etaBin1][etaBin2]);
      L0_L0bar_cosThetaProdPlane_US_LS_eta_hist[etaBin1][etaBin2]->Scale(nLLbar_eta_US_LS[etaBin1][etaBin2]/L0_L0bar_cosThetaProdPlane_US_ME_weight_eta_hist[etaBin1][etaBin2]->Integral());
      //L0_L0bar_cosThetaProdPlane_US_LS_eta_hist[etaBin1][etaBin2]->Scale(1./L0_L0bar_cosThetaProdPlane_US_ME_weight_eta_hist[etaBin1][etaBin2]->GetXaxis()->GetBinWidth(1));


      //correct signal+background
      nLLbar_eta_US[etaBin1][etaBin2] = L0_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->Integral();

      L0_L0bar_cosThetaProdPlane_US_ME_weight_eta_hist[etaBin1][etaBin2]->Sumw2();
      L0_L0bar_cosThetaProdPlane_US_ME_weight_eta_hist[etaBin1][etaBin2]->Scale(nLLbar_eta_US[etaBin1][etaBin2]/L0_L0bar_cosThetaProdPlane_US_ME_weight_eta_hist[etaBin1][etaBin2]->Integral());

      L0_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->Sumw2();
      L0_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      L0_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->SetMarkerColor(kRed);
      L0_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->SetLineColor(kRed);
      L0_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->GetXaxis()->CenterTitle();
      L0_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      L0_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->GetYaxis()->CenterTitle();
      L0_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->Divide(L0_L0bar_cosThetaProdPlane_US_ME_weight_eta_hist[etaBin1][etaBin2]);
      L0_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->Scale(nLLbar_eta_US[etaBin1][etaBin2]/L0_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->Integral());
      L0_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->Add(L0_L0bar_cosThetaProdPlane_US_LS_eta_hist[etaBin1][etaBin2], -1);
      L0_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->Scale(1./L0_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->GetXaxis()->GetBinWidth(1));
      L0_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->SetMinimum(0);
      L0_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->Draw("p e");


      TPaveText *L0_L0bar_text_eta = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
      L0_L0bar_text_eta->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
      L0_L0bar_text_eta->AddText("Minimum bias");
      L0_L0bar_text_eta->AddText("#Lambda-#bar{#Lambda}");
      L0_L0bar_text_eta->AddText("US - Background");
      L0_L0bar_text_eta->AddText("No cuts");
      L0_L0bar_text_eta->AddText(Form("%0.1f < #eta_{1} < %0.1f", eta_bins[etaBin1], eta_bins[etaBin1+1]));
      L0_L0bar_text_eta->AddText(Form("%0.1f < #eta_{2} < %0.1f", eta_bins[etaBin2], eta_bins[etaBin2+1]));
      L0_L0bar_text_eta->AddText("p_{T} integrated");
      L0_L0bar_text_eta->SetFillColorAlpha(0, 0.01);
      L0_L0bar_text_eta->Draw("same");

      L0_L0bar_cosThetaProdPlane_US_eta_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/L0_L0bar_cosThetaProdPlane_US_eta1_%i_eta2_%i.png", etaBin1, etaBin2));

      //--------------------------------------------------------------------------------

      L0_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2] = (TH1F*)inFile->Get(Form("L0_L0bar_cosThetaProdPlane_US_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));
      L0_L0bar_cosThetaProdPlane_US_LS_eta_cuts_hist[etaBin1][etaBin2] = (TH1F*)inFile->Get(Form("L0_L0bar_cosThetaProdPlane_US_LS_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));

      L0_L0bar_cosThetaProdPlane_US_ME_weight_eta_cuts_hist[etaBin1][etaBin2] = (TH1F*)inFile->Get(Form("L0_L0bar_cosThetaProdPlane_US_ME_weight_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));
      L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_cuts_hist[etaBin1][etaBin2] = (TH1F*)inFile->Get(Form("L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));


      TCanvas *L0_L0bar_cosThetaProdPlane_US_eta_cuts_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_US_eta_cuts_can_eta1_%i_eta2_%i", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_US_eta_cuts_can_eta1_%i_eta2_%i", etaBin1, etaBin2), 1200, 1000);

      //correct background first
      nLLbar_eta_US_LS_cuts[etaBin1][etaBin2] = L0_L0bar_cosThetaProdPlane_US_LS_eta_cuts_hist[etaBin1][etaBin2]->Integral();

      L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_cuts_hist[etaBin1][etaBin2]->Sumw2();
      L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_cuts_hist[etaBin1][etaBin2]->Scale(nLLbar_eta_US_LS_cuts[etaBin1][etaBin2]/L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_cuts_hist[etaBin1][etaBin2]->Integral());

      L0_L0bar_cosThetaProdPlane_US_LS_eta_cuts_hist[etaBin1][etaBin2]->Sumw2();
      L0_L0bar_cosThetaProdPlane_US_LS_eta_cuts_hist[etaBin1][etaBin2]->Divide(L0_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_cuts_hist[etaBin1][etaBin2]);
      L0_L0bar_cosThetaProdPlane_US_LS_eta_cuts_hist[etaBin1][etaBin2]->Scale(nLLbar_eta_US_LS_cuts[etaBin1][etaBin2]/L0_L0bar_cosThetaProdPlane_US_ME_weight_eta_cuts_hist[etaBin1][etaBin2]->Integral());
      //L0_L0bar_cosThetaProdPlane_US_LS_eta_cuts_hist[etaBin1][etaBin2]->Scale(1./L0_L0bar_cosThetaProdPlane_US_ME_weight_eta_cuts_hist[etaBin1][etaBin2]->GetXaxis()->GetBinWidth(1));


      //correct signal+background
      nLLbar_eta_US_cuts[etaBin1][etaBin2] = L0_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->Integral();

      L0_L0bar_cosThetaProdPlane_US_ME_weight_eta_cuts_hist[etaBin1][etaBin2]->Sumw2();
      L0_L0bar_cosThetaProdPlane_US_ME_weight_eta_cuts_hist[etaBin1][etaBin2]->Scale(nLLbar_eta_US_cuts[etaBin1][etaBin2]/L0_L0bar_cosThetaProdPlane_US_ME_weight_eta_cuts_hist[etaBin1][etaBin2]->Integral());

      L0_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->Sumw2();
      L0_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      L0_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->SetMarkerColor(kRed);
      L0_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->SetLineColor(kRed);
      L0_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->GetXaxis()->CenterTitle();
      L0_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      L0_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->GetYaxis()->CenterTitle();
      L0_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->Divide(L0_L0bar_cosThetaProdPlane_US_ME_weight_eta_cuts_hist[etaBin1][etaBin2]);
      L0_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->Scale(nLLbar_eta_US_cuts[etaBin1][etaBin2]/L0_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->Integral());
      L0_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->Add(L0_L0bar_cosThetaProdPlane_US_LS_eta_cuts_hist[etaBin1][etaBin2], -1);
      L0_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->Scale(1./L0_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->GetXaxis()->GetBinWidth(1));
      L0_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->SetMinimum(0);
      L0_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->Draw("p e");

      TPaveText *L0_L0bar_text_cuts_eta = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
      L0_L0bar_text_cuts_eta->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
      L0_L0bar_text_cuts_eta->AddText("Minimum bias");
      L0_L0bar_text_cuts_eta->AddText("#Lambda-#bar{#Lambda}");
      L0_L0bar_text_cuts_eta->AddText("US - Background");
      L0_L0bar_text_cuts_eta->AddText("Analysis cuts");
      L0_L0bar_text_cuts_eta->AddText(Form("%0.1f < #eta_{1} < %0.1f", eta_bins[etaBin1], eta_bins[etaBin1+1]));
      L0_L0bar_text_cuts_eta->AddText(Form("%0.1f < #eta_{2} < %0.1f", eta_bins[etaBin2], eta_bins[etaBin2+1]));
      L0_L0bar_text_cuts_eta->AddText("p_{T} integrated");
      L0_L0bar_text_cuts_eta->SetFillColorAlpha(0, 0.01);
      L0_L0bar_text_cuts_eta->Draw("same");

      L0_L0bar_cosThetaProdPlane_US_eta_cuts_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/L0_L0bar_cosThetaProdPlane_US_cuts_eta1_%i_eta2_%i.png", etaBin1, etaBin2));

      //_________________________________________________________________________________________________

      //L-L
      L0_L0_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2] = (TH1F*)inFile->Get(Form("L0_L0_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));
      L0_L0_cosThetaProdPlane_US_LS_eta_hist[etaBin1][etaBin2] = (TH1F*)inFile->Get(Form("L0_L0_cosThetaProdPlane_US_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));

      L0_L0_cosThetaProdPlane_US_ME_weight_eta_hist[etaBin1][etaBin2] = (TH1F*)inFile->Get(Form("L0_L0_cosThetaProdPlane_US_ME_weight_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));
      L0_L0_cosThetaProdPlane_US_LS_ME_weight_eta_hist[etaBin1][etaBin2] = (TH1F*)inFile->Get(Form("L0_L0_cosThetaProdPlane_US_LS_ME_weight_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));


      TCanvas *L0_L0_cosThetaProdPlane_US_eta_can = new TCanvas(Form("L0_L0_cosThetaProdPlane_US_eta_can_eta1_%i_eta2_%i", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_US_eta_can_eta1_%i_eta2_%i", etaBin1, etaBin2), 1200, 1000);

      //correct background first
      nLL_eta_US_LS[etaBin1][etaBin2] = L0_L0_cosThetaProdPlane_US_LS_eta_hist[etaBin1][etaBin2]->Integral();

      L0_L0_cosThetaProdPlane_US_LS_ME_weight_eta_hist[etaBin1][etaBin2]->Sumw2();
      L0_L0_cosThetaProdPlane_US_LS_ME_weight_eta_hist[etaBin1][etaBin2]->Scale(nLL_eta_US_LS[etaBin1][etaBin2]/L0_L0_cosThetaProdPlane_US_LS_ME_weight_eta_hist[etaBin1][etaBin2]->Integral());

      L0_L0_cosThetaProdPlane_US_LS_eta_hist[etaBin1][etaBin2]->Sumw2();
      L0_L0_cosThetaProdPlane_US_LS_eta_hist[etaBin1][etaBin2]->Divide(L0_L0_cosThetaProdPlane_US_LS_ME_weight_eta_hist[etaBin1][etaBin2]);
      L0_L0_cosThetaProdPlane_US_LS_eta_hist[etaBin1][etaBin2]->Scale(nLL_eta_US_LS[etaBin1][etaBin2]/L0_L0_cosThetaProdPlane_US_ME_weight_eta_hist[etaBin1][etaBin2]->Integral());
      //L0_L0_cosThetaProdPlane_US_LS_eta_hist[etaBin1][etaBin2]->Scale(1./L0_L0_cosThetaProdPlane_US_ME_weight_eta_hist[etaBin1][etaBin2]->GetXaxis()->GetBinWidth(1));


      //correct signal+background
      nLL_eta_US[etaBin1][etaBin2] = L0_L0_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->Integral();

      L0_L0_cosThetaProdPlane_US_ME_weight_eta_hist[etaBin1][etaBin2]->Sumw2();
      L0_L0_cosThetaProdPlane_US_ME_weight_eta_hist[etaBin1][etaBin2]->Scale(nLL_eta_US[etaBin1][etaBin2]/L0_L0_cosThetaProdPlane_US_ME_weight_eta_hist[etaBin1][etaBin2]->Integral());

      L0_L0_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->Sumw2();
      L0_L0_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      L0_L0_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->SetMarkerColor(kRed);
      L0_L0_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->SetLineColor(kRed);
      L0_L0_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0_L0_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->GetXaxis()->CenterTitle();
      L0_L0_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      L0_L0_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->GetYaxis()->CenterTitle();
      L0_L0_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->Divide(L0_L0_cosThetaProdPlane_US_ME_weight_eta_hist[etaBin1][etaBin2]);
      L0_L0_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->Scale(nLL_eta_US[etaBin1][etaBin2]/L0_L0_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->Integral());
      L0_L0_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->Add(L0_L0_cosThetaProdPlane_US_LS_eta_hist[etaBin1][etaBin2], -1);
      L0_L0_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->Scale(1./L0_L0_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->GetXaxis()->GetBinWidth(1));
      L0_L0_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->SetMinimum(0);
      L0_L0_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->Draw("p e");

      TPaveText *L0_L0_text_eta = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
      L0_L0_text_eta->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
      L0_L0_text_eta->AddText("Minimum bias");
      L0_L0_text_eta->AddText("#Lambda-#Lambda");
      L0_L0_text_eta->AddText("US matched to US");
      L0_L0_text_eta->AddText("No cuts");
      L0_L0_text_eta->AddText(Form("%0.1f < #eta_{1} < %0.1f", eta_bins[etaBin1], eta_bins[etaBin1+1]));
      L0_L0_text_eta->AddText(Form("%0.1f < #eta_{2} < %0.1f", eta_bins[etaBin2], eta_bins[etaBin2+1]));
      L0_L0_text_eta->AddText("p_{T} integrated");
      L0_L0_text_eta->SetFillColorAlpha(0, 0.01);
      L0_L0_text_eta->Draw("same");

      L0_L0_cosThetaProdPlane_US_eta_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_L/L0_L0_cosThetaProdPlane_US_eta1_%i_eta2_%i.png", etaBin1, etaBin2));

      //--------------------------------------------------------------------------------

      L0_L0_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2] = (TH1F*)inFile->Get(Form("L0_L0_cosThetaProdPlane_US_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));
      L0_L0_cosThetaProdPlane_US_LS_eta_cuts_hist[etaBin1][etaBin2] = (TH1F*)inFile->Get(Form("L0_L0_cosThetaProdPlane_US_LS_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));

      L0_L0_cosThetaProdPlane_US_ME_weight_eta_cuts_hist[etaBin1][etaBin2] = (TH1F*)inFile->Get(Form("L0_L0_cosThetaProdPlane_US_ME_weight_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));
      L0_L0_cosThetaProdPlane_US_LS_ME_weight_eta_cuts_hist[etaBin1][etaBin2] = (TH1F*)inFile->Get(Form("L0_L0_cosThetaProdPlane_US_LS_ME_weight_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));


      TCanvas *L0_L0_cosThetaProdPlane_US_eta_cuts_can = new TCanvas(Form("L0_L0_cosThetaProdPlane_US_eta_cuts_can_eta1_%i_eta2_%i", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_US_eta_cuts_can_eta1_%i_eta2_%i", etaBin1, etaBin2), 1200, 1000);

      //correct background first
      nLL_eta_US_LS_cuts[etaBin1][etaBin2] = L0_L0_cosThetaProdPlane_US_LS_eta_cuts_hist[etaBin1][etaBin2]->Integral();

      L0_L0_cosThetaProdPlane_US_LS_ME_weight_eta_cuts_hist[etaBin1][etaBin2]->Sumw2();
      L0_L0_cosThetaProdPlane_US_LS_ME_weight_eta_cuts_hist[etaBin1][etaBin2]->Scale(nLL_eta_US_LS_cuts[etaBin1][etaBin2]/L0_L0_cosThetaProdPlane_US_LS_ME_weight_eta_cuts_hist[etaBin1][etaBin2]->Integral());

      L0_L0_cosThetaProdPlane_US_LS_eta_cuts_hist[etaBin1][etaBin2]->Sumw2();
      L0_L0_cosThetaProdPlane_US_LS_eta_cuts_hist[etaBin1][etaBin2]->Divide(L0_L0_cosThetaProdPlane_US_LS_ME_weight_eta_cuts_hist[etaBin1][etaBin2]);
      L0_L0_cosThetaProdPlane_US_LS_eta_cuts_hist[etaBin1][etaBin2]->Scale(nLL_eta_US_LS_cuts[etaBin1][etaBin2]/L0_L0_cosThetaProdPlane_US_ME_weight_eta_cuts_hist[etaBin1][etaBin2]->Integral());
      //L0_L0_cosThetaProdPlane_US_LS_eta_cuts_hist[etaBin1][etaBin2]->Scale(1./L0_L0_cosThetaProdPlane_US_ME_weight_eta_cuts_hist[etaBin1][etaBin2]->GetXaxis()->GetBinWidth(1));


      //correct signal+background
      nLL_eta_US_cuts[etaBin1][etaBin2] = L0_L0_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->Integral();

      L0_L0_cosThetaProdPlane_US_ME_weight_eta_cuts_hist[etaBin1][etaBin2]->Sumw2();
      L0_L0_cosThetaProdPlane_US_ME_weight_eta_cuts_hist[etaBin1][etaBin2]->Scale(nLL_eta_US_cuts[etaBin1][etaBin2]/L0_L0_cosThetaProdPlane_US_ME_weight_eta_cuts_hist[etaBin1][etaBin2]->Integral());

      L0_L0_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->Sumw2();
      L0_L0_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      L0_L0_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->SetMarkerColor(kRed);
      L0_L0_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->SetLineColor(kRed);
      L0_L0_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0_L0_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->GetXaxis()->CenterTitle();
      L0_L0_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      L0_L0_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->GetYaxis()->CenterTitle();
      L0_L0_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->Divide(L0_L0_cosThetaProdPlane_US_ME_weight_eta_cuts_hist[etaBin1][etaBin2]);
      L0_L0_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->Scale(nLL_eta_US_cuts[etaBin1][etaBin2]/L0_L0_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->Integral());
      L0_L0_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->Add(L0_L0_cosThetaProdPlane_US_LS_eta_cuts_hist[etaBin1][etaBin2], -1);
      L0_L0_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->Scale(1./L0_L0_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->GetXaxis()->GetBinWidth(1));
      L0_L0_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->SetMinimum(0);
      L0_L0_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->Draw("p e");


      TPaveText *L0_L0_text_cuts_eta = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
      L0_L0_text_cuts_eta->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
      L0_L0_text_cuts_eta->AddText("Minimum bias");
      L0_L0_text_cuts_eta->AddText("#Lambda-#Lambda");
      L0_L0_text_cuts_eta->AddText("US matched to US");
      L0_L0_text_cuts_eta->AddText("Analysis cuts");
      L0_L0_text_cuts_eta->AddText(Form("%0.1f < #eta_{1} < %0.1f", eta_bins[etaBin1], eta_bins[etaBin1+1]));
      L0_L0_text_cuts_eta->AddText(Form("%0.1f < #eta_{2} < %0.1f", eta_bins[etaBin2], eta_bins[etaBin2+1]));
      L0_L0_text_cuts_eta->AddText("p_{T} integrated");
      L0_L0_text_cuts_eta->SetFillColorAlpha(0, 0.01);
      L0_L0_text_cuts_eta->Draw("same");

      L0_L0_cosThetaProdPlane_US_eta_cuts_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_L/L0_L0_cosThetaProdPlane_US_cuts_eta1_%i_eta2_%i.png", etaBin1, etaBin2));


      //_________________________________________________________________________________________________

      //Lbar-Lbar
      L0bar_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2] = (TH1F*)inFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));
      L0bar_L0bar_cosThetaProdPlane_US_LS_eta_hist[etaBin1][etaBin2] = (TH1F*)inFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_US_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));

      L0bar_L0bar_cosThetaProdPlane_US_ME_weight_eta_hist[etaBin1][etaBin2] = (TH1F*)inFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_US_ME_weight_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));
      L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_hist[etaBin1][etaBin2] = (TH1F*)inFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));


      TCanvas *L0bar_L0bar_cosThetaProdPlane_US_eta_can = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_US_eta_can_eta1_%i_eta2_%i", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_US_eta_can_eta1_%i_eta2_%i", etaBin1, etaBin2), 1200, 1000);

      //correct background first
      nLbarLbar_eta_US_LS[etaBin1][etaBin2] = L0bar_L0bar_cosThetaProdPlane_US_LS_eta_hist[etaBin1][etaBin2]->Integral();

      L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_hist[etaBin1][etaBin2]->Sumw2();
      L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_hist[etaBin1][etaBin2]->Scale(nLbarLbar_eta_US_LS[etaBin1][etaBin2]/L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_hist[etaBin1][etaBin2]->Integral());

      L0bar_L0bar_cosThetaProdPlane_US_LS_eta_hist[etaBin1][etaBin2]->Sumw2();
      L0bar_L0bar_cosThetaProdPlane_US_LS_eta_hist[etaBin1][etaBin2]->Divide(L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_hist[etaBin1][etaBin2]);
      L0bar_L0bar_cosThetaProdPlane_US_LS_eta_hist[etaBin1][etaBin2]->Scale(nLbarLbar_eta_US_LS[etaBin1][etaBin2]/L0bar_L0bar_cosThetaProdPlane_US_ME_weight_eta_hist[etaBin1][etaBin2]->Integral());
      //L0bar_L0bar_cosThetaProdPlane_US_LS_eta_hist[etaBin1][etaBin2]->Scale(1./L0bar_L0bar_cosThetaProdPlane_US_ME_weight_eta_hist[etaBin1][etaBin2]->GetXaxis()->GetBinWidth(1));


      //correct signal+background
      nLbarLbar_eta_US[etaBin1][etaBin2] = L0bar_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->Integral();

      L0bar_L0bar_cosThetaProdPlane_US_ME_weight_eta_hist[etaBin1][etaBin2]->Sumw2();
      L0bar_L0bar_cosThetaProdPlane_US_ME_weight_eta_hist[etaBin1][etaBin2]->Scale(nLbarLbar_eta_US[etaBin1][etaBin2]/L0bar_L0bar_cosThetaProdPlane_US_ME_weight_eta_hist[etaBin1][etaBin2]->Integral());

      L0bar_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->Sumw2();
      L0bar_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      L0bar_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->SetMarkerColor(kRed);
      L0bar_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->SetLineColor(kRed);
      L0bar_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0bar_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->GetXaxis()->CenterTitle();
      L0bar_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      L0bar_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->GetYaxis()->CenterTitle();
      L0bar_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->Divide(L0bar_L0bar_cosThetaProdPlane_US_ME_weight_eta_hist[etaBin1][etaBin2]);
      L0bar_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->Scale(nLbarLbar_eta_US[etaBin1][etaBin2]/L0bar_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->Integral());
      L0bar_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->Add(L0bar_L0bar_cosThetaProdPlane_US_LS_eta_hist[etaBin1][etaBin2], -1);
      L0bar_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->Scale(1./L0bar_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->GetXaxis()->GetBinWidth(1));
      L0bar_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->SetMinimum(0);
      L0bar_L0bar_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->Draw("p e");

      TPaveText *L0bar_L0bar_text_eta = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
      L0bar_L0bar_text_eta->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
      L0bar_L0bar_text_eta->AddText("Minimum bias");
      L0bar_L0bar_text_eta->AddText("#bar{#Lambda}-#bar{#Lambda}");
      L0bar_L0bar_text_eta->AddText("No cuts");
      L0bar_L0bar_text_eta->AddText(Form("%0.1f < #eta_{1} < %0.1f", eta_bins[etaBin1], eta_bins[etaBin1+1]));
      L0bar_L0bar_text_eta->AddText(Form("%0.1f < #eta_{2} < %0.1f", eta_bins[etaBin2], eta_bins[etaBin2+1]));
      L0bar_L0bar_text_eta->AddText("p_{T} integrated");
      L0bar_L0bar_text_eta->SetFillColorAlpha(0, 0.01);
      L0bar_L0bar_text_eta->Draw("same");


      L0bar_L0bar_cosThetaProdPlane_US_eta_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/Lbar_Lbar/L0bar_L0bar_cosThetaProdPlane_US_eta1_%i_eta2_%i.png", etaBin1, etaBin2));

      //--------------------------------------------------------------------------------

      L0bar_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2] = (TH1F*)inFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_US_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));
      L0bar_L0bar_cosThetaProdPlane_US_LS_eta_cuts_hist[etaBin1][etaBin2] = (TH1F*)inFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_US_LS_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));

      L0bar_L0bar_cosThetaProdPlane_US_ME_weight_eta_cuts_hist[etaBin1][etaBin2] = (TH1F*)inFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_US_ME_weight_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));
      L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_cuts_hist[etaBin1][etaBin2] = (TH1F*)inFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));


      TCanvas *L0bar_L0bar_cosThetaProdPlane_US_eta_cuts_can = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_US_eta_cuts_can_eta1_%i_eta2_%i", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_US_eta_cuts_can_eta1_%i_eta2_%i", etaBin1, etaBin2), 1200, 1000);

      //correct background first
      nLbarLbar_eta_US_LS_cuts[etaBin1][etaBin2] = L0bar_L0bar_cosThetaProdPlane_US_LS_eta_cuts_hist[etaBin1][etaBin2]->Integral();

      L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_cuts_hist[etaBin1][etaBin2]->Sumw2();
      L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_cuts_hist[etaBin1][etaBin2]->Scale(nLbarLbar_eta_US_LS_cuts[etaBin1][etaBin2]/L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_cuts_hist[etaBin1][etaBin2]->Integral());

      L0bar_L0bar_cosThetaProdPlane_US_LS_eta_cuts_hist[etaBin1][etaBin2]->Sumw2();
      L0bar_L0bar_cosThetaProdPlane_US_LS_eta_cuts_hist[etaBin1][etaBin2]->Divide(L0bar_L0bar_cosThetaProdPlane_US_LS_ME_weight_eta_cuts_hist[etaBin1][etaBin2]);
      L0bar_L0bar_cosThetaProdPlane_US_LS_eta_cuts_hist[etaBin1][etaBin2]->Scale(nLbarLbar_eta_US_LS_cuts[etaBin1][etaBin2]/L0bar_L0bar_cosThetaProdPlane_US_ME_weight_eta_cuts_hist[etaBin1][etaBin2]->Integral());
      //L0bar_L0bar_cosThetaProdPlane_US_LS_eta_cuts_hist[etaBin1][etaBin2]->Scale(1./L0bar_L0bar_cosThetaProdPlane_US_ME_weight_eta_cuts_hist[etaBin1][etaBin2]->GetXaxis()->GetBinWidth(1));


      //correct signal+background
      nLbarLbar_eta_US_cuts[etaBin1][etaBin2] = L0bar_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->Integral();

      L0bar_L0bar_cosThetaProdPlane_US_ME_weight_eta_cuts_hist[etaBin1][etaBin2]->Sumw2();
      L0bar_L0bar_cosThetaProdPlane_US_ME_weight_eta_cuts_hist[etaBin1][etaBin2]->Scale(nLbarLbar_eta_US_cuts[etaBin1][etaBin2]/L0bar_L0bar_cosThetaProdPlane_US_ME_weight_eta_cuts_hist[etaBin1][etaBin2]->Integral());

      L0bar_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->Sumw2();
      L0bar_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      L0bar_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->SetMarkerColor(kRed);
      L0bar_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->SetLineColor(kRed);
      L0bar_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0bar_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->GetXaxis()->CenterTitle();
      L0bar_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      L0bar_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->GetYaxis()->CenterTitle();
      L0bar_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->Divide(L0bar_L0bar_cosThetaProdPlane_US_ME_weight_eta_cuts_hist[etaBin1][etaBin2]);
      L0bar_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->Scale(nLbarLbar_eta_US_cuts[etaBin1][etaBin2]/L0bar_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->Integral());
      L0bar_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->Add(L0bar_L0bar_cosThetaProdPlane_US_LS_eta_cuts_hist[etaBin1][etaBin2], -1);
      L0bar_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->Scale(1./L0bar_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->GetXaxis()->GetBinWidth(1));
      L0bar_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->SetMinimum(0);
      L0bar_L0bar_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->Draw("p e");

      TPaveText *L0bar_L0bar_text_cuts_eta = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
      L0bar_L0bar_text_cuts_eta->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
      L0bar_L0bar_text_cuts_eta->AddText("Minimum bias");
      L0bar_L0bar_text_cuts_eta->AddText("#bar{#Lambda}-#bar{#Lambda}");
      L0bar_L0bar_text_cuts_eta->AddText("Analysis cuts");
      L0bar_L0bar_text_cuts_eta->AddText(Form("%0.1f < #eta_{1} < %0.1f", eta_bins[etaBin1], eta_bins[etaBin1+1]));
      L0bar_L0bar_text_cuts_eta->AddText(Form("%0.1f < #eta_{2} < %0.1f", eta_bins[etaBin2], eta_bins[etaBin2+1]));
      L0bar_L0bar_text_cuts_eta->AddText("p_{T} integrated");
      L0bar_L0bar_text_cuts_eta->SetFillColorAlpha(0, 0.01);
      L0bar_L0bar_text_cuts_eta->Draw("same");

      L0bar_L0bar_cosThetaProdPlane_US_eta_cuts_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/Lbar_Lbar/L0bar_L0bar_cosThetaProdPlane_US_cuts_eta1_%i_eta2_%i.png", etaBin1, etaBin2));

      //_________________________________________________________________________________________________

    }
  }

  SysErrSlope_hist->Write();


  inFile->Close();
  sysErrFile->Close();

  //outFile->Close();

  return;
}
