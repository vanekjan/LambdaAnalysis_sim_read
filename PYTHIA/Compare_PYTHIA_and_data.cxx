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

const int nPtBins = 8;
float const pT_bins[nPtBins+1] = { 0., 0.5, 1.,1.5, 2., 2.5, 3., 4., 5.};

void Compare_PYTHIA_and_data()
{
  //load all files
  //TFile *inFile_PYTHIA = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/output_Lambda_pp_200_MB_1B_events_hists_work.root", "READ");

  //TFile *inFile_PYTHIA = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/ME_tests_2024/output_Lambda_pp_200_MB_1B_events_hists_ME_with_SE.root", "READ");
  TFile *inFile_PYTHIA = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/2024_08_ME_weight_new/output_Lambda_pp_200_MB_1B_events_hists_Delta_phi_pi_third.root", "READ");
  //TFile *inFile_PYTHIA = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/ME_tests_2024/output_Lambda_pp_200_MB_1B_events_hists_ME_with_SE_pT_0.05.root", "READ");

  //TFile *inFile_data_Minv = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/2012/InvMass_Lambda_ana_cuts_012124.root", "read");
  TFile *inFile_data_Minv = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/2012/InvMass_Lambda_ana_cuts_work.root", "read");

  //TFile *inFile_data = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/2012/ProdPlane_Lambda_ana_cuts_standard_new_QA_hists.root", "read");
  TFile *inFile_data = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/2012/ProdPlane_Lambda_ana_cuts_work.root", "read");

  //TFile *inFile_embed = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Embedding/pp/output/RC_test/Out_file_RC_200_work.root", "read");
  TFile *inFile_embed = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Embedding/pp/output/RC_test/Out_file_RC_200_default.root", "read");
  //TFile *inFile_embed = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Embedding/pp/output/RC_test/Out_file_RC_200_straight_lines.root", "read");


  TFile *outFile_weight = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/output/weights/Weights_kinematics.root", "recreate");


  //QA histograms from PYTHIA
  //before cuts

  //daughter kinematics

  TH2F *L0_L0bar_p_pT1_vs_p_pT2_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0bar_p_pT1_vs_p_pT2_hist");
  TH2F *L0_L0bar_pi_pT1_vs_pi_pT2_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0bar_pi_pT1_vs_pi_pT2_hist");

  //------------------------------------------------------------------------------------------------------------

  //after cuts

  //L and Lbar pair kinematics

  TH2F *L0_L0bar_pT1_vs_pT2_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0bar_pT1_vs_pT2_cuts_hist");
  TH2F *L0_L0bar_eta1_vs_eta2_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0bar_eta1_vs_eta2_cuts_hist");
  TH2F *L0_L0bar_phi1_vs_phi2_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0bar_phi1_vs_phi2_cuts_hist");

  TH2F *L0_L0bar_pT1_vs_pT2_ME_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0bar_pT1_vs_pT2_ME_cuts_hist");
  TH2F *L0_L0bar_eta1_vs_eta2_ME_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0bar_eta1_vs_eta2_ME_cuts_hist");
  TH2F *L0_L0bar_phi1_vs_phi2_ME_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0bar_phi1_vs_phi2_ME_cuts_hist");

  //-------------------------------------------

  TH2F *L0_L0bar_delta_phi_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0bar_delta_phi_cuts_hist");
  TH2F *L0_L0bar_delta_eta_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0bar_delta_eta_cuts_hist");

  TH2F *L0_L0bar_delta_phi_ME_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0bar_delta_phi_ME_cuts_hist");
  TH2F *L0_L0bar_delta_eta_ME_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0bar_delta_eta_ME_cuts_hist");


  TH2F *L0_L0_pT1_vs_pT2_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0_pT1_vs_pT2_cuts_hist");
  TH2F *L0_L0_eta1_vs_eta2_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0_eta1_vs_eta2_cuts_hist");
  TH2F *L0_L0_phi1_vs_phi2_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0_phi1_vs_phi2_cuts_hist");

  TH2F *L0_L0_delta_phi_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0_delta_phi_cuts_hist");


  TH2F *L0bar_L0bar_pT1_vs_pT2_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0bar_L0bar_pT1_vs_pT2_cuts_hist");
  TH2F *L0bar_L0bar_eta1_vs_eta2_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0bar_L0bar_eta1_vs_eta2_cuts_hist");
  TH2F *L0bar_L0bar_phi1_vs_phi2_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0bar_L0bar_phi1_vs_phi2_cuts_hist");

  TH2F *L0bar_L0bar_delta_phi_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0bar_L0bar_delta_phi_cuts_hist");


  //daughter kinematics

  TH2F *L0_L0bar_p_pT1_vs_p_pT2_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0bar_p_pT1_vs_p_pT2_cuts_hist");
  TH2F *L0_L0bar_pi_pT1_vs_pi_pT2_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0bar_pi_pT1_vs_pi_pT2_cuts_hist");

  TH2F *L0_L0bar_p_eta1_vs_p_eta2_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0bar_p_eta1_vs_p_eta2_cuts_hist");
  TH2F *L0_L0bar_pi_eta1_vs_pi_eta2_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0bar_pi_eta1_vs_pi_eta2_cuts_hist");

  TH2F *L0_L0bar_p_phi1_vs_p_phi2_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0bar_p_phi1_vs_p_phi2_cuts_hist");
  TH2F *L0_L0bar_pi_phi1_vs_pi_phi2_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0bar_pi_phi1_vs_pi_phi2_cuts_hist");


  TH2F *L0_L0_p_pT1_vs_p_pT2_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0_p_pT1_vs_p_pT2_cuts_hist");
  TH2F *L0_L0_pi_pT1_vs_pi_pT2_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0_pi_pT1_vs_pi_pT2_cuts_hist");

  TH2F *L0_L0_p_eta1_vs_p_eta2_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0_p_eta1_vs_p_eta2_cuts_hist");
  TH2F *L0_L0_pi_eta1_vs_pi_eta2_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0_pi_eta1_vs_pi_eta2_cuts_hist");

  TH2F *L0_L0_p_phi1_vs_p_phi2_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0_p_phi1_vs_p_phi2_cuts_hist");
  TH2F *L0_L0_pi_phi1_vs_pi_phi2_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0_pi_phi1_vs_pi_phi2_cuts_hist");


  TH2F *L0bar_L0bar_p_pT1_vs_p_pT2_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0bar_L0bar_p_pT1_vs_p_pT2_cuts_hist");
  TH2F *L0bar_L0bar_pi_pT1_vs_pi_pT2_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0bar_L0bar_pi_pT1_vs_pi_pT2_cuts_hist");

  TH2F *L0bar_L0bar_p_eta1_vs_p_eta2_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0bar_L0bar_p_eta1_vs_p_eta2_cuts_hist");
  TH2F *L0bar_L0bar_pi_eta1_vs_pi_eta2_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0bar_L0bar_pi_eta1_vs_pi_eta2_cuts_hist");

  TH2F *L0bar_L0bar_p_phi1_vs_p_phi2_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0bar_L0bar_p_phi1_vs_p_phi2_cuts_hist");
  TH2F *L0bar_L0bar_pi_phi1_vs_pi_phi2_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0bar_L0bar_pi_phi1_vs_pi_phi2_cuts_hist");


  //compare SE and ME partner kinematics
  //L-bar and daoughters

  TH2F *L0_L0bar_Lbar_eta1_SE_vs_eta2_ME_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0bar_Lbar_eta1_SE_vs_eta2_ME_cuts_hist");
  TH2F *L0_L0bar_Lbar_phi1_SE_vs_phi2_ME_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0bar_Lbar_phi1_SE_vs_phi2_ME_cuts_hist");
  TH2F *L0_L0bar_Lbar_pT1_SE_vs_pT2_ME_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0bar_Lbar_pT1_SE_vs_pT2_ME_cuts_hist");

  //compare daughter kinemtics of SE and ME partner
  TH2F *L0_L0bar_pbar_eta1_SE_vs_eta2_ME_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0bar_pbar_eta1_SE_vs_eta2_ME_cuts_hist");
  TH2F *L0_L0bar_pbar_phi1_SE_vs_phi2_ME_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0bar_pbar_phi1_SE_vs_phi2_ME_cuts_hist");
  TH2F *L0_L0bar_pbar_pT1_SE_vs_pT2_ME_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0bar_pbar_pT1_SE_vs_pT2_ME_cuts_hist");

  TH2F *L0_L0bar_pibar_eta1_SE_vs_eta2_ME_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0bar_pibar_eta1_SE_vs_eta2_ME_cuts_hist");
  TH2F *L0_L0bar_pibar_phi1_SE_vs_phi2_ME_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0bar_pibar_phi1_SE_vs_phi2_ME_cuts_hist");
  TH2F *L0_L0bar_pibar_pT1_SE_vs_pT2_ME_cuts_hist = (TH2F*)inFile_PYTHIA->Get("L0_L0bar_pibar_pT1_SE_vs_pT2_ME_cuts_hist");


  //_______________________________________________________________________________________________________________________________________________________________________________________________

  //data histograms

  //L and Lbar pair kinematics

  TH2F *L0_L0bar_pT1_vs_pT2_US_hist = (TH2F*)inFile_data->Get("L0_L0bar_pT1_vs_pT2_US_hist");
  TH2F *L0_L0bar_eta1_vs_eta2_US_hist = (TH2F*)inFile_data->Get("L0_L0bar_eta1_vs_eta2_US_hist");
  TH2F *L0_L0bar_phi1_vs_phi2_US_hist = (TH2F*)inFile_data->Get("L0_L0bar_phi1_vs_phi2_US_hist");

  TH2F *L0_L0_pT1_vs_pT2_US_hist = (TH2F*)inFile_data->Get("L0_L0_pT1_vs_pT2_US_hist");
  TH2F *L0_L0_eta1_vs_eta2_US_hist = (TH2F*)inFile_data->Get("L0_L0_eta1_vs_eta2_US_hist");
  TH2F *L0_L0_phi1_vs_phi2_US_hist = (TH2F*)inFile_data->Get("L0_L0_phi1_vs_phi2_US_hist");

  TH2F *L0bar_L0bar_pT1_vs_pT2_US_hist = (TH2F*)inFile_data->Get("L0bar_L0bar_pT1_vs_pT2_US_hist");
  TH2F *L0bar_L0bar_eta1_vs_eta2_US_hist = (TH2F*)inFile_data->Get("L0bar_L0bar_eta1_vs_eta2_US_hist");
  TH2F *L0bar_L0bar_phi1_vs_phi2_US_hist = (TH2F*)inFile_data->Get("L0bar_L0bar_phi1_vs_phi2_US_hist");


  TH2F *L0_L0bar_pT1_vs_pT2_US_LS_hist = (TH2F*)inFile_data->Get("L0_L0bar_pT1_vs_pT2_US_LS_hist");
  TH2F *L0_L0bar_eta1_vs_eta2_US_LS_hist = (TH2F*)inFile_data->Get("L0_L0bar_eta1_vs_eta2_US_LS_hist");
  TH2F *L0_L0bar_phi1_vs_phi2_US_LS_hist = (TH2F*)inFile_data->Get("L0_L0bar_phi1_vs_phi2_US_LS_hist");

  TH2F *L0_L0_pT1_vs_pT2_US_LS_hist = (TH2F*)inFile_data->Get("L0_L0_pT1_vs_pT2_US_LS_hist");
  TH2F *L0_L0_eta1_vs_eta2_US_LS_hist = (TH2F*)inFile_data->Get("L0_L0_eta1_vs_eta2_US_LS_hist");
  TH2F *L0_L0_phi1_vs_phi2_US_LS_hist = (TH2F*)inFile_data->Get("L0_L0_phi1_vs_phi2_US_LS_hist");

  TH2F *L0bar_L0bar_pT1_vs_pT2_US_LS_hist = (TH2F*)inFile_data->Get("L0bar_L0bar_pT1_vs_pT2_US_LS_hist");
  TH2F *L0bar_L0bar_eta1_vs_eta2_US_LS_hist = (TH2F*)inFile_data->Get("L0bar_L0bar_eta1_vs_eta2_US_LS_hist");
  TH2F *L0bar_L0bar_phi1_vs_phi2_US_LS_hist = (TH2F*)inFile_data->Get("L0bar_L0bar_phi1_vs_phi2_US_LS_hist");

  //--------------------------------------

  //daughter kinematics

  TH2F *L0_L0bar_p1_pT1_vs_p2_pT2_US_hist = (TH2F*)inFile_data->Get("L0_L0bar_p1_pT1_vs_p2_pT2_US_hist");
  TH2F *L0_L0bar_pi1_pT1_vs_pi2_pT2_US_hist = (TH2F*)inFile_data->Get("L0_L0bar_pi1_pT1_vs_pi2_pT2_US_hist");

  TH2F *L0_L0bar_p1_pT1_vs_p2_pT2_US_LS_hist = (TH2F*)inFile_data->Get("L0_L0bar_p1_pT1_vs_p2_pT2_US_LS_hist");
  TH2F *L0_L0bar_pi1_pT1_vs_pi2_pT2_US_LS_hist = (TH2F*)inFile_data->Get("L0_L0bar_pi1_pT1_vs_pi2_pT2_US_LS_hist");


  TH2F *L0_L0bar_p1_eta1_vs_p2_eta2_US_hist = (TH2F*)inFile_data->Get("L0_L0bar_p1_eta1_vs_p2_eta2_US_hist");
  TH2F *L0_L0bar_pi1_eta1_vs_pi2_eta2_US_hist = (TH2F*)inFile_data->Get("L0_L0bar_pi1_eta1_vs_pi2_eta2_US_hist");

  TH2F *L0_L0bar_p1_eta1_vs_p2_eta2_US_LS_hist = (TH2F*)inFile_data->Get("L0_L0bar_p1_eta1_vs_p2_eta2_US_LS_hist");
  TH2F *L0_L0bar_pi1_eta1_vs_pi2_eta2_US_LS_hist = (TH2F*)inFile_data->Get("L0_L0bar_pi1_eta1_vs_pi2_eta2_US_LS_hist");


  TH2F *L0_L0bar_p1_phi1_vs_p2_phi2_US_hist = (TH2F*)inFile_data->Get("L0_L0bar_p1_phi1_vs_p2_phi2_US_hist");
  TH2F *L0_L0bar_pi1_phi1_vs_pi2_phi2_US_hist = (TH2F*)inFile_data->Get("L0_L0bar_pi1_phi1_vs_pi2_phi2_US_hist");

  TH2F *L0_L0bar_p1_phi1_vs_p2_phi2_US_LS_hist = (TH2F*)inFile_data->Get("L0_L0bar_p1_phi1_vs_p2_phi2_US_LS_hist");
  TH2F *L0_L0bar_pi1_phi1_vs_pi2_phi2_US_LS_hist = (TH2F*)inFile_data->Get("L0_L0bar_pi1_phi1_vs_pi2_phi2_US_LS_hist");

  //-------------------------------

  TH2F *L0_L0_p1_pT1_vs_p2_pT2_US_hist = (TH2F*)inFile_data->Get("L0_L0_p1_pT1_vs_p2_pT2_US_hist");
  TH2F *L0_L0_pi1_pT1_vs_pi2_pT2_US_hist = (TH2F*)inFile_data->Get("L0_L0_pi1_pT1_vs_pi2_pT2_US_hist");

  TH2F *L0_L0_p1_pT1_vs_p2_pT2_US_LS_hist = (TH2F*)inFile_data->Get("L0_L0_p1_pT1_vs_p2_pT2_US_LS_hist");
  TH2F *L0_L0_pi1_pT1_vs_pi2_pT2_US_LS_hist = (TH2F*)inFile_data->Get("L0_L0_pi1_pT1_vs_pi2_pT2_US_LS_hist");


  TH2F *L0_L0_p1_eta1_vs_p2_eta2_US_hist = (TH2F*)inFile_data->Get("L0_L0_p1_eta1_vs_p2_eta2_US_hist");
  TH2F *L0_L0_pi1_eta1_vs_pi2_eta2_US_hist = (TH2F*)inFile_data->Get("L0_L0_pi1_eta1_vs_pi2_eta2_US_hist");

  TH2F *L0_L0_p1_eta1_vs_p2_eta2_US_LS_hist = (TH2F*)inFile_data->Get("L0_L0_p1_eta1_vs_p2_eta2_US_LS_hist");
  TH2F *L0_L0_pi1_eta1_vs_pi2_eta2_US_LS_hist = (TH2F*)inFile_data->Get("L0_L0_pi1_eta1_vs_pi2_eta2_US_LS_hist");


  TH2F *L0_L0_p1_phi1_vs_p2_phi2_US_hist = (TH2F*)inFile_data->Get("L0_L0_p1_phi1_vs_p2_phi2_US_hist");
  TH2F *L0_L0_pi1_phi1_vs_pi2_phi2_US_hist = (TH2F*)inFile_data->Get("L0_L0_pi1_phi1_vs_pi2_phi2_US_hist");

  TH2F *L0_L0_p1_phi1_vs_p2_phi2_US_LS_hist = (TH2F*)inFile_data->Get("L0_L0_p1_phi1_vs_p2_phi2_US_LS_hist");
  TH2F *L0_L0_pi1_phi1_vs_pi2_phi2_US_LS_hist = (TH2F*)inFile_data->Get("L0_L0_pi1_phi1_vs_pi2_phi2_US_LS_hist");

  //-------------------------------

  TH2F *L0bar_L0bar_p1_pT1_vs_p2_pT2_US_hist = (TH2F*)inFile_data->Get("L0bar_L0bar_p1_pT1_vs_p2_pT2_US_hist");
  TH2F *L0bar_L0bar_pi1_pT1_vs_pi2_pT2_US_hist = (TH2F*)inFile_data->Get("L0bar_L0bar_pi1_pT1_vs_pi2_pT2_US_hist");

  TH2F *L0bar_L0bar_p1_pT1_vs_p2_pT2_US_LS_hist = (TH2F*)inFile_data->Get("L0bar_L0bar_p1_pT1_vs_p2_pT2_US_LS_hist");
  TH2F *L0bar_L0bar_pi1_pT1_vs_pi2_pT2_US_LS_hist = (TH2F*)inFile_data->Get("L0bar_L0bar_pi1_pT1_vs_pi2_pT2_US_LS_hist");


  TH2F *L0bar_L0bar_p1_eta1_vs_p2_eta2_US_hist = (TH2F*)inFile_data->Get("L0bar_L0bar_p1_eta1_vs_p2_eta2_US_hist");
  TH2F *L0bar_L0bar_pi1_eta1_vs_pi2_eta2_US_hist = (TH2F*)inFile_data->Get("L0bar_L0bar_pi1_eta1_vs_pi2_eta2_US_hist");

  TH2F *L0bar_L0bar_p1_eta1_vs_p2_eta2_US_LS_hist = (TH2F*)inFile_data->Get("L0bar_L0bar_p1_eta1_vs_p2_eta2_US_LS_hist");
  TH2F *L0bar_L0bar_pi1_eta1_vs_pi2_eta2_US_LS_hist = (TH2F*)inFile_data->Get("L0bar_L0bar_pi1_eta1_vs_pi2_eta2_US_LS_hist");


  TH2F *L0bar_L0bar_p1_phi1_vs_p2_phi2_US_hist = (TH2F*)inFile_data->Get("L0bar_L0bar_p1_phi1_vs_p2_phi2_US_hist");
  TH2F *L0bar_L0bar_pi1_phi1_vs_pi2_phi2_US_hist = (TH2F*)inFile_data->Get("L0bar_L0bar_pi1_phi1_vs_pi2_phi2_US_hist");

  TH2F *L0bar_L0bar_p1_phi1_vs_p2_phi2_US_LS_hist = (TH2F*)inFile_data->Get("L0bar_L0bar_p1_phi1_vs_p2_phi2_US_LS_hist");
  TH2F *L0bar_L0bar_pi1_phi1_vs_pi2_phi2_US_LS_hist = (TH2F*)inFile_data->Get("L0bar_L0bar_pi1_phi1_vs_pi2_phi2_US_LS_hist");

  //_______________________________________________________________________________________________________________________________________________________________________________________________

  //embedding histograms

  TH2F *L0_L0bar_p1_pT1_vs_p2_pT2_RC_US_hist = (TH2F*)inFile_embed->Get("L0_L0bar_p1_pT1_vs_p2_pT2_RC_US_hist");
  TH2F *L0_L0bar_pi1_pT1_vs_pi2_pT2_RC_US_hist = (TH2F*)inFile_embed->Get("L0_L0bar_pi1_pT1_vs_pi2_pT2_RC_US_hist");

  TH2F *L0_L0bar_p1_eta1_vs_p2_eta2_RC_US_hist = (TH2F*)inFile_embed->Get("L0_L0bar_p1_eta1_vs_p2_eta2_RC_US_hist");
  TH2F *L0_L0bar_pi1_eta1_vs_pi2_eta2_RC_US_hist = (TH2F*)inFile_embed->Get("L0_L0bar_pi1_eta1_vs_pi2_eta2_RC_US_hist");

  TH2F *L0_L0bar_p1_phi1_vs_p2_phi2_RC_US_hist = (TH2F*)inFile_embed->Get("L0_L0bar_p1_phi1_vs_p2_phi2_RC_US_hist");
  TH2F *L0_L0bar_pi1_phi1_vs_pi2_phi2_RC_US_hist = (TH2F*)inFile_embed->Get("L0_L0bar_pi1_phi1_vs_pi2_phi2_RC_US_hist");


  TH2F *L0_L0bar_p1_pT1_vs_p2_pT2_RC_US_LS_hist = (TH2F*)inFile_embed->Get("L0_L0bar_p1_pT1_vs_p2_pT2_RC_US_LS_hist");
  TH2F *L0_L0bar_pi1_pT1_vs_pi2_pT2_RC_US_LS_hist = (TH2F*)inFile_embed->Get("L0_L0bar_pi1_pT1_vs_pi2_pT2_RC_US_LS_hist");

  TH2F *L0_L0bar_p1_eta1_vs_p2_eta2_RC_US_LS_hist = (TH2F*)inFile_embed->Get("L0_L0bar_p1_eta1_vs_p2_eta2_RC_US_LS_hist");
  TH2F *L0_L0bar_pi1_eta1_vs_pi2_eta2_RC_US_LS_hist = (TH2F*)inFile_embed->Get("L0_L0bar_pi1_eta1_vs_pi2_eta2_RC_US_LS_hist");

  TH2F *L0_L0bar_p1_phi1_vs_p2_phi2_RC_US_LS_hist = (TH2F*)inFile_embed->Get("L0_L0bar_p1_phi1_vs_p2_phi2_RC_US_LS_hist");
  TH2F *L0_L0bar_pi1_phi1_vs_pi2_phi2_RC_US_LS_hist = (TH2F*)inFile_embed->Get("L0_L0bar_pi1_phi1_vs_pi2_phi2_RC_US_LS_hist");

  //----------------------------------------------------------------------------------

  TH2F *L0_L0_p1_pT1_vs_p2_pT2_RC_US_hist = (TH2F*)inFile_embed->Get("L0_L0_p1_pT1_vs_p2_pT2_RC_US_hist");
  TH2F *L0_L0_pi1_pT1_vs_pi2_pT2_RC_US_hist = (TH2F*)inFile_embed->Get("L0_L0_pi1_pT1_vs_pi2_pT2_RC_US_hist");

  TH2F *L0_L0_p1_eta1_vs_p2_eta2_RC_US_hist = (TH2F*)inFile_embed->Get("L0_L0_p1_eta1_vs_p2_eta2_RC_US_hist");
  TH2F *L0_L0_pi1_eta1_vs_pi2_eta2_RC_US_hist = (TH2F*)inFile_embed->Get("L0_L0_pi1_eta1_vs_pi2_eta2_RC_US_hist");

  TH2F *L0_L0_p1_phi1_vs_p2_phi2_RC_US_hist = (TH2F*)inFile_embed->Get("L0_L0_p1_phi1_vs_p2_phi2_RC_US_hist");
  TH2F *L0_L0_pi1_phi1_vs_pi2_phi2_RC_US_hist = (TH2F*)inFile_embed->Get("L0_L0_pi1_phi1_vs_pi2_phi2_RC_US_hist");


  TH2F *L0_L0_p1_pT1_vs_p2_pT2_RC_US_LS_hist = (TH2F*)inFile_embed->Get("L0_L0_p1_pT1_vs_p2_pT2_RC_US_LS_hist");
  TH2F *L0_L0_pi1_pT1_vs_pi2_pT2_RC_US_LS_hist = (TH2F*)inFile_embed->Get("L0_L0_pi1_pT1_vs_pi2_pT2_RC_US_LS_hist");

  TH2F *L0_L0_p1_eta1_vs_p2_eta2_RC_US_LS_hist = (TH2F*)inFile_embed->Get("L0_L0_p1_eta1_vs_p2_eta2_RC_US_LS_hist");
  TH2F *L0_L0_pi1_eta1_vs_pi2_eta2_RC_US_LS_hist = (TH2F*)inFile_embed->Get("L0_L0_pi1_eta1_vs_pi2_eta2_RC_US_LS_hist");

  TH2F *L0_L0_p1_phi1_vs_p2_phi2_RC_US_LS_hist = (TH2F*)inFile_embed->Get("L0_L0_p1_phi1_vs_p2_phi2_RC_US_LS_hist");
  TH2F *L0_L0_pi1_phi1_vs_pi2_phi2_RC_US_LS_hist = (TH2F*)inFile_embed->Get("L0_L0_pi1_phi1_vs_pi2_phi2_RC_US_LS_hist");


  //_______________________________________________________________________________________

  //draw histograms

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  //Minv
/*
  TCanvas *L_Minv_canvas = new TCanvas("L_Minv_canvas", "L_Minv_canvas", 1200, 1000);
  L_Minv_canvas->cd();

  L_Minv_US->GetXaxis()->SetTitle("M_{inv} (GeV/c^{2})");
  L_Minv_US->GetXaxis()->CenterTitle();
  L_Minv_US->Add(L_Minv_LS, -1);
  L_Minv_US->Scale(1./L_Minv_US->Integral());
  L_Minv_US->Draw("hist");

  L_Minv_hist_cuts->Scale(1./L_Minv_hist_cuts->Integral());
  L_Minv_hist_cuts->SetLineColor(kRed);
  L_Minv_hist_cuts->Draw("hist same");

  TLegend *L_Minv_leg = new TLegend(0.2, 0.4, 0.5, 0.6);
  L_Minv_leg->AddEntry(L_Minv_US, "Data (US-LS)");
  L_Minv_leg->AddEntry(L_Minv_hist_cuts, "PYTHIA after cuts");
  L_Minv_leg->SetBorderSize(0);
  L_Minv_leg->Draw("same");

  L_Minv_canvas->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/kine_QA/L_Minv.png");

*/
  //_______________________________________________________________________________________________________________________________

  //daughter pT

  //data

  TCanvas *L_Lbar_pi_pT_2D_data_can = new TCanvas("L_Lbar_pi_pT_2D_data_can", "L_Lbar_pi_pT_2D_data_can", 1200, 1000);
  L_Lbar_pi_pT_2D_data_can->cd();

  L0_L0bar_pi1_pT1_vs_pi2_pT2_US_hist->Add(L0_L0bar_pi1_pT1_vs_pi2_pT2_US_LS_hist, -1); //subtract combinatorial background
  L0_L0bar_pi1_pT1_vs_pi2_pT2_US_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L0_L0bar_pi1_pT1_vs_pi2_pT2_US_hist->GetXaxis()->CenterTitle();
  L0_L0bar_pi1_pT1_vs_pi2_pT2_US_hist->GetYaxis()->SetTitle("p_{T} GeV/c");
  L0_L0bar_pi1_pT1_vs_pi2_pT2_US_hist->GetYaxis()->CenterTitle();
  L0_L0bar_pi1_pT1_vs_pi2_pT2_US_hist->Draw("colz");

  L_Lbar_pi_pT_2D_data_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L_Lbar_pi_pT_2D_data.png");


  TCanvas *L_Lbar_pi_pT_data_can = new TCanvas("L_Lbar_pi_pT_data_can", "L_Lbar_pi_pT_data_can", 1200, 1000);
  L_Lbar_pi_pT_data_can->cd();

  //L0_L0bar_pi1_pT1_vs_pi2_pT2_US_hist->Add(L0_L0bar_pi1_pT1_vs_pi2_pT2_US_LS_hist, -1); //subtract combinatorial background

  TH1D *L_Lbar_pi_pT_data_hist = L0_L0bar_pi1_pT1_vs_pi2_pT2_US_hist->ProjectionX("L_Lbar_pi_pT_pr");
  L_Lbar_pi_pT_data_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L_Lbar_pi_pT_data_hist->GetXaxis()->CenterTitle();
  L_Lbar_pi_pT_data_hist->SetLineColor(1);
  L_Lbar_pi_pT_data_hist->Draw("hist");

  TH1D *L_Lbar_piBar_pT_data_hist = L0_L0bar_pi1_pT1_vs_pi2_pT2_US_hist->ProjectionY("L_Lbar_piBar_pT_pr");
  L_Lbar_piBar_pT_data_hist->SetLineColor(kRed);
  L_Lbar_piBar_pT_data_hist->Draw("hist same");

  TLegend *L_Lbar_pi_pT_leg = new TLegend(0.7, 0.7, 0.89, 0.89);
  L_Lbar_pi_pT_leg->AddEntry(L_Lbar_pi_pT_data_hist, "#pi^{-}");
  L_Lbar_pi_pT_leg->AddEntry(L_Lbar_piBar_pT_data_hist, "#pi^{+}");
  L_Lbar_pi_pT_leg->SetBorderSize(0);
  L_Lbar_pi_pT_leg->SetFillColorAlpha(0, 0.01);
  L_Lbar_pi_pT_leg->Draw("same");

  L_Lbar_pi_pT_data_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L_Lbar_pi_pT_data.png");

  //--------------------------------------------------

  TCanvas *L_Lbar_p_pT_2D_data_can = new TCanvas("L_Lbar_p_pT_2D_data_can", "L_Lbar_p_pT_2D_data_can", 1200, 1000);
  L_Lbar_p_pT_2D_data_can->cd();

  L0_L0bar_p1_pT1_vs_p2_pT2_US_hist->Add(L0_L0bar_p1_pT1_vs_p2_pT2_US_LS_hist, -1); //subtract combinatorial background
  L0_L0bar_p1_pT1_vs_p2_pT2_US_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L0_L0bar_p1_pT1_vs_p2_pT2_US_hist->GetXaxis()->CenterTitle();
  L0_L0bar_p1_pT1_vs_p2_pT2_US_hist->GetYaxis()->SetTitle("p_{T} GeV/c");
  L0_L0bar_p1_pT1_vs_p2_pT2_US_hist->GetYaxis()->CenterTitle();
  L0_L0bar_p1_pT1_vs_p2_pT2_US_hist->Draw("colz");

  L_Lbar_p_pT_2D_data_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L_Lbar_p_pT_2D_data.png");


  TCanvas *L_Lbar_p_pT_data_can = new TCanvas("L_Lbar_p_pT_data_can", "L_Lbar_p_pT_data_can", 1200, 1000);
  L_Lbar_p_pT_data_can->cd();

  //L0_L0bar_p1_pT1_vs_p2_pT2_US_hist->Add(L0_L0bar_p1_pT1_vs_p2_pT2_US_LS_hist, -1); //subtract combinatorial background

  TH1D *L_Lbar_p_pT_data_hist = L0_L0bar_p1_pT1_vs_p2_pT2_US_hist->ProjectionX("L_Lbar_p_pT_pr");
  L_Lbar_p_pT_data_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L_Lbar_p_pT_data_hist->GetXaxis()->CenterTitle();
  L_Lbar_p_pT_data_hist->SetLineColor(1);
  L_Lbar_p_pT_data_hist->Draw("hist");

  TH1D *L_Lbar_pBar_pT_data_hist = L0_L0bar_p1_pT1_vs_p2_pT2_US_hist->ProjectionY("L_Lbar_pBar_pT_pr");
  L_Lbar_pBar_pT_data_hist->SetLineColor(kRed);
  L_Lbar_pBar_pT_data_hist->Draw("hist same");

  TLegend *L_Lbar_p_pT_leg = new TLegend(0.7, 0.7, 0.89, 0.89);
  L_Lbar_p_pT_leg->AddEntry(L_Lbar_p_pT_data_hist, "p");
  L_Lbar_p_pT_leg->AddEntry(L_Lbar_pBar_pT_data_hist, "#bar{p}");
  L_Lbar_p_pT_leg->SetBorderSize(0);
  L_Lbar_p_pT_leg->SetFillColorAlpha(0, 0.01);
  L_Lbar_p_pT_leg->Draw("same");

  L_Lbar_p_pT_data_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L_Lbar_p_pT_data.png");


  //---------------------------------------------------------------------------------------------------------------------------------------

  TCanvas *L_Lbar_pi_eta_2D_data_can = new TCanvas("L_Lbar_pi_eta_2D_data_can", "L_Lbar_pi_eta_2D_data_can", 1200, 1000);
  L_Lbar_pi_eta_2D_data_can->cd();

  L0_L0bar_pi1_eta1_vs_pi2_eta2_US_hist->Add(L0_L0bar_pi1_eta1_vs_pi2_eta2_US_LS_hist, -1); //subtract combinatorial background
  L0_L0bar_pi1_eta1_vs_pi2_eta2_US_hist->GetXaxis()->SetTitle("#eta");
  L0_L0bar_pi1_eta1_vs_pi2_eta2_US_hist->GetXaxis()->CenterTitle();
  L0_L0bar_pi1_eta1_vs_pi2_eta2_US_hist->GetYaxis()->SetTitle("#eta");
  L0_L0bar_pi1_eta1_vs_pi2_eta2_US_hist->GetYaxis()->CenterTitle();
  L0_L0bar_pi1_eta1_vs_pi2_eta2_US_hist->Draw("colz");

  L_Lbar_pi_eta_2D_data_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L_Lbar_pi_eta_2D_data.png");


  TCanvas *L_Lbar_pi_eta_can = new TCanvas("L_Lbar_pi_eta_can", "L_Lbar_pi_eta_can", 1200, 1000);
  L_Lbar_pi_eta_can->cd();

  //L0_L0bar_pi1_eta1_vs_pi2_eta2_US_hist->Add(L0_L0bar_pi1_eta1_vs_pi2_eta2_US_LS_hist, -1); //subtract combinatorial background

  TH1D *L_Lbar_pi_eta_hist = L0_L0bar_pi1_eta1_vs_pi2_eta2_US_hist->ProjectionX("L_Lbar_pi_eta_pr");
  L_Lbar_pi_eta_hist->GetXaxis()->SetTitle("#eta");
  L_Lbar_pi_eta_hist->GetXaxis()->CenterTitle();
  L_Lbar_pi_eta_hist->SetLineColor(1);
  L_Lbar_pi_eta_hist->SetMinimum(0);
  L_Lbar_pi_eta_hist->Draw("hist");

  TH1D *L_Lbar_piBar_eta_hist = L0_L0bar_pi1_eta1_vs_pi2_eta2_US_hist->ProjectionY("L_Lbar_piBar_eta_pr");
  L_Lbar_piBar_eta_hist->SetLineColor(kRed);
  L_Lbar_piBar_eta_hist->Draw("hist same");

  TLegend *L_Lbar_pi_eta_leg = new TLegend(0.7, 0.7, 0.89, 0.89);
  L_Lbar_pi_eta_leg->AddEntry(L_Lbar_pi_eta_hist, "#pi^{-}");
  L_Lbar_pi_eta_leg->AddEntry(L_Lbar_piBar_eta_hist, "#pi^{+}");
  L_Lbar_pi_eta_leg->SetBorderSize(0);
  L_Lbar_pi_eta_leg->SetFillColorAlpha(0, 0.01);
  L_Lbar_pi_eta_leg->Draw("same");

  L_Lbar_pi_eta_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L_Lbar_pi_eta_data.png");

  //--------------------------------------------------

  TCanvas *L_Lbar_p_eta_2D_data_can = new TCanvas("L_Lbar_p_eta_2D_data_can", "L_Lbar_p_eta_2D_data_can", 1200, 1000);
  L_Lbar_p_eta_2D_data_can->cd();

  L0_L0bar_p1_eta1_vs_p2_eta2_US_hist->Add(L0_L0bar_p1_eta1_vs_p2_eta2_US_LS_hist, -1); //subtract combinatorial background
  L0_L0bar_p1_eta1_vs_p2_eta2_US_hist->GetXaxis()->SetTitle("#eta");
  L0_L0bar_p1_eta1_vs_p2_eta2_US_hist->GetXaxis()->CenterTitle();
  L0_L0bar_p1_eta1_vs_p2_eta2_US_hist->GetYaxis()->SetTitle("#eta");
  L0_L0bar_p1_eta1_vs_p2_eta2_US_hist->GetYaxis()->CenterTitle();
  L0_L0bar_p1_eta1_vs_p2_eta2_US_hist->Draw("colz");

  L_Lbar_p_eta_2D_data_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L_Lbar_p_eta_2D_data.png");


  TCanvas *L_Lbar_p_eta_can = new TCanvas("L_Lbar_p_eta_can", "L_Lbar_p_eta_can", 1200, 1000);
  L_Lbar_p_eta_can->cd();

  L0_L0bar_p1_eta1_vs_p2_eta2_US_hist->Add(L0_L0bar_p1_eta1_vs_p2_eta2_US_LS_hist, -1); //subtract combinatorial background

  TH1D *L_Lbar_p_eta_hist = L0_L0bar_p1_eta1_vs_p2_eta2_US_hist->ProjectionX("L_Lbar_p_eta_pr");
  L_Lbar_p_eta_hist->GetXaxis()->SetTitle("#eta");
  L_Lbar_p_eta_hist->GetXaxis()->CenterTitle();
  L_Lbar_p_eta_hist->SetLineColor(1);
  L_Lbar_p_eta_hist->SetMinimum(0);
  L_Lbar_p_eta_hist->Draw("hist");

  TH1D *L_Lbar_pBar_eta_hist = L0_L0bar_p1_eta1_vs_p2_eta2_US_hist->ProjectionY("L_Lbar_pBar_eta_pr");
  L_Lbar_pBar_eta_hist->SetLineColor(kRed);
  L_Lbar_pBar_eta_hist->Draw("hist same");

  TLegend *L_Lbar_p_eta_leg = new TLegend(0.7, 0.7, 0.89, 0.89);
  L_Lbar_p_eta_leg->AddEntry(L_Lbar_p_eta_hist, "p");
  L_Lbar_p_eta_leg->AddEntry(L_Lbar_pBar_eta_hist, "#bar{p}");
  L_Lbar_p_eta_leg->SetBorderSize(0);
  L_Lbar_p_eta_leg->SetFillColorAlpha(0, 0.01);
  L_Lbar_p_eta_leg->Draw("same");

  L_Lbar_p_eta_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L_Lbar_p_eta_data.png");

  //----------------------------------------------------------------------------------------------------------------------------------

  TCanvas *L_Lbar_pi_phi_2D_data_can = new TCanvas("L_Lbar_pi_phi_2D_data_can", "L_Lbar_pi_phi_2D_data_can", 1200, 1000);
  L_Lbar_pi_phi_2D_data_can->cd();

  L0_L0bar_pi1_phi1_vs_pi2_phi2_US_hist->Add(L0_L0bar_pi1_phi1_vs_pi2_phi2_US_LS_hist, -1); //subtract combinatorial background
  L0_L0bar_pi1_phi1_vs_pi2_phi2_US_hist->GetXaxis()->SetTitle("#phi");
  L0_L0bar_pi1_phi1_vs_pi2_phi2_US_hist->GetXaxis()->CenterTitle();
  L0_L0bar_pi1_phi1_vs_pi2_phi2_US_hist->GetYaxis()->SetTitle("#phi");
  L0_L0bar_pi1_phi1_vs_pi2_phi2_US_hist->GetYaxis()->CenterTitle();
  L0_L0bar_pi1_phi1_vs_pi2_phi2_US_hist->Draw("colz");

  L_Lbar_pi_phi_2D_data_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L_Lbar_pi_phi_2D_data.png");


  TCanvas *L_Lbar_pi_phi_can = new TCanvas("L_Lbar_pi_phi_can", "L_Lbar_pi_phi_can", 1200, 1000);
  L_Lbar_pi_phi_can->cd();

  //L0_L0bar_pi1_phi1_vs_pi2_phi2_US_hist->Add(L0_L0bar_pi1_phi1_vs_pi2_phi2_US_LS_hist, -1); //subtract combinatorial background

  TH1D *L_Lbar_pi_phi_hist = L0_L0bar_pi1_phi1_vs_pi2_phi2_US_hist->ProjectionX("L_Lbar_pi_phi_pr");
  L_Lbar_pi_phi_hist->GetXaxis()->SetTitle("#phi");
  L_Lbar_pi_phi_hist->GetXaxis()->CenterTitle();
  L_Lbar_pi_phi_hist->SetLineColor(1);
  L_Lbar_pi_phi_hist->SetMinimum(0);
  L_Lbar_pi_phi_hist->Draw("hist");

  TH1D *L_Lbar_piBar_phi_hist = L0_L0bar_pi1_phi1_vs_pi2_phi2_US_hist->ProjectionY("L_Lbar_piBar_phi_pr");
  L_Lbar_piBar_phi_hist->SetLineColor(kRed);
  L_Lbar_piBar_phi_hist->Draw("hist same");

  TLegend *L_Lbar_pi_phi_leg = new TLegend(0.7, 0.7, 0.89, 0.89);
  L_Lbar_pi_phi_leg->AddEntry(L_Lbar_pi_phi_hist, "#pi^{-}");
  L_Lbar_pi_phi_leg->AddEntry(L_Lbar_piBar_phi_hist, "#pi^{+}");
  L_Lbar_pi_phi_leg->SetBorderSize(0);
  L_Lbar_pi_phi_leg->SetFillColorAlpha(0, 0.01);
  L_Lbar_pi_phi_leg->Draw("same");

  L_Lbar_pi_phi_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L_Lbar_pi_phi_data.png");

  //--------------------------------------------------

  TCanvas *L_Lbar_p_phi_2D_data_can = new TCanvas("L_Lbar_p_phi_2D_data_can", "L_Lbar_p_phi_2D_data_can", 1200, 1000);
  L_Lbar_p_phi_2D_data_can->cd();

  L0_L0bar_p1_phi1_vs_p2_phi2_US_hist->Add(L0_L0bar_p1_phi1_vs_p2_phi2_US_LS_hist, -1); //subtract combinatorial background
  L0_L0bar_p1_phi1_vs_p2_phi2_US_hist->GetXaxis()->SetTitle("#phi");
  L0_L0bar_p1_phi1_vs_p2_phi2_US_hist->GetXaxis()->CenterTitle();
  L0_L0bar_p1_phi1_vs_p2_phi2_US_hist->GetYaxis()->SetTitle("#phi");
  L0_L0bar_p1_phi1_vs_p2_phi2_US_hist->GetYaxis()->CenterTitle();
  L0_L0bar_p1_phi1_vs_p2_phi2_US_hist->Draw("colz");

  L_Lbar_p_phi_2D_data_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L_Lbar_p_phi_2D_data.png");



  TCanvas *L_Lbar_p_phi_can = new TCanvas("L_Lbar_p_phi_can", "L_Lbar_p_phi_can", 1200, 1000);
  L_Lbar_p_phi_can->cd();

  //L0_L0bar_p1_phi1_vs_p2_phi2_US_hist->Add(L0_L0bar_p1_phi1_vs_p2_phi2_US_LS_hist, -1); //subtract combinatorial background

  TH1D *L_Lbar_p_phi_hist = L0_L0bar_p1_phi1_vs_p2_phi2_US_hist->ProjectionX("L_Lbar_p_phi_pr");
  L_Lbar_p_phi_hist->GetXaxis()->SetTitle("#phi");
  L_Lbar_p_phi_hist->GetXaxis()->CenterTitle();
  L_Lbar_p_phi_hist->SetLineColor(1);
  L_Lbar_p_phi_hist->SetMinimum(0);
  L_Lbar_p_phi_hist->Draw("hist");

  TH1D *L_Lbar_pBar_phi_hist = L0_L0bar_p1_phi1_vs_p2_phi2_US_hist->ProjectionY("L_Lbar_pBar_phi_pr");
  L_Lbar_pBar_phi_hist->SetLineColor(kRed);
  L_Lbar_pBar_phi_hist->Draw("hist same");

  TLegend *L_Lbar_p_phi_leg = new TLegend(0.7, 0.7, 0.89, 0.89);
  L_Lbar_p_phi_leg->AddEntry(L_Lbar_p_phi_hist, "p");
  L_Lbar_p_phi_leg->AddEntry(L_Lbar_pBar_phi_hist, "#bar{p}");
  L_Lbar_p_phi_leg->SetBorderSize(0);
  L_Lbar_p_phi_leg->SetFillColorAlpha(0, 0.01);
  L_Lbar_p_phi_leg->Draw("same");

  L_Lbar_p_phi_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L_Lbar_p_phi_data.png");

  //---------------------------------------------------------------------------------------------------------------------------------------------------------------

  TCanvas *L_L_pi_pT_2D_data_can = new TCanvas("L_L_pi_pT_2D_data_can", "L_L_pi_pT_2D_data_can", 1200, 1000);
  L_L_pi_pT_2D_data_can->cd();

  L0_L0_pi1_pT1_vs_pi2_pT2_US_hist->Add(L0_L0_pi1_pT1_vs_pi2_pT2_US_LS_hist, -1); //subtract combinatorial background
  L0_L0_pi1_pT1_vs_pi2_pT2_US_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L0_L0_pi1_pT1_vs_pi2_pT2_US_hist->GetXaxis()->CenterTitle();
  L0_L0_pi1_pT1_vs_pi2_pT2_US_hist->GetYaxis()->SetTitle("p_{T} GeV/c");
  L0_L0_pi1_pT1_vs_pi2_pT2_US_hist->GetYaxis()->CenterTitle();
  L0_L0_pi1_pT1_vs_pi2_pT2_US_hist->Draw("colz");

  L_L_pi_pT_2D_data_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_L/L_L_pi_pT_2D_data.png");


  TCanvas *L_L_pi_pT_data_can = new TCanvas("L_L_pi_pT_data_can", "L_L_pi_pT_data_can", 1200, 1000);
  L_L_pi_pT_data_can->cd();

  //L0_L0_pi1_pT1_vs_pi2_pT2_US_hist->Add(L0_L0_pi1_pT1_vs_pi2_pT2_US_LS_hist, -1); //subtract combinatorial background

  TH1D *L_L_pi_pT_data_hist = L0_L0_pi1_pT1_vs_pi2_pT2_US_hist->ProjectionX("L_L_pi_pT_pr");
  L_L_pi_pT_data_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L_L_pi_pT_data_hist->GetXaxis()->CenterTitle();
  L_L_pi_pT_data_hist->SetLineColor(1);
  L_L_pi_pT_data_hist->Draw("hist");

  TH1D *L_L_piBar_pT_data_hist = L0_L0_pi1_pT1_vs_pi2_pT2_US_hist->ProjectionY("L_L_piBar_pT_pr");
  L_L_piBar_pT_data_hist->SetLineColor(kRed);
  L_L_piBar_pT_data_hist->Draw("hist same");

  TLegend *L_L_pi_pT_leg = new TLegend(0.7, 0.7, 0.89, 0.89);
  L_L_pi_pT_leg->AddEntry(L_L_pi_pT_data_hist, "#pi^{-}_{1}");
  L_L_pi_pT_leg->AddEntry(L_L_piBar_pT_data_hist, "#pi^{-}_{2}");
  L_L_pi_pT_leg->SetBorderSize(0);
  L_L_pi_pT_leg->SetFillColorAlpha(0, 0.01);
  L_L_pi_pT_leg->Draw("same");

  L_L_pi_pT_data_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_L/L_L_pi_pT_data.png");

  //--------------------------------------------------

  TCanvas *L_L_p_pT_2D_data_can = new TCanvas("L_L_p_pT_2D_data_can", "L_L_p_pT_2D_data_can", 1200, 1000);
  L_L_p_pT_2D_data_can->cd();

  L0_L0_p1_pT1_vs_p2_pT2_US_hist->Add(L0_L0_p1_pT1_vs_p2_pT2_US_LS_hist, -1); //subtract combinatorial background
  L0_L0_p1_pT1_vs_p2_pT2_US_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L0_L0_p1_pT1_vs_p2_pT2_US_hist->GetXaxis()->CenterTitle();
  L0_L0_p1_pT1_vs_p2_pT2_US_hist->GetYaxis()->SetTitle("p_{T} GeV/c");
  L0_L0_p1_pT1_vs_p2_pT2_US_hist->GetYaxis()->CenterTitle();
  L0_L0_p1_pT1_vs_p2_pT2_US_hist->Draw("colz");

  L_L_p_pT_2D_data_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_L/L_L_p_pT_2D_data.png");


  TCanvas *L_L_p_pT_data_can = new TCanvas("L_L_p_pT_data_can", "L_L_p_pT_data_can", 1200, 1000);
  L_L_p_pT_data_can->cd();

  //L0_L0_p1_pT1_vs_p2_pT2_US_hist->Add(L0_L0_p1_pT1_vs_p2_pT2_US_LS_hist, -1); //subtract combinatorial background

  TH1D *L_L_p_pT_data_hist = L0_L0_p1_pT1_vs_p2_pT2_US_hist->ProjectionX("L_L_p_pT_pr");
  L_L_p_pT_data_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L_L_p_pT_data_hist->GetXaxis()->CenterTitle();
  L_L_p_pT_data_hist->SetLineColor(1);
  L_L_p_pT_data_hist->Draw("hist");

  TH1D *L_L_pBar_pT_data_hist = L0_L0_p1_pT1_vs_p2_pT2_US_hist->ProjectionY("L_L_pBar_pT_pr");
  L_L_pBar_pT_data_hist->SetLineColor(kRed);
  L_L_pBar_pT_data_hist->Draw("hist same");

  TLegend *L_L_p_pT_leg = new TLegend(0.7, 0.7, 0.89, 0.89);
  L_L_p_pT_leg->AddEntry(L_L_p_pT_data_hist, "p_{1}");
  L_L_p_pT_leg->AddEntry(L_L_pBar_pT_data_hist, "p_{2}");
  L_L_p_pT_leg->SetBorderSize(0);
  L_L_p_pT_leg->SetFillColorAlpha(0, 0.01);
  L_L_p_pT_leg->Draw("same");

  L_L_p_pT_data_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_L/L_L_p_pT_data.png");

  //----------------------------------------------------------------------------------------------------------------------------------

  TCanvas *L_L_pi_eta_2D_data_can = new TCanvas("L_L_pi_eta_2D_data_can", "L_L_pi_eta_2D_data_can", 1200, 1000);
  L_L_pi_eta_2D_data_can->cd();

  L0_L0_pi1_eta1_vs_pi2_eta2_US_hist->Add(L0_L0_pi1_eta1_vs_pi2_eta2_US_LS_hist, -1); //subtract combinatorial background
  L0_L0_pi1_eta1_vs_pi2_eta2_US_hist->GetXaxis()->SetTitle("#eta");
  L0_L0_pi1_eta1_vs_pi2_eta2_US_hist->GetXaxis()->CenterTitle();
  L0_L0_pi1_eta1_vs_pi2_eta2_US_hist->GetYaxis()->SetTitle("#eta");
  L0_L0_pi1_eta1_vs_pi2_eta2_US_hist->GetYaxis()->CenterTitle();
  L0_L0_pi1_eta1_vs_pi2_eta2_US_hist->Draw("colz");

  L_L_pi_eta_2D_data_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_L/L_L_pi_eta_2D_data.png");


  TCanvas *L_L_pi_eta_can = new TCanvas("L_L_pi_eta_can", "L_L_pi_eta_can", 1200, 1000);
  L_L_pi_eta_can->cd();

  //L0_L0_pi1_eta1_vs_pi2_eta2_US_hist->Add(L0_L0_pi1_eta1_vs_pi2_eta2_US_LS_hist, -1); //subtract combinatorial background

  TH1D *L_L_pi_eta_hist = L0_L0_pi1_eta1_vs_pi2_eta2_US_hist->ProjectionX("L_L_pi_eta_pr");
  L_L_pi_eta_hist->GetXaxis()->SetTitle("#eta");
  L_L_pi_eta_hist->GetXaxis()->CenterTitle();
  L_L_pi_eta_hist->SetLineColor(1);
  L_L_pi_eta_hist->SetMinimum(0);
  L_L_pi_eta_hist->Draw("hist");

  TH1D *L_L_piBar_eta_hist = L0_L0_pi1_eta1_vs_pi2_eta2_US_hist->ProjectionY("L_L_piBar_eta_pr");
  L_L_piBar_eta_hist->SetLineColor(kRed);
  L_L_piBar_eta_hist->Draw("hist same");

  TLegend *L_L_pi_eta_leg = new TLegend(0.7, 0.7, 0.89, 0.89);
  L_L_pi_eta_leg->AddEntry(L_L_pi_eta_hist, "#pi^{-}_{1}");
  L_L_pi_eta_leg->AddEntry(L_L_piBar_eta_hist, "#pi^{-}_{2}");
  L_L_pi_eta_leg->SetBorderSize(0);
  L_L_pi_eta_leg->SetFillColorAlpha(0, 0.01);
  L_L_pi_eta_leg->Draw("same");

  L_L_pi_eta_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_L/L_L_pi_eta_data.png");

  //--------------------------------------------------

  TCanvas *L_L_p_eta_2D_data_can = new TCanvas("L_L_p_eta_2D_data_can", "L_L_p_eta_2D_data_can", 1200, 1000);
  L_L_p_eta_2D_data_can->cd();

  L0_L0_p1_eta1_vs_p2_eta2_US_hist->Add(L0_L0_p1_eta1_vs_p2_eta2_US_LS_hist, -1); //subtract combinatorial background
  L0_L0_p1_eta1_vs_p2_eta2_US_hist->GetXaxis()->SetTitle("#eta");
  L0_L0_p1_eta1_vs_p2_eta2_US_hist->GetXaxis()->CenterTitle();
  L0_L0_p1_eta1_vs_p2_eta2_US_hist->GetYaxis()->SetTitle("#eta");
  L0_L0_p1_eta1_vs_p2_eta2_US_hist->GetYaxis()->CenterTitle();
  L0_L0_p1_eta1_vs_p2_eta2_US_hist->Draw("colz");

  L_L_p_eta_2D_data_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_L/L_L_p_eta_2D_data.png");


  TCanvas *L_L_p_eta_can = new TCanvas("L_L_p_eta_can", "L_L_p_eta_can", 1200, 1000);
  L_L_p_eta_can->cd();

  //L0_L0_p1_eta1_vs_p2_eta2_US_hist->Add(L0_L0_p1_eta1_vs_p2_eta2_US_LS_hist, -1); //subtract combinatorial background

  TH1D *L_L_p_eta_hist = L0_L0_p1_eta1_vs_p2_eta2_US_hist->ProjectionX("L_L_p_eta_pr");
  L_L_p_eta_hist->GetXaxis()->SetTitle("#eta");
  L_L_p_eta_hist->GetXaxis()->CenterTitle();
  L_L_p_eta_hist->SetLineColor(1);
  L_L_p_eta_hist->SetMinimum(0);
  L_L_p_eta_hist->Draw("hist");

  TH1D *L_L_pBar_eta_hist = L0_L0_p1_eta1_vs_p2_eta2_US_hist->ProjectionY("L_L_pBar_eta_pr");
  L_L_pBar_eta_hist->SetLineColor(kRed);
  L_L_pBar_eta_hist->Draw("hist same");

  TLegend *L_L_p_eta_leg = new TLegend(0.7, 0.7, 0.89, 0.89);
  L_L_p_eta_leg->AddEntry(L_L_p_eta_hist, "p_{1}");
  L_L_p_eta_leg->AddEntry(L_L_pBar_eta_hist, "p_{2}");
  L_L_p_eta_leg->SetBorderSize(0);
  L_L_p_eta_leg->SetFillColorAlpha(0, 0.01);
  L_L_p_eta_leg->Draw("same");

  L_L_p_eta_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_L/L_L_p_eta_data.png");

  //----------------------------------------------------------------------------------------------------------------------------------

  TCanvas *L_L_pi_phi_2D_data_can = new TCanvas("L_L_pi_phi_2D_data_can", "L_L_pi_phi_2D_data_can", 1200, 1000);
  L_L_pi_phi_2D_data_can->cd();

  L0_L0_pi1_phi1_vs_pi2_phi2_US_hist->Add(L0_L0_pi1_phi1_vs_pi2_phi2_US_LS_hist, -1); //subtract combinatorial background
  L0_L0_pi1_phi1_vs_pi2_phi2_US_hist->GetXaxis()->SetTitle("#phi");
  L0_L0_pi1_phi1_vs_pi2_phi2_US_hist->GetXaxis()->CenterTitle();
  L0_L0_pi1_phi1_vs_pi2_phi2_US_hist->GetYaxis()->SetTitle("#phi");
  L0_L0_pi1_phi1_vs_pi2_phi2_US_hist->GetYaxis()->CenterTitle();
  L0_L0_pi1_phi1_vs_pi2_phi2_US_hist->Draw("colz");

  L_L_pi_phi_2D_data_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_L/L_L_pi_phi_2D_data.png");


  TCanvas *L_L_pi_phi_can = new TCanvas("L_L_pi_phi_can", "L_L_pi_phi_can", 1200, 1000);
  L_L_pi_phi_can->cd();

  //L0_L0_pi1_phi1_vs_pi2_phi2_US_hist->Add(L0_L0_pi1_phi1_vs_pi2_phi2_US_LS_hist, -1); //subtract combinatorial background

  TH1D *L_L_pi_phi_hist = L0_L0_pi1_phi1_vs_pi2_phi2_US_hist->ProjectionX("L_L_pi_phi_pr");
  L_L_pi_phi_hist->GetXaxis()->SetTitle("#phi");
  L_L_pi_phi_hist->GetXaxis()->CenterTitle();
  L_L_pi_phi_hist->SetLineColor(1);
  L_L_pi_phi_hist->SetMinimum(0);
  L_L_pi_phi_hist->Draw("hist");

  TH1D *L_L_piBar_phi_hist = L0_L0_pi1_phi1_vs_pi2_phi2_US_hist->ProjectionY("L_L_piBar_phi_pr");
  L_L_piBar_phi_hist->SetLineColor(kRed);
  L_L_piBar_phi_hist->Draw("hist same");

  TLegend *L_L_pi_phi_leg = new TLegend(0.7, 0.7, 0.89, 0.89);
  L_L_pi_phi_leg->AddEntry(L_L_pi_phi_hist, "#pi^{-}_{1}");
  L_L_pi_phi_leg->AddEntry(L_L_piBar_phi_hist, "#pi^{-}_{2}");
  L_L_pi_phi_leg->SetBorderSize(0);
  L_L_pi_phi_leg->SetFillColorAlpha(0, 0.01);
  L_L_pi_phi_leg->Draw("same");

  L_L_pi_phi_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_L/L_L_pi_phi_data.png");

  //--------------------------------------------------

  TCanvas *L_L_p_phi_2D_data_can = new TCanvas("L_L_p_phi_2D_data_can", "L_L_p_phi_2D_data_can", 1200, 1000);
  L_L_p_phi_2D_data_can->cd();

  L0_L0_p1_phi1_vs_p2_phi2_US_hist->Add(L0_L0_p1_phi1_vs_p2_phi2_US_LS_hist, -1); //subtract combinatorial background
  L0_L0_p1_phi1_vs_p2_phi2_US_hist->GetXaxis()->SetTitle("#phi");
  L0_L0_p1_phi1_vs_p2_phi2_US_hist->GetXaxis()->CenterTitle();
  L0_L0_p1_phi1_vs_p2_phi2_US_hist->GetYaxis()->SetTitle("#phi");
  L0_L0_p1_phi1_vs_p2_phi2_US_hist->GetYaxis()->CenterTitle();
  L0_L0_p1_phi1_vs_p2_phi2_US_hist->Draw("colz");

  L_L_p_phi_2D_data_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_L/L_L_p_phi_2D_data.png");


  TCanvas *L_L_p_phi_can = new TCanvas("L_L_p_phi_can", "L_L_p_phi_can", 1200, 1000);
  L_L_p_phi_can->cd();

  //L0_L0_p1_phi1_vs_p2_phi2_US_hist->Add(L0_L0_p1_phi1_vs_p2_phi2_US_LS_hist, -1); //subtract combinatorial background

  TH1D *L_L_p_phi_hist = L0_L0_p1_phi1_vs_p2_phi2_US_hist->ProjectionX("L_L_p_phi_pr");
  L_L_p_phi_hist->GetXaxis()->SetTitle("#phi");
  L_L_p_phi_hist->GetXaxis()->CenterTitle();
  L_L_p_phi_hist->SetLineColor(1);
  L_L_p_phi_hist->SetMinimum(0);
  L_L_p_phi_hist->Draw("hist");

  TH1D *L_L_pBar_phi_hist = L0_L0_p1_phi1_vs_p2_phi2_US_hist->ProjectionY("L_L_pBar_phi_pr");
  L_L_pBar_phi_hist->SetLineColor(kRed);
  L_L_pBar_phi_hist->Draw("hist same");

  TLegend *L_L_p_phi_leg = new TLegend(0.7, 0.7, 0.89, 0.89);
  L_L_p_phi_leg->AddEntry(L_L_p_phi_hist, "p_{1}");
  L_L_p_phi_leg->AddEntry(L_L_pBar_phi_hist, "p_{2}");
  L_L_p_phi_leg->SetBorderSize(0);
  L_L_p_phi_leg->SetFillColorAlpha(0, 0.01);
  L_L_p_phi_leg->Draw("same");

  L_L_p_phi_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_L/L_L_p_phi_data.png");

  //----------------------------------------------------------------------------------------------------------------------------------

  TCanvas *Lbar_Lbar_pi_pT_data_can = new TCanvas("Lbar_Lbar_pi_pT_data_can", "Lbar_Lbar_pi_pT_data_can", 1200, 1000);
  Lbar_Lbar_pi_pT_data_can->cd();

  L0bar_L0bar_pi1_pT1_vs_pi2_pT2_US_hist->Add(L0bar_L0bar_pi1_pT1_vs_pi2_pT2_US_LS_hist, -1); //subtract combinatorial background

  TH1D *Lbar_Lbar_pi_pT_data_hist = L0bar_L0bar_pi1_pT1_vs_pi2_pT2_US_hist->ProjectionX("Lbar_Lbar_pi_pT_pr");
  Lbar_Lbar_pi_pT_data_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  Lbar_Lbar_pi_pT_data_hist->GetXaxis()->CenterTitle();
  Lbar_Lbar_pi_pT_data_hist->SetLineColor(1);
  Lbar_Lbar_pi_pT_data_hist->Draw("hist");

  TH1D *Lbar_Lbar_piBar_pT_data_hist = L0bar_L0bar_pi1_pT1_vs_pi2_pT2_US_hist->ProjectionY("Lbar_Lbar_piBar_pT_pr");
  Lbar_Lbar_piBar_pT_data_hist->SetLineColor(kRed);
  Lbar_Lbar_piBar_pT_data_hist->Draw("hist same");

  TLegend *Lbar_Lbar_pi_pT_leg = new TLegend(0.7, 0.7, 0.89, 0.89);
  Lbar_Lbar_pi_pT_leg->AddEntry(Lbar_Lbar_pi_pT_data_hist, "#pi^{+}_{1}");
  Lbar_Lbar_pi_pT_leg->AddEntry(Lbar_Lbar_piBar_pT_data_hist, "#pi^{+}_{2}");
  Lbar_Lbar_pi_pT_leg->SetBorderSize(0);
  Lbar_Lbar_pi_pT_leg->SetFillColorAlpha(0, 0.01);
  Lbar_Lbar_pi_pT_leg->Draw("same");

  Lbar_Lbar_pi_pT_data_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/Lbar_Lbar/Lbar_Lbar_pi_pT_data.png");

  //--------------------------------------------------


  TCanvas *Lbar_Lbar_p_pT_data_can = new TCanvas("Lbar_Lbar_p_pT_data_can", "Lbar_Lbar_p_pT_data_can", 1200, 1000);
  Lbar_Lbar_p_pT_data_can->cd();

  L0bar_L0bar_p1_pT1_vs_p2_pT2_US_hist->Add(L0bar_L0bar_p1_pT1_vs_p2_pT2_US_LS_hist, -1); //subtract combinatorial background

  TH1D *Lbar_Lbar_p_pT_data_hist = L0bar_L0bar_p1_pT1_vs_p2_pT2_US_hist->ProjectionX("Lbar_Lbar_p_pT_pr");
  Lbar_Lbar_p_pT_data_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  Lbar_Lbar_p_pT_data_hist->GetXaxis()->CenterTitle();
  Lbar_Lbar_p_pT_data_hist->SetLineColor(1);
  Lbar_Lbar_p_pT_data_hist->Draw("hist");

  TH1D *Lbar_Lbar_pBar_pT_data_hist = L0bar_L0bar_p1_pT1_vs_p2_pT2_US_hist->ProjectionY("Lbar_Lbar_pBar_pT_pr");
  Lbar_Lbar_pBar_pT_data_hist->SetLineColor(kRed);
  Lbar_Lbar_pBar_pT_data_hist->Draw("hist same");

  TLegend *Lbar_Lbar_p_pT_leg = new TLegend(0.7, 0.7, 0.89, 0.89);
  Lbar_Lbar_p_pT_leg->AddEntry(Lbar_Lbar_p_pT_data_hist, "#bar{p}_{1}");
  Lbar_Lbar_p_pT_leg->AddEntry(Lbar_Lbar_pBar_pT_data_hist, "#bar{p}_{2}");
  Lbar_Lbar_p_pT_leg->SetBorderSize(0);
  Lbar_Lbar_p_pT_leg->SetFillColorAlpha(0, 0.01);
  Lbar_Lbar_p_pT_leg->Draw("same");

  Lbar_Lbar_p_pT_data_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/Lbar_Lbar/Lbar_Lbar_p_pT_data.png");

  //-------------------------------------------------------------------------------------------------------------------------------------------------------------------

  TCanvas *Lbar_Lbar_pi_eta_can = new TCanvas("Lbar_Lbar_pi_eta_can", "Lbar_Lbar_pi_eta_can", 1200, 1000);
  Lbar_Lbar_pi_eta_can->cd();

  L0bar_L0bar_pi1_eta1_vs_pi2_eta2_US_hist->Add(L0bar_L0bar_pi1_eta1_vs_pi2_eta2_US_LS_hist, -1); //subtract combinatorial background

  TH1D *Lbar_Lbar_pi_eta_hist = L0bar_L0bar_pi1_eta1_vs_pi2_eta2_US_hist->ProjectionX("Lbar_Lbar_pi_eta_pr");
  Lbar_Lbar_pi_eta_hist->GetXaxis()->SetTitle("#eta");
  Lbar_Lbar_pi_eta_hist->GetXaxis()->CenterTitle();
  Lbar_Lbar_pi_eta_hist->SetLineColor(1);
  Lbar_Lbar_pi_eta_hist->SetMinimum(0);
  Lbar_Lbar_pi_eta_hist->Draw("hist");

  TH1D *Lbar_Lbar_piBar_eta_hist = L0bar_L0bar_pi1_eta1_vs_pi2_eta2_US_hist->ProjectionY("Lbar_Lbar_piBar_eta_pr");
  Lbar_Lbar_piBar_eta_hist->SetLineColor(kRed);
  Lbar_Lbar_piBar_eta_hist->Draw("hist same");

  TLegend *Lbar_Lbar_pi_eta_leg = new TLegend(0.7, 0.7, 0.89, 0.89);
  Lbar_Lbar_pi_eta_leg->AddEntry(Lbar_Lbar_pi_eta_hist, "#pi^{+}_{1}");
  Lbar_Lbar_pi_eta_leg->AddEntry(Lbar_Lbar_piBar_eta_hist, "#pi^{+}_{2}");
  Lbar_Lbar_pi_eta_leg->SetBorderSize(0);
  Lbar_Lbar_pi_eta_leg->SetFillColorAlpha(0, 0.01);
  Lbar_Lbar_pi_eta_leg->Draw("same");

  Lbar_Lbar_pi_eta_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/Lbar_Lbar/Lbar_Lbar_pi_eta_data.png");

  //--------------------------------------------------


  TCanvas *Lbar_Lbar_p_eta_can = new TCanvas("Lbar_Lbar_p_eta_can", "Lbar_Lbar_p_eta_can", 1200, 1000);
  Lbar_Lbar_p_eta_can->cd();

  L0bar_L0bar_p1_eta1_vs_p2_eta2_US_hist->Add(L0bar_L0bar_p1_eta1_vs_p2_eta2_US_LS_hist, -1); //subtract combinatorial background

  TH1D *Lbar_Lbar_p_eta_hist = L0bar_L0bar_p1_eta1_vs_p2_eta2_US_hist->ProjectionX("Lbar_Lbar_p_eta_pr");
  Lbar_Lbar_p_eta_hist->GetXaxis()->SetTitle("#eta");
  Lbar_Lbar_p_eta_hist->GetXaxis()->CenterTitle();
  Lbar_Lbar_p_eta_hist->SetLineColor(1);
  Lbar_Lbar_p_eta_hist->SetMinimum(0);
  Lbar_Lbar_p_eta_hist->Draw("hist");

  TH1D *Lbar_Lbar_pBar_eta_hist = L0bar_L0bar_p1_eta1_vs_p2_eta2_US_hist->ProjectionY("Lbar_Lbar_pBar_eta_pr");
  Lbar_Lbar_pBar_eta_hist->SetLineColor(kRed);
  Lbar_Lbar_pBar_eta_hist->Draw("hist same");

  TLegend *Lbar_Lbar_p_eta_leg = new TLegend(0.7, 0.7, 0.89, 0.89);
  Lbar_Lbar_p_eta_leg->AddEntry(Lbar_Lbar_p_eta_hist, "#bar{p}_{1}");
  Lbar_Lbar_p_eta_leg->AddEntry(Lbar_Lbar_pBar_eta_hist, "#bar{p}_{2}");
  Lbar_Lbar_p_eta_leg->SetBorderSize(0);
  Lbar_Lbar_p_eta_leg->SetFillColorAlpha(0, 0.01);
  Lbar_Lbar_p_eta_leg->Draw("same");

  Lbar_Lbar_p_eta_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/Lbar_Lbar/Lbar_Lbar_p_eta_data.png");

  //----------------------------------------------------------------------------------------------------------------------------------

  TCanvas *Lbar_Lbar_pi_phi_can = new TCanvas("Lbar_Lbar_pi_phi_can", "Lbar_Lbar_pi_phi_can", 1200, 1000);
  Lbar_Lbar_pi_phi_can->cd();

  L0bar_L0bar_pi1_phi1_vs_pi2_phi2_US_hist->Add(L0bar_L0bar_pi1_phi1_vs_pi2_phi2_US_LS_hist, -1); //subtract combinatorial background

  TH1D *Lbar_Lbar_pi_phi_hist = L0bar_L0bar_pi1_phi1_vs_pi2_phi2_US_hist->ProjectionX("Lbar_Lbar_pi_phi_pr");
  Lbar_Lbar_pi_phi_hist->GetXaxis()->SetTitle("#phi");
  Lbar_Lbar_pi_phi_hist->GetXaxis()->CenterTitle();
  Lbar_Lbar_pi_phi_hist->SetLineColor(1);
  Lbar_Lbar_pi_phi_hist->SetMinimum(0);
  Lbar_Lbar_pi_phi_hist->Draw("hist");

  TH1D *Lbar_Lbar_piBar_phi_hist = L0bar_L0bar_pi1_phi1_vs_pi2_phi2_US_hist->ProjectionY("Lbar_Lbar_piBar_phi_pr");
  Lbar_Lbar_piBar_phi_hist->SetLineColor(kRed);
  Lbar_Lbar_piBar_phi_hist->Draw("hist same");

  TLegend *Lbar_Lbar_pi_phi_leg = new TLegend(0.7, 0.7, 0.89, 0.89);
  Lbar_Lbar_pi_phi_leg->AddEntry(Lbar_Lbar_pi_phi_hist, "#pi^{+}_{1}");
  Lbar_Lbar_pi_phi_leg->AddEntry(Lbar_Lbar_piBar_phi_hist, "#pi^{+}_{2}");
  Lbar_Lbar_pi_phi_leg->SetBorderSize(0);
  Lbar_Lbar_pi_phi_leg->SetFillColorAlpha(0, 0.01);
  Lbar_Lbar_pi_phi_leg->Draw("same");

  Lbar_Lbar_pi_phi_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/Lbar_Lbar/Lbar_Lbar_pi_phi_data.png");

  //--------------------------------------------------


  TCanvas *Lbar_Lbar_p_phi_can = new TCanvas("Lbar_Lbar_p_phi_can", "Lbar_Lbar_p_phi_can", 1200, 1000);
  Lbar_Lbar_p_phi_can->cd();

  L0bar_L0bar_p1_phi1_vs_p2_phi2_US_hist->Add(L0bar_L0bar_p1_phi1_vs_p2_phi2_US_LS_hist, -1); //subtract combinatorial background

  TH1D *Lbar_Lbar_p_phi_hist = L0bar_L0bar_p1_phi1_vs_p2_phi2_US_hist->ProjectionX("Lbar_Lbar_p_phi_pr");
  Lbar_Lbar_p_phi_hist->GetXaxis()->SetTitle("#phi");
  Lbar_Lbar_p_phi_hist->GetXaxis()->CenterTitle();
  Lbar_Lbar_p_phi_hist->SetLineColor(1);
  Lbar_Lbar_p_phi_hist->SetMinimum(0);
  Lbar_Lbar_p_phi_hist->Draw("hist");

  TH1D *Lbar_Lbar_pBar_phi_hist = L0bar_L0bar_p1_phi1_vs_p2_phi2_US_hist->ProjectionY("Lbar_Lbar_pBar_phi_pr");
  Lbar_Lbar_pBar_phi_hist->SetLineColor(kRed);
  Lbar_Lbar_pBar_phi_hist->Draw("hist same");

  TLegend *Lbar_Lbar_p_phi_leg = new TLegend(0.7, 0.7, 0.89, 0.89);
  Lbar_Lbar_p_phi_leg->AddEntry(Lbar_Lbar_p_phi_hist, "#bar{p}_{1}");
  Lbar_Lbar_p_phi_leg->AddEntry(Lbar_Lbar_pBar_phi_hist, "#bar{p}_{2}");
  Lbar_Lbar_p_phi_leg->SetBorderSize(0);
  Lbar_Lbar_p_phi_leg->SetFillColorAlpha(0, 0.01);
  Lbar_Lbar_p_phi_leg->Draw("same");

  Lbar_Lbar_p_phi_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/Lbar_Lbar/Lbar_Lbar_p_phi_data.png");

  //________________________________________________________________________________________________________________________________________________________________________________________________

  //PYTHIA

  TCanvas *L_Lbar_pi_pT_2D_PYTHIA_can = new TCanvas("L_Lbar_pi_pT_2D_PYTHIA_can", "L_Lbar_pi_pT_2D_PYTHIA_can", 1200, 1000);
  L_Lbar_pi_pT_2D_PYTHIA_can->cd();

  L0_L0bar_pi_pT1_vs_pi_pT2_cuts_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L0_L0bar_pi_pT1_vs_pi_pT2_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0bar_pi_pT1_vs_pi_pT2_cuts_hist->GetYaxis()->SetTitle("p_{T} GeV/c");
  L0_L0bar_pi_pT1_vs_pi_pT2_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0bar_pi_pT1_vs_pi_pT2_cuts_hist->Draw("colz");

  L_Lbar_pi_pT_2D_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L_Lbar_pi_pT_2D_PYTHIA.png");


  TCanvas *L_Lbar_pi_pT_PYTHIA_can = new TCanvas("L_Lbar_pi_pT_PYTHIA_can", "L_Lbar_pi_pT_PYTHIA_can", 1200, 1000);
  L_Lbar_pi_pT_PYTHIA_can->cd();

  TH1D *L_Lbar_pi_pT_PYTHIA_hist = L0_L0bar_pi_pT1_vs_pi_pT2_cuts_hist->ProjectionX("L_Lbar_pi_pT_PYTHIA_pr");
  L_Lbar_pi_pT_PYTHIA_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L_Lbar_pi_pT_PYTHIA_hist->GetXaxis()->CenterTitle();
  L_Lbar_pi_pT_PYTHIA_hist->SetLineColor(1);
  L_Lbar_pi_pT_PYTHIA_hist->Draw("hist");

  TH1D *L_Lbar_piBar_pT_PYTHIA_hist = L0_L0bar_pi_pT1_vs_pi_pT2_cuts_hist->ProjectionY("L_Lbar_piBar_pT_PYTHIA_pr");
  L_Lbar_piBar_pT_PYTHIA_hist->SetLineColor(kRed);
  L_Lbar_piBar_pT_PYTHIA_hist->Draw("hist same");

  L_Lbar_pi_pT_leg->Draw("same");

  L_Lbar_pi_pT_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L_Lbar_pi_pT_PYTHIA.png");

  //--------------------------------------------------

  TCanvas *L_Lbar_p_pT_2D_PYTHIA_can = new TCanvas("L_Lbar_p_pT_2D_PYTHIA_can", "L_Lbar_p_pT_2D_PYTHIA_can", 1200, 1000);
  L_Lbar_p_pT_2D_PYTHIA_can->cd();

  L0_L0bar_p_pT1_vs_p_pT2_cuts_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L0_L0bar_p_pT1_vs_p_pT2_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0bar_p_pT1_vs_p_pT2_cuts_hist->GetYaxis()->SetTitle("p_{T} GeV/c");
  L0_L0bar_p_pT1_vs_p_pT2_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0bar_p_pT1_vs_p_pT2_cuts_hist->Draw("colz");

  L_Lbar_p_pT_2D_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L_Lbar_p_pT_2D_PYTHIA.png");


  TCanvas *L_Lbar_p_pT_PYTHIA_can = new TCanvas("L_Lbar_p_pT_PYTHIA_can", "L_Lbar_p_pT_PYTHIA_can", 1200, 1000);
  L_Lbar_p_pT_PYTHIA_can->cd();

  TH1D *L_Lbar_p_pT_PYTHIA_hist = L0_L0bar_p_pT1_vs_p_pT2_cuts_hist->ProjectionX("L_Lbar_p_pT_PYTHIA_pr");
  L_Lbar_p_pT_PYTHIA_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L_Lbar_p_pT_PYTHIA_hist->GetXaxis()->CenterTitle();
  L_Lbar_p_pT_PYTHIA_hist->SetLineColor(1);
  L_Lbar_p_pT_PYTHIA_hist->Draw("hist");

  TH1D *L_Lbar_pBar_pT_PYTHIA_hist = L0_L0bar_p_pT1_vs_p_pT2_cuts_hist->ProjectionY("L_Lbar_pBar_pT_PYTHIA_pr");
  L_Lbar_pBar_pT_PYTHIA_hist->SetLineColor(kRed);
  L_Lbar_pBar_pT_PYTHIA_hist->Draw("hist same");

  L_Lbar_p_pT_leg->Draw("same");

  L_Lbar_p_pT_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L_Lbar_p_pT_PYTHIA.png");

  //-------------------------------------------------------------------------------------------------------------------------------------------------

  TCanvas *L_Lbar_pi_eta_2D_PYTHIA_can = new TCanvas("L_Lbar_pi_eta_2D_PYTHIA_can", "L_Lbar_pi_eta_2D_PYTHIA_can", 1200, 1000);
  L_Lbar_pi_eta_2D_PYTHIA_can->cd();

  L0_L0bar_pi_eta1_vs_pi_eta2_cuts_hist->GetXaxis()->SetTitle("#eta");
  L0_L0bar_pi_eta1_vs_pi_eta2_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0bar_pi_eta1_vs_pi_eta2_cuts_hist->GetYaxis()->SetTitle("#eta");
  L0_L0bar_pi_eta1_vs_pi_eta2_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0bar_pi_eta1_vs_pi_eta2_cuts_hist->Draw("colz");

  L_Lbar_pi_eta_2D_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L_Lbar_pi_eta_2D_PYTHIA.png");


  TCanvas *L_Lbar_pi_eta_PYTHIA_can = new TCanvas("L_Lbar_pi_eta_PYTHIA_can", "L_Lbar_pi_eta_PYTHIA_can", 1200, 1000);
  L_Lbar_pi_eta_PYTHIA_can->cd();

  TH1D *L_Lbar_pi_eta_PYTHIA_hist = L0_L0bar_pi_eta1_vs_pi_eta2_cuts_hist->ProjectionX("L_Lbar_pi_eta_PYTHIA_pr");
  L_Lbar_pi_eta_PYTHIA_hist->GetXaxis()->SetTitle("#eta");
  L_Lbar_pi_eta_PYTHIA_hist->GetXaxis()->CenterTitle();
  L_Lbar_pi_eta_PYTHIA_hist->SetLineColor(1);
  L_Lbar_pi_eta_PYTHIA_hist->SetMinimum(0);
  L_Lbar_pi_eta_PYTHIA_hist->Draw("hist");

  TH1D *L_Lbar_piBar_eta_PYTHIA_hist = L0_L0bar_pi_eta1_vs_pi_eta2_cuts_hist->ProjectionY("L_Lbar_piBar_eta_PYTHIA_pr");
  L_Lbar_piBar_eta_PYTHIA_hist->SetLineColor(kRed);
  L_Lbar_piBar_eta_PYTHIA_hist->Draw("hist same");

  L_Lbar_pi_eta_leg->Draw("same");

  L_Lbar_pi_eta_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L_Lbar_pi_eta_PYTHIA.png");

  //--------------------------------------------------

  TCanvas *L_Lbar_p_eta_2D_PYTHIA_can = new TCanvas("L_Lbar_p_eta_2D_PYTHIA_can", "L_Lbar_p_eta_2D_PYTHIA_can", 1200, 1000);
  L_Lbar_p_eta_2D_PYTHIA_can->cd();

  L0_L0bar_p_eta1_vs_p_eta2_cuts_hist->GetXaxis()->SetTitle("#eta");
  L0_L0bar_p_eta1_vs_p_eta2_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0bar_p_eta1_vs_p_eta2_cuts_hist->GetYaxis()->SetTitle("#eta");
  L0_L0bar_p_eta1_vs_p_eta2_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0bar_p_eta1_vs_p_eta2_cuts_hist->Draw("colz");

  L_Lbar_p_eta_2D_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L_Lbar_p_eta_2D_PYTHIA.png");


  TCanvas *L_Lbar_p_eta_PYTHIA_can = new TCanvas("L_Lbar_p_eta_PYTHIA_can", "L_Lbar_p_eta_PYTHIA_can", 1200, 1000);
  L_Lbar_p_eta_PYTHIA_can->cd();

  TH1D *L_Lbar_p_eta_PYTHIA_hist = L0_L0bar_p_eta1_vs_p_eta2_cuts_hist->ProjectionX("L_Lbar_p_eta_PYTHIA_pr");
  L_Lbar_p_eta_PYTHIA_hist->GetXaxis()->SetTitle("#eta");
  L_Lbar_p_eta_PYTHIA_hist->GetXaxis()->CenterTitle();
  L_Lbar_p_eta_PYTHIA_hist->SetLineColor(1);
  L_Lbar_p_eta_PYTHIA_hist->SetMinimum(0);
  L_Lbar_p_eta_PYTHIA_hist->Draw("hist");

  TH1D *L_Lbar_pBar_eta_PYTHIA_hist = L0_L0bar_p_eta1_vs_p_eta2_cuts_hist->ProjectionY("L_Lbar_pBar_eta_PYTHIA_pr");
  L_Lbar_pBar_eta_PYTHIA_hist->SetLineColor(kRed);
  L_Lbar_pBar_eta_PYTHIA_hist->Draw("hist same");

  L_Lbar_p_eta_leg->Draw("same");

  L_Lbar_p_eta_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L_Lbar_p_eta_PYTHIA.png");

  //-------------------------------------------------------------------------------------------------------------------------------------------------

  TCanvas *L_Lbar_pi_phi_2D_PYTHIA_can = new TCanvas("L_Lbar_pi_phi_2D_PYTHIA_can", "L_Lbar_pi_phi_2D_PYTHIA_can", 1200, 1000);
  L_Lbar_pi_phi_2D_PYTHIA_can->cd();

  L0_L0bar_pi_phi1_vs_pi_phi2_cuts_hist->GetXaxis()->SetTitle("#phi");
  L0_L0bar_pi_phi1_vs_pi_phi2_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0bar_pi_phi1_vs_pi_phi2_cuts_hist->GetYaxis()->SetTitle("#phi");
  L0_L0bar_pi_phi1_vs_pi_phi2_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0bar_pi_phi1_vs_pi_phi2_cuts_hist->Draw("colz");

  L_Lbar_pi_phi_2D_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L_Lbar_pi_phi_2D_PYTHIA.png");


  TCanvas *L_Lbar_pi_phi_PYTHIA_can = new TCanvas("L_Lbar_pi_phi_PYTHIA_can", "L_Lbar_pi_phi_PYTHIA_can", 1200, 1000);
  L_Lbar_pi_phi_PYTHIA_can->cd();

  TH1D *L_Lbar_pi_phi_PYTHIA_hist = L0_L0bar_pi_phi1_vs_pi_phi2_cuts_hist->ProjectionX("L_Lbar_pi_phi_PYTHIA_pr");
  L_Lbar_pi_phi_PYTHIA_hist->GetXaxis()->SetTitle("#phi");
  L_Lbar_pi_phi_PYTHIA_hist->GetXaxis()->CenterTitle();
  L_Lbar_pi_phi_PYTHIA_hist->SetLineColor(1);
  L_Lbar_pi_phi_PYTHIA_hist->SetMinimum(0);
  L_Lbar_pi_phi_PYTHIA_hist->Draw("hist");

  TH1D *L_Lbar_piBar_phi_PYTHIA_hist = L0_L0bar_pi_phi1_vs_pi_phi2_cuts_hist->ProjectionY("L_Lbar_piBar_phi_PYTHIA_pr");
  L_Lbar_piBar_phi_PYTHIA_hist->SetLineColor(kRed);
  L_Lbar_piBar_phi_PYTHIA_hist->Draw("hist same");

  L_Lbar_pi_phi_leg->Draw("same");

  L_Lbar_pi_phi_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L_Lbar_pi_phi_PYTHIA.png");

  //--------------------------------------------------

  TCanvas *L_Lbar_p_phi_2D_PYTHIA_can = new TCanvas("L_Lbar_p_phi_2D_PYTHIA_can", "L_Lbar_p_phi_2D_PYTHIA_can", 1200, 1000);
  L_Lbar_p_phi_2D_PYTHIA_can->cd();

  L0_L0bar_p_phi1_vs_p_phi2_cuts_hist->GetXaxis()->SetTitle("#phi");
  L0_L0bar_p_phi1_vs_p_phi2_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0bar_p_phi1_vs_p_phi2_cuts_hist->GetYaxis()->SetTitle("#phi");
  L0_L0bar_p_phi1_vs_p_phi2_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0bar_p_phi1_vs_p_phi2_cuts_hist->Draw("colz");

  L_Lbar_p_phi_2D_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L_Lbar_p_phi_2D_PYTHIA.png");


  TCanvas *L_Lbar_p_phi_PYTHIA_can = new TCanvas("L_Lbar_p_phi_PYTHIA_can", "L_Lbar_p_phi_PYTHIA_can", 1200, 1000);
  L_Lbar_p_phi_PYTHIA_can->cd();

  TH1D *L_Lbar_p_phi_PYTHIA_hist = L0_L0bar_p_phi1_vs_p_phi2_cuts_hist->ProjectionX("L_Lbar_p_phi_PYTHIA_pr");
  L_Lbar_p_phi_PYTHIA_hist->GetXaxis()->SetTitle("#phi");
  L_Lbar_p_phi_PYTHIA_hist->GetXaxis()->CenterTitle();
  L_Lbar_p_phi_PYTHIA_hist->SetLineColor(1);
  L_Lbar_p_phi_PYTHIA_hist->SetMinimum(0);
  L_Lbar_p_phi_PYTHIA_hist->Draw("hist");

  TH1D *L_Lbar_pBar_phi_PYTHIA_hist = L0_L0bar_p_phi1_vs_p_phi2_cuts_hist->ProjectionY("L_Lbar_pBar_phi_PYTHIA_pr");
  L_Lbar_pBar_phi_PYTHIA_hist->SetLineColor(kRed);
  L_Lbar_pBar_phi_PYTHIA_hist->Draw("hist same");

  L_Lbar_p_phi_leg->Draw("same");

  L_Lbar_p_phi_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L_Lbar_p_phi_PYTHIA.png");

  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  TCanvas *L_L_pi_pT_2D_PYTHIA_can = new TCanvas("L_L_pi_pT_2D_PYTHIA_can", "L_L_pi_pT_2D_PYTHIA_can", 1200, 1000);
  L_L_pi_pT_2D_PYTHIA_can->cd();

  L0_L0_pi_pT1_vs_pi_pT2_cuts_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L0_L0_pi_pT1_vs_pi_pT2_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0_pi_pT1_vs_pi_pT2_cuts_hist->GetYaxis()->SetTitle("p_{T} GeV/c");
  L0_L0_pi_pT1_vs_pi_pT2_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0_pi_pT1_vs_pi_pT2_cuts_hist->Draw("colz");

  L_L_pi_pT_2D_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_L/L_L_pi_pT_2D_PYTHIA.png");


  TCanvas *L_L_pi_pT_PYTHIA_can = new TCanvas("L_L_pi_pT_PYTHIA_can", "L_L_pi_pT_PYTHIA_can", 1200, 1000);
  L_L_pi_pT_PYTHIA_can->cd();

  TH1D *L_L_pi_pT_PYTHIA_hist = L0_L0_pi_pT1_vs_pi_pT2_cuts_hist->ProjectionX("L_L_pi_pT_PYTHIA_pr");
  L_L_pi_pT_PYTHIA_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L_L_pi_pT_PYTHIA_hist->GetXaxis()->CenterTitle();
  L_L_pi_pT_PYTHIA_hist->SetLineColor(1);
  L_L_pi_pT_PYTHIA_hist->Draw("hist");

  TH1D *L_L_piBar_pT_PYTHIA_hist = L0_L0_pi_pT1_vs_pi_pT2_cuts_hist->ProjectionY("L_L_piBar_pT_PYTHIA_pr");
  L_L_piBar_pT_PYTHIA_hist->SetLineColor(kRed);
  L_L_piBar_pT_PYTHIA_hist->Draw("hist same");

  L_L_pi_pT_leg->Draw("same");

  L_L_pi_pT_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_L/L_L_pi_pT_PYTHIA.png");

  //--------------------------------------------------

  TCanvas *L_L_p_pT_2D_PYTHIA_can = new TCanvas("L_L_p_pT_2D_PYTHIA_can", "L_L_p_pT_2D_PYTHIA_can", 1200, 1000);
  L_L_p_pT_2D_PYTHIA_can->cd();

  L0_L0_p_pT1_vs_p_pT2_cuts_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L0_L0_p_pT1_vs_p_pT2_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0_p_pT1_vs_p_pT2_cuts_hist->GetYaxis()->SetTitle("p_{T} GeV/c");
  L0_L0_p_pT1_vs_p_pT2_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0_p_pT1_vs_p_pT2_cuts_hist->Draw("colz");

  L_L_p_pT_2D_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_L/L_L_p_pT_2D_PYTHIA.png");


  TCanvas *L_L_p_pT_PYTHIA_can = new TCanvas("L_L_p_pT_PYTHIA_can", "L_L_p_pT_PYTHIA_can", 1200, 1000);
  L_L_p_pT_PYTHIA_can->cd();

  TH1D *L_L_p_pT_PYTHIA_hist = L0_L0_p_pT1_vs_p_pT2_cuts_hist->ProjectionX("L_L_p_pT_PYTHIA_pr");
  L_L_p_pT_PYTHIA_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L_L_p_pT_PYTHIA_hist->GetXaxis()->CenterTitle();
  L_L_p_pT_PYTHIA_hist->SetLineColor(1);
  L_L_p_pT_PYTHIA_hist->Draw("hist");

  TH1D *L_L_pBar_pT_PYTHIA_hist = L0_L0_p_pT1_vs_p_pT2_cuts_hist->ProjectionY("L_L_pBar_pT_PYTHIA_pr");
  L_L_pBar_pT_PYTHIA_hist->SetLineColor(kRed);
  L_L_pBar_pT_PYTHIA_hist->Draw("hist same");

  L_L_p_pT_leg->Draw("same");

  L_L_p_pT_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_L/L_L_p_pT_PYTHIA.png");

   //-------------------------------------------------------------------------------------------------------------------------------------------------

  TCanvas *L_L_pi_eta_2D_PYTHIA_can = new TCanvas("L_L_pi_eta_2D_PYTHIA_can", "L_L_pi_eta_2D_PYTHIA_can", 1200, 1000);
  L_L_pi_eta_2D_PYTHIA_can->cd();

  L0_L0_pi_eta1_vs_pi_eta2_cuts_hist->GetXaxis()->SetTitle("#eta");
  L0_L0_pi_eta1_vs_pi_eta2_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0_pi_eta1_vs_pi_eta2_cuts_hist->GetYaxis()->SetTitle("#eta");
  L0_L0_pi_eta1_vs_pi_eta2_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0_pi_eta1_vs_pi_eta2_cuts_hist->Draw("colz");

  L_L_pi_eta_2D_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_L/L_L_pi_eta_2D_PYTHIA.png");


  TCanvas *L_L_pi_eta_PYTHIA_can = new TCanvas("L_L_pi_eta_PYTHIA_can", "L_L_pi_eta_PYTHIA_can", 1200, 1000);
  L_L_pi_eta_PYTHIA_can->cd();

  TH1D *L_L_pi_eta_PYTHIA_hist = L0_L0_pi_eta1_vs_pi_eta2_cuts_hist->ProjectionX("L_L_pi_eta_PYTHIA_pr");
  L_L_pi_eta_PYTHIA_hist->GetXaxis()->SetTitle("#eta");
  L_L_pi_eta_PYTHIA_hist->GetXaxis()->CenterTitle();
  L_L_pi_eta_PYTHIA_hist->SetLineColor(1);
  L_L_pi_eta_PYTHIA_hist->SetMinimum(0);
  L_L_pi_eta_PYTHIA_hist->Draw("hist");

  TH1D *L_L_piBar_eta_PYTHIA_hist = L0_L0_pi_eta1_vs_pi_eta2_cuts_hist->ProjectionY("L_L_piBar_eta_PYTHIA_pr");
  L_L_piBar_eta_PYTHIA_hist->SetLineColor(kRed);
  L_L_piBar_eta_PYTHIA_hist->Draw("hist same");

  L_L_pi_eta_leg->Draw("same");

  L_L_pi_eta_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_L/L_L_pi_eta_PYTHIA.png");

  //--------------------------------------------------

  TCanvas *L_L_p_eta_2D_PYTHIA_can = new TCanvas("L_L_p_eta_2D_PYTHIA_can", "L_L_p_eta_2D_PYTHIA_can", 1200, 1000);
  L_L_p_eta_2D_PYTHIA_can->cd();

  L0_L0_p_eta1_vs_p_eta2_cuts_hist->GetXaxis()->SetTitle("#eta");
  L0_L0_p_eta1_vs_p_eta2_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0_p_eta1_vs_p_eta2_cuts_hist->GetYaxis()->SetTitle("#eta");
  L0_L0_p_eta1_vs_p_eta2_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0_p_eta1_vs_p_eta2_cuts_hist->Draw("colz");

  L_L_p_eta_2D_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_L/L_L_p_eta_2D_PYTHIA.png");


  TCanvas *L_L_p_eta_PYTHIA_can = new TCanvas("L_L_p_eta_PYTHIA_can", "L_L_p_eta_PYTHIA_can", 1200, 1000);
  L_L_p_eta_PYTHIA_can->cd();

  TH1D *L_L_p_eta_PYTHIA_hist = L0_L0_p_eta1_vs_p_eta2_cuts_hist->ProjectionX("L_L_p_eta_PYTHIA_pr");
  L_L_p_eta_PYTHIA_hist->GetXaxis()->SetTitle("#eta");
  L_L_p_eta_PYTHIA_hist->GetXaxis()->CenterTitle();
  L_L_p_eta_PYTHIA_hist->SetLineColor(1);
  L_L_p_eta_PYTHIA_hist->SetMinimum(0);
  L_L_p_eta_PYTHIA_hist->Draw("hist");

  TH1D *L_L_pBar_eta_PYTHIA_hist = L0_L0_p_eta1_vs_p_eta2_cuts_hist->ProjectionY("L_L_pBar_eta_PYTHIA_pr");
  L_L_pBar_eta_PYTHIA_hist->SetLineColor(kRed);
  L_L_pBar_eta_PYTHIA_hist->Draw("hist same");

  L_L_p_eta_leg->Draw("same");

  L_L_p_eta_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_L/L_L_p_eta_PYTHIA.png");

  //-------------------------------------------------------------------------------------------------------------------------------------------------

  TCanvas *L_L_pi_phi_2D_PYTHIA_can = new TCanvas("L_L_pi_phi_2D_PYTHIA_can", "L_L_pi_phi_2D_PYTHIA_can", 1200, 1000);
  L_L_pi_phi_2D_PYTHIA_can->cd();

  L0_L0_pi_phi1_vs_pi_phi2_cuts_hist->GetXaxis()->SetTitle("#phi");
  L0_L0_pi_phi1_vs_pi_phi2_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0_pi_phi1_vs_pi_phi2_cuts_hist->GetYaxis()->SetTitle("#phi");
  L0_L0_pi_phi1_vs_pi_phi2_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0_pi_phi1_vs_pi_phi2_cuts_hist->Draw("colz");

  L_L_pi_phi_2D_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_L/L_L_pi_phi_2D_PYTHIA.png");


  TCanvas *L_L_pi_phi_PYTHIA_can = new TCanvas("L_L_pi_phi_PYTHIA_can", "L_L_pi_phi_PYTHIA_can", 1200, 1000);
  L_L_pi_phi_PYTHIA_can->cd();

  TH1D *L_L_pi_phi_PYTHIA_hist = L0_L0_pi_phi1_vs_pi_phi2_cuts_hist->ProjectionX("L_L_pi_phi_PYTHIA_pr");
  L_L_pi_phi_PYTHIA_hist->GetXaxis()->SetTitle("#phi");
  L_L_pi_phi_PYTHIA_hist->GetXaxis()->CenterTitle();
  L_L_pi_phi_PYTHIA_hist->SetLineColor(1);
  L_L_pi_phi_PYTHIA_hist->SetMinimum(0);
  L_L_pi_phi_PYTHIA_hist->Draw("hist");

  TH1D *L_L_piBar_phi_PYTHIA_hist = L0_L0_pi_phi1_vs_pi_phi2_cuts_hist->ProjectionY("L_L_piBar_phi_PYTHIA_pr");
  L_L_piBar_phi_PYTHIA_hist->SetLineColor(kRed);
  L_L_piBar_phi_PYTHIA_hist->Draw("hist same");

  L_L_pi_phi_leg->Draw("same");

  L_L_pi_phi_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_L/L_L_pi_phi_PYTHIA.png");

  //--------------------------------------------------

  TCanvas *L_L_p_phi_2D_PYTHIA_can = new TCanvas("L_L_p_phi_2D_PYTHIA_can", "L_L_p_phi_2D_PYTHIA_can", 1200, 1000);
  L_L_p_phi_2D_PYTHIA_can->cd();

  L0_L0_p_phi1_vs_p_phi2_cuts_hist->GetXaxis()->SetTitle("#phi");
  L0_L0_p_phi1_vs_p_phi2_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0_p_phi1_vs_p_phi2_cuts_hist->GetYaxis()->SetTitle("#phi");
  L0_L0_p_phi1_vs_p_phi2_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0_p_phi1_vs_p_phi2_cuts_hist->Draw("colz");

  L_L_p_phi_2D_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_L/L_L_p_phi_2D_PYTHIA.png");


  TCanvas *L_L_p_phi_PYTHIA_can = new TCanvas("L_L_p_phi_PYTHIA_can", "L_L_p_phi_PYTHIA_can", 1200, 1000);
  L_L_p_phi_PYTHIA_can->cd();

  TH1D *L_L_p_phi_PYTHIA_hist = L0_L0_p_phi1_vs_p_phi2_cuts_hist->ProjectionX("L_L_p_phi_PYTHIA_pr");
  L_L_p_phi_PYTHIA_hist->GetXaxis()->SetTitle("#phi");
  L_L_p_phi_PYTHIA_hist->GetXaxis()->CenterTitle();
  L_L_p_phi_PYTHIA_hist->SetLineColor(1);
  L_L_p_phi_PYTHIA_hist->SetMinimum(0);
  L_L_p_phi_PYTHIA_hist->Draw("hist");

  TH1D *L_L_pBar_phi_PYTHIA_hist = L0_L0_p_phi1_vs_p_phi2_cuts_hist->ProjectionY("L_L_pBar_phi_PYTHIA_pr");
  L_L_pBar_phi_PYTHIA_hist->SetLineColor(kRed);
  L_L_pBar_phi_PYTHIA_hist->Draw("hist same");

  L_L_p_phi_leg->Draw("same");

  L_L_p_phi_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_L/L_L_p_phi_PYTHIA.png");

  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  TCanvas *Lbar_Lbar_pi_pT_2D_PYTHIA_can = new TCanvas("Lbar_Lbar_pi_pT_2D_PYTHIA_can", "Lbar_Lbar_pi_pT_2D_PYTHIA_can", 1200, 1000);
  Lbar_Lbar_pi_pT_2D_PYTHIA_can->cd();

  L0bar_L0bar_pi_pT1_vs_pi_pT2_cuts_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L0bar_L0bar_pi_pT1_vs_pi_pT2_cuts_hist->GetXaxis()->CenterTitle();
  L0bar_L0bar_pi_pT1_vs_pi_pT2_cuts_hist->GetYaxis()->SetTitle("p_{T} GeV/c");
  L0bar_L0bar_pi_pT1_vs_pi_pT2_cuts_hist->GetYaxis()->CenterTitle();
  L0bar_L0bar_pi_pT1_vs_pi_pT2_cuts_hist->Draw("colz");

  Lbar_Lbar_pi_pT_2D_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/Lbar_Lbar/Lbar_Lbar_pi_pT_2D_PYTHIA.png");


  TCanvas *Lbar_Lbar_pi_pT_PYTHIA_can = new TCanvas("Lbar_Lbar_pi_pT_PYTHIA_can", "Lbar_Lbar_pi_pT_PYTHIA_can", 1200, 1000);
  Lbar_Lbar_pi_pT_PYTHIA_can->cd();

  TH1D *Lbar_Lbar_pi_pT_PYTHIA_hist = L0bar_L0bar_pi_pT1_vs_pi_pT2_cuts_hist->ProjectionX("Lbar_Lbar_pi_pT_PYTHIA_pr");
  Lbar_Lbar_pi_pT_PYTHIA_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  Lbar_Lbar_pi_pT_PYTHIA_hist->GetXaxis()->CenterTitle();
  Lbar_Lbar_pi_pT_PYTHIA_hist->SetLineColor(1);
  Lbar_Lbar_pi_pT_PYTHIA_hist->Draw("hist");

  TH1D *Lbar_Lbar_piBar_pT_PYTHIA_hist = L0bar_L0bar_pi_pT1_vs_pi_pT2_cuts_hist->ProjectionY("Lbar_Lbar_piBar_pT_PYTHIA_pr");
  Lbar_Lbar_piBar_pT_PYTHIA_hist->SetLineColor(kRed);
  Lbar_Lbar_piBar_pT_PYTHIA_hist->Draw("hist same");

  Lbar_Lbar_pi_pT_leg->Draw("same");

  Lbar_Lbar_pi_pT_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/Lbar_Lbar/Lbar_Lbar_pi_pT_PYTHIA.png");

  //--------------------------------------------------

  TCanvas *Lbar_Lbar_p_pT_2D_PYTHIA_can = new TCanvas("Lbar_Lbar_p_pT_2D_PYTHIA_can", "Lbar_Lbar_p_pT_2D_PYTHIA_can", 1200, 1000);
  Lbar_Lbar_p_pT_2D_PYTHIA_can->cd();

  L0bar_L0bar_p_pT1_vs_p_pT2_cuts_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L0bar_L0bar_p_pT1_vs_p_pT2_cuts_hist->GetXaxis()->CenterTitle();
  L0bar_L0bar_p_pT1_vs_p_pT2_cuts_hist->GetYaxis()->SetTitle("p_{T} GeV/c");
  L0bar_L0bar_p_pT1_vs_p_pT2_cuts_hist->GetYaxis()->CenterTitle();
  L0bar_L0bar_p_pT1_vs_p_pT2_cuts_hist->Draw("colz");

  Lbar_Lbar_p_pT_2D_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/Lbar_Lbar/Lbar_Lbar_p_pT_2D_PYTHIA.png");


  TCanvas *Lbar_Lbar_p_pT_PYTHIA_can = new TCanvas("Lbar_Lbar_p_pT_PYTHIA_can", "Lbar_Lbar_p_pT_PYTHIA_can", 1200, 1000);
  Lbar_Lbar_p_pT_PYTHIA_can->cd();

  TH1D *Lbar_Lbar_p_pT_PYTHIA_hist = L0bar_L0bar_p_pT1_vs_p_pT2_cuts_hist->ProjectionX("Lbar_Lbar_p_pT_PYTHIA_pr");
  Lbar_Lbar_p_pT_PYTHIA_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  Lbar_Lbar_p_pT_PYTHIA_hist->GetXaxis()->CenterTitle();
  Lbar_Lbar_p_pT_PYTHIA_hist->SetLineColor(1);
  Lbar_Lbar_p_pT_PYTHIA_hist->Draw("hist");

  TH1D *Lbar_Lbar_pBar_pT_PYTHIA_hist = L0bar_L0bar_p_pT1_vs_p_pT2_cuts_hist->ProjectionY("Lbar_Lbar_pBar_pT_PYTHIA_pr");
  Lbar_Lbar_pBar_pT_PYTHIA_hist->SetLineColor(kRed);
  Lbar_Lbar_pBar_pT_PYTHIA_hist->Draw("hist same");

  Lbar_Lbar_p_pT_leg->Draw("same");

  Lbar_Lbar_p_pT_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/Lbar_Lbar/Lbar_Lbar_p_pT_PYTHIA.png");

  //--------------------------------------------------------------------------------------------------------------------------------

  //compare SE and ME partner kinematics
  //L-bar and daoughters

  TCanvas *L0_L0bar_Lbar_eta1_SE_vs_eta2_ME_cuts_can = new TCanvas("L0_L0bar_Lbar_eta1_SE_vs_eta2_ME_cuts_can", "L0_L0bar_Lbar_eta1_SE_vs_eta2_ME_cuts_can", 1200, 1000);
  L0_L0bar_Lbar_eta1_SE_vs_eta2_ME_cuts_can->cd();

  L0_L0bar_Lbar_eta1_SE_vs_eta2_ME_cuts_hist->GetXaxis()->SetTitle("y");
  L0_L0bar_Lbar_eta1_SE_vs_eta2_ME_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0bar_Lbar_eta1_SE_vs_eta2_ME_cuts_hist->GetYaxis()->SetTitle("y");
  L0_L0bar_Lbar_eta1_SE_vs_eta2_ME_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0bar_Lbar_eta1_SE_vs_eta2_ME_cuts_hist->Draw("colz");

  L0_L0bar_Lbar_eta1_SE_vs_eta2_ME_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L0_L0bar_Lbar_eta1_SE_vs_eta2_ME_cuts.png");


  TCanvas *L0_L0bar_Lbar_phi1_SE_vs_phi2_ME_cuts_can = new TCanvas("L0_L0bar_Lbar_phi1_SE_vs_phi2_ME_cuts_can", "L0_L0bar_Lbar_phi1_SE_vs_phi2_ME_cuts_can", 1200, 1000);
  L0_L0bar_Lbar_phi1_SE_vs_phi2_ME_cuts_can->cd();

  L0_L0bar_Lbar_phi1_SE_vs_phi2_ME_cuts_hist->GetXaxis()->SetTitle("#phi");
  L0_L0bar_Lbar_phi1_SE_vs_phi2_ME_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0bar_Lbar_phi1_SE_vs_phi2_ME_cuts_hist->GetYaxis()->SetTitle("#phi");
  L0_L0bar_Lbar_phi1_SE_vs_phi2_ME_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0bar_Lbar_phi1_SE_vs_phi2_ME_cuts_hist->Draw("colz");

  L0_L0bar_Lbar_phi1_SE_vs_phi2_ME_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L0_L0bar_Lbar_phi1_SE_vs_phi2_ME_cuts.png");


  TCanvas *L0_L0bar_Lbar_pT1_SE_vs_pT2_ME_cuts_can = new TCanvas("L0_L0bar_Lbar_pT1_SE_vs_pT2_ME_cuts_can", "L0_L0bar_Lbar_pT1_SE_vs_pT2_ME_cuts_can", 1200, 1000);
  L0_L0bar_Lbar_pT1_SE_vs_pT2_ME_cuts_can->cd();

  L0_L0bar_Lbar_pT1_SE_vs_pT2_ME_cuts_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L0_L0bar_Lbar_pT1_SE_vs_pT2_ME_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0bar_Lbar_pT1_SE_vs_pT2_ME_cuts_hist->GetYaxis()->SetTitle("p_{T} GeV/c");
  L0_L0bar_Lbar_pT1_SE_vs_pT2_ME_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0bar_Lbar_pT1_SE_vs_pT2_ME_cuts_hist->Draw("colz");

  L0_L0bar_Lbar_pT1_SE_vs_pT2_ME_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L0_L0bar_Lbar_pT1_SE_vs_pT2_ME_cuts.png");

  //------------------------------------------------------

  //compare daughter kinemtics of SE and ME partner

  TCanvas *L0_L0bar_pbar_eta1_SE_vs_eta2_ME_cuts_can = new TCanvas("L0_L0bar_pbar_eta1_SE_vs_eta2_ME_cuts_can", "L0_L0bar_pbar_eta1_SE_vs_eta2_ME_cuts_can", 1200, 1000);
  L0_L0bar_pbar_eta1_SE_vs_eta2_ME_cuts_can->cd();

  L0_L0bar_pbar_eta1_SE_vs_eta2_ME_cuts_hist->GetXaxis()->SetTitle("#eta");
  L0_L0bar_pbar_eta1_SE_vs_eta2_ME_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0bar_pbar_eta1_SE_vs_eta2_ME_cuts_hist->GetYaxis()->SetTitle("#eta");
  L0_L0bar_pbar_eta1_SE_vs_eta2_ME_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0bar_pbar_eta1_SE_vs_eta2_ME_cuts_hist->Draw("colz");

  L0_L0bar_pbar_eta1_SE_vs_eta2_ME_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L0_L0bar_pbar_eta1_SE_vs_eta2_ME_cuts.png");


  TCanvas *L0_L0bar_pbar_phi1_SE_vs_phi2_ME_cuts_can = new TCanvas("L0_L0bar_pbar_phi1_SE_vs_phi2_ME_cuts_can", "L0_L0bar_pbar_phi1_SE_vs_phi2_ME_cuts_can", 1200, 1000);
  L0_L0bar_pbar_phi1_SE_vs_phi2_ME_cuts_can->cd();

  L0_L0bar_pbar_phi1_SE_vs_phi2_ME_cuts_hist->GetXaxis()->SetTitle("#phi");
  L0_L0bar_pbar_phi1_SE_vs_phi2_ME_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0bar_pbar_phi1_SE_vs_phi2_ME_cuts_hist->GetYaxis()->SetTitle("#phi");
  L0_L0bar_pbar_phi1_SE_vs_phi2_ME_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0bar_pbar_phi1_SE_vs_phi2_ME_cuts_hist->Draw("colz");

  L0_L0bar_pbar_phi1_SE_vs_phi2_ME_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L0_L0bar_pbar_phi1_SE_vs_phi2_ME_cuts.png");


  TCanvas *L0_L0bar_pbar_pT1_SE_vs_pT2_ME_cuts_can = new TCanvas("L0_L0bar_pbar_pT1_SE_vs_pT2_ME_cuts_can", "L0_L0bar_pbar_pT1_SE_vs_pT2_ME_cuts_can", 1200, 1000);
  L0_L0bar_pbar_pT1_SE_vs_pT2_ME_cuts_can->cd();

  L0_L0bar_pbar_pT1_SE_vs_pT2_ME_cuts_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L0_L0bar_pbar_pT1_SE_vs_pT2_ME_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0bar_pbar_pT1_SE_vs_pT2_ME_cuts_hist->GetYaxis()->SetTitle("p_{T} GeV/c");
  L0_L0bar_pbar_pT1_SE_vs_pT2_ME_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0bar_pbar_pT1_SE_vs_pT2_ME_cuts_hist->Draw("colz");

  L0_L0bar_pbar_pT1_SE_vs_pT2_ME_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L0_L0bar_pbar_pT1_SE_vs_pT2_ME_cuts.png");

  //---------------------------------------------------

  TCanvas *L0_L0bar_pibar_eta1_SE_vs_eta2_ME_cuts_can = new TCanvas("L0_L0bar_pibar_eta1_SE_vs_eta2_ME_cuts_can", "L0_L0bar_pibar_eta1_SE_vs_eta2_ME_cuts_can", 1200, 1000);
  L0_L0bar_pibar_eta1_SE_vs_eta2_ME_cuts_can->cd();

  L0_L0bar_pibar_eta1_SE_vs_eta2_ME_cuts_hist->GetXaxis()->SetTitle("#eta");
  L0_L0bar_pibar_eta1_SE_vs_eta2_ME_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0bar_pibar_eta1_SE_vs_eta2_ME_cuts_hist->GetYaxis()->SetTitle("#eta");
  L0_L0bar_pibar_eta1_SE_vs_eta2_ME_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0bar_pibar_eta1_SE_vs_eta2_ME_cuts_hist->Draw("colz");

  L0_L0bar_pibar_eta1_SE_vs_eta2_ME_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L0_L0bar_pibar_eta1_SE_vs_eta2_ME_cuts.png");


  TCanvas *L0_L0bar_pibar_phi1_SE_vs_phi2_ME_cuts_can = new TCanvas("L0_L0bar_pibar_phi1_SE_vs_phi2_ME_cuts_can", "L0_L0bar_pibar_phi1_SE_vs_phi2_ME_cuts_can", 1200, 1000);
  L0_L0bar_pibar_phi1_SE_vs_phi2_ME_cuts_can->cd();

  L0_L0bar_pibar_phi1_SE_vs_phi2_ME_cuts_hist->GetXaxis()->SetTitle("#phi");
  L0_L0bar_pibar_phi1_SE_vs_phi2_ME_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0bar_pibar_phi1_SE_vs_phi2_ME_cuts_hist->GetYaxis()->SetTitle("#phi");
  L0_L0bar_pibar_phi1_SE_vs_phi2_ME_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0bar_pibar_phi1_SE_vs_phi2_ME_cuts_hist->Draw("colz");

  L0_L0bar_pibar_phi1_SE_vs_phi2_ME_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L0_L0bar_pibar_phi1_SE_vs_phi2_ME_cuts.png");


  TCanvas *L0_L0bar_pibar_pT1_SE_vs_pT2_ME_cuts_can = new TCanvas("L0_L0bar_pibar_pT1_SE_vs_pT2_ME_cuts_can", "L0_L0bar_pibar_pT1_SE_vs_pT2_ME_cuts_can", 1200, 1000);
  L0_L0bar_pibar_pT1_SE_vs_pT2_ME_cuts_can->cd();

  L0_L0bar_pibar_pT1_SE_vs_pT2_ME_cuts_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L0_L0bar_pibar_pT1_SE_vs_pT2_ME_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0bar_pibar_pT1_SE_vs_pT2_ME_cuts_hist->GetYaxis()->SetTitle("p_{T} GeV/c");
  L0_L0bar_pibar_pT1_SE_vs_pT2_ME_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0bar_pibar_pT1_SE_vs_pT2_ME_cuts_hist->Draw("colz");

  L0_L0bar_pibar_pT1_SE_vs_pT2_ME_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L0_L0bar_pibar_pT1_SE_vs_pT2_ME_cuts.png");


  //________________________________________________________________________________________________________________________________________________________________________________________________

  //embedding
  TCanvas *L_Lbar_pi_pT_2D_embedding_can = new TCanvas("L_Lbar_pi_pT_2D_embedding_can", "L_Lbar_pi_pT_2D_embedding_can", 1200, 1000);
  L_Lbar_pi_pT_2D_embedding_can->cd();

  L0_L0bar_pi1_pT1_vs_pi2_pT2_RC_US_hist->Add(L0_L0bar_pi1_pT1_vs_pi2_pT2_RC_US_LS_hist, -1); //subtract combinatorial background
  L0_L0bar_pi1_pT1_vs_pi2_pT2_RC_US_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L0_L0bar_pi1_pT1_vs_pi2_pT2_RC_US_hist->GetXaxis()->CenterTitle();
  L0_L0bar_pi1_pT1_vs_pi2_pT2_RC_US_hist->GetYaxis()->SetTitle("p_{T} GeV/c");
  L0_L0bar_pi1_pT1_vs_pi2_pT2_RC_US_hist->GetYaxis()->CenterTitle();
  L0_L0bar_pi1_pT1_vs_pi2_pT2_RC_US_hist->Draw("colz");

  L_Lbar_pi_pT_2D_embedding_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L_Lbar_pi_pT_2D_embedding.png");


  TCanvas *L_Lbar_pi_pT_embedding_can = new TCanvas("L_Lbar_pi_pT_embedding_can", "L_Lbar_pi_pT_embedding_can", 1200, 1000);
  L_Lbar_pi_pT_embedding_can->cd();

  //L0_L0bar_pi1_pT1_vs_pi2_pT2_RC_US_hist->Add(L0_L0bar_pi1_pT1_vs_pi2_pT2_RC_US_LS_hist, -1); //subtract combinatorial background

  TH1D *L_Lbar_pi_pT_embedding_hist = L0_L0bar_pi1_pT1_vs_pi2_pT2_RC_US_hist->ProjectionX("L_Lbar_pi_pT_embed_pr");
  L_Lbar_pi_pT_embedding_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L_Lbar_pi_pT_embedding_hist->GetXaxis()->CenterTitle();
  L_Lbar_pi_pT_embedding_hist->SetLineColor(1);
  L_Lbar_pi_pT_embedding_hist->Draw("hist");

  TH1D *L_Lbar_piBar_pT_embedding_hist = L0_L0bar_pi1_pT1_vs_pi2_pT2_RC_US_hist->ProjectionY("L_Lbar_piBar_pT_embed_pr");
  L_Lbar_piBar_pT_embedding_hist->SetLineColor(kRed);
  L_Lbar_piBar_pT_embedding_hist->Draw("hist same");

  L_Lbar_pi_pT_leg->Draw("same");

  L_Lbar_pi_pT_embedding_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L_Lbar_pi_pT_embedding.png");

  //--------------------------------------------------

  TCanvas *L_Lbar_pi_eta_2D_embedding_can = new TCanvas("L_Lbar_pi_eta_2D_embedding_can", "L_Lbar_pi_eta_2D_embedding_can", 1200, 1000);
  L_Lbar_pi_eta_2D_embedding_can->cd();

  L0_L0bar_pi1_eta1_vs_pi2_eta2_RC_US_hist->Add(L0_L0bar_pi1_eta1_vs_pi2_eta2_RC_US_LS_hist, -1); //subtract combinatorial background
  L0_L0bar_pi1_eta1_vs_pi2_eta2_RC_US_hist->GetXaxis()->SetTitle("#eta");
  L0_L0bar_pi1_eta1_vs_pi2_eta2_RC_US_hist->GetXaxis()->CenterTitle();
  L0_L0bar_pi1_eta1_vs_pi2_eta2_RC_US_hist->GetYaxis()->SetTitle("#eta");
  L0_L0bar_pi1_eta1_vs_pi2_eta2_RC_US_hist->GetYaxis()->CenterTitle();
  L0_L0bar_pi1_eta1_vs_pi2_eta2_RC_US_hist->Draw("colz");

  L_Lbar_pi_eta_2D_embedding_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L_Lbar_pi_eta_2D_embedding.png");


  TCanvas *L_Lbar_pi_eta_embedding_can = new TCanvas("L_Lbar_pi_eta_embedding_can", "L_Lbar_pi_eta_embedding_can", 1200, 1000);
  L_Lbar_pi_eta_embedding_can->cd();

  //L0_L0bar_pi1_eta1_vs_pi2_eta2_RC_US_hist->Add(L0_L0bar_pi1_eta1_vs_pi2_eta2_RC_US_LS_hist, -1); //subtract combinatorial background

  TH1D *L_Lbar_pi_eta_embedding_hist = L0_L0bar_pi1_eta1_vs_pi2_eta2_RC_US_hist->ProjectionX("L_Lbar_pi_eta_embed_pr");
  L_Lbar_pi_eta_embedding_hist->GetXaxis()->SetTitle("#eta");
  L_Lbar_pi_eta_embedding_hist->GetXaxis()->CenterTitle();
  L_Lbar_pi_eta_embedding_hist->SetLineColor(1);
  L_Lbar_pi_eta_embedding_hist->SetMinimum(0);
  L_Lbar_pi_eta_embedding_hist->Draw("hist");

  TH1D *L_Lbar_piBar_eta_embedding_hist = L0_L0bar_pi1_eta1_vs_pi2_eta2_RC_US_hist->ProjectionY("L_Lbar_piBar_eta_embed_pr");
  L_Lbar_piBar_eta_embedding_hist->SetLineColor(kRed);
  L_Lbar_piBar_eta_embedding_hist->Draw("hist same");

  L_Lbar_pi_pT_leg->Draw("same");

  L_Lbar_pi_eta_embedding_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L_Lbar_pi_eta_embedding.png");

  //--------------------------------------------------

  TCanvas *L_Lbar_pi_phi_2D_embedding_can = new TCanvas("L_Lbar_pi_phi_2D_embedding_can", "L_Lbar_pi_phi_2D_embedding_can", 1200, 1000);
  L_Lbar_pi_phi_2D_embedding_can->cd();

  L0_L0bar_pi1_phi1_vs_pi2_phi2_RC_US_hist->Add(L0_L0bar_pi1_phi1_vs_pi2_phi2_RC_US_LS_hist, -1); //subtract combinatorial background
  L0_L0bar_pi1_phi1_vs_pi2_phi2_RC_US_hist->GetXaxis()->SetTitle("#phi");
  L0_L0bar_pi1_phi1_vs_pi2_phi2_RC_US_hist->GetXaxis()->CenterTitle();
  L0_L0bar_pi1_phi1_vs_pi2_phi2_RC_US_hist->GetYaxis()->SetTitle("#phi");
  L0_L0bar_pi1_phi1_vs_pi2_phi2_RC_US_hist->GetYaxis()->CenterTitle();
  L0_L0bar_pi1_phi1_vs_pi2_phi2_RC_US_hist->Draw("colz");

  L_Lbar_pi_phi_2D_embedding_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L_Lbar_pi_phi_2D_embedding.png");


  TCanvas *L_Lbar_pi_phi_embedding_can = new TCanvas("L_Lbar_pi_phi_embedding_can", "L_Lbar_pi_phi_embedding_can", 1200, 1000);
  L_Lbar_pi_phi_embedding_can->cd();

  //L0_L0bar_pi1_phi1_vs_pi2_phi2_RC_US_hist->Add(L0_L0bar_pi1_phi1_vs_pi2_phi2_RC_US_LS_hist, -1); //subtract combinatorial background

  TH1D *L_Lbar_pi_phi_embedding_hist = L0_L0bar_pi1_phi1_vs_pi2_phi2_RC_US_hist->ProjectionX("L_Lbar_pi_phi_embed_pr");
  L_Lbar_pi_phi_embedding_hist->GetXaxis()->SetTitle("#phi");
  L_Lbar_pi_phi_embedding_hist->GetXaxis()->CenterTitle();
  L_Lbar_pi_phi_embedding_hist->SetLineColor(1);
  L_Lbar_pi_phi_embedding_hist->SetMinimum(0);
  L_Lbar_pi_phi_embedding_hist->Draw("hist");

  TH1D *L_Lbar_piBar_phi_embedding_hist = L0_L0bar_pi1_phi1_vs_pi2_phi2_RC_US_hist->ProjectionY("L_Lbar_piBar_phi_embed_pr");
  L_Lbar_piBar_phi_embedding_hist->SetLineColor(kRed);
  L_Lbar_piBar_phi_embedding_hist->Draw("hist same");

  L_Lbar_pi_pT_leg->Draw("same");

  L_Lbar_pi_phi_embedding_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L_Lbar_pi_phi_embedding.png");

  //--------------------------------------------------------------------------------------------------------------------------------------


  TCanvas *L_Lbar_p_pT_2D_embedding_can = new TCanvas("L_Lbar_p_pT_2D_embedding_can", "L_Lbar_p_pT_2D_embedding_can", 1200, 1000);
  L_Lbar_p_pT_2D_embedding_can->cd();

  L0_L0bar_p1_pT1_vs_p2_pT2_RC_US_hist->Add(L0_L0bar_p1_pT1_vs_p2_pT2_RC_US_LS_hist, -1); //subtract combinatorial background
  L0_L0bar_p1_pT1_vs_p2_pT2_RC_US_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L0_L0bar_p1_pT1_vs_p2_pT2_RC_US_hist->GetXaxis()->CenterTitle();
  L0_L0bar_p1_pT1_vs_p2_pT2_RC_US_hist->GetYaxis()->SetTitle("p_{T} GeV/c");
  L0_L0bar_p1_pT1_vs_p2_pT2_RC_US_hist->GetYaxis()->CenterTitle();
  L0_L0bar_p1_pT1_vs_p2_pT2_RC_US_hist->Draw("colz");

  L_Lbar_p_pT_2D_embedding_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L_Lbar_p_pT_2D_embedding.png");

  TCanvas *L_Lbar_p_pT_embedding_can = new TCanvas("L_Lbar_p_pT_embedding_can", "L_Lbar_p_pT_embedding_can", 1200, 1000);
  L_Lbar_p_pT_embedding_can->cd();

  TH1D *L_Lbar_p_pT_embedding_hist = L0_L0bar_p1_pT1_vs_p2_pT2_RC_US_hist->ProjectionX("L_Lbar_p_pT_embed_pr");
  L_Lbar_p_pT_embedding_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L_Lbar_p_pT_embedding_hist->GetXaxis()->CenterTitle();
  L_Lbar_p_pT_embedding_hist->SetLineColor(1);
  L_Lbar_p_pT_embedding_hist->Draw("hist");

  TH1D *L_Lbar_pBar_pT_embedding_hist = L0_L0bar_p1_pT1_vs_p2_pT2_RC_US_hist->ProjectionY("L_Lbar_pBar_pT_embed_pr");
  L_Lbar_pBar_pT_embedding_hist->SetLineColor(kRed);
  L_Lbar_pBar_pT_embedding_hist->Draw("hist same");

  L_Lbar_p_pT_leg->Draw("same");

  L_Lbar_p_pT_embedding_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L_Lbar_p_pT_embedding.png");

  //--------------------------------------------------

  TCanvas *L_Lbar_p_eta_2D_embedding_can = new TCanvas("L_Lbar_p_eta_2D_embedding_can", "L_Lbar_p_eta_2D_embedding_can", 1200, 1000);
  L_Lbar_p_eta_2D_embedding_can->cd();

  L0_L0bar_p1_eta1_vs_p2_eta2_RC_US_hist->Add(L0_L0bar_p1_eta1_vs_p2_eta2_RC_US_LS_hist, -1); //subtract combinatorial background
  L0_L0bar_p1_eta1_vs_p2_eta2_RC_US_hist->GetXaxis()->SetTitle("#eta");
  L0_L0bar_p1_eta1_vs_p2_eta2_RC_US_hist->GetXaxis()->CenterTitle();
  L0_L0bar_p1_eta1_vs_p2_eta2_RC_US_hist->GetYaxis()->SetTitle("#eta");
  L0_L0bar_p1_eta1_vs_p2_eta2_RC_US_hist->GetYaxis()->CenterTitle();
  L0_L0bar_p1_eta1_vs_p2_eta2_RC_US_hist->Draw("colz");

  L_Lbar_p_eta_2D_embedding_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L_Lbar_p_eta_2D_embedding.png");

  TCanvas *L_Lbar_p_eta_embedding_can = new TCanvas("L_Lbar_p_eta_embedding_can", "L_Lbar_p_eta_embedding_can", 1200, 1000);
  L_Lbar_p_eta_embedding_can->cd();

  TH1D *L_Lbar_p_eta_embedding_hist = L0_L0bar_p1_eta1_vs_p2_eta2_RC_US_hist->ProjectionX("L_Lbar_p_eta_embed_pr");
  L_Lbar_p_eta_embedding_hist->GetXaxis()->SetTitle("#eta");
  L_Lbar_p_eta_embedding_hist->GetXaxis()->CenterTitle();
  L_Lbar_p_eta_embedding_hist->SetLineColor(1);
  L_Lbar_p_eta_embedding_hist->SetMinimum(0);
  L_Lbar_p_eta_embedding_hist->Draw("hist");

  TH1D *L_Lbar_pBar_eta_embedding_hist = L0_L0bar_p1_eta1_vs_p2_eta2_RC_US_hist->ProjectionY("L_Lbar_pBar_eta_embed_pr");
  L_Lbar_pBar_eta_embedding_hist->SetLineColor(kRed);
  L_Lbar_pBar_eta_embedding_hist->Draw("hist same");

  L_Lbar_p_pT_leg->Draw("same");

  L_Lbar_p_eta_embedding_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L_Lbar_p_eta_embedding.png");

  //--------------------------------------------------

  TCanvas *L_Lbar_p_phi_2D_embedding_can = new TCanvas("L_Lbar_p_phi_2D_embedding_can", "L_Lbar_p_phi_2D_embedding_can", 1200, 1000);
  L_Lbar_p_phi_2D_embedding_can->cd();

  L0_L0bar_p1_phi1_vs_p2_phi2_RC_US_hist->Add(L0_L0bar_p1_phi1_vs_p2_phi2_RC_US_LS_hist, -1); //subtract combinatorial background
  L0_L0bar_p1_phi1_vs_p2_phi2_RC_US_hist->GetXaxis()->SetTitle("#phi");
  L0_L0bar_p1_phi1_vs_p2_phi2_RC_US_hist->GetXaxis()->CenterTitle();
  L0_L0bar_p1_phi1_vs_p2_phi2_RC_US_hist->GetYaxis()->SetTitle("#phi");
  L0_L0bar_p1_phi1_vs_p2_phi2_RC_US_hist->GetYaxis()->CenterTitle();
  L0_L0bar_p1_phi1_vs_p2_phi2_RC_US_hist->Draw("colz");

  L_Lbar_p_phi_2D_embedding_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L_Lbar_p_phi_2D_embedding.png");



  TCanvas *L_Lbar_p_phi_embedding_can = new TCanvas("L_Lbar_p_phi_embedding_can", "L_Lbar_p_phi_embedding_can", 1200, 1000);
  L_Lbar_p_phi_embedding_can->cd();

  TH1D *L_Lbar_p_phi_embedding_hist = L0_L0bar_p1_phi1_vs_p2_phi2_RC_US_hist->ProjectionX("L_Lbar_p_phi_embed_pr");
  L_Lbar_p_phi_embedding_hist->GetXaxis()->SetTitle("#phi");
  L_Lbar_p_phi_embedding_hist->GetXaxis()->CenterTitle();
  L_Lbar_p_phi_embedding_hist->SetLineColor(1);
  L_Lbar_p_phi_embedding_hist->SetMinimum(0);
  L_Lbar_p_phi_embedding_hist->Draw("hist");

  TH1D *L_Lbar_pBar_phi_embedding_hist = L0_L0bar_p1_phi1_vs_p2_phi2_RC_US_hist->ProjectionY("L_Lbar_pBar_phi_embed_pr");
  L_Lbar_pBar_phi_embedding_hist->SetLineColor(kRed);
  L_Lbar_pBar_phi_embedding_hist->Draw("hist same");

  L_Lbar_p_pT_leg->Draw("same");

  L_Lbar_p_phi_embedding_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_Lbar/L_Lbar_p_phi_embedding.png");

  //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  TCanvas *L_L_pi_pT_2D_embedding_can = new TCanvas("L_L_pi_pT_2D_embedding_can", "L_L_pi_pT_2D_embedding_can", 1200, 1000);
  L_L_pi_pT_2D_embedding_can->cd();

  L0_L0_pi1_pT1_vs_pi2_pT2_RC_US_hist->Add(L0_L0_pi1_pT1_vs_pi2_pT2_RC_US_LS_hist, -1); //subtract combinatorial background
  L0_L0_pi1_pT1_vs_pi2_pT2_RC_US_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L0_L0_pi1_pT1_vs_pi2_pT2_RC_US_hist->GetXaxis()->CenterTitle();
  L0_L0_pi1_pT1_vs_pi2_pT2_RC_US_hist->GetYaxis()->SetTitle("p_{T} GeV/c");
  L0_L0_pi1_pT1_vs_pi2_pT2_RC_US_hist->GetYaxis()->CenterTitle();
  L0_L0_pi1_pT1_vs_pi2_pT2_RC_US_hist->Draw("colz");

  L_L_pi_pT_2D_embedding_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_L/L_L_pi_pT_2D_embedding.png");


  TCanvas *L_L_pi_pT_embedding_can = new TCanvas("L_L_pi_pT_embedding_can", "L_L_pi_pT_embedding_can", 1200, 1000);
  L_L_pi_pT_embedding_can->cd();

  //L0_L0_pi1_pT1_vs_pi2_pT2_RC_US_hist->Add(L0_L0_pi1_pT1_vs_pi2_pT2_RC_US_LS_hist, -1); //subtract combinatorial background

  TH1D *L_L_pi_pT_embedding_hist = L0_L0_pi1_pT1_vs_pi2_pT2_RC_US_hist->ProjectionX("L_L_pi_pT_embed_pr");
  L_L_pi_pT_embedding_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L_L_pi_pT_embedding_hist->GetXaxis()->CenterTitle();
  L_L_pi_pT_embedding_hist->SetLineColor(1);
  L_L_pi_pT_embedding_hist->Draw("hist");

  TH1D *L_L_piBar_pT_embedding_hist = L0_L0_pi1_pT1_vs_pi2_pT2_RC_US_hist->ProjectionY("L_L_piBar_pT_embed_pr");
  L_L_piBar_pT_embedding_hist->SetLineColor(kRed);
  L_L_piBar_pT_embedding_hist->Draw("hist same");

  L_L_pi_pT_leg->Draw("same");

  L_L_pi_pT_embedding_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_L/L_L_pi_pT_embedding.png");

  //--------------------------------------------------

  TCanvas *L_L_pi_eta_2D_embedding_can = new TCanvas("L_L_pi_eta_2D_embedding_can", "L_L_pi_eta_2D_embedding_can", 1200, 1000);
  L_L_pi_eta_2D_embedding_can->cd();

  L0_L0_pi1_eta1_vs_pi2_eta2_RC_US_hist->Add(L0_L0_pi1_eta1_vs_pi2_eta2_RC_US_LS_hist, -1); //subtract combinatorial background
  L0_L0_pi1_eta1_vs_pi2_eta2_RC_US_hist->GetXaxis()->SetTitle("#eta");
  L0_L0_pi1_eta1_vs_pi2_eta2_RC_US_hist->GetXaxis()->CenterTitle();
  L0_L0_pi1_eta1_vs_pi2_eta2_RC_US_hist->GetYaxis()->SetTitle("#eta");
  L0_L0_pi1_eta1_vs_pi2_eta2_RC_US_hist->GetYaxis()->CenterTitle();
  L0_L0_pi1_eta1_vs_pi2_eta2_RC_US_hist->Draw("colz");

  L_L_pi_eta_2D_embedding_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_L/L_L_pi_eta_2D_embedding.png");


  TCanvas *L_L_pi_eta_embedding_can = new TCanvas("L_L_pi_eta_embedding_can", "L_L_pi_eta_embedding_can", 1200, 1000);
  L_L_pi_eta_embedding_can->cd();

  //L0_L0_pi1_eta1_vs_pi2_eta2_RC_US_hist->Add(L0_L0_pi1_eta1_vs_pi2_eta2_RC_US_LS_hist, -1); //subtract combinatorial background

  TH1D *L_L_pi_eta_embedding_hist = L0_L0_pi1_eta1_vs_pi2_eta2_RC_US_hist->ProjectionX("L_L_pi_eta_embed_pr");
  L_L_pi_eta_embedding_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L_L_pi_eta_embedding_hist->GetXaxis()->CenterTitle();
  L_L_pi_eta_embedding_hist->SetLineColor(1);
  L_L_pi_eta_embedding_hist->SetMinimum(0);
  L_L_pi_eta_embedding_hist->Draw("hist");

  TH1D *L_L_piBar_eta_embedding_hist = L0_L0_pi1_eta1_vs_pi2_eta2_RC_US_hist->ProjectionY("L_L_piBar_eta_embed_pr");
  L_L_piBar_eta_embedding_hist->SetLineColor(kRed);
  L_L_piBar_eta_embedding_hist->Draw("hist same");

  L_L_pi_pT_leg->Draw("same");

  L_L_pi_eta_embedding_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_L/L_L_pi_eta_embedding.png");

  //--------------------------------------------------

  TCanvas *L_L_pi_phi_2D_embedding_can = new TCanvas("L_L_pi_phi_2D_embedding_can", "L_L_pi_phi_2D_embedding_can", 1200, 1000);
  L_L_pi_phi_2D_embedding_can->cd();

  L0_L0_pi1_phi1_vs_pi2_phi2_RC_US_hist->Add(L0_L0_pi1_phi1_vs_pi2_phi2_RC_US_LS_hist, -1); //subtract combinatorial background
  L0_L0_pi1_phi1_vs_pi2_phi2_RC_US_hist->GetXaxis()->SetTitle("#phi");
  L0_L0_pi1_phi1_vs_pi2_phi2_RC_US_hist->GetXaxis()->CenterTitle();
  L0_L0_pi1_phi1_vs_pi2_phi2_RC_US_hist->GetYaxis()->SetTitle("#phi");
  L0_L0_pi1_phi1_vs_pi2_phi2_RC_US_hist->GetYaxis()->CenterTitle();
  L0_L0_pi1_phi1_vs_pi2_phi2_RC_US_hist->Draw("colz");

  L_L_pi_phi_2D_embedding_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_L/L_L_pi_phi_2D_embedding.png");


  TCanvas *L_L_pi_phi_embedding_can = new TCanvas("L_L_pi_phi_embedding_can", "L_L_pi_phi_embedding_can", 1200, 1000);
  L_L_pi_phi_embedding_can->cd();

  //L0_L0_pi1_phi1_vs_pi2_phi2_RC_US_hist->Add(L0_L0_pi1_phi1_vs_pi2_phi2_RC_US_LS_hist, -1); //subtract combinatorial background

  TH1D *L_L_pi_phi_embedding_hist = L0_L0_pi1_phi1_vs_pi2_phi2_RC_US_hist->ProjectionX("L_L_pi_phi_embed_pr");
  L_L_pi_phi_embedding_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L_L_pi_phi_embedding_hist->GetXaxis()->CenterTitle();
  L_L_pi_phi_embedding_hist->SetLineColor(1);
  L_L_pi_phi_embedding_hist->SetMinimum(0);
  L_L_pi_phi_embedding_hist->Draw("hist");

  TH1D *L_L_piBar_phi_embedding_hist = L0_L0_pi1_phi1_vs_pi2_phi2_RC_US_hist->ProjectionY("L_L_piBar_phi_embed_pr");
  L_L_piBar_phi_embedding_hist->SetLineColor(kRed);
  L_L_piBar_phi_embedding_hist->Draw("hist same");

  L_L_pi_pT_leg->Draw("same");

  L_L_pi_phi_embedding_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_L/L_L_pi_phi_embedding.png");

  //--------------------------------------------------------------------------------------------------------------------------------------


  TCanvas *L_L_p_pT_2D_embedding_can = new TCanvas("L_L_p_pT_2D_embedding_can", "L_L_p_pT_2D_embedding_can", 1200, 1000);
  L_L_p_pT_2D_embedding_can->cd();

  L0_L0_p1_pT1_vs_p2_pT2_RC_US_hist->Add(L0_L0_p1_pT1_vs_p2_pT2_RC_US_LS_hist, -1); //subtract combinatorial background
  L0_L0_p1_pT1_vs_p2_pT2_RC_US_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L0_L0_p1_pT1_vs_p2_pT2_RC_US_hist->GetXaxis()->CenterTitle();
  L0_L0_p1_pT1_vs_p2_pT2_RC_US_hist->GetYaxis()->SetTitle("p_{T} GeV/c");
  L0_L0_p1_pT1_vs_p2_pT2_RC_US_hist->GetYaxis()->CenterTitle();
  L0_L0_p1_pT1_vs_p2_pT2_RC_US_hist->Draw("colz");

  L_L_p_pT_2D_embedding_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_L/L_L_p_pT_2D_embedding.png");

  TCanvas *L_L_p_pT_embedding_can = new TCanvas("L_L_p_pT_embedding_can", "L_L_p_pT_embedding_can", 1200, 1000);
  L_L_p_pT_embedding_can->cd();

  TH1D *L_L_p_pT_embedding_hist = L0_L0_p1_pT1_vs_p2_pT2_RC_US_hist->ProjectionX("L_L_p_pT_embed_pr");
  L_L_p_pT_embedding_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L_L_p_pT_embedding_hist->GetXaxis()->CenterTitle();
  L_L_p_pT_embedding_hist->SetLineColor(1);
  L_L_p_pT_embedding_hist->Draw("hist");

  TH1D *L_L_pBar_pT_embedding_hist = L0_L0_p1_pT1_vs_p2_pT2_RC_US_hist->ProjectionY("L_L_pBar_pT_embed_pr");
  L_L_pBar_pT_embedding_hist->SetLineColor(kRed);
  L_L_pBar_pT_embedding_hist->Draw("hist same");

  L_L_p_pT_leg->Draw("same");

  L_L_p_pT_embedding_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_L/L_L_p_pT_embedding.png");

  //--------------------------------------------------

  TCanvas *L_L_p_eta_2D_embedding_can = new TCanvas("L_L_p_eta_2D_embedding_can", "L_L_p_eta_2D_embedding_can", 1200, 1000);
  L_L_p_eta_2D_embedding_can->cd();

  L0_L0_p1_eta1_vs_p2_eta2_RC_US_hist->Add(L0_L0_p1_eta1_vs_p2_eta2_RC_US_LS_hist, -1); //subtract combinatorial background
  L0_L0_p1_eta1_vs_p2_eta2_RC_US_hist->GetXaxis()->SetTitle("#eta");
  L0_L0_p1_eta1_vs_p2_eta2_RC_US_hist->GetXaxis()->CenterTitle();
  L0_L0_p1_eta1_vs_p2_eta2_RC_US_hist->GetYaxis()->SetTitle("#eta");
  L0_L0_p1_eta1_vs_p2_eta2_RC_US_hist->GetYaxis()->CenterTitle();
  L0_L0_p1_eta1_vs_p2_eta2_RC_US_hist->Draw("colz");

  L_L_p_eta_2D_embedding_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_L/L_L_p_eta_2D_embedding.png");

  TCanvas *L_L_p_eta_embedding_can = new TCanvas("L_L_p_eta_embedding_can", "L_L_p_eta_embedding_can", 1200, 1000);
  L_L_p_eta_embedding_can->cd();

  TH1D *L_L_p_eta_embedding_hist = L0_L0_p1_eta1_vs_p2_eta2_RC_US_hist->ProjectionX("L_L_p_eta_embed_pr");
  L_L_p_eta_embedding_hist->GetXaxis()->SetTitle("#eta");
  L_L_p_eta_embedding_hist->GetXaxis()->CenterTitle();
  L_L_p_eta_embedding_hist->SetLineColor(1);
  L_L_p_eta_embedding_hist->SetMinimum(0);
  L_L_p_eta_embedding_hist->Draw("hist");

  TH1D *L_L_pBar_eta_embedding_hist = L0_L0_p1_eta1_vs_p2_eta2_RC_US_hist->ProjectionY("L_L_pBar_eta_embed_pr");
  L_L_pBar_eta_embedding_hist->SetLineColor(kRed);
  L_L_pBar_eta_embedding_hist->Draw("hist same");

  L_L_p_pT_leg->Draw("same");

  L_L_p_eta_embedding_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_L/L_L_p_eta_embedding.png");

  //--------------------------------------------------

  TCanvas *L_L_p_phi_2D_embedding_can = new TCanvas("L_L_p_phi_2D_embedding_can", "L_L_p_phi_2D_embedding_can", 1200, 1000);
  L_L_p_phi_2D_embedding_can->cd();

  L0_L0_p1_phi1_vs_p2_phi2_RC_US_hist->Add(L0_L0_p1_phi1_vs_p2_phi2_RC_US_LS_hist, -1); //subtract combinatorial background
  L0_L0_p1_phi1_vs_p2_phi2_RC_US_hist->GetXaxis()->SetTitle("#phi");
  L0_L0_p1_phi1_vs_p2_phi2_RC_US_hist->GetXaxis()->CenterTitle();
  L0_L0_p1_phi1_vs_p2_phi2_RC_US_hist->GetYaxis()->SetTitle("#phi");
  L0_L0_p1_phi1_vs_p2_phi2_RC_US_hist->GetYaxis()->CenterTitle();
  L0_L0_p1_phi1_vs_p2_phi2_RC_US_hist->Draw("colz");

  L_L_p_phi_2D_embedding_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_L/L_L_p_phi_2D_embedding.png");



  TCanvas *L_L_p_phi_embedding_can = new TCanvas("L_L_p_phi_embedding_can", "L_L_p_phi_embedding_can", 1200, 1000);
  L_L_p_phi_embedding_can->cd();

  TH1D *L_L_p_phi_embedding_hist = L0_L0_p1_phi1_vs_p2_phi2_RC_US_hist->ProjectionX("L_L_p_phi_embed_pr");
  L_L_p_phi_embedding_hist->GetXaxis()->SetTitle("#phi");
  L_L_p_phi_embedding_hist->GetXaxis()->CenterTitle();
  L_L_p_phi_embedding_hist->SetLineColor(1);
  L_L_p_phi_embedding_hist->SetMinimum(0);
  L_L_p_phi_embedding_hist->Draw("hist");

  TH1D *L_L_pBar_phi_embedding_hist = L0_L0_p1_phi1_vs_p2_phi2_RC_US_hist->ProjectionY("L_L_pBar_phi_embed_pr");
  L_L_pBar_phi_embedding_hist->SetLineColor(kRed);
  L_L_pBar_phi_embedding_hist->Draw("hist same");

  L_L_p_pT_leg->Draw("same");

  L_L_p_phi_embedding_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L_L/L_L_p_phi_embedding.png");

  //________________________________________________________________________________________________________________________________________________________________________________________________

  //calculate and save weight histograms

  //L-Lbar kinematics

  TCanvas *L_Lbar_pT_data_can = new TCanvas("L_Lbar_pT_data_can", "L_Lbar_pT_data_can", 1200, 1000);
  L_Lbar_pT_data_can->cd();

  L0_L0bar_pT1_vs_pT2_US_hist->Add(L0_L0bar_pT1_vs_pT2_US_LS_hist, -1);
  L0_L0bar_pT1_vs_pT2_US_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L0_L0bar_pT1_vs_pT2_US_hist->GetXaxis()->CenterTitle();
  L0_L0bar_pT1_vs_pT2_US_hist->GetYaxis()->SetTitle("p_{T} GeV/c");
  L0_L0bar_pT1_vs_pT2_US_hist->GetYaxis()->CenterTitle();
  L0_L0bar_pT1_vs_pT2_US_hist->Draw("colz");

  L_Lbar_pT_data_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine_new/L_Lbar/L_Lbar_pT_data.png");


  TCanvas *L_Lbar_pT_PYTHIA_can = new TCanvas("L_Lbar_pT_PYTHIA_can", "L_Lbar_pT_PYTHIA_can", 1200, 1000);
  L_Lbar_pT_PYTHIA_can->cd();

  L0_L0bar_pT1_vs_pT2_cuts_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L0_L0bar_pT1_vs_pT2_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0bar_pT1_vs_pT2_cuts_hist->GetYaxis()->SetTitle("p_{T} GeV/c");
  L0_L0bar_pT1_vs_pT2_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0bar_pT1_vs_pT2_cuts_hist->Draw("colz");

  L_Lbar_pT_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine_new/L_Lbar/L_Lbar_pT_PYTHIA.png");


  TCanvas *L_Lbar_pT_PYTHIA_ME_can = new TCanvas("L_Lbar_pT_PYTHIA_ME_can", "L_Lbar_pT_PYTHIA_ME_can", 1200, 1000);
  L_Lbar_pT_PYTHIA_ME_can->cd();

  L0_L0bar_pT1_vs_pT2_ME_cuts_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L0_L0bar_pT1_vs_pT2_ME_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0bar_pT1_vs_pT2_ME_cuts_hist->GetYaxis()->SetTitle("p_{T} GeV/c");
  L0_L0bar_pT1_vs_pT2_ME_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0bar_pT1_vs_pT2_ME_cuts_hist->Draw("colz");

  L_Lbar_pT_PYTHIA_ME_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine_new/L_Lbar/L_Lbar_pT_PYTHIA_ME.png");


  TCanvas *L_Lbar_pT_weight_can = new TCanvas("L_Lbar_pT_weight_can", "L_Lbar_pT_weight_can", 1200, 1000);
  L_Lbar_pT_weight_can->cd();

  L0_L0bar_pT1_vs_pT2_cuts_hist->Scale(1./L0_L0bar_pT1_vs_pT2_cuts_hist->Integral());

  TH2F *L0_L0bar_pT1_vs_pT2_weight_hist = (TH2F*)L0_L0bar_pT1_vs_pT2_US_hist->Clone("L0_L0bar_pT1_vs_pT2_weight_hist");
  L0_L0bar_pT1_vs_pT2_weight_hist->Scale(1./L0_L0bar_pT1_vs_pT2_weight_hist->Integral());
  L0_L0bar_pT1_vs_pT2_weight_hist->Divide(L0_L0bar_pT1_vs_pT2_cuts_hist);
  L0_L0bar_pT1_vs_pT2_weight_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L0_L0bar_pT1_vs_pT2_weight_hist->GetXaxis()->CenterTitle();
  L0_L0bar_pT1_vs_pT2_weight_hist->GetYaxis()->SetTitle("p_{T} GeV/c");
  L0_L0bar_pT1_vs_pT2_weight_hist->GetYaxis()->CenterTitle();
  L0_L0bar_pT1_vs_pT2_weight_hist->Draw("colz");
  L0_L0bar_pT1_vs_pT2_weight_hist->Write();

  L_Lbar_pT_weight_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine_new/L_Lbar/L_Lbar_pT_weight.png");

  //---------------------------------------------------------------

  TCanvas *L_Lbar_eta_data_can = new TCanvas("L_Lbar_eta_data_can", "L_Lbar_eta_data_can", 1200, 1000);
  L_Lbar_eta_data_can->cd();

  L0_L0bar_eta1_vs_eta2_US_hist->Add(L0_L0bar_eta1_vs_eta2_US_LS_hist, -1);
  L0_L0bar_eta1_vs_eta2_US_hist->GetXaxis()->SetTitle("#eta");
  L0_L0bar_eta1_vs_eta2_US_hist->GetXaxis()->CenterTitle();
  L0_L0bar_eta1_vs_eta2_US_hist->GetYaxis()->SetTitle("#eta");
  L0_L0bar_eta1_vs_eta2_US_hist->GetYaxis()->CenterTitle();
  L0_L0bar_eta1_vs_eta2_US_hist->Draw("colz");

  L_Lbar_eta_data_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine_new/L_Lbar/L_Lbar_eta_data.png");


  TCanvas *L_Lbar_eta_PYTHIA_can = new TCanvas("L_Lbar_eta_PYTHIA_can", "L_Lbar_eta_PYTHIA_can", 1200, 1000);
  L_Lbar_eta_PYTHIA_can->cd();

  L0_L0bar_eta1_vs_eta2_cuts_hist->GetXaxis()->SetTitle("y");
  L0_L0bar_eta1_vs_eta2_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0bar_eta1_vs_eta2_cuts_hist->GetYaxis()->SetTitle("y");
  L0_L0bar_eta1_vs_eta2_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0bar_eta1_vs_eta2_cuts_hist->Draw("colz");

  L_Lbar_eta_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine_new/L_Lbar/L_Lbar_eta_PYTHIA.png");


  TCanvas *L_Lbar_eta_PYTHIA_ME_can = new TCanvas("L_Lbar_eta_PYTHIA_ME_can", "L_Lbar_eta_PYTHIA_ME_can", 1200, 1000);
  L_Lbar_eta_PYTHIA_ME_can->cd();

  L0_L0bar_eta1_vs_eta2_ME_cuts_hist->GetXaxis()->SetTitle("y");
  L0_L0bar_eta1_vs_eta2_ME_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0bar_eta1_vs_eta2_ME_cuts_hist->GetYaxis()->SetTitle("y");
  L0_L0bar_eta1_vs_eta2_ME_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0bar_eta1_vs_eta2_ME_cuts_hist->Draw("colz");

  L_Lbar_eta_PYTHIA_ME_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine_new/L_Lbar/L_Lbar_eta_PYTHIA_ME.png");


  TCanvas *L_Lbar_eta_weight_can = new TCanvas("L_Lbar_eta_weight_can", "L_Lbar_eta_weight_can", 1200, 1000);
  L_Lbar_eta_weight_can->cd();

  L0_L0bar_eta1_vs_eta2_cuts_hist->Scale(1./L0_L0bar_eta1_vs_eta2_cuts_hist->Integral());

  TH2F *L0_L0bar_eta1_vs_eta2_weight_hist = (TH2F*)L0_L0bar_eta1_vs_eta2_US_hist->Clone("L0_L0bar_eta1_vs_eta2_weight_hist");
  L0_L0bar_eta1_vs_eta2_weight_hist->Scale(1./L0_L0bar_eta1_vs_eta2_weight_hist->Integral());
  L0_L0bar_eta1_vs_eta2_weight_hist->Divide(L0_L0bar_eta1_vs_eta2_cuts_hist);
  L0_L0bar_eta1_vs_eta2_weight_hist->GetXaxis()->SetTitle("#eta");
  L0_L0bar_eta1_vs_eta2_weight_hist->GetXaxis()->CenterTitle();
  L0_L0bar_eta1_vs_eta2_weight_hist->GetYaxis()->SetTitle("#eta");
  L0_L0bar_eta1_vs_eta2_weight_hist->GetYaxis()->CenterTitle();
  L0_L0bar_eta1_vs_eta2_weight_hist->Draw("colz");
  L0_L0bar_eta1_vs_eta2_weight_hist->Write();

  L_Lbar_eta_weight_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine_new/L_Lbar/L_Lbar_eta_weight.png");

  //---------------------------------------------------------------

  TCanvas *L_Lbar_phi_data_can = new TCanvas("L_Lbar_phi_data_can", "L_Lbar_phi_data_can", 1200, 1000);
  L_Lbar_phi_data_can->cd();

  L0_L0bar_phi1_vs_phi2_US_hist->Add(L0_L0bar_phi1_vs_phi2_US_LS_hist, -1);
  L0_L0bar_phi1_vs_phi2_US_hist->GetXaxis()->SetTitle("#phi");
  L0_L0bar_phi1_vs_phi2_US_hist->GetXaxis()->CenterTitle();
  L0_L0bar_phi1_vs_phi2_US_hist->GetYaxis()->SetTitle("#phi");
  L0_L0bar_phi1_vs_phi2_US_hist->GetYaxis()->CenterTitle();
  L0_L0bar_phi1_vs_phi2_US_hist->Draw("colz");

  L_Lbar_phi_data_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine_new/L_Lbar/L_Lbar_phi_data.png");


  TCanvas *L_Lbar_phi_PYTHIA_can = new TCanvas("L_Lbar_phi_PYTHIA_can", "L_Lbar_phi_PYTHIA_can", 1200, 1000);
  L_Lbar_phi_PYTHIA_can->cd();

  L0_L0bar_phi1_vs_phi2_cuts_hist->GetXaxis()->SetTitle("#phi");
  L0_L0bar_phi1_vs_phi2_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0bar_phi1_vs_phi2_cuts_hist->GetYaxis()->SetTitle("#phi");
  L0_L0bar_phi1_vs_phi2_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0bar_phi1_vs_phi2_cuts_hist->Draw("colz");

  L_Lbar_phi_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine_new/L_Lbar/L_Lbar_phi_PYTHIA.png");


  TCanvas *L_Lbar_phi_PYTHIA_ME_can = new TCanvas("L_Lbar_phi_PYTHIA_ME_can", "L_Lbar_phi_PYTHIA_ME_can", 1200, 1000);
  L_Lbar_phi_PYTHIA_ME_can->cd();

  L0_L0bar_phi1_vs_phi2_ME_cuts_hist->GetXaxis()->SetTitle("y");
  L0_L0bar_phi1_vs_phi2_ME_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0bar_phi1_vs_phi2_ME_cuts_hist->GetYaxis()->SetTitle("y");
  L0_L0bar_phi1_vs_phi2_ME_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0bar_phi1_vs_phi2_ME_cuts_hist->Draw("colz");

  L_Lbar_phi_PYTHIA_ME_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine_new/L_Lbar/L_Lbar_phi_PYTHIA_ME.png");


  TCanvas *L_Lbar_phi_weight_can = new TCanvas("L_Lbar_phi_weight_can", "L_Lbar_phi_weight_can", 1200, 1000);
  L_Lbar_phi_weight_can->cd();

  L0_L0bar_phi1_vs_phi2_cuts_hist->Scale(1./L0_L0bar_phi1_vs_phi2_cuts_hist->Integral());

  TH2F *L0_L0bar_phi1_vs_phi2_weight_hist = (TH2F*)L0_L0bar_phi1_vs_phi2_US_hist->Clone("L0_L0bar_phi1_vs_phi2_weight_hist");
  L0_L0bar_phi1_vs_phi2_weight_hist->Scale(1./L0_L0bar_phi1_vs_phi2_weight_hist->Integral());
  L0_L0bar_phi1_vs_phi2_weight_hist->Divide(L0_L0bar_phi1_vs_phi2_cuts_hist);
  L0_L0bar_phi1_vs_phi2_weight_hist->GetXaxis()->SetTitle("#phi");
  L0_L0bar_phi1_vs_phi2_weight_hist->GetXaxis()->CenterTitle();
  L0_L0bar_phi1_vs_phi2_weight_hist->GetYaxis()->SetTitle("#phi");
  L0_L0bar_phi1_vs_phi2_weight_hist->GetYaxis()->CenterTitle();
  L0_L0bar_phi1_vs_phi2_weight_hist->Draw("colz");
  L0_L0bar_phi1_vs_phi2_weight_hist->Write();

  L_Lbar_phi_weight_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine_new/L_Lbar/L_Lbar_phi_weight.png");

  //--------------------

  TCanvas *L_Lbar_delta_phi_can = new TCanvas("L_Lbar_delta_phi_can", "L_Lbar_delta_phi_can", 1200, 1000);
  L_Lbar_delta_phi_can->cd();

  L0_L0bar_delta_phi_cuts_hist->SetLineColor(1);
  L0_L0bar_delta_phi_cuts_hist->GetXaxis()->SetTitle("#Delta#phi");
  L0_L0bar_delta_phi_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0bar_delta_phi_cuts_hist->SetMinimum(0);
  L0_L0bar_delta_phi_cuts_hist->Draw("hist");

  TPaveText *text_LLbar = new TPaveText(0.4, 0.6, 0.8, 0.8, "NDC");
  //cent_text->AddText("STAR preliminary");
  //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
  text_LLbar->SetTextFont(43);
  text_LLbar->SetTextSize(33);
  text_LLbar->AddText(Form("MC p+p #sqrt{s} = %i GeV", 200));
  text_LLbar->AddText("Minimum bias");
  text_LLbar->AddText("#Lambda#bar{#Lambda}");
  text_LLbar->SetFillColorAlpha(0, 0.01);
  text_LLbar->Draw("same");

  L_Lbar_delta_phi_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine_new/L_Lbar/L_Lbar_delta_phi.png");


  TCanvas *L_Lbar_delta_phi_ME_can = new TCanvas("L_Lbar_delta_phi_ME_can", "L_Lbar_delta_phi_ME_can", 1200, 1000);
  L_Lbar_delta_phi_ME_can->cd();

  L0_L0bar_delta_phi_ME_cuts_hist->SetLineColor(1);
  L0_L0bar_delta_phi_ME_cuts_hist->GetXaxis()->SetTitle("#Delta#phi");
  L0_L0bar_delta_phi_ME_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0bar_delta_phi_ME_cuts_hist->SetMinimum(0);
  L0_L0bar_delta_phi_ME_cuts_hist->Draw("hist");

  L_Lbar_delta_phi_ME_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine_new/L_Lbar/L_Lbar_delta_phi_ME.png");

  //--------------------

  TCanvas *L_Lbar_delta_eta_can = new TCanvas("L_Lbar_delta_eta_can", "L_Lbar_delta_eta_can", 1200, 1000);
  L_Lbar_delta_eta_can->cd();

  L0_L0bar_delta_eta_cuts_hist->SetLineColor(1);
  L0_L0bar_delta_eta_cuts_hist->GetXaxis()->SetTitle("#Deltay");
  L0_L0bar_delta_eta_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0bar_delta_eta_cuts_hist->SetMinimum(0);
  L0_L0bar_delta_eta_cuts_hist->Draw("hist");

  text_LLbar->Draw("same");

  L_Lbar_delta_eta_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine_new/L_Lbar/L_Lbar_delta_eta.png");


  TCanvas *L_Lbar_delta_eta_ME_can = new TCanvas("L_Lbar_delta_eta_ME_can", "L_Lbar_delta_eta_ME_can", 1200, 1000);
  L_Lbar_delta_eta_ME_can->cd();

  L0_L0bar_delta_eta_ME_cuts_hist->SetLineColor(1);
  L0_L0bar_delta_eta_ME_cuts_hist->GetXaxis()->SetTitle("#Deltay");
  L0_L0bar_delta_eta_ME_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0bar_delta_eta_ME_cuts_hist->SetMinimum(0);
  L0_L0bar_delta_eta_ME_cuts_hist->Draw("hist");

  L_Lbar_delta_eta_ME_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine_new/L_Lbar/L_Lbar_delta_eta_ME.png");




  //--------------------------------------------------------------------------------------------------------------------------------------

  //L-L kinematics

  TCanvas *L_L_pT_data_can = new TCanvas("L_L_pT_data_can", "L_L_pT_data_can", 1200, 1000);
  L_L_pT_data_can->cd();

  L0_L0_pT1_vs_pT2_US_hist->Add(L0_L0_pT1_vs_pT2_US_LS_hist, -1);
  L0_L0_pT1_vs_pT2_US_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L0_L0_pT1_vs_pT2_US_hist->GetXaxis()->CenterTitle();
  L0_L0_pT1_vs_pT2_US_hist->GetYaxis()->SetTitle("p_{T} GeV/c");
  L0_L0_pT1_vs_pT2_US_hist->GetYaxis()->CenterTitle();
  L0_L0_pT1_vs_pT2_US_hist->Draw("colz");

  L_L_pT_data_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine_new/L_L/L_L_pT_data.png");


  TCanvas *L_L_pT_PYTHIA_can = new TCanvas("L_L_pT_PYTHIA_can", "L_L_pT_PYTHIA_can", 1200, 1000);
  L_L_pT_PYTHIA_can->cd();

  L0_L0_pT1_vs_pT2_cuts_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L0_L0_pT1_vs_pT2_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0_pT1_vs_pT2_cuts_hist->GetYaxis()->SetTitle("p_{T} GeV/c");
  L0_L0_pT1_vs_pT2_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0_pT1_vs_pT2_cuts_hist->Draw("colz");

  L_L_pT_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine_new/L_L/L_L_pT_PYTHIA.png");


  TCanvas *L_L_pT_weight_can = new TCanvas("L_L_pT_weight_can", "L_L_pT_weight_can", 1200, 1000);
  L_L_pT_weight_can->cd();

  L0_L0_pT1_vs_pT2_cuts_hist->Scale(1./L0_L0_pT1_vs_pT2_cuts_hist->Integral());

  TH2F *L0_L0_pT1_vs_pT2_weight_hist = (TH2F*)L0_L0_pT1_vs_pT2_US_hist->Clone("L0_L0_pT1_vs_pT2_weight_hist");
  L0_L0_pT1_vs_pT2_weight_hist->Scale(1./L0_L0_pT1_vs_pT2_weight_hist->Integral());
  L0_L0_pT1_vs_pT2_weight_hist->Divide(L0_L0_pT1_vs_pT2_cuts_hist);
  L0_L0_pT1_vs_pT2_weight_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L0_L0_pT1_vs_pT2_weight_hist->GetXaxis()->CenterTitle();
  L0_L0_pT1_vs_pT2_weight_hist->GetYaxis()->SetTitle("p_{T} GeV/c");
  L0_L0_pT1_vs_pT2_weight_hist->GetYaxis()->CenterTitle();
  L0_L0_pT1_vs_pT2_weight_hist->Draw("colz");
  L0_L0_pT1_vs_pT2_weight_hist->Write();

  L_L_pT_weight_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine_new/L_L/L_L_pT_weight.png");

  //---------------------------------------------------------------

  TCanvas *L_L_eta_data_can = new TCanvas("L_L_eta_data_can", "L_L_eta_data_can", 1200, 1000);
  L_L_eta_data_can->cd();

  L0_L0_eta1_vs_eta2_US_hist->Add(L0_L0_eta1_vs_eta2_US_LS_hist, -1);
  L0_L0_eta1_vs_eta2_US_hist->GetXaxis()->SetTitle("#eta");
  L0_L0_eta1_vs_eta2_US_hist->GetXaxis()->CenterTitle();
  L0_L0_eta1_vs_eta2_US_hist->GetYaxis()->SetTitle("#eta");
  L0_L0_eta1_vs_eta2_US_hist->GetYaxis()->CenterTitle();
  L0_L0_eta1_vs_eta2_US_hist->Draw("colz");

  L_L_eta_data_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine_new/L_L/L_L_eta_data.png");


  TCanvas *L_L_eta_PYTHIA_can = new TCanvas("L_L_eta_PYTHIA_can", "L_L_eta_PYTHIA_can", 1200, 1000);
  L_L_eta_PYTHIA_can->cd();

  L0_L0_eta1_vs_eta2_cuts_hist->GetXaxis()->SetTitle("#eta");
  L0_L0_eta1_vs_eta2_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0_eta1_vs_eta2_cuts_hist->GetYaxis()->SetTitle("#eta");
  L0_L0_eta1_vs_eta2_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0_eta1_vs_eta2_cuts_hist->Draw("colz");

  L_L_eta_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine_new/L_L/L_L_eta_PYTHIA.png");


  TCanvas *L_L_eta_weight_can = new TCanvas("L_L_eta_weight_can", "L_L_eta_weight_can", 1200, 1000);
  L_L_eta_weight_can->cd();

  L0_L0_eta1_vs_eta2_cuts_hist->Scale(1./L0_L0_eta1_vs_eta2_cuts_hist->Integral());

  TH2F *L0_L0_eta1_vs_eta2_weight_hist = (TH2F*)L0_L0_eta1_vs_eta2_US_hist->Clone("L0_L0_eta1_vs_eta2_weight_hist");
  L0_L0_eta1_vs_eta2_weight_hist->Scale(1./L0_L0_eta1_vs_eta2_weight_hist->Integral());
  L0_L0_eta1_vs_eta2_weight_hist->Divide(L0_L0_eta1_vs_eta2_cuts_hist);
  L0_L0_eta1_vs_eta2_weight_hist->GetXaxis()->SetTitle("#eta");
  L0_L0_eta1_vs_eta2_weight_hist->GetXaxis()->CenterTitle();
  L0_L0_eta1_vs_eta2_weight_hist->GetYaxis()->SetTitle("#eta");
  L0_L0_eta1_vs_eta2_weight_hist->GetYaxis()->CenterTitle();
  L0_L0_eta1_vs_eta2_weight_hist->Draw("colz");
  L0_L0_eta1_vs_eta2_weight_hist->Write();

  L_L_eta_weight_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine_new/L_L/L_L_eta_weight.png");

  //---------------------------------------------------------------

  TCanvas *L_L_phi_data_can = new TCanvas("L_L_phi_data_can", "L_L_phi_data_can", 1200, 1000);
  L_L_phi_data_can->cd();

  L0_L0_phi1_vs_phi2_US_hist->Add(L0_L0_phi1_vs_phi2_US_LS_hist, -1);
  L0_L0_phi1_vs_phi2_US_hist->GetXaxis()->SetTitle("#phi");
  L0_L0_phi1_vs_phi2_US_hist->GetXaxis()->CenterTitle();
  L0_L0_phi1_vs_phi2_US_hist->GetYaxis()->SetTitle("#phi");
  L0_L0_phi1_vs_phi2_US_hist->GetYaxis()->CenterTitle();
  L0_L0_phi1_vs_phi2_US_hist->Draw("colz");

  L_L_phi_data_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine_new/L_L/L_L_phi_data.png");


  TCanvas *L_L_phi_PYTHIA_can = new TCanvas("L_L_phi_PYTHIA_can", "L_L_phi_PYTHIA_can", 1200, 1000);
  L_L_phi_PYTHIA_can->cd();

  L0_L0_phi1_vs_phi2_cuts_hist->GetXaxis()->SetTitle("#phi");
  L0_L0_phi1_vs_phi2_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0_phi1_vs_phi2_cuts_hist->GetYaxis()->SetTitle("#phi");
  L0_L0_phi1_vs_phi2_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0_phi1_vs_phi2_cuts_hist->Draw("colz");

  L_L_phi_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine_new/L_L/L_L_phi_PYTHIA.png");


  TCanvas *L_L_phi_weight_can = new TCanvas("L_L_phi_weight_can", "L_L_phi_weight_can", 1200, 1000);
  L_L_phi_weight_can->cd();

  L0_L0_phi1_vs_phi2_cuts_hist->Scale(1./L0_L0_phi1_vs_phi2_cuts_hist->Integral());

  TH2F *L0_L0_phi1_vs_phi2_weight_hist = (TH2F*)L0_L0_phi1_vs_phi2_US_hist->Clone("L0_L0_phi1_vs_phi2_weight_hist");
  L0_L0_phi1_vs_phi2_weight_hist->Scale(1./L0_L0_phi1_vs_phi2_weight_hist->Integral());
  L0_L0_phi1_vs_phi2_weight_hist->Divide(L0_L0_phi1_vs_phi2_cuts_hist);
  L0_L0_phi1_vs_phi2_weight_hist->GetXaxis()->SetTitle("#phi");
  L0_L0_phi1_vs_phi2_weight_hist->GetXaxis()->CenterTitle();
  L0_L0_phi1_vs_phi2_weight_hist->GetYaxis()->SetTitle("#phi");
  L0_L0_phi1_vs_phi2_weight_hist->GetYaxis()->CenterTitle();
  L0_L0_phi1_vs_phi2_weight_hist->Draw("colz");
  L0_L0_phi1_vs_phi2_weight_hist->Write();

  L_L_phi_weight_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine_new/L_L/L_L_phi_weight.png");


  TCanvas *L_L_delta_phi_can = new TCanvas("L_L_delta_phi_can", "L_L_delta_phi_can", 1200, 1000);
  L_L_delta_phi_can->cd();

  L0_L0_delta_phi_cuts_hist->SetLineColor(1);
  L0_L0_delta_phi_cuts_hist->GetXaxis()->SetTitle("#Delta#phi");
  L0_L0_delta_phi_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0_delta_phi_cuts_hist->SetMinimum(0);
  L0_L0_delta_phi_cuts_hist->Draw("hist");

  L_L_delta_phi_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine_new/L_L/L_L_delta_phi.png");

  //--------------------------------------------------------------------------------------------------------------------------------------

  //Lbar-Lbar kinematics

  TCanvas *Lbar_Lbar_pT_data_can = new TCanvas("Lbar_Lbar_pT_data_can", "Lbar_Lbar_pT_data_can", 1200, 1000);
  Lbar_Lbar_pT_data_can->cd();

  L0bar_L0bar_pT1_vs_pT2_US_hist->Add(L0bar_L0bar_pT1_vs_pT2_US_LS_hist, -1);
  L0bar_L0bar_pT1_vs_pT2_US_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L0bar_L0bar_pT1_vs_pT2_US_hist->GetXaxis()->CenterTitle();
  L0bar_L0bar_pT1_vs_pT2_US_hist->GetYaxis()->SetTitle("p_{T} GeV/c");
  L0bar_L0bar_pT1_vs_pT2_US_hist->GetYaxis()->CenterTitle();
  L0bar_L0bar_pT1_vs_pT2_US_hist->Draw("colz");

  Lbar_Lbar_pT_data_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine_new/Lbar_Lbar/Lbar_Lbar_pT_data.png");


  TCanvas *Lbar_Lbar_pT_PYTHIA_can = new TCanvas("Lbar_Lbar_pT_PYTHIA_can", "Lbar_Lbar_pT_PYTHIA_can", 1200, 1000);
  Lbar_Lbar_pT_PYTHIA_can->cd();

  L0bar_L0bar_pT1_vs_pT2_cuts_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L0bar_L0bar_pT1_vs_pT2_cuts_hist->GetXaxis()->CenterTitle();
  L0bar_L0bar_pT1_vs_pT2_cuts_hist->GetYaxis()->SetTitle("p_{T} GeV/c");
  L0bar_L0bar_pT1_vs_pT2_cuts_hist->GetYaxis()->CenterTitle();
  L0bar_L0bar_pT1_vs_pT2_cuts_hist->Draw("colz");

  Lbar_Lbar_pT_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine_new/Lbar_Lbar/Lbar_Lbar_pT_PYTHIA.png");


  TCanvas *Lbar_Lbar_pT_weight_can = new TCanvas("Lbar_Lbar_pT_weight_can", "Lbar_Lbar_pT_weight_can", 1200, 1000);
  Lbar_Lbar_pT_weight_can->cd();

  L0bar_L0bar_pT1_vs_pT2_cuts_hist->Scale(1./L0bar_L0bar_pT1_vs_pT2_cuts_hist->Integral());

  TH2F *L0bar_L0bar_pT1_vs_pT2_weight_hist = (TH2F*)L0bar_L0bar_pT1_vs_pT2_US_hist->Clone("L0bar_L0bar_pT1_vs_pT2_weight_hist");
  L0bar_L0bar_pT1_vs_pT2_weight_hist->Scale(1./L0bar_L0bar_pT1_vs_pT2_weight_hist->Integral());
  L0bar_L0bar_pT1_vs_pT2_weight_hist->Divide(L0bar_L0bar_pT1_vs_pT2_cuts_hist);
  L0bar_L0bar_pT1_vs_pT2_weight_hist->GetXaxis()->SetTitle("p_{T} GeV/c");
  L0bar_L0bar_pT1_vs_pT2_weight_hist->GetXaxis()->CenterTitle();
  L0bar_L0bar_pT1_vs_pT2_weight_hist->GetYaxis()->SetTitle("p_{T} GeV/c");
  L0bar_L0bar_pT1_vs_pT2_weight_hist->GetYaxis()->CenterTitle();
  L0bar_L0bar_pT1_vs_pT2_weight_hist->Draw("colz");
  L0bar_L0bar_pT1_vs_pT2_weight_hist->Write();

  Lbar_Lbar_pT_weight_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine_new/Lbar_Lbar/Lbar_Lbar_pT_weight.png");

  //---------------------------------------------------------------

  TCanvas *Lbar_Lbar_eta_data_can = new TCanvas("Lbar_Lbar_eta_data_can", "Lbar_Lbar_eta_data_can", 1200, 1000);
  Lbar_Lbar_eta_data_can->cd();

  L0bar_L0bar_eta1_vs_eta2_US_hist->Add(L0bar_L0bar_eta1_vs_eta2_US_LS_hist, -1);
  L0bar_L0bar_eta1_vs_eta2_US_hist->GetXaxis()->SetTitle("#eta");
  L0bar_L0bar_eta1_vs_eta2_US_hist->GetXaxis()->CenterTitle();
  L0bar_L0bar_eta1_vs_eta2_US_hist->GetYaxis()->SetTitle("#eta");
  L0bar_L0bar_eta1_vs_eta2_US_hist->GetYaxis()->CenterTitle();
  L0bar_L0bar_eta1_vs_eta2_US_hist->Draw("colz");

  Lbar_Lbar_eta_data_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine_new/Lbar_Lbar/Lbar_Lbar_eta_data.png");


  TCanvas *Lbar_Lbar_eta_PYTHIA_can = new TCanvas("Lbar_Lbar_eta_PYTHIA_can", "Lbar_Lbar_eta_PYTHIA_can", 1200, 1000);
  Lbar_Lbar_eta_PYTHIA_can->cd();

  L0bar_L0bar_eta1_vs_eta2_cuts_hist->GetXaxis()->SetTitle("#eta");
  L0bar_L0bar_eta1_vs_eta2_cuts_hist->GetXaxis()->CenterTitle();
  L0bar_L0bar_eta1_vs_eta2_cuts_hist->GetYaxis()->SetTitle("#eta");
  L0bar_L0bar_eta1_vs_eta2_cuts_hist->GetYaxis()->CenterTitle();
  L0bar_L0bar_eta1_vs_eta2_cuts_hist->Draw("colz");

  Lbar_Lbar_eta_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine_new/Lbar_Lbar/Lbar_Lbar_eta_PYTHIA.png");


  TCanvas *Lbar_Lbar_eta_weight_can = new TCanvas("Lbar_Lbar_eta_weight_can", "Lbar_Lbar_eta_weight_can", 1200, 1000);
  Lbar_Lbar_eta_weight_can->cd();

  L0bar_L0bar_eta1_vs_eta2_cuts_hist->Scale(1./L0bar_L0bar_eta1_vs_eta2_cuts_hist->Integral());

  TH2F *L0bar_L0bar_eta1_vs_eta2_weight_hist = (TH2F*)L0bar_L0bar_eta1_vs_eta2_US_hist->Clone("L0bar_L0bar_eta1_vs_eta2_weight_hist");
  L0bar_L0bar_eta1_vs_eta2_weight_hist->Scale(1./L0bar_L0bar_eta1_vs_eta2_weight_hist->Integral());
  L0bar_L0bar_eta1_vs_eta2_weight_hist->Divide(L0bar_L0bar_eta1_vs_eta2_cuts_hist);
  L0bar_L0bar_eta1_vs_eta2_weight_hist->GetXaxis()->SetTitle("#eta");
  L0bar_L0bar_eta1_vs_eta2_weight_hist->GetXaxis()->CenterTitle();
  L0bar_L0bar_eta1_vs_eta2_weight_hist->GetYaxis()->SetTitle("#eta");
  L0bar_L0bar_eta1_vs_eta2_weight_hist->GetYaxis()->CenterTitle();
  L0bar_L0bar_eta1_vs_eta2_weight_hist->Draw("colz");
  L0bar_L0bar_eta1_vs_eta2_weight_hist->Write();

  Lbar_Lbar_eta_weight_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine_new/Lbar_Lbar/Lbar_Lbar_eta_weight.png");

  //---------------------------------------------------------------

  TCanvas *Lbar_Lbar_phi_data_can = new TCanvas("Lbar_Lbar_phi_data_can", "Lbar_Lbar_phi_data_can", 1200, 1000);
  Lbar_Lbar_phi_data_can->cd();

  L0bar_L0bar_phi1_vs_phi2_US_hist->Add(L0bar_L0bar_phi1_vs_phi2_US_LS_hist, -1);
  L0bar_L0bar_phi1_vs_phi2_US_hist->GetXaxis()->SetTitle("#phi");
  L0bar_L0bar_phi1_vs_phi2_US_hist->GetXaxis()->CenterTitle();
  L0bar_L0bar_phi1_vs_phi2_US_hist->GetYaxis()->SetTitle("#phi");
  L0bar_L0bar_phi1_vs_phi2_US_hist->GetYaxis()->CenterTitle();
  L0bar_L0bar_phi1_vs_phi2_US_hist->Draw("colz");

  Lbar_Lbar_phi_data_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine_new/Lbar_Lbar/Lbar_Lbar_phi_data.png");


  TCanvas *Lbar_Lbar_phi_PYTHIA_can = new TCanvas("Lbar_Lbar_phi_PYTHIA_can", "Lbar_Lbar_phi_PYTHIA_can", 1200, 1000);
  Lbar_Lbar_phi_PYTHIA_can->cd();

  L0bar_L0bar_phi1_vs_phi2_cuts_hist->GetXaxis()->SetTitle("#phi");
  L0bar_L0bar_phi1_vs_phi2_cuts_hist->GetXaxis()->CenterTitle();
  L0bar_L0bar_phi1_vs_phi2_cuts_hist->GetYaxis()->SetTitle("#phi");
  L0bar_L0bar_phi1_vs_phi2_cuts_hist->GetYaxis()->CenterTitle();
  L0bar_L0bar_phi1_vs_phi2_cuts_hist->Draw("colz");

  Lbar_Lbar_phi_PYTHIA_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine_new/Lbar_Lbar/Lbar_Lbar_phi_PYTHIA.png");


  TCanvas *Lbar_Lbar_phi_weight_can = new TCanvas("Lbar_Lbar_phi_weight_can", "Lbar_Lbar_phi_weight_can", 1200, 1000);
  Lbar_Lbar_phi_weight_can->cd();

  L0bar_L0bar_phi1_vs_phi2_cuts_hist->Scale(1./L0bar_L0bar_phi1_vs_phi2_cuts_hist->Integral());

  TH2F *L0bar_L0bar_phi1_vs_phi2_weight_hist = (TH2F*)L0bar_L0bar_phi1_vs_phi2_US_hist->Clone("L0bar_L0bar_phi1_vs_phi2_weight_hist");
  L0bar_L0bar_phi1_vs_phi2_weight_hist->Scale(1./L0bar_L0bar_phi1_vs_phi2_weight_hist->Integral());
  L0bar_L0bar_phi1_vs_phi2_weight_hist->Divide(L0bar_L0bar_phi1_vs_phi2_cuts_hist);
  L0bar_L0bar_phi1_vs_phi2_weight_hist->GetXaxis()->SetTitle("#phi");
  L0bar_L0bar_phi1_vs_phi2_weight_hist->GetXaxis()->CenterTitle();
  L0bar_L0bar_phi1_vs_phi2_weight_hist->GetYaxis()->SetTitle("#phi");
  L0bar_L0bar_phi1_vs_phi2_weight_hist->GetYaxis()->CenterTitle();
  L0bar_L0bar_phi1_vs_phi2_weight_hist->Draw("colz");
  L0bar_L0bar_phi1_vs_phi2_weight_hist->Write();

  Lbar_Lbar_phi_weight_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine_new/Lbar_Lbar/Lbar_Lbar_phi_weight.png");


  TCanvas *Lbar_Lbar_delta_phi_can = new TCanvas("Lbar_Lbar_delta_phi_can", "Lbar_Lbar_delta_phi_can", 1200, 1000);
  Lbar_Lbar_delta_phi_can->cd();

  L0bar_L0bar_delta_phi_cuts_hist->SetLineColor(1);
  L0bar_L0bar_delta_phi_cuts_hist->GetXaxis()->SetTitle("#Delta#phi");
  L0bar_L0bar_delta_phi_cuts_hist->GetXaxis()->CenterTitle();
  L0bar_L0bar_delta_phi_cuts_hist->SetMinimum(0);
  L0bar_L0bar_delta_phi_cuts_hist->Draw("hist");

  Lbar_Lbar_delta_phi_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine_new/Lbar_Lbar/Lbar_Lbar_delta_phi.png");

  //________________________________________________________________________________________________________________________________________________________________________________________________


  inFile_PYTHIA->Close();
  inFile_data->Close();
  inFile_data_Minv->Close();

  outFile_weight->Close();

  return;

}
