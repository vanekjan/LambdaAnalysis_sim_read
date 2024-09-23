#include<iostream>
#include<fstream>
#include<vector>
#include"TH1.h"
#include"TH2.h"
#include"TF1.h"
#include"TGraphErrors.h"
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
//const int corr_err = 0 - no stat. err. corection, 1 - correct sys. err. for statistical precision
void Closure_test_new_ME(const int energy = 200, const int corr_err = 0)
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

    //inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/output_Lambda_pp_200_MB_1B_events_hists_work.root", "READ");

    //inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/2024_07_ME_weight/output_Lambda_pp_200_MB_1B_events_hists_ME_weight_1.root", "READ");

    inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/2024_08_ME_weight_new/output_Lambda_pp_200_MB_1B_events_hists_Delta_phi_pi_third.root", "READ");
    //inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/2024_08_ME_weight_new/output_Lambda_pp_200_MB_1B_events_hists_Delta_phi_pi_quater.root", "READ");

    //momentum smearing
    //inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/mom_smear_new/output_Lambda_pp_200_MB_1B_events_hists_pair_weight_mom_smear_full.root", "READ");
    //inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/mom_smear_new/output_Lambda_pp_200_MB_1B_events_hists_pair_weight_no_smear_full.root", "READ");

  }
  else
  {
    cout<<"Not a valid collision energy! Abborting!"<<endl;
    return;
  }


  TFile *sysErrFile;

  if(corr_err == 0) sysErrFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/output/sys_err/SysErrSlope_nocorr.root", "recreate");
  else if(corr_err == 1) sysErrFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/output/sys_err/SysErrSlope_corr.root", "recreate");
  else
  {
    cout<<"Wrong second argument"<<endl;
    return;
  }


  //histogram to store systematic error from slope difference in closure test
  //each bin is for one L charge combination: bin 1 - LLbar, bin 2 - LL, bin 3 - LbarLbar
  TH1F *SysErrSlope_hist = new TH1F("SysErrSlope_hist", "SysErrSlope_hist", 3, 0, 3);

  TH1F *SysErrSlope_delta_eta_hist[2];
  TH1F *ResidualPolarization_delta_eta_hist[2];

  TH1F *SysErrSlope_delta_phi_hist[2];
  TH1F *ResidualPolarization_delta_phi_hist[2];

  TH1F *SysErrSlope_delta_eta_delta_phi_hist[2];
  TH1F *ResidualPolarization_delta_eta_delta_phi_hist[2];

  for( unsigned int delta_eta_bin = 0; delta_eta_bin < 2; delta_eta_bin++ )
  {
    SysErrSlope_delta_eta_hist[delta_eta_bin] = new TH1F(Form("SysErrSlope_delta_eta_hist_%i", delta_eta_bin), Form("SysErrSlope_delta_eta_hist_%i", delta_eta_bin), 3, 0, 3);
    ResidualPolarization_delta_eta_hist[delta_eta_bin] = new TH1F(Form("ResidualPolarization_delta_eta_hist_%i", delta_eta_bin), Form("ResidualPolarization_delta_eta_hist_%i", delta_eta_bin), 3, 0, 3);

    SysErrSlope_delta_phi_hist[delta_eta_bin] = new TH1F(Form("SysErrSlope_delta_phi_hist_%i", delta_eta_bin), Form("SysErrSlope_delta_phi_hist_%i", delta_eta_bin), 3, 0, 3);
    ResidualPolarization_delta_phi_hist[delta_eta_bin] = new TH1F(Form("ResidualPolarization_delta_phi_hist_%i", delta_eta_bin), Form("ResidualPolarization_delta_phi_hist_%i", delta_eta_bin), 3, 0, 3);

    SysErrSlope_delta_eta_delta_phi_hist[delta_eta_bin] = new TH1F(Form("SysErrSlope_delta_eta_delta_phi_hist_%i", delta_eta_bin), Form("SysErrSlope_delta_eta_delta_phi_hist_%i", delta_eta_bin), 3, 0, 3);
    ResidualPolarization_delta_eta_delta_phi_hist[delta_eta_bin] = new TH1F(Form("ResidualPolarization_delta_eta_delta_phi_hist_%i", delta_eta_bin), Form("ResidualPolarization_delta_eta_delta_phi_hist_%i", delta_eta_bin), 3, 0, 3);
  }

  TH1F *SysErrSlope_delta_phi_less_hist[3];
  TH1F *SysErrSlope_delta_phi_more_hist[3];

  TH1F *ResidualPolarization_delta_phi_less_hist[3];
  TH1F *ResidualPolarization_delta_phi_more_hist[3];

  for(unsigned int delta_eta_bin = 0; delta_eta_bin < 3; delta_eta_bin++)
  {
    SysErrSlope_delta_phi_less_hist[delta_eta_bin] = new TH1F(Form("SysErrSlope_delta_phi_less_hist_%i", delta_eta_bin), Form("SysErrSlope_delta_phi_less_hist_%i", delta_eta_bin), 3, 0, 3);
    SysErrSlope_delta_phi_more_hist[delta_eta_bin] = new TH1F(Form("SysErrSlope_delta_phi_more_hist_%i", delta_eta_bin), Form("SysErrSlope_delta_phi_more_hist_%i", delta_eta_bin), 3, 0, 3);

    ResidualPolarization_delta_phi_less_hist[delta_eta_bin] = new TH1F(Form("ResidualPolarization_delta_phi_less_hist_%i", delta_eta_bin), Form("ResidualPolarization_delta_phi_less_hist_%i", delta_eta_bin), 3, 0, 3);
    ResidualPolarization_delta_phi_more_hist[delta_eta_bin] = new TH1F(Form("ResidualPolarization_delta_phi_more_hist_%i", delta_eta_bin), Form("ResidualPolarization_delta_phi_more_hist_%i", delta_eta_bin), 3, 0, 3);

  }


  TFile *out_file_polarization = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/output/polarization/Polarization_PYTHIA.root", "recreate");


  //histograms

  //before cuts
  //True MC
  TH1F *L0_L0bar_cosThetaProdPlane = (TH1F*)inFile->Get("L0_L0bar_cosThetaProdPlane");
  TH1F *L0_L0bar_cosThetaProdPlane_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0bar_cosThetaProdPlane_eta_hist[nEtaBins][nEtaBins];

  TH2F *L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_hist = (TH2F*)inFile->Get("L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_hist");
  TH2F *L0_L0bar_cos_theta_star_vs_delta_eta_hist = (TH2F*)inFile->Get("L0_L0bar_cos_theta_star_vs_delta_eta_hist");
  TH2F *L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_hist = (TH2F*)inFile->Get("L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_hist");



  TH1F *L0_L0_cosThetaProdPlane = (TH1F*)inFile->Get("L0_L0_cosThetaProdPlane");
  TH1F *L0_L0_cosThetaProdPlane_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0_cosThetaProdPlane_eta_hist[nEtaBins][nEtaBins];

  TH2F *L0_L0_cos_theta_star_vs_delta_phi_for_corr_hist = (TH2F*)inFile->Get("L0_L0_cos_theta_star_vs_delta_phi_for_corr_hist");
  TH2F *L0_L0_cos_theta_star_vs_delta_eta_hist = (TH2F*)inFile->Get("L0_L0_cos_theta_star_vs_delta_eta_hist");
  TH2F *L0_L0_cos_theta_star_vs_delta_eta_delta_phi_hist = (TH2F*)inFile->Get("L0_L0_cos_theta_star_vs_delta_eta_delta_phi_hist");



  TH1F *L0bar_L0bar_cosThetaProdPlane = (TH1F*)inFile->Get("L0bar_L0bar_cosThetaProdPlane");
  TH1F *L0bar_L0bar_cosThetaProdPlane_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar_L0bar_cosThetaProdPlane_eta_hist[nEtaBins][nEtaBins];

  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_phi_for_corr_hist = (TH2F*)inFile->Get("L0bar_L0bar_cos_theta_star_vs_delta_phi_for_corr_hist");
  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_eta_hist = (TH2F*)inFile->Get("L0bar_L0bar_cos_theta_star_vs_delta_eta_hist");
  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_hist = (TH2F*)inFile->Get("L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_hist");

  //pi kinematics
  //TH2F *L0_L0bar_pi_pT1_vs_pi_pT2_hist = (TH2F*)inFile->Get("L0_L0bar_pi_pT1_vs_pi_pT2_hist");

  //---------------------------------------------------------------------------------------------------------
  //mixed event
  //True MC
  TH1F *L0_L0bar_cosThetaProdPlane_ME = (TH1F*)inFile->Get("L0_L0bar_cosThetaProdPlane_ME");
  TH1F *L0_L0bar_cosThetaProdPlane_ME_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0bar_cosThetaProdPlane_ME_eta_hist[nEtaBins][nEtaBins];


  TH1F *L0_L0_cosThetaProdPlane_ME = (TH1F*)inFile->Get("L0_L0_cosThetaProdPlane_ME");
  TH1F *L0_L0_cosThetaProdPlane_ME_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0_cosThetaProdPlane_ME_eta_hist[nEtaBins][nEtaBins];


  TH1F *L0bar_L0bar_cosThetaProdPlane_ME = (TH1F*)inFile->Get("L0bar_L0bar_cosThetaProdPlane_ME");
  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_eta_hist[nEtaBins][nEtaBins];

  //pi kinematics
  //TH2F *L0_L0bar_pi_pT1_vs_pi_pT2_ME_hist = (TH2F*)inFile->Get("L0_L0bar_pi_pT1_vs_pi_pT2_ME_hist");


  //_____________________________________________________________________________________________________________________________

  //after cuts

  //True MC
  TH1F *L0_L0bar_cosThetaProdPlane_cuts = (TH1F*)inFile->Get("L0_L0bar_cosThetaProdPlane_cuts");
  TH1F *L0_L0bar_cosThetaProdPlane_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0bar_cosThetaProdPlane_eta_cuts_hist[nEtaBins][nEtaBins];

  TH2F *L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_cuts_hist = (TH2F*)inFile->Get("L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_cuts_hist");
  TH2F *L0_L0bar_cos_theta_star_vs_delta_eta_cuts_hist = (TH2F*)inFile->Get("L0_L0bar_cos_theta_star_vs_delta_eta_cuts_hist");
  TH2F *L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_cuts_hist = (TH2F*)inFile->Get("L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_cuts_hist");

  //QA
  TH2F *L0_L0bar_delta_phi_cuts_hist = (TH2F*)inFile->Get("L0_L0bar_delta_phi_cuts_hist");
  TH2F *L0_L0bar_delta_eta_cuts_hist = (TH2F*)inFile->Get("L0_L0bar_delta_eta_cuts_hist");


  TH1F *L0_L0_cosThetaProdPlane_cuts = (TH1F*)inFile->Get("L0_L0_cosThetaProdPlane_cuts");
  TH1F *L0_L0_cosThetaProdPlane_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0_cosThetaProdPlane_eta_cuts_hist[nEtaBins][nEtaBins];

  TH2F *L0_L0_cos_theta_star_vs_delta_phi_for_corr_cuts_hist = (TH2F*)inFile->Get("L0_L0_cos_theta_star_vs_delta_phi_for_corr_cuts_hist");

  TH2F *L0_L0_cos_theta_star_vs_delta_eta_cuts_hist = (TH2F*)inFile->Get("L0_L0_cos_theta_star_vs_delta_eta_cuts_hist");

  TH2F *L0_L0_cos_theta_star_vs_delta_eta_delta_phi_cuts_hist = (TH2F*)inFile->Get("L0_L0_cos_theta_star_vs_delta_eta_delta_phi_cuts_hist");


  TH1F *L0bar_L0bar_cosThetaProdPlane_cuts = (TH1F*)inFile->Get("L0bar_L0bar_cosThetaProdPlane_cuts");
  TH1F *L0bar_L0bar_cosThetaProdPlane_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar_L0bar_cosThetaProdPlane_eta_cuts_hist[nEtaBins][nEtaBins];

  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_phi_for_corr_cuts_hist = (TH2F*)inFile->Get("L0bar_L0bar_cos_theta_star_vs_delta_phi_for_corr_cuts_hist");

  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_eta_cuts_hist = (TH2F*)inFile->Get("L0bar_L0bar_cos_theta_star_vs_delta_eta_cuts_hist");

  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_cuts_hist = (TH2F*)inFile->Get("L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_cuts_hist");

  //pi kinematics
  //TH2F *L0_L0bar_pi_pT1_vs_pi_pT2_cuts_hist = (TH2F*)inFile->Get("L0_L0bar_pi_pT1_vs_pi_pT2_cuts_hist");

  //---------------------------------------------------------------------------------------------------------


  //mixed event
  //true MC
  //weight from distributions after cuts
  TH1F *L0_L0bar_cosThetaProdPlane_ME_cuts = (TH1F*)inFile->Get("L0_L0bar_cosThetaProdPlane_ME_cuts");
  TH1F *L0_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[nEtaBins][nEtaBins];

  TH2F *L0_L0bar_cos_theta_star_vs_delta_phi_ME_cuts_hist = (TH2F*)inFile->Get("L0_L0bar_cos_theta_star_vs_delta_phi_ME_cuts_hist");
  TH2F *L0_L0bar_cos_theta_star_vs_delta_phi_ME_cuts_2_hist = (TH2F*)inFile->Get("L0_L0bar_cos_theta_star_vs_delta_phi_ME_cuts_2_hist");

  TH2F *L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_ME_cuts_hist = (TH2F*)inFile->Get("L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_ME_cuts_hist");

  TH2F *L0_L0bar_cos_theta_star_vs_delta_eta_ME_cuts_hist = (TH2F*)inFile->Get("L0_L0bar_cos_theta_star_vs_delta_eta_ME_cuts_hist");
  TH2F *L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_ME_cuts_hist = (TH2F*)inFile->Get("L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_ME_cuts_hist");

  //QA
  TH2F *L0_L0bar_delta_phi_ME_cuts_hist = (TH2F*)inFile->Get("L0_L0bar_delta_phi_ME_cuts_hist");
  TH2F *L0_L0bar_delta_eta_ME_cuts_hist = (TH2F*)inFile->Get("L0_L0bar_delta_eta_ME_cuts_hist");


  TH1F *L0_L0_cosThetaProdPlane_ME_cuts = (TH1F*)inFile->Get("L0_L0_cosThetaProdPlane_ME_cuts");
  TH1F *L0_L0_cosThetaProdPlane_ME_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0_cosThetaProdPlane_ME_eta_cuts_hist[nEtaBins][nEtaBins];

  TH2F *L0_L0_cos_theta_star_vs_delta_phi_ME_cuts_hist = (TH2F*)inFile->Get("L0_L0_cos_theta_star_vs_delta_phi_ME_cuts_hist");
  TH2F *L0_L0_cos_theta_star_vs_delta_phi_ME_cuts_2_hist = (TH2F*)inFile->Get("L0_L0_cos_theta_star_vs_delta_phi_ME_cuts_2_hist");

  TH2F *L0_L0_cos_theta_star_vs_delta_phi_for_corr_ME_cuts_hist = (TH2F*)inFile->Get("L0_L0_cos_theta_star_vs_delta_phi_for_corr_ME_cuts_hist");

  TH2F *L0_L0_cos_theta_star_vs_delta_eta_ME_cuts_hist = (TH2F*)inFile->Get("L0_L0_cos_theta_star_vs_delta_eta_ME_cuts_hist");
  TH2F *L0_L0_cos_theta_star_vs_delta_eta_delta_phi_ME_cuts_hist = (TH2F*)inFile->Get("L0_L0_cos_theta_star_vs_delta_eta_delta_phi_ME_cuts_hist");


  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_cuts = (TH1F*)inFile->Get("L0bar_L0bar_cosThetaProdPlane_ME_cuts");
  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[nEtaBins][nEtaBins];

  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_phi_ME_cuts_hist = (TH2F*)inFile->Get("L0bar_L0bar_cos_theta_star_vs_delta_phi_ME_cuts_hist");
  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_phi_ME_cuts_2_hist = (TH2F*)inFile->Get("L0bar_L0bar_cos_theta_star_vs_delta_phi_ME_cuts_2_hist");

  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_phi_for_corr_ME_cuts_hist = (TH2F*)inFile->Get("L0bar_L0bar_cos_theta_star_vs_delta_phi_for_corr_ME_cuts_hist");

  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_eta_ME_cuts_hist = (TH2F*)inFile->Get("L0bar_L0bar_cos_theta_star_vs_delta_eta_ME_cuts_hist");
  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_ME_cuts_hist = (TH2F*)inFile->Get("L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_ME_cuts_hist");

  //pi kinematics
  //TH2F *L0_L0bar_pi_pT1_vs_pi_pT2_ME_cuts_hist = (TH2F*)inFile->Get("L0_L0bar_pi_pT1_vs_pi_pT2_ME_cuts_hist");

  //_________________________________________________________________________________________________________________________________________________________________________________


  //True MC before cuts
  TCanvas *L0_L0bar_cosThetaProdPlane_can = new TCanvas("L0_L0bar_cosThetaProdPlane_can", "L0_L0bar_cosThetaProdPlane_can", 1200, 1000);
  TCanvas *L0_L0_cosThetaProdPlane_can = new TCanvas("L0_L0_cosThetaProdPlane_can", "L0_L0_cosThetaProdPlane_can", 1200, 1000);
  TCanvas *L0bar_L0bar_cosThetaProdPlane_can = new TCanvas("L0bar_L0bar_cosThetaProdPlane_can", "L0bar_L0bar_cosThetaProdPlane_can", 1200, 1000);

  //True MC after cuts
  TCanvas *L0_L0bar_cosThetaProdPlane_cuts_can = new TCanvas("L0_L0bar_cosThetaProdPlane_cuts_can", "L0_L0bar_cosThetaProdPlane_cuts_can", 1200, 1000);
  TCanvas *L0_L0_cosThetaProdPlane_cuts_can = new TCanvas("L0_L0_cosThetaProdPlane_cuts_can", "L0_L0_cosThetaProdPlane_cuts_can", 1200, 1000);
  TCanvas *L0bar_L0bar_cosThetaProdPlane_cuts_can = new TCanvas("L0bar_L0bar_cosThetaProdPlane_cuts_can", "L0bar_L0bar_cosThetaProdPlane_cuts_can", 1200, 1000);


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

  L0_L0bar_cosThetaProdPlane_ME->Sumw2();
  //L0_L0bar_cosThetaProdPlane_ME->Add(L0_L0bar_cosThetaProdPlane_LS_ME);
  L0_L0bar_cosThetaProdPlane_ME->Scale(nL0L0bar/L0_L0bar_cosThetaProdPlane_ME->Integral());
  //L0_L0bar_cosThetaProdPlane_ME->Draw("p e");

  L0_L0bar_cosThetaProdPlane->Sumw2();
  L0_L0bar_cosThetaProdPlane->SetMarkerStyle(20);
  L0_L0bar_cosThetaProdPlane->SetMarkerColor(kRed);
  L0_L0bar_cosThetaProdPlane->SetLineColor(kRed);
  L0_L0bar_cosThetaProdPlane->GetXaxis()->SetTitle("cos(#theta*)");
  L0_L0bar_cosThetaProdPlane->GetXaxis()->CenterTitle();
  L0_L0bar_cosThetaProdPlane->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0_L0bar_cosThetaProdPlane->GetYaxis()->CenterTitle();
  L0_L0bar_cosThetaProdPlane->Divide(L0_L0bar_cosThetaProdPlane_ME); //correct using ME
  L0_L0bar_cosThetaProdPlane->Scale(nL0L0bar/L0_L0bar_cosThetaProdPlane->Integral()); //scale back
  L0_L0bar_cosThetaProdPlane->Scale(1./L0_L0bar_cosThetaProdPlane->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
  L0_L0bar_cosThetaProdPlane->Fit(fit_L0_L0bar_before_cuts, "i 0 r");
  L0_L0bar_cosThetaProdPlane->SetMinimum(0);
  //L0_L0bar_cosThetaProdPlane->GetYaxis()->SetRangeUser(-20e9, 0);
  L0_L0bar_cosThetaProdPlane->Draw("p e");

  //L0_L0bar_cosThetaProdPlane_ME->Draw("p e same");

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


  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  //after cuts
  TF1 *fit_L0_L0bar_after_cuts = new TF1("fit_L0_L0bar_after_cuts", "[0]*(1+[1]*x)", -1, 1);
  fit_L0_L0bar_after_cuts->SetParameters(1000, 0.10);

  TF1 *fit_L0_L0bar_after_cuts_ME = new TF1("fit_L0_L0bar_after_cuts_ME", "[0]*(1+[1]*x)", -1, 1);
  fit_L0_L0bar_after_cuts_ME->SetParameters(1000, 0.10);

  L0_L0bar_cosThetaProdPlane_cuts_can->cd();

  float nL0L0bar_cuts = L0_L0bar_cosThetaProdPlane_cuts->Integral();

  L0_L0bar_cosThetaProdPlane_ME_cuts->Sumw2();
  L0_L0bar_cosThetaProdPlane_ME_cuts->Scale(nL0L0bar_cuts/L0_L0bar_cosThetaProdPlane_ME_cuts->Integral());
  L0_L0bar_cosThetaProdPlane_ME_cuts->Scale(1./L0_L0bar_cosThetaProdPlane_ME_cuts->GetXaxis()->GetBinWidth(1));
  L0_L0bar_cosThetaProdPlane_ME_cuts->Fit(fit_L0_L0bar_after_cuts_ME, "i 0 r");
  //L0_L0bar_cosThetaProdPlane_ME->Draw("p e");

  L0_L0bar_cosThetaProdPlane_cuts->Sumw2();
  L0_L0bar_cosThetaProdPlane_cuts->SetMarkerStyle(20);
  L0_L0bar_cosThetaProdPlane_cuts->SetMarkerColor(kRed);
  L0_L0bar_cosThetaProdPlane_cuts->SetLineColor(kRed);
  L0_L0bar_cosThetaProdPlane_cuts->GetXaxis()->SetTitle("cos(#theta*)");
  L0_L0bar_cosThetaProdPlane_cuts->GetXaxis()->CenterTitle();
  L0_L0bar_cosThetaProdPlane_cuts->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0_L0bar_cosThetaProdPlane_cuts->GetYaxis()->CenterTitle();
  //L0_L0bar_cosThetaProdPlane_cuts->Divide(L0_L0bar_cosThetaProdPlane_ME_cuts); //correct using ME
  //L0_L0bar_cosThetaProdPlane_cuts->Scale(nL0L0bar_cuts/L0_L0bar_cosThetaProdPlane_cuts->Integral()); //scale back
  L0_L0bar_cosThetaProdPlane_cuts->Scale(1./L0_L0bar_cosThetaProdPlane_cuts->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
  L0_L0bar_cosThetaProdPlane_cuts->Fit(fit_L0_L0bar_after_cuts, "i 0 r");
  L0_L0bar_cosThetaProdPlane_cuts->SetMinimum(0);
  //L0_L0bar_cosThetaProdPlane_cuts->GetYaxis()->SetRangeUser(-20e9, 0);
  L0_L0bar_cosThetaProdPlane_cuts->Draw("p e");

  L0_L0bar_cosThetaProdPlane_ME_cuts->Draw("p e same");

  fit_L0_L0bar_after_cuts->SetLineColor(1);
  fit_L0_L0bar_after_cuts->Draw("same");

  fit_L0_L0bar_after_cuts_ME->SetLineColor(kBlue);
  fit_L0_L0bar_after_cuts_ME->Draw("same");


  float L0_L0bar_slope_cuts = fit_L0_L0bar_after_cuts->GetParameter(1)/(L0_alpha*L0bar_alpha);
  float L0_L0bar_slope_cuts_err = fit_L0_L0bar_after_cuts->GetParError(1)/(L0_alpha*L0bar_alpha);


  //statistical error correction for systematic error
  float SysErrSlope_L0_L0bar_corr = sqrt( fabs( fit_L0_L0bar_after_cuts->GetParError(1)*fit_L0_L0bar_after_cuts->GetParError(1) - fit_L0_L0bar_after_cuts_ME->GetParError(1)*fit_L0_L0bar_after_cuts_ME->GetParError(1) ) );

  float SysErrSlope_L0_L0bar_work = ( fabs( fit_L0_L0bar_after_cuts->GetParameter(1) - fit_L0_L0bar_after_cuts_ME->GetParameter(1) ) - SysErrSlope_L0_L0bar_corr )/fabs(fit_L0_L0bar_after_cuts->GetParameter(1));

  float L0_L0bar_slope_cuts_sys_err = 0;

  if( SysErrSlope_L0_L0bar_work > 0 ) L0_L0bar_slope_cuts_sys_err = SysErrSlope_L0_L0bar_work; //store sys. err. only if it is larger than statistical fluctuations

  float L0_L0bar_slope_cuts_sys_err_no_corr = fabs( fit_L0_L0bar_after_cuts->GetParameter(1) - fit_L0_L0bar_after_cuts_ME->GetParameter(1) )/fit_L0_L0bar_after_cuts->GetParameter(1);

  if(corr_err == 1) SysErrSlope_hist->SetBinContent(1, L0_L0bar_slope_cuts_sys_err);
  if(corr_err == 0) SysErrSlope_hist->SetBinContent(1, L0_L0bar_slope_cuts_sys_err_no_corr);

  TPaveText *L0_L0bar_text_MC_cuts = new TPaveText(0.5, 0.2, 0.8, 0.5, "NDC");
  L0_L0bar_text_MC_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  L0_L0bar_text_MC_cuts->AddText("Minimum bias");
  L0_L0bar_text_MC_cuts->AddText("#Lambda-#bar{#Lambda}");
  L0_L0bar_text_MC_cuts->AddText("True MC");
  L0_L0bar_text_MC_cuts->AddText("Analysis cuts");
  L0_L0bar_text_MC_cuts->AddText(Form("P = %.3f #pm %.3f", L0_L0bar_slope_cuts, L0_L0bar_slope_cuts_err ));
  L0_L0bar_text_MC_cuts->AddText(Form("P_{ME} = %.3f", fit_L0_L0bar_after_cuts_ME->GetParameter(1)/(L0_alpha*L0bar_alpha) ));
  if(corr_err == 1) L0_L0bar_text_MC_cuts->AddText(Form("Slope sys. err: %.1f %%", L0_L0bar_slope_cuts_sys_err*100 ));
  if(corr_err == 0) L0_L0bar_text_MC_cuts->AddText(Form("Slope sys. err: %.1f %%", L0_L0bar_slope_cuts_sys_err_no_corr*100 ));
  L0_L0bar_text_MC_cuts->SetFillColorAlpha(0, 0.01);
  L0_L0bar_text_MC_cuts->Draw("same");

  L0_L0bar_cosThetaProdPlane_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/L0_L0bar_cosThetaProdPlane_cuts.png");


  //___________________________________________________________________________________________________________________________________________________________________________________________________________________________

  //L-L
  //before cuts
  TF1 *fit_L0_L0_before_cuts = new TF1("fit_L0_L0_before_cuts", "[0]*(1+[1]*x)", -1, 1);
  fit_L0_L0_before_cuts->SetParameters(1000, 0.10);

  L0_L0_cosThetaProdPlane_can->cd();

  float nL0L0 = L0_L0_cosThetaProdPlane->Integral();

  L0_L0_cosThetaProdPlane_ME->Sumw2();
  //L0_L0_cosThetaProdPlane_ME->Add(L0_L0_cosThetaProdPlane_LS_ME);
  L0_L0_cosThetaProdPlane_ME->Scale(nL0L0/L0_L0_cosThetaProdPlane_ME->Integral());
  //L0_L0_cosThetaProdPlane_ME->Draw("p e");

  L0_L0_cosThetaProdPlane->Sumw2();
  L0_L0_cosThetaProdPlane->SetMarkerStyle(20);
  L0_L0_cosThetaProdPlane->SetMarkerColor(kRed);
  L0_L0_cosThetaProdPlane->SetLineColor(kRed);
  L0_L0_cosThetaProdPlane->GetXaxis()->SetTitle("cos(#theta*)");
  L0_L0_cosThetaProdPlane->GetXaxis()->CenterTitle();
  L0_L0_cosThetaProdPlane->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0_L0_cosThetaProdPlane->GetYaxis()->CenterTitle();
  L0_L0_cosThetaProdPlane->Divide(L0_L0_cosThetaProdPlane_ME); //correct using ME
  L0_L0_cosThetaProdPlane->Scale(nL0L0/L0_L0_cosThetaProdPlane->Integral()); //scale back
  //L0_L0_cosThetaProdPlane->Scale(1./L0_L0_cosThetaProdPlane->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
  L0_L0_cosThetaProdPlane->Fit(fit_L0_L0_before_cuts, "i 0 r");
  L0_L0_cosThetaProdPlane->SetMinimum(0);
  //L0_L0_cosThetaProdPlane->GetYaxis()->SetRangeUser(-20e9, 0);
  L0_L0_cosThetaProdPlane->Draw("p e");

  //L0_L0_cosThetaProdPlane_ME->Draw("p e same");

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

  L0_L0_cosThetaProdPlane_ME_cuts->Sumw2();
  L0_L0_cosThetaProdPlane_ME_cuts->Scale(nL0L0_cuts/L0_L0_cosThetaProdPlane_ME_cuts->Integral());
  L0_L0_cosThetaProdPlane_ME_cuts->Fit(fit_L0_L0_after_cuts_ME, "i 0 r");
  //L0_L0_cosThetaProdPlane_ME->Draw("p e");

  L0_L0_cosThetaProdPlane_cuts->Sumw2();
  L0_L0_cosThetaProdPlane_cuts->SetMarkerStyle(20);
  L0_L0_cosThetaProdPlane_cuts->SetMarkerColor(kRed);
  L0_L0_cosThetaProdPlane_cuts->SetLineColor(kRed);
  L0_L0_cosThetaProdPlane_cuts->GetXaxis()->SetTitle("cos(#theta*)");
  L0_L0_cosThetaProdPlane_cuts->GetXaxis()->CenterTitle();
  L0_L0_cosThetaProdPlane_cuts->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0_L0_cosThetaProdPlane_cuts->GetYaxis()->CenterTitle();
  //L0_L0_cosThetaProdPlane_cuts->Divide(L0_L0_cosThetaProdPlane_ME_cuts); //correct using ME
  //L0_L0_cosThetaProdPlane_cuts->Scale(nL0L0_cuts/L0_L0_cosThetaProdPlane_cuts->Integral()); //scale back
  //L0_L0_cosThetaProdPlane_cuts->Scale(1./L0_L0_cosThetaProdPlane_cuts->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
  L0_L0_cosThetaProdPlane_cuts->Fit(fit_L0_L0_after_cuts, "i 0 r");
  L0_L0_cosThetaProdPlane_cuts->SetMinimum(0);
  //L0_L0_cosThetaProdPlane_cuts->GetYaxis()->SetRangeUser(-20e9, 0);
  L0_L0_cosThetaProdPlane_cuts->Draw("p e");

  L0_L0_cosThetaProdPlane_ME_cuts->Draw("p e same");

  fit_L0_L0_after_cuts->SetLineColor(1);
  fit_L0_L0_after_cuts->Draw("same");

  float L0_L0_slope_cuts = fit_L0_L0_after_cuts->GetParameter(1)/(L0_alpha*L0_alpha);
  float L0_L0_slope_cuts_err = fit_L0_L0_after_cuts->GetParError(1)/(L0_alpha*L0_alpha);

  //statistical error correction for systematic error
  float SysErrSlope_L0_L0_corr = sqrt( fabs( fit_L0_L0_after_cuts->GetParError(1)*fit_L0_L0_after_cuts->GetParError(1) - fit_L0_L0_after_cuts_ME->GetParError(1)*fit_L0_L0_after_cuts_ME->GetParError(1) ) );

  float SysErrSlope_L0_L0_work = ( fabs( fit_L0_L0_after_cuts->GetParameter(1) - fit_L0_L0_after_cuts_ME->GetParameter(1) ) - SysErrSlope_L0_L0_corr )/fabs(fit_L0_L0_after_cuts->GetParameter(1));

  float L0_L0_slope_cuts_sys_err = 0;

  if( SysErrSlope_L0_L0_work > 0 ) L0_L0_slope_cuts_sys_err = SysErrSlope_L0_L0_work; //store sys. err. only if it is larger than statistical fluctuations

  float L0_L0_slope_cuts_sys_err_no_corr = fabs( fit_L0_L0_after_cuts->GetParameter(1) - fit_L0_L0_after_cuts_ME->GetParameter(1) )/fabs(fit_L0_L0_after_cuts->GetParameter(1));

  if(corr_err == 1) SysErrSlope_hist->SetBinContent(2, L0_L0_slope_cuts_sys_err);
  if(corr_err == 0) SysErrSlope_hist->SetBinContent(2, L0_L0_slope_cuts_sys_err_no_corr);

  TPaveText *L0_L0_text_MC_cuts = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
  L0_L0_text_MC_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  L0_L0_text_MC_cuts->AddText("Minimum bias");
  L0_L0_text_MC_cuts->AddText("#Lambda-#Lambda");
  L0_L0_text_MC_cuts->AddText("True MC");
  L0_L0_text_MC_cuts->AddText("Analysis cuts");
  L0_L0_text_MC_cuts->AddText(Form("P = %.3f #pm %.3f", L0_L0_slope_cuts, L0_L0_slope_cuts_err));
  if(corr_err == 1) L0_L0_text_MC_cuts->AddText(Form("Sys. err: %.1f %%", L0_L0_slope_cuts_sys_err*100 ));
  if(corr_err == 0) L0_L0_text_MC_cuts->AddText(Form("Sys. err: %.1f %%", L0_L0_slope_cuts_sys_err_no_corr*100 ));
  L0_L0_text_MC_cuts->SetFillColorAlpha(0, 0.01);
  L0_L0_text_MC_cuts->Draw("same");

  L0_L0_cosThetaProdPlane_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_L/L0_L0_cosThetaProdPlane_cuts.png");


  //_________________________________________________________________________

  //Lbar-Lbar
  //before cuts
  TF1 *fit_L0bar_L0bar_before_cuts = new TF1("fit_L0bar_L0bar_before_cuts", "[0]*(1+[1]*x)", -1, 1);
  fit_L0bar_L0bar_before_cuts->SetParameters(1000, 0.10);

  L0bar_L0bar_cosThetaProdPlane_can->cd();

  float nL0barL0bar = L0bar_L0bar_cosThetaProdPlane->Integral();

  L0bar_L0bar_cosThetaProdPlane_ME->Sumw2();
  //L0bar_L0bar_cosThetaProdPlane_ME->Add(L0bar_L0bar_cosThetaProdPlane_LS_ME);
  L0bar_L0bar_cosThetaProdPlane_ME->Scale(nL0barL0bar/L0bar_L0bar_cosThetaProdPlane_ME->Integral());
  //L0bar_L0bar_cosThetaProdPlane_ME->Draw("p e");

  L0bar_L0bar_cosThetaProdPlane->Sumw2();
  L0bar_L0bar_cosThetaProdPlane->SetMarkerStyle(20);
  L0bar_L0bar_cosThetaProdPlane->SetMarkerColor(kRed);
  L0bar_L0bar_cosThetaProdPlane->SetLineColor(kRed);
  L0bar_L0bar_cosThetaProdPlane->GetXaxis()->SetTitle("cos(#theta*)");
  L0bar_L0bar_cosThetaProdPlane->GetXaxis()->CenterTitle();
  L0bar_L0bar_cosThetaProdPlane->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0bar_L0bar_cosThetaProdPlane->GetYaxis()->CenterTitle();
  L0bar_L0bar_cosThetaProdPlane->Divide(L0bar_L0bar_cosThetaProdPlane_ME); //correct using ME
  L0bar_L0bar_cosThetaProdPlane->Scale(nL0barL0bar/L0bar_L0bar_cosThetaProdPlane->Integral()); //scale back
  //L0bar_L0bar_cosThetaProdPlane->Scale(1./L0bar_L0bar_cosThetaProdPlane->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
  L0bar_L0bar_cosThetaProdPlane->Fit(fit_L0bar_L0bar_before_cuts, "i 0 r");
  L0bar_L0bar_cosThetaProdPlane->SetMinimum(0);
  //L0bar_L0bar_cosThetaProdPlane->GetYaxis()->SetRangeUser(-20e9, 0);
  L0bar_L0bar_cosThetaProdPlane->Draw("p e");

  //L0bar_L0bar_cosThetaProdPlane_ME->Draw("p e same");

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

  L0bar_L0bar_cosThetaProdPlane_ME_cuts->Sumw2();
  L0bar_L0bar_cosThetaProdPlane_ME_cuts->Scale(nL0barL0bar_cuts/L0bar_L0bar_cosThetaProdPlane_ME_cuts->Integral());
  L0bar_L0bar_cosThetaProdPlane_ME_cuts->Fit(fit_L0bar_L0bar_after_cuts_ME, "i 0 r");
  //L0bar_L0bar_cosThetaProdPlane_ME->Draw("p e");

  L0bar_L0bar_cosThetaProdPlane_cuts->Sumw2();
  L0bar_L0bar_cosThetaProdPlane_cuts->SetMarkerStyle(20);
  L0bar_L0bar_cosThetaProdPlane_cuts->SetMarkerColor(kRed);
  L0bar_L0bar_cosThetaProdPlane_cuts->SetLineColor(kRed);
  L0bar_L0bar_cosThetaProdPlane_cuts->GetXaxis()->SetTitle("cos(#theta*)");
  L0bar_L0bar_cosThetaProdPlane_cuts->GetXaxis()->CenterTitle();
  L0bar_L0bar_cosThetaProdPlane_cuts->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0bar_L0bar_cosThetaProdPlane_cuts->GetYaxis()->CenterTitle();
  //L0bar_L0bar_cosThetaProdPlane_cuts->Divide(L0bar_L0bar_cosThetaProdPlane_ME_cuts); //correct using ME
  //L0bar_L0bar_cosThetaProdPlane_cuts->Scale(nL0barL0bar_cuts/L0bar_L0bar_cosThetaProdPlane_cuts->Integral()); //scale back
  //L0bar_L0bar_cosThetaProdPlane_cuts->Scale(1./L0bar_L0bar_cosThetaProdPlane_cuts->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
  L0bar_L0bar_cosThetaProdPlane_cuts->Fit(fit_L0bar_L0bar_after_cuts, "i 0 r");
  L0bar_L0bar_cosThetaProdPlane_cuts->SetMinimum(0);
  //L0bar_L0bar_cosThetaProdPlane_cuts->GetYaxis()->SetRangeUser(-20e9, 0);
  L0bar_L0bar_cosThetaProdPlane_cuts->Draw("p e");

  L0bar_L0bar_cosThetaProdPlane_ME_cuts->Draw("p e same");

  fit_L0bar_L0bar_after_cuts->SetLineColor(1);
  fit_L0bar_L0bar_after_cuts->Draw("same");

  float L0bar_L0bar_slope_cuts = fit_L0bar_L0bar_after_cuts->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
  float L0bar_L0bar_slope_cuts_err = fit_L0bar_L0bar_after_cuts->GetParError(1)/(L0bar_alpha*L0bar_alpha);

  //statistical error correction for systematic error
  float SysErrSlope_L0bar_L0bar_corr = sqrt( fabs( fit_L0bar_L0bar_after_cuts->GetParError(1)*fit_L0bar_L0bar_after_cuts->GetParError(1) - fit_L0bar_L0bar_after_cuts_ME->GetParError(1)*fit_L0bar_L0bar_after_cuts_ME->GetParError(1) ) );

  float SysErrSlope_L0bar_L0bar_work = ( fabs( fit_L0bar_L0bar_after_cuts->GetParameter(1) - fit_L0bar_L0bar_after_cuts_ME->GetParameter(1) ) - SysErrSlope_L0bar_L0bar_corr )/fabs(fit_L0bar_L0bar_after_cuts->GetParameter(1));

  float L0bar_L0bar_slope_cuts_sys_err = 0;

  if( SysErrSlope_L0bar_L0bar_work > 0 ) L0bar_L0bar_slope_cuts_sys_err = SysErrSlope_L0bar_L0bar_work; //store sys. err. only if it is larger than statistical fluctuations

  float L0bar_L0bar_slope_cuts_sys_err_no_corr = fabs( fit_L0bar_L0bar_after_cuts->GetParameter(1) - fit_L0bar_L0bar_after_cuts_ME->GetParameter(1) )/fabs(fit_L0bar_L0bar_after_cuts->GetParameter(1));

  if(corr_err == 1) SysErrSlope_hist->SetBinContent(3, L0bar_L0bar_slope_cuts_sys_err);
  if(corr_err == 0) SysErrSlope_hist->SetBinContent(3, L0bar_L0bar_slope_cuts_sys_err_no_corr);

  TPaveText *L0bar_L0bar_text_MC_cuts = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
  L0bar_L0bar_text_MC_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  L0bar_L0bar_text_MC_cuts->AddText("Minimum bias");
  L0bar_L0bar_text_MC_cuts->AddText("#bar{#Lambda}-#bar{#Lambda}");
  L0bar_L0bar_text_MC_cuts->AddText("True MC");
  L0bar_L0bar_text_MC_cuts->AddText("Analysis cuts");
  L0bar_L0bar_text_MC_cuts->AddText(Form("P = %.3f #pm %.3f", L0bar_L0bar_slope_cuts, L0bar_L0bar_slope_cuts_err));
  if(corr_err == 1) L0bar_L0bar_text_MC_cuts->AddText(Form("Slope sys. err: %.1f %%", L0bar_L0bar_slope_cuts_sys_err*100 ));
  if(corr_err == 0) L0bar_L0bar_text_MC_cuts->AddText(Form("Slope sys. err: %.1f %%", L0bar_L0bar_slope_cuts_sys_err_no_corr*100 ));
  L0bar_L0bar_text_MC_cuts->SetFillColorAlpha(0, 0.01);
  L0bar_L0bar_text_MC_cuts->Draw("same");

  L0bar_L0bar_cosThetaProdPlane_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/Lbar_Lbar/L0bar_L0bar_cosThetaProdPlane_cuts.png");

  //---------------------------

  sysErrFile->cd();

  SysErrSlope_hist->Write();

  //___________________________________________________________________________________________________________________________________________________________________________________________________________

  TGraphErrors *PolarizationGraph_delta_eta_PYTHIA[2];
  TGraphErrors *PolarizationGraph_delta_phi_PYTHIA[2];
  TGraphErrors *PolarizationGraph_delta_eta_delta_phi_PYTHIA[2];

  //-----------------------

  //for Delta phi projections

  int bin_min = 1;
  int bin_middle = 20;
  int bin_max = 60;


  for( unsigned int delta_eta_bin = 1; delta_eta_bin < 3; delta_eta_bin++)
  {
    PolarizationGraph_delta_eta_PYTHIA[delta_eta_bin-1] = new TGraphErrors(3);
    PolarizationGraph_delta_phi_PYTHIA[delta_eta_bin-1] = new TGraphErrors(3);
    PolarizationGraph_delta_eta_delta_phi_PYTHIA[delta_eta_bin-1] = new TGraphErrors(3);

    //---------------------------------------------------------------------------------------

    //for Delta phi projections

    int bin_min_loop = 0;
    int bin_max_loop = 0;

    float bin_edge = 0;

    if(delta_eta_bin == 1)
    {
      bin_min_loop = bin_min;
      bin_max_loop = bin_middle;
    }
    if(delta_eta_bin == 2)
    {
      bin_min_loop = bin_middle+1;
      bin_max_loop = bin_max;
    }

    //_______________________________________________________________________________________

    out_file_polarization->cd();

    //before cuts
    //L-Lbar
    TCanvas *L0_L0bar_cosThetaProdPlane_delta_eta_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_delta_eta_can_%i", delta_eta_bin), Form("L0_L0bar_cosThetaProdPlane_delta_eta_can_%i", delta_eta_bin), 1200, 1000);
    L0_L0bar_cosThetaProdPlane_delta_eta_can->cd();

    TF1 *fit_L0_L0bar_after_delta_eta = new TF1("fit_L0_L0bar_after_delta_eta", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0bar_after_delta_eta->SetParameters(1000, 0.10);

    TH1D *L0_L0bar_cosThetaProdPlane_delta_eta = L0_L0bar_cos_theta_star_vs_delta_eta_hist->ProjectionX(Form("L_Lbar_proj_eta_SE_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
    float nL0L0bar_delta_eta = L0_L0bar_cosThetaProdPlane_delta_eta->Integral();
    L0_L0bar_cosThetaProdPlane_delta_eta->Sumw2();
    L0_L0bar_cosThetaProdPlane_delta_eta->SetMarkerStyle(20);
    L0_L0bar_cosThetaProdPlane_delta_eta->SetMarkerColor(kRed);
    L0_L0bar_cosThetaProdPlane_delta_eta->SetLineColor(kRed);
    L0_L0bar_cosThetaProdPlane_delta_eta->GetXaxis()->SetTitle("cos(#theta*)");
    L0_L0bar_cosThetaProdPlane_delta_eta->GetXaxis()->CenterTitle();
    L0_L0bar_cosThetaProdPlane_delta_eta->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    L0_L0bar_cosThetaProdPlane_delta_eta->GetYaxis()->CenterTitle();
    //L0_L0bar_cosThetaProdPlane_delta_eta->Divide(L0_L0bar_cosThetaProdPlane_ME_delta_eta); //correct using ME
    //L0_L0bar_cosThetaProdPlane_delta_eta->Scale(nL0L0bar_delta_eta/L0_L0bar_cosThetaProdPlane_delta_eta->Integral()); //scale back
    L0_L0bar_cosThetaProdPlane_delta_eta->Scale(1./L0_L0bar_cosThetaProdPlane_delta_eta->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
    L0_L0bar_cosThetaProdPlane_delta_eta->Fit(fit_L0_L0bar_after_delta_eta, "i 0 r");
    L0_L0bar_cosThetaProdPlane_delta_eta->SetMinimum(0);
    //L0_L0bar_cosThetaProdPlane_delta_eta->GetYaxis()->SetRangeUser(-20e9, 0);
    L0_L0bar_cosThetaProdPlane_delta_eta->Draw("p e");

    fit_L0_L0bar_after_delta_eta->SetLineColor(1);
    fit_L0_L0bar_after_delta_eta->Draw("same");


    float L0_L0bar_slope_delta_eta = fit_L0_L0bar_after_delta_eta->GetParameter(1)/(L0_alpha*L0bar_alpha);
    float L0_L0bar_slope_delta_eta_err = fit_L0_L0bar_after_delta_eta->GetParError(1)/(L0_alpha*L0bar_alpha);

    TPaveText *L0_L0bar_text_MC_delta_eta = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    L0_L0bar_text_MC_delta_eta->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0_L0bar_text_MC_delta_eta->AddText("Minimum bias");
    L0_L0bar_text_MC_delta_eta->AddText("#Lambda-#bar{#Lambda}");
    L0_L0bar_text_MC_delta_eta->AddText("True MC");
    L0_L0bar_text_MC_delta_eta->AddText("No cuts");
    L0_L0bar_text_MC_delta_eta->AddText(Form("P = %.3f #pm %.3f", L0_L0bar_slope_delta_eta, L0_L0bar_slope_delta_eta_err ));
    L0_L0bar_text_MC_delta_eta->SetFillColorAlpha(0, 0.01);
    L0_L0bar_text_MC_delta_eta->Draw("same");

    PolarizationGraph_delta_eta_PYTHIA[delta_eta_bin-1]->SetPoint(1, L0_L0bar_slope_delta_eta, 1.15);
    PolarizationGraph_delta_eta_PYTHIA[delta_eta_bin-1]->SetPointError(1, fabs(L0_L0bar_slope_delta_eta_err), 0);

    L0_L0bar_cosThetaProdPlane_delta_eta_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/L0_L0bar_cosThetaProdPlane_delta_eta_%i.png", delta_eta_bin));

    //-----------------------------------------------------------

    //L-L
    TCanvas *L0_L0_cosThetaProdPlane_delta_eta_can = new TCanvas(Form("L0_L0_cosThetaProdPlane_delta_eta_can_%i", delta_eta_bin), Form("L0_L0_cosThetaProdPlane_delta_eta_can_%i", delta_eta_bin), 1200, 1000);
    L0_L0_cosThetaProdPlane_delta_eta_can->cd();

    TF1 *fit_L0_L0_after_delta_eta = new TF1("fit_L0_L0_after_delta_eta", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0_after_delta_eta->SetParameters(1000, 0.10);

    TH1D *L0_L0_cosThetaProdPlane_delta_eta = L0_L0_cos_theta_star_vs_delta_eta_hist->ProjectionX(Form("LL_proj_eta_SE_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
    float nL0L0_delta_eta = L0_L0_cosThetaProdPlane_delta_eta->Integral();
    L0_L0_cosThetaProdPlane_delta_eta->Sumw2();
    L0_L0_cosThetaProdPlane_delta_eta->SetMarkerStyle(20);
    L0_L0_cosThetaProdPlane_delta_eta->SetMarkerColor(kRed);
    L0_L0_cosThetaProdPlane_delta_eta->SetLineColor(kRed);
    L0_L0_cosThetaProdPlane_delta_eta->GetXaxis()->SetTitle("cos(#theta*)");
    L0_L0_cosThetaProdPlane_delta_eta->GetXaxis()->CenterTitle();
    L0_L0_cosThetaProdPlane_delta_eta->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    L0_L0_cosThetaProdPlane_delta_eta->GetYaxis()->CenterTitle();
    //L0_L0_cosThetaProdPlane_delta_eta->Divide(L0_L0_cosThetaProdPlane_ME_delta_eta); //correct using ME
    //L0_L0_cosThetaProdPlane_delta_eta->Scale(nL0L0_delta_eta/L0_L0_cosThetaProdPlane_delta_eta->Integral()); //scale back
    L0_L0_cosThetaProdPlane_delta_eta->Scale(1./L0_L0_cosThetaProdPlane_delta_eta->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
    L0_L0_cosThetaProdPlane_delta_eta->Fit(fit_L0_L0_after_delta_eta, "i 0 r");
    L0_L0_cosThetaProdPlane_delta_eta->SetMinimum(0);
    //L0_L0_cosThetaProdPlane_delta_eta->GetYaxis()->SetRangeUser(-20e9, 0);
    L0_L0_cosThetaProdPlane_delta_eta->Draw("p e");

    fit_L0_L0_after_delta_eta->SetLineColor(1);
    fit_L0_L0_after_delta_eta->Draw("same");

    float L0_L0_slope_delta_eta = fit_L0_L0_after_delta_eta->GetParameter(1)/(L0_alpha*L0_alpha);
    float L0_L0_slope_delta_eta_err = fit_L0_L0_after_delta_eta->GetParError(1)/(L0_alpha*L0_alpha);


    TPaveText *L0_L0_text_MC_delta_eta = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    L0_L0_text_MC_delta_eta->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0_L0_text_MC_delta_eta->AddText("Minimum bias");
    L0_L0_text_MC_delta_eta->AddText("#Lambda-#Lambda");
    L0_L0_text_MC_delta_eta->AddText("True MC");
    L0_L0_text_MC_delta_eta->AddText("No cuts");
    L0_L0_text_MC_delta_eta->AddText(Form("P = %.3f #pm %.3f", L0_L0_slope_delta_eta, L0_L0_slope_delta_eta_err ));
    L0_L0_text_MC_delta_eta->SetFillColorAlpha(0, 0.01);
    L0_L0_text_MC_delta_eta->Draw("same");

    PolarizationGraph_delta_eta_PYTHIA[delta_eta_bin-1]->SetPoint(2, L0_L0_slope_delta_eta, 2.15);
    PolarizationGraph_delta_eta_PYTHIA[delta_eta_bin-1]->SetPointError(2, fabs(L0_L0_slope_delta_eta_err), 0);

    L0_L0_cosThetaProdPlane_delta_eta_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_L/L0_L0_cosThetaProdPlane_delta_eta_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------

    //Lbar-Lbar
    TCanvas *L0bar_L0bar_cosThetaProdPlane_delta_eta_can = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_delta_eta_can_%i", delta_eta_bin), Form("L0bar_L0bar_cosThetaProdPlane_delta_eta_can_%i", delta_eta_bin), 1200, 1000);
    L0bar_L0bar_cosThetaProdPlane_delta_eta_can->cd();

    TF1 *fit_L0bar_L0bar_after_delta_eta = new TF1("fit_L0bar_L0bar_after_delta_eta", "[0]*(1+[1]*x)", -1, 1);
    fit_L0bar_L0bar_after_delta_eta->SetParameters(1000, 0.10);

    TH1D *L0bar_L0bar_cosThetaProdPlane_delta_eta = L0bar_L0bar_cos_theta_star_vs_delta_eta_hist->ProjectionX(Form("LbaLbar_proj_eta_SE_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
    float nL0barL0bar_delta_eta = L0bar_L0bar_cosThetaProdPlane_delta_eta->Integral();
    L0bar_L0bar_cosThetaProdPlane_delta_eta->Sumw2();
    L0bar_L0bar_cosThetaProdPlane_delta_eta->SetMarkerStyle(20);
    L0bar_L0bar_cosThetaProdPlane_delta_eta->SetMarkerColor(kRed);
    L0bar_L0bar_cosThetaProdPlane_delta_eta->SetLineColor(kRed);
    L0bar_L0bar_cosThetaProdPlane_delta_eta->GetXaxis()->SetTitle("cos(#theta*)");
    L0bar_L0bar_cosThetaProdPlane_delta_eta->GetXaxis()->CenterTitle();
    L0bar_L0bar_cosThetaProdPlane_delta_eta->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    L0bar_L0bar_cosThetaProdPlane_delta_eta->GetYaxis()->CenterTitle();
    //L0bar_L0bar_cosThetaProdPlane_delta_eta->Divide(L0bar_L0bar_cosThetaProdPlane_ME_delta_eta); //correct using ME
    //L0bar_L0bar_cosThetaProdPlane_delta_eta->Scale(nL0barL0bar_delta_eta/L0bar_L0bar_cosThetaProdPlane_delta_eta->Integral()); //scale back
    L0bar_L0bar_cosThetaProdPlane_delta_eta->Scale(1./L0bar_L0bar_cosThetaProdPlane_delta_eta->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
    L0bar_L0bar_cosThetaProdPlane_delta_eta->Fit(fit_L0bar_L0bar_after_delta_eta, "i 0 r");
    L0bar_L0bar_cosThetaProdPlane_delta_eta->SetMinimum(0);
    //L0bar_L0bar_cosThetaProdPlane_delta_eta->GetYaxis()->SetRangeUser(-20e9, 0);
    L0bar_L0bar_cosThetaProdPlane_delta_eta->Draw("p e");

    fit_L0bar_L0bar_after_delta_eta->SetLineColor(1);
    fit_L0bar_L0bar_after_delta_eta->Draw("same");

    float L0bar_L0bar_slope_delta_eta = fit_L0bar_L0bar_after_delta_eta->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
    float L0bar_L0bar_slope_delta_eta_err = fit_L0bar_L0bar_after_delta_eta->GetParError(1)/(L0bar_alpha*L0bar_alpha);


    TPaveText *L0bar_L0bar_text_MC_delta_eta = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    L0bar_L0bar_text_MC_delta_eta->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0bar_L0bar_text_MC_delta_eta->AddText("Minimum bias");
    L0bar_L0bar_text_MC_delta_eta->AddText("#bar{#Lambda}-#bar{#Lambda}");
    L0bar_L0bar_text_MC_delta_eta->AddText("True MC");
    L0bar_L0bar_text_MC_delta_eta->AddText("No cuts");
    L0bar_L0bar_text_MC_delta_eta->AddText(Form("P = %.3f #pm %.3f", L0bar_L0bar_slope_delta_eta, L0bar_L0bar_slope_delta_eta_err ));
    L0bar_L0bar_text_MC_delta_eta->SetFillColorAlpha(0, 0.01);
    L0bar_L0bar_text_MC_delta_eta->Draw("same");

    PolarizationGraph_delta_eta_PYTHIA[delta_eta_bin-1]->SetPoint(3, L0bar_L0bar_slope_delta_eta, 3.15);
    PolarizationGraph_delta_eta_PYTHIA[delta_eta_bin-1]->SetPointError(3, fabs(L0bar_L0bar_slope_delta_eta_err), 0);

    L0bar_L0bar_cosThetaProdPlane_delta_eta_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/Lbar_Lbar/L0bar_L0bar_cosThetaProdPlane_delta_eta_%i.png", delta_eta_bin));

    PolarizationGraph_delta_eta_PYTHIA[delta_eta_bin-1]->Write(Form("PolarizationGraph_delta_eta_PYTHIA_%i", delta_eta_bin-1));

    //-----------------------------------------------------------------------------------------------------------------------------

    //before cuts
    //L-Lbar
    TCanvas *L0_L0bar_cosThetaProdPlane_delta_phi_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_delta_phi_can_%i", delta_eta_bin), Form("L0_L0bar_cosThetaProdPlane_delta_phi_can_%i", delta_eta_bin), 1200, 1000);
    L0_L0bar_cosThetaProdPlane_delta_phi_can->cd();

    TF1 *fit_L0_L0bar_after_delta_phi = new TF1("fit_L0_L0bar_after_delta_phi", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0bar_after_delta_phi->SetParameters(1000, 0.10);

    TH1D *L0_L0bar_cosThetaProdPlane_delta_phi = L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_hist->ProjectionX(Form("L_Lbar_proj_phi_SE_%i", delta_eta_bin), bin_min_loop, bin_max_loop);
    float nL0L0bar_delta_phi = L0_L0bar_cosThetaProdPlane_delta_phi->Integral();
    L0_L0bar_cosThetaProdPlane_delta_phi->Sumw2();
    L0_L0bar_cosThetaProdPlane_delta_phi->SetMarkerStyle(20);
    L0_L0bar_cosThetaProdPlane_delta_phi->SetMarkerColor(kRed);
    L0_L0bar_cosThetaProdPlane_delta_phi->SetLineColor(kRed);
    L0_L0bar_cosThetaProdPlane_delta_phi->GetXaxis()->SetTitle("cos(#theta*)");
    L0_L0bar_cosThetaProdPlane_delta_phi->GetXaxis()->CenterTitle();
    L0_L0bar_cosThetaProdPlane_delta_phi->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    L0_L0bar_cosThetaProdPlane_delta_phi->GetYaxis()->CenterTitle();
    //L0_L0bar_cosThetaProdPlane_delta_phi->Divide(L0_L0bar_cosThetaProdPlane_ME_delta_phi); //correct using ME
    //L0_L0bar_cosThetaProdPlane_delta_phi->Scale(nL0L0bar_delta_phi/L0_L0bar_cosThetaProdPlane_delta_phi->Integral()); //scale back
    L0_L0bar_cosThetaProdPlane_delta_phi->Scale(1./L0_L0bar_cosThetaProdPlane_delta_phi->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
    L0_L0bar_cosThetaProdPlane_delta_phi->Fit(fit_L0_L0bar_after_delta_phi, "i 0 r");
    L0_L0bar_cosThetaProdPlane_delta_phi->SetMinimum(0);
    //L0_L0bar_cosThetaProdPlane_delta_phi->GetYaxis()->SetRangeUser(-20e9, 0);
    L0_L0bar_cosThetaProdPlane_delta_phi->Draw("p e");

    fit_L0_L0bar_after_delta_phi->SetLineColor(1);
    fit_L0_L0bar_after_delta_phi->Draw("same");


    float L0_L0bar_slope_delta_phi = fit_L0_L0bar_after_delta_phi->GetParameter(1)/(L0_alpha*L0bar_alpha);
    float L0_L0bar_slope_delta_phi_err = fit_L0_L0bar_after_delta_phi->GetParError(1)/(L0_alpha*L0bar_alpha);

    TPaveText *L0_L0bar_text_MC_delta_phi = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    L0_L0bar_text_MC_delta_phi->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0_L0bar_text_MC_delta_phi->AddText("Minimum bias");
    L0_L0bar_text_MC_delta_phi->AddText("#Lambda-#bar{#Lambda}");
    L0_L0bar_text_MC_delta_phi->AddText("True MC");
    L0_L0bar_text_MC_delta_phi->AddText("No cuts");
    L0_L0bar_text_MC_delta_phi->AddText(Form("P = %.3f #pm %.3f", L0_L0bar_slope_delta_phi, L0_L0bar_slope_delta_phi_err ));
    L0_L0bar_text_MC_delta_phi->SetFillColorAlpha(0, 0.01);
    L0_L0bar_text_MC_delta_phi->Draw("same");

    PolarizationGraph_delta_phi_PYTHIA[delta_eta_bin-1]->SetPoint(1, L0_L0bar_slope_delta_phi, 1.15);
    PolarizationGraph_delta_phi_PYTHIA[delta_eta_bin-1]->SetPointError(1, fabs(L0_L0bar_slope_delta_phi_err), 0);

    L0_L0bar_cosThetaProdPlane_delta_phi_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/L0_L0bar_cosThetaProdPlane_delta_phi_%i.png", delta_eta_bin));

    //-----------------------------------------------------------

    //L-L
    TCanvas *L0_L0_cosThetaProdPlane_delta_phi_can = new TCanvas(Form("L0_L0_cosThetaProdPlane_delta_phi_can_%i", delta_eta_bin), Form("L0_L0_cosThetaProdPlane_delta_phi_can_%i", delta_eta_bin), 1200, 1000);
    L0_L0_cosThetaProdPlane_delta_phi_can->cd();

    TF1 *fit_L0_L0_after_delta_phi = new TF1("fit_L0_L0_after_delta_phi", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0_after_delta_phi->SetParameters(1000, 0.10);

    TH1D *L0_L0_cosThetaProdPlane_delta_phi = L0_L0_cos_theta_star_vs_delta_phi_for_corr_hist->ProjectionX(Form("LL_proj_phi_SE_%i", delta_eta_bin), bin_min_loop, bin_max_loop);
    float nL0L0_delta_phi = L0_L0_cosThetaProdPlane_delta_phi->Integral();
    L0_L0_cosThetaProdPlane_delta_phi->Sumw2();
    L0_L0_cosThetaProdPlane_delta_phi->SetMarkerStyle(20);
    L0_L0_cosThetaProdPlane_delta_phi->SetMarkerColor(kRed);
    L0_L0_cosThetaProdPlane_delta_phi->SetLineColor(kRed);
    L0_L0_cosThetaProdPlane_delta_phi->GetXaxis()->SetTitle("cos(#theta*)");
    L0_L0_cosThetaProdPlane_delta_phi->GetXaxis()->CenterTitle();
    L0_L0_cosThetaProdPlane_delta_phi->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    L0_L0_cosThetaProdPlane_delta_phi->GetYaxis()->CenterTitle();
    //L0_L0_cosThetaProdPlane_delta_phi->Divide(L0_L0_cosThetaProdPlane_ME_delta_phi); //correct using ME
    //L0_L0_cosThetaProdPlane_delta_phi->Scale(nL0L0_delta_phi/L0_L0_cosThetaProdPlane_delta_phi->Integral()); //scale back
    L0_L0_cosThetaProdPlane_delta_phi->Scale(1./L0_L0_cosThetaProdPlane_delta_phi->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
    L0_L0_cosThetaProdPlane_delta_phi->Fit(fit_L0_L0_after_delta_phi, "i 0 r");
    L0_L0_cosThetaProdPlane_delta_phi->SetMinimum(0);
    //L0_L0_cosThetaProdPlane_delta_phi->GetYaxis()->SetRangeUser(-20e9, 0);
    L0_L0_cosThetaProdPlane_delta_phi->Draw("p e");

    fit_L0_L0_after_delta_phi->SetLineColor(1);
    fit_L0_L0_after_delta_phi->Draw("same");

    float L0_L0_slope_delta_phi = fit_L0_L0_after_delta_phi->GetParameter(1)/(L0_alpha*L0_alpha);
    float L0_L0_slope_delta_phi_err = fit_L0_L0_after_delta_phi->GetParError(1)/(L0_alpha*L0_alpha);


    TPaveText *L0_L0_text_MC_delta_phi = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    L0_L0_text_MC_delta_phi->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0_L0_text_MC_delta_phi->AddText("Minimum bias");
    L0_L0_text_MC_delta_phi->AddText("#Lambda-#Lambda");
    L0_L0_text_MC_delta_phi->AddText("True MC");
    L0_L0_text_MC_delta_phi->AddText("No cuts");
    L0_L0_text_MC_delta_phi->AddText(Form("P = %.3f #pm %.3f", L0_L0_slope_delta_phi, L0_L0_slope_delta_phi_err ));
    L0_L0_text_MC_delta_phi->SetFillColorAlpha(0, 0.01);
    L0_L0_text_MC_delta_phi->Draw("same");

    PolarizationGraph_delta_phi_PYTHIA[delta_eta_bin-1]->SetPoint(2, L0_L0_slope_delta_phi, 2.15);
    PolarizationGraph_delta_phi_PYTHIA[delta_eta_bin-1]->SetPointError(2, fabs(L0_L0_slope_delta_phi_err), 0);

    L0_L0_cosThetaProdPlane_delta_phi_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_L/L0_L0_cosThetaProdPlane_delta_phi_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------

    //Lbar-Lbar
    TCanvas *L0bar_L0bar_cosThetaProdPlane_delta_phi_can = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_delta_phi_can_%i", delta_eta_bin), Form("L0bar_L0bar_cosThetaProdPlane_delta_phi_can_%i", delta_eta_bin), 1200, 1000);
    L0bar_L0bar_cosThetaProdPlane_delta_phi_can->cd();

    TF1 *fit_L0bar_L0bar_after_delta_phi = new TF1("fit_L0bar_L0bar_after_delta_phi", "[0]*(1+[1]*x)", -1, 1);
    fit_L0bar_L0bar_after_delta_phi->SetParameters(1000, 0.10);

    TH1D *L0bar_L0bar_cosThetaProdPlane_delta_phi = L0bar_L0bar_cos_theta_star_vs_delta_phi_for_corr_hist->ProjectionX(Form("LbaLbar_proj_phi_SE_%i", delta_eta_bin), bin_min_loop, bin_max_loop);
    float nL0barL0bar_delta_phi = L0bar_L0bar_cosThetaProdPlane_delta_phi->Integral();
    L0bar_L0bar_cosThetaProdPlane_delta_phi->Sumw2();
    L0bar_L0bar_cosThetaProdPlane_delta_phi->SetMarkerStyle(20);
    L0bar_L0bar_cosThetaProdPlane_delta_phi->SetMarkerColor(kRed);
    L0bar_L0bar_cosThetaProdPlane_delta_phi->SetLineColor(kRed);
    L0bar_L0bar_cosThetaProdPlane_delta_phi->GetXaxis()->SetTitle("cos(#theta*)");
    L0bar_L0bar_cosThetaProdPlane_delta_phi->GetXaxis()->CenterTitle();
    L0bar_L0bar_cosThetaProdPlane_delta_phi->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    L0bar_L0bar_cosThetaProdPlane_delta_phi->GetYaxis()->CenterTitle();
    //L0bar_L0bar_cosThetaProdPlane_delta_phi->Divide(L0bar_L0bar_cosThetaProdPlane_ME_delta_phi); //correct using ME
    //L0bar_L0bar_cosThetaProdPlane_delta_phi->Scale(nL0barL0bar_delta_phi/L0bar_L0bar_cosThetaProdPlane_delta_phi->Integral()); //scale back
    L0bar_L0bar_cosThetaProdPlane_delta_phi->Scale(1./L0bar_L0bar_cosThetaProdPlane_delta_phi->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
    L0bar_L0bar_cosThetaProdPlane_delta_phi->Fit(fit_L0bar_L0bar_after_delta_phi, "i 0 r");
    L0bar_L0bar_cosThetaProdPlane_delta_phi->SetMinimum(0);
    //L0bar_L0bar_cosThetaProdPlane_delta_phi->GetYaxis()->SetRangeUser(-20e9, 0);
    L0bar_L0bar_cosThetaProdPlane_delta_phi->Draw("p e");

    fit_L0bar_L0bar_after_delta_phi->SetLineColor(1);
    fit_L0bar_L0bar_after_delta_phi->Draw("same");

    float L0bar_L0bar_slope_delta_phi = fit_L0bar_L0bar_after_delta_phi->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
    float L0bar_L0bar_slope_delta_phi_err = fit_L0bar_L0bar_after_delta_phi->GetParError(1)/(L0bar_alpha*L0bar_alpha);


    TPaveText *L0bar_L0bar_text_MC_delta_phi = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    L0bar_L0bar_text_MC_delta_phi->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0bar_L0bar_text_MC_delta_phi->AddText("Minimum bias");
    L0bar_L0bar_text_MC_delta_phi->AddText("#bar{#Lambda}-#bar{#Lambda}");
    L0bar_L0bar_text_MC_delta_phi->AddText("True MC");
    L0bar_L0bar_text_MC_delta_phi->AddText("No cuts");
    L0bar_L0bar_text_MC_delta_phi->AddText(Form("P = %.3f #pm %.3f", L0bar_L0bar_slope_delta_phi, L0bar_L0bar_slope_delta_phi_err ));
    L0bar_L0bar_text_MC_delta_phi->SetFillColorAlpha(0, 0.01);
    L0bar_L0bar_text_MC_delta_phi->Draw("same");

    PolarizationGraph_delta_phi_PYTHIA[delta_eta_bin-1]->SetPoint(3, L0bar_L0bar_slope_delta_phi, 3.15);
    PolarizationGraph_delta_phi_PYTHIA[delta_eta_bin-1]->SetPointError(3, fabs(L0bar_L0bar_slope_delta_phi_err), 0);

    L0bar_L0bar_cosThetaProdPlane_delta_phi_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/Lbar_Lbar/L0bar_L0bar_cosThetaProdPlane_delta_phi_%i.png", delta_eta_bin));

    PolarizationGraph_delta_phi_PYTHIA[delta_eta_bin-1]->Write(Form("PolarizationGraph_delta_phi_PYTHIA_%i", delta_eta_bin-1));

    //-----------------------------------------------------------------------------------------------------------------------------

    //before cuts
    //L-Lbar
    TCanvas *L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_can_%i", delta_eta_bin), Form("L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_can_%i", delta_eta_bin), 1200, 1000);
    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_can->cd();

    TF1 *fit_L0_L0bar_after_delta_eta_delta_phi = new TF1("fit_L0_L0bar_after_delta_eta_delta_phi", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0bar_after_delta_eta_delta_phi->SetParameters(1000, 0.10);

    TH1D *L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi = L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_hist->ProjectionX(Form("LLbar_proj_eta_phi_SE_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
    float nL0L0bar_delta_eta_delta_phi = L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi->Integral();
    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi->Sumw2();
    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi->SetMarkerStyle(20);
    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi->SetMarkerSize(1.5);
    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi->SetMarkerColor(kRed);
    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi->SetLineColor(kRed);
    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi->GetXaxis()->SetTitle("cos(#theta*)");
    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi->GetXaxis()->CenterTitle();
    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi->GetYaxis()->CenterTitle();
    //L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi->Divide(L0_L0bar_cosThetaProdPlane_ME_delta_eta_delta_phi); //correct using ME
    //L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi->Scale(nL0L0bar_delta_eta_delta_phi/L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi->Integral()); //scale back
    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi->Scale(1./L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi->Fit(fit_L0_L0bar_after_delta_eta_delta_phi, "i 0 r");
    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi->SetMinimum(0);
    //L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi->GetYaxis()->SetRangeUser(-20e9, 0);
    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi->Draw("p e");

    fit_L0_L0bar_after_delta_eta_delta_phi->SetLineColor(1);
    fit_L0_L0bar_after_delta_eta_delta_phi->Draw("same");


    float L0_L0bar_slope_delta_eta_delta_phi = fit_L0_L0bar_after_delta_eta_delta_phi->GetParameter(1)/(L0_alpha*L0bar_alpha);
    float L0_L0bar_slope_delta_eta_delta_phi_err = fit_L0_L0bar_after_delta_eta_delta_phi->GetParError(1)/(L0_alpha*L0bar_alpha);

    TPaveText *L0_L0bar_text_MC_delta_eta_delta_phi = new TPaveText(0.5, 0.2, 0.8, 0.5, "NDC");
    L0_L0bar_text_MC_delta_eta_delta_phi->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0_L0bar_text_MC_delta_eta_delta_phi->AddText("Minimum bias");
    L0_L0bar_text_MC_delta_eta_delta_phi->AddText("#Lambda-#bar{#Lambda}");
    L0_L0bar_text_MC_delta_eta_delta_phi->AddText("True MC");
    L0_L0bar_text_MC_delta_eta_delta_phi->AddText("No cuts");
    L0_L0bar_text_MC_delta_eta_delta_phi->AddText(Form("P = %.3f #pm %.3f", L0_L0bar_slope_delta_eta_delta_phi, L0_L0bar_slope_delta_eta_delta_phi_err ));
    if( (delta_eta_bin-1) == 0 )L0_L0bar_text_MC_delta_eta_delta_phi->AddText("|#Delta#it{y}| < 0.5, |#Delta#phi| < #pi/3");
    if( (delta_eta_bin-1) == 1 )L0_L0bar_text_MC_delta_eta_delta_phi->AddText("0.5 < |#Delta#it{y}| < 2.0 or |#Delta#phi| > #pi/3");
    L0_L0bar_text_MC_delta_eta_delta_phi->SetFillColorAlpha(0, 0.01);
    L0_L0bar_text_MC_delta_eta_delta_phi->Draw("same");

    PolarizationGraph_delta_eta_delta_phi_PYTHIA[delta_eta_bin-1]->SetPoint(1, L0_L0bar_slope_delta_eta_delta_phi, 1.15);
    PolarizationGraph_delta_eta_delta_phi_PYTHIA[delta_eta_bin-1]->SetPointError(1, fabs(L0_L0bar_slope_delta_eta_delta_phi_err), 0);

    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------

    //L-L
    TCanvas *L0_L0_cosThetaProdPlane_delta_eta_delta_phi_can = new TCanvas(Form("L0_L0_cosThetaProdPlane_delta_eta_delta_phi_can_%i", delta_eta_bin), Form("L0_L0_cosThetaProdPlane_delta_eta_delta_phi_can_%i", delta_eta_bin), 1200, 1000);
    L0_L0_cosThetaProdPlane_delta_eta_delta_phi_can->cd();

    TF1 *fit_L0_L0_after_delta_eta_delta_phi = new TF1("fit_L0_L0_after_delta_eta_delta_phi", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0_after_delta_eta_delta_phi->SetParameters(1000, 0.10);

    TH1D *L0_L0_cosThetaProdPlane_delta_eta_delta_phi = L0_L0_cos_theta_star_vs_delta_eta_delta_phi_hist->ProjectionX(Form("LL_proj_eta_phi_SE_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
    float nL0L0_delta_eta_delta_phi = L0_L0_cosThetaProdPlane_delta_eta_delta_phi->Integral();
    L0_L0_cosThetaProdPlane_delta_eta_delta_phi->Sumw2();
    L0_L0_cosThetaProdPlane_delta_eta_delta_phi->SetMarkerStyle(20);
    L0_L0_cosThetaProdPlane_delta_eta_delta_phi->SetMarkerColor(kRed);
    L0_L0_cosThetaProdPlane_delta_eta_delta_phi->SetLineColor(kRed);
    L0_L0_cosThetaProdPlane_delta_eta_delta_phi->GetXaxis()->SetTitle("cos(#theta*)");
    L0_L0_cosThetaProdPlane_delta_eta_delta_phi->GetXaxis()->CenterTitle();
    L0_L0_cosThetaProdPlane_delta_eta_delta_phi->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    L0_L0_cosThetaProdPlane_delta_eta_delta_phi->GetYaxis()->CenterTitle();
    //L0_L0_cosThetaProdPlane_delta_eta_delta_phi->Divide(L0_L0_cosThetaProdPlane_ME_delta_eta_delta_phi); //correct using ME
    //L0_L0_cosThetaProdPlane_delta_eta_delta_phi->Scale(nL0L0_delta_eta_delta_phi/L0_L0_cosThetaProdPlane_delta_eta_delta_phi->Integral()); //scale back
    L0_L0_cosThetaProdPlane_delta_eta_delta_phi->Scale(1./L0_L0_cosThetaProdPlane_delta_eta_delta_phi->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
    L0_L0_cosThetaProdPlane_delta_eta_delta_phi->Fit(fit_L0_L0_after_delta_eta_delta_phi, "i 0 r");
    L0_L0_cosThetaProdPlane_delta_eta_delta_phi->SetMinimum(0);
    //L0_L0_cosThetaProdPlane_delta_eta_delta_phi->GetYaxis()->SetRangeUser(-20e9, 0);
    L0_L0_cosThetaProdPlane_delta_eta_delta_phi->Draw("p e");

    fit_L0_L0_after_delta_eta_delta_phi->SetLineColor(1);
    fit_L0_L0_after_delta_eta_delta_phi->Draw("same");

    float L0_L0_slope_delta_eta_delta_phi = fit_L0_L0_after_delta_eta_delta_phi->GetParameter(1)/(L0_alpha*L0_alpha);
    float L0_L0_slope_delta_eta_delta_phi_err = fit_L0_L0_after_delta_eta_delta_phi->GetParError(1)/(L0_alpha*L0_alpha);


    TPaveText *L0_L0_text_MC_delta_eta_delta_phi = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    L0_L0_text_MC_delta_eta_delta_phi->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0_L0_text_MC_delta_eta_delta_phi->AddText("Minimum bias");
    L0_L0_text_MC_delta_eta_delta_phi->AddText("#Lambda-#Lambda");
    L0_L0_text_MC_delta_eta_delta_phi->AddText("True MC");
    L0_L0_text_MC_delta_eta_delta_phi->AddText("No cuts");
    L0_L0_text_MC_delta_eta_delta_phi->AddText(Form("P = %.3f #pm %.3f", L0_L0_slope_delta_eta_delta_phi, L0_L0_slope_delta_eta_delta_phi_err ));
    L0_L0_text_MC_delta_eta_delta_phi->SetFillColorAlpha(0, 0.01);
    L0_L0_text_MC_delta_eta_delta_phi->Draw("same");

    PolarizationGraph_delta_eta_delta_phi_PYTHIA[delta_eta_bin-1]->SetPoint(2, L0_L0_slope_delta_eta_delta_phi, 2.15);
    PolarizationGraph_delta_eta_delta_phi_PYTHIA[delta_eta_bin-1]->SetPointError(2, fabs(L0_L0_slope_delta_eta_delta_phi_err), 0);

    L0_L0_cosThetaProdPlane_delta_eta_delta_phi_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_L/L0_L0_cosThetaProdPlane_delta_eta_delta_phi_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------

    //Lbar-Lbar
    TCanvas *L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_can = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_can_%i", delta_eta_bin), Form("L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_can_%i", delta_eta_bin), 1200, 1000);
    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_can->cd();

    TF1 *fit_L0bar_L0bar_after_delta_eta_delta_phi = new TF1("fit_L0bar_L0bar_after_delta_eta_delta_phi", "[0]*(1+[1]*x)", -1, 1);
    fit_L0bar_L0bar_after_delta_eta_delta_phi->SetParameters(1000, 0.10);

    TH1D *L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi = L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_hist->ProjectionX(Form("LbarLbar_proj_eta_phi_SE_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
    float nL0barL0bar_delta_eta_delta_phi = L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi->Integral();
    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi->Sumw2();
    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi->SetMarkerStyle(20);
    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi->SetMarkerColor(kRed);
    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi->SetLineColor(kRed);
    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi->GetXaxis()->SetTitle("cos(#theta*)");
    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi->GetXaxis()->CenterTitle();
    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi->GetYaxis()->CenterTitle();
    //L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi->Divide(L0bar_L0bar_cosThetaProdPlane_ME_delta_eta_delta_phi); //correct using ME
    //L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi->Scale(nL0barL0bar_delta_eta_delta_phi/L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi->Integral()); //scale back
    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi->Scale(1./L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi->Fit(fit_L0bar_L0bar_after_delta_eta_delta_phi, "i 0 r");
    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi->SetMinimum(0);
    //L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi->GetYaxis()->SetRangeUser(-20e9, 0);
    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi->Draw("p e");

    fit_L0bar_L0bar_after_delta_eta_delta_phi->SetLineColor(1);
    fit_L0bar_L0bar_after_delta_eta_delta_phi->Draw("same");

    float L0bar_L0bar_slope_delta_eta_delta_phi = fit_L0bar_L0bar_after_delta_eta_delta_phi->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
    float L0bar_L0bar_slope_delta_eta_delta_phi_err = fit_L0bar_L0bar_after_delta_eta_delta_phi->GetParError(1)/(L0bar_alpha*L0bar_alpha);


    TPaveText *L0bar_L0bar_text_MC_delta_eta_delta_phi = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    L0bar_L0bar_text_MC_delta_eta_delta_phi->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0bar_L0bar_text_MC_delta_eta_delta_phi->AddText("Minimum bias");
    L0bar_L0bar_text_MC_delta_eta_delta_phi->AddText("#bar{#Lambda}-#bar{#Lambda}");
    L0bar_L0bar_text_MC_delta_eta_delta_phi->AddText("True MC");
    L0bar_L0bar_text_MC_delta_eta_delta_phi->AddText("No cuts");
    L0bar_L0bar_text_MC_delta_eta_delta_phi->AddText(Form("P = %.3f #pm %.3f", L0bar_L0bar_slope_delta_eta_delta_phi, L0bar_L0bar_slope_delta_eta_delta_phi_err ));
    L0bar_L0bar_text_MC_delta_eta_delta_phi->SetFillColorAlpha(0, 0.01);
    L0bar_L0bar_text_MC_delta_eta_delta_phi->Draw("same");

    PolarizationGraph_delta_eta_delta_phi_PYTHIA[delta_eta_bin-1]->SetPoint(3, L0bar_L0bar_slope_delta_eta_delta_phi, 3.15);
    PolarizationGraph_delta_eta_delta_phi_PYTHIA[delta_eta_bin-1]->SetPointError(3, fabs(L0bar_L0bar_slope_delta_eta_delta_phi_err), 0);

    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/Lbar_Lbar/L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_%i.png", delta_eta_bin));

    PolarizationGraph_delta_eta_delta_phi_PYTHIA[delta_eta_bin-1]->Write(Form("PolarizationGraph_delta_eta_delta_phi_PYTHIA_%i", delta_eta_bin-1));


    //_______________________________________________________________________________________________________________________________________________________________________________________________________________________

    sysErrFile->cd();

    //after cuts
    //L-Lbar
    TCanvas *L0_L0bar_cosThetaProdPlane_delta_eta_cuts_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_delta_eta_cuts_can_%i", delta_eta_bin), Form("L0_L0bar_cosThetaProdPlane_delta_eta_cuts_can_%i", delta_eta_bin), 1200, 1000);
    L0_L0bar_cosThetaProdPlane_delta_eta_cuts_can->cd();

    TF1 *fit_L0_L0bar_after_delta_eta_cuts = new TF1("fit_L0_L0bar_after_delta_eta_cuts", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0bar_after_delta_eta_cuts->SetParameters(1000, 0.10);

    TF1 *fit_L0_L0bar_after_delta_eta_cuts_ME = new TF1("fit_L0_L0bar_after_delta_eta_cuts_ME", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0bar_after_delta_eta_cuts_ME->SetParameters(1000, 0.10);

    TH1D *L0_L0bar_cosThetaProdPlane_delta_eta_cuts = L0_L0bar_cos_theta_star_vs_delta_eta_cuts_hist->ProjectionX(Form("proj_eta_SE_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
    float nL0L0bar_delta_eta_cuts = L0_L0bar_cosThetaProdPlane_delta_eta_cuts->Integral();

    TH1D *L0_L0bar_cosThetaProdPlane_ME_delta_eta_cuts = L0_L0bar_cos_theta_star_vs_delta_eta_ME_cuts_hist->ProjectionX(Form("proj_eta_ME_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
    L0_L0bar_cosThetaProdPlane_ME_delta_eta_cuts->Sumw2();
    L0_L0bar_cosThetaProdPlane_ME_delta_eta_cuts->SetMarkerStyle(24);
    L0_L0bar_cosThetaProdPlane_ME_delta_eta_cuts->SetMarkerSize(1.5);
    L0_L0bar_cosThetaProdPlane_ME_delta_eta_cuts->SetMarkerColor(1);
    L0_L0bar_cosThetaProdPlane_ME_delta_eta_cuts->SetLineColor(1);
    L0_L0bar_cosThetaProdPlane_ME_delta_eta_cuts->Scale(nL0L0bar_delta_eta_cuts/L0_L0bar_cosThetaProdPlane_ME_delta_eta_cuts->Integral());
    L0_L0bar_cosThetaProdPlane_ME_delta_eta_cuts->Scale(1./L0_L0bar_cosThetaProdPlane_ME_delta_eta_cuts->GetXaxis()->GetBinWidth(1));
    L0_L0bar_cosThetaProdPlane_ME_delta_eta_cuts->Fit(fit_L0_L0bar_after_delta_eta_cuts_ME, "i 0 r");
    //L0_L0bar_cosThetaProdPlane_ME->Draw("p e");


    L0_L0bar_cosThetaProdPlane_delta_eta_cuts->Sumw2();
    L0_L0bar_cosThetaProdPlane_delta_eta_cuts->SetMarkerStyle(20);
    L0_L0bar_cosThetaProdPlane_delta_eta_cuts->SetMarkerSize(1.5);
    L0_L0bar_cosThetaProdPlane_delta_eta_cuts->SetMarkerColor(kRed);
    L0_L0bar_cosThetaProdPlane_delta_eta_cuts->SetLineColor(kRed);
    L0_L0bar_cosThetaProdPlane_delta_eta_cuts->GetXaxis()->SetTitle("cos(#theta*)");
    L0_L0bar_cosThetaProdPlane_delta_eta_cuts->GetXaxis()->CenterTitle();
    L0_L0bar_cosThetaProdPlane_delta_eta_cuts->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    L0_L0bar_cosThetaProdPlane_delta_eta_cuts->GetYaxis()->CenterTitle();
    //L0_L0bar_cosThetaProdPlane_delta_eta_cuts->Divide(L0_L0bar_cosThetaProdPlane_ME_delta_eta_cuts); //correct using ME
    //L0_L0bar_cosThetaProdPlane_delta_eta_cuts->Scale(nL0L0bar_delta_eta_cuts/L0_L0bar_cosThetaProdPlane_delta_eta_cuts->Integral()); //scale back
    L0_L0bar_cosThetaProdPlane_delta_eta_cuts->Scale(1./L0_L0bar_cosThetaProdPlane_delta_eta_cuts->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
    L0_L0bar_cosThetaProdPlane_delta_eta_cuts->Fit(fit_L0_L0bar_after_delta_eta_cuts, "i 0 r");
    L0_L0bar_cosThetaProdPlane_delta_eta_cuts->SetMinimum(0);
    //L0_L0bar_cosThetaProdPlane_delta_eta_cuts->GetYaxis()->SetRangeUser(-20e9, 0);
    L0_L0bar_cosThetaProdPlane_delta_eta_cuts->Draw("p e");

    L0_L0bar_cosThetaProdPlane_ME_delta_eta_cuts->Draw("p e same");

    fit_L0_L0bar_after_delta_eta_cuts->SetLineColor(1);
    fit_L0_L0bar_after_delta_eta_cuts->Draw("same");

    fit_L0_L0bar_after_delta_eta_cuts_ME->SetLineColor(1);
    fit_L0_L0bar_after_delta_eta_cuts_ME->SetLineStyle(7);
    fit_L0_L0bar_after_delta_eta_cuts_ME->Draw("same");


    float L0_L0bar_slope_delta_eta_cuts = fit_L0_L0bar_after_delta_eta_cuts->GetParameter(1)/(L0_alpha*L0bar_alpha);
    float L0_L0bar_slope_delta_eta_cuts_err = fit_L0_L0bar_after_delta_eta_cuts->GetParError(1)/(L0_alpha*L0bar_alpha);

    //statistical error correction for systematic error
    float SysErrSlope_L0_L0bar_delta_eta_corr = sqrt( fabs( fit_L0_L0bar_after_delta_eta_cuts->GetParError(1)*fit_L0_L0bar_after_delta_eta_cuts->GetParError(1) - fit_L0_L0bar_after_delta_eta_cuts_ME->GetParError(1)*fit_L0_L0bar_after_delta_eta_cuts_ME->GetParError(1) ) );

    float SysErrSlope_L0_L0bar_delta_eta_work = ( fabs( fit_L0_L0bar_after_delta_eta_cuts->GetParameter(1) - fit_L0_L0bar_after_delta_eta_cuts_ME->GetParameter(1) ) - SysErrSlope_L0_L0bar_delta_eta_corr )/fabs(fit_L0_L0bar_after_delta_eta_cuts->GetParameter(1));

    float L0_L0bar_slope_delta_eta_cuts_sys_err = 0;

    if( SysErrSlope_L0_L0bar_delta_eta_work > 0 ) L0_L0bar_slope_delta_eta_cuts_sys_err = SysErrSlope_L0_L0bar_delta_eta_work; //store sys. err. only if it is larger than statistical fluctuations

    float L0_L0bar_slope_delta_eta_cuts_sys_err_no_corr = fabs( fit_L0_L0bar_after_delta_eta_cuts->GetParameter(1) - fit_L0_L0bar_after_delta_eta_cuts_ME->GetParameter(1) )/fabs(fit_L0_L0bar_after_delta_eta_cuts->GetParameter(1));

    if(corr_err == 1) SysErrSlope_delta_eta_hist[delta_eta_bin-1]->SetBinContent(1, L0_L0bar_slope_delta_eta_cuts_sys_err);
    if(corr_err == 0) SysErrSlope_delta_eta_hist[delta_eta_bin-1]->SetBinContent(1, L0_L0bar_slope_delta_eta_cuts_sys_err_no_corr);

    TPaveText *L0_L0bar_text_MC_delta_eta_cuts = new TPaveText(0.5, 0.2, 0.8, 0.45, "NDC");
    L0_L0bar_text_MC_delta_eta_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0_L0bar_text_MC_delta_eta_cuts->AddText("Minimum bias");
    L0_L0bar_text_MC_delta_eta_cuts->AddText("#Lambda-#bar{#Lambda}");
    L0_L0bar_text_MC_delta_eta_cuts->AddText("True MC");
    L0_L0bar_text_MC_delta_eta_cuts->AddText("Analysis cuts");
    L0_L0bar_text_MC_delta_eta_cuts->AddText(Form("P = %.3f #pm %.3f", L0_L0bar_slope_delta_eta_cuts, L0_L0bar_slope_delta_eta_cuts_err ));
    L0_L0bar_text_MC_delta_eta_cuts->AddText(Form("P_{ME} = %.3f", fit_L0_L0bar_after_delta_eta_cuts_ME->GetParameter(1)/(L0_alpha*L0bar_alpha) ));
    if( (delta_eta_bin-1) == 0 )L0_L0bar_text_MC_delta_eta_cuts->AddText("|#Delta#it{y}| < 0.5");
    if( (delta_eta_bin-1) == 1 )L0_L0bar_text_MC_delta_eta_cuts->AddText("0.5 < |#Delta#it{y}| < 2.0");
    //if(corr_err == 1) L0_L0bar_text_MC_delta_eta_cuts->AddText(Form("Slope sys. err: %.1f %%", L0_L0bar_slope_delta_eta_cuts_sys_err*100 ));
    //if(corr_err == 0) L0_L0bar_text_MC_delta_eta_cuts->AddText(Form("Slope sys. err: %.1f %%", L0_L0bar_slope_delta_eta_cuts_sys_err_no_corr*100 ));
    L0_L0bar_text_MC_delta_eta_cuts->SetFillColorAlpha(0, 0.01);
    L0_L0bar_text_MC_delta_eta_cuts->Draw("same");

    L0_L0bar_cosThetaProdPlane_delta_eta_cuts_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/L0_L0bar_cosThetaProdPlane_delta_eta_cuts_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------

    TCanvas *L0_L0bar_cosThetaProdPlane_delta_eta_cuts_corr_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_delta_eta_cuts_corr_can_%i", delta_eta_bin), Form("L0_L0bar_cosThetaProdPlane_delta_eta_cuts_corr_can_%i", delta_eta_bin), 1200, 1000);
    L0_L0bar_cosThetaProdPlane_delta_eta_cuts_corr_can->cd();

    TF1 *fit_L0_L0bar_after_delta_eta_cuts_corr = new TF1("fit_L0_L0bar_after_delta_eta_cuts_corr", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0bar_after_delta_eta_cuts_corr->SetParameters(1000, 0.10);

    TH1D *L0_L0bar_cosThetaProdPlane_delta_eta_cuts_corr = (TH1D*)L0_L0bar_cosThetaProdPlane_delta_eta_cuts->Clone(Form("L0_L0bar_cosThetaProdPlane_delta_eta_cuts_corr_%i", delta_eta_bin));
    L0_L0bar_cosThetaProdPlane_delta_eta_cuts_corr->Divide(L0_L0bar_cosThetaProdPlane_ME_delta_eta_cuts); //correct
    L0_L0bar_cosThetaProdPlane_delta_eta_cuts_corr->Scale(nL0L0bar_delta_eta_cuts/L0_L0bar_cosThetaProdPlane_delta_eta_cuts_corr->Integral()); //scale back
    L0_L0bar_cosThetaProdPlane_delta_eta_cuts_corr->Scale(1./L0_L0bar_cosThetaProdPlane_delta_eta_cuts_corr->GetXaxis()->GetBinWidth(1)); //divide by bin width
    //L0_L0bar_cosThetaProdPlane_delta_eta_cuts_corr->Scale(nL0L0bar_delta_eta_cuts);
    L0_L0bar_cosThetaProdPlane_delta_eta_cuts_corr->Fit(fit_L0_L0bar_after_delta_eta_cuts_corr, "i 0 r");
    L0_L0bar_cosThetaProdPlane_delta_eta_cuts_corr->SetMinimum(0);
    L0_L0bar_cosThetaProdPlane_delta_eta_cuts_corr->Draw("p e");

    fit_L0_L0bar_after_delta_eta_cuts_corr->SetLineColor(1);
    fit_L0_L0bar_after_delta_eta_cuts_corr->Draw("same");

    float L0_L0bar_slope_delta_eta_cuts_corr = fit_L0_L0bar_after_delta_eta_cuts_corr->GetParameter(1)/(L0_alpha*L0bar_alpha);
    float L0_L0bar_slope_delta_eta_cuts_corr_err = fit_L0_L0bar_after_delta_eta_cuts_corr->GetParError(1)/(L0_alpha*L0bar_alpha);

    ResidualPolarization_delta_eta_hist[delta_eta_bin-1]->SetBinContent(1, L0_L0bar_slope_delta_eta_cuts_corr);
    ResidualPolarization_delta_eta_hist[delta_eta_bin-1]->SetBinError(1, L0_L0bar_slope_delta_eta_cuts_corr_err);


    TPaveText *L0_L0bar_text_MC_delta_eta_cuts_corr = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    L0_L0bar_text_MC_delta_eta_cuts_corr->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0_L0bar_text_MC_delta_eta_cuts_corr->AddText("Minimum bias");
    L0_L0bar_text_MC_delta_eta_cuts_corr->AddText("#Lambda-#bar{#Lambda}");
    L0_L0bar_text_MC_delta_eta_cuts_corr->AddText("True MC");
    L0_L0bar_text_MC_delta_eta_cuts_corr->AddText("Analysis cuts");
    L0_L0bar_text_MC_delta_eta_cuts_corr->AddText(Form("P = %.3f #pm %.3f", L0_L0bar_slope_delta_eta_cuts_corr, L0_L0bar_slope_delta_eta_cuts_corr_err ));
    //if(corr_err == 1) L0_L0bar_text_MC_delta_eta_cuts_corr->AddText(Form("Slope sys. err: %.1f %%", L0_L0bar_slope_delta_eta_cuts_corr_sys_err*100 ));
    //if(corr_err == 0) L0_L0bar_text_MC_delta_eta_cuts_corr->AddText(Form("Slope sys. err: %.1f %%", L0_L0bar_slope_delta_eta_cuts_corr_sys_err_no_corr*100 ));
    L0_L0bar_text_MC_delta_eta_cuts_corr->SetFillColorAlpha(0, 0.01);
    L0_L0bar_text_MC_delta_eta_cuts_corr->Draw("same");


    L0_L0bar_cosThetaProdPlane_delta_eta_cuts_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/L0_L0bar_cosThetaProdPlane_delta_eta_cuts_corr_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------

    //L-L
    TCanvas *L0_L0_cosThetaProdPlane_delta_eta_cuts_can = new TCanvas(Form("L0_L0_cosThetaProdPlane_delta_eta_cuts_can_%i", delta_eta_bin), Form("L0_L0_cosThetaProdPlane_delta_eta_cuts_can_%i", delta_eta_bin), 1200, 1000);
    L0_L0_cosThetaProdPlane_delta_eta_cuts_can->cd();

    TF1 *fit_L0_L0_after_delta_eta_cuts = new TF1("fit_L0_L0_after_delta_eta_cuts", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0_after_delta_eta_cuts->SetParameters(1000, 0.10);

    TF1 *fit_L0_L0_after_delta_eta_cuts_ME = new TF1("fit_L0_L0_after_delta_eta_cuts_ME", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0_after_delta_eta_cuts_ME->SetParameters(1000, 0.10);

    TH1D *L0_L0_cosThetaProdPlane_delta_eta_cuts = L0_L0_cos_theta_star_vs_delta_eta_cuts_hist->ProjectionX(Form("proj_eta_SE_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
    float nL0L0_delta_eta_cuts = L0_L0_cosThetaProdPlane_delta_eta_cuts->Integral();

    TH1D *L0_L0_cosThetaProdPlane_ME_delta_eta_cuts = L0_L0_cos_theta_star_vs_delta_eta_ME_cuts_hist->ProjectionX(Form("proj_eta_ME_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
    L0_L0_cosThetaProdPlane_ME_delta_eta_cuts->Sumw2();
    L0_L0_cosThetaProdPlane_ME_delta_eta_cuts->Scale(nL0L0_delta_eta_cuts/L0_L0_cosThetaProdPlane_ME_delta_eta_cuts->Integral());
    L0_L0_cosThetaProdPlane_ME_delta_eta_cuts->Scale(1./L0_L0_cosThetaProdPlane_ME_delta_eta_cuts->GetXaxis()->GetBinWidth(1));
    L0_L0_cosThetaProdPlane_ME_delta_eta_cuts->Fit(fit_L0_L0_after_delta_eta_cuts_ME, "i 0 r");
    //L0_L0_cosThetaProdPlane_ME->Draw("p e");


    L0_L0_cosThetaProdPlane_delta_eta_cuts->Sumw2();
    L0_L0_cosThetaProdPlane_delta_eta_cuts->SetMarkerStyle(20);
    L0_L0_cosThetaProdPlane_delta_eta_cuts->SetMarkerColor(kRed);
    L0_L0_cosThetaProdPlane_delta_eta_cuts->SetLineColor(kRed);
    L0_L0_cosThetaProdPlane_delta_eta_cuts->GetXaxis()->SetTitle("cos(#theta*)");
    L0_L0_cosThetaProdPlane_delta_eta_cuts->GetXaxis()->CenterTitle();
    L0_L0_cosThetaProdPlane_delta_eta_cuts->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    L0_L0_cosThetaProdPlane_delta_eta_cuts->GetYaxis()->CenterTitle();
    //L0_L0_cosThetaProdPlane_delta_eta_cuts->Divide(L0_L0_cosThetaProdPlane_ME_delta_eta_cuts); //correct using ME
    //L0_L0_cosThetaProdPlane_delta_eta_cuts->Scale(nL0L0_delta_eta_cuts/L0_L0_cosThetaProdPlane_delta_eta_cuts->Integral()); //scale back
    L0_L0_cosThetaProdPlane_delta_eta_cuts->Scale(1./L0_L0_cosThetaProdPlane_delta_eta_cuts->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
    L0_L0_cosThetaProdPlane_delta_eta_cuts->Fit(fit_L0_L0_after_delta_eta_cuts, "i 0 r");
    L0_L0_cosThetaProdPlane_delta_eta_cuts->SetMinimum(0);
    //L0_L0_cosThetaProdPlane_delta_eta_cuts->GetYaxis()->SetRangeUser(-20e9, 0);
    L0_L0_cosThetaProdPlane_delta_eta_cuts->Draw("p e");

    L0_L0_cosThetaProdPlane_ME_delta_eta_cuts->Draw("p e same");

    fit_L0_L0_after_delta_eta_cuts->SetLineColor(1);
    fit_L0_L0_after_delta_eta_cuts->Draw("same");

    fit_L0_L0_after_delta_eta_cuts_ME->SetLineColor(kBlue);
    fit_L0_L0_after_delta_eta_cuts_ME->Draw("same");


    float L0_L0_slope_delta_eta_cuts = fit_L0_L0_after_delta_eta_cuts->GetParameter(1)/(L0_alpha*L0_alpha);
    float L0_L0_slope_delta_eta_cuts_err = fit_L0_L0_after_delta_eta_cuts->GetParError(1)/(L0_alpha*L0_alpha);


    //statistical error correction for systematic error
    float SysErrSlope_L0_L0_delta_eta_corr = sqrt( fabs( fit_L0_L0_after_delta_eta_cuts->GetParError(1)*fit_L0_L0_after_delta_eta_cuts->GetParError(1) - fit_L0_L0_after_delta_eta_cuts_ME->GetParError(1)*fit_L0_L0_after_delta_eta_cuts_ME->GetParError(1) ) );

    float SysErrSlope_L0_L0_delta_eta_work = ( fabs( fit_L0_L0_after_delta_eta_cuts->GetParameter(1) - fit_L0_L0_after_delta_eta_cuts_ME->GetParameter(1) ) - SysErrSlope_L0_L0_delta_eta_corr )/fabs(fit_L0_L0_after_delta_eta_cuts->GetParameter(1));

    float L0_L0_slope_delta_eta_cuts_sys_err = 0;

    if( SysErrSlope_L0_L0_delta_eta_work > 0 ) L0_L0_slope_delta_eta_cuts_sys_err = SysErrSlope_L0_L0_delta_eta_work; //store sys. err. only if it is larger than statistical fluctuations

    float L0_L0_slope_delta_eta_cuts_sys_err_no_corr = fabs( fit_L0_L0_after_delta_eta_cuts->GetParameter(1) - fit_L0_L0_after_delta_eta_cuts_ME->GetParameter(1) )/fabs(fit_L0_L0_after_delta_eta_cuts->GetParameter(1));

    if(corr_err == 1) SysErrSlope_delta_eta_hist[delta_eta_bin-1]->SetBinContent(2, L0_L0_slope_delta_eta_cuts_sys_err);
    if(corr_err == 0) SysErrSlope_delta_eta_hist[delta_eta_bin-1]->SetBinContent(2, L0_L0_slope_delta_eta_cuts_sys_err_no_corr);

    TPaveText *L0_L0_text_MC_delta_eta_cuts = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    L0_L0_text_MC_delta_eta_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0_L0_text_MC_delta_eta_cuts->AddText("Minimum bias");
    L0_L0_text_MC_delta_eta_cuts->AddText("#Lambda-#Lambda");
    L0_L0_text_MC_delta_eta_cuts->AddText("True MC");
    L0_L0_text_MC_delta_eta_cuts->AddText("Analysis cuts");
    L0_L0_text_MC_delta_eta_cuts->AddText(Form("P = %.3f #pm %.3f", L0_L0_slope_delta_eta_cuts, L0_L0_slope_delta_eta_cuts_err ));
    L0_L0_text_MC_delta_eta_cuts->AddText(Form("P_{ME} = %.3f", fit_L0_L0_after_delta_eta_cuts_ME->GetParameter(1)/(L0_alpha*L0_alpha) ));
    if( (delta_eta_bin-1) == 0 )L0_L0_text_MC_delta_eta_cuts->AddText("|#Delta#it{y}| < 0.5");
    if( (delta_eta_bin-1) == 1 )L0_L0_text_MC_delta_eta_cuts->AddText("0.5 < |#Delta#it{y}| < 2.0");
    if(corr_err == 1) L0_L0_text_MC_delta_eta_cuts->AddText(Form("Slope sys. err: %.1f %%", L0_L0_slope_delta_eta_cuts_sys_err*100 ));
    if(corr_err == 0) L0_L0_text_MC_delta_eta_cuts->AddText(Form("Slope sys. err: %.1f %%", L0_L0_slope_delta_eta_cuts_sys_err_no_corr*100 ));
    L0_L0_text_MC_delta_eta_cuts->SetFillColorAlpha(0, 0.01);
    L0_L0_text_MC_delta_eta_cuts->Draw("same");

    L0_L0_cosThetaProdPlane_delta_eta_cuts_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_L/L0_L0_cosThetaProdPlane_delta_eta_cuts_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------

    TCanvas *L0_L0_cosThetaProdPlane_delta_eta_cuts_corr_can = new TCanvas(Form("L0_L0_cosThetaProdPlane_delta_eta_cuts_corr_can_%i", delta_eta_bin), Form("L0_L0_cosThetaProdPlane_delta_eta_cuts_corr_can_%i", delta_eta_bin), 1200, 1000);
    L0_L0_cosThetaProdPlane_delta_eta_cuts_corr_can->cd();

    TF1 *fit_L0_L0_after_delta_eta_cuts_corr = new TF1("fit_L0_L0_after_delta_eta_cuts_corr", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0_after_delta_eta_cuts_corr->SetParameters(1000, 0.10);

    TH1D *L0_L0_cosThetaProdPlane_delta_eta_cuts_corr = (TH1D*)L0_L0_cosThetaProdPlane_delta_eta_cuts->Clone(Form("L0_L0_cosThetaProdPlane_delta_eta_cuts_corr_%i", delta_eta_bin));
    L0_L0_cosThetaProdPlane_delta_eta_cuts_corr->Divide(L0_L0_cosThetaProdPlane_ME_delta_eta_cuts); //correct
    L0_L0_cosThetaProdPlane_delta_eta_cuts_corr->Scale(nL0L0_delta_eta_cuts/L0_L0_cosThetaProdPlane_delta_eta_cuts_corr->Integral()); //scale back
    L0_L0_cosThetaProdPlane_delta_eta_cuts_corr->Scale(1./L0_L0_cosThetaProdPlane_delta_eta_cuts_corr->GetXaxis()->GetBinWidth(1)); //divide by bin width
    //L0_L0_cosThetaProdPlane_delta_eta_cuts_corr->Scale(nL0L0_delta_eta_cuts);
    L0_L0_cosThetaProdPlane_delta_eta_cuts_corr->Fit(fit_L0_L0_after_delta_eta_cuts_corr, "i 0 r");
    L0_L0_cosThetaProdPlane_delta_eta_cuts_corr->SetMinimum(0);
    L0_L0_cosThetaProdPlane_delta_eta_cuts_corr->Draw("p e");

    fit_L0_L0_after_delta_eta_cuts_corr->SetLineColor(1);
    fit_L0_L0_after_delta_eta_cuts_corr->Draw("same");

    float L0_L0_slope_delta_eta_cuts_corr = fit_L0_L0_after_delta_eta_cuts_corr->GetParameter(1)/(L0_alpha*L0_alpha);
    float L0_L0_slope_delta_eta_cuts_corr_err = fit_L0_L0_after_delta_eta_cuts_corr->GetParError(1)/(L0_alpha*L0_alpha);

    ResidualPolarization_delta_eta_hist[delta_eta_bin-1]->SetBinContent(2, L0_L0_slope_delta_eta_cuts_corr);
    ResidualPolarization_delta_eta_hist[delta_eta_bin-1]->SetBinError(2, L0_L0_slope_delta_eta_cuts_corr_err);


    TPaveText *L0_L0_text_MC_delta_eta_cuts_corr = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    L0_L0_text_MC_delta_eta_cuts_corr->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0_L0_text_MC_delta_eta_cuts_corr->AddText("Minimum bias");
    L0_L0_text_MC_delta_eta_cuts_corr->AddText("#Lambda-#Lambda");
    L0_L0_text_MC_delta_eta_cuts_corr->AddText("True MC");
    L0_L0_text_MC_delta_eta_cuts_corr->AddText("Analysis cuts");
    L0_L0_text_MC_delta_eta_cuts_corr->AddText(Form("P = %.3f #pm %.3f", L0_L0_slope_delta_eta_cuts_corr, L0_L0_slope_delta_eta_cuts_corr_err ));
    //if(corr_err == 1) L0_L0_text_MC_delta_eta_cuts_corr->AddText(Form("Slope sys. err: %.1f %%", L0_L0_slope_delta_eta_cuts_corr_sys_err*100 ));
    //if(corr_err == 0) L0_L0_text_MC_delta_eta_cuts_corr->AddText(Form("Slope sys. err: %.1f %%", L0_L0_slope_delta_eta_cuts_corr_sys_err_no_corr*100 ));
    L0_L0_text_MC_delta_eta_cuts_corr->SetFillColorAlpha(0, 0.01);
    L0_L0_text_MC_delta_eta_cuts_corr->Draw("same");


    L0_L0_cosThetaProdPlane_delta_eta_cuts_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_L/L0_L0_cosThetaProdPlane_delta_eta_cuts_corr_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------

    //Lbar-Lbar
    TCanvas *L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts_can = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts_can_%i", delta_eta_bin), Form("L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts_can_%i", delta_eta_bin), 1200, 1000);
    L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts_can->cd();

    TF1 *fit_L0bar_L0bar_after_delta_eta_cuts = new TF1("fit_L0bar_L0bar_after_delta_eta_cuts", "[0]*(1+[1]*x)", -1, 1);
    fit_L0bar_L0bar_after_delta_eta_cuts->SetParameters(1000, 0.10);

    TF1 *fit_L0bar_L0bar_after_delta_eta_cuts_ME = new TF1("fit_L0bar_L0bar_after_delta_eta_cuts_ME", "[0]*(1+[1]*x)", -1, 1);
    fit_L0bar_L0bar_after_delta_eta_cuts_ME->SetParameters(1000, 0.10);

    TH1D *L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts = L0bar_L0bar_cos_theta_star_vs_delta_eta_cuts_hist->ProjectionX(Form("proj_eta_SE_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
    float nL0barL0bar_delta_eta_cuts = L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts->Integral();

    TH1D *L0bar_L0bar_cosThetaProdPlane_ME_delta_eta_cuts = L0bar_L0bar_cos_theta_star_vs_delta_eta_ME_cuts_hist->ProjectionX(Form("proj_eta_ME_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
    L0bar_L0bar_cosThetaProdPlane_ME_delta_eta_cuts->Sumw2();
    L0bar_L0bar_cosThetaProdPlane_ME_delta_eta_cuts->Scale(nL0barL0bar_delta_eta_cuts/L0bar_L0bar_cosThetaProdPlane_ME_delta_eta_cuts->Integral());
    L0bar_L0bar_cosThetaProdPlane_ME_delta_eta_cuts->Scale(1./L0bar_L0bar_cosThetaProdPlane_ME_delta_eta_cuts->GetXaxis()->GetBinWidth(1));
    L0bar_L0bar_cosThetaProdPlane_ME_delta_eta_cuts->Fit(fit_L0bar_L0bar_after_delta_eta_cuts_ME, "i 0 r");
    //L0bar_L0bar_cosThetaProdPlane_ME->Draw("p e");


    L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts->Sumw2();
    L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts->SetMarkerStyle(20);
    L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts->SetMarkerColor(kRed);
    L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts->SetLineColor(kRed);
    L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts->GetXaxis()->SetTitle("cos(#theta*)");
    L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts->GetXaxis()->CenterTitle();
    L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts->GetYaxis()->CenterTitle();
    //L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts->Divide(L0bar_L0bar_cosThetaProdPlane_ME_delta_eta_cuts); //correct using ME
    //L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts->Scale(nL0barL0bar_delta_eta_cuts/L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts->Integral()); //scale back
    L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts->Scale(1./L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
    L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts->Fit(fit_L0bar_L0bar_after_delta_eta_cuts, "i 0 r");
    L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts->SetMinimum(0);
    //L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts->GetYaxis()->SetRangeUser(-20e9, 0);
    L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts->Draw("p e");

    L0bar_L0bar_cosThetaProdPlane_ME_delta_eta_cuts->Draw("p e same");

    fit_L0bar_L0bar_after_delta_eta_cuts->SetLineColor(1);
    fit_L0bar_L0bar_after_delta_eta_cuts->Draw("same");

    fit_L0bar_L0bar_after_delta_eta_cuts_ME->SetLineColor(kBlue);
    fit_L0bar_L0bar_after_delta_eta_cuts_ME->Draw("same");


    float L0bar_L0bar_slope_delta_eta_cuts = fit_L0bar_L0bar_after_delta_eta_cuts->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
    float L0bar_L0bar_slope_delta_eta_cuts_err = fit_L0bar_L0bar_after_delta_eta_cuts->GetParError(1)/(L0bar_alpha*L0bar_alpha);


    //statistical error correction for systematic error
    float SysErrSlope_L0bar_L0bar_delta_eta_corr = sqrt( fabs( fit_L0bar_L0bar_after_delta_eta_cuts->GetParError(1)*fit_L0bar_L0bar_after_delta_eta_cuts->GetParError(1) - fit_L0bar_L0bar_after_delta_eta_cuts_ME->GetParError(1)*fit_L0bar_L0bar_after_delta_eta_cuts_ME->GetParError(1) ) );

    float SysErrSlope_L0bar_L0bar_delta_eta_work = ( fabs( fit_L0bar_L0bar_after_delta_eta_cuts->GetParameter(1) - fit_L0bar_L0bar_after_delta_eta_cuts_ME->GetParameter(1) ) - SysErrSlope_L0bar_L0bar_delta_eta_corr )/fabs(fit_L0bar_L0bar_after_delta_eta_cuts->GetParameter(1));

    float L0bar_L0bar_slope_delta_eta_cuts_sys_err = 0;

    if( SysErrSlope_L0bar_L0bar_delta_eta_work > 0 ) L0bar_L0bar_slope_delta_eta_cuts_sys_err = SysErrSlope_L0bar_L0bar_delta_eta_work; //store sys. err. only if it is larger than statistical fluctuations

    float L0bar_L0bar_slope_delta_eta_cuts_sys_err_no_corr = fabs( fit_L0bar_L0bar_after_delta_eta_cuts->GetParameter(1) - fit_L0bar_L0bar_after_delta_eta_cuts_ME->GetParameter(1) )/fabs(fit_L0bar_L0bar_after_delta_eta_cuts->GetParameter(1));

    if(corr_err == 1) SysErrSlope_delta_eta_hist[delta_eta_bin-1]->SetBinContent(3, L0bar_L0bar_slope_delta_eta_cuts_sys_err);
    if(corr_err == 0) SysErrSlope_delta_eta_hist[delta_eta_bin-1]->SetBinContent(3, L0bar_L0bar_slope_delta_eta_cuts_sys_err_no_corr);

    TPaveText *L0bar_L0bar_text_MC_delta_eta_cuts = new TPaveText(0.5, 0.2, 0.8, 0.45, "NDC");
    L0bar_L0bar_text_MC_delta_eta_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0bar_L0bar_text_MC_delta_eta_cuts->AddText("Minimum bias");
    L0bar_L0bar_text_MC_delta_eta_cuts->AddText("#bar{#Lambda}-#bar{#Lambda}");
    L0bar_L0bar_text_MC_delta_eta_cuts->AddText("True MC");
    L0bar_L0bar_text_MC_delta_eta_cuts->AddText("Analysis cuts");
    L0bar_L0bar_text_MC_delta_eta_cuts->AddText(Form("P = %.3f #pm %.3f", L0bar_L0bar_slope_delta_eta_cuts, L0bar_L0bar_slope_delta_eta_cuts_err ));
    L0bar_L0bar_text_MC_delta_eta_cuts->AddText(Form("P_{ME} = %.3f", fit_L0bar_L0bar_after_delta_eta_cuts_ME->GetParameter(1)/(L0bar_alpha*L0bar_alpha) ));
    if( (delta_eta_bin-1) == 0 )L0bar_L0bar_text_MC_delta_eta_cuts->AddText("|#Delta#it{y}| < 0.5");
    if( (delta_eta_bin-1) == 1 )L0bar_L0bar_text_MC_delta_eta_cuts->AddText("0.5 < |#Delta#it{y}| < 2.0");
    if(corr_err == 1) L0bar_L0bar_text_MC_delta_eta_cuts->AddText(Form("Slope sys. err: %.1f %%", L0bar_L0bar_slope_delta_eta_cuts_sys_err*100 ));
    if(corr_err == 0) L0bar_L0bar_text_MC_delta_eta_cuts->AddText(Form("Slope sys. err: %.1f %%", L0bar_L0bar_slope_delta_eta_cuts_sys_err_no_corr*100 ));
    L0bar_L0bar_text_MC_delta_eta_cuts->SetFillColorAlpha(0, 0.01);
    L0bar_L0bar_text_MC_delta_eta_cuts->Draw("same");

    L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/Lbar_Lbar/L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------

    TCanvas *L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts_corr_can = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts_corr_can_%i", delta_eta_bin), Form("L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts_corr_can_%i", delta_eta_bin), 1200, 1000);
    L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts_corr_can->cd();

    TF1 *fit_L0bar_L0bar_after_delta_eta_cuts_corr = new TF1("fit_L0bar_L0bar_after_delta_eta_cuts_corr", "[0]*(1+[1]*x)", -1, 1);
    fit_L0bar_L0bar_after_delta_eta_cuts_corr->SetParameters(1000, 0.10);

    TH1D *L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts_corr = (TH1D*)L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts->Clone(Form("L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts_corr_%i", delta_eta_bin));
    L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts_corr->Divide(L0bar_L0bar_cosThetaProdPlane_ME_delta_eta_cuts); //correct
    L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts_corr->Scale(nL0barL0bar_delta_eta_cuts/L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts_corr->Integral()); //scale back
    L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts_corr->Scale(1./L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts_corr->GetXaxis()->GetBinWidth(1)); //divide by bin width
    //L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts_corr->Scale(nL0barL0bar_delta_eta_cuts);
    L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts_corr->Fit(fit_L0bar_L0bar_after_delta_eta_cuts_corr, "i 0 r");
    L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts_corr->SetMinimum(0);
    L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts_corr->Draw("p e");

    fit_L0bar_L0bar_after_delta_eta_cuts_corr->SetLineColor(1);
    fit_L0bar_L0bar_after_delta_eta_cuts_corr->Draw("same");

    float L0bar_L0bar_slope_delta_eta_cuts_corr = fit_L0bar_L0bar_after_delta_eta_cuts_corr->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
    float L0bar_L0bar_slope_delta_eta_cuts_corr_err = fit_L0bar_L0bar_after_delta_eta_cuts_corr->GetParError(1)/(L0bar_alpha*L0bar_alpha);

    ResidualPolarization_delta_eta_hist[delta_eta_bin-1]->SetBinContent(3, L0bar_L0bar_slope_delta_eta_cuts_corr);
    ResidualPolarization_delta_eta_hist[delta_eta_bin-1]->SetBinError(3, L0bar_L0bar_slope_delta_eta_cuts_corr_err);


    TPaveText *L0bar_L0bar_text_MC_delta_eta_cuts_corr = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    L0bar_L0bar_text_MC_delta_eta_cuts_corr->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0bar_L0bar_text_MC_delta_eta_cuts_corr->AddText("Minimum bias");
    L0bar_L0bar_text_MC_delta_eta_cuts_corr->AddText("#bar{#Lambda}-#bar{#Lambda}");
    L0bar_L0bar_text_MC_delta_eta_cuts_corr->AddText("True MC");
    L0bar_L0bar_text_MC_delta_eta_cuts_corr->AddText("Analysis cuts");
    L0bar_L0bar_text_MC_delta_eta_cuts_corr->AddText(Form("P = %.3f #pm %.3f", L0bar_L0bar_slope_delta_eta_cuts_corr, L0bar_L0bar_slope_delta_eta_cuts_corr_err ));
    //if(corr_err == 1) L0bar_L0bar_text_MC_delta_eta_cuts_corr->AddText(Form("Slope sys. err: %.1f %%", L0bar_L0bar_slope_delta_eta_cuts_corr_sys_err*100 ));
    //if(corr_err == 0) L0bar_L0bar_text_MC_delta_eta_cuts_corr->AddText(Form("Slope sys. err: %.1f %%", L0bar_L0bar_slope_delta_eta_cuts_corr_sys_err_no_corr*100 ));
    L0bar_L0bar_text_MC_delta_eta_cuts_corr->SetFillColorAlpha(0, 0.01);
    L0bar_L0bar_text_MC_delta_eta_cuts_corr->Draw("same");


    L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/Lbar_Lbar/L0bar_L0bar_cosThetaProdPlane_delta_eta_cuts_corr_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------

    SysErrSlope_delta_eta_hist[delta_eta_bin-1]->Write();
    ResidualPolarization_delta_eta_hist[delta_eta_bin-1]->Write();

    //_________________________________________________________________________________________________________________________________________________________________________________________________

    //after cuts
    //L-Lbar
    TCanvas *L0_L0bar_cosThetaProdPlane_delta_phi_cuts_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_delta_phi_cuts_can_%i", delta_eta_bin), Form("L0_L0bar_cosThetaProdPlane_delta_phi_cuts_can_%i", delta_eta_bin), 1200, 1000);
    L0_L0bar_cosThetaProdPlane_delta_phi_cuts_can->cd();

    TF1 *fit_L0_L0bar_after_delta_phi_cuts = new TF1("fit_L0_L0bar_after_delta_phi_cuts", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0bar_after_delta_phi_cuts->SetParameters(1000, 0.10);

    TF1 *fit_L0_L0bar_after_delta_phi_cuts_ME = new TF1("fit_L0_L0bar_after_delta_phi_cuts_ME", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0bar_after_delta_phi_cuts_ME->SetParameters(1000, 0.10);

    TH1D *L0_L0bar_cosThetaProdPlane_delta_phi_cuts = L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_cuts_hist->ProjectionX(Form("proj_phi_SE_%i", delta_eta_bin), bin_min_loop, bin_max_loop);
    float nL0L0bar_delta_phi_cuts = L0_L0bar_cosThetaProdPlane_delta_phi_cuts->Integral();

    TH1D *L0_L0bar_cosThetaProdPlane_ME_delta_phi_cuts = L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_ME_cuts_hist->ProjectionX(Form("proj_phi_ME_%i", delta_eta_bin), bin_min_loop, bin_max_loop);
    L0_L0bar_cosThetaProdPlane_ME_delta_phi_cuts->Sumw2();
    L0_L0bar_cosThetaProdPlane_ME_delta_phi_cuts->Scale(nL0L0bar_delta_phi_cuts/L0_L0bar_cosThetaProdPlane_ME_delta_phi_cuts->Integral());
    L0_L0bar_cosThetaProdPlane_ME_delta_phi_cuts->Scale(1./L0_L0bar_cosThetaProdPlane_ME_delta_phi_cuts->GetXaxis()->GetBinWidth(1));
    L0_L0bar_cosThetaProdPlane_ME_delta_phi_cuts->Fit(fit_L0_L0bar_after_delta_phi_cuts_ME, "i 0 r");
    //L0_L0bar_cosThetaProdPlane_ME->Draw("p e");

    L0_L0bar_cosThetaProdPlane_delta_phi_cuts->Sumw2();
    L0_L0bar_cosThetaProdPlane_delta_phi_cuts->SetMarkerStyle(20);
    L0_L0bar_cosThetaProdPlane_delta_phi_cuts->SetMarkerColor(kRed);
    L0_L0bar_cosThetaProdPlane_delta_phi_cuts->SetLineColor(kRed);
    L0_L0bar_cosThetaProdPlane_delta_phi_cuts->GetXaxis()->SetTitle("cos(#theta*)");
    L0_L0bar_cosThetaProdPlane_delta_phi_cuts->GetXaxis()->CenterTitle();
    L0_L0bar_cosThetaProdPlane_delta_phi_cuts->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    L0_L0bar_cosThetaProdPlane_delta_phi_cuts->GetYaxis()->CenterTitle();
    //L0_L0bar_cosThetaProdPlane_delta_phi_cuts->Divide(L0_L0bar_cosThetaProdPlane_ME_delta_phi_cuts); //correct using ME
    //L0_L0bar_cosThetaProdPlane_delta_phi_cuts->Scale(nL0L0bar_delta_phi_cuts/L0_L0bar_cosThetaProdPlane_delta_phi_cuts->Integral()); //scale back
    L0_L0bar_cosThetaProdPlane_delta_phi_cuts->Scale(1./L0_L0bar_cosThetaProdPlane_delta_phi_cuts->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
    L0_L0bar_cosThetaProdPlane_delta_phi_cuts->Fit(fit_L0_L0bar_after_delta_phi_cuts, "i 0 r");
    L0_L0bar_cosThetaProdPlane_delta_phi_cuts->SetMinimum(0);
    //L0_L0bar_cosThetaProdPlane_delta_phi_cuts->GetYaxis()->SetRangeUser(-20e9, 0);
    L0_L0bar_cosThetaProdPlane_delta_phi_cuts->Draw("p e");

    L0_L0bar_cosThetaProdPlane_ME_delta_phi_cuts->Draw("p e same");

    fit_L0_L0bar_after_delta_phi_cuts->SetLineColor(1);
    fit_L0_L0bar_after_delta_phi_cuts->Draw("same");

    fit_L0_L0bar_after_delta_phi_cuts_ME->SetLineColor(kBlue);
    fit_L0_L0bar_after_delta_phi_cuts_ME->Draw("same");

    float L0_L0bar_slope_delta_phi_cuts = fit_L0_L0bar_after_delta_phi_cuts->GetParameter(1)/(L0_alpha*L0bar_alpha);
    float L0_L0bar_slope_delta_phi_cuts_err = fit_L0_L0bar_after_delta_phi_cuts->GetParError(1)/(L0_alpha*L0bar_alpha);

    //statistical error correction for systematic error
    float SysErrSlope_L0_L0bar_delta_phi_corr = sqrt( fabs( fit_L0_L0bar_after_delta_phi_cuts->GetParError(1)*fit_L0_L0bar_after_delta_phi_cuts->GetParError(1) - fit_L0_L0bar_after_delta_phi_cuts_ME->GetParError(1)*fit_L0_L0bar_after_delta_phi_cuts_ME->GetParError(1) ) );

    float SysErrSlope_L0_L0bar_delta_phi_work = ( fabs( fit_L0_L0bar_after_delta_phi_cuts->GetParameter(1) - fit_L0_L0bar_after_delta_phi_cuts_ME->GetParameter(1) ) - SysErrSlope_L0_L0bar_delta_phi_corr )/fabs(fit_L0_L0bar_after_delta_phi_cuts->GetParameter(1));

    float L0_L0bar_slope_delta_phi_cuts_sys_err = 0;

    if( SysErrSlope_L0_L0bar_delta_phi_work > 0 ) L0_L0bar_slope_delta_phi_cuts_sys_err = SysErrSlope_L0_L0bar_delta_phi_work; //store sys. err. only if it is larger than statistical fluctuations

    float L0_L0bar_slope_delta_phi_cuts_sys_err_no_corr = fabs( fit_L0_L0bar_after_delta_phi_cuts->GetParameter(1) - fit_L0_L0bar_after_delta_phi_cuts_ME->GetParameter(1) )/fabs(fit_L0_L0bar_after_delta_phi_cuts->GetParameter(1));

    if(corr_err == 1) SysErrSlope_delta_phi_hist[delta_eta_bin-1]->SetBinContent(1, L0_L0bar_slope_delta_phi_cuts_sys_err);
    if(corr_err == 0) SysErrSlope_delta_phi_hist[delta_eta_bin-1]->SetBinContent(1, L0_L0bar_slope_delta_phi_cuts_sys_err_no_corr);

    TPaveText *L0_L0bar_text_MC_delta_phi_cuts = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    L0_L0bar_text_MC_delta_phi_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0_L0bar_text_MC_delta_phi_cuts->AddText("Minimum bias");
    L0_L0bar_text_MC_delta_phi_cuts->AddText("#Lambda-#bar{#Lambda}");
    L0_L0bar_text_MC_delta_phi_cuts->AddText("True MC");
    L0_L0bar_text_MC_delta_phi_cuts->AddText("Analysis cuts");
    L0_L0bar_text_MC_delta_phi_cuts->AddText(Form("P = %.3f #pm %.3f", L0_L0bar_slope_delta_phi_cuts, L0_L0bar_slope_delta_phi_cuts_err ));
    L0_L0bar_text_MC_delta_phi_cuts->AddText(Form("P_{ME} = %.3f", fit_L0_L0bar_after_delta_phi_cuts_ME->GetParameter(1)/(L0_alpha*L0bar_alpha) ));
    if(corr_err == 1) L0_L0bar_text_MC_delta_phi_cuts->AddText(Form("Slope sys. err: %.1f %%", L0_L0bar_slope_delta_phi_cuts_sys_err*100 ));
    if(corr_err == 0) L0_L0bar_text_MC_delta_phi_cuts->AddText(Form("Slope sys. err: %.1f %%", L0_L0bar_slope_delta_phi_cuts_sys_err_no_corr*100 ));
    L0_L0bar_text_MC_delta_phi_cuts->SetFillColorAlpha(0, 0.01);
    L0_L0bar_text_MC_delta_phi_cuts->Draw("same");

    L0_L0bar_cosThetaProdPlane_delta_phi_cuts_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/L0_L0bar_cosThetaProdPlane_delta_phi_cuts_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------

    TCanvas *L0_L0bar_cosThetaProdPlane_delta_phi_cuts_corr_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_delta_phi_cuts_corr_can_%i", delta_eta_bin), Form("L0_L0bar_cosThetaProdPlane_delta_phi_cuts_corr_can_%i", delta_eta_bin), 1200, 1000);
    L0_L0bar_cosThetaProdPlane_delta_phi_cuts_corr_can->cd();

    TF1 *fit_L0_L0bar_after_delta_phi_cuts_corr = new TF1("fit_L0_L0bar_after_delta_phi_cuts_corr", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0bar_after_delta_phi_cuts_corr->SetParameters(1000, 0.10);

    TH1D *L0_L0bar_cosThetaProdPlane_delta_phi_cuts_corr = (TH1D*)L0_L0bar_cosThetaProdPlane_delta_phi_cuts->Clone(Form("L0_L0bar_cosThetaProdPlane_delta_phi_cuts_corr_%i", delta_eta_bin));
    L0_L0bar_cosThetaProdPlane_delta_phi_cuts_corr->Divide(L0_L0bar_cosThetaProdPlane_ME_delta_phi_cuts); //correct
    L0_L0bar_cosThetaProdPlane_delta_phi_cuts_corr->Scale(nL0L0bar_delta_phi_cuts/L0_L0bar_cosThetaProdPlane_delta_phi_cuts_corr->Integral()); //scale back
    L0_L0bar_cosThetaProdPlane_delta_phi_cuts_corr->Scale(1./L0_L0bar_cosThetaProdPlane_delta_phi_cuts_corr->GetXaxis()->GetBinWidth(1)); //divide by bin width
    //L0_L0bar_cosThetaProdPlane_delta_phi_cuts_corr->Scale(nL0L0bar_delta_phi_cuts);
    L0_L0bar_cosThetaProdPlane_delta_phi_cuts_corr->Fit(fit_L0_L0bar_after_delta_phi_cuts_corr, "i 0 r");
    L0_L0bar_cosThetaProdPlane_delta_phi_cuts_corr->SetMinimum(0);
    L0_L0bar_cosThetaProdPlane_delta_phi_cuts_corr->Draw("p e");

    fit_L0_L0bar_after_delta_phi_cuts_corr->SetLineColor(1);
    fit_L0_L0bar_after_delta_phi_cuts_corr->Draw("same");

    float L0_L0bar_slope_delta_phi_cuts_corr = fit_L0_L0bar_after_delta_phi_cuts_corr->GetParameter(1)/(L0_alpha*L0bar_alpha);
    float L0_L0bar_slope_delta_phi_cuts_corr_err = fit_L0_L0bar_after_delta_phi_cuts_corr->GetParError(1)/(L0_alpha*L0bar_alpha);

    ResidualPolarization_delta_phi_hist[delta_eta_bin-1]->SetBinContent(1, L0_L0bar_slope_delta_phi_cuts_corr);
    ResidualPolarization_delta_phi_hist[delta_eta_bin-1]->SetBinError(1, L0_L0bar_slope_delta_phi_cuts_corr_err);


    TPaveText *L0_L0bar_text_MC_delta_phi_cuts_corr = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    L0_L0bar_text_MC_delta_phi_cuts_corr->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0_L0bar_text_MC_delta_phi_cuts_corr->AddText("Minimum bias");
    L0_L0bar_text_MC_delta_phi_cuts_corr->AddText("#Lambda-#bar{#Lambda}");
    L0_L0bar_text_MC_delta_phi_cuts_corr->AddText("True MC");
    L0_L0bar_text_MC_delta_phi_cuts_corr->AddText("Analysis cuts");
    L0_L0bar_text_MC_delta_phi_cuts_corr->AddText(Form("P = %.3f #pm %.3f", L0_L0bar_slope_delta_phi_cuts_corr, L0_L0bar_slope_delta_phi_cuts_corr_err ));
    //if(corr_err == 1) L0_L0bar_text_MC_delta_phi_cuts_corr->AddText(Form("Slope sys. err: %.1f %%", L0_L0bar_slope_delta_phi_cuts_corr_sys_err*100 ));
    //if(corr_err == 0) L0_L0bar_text_MC_delta_phi_cuts_corr->AddText(Form("Slope sys. err: %.1f %%", L0_L0bar_slope_delta_phi_cuts_corr_sys_err_no_corr*100 ));
    L0_L0bar_text_MC_delta_phi_cuts_corr->SetFillColorAlpha(0, 0.01);
    L0_L0bar_text_MC_delta_phi_cuts_corr->Draw("same");


    L0_L0bar_cosThetaProdPlane_delta_phi_cuts_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/L0_L0bar_cosThetaProdPlane_delta_phi_cuts_corr_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------

    //L-L
    TCanvas *L0_L0_cosThetaProdPlane_delta_phi_cuts_can = new TCanvas(Form("L0_L0_cosThetaProdPlane_delta_phi_cuts_can_%i", delta_eta_bin), Form("L0_L0_cosThetaProdPlane_delta_phi_cuts_can_%i", delta_eta_bin), 1200, 1000);
    L0_L0_cosThetaProdPlane_delta_phi_cuts_can->cd();

    TF1 *fit_L0_L0_after_delta_phi_cuts = new TF1("fit_L0_L0_after_delta_phi_cuts", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0_after_delta_phi_cuts->SetParameters(1000, 0.10);

    TF1 *fit_L0_L0_after_delta_phi_cuts_ME = new TF1("fit_L0_L0_after_delta_phi_cuts_ME", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0_after_delta_phi_cuts_ME->SetParameters(1000, 0.10);

    TH1D *L0_L0_cosThetaProdPlane_delta_phi_cuts = L0_L0_cos_theta_star_vs_delta_phi_for_corr_cuts_hist->ProjectionX(Form("proj_phi_SE_%i", delta_eta_bin), bin_min_loop, bin_max_loop);
    float nL0L0_delta_phi_cuts = L0_L0_cosThetaProdPlane_delta_phi_cuts->Integral();

    TH1D *L0_L0_cosThetaProdPlane_ME_delta_phi_cuts = L0_L0_cos_theta_star_vs_delta_phi_for_corr_ME_cuts_hist->ProjectionX(Form("proj_phi_ME_%i", delta_eta_bin), bin_min_loop, bin_max_loop);
    L0_L0_cosThetaProdPlane_ME_delta_phi_cuts->Sumw2();
    L0_L0_cosThetaProdPlane_ME_delta_phi_cuts->Scale(nL0L0_delta_phi_cuts/L0_L0_cosThetaProdPlane_ME_delta_phi_cuts->Integral());
    L0_L0_cosThetaProdPlane_ME_delta_phi_cuts->Scale(1./L0_L0_cosThetaProdPlane_ME_delta_phi_cuts->GetXaxis()->GetBinWidth(1));
    L0_L0_cosThetaProdPlane_ME_delta_phi_cuts->Fit(fit_L0_L0_after_delta_phi_cuts_ME, "i 0 r");
    //L0_L0_cosThetaProdPlane_ME->Draw("p e");


    L0_L0_cosThetaProdPlane_delta_phi_cuts->Sumw2();
    L0_L0_cosThetaProdPlane_delta_phi_cuts->SetMarkerStyle(20);
    L0_L0_cosThetaProdPlane_delta_phi_cuts->SetMarkerColor(kRed);
    L0_L0_cosThetaProdPlane_delta_phi_cuts->SetLineColor(kRed);
    L0_L0_cosThetaProdPlane_delta_phi_cuts->GetXaxis()->SetTitle("cos(#theta*)");
    L0_L0_cosThetaProdPlane_delta_phi_cuts->GetXaxis()->CenterTitle();
    L0_L0_cosThetaProdPlane_delta_phi_cuts->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    L0_L0_cosThetaProdPlane_delta_phi_cuts->GetYaxis()->CenterTitle();
    //L0_L0_cosThetaProdPlane_delta_phi_cuts->Divide(L0_L0_cosThetaProdPlane_ME_delta_phi_cuts); //correct using ME
    //L0_L0_cosThetaProdPlane_delta_phi_cuts->Scale(nL0L0_delta_phi_cuts/L0_L0_cosThetaProdPlane_delta_phi_cuts->Integral()); //scale back
    L0_L0_cosThetaProdPlane_delta_phi_cuts->Scale(1./L0_L0_cosThetaProdPlane_delta_phi_cuts->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
    L0_L0_cosThetaProdPlane_delta_phi_cuts->Fit(fit_L0_L0_after_delta_phi_cuts, "i 0 r");
    L0_L0_cosThetaProdPlane_delta_phi_cuts->SetMinimum(0);
    //L0_L0_cosThetaProdPlane_delta_phi_cuts->GetYaxis()->SetRangeUser(-20e9, 0);
    L0_L0_cosThetaProdPlane_delta_phi_cuts->Draw("p e");

    L0_L0_cosThetaProdPlane_ME_delta_phi_cuts->Draw("p e same");

    fit_L0_L0_after_delta_phi_cuts->SetLineColor(1);
    fit_L0_L0_after_delta_phi_cuts->Draw("same");

    fit_L0_L0_after_delta_phi_cuts_ME->SetLineColor(kBlue);
    fit_L0_L0_after_delta_phi_cuts_ME->Draw("same");


    float L0_L0_slope_delta_phi_cuts = fit_L0_L0_after_delta_phi_cuts->GetParameter(1)/(L0_alpha*L0_alpha);
    float L0_L0_slope_delta_phi_cuts_err = fit_L0_L0_after_delta_phi_cuts->GetParError(1)/(L0_alpha*L0_alpha);


    //statistical error correction for systematic error
    float SysErrSlope_L0_L0_delta_phi_corr = sqrt( fabs( fit_L0_L0_after_delta_phi_cuts->GetParError(1)*fit_L0_L0_after_delta_phi_cuts->GetParError(1) - fit_L0_L0_after_delta_phi_cuts_ME->GetParError(1)*fit_L0_L0_after_delta_phi_cuts_ME->GetParError(1) ) );

    float SysErrSlope_L0_L0_delta_phi_work = ( fabs( fit_L0_L0_after_delta_phi_cuts->GetParameter(1) - fit_L0_L0_after_delta_phi_cuts_ME->GetParameter(1) ) - SysErrSlope_L0_L0_delta_phi_corr )/fabs(fit_L0_L0_after_delta_phi_cuts->GetParameter(1));

    float L0_L0_slope_delta_phi_cuts_sys_err = 0;

    if( SysErrSlope_L0_L0_delta_phi_work > 0 ) L0_L0_slope_delta_phi_cuts_sys_err = SysErrSlope_L0_L0_delta_phi_work; //store sys. err. only if it is larger than statistical fluctuations

    float L0_L0_slope_delta_phi_cuts_sys_err_no_corr = fabs( fit_L0_L0_after_delta_phi_cuts->GetParameter(1) - fit_L0_L0_after_delta_phi_cuts_ME->GetParameter(1) )/fabs(fit_L0_L0_after_delta_phi_cuts->GetParameter(1));

    if(corr_err == 1) SysErrSlope_delta_phi_hist[delta_eta_bin-1]->SetBinContent(2, L0_L0_slope_delta_phi_cuts_sys_err);
    if(corr_err == 0) SysErrSlope_delta_phi_hist[delta_eta_bin-1]->SetBinContent(2, L0_L0_slope_delta_phi_cuts_sys_err_no_corr);

    TPaveText *L0_L0_text_MC_delta_phi_cuts = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    L0_L0_text_MC_delta_phi_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0_L0_text_MC_delta_phi_cuts->AddText("Minimum bias");
    L0_L0_text_MC_delta_phi_cuts->AddText("#Lambda-#Lambda");
    L0_L0_text_MC_delta_phi_cuts->AddText("True MC");
    L0_L0_text_MC_delta_phi_cuts->AddText("Analysis cuts");
    L0_L0_text_MC_delta_phi_cuts->AddText(Form("P = %.3f #pm %.3f", L0_L0_slope_delta_phi_cuts, L0_L0_slope_delta_phi_cuts_err ));
    L0_L0_text_MC_delta_phi_cuts->AddText(Form("P_{ME} = %.3f", fit_L0_L0_after_delta_phi_cuts_ME->GetParameter(1)/(L0_alpha*L0_alpha) ));
    if(corr_err == 1) L0_L0_text_MC_delta_phi_cuts->AddText(Form("Slope sys. err: %.1f %%", L0_L0_slope_delta_phi_cuts_sys_err*100 ));
    if(corr_err == 0) L0_L0_text_MC_delta_phi_cuts->AddText(Form("Slope sys. err: %.1f %%", L0_L0_slope_delta_phi_cuts_sys_err_no_corr*100 ));
    L0_L0_text_MC_delta_phi_cuts->SetFillColorAlpha(0, 0.01);
    L0_L0_text_MC_delta_phi_cuts->Draw("same");

    L0_L0_cosThetaProdPlane_delta_phi_cuts_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_L/L0_L0_cosThetaProdPlane_delta_phi_cuts_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------

    TCanvas *L0_L0_cosThetaProdPlane_delta_phi_cuts_corr_can = new TCanvas(Form("L0_L0_cosThetaProdPlane_delta_phi_cuts_corr_can_%i", delta_eta_bin), Form("L0_L0_cosThetaProdPlane_delta_phi_cuts_corr_can_%i", delta_eta_bin), 1200, 1000);
    L0_L0_cosThetaProdPlane_delta_phi_cuts_corr_can->cd();

    TF1 *fit_L0_L0_after_delta_phi_cuts_corr = new TF1("fit_L0_L0_after_delta_phi_cuts_corr", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0_after_delta_phi_cuts_corr->SetParameters(1000, 0.10);

    TH1D *L0_L0_cosThetaProdPlane_delta_phi_cuts_corr = (TH1D*)L0_L0_cosThetaProdPlane_delta_phi_cuts->Clone(Form("L0_L0_cosThetaProdPlane_delta_phi_cuts_corr_%i", delta_eta_bin));
    L0_L0_cosThetaProdPlane_delta_phi_cuts_corr->Divide(L0_L0_cosThetaProdPlane_ME_delta_phi_cuts); //correct
    L0_L0_cosThetaProdPlane_delta_phi_cuts_corr->Scale(nL0L0_delta_phi_cuts/L0_L0_cosThetaProdPlane_delta_phi_cuts_corr->Integral()); //scale back
    L0_L0_cosThetaProdPlane_delta_phi_cuts_corr->Scale(1./L0_L0_cosThetaProdPlane_delta_phi_cuts_corr->GetXaxis()->GetBinWidth(1)); //divide by bin width
    //L0_L0_cosThetaProdPlane_delta_phi_cuts_corr->Scale(nL0L0_delta_phi_cuts);
    L0_L0_cosThetaProdPlane_delta_phi_cuts_corr->Fit(fit_L0_L0_after_delta_phi_cuts_corr, "i 0 r");
    L0_L0_cosThetaProdPlane_delta_phi_cuts_corr->SetMinimum(0);
    L0_L0_cosThetaProdPlane_delta_phi_cuts_corr->Draw("p e");

    fit_L0_L0_after_delta_phi_cuts_corr->SetLineColor(1);
    fit_L0_L0_after_delta_phi_cuts_corr->Draw("same");

    float L0_L0_slope_delta_phi_cuts_corr = fit_L0_L0_after_delta_phi_cuts_corr->GetParameter(1)/(L0_alpha*L0_alpha);
    float L0_L0_slope_delta_phi_cuts_corr_err = fit_L0_L0_after_delta_phi_cuts_corr->GetParError(1)/(L0_alpha*L0_alpha);

    ResidualPolarization_delta_phi_hist[delta_eta_bin-1]->SetBinContent(2, L0_L0_slope_delta_phi_cuts_corr);
    ResidualPolarization_delta_phi_hist[delta_eta_bin-1]->SetBinError(2, L0_L0_slope_delta_phi_cuts_corr_err);


    TPaveText *L0_L0_text_MC_delta_phi_cuts_corr = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    L0_L0_text_MC_delta_phi_cuts_corr->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0_L0_text_MC_delta_phi_cuts_corr->AddText("Minimum bias");
    L0_L0_text_MC_delta_phi_cuts_corr->AddText("#Lambda-#Lambda");
    L0_L0_text_MC_delta_phi_cuts_corr->AddText("True MC");
    L0_L0_text_MC_delta_phi_cuts_corr->AddText("Analysis cuts");
    L0_L0_text_MC_delta_phi_cuts_corr->AddText(Form("P = %.3f #pm %.3f", L0_L0_slope_delta_phi_cuts_corr, L0_L0_slope_delta_phi_cuts_corr_err ));
    //if(corr_err == 1) L0_L0_text_MC_delta_phi_cuts_corr->AddText(Form("Slope sys. err: %.1f %%", L0_L0_slope_delta_phi_cuts_corr_sys_err*100 ));
    //if(corr_err == 0) L0_L0_text_MC_delta_phi_cuts_corr->AddText(Form("Slope sys. err: %.1f %%", L0_L0_slope_delta_phi_cuts_corr_sys_err_no_corr*100 ));
    L0_L0_text_MC_delta_phi_cuts_corr->SetFillColorAlpha(0, 0.01);
    L0_L0_text_MC_delta_phi_cuts_corr->Draw("same");


    L0_L0_cosThetaProdPlane_delta_phi_cuts_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_L/L0_L0_cosThetaProdPlane_delta_phi_cuts_corr_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------

    //Lbar-Lbar
    TCanvas *L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts_can = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts_can_%i", delta_eta_bin), Form("L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts_can_%i", delta_eta_bin), 1200, 1000);
    L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts_can->cd();

    TF1 *fit_L0bar_L0bar_after_delta_phi_cuts = new TF1("fit_L0bar_L0bar_after_delta_phi_cuts", "[0]*(1+[1]*x)", -1, 1);
    fit_L0bar_L0bar_after_delta_phi_cuts->SetParameters(1000, 0.10);

    TF1 *fit_L0bar_L0bar_after_delta_phi_cuts_ME = new TF1("fit_L0bar_L0bar_after_delta_phi_cuts_ME", "[0]*(1+[1]*x)", -1, 1);
    fit_L0bar_L0bar_after_delta_phi_cuts_ME->SetParameters(1000, 0.10);

    TH1D *L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts = L0bar_L0bar_cos_theta_star_vs_delta_phi_for_corr_cuts_hist->ProjectionX(Form("proj_phi_SE_%i", delta_eta_bin), bin_min_loop, bin_max_loop);
    float nL0barL0bar_delta_phi_cuts = L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts->Integral();

    TH1D *L0bar_L0bar_cosThetaProdPlane_ME_delta_phi_cuts = L0bar_L0bar_cos_theta_star_vs_delta_phi_for_corr_ME_cuts_hist->ProjectionX(Form("proj_phi_ME_%i", delta_eta_bin), bin_min_loop, bin_max_loop);
    L0bar_L0bar_cosThetaProdPlane_ME_delta_phi_cuts->Sumw2();
    L0bar_L0bar_cosThetaProdPlane_ME_delta_phi_cuts->Scale(nL0barL0bar_delta_phi_cuts/L0bar_L0bar_cosThetaProdPlane_ME_delta_phi_cuts->Integral());
    L0bar_L0bar_cosThetaProdPlane_ME_delta_phi_cuts->Scale(1./L0bar_L0bar_cosThetaProdPlane_ME_delta_phi_cuts->GetXaxis()->GetBinWidth(1));
    L0bar_L0bar_cosThetaProdPlane_ME_delta_phi_cuts->Fit(fit_L0bar_L0bar_after_delta_phi_cuts_ME, "i 0 r");
    //L0bar_L0bar_cosThetaProdPlane_ME->Draw("p e");


    L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts->Sumw2();
    L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts->SetMarkerStyle(20);
    L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts->SetMarkerColor(kRed);
    L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts->SetLineColor(kRed);
    L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts->GetXaxis()->SetTitle("cos(#theta*)");
    L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts->GetXaxis()->CenterTitle();
    L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts->GetYaxis()->CenterTitle();
    //L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts->Divide(L0bar_L0bar_cosThetaProdPlane_ME_delta_phi_cuts); //correct using ME
    //L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts->Scale(nL0barL0bar_delta_phi_cuts/L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts->Integral()); //scale back
    L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts->Scale(1./L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
    L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts->Fit(fit_L0bar_L0bar_after_delta_phi_cuts, "i 0 r");
    L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts->SetMinimum(0);
    //L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts->GetYaxis()->SetRangeUser(-20e9, 0);
    L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts->Draw("p e");

    L0bar_L0bar_cosThetaProdPlane_ME_delta_phi_cuts->Draw("p e same");

    fit_L0bar_L0bar_after_delta_phi_cuts->SetLineColor(1);
    fit_L0bar_L0bar_after_delta_phi_cuts->Draw("same");

    fit_L0bar_L0bar_after_delta_phi_cuts_ME->SetLineColor(kBlue);
    fit_L0bar_L0bar_after_delta_phi_cuts_ME->Draw("same");


    float L0bar_L0bar_slope_delta_phi_cuts = fit_L0bar_L0bar_after_delta_phi_cuts->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
    float L0bar_L0bar_slope_delta_phi_cuts_err = fit_L0bar_L0bar_after_delta_phi_cuts->GetParError(1)/(L0bar_alpha*L0bar_alpha);


    //statistical error correction for systematic error
    float SysErrSlope_L0bar_L0bar_delta_phi_corr = sqrt( fabs( fit_L0bar_L0bar_after_delta_phi_cuts->GetParError(1)*fit_L0bar_L0bar_after_delta_phi_cuts->GetParError(1) - fit_L0bar_L0bar_after_delta_phi_cuts_ME->GetParError(1)*fit_L0bar_L0bar_after_delta_phi_cuts_ME->GetParError(1) ) );

    float SysErrSlope_L0bar_L0bar_delta_phi_work = ( fabs( fit_L0bar_L0bar_after_delta_phi_cuts->GetParameter(1) - fit_L0bar_L0bar_after_delta_phi_cuts_ME->GetParameter(1) ) - SysErrSlope_L0bar_L0bar_delta_phi_corr )/fabs(fit_L0bar_L0bar_after_delta_phi_cuts->GetParameter(1));

    float L0bar_L0bar_slope_delta_phi_cuts_sys_err = 0;

    if( SysErrSlope_L0bar_L0bar_delta_phi_work > 0 ) L0bar_L0bar_slope_delta_phi_cuts_sys_err = SysErrSlope_L0bar_L0bar_delta_phi_work; //store sys. err. only if it is larger than statistical fluctuations

    float L0bar_L0bar_slope_delta_phi_cuts_sys_err_no_corr = fabs( fit_L0bar_L0bar_after_delta_phi_cuts->GetParameter(1) - fit_L0bar_L0bar_after_delta_phi_cuts_ME->GetParameter(1) )/fabs(fit_L0bar_L0bar_after_delta_phi_cuts->GetParameter(1));

    if(corr_err == 1) SysErrSlope_delta_phi_hist[delta_eta_bin-1]->SetBinContent(3, L0bar_L0bar_slope_delta_phi_cuts_sys_err);
    if(corr_err == 0) SysErrSlope_delta_phi_hist[delta_eta_bin-1]->SetBinContent(3, L0bar_L0bar_slope_delta_phi_cuts_sys_err_no_corr);

    TPaveText *L0bar_L0bar_text_MC_delta_phi_cuts = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    L0bar_L0bar_text_MC_delta_phi_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0bar_L0bar_text_MC_delta_phi_cuts->AddText("Minimum bias");
    L0bar_L0bar_text_MC_delta_phi_cuts->AddText("#bar{#Lambda}-#bar{#Lambda}");
    L0bar_L0bar_text_MC_delta_phi_cuts->AddText("True MC");
    L0bar_L0bar_text_MC_delta_phi_cuts->AddText("Analysis cuts");
    L0bar_L0bar_text_MC_delta_phi_cuts->AddText(Form("P = %.3f #pm %.3f", L0bar_L0bar_slope_delta_phi_cuts, L0bar_L0bar_slope_delta_phi_cuts_err ));
    L0bar_L0bar_text_MC_delta_phi_cuts->AddText(Form("P_{ME} = %.3f", fit_L0bar_L0bar_after_delta_phi_cuts_ME->GetParameter(1)/(L0bar_alpha*L0bar_alpha) ));
    if(corr_err == 1) L0bar_L0bar_text_MC_delta_phi_cuts->AddText(Form("Slope sys. err: %.1f %%", L0bar_L0bar_slope_delta_phi_cuts_sys_err*100 ));
    if(corr_err == 0) L0bar_L0bar_text_MC_delta_phi_cuts->AddText(Form("Slope sys. err: %.1f %%", L0bar_L0bar_slope_delta_phi_cuts_sys_err_no_corr*100 ));
    L0bar_L0bar_text_MC_delta_phi_cuts->SetFillColorAlpha(0, 0.01);
    L0bar_L0bar_text_MC_delta_phi_cuts->Draw("same");

    L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/Lbar_Lbar/L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------

    TCanvas *L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts_corr_can = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts_corr_can_%i", delta_eta_bin), Form("L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts_corr_can_%i", delta_eta_bin), 1200, 1000);
    L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts_corr_can->cd();

    TF1 *fit_L0bar_L0bar_after_delta_phi_cuts_corr = new TF1("fit_L0bar_L0bar_after_delta_phi_cuts_corr", "[0]*(1+[1]*x)", -1, 1);
    fit_L0bar_L0bar_after_delta_phi_cuts_corr->SetParameters(1000, 0.10);

    TH1D *L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts_corr = (TH1D*)L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts->Clone(Form("L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts_corr_%i", delta_eta_bin));
    L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts_corr->Divide(L0bar_L0bar_cosThetaProdPlane_ME_delta_phi_cuts); //correct
    L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts_corr->Scale(nL0barL0bar_delta_phi_cuts/L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts_corr->Integral()); //scale back
    L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts_corr->Scale(1./L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts_corr->GetXaxis()->GetBinWidth(1)); //divide by bin width
    //L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts_corr->Scale(nL0barL0bar_delta_phi_cuts);
    L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts_corr->Fit(fit_L0bar_L0bar_after_delta_phi_cuts_corr, "i 0 r");
    L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts_corr->SetMinimum(0);
    L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts_corr->Draw("p e");

    fit_L0bar_L0bar_after_delta_phi_cuts_corr->SetLineColor(1);
    fit_L0bar_L0bar_after_delta_phi_cuts_corr->Draw("same");

    float L0bar_L0bar_slope_delta_phi_cuts_corr = fit_L0bar_L0bar_after_delta_phi_cuts_corr->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
    float L0bar_L0bar_slope_delta_phi_cuts_corr_err = fit_L0bar_L0bar_after_delta_phi_cuts_corr->GetParError(1)/(L0bar_alpha*L0bar_alpha);

    ResidualPolarization_delta_phi_hist[delta_eta_bin-1]->SetBinContent(3, L0bar_L0bar_slope_delta_phi_cuts_corr);
    ResidualPolarization_delta_phi_hist[delta_eta_bin-1]->SetBinError(3, L0bar_L0bar_slope_delta_phi_cuts_corr_err);


    TPaveText *L0bar_L0bar_text_MC_delta_phi_cuts_corr = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    L0bar_L0bar_text_MC_delta_phi_cuts_corr->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0bar_L0bar_text_MC_delta_phi_cuts_corr->AddText("Minimum bias");
    L0bar_L0bar_text_MC_delta_phi_cuts_corr->AddText("#bar{#Lambda}-#bar{#Lambda}");
    L0bar_L0bar_text_MC_delta_phi_cuts_corr->AddText("True MC");
    L0bar_L0bar_text_MC_delta_phi_cuts_corr->AddText("Analysis cuts");
    L0bar_L0bar_text_MC_delta_phi_cuts_corr->AddText(Form("P = %.3f #pm %.3f", L0bar_L0bar_slope_delta_phi_cuts_corr, L0bar_L0bar_slope_delta_phi_cuts_corr_err ));
    //if(corr_err == 1) L0bar_L0bar_text_MC_delta_phi_cuts_corr->AddText(Form("Slope sys. err: %.1f %%", L0bar_L0bar_slope_delta_phi_cuts_corr_sys_err*100 ));
    //if(corr_err == 0) L0bar_L0bar_text_MC_delta_phi_cuts_corr->AddText(Form("Slope sys. err: %.1f %%", L0bar_L0bar_slope_delta_phi_cuts_corr_sys_err_no_corr*100 ));
    L0bar_L0bar_text_MC_delta_phi_cuts_corr->SetFillColorAlpha(0, 0.01);
    L0bar_L0bar_text_MC_delta_phi_cuts_corr->Draw("same");


    L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/Lbar_Lbar/L0bar_L0bar_cosThetaProdPlane_delta_phi_cuts_corr_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------

    SysErrSlope_delta_phi_hist[delta_eta_bin-1]->Write();
    ResidualPolarization_delta_phi_hist[delta_eta_bin-1]->Write();

    //_________________________________________________________________________________________________________________________________________________________________________________________________

    sysErrFile->cd();

    //after cuts
    //L-Lbar
    TCanvas *L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_can_%i", delta_eta_bin), Form("L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_can_%i", delta_eta_bin), 1200, 1000);
    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_can->cd();

    TF1 *fit_L0_L0bar_after_delta_eta_delta_phi_cuts = new TF1("fit_L0_L0bar_after_delta_eta_delta_phi_cuts", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0bar_after_delta_eta_delta_phi_cuts->SetParameters(1000, 0.10);

    TF1 *fit_L0_L0bar_after_delta_eta_delta_phi_cuts_ME = new TF1("fit_L0_L0bar_after_delta_eta_delta_phi_cuts_ME", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0bar_after_delta_eta_delta_phi_cuts_ME->SetParameters(1000, 0.10);

    TH1D *L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts = L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_cuts_hist->ProjectionX(Form("proj_eta_SE_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
    float nL0L0bar_delta_eta_delta_phi_cuts = L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->Integral();

    TH1D *L0_L0bar_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts = L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_ME_cuts_hist->ProjectionX(Form("proj_eta_ME_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
    L0_L0bar_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts->Sumw2();
    L0_L0bar_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts->SetMarkerStyle(24);
    L0_L0bar_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts->SetMarkerSize(1.5);
    L0_L0bar_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts->SetMarkerColor(1);
    L0_L0bar_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts->SetLineColor(1);
    L0_L0bar_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts->Scale(nL0L0bar_delta_eta_delta_phi_cuts/L0_L0bar_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts->Integral());
    L0_L0bar_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts->Scale(1./L0_L0bar_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts->GetXaxis()->GetBinWidth(1));
    L0_L0bar_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts->Fit(fit_L0_L0bar_after_delta_eta_delta_phi_cuts_ME, "i 0 r");
    //L0_L0bar_cosThetaProdPlane_ME->Draw("p e");


    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->Sumw2();
    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->SetMarkerStyle(20);
    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->SetMarkerSize(1.5);
    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->SetMarkerColor(kRed);
    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->SetLineColor(kRed);
    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->GetXaxis()->SetTitle("cos(#theta*)");
    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->GetXaxis()->CenterTitle();
    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->GetYaxis()->CenterTitle();
    //L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->Divide(L0_L0bar_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts); //correct using ME
    //L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->Scale(nL0L0bar_delta_eta_delta_phi_cuts/L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->Integral()); //scale back
    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->Scale(1./L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->Fit(fit_L0_L0bar_after_delta_eta_delta_phi_cuts, "i 0 r");
    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->SetMinimum(0);
    //L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->GetYaxis()->SetRangeUser(-20e9, 0);
    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->Draw("p e");

    L0_L0bar_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts->Draw("p e same");

    fit_L0_L0bar_after_delta_eta_delta_phi_cuts->SetLineColor(1);
    fit_L0_L0bar_after_delta_eta_delta_phi_cuts->Draw("same");

    fit_L0_L0bar_after_delta_eta_delta_phi_cuts_ME->SetLineColor(1);
    fit_L0_L0bar_after_delta_eta_delta_phi_cuts_ME->SetLineStyle(7);
    fit_L0_L0bar_after_delta_eta_delta_phi_cuts_ME->Draw("same");

    TLegend *L0_L0bar_leg = new TLegend(0.15, 0.2, 0.45, 0.45);
    L0_L0bar_leg->AddEntry(L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts, "Same-event");
    L0_L0bar_leg->AddEntry(fit_L0_L0bar_after_delta_eta_delta_phi_cuts, "Same-event fit");
    L0_L0bar_leg->AddEntry(L0_L0bar_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts, "Mixed-event");
    L0_L0bar_leg->AddEntry(fit_L0_L0bar_after_delta_eta_delta_phi_cuts_ME, "Mixed-event fit");
    //L0_L0bar_leg->AddEntry(fitL0_L0bar_US_ThetaStar_no_corr, "Linear fit to US");
    L0_L0bar_leg->SetBorderSize(0);
    L0_L0bar_leg->SetFillColorAlpha(0, 0.01);
    L0_L0bar_leg->Draw("same");


    float L0_L0bar_slope_delta_eta_delta_phi_cuts = fit_L0_L0bar_after_delta_eta_delta_phi_cuts->GetParameter(1)/(L0_alpha*L0bar_alpha);
    float L0_L0bar_slope_delta_eta_delta_phi_cuts_err = fit_L0_L0bar_after_delta_eta_delta_phi_cuts->GetParError(1)/(L0_alpha*L0bar_alpha);

    //statistical error correction for systematic error
    float SysErrSlope_L0_L0bar_delta_eta_delta_phi_corr = sqrt( fabs( fit_L0_L0bar_after_delta_eta_delta_phi_cuts->GetParError(1)*fit_L0_L0bar_after_delta_eta_delta_phi_cuts->GetParError(1) - fit_L0_L0bar_after_delta_eta_delta_phi_cuts_ME->GetParError(1)*fit_L0_L0bar_after_delta_eta_delta_phi_cuts_ME->GetParError(1) ) );

    float SysErrSlope_L0_L0bar_delta_eta_delta_phi_work = ( fabs( fit_L0_L0bar_after_delta_eta_delta_phi_cuts->GetParameter(1) - fit_L0_L0bar_after_delta_eta_delta_phi_cuts_ME->GetParameter(1) ) - SysErrSlope_L0_L0bar_delta_eta_delta_phi_corr )/fabs(fit_L0_L0bar_after_delta_eta_delta_phi_cuts->GetParameter(1));

    float L0_L0bar_slope_delta_eta_delta_phi_cuts_sys_err = 0;

    if( SysErrSlope_L0_L0bar_delta_eta_delta_phi_work > 0 ) L0_L0bar_slope_delta_eta_delta_phi_cuts_sys_err = SysErrSlope_L0_L0bar_delta_eta_delta_phi_work; //store sys. err. only if it is larger than statistical fluctuations

    float L0_L0bar_slope_delta_eta_delta_phi_cuts_sys_err_no_corr = fabs( fit_L0_L0bar_after_delta_eta_delta_phi_cuts->GetParameter(1) - fit_L0_L0bar_after_delta_eta_delta_phi_cuts_ME->GetParameter(1) )/fabs(fit_L0_L0bar_after_delta_eta_delta_phi_cuts->GetParameter(1));

    if(corr_err == 1) SysErrSlope_delta_eta_delta_phi_hist[delta_eta_bin-1]->SetBinContent(1, L0_L0bar_slope_delta_eta_delta_phi_cuts_sys_err);
    if(corr_err == 0) SysErrSlope_delta_eta_delta_phi_hist[delta_eta_bin-1]->SetBinContent(1, L0_L0bar_slope_delta_eta_delta_phi_cuts_sys_err_no_corr);

    TPaveText *L0_L0bar_text_MC_delta_eta_delta_phi_cuts = new TPaveText(0.5, 0.2, 0.8, 0.5, "NDC");
    L0_L0bar_text_MC_delta_eta_delta_phi_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0_L0bar_text_MC_delta_eta_delta_phi_cuts->AddText("Minimum bias");
    L0_L0bar_text_MC_delta_eta_delta_phi_cuts->AddText("#Lambda-#bar{#Lambda}");
    L0_L0bar_text_MC_delta_eta_delta_phi_cuts->AddText("True MC");
    L0_L0bar_text_MC_delta_eta_delta_phi_cuts->AddText("Analysis cuts");
    L0_L0bar_text_MC_delta_eta_delta_phi_cuts->AddText(Form("P = %.3f #pm %.3f", L0_L0bar_slope_delta_eta_delta_phi_cuts, L0_L0bar_slope_delta_eta_delta_phi_cuts_err ));
    L0_L0bar_text_MC_delta_eta_delta_phi_cuts->AddText(Form("P_{ME} = %.3f", fit_L0_L0bar_after_delta_eta_delta_phi_cuts_ME->GetParameter(1)/(L0_alpha*L0bar_alpha) ));
    if( (delta_eta_bin-1) == 0 )L0_L0bar_text_MC_delta_eta_delta_phi_cuts->AddText("|#Delta#it{y}| < 0.5, |#Delta#phi| < #pi/3");
    if( (delta_eta_bin-1) == 1 )L0_L0bar_text_MC_delta_eta_delta_phi_cuts->AddText("0.5 < |#Delta#it{y}| < 2.0 or |#Delta#phi| > #pi/3");
    //if(corr_err == 1) L0_L0bar_text_MC_delta_eta_delta_phi_cuts->AddText(Form("Slope sys. err: %.1f %%", L0_L0bar_slope_delta_eta_delta_phi_cuts_sys_err*100 ));
    //if(corr_err == 0) L0_L0bar_text_MC_delta_eta_delta_phi_cuts->AddText(Form("Slope sys. err: %.1f %%", L0_L0bar_slope_delta_eta_delta_phi_cuts_sys_err_no_corr*100 ));
    L0_L0bar_text_MC_delta_eta_delta_phi_cuts->SetFillColorAlpha(0, 0.01);
    L0_L0bar_text_MC_delta_eta_delta_phi_cuts->Draw("same");

    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_%i.png", delta_eta_bin));

    //-----------------------------------------------------------

    TCanvas *L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr_can_%i", delta_eta_bin), Form("L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr_can_%i", delta_eta_bin), 1200, 1000);
    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr_can->cd();

    TF1 *fit_L0_L0bar_after_delta_eta_delta_phi_cuts_corr = new TF1("fit_L0_L0bar_after_delta_eta_delta_phi_cuts_corr", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0bar_after_delta_eta_delta_phi_cuts_corr->SetParameters(1000, 0.10);

    TH1D *L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr = (TH1D*)L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->Clone(Form("L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr_%i", delta_eta_bin));
    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr->Divide(L0_L0bar_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts); //correct
    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr->Scale(nL0L0bar_delta_eta_delta_phi_cuts/L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr->Integral()); //scale back
    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr->Scale(1./L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr->GetXaxis()->GetBinWidth(1)); //divide by bin width
    //L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr->Scale(nL0L0bar_delta_eta_delta_phi_cuts);
    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr->Fit(fit_L0_L0bar_after_delta_eta_delta_phi_cuts_corr, "i 0 r");
    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr->SetMinimum(0);
    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr->Draw("p e");

    fit_L0_L0bar_after_delta_eta_delta_phi_cuts_corr->SetLineColor(1);
    fit_L0_L0bar_after_delta_eta_delta_phi_cuts_corr->Draw("same");

    float L0_L0bar_slope_delta_eta_delta_phi_cuts_corr = fit_L0_L0bar_after_delta_eta_delta_phi_cuts_corr->GetParameter(1)/(L0_alpha*L0bar_alpha);
    float L0_L0bar_slope_delta_eta_delta_phi_cuts_corr_err = fit_L0_L0bar_after_delta_eta_delta_phi_cuts_corr->GetParError(1)/(L0_alpha*L0bar_alpha);

    ResidualPolarization_delta_eta_delta_phi_hist[delta_eta_bin-1]->SetBinContent(1, L0_L0bar_slope_delta_eta_delta_phi_cuts_corr);
    ResidualPolarization_delta_eta_delta_phi_hist[delta_eta_bin-1]->SetBinError(1, L0_L0bar_slope_delta_eta_delta_phi_cuts_corr_err);


    TPaveText *L0_L0bar_text_MC_delta_eta_delta_phi_cuts_corr = new TPaveText(0.5, 0.2, 0.8, 0.5, "NDC");
    L0_L0bar_text_MC_delta_eta_delta_phi_cuts_corr->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0_L0bar_text_MC_delta_eta_delta_phi_cuts_corr->AddText("Minimum bias");
    L0_L0bar_text_MC_delta_eta_delta_phi_cuts_corr->AddText("#Lambda-#bar{#Lambda}");
    L0_L0bar_text_MC_delta_eta_delta_phi_cuts_corr->AddText("True MC");
    L0_L0bar_text_MC_delta_eta_delta_phi_cuts_corr->AddText("Analysis cuts");
    L0_L0bar_text_MC_delta_eta_delta_phi_cuts_corr->AddText(Form("P = %.3f #pm %.3f", L0_L0bar_slope_delta_eta_delta_phi_cuts_corr, L0_L0bar_slope_delta_eta_delta_phi_cuts_corr_err ));
    if( (delta_eta_bin-1) == 0 )L0_L0bar_text_MC_delta_eta_delta_phi_cuts_corr->AddText("|#Delta#it{y}| < 0.5, |#Delta#phi| < #pi/3");
    if( (delta_eta_bin-1) == 1 )L0_L0bar_text_MC_delta_eta_delta_phi_cuts_corr->AddText("0.5 < |#Delta#it{y}| < 2.0 or |#Delta#phi| > #pi/3");
    //if(corr_err == 1) L0_L0bar_text_MC_delta_eta_delta_phi_cuts_corr->AddText(Form("Slope sys. err: %.1f %%", L0_L0bar_slope_delta_eta_delta_phi_cuts_corr_sys_err*100 ));
    //if(corr_err == 0) L0_L0bar_text_MC_delta_eta_delta_phi_cuts_corr->AddText(Form("Slope sys. err: %.1f %%", L0_L0bar_slope_delta_eta_delta_phi_cuts_corr_sys_err_no_corr*100 ));
    L0_L0bar_text_MC_delta_eta_delta_phi_cuts_corr->SetFillColorAlpha(0, 0.01);
    L0_L0bar_text_MC_delta_eta_delta_phi_cuts_corr->Draw("same");


    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------

    //L-L
    TCanvas *L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts_can = new TCanvas(Form("L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts_can_%i", delta_eta_bin), Form("L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts_can_%i", delta_eta_bin), 1200, 1000);
    L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts_can->cd();

    TF1 *fit_L0_L0_after_delta_eta_delta_phi_cuts = new TF1("fit_L0_L0_after_delta_eta_delta_phi_cuts", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0_after_delta_eta_delta_phi_cuts->SetParameters(1000, 0.10);

    TF1 *fit_L0_L0_after_delta_eta_delta_phi_cuts_ME = new TF1("fit_L0_L0_after_delta_eta_delta_phi_cuts_ME", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0_after_delta_eta_delta_phi_cuts_ME->SetParameters(1000, 0.10);

    TH1D *L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts = L0_L0_cos_theta_star_vs_delta_eta_delta_phi_cuts_hist->ProjectionX(Form("proj_eta_SE_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
    float nL0L0_delta_eta_delta_phi_cuts = L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts->Integral();

    TH1D *L0_L0_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts = L0_L0_cos_theta_star_vs_delta_eta_delta_phi_ME_cuts_hist->ProjectionX(Form("proj_eta_ME_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
    L0_L0_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts->Sumw2();
    L0_L0_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts->Scale(nL0L0_delta_eta_delta_phi_cuts/L0_L0_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts->Integral());
    L0_L0_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts->Scale(1./L0_L0_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts->GetXaxis()->GetBinWidth(1));
    L0_L0_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts->Fit(fit_L0_L0_after_delta_eta_delta_phi_cuts_ME, "i 0 r");
    //L0_L0_cosThetaProdPlane_ME->Draw("p e");


    L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts->Sumw2();
    L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts->SetMarkerStyle(20);
    L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts->SetMarkerColor(kRed);
    L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts->SetLineColor(kRed);
    L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts->GetXaxis()->SetTitle("cos(#theta*)");
    L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts->GetXaxis()->CenterTitle();
    L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts->GetYaxis()->CenterTitle();
    //L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts->Divide(L0_L0_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts); //correct using ME
    //L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts->Scale(nL0L0_delta_eta_delta_phi_cuts/L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts->Integral()); //scale back
    L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts->Scale(1./L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
    L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts->Fit(fit_L0_L0_after_delta_eta_delta_phi_cuts, "i 0 r");
    L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts->SetMinimum(0);
    //L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts->GetYaxis()->SetRangeUser(-20e9, 0);
    L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts->Draw("p e");

    L0_L0_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts->Draw("p e same");

    fit_L0_L0_after_delta_eta_delta_phi_cuts->SetLineColor(1);
    fit_L0_L0_after_delta_eta_delta_phi_cuts->Draw("same");

    fit_L0_L0_after_delta_eta_delta_phi_cuts_ME->SetLineColor(kBlue);
    fit_L0_L0_after_delta_eta_delta_phi_cuts_ME->Draw("same");


    float L0_L0_slope_delta_eta_delta_phi_cuts = fit_L0_L0_after_delta_eta_delta_phi_cuts->GetParameter(1)/(L0_alpha*L0_alpha);
    float L0_L0_slope_delta_eta_delta_phi_cuts_err = fit_L0_L0_after_delta_eta_delta_phi_cuts->GetParError(1)/(L0_alpha*L0_alpha);


    //statistical error correction for systematic error
    float SysErrSlope_L0_L0_delta_eta_delta_phi_corr = sqrt( fabs( fit_L0_L0_after_delta_eta_delta_phi_cuts->GetParError(1)*fit_L0_L0_after_delta_eta_delta_phi_cuts->GetParError(1) - fit_L0_L0_after_delta_eta_delta_phi_cuts_ME->GetParError(1)*fit_L0_L0_after_delta_eta_delta_phi_cuts_ME->GetParError(1) ) );

    float SysErrSlope_L0_L0_delta_eta_delta_phi_work = ( fabs( fit_L0_L0_after_delta_eta_delta_phi_cuts->GetParameter(1) - fit_L0_L0_after_delta_eta_delta_phi_cuts_ME->GetParameter(1) ) - SysErrSlope_L0_L0_delta_eta_delta_phi_corr )/fabs(fit_L0_L0_after_delta_eta_delta_phi_cuts->GetParameter(1));

    float L0_L0_slope_delta_eta_delta_phi_cuts_sys_err = 0;

    if( SysErrSlope_L0_L0_delta_eta_delta_phi_work > 0 ) L0_L0_slope_delta_eta_delta_phi_cuts_sys_err = SysErrSlope_L0_L0_delta_eta_delta_phi_work; //store sys. err. only if it is larger than statistical fluctuations

    float L0_L0_slope_delta_eta_delta_phi_cuts_sys_err_no_corr = fabs( fit_L0_L0_after_delta_eta_delta_phi_cuts->GetParameter(1) - fit_L0_L0_after_delta_eta_delta_phi_cuts_ME->GetParameter(1) )/fabs(fit_L0_L0_after_delta_eta_delta_phi_cuts->GetParameter(1));

    if(corr_err == 1) SysErrSlope_delta_eta_delta_phi_hist[delta_eta_bin-1]->SetBinContent(2, L0_L0_slope_delta_eta_delta_phi_cuts_sys_err);
    if(corr_err == 0) SysErrSlope_delta_eta_delta_phi_hist[delta_eta_bin-1]->SetBinContent(2, L0_L0_slope_delta_eta_delta_phi_cuts_sys_err_no_corr);

    TPaveText *L0_L0_text_MC_delta_eta_delta_phi_cuts = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    L0_L0_text_MC_delta_eta_delta_phi_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0_L0_text_MC_delta_eta_delta_phi_cuts->AddText("Minimum bias");
    L0_L0_text_MC_delta_eta_delta_phi_cuts->AddText("#Lambda-#Lambda");
    L0_L0_text_MC_delta_eta_delta_phi_cuts->AddText("True MC");
    L0_L0_text_MC_delta_eta_delta_phi_cuts->AddText("Analysis cuts");
    L0_L0_text_MC_delta_eta_delta_phi_cuts->AddText(Form("P = %.3f #pm %.3f", L0_L0_slope_delta_eta_delta_phi_cuts, L0_L0_slope_delta_eta_delta_phi_cuts_err ));
    L0_L0_text_MC_delta_eta_delta_phi_cuts->AddText(Form("P_{ME} = %.3f", fit_L0_L0_after_delta_eta_delta_phi_cuts_ME->GetParameter(1)/(L0_alpha*L0_alpha) ));
    if(corr_err == 1) L0_L0_text_MC_delta_eta_delta_phi_cuts->AddText(Form("Slope sys. err: %.1f %%", L0_L0_slope_delta_eta_delta_phi_cuts_sys_err*100 ));
    if(corr_err == 0) L0_L0_text_MC_delta_eta_delta_phi_cuts->AddText(Form("Slope sys. err: %.1f %%", L0_L0_slope_delta_eta_delta_phi_cuts_sys_err_no_corr*100 ));
    L0_L0_text_MC_delta_eta_delta_phi_cuts->SetFillColorAlpha(0, 0.01);
    L0_L0_text_MC_delta_eta_delta_phi_cuts->Draw("same");

    L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_L/L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts_%i.png", delta_eta_bin));

    //-----------------------------------------------------------

    TCanvas *L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr_can = new TCanvas(Form("L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr_can_%i", delta_eta_bin), Form("L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr_can_%i", delta_eta_bin), 1200, 1000);
    L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr_can->cd();

    TF1 *fit_L0_L0_after_delta_eta_delta_phi_cuts_corr = new TF1("fit_L0_L0_after_delta_eta_delta_phi_cuts_corr", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0_after_delta_eta_delta_phi_cuts_corr->SetParameters(1000, 0.10);

    TH1D *L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr = (TH1D*)L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts->Clone(Form("L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr_%i", delta_eta_bin));
    L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr->Divide(L0_L0_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts); //correct
    L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr->Scale(nL0L0_delta_eta_delta_phi_cuts/L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr->Integral()); //scale back
    L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr->Scale(1./L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr->GetXaxis()->GetBinWidth(1)); //divide by bin width
    //L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr->Scale(nL0L0_delta_eta_delta_phi_cuts);
    L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr->Fit(fit_L0_L0_after_delta_eta_delta_phi_cuts_corr, "i 0 r");
    L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr->SetMinimum(0);
    L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr->Draw("p e");

    fit_L0_L0_after_delta_eta_delta_phi_cuts_corr->SetLineColor(1);
    fit_L0_L0_after_delta_eta_delta_phi_cuts_corr->Draw("same");

    float L0_L0_slope_delta_eta_delta_phi_cuts_corr = fit_L0_L0_after_delta_eta_delta_phi_cuts_corr->GetParameter(1)/(L0_alpha*L0_alpha);
    float L0_L0_slope_delta_eta_delta_phi_cuts_corr_err = fit_L0_L0_after_delta_eta_delta_phi_cuts_corr->GetParError(1)/(L0_alpha*L0_alpha);

    ResidualPolarization_delta_eta_delta_phi_hist[delta_eta_bin-1]->SetBinContent(2, L0_L0_slope_delta_eta_delta_phi_cuts_corr);
    ResidualPolarization_delta_eta_delta_phi_hist[delta_eta_bin-1]->SetBinError(2, L0_L0_slope_delta_eta_delta_phi_cuts_corr_err);


    TPaveText *L0_L0_text_MC_delta_eta_delta_phi_cuts_corr = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    L0_L0_text_MC_delta_eta_delta_phi_cuts_corr->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0_L0_text_MC_delta_eta_delta_phi_cuts_corr->AddText("Minimum bias");
    L0_L0_text_MC_delta_eta_delta_phi_cuts_corr->AddText("#Lambda-#Lambda");
    L0_L0_text_MC_delta_eta_delta_phi_cuts_corr->AddText("True MC");
    L0_L0_text_MC_delta_eta_delta_phi_cuts_corr->AddText("Analysis cuts");
    L0_L0_text_MC_delta_eta_delta_phi_cuts_corr->AddText(Form("P = %.3f #pm %.3f", L0_L0_slope_delta_eta_delta_phi_cuts_corr, L0_L0_slope_delta_eta_delta_phi_cuts_corr_err ));
    //if(corr_err == 1) L0_L0_text_MC_delta_eta_delta_phi_cuts_corr->AddText(Form("Slope sys. err: %.1f %%", L0_L0_slope_delta_eta_delta_phi_cuts_corr_sys_err*100 ));
    //if(corr_err == 0) L0_L0_text_MC_delta_eta_delta_phi_cuts_corr->AddText(Form("Slope sys. err: %.1f %%", L0_L0_slope_delta_eta_delta_phi_cuts_corr_sys_err_no_corr*100 ));
    L0_L0_text_MC_delta_eta_delta_phi_cuts_corr->SetFillColorAlpha(0, 0.01);
    L0_L0_text_MC_delta_eta_delta_phi_cuts_corr->Draw("same");


    L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_L/L0_L0_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------

    //Lbar-Lbar
    TCanvas *L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_can = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_can_%i", delta_eta_bin), Form("L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_can_%i", delta_eta_bin), 1200, 1000);
    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_can->cd();

    TF1 *fit_L0bar_L0bar_after_delta_eta_delta_phi_cuts = new TF1("fit_L0bar_L0bar_after_delta_eta_delta_phi_cuts", "[0]*(1+[1]*x)", -1, 1);
    fit_L0bar_L0bar_after_delta_eta_delta_phi_cuts->SetParameters(1000, 0.10);

    TF1 *fit_L0bar_L0bar_after_delta_eta_delta_phi_cuts_ME = new TF1("fit_L0bar_L0bar_after_delta_eta_delta_phi_cuts_ME", "[0]*(1+[1]*x)", -1, 1);
    fit_L0bar_L0bar_after_delta_eta_delta_phi_cuts_ME->SetParameters(1000, 0.10);

    TH1D *L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts = L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_cuts_hist->ProjectionX(Form("proj_eta_SE_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
    float nL0barL0bar_delta_eta_delta_phi_cuts = L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->Integral();

    TH1D *L0bar_L0bar_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts = L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_ME_cuts_hist->ProjectionX(Form("proj_eta_ME_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
    L0bar_L0bar_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts->Sumw2();
    L0bar_L0bar_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts->Scale(nL0barL0bar_delta_eta_delta_phi_cuts/L0bar_L0bar_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts->Integral());
    L0bar_L0bar_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts->Scale(1./L0bar_L0bar_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts->GetXaxis()->GetBinWidth(1));
    L0bar_L0bar_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts->Fit(fit_L0bar_L0bar_after_delta_eta_delta_phi_cuts_ME, "i 0 r");
    //L0bar_L0bar_cosThetaProdPlane_ME->Draw("p e");


    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->Sumw2();
    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->SetMarkerStyle(20);
    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->SetMarkerColor(kRed);
    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->SetLineColor(kRed);
    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->GetXaxis()->SetTitle("cos(#theta*)");
    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->GetXaxis()->CenterTitle();
    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->GetYaxis()->CenterTitle();
    //L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->Divide(L0bar_L0bar_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts); //correct using ME
    //L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->Scale(nL0barL0bar_delta_eta_delta_phi_cuts/L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->Integral()); //scale back
    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->Scale(1./L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->Fit(fit_L0bar_L0bar_after_delta_eta_delta_phi_cuts, "i 0 r");
    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->SetMinimum(0);
    //L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->GetYaxis()->SetRangeUser(-20e9, 0);
    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->Draw("p e");

    L0bar_L0bar_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts->Draw("p e same");

    fit_L0bar_L0bar_after_delta_eta_delta_phi_cuts->SetLineColor(1);
    fit_L0bar_L0bar_after_delta_eta_delta_phi_cuts->Draw("same");

    fit_L0bar_L0bar_after_delta_eta_delta_phi_cuts_ME->SetLineColor(kBlue);
    fit_L0bar_L0bar_after_delta_eta_delta_phi_cuts_ME->Draw("same");


    float L0bar_L0bar_slope_delta_eta_delta_phi_cuts = fit_L0bar_L0bar_after_delta_eta_delta_phi_cuts->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
    float L0bar_L0bar_slope_delta_eta_delta_phi_cuts_err = fit_L0bar_L0bar_after_delta_eta_delta_phi_cuts->GetParError(1)/(L0bar_alpha*L0bar_alpha);


    //statistical error correction for systematic error
    float SysErrSlope_L0bar_L0bar_delta_eta_delta_phi_corr = sqrt( fabs( fit_L0bar_L0bar_after_delta_eta_delta_phi_cuts->GetParError(1)*fit_L0bar_L0bar_after_delta_eta_delta_phi_cuts->GetParError(1) - fit_L0bar_L0bar_after_delta_eta_delta_phi_cuts_ME->GetParError(1)*fit_L0bar_L0bar_after_delta_eta_delta_phi_cuts_ME->GetParError(1) ) );

    float SysErrSlope_L0bar_L0bar_delta_eta_delta_phi_work = ( fabs( fit_L0bar_L0bar_after_delta_eta_delta_phi_cuts->GetParameter(1) - fit_L0bar_L0bar_after_delta_eta_delta_phi_cuts_ME->GetParameter(1) ) - SysErrSlope_L0bar_L0bar_delta_eta_delta_phi_corr )/fabs(fit_L0bar_L0bar_after_delta_eta_delta_phi_cuts->GetParameter(1));

    float L0bar_L0bar_slope_delta_eta_delta_phi_cuts_sys_err = 0;

    if( SysErrSlope_L0bar_L0bar_delta_eta_delta_phi_work > 0 ) L0bar_L0bar_slope_delta_eta_delta_phi_cuts_sys_err = SysErrSlope_L0bar_L0bar_delta_eta_delta_phi_work; //store sys. err. only if it is larger than statistical fluctuations

    float L0bar_L0bar_slope_delta_eta_delta_phi_cuts_sys_err_no_corr = fabs( fit_L0bar_L0bar_after_delta_eta_delta_phi_cuts->GetParameter(1) - fit_L0bar_L0bar_after_delta_eta_delta_phi_cuts_ME->GetParameter(1) )/fabs(fit_L0bar_L0bar_after_delta_eta_delta_phi_cuts->GetParameter(1));

    if(corr_err == 1) SysErrSlope_delta_eta_delta_phi_hist[delta_eta_bin-1]->SetBinContent(3, L0bar_L0bar_slope_delta_eta_delta_phi_cuts_sys_err);
    if(corr_err == 0) SysErrSlope_delta_eta_delta_phi_hist[delta_eta_bin-1]->SetBinContent(3, L0bar_L0bar_slope_delta_eta_delta_phi_cuts_sys_err_no_corr);

    TPaveText *L0bar_L0bar_text_MC_delta_eta_delta_phi_cuts = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    L0bar_L0bar_text_MC_delta_eta_delta_phi_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0bar_L0bar_text_MC_delta_eta_delta_phi_cuts->AddText("Minimum bias");
    L0bar_L0bar_text_MC_delta_eta_delta_phi_cuts->AddText("#bar{#Lambda}-#bar{#Lambda}");
    L0bar_L0bar_text_MC_delta_eta_delta_phi_cuts->AddText("True MC");
    L0bar_L0bar_text_MC_delta_eta_delta_phi_cuts->AddText("Analysis cuts");
    L0bar_L0bar_text_MC_delta_eta_delta_phi_cuts->AddText(Form("P = %.3f #pm %.3f", L0bar_L0bar_slope_delta_eta_delta_phi_cuts, L0bar_L0bar_slope_delta_eta_delta_phi_cuts_err ));
    L0bar_L0bar_text_MC_delta_eta_delta_phi_cuts->AddText(Form("P_{ME} = %.3f", fit_L0bar_L0bar_after_delta_eta_delta_phi_cuts_ME->GetParameter(1)/(L0bar_alpha*L0bar_alpha) ));
    if(corr_err == 1) L0bar_L0bar_text_MC_delta_eta_delta_phi_cuts->AddText(Form("Slope sys. err: %.1f %%", L0bar_L0bar_slope_delta_eta_delta_phi_cuts_sys_err*100 ));
    if(corr_err == 0) L0bar_L0bar_text_MC_delta_eta_delta_phi_cuts->AddText(Form("Slope sys. err: %.1f %%", L0bar_L0bar_slope_delta_eta_delta_phi_cuts_sys_err_no_corr*100 ));
    L0bar_L0bar_text_MC_delta_eta_delta_phi_cuts->SetFillColorAlpha(0, 0.01);
    L0bar_L0bar_text_MC_delta_eta_delta_phi_cuts->Draw("same");

    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/Lbar_Lbar/L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------

    TCanvas *L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr_can = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr_can_%i", delta_eta_bin), Form("L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr_can_%i", delta_eta_bin), 1200, 1000);
    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr_can->cd();

    TF1 *fit_L0bar_L0bar_after_delta_eta_delta_phi_cuts_corr = new TF1("fit_L0bar_L0bar_after_delta_eta_delta_phi_cuts_corr", "[0]*(1+[1]*x)", -1, 1);
    fit_L0bar_L0bar_after_delta_eta_delta_phi_cuts_corr->SetParameters(1000, 0.10);

    TH1D *L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr = (TH1D*)L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts->Clone(Form("L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr_%i", delta_eta_bin));
    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr->Divide(L0bar_L0bar_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts); //correct
    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr->Scale(nL0barL0bar_delta_eta_delta_phi_cuts/L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr->Integral()); //scale back
    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr->Scale(1./L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr->GetXaxis()->GetBinWidth(1)); //divide by bin width
    //L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr->Scale(nL0barL0bar_delta_eta_delta_phi_cuts);
    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr->Fit(fit_L0bar_L0bar_after_delta_eta_delta_phi_cuts_corr, "i 0 r");
    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr->SetMinimum(0);
    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr->Draw("p e");

    fit_L0bar_L0bar_after_delta_eta_delta_phi_cuts_corr->SetLineColor(1);
    fit_L0bar_L0bar_after_delta_eta_delta_phi_cuts_corr->Draw("same");

    float L0bar_L0bar_slope_delta_eta_delta_phi_cuts_corr = fit_L0bar_L0bar_after_delta_eta_delta_phi_cuts_corr->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
    float L0bar_L0bar_slope_delta_eta_delta_phi_cuts_corr_err = fit_L0bar_L0bar_after_delta_eta_delta_phi_cuts_corr->GetParError(1)/(L0bar_alpha*L0bar_alpha);

    ResidualPolarization_delta_eta_delta_phi_hist[delta_eta_bin-1]->SetBinContent(3, L0bar_L0bar_slope_delta_eta_delta_phi_cuts_corr);
    ResidualPolarization_delta_eta_delta_phi_hist[delta_eta_bin-1]->SetBinError(3, L0bar_L0bar_slope_delta_eta_delta_phi_cuts_corr_err);


    TPaveText *L0bar_L0bar_text_MC_delta_eta_delta_phi_cuts_corr = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    L0bar_L0bar_text_MC_delta_eta_delta_phi_cuts_corr->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0bar_L0bar_text_MC_delta_eta_delta_phi_cuts_corr->AddText("Minimum bias");
    L0bar_L0bar_text_MC_delta_eta_delta_phi_cuts_corr->AddText("#bar{#Lambda}-#bar{#Lambda}");
    L0bar_L0bar_text_MC_delta_eta_delta_phi_cuts_corr->AddText("True MC");
    L0bar_L0bar_text_MC_delta_eta_delta_phi_cuts_corr->AddText("Analysis cuts");
    L0bar_L0bar_text_MC_delta_eta_delta_phi_cuts_corr->AddText(Form("P = %.3f #pm %.3f", L0bar_L0bar_slope_delta_eta_delta_phi_cuts_corr, L0bar_L0bar_slope_delta_eta_delta_phi_cuts_corr_err ));
    //if(corr_err == 1) L0bar_L0bar_text_MC_delta_eta_delta_phi_cuts_corr->AddText(Form("Slope sys. err: %.1f %%", L0bar_L0bar_slope_delta_eta_delta_phi_cuts_corr_sys_err*100 ));
    //if(corr_err == 0) L0bar_L0bar_text_MC_delta_eta_delta_phi_cuts_corr->AddText(Form("Slope sys. err: %.1f %%", L0bar_L0bar_slope_delta_eta_delta_phi_cuts_corr_sys_err_no_corr*100 ));
    L0bar_L0bar_text_MC_delta_eta_delta_phi_cuts_corr->SetFillColorAlpha(0, 0.01);
    L0bar_L0bar_text_MC_delta_eta_delta_phi_cuts_corr->Draw("same");


    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/Lbar_Lbar/L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------

    SysErrSlope_delta_eta_delta_phi_hist[delta_eta_bin-1]->Write();
    ResidualPolarization_delta_eta_delta_phi_hist[delta_eta_bin-1]->Write();

  }


  //_________________________________________________________________________________________________________________________________________________________________________

  //with Delta phi

  //for( unsigned int delta_eta_bin = 1; delta_eta_bin < 4; delta_eta_bin++)
  for( unsigned int delta_eta_bin = 0; delta_eta_bin < 3; delta_eta_bin++)
  {
    //after cuts
    TCanvas *L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts_can_%i", delta_eta_bin), Form("L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts_can_%i", delta_eta_bin), 1200, 1000);
    L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts_can->cd();

    TF1 *fit_L0_L0bar_after_less_delta_phi_cuts = new TF1("fit_L0_L0bar_after_less_delta_phi_cuts", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0bar_after_less_delta_phi_cuts->SetParameters(1000, 0.10);

    TF1 *fit_L0_L0bar_after_less_delta_phi_cuts_ME = new TF1("fit_L0_L0bar_after_less_delta_phi_cuts_ME", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0bar_after_less_delta_phi_cuts_ME->SetParameters(1000, 0.10);

    //TH1D *L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts = L0_L0bar_cos_theta_star_vs_delta_phi_cuts_hist->ProjectionX(Form("proj_SE_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
    TH1D *L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts = L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_cuts_hist->ProjectionX(Form("proj_SE_%i", delta_eta_bin), 1, 5*delta_eta_bin+5);
    //TH1D *L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts = L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_cuts_hist->ProjectionX(Form("proj_SE_%i", delta_eta_bin), 1, delta_eta_bin+1);
    float nL0L0bar_less_delta_phi_cuts = L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts->Integral();

    //TH1D *L0_L0bar_cosThetaProdPlane_ME_less_delta_phi_cuts = L0_L0bar_cos_theta_star_vs_delta_phi_ME_cuts_hist->ProjectionX(Form("proj_ME_%i", delta_eta_bin), delta_eta_bin+1, delta_eta_bin);
    TH1D *L0_L0bar_cosThetaProdPlane_ME_less_delta_phi_cuts = L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_ME_cuts_hist->ProjectionX(Form("proj_ME_%i", delta_eta_bin), 1, 5*delta_eta_bin+5);
    //TH1D *L0_L0bar_cosThetaProdPlane_ME_less_delta_phi_cuts = L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_ME_cuts_hist->ProjectionX(Form("proj_ME_%i", delta_eta_bin), 1, delta_eta_bin+1);
    L0_L0bar_cosThetaProdPlane_ME_less_delta_phi_cuts->Sumw2();
    L0_L0bar_cosThetaProdPlane_ME_less_delta_phi_cuts->Scale(nL0L0bar_less_delta_phi_cuts/L0_L0bar_cosThetaProdPlane_ME_less_delta_phi_cuts->Integral());
    L0_L0bar_cosThetaProdPlane_ME_less_delta_phi_cuts->Scale(1./L0_L0bar_cosThetaProdPlane_ME_less_delta_phi_cuts->GetXaxis()->GetBinWidth(1));
    L0_L0bar_cosThetaProdPlane_ME_less_delta_phi_cuts->Fit(fit_L0_L0bar_after_less_delta_phi_cuts_ME, "i 0 r");
    //L0_L0bar_cosThetaProdPlane_ME->Draw("p e");


    L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts->Sumw2();
    L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts->SetMarkerStyle(20);
    L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts->SetMarkerColor(kRed);
    L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts->SetLineColor(kRed);
    L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts->GetXaxis()->SetTitle("cos(#theta*)");
    L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts->GetXaxis()->CenterTitle();
    L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts->GetYaxis()->CenterTitle();
    //L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts->Divide(L0_L0bar_cosThetaProdPlane_ME_less_delta_phi_cuts); //correct using ME
    //L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts->Scale(nL0L0bar_less_delta_phi_cuts/L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts->Integral()); //scale back
    L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts->Scale(1./L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
    L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts->Fit(fit_L0_L0bar_after_less_delta_phi_cuts, "i 0 r");
    L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts->SetMinimum(0);
    //L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts->GetYaxis()->SetRangeUser(-20e9, 0);
    L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts->Draw("p e");

    L0_L0bar_cosThetaProdPlane_ME_less_delta_phi_cuts->Draw("p e same");

    fit_L0_L0bar_after_less_delta_phi_cuts->SetLineColor(1);
    fit_L0_L0bar_after_less_delta_phi_cuts->Draw("same");

    fit_L0_L0bar_after_less_delta_phi_cuts_ME->SetLineColor(kBlue);
    fit_L0_L0bar_after_less_delta_phi_cuts_ME->Draw("same");


    float L0_L0bar_slope_less_delta_phi_cuts = fit_L0_L0bar_after_less_delta_phi_cuts->GetParameter(1)/(L0_alpha*L0bar_alpha);
    float L0_L0bar_slope_less_delta_phi_cuts_err = fit_L0_L0bar_after_less_delta_phi_cuts->GetParError(1)/(L0_alpha*L0bar_alpha);

    //statistical error correction for systematic error
    float SysErrSlope_L0_L0bar_less_delta_phi_corr = sqrt( fabs( fit_L0_L0bar_after_less_delta_phi_cuts->GetParError(1)*fit_L0_L0bar_after_less_delta_phi_cuts->GetParError(1) - fit_L0_L0bar_after_less_delta_phi_cuts_ME->GetParError(1)*fit_L0_L0bar_after_less_delta_phi_cuts_ME->GetParError(1) ) );

    float SysErrSlope_L0_L0bar_less_delta_phi_work = ( fabs( fit_L0_L0bar_after_less_delta_phi_cuts->GetParameter(1) - fit_L0_L0bar_after_less_delta_phi_cuts_ME->GetParameter(1) ) - SysErrSlope_L0_L0bar_less_delta_phi_corr )/fabs(fit_L0_L0bar_after_less_delta_phi_cuts->GetParameter(1));

    float L0_L0bar_slope_less_delta_phi_cuts_sys_err = 0;

    if( SysErrSlope_L0_L0bar_less_delta_phi_work > 0 ) L0_L0bar_slope_less_delta_phi_cuts_sys_err = SysErrSlope_L0_L0bar_less_delta_phi_work; //store sys. err. only if it is larger than statistical fluctuations

    float L0_L0bar_slope_less_delta_phi_cuts_sys_err_no_corr = fabs( fit_L0_L0bar_after_less_delta_phi_cuts->GetParameter(1) - fit_L0_L0bar_after_less_delta_phi_cuts_ME->GetParameter(1) )/fabs(fit_L0_L0bar_after_less_delta_phi_cuts->GetParameter(1));

    if(corr_err == 1) SysErrSlope_delta_phi_less_hist[delta_eta_bin]->SetBinContent(1, L0_L0bar_slope_less_delta_phi_cuts_sys_err);
    if(corr_err == 0) SysErrSlope_delta_phi_less_hist[delta_eta_bin]->SetBinContent(1, L0_L0bar_slope_less_delta_phi_cuts_sys_err_no_corr);

    TPaveText *L0_L0bar_text_MC_less_delta_phi_cuts = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    L0_L0bar_text_MC_less_delta_phi_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0_L0bar_text_MC_less_delta_phi_cuts->AddText("Minimum bias");
    L0_L0bar_text_MC_less_delta_phi_cuts->AddText("#Lambda-#bar{#Lambda}");
    L0_L0bar_text_MC_less_delta_phi_cuts->AddText("True MC");
    L0_L0bar_text_MC_less_delta_phi_cuts->AddText("Analysis cuts");
    L0_L0bar_text_MC_less_delta_phi_cuts->AddText(Form("P = %.3f #pm %.3f", L0_L0bar_slope_less_delta_phi_cuts, L0_L0bar_slope_less_delta_phi_cuts_err ));
    L0_L0bar_text_MC_less_delta_phi_cuts->AddText(Form("P_{ME} = %.3f", fit_L0_L0bar_after_less_delta_phi_cuts_ME->GetParameter(1)/(L0_alpha*L0bar_alpha) ));
    if(corr_err == 1) L0_L0bar_text_MC_less_delta_phi_cuts->AddText(Form("Slope sys. err: %.1f %%", L0_L0bar_slope_less_delta_phi_cuts_sys_err*100 ));
    if(corr_err == 0) L0_L0bar_text_MC_less_delta_phi_cuts->AddText(Form("Slope sys. err: %.1f %%", L0_L0bar_slope_less_delta_phi_cuts_sys_err_no_corr*100 ));
    L0_L0bar_text_MC_less_delta_phi_cuts->SetFillColorAlpha(0, 0.01);
    L0_L0bar_text_MC_less_delta_phi_cuts->Draw("same");

    L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------

    TCanvas *L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts_corr_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts_corr_can_%i", delta_eta_bin), Form("L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts_corr_can_%i", delta_eta_bin), 1200, 1000);
    L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts_corr_can->cd();

    TF1 *fit_L0_L0bar_after_less_delta_phi_cuts_corr = new TF1("fit_L0_L0bar_after_less_delta_phi_cuts_corr", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0bar_after_less_delta_phi_cuts_corr->SetParameters(1000, 0.10);

    TH1D *L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts_corr = (TH1D*)L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts->Clone(Form("L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts_corr_%i", delta_eta_bin));
    L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts_corr->Divide(L0_L0bar_cosThetaProdPlane_ME_less_delta_phi_cuts); //correct
    L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts_corr->Scale(nL0L0bar_less_delta_phi_cuts/L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts_corr->Integral()); //scale back
    L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts_corr->Scale(1./L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts_corr->GetXaxis()->GetBinWidth(1)); //divide by bin width
    //L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts_corr->Scale(nL0L0bar_less_delta_phi_cuts);
    L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts_corr->Fit(fit_L0_L0bar_after_less_delta_phi_cuts_corr, "i 0 r");
    L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts_corr->SetMinimum(0);
    L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts_corr->Draw("p e");

    fit_L0_L0bar_after_less_delta_phi_cuts_corr->SetLineColor(1);
    fit_L0_L0bar_after_less_delta_phi_cuts_corr->Draw("same");

    float L0_L0bar_slope_less_delta_phi_cuts_corr = fit_L0_L0bar_after_less_delta_phi_cuts_corr->GetParameter(1)/(L0_alpha*L0bar_alpha);
    float L0_L0bar_slope_less_delta_phi_cuts_corr_err = fit_L0_L0bar_after_less_delta_phi_cuts_corr->GetParError(1)/(L0_alpha*L0bar_alpha);

    ResidualPolarization_delta_phi_less_hist[delta_eta_bin]->SetBinContent(1, L0_L0bar_slope_less_delta_phi_cuts_corr);
    ResidualPolarization_delta_phi_less_hist[delta_eta_bin]->SetBinError(1, L0_L0bar_slope_less_delta_phi_cuts_corr_err);


    TPaveText *L0_L0bar_text_MC_less_delta_phi_cuts_corr = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    L0_L0bar_text_MC_less_delta_phi_cuts_corr->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0_L0bar_text_MC_less_delta_phi_cuts_corr->AddText("Minimum bias");
    L0_L0bar_text_MC_less_delta_phi_cuts_corr->AddText("#Lambda-#bar{#Lambda}");
    L0_L0bar_text_MC_less_delta_phi_cuts_corr->AddText("True MC");
    L0_L0bar_text_MC_less_delta_phi_cuts_corr->AddText("Analysis cuts");
    L0_L0bar_text_MC_less_delta_phi_cuts_corr->AddText(Form("P = %.3f #pm %.3f", L0_L0bar_slope_less_delta_phi_cuts_corr, L0_L0bar_slope_less_delta_phi_cuts_corr_err ));
    //if(corr_err == 1) L0_L0bar_text_MC_less_delta_phi_cuts_corr->AddText(Form("Slope sys. err: %.1f %%", L0_L0bar_slope_less_delta_phi_cuts_corr_sys_err*100 ));
    //if(corr_err == 0) L0_L0bar_text_MC_less_delta_phi_cuts_corr->AddText(Form("Slope sys. err: %.1f %%", L0_L0bar_slope_less_delta_phi_cuts_corr_sys_err_no_corr*100 ));
    L0_L0bar_text_MC_less_delta_phi_cuts_corr->SetFillColorAlpha(0, 0.01);
    L0_L0bar_text_MC_less_delta_phi_cuts_corr->Draw("same");


    L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/L0_L0bar_cosThetaProdPlane_less_delta_phi_cuts_corr_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------

    TCanvas *L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts_can_%i", delta_eta_bin), Form("L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts_can_%i", delta_eta_bin), 1200, 1000);
    L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts_can->cd();

    TF1 *fit_L0_L0bar_after_more_delta_phi_cuts = new TF1("fit_L0_L0bar_after_more_delta_phi_cuts", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0bar_after_more_delta_phi_cuts->SetParameters(1000, 0.10);

    TF1 *fit_L0_L0bar_after_more_delta_phi_cuts_ME = new TF1("fit_L0_L0bar_after_more_delta_phi_cuts_ME", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0bar_after_more_delta_phi_cuts_ME->SetParameters(1000, 0.10);


    //TH1D *L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts = L0_L0bar_cos_theta_star_vs_delta_phi_cuts_2_hist->ProjectionX(Form("proj_SE_2_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
    TH1D *L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts = L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_cuts_hist->ProjectionX(Form("proj_SE_2_%i", delta_eta_bin), 5*delta_eta_bin+1, 20);
    float nL0L0bar_more_delta_phi_cuts = L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts->Integral();

    //TH1D *L0_L0bar_cosThetaProdPlane_ME_more_delta_phi_cuts = L0_L0bar_cos_theta_star_vs_delta_phi_ME_cuts_2_hist->ProjectionX(Form("proj_ME_2_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
    TH1D *L0_L0bar_cosThetaProdPlane_ME_more_delta_phi_cuts = L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_ME_cuts_hist->ProjectionX(Form("proj_ME_2_%i", delta_eta_bin), 5*delta_eta_bin+1, 20);
    L0_L0bar_cosThetaProdPlane_ME_more_delta_phi_cuts->Sumw2();
    L0_L0bar_cosThetaProdPlane_ME_more_delta_phi_cuts->Scale(nL0L0bar_more_delta_phi_cuts/L0_L0bar_cosThetaProdPlane_ME_more_delta_phi_cuts->Integral());
    L0_L0bar_cosThetaProdPlane_ME_more_delta_phi_cuts->Scale(1./L0_L0bar_cosThetaProdPlane_ME_more_delta_phi_cuts->GetXaxis()->GetBinWidth(1));
    L0_L0bar_cosThetaProdPlane_ME_more_delta_phi_cuts->Fit(fit_L0_L0bar_after_more_delta_phi_cuts_ME, "i 0 r");
    //L0_L0bar_cosThetaProdPlane_ME->Draw("p e");


    L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts->Sumw2();
    L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts->SetMarkerStyle(20);
    L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts->SetMarkerColor(kRed);
    L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts->SetLineColor(kRed);
    L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts->GetXaxis()->SetTitle("cos(#theta*)");
    L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts->GetXaxis()->CenterTitle();
    L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts->GetYaxis()->CenterTitle();
    //L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts->Divide(L0_L0bar_cosThetaProdPlane_ME_more_delta_phi_cuts); //correct using ME
    //L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts->Scale(nL0L0bar_more_delta_phi_cuts/L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts->Integral()); //scale back
    L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts->Scale(1./L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
    L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts->Fit(fit_L0_L0bar_after_more_delta_phi_cuts, "i 0 r");
    L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts->SetMinimum(0);
    //L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts->GetYaxis()->SetRangeUser(-20e9, 0);
    L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts->Draw("p e");

    L0_L0bar_cosThetaProdPlane_ME_more_delta_phi_cuts->Draw("p e same");

    fit_L0_L0bar_after_more_delta_phi_cuts->SetLineColor(1);
    fit_L0_L0bar_after_more_delta_phi_cuts->Draw("same");

    fit_L0_L0bar_after_more_delta_phi_cuts_ME->SetLineColor(kBlue);
    fit_L0_L0bar_after_more_delta_phi_cuts_ME->Draw("same");


    float L0_L0bar_slope_more_delta_phi_cuts = fit_L0_L0bar_after_more_delta_phi_cuts->GetParameter(1)/(L0_alpha*L0bar_alpha);
    float L0_L0bar_slope_more_delta_phi_cuts_err = fit_L0_L0bar_after_more_delta_phi_cuts->GetParError(1)/(L0_alpha*L0bar_alpha);

    //statistical error correction for systematic error
    float SysErrSlope_L0_L0bar_more_delta_phi_corr = sqrt( fabs( fit_L0_L0bar_after_more_delta_phi_cuts->GetParError(1)*fit_L0_L0bar_after_more_delta_phi_cuts->GetParError(1) - fit_L0_L0bar_after_more_delta_phi_cuts_ME->GetParError(1)*fit_L0_L0bar_after_more_delta_phi_cuts_ME->GetParError(1) ) );

    float SysErrSlope_L0_L0bar_more_delta_phi_work = ( fabs( fit_L0_L0bar_after_more_delta_phi_cuts->GetParameter(1) - fit_L0_L0bar_after_more_delta_phi_cuts_ME->GetParameter(1) ) - SysErrSlope_L0_L0bar_more_delta_phi_corr )/fabs(fit_L0_L0bar_after_more_delta_phi_cuts->GetParameter(1));

    float L0_L0bar_slope_more_delta_phi_cuts_sys_err = 0;

    if( SysErrSlope_L0_L0bar_more_delta_phi_work > 0 ) L0_L0bar_slope_more_delta_phi_cuts_sys_err = SysErrSlope_L0_L0bar_more_delta_phi_work; //store sys. err. only if it is larger than statistical fluctuations

    float L0_L0bar_slope_more_delta_phi_cuts_sys_err_no_corr = fabs( fit_L0_L0bar_after_more_delta_phi_cuts->GetParameter(1) - fit_L0_L0bar_after_more_delta_phi_cuts_ME->GetParameter(1) )/fabs(fit_L0_L0bar_after_more_delta_phi_cuts->GetParameter(1));

    if(corr_err == 1) SysErrSlope_delta_phi_more_hist[delta_eta_bin]->SetBinContent(1, L0_L0bar_slope_more_delta_phi_cuts_sys_err);
    if(corr_err == 0) SysErrSlope_delta_phi_more_hist[delta_eta_bin]->SetBinContent(1, L0_L0bar_slope_more_delta_phi_cuts_sys_err_no_corr);

    TPaveText *L0_L0bar_text_MC_more_delta_phi_cuts = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    L0_L0bar_text_MC_more_delta_phi_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0_L0bar_text_MC_more_delta_phi_cuts->AddText("Minimum bias");
    L0_L0bar_text_MC_more_delta_phi_cuts->AddText("#Lambda-#bar{#Lambda}");
    L0_L0bar_text_MC_more_delta_phi_cuts->AddText("True MC");
    L0_L0bar_text_MC_more_delta_phi_cuts->AddText("Analysis cuts");
    L0_L0bar_text_MC_more_delta_phi_cuts->AddText(Form("P = %.3f #pm %.3f", L0_L0bar_slope_more_delta_phi_cuts, L0_L0bar_slope_more_delta_phi_cuts_err ));
    L0_L0bar_text_MC_more_delta_phi_cuts->AddText(Form("P_{ME} = %.3f", fit_L0_L0bar_after_more_delta_phi_cuts_ME->GetParameter(1)/(L0_alpha*L0bar_alpha) ));
    if(corr_err == 1) L0_L0bar_text_MC_more_delta_phi_cuts->AddText(Form("Slope sys. err: %.1f %%", L0_L0bar_slope_more_delta_phi_cuts_sys_err*100 ));
    if(corr_err == 0) L0_L0bar_text_MC_more_delta_phi_cuts->AddText(Form("Slope sys. err: %.1f %%", L0_L0bar_slope_more_delta_phi_cuts_sys_err_no_corr*100 ));
    L0_L0bar_text_MC_more_delta_phi_cuts->SetFillColorAlpha(0, 0.01);
    L0_L0bar_text_MC_more_delta_phi_cuts->Draw("same");

    L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts_%i.png", delta_eta_bin));

    //--------------------------------------------------------------------------------------------

    TCanvas *L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts_corr_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts_corr_can_%i", delta_eta_bin), Form("L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts_corr_can_%i", delta_eta_bin), 1200, 1000);
    L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts_corr_can->cd();

    TF1 *fit_L0_L0bar_after_more_delta_phi_cuts_corr = new TF1("fit_L0_L0bar_after_more_delta_phi_cuts_corr", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0bar_after_more_delta_phi_cuts_corr->SetParameters(1000, 0.10);

    TH1D *L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts_corr = (TH1D*)L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts->Clone(Form("L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts_corr_%i", delta_eta_bin));
    L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts_corr->Divide(L0_L0bar_cosThetaProdPlane_ME_more_delta_phi_cuts); //correct
    L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts_corr->Scale(nL0L0bar_more_delta_phi_cuts/L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts_corr->Integral()); //scale back
    L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts_corr->Scale(1./L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts_corr->GetXaxis()->GetBinWidth(1)); //divide by bin width
    //L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts_corr->Scale(nL0L0bar_more_delta_phi_cuts);
    L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts_corr->Fit(fit_L0_L0bar_after_more_delta_phi_cuts_corr, "i 0 r");
    L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts_corr->SetMinimum(0);
    L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts_corr->Draw("p e");

    fit_L0_L0bar_after_more_delta_phi_cuts_corr->SetLineColor(1);
    fit_L0_L0bar_after_more_delta_phi_cuts_corr->Draw("same");

    float L0_L0bar_slope_more_delta_phi_cuts_corr = fit_L0_L0bar_after_more_delta_phi_cuts_corr->GetParameter(1)/(L0_alpha*L0bar_alpha);
    float L0_L0bar_slope_more_delta_phi_cuts_corr_err = fit_L0_L0bar_after_more_delta_phi_cuts_corr->GetParError(1)/(L0_alpha*L0bar_alpha);

    ResidualPolarization_delta_phi_more_hist[delta_eta_bin]->SetBinContent(1, L0_L0bar_slope_more_delta_phi_cuts_corr);
    ResidualPolarization_delta_phi_more_hist[delta_eta_bin]->SetBinError(1, L0_L0bar_slope_more_delta_phi_cuts_corr_err);


    TPaveText *L0_L0bar_text_MC_more_delta_phi_cuts_corr = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    L0_L0bar_text_MC_more_delta_phi_cuts_corr->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0_L0bar_text_MC_more_delta_phi_cuts_corr->AddText("Minimum bias");
    L0_L0bar_text_MC_more_delta_phi_cuts_corr->AddText("#Lambda-#bar{#Lambda}");
    L0_L0bar_text_MC_more_delta_phi_cuts_corr->AddText("True MC");
    L0_L0bar_text_MC_more_delta_phi_cuts_corr->AddText("Analysis cuts");
    L0_L0bar_text_MC_more_delta_phi_cuts_corr->AddText(Form("P = %.3f #pm %.3f", L0_L0bar_slope_more_delta_phi_cuts_corr, L0_L0bar_slope_more_delta_phi_cuts_corr_err ));
    //if(corr_err == 1) L0_L0bar_text_MC_more_delta_phi_cuts_corr->AddText(Form("Slope sys. err: %.1f %%", L0_L0bar_slope_more_delta_phi_cuts_corr_sys_err*100 ));
    //if(corr_err == 0) L0_L0bar_text_MC_more_delta_phi_cuts_corr->AddText(Form("Slope sys. err: %.1f %%", L0_L0bar_slope_more_delta_phi_cuts_corr_sys_err_no_corr*100 ));
    L0_L0bar_text_MC_more_delta_phi_cuts_corr->SetFillColorAlpha(0, 0.01);
    L0_L0bar_text_MC_more_delta_phi_cuts_corr->Draw("same");


    L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/L0_L0bar_cosThetaProdPlane_more_delta_phi_cuts_corr_%i.png", delta_eta_bin));

    //________________________________________________________________________________________________________________________________________

    TCanvas *L0_L0_cosThetaProdPlane_less_delta_phi_cuts_can = new TCanvas(Form("L0_L0_cosThetaProdPlane_less_delta_phi_cuts_can_%i", delta_eta_bin), Form("L0_L0_cosThetaProdPlane_less_delta_phi_cuts_can_%i", delta_eta_bin), 1200, 1000);
    L0_L0_cosThetaProdPlane_less_delta_phi_cuts_can->cd();

    TF1 *fit_L0_L0_after_less_delta_phi_cuts = new TF1("fit_L0_L0_after_less_delta_phi_cuts", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0_after_less_delta_phi_cuts->SetParameters(1000, 0.10);

    TF1 *fit_L0_L0_after_less_delta_phi_cuts_ME = new TF1("fit_L0_L0_after_less_delta_phi_cuts_ME", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0_after_less_delta_phi_cuts_ME->SetParameters(1000, 0.10);

    //TH1D *L0_L0_cosThetaProdPlane_less_delta_phi_cuts = L0_L0_cos_theta_star_vs_delta_phi_cuts_hist->ProjectionX(Form("proj_SE_L_L_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
    TH1D *L0_L0_cosThetaProdPlane_less_delta_phi_cuts = L0_L0_cos_theta_star_vs_delta_phi_for_corr_cuts_hist->ProjectionX(Form("proj_SE_L_L_%i", delta_eta_bin), 1, 5*delta_eta_bin+5);
    float nL0L0_less_delta_phi_cuts = L0_L0_cosThetaProdPlane_less_delta_phi_cuts->Integral();

    //TH1D *L0_L0_cosThetaProdPlane_ME_less_delta_phi_cuts = L0_L0_cos_theta_star_vs_delta_phi_ME_cuts_hist->ProjectionX(Form("proj_ME_L_L_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
    TH1D *L0_L0_cosThetaProdPlane_ME_less_delta_phi_cuts = L0_L0_cos_theta_star_vs_delta_phi_for_corr_ME_cuts_hist->ProjectionX(Form("proj_ME_L_L_%i", delta_eta_bin), 1, 5*delta_eta_bin+5);
    L0_L0_cosThetaProdPlane_ME_less_delta_phi_cuts->Sumw2();
    L0_L0_cosThetaProdPlane_ME_less_delta_phi_cuts->Scale(nL0L0_less_delta_phi_cuts/L0_L0_cosThetaProdPlane_ME_less_delta_phi_cuts->Integral());
    L0_L0_cosThetaProdPlane_ME_less_delta_phi_cuts->Scale(1./L0_L0_cosThetaProdPlane_ME_less_delta_phi_cuts->GetXaxis()->GetBinWidth(1));
    L0_L0_cosThetaProdPlane_ME_less_delta_phi_cuts->Fit(fit_L0_L0_after_less_delta_phi_cuts_ME, "i 0 r");
    //L0_L0_cosThetaProdPlane_ME->Draw("p e");


    L0_L0_cosThetaProdPlane_less_delta_phi_cuts->Sumw2();
    L0_L0_cosThetaProdPlane_less_delta_phi_cuts->SetMarkerStyle(20);
    L0_L0_cosThetaProdPlane_less_delta_phi_cuts->SetMarkerColor(kRed);
    L0_L0_cosThetaProdPlane_less_delta_phi_cuts->SetLineColor(kRed);
    L0_L0_cosThetaProdPlane_less_delta_phi_cuts->GetXaxis()->SetTitle("cos(#theta*)");
    L0_L0_cosThetaProdPlane_less_delta_phi_cuts->GetXaxis()->CenterTitle();
    L0_L0_cosThetaProdPlane_less_delta_phi_cuts->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    L0_L0_cosThetaProdPlane_less_delta_phi_cuts->GetYaxis()->CenterTitle();
    //L0_L0_cosThetaProdPlane_less_delta_phi_cuts->Divide(L0_L0_cosThetaProdPlane_ME_less_delta_phi_cuts); //correct using ME
    //L0_L0_cosThetaProdPlane_less_delta_phi_cuts->Scale(nL0L0_less_delta_phi_cuts/L0_L0_cosThetaProdPlane_less_delta_phi_cuts->Integral()); //scale back
    L0_L0_cosThetaProdPlane_less_delta_phi_cuts->Scale(1./L0_L0_cosThetaProdPlane_less_delta_phi_cuts->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
    L0_L0_cosThetaProdPlane_less_delta_phi_cuts->Fit(fit_L0_L0_after_less_delta_phi_cuts, "i 0 r");
    L0_L0_cosThetaProdPlane_less_delta_phi_cuts->SetMinimum(0);
    //L0_L0_cosThetaProdPlane_less_delta_phi_cuts->GetYaxis()->SetRangeUser(-20e9, 0);
    L0_L0_cosThetaProdPlane_less_delta_phi_cuts->Draw("p e");

    L0_L0_cosThetaProdPlane_ME_less_delta_phi_cuts->Draw("p e same");

    fit_L0_L0_after_less_delta_phi_cuts->SetLineColor(1);
    fit_L0_L0_after_less_delta_phi_cuts->Draw("same");

    fit_L0_L0_after_less_delta_phi_cuts_ME->SetLineColor(kBlue);
    fit_L0_L0_after_less_delta_phi_cuts_ME->Draw("same");


    float L0_L0_slope_less_delta_phi_cuts = fit_L0_L0_after_less_delta_phi_cuts->GetParameter(1)/(L0_alpha*L0_alpha);
    float L0_L0_slope_less_delta_phi_cuts_err = fit_L0_L0_after_less_delta_phi_cuts->GetParError(1)/(L0_alpha*L0_alpha);

    //statistical error correction for systematic error
    float SysErrSlope_L0_L0_less_delta_phi_corr = sqrt( fabs( fit_L0_L0_after_less_delta_phi_cuts->GetParError(1)*fit_L0_L0_after_less_delta_phi_cuts->GetParError(1) - fit_L0_L0_after_less_delta_phi_cuts_ME->GetParError(1)*fit_L0_L0_after_less_delta_phi_cuts_ME->GetParError(1) ) );

    float SysErrSlope_L0_L0_less_delta_phi_work = ( fabs( fit_L0_L0_after_less_delta_phi_cuts->GetParameter(1) - fit_L0_L0_after_less_delta_phi_cuts_ME->GetParameter(1) ) - SysErrSlope_L0_L0_less_delta_phi_corr )/fabs(fit_L0_L0_after_less_delta_phi_cuts->GetParameter(1));

    float L0_L0_slope_less_delta_phi_cuts_sys_err = 0;

    if( SysErrSlope_L0_L0_less_delta_phi_work > 0 ) L0_L0_slope_less_delta_phi_cuts_sys_err = SysErrSlope_L0_L0_less_delta_phi_work; //store sys. err. only if it is larger than statistical fluctuations

    float L0_L0_slope_less_delta_phi_cuts_sys_err_no_corr = fabs( fit_L0_L0_after_less_delta_phi_cuts->GetParameter(1) - fit_L0_L0_after_less_delta_phi_cuts_ME->GetParameter(1) )/fabs(fit_L0_L0_after_less_delta_phi_cuts->GetParameter(1));

    if(corr_err == 1) SysErrSlope_delta_phi_less_hist[delta_eta_bin]->SetBinContent(2, L0_L0_slope_less_delta_phi_cuts_sys_err);
    if(corr_err == 0) SysErrSlope_delta_phi_less_hist[delta_eta_bin]->SetBinContent(2, L0_L0_slope_less_delta_phi_cuts_sys_err_no_corr);

    TPaveText *L0_L0_text_MC_less_delta_phi_cuts = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    L0_L0_text_MC_less_delta_phi_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0_L0_text_MC_less_delta_phi_cuts->AddText("Minimum bias");
    L0_L0_text_MC_less_delta_phi_cuts->AddText("#Lambda-#Lambda");
    L0_L0_text_MC_less_delta_phi_cuts->AddText("True MC");
    L0_L0_text_MC_less_delta_phi_cuts->AddText("Analysis cuts");
    L0_L0_text_MC_less_delta_phi_cuts->AddText(Form("P = %.3f #pm %.3f", L0_L0_slope_less_delta_phi_cuts, L0_L0_slope_less_delta_phi_cuts_err ));
    L0_L0_text_MC_less_delta_phi_cuts->AddText(Form("P_{ME} = %.3f", fit_L0_L0_after_less_delta_phi_cuts_ME->GetParameter(1)/(L0_alpha*L0_alpha) ));
    if(corr_err == 1) L0_L0_text_MC_less_delta_phi_cuts->AddText(Form("Slope sys. err: %.1f %%", L0_L0_slope_less_delta_phi_cuts_sys_err*100 ));
    if(corr_err == 0) L0_L0_text_MC_less_delta_phi_cuts->AddText(Form("Slope sys. err: %.1f %%", L0_L0_slope_less_delta_phi_cuts_sys_err_no_corr*100 ));
    L0_L0_text_MC_less_delta_phi_cuts->SetFillColorAlpha(0, 0.01);
    L0_L0_text_MC_less_delta_phi_cuts->Draw("same");

    L0_L0_cosThetaProdPlane_less_delta_phi_cuts_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_L/L0_L0_cosThetaProdPlane_less_delta_phi_cuts_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------

    TCanvas *L0_L0_cosThetaProdPlane_less_delta_phi_cuts_corr_can = new TCanvas(Form("L0_L0_cosThetaProdPlane_less_delta_phi_cuts_corr_can_%i", delta_eta_bin), Form("L0_L0_cosThetaProdPlane_less_delta_phi_cuts_corr_can_%i", delta_eta_bin), 1200, 1000);
    L0_L0_cosThetaProdPlane_less_delta_phi_cuts_corr_can->cd();

    TF1 *fit_L0_L0_after_less_delta_phi_cuts_corr = new TF1("fit_L0_L0_after_less_delta_phi_cuts_corr", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0_after_less_delta_phi_cuts_corr->SetParameters(1000, 0.10);

    TH1D *L0_L0_cosThetaProdPlane_less_delta_phi_cuts_corr = (TH1D*)L0_L0_cosThetaProdPlane_less_delta_phi_cuts->Clone(Form("L0_L0_cosThetaProdPlane_less_delta_phi_cuts_corr_%i", delta_eta_bin));
    L0_L0_cosThetaProdPlane_less_delta_phi_cuts_corr->Divide(L0_L0_cosThetaProdPlane_ME_less_delta_phi_cuts); //correct
    L0_L0_cosThetaProdPlane_less_delta_phi_cuts_corr->Scale(nL0L0_less_delta_phi_cuts/L0_L0_cosThetaProdPlane_less_delta_phi_cuts_corr->Integral()); //scale back
    L0_L0_cosThetaProdPlane_less_delta_phi_cuts_corr->Scale(1./L0_L0_cosThetaProdPlane_less_delta_phi_cuts_corr->GetXaxis()->GetBinWidth(1)); //divide by bin width
    //L0_L0_cosThetaProdPlane_less_delta_phi_cuts_corr->Scale(nL0L0_less_delta_phi_cuts);
    L0_L0_cosThetaProdPlane_less_delta_phi_cuts_corr->Fit(fit_L0_L0_after_less_delta_phi_cuts_corr, "i 0 r");
    L0_L0_cosThetaProdPlane_less_delta_phi_cuts_corr->SetMinimum(0);
    L0_L0_cosThetaProdPlane_less_delta_phi_cuts_corr->Draw("p e");

    fit_L0_L0_after_less_delta_phi_cuts_corr->SetLineColor(1);
    fit_L0_L0_after_less_delta_phi_cuts_corr->Draw("same");

    float L0_L0_slope_less_delta_phi_cuts_corr = fit_L0_L0_after_less_delta_phi_cuts_corr->GetParameter(1)/(L0_alpha*L0_alpha);
    float L0_L0_slope_less_delta_phi_cuts_corr_err = fit_L0_L0_after_less_delta_phi_cuts_corr->GetParError(1)/(L0_alpha*L0_alpha);

    ResidualPolarization_delta_phi_less_hist[delta_eta_bin]->SetBinContent(2, L0_L0_slope_less_delta_phi_cuts_corr);
    ResidualPolarization_delta_phi_less_hist[delta_eta_bin]->SetBinError(2, L0_L0_slope_less_delta_phi_cuts_corr_err);


    TPaveText *L0_L0_text_MC_less_delta_phi_cuts_corr = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    L0_L0_text_MC_less_delta_phi_cuts_corr->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0_L0_text_MC_less_delta_phi_cuts_corr->AddText("Minimum bias");
    L0_L0_text_MC_less_delta_phi_cuts_corr->AddText("#Lambda-#Lambda");
    L0_L0_text_MC_less_delta_phi_cuts_corr->AddText("True MC");
    L0_L0_text_MC_less_delta_phi_cuts_corr->AddText("Analysis cuts");
    L0_L0_text_MC_less_delta_phi_cuts_corr->AddText(Form("P = %.3f #pm %.3f", L0_L0_slope_less_delta_phi_cuts_corr, L0_L0_slope_less_delta_phi_cuts_corr_err ));
    //if(corr_err == 1) L0_L0_text_MC_less_delta_phi_cuts_corr->AddText(Form("Slope sys. err: %.1f %%", L0_L0_slope_less_delta_phi_cuts_corr_sys_err*100 ));
    //if(corr_err == 0) L0_L0_text_MC_less_delta_phi_cuts_corr->AddText(Form("Slope sys. err: %.1f %%", L0_L0_slope_less_delta_phi_cuts_corr_sys_err_no_corr*100 ));
    L0_L0_text_MC_less_delta_phi_cuts_corr->SetFillColorAlpha(0, 0.01);
    L0_L0_text_MC_less_delta_phi_cuts_corr->Draw("same");


    L0_L0_cosThetaProdPlane_less_delta_phi_cuts_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_L/L0_L0_cosThetaProdPlane_less_delta_phi_cuts_corr_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------

    TCanvas *L0_L0_cosThetaProdPlane_more_delta_phi_cuts_can = new TCanvas(Form("L0_L0_cosThetaProdPlane_more_delta_phi_cuts_can_%i", delta_eta_bin), Form("L0_L0_cosThetaProdPlane_more_delta_phi_cuts_can_%i", delta_eta_bin), 1200, 1000);
    L0_L0_cosThetaProdPlane_more_delta_phi_cuts_can->cd();

    TF1 *fit_L0_L0_after_more_delta_phi_cuts = new TF1("fit_L0_L0_after_more_delta_phi_cuts", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0_after_more_delta_phi_cuts->SetParameters(1000, 0.10);

    TF1 *fit_L0_L0_after_more_delta_phi_cuts_ME = new TF1("fit_L0_L0_after_more_delta_phi_cuts_ME", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0_after_more_delta_phi_cuts_ME->SetParameters(1000, 0.10);


    //TH1D *L0_L0_cosThetaProdPlane_more_delta_phi_cuts = L0_L0_cos_theta_star_vs_delta_phi_cuts_2_hist->ProjectionX(Form("proj_SE_2_L_L_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
    TH1D *L0_L0_cosThetaProdPlane_more_delta_phi_cuts = L0_L0_cos_theta_star_vs_delta_phi_for_corr_cuts_hist->ProjectionX(Form("proj_SE_2_L_L_%i", delta_eta_bin), 5*delta_eta_bin+1, 20);
    float nL0L0_more_delta_phi_cuts = L0_L0_cosThetaProdPlane_more_delta_phi_cuts->Integral();

    //TH1D *L0_L0_cosThetaProdPlane_ME_more_delta_phi_cuts = L0_L0_cos_theta_star_vs_delta_phi_ME_cuts_2_hist->ProjectionX(Form("proj_ME_2_L_L_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
    TH1D *L0_L0_cosThetaProdPlane_ME_more_delta_phi_cuts = L0_L0_cos_theta_star_vs_delta_phi_for_corr_ME_cuts_hist->ProjectionX(Form("proj_ME_2_L_L_%i", delta_eta_bin), 5*delta_eta_bin+1, 20);
    L0_L0_cosThetaProdPlane_ME_more_delta_phi_cuts->Sumw2();
    L0_L0_cosThetaProdPlane_ME_more_delta_phi_cuts->Scale(nL0L0_more_delta_phi_cuts/L0_L0_cosThetaProdPlane_ME_more_delta_phi_cuts->Integral());
    L0_L0_cosThetaProdPlane_ME_more_delta_phi_cuts->Scale(1./L0_L0_cosThetaProdPlane_ME_more_delta_phi_cuts->GetXaxis()->GetBinWidth(1));
    L0_L0_cosThetaProdPlane_ME_more_delta_phi_cuts->Fit(fit_L0_L0_after_more_delta_phi_cuts_ME, "i 0 r");
    //L0_L0_cosThetaProdPlane_ME->Draw("p e");


    L0_L0_cosThetaProdPlane_more_delta_phi_cuts->Sumw2();
    L0_L0_cosThetaProdPlane_more_delta_phi_cuts->SetMarkerStyle(20);
    L0_L0_cosThetaProdPlane_more_delta_phi_cuts->SetMarkerColor(kRed);
    L0_L0_cosThetaProdPlane_more_delta_phi_cuts->SetLineColor(kRed);
    L0_L0_cosThetaProdPlane_more_delta_phi_cuts->GetXaxis()->SetTitle("cos(#theta*)");
    L0_L0_cosThetaProdPlane_more_delta_phi_cuts->GetXaxis()->CenterTitle();
    L0_L0_cosThetaProdPlane_more_delta_phi_cuts->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    L0_L0_cosThetaProdPlane_more_delta_phi_cuts->GetYaxis()->CenterTitle();
    //L0_L0_cosThetaProdPlane_more_delta_phi_cuts->Divide(L0_L0_cosThetaProdPlane_ME_more_delta_phi_cuts); //correct using ME
    //L0_L0_cosThetaProdPlane_more_delta_phi_cuts->Scale(nL0L0_more_delta_phi_cuts/L0_L0_cosThetaProdPlane_more_delta_phi_cuts->Integral()); //scale back
    L0_L0_cosThetaProdPlane_more_delta_phi_cuts->Scale(1./L0_L0_cosThetaProdPlane_more_delta_phi_cuts->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
    L0_L0_cosThetaProdPlane_more_delta_phi_cuts->Fit(fit_L0_L0_after_more_delta_phi_cuts, "i 0 r");
    L0_L0_cosThetaProdPlane_more_delta_phi_cuts->SetMinimum(0);
    //L0_L0_cosThetaProdPlane_more_delta_phi_cuts->GetYaxis()->SetRangeUser(-20e9, 0);
    L0_L0_cosThetaProdPlane_more_delta_phi_cuts->Draw("p e");

    L0_L0_cosThetaProdPlane_ME_more_delta_phi_cuts->Draw("p e same");

    fit_L0_L0_after_more_delta_phi_cuts->SetLineColor(1);
    fit_L0_L0_after_more_delta_phi_cuts->Draw("same");

    fit_L0_L0_after_more_delta_phi_cuts_ME->SetLineColor(kBlue);
    fit_L0_L0_after_more_delta_phi_cuts_ME->Draw("same");


    float L0_L0_slope_more_delta_phi_cuts = fit_L0_L0_after_more_delta_phi_cuts->GetParameter(1)/(L0_alpha*L0_alpha);
    float L0_L0_slope_more_delta_phi_cuts_err = fit_L0_L0_after_more_delta_phi_cuts->GetParError(1)/(L0_alpha*L0_alpha);

    //statistical error correction for systematic error
    float SysErrSlope_L0_L0_more_delta_phi_corr = sqrt( fabs( fit_L0_L0_after_more_delta_phi_cuts->GetParError(1)*fit_L0_L0_after_more_delta_phi_cuts->GetParError(1) - fit_L0_L0_after_more_delta_phi_cuts_ME->GetParError(1)*fit_L0_L0_after_more_delta_phi_cuts_ME->GetParError(1) ) );

    float SysErrSlope_L0_L0_more_delta_phi_work = ( fabs( fit_L0_L0_after_more_delta_phi_cuts->GetParameter(1) - fit_L0_L0_after_more_delta_phi_cuts_ME->GetParameter(1) ) - SysErrSlope_L0_L0_more_delta_phi_corr )/fabs(fit_L0_L0_after_more_delta_phi_cuts->GetParameter(1));

    float L0_L0_slope_more_delta_phi_cuts_sys_err = 0;

    if( SysErrSlope_L0_L0_more_delta_phi_work > 0 ) L0_L0_slope_more_delta_phi_cuts_sys_err = SysErrSlope_L0_L0_more_delta_phi_work; //store sys. err. only if it is larger than statistical fluctuations

    float L0_L0_slope_more_delta_phi_cuts_sys_err_no_corr = fabs( fit_L0_L0_after_more_delta_phi_cuts->GetParameter(1) - fit_L0_L0_after_more_delta_phi_cuts_ME->GetParameter(1) )/fabs(fit_L0_L0_after_more_delta_phi_cuts->GetParameter(1));

    if(corr_err == 1) SysErrSlope_delta_phi_more_hist[delta_eta_bin]->SetBinContent(2, L0_L0_slope_more_delta_phi_cuts_sys_err);
    if(corr_err == 0) SysErrSlope_delta_phi_more_hist[delta_eta_bin]->SetBinContent(2, L0_L0_slope_more_delta_phi_cuts_sys_err_no_corr);

    TPaveText *L0_L0_text_MC_more_delta_phi_cuts = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    L0_L0_text_MC_more_delta_phi_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0_L0_text_MC_more_delta_phi_cuts->AddText("Minimum bias");
    L0_L0_text_MC_more_delta_phi_cuts->AddText("#Lambda-#Lambda");
    L0_L0_text_MC_more_delta_phi_cuts->AddText("True MC");
    L0_L0_text_MC_more_delta_phi_cuts->AddText("Analysis cuts");
    L0_L0_text_MC_more_delta_phi_cuts->AddText(Form("P = %.3f #pm %.3f", L0_L0_slope_more_delta_phi_cuts, L0_L0_slope_more_delta_phi_cuts_err ));
    L0_L0_text_MC_more_delta_phi_cuts->AddText(Form("P_{ME} = %.3f", fit_L0_L0_after_more_delta_phi_cuts_ME->GetParameter(1)/(L0_alpha*L0_alpha) ));
    if(corr_err == 1) L0_L0_text_MC_more_delta_phi_cuts->AddText(Form("Slope sys. err: %.1f %%", L0_L0_slope_more_delta_phi_cuts_sys_err*100 ));
    if(corr_err == 0) L0_L0_text_MC_more_delta_phi_cuts->AddText(Form("Slope sys. err: %.1f %%", L0_L0_slope_more_delta_phi_cuts_sys_err_no_corr*100 ));
    L0_L0_text_MC_more_delta_phi_cuts->SetFillColorAlpha(0, 0.01);
    L0_L0_text_MC_more_delta_phi_cuts->Draw("same");

    L0_L0_cosThetaProdPlane_more_delta_phi_cuts_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_L/L0_L0_cosThetaProdPlane_more_delta_phi_cuts_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------

    TCanvas *L0_L0_cosThetaProdPlane_more_delta_phi_cuts_corr_can = new TCanvas(Form("L0_L0_cosThetaProdPlane_more_delta_phi_cuts_corr_can_%i", delta_eta_bin), Form("L0_L0_cosThetaProdPlane_more_delta_phi_cuts_corr_can_%i", delta_eta_bin), 1200, 1000);
    L0_L0_cosThetaProdPlane_more_delta_phi_cuts_corr_can->cd();

    TF1 *fit_L0_L0_after_more_delta_phi_cuts_corr = new TF1("fit_L0_L0_after_more_delta_phi_cuts_corr", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0_after_more_delta_phi_cuts_corr->SetParameters(1000, 0.10);

    TH1D *L0_L0_cosThetaProdPlane_more_delta_phi_cuts_corr = (TH1D*)L0_L0_cosThetaProdPlane_more_delta_phi_cuts->Clone(Form("L0_L0_cosThetaProdPlane_more_delta_phi_cuts_corr_%i", delta_eta_bin));
    L0_L0_cosThetaProdPlane_more_delta_phi_cuts_corr->Divide(L0_L0_cosThetaProdPlane_ME_more_delta_phi_cuts); //correct
    L0_L0_cosThetaProdPlane_more_delta_phi_cuts_corr->Scale(nL0L0_more_delta_phi_cuts/L0_L0_cosThetaProdPlane_more_delta_phi_cuts_corr->Integral()); //scale back
    L0_L0_cosThetaProdPlane_more_delta_phi_cuts_corr->Scale(1./L0_L0_cosThetaProdPlane_more_delta_phi_cuts_corr->GetXaxis()->GetBinWidth(1)); //divide by bin width
    //L0_L0_cosThetaProdPlane_more_delta_phi_cuts_corr->Scale(nL0L0_more_delta_phi_cuts);
    L0_L0_cosThetaProdPlane_more_delta_phi_cuts_corr->Fit(fit_L0_L0_after_more_delta_phi_cuts_corr, "i 0 r");
    L0_L0_cosThetaProdPlane_more_delta_phi_cuts_corr->SetMinimum(0);
    L0_L0_cosThetaProdPlane_more_delta_phi_cuts_corr->Draw("p e");

    fit_L0_L0_after_more_delta_phi_cuts_corr->SetLineColor(1);
    fit_L0_L0_after_more_delta_phi_cuts_corr->Draw("same");

    float L0_L0_slope_more_delta_phi_cuts_corr = fit_L0_L0_after_more_delta_phi_cuts_corr->GetParameter(1)/(L0_alpha*L0_alpha);
    float L0_L0_slope_more_delta_phi_cuts_corr_err = fit_L0_L0_after_more_delta_phi_cuts_corr->GetParError(1)/(L0_alpha*L0_alpha);

    ResidualPolarization_delta_phi_more_hist[delta_eta_bin]->SetBinContent(2, L0_L0_slope_more_delta_phi_cuts_corr);
    ResidualPolarization_delta_phi_more_hist[delta_eta_bin]->SetBinError(2, L0_L0_slope_more_delta_phi_cuts_corr_err);


    TPaveText *L0_L0_text_MC_more_delta_phi_cuts_corr = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    L0_L0_text_MC_more_delta_phi_cuts_corr->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0_L0_text_MC_more_delta_phi_cuts_corr->AddText("Minimum bias");
    L0_L0_text_MC_more_delta_phi_cuts_corr->AddText("#Lambda-#Lambda");
    L0_L0_text_MC_more_delta_phi_cuts_corr->AddText("True MC");
    L0_L0_text_MC_more_delta_phi_cuts_corr->AddText("Analysis cuts");
    L0_L0_text_MC_more_delta_phi_cuts_corr->AddText(Form("P = %.3f #pm %.3f", L0_L0_slope_more_delta_phi_cuts_corr, L0_L0_slope_more_delta_phi_cuts_corr_err ));
    //if(corr_err == 1) L0_L0_text_MC_more_delta_phi_cuts_corr->AddText(Form("Slope sys. err: %.1f %%", L0_L0_slope_more_delta_phi_cuts_corr_sys_err*100 ));
    //if(corr_err == 0) L0_L0_text_MC_more_delta_phi_cuts_corr->AddText(Form("Slope sys. err: %.1f %%", L0_L0_slope_more_delta_phi_cuts_corr_sys_err_no_corr*100 ));
    L0_L0_text_MC_more_delta_phi_cuts_corr->SetFillColorAlpha(0, 0.01);
    L0_L0_text_MC_more_delta_phi_cuts_corr->Draw("same");


    L0_L0_cosThetaProdPlane_more_delta_phi_cuts_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_L/L0_L0_cosThetaProdPlane_more_delta_phi_cuts_corr_%i.png", delta_eta_bin));

    //________________________________________________________________________________________________________________________________________

    TCanvas *L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts_can = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts_can_%i", delta_eta_bin), Form("L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts_can_%i", delta_eta_bin), 1200, 1000);
    L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts_can->cd();

    TF1 *fit_L0bar_L0bar_after_less_delta_phi_cuts = new TF1("fit_L0bar_L0bar_after_less_delta_phi_cuts", "[0]*(1+[1]*x)", -1, 1);
    fit_L0bar_L0bar_after_less_delta_phi_cuts->SetParameters(1000, 0.10);

    TF1 *fit_L0bar_L0bar_after_less_delta_phi_cuts_ME = new TF1("fit_L0bar_L0bar_after_less_delta_phi_cuts_ME", "[0]*(1+[1]*x)", -1, 1);
    fit_L0bar_L0bar_after_less_delta_phi_cuts_ME->SetParameters(1000, 0.10);

    //TH1D *L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts = L0bar_L0bar_cos_theta_star_vs_delta_phi_cuts_hist->ProjectionX(Form("proj_SE_Lbar_Lbar_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
    TH1D *L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts = L0bar_L0bar_cos_theta_star_vs_delta_phi_for_corr_cuts_hist->ProjectionX(Form("proj_SE_Lbar_Lbar_%i", delta_eta_bin), 1, 5*delta_eta_bin+5);
    float nL0barL0bar_less_delta_phi_cuts = L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts->Integral();

    //TH1D *L0bar_L0bar_cosThetaProdPlane_ME_less_delta_phi_cuts = L0bar_L0bar_cos_theta_star_vs_delta_phi_ME_cuts_hist->ProjectionX(Form("proj_ME_Lbar_Lbar_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
    TH1D *L0bar_L0bar_cosThetaProdPlane_ME_less_delta_phi_cuts = L0bar_L0bar_cos_theta_star_vs_delta_phi_for_corr_ME_cuts_hist->ProjectionX(Form("proj_ME_Lbar_Lbar_%i", delta_eta_bin), 1, 5*delta_eta_bin+5);
    L0bar_L0bar_cosThetaProdPlane_ME_less_delta_phi_cuts->Sumw2();
    L0bar_L0bar_cosThetaProdPlane_ME_less_delta_phi_cuts->Scale(nL0barL0bar_less_delta_phi_cuts/L0bar_L0bar_cosThetaProdPlane_ME_less_delta_phi_cuts->Integral());
    L0bar_L0bar_cosThetaProdPlane_ME_less_delta_phi_cuts->Scale(1./L0bar_L0bar_cosThetaProdPlane_ME_less_delta_phi_cuts->GetXaxis()->GetBinWidth(1));
    L0bar_L0bar_cosThetaProdPlane_ME_less_delta_phi_cuts->Fit(fit_L0bar_L0bar_after_less_delta_phi_cuts_ME, "i 0 r");
    //L0bar_L0bar_cosThetaProdPlane_ME->Draw("p e");


    L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts->Sumw2();
    L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts->SetMarkerStyle(20);
    L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts->SetMarkerColor(kRed);
    L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts->SetLineColor(kRed);
    L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts->GetXaxis()->SetTitle("cos(#theta*)");
    L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts->GetXaxis()->CenterTitle();
    L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts->GetYaxis()->CenterTitle();
    //L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts->Divide(L0bar_L0bar_cosThetaProdPlane_ME_less_delta_phi_cuts); //correct using ME
    //L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts->Scale(nL0barL0bar_less_delta_phi_cuts/L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts->Integral()); //scale back
    L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts->Scale(1./L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
    L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts->Fit(fit_L0bar_L0bar_after_less_delta_phi_cuts, "i 0 r");
    L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts->SetMinimum(0);
    //L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts->GetYaxis()->SetRangeUser(-20e9, 0);
    L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts->Draw("p e");

    L0bar_L0bar_cosThetaProdPlane_ME_less_delta_phi_cuts->Draw("p e same");

    fit_L0bar_L0bar_after_less_delta_phi_cuts->SetLineColor(1);
    fit_L0bar_L0bar_after_less_delta_phi_cuts->Draw("same");

    fit_L0bar_L0bar_after_less_delta_phi_cuts_ME->SetLineColor(kBlue);
    fit_L0bar_L0bar_after_less_delta_phi_cuts_ME->Draw("same");


    float L0bar_L0bar_slope_less_delta_phi_cuts = fit_L0bar_L0bar_after_less_delta_phi_cuts->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
    float L0bar_L0bar_slope_less_delta_phi_cuts_err = fit_L0bar_L0bar_after_less_delta_phi_cuts->GetParError(1)/(L0bar_alpha*L0bar_alpha);

    //statistical error correction for systematic error
    float SysErrSlope_L0bar_L0bar_less_delta_phi_corr = sqrt( fabs( fit_L0bar_L0bar_after_less_delta_phi_cuts->GetParError(1)*fit_L0bar_L0bar_after_less_delta_phi_cuts->GetParError(1) - fit_L0bar_L0bar_after_less_delta_phi_cuts_ME->GetParError(1)*fit_L0bar_L0bar_after_less_delta_phi_cuts_ME->GetParError(1) ) );

    float SysErrSlope_L0bar_L0bar_less_delta_phi_work = ( fabs( fit_L0bar_L0bar_after_less_delta_phi_cuts->GetParameter(1) - fit_L0bar_L0bar_after_less_delta_phi_cuts_ME->GetParameter(1) ) - SysErrSlope_L0bar_L0bar_less_delta_phi_corr )/fabs(fit_L0bar_L0bar_after_less_delta_phi_cuts->GetParameter(1));

    float L0bar_L0bar_slope_less_delta_phi_cuts_sys_err = 0;

    if( SysErrSlope_L0bar_L0bar_less_delta_phi_work > 0 ) L0bar_L0bar_slope_less_delta_phi_cuts_sys_err = SysErrSlope_L0bar_L0bar_less_delta_phi_work; //store sys. err. only if it is larger than statistical fluctuations

    float L0bar_L0bar_slope_less_delta_phi_cuts_sys_err_no_corr = fabs( fit_L0bar_L0bar_after_less_delta_phi_cuts->GetParameter(1) - fit_L0bar_L0bar_after_less_delta_phi_cuts_ME->GetParameter(1) )/fabs(fit_L0bar_L0bar_after_less_delta_phi_cuts->GetParameter(1));

    if(corr_err == 1) SysErrSlope_delta_phi_less_hist[delta_eta_bin]->SetBinContent(3, L0bar_L0bar_slope_less_delta_phi_cuts_sys_err);
    if(corr_err == 0) SysErrSlope_delta_phi_less_hist[delta_eta_bin]->SetBinContent(3, L0bar_L0bar_slope_less_delta_phi_cuts_sys_err_no_corr);

    TPaveText *L0bar_L0bar_text_MC_less_delta_phi_cuts = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    L0bar_L0bar_text_MC_less_delta_phi_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0bar_L0bar_text_MC_less_delta_phi_cuts->AddText("Minimum bias");
    L0bar_L0bar_text_MC_less_delta_phi_cuts->AddText("#bar{#Lambda}-#bar{#Lambda}");
    L0bar_L0bar_text_MC_less_delta_phi_cuts->AddText("True MC");
    L0bar_L0bar_text_MC_less_delta_phi_cuts->AddText("Analysis cuts");
    L0bar_L0bar_text_MC_less_delta_phi_cuts->AddText(Form("P = %.3f #pm %.3f", L0bar_L0bar_slope_less_delta_phi_cuts, L0bar_L0bar_slope_less_delta_phi_cuts_err ));
    L0bar_L0bar_text_MC_less_delta_phi_cuts->AddText(Form("P_{ME} = %.3f", fit_L0bar_L0bar_after_less_delta_phi_cuts_ME->GetParameter(1)/(L0bar_alpha*L0bar_alpha) ));
    if(corr_err == 1) L0bar_L0bar_text_MC_less_delta_phi_cuts->AddText(Form("Slope sys. err: %.1f %%", L0bar_L0bar_slope_less_delta_phi_cuts_sys_err*100 ));
    if(corr_err == 0) L0bar_L0bar_text_MC_less_delta_phi_cuts->AddText(Form("Slope sys. err: %.1f %%", L0bar_L0bar_slope_less_delta_phi_cuts_sys_err_no_corr*100 ));
    L0bar_L0bar_text_MC_less_delta_phi_cuts->SetFillColorAlpha(0, 0.01);
    L0bar_L0bar_text_MC_less_delta_phi_cuts->Draw("same");

    L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/Lbar_Lbar/L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------

    TCanvas *L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts_corr_can = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts_corr_can_%i", delta_eta_bin), Form("L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts_corr_can_%i", delta_eta_bin), 1200, 1000);
    L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts_corr_can->cd();

    TF1 *fit_L0bar_L0bar_after_less_delta_phi_cuts_corr = new TF1("fit_L0bar_L0bar_after_less_delta_phi_cuts_corr", "[0]*(1+[1]*x)", -1, 1);
    fit_L0bar_L0bar_after_less_delta_phi_cuts_corr->SetParameters(1000, 0.10);

    TH1D *L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts_corr = (TH1D*)L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts->Clone(Form("L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts_corr_%i", delta_eta_bin));
    L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts_corr->Divide(L0bar_L0bar_cosThetaProdPlane_ME_less_delta_phi_cuts); //correct
    L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts_corr->Scale(nL0barL0bar_less_delta_phi_cuts/L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts_corr->Integral()); //scale back
    L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts_corr->Scale(1./L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts_corr->GetXaxis()->GetBinWidth(1)); //divide by bin width
    //L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts_corr->Scale(nL0barL0bar_less_delta_phi_cuts);
    L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts_corr->Fit(fit_L0bar_L0bar_after_less_delta_phi_cuts_corr, "i 0 r");
    L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts_corr->SetMinimum(0);
    L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts_corr->Draw("p e");

    fit_L0bar_L0bar_after_less_delta_phi_cuts_corr->SetLineColor(1);
    fit_L0bar_L0bar_after_less_delta_phi_cuts_corr->Draw("same");

    float L0bar_L0bar_slope_less_delta_phi_cuts_corr = fit_L0bar_L0bar_after_less_delta_phi_cuts_corr->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
    float L0bar_L0bar_slope_less_delta_phi_cuts_corr_err = fit_L0bar_L0bar_after_less_delta_phi_cuts_corr->GetParError(1)/(L0bar_alpha*L0bar_alpha);

    ResidualPolarization_delta_phi_less_hist[delta_eta_bin]->SetBinContent(3, L0bar_L0bar_slope_less_delta_phi_cuts_corr);
    ResidualPolarization_delta_phi_less_hist[delta_eta_bin]->SetBinError(3, L0bar_L0bar_slope_less_delta_phi_cuts_corr_err);


    TPaveText *L0bar_L0bar_text_MC_less_delta_phi_cuts_corr = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    L0bar_L0bar_text_MC_less_delta_phi_cuts_corr->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0bar_L0bar_text_MC_less_delta_phi_cuts_corr->AddText("Minimum bias");
    L0bar_L0bar_text_MC_less_delta_phi_cuts_corr->AddText("#bar{#Lambda}-#bar{#Lambda}");
    L0bar_L0bar_text_MC_less_delta_phi_cuts_corr->AddText("True MC");
    L0bar_L0bar_text_MC_less_delta_phi_cuts_corr->AddText("Analysis cuts");
    L0bar_L0bar_text_MC_less_delta_phi_cuts_corr->AddText(Form("P = %.3f #pm %.3f", L0bar_L0bar_slope_less_delta_phi_cuts_corr, L0bar_L0bar_slope_less_delta_phi_cuts_corr_err ));
    //if(corr_err == 1) L0bar_L0bar_text_MC_less_delta_phi_cuts_corr->AddText(Form("Slope sys. err: %.1f %%", L0bar_L0bar_slope_less_delta_phi_cuts_corr_sys_err*100 ));
    //if(corr_err == 0) L0bar_L0bar_text_MC_less_delta_phi_cuts_corr->AddText(Form("Slope sys. err: %.1f %%", L0bar_L0bar_slope_less_delta_phi_cuts_corr_sys_err_no_corr*100 ));
    L0bar_L0bar_text_MC_less_delta_phi_cuts_corr->SetFillColorAlpha(0, 0.01);
    L0bar_L0bar_text_MC_less_delta_phi_cuts_corr->Draw("same");


    L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/Lbar_Lbar/L0bar_L0bar_cosThetaProdPlane_less_delta_phi_cuts_corr_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------

    TCanvas *L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts_can = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts_can_%i", delta_eta_bin), Form("L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts_can_%i", delta_eta_bin), 1200, 1000);
    L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts_can->cd();

    TF1 *fit_L0bar_L0bar_after_more_delta_phi_cuts = new TF1("fit_L0bar_L0bar_after_more_delta_phi_cuts", "[0]*(1+[1]*x)", -1, 1);
    fit_L0bar_L0bar_after_more_delta_phi_cuts->SetParameters(1000, 0.10);

    TF1 *fit_L0bar_L0bar_after_more_delta_phi_cuts_ME = new TF1("fit_L0bar_L0bar_after_more_delta_phi_cuts_ME", "[0]*(1+[1]*x)", -1, 1);
    fit_L0bar_L0bar_after_more_delta_phi_cuts_ME->SetParameters(1000, 0.10);


    //TH1D *L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts = L0bar_L0bar_cos_theta_star_vs_delta_phi_cuts_2_hist->ProjectionX(Form("proj_SE_2_Lbar_Lbar_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
    TH1D *L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts = L0bar_L0bar_cos_theta_star_vs_delta_phi_for_corr_cuts_hist->ProjectionX(Form("proj_SE_2_Lbar_Lbar_%i", delta_eta_bin), 5*delta_eta_bin+1, 20);
    float nL0barL0bar_more_delta_phi_cuts = L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts->Integral();

    //TH1D *L0bar_L0bar_cosThetaProdPlane_ME_more_delta_phi_cuts = L0bar_L0bar_cos_theta_star_vs_delta_phi_ME_cuts_2_hist->ProjectionX(Form("proj_ME_2_Lbar_Lbar_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
    TH1D *L0bar_L0bar_cosThetaProdPlane_ME_more_delta_phi_cuts = L0bar_L0bar_cos_theta_star_vs_delta_phi_for_corr_ME_cuts_hist->ProjectionX(Form("proj_ME_2_Lbar_Lbar_%i", delta_eta_bin), 5*delta_eta_bin+1, 20);
    L0bar_L0bar_cosThetaProdPlane_ME_more_delta_phi_cuts->Sumw2();
    L0bar_L0bar_cosThetaProdPlane_ME_more_delta_phi_cuts->Scale(nL0barL0bar_more_delta_phi_cuts/L0bar_L0bar_cosThetaProdPlane_ME_more_delta_phi_cuts->Integral());
    L0bar_L0bar_cosThetaProdPlane_ME_more_delta_phi_cuts->Scale(1./L0bar_L0bar_cosThetaProdPlane_ME_more_delta_phi_cuts->GetXaxis()->GetBinWidth(1));
    L0bar_L0bar_cosThetaProdPlane_ME_more_delta_phi_cuts->Fit(fit_L0bar_L0bar_after_more_delta_phi_cuts_ME, "i 0 r");
    //L0bar_L0bar_cosThetaProdPlane_ME->Draw("p e");


    L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts->Sumw2();
    L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts->SetMarkerStyle(20);
    L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts->SetMarkerColor(kRed);
    L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts->SetLineColor(kRed);
    L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts->GetXaxis()->SetTitle("cos(#theta*)");
    L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts->GetXaxis()->CenterTitle();
    L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts->GetYaxis()->CenterTitle();
    //L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts->Divide(L0bar_L0bar_cosThetaProdPlane_ME_more_delta_phi_cuts); //correct using ME
    //L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts->Scale(nL0barL0bar_more_delta_phi_cuts/L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts->Integral()); //scale back
    L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts->Scale(1./L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
    L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts->Fit(fit_L0bar_L0bar_after_more_delta_phi_cuts, "i 0 r");
    L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts->SetMinimum(0);
    //L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts->GetYaxis()->SetRangeUser(-20e9, 0);
    L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts->Draw("p e");

    L0bar_L0bar_cosThetaProdPlane_ME_more_delta_phi_cuts->Draw("p e same");

    fit_L0bar_L0bar_after_more_delta_phi_cuts->SetLineColor(1);
    fit_L0bar_L0bar_after_more_delta_phi_cuts->Draw("same");

    fit_L0bar_L0bar_after_more_delta_phi_cuts_ME->SetLineColor(kBlue);
    fit_L0bar_L0bar_after_more_delta_phi_cuts_ME->Draw("same");


    float L0bar_L0bar_slope_more_delta_phi_cuts = fit_L0bar_L0bar_after_more_delta_phi_cuts->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
    float L0bar_L0bar_slope_more_delta_phi_cuts_err = fit_L0bar_L0bar_after_more_delta_phi_cuts->GetParError(1)/(L0bar_alpha*L0bar_alpha);

    //statistical error correction for systematic error
    float SysErrSlope_L0bar_L0bar_more_delta_phi_corr = sqrt( fabs( fit_L0bar_L0bar_after_more_delta_phi_cuts->GetParError(1)*fit_L0bar_L0bar_after_more_delta_phi_cuts->GetParError(1) - fit_L0bar_L0bar_after_more_delta_phi_cuts_ME->GetParError(1)*fit_L0bar_L0bar_after_more_delta_phi_cuts_ME->GetParError(1) ) );

    float SysErrSlope_L0bar_L0bar_more_delta_phi_work = ( fabs( fit_L0bar_L0bar_after_more_delta_phi_cuts->GetParameter(1) - fit_L0bar_L0bar_after_more_delta_phi_cuts_ME->GetParameter(1) ) - SysErrSlope_L0bar_L0bar_more_delta_phi_corr )/fabs(fit_L0bar_L0bar_after_more_delta_phi_cuts->GetParameter(1));

    float L0bar_L0bar_slope_more_delta_phi_cuts_sys_err = 0;

    if( SysErrSlope_L0bar_L0bar_more_delta_phi_work > 0 ) L0bar_L0bar_slope_more_delta_phi_cuts_sys_err = SysErrSlope_L0bar_L0bar_more_delta_phi_work; //store sys. err. only if it is larger than statistical fluctuations

    float L0bar_L0bar_slope_more_delta_phi_cuts_sys_err_no_corr = fabs( fit_L0bar_L0bar_after_more_delta_phi_cuts->GetParameter(1) - fit_L0bar_L0bar_after_more_delta_phi_cuts_ME->GetParameter(1) )/fabs(fit_L0bar_L0bar_after_more_delta_phi_cuts->GetParameter(1));

    if(corr_err == 1) SysErrSlope_delta_phi_more_hist[delta_eta_bin]->SetBinContent(3, L0bar_L0bar_slope_more_delta_phi_cuts_sys_err);
    if(corr_err == 0) SysErrSlope_delta_phi_more_hist[delta_eta_bin]->SetBinContent(3, L0bar_L0bar_slope_more_delta_phi_cuts_sys_err_no_corr);

    TPaveText *L0bar_L0bar_text_MC_more_delta_phi_cuts = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    L0bar_L0bar_text_MC_more_delta_phi_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0bar_L0bar_text_MC_more_delta_phi_cuts->AddText("Minimum bias");
    L0bar_L0bar_text_MC_more_delta_phi_cuts->AddText("#bar{#Lambda}-#bar{#Lambda}");
    L0bar_L0bar_text_MC_more_delta_phi_cuts->AddText("True MC");
    L0bar_L0bar_text_MC_more_delta_phi_cuts->AddText("Analysis cuts");
    L0bar_L0bar_text_MC_more_delta_phi_cuts->AddText(Form("P = %.3f #pm %.3f", L0bar_L0bar_slope_more_delta_phi_cuts, L0bar_L0bar_slope_more_delta_phi_cuts_err ));
    L0bar_L0bar_text_MC_more_delta_phi_cuts->AddText(Form("P_{ME} = %.3f", fit_L0bar_L0bar_after_more_delta_phi_cuts_ME->GetParameter(1)/(L0bar_alpha*L0bar_alpha) ));
    if(corr_err == 1) L0bar_L0bar_text_MC_more_delta_phi_cuts->AddText(Form("Slope sys. err: %.1f %%", L0bar_L0bar_slope_more_delta_phi_cuts_sys_err*100 ));
    if(corr_err == 0) L0bar_L0bar_text_MC_more_delta_phi_cuts->AddText(Form("Slope sys. err: %.1f %%", L0bar_L0bar_slope_more_delta_phi_cuts_sys_err_no_corr*100 ));
    L0bar_L0bar_text_MC_more_delta_phi_cuts->SetFillColorAlpha(0, 0.01);
    L0bar_L0bar_text_MC_more_delta_phi_cuts->Draw("same");

    L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/Lbar_Lbar/L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------

    TCanvas *L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts_corr_can = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts_corr_can_%i", delta_eta_bin), Form("L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts_corr_can_%i", delta_eta_bin), 1200, 1000);
    L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts_corr_can->cd();

    TF1 *fit_L0bar_L0bar_after_more_delta_phi_cuts_corr = new TF1("fit_L0bar_L0bar_after_more_delta_phi_cuts_corr", "[0]*(1+[1]*x)", -1, 1);
    fit_L0bar_L0bar_after_more_delta_phi_cuts_corr->SetParameters(1000, 0.10);

    TH1D *L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts_corr = (TH1D*)L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts->Clone(Form("L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts_corr_%i", delta_eta_bin));
    L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts_corr->Divide(L0bar_L0bar_cosThetaProdPlane_ME_more_delta_phi_cuts); //correct
    L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts_corr->Scale(nL0barL0bar_more_delta_phi_cuts/L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts_corr->Integral()); //scale back
    L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts_corr->Scale(1./L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts_corr->GetXaxis()->GetBinWidth(1)); //divide by bin width
    //L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts_corr->Scale(nL0barL0bar_more_delta_phi_cuts);
    L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts_corr->Fit(fit_L0bar_L0bar_after_more_delta_phi_cuts_corr, "i 0 r");
    L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts_corr->SetMinimum(0);
    L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts_corr->Draw("p e");

    fit_L0bar_L0bar_after_more_delta_phi_cuts_corr->SetLineColor(1);
    fit_L0bar_L0bar_after_more_delta_phi_cuts_corr->Draw("same");

    float L0bar_L0bar_slope_more_delta_phi_cuts_corr = fit_L0bar_L0bar_after_more_delta_phi_cuts_corr->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
    float L0bar_L0bar_slope_more_delta_phi_cuts_corr_err = fit_L0bar_L0bar_after_more_delta_phi_cuts_corr->GetParError(1)/(L0bar_alpha*L0bar_alpha);

    ResidualPolarization_delta_phi_more_hist[delta_eta_bin]->SetBinContent(3, L0bar_L0bar_slope_more_delta_phi_cuts_corr);
    ResidualPolarization_delta_phi_more_hist[delta_eta_bin]->SetBinError(3, L0bar_L0bar_slope_more_delta_phi_cuts_corr_err);


    TPaveText *L0bar_L0bar_text_MC_more_delta_phi_cuts_corr = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    L0bar_L0bar_text_MC_more_delta_phi_cuts_corr->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0bar_L0bar_text_MC_more_delta_phi_cuts_corr->AddText("Minimum bias");
    L0bar_L0bar_text_MC_more_delta_phi_cuts_corr->AddText("#bar{#Lambda}-#bar{#Lambda}");
    L0bar_L0bar_text_MC_more_delta_phi_cuts_corr->AddText("True MC");
    L0bar_L0bar_text_MC_more_delta_phi_cuts_corr->AddText("Analysis cuts");
    L0bar_L0bar_text_MC_more_delta_phi_cuts_corr->AddText(Form("P = %.3f #pm %.3f", L0bar_L0bar_slope_more_delta_phi_cuts_corr, L0bar_L0bar_slope_more_delta_phi_cuts_corr_err ));
    //if(corr_err == 1) L0bar_L0bar_text_MC_more_delta_phi_cuts_corr->AddText(Form("Slope sys. err: %.1f %%", L0bar_L0bar_slope_more_delta_phi_cuts_corr_sys_err*100 ));
    //if(corr_err == 0) L0bar_L0bar_text_MC_more_delta_phi_cuts_corr->AddText(Form("Slope sys. err: %.1f %%", L0bar_L0bar_slope_more_delta_phi_cuts_corr_sys_err_no_corr*100 ));
    L0bar_L0bar_text_MC_more_delta_phi_cuts_corr->SetFillColorAlpha(0, 0.01);
    L0bar_L0bar_text_MC_more_delta_phi_cuts_corr->Draw("same");


    L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/Lbar_Lbar/L0bar_L0bar_cosThetaProdPlane_more_delta_phi_cuts_corr_%i.png", delta_eta_bin));


    //---------------------------------------------

    SysErrSlope_delta_phi_less_hist[delta_eta_bin]->Write();
    SysErrSlope_delta_phi_more_hist[delta_eta_bin]->Write();

    ResidualPolarization_delta_phi_less_hist[delta_eta_bin]->Write();
    ResidualPolarization_delta_phi_more_hist[delta_eta_bin]->Write();

  }


  //couple of QA histograms

  TCanvas *L_Lbar_delta_phi_can = new TCanvas("L_Lbar_delta_phi_can", "L_Lbar_delta_phi_can", 1200, 1000);
  L_Lbar_delta_phi_can->cd();

  L0_L0bar_delta_phi_cuts_hist->SetLineColor(1);
  L0_L0bar_delta_phi_cuts_hist->GetXaxis()->SetTitle("#Delta#phi");
  L0_L0bar_delta_phi_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0bar_delta_phi_cuts_hist->SetMinimum(0);
  L0_L0bar_delta_phi_cuts_hist->Draw("hist");


  L0_L0bar_delta_phi_ME_cuts_hist->SetLineColor(kRed);
  L0_L0bar_delta_phi_ME_cuts_hist->GetXaxis()->SetTitle("#Delta#phi");
  L0_L0bar_delta_phi_ME_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0bar_delta_phi_ME_cuts_hist->Scale(L0_L0bar_delta_phi_cuts_hist->Integral()/L0_L0bar_delta_phi_ME_cuts_hist->Integral());
  L0_L0bar_delta_phi_ME_cuts_hist->SetMinimum(0);
  L0_L0bar_delta_phi_ME_cuts_hist->Draw("hist same");

  L_Lbar_delta_phi_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/L_Lbar_delta_phi.png");


  TFile *Delta_phi_file = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/output/Delta_phi_weight/Delta_phi_reweight.root", "recreate");
  Delta_phi_file->cd();

  TH1F *Delta_phi_reweight = (TH1F*)L0_L0bar_delta_phi_cuts_hist->Clone("Delta_phi_reweight");
  Delta_phi_reweight->Divide(L0_L0bar_delta_phi_ME_cuts_hist);
  Delta_phi_reweight->Write("Delta_phi_reweight");



  //---------------------------------------

  TCanvas *L_Lbar_delta_eta_can = new TCanvas("L_Lbar_delta_eta_can", "L_Lbar_delta_eta_can", 1200, 1000);
  L_Lbar_delta_eta_can->cd();

  L0_L0bar_delta_eta_cuts_hist->SetLineColor(1);
  L0_L0bar_delta_eta_cuts_hist->GetXaxis()->SetTitle("#Delta#eta");
  L0_L0bar_delta_eta_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0bar_delta_eta_cuts_hist->SetMinimum(0);
  L0_L0bar_delta_eta_cuts_hist->Draw("hist");


  L0_L0bar_delta_eta_ME_cuts_hist->SetLineColor(kRed);
  L0_L0bar_delta_eta_ME_cuts_hist->GetXaxis()->SetTitle("#Delta#eta");
  L0_L0bar_delta_eta_ME_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0bar_delta_eta_ME_cuts_hist->Scale(L0_L0bar_delta_eta_cuts_hist->Integral()/L0_L0bar_delta_eta_ME_cuts_hist->Integral());
  L0_L0bar_delta_eta_ME_cuts_hist->SetMinimum(0);
  L0_L0bar_delta_eta_ME_cuts_hist->Draw("hist same");

  L_Lbar_delta_eta_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/L_Lbar_delta_eta.png");

  //_______________________________________________________________________________________________________


  inFile->Close();
  sysErrFile->Close();
  Delta_phi_file->Close();

  //outFile->Close();

  return;
}
