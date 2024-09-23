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
void Read_PYTHIA_hists(const int energy = 510)
{
  const float L0_alpha = 0.732; //decay parameter of L0
  const float L0bar_alpha = -0.758; //decay paramteter of L0bar

  const int nEvents = 1e9;

  //load all files
  TFile *inFile;

  if(energy == 510)
  {
    inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run17/output_Lambda_pp_510_MB_1B_events.root", "READ");
    //inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/output_Lambda_pp_510_MB_1M_events.root", "READ");
  }

  if(energy == 200)
  {
    inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/output_Lambda_pp_200_MB_1B_events_hists_work.root", "READ");

    //inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/ME_tests_Delta_phi_06_2024/N_SE_in_ME/output_Lambda_pp_200_MB_1B_events_hists_1_SE_in_ME.root", "READ");
    //inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/ME_tests_Delta_phi_06_2024/N_SE_in_ME/output_Lambda_pp_200_MB_1B_events_hists_10_SE_in_ME.root", "READ");
    //inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/ME_tests_Delta_phi_06_2024/N_SE_in_ME/output_Lambda_pp_200_MB_1B_events_hists_100_SE_in_ME.root", "READ");
  }



  const int nPtBins = 8;
  float const pT_bins[nPtBins+1] = { 0., 0.5, 1.,1.5, 2., 2.5, 3., 4., 5.};

  const int nEtaBins = 3;
  float const eta_bins[nEtaBins+1] = { -1, -0.4, 0.4, 1 };
  //______________________________________________________________________________________________

  //histograms
  TH2F *L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_cuts_hist = (TH2F*)inFile->Get("L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_cuts_hist");

  TH2F *L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_ME_cuts_hist = (TH2F*)inFile->Get("L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_ME_cuts_hist");

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  //projections to Delta phi in individual dN/d cos(theta*) bins
  for(unsigned int bin = 0; bin < 10; bin++)
  {
    TH1D* Delta_phi_proj = L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_cuts_hist->ProjectionY(Form("Proj_%i", bin), bin+1, bin+1);

    TH1D* Delta_phi_ME_proj = L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_ME_cuts_hist->ProjectionY(Form("Proj_ME_%i", bin), bin+1, bin+1);

    TCanvas *Delta_phi_proj_can = new TCanvas(Form("Delta_phi_proj_can_%i", bin), Form("Delta_phi_proj_can_%i", bin), 1200, 1000);
    Delta_phi_proj_can->cd();

    Delta_phi_proj->GetXaxis()->SetTitle("#Delta#phi");
    Delta_phi_proj->GetXaxis()->CenterTitle();
    Delta_phi_proj->SetLineColor(1);
    Delta_phi_proj->Draw("hist");

    Delta_phi_ME_proj->SetLineColor(kRed);
    Delta_phi_ME_proj->Scale(Delta_phi_proj->Integral()/Delta_phi_ME_proj->Integral());
    Delta_phi_ME_proj->Draw("hist same");

    TPaveText *text = new TPaveText(0.5, 0.7, 0.8, 0.9, "NDC");
    text->AddText(Form("%.1f < dN/dcos(#theta*) < %.1f", -1+bin*0.2, -0.8+bin*0.2));
    text->SetFillColorAlpha(0, 0.01);
    text->Draw("same");

    Delta_phi_proj_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/Delta_phi_QA/L_Lbar/L0_L0bar_Delta_phi_cuts_%i.png", bin));

  }


  inFile->Close();

  return;

}
