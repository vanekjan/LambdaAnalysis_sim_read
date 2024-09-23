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
void Closure_test_K0s_new_ME(const int energy = 510, const int corr_err = 0)
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
  float const eta_bins[nEtaBins+1] = { -1, -0.2, 0.2, 1 };

  //TFile *outFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/output/K0s_cosThetaStar_eff_new_US_LS_work.root", "recreate");

  //load all files
  TFile *inFile;

  if(energy == 510)
  {
    inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run17/output_K0s_pp_510_MB_1B_events_new_hists.root", "READ");
  }
  else if(energy == 200)
  {
    //inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/output_K0s_pp_200_MB_1B_events_hists_work.root", "READ");


    //inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/2024_08_ME_weight_new/output_K0s_pp_200_MB_1B_events_hists_Delta_phi_quater.root", "READ");
    inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/2024_08_ME_weight_new/output_K0s_pp_200_MB_1B_events_hists_Delta_phi_third.root", "READ");
  }
  else
  {
    cout<<"Not a valid collision energy! Abborting!"<<endl;
    return;
  }

  //----------------------------------------------------------------

  TFile *sysErrFile;

  if(corr_err == 0) sysErrFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/output/sys_err/SysErrSlope_K0s_nocorr.root", "recreate");
  else if(corr_err == 1) sysErrFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/output/sys_err/SysErrSlope_K0s_corr.root", "recreate");
  else
  {
    cout<<"Wrong second argument"<<endl;
    return;
  }


  //histogram to store systematic error from slope difference in closure test
  //each bin is for one L charge combination: bin 1 - LLbar, bin 2 - LL, bin 3 - LbarLbar
  TH1F *SysErrSlope_hist = new TH1F("SysErrSlope_hist", "SysErrSlope_hist", 1, 0, 1);

  TH1F *SysErrSlope_delta_eta_hist[2];
  TH1F *ResidualPolarization_delta_eta_hist[2];

  TH1F *SysErrSlope_delta_phi_hist[2];
  TH1F *ResidualPolarization_delta_phi_hist[2];

  TH1F *SysErrSlope_delta_eta_delta_phi_hist[2];
  TH1F *ResidualPolarization_delta_eta_delta_phi_hist[2];

  for( unsigned int delta_eta_bin = 0; delta_eta_bin < 2; delta_eta_bin++ )
  {
    SysErrSlope_delta_eta_hist[delta_eta_bin] = new TH1F(Form("SysErrSlope_K0s_delta_eta_hist_%i", delta_eta_bin), Form("SysErrSlope_K0s_delta_eta_hist_%i", delta_eta_bin), 1, 0, 1);
    ResidualPolarization_delta_eta_hist[delta_eta_bin] = new TH1F(Form("ResidualPolarization_K0s_delta_eta_hist_%i", delta_eta_bin), Form("ResidualPolarization_K0s_delta_eta_hist_%i", delta_eta_bin), 1, 0, 1);

    SysErrSlope_delta_phi_hist[delta_eta_bin] = new TH1F(Form("SysErrSlope_K0s_delta_phi_hist_%i", delta_eta_bin), Form("SysErrSlope_K0s_delta_phi_hist_%i", delta_eta_bin), 1, 0, 1);
    ResidualPolarization_delta_phi_hist[delta_eta_bin] = new TH1F(Form("ResidualPolarization_K0s_delta_phi_hist_%i", delta_eta_bin), Form("ResidualPolarization_K0s_delta_phi_hist_%i", delta_eta_bin), 1, 0, 1);

    SysErrSlope_delta_eta_delta_phi_hist[delta_eta_bin] = new TH1F(Form("SysErrSlope_K0s_delta_eta_delta_phi_hist_%i", delta_eta_bin), Form("SysErrSlope_K0s_delta_eta_delta_phi_hist_%i", delta_eta_bin), 1, 0, 1);
    ResidualPolarization_delta_eta_delta_phi_hist[delta_eta_bin] = new TH1F(Form("ResidualPolarization_K0s_delta_eta_delta_phi_hist_%i", delta_eta_bin), Form("ResidualPolarization_K0s_delta_eta_delta_phi_hist_%i", delta_eta_bin), 1, 0, 1);
  }

  //----------------------------------------------------------------

  TFile *out_file_polarization = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/output/polarization/Polarization_K0s_PYTHIA.root", "recreate");

  //histograms

  //true MC
  TH1D *K0s_K0s_cosThetaProdPlane_hist = (TH1D*)inFile->Get("K0s_K0s_cosThetaProdPlane");

  TH1D *K0s_K0s_cosThetaProdPlane_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1D *K0s_K0s_cosThetaProdPlane_eta_hist[nEtaBins][nEtaBins];

  TH2F *K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_hist = (TH2F*)inFile->Get("K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_hist");
  TH2F *K0s_K0s_cos_theta_star_vs_delta_eta_hist = (TH2F*)inFile->Get("K0s_K0s_cos_theta_star_vs_delta_eta_hist");
  TH2F *K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_hist = (TH2F*)inFile->Get("K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_hist");


  TH1D *K0s_K0s_cosThetaProdPlane_cuts_hist = (TH1D*)inFile->Get("K0s_K0s_cosThetaProdPlane_cuts");

  TH1D *K0s_K0s_cosThetaProdPlane_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1D *K0s_K0s_cosThetaProdPlane_eta_cuts_hist[nEtaBins][nEtaBins];

  TH2F *K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_cuts_hist = (TH2F*)inFile->Get("K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_cuts_hist");
  TH2F *K0s_K0s_cos_theta_star_vs_delta_eta_cuts_hist = (TH2F*)inFile->Get("K0s_K0s_cos_theta_star_vs_delta_eta_cuts_hist");
  TH2F *K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_cuts_hist = (TH2F*)inFile->Get("K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_cuts_hist");

  //mixed event
  TH1D *K0s_K0s_cosThetaProdPlane_ME_hist = (TH1D*)inFile->Get("K0s_K0s_cosThetaProdPlane_ME");

  TH1D *K0s_K0s_cosThetaProdPlane_ME_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1D *K0s_K0s_cosThetaProdPlane_ME_eta_hist[nEtaBins][nEtaBins];


  TH1D *K0s_K0s_cosThetaProdPlane_ME_cuts_hist = (TH1D*)inFile->Get("K0s_K0s_cosThetaProdPlane_ME_cuts");

  TH1D *K0s_K0s_cosThetaProdPlane_ME_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1D *K0s_K0s_cosThetaProdPlane_ME_eta_cuts_hist[nEtaBins][nEtaBins];

  TH2F *K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_ME_cuts_hist = (TH2F*)inFile->Get("K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_ME_cuts_hist");
  TH2F *K0s_K0s_cos_theta_star_vs_delta_eta_ME_cuts_hist = (TH2F*)inFile->Get("K0s_K0s_cos_theta_star_vs_delta_eta_ME_cuts_hist");
  TH2F *K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_ME_cuts_hist = (TH2F*)inFile->Get("K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_ME_cuts_hist");

  //_________________________________________________________________________________________________________________________


  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  //true MC
  //before cuts
  TCanvas *K0s_K0s_cosThetaProdPlane_can = new TCanvas("K0s_K0s_cosThetaProdPlane_can", "K0s_K0s_cosThetaProdPlane_can", 1200, 1000);

  //signal+background
  float nK0sK0s = K0s_K0s_cosThetaProdPlane_hist->Integral();

  K0s_K0s_cosThetaProdPlane_ME_hist->Sumw2();
  K0s_K0s_cosThetaProdPlane_ME_hist->Scale(nK0sK0s/K0s_K0s_cosThetaProdPlane_ME_hist->Integral());
  K0s_K0s_cosThetaProdPlane_ME_hist->Scale(1./K0s_K0s_cosThetaProdPlane_ME_hist->GetXaxis()->GetBinWidth(1));

  K0s_K0s_cosThetaProdPlane_hist->SetMarkerStyle(20);
  K0s_K0s_cosThetaProdPlane_hist->SetMarkerColor(kRed);
  K0s_K0s_cosThetaProdPlane_hist->SetLineColor(kRed);
  K0s_K0s_cosThetaProdPlane_hist->GetXaxis()->SetTitle("cos(#theta*)");
  K0s_K0s_cosThetaProdPlane_hist->GetXaxis()->CenterTitle();
  K0s_K0s_cosThetaProdPlane_hist->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  K0s_K0s_cosThetaProdPlane_hist->GetYaxis()->CenterTitle();
  K0s_K0s_cosThetaProdPlane_hist->Sumw2();
  //K0s_K0s_cosThetaProdPlane_hist->Divide(K0s_K0s_cosThetaProdPlane_ME_hist);
  K0s_K0s_cosThetaProdPlane_hist->Scale(nK0sK0s/K0s_K0s_cosThetaProdPlane_hist->Integral());
  K0s_K0s_cosThetaProdPlane_hist->Scale(1./K0s_K0s_cosThetaProdPlane_hist->GetXaxis()->GetBinWidth(1));
  K0s_K0s_cosThetaProdPlane_hist->SetMinimum(0);
  K0s_K0s_cosThetaProdPlane_hist->Draw("p e");

  K0s_K0s_cosThetaProdPlane_ME_hist->Draw("p e same");


  TPaveText *K0s_K0s_MC_text = new TPaveText(0.4, 0.2, 0.8, 0.4, "NDC");
  //K0s_K0s_MC_text->AddText("STAR preliminary");
  //((TText*)K0s_K0s_MC_text->GetListOfLines()->Last())->SetTextColor(2);
  K0s_K0s_MC_text->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  K0s_K0s_MC_text->AddText("Minimum bias");
  K0s_K0s_MC_text->AddText("K_{s}^{0}-K_{s}^{0}");
  K0s_K0s_MC_text->AddText("True MC");
  K0s_K0s_MC_text->AddText("No cuts");
  K0s_K0s_MC_text->SetFillColorAlpha(0, 0.01);
  K0s_K0s_MC_text->Draw("same");

  K0s_K0s_cosThetaProdPlane_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/K0s/K0s_K0s_cosThetaProdPlane.png");

  //----------------------------------------------------------------------------------------------------------

  //_____________________________________________________________________________________________________________________________________________________________________________________

  //after analysis cuts
  TCanvas *K0s_K0s_cosThetaProdPlane_cuts_can = new TCanvas("K0s_K0s_cosThetaProdPlane_cuts_can", "K0s_K0s_cosThetaProdPlane_cuts_can", 1200, 1000);

  float nK0sK0s_cuts = K0s_K0s_cosThetaProdPlane_cuts_hist->Integral();

  K0s_K0s_cosThetaProdPlane_ME_cuts_hist->Sumw2();
  K0s_K0s_cosThetaProdPlane_ME_cuts_hist->Scale(nK0sK0s_cuts/K0s_K0s_cosThetaProdPlane_ME_cuts_hist->Integral());
  K0s_K0s_cosThetaProdPlane_ME_cuts_hist->Scale(1./K0s_K0s_cosThetaProdPlane_ME_cuts_hist->GetXaxis()->GetBinWidth(1));

  K0s_K0s_cosThetaProdPlane_cuts_hist->SetMarkerStyle(20);
  K0s_K0s_cosThetaProdPlane_cuts_hist->SetMarkerColor(kRed);
  K0s_K0s_cosThetaProdPlane_cuts_hist->SetLineColor(kRed);
  K0s_K0s_cosThetaProdPlane_cuts_hist->GetXaxis()->SetTitle("cos(#theta*)");
  K0s_K0s_cosThetaProdPlane_cuts_hist->GetXaxis()->CenterTitle();
  K0s_K0s_cosThetaProdPlane_cuts_hist->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  K0s_K0s_cosThetaProdPlane_cuts_hist->GetYaxis()->CenterTitle();
  K0s_K0s_cosThetaProdPlane_cuts_hist->Sumw2();
  //K0s_K0s_cosThetaProdPlane_cuts_hist->Divide(K0s_K0s_cosThetaProdPlane_ME_cuts_hist); //correct using ME
  K0s_K0s_cosThetaProdPlane_cuts_hist->Scale(nK0sK0s_cuts/K0s_K0s_cosThetaProdPlane_cuts_hist->Integral()); //scale back
  K0s_K0s_cosThetaProdPlane_cuts_hist->Scale(1./K0s_K0s_cosThetaProdPlane_cuts_hist->GetXaxis()->GetBinWidth(1)); //bin width
  K0s_K0s_cosThetaProdPlane_cuts_hist->SetMinimum(0);
  K0s_K0s_cosThetaProdPlane_cuts_hist->Draw("p e");

  K0s_K0s_cosThetaProdPlane_ME_cuts_hist->Draw("p e same");

  TPaveText *K0s_K0s_MC_text_cuts = new TPaveText(0.4, 0.2, 0.8, 0.4, "NDC");
  //K0s_K0s_MC_text_cuts->AddText("STAR preliminary");
  //((TText*)K0s_K0s_MC_text_cuts->GetListOfLines()->Last())->SetTextColor(2);
  K0s_K0s_MC_text_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  K0s_K0s_MC_text_cuts->AddText("Minimum bias");
  K0s_K0s_MC_text_cuts->AddText("K_{s}^{0}-K_{s}^{0}");
  K0s_K0s_MC_text_cuts->AddText("True MC");
  K0s_K0s_MC_text_cuts->AddText("Analysis cuts");
  K0s_K0s_MC_text_cuts->SetFillColorAlpha(0, 0.01);
  K0s_K0s_MC_text_cuts->Draw("same");

  K0s_K0s_cosThetaProdPlane_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/K0s/K0s_K0s_cosThetaProdPlane_cuts.png");

  //----------------------------------------------------------------------------------------------------------

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
    PolarizationGraph_delta_eta_PYTHIA[delta_eta_bin-1] = new TGraphErrors(1);
    PolarizationGraph_delta_phi_PYTHIA[delta_eta_bin-1] = new TGraphErrors(1);
    PolarizationGraph_delta_eta_delta_phi_PYTHIA[delta_eta_bin-1] = new TGraphErrors(1);

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

    //delta eta, before cuts
    TCanvas *K0s_K0s_cosThetaProdPlane_delta_eta_can = new TCanvas(Form("K0s_K0s_cosThetaProdPlane_delta_eta_can_%i", delta_eta_bin), Form("K0s_K0s_cosThetaProdPlane_delta_eta_can_%i", delta_eta_bin), 1200, 1000);
    K0s_K0s_cosThetaProdPlane_delta_eta_can->cd();

    TF1 *fit_K0s_K0s_after_delta_eta = new TF1("fit_K0s_K0s_after_delta_eta", "[0]*(1+[1]*x)", -1, 1);
    fit_K0s_K0s_after_delta_eta->SetParameters(1000, 0.10);

    TH1D *K0s_K0s_cosThetaProdPlane_delta_eta = K0s_K0s_cos_theta_star_vs_delta_eta_hist->ProjectionX(Form("proj_eta_SE_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
    float nK0sK0s_delta_eta = K0s_K0s_cosThetaProdPlane_delta_eta->Integral();
    K0s_K0s_cosThetaProdPlane_delta_eta->Sumw2();
    K0s_K0s_cosThetaProdPlane_delta_eta->SetMarkerStyle(20);
    K0s_K0s_cosThetaProdPlane_delta_eta->SetMarkerColor(kRed);
    K0s_K0s_cosThetaProdPlane_delta_eta->SetLineColor(kRed);
    K0s_K0s_cosThetaProdPlane_delta_eta->GetXaxis()->SetTitle("cos(#theta*)");
    K0s_K0s_cosThetaProdPlane_delta_eta->GetXaxis()->CenterTitle();
    K0s_K0s_cosThetaProdPlane_delta_eta->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    K0s_K0s_cosThetaProdPlane_delta_eta->GetYaxis()->CenterTitle();
    //K0s_K0s_cosThetaProdPlane_delta_eta->Divide(K0s_K0s_cosThetaProdPlane_ME_delta_eta); //correct using ME
    //K0s_K0s_cosThetaProdPlane_delta_eta->Scale(nK0sK0s_delta_eta/K0s_K0s_cosThetaProdPlane_delta_eta->Integral()); //scale back
    K0s_K0s_cosThetaProdPlane_delta_eta->Scale(1./K0s_K0s_cosThetaProdPlane_delta_eta->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
    K0s_K0s_cosThetaProdPlane_delta_eta->Fit(fit_K0s_K0s_after_delta_eta, "i 0 r");
    K0s_K0s_cosThetaProdPlane_delta_eta->SetMinimum(0);
    //K0s_K0s_cosThetaProdPlane_delta_eta->GetYaxis()->SetRangeUser(-20e9, 0);
    K0s_K0s_cosThetaProdPlane_delta_eta->Draw("p e");

    fit_K0s_K0s_after_delta_eta->SetLineColor(1);
    fit_K0s_K0s_after_delta_eta->Draw("same");

    float K0s_K0s_slope_delta_eta = fit_K0s_K0s_after_delta_eta->GetParameter(1)/(L0_alpha*L0_alpha);
    float K0s_K0s_slope_delta_eta_err = fit_K0s_K0s_after_delta_eta->GetParError(1)/(L0_alpha*L0_alpha);

    TPaveText *K0s_K0s_text_MC_delta_eta = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    K0s_K0s_text_MC_delta_eta->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    K0s_K0s_text_MC_delta_eta->AddText("Minimum bias");
    K0s_K0s_text_MC_delta_eta->AddText("K_{s}^{0}-K_{s}^{0}");
    K0s_K0s_text_MC_delta_eta->AddText("True MC");
    K0s_K0s_text_MC_delta_eta->AddText("Analysis cuts");
    K0s_K0s_text_MC_delta_eta->AddText(Form("P = %.3f #pm %.3f", K0s_K0s_slope_delta_eta, K0s_K0s_slope_delta_eta_err ));
    K0s_K0s_text_MC_delta_eta->SetFillColorAlpha(0, 0.01);
    K0s_K0s_text_MC_delta_eta->Draw("same");

    PolarizationGraph_delta_eta_PYTHIA[delta_eta_bin-1]->SetPoint(1, K0s_K0s_slope_delta_eta, 4.15);
    PolarizationGraph_delta_eta_PYTHIA[delta_eta_bin-1]->SetPointError(1, fabs(K0s_K0s_slope_delta_eta_err), 0);

    K0s_K0s_cosThetaProdPlane_delta_eta_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/K0s/K0s_K0s_cosThetaProdPlane_delta_eta_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------

    //delta phi, before cuts
    TCanvas *K0s_K0s_cosThetaProdPlane_delta_phi_can = new TCanvas(Form("K0s_K0s_cosThetaProdPlane_delta_phi_can_%i", delta_eta_bin), Form("K0s_K0s_cosThetaProdPlane_delta_phi_can_%i", delta_eta_bin), 1200, 1000);
    K0s_K0s_cosThetaProdPlane_delta_phi_can->cd();

    TF1 *fit_K0s_K0s_after_delta_phi = new TF1("fit_K0s_K0s_after_delta_phi", "[0]*(1+[1]*x)", -1, 1);
    fit_K0s_K0s_after_delta_phi->SetParameters(1000, 0.10);

    TH1D *K0s_K0s_cosThetaProdPlane_delta_phi = K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_hist->ProjectionX(Form("proj_phi_SE_%i", delta_eta_bin), bin_min_loop, bin_max_loop);
    float nK0sK0s_delta_phi = K0s_K0s_cosThetaProdPlane_delta_phi->Integral();
    K0s_K0s_cosThetaProdPlane_delta_phi->Sumw2();
    K0s_K0s_cosThetaProdPlane_delta_phi->SetMarkerStyle(20);
    K0s_K0s_cosThetaProdPlane_delta_phi->SetMarkerColor(kRed);
    K0s_K0s_cosThetaProdPlane_delta_phi->SetLineColor(kRed);
    K0s_K0s_cosThetaProdPlane_delta_phi->GetXaxis()->SetTitle("cos(#theta*)");
    K0s_K0s_cosThetaProdPlane_delta_phi->GetXaxis()->CenterTitle();
    K0s_K0s_cosThetaProdPlane_delta_phi->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    K0s_K0s_cosThetaProdPlane_delta_phi->GetYaxis()->CenterTitle();
    //K0s_K0s_cosThetaProdPlane_delta_phi->Divide(K0s_K0s_cosThetaProdPlane_ME_delta_phi); //correct using ME
    //K0s_K0s_cosThetaProdPlane_delta_phi->Scale(nK0sK0s_delta_phi/K0s_K0s_cosThetaProdPlane_delta_phi->Integral()); //scale back
    K0s_K0s_cosThetaProdPlane_delta_phi->Scale(1./K0s_K0s_cosThetaProdPlane_delta_phi->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
    K0s_K0s_cosThetaProdPlane_delta_phi->Fit(fit_K0s_K0s_after_delta_phi, "i 0 r");
    K0s_K0s_cosThetaProdPlane_delta_phi->SetMinimum(0);
    //K0s_K0s_cosThetaProdPlane_delta_phi->GetYaxis()->SetRangeUser(-20e9, 0);
    K0s_K0s_cosThetaProdPlane_delta_phi->Draw("p e");

    fit_K0s_K0s_after_delta_phi->SetLineColor(1);
    fit_K0s_K0s_after_delta_phi->Draw("same");

    float K0s_K0s_slope_delta_phi = fit_K0s_K0s_after_delta_phi->GetParameter(1)/(L0_alpha*L0_alpha);
    float K0s_K0s_slope_delta_phi_err = fit_K0s_K0s_after_delta_phi->GetParError(1)/(L0_alpha*L0_alpha);

    TPaveText *K0s_K0s_text_MC_delta_phi = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    K0s_K0s_text_MC_delta_phi->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    K0s_K0s_text_MC_delta_phi->AddText("Minimum bias");
    K0s_K0s_text_MC_delta_phi->AddText("K_{s}^{0}-K_{s}^{0}");
    K0s_K0s_text_MC_delta_phi->AddText("True MC");
    K0s_K0s_text_MC_delta_phi->AddText("Analysis cuts");
    K0s_K0s_text_MC_delta_phi->AddText(Form("P = %.3f #pm %.3f", K0s_K0s_slope_delta_phi, K0s_K0s_slope_delta_phi_err ));
    K0s_K0s_text_MC_delta_phi->SetFillColorAlpha(0, 0.01);
    K0s_K0s_text_MC_delta_phi->Draw("same");

    PolarizationGraph_delta_phi_PYTHIA[delta_eta_bin-1]->SetPoint(1, K0s_K0s_slope_delta_phi, 4.15);
    PolarizationGraph_delta_phi_PYTHIA[delta_eta_bin-1]->SetPointError(1, fabs(K0s_K0s_slope_delta_phi_err), 0);

    K0s_K0s_cosThetaProdPlane_delta_phi_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/K0s/K0s_K0s_cosThetaProdPlane_delta_phi_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------


    //delta eta delta phi, before cuts
    TCanvas *K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_can = new TCanvas(Form("K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_can_%i", delta_eta_bin), Form("K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_can_%i", delta_eta_bin), 1200, 1000);
    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_can->cd();

    TF1 *fit_K0s_K0s_after_delta_eta_delta_phi = new TF1("fit_K0s_K0s_after_delta_eta_delta_phi", "[0]*(1+[1]*x)", -1, 1);
    fit_K0s_K0s_after_delta_eta_delta_phi->SetParameters(1000, 0.10);



    TH1D *K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi = K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_hist->ProjectionX(Form("proj_eta_SE_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
    float nK0sK0s_delta_eta_delta_phi = K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi->Integral();
    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi->Sumw2();
    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi->SetMarkerStyle(20);
    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi->SetMarkerSize(1.5);
    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi->SetMarkerColor(kRed);
    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi->SetLineColor(kRed);
    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi->GetXaxis()->SetTitle("cos(#theta*)");
    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi->GetXaxis()->CenterTitle();
    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi->GetYaxis()->CenterTitle();
    //K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi->Divide(K0s_K0s_cosThetaProdPlane_ME_delta_eta_delta_phi); //correct using ME
    //K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi->Scale(nK0sK0s_delta_eta_delta_phi/K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi->Integral()); //scale back
    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi->Scale(1./K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi->Fit(fit_K0s_K0s_after_delta_eta_delta_phi, "i 0 r");
    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi->SetMinimum(0);
    //K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi->GetYaxis()->SetRangeUser(-20e9, 0);
    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi->Draw("p e");

    fit_K0s_K0s_after_delta_eta_delta_phi->SetLineColor(1);
    fit_K0s_K0s_after_delta_eta_delta_phi->Draw("same");

    float K0s_K0s_slope_delta_eta_delta_phi = fit_K0s_K0s_after_delta_eta_delta_phi->GetParameter(1)/(L0_alpha*L0_alpha);
    float K0s_K0s_slope_delta_eta_delta_phi_err = fit_K0s_K0s_after_delta_eta_delta_phi->GetParError(1)/(L0_alpha*L0_alpha);

    TPaveText *K0s_K0s_text_MC_delta_eta_delta_phi = new TPaveText(0.5, 0.2, 0.8, 0.4, "NDC");
    K0s_K0s_text_MC_delta_eta_delta_phi->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    K0s_K0s_text_MC_delta_eta_delta_phi->AddText("Minimum bias");
    K0s_K0s_text_MC_delta_eta_delta_phi->AddText("K_{s}^{0}-K_{s}^{0}");
    K0s_K0s_text_MC_delta_eta_delta_phi->AddText("True MC");
    K0s_K0s_text_MC_delta_eta_delta_phi->AddText("Analysis cuts");
    K0s_K0s_text_MC_delta_eta_delta_phi->AddText(Form("P = %.3f #pm %.3f", K0s_K0s_slope_delta_eta_delta_phi, K0s_K0s_slope_delta_eta_delta_phi_err ));
    if( (delta_eta_bin-1) == 0 )K0s_K0s_text_MC_delta_eta_delta_phi->AddText("|#Delta#it{y}| < 0.5, |#Delta#phi| < #pi/3");
    if( (delta_eta_bin-1) == 1 )K0s_K0s_text_MC_delta_eta_delta_phi->AddText("0.5 < |#Delta#it{y}| < 2.0 or |#Delta#phi| > #pi/3");
    K0s_K0s_text_MC_delta_eta_delta_phi->SetFillColorAlpha(0, 0.01);
    K0s_K0s_text_MC_delta_eta_delta_phi->Draw("same");

    PolarizationGraph_delta_eta_delta_phi_PYTHIA[delta_eta_bin-1]->SetPoint(1, K0s_K0s_slope_delta_eta_delta_phi, 4.15);
    PolarizationGraph_delta_eta_delta_phi_PYTHIA[delta_eta_bin-1]->SetPointError(1, fabs(K0s_K0s_slope_delta_eta_delta_phi_err), 0);

    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/K0s/K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_%i.png", delta_eta_bin));

    //--------------------------------------------------------------------------------

    PolarizationGraph_delta_eta_PYTHIA[delta_eta_bin-1]->Write(Form("PolarizationGraph_K0s_delta_eta_PYTHIA_%i", delta_eta_bin-1));
    PolarizationGraph_delta_phi_PYTHIA[delta_eta_bin-1]->Write(Form("PolarizationGraph_K0s_delta_phi_PYTHIA_%i", delta_eta_bin-1));
    PolarizationGraph_delta_eta_delta_phi_PYTHIA[delta_eta_bin-1]->Write(Form("PolarizationGraph_K0s_delta_eta_delta_phi_PYTHIA_%i", delta_eta_bin-1));

    //________________________________________________________________________________________________________________________________________________________________________________________________________________________________


    sysErrFile->cd();

    //delta eta, after cuts
    TCanvas *K0s_K0s_cosThetaProdPlane_delta_eta_cuts_can = new TCanvas(Form("K0s_K0s_cosThetaProdPlane_delta_eta_cuts_can_%i", delta_eta_bin), Form("K0s_K0s_cosThetaProdPlane_delta_eta_cuts_can_%i", delta_eta_bin), 1200, 1000);
    K0s_K0s_cosThetaProdPlane_delta_eta_cuts_can->cd();

    TF1 *fit_K0s_K0s_after_delta_eta_cuts = new TF1("fit_K0s_K0s_after_delta_eta_cuts", "[0]*(1+[1]*x)", -1, 1);
    fit_K0s_K0s_after_delta_eta_cuts->SetParameters(1000, 0.10);

    TF1 *fit_K0s_K0s_after_delta_eta_cuts_ME = new TF1("fit_K0s_K0s_after_delta_eta_cuts_ME", "[0]*(1+[1]*x)", -1, 1);
    fit_K0s_K0s_after_delta_eta_cuts_ME->SetParameters(1000, 0.10);

    TH1D *K0s_K0s_cosThetaProdPlane_delta_eta_cuts = K0s_K0s_cos_theta_star_vs_delta_eta_cuts_hist->ProjectionX(Form("proj_eta_SE_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
    float nK0sK0s_delta_eta_cuts = K0s_K0s_cosThetaProdPlane_delta_eta_cuts->Integral();

    TH1D *K0s_K0s_cosThetaProdPlane_ME_delta_eta_cuts = K0s_K0s_cos_theta_star_vs_delta_eta_ME_cuts_hist->ProjectionX(Form("proj_eta_ME_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
    K0s_K0s_cosThetaProdPlane_ME_delta_eta_cuts->Sumw2();
    K0s_K0s_cosThetaProdPlane_ME_delta_eta_cuts->SetMarkerStyle(24);
    K0s_K0s_cosThetaProdPlane_ME_delta_eta_cuts->SetMarkerSize(1.5);
    K0s_K0s_cosThetaProdPlane_ME_delta_eta_cuts->SetMarkerColor(1);
    K0s_K0s_cosThetaProdPlane_ME_delta_eta_cuts->SetLineColor(1);
    K0s_K0s_cosThetaProdPlane_ME_delta_eta_cuts->Scale(nK0sK0s_delta_eta_cuts/K0s_K0s_cosThetaProdPlane_ME_delta_eta_cuts->Integral());
    K0s_K0s_cosThetaProdPlane_ME_delta_eta_cuts->Scale(1./K0s_K0s_cosThetaProdPlane_ME_delta_eta_cuts->GetXaxis()->GetBinWidth(1));
    K0s_K0s_cosThetaProdPlane_ME_delta_eta_cuts->Fit(fit_K0s_K0s_after_delta_eta_cuts_ME, "i 0 r");
    //K0s_K0s_cosThetaProdPlane_ME->Draw("p e");


    K0s_K0s_cosThetaProdPlane_delta_eta_cuts->Sumw2();
    K0s_K0s_cosThetaProdPlane_delta_eta_cuts->SetMarkerStyle(20);
    K0s_K0s_cosThetaProdPlane_delta_eta_cuts->SetMarkerColor(kRed);
    K0s_K0s_cosThetaProdPlane_delta_eta_cuts->SetLineColor(kRed);
    K0s_K0s_cosThetaProdPlane_delta_eta_cuts->GetXaxis()->SetTitle("cos(#theta*)");
    K0s_K0s_cosThetaProdPlane_delta_eta_cuts->GetXaxis()->CenterTitle();
    K0s_K0s_cosThetaProdPlane_delta_eta_cuts->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    K0s_K0s_cosThetaProdPlane_delta_eta_cuts->GetYaxis()->CenterTitle();
    //K0s_K0s_cosThetaProdPlane_delta_eta_cuts->Divide(K0s_K0s_cosThetaProdPlane_ME_delta_eta_cuts); //correct using ME
    //K0s_K0s_cosThetaProdPlane_delta_eta_cuts->Scale(nK0sK0s_delta_eta_cuts/K0s_K0s_cosThetaProdPlane_delta_eta_cuts->Integral()); //scale back
    K0s_K0s_cosThetaProdPlane_delta_eta_cuts->Scale(1./K0s_K0s_cosThetaProdPlane_delta_eta_cuts->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
    K0s_K0s_cosThetaProdPlane_delta_eta_cuts->Fit(fit_K0s_K0s_after_delta_eta_cuts, "i 0 r");
    K0s_K0s_cosThetaProdPlane_delta_eta_cuts->SetMinimum(0);
    //K0s_K0s_cosThetaProdPlane_delta_eta_cuts->GetYaxis()->SetRangeUser(-20e9, 0);
    K0s_K0s_cosThetaProdPlane_delta_eta_cuts->Draw("p e");

    K0s_K0s_cosThetaProdPlane_ME_delta_eta_cuts->Draw("p e same");

    fit_K0s_K0s_after_delta_eta_cuts->SetLineColor(1);
    fit_K0s_K0s_after_delta_eta_cuts->Draw("same");

    fit_K0s_K0s_after_delta_eta_cuts_ME->SetLineColor(kBlue);
    fit_K0s_K0s_after_delta_eta_cuts_ME->Draw("same");


    float K0s_K0s_slope_delta_eta_cuts = fit_K0s_K0s_after_delta_eta_cuts->GetParameter(1)/(L0_alpha*L0_alpha);
    float K0s_K0s_slope_delta_eta_cuts_err = fit_K0s_K0s_after_delta_eta_cuts->GetParError(1)/(L0_alpha*L0_alpha);


    //statistical error correction for systematic error
    float SysErrSlope_K0s_K0s_delta_eta_corr = sqrt( fabs( fit_K0s_K0s_after_delta_eta_cuts->GetParError(1)*fit_K0s_K0s_after_delta_eta_cuts->GetParError(1) - fit_K0s_K0s_after_delta_eta_cuts_ME->GetParError(1)*fit_K0s_K0s_after_delta_eta_cuts_ME->GetParError(1) ) );

    float SysErrSlope_K0s_K0s_delta_eta_work = ( fabs( fit_K0s_K0s_after_delta_eta_cuts->GetParameter(1) - fit_K0s_K0s_after_delta_eta_cuts_ME->GetParameter(1) ) - SysErrSlope_K0s_K0s_delta_eta_corr )/fabs(fit_K0s_K0s_after_delta_eta_cuts->GetParameter(1));

    float K0s_K0s_slope_delta_eta_cuts_sys_err = 0;

    if( SysErrSlope_K0s_K0s_delta_eta_work > 0 ) K0s_K0s_slope_delta_eta_cuts_sys_err = SysErrSlope_K0s_K0s_delta_eta_work; //store sys. err. only if it is larger than statistical fluctuations

    float K0s_K0s_slope_delta_eta_cuts_sys_err_no_corr = fabs( fit_K0s_K0s_after_delta_eta_cuts->GetParameter(1) - fit_K0s_K0s_after_delta_eta_cuts_ME->GetParameter(1) )/fabs(fit_K0s_K0s_after_delta_eta_cuts->GetParameter(1));

    if(corr_err == 1) SysErrSlope_delta_eta_hist[delta_eta_bin-1]->SetBinContent(1, K0s_K0s_slope_delta_eta_cuts_sys_err);
    if(corr_err == 0) SysErrSlope_delta_eta_hist[delta_eta_bin-1]->SetBinContent(1, K0s_K0s_slope_delta_eta_cuts_sys_err_no_corr);

    TPaveText *K0s_K0s_text_MC_delta_eta_cuts = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    K0s_K0s_text_MC_delta_eta_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    K0s_K0s_text_MC_delta_eta_cuts->AddText("Minimum bias");
    K0s_K0s_text_MC_delta_eta_cuts->AddText("K_{s}^{0}-K_{s}^{0}");
    K0s_K0s_text_MC_delta_eta_cuts->AddText("True MC");
    K0s_K0s_text_MC_delta_eta_cuts->AddText("Analysis cuts");
    K0s_K0s_text_MC_delta_eta_cuts->AddText(Form("P = %.3f #pm %.3f", K0s_K0s_slope_delta_eta_cuts, K0s_K0s_slope_delta_eta_cuts_err ));
    K0s_K0s_text_MC_delta_eta_cuts->AddText(Form("P_{ME} = %.3f", fit_K0s_K0s_after_delta_eta_cuts_ME->GetParameter(1)/(L0_alpha*L0_alpha) ));
    //if(corr_err == 1) K0s_K0s_text_MC_delta_eta_cuts->AddText(Form("Slope sys. err: %.1f %%", K0s_K0s_slope_delta_eta_cuts_sys_err*100 ));
    //if(corr_err == 0) K0s_K0s_text_MC_delta_eta_cuts->AddText(Form("Slope sys. err: %.1f %%", K0s_K0s_slope_delta_eta_cuts_sys_err_no_corr*100 ));
    K0s_K0s_text_MC_delta_eta_cuts->SetFillColorAlpha(0, 0.01);
    K0s_K0s_text_MC_delta_eta_cuts->Draw("same");

    K0s_K0s_cosThetaProdPlane_delta_eta_cuts_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/K0s/K0s_K0s_cosThetaProdPlane_delta_eta_cuts_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------

    TCanvas *K0s_K0s_cosThetaProdPlane_delta_eta_cuts_corr_can = new TCanvas(Form("K0s_K0s_cosThetaProdPlane_delta_eta_cuts_corr_can_%i", delta_eta_bin), Form("K0s_K0s_cosThetaProdPlane_delta_eta_cuts_corr_can_%i", delta_eta_bin), 1200, 1000);
    K0s_K0s_cosThetaProdPlane_delta_eta_cuts_corr_can->cd();

    TF1 *fit_K0s_K0s_after_delta_eta_cuts_corr = new TF1("fit_K0s_K0s_after_delta_eta_cuts_corr", "[0]*(1+[1]*x)", -1, 1);
    fit_K0s_K0s_after_delta_eta_cuts_corr->SetParameters(1000, 0.10);

    TH1D *K0s_K0s_cosThetaProdPlane_delta_eta_cuts_corr = (TH1D*)K0s_K0s_cosThetaProdPlane_delta_eta_cuts->Clone(Form("K0s_K0s_cosThetaProdPlane_delta_eta_cuts_corr_%i", delta_eta_bin));
    K0s_K0s_cosThetaProdPlane_delta_eta_cuts_corr->Divide(K0s_K0s_cosThetaProdPlane_ME_delta_eta_cuts); //correct
    K0s_K0s_cosThetaProdPlane_delta_eta_cuts_corr->Scale(nK0sK0s_delta_eta_cuts/K0s_K0s_cosThetaProdPlane_delta_eta_cuts_corr->Integral()); //scale back
    K0s_K0s_cosThetaProdPlane_delta_eta_cuts_corr->Scale(1./K0s_K0s_cosThetaProdPlane_delta_eta_cuts_corr->GetXaxis()->GetBinWidth(1)); //divide by bin width
    //K0s_K0s_cosThetaProdPlane_delta_eta_cuts_corr->Scale(nK0sK0s_delta_eta_cuts);
    K0s_K0s_cosThetaProdPlane_delta_eta_cuts_corr->Fit(fit_K0s_K0s_after_delta_eta_cuts_corr, "i 0 r");
    K0s_K0s_cosThetaProdPlane_delta_eta_cuts_corr->SetMinimum(0);
    K0s_K0s_cosThetaProdPlane_delta_eta_cuts_corr->Draw("p e");

    fit_K0s_K0s_after_delta_eta_cuts_corr->SetLineColor(1);
    fit_K0s_K0s_after_delta_eta_cuts_corr->Draw("same");

    float K0s_K0s_slope_delta_eta_cuts_corr = fit_K0s_K0s_after_delta_eta_cuts_corr->GetParameter(1)/(L0_alpha*L0_alpha);
    float K0s_K0s_slope_delta_eta_cuts_corr_err = fit_K0s_K0s_after_delta_eta_cuts_corr->GetParError(1)/(L0_alpha*L0_alpha);

    ResidualPolarization_delta_eta_hist[delta_eta_bin-1]->SetBinContent(1, K0s_K0s_slope_delta_eta_cuts_corr);
    ResidualPolarization_delta_eta_hist[delta_eta_bin-1]->SetBinError(1, K0s_K0s_slope_delta_eta_cuts_corr_err);


    TPaveText *K0s_K0s_text_MC_delta_eta_cuts_corr = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    K0s_K0s_text_MC_delta_eta_cuts_corr->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    K0s_K0s_text_MC_delta_eta_cuts_corr->AddText("Minimum bias");
    K0s_K0s_text_MC_delta_eta_cuts_corr->AddText("K_{s}^{0}-K_{s}^{0}");
    K0s_K0s_text_MC_delta_eta_cuts_corr->AddText("True MC");
    K0s_K0s_text_MC_delta_eta_cuts_corr->AddText("Analysis cuts");
    K0s_K0s_text_MC_delta_eta_cuts_corr->AddText(Form("P = %.3f #pm %.3f", K0s_K0s_slope_delta_eta_cuts_corr, K0s_K0s_slope_delta_eta_cuts_corr_err ));
    //if(corr_err == 1) K0s_K0s_text_MC_delta_eta_cuts_corr->AddText(Form("Slope sys. err: %.1f %%", K0s_K0s_slope_delta_eta_cuts_corr_sys_err*100 ));
    //if(corr_err == 0) K0s_K0s_text_MC_delta_eta_cuts_corr->AddText(Form("Slope sys. err: %.1f %%", K0s_K0s_slope_delta_eta_cuts_corr_sys_err_no_corr*100 ));
    K0s_K0s_text_MC_delta_eta_cuts_corr->SetFillColorAlpha(0, 0.01);
    K0s_K0s_text_MC_delta_eta_cuts_corr->Draw("same");


    K0s_K0s_cosThetaProdPlane_delta_eta_cuts_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/K0s/K0s_K0s_cosThetaProdPlane_delta_eta_cuts_corr_%i.png", delta_eta_bin));

    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    //delta phi, after cuts
    TCanvas *K0s_K0s_cosThetaProdPlane_delta_phi_cuts_can = new TCanvas(Form("K0s_K0s_cosThetaProdPlane_delta_phi_cuts_can_%i", delta_eta_bin), Form("K0s_K0s_cosThetaProdPlane_delta_phi_cuts_can_%i", delta_eta_bin), 1200, 1000);
    K0s_K0s_cosThetaProdPlane_delta_phi_cuts_can->cd();

    TF1 *fit_K0s_K0s_after_delta_phi_cuts = new TF1("fit_K0s_K0s_after_delta_phi_cuts", "[0]*(1+[1]*x)", -1, 1);
    fit_K0s_K0s_after_delta_phi_cuts->SetParameters(1000, 0.10);

    TF1 *fit_K0s_K0s_after_delta_phi_cuts_ME = new TF1("fit_K0s_K0s_after_delta_phi_cuts_ME", "[0]*(1+[1]*x)", -1, 1);
    fit_K0s_K0s_after_delta_phi_cuts_ME->SetParameters(1000, 0.10);

    TH1D *K0s_K0s_cosThetaProdPlane_delta_phi_cuts = K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_cuts_hist->ProjectionX(Form("proj_phi_SE_%i", delta_eta_bin), bin_min_loop, bin_max_loop);
    float nK0sK0s_delta_phi_cuts = K0s_K0s_cosThetaProdPlane_delta_phi_cuts->Integral();

    TH1D *K0s_K0s_cosThetaProdPlane_ME_delta_phi_cuts = K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_ME_cuts_hist->ProjectionX(Form("proj_phi_ME_%i", delta_eta_bin), bin_min_loop, bin_max_loop);
    K0s_K0s_cosThetaProdPlane_ME_delta_phi_cuts->Sumw2();
    K0s_K0s_cosThetaProdPlane_ME_delta_phi_cuts->SetMarkerStyle(24);
    K0s_K0s_cosThetaProdPlane_ME_delta_phi_cuts->SetMarkerSize(1.5);
    K0s_K0s_cosThetaProdPlane_ME_delta_phi_cuts->SetMarkerColor(1);
    K0s_K0s_cosThetaProdPlane_ME_delta_phi_cuts->SetLineColor(1);
    K0s_K0s_cosThetaProdPlane_ME_delta_phi_cuts->Scale(nK0sK0s_delta_phi_cuts/K0s_K0s_cosThetaProdPlane_ME_delta_phi_cuts->Integral());
    K0s_K0s_cosThetaProdPlane_ME_delta_phi_cuts->Scale(1./K0s_K0s_cosThetaProdPlane_ME_delta_phi_cuts->GetXaxis()->GetBinWidth(1));
    K0s_K0s_cosThetaProdPlane_ME_delta_phi_cuts->Fit(fit_K0s_K0s_after_delta_phi_cuts_ME, "i 0 r");
    //K0s_K0s_cosThetaProdPlane_ME->Draw("p e");


    K0s_K0s_cosThetaProdPlane_delta_phi_cuts->Sumw2();
    K0s_K0s_cosThetaProdPlane_delta_phi_cuts->SetMarkerStyle(20);
    K0s_K0s_cosThetaProdPlane_delta_phi_cuts->SetMarkerColor(kRed);
    K0s_K0s_cosThetaProdPlane_delta_phi_cuts->SetLineColor(kRed);
    K0s_K0s_cosThetaProdPlane_delta_phi_cuts->GetXaxis()->SetTitle("cos(#theta*)");
    K0s_K0s_cosThetaProdPlane_delta_phi_cuts->GetXaxis()->CenterTitle();
    K0s_K0s_cosThetaProdPlane_delta_phi_cuts->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    K0s_K0s_cosThetaProdPlane_delta_phi_cuts->GetYaxis()->CenterTitle();
    //K0s_K0s_cosThetaProdPlane_delta_phi_cuts->Divide(K0s_K0s_cosThetaProdPlane_ME_delta_phi_cuts); //correct using ME
    //K0s_K0s_cosThetaProdPlane_delta_phi_cuts->Scale(nK0sK0s_delta_phi_cuts/K0s_K0s_cosThetaProdPlane_delta_phi_cuts->Integral()); //scale back
    K0s_K0s_cosThetaProdPlane_delta_phi_cuts->Scale(1./K0s_K0s_cosThetaProdPlane_delta_phi_cuts->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
    K0s_K0s_cosThetaProdPlane_delta_phi_cuts->Fit(fit_K0s_K0s_after_delta_phi_cuts, "i 0 r");
    K0s_K0s_cosThetaProdPlane_delta_phi_cuts->SetMinimum(0);
    //K0s_K0s_cosThetaProdPlane_delta_phi_cuts->GetYaxis()->SetRangeUser(-20e9, 0);
    K0s_K0s_cosThetaProdPlane_delta_phi_cuts->Draw("p e");

    K0s_K0s_cosThetaProdPlane_ME_delta_phi_cuts->Draw("p e same");

    fit_K0s_K0s_after_delta_phi_cuts->SetLineColor(1);
    fit_K0s_K0s_after_delta_phi_cuts->Draw("same");

    fit_K0s_K0s_after_delta_phi_cuts_ME->SetLineColor(kBlue);
    fit_K0s_K0s_after_delta_phi_cuts_ME->Draw("same");


    float K0s_K0s_slope_delta_phi_cuts = fit_K0s_K0s_after_delta_phi_cuts->GetParameter(1)/(L0_alpha*L0_alpha);
    float K0s_K0s_slope_delta_phi_cuts_err = fit_K0s_K0s_after_delta_phi_cuts->GetParError(1)/(L0_alpha*L0_alpha);


    //statistical error correction for systematic error
    float SysErrSlope_K0s_K0s_delta_phi_corr = sqrt( fabs( fit_K0s_K0s_after_delta_phi_cuts->GetParError(1)*fit_K0s_K0s_after_delta_phi_cuts->GetParError(1) - fit_K0s_K0s_after_delta_phi_cuts_ME->GetParError(1)*fit_K0s_K0s_after_delta_phi_cuts_ME->GetParError(1) ) );

    float SysErrSlope_K0s_K0s_delta_phi_work = ( fabs( fit_K0s_K0s_after_delta_phi_cuts->GetParameter(1) - fit_K0s_K0s_after_delta_phi_cuts_ME->GetParameter(1) ) - SysErrSlope_K0s_K0s_delta_phi_corr )/fabs(fit_K0s_K0s_after_delta_phi_cuts->GetParameter(1));

    float K0s_K0s_slope_delta_phi_cuts_sys_err = 0;

    if( SysErrSlope_K0s_K0s_delta_phi_work > 0 ) K0s_K0s_slope_delta_phi_cuts_sys_err = SysErrSlope_K0s_K0s_delta_phi_work; //store sys. err. only if it is larger than statistical fluctuations

    float K0s_K0s_slope_delta_phi_cuts_sys_err_no_corr = fabs( fit_K0s_K0s_after_delta_phi_cuts->GetParameter(1) - fit_K0s_K0s_after_delta_phi_cuts_ME->GetParameter(1) )/fabs(fit_K0s_K0s_after_delta_phi_cuts->GetParameter(1));

    if(corr_err == 1) SysErrSlope_delta_phi_hist[delta_eta_bin-1]->SetBinContent(1, K0s_K0s_slope_delta_phi_cuts_sys_err);
    if(corr_err == 0) SysErrSlope_delta_phi_hist[delta_eta_bin-1]->SetBinContent(1, K0s_K0s_slope_delta_phi_cuts_sys_err_no_corr);

    TPaveText *K0s_K0s_text_MC_delta_phi_cuts = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    K0s_K0s_text_MC_delta_phi_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    K0s_K0s_text_MC_delta_phi_cuts->AddText("Minimum bias");
    K0s_K0s_text_MC_delta_phi_cuts->AddText("K_{s}^{0}-K_{s}^{0}");
    K0s_K0s_text_MC_delta_phi_cuts->AddText("True MC");
    K0s_K0s_text_MC_delta_phi_cuts->AddText("Analysis cuts");
    K0s_K0s_text_MC_delta_phi_cuts->AddText(Form("P = %.3f #pm %.3f", K0s_K0s_slope_delta_phi_cuts, K0s_K0s_slope_delta_phi_cuts_err ));
    K0s_K0s_text_MC_delta_phi_cuts->AddText(Form("P_{ME} = %.3f", fit_K0s_K0s_after_delta_phi_cuts_ME->GetParameter(1)/(L0_alpha*L0_alpha) ));
    //if(corr_err == 1) K0s_K0s_text_MC_delta_phi_cuts->AddText(Form("Slope sys. err: %.1f %%", K0s_K0s_slope_delta_phi_cuts_sys_err*100 ));
    //if(corr_err == 0) K0s_K0s_text_MC_delta_phi_cuts->AddText(Form("Slope sys. err: %.1f %%", K0s_K0s_slope_delta_phi_cuts_sys_err_no_corr*100 ));
    K0s_K0s_text_MC_delta_phi_cuts->SetFillColorAlpha(0, 0.01);
    K0s_K0s_text_MC_delta_phi_cuts->Draw("same");

    K0s_K0s_cosThetaProdPlane_delta_phi_cuts_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/K0s/K0s_K0s_cosThetaProdPlane_delta_phi_cuts_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------

    TCanvas *K0s_K0s_cosThetaProdPlane_delta_phi_cuts_corr_can = new TCanvas(Form("K0s_K0s_cosThetaProdPlane_delta_phi_cuts_corr_can_%i", delta_eta_bin), Form("K0s_K0s_cosThetaProdPlane_delta_phi_cuts_corr_can_%i", delta_eta_bin), 1200, 1000);
    K0s_K0s_cosThetaProdPlane_delta_phi_cuts_corr_can->cd();

    TF1 *fit_K0s_K0s_after_delta_phi_cuts_corr = new TF1("fit_K0s_K0s_after_delta_phi_cuts_corr", "[0]*(1+[1]*x)", -1, 1);
    fit_K0s_K0s_after_delta_phi_cuts_corr->SetParameters(1000, 0.10);

    TH1D *K0s_K0s_cosThetaProdPlane_delta_phi_cuts_corr = (TH1D*)K0s_K0s_cosThetaProdPlane_delta_phi_cuts->Clone(Form("K0s_K0s_cosThetaProdPlane_delta_phi_cuts_corr_%i", delta_eta_bin));
    K0s_K0s_cosThetaProdPlane_delta_phi_cuts_corr->Divide(K0s_K0s_cosThetaProdPlane_ME_delta_phi_cuts); //correct
    K0s_K0s_cosThetaProdPlane_delta_phi_cuts_corr->Scale(nK0sK0s_delta_phi_cuts/K0s_K0s_cosThetaProdPlane_delta_phi_cuts_corr->Integral()); //scale back
    K0s_K0s_cosThetaProdPlane_delta_phi_cuts_corr->Scale(1./K0s_K0s_cosThetaProdPlane_delta_phi_cuts_corr->GetXaxis()->GetBinWidth(1)); //divide by bin width
    //K0s_K0s_cosThetaProdPlane_delta_phi_cuts_corr->Scale(nK0sK0s_delta_phi_cuts);
    K0s_K0s_cosThetaProdPlane_delta_phi_cuts_corr->Fit(fit_K0s_K0s_after_delta_phi_cuts_corr, "i 0 r");
    K0s_K0s_cosThetaProdPlane_delta_phi_cuts_corr->SetMinimum(0);
    K0s_K0s_cosThetaProdPlane_delta_phi_cuts_corr->Draw("p e");

    fit_K0s_K0s_after_delta_phi_cuts_corr->SetLineColor(1);
    fit_K0s_K0s_after_delta_phi_cuts_corr->Draw("same");

    float K0s_K0s_slope_delta_phi_cuts_corr = fit_K0s_K0s_after_delta_phi_cuts_corr->GetParameter(1)/(L0_alpha*L0_alpha);
    float K0s_K0s_slope_delta_phi_cuts_corr_err = fit_K0s_K0s_after_delta_phi_cuts_corr->GetParError(1)/(L0_alpha*L0_alpha);

    ResidualPolarization_delta_phi_hist[delta_eta_bin-1]->SetBinContent(1, K0s_K0s_slope_delta_phi_cuts_corr);
    ResidualPolarization_delta_phi_hist[delta_eta_bin-1]->SetBinError(1, K0s_K0s_slope_delta_phi_cuts_corr_err);


    TPaveText *K0s_K0s_text_MC_delta_phi_cuts_corr = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
    K0s_K0s_text_MC_delta_phi_cuts_corr->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    K0s_K0s_text_MC_delta_phi_cuts_corr->AddText("Minimum bias");
    K0s_K0s_text_MC_delta_phi_cuts_corr->AddText("K_{s}^{0}-K_{s}^{0}");
    K0s_K0s_text_MC_delta_phi_cuts_corr->AddText("True MC");
    K0s_K0s_text_MC_delta_phi_cuts_corr->AddText("Analysis cuts");
    K0s_K0s_text_MC_delta_phi_cuts_corr->AddText(Form("P = %.3f #pm %.3f", K0s_K0s_slope_delta_phi_cuts_corr, K0s_K0s_slope_delta_phi_cuts_corr_err ));
    //if(corr_err == 1) K0s_K0s_text_MC_delta_phi_cuts_corr->AddText(Form("Slope sys. err: %.1f %%", K0s_K0s_slope_delta_phi_cuts_corr_sys_err*100 ));
    //if(corr_err == 0) K0s_K0s_text_MC_delta_phi_cuts_corr->AddText(Form("Slope sys. err: %.1f %%", K0s_K0s_slope_delta_phi_cuts_corr_sys_err_no_corr*100 ));
    K0s_K0s_text_MC_delta_phi_cuts_corr->SetFillColorAlpha(0, 0.01);
    K0s_K0s_text_MC_delta_phi_cuts_corr->Draw("same");


    K0s_K0s_cosThetaProdPlane_delta_phi_cuts_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/K0s/K0s_K0s_cosThetaProdPlane_delta_phi_cuts_corr_%i.png", delta_eta_bin));

    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


    //delta eta delta phi, after cuts
    TCanvas *K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts_can = new TCanvas(Form("K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts_can_%i", delta_eta_bin), Form("K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts_can_%i", delta_eta_bin), 1200, 1000);
    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts_can->cd();

    TF1 *fit_K0s_K0s_after_delta_eta_delta_phi_cuts = new TF1("fit_K0s_K0s_after_delta_eta_delta_phi_cuts", "[0]*(1+[1]*x)", -1, 1);
    fit_K0s_K0s_after_delta_eta_delta_phi_cuts->SetParameters(1000, 0.10);

    TF1 *fit_K0s_K0s_after_delta_eta_delta_phi_cuts_ME = new TF1("fit_K0s_K0s_after_delta_eta_delta_phi_cuts_ME", "[0]*(1+[1]*x)", -1, 1);
    fit_K0s_K0s_after_delta_eta_delta_phi_cuts_ME->SetParameters(1000, 0.10);

    TH1D *K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts = K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_cuts_hist->ProjectionX(Form("proj_eta_SE_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
    float nK0sK0s_delta_eta_delta_phi_cuts = K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts->Integral();

    TH1D *K0s_K0s_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts = K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_ME_cuts_hist->ProjectionX(Form("proj_eta_ME_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
    K0s_K0s_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts->Sumw2();
    K0s_K0s_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts->SetMarkerStyle(24);
    K0s_K0s_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts->SetMarkerSize(1.5);
    K0s_K0s_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts->SetMarkerColor(1);
    K0s_K0s_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts->SetLineColor(1);
    K0s_K0s_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts->Scale(nK0sK0s_delta_eta_delta_phi_cuts/K0s_K0s_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts->Integral());
    K0s_K0s_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts->Scale(1./K0s_K0s_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts->GetXaxis()->GetBinWidth(1));
    K0s_K0s_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts->Fit(fit_K0s_K0s_after_delta_eta_delta_phi_cuts_ME, "i 0 r");
    //K0s_K0s_cosThetaProdPlane_ME->Draw("p e");


    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts->Sumw2();
    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts->SetMarkerStyle(20);
    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts->SetMarkerSize(1.5);
    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts->SetMarkerColor(kRed);
    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts->SetLineColor(kRed);
    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts->GetXaxis()->SetTitle("cos(#theta*)");
    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts->GetXaxis()->CenterTitle();
    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts->GetYaxis()->CenterTitle();
    //K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts->Divide(K0s_K0s_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts); //correct using ME
    //K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts->Scale(nK0sK0s_delta_eta_delta_phi_cuts/K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts->Integral()); //scale back
    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts->Scale(1./K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts->Fit(fit_K0s_K0s_after_delta_eta_delta_phi_cuts, "i 0 r");
    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts->SetMinimum(0);
    //K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts->GetYaxis()->SetRangeUser(-20e9, 0);
    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts->Draw("p e");

    K0s_K0s_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts->Draw("p e same");

    fit_K0s_K0s_after_delta_eta_delta_phi_cuts->SetLineColor(1);
    fit_K0s_K0s_after_delta_eta_delta_phi_cuts->Draw("same");

    fit_K0s_K0s_after_delta_eta_delta_phi_cuts_ME->SetLineColor(1);
    fit_K0s_K0s_after_delta_eta_delta_phi_cuts_ME->SetLineStyle(7);
    fit_K0s_K0s_after_delta_eta_delta_phi_cuts_ME->Draw("same");

    TLegend *K0s_K0s_leg = new TLegend(0.15, 0.2, 0.45, 0.45);
    K0s_K0s_leg->AddEntry(K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts, "Same-event");
    K0s_K0s_leg->AddEntry(fit_K0s_K0s_after_delta_eta_delta_phi_cuts, "Same-event fit");
    K0s_K0s_leg->AddEntry(K0s_K0s_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts, "Mixed-event");
    K0s_K0s_leg->AddEntry(fit_K0s_K0s_after_delta_eta_delta_phi_cuts_ME, "Mixed-event fit");
    //K0s_K0s_leg->AddEntry(fitK0s_K0s_US_ThetaStar_no_corr, "Linear fit to US");
    K0s_K0s_leg->SetBorderSize(0);
    K0s_K0s_leg->SetFillColorAlpha(0, 0.01);
    K0s_K0s_leg->Draw("same");


    float K0s_K0s_slope_delta_eta_delta_phi_cuts = fit_K0s_K0s_after_delta_eta_delta_phi_cuts->GetParameter(1)/(L0_alpha*L0_alpha);
    float K0s_K0s_slope_delta_eta_delta_phi_cuts_err = fit_K0s_K0s_after_delta_eta_delta_phi_cuts->GetParError(1)/(L0_alpha*L0_alpha);


    //statistical error correction for systematic error
    float SysErrSlope_K0s_K0s_delta_eta_delta_phi_corr = sqrt( fabs( fit_K0s_K0s_after_delta_eta_delta_phi_cuts->GetParError(1)*fit_K0s_K0s_after_delta_eta_delta_phi_cuts->GetParError(1) - fit_K0s_K0s_after_delta_eta_delta_phi_cuts_ME->GetParError(1)*fit_K0s_K0s_after_delta_eta_delta_phi_cuts_ME->GetParError(1) ) );

    float SysErrSlope_K0s_K0s_delta_eta_delta_phi_work = ( fabs( fit_K0s_K0s_after_delta_eta_delta_phi_cuts->GetParameter(1) - fit_K0s_K0s_after_delta_eta_delta_phi_cuts_ME->GetParameter(1) ) - SysErrSlope_K0s_K0s_delta_eta_delta_phi_corr )/fabs(fit_K0s_K0s_after_delta_eta_delta_phi_cuts->GetParameter(1));

    float K0s_K0s_slope_delta_eta_delta_phi_cuts_sys_err = 0;

    if( SysErrSlope_K0s_K0s_delta_eta_delta_phi_work > 0 ) K0s_K0s_slope_delta_eta_delta_phi_cuts_sys_err = SysErrSlope_K0s_K0s_delta_eta_delta_phi_work; //store sys. err. only if it is larger than statistical fluctuations

    float K0s_K0s_slope_delta_eta_delta_phi_cuts_sys_err_no_corr = fabs( fit_K0s_K0s_after_delta_eta_delta_phi_cuts->GetParameter(1) - fit_K0s_K0s_after_delta_eta_delta_phi_cuts_ME->GetParameter(1) )/fabs(fit_K0s_K0s_after_delta_eta_delta_phi_cuts->GetParameter(1));

    if(corr_err == 1) SysErrSlope_delta_eta_delta_phi_hist[delta_eta_bin-1]->SetBinContent(1, K0s_K0s_slope_delta_eta_delta_phi_cuts_sys_err);
    if(corr_err == 0) SysErrSlope_delta_eta_delta_phi_hist[delta_eta_bin-1]->SetBinContent(1, K0s_K0s_slope_delta_eta_delta_phi_cuts_sys_err_no_corr);

    TPaveText *K0s_K0s_text_MC_delta_eta_delta_phi_cuts = new TPaveText(0.5, 0.2, 0.9, 0.45, "NDC");
    K0s_K0s_text_MC_delta_eta_delta_phi_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    K0s_K0s_text_MC_delta_eta_delta_phi_cuts->AddText("Minimum bias");
    K0s_K0s_text_MC_delta_eta_delta_phi_cuts->AddText("K_{s}^{0}-K_{s}^{0}");
    K0s_K0s_text_MC_delta_eta_delta_phi_cuts->AddText("True MC");
    K0s_K0s_text_MC_delta_eta_delta_phi_cuts->AddText("Analysis cuts");
    K0s_K0s_text_MC_delta_eta_delta_phi_cuts->AddText(Form("P = %.3f #pm %.3f", K0s_K0s_slope_delta_eta_delta_phi_cuts, K0s_K0s_slope_delta_eta_delta_phi_cuts_err ));
    K0s_K0s_text_MC_delta_eta_delta_phi_cuts->AddText(Form("P_{ME} = %.3f", fit_K0s_K0s_after_delta_eta_delta_phi_cuts_ME->GetParameter(1)/(L0_alpha*L0_alpha) ));
    if( (delta_eta_bin-1) == 0 )K0s_K0s_text_MC_delta_eta_delta_phi_cuts->AddText("|#Delta#it{y}| < 0.5, |#Delta#phi| < #pi/3");
    if( (delta_eta_bin-1) == 1 )K0s_K0s_text_MC_delta_eta_delta_phi_cuts->AddText("0.5 < |#Delta#it{y}| < 2.0 or |#Delta#phi| > #pi/3");
    //if(corr_err == 1) K0s_K0s_text_MC_delta_eta_delta_phi_cuts->AddText(Form("Slope sys. err: %.1f %%", K0s_K0s_slope_delta_eta_delta_phi_cuts_sys_err*100 ));
    //if(corr_err == 0) K0s_K0s_text_MC_delta_eta_delta_phi_cuts->AddText(Form("Slope sys. err: %.1f %%", K0s_K0s_slope_delta_eta_delta_phi_cuts_sys_err_no_corr*100 ));
    K0s_K0s_text_MC_delta_eta_delta_phi_cuts->SetFillColorAlpha(0, 0.01);
    K0s_K0s_text_MC_delta_eta_delta_phi_cuts->Draw("same");

    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/K0s/K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------

    TCanvas *K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr_can = new TCanvas(Form("K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr_can_%i", delta_eta_bin), Form("K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr_can_%i", delta_eta_bin), 1200, 1000);
    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr_can->cd();

    TF1 *fit_K0s_K0s_after_delta_eta_delta_phi_cuts_corr = new TF1("fit_K0s_K0s_after_delta_eta_delta_phi_cuts_corr", "[0]*(1+[1]*x)", -1, 1);
    fit_K0s_K0s_after_delta_eta_delta_phi_cuts_corr->SetParameters(1000, 0.10);

    TH1D *K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr = (TH1D*)K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts->Clone(Form("K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr_%i", delta_eta_bin));
    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr->Divide(K0s_K0s_cosThetaProdPlane_ME_delta_eta_delta_phi_cuts); //correct
    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr->Scale(nK0sK0s_delta_eta_delta_phi_cuts/K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr->Integral()); //scale back
    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr->Scale(1./K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr->GetXaxis()->GetBinWidth(1)); //divide by bin width
    //K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr->Scale(nK0sK0s_delta_eta_delta_phi_cuts);
    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr->Fit(fit_K0s_K0s_after_delta_eta_delta_phi_cuts_corr, "i 0 r");
    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr->SetMinimum(0);
    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr->Draw("p e");

    fit_K0s_K0s_after_delta_eta_delta_phi_cuts_corr->SetLineColor(1);
    fit_K0s_K0s_after_delta_eta_delta_phi_cuts_corr->Draw("same");

    float K0s_K0s_slope_delta_eta_delta_phi_cuts_corr = fit_K0s_K0s_after_delta_eta_delta_phi_cuts_corr->GetParameter(1)/(L0_alpha*L0_alpha);
    float K0s_K0s_slope_delta_eta_delta_phi_cuts_corr_err = fit_K0s_K0s_after_delta_eta_delta_phi_cuts_corr->GetParError(1)/(L0_alpha*L0_alpha);

    ResidualPolarization_delta_eta_delta_phi_hist[delta_eta_bin-1]->SetBinContent(1, K0s_K0s_slope_delta_eta_delta_phi_cuts_corr);
    ResidualPolarization_delta_eta_delta_phi_hist[delta_eta_bin-1]->SetBinError(1, K0s_K0s_slope_delta_eta_delta_phi_cuts_corr_err);


    TPaveText *K0s_K0s_text_MC_delta_eta_delta_phi_cuts_corr = new TPaveText(0.5, 0.2, 0.9, 0.45, "NDC");
    K0s_K0s_text_MC_delta_eta_delta_phi_cuts_corr->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    K0s_K0s_text_MC_delta_eta_delta_phi_cuts_corr->AddText("Minimum bias");
    K0s_K0s_text_MC_delta_eta_delta_phi_cuts_corr->AddText("K_{s}^{0}-K_{s}^{0}");
    K0s_K0s_text_MC_delta_eta_delta_phi_cuts_corr->AddText("True MC");
    K0s_K0s_text_MC_delta_eta_delta_phi_cuts_corr->AddText("Analysis cuts");
    K0s_K0s_text_MC_delta_eta_delta_phi_cuts_corr->AddText(Form("P = %.3f #pm %.3f", K0s_K0s_slope_delta_eta_delta_phi_cuts_corr, K0s_K0s_slope_delta_eta_delta_phi_cuts_corr_err ));
    if( (delta_eta_bin-1) == 0 )K0s_K0s_text_MC_delta_eta_delta_phi_cuts_corr->AddText("|#Delta#it{y}| < 0.5, |#Delta#phi| < #pi/3");
    if( (delta_eta_bin-1) == 1 )K0s_K0s_text_MC_delta_eta_delta_phi_cuts_corr->AddText("0.5 < |#Delta#it{y}| < 2.0 or |#Delta#phi| > #pi/3");
    //if(corr_err == 1) K0s_K0s_text_MC_delta_eta_delta_phi_cuts_corr->AddText(Form("Slope sys. err: %.1f %%", K0s_K0s_slope_delta_eta_delta_phi_cuts_corr_sys_err*100 ));
    //if(corr_err == 0) K0s_K0s_text_MC_delta_eta_delta_phi_cuts_corr->AddText(Form("Slope sys. err: %.1f %%", K0s_K0s_slope_delta_eta_delta_phi_cuts_corr_sys_err_no_corr*100 ));
    K0s_K0s_text_MC_delta_eta_delta_phi_cuts_corr->SetFillColorAlpha(0, 0.01);
    K0s_K0s_text_MC_delta_eta_delta_phi_cuts_corr->Draw("same");


    K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/K0s/K0s_K0s_cosThetaProdPlane_delta_eta_delta_phi_cuts_corr_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------


    SysErrSlope_delta_eta_hist[delta_eta_bin-1]->Write();
    ResidualPolarization_delta_eta_hist[delta_eta_bin-1]->Write();

    SysErrSlope_delta_phi_hist[delta_eta_bin-1]->Write();
    ResidualPolarization_delta_phi_hist[delta_eta_bin-1]->Write();

    SysErrSlope_delta_eta_delta_phi_hist[delta_eta_bin-1]->Write();
    ResidualPolarization_delta_eta_delta_phi_hist[delta_eta_bin-1]->Write();

  }


  //_____________________________________________________________________________________________________________________________________________________________________________________


  inFile->Close();

  sysErrFile->Close();
  out_file_polarization->Close();


  return;

}
