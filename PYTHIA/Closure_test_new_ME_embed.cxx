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
void Closure_test_new_ME_embed(const int energy = 200, const int corr_err = 0)
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

    //inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/2024_08_ME_weight_new/output_Lambda_pp_200_MB_1B_events_hists_Delta_phi_pi_quater.root", "READ");
    //inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/2024_08_ME_weight_new/output_Lambda_pp_200_MB_1B_events_hists_Delta_phi_pi_third.root", "READ");

    inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/Embed_weight/output_Lambda_pp_200_MB_1B_events_hists_work.root", "READ");

    //inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/ME_tests_Delta_phi_06_2024/N_SE_in_ME/output_Lambda_pp_200_MB_1B_events_hists_1_SE_in_ME.root", "READ");
    //inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/ME_tests_Delta_phi_06_2024/N_SE_in_ME/output_Lambda_pp_200_MB_1B_events_hists_10_SE_in_ME.root", "READ");
    //inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/ME_tests_Delta_phi_06_2024/N_SE_in_ME/output_Lambda_pp_200_MB_1B_events_hists_100_SE_in_ME.root", "READ");

  }
  else
  {
    cout<<"Not a valid collision energy! Abborting!"<<endl;
    return;
  }


  //histograms

  //before cuts
  //True MC
  TH1F *L0_L0bar_cosThetaProdPlane_embed = (TH1F*)inFile->Get("L0_L0bar_cosThetaProdPlane_embed");

  TH2F *L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_embed_hist = (TH2F*)inFile->Get("L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_embed_hist");
  TH2F *L0_L0bar_cos_theta_star_vs_delta_eta_embed_hist = (TH2F*)inFile->Get("L0_L0bar_cos_theta_star_vs_delta_eta_embed_hist");
  TH2F *L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_embed_hist = (TH2F*)inFile->Get("L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_embed_hist");

  //----------------------------------

  TH1F *L0_L0_cosThetaProdPlane_embed = (TH1F*)inFile->Get("L0_L0_cosThetaProdPlane_embed");

  TH2F *L0_L0_cos_theta_star_vs_delta_phi_for_corr_embed_hist = (TH2F*)inFile->Get("L0_L0_cos_theta_star_vs_delta_phi_for_corr_embed_hist");
  TH2F *L0_L0_cos_theta_star_vs_delta_eta_embed_hist = (TH2F*)inFile->Get("L0_L0_cos_theta_star_vs_delta_eta_embed_hist");
  TH2F *L0_L0_cos_theta_star_vs_delta_eta_delta_phi_embed_hist = (TH2F*)inFile->Get("L0_L0_cos_theta_star_vs_delta_eta_delta_phi_embed_hist");

  //----------------------------------

  TH1F *L0bar_L0bar_cosThetaProdPlane_embed = (TH1F*)inFile->Get("L0bar_L0bar_cosThetaProdPlane_embed");

  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_phi_for_corr_embed_hist = (TH2F*)inFile->Get("L0bar_L0bar_cos_theta_star_vs_delta_phi_for_corr_embed_hist");
  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_eta_embed_hist = (TH2F*)inFile->Get("L0bar_L0bar_cos_theta_star_vs_delta_eta_embed_hist");
  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_embed_hist = (TH2F*)inFile->Get("L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_embed_hist");


  //_____________________________________________________________________________________________________________________________


  //True MC before cuts
  TCanvas *L0_L0bar_cosThetaProdPlane_embed_can = new TCanvas("L0_L0bar_cosThetaProdPlane_embed_can", "L0_L0bar_cosThetaProdPlane_embed_can", 1200, 1000);
  TCanvas *L0_L0_cosThetaProdPlane_embed_can = new TCanvas("L0_L0_cosThetaProdPlane_embed_can", "L0_L0_cosThetaProdPlane_embed_can", 1200, 1000);
  TCanvas *L0bar_L0bar_cosThetaProdPlane_embed_can = new TCanvas("L0bar_L0bar_cosThetaProdPlane_embed_can", "L0bar_L0bar_cosThetaProdPlane_embed_can", 1200, 1000);


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

  L0_L0bar_cosThetaProdPlane_embed_can->cd();

  float nL0L0bar = L0_L0bar_cosThetaProdPlane_embed->Integral();

  //L0_L0bar_cosThetaProdPlane_embed_ME->Sumw2();
  //L0_L0bar_cosThetaProdPlane_embed_ME->Add(L0_L0bar_cosThetaProdPlane_embed_LS_ME);
  //L0_L0bar_cosThetaProdPlane_embed_ME->Scale(nL0L0bar/L0_L0bar_cosThetaProdPlane_embed_ME->Integral());
  //L0_L0bar_cosThetaProdPlane_embed_ME->Draw("p e");

  L0_L0bar_cosThetaProdPlane_embed->Sumw2();
  L0_L0bar_cosThetaProdPlane_embed->SetMarkerStyle(20);
  L0_L0bar_cosThetaProdPlane_embed->SetMarkerColor(kRed);
  L0_L0bar_cosThetaProdPlane_embed->SetLineColor(kRed);
  L0_L0bar_cosThetaProdPlane_embed->GetXaxis()->SetTitle("cos(#theta*)");
  L0_L0bar_cosThetaProdPlane_embed->GetXaxis()->CenterTitle();
  L0_L0bar_cosThetaProdPlane_embed->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0_L0bar_cosThetaProdPlane_embed->GetYaxis()->CenterTitle();
  //L0_L0bar_cosThetaProdPlane_embed->Divide(L0_L0bar_cosThetaProdPlane_embed_ME); //correct using ME
  //L0_L0bar_cosThetaProdPlane_embed->Scale(nL0L0bar/L0_L0bar_cosThetaProdPlane_embed->Integral()); //scale back
  L0_L0bar_cosThetaProdPlane_embed->Scale(1./L0_L0bar_cosThetaProdPlane_embed->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
  L0_L0bar_cosThetaProdPlane_embed->Fit(fit_L0_L0bar_before_cuts, "i 0 r");
  L0_L0bar_cosThetaProdPlane_embed->SetMinimum(0);
  //L0_L0bar_cosThetaProdPlane_embed->GetYaxis()->SetRangeUser(-20e9, 0);
  L0_L0bar_cosThetaProdPlane_embed->Draw("p e");

  //L0_L0bar_cosThetaProdPlane_embed_ME->Draw("p e same");

  fit_L0_L0bar_before_cuts->SetLineColor(1);
  fit_L0_L0bar_before_cuts->Draw("same");


  float L0_L0bar_slope = fit_L0_L0bar_before_cuts->GetParameter(1)/(L0_alpha*L0bar_alpha);
  float L0_L0bar_slope_err = fit_L0_L0bar_before_cuts->GetParError(1)/(L0_alpha*L0bar_alpha);


  TPaveText *L0_L0bar_text_MC = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
  L0_L0bar_text_MC->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  L0_L0bar_text_MC->AddText("Minimum bias");
  L0_L0bar_text_MC->AddText("#Lambda-#bar{#Lambda}");
  L0_L0bar_text_MC->AddText("True MC");
  L0_L0bar_text_MC->AddText("Embed. eff.");
  L0_L0bar_text_MC->AddText(Form("Fit slope: %.3f #pm %.3f", L0_L0bar_slope, L0_L0bar_slope_err));
  L0_L0bar_text_MC->SetFillColorAlpha(0, 0.01);
  L0_L0bar_text_MC->Draw("same");

  L0_L0bar_cosThetaProdPlane_embed_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/L0_L0bar_cosThetaProdPlane_embed.png");


  //___________________________________________________________________________________________________________________________________________________________________________________________________________________________

  //L-L
  //before cuts
  TF1 *fit_L0_L0_before_cuts = new TF1("fit_L0_L0_before_cuts", "[0]*(1+[1]*x)", -1, 1);
  fit_L0_L0_before_cuts->SetParameters(1000, 0.10);

  L0_L0_cosThetaProdPlane_embed_can->cd();

  float nL0L0 = L0_L0_cosThetaProdPlane_embed->Integral();

  //L0_L0_cosThetaProdPlane_embed_ME->Sumw2();
  //L0_L0_cosThetaProdPlane_embed_ME->Add(L0_L0_cosThetaProdPlane_embed_LS_ME);
  //L0_L0_cosThetaProdPlane_embed_ME->Scale(nL0L0/L0_L0_cosThetaProdPlane_embed_ME->Integral());
  //L0_L0_cosThetaProdPlane_embed_ME->Draw("p e");

  L0_L0_cosThetaProdPlane_embed->Sumw2();
  L0_L0_cosThetaProdPlane_embed->SetMarkerStyle(20);
  L0_L0_cosThetaProdPlane_embed->SetMarkerColor(kRed);
  L0_L0_cosThetaProdPlane_embed->SetLineColor(kRed);
  L0_L0_cosThetaProdPlane_embed->GetXaxis()->SetTitle("cos(#theta*)");
  L0_L0_cosThetaProdPlane_embed->GetXaxis()->CenterTitle();
  L0_L0_cosThetaProdPlane_embed->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0_L0_cosThetaProdPlane_embed->GetYaxis()->CenterTitle();
  //L0_L0_cosThetaProdPlane_embed->Divide(L0_L0_cosThetaProdPlane_embed_ME); //correct using ME
  //L0_L0_cosThetaProdPlane_embed->Scale(nL0L0/L0_L0_cosThetaProdPlane_embed->Integral()); //scale back
  L0_L0_cosThetaProdPlane_embed->Scale(1./L0_L0_cosThetaProdPlane_embed->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
  L0_L0_cosThetaProdPlane_embed->Fit(fit_L0_L0_before_cuts, "i 0 r");
  L0_L0_cosThetaProdPlane_embed->SetMinimum(0);
  //L0_L0_cosThetaProdPlane_embed->GetYaxis()->SetRangeUser(-20e9, 0);
  L0_L0_cosThetaProdPlane_embed->Draw("p e");

  //L0_L0_cosThetaProdPlane_embed_ME->Draw("p e same");

  fit_L0_L0_before_cuts->SetLineColor(1);
  fit_L0_L0_before_cuts->Draw("same");

  float L0_L0_slope = fit_L0_L0_before_cuts->GetParameter(1)/(L0_alpha*L0_alpha);
  float L0_L0_slope_err = fit_L0_L0_before_cuts->GetParError(1)/(L0_alpha*L0_alpha);

  TPaveText *L0_L0_text_MC = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
  L0_L0_text_MC->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  L0_L0_text_MC->AddText("Minimum bias");
  L0_L0_text_MC->AddText("#Lambda-#Lambda");
  L0_L0_text_MC->AddText("True MC");
  L0_L0_text_MC->AddText("Embed. eff.");
  L0_L0_text_MC->AddText(Form("P = %.3f #pm %.3f", L0_L0_slope, L0_L0_slope_err));
  L0_L0_text_MC->SetFillColorAlpha(0, 0.01);
  L0_L0_text_MC->Draw("same");

  L0_L0_cosThetaProdPlane_embed_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_L/L0_L0_cosThetaProdPlane_embed.png");


  //_________________________________________________________________________

  //Lbar-Lbar
  //before cuts
  TF1 *fit_L0bar_L0bar_before_cuts = new TF1("fit_L0bar_L0bar_before_cuts", "[0]*(1+[1]*x)", -1, 1);
  fit_L0bar_L0bar_before_cuts->SetParameters(1000, 0.10);

  L0bar_L0bar_cosThetaProdPlane_embed_can->cd();

  float nL0barL0bar = L0bar_L0bar_cosThetaProdPlane_embed->Integral();

  //L0bar_L0bar_cosThetaProdPlane_embed_ME->Sumw2();
  //L0bar_L0bar_cosThetaProdPlane_embed_ME->Add(L0bar_L0bar_cosThetaProdPlane_embed_LS_ME);
  //L0bar_L0bar_cosThetaProdPlane_embed_ME->Scale(nL0barL0bar/L0bar_L0bar_cosThetaProdPlane_embed_ME->Integral());
  //L0bar_L0bar_cosThetaProdPlane_embed_ME->Draw("p e");

  L0bar_L0bar_cosThetaProdPlane_embed->Sumw2();
  L0bar_L0bar_cosThetaProdPlane_embed->SetMarkerStyle(20);
  L0bar_L0bar_cosThetaProdPlane_embed->SetMarkerColor(kRed);
  L0bar_L0bar_cosThetaProdPlane_embed->SetLineColor(kRed);
  L0bar_L0bar_cosThetaProdPlane_embed->GetXaxis()->SetTitle("cos(#theta*)");
  L0bar_L0bar_cosThetaProdPlane_embed->GetXaxis()->CenterTitle();
  L0bar_L0bar_cosThetaProdPlane_embed->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0bar_L0bar_cosThetaProdPlane_embed->GetYaxis()->CenterTitle();
  //L0bar_L0bar_cosThetaProdPlane_embed->Divide(L0bar_L0bar_cosThetaProdPlane_embed_ME); //correct using ME
  //L0bar_L0bar_cosThetaProdPlane_embed->Scale(nL0barL0bar/L0bar_L0bar_cosThetaProdPlane_embed->Integral()); //scale back
  L0bar_L0bar_cosThetaProdPlane_embed->Scale(1./L0bar_L0bar_cosThetaProdPlane_embed->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
  L0bar_L0bar_cosThetaProdPlane_embed->Fit(fit_L0bar_L0bar_before_cuts, "i 0 r");
  L0bar_L0bar_cosThetaProdPlane_embed->SetMinimum(0);
  //L0bar_L0bar_cosThetaProdPlane_embed->GetYaxis()->SetRangeUser(-20e9, 0);
  L0bar_L0bar_cosThetaProdPlane_embed->Draw("p e");

  //L0bar_L0bar_cosThetaProdPlane_embed_ME->Draw("p e same");

  fit_L0bar_L0bar_before_cuts->SetLineColor(1);
  fit_L0bar_L0bar_before_cuts->Draw("same");


  float L0bar_L0bar_slope = fit_L0bar_L0bar_before_cuts->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
  float L0bar_L0bar_slope_err = fit_L0bar_L0bar_before_cuts->GetParError(1)/(L0bar_alpha*L0bar_alpha);

  TPaveText *L0bar_L0bar_text_MC = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
  L0bar_L0bar_text_MC->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  L0bar_L0bar_text_MC->AddText("Minimum bias");
  L0bar_L0bar_text_MC->AddText("#bar{#Lambda}-#bar{#Lambda}");
  L0bar_L0bar_text_MC->AddText("True MC");
  L0bar_L0bar_text_MC->AddText("Embed. eff.");
  L0bar_L0bar_text_MC->AddText(Form("P = %.3f #pm %.3f", L0bar_L0bar_slope, L0bar_L0bar_slope_err));
  L0bar_L0bar_text_MC->SetFillColorAlpha(0, 0.01);
  L0bar_L0bar_text_MC->Draw("same");

  L0bar_L0bar_cosThetaProdPlane_embed_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/Lbar_Lbar/L0bar_L0bar_cosThetaProdPlane_embed.png");


  //_________________________________________________________________________________________________________________________________________________________________________

  //for Delta phi projections

  int bin_min = 1;
  int bin_middle = 20;
  int bin_max = 60;

  for( unsigned int delta_eta_bin = 1; delta_eta_bin < 3; delta_eta_bin++)
  {
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

    //before cuts
    //L-Lbar
    TCanvas *L0_L0bar_cosThetaProdPlane_delta_eta_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_delta_eta_can_%i", delta_eta_bin), Form("L0_L0bar_cosThetaProdPlane_delta_eta_can_%i", delta_eta_bin), 1200, 1000);
    L0_L0bar_cosThetaProdPlane_delta_eta_can->cd();

    TF1 *fit_L0_L0bar_after_delta_eta = new TF1("fit_L0_L0bar_after_delta_eta", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0bar_after_delta_eta->SetParameters(1000, 0.10);

    TH1D *L0_L0bar_cosThetaProdPlane_delta_eta = L0_L0bar_cos_theta_star_vs_delta_eta_embed_hist->ProjectionX(Form("L_Lbar_proj_eta_SE_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
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
    L0_L0bar_text_MC_delta_eta->AddText("Embed. eff.");
    L0_L0bar_text_MC_delta_eta->AddText(Form("P = %.3f #pm %.3f", L0_L0bar_slope_delta_eta, L0_L0bar_slope_delta_eta_err ));
    L0_L0bar_text_MC_delta_eta->SetFillColorAlpha(0, 0.01);
    L0_L0bar_text_MC_delta_eta->Draw("same");

    L0_L0bar_cosThetaProdPlane_delta_eta_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/L0_L0bar_cosThetaProdPlane_delta_eta_embed_%i.png", delta_eta_bin));

    //-----------------------------------------------------------

    //L-L
    TCanvas *L0_L0_cosThetaProdPlane_delta_eta_can = new TCanvas(Form("L0_L0_cosThetaProdPlane_delta_eta_can_%i", delta_eta_bin), Form("L0_L0_cosThetaProdPlane_delta_eta_can_%i", delta_eta_bin), 1200, 1000);
    L0_L0_cosThetaProdPlane_delta_eta_can->cd();

    TF1 *fit_L0_L0_after_delta_eta = new TF1("fit_L0_L0_after_delta_eta", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0_after_delta_eta->SetParameters(1000, 0.10);

    TH1D *L0_L0_cosThetaProdPlane_delta_eta = L0_L0_cos_theta_star_vs_delta_eta_embed_hist->ProjectionX(Form("LL_proj_eta_SE_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
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
    L0_L0_text_MC_delta_eta->AddText("Embed. eff.");
    L0_L0_text_MC_delta_eta->AddText(Form("P = %.3f #pm %.3f", L0_L0_slope_delta_eta, L0_L0_slope_delta_eta_err ));
    L0_L0_text_MC_delta_eta->SetFillColorAlpha(0, 0.01);
    L0_L0_text_MC_delta_eta->Draw("same");

    L0_L0_cosThetaProdPlane_delta_eta_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_L/L0_L0_cosThetaProdPlane_delta_eta_embed_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------

    //Lbar-Lbar
    TCanvas *L0bar_L0bar_cosThetaProdPlane_delta_eta_can = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_delta_eta_can_%i", delta_eta_bin), Form("L0bar_L0bar_cosThetaProdPlane_delta_eta_can_%i", delta_eta_bin), 1200, 1000);
    L0bar_L0bar_cosThetaProdPlane_delta_eta_can->cd();

    TF1 *fit_L0bar_L0bar_after_delta_eta = new TF1("fit_L0bar_L0bar_after_delta_eta", "[0]*(1+[1]*x)", -1, 1);
    fit_L0bar_L0bar_after_delta_eta->SetParameters(1000, 0.10);

    TH1D *L0bar_L0bar_cosThetaProdPlane_delta_eta = L0bar_L0bar_cos_theta_star_vs_delta_eta_embed_hist->ProjectionX(Form("LbaLbar_proj_eta_SE_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
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
    L0bar_L0bar_text_MC_delta_eta->AddText("Embed. eff.");
    L0bar_L0bar_text_MC_delta_eta->AddText(Form("P = %.3f #pm %.3f", L0bar_L0bar_slope_delta_eta, L0bar_L0bar_slope_delta_eta_err ));
    L0bar_L0bar_text_MC_delta_eta->SetFillColorAlpha(0, 0.01);
    L0bar_L0bar_text_MC_delta_eta->Draw("same");

    L0bar_L0bar_cosThetaProdPlane_delta_eta_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/Lbar_Lbar/L0bar_L0bar_cosThetaProdPlane_delta_eta_embed_%i.png", delta_eta_bin));

    //-----------------------------------------------------------------------------------------------------------------------------

    //before cuts
    //L-Lbar
    TCanvas *L0_L0bar_cosThetaProdPlane_delta_phi_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_delta_phi_can_%i", delta_eta_bin), Form("L0_L0bar_cosThetaProdPlane_delta_phi_can_%i", delta_eta_bin), 1200, 1000);
    L0_L0bar_cosThetaProdPlane_delta_phi_can->cd();

    TF1 *fit_L0_L0bar_after_delta_phi = new TF1("fit_L0_L0bar_after_delta_phi", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0bar_after_delta_phi->SetParameters(1000, 0.10);

    TH1D *L0_L0bar_cosThetaProdPlane_delta_phi = L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_embed_hist->ProjectionX(Form("L_Lbar_proj_phi_SE_%i", delta_eta_bin), bin_min_loop, bin_max_loop);
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
    L0_L0bar_text_MC_delta_phi->AddText("Embed. eff.");
    L0_L0bar_text_MC_delta_phi->AddText(Form("P = %.3f #pm %.3f", L0_L0bar_slope_delta_phi, L0_L0bar_slope_delta_phi_err ));
    L0_L0bar_text_MC_delta_phi->SetFillColorAlpha(0, 0.01);
    L0_L0bar_text_MC_delta_phi->Draw("same");

    L0_L0bar_cosThetaProdPlane_delta_phi_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/L0_L0bar_cosThetaProdPlane_delta_phi_embed_%i.png", delta_eta_bin));

    //-----------------------------------------------------------

    //L-L
    TCanvas *L0_L0_cosThetaProdPlane_delta_phi_can = new TCanvas(Form("L0_L0_cosThetaProdPlane_delta_phi_can_%i", delta_eta_bin), Form("L0_L0_cosThetaProdPlane_delta_phi_can_%i", delta_eta_bin), 1200, 1000);
    L0_L0_cosThetaProdPlane_delta_phi_can->cd();

    TF1 *fit_L0_L0_after_delta_phi = new TF1("fit_L0_L0_after_delta_phi", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0_after_delta_phi->SetParameters(1000, 0.10);

    TH1D *L0_L0_cosThetaProdPlane_delta_phi = L0_L0_cos_theta_star_vs_delta_phi_for_corr_embed_hist->ProjectionX(Form("LL_proj_phi_SE_%i", delta_eta_bin), bin_min_loop, bin_max_loop);
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
    L0_L0_text_MC_delta_phi->AddText("Embed. eff.");
    L0_L0_text_MC_delta_phi->AddText(Form("P = %.3f #pm %.3f", L0_L0_slope_delta_phi, L0_L0_slope_delta_phi_err ));
    L0_L0_text_MC_delta_phi->SetFillColorAlpha(0, 0.01);
    L0_L0_text_MC_delta_phi->Draw("same");

    L0_L0_cosThetaProdPlane_delta_phi_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_L/L0_L0_cosThetaProdPlane_delta_phi_embed_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------

    //Lbar-Lbar
    TCanvas *L0bar_L0bar_cosThetaProdPlane_delta_phi_can = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_delta_phi_can_%i", delta_eta_bin), Form("L0bar_L0bar_cosThetaProdPlane_delta_phi_can_%i", delta_eta_bin), 1200, 1000);
    L0bar_L0bar_cosThetaProdPlane_delta_phi_can->cd();

    TF1 *fit_L0bar_L0bar_after_delta_phi = new TF1("fit_L0bar_L0bar_after_delta_phi", "[0]*(1+[1]*x)", -1, 1);
    fit_L0bar_L0bar_after_delta_phi->SetParameters(1000, 0.10);

    TH1D *L0bar_L0bar_cosThetaProdPlane_delta_phi = L0bar_L0bar_cos_theta_star_vs_delta_phi_for_corr_embed_hist->ProjectionX(Form("LbaLbar_proj_phi_SE_%i", delta_eta_bin), bin_min_loop, bin_max_loop);
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
    L0bar_L0bar_text_MC_delta_phi->AddText("Embed. eff.");
    L0bar_L0bar_text_MC_delta_phi->AddText(Form("P = %.3f #pm %.3f", L0bar_L0bar_slope_delta_phi, L0bar_L0bar_slope_delta_phi_err ));
    L0bar_L0bar_text_MC_delta_phi->SetFillColorAlpha(0, 0.01);
    L0bar_L0bar_text_MC_delta_phi->Draw("same");

    L0bar_L0bar_cosThetaProdPlane_delta_phi_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/Lbar_Lbar/L0bar_L0bar_cosThetaProdPlane_delta_phi_embed_%i.png", delta_eta_bin));

    //-----------------------------------------------------------------------------------------------------------------------------

    //before cuts
    //L-Lbar
    TCanvas *L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_can_%i", delta_eta_bin), Form("L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_can_%i", delta_eta_bin), 1200, 1000);
    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_can->cd();

    TF1 *fit_L0_L0bar_after_delta_eta_delta_phi = new TF1("fit_L0_L0bar_after_delta_eta_delta_phi", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0bar_after_delta_eta_delta_phi->SetParameters(1000, 0.10);

    TH1D *L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi = L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_embed_hist->ProjectionX(Form("LLbar_proj_eta_phi_SE_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
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
    L0_L0bar_text_MC_delta_eta_delta_phi->AddText("Embed. eff.");
    L0_L0bar_text_MC_delta_eta_delta_phi->AddText(Form("P = %.3f #pm %.3f", L0_L0bar_slope_delta_eta_delta_phi, L0_L0bar_slope_delta_eta_delta_phi_err ));
    if( (delta_eta_bin-1) == 0 )L0_L0bar_text_MC_delta_eta_delta_phi->AddText("|#Delta#it{y}| < 0.5, |#Delta#phi| < #pi/3");
    if( (delta_eta_bin-1) == 1 )L0_L0bar_text_MC_delta_eta_delta_phi->AddText("0.5 < |#Delta#it{y}| < 2.0 or |#Delta#phi| > #pi/3");
    L0_L0bar_text_MC_delta_eta_delta_phi->SetFillColorAlpha(0, 0.01);
    L0_L0bar_text_MC_delta_eta_delta_phi->Draw("same");

    L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_Lbar/L0_L0bar_cosThetaProdPlane_delta_eta_delta_phi_embed_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------

    //L-L
    TCanvas *L0_L0_cosThetaProdPlane_delta_eta_delta_phi_can = new TCanvas(Form("L0_L0_cosThetaProdPlane_delta_eta_delta_phi_can_%i", delta_eta_bin), Form("L0_L0_cosThetaProdPlane_delta_eta_delta_phi_can_%i", delta_eta_bin), 1200, 1000);
    L0_L0_cosThetaProdPlane_delta_eta_delta_phi_can->cd();

    TF1 *fit_L0_L0_after_delta_eta_delta_phi = new TF1("fit_L0_L0_after_delta_eta_delta_phi", "[0]*(1+[1]*x)", -1, 1);
    fit_L0_L0_after_delta_eta_delta_phi->SetParameters(1000, 0.10);

    TH1D *L0_L0_cosThetaProdPlane_delta_eta_delta_phi = L0_L0_cos_theta_star_vs_delta_eta_delta_phi_embed_hist->ProjectionX(Form("LL_proj_eta_phi_SE_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
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
    L0_L0_text_MC_delta_eta_delta_phi->AddText("Embed. eff.");
    L0_L0_text_MC_delta_eta_delta_phi->AddText(Form("P = %.3f #pm %.3f", L0_L0_slope_delta_eta_delta_phi, L0_L0_slope_delta_eta_delta_phi_err ));
    L0_L0_text_MC_delta_eta_delta_phi->SetFillColorAlpha(0, 0.01);
    L0_L0_text_MC_delta_eta_delta_phi->Draw("same");

    L0_L0_cosThetaProdPlane_delta_eta_delta_phi_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/L_L/L0_L0_cosThetaProdPlane_delta_eta_delta_phi_embed_%i.png", delta_eta_bin));

    //-------------------------------------------------------------------------------------------------------

    //Lbar-Lbar
    TCanvas *L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_can = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_can_%i", delta_eta_bin), Form("L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_can_%i", delta_eta_bin), 1200, 1000);
    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_can->cd();

    TF1 *fit_L0bar_L0bar_after_delta_eta_delta_phi = new TF1("fit_L0bar_L0bar_after_delta_eta_delta_phi", "[0]*(1+[1]*x)", -1, 1);
    fit_L0bar_L0bar_after_delta_eta_delta_phi->SetParameters(1000, 0.10);

    TH1D *L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi = L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_embed_hist->ProjectionX(Form("LbarLbar_proj_eta_phi_SE_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);
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
    L0bar_L0bar_text_MC_delta_eta_delta_phi->AddText("Embed. eff.");
    L0bar_L0bar_text_MC_delta_eta_delta_phi->AddText(Form("P = %.3f #pm %.3f", L0bar_L0bar_slope_delta_eta_delta_phi, L0bar_L0bar_slope_delta_eta_delta_phi_err ));
    L0bar_L0bar_text_MC_delta_eta_delta_phi->SetFillColorAlpha(0, 0.01);
    L0bar_L0bar_text_MC_delta_eta_delta_phi->Draw("same");

    L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda/Lbar_Lbar/L0bar_L0bar_cosThetaProdPlane_delta_eta_delta_phi_embed_%i.png", delta_eta_bin));


    //_______________________________________________________________________________________________________________________________________________________________________________________________________________________
  }


  inFile->Close();


  //outFile->Close();

  return;
}
