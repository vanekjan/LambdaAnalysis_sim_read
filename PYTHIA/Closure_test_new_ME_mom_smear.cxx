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
void Closure_test_new_ME_mom_smear(const int energy = 200, const int corr_err = 0)
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

    //inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/2024_08_ME_weight_new/output_Lambda_pp_200_MB_1B_events_hists_Delta_phi_pi_third.root", "READ");
    //inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/2024_08_ME_weight_new/output_Lambda_pp_200_MB_1B_events_hists_Delta_phi_pi_quater.root", "READ");

    //momentum smearing
    //inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/mom_smear_new/output_Lambda_pp_200_MB_1B_events_hists_pair_weight_mom_smear_full.root", "READ");
    inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/mom_smear_new/output_Lambda_pp_200_MB_1B_events_hists_pair_weight_no_smear_full.root", "READ");

    //inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/mom_smear_pT/output_Lambda_pp_200_MB_1B_events_hists_pT.root", "READ");
    //inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/mom_smear_pT/output_Lambda_pp_200_MB_1B_events_hists_pT_times_20.root", "READ");

  }
  else
  {
    cout<<"Not a valid collision energy! Abborting!"<<endl;
    return;
  }






  //histograms

  //before cuts
  //True MC
  TH1F *L0_L0bar_cosThetaProdPlane = (TH1F*)inFile->Get("L0_L0bar_cosThetaProdPlane");

  TH1F *L0_L0_cosThetaProdPlane = (TH1F*)inFile->Get("L0_L0_cosThetaProdPlane");

  TH1F *L0bar_L0bar_cosThetaProdPlane = (TH1F*)inFile->Get("L0bar_L0bar_cosThetaProdPlane");

  //_____________________________________________________________________________________________________________________________

  //after cuts

  //True MC
  TH1F *L0_L0bar_cosThetaProdPlane_cuts = (TH1F*)inFile->Get("L0_L0bar_cosThetaProdPlane_cuts");

  TH1F *L0_L0_cosThetaProdPlane_cuts = (TH1F*)inFile->Get("L0_L0_cosThetaProdPlane_cuts");

  TH1F *L0bar_L0bar_cosThetaProdPlane_cuts = (TH1F*)inFile->Get("L0bar_L0bar_cosThetaProdPlane_cuts");

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



  L0_L0bar_cosThetaProdPlane->Sumw2();
  L0_L0bar_cosThetaProdPlane->SetMarkerStyle(20);
  L0_L0bar_cosThetaProdPlane->SetMarkerSize(1.5);
  L0_L0bar_cosThetaProdPlane->SetMarkerColor(kRed);
  L0_L0bar_cosThetaProdPlane->SetLineColor(kRed);
  L0_L0bar_cosThetaProdPlane->GetXaxis()->SetTitle("cos(#theta*)");
  L0_L0bar_cosThetaProdPlane->GetXaxis()->CenterTitle();
  L0_L0bar_cosThetaProdPlane->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0_L0bar_cosThetaProdPlane->GetYaxis()->CenterTitle();
  L0_L0bar_cosThetaProdPlane->Scale(nL0L0bar/L0_L0bar_cosThetaProdPlane->Integral()); //scale back
  L0_L0bar_cosThetaProdPlane->Scale(1./L0_L0bar_cosThetaProdPlane->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
  L0_L0bar_cosThetaProdPlane->Fit(fit_L0_L0bar_before_cuts, "i 0 r");
  L0_L0bar_cosThetaProdPlane->SetMinimum(0);
  //L0_L0bar_cosThetaProdPlane->GetYaxis()->SetRangeUser(-20e9, 0);
  L0_L0bar_cosThetaProdPlane->Draw("p e");


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

  L0_L0bar_cosThetaProdPlane_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda_mom_smear/L0_L0bar_cosThetaProdPlane.png");


  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  //after cuts
  TF1 *fit_L0_L0bar_after_cuts = new TF1("fit_L0_L0bar_after_cuts", "[0]*(1+[1]*x)", -1, 1);
  fit_L0_L0bar_after_cuts->SetParameters(1000, 0.10);

  L0_L0bar_cosThetaProdPlane_cuts_can->cd();

  float nL0L0bar_cuts = L0_L0bar_cosThetaProdPlane_cuts->Integral();


  L0_L0bar_cosThetaProdPlane_cuts->Sumw2();
  L0_L0bar_cosThetaProdPlane_cuts->SetMarkerStyle(20);
  L0_L0bar_cosThetaProdPlane_cuts->SetMarkerSize(1.5);
  L0_L0bar_cosThetaProdPlane_cuts->SetMarkerColor(kRed);
  L0_L0bar_cosThetaProdPlane_cuts->SetLineColor(kRed);
  L0_L0bar_cosThetaProdPlane_cuts->GetXaxis()->SetTitle("cos(#theta*)");
  L0_L0bar_cosThetaProdPlane_cuts->GetXaxis()->CenterTitle();
  L0_L0bar_cosThetaProdPlane_cuts->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0_L0bar_cosThetaProdPlane_cuts->GetYaxis()->CenterTitle();
  //L0_L0bar_cosThetaProdPlane_cuts->Scale(nL0L0bar_cuts/L0_L0bar_cosThetaProdPlane_cuts->Integral()); //scale back
  L0_L0bar_cosThetaProdPlane_cuts->Scale(1./L0_L0bar_cosThetaProdPlane_cuts->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
  L0_L0bar_cosThetaProdPlane_cuts->Fit(fit_L0_L0bar_after_cuts, "i 0 r");
  L0_L0bar_cosThetaProdPlane_cuts->SetMinimum(0);
  //L0_L0bar_cosThetaProdPlane_cuts->GetYaxis()->SetRangeUser(-20e9, 0);
  L0_L0bar_cosThetaProdPlane_cuts->Draw("p e");


  fit_L0_L0bar_after_cuts->SetLineColor(1);
  fit_L0_L0bar_after_cuts->Draw("same");


  float L0_L0bar_slope_cuts = fit_L0_L0bar_after_cuts->GetParameter(1)/(L0_alpha*L0bar_alpha);
  float L0_L0bar_slope_cuts_err = fit_L0_L0bar_after_cuts->GetParError(1)/(L0_alpha*L0bar_alpha);


  TPaveText *L0_L0bar_text_MC_cuts = new TPaveText(0.5, 0.2, 0.8, 0.5, "NDC");
  L0_L0bar_text_MC_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  L0_L0bar_text_MC_cuts->AddText("Minimum bias");
  L0_L0bar_text_MC_cuts->AddText("#Lambda-#bar{#Lambda}");
  L0_L0bar_text_MC_cuts->AddText("True MC");
  L0_L0bar_text_MC_cuts->AddText("Analysis cuts");
  L0_L0bar_text_MC_cuts->AddText("No smearing");
  //L0_L0bar_text_MC_cuts->AddText("#sigma_{p_{T}} smearing");
  //L0_L0bar_text_MC_cuts->AddText("20#sigma_{p_{T}} smearing");
  L0_L0bar_text_MC_cuts->AddText(Form("P = %.3f #pm %.3f", L0_L0bar_slope_cuts, L0_L0bar_slope_cuts_err ));
  L0_L0bar_text_MC_cuts->SetFillColorAlpha(0, 0.01);
  L0_L0bar_text_MC_cuts->Draw("same");

  L0_L0bar_cosThetaProdPlane_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda_mom_smear/L0_L0bar_cosThetaProdPlane_cuts.png");


  //___________________________________________________________________________________________________________________________________________________________________________________________________________________________

  //L-L
  //before cuts
  TF1 *fit_L0_L0_before_cuts = new TF1("fit_L0_L0_before_cuts", "[0]*(1+[1]*x)", -1, 1);
  fit_L0_L0_before_cuts->SetParameters(1000, 0.10);

  L0_L0_cosThetaProdPlane_can->cd();

  float nL0L0 = L0_L0_cosThetaProdPlane->Integral();


  L0_L0_cosThetaProdPlane->Sumw2();
  L0_L0_cosThetaProdPlane->SetMarkerStyle(20);
  L0_L0_cosThetaProdPlane->SetMarkerColor(kRed);
  L0_L0_cosThetaProdPlane->SetLineColor(kRed);
  L0_L0_cosThetaProdPlane->GetXaxis()->SetTitle("cos(#theta*)");
  L0_L0_cosThetaProdPlane->GetXaxis()->CenterTitle();
  L0_L0_cosThetaProdPlane->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0_L0_cosThetaProdPlane->GetYaxis()->CenterTitle();
  L0_L0_cosThetaProdPlane->Scale(nL0L0/L0_L0_cosThetaProdPlane->Integral()); //scale back
  //L0_L0_cosThetaProdPlane->Scale(1./L0_L0_cosThetaProdPlane->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
  L0_L0_cosThetaProdPlane->Fit(fit_L0_L0_before_cuts, "i 0 r");
  L0_L0_cosThetaProdPlane->SetMinimum(0);
  //L0_L0_cosThetaProdPlane->GetYaxis()->SetRangeUser(-20e9, 0);
  L0_L0_cosThetaProdPlane->Draw("p e");

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

  L0_L0_cosThetaProdPlane_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda_mom_smear/L0_L0_cosThetaProdPlane.png");

  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  //after cuts
  TF1 *fit_L0_L0_after_cuts = new TF1("fit_L0_L0_after_cuts", "[0]*(1+[1]*x)", -1, 1);
  fit_L0_L0_after_cuts->SetParameters(1000, 0.10);


  L0_L0_cosThetaProdPlane_cuts_can->cd();

  float nL0L0_cuts = L0_L0_cosThetaProdPlane_cuts->Integral();

  L0_L0_cosThetaProdPlane_cuts->Sumw2();
  L0_L0_cosThetaProdPlane_cuts->SetMarkerStyle(20);
  L0_L0_cosThetaProdPlane_cuts->SetMarkerColor(kRed);
  L0_L0_cosThetaProdPlane_cuts->SetLineColor(kRed);
  L0_L0_cosThetaProdPlane_cuts->GetXaxis()->SetTitle("cos(#theta*)");
  L0_L0_cosThetaProdPlane_cuts->GetXaxis()->CenterTitle();
  L0_L0_cosThetaProdPlane_cuts->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0_L0_cosThetaProdPlane_cuts->GetYaxis()->CenterTitle();
  //L0_L0_cosThetaProdPlane_cuts->Scale(nL0L0_cuts/L0_L0_cosThetaProdPlane_cuts->Integral()); //scale back
  //L0_L0_cosThetaProdPlane_cuts->Scale(1./L0_L0_cosThetaProdPlane_cuts->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
  L0_L0_cosThetaProdPlane_cuts->Fit(fit_L0_L0_after_cuts, "i 0 r");
  L0_L0_cosThetaProdPlane_cuts->SetMinimum(0);
  //L0_L0_cosThetaProdPlane_cuts->GetYaxis()->SetRangeUser(-20e9, 0);
  L0_L0_cosThetaProdPlane_cuts->Draw("p e");

  fit_L0_L0_after_cuts->SetLineColor(1);
  fit_L0_L0_after_cuts->Draw("same");

  float L0_L0_slope_cuts = fit_L0_L0_after_cuts->GetParameter(1)/(L0_alpha*L0_alpha);
  float L0_L0_slope_cuts_err = fit_L0_L0_after_cuts->GetParError(1)/(L0_alpha*L0_alpha);


  TPaveText *L0_L0_text_MC_cuts = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
  L0_L0_text_MC_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  L0_L0_text_MC_cuts->AddText("Minimum bias");
  L0_L0_text_MC_cuts->AddText("#Lambda-#Lambda");
  L0_L0_text_MC_cuts->AddText("True MC");
  L0_L0_text_MC_cuts->AddText("Analysis cuts");
  L0_L0_text_MC_cuts->AddText(Form("P = %.3f #pm %.3f", L0_L0_slope_cuts, L0_L0_slope_cuts_err));
  L0_L0_text_MC_cuts->SetFillColorAlpha(0, 0.01);
  L0_L0_text_MC_cuts->Draw("same");

  L0_L0_cosThetaProdPlane_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda_mom_smear/L0_L0_cosThetaProdPlane_cuts.png");


  //_________________________________________________________________________

  //Lbar-Lbar
  //before cuts
  TF1 *fit_L0bar_L0bar_before_cuts = new TF1("fit_L0bar_L0bar_before_cuts", "[0]*(1+[1]*x)", -1, 1);
  fit_L0bar_L0bar_before_cuts->SetParameters(1000, 0.10);

  L0bar_L0bar_cosThetaProdPlane_can->cd();

  float nL0barL0bar = L0bar_L0bar_cosThetaProdPlane->Integral();


  L0bar_L0bar_cosThetaProdPlane->Sumw2();
  L0bar_L0bar_cosThetaProdPlane->SetMarkerStyle(20);
  L0bar_L0bar_cosThetaProdPlane->SetMarkerColor(kRed);
  L0bar_L0bar_cosThetaProdPlane->SetLineColor(kRed);
  L0bar_L0bar_cosThetaProdPlane->GetXaxis()->SetTitle("cos(#theta*)");
  L0bar_L0bar_cosThetaProdPlane->GetXaxis()->CenterTitle();
  L0bar_L0bar_cosThetaProdPlane->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0bar_L0bar_cosThetaProdPlane->GetYaxis()->CenterTitle();
  L0bar_L0bar_cosThetaProdPlane->Scale(nL0barL0bar/L0bar_L0bar_cosThetaProdPlane->Integral()); //scale back
  //L0bar_L0bar_cosThetaProdPlane->Scale(1./L0bar_L0bar_cosThetaProdPlane->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
  L0bar_L0bar_cosThetaProdPlane->Fit(fit_L0bar_L0bar_before_cuts, "i 0 r");
  L0bar_L0bar_cosThetaProdPlane->SetMinimum(0);
  //L0bar_L0bar_cosThetaProdPlane->GetYaxis()->SetRangeUser(-20e9, 0);
  L0bar_L0bar_cosThetaProdPlane->Draw("p e");

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

  L0bar_L0bar_cosThetaProdPlane_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda_mom_smear/L0bar_L0bar_cosThetaProdPlane.png");

  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  //after cuts
  TF1 *fit_L0bar_L0bar_after_cuts = new TF1("fit_L0bar_L0bar_after_cuts", "[0]*(1+[1]*x)", -1, 1);
  fit_L0bar_L0bar_after_cuts->SetParameters(1000, 0.10);


  L0bar_L0bar_cosThetaProdPlane_cuts_can->cd();

  float nL0barL0bar_cuts = L0bar_L0bar_cosThetaProdPlane_cuts->Integral();

  L0bar_L0bar_cosThetaProdPlane_cuts->Sumw2();
  L0bar_L0bar_cosThetaProdPlane_cuts->SetMarkerStyle(20);
  L0bar_L0bar_cosThetaProdPlane_cuts->SetMarkerColor(kRed);
  L0bar_L0bar_cosThetaProdPlane_cuts->SetLineColor(kRed);
  L0bar_L0bar_cosThetaProdPlane_cuts->GetXaxis()->SetTitle("cos(#theta*)");
  L0bar_L0bar_cosThetaProdPlane_cuts->GetXaxis()->CenterTitle();
  L0bar_L0bar_cosThetaProdPlane_cuts->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0bar_L0bar_cosThetaProdPlane_cuts->GetYaxis()->CenterTitle();
  //L0bar_L0bar_cosThetaProdPlane_cuts->Scale(nL0barL0bar_cuts/L0bar_L0bar_cosThetaProdPlane_cuts->Integral()); //scale back
  //L0bar_L0bar_cosThetaProdPlane_cuts->Scale(1./L0bar_L0bar_cosThetaProdPlane_cuts->GetXaxis()->GetBinWidth(1)); //scale by bin widhth
  L0bar_L0bar_cosThetaProdPlane_cuts->Fit(fit_L0bar_L0bar_after_cuts, "i 0 r");
  L0bar_L0bar_cosThetaProdPlane_cuts->SetMinimum(0);
  //L0bar_L0bar_cosThetaProdPlane_cuts->GetYaxis()->SetRangeUser(-20e9, 0);
  L0bar_L0bar_cosThetaProdPlane_cuts->Draw("p e");

  fit_L0bar_L0bar_after_cuts->SetLineColor(1);
  fit_L0bar_L0bar_after_cuts->Draw("same");

  float L0bar_L0bar_slope_cuts = fit_L0bar_L0bar_after_cuts->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
  float L0bar_L0bar_slope_cuts_err = fit_L0bar_L0bar_after_cuts->GetParError(1)/(L0bar_alpha*L0bar_alpha);


  TPaveText *L0bar_L0bar_text_MC_cuts = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
  L0bar_L0bar_text_MC_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  L0bar_L0bar_text_MC_cuts->AddText("Minimum bias");
  L0bar_L0bar_text_MC_cuts->AddText("#bar{#Lambda}-#bar{#Lambda}");
  L0bar_L0bar_text_MC_cuts->AddText("True MC");
  L0bar_L0bar_text_MC_cuts->AddText("Analysis cuts");
  L0bar_L0bar_text_MC_cuts->AddText(Form("P = %.3f #pm %.3f", L0bar_L0bar_slope_cuts, L0bar_L0bar_slope_cuts_err));
  L0bar_L0bar_text_MC_cuts->SetFillColorAlpha(0, 0.01);
  L0bar_L0bar_text_MC_cuts->Draw("same");

  L0bar_L0bar_cosThetaProdPlane_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/Lambda_mom_smear/L0bar_L0bar_cosThetaProdPlane_cuts.png");

  //_____________________________________________________________________________________

  inFile->Close();

  return;
}
