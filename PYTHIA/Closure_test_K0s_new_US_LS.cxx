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
void Closure_test_K0s_new_US_LS(const int energy = 510)
{

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
    inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/output_K0s_pp_200_MB_1B_events_hists_work.root", "READ");
  }
  else
  {
    cout<<"Not a valid collision energy! Abborting!"<<endl;
    return;
  }



  //histograms

  //true MC
  TH1D *K0s_K0s_cosThetaProdPlane_hist = (TH1D*)inFile->Get("K0s_K0s_cosThetaProdPlane");

  TH1D *K0s_K0s_cosThetaProdPlane_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1D *K0s_K0s_cosThetaProdPlane_eta_hist[nEtaBins][nEtaBins];


  TH1D *K0s_K0s_cosThetaProdPlane_cuts_hist = (TH1D*)inFile->Get("K0s_K0s_cosThetaProdPlane_cuts");

  TH1D *K0s_K0s_cosThetaProdPlane_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1D *K0s_K0s_cosThetaProdPlane_eta_cuts_hist[nEtaBins][nEtaBins];

  //mixed event
  TH1D *K0s_K0s_cosThetaProdPlane_ME_hist = (TH1D*)inFile->Get("K0s_K0s_cosThetaProdPlane_ME");

  TH1D *K0s_K0s_cosThetaProdPlane_ME_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1D *K0s_K0s_cosThetaProdPlane_ME_eta_hist[nEtaBins][nEtaBins];


  TH1D *K0s_K0s_cosThetaProdPlane_ME_cuts_hist = (TH1D*)inFile->Get("K0s_K0s_cosThetaProdPlane_ME_cuts");

  TH1D *K0s_K0s_cosThetaProdPlane_ME_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1D *K0s_K0s_cosThetaProdPlane_ME_eta_cuts_hist[nEtaBins][nEtaBins];

  //---------------------------------------------------------------------------------------------------------

  //need to update hist names from file with new production
  //US matched to US
  TH1D *K0s_K0s_cosThetaProdPlane_US_hist = (TH1D*)inFile->Get("K0s_K0s_cosThetaProdPlane_US");

  TH1D *K0s_K0s_cosThetaProdPlane_US_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1D *K0s_K0s_cosThetaProdPlane_US_eta_hist[nEtaBins][nEtaBins];


  TH1D *K0s_K0s_cosThetaProdPlane_US_cuts_hist = (TH1D*)inFile->Get("K0s_K0s_cosThetaProdPlane_US_cuts");

  TH1D *K0s_K0s_cosThetaProdPlane_US_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1D *K0s_K0s_cosThetaProdPlane_US_eta_cuts_hist[nEtaBins][nEtaBins];

  //mixed event
  TH1D *K0s_K0s_cosThetaProdPlane_US_ME_hist = (TH1D*)inFile->Get("K0s_K0s_cosThetaProdPlane_US_ME");

  TH1D *K0s_K0s_cosThetaProdPlane_US_ME_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1D *K0s_K0s_cosThetaProdPlane_US_ME_eta_hist[nEtaBins][nEtaBins];


  TH1D *K0s_K0s_cosThetaProdPlane_US_ME_cuts_hist = (TH1D*)inFile->Get("K0s_K0s_cosThetaProdPlane_US_ME_cuts");

  TH1D *K0s_K0s_cosThetaProdPlane_US_ME_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1D *K0s_K0s_cosThetaProdPlane_US_ME_eta_cuts_hist[nEtaBins][nEtaBins];

  //---------------------------------------------------------------------------------------------------------

  //US mathced to LS
  TH1D *K0s_K0s_cosThetaProdPlane_US_LS_hist = (TH1D*)inFile->Get("K0s_K0s_cosThetaProdPlane_US_LS");

  TH1D *K0s_K0s_cosThetaProdPlane_US_LS_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1D *K0s_K0s_cosThetaProdPlane_US_LS_eta_hist[nEtaBins][nEtaBins];


  TH1D *K0s_K0s_cosThetaProdPlane_US_LS_cuts_hist = (TH1D*)inFile->Get("K0s_K0s_cosThetaProdPlane_US_LS_cuts");

  TH1D *K0s_K0s_cosThetaProdPlane_US_LS_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1D *K0s_K0s_cosThetaProdPlane_US_LS_eta_cuts_hist[nEtaBins][nEtaBins];

  //mixed event
  TH1D *K0s_K0s_cosThetaProdPlane_US_LS_ME_hist = (TH1D*)inFile->Get("K0s_K0s_cosThetaProdPlane_US_LS_ME");

  TH1D *K0s_K0s_cosThetaProdPlane_US_LS_ME_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1D *K0s_K0s_cosThetaProdPlane_US_LS_ME_eta_hist[nEtaBins][nEtaBins];


  TH1D *K0s_K0s_cosThetaProdPlane_US_LS_ME_cuts_hist = (TH1D*)inFile->Get("K0s_K0s_cosThetaProdPlane_US_LS_ME_cuts");

  TH1D *K0s_K0s_cosThetaProdPlane_US_LS_ME_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1D *K0s_K0s_cosThetaProdPlane_US_LS_ME_eta_cuts_hist[nEtaBins][nEtaBins];

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

  K0s_K0s_cosThetaProdPlane_hist->SetMarkerStyle(20);
  K0s_K0s_cosThetaProdPlane_hist->SetMarkerColor(kRed);
  K0s_K0s_cosThetaProdPlane_hist->SetLineColor(kRed);
  K0s_K0s_cosThetaProdPlane_hist->GetXaxis()->SetTitle("cos(#theta*)");
  K0s_K0s_cosThetaProdPlane_hist->GetXaxis()->CenterTitle();
  K0s_K0s_cosThetaProdPlane_hist->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  K0s_K0s_cosThetaProdPlane_hist->GetYaxis()->CenterTitle();
  K0s_K0s_cosThetaProdPlane_hist->Sumw2();
  K0s_K0s_cosThetaProdPlane_hist->Divide(K0s_K0s_cosThetaProdPlane_ME_hist);
  K0s_K0s_cosThetaProdPlane_hist->Scale(nK0sK0s/K0s_K0s_cosThetaProdPlane_hist->Integral());
  //K0s_K0s_cosThetaProdPlane_hist->Scale(1./K0s_K0s_cosThetaProdPlane_hist->GetXaxis()->GetBinWidth(1));
  K0s_K0s_cosThetaProdPlane_hist->SetMinimum(0);
  K0s_K0s_cosThetaProdPlane_hist->Draw("p e");


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

  //after analysis cuts
  TCanvas *K0s_K0s_cosThetaProdPlane_cuts_can = new TCanvas("K0s_K0s_cosThetaProdPlane_cuts_can", "K0s_K0s_cosThetaProdPlane_cuts_can", 1200, 1000);

  float nK0sK0s_cuts = K0s_K0s_cosThetaProdPlane_cuts_hist->Integral();

  K0s_K0s_cosThetaProdPlane_ME_cuts_hist->Sumw2();
  K0s_K0s_cosThetaProdPlane_ME_cuts_hist->Scale(nK0sK0s_cuts/K0s_K0s_cosThetaProdPlane_ME_cuts_hist->Integral());

  K0s_K0s_cosThetaProdPlane_cuts_hist->SetMarkerStyle(20);
  K0s_K0s_cosThetaProdPlane_cuts_hist->SetMarkerColor(kRed);
  K0s_K0s_cosThetaProdPlane_cuts_hist->SetLineColor(kRed);
  K0s_K0s_cosThetaProdPlane_cuts_hist->GetXaxis()->SetTitle("cos(#theta*)");
  K0s_K0s_cosThetaProdPlane_cuts_hist->GetXaxis()->CenterTitle();
  K0s_K0s_cosThetaProdPlane_cuts_hist->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  K0s_K0s_cosThetaProdPlane_cuts_hist->GetYaxis()->CenterTitle();
  K0s_K0s_cosThetaProdPlane_cuts_hist->Sumw2();
  K0s_K0s_cosThetaProdPlane_cuts_hist->Divide(K0s_K0s_cosThetaProdPlane_ME_cuts_hist); //correct using ME
  K0s_K0s_cosThetaProdPlane_cuts_hist->Scale(nK0sK0s_cuts/K0s_K0s_cosThetaProdPlane_cuts_hist->Integral()); //scale back
  //K0s_K0s_cosThetaProdPlane_cuts_hist->Scale(1./K0s_K0s_cosThetaProdPlane_cuts_hist->GetXaxis()->GetBinWidth(1)); //bin width
  K0s_K0s_cosThetaProdPlane_cuts_hist->SetMinimum(0);
  K0s_K0s_cosThetaProdPlane_cuts_hist->Draw("p e");

  //K0s_K0s_cosThetaProdPlane_ME_cuts_hist->Draw("same");

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

  //_____________________________________________________________________________________________________________________________________________________________________________________

  //before cuts
  TCanvas *K0s_K0s_cosThetaProdPlane_US_can = new TCanvas("K0s_K0s_cosThetaProdPlane_US_can", "K0s_K0s_cosThetaProdPlane_US_can", 1200, 1000);

  //background
  float nK0sK0s_US_LS = K0s_K0s_cosThetaProdPlane_US_LS_hist->Integral();

  K0s_K0s_cosThetaProdPlane_US_LS_ME_hist->Sumw2();
  K0s_K0s_cosThetaProdPlane_US_LS_ME_hist->Scale(nK0sK0s_US_LS/K0s_K0s_cosThetaProdPlane_US_LS_ME_hist->Integral()); //scale ME to mach data

  K0s_K0s_cosThetaProdPlane_US_LS_hist->Sumw2();
  K0s_K0s_cosThetaProdPlane_US_LS_hist->Divide(K0s_K0s_cosThetaProdPlane_US_LS_ME_hist); //correction using ME
  K0s_K0s_cosThetaProdPlane_US_LS_hist->Scale(nK0sK0s_US_LS/K0s_K0s_cosThetaProdPlane_US_LS_hist->Integral()); //scale back

  //signal+background
  float nK0sK0s_US = K0s_K0s_cosThetaProdPlane_US_hist->Integral();

  K0s_K0s_cosThetaProdPlane_US_ME_hist->Sumw2();
  K0s_K0s_cosThetaProdPlane_US_ME_hist->Scale(nK0sK0s_US/K0s_K0s_cosThetaProdPlane_US_ME_hist->Integral());

  K0s_K0s_cosThetaProdPlane_US_hist->SetMarkerStyle(20);
  K0s_K0s_cosThetaProdPlane_US_hist->SetMarkerColor(kRed);
  K0s_K0s_cosThetaProdPlane_US_hist->SetLineColor(kRed);
  K0s_K0s_cosThetaProdPlane_US_hist->GetXaxis()->SetTitle("cos(#theta*)");
  K0s_K0s_cosThetaProdPlane_US_hist->GetXaxis()->CenterTitle();
  K0s_K0s_cosThetaProdPlane_US_hist->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  K0s_K0s_cosThetaProdPlane_US_hist->GetYaxis()->CenterTitle();
  K0s_K0s_cosThetaProdPlane_US_hist->Sumw2();
  K0s_K0s_cosThetaProdPlane_US_hist->Divide(K0s_K0s_cosThetaProdPlane_US_ME_hist);
  K0s_K0s_cosThetaProdPlane_US_hist->Scale(nK0sK0s_US/K0s_K0s_cosThetaProdPlane_US_hist->Integral());
  K0s_K0s_cosThetaProdPlane_US_hist->Add(K0s_K0s_cosThetaProdPlane_US_LS_hist, -1); //subtract background
  K0s_K0s_cosThetaProdPlane_US_hist->Scale(1./K0s_K0s_cosThetaProdPlane_US_hist->GetXaxis()->GetBinWidth(1));
  K0s_K0s_cosThetaProdPlane_US_hist->SetMinimum(0);
  K0s_K0s_cosThetaProdPlane_US_hist->Draw("p e");


  TPaveText *K0s_K0s_text = new TPaveText(0.4, 0.2, 0.8, 0.4, "NDC");
  //K0s_K0s_text->AddText("STAR preliminary");
  //((TText*)K0s_K0s_text->GetListOfLines()->Last())->SetTextColor(2);
  K0s_K0s_text->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  K0s_K0s_text->AddText("Minimum bias");
  K0s_K0s_text->AddText("K_{s}^{0}-K_{s}^{0}");
  K0s_K0s_text->AddText("US - Background");
  K0s_K0s_text->AddText("No cuts");
  K0s_K0s_text->SetFillColorAlpha(0, 0.01);
  K0s_K0s_text->Draw("same");

  K0s_K0s_cosThetaProdPlane_US_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/K0s/K0s_K0s_cosThetaProdPlane_US.png");

  //----------------------------------------------------------------------------------------------------------

  //after analysis cuts
  TCanvas *K0s_K0s_cosThetaProdPlane_US_cuts_can = new TCanvas("K0s_K0s_cosThetaProdPlane_US_cuts_can", "K0s_K0s_cosThetaProdPlane_US_cuts_can", 1200, 1000);

  //background
  float nK0sK0s_US_LS_cuts = K0s_K0s_cosThetaProdPlane_US_LS_cuts_hist->Integral();

  K0s_K0s_cosThetaProdPlane_US_LS_ME_cuts_hist->Sumw2();
  K0s_K0s_cosThetaProdPlane_US_LS_ME_cuts_hist->Scale(nK0sK0s_US_LS_cuts/K0s_K0s_cosThetaProdPlane_US_LS_ME_cuts_hist->Integral());

  K0s_K0s_cosThetaProdPlane_US_LS_cuts_hist->Sumw2();
  K0s_K0s_cosThetaProdPlane_US_LS_cuts_hist->Divide(K0s_K0s_cosThetaProdPlane_US_LS_ME_cuts_hist);
  K0s_K0s_cosThetaProdPlane_US_LS_cuts_hist->Scale(nK0sK0s_US_LS_cuts/K0s_K0s_cosThetaProdPlane_US_LS_cuts_hist->Integral());


  //signal+background
  float nK0sK0s_US_cuts = K0s_K0s_cosThetaProdPlane_US_cuts_hist->Integral();

  K0s_K0s_cosThetaProdPlane_US_ME_cuts_hist->Sumw2();
  K0s_K0s_cosThetaProdPlane_US_ME_cuts_hist->Scale(nK0sK0s_US_cuts/K0s_K0s_cosThetaProdPlane_US_ME_cuts_hist->Integral());

  K0s_K0s_cosThetaProdPlane_US_cuts_hist->SetMarkerStyle(20);
  K0s_K0s_cosThetaProdPlane_US_cuts_hist->SetMarkerColor(kRed);
  K0s_K0s_cosThetaProdPlane_US_cuts_hist->SetLineColor(kRed);
  K0s_K0s_cosThetaProdPlane_US_cuts_hist->GetXaxis()->SetTitle("cos(#theta*)");
  K0s_K0s_cosThetaProdPlane_US_cuts_hist->GetXaxis()->CenterTitle();
  K0s_K0s_cosThetaProdPlane_US_cuts_hist->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  K0s_K0s_cosThetaProdPlane_US_cuts_hist->GetYaxis()->CenterTitle();
  K0s_K0s_cosThetaProdPlane_US_cuts_hist->Sumw2();
  K0s_K0s_cosThetaProdPlane_US_cuts_hist->Divide(K0s_K0s_cosThetaProdPlane_US_ME_cuts_hist); //correct using ME
  K0s_K0s_cosThetaProdPlane_US_cuts_hist->Scale(nK0sK0s_US_cuts/K0s_K0s_cosThetaProdPlane_US_cuts_hist->Integral()); //scale back
  K0s_K0s_cosThetaProdPlane_US_cuts_hist->Add(K0s_K0s_cosThetaProdPlane_US_LS_cuts_hist, -1); //subtract background
  //K0s_K0s_cosThetaProdPlane_US_cuts_hist->Scale(1./K0s_K0s_cosThetaProdPlane_US_cuts_hist->GetXaxis()->GetBinWidth(1)); //bin width
  K0s_K0s_cosThetaProdPlane_US_cuts_hist->SetMinimum(0);
  K0s_K0s_cosThetaProdPlane_US_cuts_hist->Draw("p e");

  //K0s_K0s_cosThetaProdPlane_US_ME_cuts_hist->Draw("same");

  TPaveText *K0s_K0s_text_cuts = new TPaveText(0.4, 0.2, 0.8, 0.4, "NDC");
  //K0s_K0s_text_cuts->AddText("STAR preliminary");
  //((TText*)K0s_K0s_text_cuts->GetListOfLines()->Last())->SetTextColor(2);
  K0s_K0s_text_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  K0s_K0s_text_cuts->AddText("Minimum bias");
  K0s_K0s_text_cuts->AddText("K_{s}^{0}-K_{s}^{0}");
  K0s_K0s_text_cuts->AddText("US - Background");
  K0s_K0s_text_cuts->AddText("Analysis cuts");
  K0s_K0s_text_cuts->SetFillColorAlpha(0, 0.01);
  K0s_K0s_text_cuts->Draw("same");

  K0s_K0s_cosThetaProdPlane_US_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/K0s/K0s_K0s_cosThetaProdPlane_US_cuts.png");

  //_____________________________________________________________________________________________________________________________________________________________________________________


  float nK0sK0s_US_pT[nPtBins_corr][nPtBins_corr];
  float nK0sK0s_US_LS_pT[nPtBins_corr][nPtBins_corr];

  float nK0sK0s_US_pT_cuts[nPtBins_corr][nPtBins_corr];
  float nK0sK0s_US_LS_pT_cuts[nPtBins_corr][nPtBins_corr];

  //K0s-K0s correlation histograms in bins
  for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
  {
    for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
    {
      //before cuts
      K0s_K0s_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2] = (TH1D*)inFile->Get(Form("K0s_K0s_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
      K0s_K0s_cosThetaProdPlane_US_ME_pT_hist[pTbin1][pTbin2] = (TH1D*)inFile->Get(Form("K0s_K0s_cosThetaProdPlane_US_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));

      K0s_K0s_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2] = (TH1D*)inFile->Get(Form("K0s_K0s_cosThetaProdPlane_US_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
      K0s_K0s_cosThetaProdPlane_US_LS_ME_pT_hist[pTbin1][pTbin2] = (TH1D*)inFile->Get(Form("K0s_K0s_cosThetaProdPlane_US_LS_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));


      TCanvas *K0s_K0s_cosThetaProdPlane_US_pT_can = new TCanvas(Form("K0s_K0s_cosThetaProdPlane_US_pT_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("K0s_K0s_cosThetaProdPlane_US_pT_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);

      //background
      nK0sK0s_US_pT[pTbin1][pTbin2] = K0s_K0s_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2]->Integral();

      K0s_K0s_cosThetaProdPlane_US_LS_ME_pT_hist[pTbin1][pTbin2]->Sumw2();
      K0s_K0s_cosThetaProdPlane_US_LS_ME_pT_hist[pTbin1][pTbin2]->Scale(nK0sK0s_US_pT[pTbin1][pTbin2]/K0s_K0s_cosThetaProdPlane_US_LS_ME_pT_hist[pTbin1][pTbin2]->Integral());

      K0s_K0s_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2]->Sumw2();
      K0s_K0s_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2]->Divide(K0s_K0s_cosThetaProdPlane_US_LS_ME_pT_hist[pTbin1][pTbin2]);
      K0s_K0s_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2]->Scale(nK0sK0s_US_pT[pTbin1][pTbin2]/K0s_K0s_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2]->Integral());

      //signal+background
      nK0sK0s_US_LS_pT[pTbin1][pTbin2] = K0s_K0s_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->Integral();

      K0s_K0s_cosThetaProdPlane_US_ME_pT_hist[pTbin1][pTbin2]->Sumw2();
      K0s_K0s_cosThetaProdPlane_US_ME_pT_hist[pTbin1][pTbin2]->Scale(nK0sK0s_US_LS_pT[pTbin1][pTbin2]/K0s_K0s_cosThetaProdPlane_US_ME_pT_hist[pTbin1][pTbin2]->Integral());

      K0s_K0s_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->Sumw2();
      K0s_K0s_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      K0s_K0s_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->SetMarkerColor(kRed);
      K0s_K0s_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->SetLineColor(kRed);
      K0s_K0s_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->GetXaxis()->SetTitle("cos(#theta*)");
      K0s_K0s_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      K0s_K0s_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      K0s_K0s_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      K0s_K0s_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->Divide(K0s_K0s_cosThetaProdPlane_US_ME_pT_hist[pTbin1][pTbin2]);
      K0s_K0s_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->Scale(nK0sK0s_US_LS_pT[pTbin1][pTbin2]/K0s_K0s_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->Integral());
      K0s_K0s_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->Add(K0s_K0s_cosThetaProdPlane_US_LS_pT_hist[pTbin1][pTbin2], -1);
      K0s_K0s_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->Scale(1./K0s_K0s_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->GetXaxis()->GetBinWidth(1));
      K0s_K0s_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->SetMinimum(0);
      K0s_K0s_cosThetaProdPlane_US_pT_hist[pTbin1][pTbin2]->Draw("p e");


      TPaveText *K0s_K0s_text_pT = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
      K0s_K0s_text_pT->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
      K0s_K0s_text_pT->AddText("Minimum bias");
      K0s_K0s_text_pT->AddText("K_{s}^{0}-K_{s}^{0}");
      K0s_K0s_text_pT->AddText("US - Background");
      K0s_K0s_text_pT->AddText("No cuts");
      K0s_K0s_text_pT->AddText("|#eta| < 1");
      K0s_K0s_text_pT->AddText(Form("%0.1f < p_{T}^{1} < %0.1f", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      K0s_K0s_text_pT->AddText(Form("%0.1f < p_{T}^{2} < %0.1f", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      K0s_K0s_text_pT->SetFillColorAlpha(0, 0.01);
      K0s_K0s_text_pT->Draw("same");

      K0s_K0s_cosThetaProdPlane_US_pT_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/K0s/K0s_K0s_cosThetaProdPlane_US_pT1_%i_pT2_%i.png", pTbin1, pTbin2));

      //----------------------------------------------------------------------------------------------------------

      //after cuts
      K0s_K0s_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2] = (TH1D*)inFile->Get(Form("K0s_K0s_cosThetaProdPlane_US_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
      K0s_K0s_cosThetaProdPlane_US_ME_pT_cuts_hist[pTbin1][pTbin2] = (TH1D*)inFile->Get(Form("K0s_K0s_cosThetaProdPlane_US_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));

      K0s_K0s_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2] = (TH1D*)inFile->Get(Form("K0s_K0s_cosThetaProdPlane_US_LS_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
      K0s_K0s_cosThetaProdPlane_US_LS_ME_pT_cuts_hist[pTbin1][pTbin2] = (TH1D*)inFile->Get(Form("K0s_K0s_cosThetaProdPlane_US_LS_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));

      TCanvas *K0s_K0s_cosThetaProdPlane_US_pT_cuts_can = new TCanvas(Form("K0s_K0s_cosThetaProdPlane_US_pT_cuts_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("K0s_K0s_cosThetaProdPlane_US_pT_cuts_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);

      //background
      nK0sK0s_US_LS_pT_cuts[pTbin1][pTbin2] = K0s_K0s_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2]->Integral();

      K0s_K0s_cosThetaProdPlane_US_LS_ME_pT_cuts_hist[pTbin1][pTbin2]->Sumw2();
      K0s_K0s_cosThetaProdPlane_US_LS_ME_pT_cuts_hist[pTbin1][pTbin2]->Scale(nK0sK0s_US_LS_pT_cuts[pTbin1][pTbin2]/K0s_K0s_cosThetaProdPlane_US_LS_ME_pT_cuts_hist[pTbin1][pTbin2]->Integral());

      K0s_K0s_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2]->Sumw2();
      K0s_K0s_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2]->Divide(K0s_K0s_cosThetaProdPlane_US_LS_ME_pT_cuts_hist[pTbin1][pTbin2]);
      K0s_K0s_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2]->Scale(1./K0s_K0s_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2]->Integral());

      //signal+background
      nK0sK0s_US_pT_cuts[pTbin1][pTbin2] = K0s_K0s_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->Integral();

      K0s_K0s_cosThetaProdPlane_US_ME_pT_cuts_hist[pTbin1][pTbin2]->Sumw2();
      K0s_K0s_cosThetaProdPlane_US_ME_pT_cuts_hist[pTbin1][pTbin2]->Scale(nK0sK0s_US_pT_cuts[pTbin1][pTbin2]/K0s_K0s_cosThetaProdPlane_US_ME_pT_cuts_hist[pTbin1][pTbin2]->Integral());

      K0s_K0s_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      K0s_K0s_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->SetMarkerColor(kRed);
      K0s_K0s_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->SetLineColor(kRed);
      K0s_K0s_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->GetXaxis()->SetTitle("cos(#theta*)");
      K0s_K0s_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      K0s_K0s_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      K0s_K0s_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      K0s_K0s_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->Divide(K0s_K0s_cosThetaProdPlane_US_ME_pT_cuts_hist[pTbin1][pTbin2]);
      K0s_K0s_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->Scale(1./K0s_K0s_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->Integral());
      K0s_K0s_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->Add(K0s_K0s_cosThetaProdPlane_US_LS_pT_cuts_hist[pTbin1][pTbin2], -1);
      K0s_K0s_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->Scale(1./K0s_K0s_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->GetXaxis()->GetBinWidth(1));
      K0s_K0s_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->SetMinimum(0);
      K0s_K0s_cosThetaProdPlane_US_pT_cuts_hist[pTbin1][pTbin2]->Draw("p e");


      TPaveText *K0s_K0s_text_pT_cuts = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
      K0s_K0s_text_pT_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
      K0s_K0s_text_pT_cuts->AddText("Minimum bias");
      K0s_K0s_text_pT_cuts->AddText("K_{s}^{0}-K_{s}^{0}");
      K0s_K0s_text_pT_cuts->AddText("US - Background");
      K0s_K0s_text_pT_cuts->AddText("Analysis cuts");
      K0s_K0s_text_pT_cuts->AddText("|#eta| < 1");
      K0s_K0s_text_pT_cuts->AddText(Form("%0.1f < p_{T}^{1} < %0.1f", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      K0s_K0s_text_pT_cuts->AddText(Form("%0.1f < p_{T}^{2} < %0.1f", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      K0s_K0s_text_pT_cuts->SetFillColorAlpha(0, 0.01);
      K0s_K0s_text_pT_cuts->Draw("same");

      K0s_K0s_cosThetaProdPlane_US_pT_cuts_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/K0s/K0s_K0s_cosThetaProdPlane_US_cuts_pT1_%i_pT2_%i_eff.png", pTbin1, pTbin2));
      //_________________________________________________________________________________________________

    }
  }


  float nK0sK0s_US_eta[nEtaBins][nEtaBins];
  float nK0sK0s_US_LS_eta[nEtaBins][nEtaBins];

  float nK0sK0s_US_eta_cuts[nEtaBins][nEtaBins];
  float nK0sK0s_US_LS_eta_cuts[nEtaBins][nEtaBins];

  for(unsigned int etaBin1 = 0; etaBin1 < nEtaBins; etaBin1++)
  {
    for(unsigned int etaBin2 = 0; etaBin2 < nEtaBins; etaBin2++)
    {
      //before cuts
      K0s_K0s_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2] = (TH1D*)inFile->Get(Form("K0s_K0s_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));
      K0s_K0s_cosThetaProdPlane_US_ME_eta_hist[etaBin1][etaBin2] = (TH1D*)inFile->Get(Form("K0s_K0s_cosThetaProdPlane_US_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));

      K0s_K0s_cosThetaProdPlane_US_LS_eta_hist[etaBin1][etaBin2] = (TH1D*)inFile->Get(Form("K0s_K0s_cosThetaProdPlane_US_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));
      K0s_K0s_cosThetaProdPlane_US_LS_ME_eta_hist[etaBin1][etaBin2] = (TH1D*)inFile->Get(Form("K0s_K0s_cosThetaProdPlane_US_LS_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));


      TCanvas *K0s_K0s_cosThetaProdPlane_US_eta_can = new TCanvas(Form("K0s_K0s_cosThetaProdPlane_US_eta_can_eta1_%i_eta2_%i", etaBin1, etaBin2), Form("K0s_K0s_cosThetaProdPlane_US_eta_can_eta1_%i_eta2_%i", etaBin1, etaBin2), 1200, 1000);

      //background
      nK0sK0s_US_LS_eta[etaBin1][etaBin2] = K0s_K0s_cosThetaProdPlane_US_LS_eta_hist[etaBin1][etaBin2]->Integral();

      K0s_K0s_cosThetaProdPlane_US_LS_ME_eta_hist[etaBin1][etaBin2]->Sumw2();
      K0s_K0s_cosThetaProdPlane_US_LS_ME_eta_hist[etaBin1][etaBin2]->Scale(nK0sK0s_US_LS_eta[etaBin1][etaBin2]/K0s_K0s_cosThetaProdPlane_US_LS_ME_eta_hist[etaBin1][etaBin2]->Integral());

      K0s_K0s_cosThetaProdPlane_US_LS_eta_hist[etaBin1][etaBin2]->Sumw2();
      K0s_K0s_cosThetaProdPlane_US_LS_eta_hist[etaBin1][etaBin2]->Divide(K0s_K0s_cosThetaProdPlane_US_LS_ME_eta_hist[etaBin1][etaBin2]);
      K0s_K0s_cosThetaProdPlane_US_LS_eta_hist[etaBin1][etaBin2]->Scale(nK0sK0s_US_LS_eta[etaBin1][etaBin2]/K0s_K0s_cosThetaProdPlane_US_LS_eta_hist[etaBin1][etaBin2]->Integral());

      //signal+background
      nK0sK0s_US_eta[etaBin1][etaBin2] = K0s_K0s_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->Integral();

      K0s_K0s_cosThetaProdPlane_US_ME_eta_hist[etaBin1][etaBin2]->Sumw2();
      K0s_K0s_cosThetaProdPlane_US_ME_eta_hist[etaBin1][etaBin2]->Scale(nK0sK0s_US_eta[etaBin1][etaBin2]/K0s_K0s_cosThetaProdPlane_US_ME_eta_hist[etaBin1][etaBin2]->Integral());

      K0s_K0s_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->Sumw2();
      K0s_K0s_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      K0s_K0s_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->SetMarkerColor(kRed);
      K0s_K0s_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->SetLineColor(kRed);
      K0s_K0s_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->GetXaxis()->SetTitle("cos(#theta*)");
      K0s_K0s_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->GetXaxis()->CenterTitle();
      K0s_K0s_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      K0s_K0s_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->GetYaxis()->CenterTitle();
      K0s_K0s_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->Divide(K0s_K0s_cosThetaProdPlane_US_ME_eta_hist[etaBin1][etaBin2]);
      K0s_K0s_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->Scale(nK0sK0s_US_eta[etaBin1][etaBin2]/K0s_K0s_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->Integral());
      K0s_K0s_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->Add(K0s_K0s_cosThetaProdPlane_US_LS_eta_hist[etaBin1][etaBin2]);
      K0s_K0s_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->Scale(1./K0s_K0s_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->GetXaxis()->GetBinWidth(1));
      K0s_K0s_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->SetMinimum(0);
      K0s_K0s_cosThetaProdPlane_US_eta_hist[etaBin1][etaBin2]->Draw("p e");



      TPaveText *K0s_K0s_text_eta = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
      K0s_K0s_text_eta->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
      K0s_K0s_text_eta->AddText("Minimum bias");
      K0s_K0s_text_eta->AddText("K_{s}^{0}-K_{s}^{0}");
      K0s_K0s_text_eta->AddText("US - Background");
      K0s_K0s_text_eta->AddText("No cuts");
      K0s_K0s_text_eta->AddText(Form("%0.1f < #eta_{1} < %0.1f", eta_bins[etaBin1], eta_bins[etaBin1+1]));
      K0s_K0s_text_eta->AddText(Form("%0.1f < #eta_{2} < %0.1f", eta_bins[etaBin2], eta_bins[etaBin2+1]));
      K0s_K0s_text_eta->AddText("p_{T} integrated");
      K0s_K0s_text_eta->SetFillColorAlpha(0, 0.01);
      K0s_K0s_text_eta->Draw("same");

      K0s_K0s_cosThetaProdPlane_US_eta_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/K0s/K0s_K0s_cosThetaProdPlane_US_eta1_%i_eta2_%i.png", etaBin1, etaBin2));

      //----------------------------------------------------------------------------------------------------------------------------

      //after cuts
      K0s_K0s_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2] = (TH1D*)inFile->Get(Form("K0s_K0s_cosThetaProdPlane_US_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));
      K0s_K0s_cosThetaProdPlane_US_ME_eta_cuts_hist[etaBin1][etaBin2] = (TH1D*)inFile->Get(Form("K0s_K0s_cosThetaProdPlane_US_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));

      K0s_K0s_cosThetaProdPlane_US_LS_eta_cuts_hist[etaBin1][etaBin2] = (TH1D*)inFile->Get(Form("K0s_K0s_cosThetaProdPlane_US_LS_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));
      K0s_K0s_cosThetaProdPlane_US_LS_ME_eta_cuts_hist[etaBin1][etaBin2] = (TH1D*)inFile->Get(Form("K0s_K0s_cosThetaProdPlane_US_LS_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));


      TCanvas *K0s_K0s_cosThetaProdPlane_US_eta_corr_can = new TCanvas(Form("K0s_K0s_cosThetaProdPlane_US_eta_corr_can_eta1_%i_eta2_%i", etaBin1, etaBin2), Form("K0s_K0s_cosThetaProdPlane_US_eta_corr_can_eta1_%i_eta2_%i", etaBin1, etaBin2), 1200, 1000);

      //background
      nK0sK0s_US_LS_eta_cuts[etaBin1][etaBin2] = K0s_K0s_cosThetaProdPlane_US_LS_eta_cuts_hist[etaBin1][etaBin2]->Integral();

      K0s_K0s_cosThetaProdPlane_US_LS_ME_eta_cuts_hist[etaBin1][etaBin2]->Sumw2();
      K0s_K0s_cosThetaProdPlane_US_LS_ME_eta_cuts_hist[etaBin1][etaBin2]->Scale(nK0sK0s_US_LS_eta_cuts[etaBin1][etaBin2]/K0s_K0s_cosThetaProdPlane_US_LS_ME_eta_cuts_hist[etaBin1][etaBin2]->Integral());


      K0s_K0s_cosThetaProdPlane_US_LS_eta_cuts_hist[etaBin1][etaBin2]->Sumw2();
      K0s_K0s_cosThetaProdPlane_US_LS_eta_cuts_hist[etaBin1][etaBin2]->Divide(K0s_K0s_cosThetaProdPlane_US_LS_ME_eta_cuts_hist[etaBin1][etaBin2]);
      K0s_K0s_cosThetaProdPlane_US_LS_eta_cuts_hist[etaBin1][etaBin2]->Scale(nK0sK0s_US_LS_eta_cuts[etaBin1][etaBin2]/K0s_K0s_cosThetaProdPlane_US_LS_eta_cuts_hist[etaBin1][etaBin2]->Integral());


      //signal+background
      nK0sK0s_US_eta_cuts[etaBin1][etaBin2] = K0s_K0s_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->Integral();

      K0s_K0s_cosThetaProdPlane_US_ME_eta_cuts_hist[etaBin1][etaBin2]->Sumw2();
      K0s_K0s_cosThetaProdPlane_US_ME_eta_cuts_hist[etaBin1][etaBin2]->Scale(nK0sK0s_US_eta_cuts[etaBin1][etaBin2]/K0s_K0s_cosThetaProdPlane_US_ME_eta_cuts_hist[etaBin1][etaBin2]->Integral());

      K0s_K0s_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->Sumw2();
      K0s_K0s_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      K0s_K0s_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->SetMarkerColor(kRed);
      K0s_K0s_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->SetLineColor(kRed);
      K0s_K0s_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->GetXaxis()->SetTitle("cos(#theta*)");
      K0s_K0s_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->GetXaxis()->CenterTitle();
      K0s_K0s_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      K0s_K0s_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->GetYaxis()->CenterTitle();
      K0s_K0s_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->Divide(K0s_K0s_cosThetaProdPlane_US_ME_eta_cuts_hist[etaBin1][etaBin2]);
      K0s_K0s_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->Scale(nK0sK0s_US_eta_cuts[etaBin1][etaBin2]/K0s_K0s_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->Integral());
      K0s_K0s_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->Add(K0s_K0s_cosThetaProdPlane_US_LS_eta_cuts_hist[etaBin1][etaBin2], -1);
      K0s_K0s_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->Scale(1./K0s_K0s_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->GetXaxis()->GetBinWidth(1));
      K0s_K0s_cosThetaProdPlane_US_eta_cuts_hist[etaBin1][etaBin2]->SetMinimum(0);


      TPaveText *K0s_K0s_text_eta_cuts = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
      K0s_K0s_text_eta_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
      K0s_K0s_text_eta_cuts->AddText("Minimum bias");
      K0s_K0s_text_eta_cuts->AddText("K_{s}^{0}-K_{s}^{0}");
      K0s_K0s_text_eta_cuts->AddText("US - Background");
      K0s_K0s_text_eta_cuts->AddText("Analysis cuts");
      K0s_K0s_text_eta_cuts->AddText(Form("%0.1f < #eta_{1} < %0.1f", eta_bins[etaBin1], eta_bins[etaBin1+1]));
      K0s_K0s_text_eta_cuts->AddText(Form("%0.1f < #eta_{2} < %0.1f", eta_bins[etaBin2], eta_bins[etaBin2+1]));
      K0s_K0s_text_eta_cuts->AddText("p_{T} integrated");
      K0s_K0s_text_eta_cuts->SetFillColorAlpha(0, 0.01);
      K0s_K0s_text_eta_cuts->Draw("same");

      K0s_K0s_cosThetaProdPlane_US_eta_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/closure_test/K0s/K0s_K0s_cosThetaProdPlane_US_cuts_eta1_%i_eta2_%i_eff.png", etaBin1, etaBin2));
      //_________________________________________________________________________________________________

    }
  }


  inFile->Close();
  //outFile->Close();

  return;

}
