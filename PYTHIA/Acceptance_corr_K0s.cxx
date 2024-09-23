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
void Acceptance_corr_K0s(const int energy = 510)
{

  const int nPtBins = 8;
  float const pT_bins[nPtBins+1] = { 0., 0.5, 1.,1.5, 2., 2.5, 3., 4., 5.};

  const int nPtBins_corr = 2;
  float const pT_bins_corr[nPtBins_corr+1] = { 0.5, 1.5, 5.};

  const int nEtaBins = 3;
  float const eta_bins[nEtaBins+1] = { -1, -0.2, 0.2, 1 };

  TFile *outFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/output/K0s_cosThetaStar_eff_work.root", "recreate");

  //load all files
  TFile *inFile;

  if(energy == 510)
  {
    //inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run17/output_K0s_pp_510_MB_1B_events.root", "READ");
    inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run17/output_K0s_pp_510_MB_1B_events_tight_eta.root", "READ");
  }
  else if(energy == 200)
  {
    //inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/output_K0s_pp_200_MB_1B_events.root", "READ");
    inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/output_K0s_pp_200_MB_1B_events_tight_eta.root", "READ");
  }
  else
  {
    cout<<"Not a valid collision energy! Abborting!"<<endl;
  }



  //histograms
  TH1D *K0s_cosThetaProdPlane_open_hist[nPtBins+1][nEtaBins+1];
  TH1D *K0s_cosThetaProdPlane_ana_hist[nPtBins+1][nEtaBins+1];

  //need tu update hist names from file with new production
  TH1D *K0s_K0s_cosThetaProdPlane_hist = (TH1D*)inFile->Get("K0s_K0s_cosThetaProdPlane");

  TH1D *K0s_K0s_cosThetaProdPlane_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1D *K0s_K0s_cosThetaProdPlane_pT_hist_tight_eta[nPtBins_corr][nPtBins_corr];

  TH1D *K0s_K0s_cosThetaProdPlane_eta_hist[nEtaBins][nEtaBins];


  TH1D *K0s_K0s_cosThetaProdPlane_cuts_hist = (TH1D*)inFile->Get("K0s_K0s_cosThetaProdPlane_cuts");

  TH1D *K0s_K0s_cosThetaProdPlane_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1D *K0s_K0s_cosThetaProdPlane_pT_cuts_hist_tight_eta[nPtBins_corr][nPtBins_corr];

  TH1D *K0s_K0s_cosThetaProdPlane_eta_cuts_hist[nEtaBins][nEtaBins];


  TCanvas *K0s_cosThetaProdPlane_corr_can[nPtBins+1][nEtaBins+1];

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  outFile->cd();

  for(unsigned int pTbin = 0; pTbin < nPtBins+1; pTbin++)
  {

    TString *pT_range = new TString();
    if(pTbin < nPtBins) pT_range->Form("%0.1f < #it{p}_{T} < %0.1f GeV/#it{c}", pT_bins[pTbin], pT_bins[pTbin+1]);
    else pT_range->Form("p_{T} integrated");


    for(unsigned int etaBin = 0; etaBin < nEtaBins+1; etaBin++)
    {



      TString *eta_range = new TString();
      if(etaBin < nEtaBins) eta_range->Form("%0.1f < #eta < %0.1f", eta_bins[etaBin], eta_bins[etaBin+1]);
      else eta_range->Form("-1 < #eta < 1");


      K0s_cosThetaProdPlane_corr_can[pTbin][etaBin] = new TCanvas(Form("K0s_cosThetaProdPlane_corr_can_pT_%i_eta_%i", pTbin, etaBin), Form("K0s_cosThetaProdPlane_corr_can_pT_%i_eta_%i", pTbin, etaBin), 1200, 1000);


      TPaveText *cent_text_2 = new TPaveText(0.4, 0.2, 0.8, 0.4, "NDC");
      //cent_text->AddText("STAR preliminary");
      //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
      cent_text_2->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
      cent_text_2->AddText("Minimum bias");
      cent_text_2->AddText("K_{s}^{0}");
      cent_text_2->AddText(eta_range->Data());
      cent_text_2->AddText(pT_range->Data());
      cent_text_2->SetFillColorAlpha(0, 0.01);
      //cent_text_2->Draw("same");


      K0s_cosThetaProdPlane_corr_can[pTbin][etaBin]->cd();

      K0s_cosThetaProdPlane_open_hist[pTbin][etaBin] = (TH1D*)inFile->Get(Form("K0s_cosThetaProdPlane_pT_%i_eta_%i", pTbin, etaBin));
      K0s_cosThetaProdPlane_open_hist[pTbin][etaBin]->SetNameTitle(Form("K0s_cosThetaProdPlane_open_pT_%i_eta_%i", pTbin, etaBin), Form("K0s_cosThetaProdPlane_open_pT_%i_eta_%i", pTbin, etaBin));

      K0s_cosThetaProdPlane_ana_hist[pTbin][etaBin] = (TH1D*)inFile->Get(Form("K0s_cosThetaProdPlane_cuts_pT_%i_eta_%i", pTbin, etaBin));
      K0s_cosThetaProdPlane_ana_hist[pTbin][etaBin]->SetNameTitle(Form("K0s_cosThetaProdPlane_eff_pT_%i_eta_%i", pTbin, etaBin), Form("K0s_cosThetaProdPlane_eff_pT_%i_eta_%i", pTbin, etaBin));
      K0s_cosThetaProdPlane_ana_hist[pTbin][etaBin]->SetMarkerStyle(20);
      K0s_cosThetaProdPlane_ana_hist[pTbin][etaBin]->SetMarkerColor(kRed);
      K0s_cosThetaProdPlane_ana_hist[pTbin][etaBin]->SetLineColor(kRed);
      K0s_cosThetaProdPlane_ana_hist[pTbin][etaBin]->GetXaxis()->SetTitle("cos(#theta*)");
      K0s_cosThetaProdPlane_ana_hist[pTbin][etaBin]->GetXaxis()->CenterTitle();
      K0s_cosThetaProdPlane_ana_hist[pTbin][etaBin]->GetYaxis()->SetTitle("Acceptance");
      K0s_cosThetaProdPlane_ana_hist[pTbin][etaBin]->GetYaxis()->CenterTitle();
      K0s_cosThetaProdPlane_ana_hist[pTbin][etaBin]->Sumw2();
      K0s_cosThetaProdPlane_ana_hist[pTbin][etaBin]->Divide(K0s_cosThetaProdPlane_ana_hist[pTbin][etaBin], K0s_cosThetaProdPlane_open_hist[pTbin][etaBin], 1, 1, "b"); //binomial errors
      K0s_cosThetaProdPlane_ana_hist[pTbin][etaBin]->SetMinimum(0);
      K0s_cosThetaProdPlane_ana_hist[pTbin][etaBin]->Draw("p e");
      K0s_cosThetaProdPlane_ana_hist[pTbin][etaBin]->Write();


      TLegend *PYTHIA_legend = new TLegend(0.2, 0.4, 0.4, 0.5 );
      PYTHIA_legend->AddEntry(K0s_cosThetaProdPlane_ana_hist[pTbin][etaBin], "PYTHIA 8.162");
      PYTHIA_legend->SetBorderSize(0);
      PYTHIA_legend->SetFillColorAlpha(0, 0.01);
      PYTHIA_legend->Draw("same");


      cent_text_2->Draw("same");


      //_____________________________________________________________________________________________________________________________
      TPaveText *cent_text_3 = new TPaveText(0.4, 0.2, 0.8, 0.4, "NDC");
      //cent_text->AddText("STAR preliminary");
      //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
      cent_text_3->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
      cent_text_3->AddText("Minimum bias");
      cent_text_3->AddText("#bar{#Lambda^{0}}");
      cent_text_3->AddText(eta_range->Data());
      cent_text_3->AddText(pT_range->Data());
      cent_text_3->SetFillColorAlpha(0, 0.01);
      //cent_text_3->Draw("same")


      K0s_cosThetaProdPlane_corr_can[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/efficiency/K0s/K0s_cosThetaProdPlane_eff_pT_%i_eta_%i.png", pTbin, etaBin));


      //save histograms to multi-page PDF
      if(pTbin == 0 && etaBin == 0)
      {
        K0s_cosThetaProdPlane_corr_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/efficiency/K0s/K0s_cosThetaProdPlane_eff.pdf(", "pdf");
      }
      else if( pTbin == nPtBins && etaBin == nEtaBins)
      {
        K0s_cosThetaProdPlane_corr_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/efficiency/K0s/K0s_cosThetaProdPlane_eff.pdf)", "pdf");

      }
      else
      {
        K0s_cosThetaProdPlane_corr_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/efficiency/K0s/K0s_cosThetaProdPlane_eff.pdf", "pdf");
      }

    }
  }

  TCanvas *K0s_K0s_cosThetaProdPlane_can = new TCanvas("K0s_K0s_cosThetaProdPlane_can", "K0s_K0s_cosThetaProdPlane_can", 1200, 1000);

  K0s_K0s_cosThetaProdPlane_cuts_hist->SetNameTitle("K0s_K0s_cosThetaProdPlane_eff", "K0s_K0s_cosThetaProdPlane_eff");
  K0s_K0s_cosThetaProdPlane_cuts_hist->SetMarkerStyle(20);
  K0s_K0s_cosThetaProdPlane_cuts_hist->SetMarkerColor(kRed);
  K0s_K0s_cosThetaProdPlane_cuts_hist->SetLineColor(kRed);
  K0s_K0s_cosThetaProdPlane_cuts_hist->GetXaxis()->SetTitle("cos(#theta*)");
  K0s_K0s_cosThetaProdPlane_cuts_hist->GetXaxis()->CenterTitle();
  K0s_K0s_cosThetaProdPlane_cuts_hist->GetYaxis()->SetTitle("Acceptance");
  K0s_K0s_cosThetaProdPlane_cuts_hist->GetYaxis()->CenterTitle();
  K0s_K0s_cosThetaProdPlane_cuts_hist->Sumw2();
  K0s_K0s_cosThetaProdPlane_cuts_hist->Divide(K0s_K0s_cosThetaProdPlane_cuts_hist, K0s_K0s_cosThetaProdPlane_hist, 1, 1, "b"); //binomial errors
  K0s_K0s_cosThetaProdPlane_cuts_hist->SetMinimum(0);
  K0s_K0s_cosThetaProdPlane_cuts_hist->Draw("p e");
  K0s_K0s_cosThetaProdPlane_cuts_hist->Write();

  TPaveText *K0s_K0s_text = new TPaveText(0.4, 0.2, 0.8, 0.4, "NDC");
  //K0s_K0s_text->AddText("STAR preliminary");
  //((TText*)K0s_K0s_text->GetListOfLines()->Last())->SetTextColor(2);
  K0s_K0s_text->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  K0s_K0s_text->AddText("Minimum bias");
  K0s_K0s_text->AddText("K_{s}^{0}-K_{s}^{0}");
  K0s_K0s_text->AddText("|#eta| < 1");
  K0s_K0s_text->AddText("#it{p}_{T} integrated");
  K0s_K0s_text->SetFillColorAlpha(0, 0.01);
  K0s_K0s_text->Draw("same");

  K0s_K0s_cosThetaProdPlane_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/efficiency/K0s/correlations/K0s_K0s_cosThetaProdPlane_eff.png");



  //K0s-K0s correlation histograms in bins
  for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
  {
    for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
    {
      TCanvas *K0s_K0s_cosThetaProdPlane_pT_corr_can = new TCanvas(Form("K0s_K0s_cosThetaProdPlane_pT_corr_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("K0s_K0s_cosThetaProdPlane_pT_corr_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);

      K0s_K0s_cosThetaProdPlane_pT_hist[pTbin1][pTbin2] = (TH1D*)inFile->Get(Form("K0s_K0s_cosThetaProdPlane_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
      K0s_K0s_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2] = (TH1D*)inFile->Get(Form("K0s_K0s_cosThetaProdPlane_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));

      K0s_K0s_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->Sumw2();
      K0s_K0s_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->SetNameTitle(Form("K0s_K0s_cosThetaProdPlane_pT1_%i_pT2_%i_eff", pTbin1, pTbin2), Form("K0s_K0s_cosThetaProdPlane_pT1_%i_pT2_%i_eff", pTbin1, pTbin2));
      K0s_K0s_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->Divide(K0s_K0s_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2], K0s_K0s_cosThetaProdPlane_pT_hist[pTbin1][pTbin2], 1, 1, "b");
      K0s_K0s_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      K0s_K0s_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->SetMarkerColor(kRed);
      K0s_K0s_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->SetLineColor(kRed);
      K0s_K0s_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->GetXaxis()->SetTitle("cos(#theta*)");
      K0s_K0s_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      K0s_K0s_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->GetYaxis()->SetTitle("Acceptance");
      K0s_K0s_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      K0s_K0s_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->SetMinimum(0);
      K0s_K0s_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->Draw("p e");
      K0s_K0s_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->Write();

      TPaveText *K0s_K0s_text_pT = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
      K0s_K0s_text_pT->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
      K0s_K0s_text_pT->AddText("Minimum bias");
      K0s_K0s_text_pT->AddText("K_{s}^{0}-K_{s}^{0}");
      //K0s_K0s_text_pt->AddText("No cuts_hist");
      K0s_K0s_text_pT->AddText("|#eta| < 1");
      K0s_K0s_text_pT->AddText(Form("%0.1f < p_{T}^{1} < %0.1f", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      K0s_K0s_text_pT->AddText(Form("%0.1f < p_{T}^{2} < %0.1f", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      K0s_K0s_text_pT->SetFillColorAlpha(0, 0.01);
      K0s_K0s_text_pT->Draw("same");

      K0s_K0s_cosThetaProdPlane_pT_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/efficiency/K0s/correlations/K0s_K0s_cosThetaProdPlane_pT1_%i_pT2_%i_eff.png", pTbin1, pTbin2));
      //_________________________________________________________________________________________________


      TCanvas *K0s_K0s_cosThetaProdPlane_tight_eta_pT_corr_can = new TCanvas(Form("K0s_K0s_cosThetaProdPlane_pT_corr_can_tight_eta_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("K0s_K0s_cosThetaProdPlane_pT_corr_can_tight_eta_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);

      K0s_K0s_cosThetaProdPlane_pT_hist_tight_eta[pTbin1][pTbin2] = (TH1D*)inFile->Get(Form("K0s_K0s_cosThetaProdPlane_tight_eta_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
      K0s_K0s_cosThetaProdPlane_pT_cuts_hist_tight_eta[pTbin1][pTbin2] = (TH1D*)inFile->Get(Form("K0s_K0s_cosThetaProdPlane_cuts_tight_eta_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));

      K0s_K0s_cosThetaProdPlane_pT_cuts_hist_tight_eta[pTbin1][pTbin2]->Sumw2();
      K0s_K0s_cosThetaProdPlane_pT_cuts_hist_tight_eta[pTbin1][pTbin2]->SetNameTitle(Form("K0s_K0s_cosThetaProdPlane_tight_eta_pT1_%i_pT2_%i_eff", pTbin1, pTbin2), Form("K0s_K0s_cosThetaProdPlane_tight_eta_pT1_%i_pT2_%i_eff", pTbin1, pTbin2));
      K0s_K0s_cosThetaProdPlane_pT_cuts_hist_tight_eta[pTbin1][pTbin2]->Divide(K0s_K0s_cosThetaProdPlane_pT_cuts_hist_tight_eta[pTbin1][pTbin2], K0s_K0s_cosThetaProdPlane_pT_hist_tight_eta[pTbin1][pTbin2], 1, 1, "b");
      K0s_K0s_cosThetaProdPlane_pT_cuts_hist_tight_eta[pTbin1][pTbin2]->SetMarkerStyle(20);
      K0s_K0s_cosThetaProdPlane_pT_cuts_hist_tight_eta[pTbin1][pTbin2]->SetMarkerColor(kRed);
      K0s_K0s_cosThetaProdPlane_pT_cuts_hist_tight_eta[pTbin1][pTbin2]->SetLineColor(kRed);
      K0s_K0s_cosThetaProdPlane_pT_cuts_hist_tight_eta[pTbin1][pTbin2]->GetXaxis()->SetTitle("cos(#theta*)");
      K0s_K0s_cosThetaProdPlane_pT_cuts_hist_tight_eta[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      K0s_K0s_cosThetaProdPlane_pT_cuts_hist_tight_eta[pTbin1][pTbin2]->GetYaxis()->SetTitle("Acceptance");
      K0s_K0s_cosThetaProdPlane_pT_cuts_hist_tight_eta[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      K0s_K0s_cosThetaProdPlane_pT_cuts_hist_tight_eta[pTbin1][pTbin2]->SetMinimum(0);
      K0s_K0s_cosThetaProdPlane_pT_cuts_hist_tight_eta[pTbin1][pTbin2]->Draw("p e");
      K0s_K0s_cosThetaProdPlane_pT_cuts_hist_tight_eta[pTbin1][pTbin2]->Write();

      TPaveText *K0s_K0s_text_pT_tight_eta = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
      K0s_K0s_text_pT_tight_eta->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
      K0s_K0s_text_pT_tight_eta->AddText("Minimum bias");
      K0s_K0s_text_pT_tight_eta->AddText("K_{s}^{0}-K_{s}^{0}");
      //K0s_K0s_text_pT_tight_eta->AddText("No cuts_hist_tight_eta");
      K0s_K0s_text_pT_tight_eta->AddText("|#eta| < 0.2");
      K0s_K0s_text_pT_tight_eta->AddText(Form("%0.1f < p_{T}^{1} < %0.1f", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      K0s_K0s_text_pT_tight_eta->AddText(Form("%0.1f < p_{T}^{2} < %0.1f", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      K0s_K0s_text_pT_tight_eta->SetFillColorAlpha(0, 0.01);
      K0s_K0s_text_pT_tight_eta->Draw("same");

      K0s_K0s_cosThetaProdPlane_tight_eta_pT_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/efficiency/K0s/correlations/K0s_K0s_cosThetaProdPlane_tight_eta_pT1_%i_pT2_%i_eff.png", pTbin1, pTbin2));
      //_________________________________________________________________________________________________



    }
  }

  for(unsigned int etaBin1 = 0; etaBin1 < nEtaBins; etaBin1++)
  {
    for(unsigned int etaBin2 = 0; etaBin2 < nEtaBins; etaBin2++)
    {

      TCanvas *K0s_K0s_cosThetaProdPlane_eta_corr_can = new TCanvas(Form("K0s_K0s_cosThetaProdPlane_eta_corr_can_eta1_%i_eta2_%i", etaBin1, etaBin2), Form("K0s_K0s_cosThetaProdPlane_eta_corr_can_eta1_%i_eta2_%i", etaBin1, etaBin2), 1200, 1000);

      K0s_K0s_cosThetaProdPlane_eta_hist[etaBin1][etaBin2] = (TH1D*)inFile->Get(Form("K0s_K0s_cosThetaProdPlane_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));
      K0s_K0s_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2] = (TH1D*)inFile->Get(Form("K0s_K0s_cosThetaProdPlane_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));

      K0s_K0s_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->Sumw2();
      K0s_K0s_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->SetNameTitle(Form("K0s_K0s_cosThetaProdPlane_eta1_%i_eta2_%i_eff", etaBin1, etaBin2), Form("K0s_K0s_cosThetaProdPlane_eta1_%i_eta2_%i_eff", etaBin1, etaBin2));
      K0s_K0s_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->Divide(K0s_K0s_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2], K0s_K0s_cosThetaProdPlane_eta_hist[etaBin1][etaBin2], 1, 1, "b");
      K0s_K0s_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      K0s_K0s_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->SetMarkerColor(kRed);
      K0s_K0s_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->SetLineColor(kRed);
      K0s_K0s_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->GetXaxis()->SetTitle("cos(#theta*)");
      K0s_K0s_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->GetXaxis()->CenterTitle();
      K0s_K0s_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->GetYaxis()->SetTitle("Acceptance");
      K0s_K0s_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->GetYaxis()->CenterTitle();
      K0s_K0s_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->SetMinimum(0);
      K0s_K0s_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->Draw("p e");
      K0s_K0s_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->Write();

      TPaveText *K0s_K0s_text_eta = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
      K0s_K0s_text_eta->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
      K0s_K0s_text_eta->AddText("Minimum bias");
      K0s_K0s_text_eta->AddText("K_{s}^{0}-K_{s}^{0}");
      //K0s_K0s_text_pt->AddText("No cuts_hist");
      K0s_K0s_text_eta->AddText(Form("%0.1f < #eta_{1} < %0.1f", eta_bins[etaBin1], eta_bins[etaBin1+1]));
      K0s_K0s_text_eta->AddText(Form("%0.1f < #eta_{2} < %0.1f", eta_bins[etaBin2], eta_bins[etaBin2+1]));
      K0s_K0s_text_eta->AddText("p_{T} integrated");
      K0s_K0s_text_eta->SetFillColorAlpha(0, 0.01);
      K0s_K0s_text_eta->Draw("same");

      K0s_K0s_cosThetaProdPlane_eta_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/efficiency/K0s/correlations/K0s_K0s_cosThetaProdPlane_eta1_%i_eta2_%i_eff.png", etaBin1, etaBin2));
      //_________________________________________________________________________________________________


    }
  }


  inFile->Close();
  outFile->Close();

  return;

}
