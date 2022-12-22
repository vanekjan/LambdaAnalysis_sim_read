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
void Acceptance_corr(const int energy = 510)
{

  const int nPtBins = 8;
  float const pT_bins[nPtBins+1] = { 0., 0.5, 1.,1.5, 2., 2.5, 3., 4., 5.};

  const int nPtBins_corr = 2;
  float const pT_bins_corr[nPtBins_corr+1] = { 0.5, 1.5, 5.};

  const int nEtaBins = 3;
  float const eta_bins[nEtaBins+1] = { -1, -0.4, 0.4, 1 };

  TFile *outFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/output/L_cosThetaStar_eff_work.root", "recreate");

  //load all files
  TFile *inFile;


  if(energy == 510)
  {
    inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run17/output_Lambda_pp_510_MB_1B_events_new.root", "READ");
  }
  else if(energy == 200)
  {
    inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/output_Lambda_pp_200_MB_1B_events_new.root", "READ");
  }
  else
  {
    cout<<"Not a valid collision energy! Abborting!"<<endl;
  }



  //histograms

  TH1D *L0_cosThetaProdPlane_open_hist[nPtBins+1][nEtaBins+1];
  TH1D *L0_cosThetaProdPlane_ana_hist[nPtBins+1][nEtaBins+1];

  TH1D *L0bar_cosThetaProdPlane_open_hist[nPtBins+1][nEtaBins+1];
  TH1D *L0bar_cosThetaProdPlane_ana_hist[nPtBins+1][nEtaBins+1];


  TH1D *L0_L0bar_cosThetaProdPlane = (TH1D*)inFile->Get("L0_L0bar_cosThetaProdPlane");
  TH1D *L0_L0_cosThetaProdPlane = (TH1D*)inFile->Get("L0_L0_cosThetaProdPlane");
  TH1D *L0bar_L0bar_cosThetaProdPlane = (TH1D*)inFile->Get("L0bar_L0bar_cosThetaProdPlane");

  TH1D *L0_L0bar_cosThetaProdPlane_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1D *L0_L0bar_cosThetaProdPlane_eta_hist[nEtaBins][nEtaBins];

  TH1D *L0_L0_cosThetaProdPlane_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1D *L0_L0_cosThetaProdPlane_eta_hist[nEtaBins][nEtaBins];

  TH1D *L0bar_L0bar_cosThetaProdPlane_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1D *L0bar_L0bar_cosThetaProdPlane_eta_hist[nEtaBins][nEtaBins];


  TH1D *L0_L0bar_cosThetaProdPlane_cuts = (TH1D*)inFile->Get("L0_L0bar_cosThetaProdPlane_cuts");
  TH1D *L0_L0_cosThetaProdPlane_cuts = (TH1D*)inFile->Get("L0_L0_cosThetaProdPlane_cuts");
  TH1D *L0bar_L0bar_cosThetaProdPlane_cuts = (TH1D*)inFile->Get("L0bar_L0bar_cosThetaProdPlane_cuts");

  TH1D *L0_L0bar_cosThetaProdPlane_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1D *L0_L0bar_cosThetaProdPlane_eta_cuts_hist[nEtaBins][nEtaBins];

  TH1D *L0_L0_cosThetaProdPlane_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1D *L0_L0_cosThetaProdPlane_eta_cuts_hist[nEtaBins][nEtaBins];

  TH1D *L0bar_L0bar_cosThetaProdPlane_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1D *L0bar_L0bar_cosThetaProdPlane_eta_cuts_hist[nEtaBins][nEtaBins];


  //mixed event histograms
  TH1D *L0_L0bar_cosThetaProdPlane_ME = (TH1D*)inFile->Get("L0_L0bar_cosThetaProdPlane_ME");
  TH1D *L0_L0_cosThetaProdPlane_ME = (TH1D*)inFile->Get("L0_L0_cosThetaProdPlane_ME");
  TH1D *L0bar_L0bar_cosThetaProdPlane_ME = (TH1D*)inFile->Get("L0bar_L0bar_cosThetaProdPlane_ME");

  TH1D *L0_L0bar_cosThetaProdPlane_ME_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1D *L0_L0bar_cosThetaProdPlane_ME_eta_hist[nEtaBins][nEtaBins];

  TH1D *L0_L0_cosThetaProdPlane_ME_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1D *L0_L0_cosThetaProdPlane_ME_eta_hist[nEtaBins][nEtaBins];

  TH1D *L0bar_L0bar_cosThetaProdPlane_ME_pT_hist[nPtBins_corr][nPtBins_corr];
  TH1D *L0bar_L0bar_cosThetaProdPlane_ME_eta_hist[nEtaBins][nEtaBins];


  TH1D *L0_L0bar_cosThetaProdPlane_ME_cuts = (TH1D*)inFile->Get("L0_L0bar_cosThetaProdPlane_ME_cuts");
  TH1D *L0_L0_cosThetaProdPlane_ME_cuts = (TH1D*)inFile->Get("L0_L0_cosThetaProdPlane_ME_cuts");
  TH1D *L0bar_L0bar_cosThetaProdPlane_ME_cuts = (TH1D*)inFile->Get("L0bar_L0bar_cosThetaProdPlane_ME_cuts");

  TH1D *L0_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1D *L0_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[nEtaBins][nEtaBins];

  TH1D *L0_L0_cosThetaProdPlane_ME_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1D *L0_L0_cosThetaProdPlane_ME_eta_cuts_hist[nEtaBins][nEtaBins];

  TH1D *L0bar_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[nPtBins_corr][nPtBins_corr];
  TH1D *L0bar_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[nEtaBins][nEtaBins];


  //_____________________________________________________________________________

  TCanvas *L0_cosThetaProdPlane_corr_can[nPtBins+1][nEtaBins+1];

  TCanvas *L0bar_cosThetaProdPlane_corr_can[nPtBins+1][nEtaBins+1];

  TCanvas *L0_L0bar_cosThetaProdPlane_corr_can = new TCanvas("L0_L0bar_cosThetaProdPlane_corr_can", "L0_L0bar_cosThetaProdPlane_corr_can", 1200, 1000);
  TCanvas *L0_L0_cosThetaProdPlane_corr_can = new TCanvas("L0_L0_cosThetaProdPlane_corr_can", "L0_L0_cosThetaProdPlane_corr_can", 1200, 1000);
  TCanvas *L0bar_L0bar_cosThetaProdPlane_corr_can = new TCanvas("L0bar_L0bar_cosThetaProdPlane_corr_can", "L0bar_L0bar_cosThetaProdPlane_corr_can", 1200, 1000);

  TCanvas *L0_L0bar_cosThetaProdPlane_ME_corr_can = new TCanvas("L0_L0bar_cosThetaProdPlane_ME_corr_can", "L0_L0bar_cosThetaProdPlane_ME_corr_can", 1200, 1000);
  TCanvas *L0_L0_cosThetaProdPlane_ME_corr_can = new TCanvas("L0_L0_cosThetaProdPlane_ME_corr_can", "L0_L0_cosThetaProdPlane_ME_corr_can", 1200, 1000);
  TCanvas *L0bar_L0bar_cosThetaProdPlane_ME_corr_can = new TCanvas("L0bar_L0bar_cosThetaProdPlane_ME_corr_can", "L0bar_L0bar_cosThetaProdPlane_ME_corr_can", 1200, 1000);


  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);



  outFile->cd();

  //single L efficiencies
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


      L0_cosThetaProdPlane_corr_can[pTbin][etaBin] = new TCanvas(Form("L0_cosThetaProdPlane_corr_can_pT_%i_eta_%i", pTbin, etaBin), Form("L0_cosThetaProdPlane_corr_can_pT_%i_eta_%i", pTbin, etaBin), 1200, 1000);

      L0bar_cosThetaProdPlane_corr_can[pTbin][etaBin] = new TCanvas(Form("L0bar_cosThetaProdPlane_corr_can_pT_%i_eta_%i", pTbin, etaBin), Form("L0bar_cosThetaProdPlane_corr_can_pT_%i_eta_%i", pTbin, etaBin), 1200, 1000);



      TPaveText *cent_text_2 = new TPaveText(0.4, 0.2, 0.8, 0.4, "NDC");
      //cent_text->AddText("STAR preliminary");
      //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
      cent_text_2->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
      cent_text_2->AddText("Minimum bias");
      cent_text_2->AddText("#Lambda^{0}");
      cent_text_2->AddText(eta_range->Data());
      cent_text_2->AddText(pT_range->Data());
      cent_text_2->SetFillColorAlpha(0, 0.01);
      //cent_text_2->Draw("same");


      L0_cosThetaProdPlane_corr_can[pTbin][etaBin]->cd();

      L0_cosThetaProdPlane_open_hist[pTbin][etaBin] = (TH1D*)inFile->Get(Form("L0_cosThetaProdPlane_pT_%i_eta_%i", pTbin, etaBin));
      L0_cosThetaProdPlane_open_hist[pTbin][etaBin]->SetNameTitle(Form("L0_cosThetaProdPlane_open_pT_%i_eta_%i", pTbin, etaBin), Form("L0_cosThetaProdPlane_open_pT_%i_eta_%i", pTbin, etaBin));

      L0_cosThetaProdPlane_ana_hist[pTbin][etaBin] = (TH1D*)inFile->Get(Form("L0_cosThetaProdPlane_cuts_pT_%i_eta_%i", pTbin, etaBin));
      L0_cosThetaProdPlane_ana_hist[pTbin][etaBin]->SetNameTitle(Form("L0_cosThetaProdPlane_eff_pT_%i_eta_%i", pTbin, etaBin), Form("L0_cosThetaProdPlane_eff_pT_%i_eta_%i", pTbin, etaBin));
      L0_cosThetaProdPlane_ana_hist[pTbin][etaBin]->SetMarkerStyle(20);
      L0_cosThetaProdPlane_ana_hist[pTbin][etaBin]->SetMarkerColor(kRed);
      L0_cosThetaProdPlane_ana_hist[pTbin][etaBin]->SetLineColor(kRed);
      L0_cosThetaProdPlane_ana_hist[pTbin][etaBin]->GetXaxis()->SetTitle("cos(#theta*)");
      L0_cosThetaProdPlane_ana_hist[pTbin][etaBin]->GetXaxis()->CenterTitle();
      L0_cosThetaProdPlane_ana_hist[pTbin][etaBin]->GetYaxis()->SetTitle("Efficiency");
      L0_cosThetaProdPlane_ana_hist[pTbin][etaBin]->GetYaxis()->CenterTitle();
      L0_cosThetaProdPlane_ana_hist[pTbin][etaBin]->Sumw2();
      L0_cosThetaProdPlane_ana_hist[pTbin][etaBin]->Divide(L0_cosThetaProdPlane_ana_hist[pTbin][etaBin] ,L0_cosThetaProdPlane_open_hist[pTbin][etaBin], 1, 1, "b"); //binomial errors for efficiency
      L0_cosThetaProdPlane_ana_hist[pTbin][etaBin]->SetMinimum(0);
      L0_cosThetaProdPlane_ana_hist[pTbin][etaBin]->Draw("p e");
      L0_cosThetaProdPlane_ana_hist[pTbin][etaBin]->Write();


      TLegend *PYTHIA_legend = new TLegend(0.2, 0.4, 0.4, 0.5 );
      PYTHIA_legend->AddEntry(L0_cosThetaProdPlane_ana_hist[pTbin][etaBin], "PYTHIA 8.307");
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




      L0bar_cosThetaProdPlane_corr_can[pTbin][etaBin]->cd();

      L0bar_cosThetaProdPlane_open_hist[pTbin][etaBin] = (TH1D*)inFile->Get(Form("L0bar_cosThetaProdPlane_pT_%i_eta_%i", pTbin, etaBin));
      L0bar_cosThetaProdPlane_open_hist[pTbin][etaBin]->SetNameTitle(Form("L0bar_cosThetaProdPlane_open_pT_%i_eta_%i", pTbin, etaBin), Form("L0bar_cosThetaProdPlane_open_pT_%i_eta_%i", pTbin, etaBin));

      L0bar_cosThetaProdPlane_ana_hist[pTbin][etaBin] = (TH1D*)inFile->Get(Form("L0bar_cosThetaProdPlane_cuts_pT_%i_eta_%i", pTbin, etaBin));
      L0bar_cosThetaProdPlane_ana_hist[pTbin][etaBin]->SetNameTitle(Form("L0bar_cosThetaProdPlane_eff_pT_%i_eta_%i", pTbin, etaBin), Form("L0bar_cosThetaProdPlane_eff_pT_%i_eta_%i", pTbin, etaBin));
      L0bar_cosThetaProdPlane_ana_hist[pTbin][etaBin]->SetMarkerStyle(20);
      L0bar_cosThetaProdPlane_ana_hist[pTbin][etaBin]->SetMarkerColor(kRed);
      L0bar_cosThetaProdPlane_ana_hist[pTbin][etaBin]->SetLineColor(kRed);
      L0bar_cosThetaProdPlane_ana_hist[pTbin][etaBin]->GetXaxis()->SetTitle("cos(#theta*)");
      L0bar_cosThetaProdPlane_ana_hist[pTbin][etaBin]->GetXaxis()->CenterTitle();
      L0bar_cosThetaProdPlane_ana_hist[pTbin][etaBin]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      L0bar_cosThetaProdPlane_ana_hist[pTbin][etaBin]->GetYaxis()->CenterTitle();
      L0bar_cosThetaProdPlane_ana_hist[pTbin][etaBin]->Sumw2();
      L0bar_cosThetaProdPlane_ana_hist[pTbin][etaBin]->Divide(L0bar_cosThetaProdPlane_ana_hist[pTbin][etaBin], L0bar_cosThetaProdPlane_open_hist[pTbin][etaBin], 1, 1, "b"); //scale by vin width to obtain d N/d cos(theta*)
      L0bar_cosThetaProdPlane_ana_hist[pTbin][etaBin]->SetMinimum(0);
      L0bar_cosThetaProdPlane_ana_hist[pTbin][etaBin]->Draw("p e");
      L0bar_cosThetaProdPlane_ana_hist[pTbin][etaBin]->Write();


      PYTHIA_legend->Draw("same");

      cent_text_3->Draw("same");


      L0_cosThetaProdPlane_corr_can[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/efficiency/Lambda/L0_cosThetaProdPlane_eff_pT_%i_eta_%i.png", pTbin, etaBin));

      L0bar_cosThetaProdPlane_corr_can[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/efficiency/Lambda/L0bar_cosThetaProdPlane_eff_pT_%i_eta_%i.png", pTbin, etaBin));



      //save histograms to multi-page PDF
      if(pTbin == 0 && etaBin == 0)
      {
        L0_cosThetaProdPlane_corr_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/efficiency/Lambda/L0_cosThetaProdPlane_eff.pdf(", "pdf");

        L0bar_cosThetaProdPlane_corr_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/efficiency/Lambda/L0bar_cosThetaProdPlane_eff.pdf(", "pdf");


      }
      else if( pTbin == nPtBins && etaBin == nEtaBins)
      {
        L0_cosThetaProdPlane_corr_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/efficiency/Lambda/L0_cosThetaProdPlane_eff.pdf)", "pdf");

        L0bar_cosThetaProdPlane_corr_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/efficiency/Lambda/L0bar_cosThetaProdPlane_eff.pdf)", "pdf");

      }
      else
      {
        L0_cosThetaProdPlane_corr_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/efficiency/Lambda/L0_cosThetaProdPlane_eff.pdf", "pdf");

        L0bar_cosThetaProdPlane_corr_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/efficiency/Lambda/L0bar_cosThetaProdPlane_eff.pdf", "pdf");
      }

    }
  }
  //__________________________________________________________________________________________________________________________________


  //L-L correlations efficiencies
  L0_L0bar_cosThetaProdPlane_corr_can->cd();

  L0_L0bar_cosThetaProdPlane_cuts->Sumw2();
  L0_L0bar_cosThetaProdPlane_cuts->SetNameTitle("L0_L0bar_cosThetaProdPlane_eff", "L0_L0bar_cosThetaProdPlane_eff");
  L0_L0bar_cosThetaProdPlane_cuts->Divide(L0_L0bar_cosThetaProdPlane_cuts, L0_L0bar_cosThetaProdPlane, 1, 1, "b");
  L0_L0bar_cosThetaProdPlane_cuts->SetMarkerStyle(20);
  L0_L0bar_cosThetaProdPlane_cuts->SetMarkerColor(kRed);
  L0_L0bar_cosThetaProdPlane_cuts->SetLineColor(kRed);
  L0_L0bar_cosThetaProdPlane_cuts->GetXaxis()->SetTitle("cos(#theta*)");
  L0_L0bar_cosThetaProdPlane_cuts->GetXaxis()->CenterTitle();
  L0_L0bar_cosThetaProdPlane_cuts->GetYaxis()->SetTitle("Acceptance");
  L0_L0bar_cosThetaProdPlane_cuts->GetYaxis()->CenterTitle();
  L0_L0bar_cosThetaProdPlane_cuts->SetMinimum(0);
  L0_L0bar_cosThetaProdPlane_cuts->Draw("p e");
  L0_L0bar_cosThetaProdPlane_cuts->Write();

  TPaveText *L0_L0bar_text_pt = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
  L0_L0bar_text_pt->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  L0_L0bar_text_pt->AddText("Minimum bias");
  L0_L0bar_text_pt->AddText("#Lambda-#bar{#Lambda}");
  //L0_L0bar_text_pt->AddText("No cuts");
  //L0bar_text_pt->AddText(pT_range->Data());
  L0_L0bar_text_pt->SetFillColorAlpha(0, 0.01);
  L0_L0bar_text_pt->Draw("same");

  L0_L0bar_cosThetaProdPlane_corr_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/efficiency/Lambda/correlations/L0_L0bar_cosThetaProdPlane_corr_eff.png");



  L0_L0_cosThetaProdPlane_corr_can->cd();

  L0_L0_cosThetaProdPlane_cuts->Sumw2();
  L0_L0_cosThetaProdPlane_cuts->SetNameTitle("L0_L0_cosThetaProdPlane_eff", "L0_L0_cosThetaProdPlane_eff");
  L0_L0_cosThetaProdPlane_cuts->Divide(L0_L0_cosThetaProdPlane_cuts, L0_L0_cosThetaProdPlane, 1, 1, "b");
  L0_L0_cosThetaProdPlane_cuts->SetMarkerStyle(20);
  L0_L0_cosThetaProdPlane_cuts->SetMarkerColor(kRed);
  L0_L0_cosThetaProdPlane_cuts->SetLineColor(kRed);
  L0_L0_cosThetaProdPlane_cuts->GetXaxis()->SetTitle("cos(#theta*)");
  L0_L0_cosThetaProdPlane_cuts->GetXaxis()->CenterTitle();
  L0_L0_cosThetaProdPlane_cuts->GetYaxis()->SetTitle("Acceptance");
  L0_L0_cosThetaProdPlane_cuts->GetYaxis()->CenterTitle();
  L0_L0_cosThetaProdPlane_cuts->SetMinimum(0);
  L0_L0_cosThetaProdPlane_cuts->Draw("p e");
  L0_L0_cosThetaProdPlane_cuts->Write();

  TPaveText *L0_L0_text_pt = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
  L0_L0_text_pt->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  L0_L0_text_pt->AddText("Minimum bias");
  L0_L0_text_pt->AddText("#Lambda-#Lambda");
  //L0_L0_text_pt->AddText("No cuts");
  //L0bar_text_pt->AddText(pT_range->Data());
  L0_L0_text_pt->SetFillColorAlpha(0, 0.01);
  L0_L0_text_pt->Draw("same");

  L0_L0_cosThetaProdPlane_corr_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/efficiency/Lambda/correlations/L0_L0_cosThetaProdPlane_corr_eff.png");



  L0bar_L0bar_cosThetaProdPlane_corr_can->cd();

  L0bar_L0bar_cosThetaProdPlane_cuts->Sumw2();
  L0bar_L0bar_cosThetaProdPlane_cuts->SetNameTitle("L0bar_L0bar_cosThetaProdPlane_eff", "L0bar_L0bar_cosThetaProdPlane_eff");
  L0bar_L0bar_cosThetaProdPlane_cuts->Divide(L0bar_L0bar_cosThetaProdPlane_cuts, L0bar_L0bar_cosThetaProdPlane, 1, 1, "b");
  L0bar_L0bar_cosThetaProdPlane_cuts->SetMarkerStyle(20);
  L0bar_L0bar_cosThetaProdPlane_cuts->SetMarkerColor(kRed);
  L0bar_L0bar_cosThetaProdPlane_cuts->SetLineColor(kRed);
  L0bar_L0bar_cosThetaProdPlane_cuts->GetXaxis()->SetTitle("cos(#theta*)");
  L0bar_L0bar_cosThetaProdPlane_cuts->GetXaxis()->CenterTitle();
  L0bar_L0bar_cosThetaProdPlane_cuts->GetYaxis()->SetTitle("Acceptance");
  L0bar_L0bar_cosThetaProdPlane_cuts->GetYaxis()->CenterTitle();
  L0bar_L0bar_cosThetaProdPlane_cuts->SetMinimum(0);
  L0bar_L0bar_cosThetaProdPlane_cuts->Draw("p e");
  L0bar_L0bar_cosThetaProdPlane_cuts->Write();

  TPaveText *L0bar_L0bar_text_pt = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
  L0bar_L0bar_text_pt->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  L0bar_L0bar_text_pt->AddText("Minimum bias");
  L0bar_L0bar_text_pt->AddText("#bar{#Lambda}-#bar{#Lambda}");
  //L0bar_L0bar_text_pt->AddText("No cuts");
  //L0bar_text_pt->AddText(pT_range->Data());
  L0bar_L0bar_text_pt->SetFillColorAlpha(0, 0.01);
  L0bar_L0bar_text_pt->Draw("same");

  L0bar_L0bar_cosThetaProdPlane_corr_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/efficiency/Lambda/correlations/L0bar_L0bar_cosThetaProdPlane_corr_eff.png");


  //L-L correlations efficiencies - mixed event
  L0_L0bar_cosThetaProdPlane_ME_corr_can->cd();

  L0_L0bar_cosThetaProdPlane_ME_cuts->Sumw2();
  L0_L0bar_cosThetaProdPlane_ME_cuts->SetNameTitle("L0_L0bar_cosThetaProdPlane_ME_eff", "L0_L0bar_cosThetaProdPlane_ME_eff");
  L0_L0bar_cosThetaProdPlane_ME_cuts->Divide(L0_L0bar_cosThetaProdPlane_ME_cuts, L0_L0bar_cosThetaProdPlane_ME, 1, 1, "b");
  L0_L0bar_cosThetaProdPlane_ME_cuts->SetMarkerStyle(20);
  L0_L0bar_cosThetaProdPlane_ME_cuts->SetMarkerColor(kRed);
  L0_L0bar_cosThetaProdPlane_ME_cuts->SetLineColor(kRed);
  L0_L0bar_cosThetaProdPlane_ME_cuts->GetXaxis()->SetTitle("cos(#theta*)");
  L0_L0bar_cosThetaProdPlane_ME_cuts->GetXaxis()->CenterTitle();
  L0_L0bar_cosThetaProdPlane_ME_cuts->GetYaxis()->SetTitle("Acceptance");
  L0_L0bar_cosThetaProdPlane_ME_cuts->GetYaxis()->CenterTitle();
  L0_L0bar_cosThetaProdPlane_ME_cuts->SetMinimum(0);
  L0_L0bar_cosThetaProdPlane_ME_cuts->Draw("p e");
  L0_L0bar_cosThetaProdPlane_ME_cuts->Write();

  L0_L0bar_text_pt->Draw("same");

  L0_L0bar_cosThetaProdPlane_ME_corr_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/efficiency/Lambda/correlations/L0_L0bar_cosThetaProdPlane_ME_corr_eff.png");



  L0_L0_cosThetaProdPlane_ME_corr_can->cd();

  L0_L0_cosThetaProdPlane_ME_cuts->Sumw2();
  L0_L0_cosThetaProdPlane_ME_cuts->SetNameTitle("L0_L0_cosThetaProdPlane_ME_eff", "L0_L0_cosThetaProdPlane_ME_eff");
  L0_L0_cosThetaProdPlane_ME_cuts->Divide(L0_L0_cosThetaProdPlane_ME_cuts, L0_L0_cosThetaProdPlane_ME, 1, 1, "b");
  L0_L0_cosThetaProdPlane_ME_cuts->SetMarkerStyle(20);
  L0_L0_cosThetaProdPlane_ME_cuts->SetMarkerColor(kRed);
  L0_L0_cosThetaProdPlane_ME_cuts->SetLineColor(kRed);
  L0_L0_cosThetaProdPlane_ME_cuts->GetXaxis()->SetTitle("cos(#theta*)");
  L0_L0_cosThetaProdPlane_ME_cuts->GetXaxis()->CenterTitle();
  L0_L0_cosThetaProdPlane_ME_cuts->GetYaxis()->SetTitle("Acceptance");
  L0_L0_cosThetaProdPlane_ME_cuts->GetYaxis()->CenterTitle();
  L0_L0_cosThetaProdPlane_ME_cuts->SetMinimum(0);
  L0_L0_cosThetaProdPlane_ME_cuts->Draw("p e");
  L0_L0_cosThetaProdPlane_ME_cuts->Write();

  L0_L0_text_pt->Draw("same");

  L0_L0_cosThetaProdPlane_ME_corr_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/efficiency/Lambda/correlations/L0_L0_cosThetaProdPlane_ME_corr_eff.png");



  L0bar_L0bar_cosThetaProdPlane_ME_corr_can->cd();

  L0bar_L0bar_cosThetaProdPlane_ME_cuts->Sumw2();
  L0bar_L0bar_cosThetaProdPlane_ME_cuts->SetNameTitle("L0bar_L0bar_cosThetaProdPlane_ME_eff", "L0bar_L0bar_cosThetaProdPlane_ME_eff");
  L0bar_L0bar_cosThetaProdPlane_ME_cuts->Divide(L0bar_L0bar_cosThetaProdPlane_ME_cuts, L0bar_L0bar_cosThetaProdPlane_ME, 1, 1, "b");
  L0bar_L0bar_cosThetaProdPlane_ME_cuts->SetMarkerStyle(20);
  L0bar_L0bar_cosThetaProdPlane_ME_cuts->SetMarkerColor(kRed);
  L0bar_L0bar_cosThetaProdPlane_ME_cuts->SetLineColor(kRed);
  L0bar_L0bar_cosThetaProdPlane_ME_cuts->GetXaxis()->SetTitle("cos(#theta*)");
  L0bar_L0bar_cosThetaProdPlane_ME_cuts->GetXaxis()->CenterTitle();
  L0bar_L0bar_cosThetaProdPlane_ME_cuts->GetYaxis()->SetTitle("Acceptance");
  L0bar_L0bar_cosThetaProdPlane_ME_cuts->GetYaxis()->CenterTitle();
  L0bar_L0bar_cosThetaProdPlane_ME_cuts->SetMinimum(0);
  L0bar_L0bar_cosThetaProdPlane_ME_cuts->Draw("p e");
  L0bar_L0bar_cosThetaProdPlane_ME_cuts->Write();

  L0bar_L0bar_text_pt->Draw("same");

  L0bar_L0bar_cosThetaProdPlane_ME_corr_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/efficiency/Lambda/correlations/L0bar_L0bar_cosThetaProdPlane_ME_corr_eff.png");


  //L-L correlation histograms in bins
  for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
  {
    for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
    {

      //signal efficiencies
      TCanvas *L0_L0bar_cosThetaProdPlane_pT_corr_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_pT_corr_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_pT_corr_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);

      L0_L0bar_cosThetaProdPlane_pT_hist[pTbin1][pTbin2] = (TH1D*)inFile->Get(Form("L0_L0bar_cosThetaProdPlane_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
      L0_L0bar_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2] = (TH1D*)inFile->Get(Form("L0_L0bar_cosThetaProdPlane_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));

      L0_L0bar_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->Sumw2();
      L0_L0bar_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->SetNameTitle(Form("L0_L0bar_cosThetaProdPlane_pT1_%i_pT2_%i_eff", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_pT1_%i_pT2_%i_eff", pTbin1, pTbin2));
      L0_L0bar_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->Divide(L0_L0bar_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2], L0_L0bar_cosThetaProdPlane_pT_hist[pTbin1][pTbin2], 1, 1, "b");
      L0_L0bar_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0_L0bar_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0_L0bar_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->SetLineColor(kRed);
      L0_L0bar_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0_L0bar_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_L0bar_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->GetYaxis()->SetTitle("Acceptance");
      L0_L0bar_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_L0bar_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->SetMinimum(0);
      L0_L0bar_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->Draw("p e");
      L0_L0bar_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->Write();

      TPaveText *L0_L0bar_text_pT = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
      L0_L0bar_text_pT->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
      L0_L0bar_text_pT->AddText("Minimum bias");
      L0_L0bar_text_pT->AddText("#Lambda-#bar{#Lambda}");
      //L0_L0bar_text_pt->AddText("No cuts_hist");
      L0_L0bar_text_pT->AddText("|#eta| < 1");
      L0_L0bar_text_pT->AddText(Form("%0.1f < p_{T}^{1} < %0.1f", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0_L0bar_text_pT->AddText(Form("%0.1f < p_{T}^{2} < %0.1f", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0_L0bar_text_pT->SetFillColorAlpha(0, 0.01);
      L0_L0bar_text_pT->Draw("same");

      L0_L0bar_cosThetaProdPlane_pT_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/efficiency/Lambda/correlations/L0_L0bar_cosThetaProdPlane_pT1_%i_pT2_%i_eff.png", pTbin1, pTbin2));
      //_________________________________________________________________________________________________


      TCanvas *L0_L0_cosThetaProdPlane_pT_corr_can = new TCanvas(Form("L0_L0_cosThetaProdPlane_pT_corr_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_pT_corr_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);

      L0_L0_cosThetaProdPlane_pT_hist[pTbin1][pTbin2] = (TH1D*)inFile->Get(Form("L0_L0_cosThetaProdPlane_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
      L0_L0_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2] = (TH1D*)inFile->Get(Form("L0_L0_cosThetaProdPlane_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));

      L0_L0_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->Sumw2();
      L0_L0_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->SetNameTitle(Form("L0_L0_cosThetaProdPlane_pT1_%i_pT2_%i_eff", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_pT1_%i_pT2_%i_eff", pTbin1, pTbin2));
      L0_L0_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->Divide(L0_L0_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2], L0_L0_cosThetaProdPlane_pT_hist[pTbin1][pTbin2], 1, 1, "b");
      L0_L0_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0_L0_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0_L0_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->SetLineColor(kRed);
      L0_L0_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0_L0_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_L0_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->GetYaxis()->SetTitle("Acceptance");
      L0_L0_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_L0_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->SetMinimum(0);
      L0_L0_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->Draw("p e");
      L0_L0_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->Write();

      TPaveText *L0_L0_text_pT = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
      L0_L0_text_pT->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
      L0_L0_text_pT->AddText("Minimum bias");
      L0_L0_text_pT->AddText("#Lambda-#Lambda");
      //L0_L0_text_pt->AddText("No cuts_hist");
      L0_L0_text_pT->AddText("|#eta| < 1");
      L0_L0_text_pT->AddText(Form("%0.1f < p_{T}^{1} < %0.1f", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0_L0_text_pT->AddText(Form("%0.1f < p_{T}^{2} < %0.1f", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0_L0_text_pT->SetFillColorAlpha(0, 0.01);
      L0_L0_text_pT->Draw("same");

      L0_L0_cosThetaProdPlane_pT_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/efficiency/Lambda/correlations/L0_L0_cosThetaProdPlane_pT1_%i_pT2_%i_eff.png", pTbin1, pTbin2));
      //_________________________________________________________________________________________________

      TCanvas *L0bar_L0bar_cosThetaProdPlane_pT_corr_can = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_pT_corr_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_pT_corr_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);

      L0bar_L0bar_cosThetaProdPlane_pT_hist[pTbin1][pTbin2] = (TH1D*)inFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
      L0bar_L0bar_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2] = (TH1D*)inFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));

      L0bar_L0bar_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->Sumw2();
      L0bar_L0bar_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->SetNameTitle(Form("L0bar_L0bar_cosThetaProdPlane_pT1_%i_pT2_%i_eff", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_pT1_%i_pT2_%i_eff", pTbin1, pTbin2));
      L0bar_L0bar_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->Divide(L0bar_L0bar_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2], L0bar_L0bar_cosThetaProdPlane_pT_hist[pTbin1][pTbin2], 1, 1, "b");
      L0bar_L0bar_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0bar_L0bar_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0bar_L0bar_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->SetLineColor(kRed);
      L0bar_L0bar_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0bar_L0bar_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0bar_L0bar_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->GetYaxis()->SetTitle("Acceptance");
      L0bar_L0bar_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0bar_L0bar_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->SetMinimum(0);
      L0bar_L0bar_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->Draw("p e");
      L0bar_L0bar_cosThetaProdPlane_pT_cuts_hist[pTbin1][pTbin2]->Write();

      TPaveText *L0bar_L0bar_text_pT = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
      L0bar_L0bar_text_pT->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
      L0bar_L0bar_text_pT->AddText("Minimum bias");
      L0bar_L0bar_text_pT->AddText("#bar{#Lambda}-#bar{#Lambda}");
      L0bar_L0bar_text_pT->AddText("|#eta| < 1");
      //L0bar_L0bar_text_pt->AddText("No cuts_hist");
      L0bar_L0bar_text_pT->AddText(Form("%0.1f < p_{T}^{1} < %0.1f", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0bar_L0bar_text_pT->AddText(Form("%0.1f < p_{T}^{2} < %0.1f", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0bar_L0bar_text_pT->SetFillColorAlpha(0, 0.01);
      L0bar_L0bar_text_pT->Draw("same");

      L0bar_L0bar_cosThetaProdPlane_pT_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/efficiency/Lambda/correlations/L0bar_L0bar_cosThetaProdPlane_pT1_%i_pT2_%i_eff.png", pTbin1, pTbin2));
      //_________________________________________________________________________________________________

      //mixed event efficiencies
      TCanvas *L0_L0bar_cosThetaProdPlane_ME_pT_corr_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_ME_pT_corr_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_ME_pT_corr_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);

      L0_L0bar_cosThetaProdPlane_ME_pT_hist[pTbin1][pTbin2] = (TH1D*)inFile->Get(Form("L0_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
      L0_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2] = (TH1D*)inFile->Get(Form("L0_L0bar_cosThetaProdPlane_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));

      L0_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->Sumw2();
      L0_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->SetNameTitle(Form("L0_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_eff", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_eff", pTbin1, pTbin2));
      L0_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->Divide(L0_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2], L0_L0bar_cosThetaProdPlane_ME_pT_hist[pTbin1][pTbin2], 1, 1, "b");
      L0_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->SetLineColor(kRed);
      L0_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->GetYaxis()->SetTitle("Acceptance");
      L0_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->SetMinimum(0);
      L0_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->Draw("p e");
      L0_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->Write();

      L0_L0bar_text_pT->Draw("same");

      L0_L0bar_cosThetaProdPlane_ME_pT_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/efficiency/Lambda/correlations/L0_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_eff.png", pTbin1, pTbin2));
      //_________________________________________________________________________________________________


      TCanvas *L0_L0_cosThetaProdPlane_ME_pT_corr_can = new TCanvas(Form("L0_L0_cosThetaProdPlane_ME_pT_corr_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_ME_pT_corr_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);

      L0_L0_cosThetaProdPlane_ME_pT_hist[pTbin1][pTbin2] = (TH1D*)inFile->Get(Form("L0_L0_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
      L0_L0_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2] = (TH1D*)inFile->Get(Form("L0_L0_cosThetaProdPlane_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));

      L0_L0_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->Sumw2();
      L0_L0_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->SetNameTitle(Form("L0_L0_cosThetaProdPlane_ME_pT1_%i_pT2_%i_eff", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_ME_pT1_%i_pT2_%i_eff", pTbin1, pTbin2));
      L0_L0_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->Divide(L0_L0_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2], L0_L0_cosThetaProdPlane_ME_pT_hist[pTbin1][pTbin2], 1, 1, "b");
      L0_L0_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0_L0_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0_L0_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->SetLineColor(kRed);
      L0_L0_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0_L0_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_L0_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->GetYaxis()->SetTitle("Acceptance");
      L0_L0_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_L0_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->SetMinimum(0);
      L0_L0_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->Draw("p e");
      L0_L0_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->Write();

      L0_L0_text_pT->Draw("same");

      L0_L0_cosThetaProdPlane_ME_pT_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/efficiency/Lambda/correlations/L0_L0_cosThetaProdPlane_ME_pT1_%i_pT2_%i_eff.png", pTbin1, pTbin2));
      //_________________________________________________________________________________________________

      TCanvas *L0bar_L0bar_cosThetaProdPlane_ME_pT_corr_can = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_ME_pT_corr_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_ME_pT_corr_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);

      L0bar_L0bar_cosThetaProdPlane_ME_pT_hist[pTbin1][pTbin2] = (TH1D*)inFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
      L0bar_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2] = (TH1D*)inFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_ME_cuts_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));

      L0bar_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->Sumw2();
      L0bar_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->SetNameTitle(Form("L0bar_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_eff", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_eff", pTbin1, pTbin2));
      L0bar_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->Divide(L0bar_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2], L0bar_L0bar_cosThetaProdPlane_ME_pT_hist[pTbin1][pTbin2], 1, 1, "b");
      L0bar_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0bar_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0bar_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->SetLineColor(kRed);
      L0bar_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0bar_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0bar_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->GetYaxis()->SetTitle("Acceptance");
      L0bar_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0bar_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->SetMinimum(0);
      L0bar_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->Draw("p e");
      L0bar_L0bar_cosThetaProdPlane_ME_pT_cuts_hist[pTbin1][pTbin2]->Write();

      L0bar_L0bar_text_pT->Draw("same");

      L0bar_L0bar_cosThetaProdPlane_ME_pT_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/efficiency/Lambda/correlations/L0bar_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_eff.png", pTbin1, pTbin2));
      //_________________________________________________________________________________________________

    }
  }

  for(unsigned int etaBin1 = 0; etaBin1 < nEtaBins; etaBin1++)
  {
    for(unsigned int etaBin2 = 0; etaBin2 < nEtaBins; etaBin2++)
    {
      //signal distributions
      TCanvas *L0_L0bar_cosThetaProdPlane_eta_corr_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_eta_corr_can_eta1_%i_eta2_%i", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_eta_corr_can_eta1_%i_eta2_%i", etaBin1, etaBin2), 1200, 1000);

      L0_L0bar_cosThetaProdPlane_eta_hist[etaBin1][etaBin2] = (TH1D*)inFile->Get(Form("L0_L0bar_cosThetaProdPlane_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));
      L0_L0bar_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2] = (TH1D*)inFile->Get(Form("L0_L0bar_cosThetaProdPlane_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));

      L0_L0bar_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->Sumw2();
      L0_L0bar_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->SetNameTitle(Form("L0_L0bar_cosThetaProdPlane_eta1_%i_eta2_%i_eff", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_eta1_%i_eta2_%i_eff", etaBin1, etaBin2));
      L0_L0bar_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->Divide(L0_L0bar_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2], L0_L0bar_cosThetaProdPlane_eta_hist[etaBin1][etaBin2], 1, 1, "b");
      L0_L0bar_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      L0_L0bar_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->SetMarkerColor(kRed);
      L0_L0bar_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->SetLineColor(kRed);
      L0_L0bar_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0_L0bar_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->GetXaxis()->CenterTitle();
      L0_L0bar_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->GetYaxis()->SetTitle("Acceptance");
      L0_L0bar_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->GetYaxis()->CenterTitle();
      L0_L0bar_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->SetMinimum(0);
      L0_L0bar_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->Draw("p e");
      L0_L0bar_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->Write();

      TPaveText *L0_L0bar_text_eta = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
      L0_L0bar_text_eta->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
      L0_L0bar_text_eta->AddText("Minimum bias");
      L0_L0bar_text_eta->AddText("#Lambda-#bar{#Lambda}");
      //L0_L0bar_text_pt->AddText("No cuts_hist");
      L0_L0bar_text_eta->AddText(Form("%0.1f < #eta_{1} < %0.1f", eta_bins[etaBin1], eta_bins[etaBin1+1]));
      L0_L0bar_text_eta->AddText(Form("%0.1f < #eta_{2} < %0.1f", eta_bins[etaBin2], eta_bins[etaBin2+1]));
      L0_L0bar_text_eta->AddText("p_{T} integrated");
      L0_L0bar_text_eta->SetFillColorAlpha(0, 0.01);
      L0_L0bar_text_eta->Draw("same");

      L0_L0bar_cosThetaProdPlane_eta_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/efficiency/Lambda/correlations/L0_L0bar_cosThetaProdPlane_eta1_%i_eta2_%i_eff.png", etaBin1, etaBin2));
      //_________________________________________________________________________________________________


      TCanvas *L0_L0_cosThetaProdPlane_eta_corr_can = new TCanvas(Form("L0_L0_cosThetaProdPlane_eta_corr_can_eta1_%i_eta2_%i", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_eta_corr_can_eta1_%i_eta2_%i", etaBin1, etaBin2), 1200, 1000);

      L0_L0_cosThetaProdPlane_eta_hist[etaBin1][etaBin2] = (TH1D*)inFile->Get(Form("L0_L0_cosThetaProdPlane_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));
      L0_L0_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2] = (TH1D*)inFile->Get(Form("L0_L0_cosThetaProdPlane_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));

      L0_L0_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->Sumw2();
      L0_L0_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->SetNameTitle(Form("L0_L0_cosThetaProdPlane_eta1_%i_eta2_%i_eff", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_eta1_%i_eta2_%i_eff", etaBin1, etaBin2));
      L0_L0_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->Divide(L0_L0_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2], L0_L0_cosThetaProdPlane_eta_hist[etaBin1][etaBin2], 1, 1, "b");
      L0_L0_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      L0_L0_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->SetMarkerColor(kRed);
      L0_L0_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->SetLineColor(kRed);
      L0_L0_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0_L0_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->GetXaxis()->CenterTitle();
      L0_L0_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->GetYaxis()->SetTitle("Acceptance");
      L0_L0_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->GetYaxis()->CenterTitle();
      L0_L0_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->SetMinimum(0);
      L0_L0_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->Draw("p e");
      L0_L0_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->Write();

      TPaveText *L0_L0_text_eta = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
      L0_L0_text_eta->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
      L0_L0_text_eta->AddText("Minimum bias");
      L0_L0_text_eta->AddText("#Lambda-#Lambda");
      //L0_L0_text_pt->AddText("No cuts_hist");
      L0_L0_text_eta->AddText(Form("%0.1f < #eta_{1} < %0.1f", eta_bins[etaBin1], eta_bins[etaBin1+1]));
      L0_L0_text_eta->AddText(Form("%0.1f < #eta_{2} < %0.1f", eta_bins[etaBin2], eta_bins[etaBin2+1]));
      L0_L0_text_eta->AddText("p_{T} integrated");
      L0_L0_text_eta->SetFillColorAlpha(0, 0.01);
      L0_L0_text_eta->Draw("same");

      L0_L0_cosThetaProdPlane_eta_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/efficiency/Lambda/correlations/L0_L0_cosThetaProdPlane_eta1_%i_eta2_%i_eff.png", etaBin1, etaBin2));
      //_________________________________________________________________________________________________

      TCanvas *L0bar_L0bar_cosThetaProdPlane_eta_corr_can = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_eta_corr_can_eta1_%i_eta2_%i", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_eta_corr_can_eta1_%i_eta2_%i", etaBin1, etaBin2), 1200, 1000);

      L0bar_L0bar_cosThetaProdPlane_eta_hist[etaBin1][etaBin2] = (TH1D*)inFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));
      L0bar_L0bar_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2] = (TH1D*)inFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));

      L0bar_L0bar_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->Sumw2();
      L0bar_L0bar_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->SetNameTitle(Form("L0bar_L0bar_cosThetaProdPlane_eta1_%i_eta2_%i_eff", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_eta1_%i_eta2_%i_eff", etaBin1, etaBin2));
      L0bar_L0bar_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->Divide(L0bar_L0bar_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2], L0bar_L0bar_cosThetaProdPlane_eta_hist[etaBin1][etaBin2], 1, 1, "b");
      L0bar_L0bar_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      L0bar_L0bar_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->SetMarkerColor(kRed);
      L0bar_L0bar_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->SetLineColor(kRed);
      L0bar_L0bar_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0bar_L0bar_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->GetXaxis()->CenterTitle();
      L0bar_L0bar_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->GetYaxis()->SetTitle("Acceptance");
      L0bar_L0bar_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->GetYaxis()->CenterTitle();
      L0bar_L0bar_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->SetMinimum(0);
      L0bar_L0bar_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->Draw("p e");
      L0bar_L0bar_cosThetaProdPlane_eta_cuts_hist[etaBin1][etaBin2]->Write();

      TPaveText *L0bar_L0bar_text_eta = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
      L0bar_L0bar_text_eta->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
      L0bar_L0bar_text_eta->AddText("Minimum bias");
      L0bar_L0bar_text_eta->AddText("#bar{#Lambda}-#bar{#Lambda}");
      //L0bar_L0bar_text_pt->AddText("No cuts_hist");
      L0bar_L0bar_text_eta->AddText(Form("%0.1f < #eta_{1} < %0.1f", eta_bins[etaBin1], eta_bins[etaBin1+1]));
      L0bar_L0bar_text_eta->AddText(Form("%0.1f < #eta_{2} < %0.1f", eta_bins[etaBin2], eta_bins[etaBin2+1]));
      L0bar_L0bar_text_eta->AddText("p_{T} integrated");
      L0bar_L0bar_text_eta->SetFillColorAlpha(0, 0.01);
      L0bar_L0bar_text_eta->Draw("same");

      L0bar_L0bar_cosThetaProdPlane_eta_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/efficiency/Lambda/correlations/L0bar_L0bar_cosThetaProdPlane_eta1_%i_eta2_%i_eff.png", etaBin1, etaBin2));
      //_________________________________________________________________________________________________

      //mixed event distributins
      TCanvas *L0_L0bar_cosThetaProdPlane_ME_eta_corr_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_ME_eta_corr_can_eta1_%i_eta2_%i", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_ME_eta_corr_can_eta1_%i_eta2_%i", etaBin1, etaBin2), 1200, 1000);

      L0_L0bar_cosThetaProdPlane_ME_eta_hist[etaBin1][etaBin2] = (TH1D*)inFile->Get(Form("L0_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));
      L0_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2] = (TH1D*)inFile->Get(Form("L0_L0bar_cosThetaProdPlane_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));

      L0_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->Sumw2();
      L0_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->SetNameTitle(Form("L0_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_eff", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_eff", etaBin1, etaBin2));
      L0_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->Divide(L0_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2], L0_L0bar_cosThetaProdPlane_ME_eta_hist[etaBin1][etaBin2], 1, 1, "b");
      L0_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      L0_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->SetMarkerColor(kRed);
      L0_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->SetLineColor(kRed);
      L0_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->GetXaxis()->CenterTitle();
      L0_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->GetYaxis()->SetTitle("Acceptance");
      L0_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->GetYaxis()->CenterTitle();
      L0_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->SetMinimum(0);
      L0_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->Draw("p e");
      L0_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->Write();

      L0_L0bar_text_eta->Draw("same");

      L0_L0bar_cosThetaProdPlane_ME_eta_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/efficiency/Lambda/correlations/L0_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_eff.png", etaBin1, etaBin2));
      //_________________________________________________________________________________________________


      TCanvas *L0_L0_cosThetaProdPlane_ME_eta_corr_can = new TCanvas(Form("L0_L0_cosThetaProdPlane_ME_eta_corr_can_eta1_%i_eta2_%i", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_ME_eta_corr_can_eta1_%i_eta2_%i", etaBin1, etaBin2), 1200, 1000);

      L0_L0_cosThetaProdPlane_ME_eta_hist[etaBin1][etaBin2] = (TH1D*)inFile->Get(Form("L0_L0_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));
      L0_L0_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2] = (TH1D*)inFile->Get(Form("L0_L0_cosThetaProdPlane_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));

      L0_L0_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->Sumw2();
      L0_L0_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->SetNameTitle(Form("L0_L0_cosThetaProdPlane_ME_eta1_%i_eta2_%i_eff", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_ME_eta1_%i_eta2_%i_eff", etaBin1, etaBin2));
      L0_L0_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->Divide(L0_L0_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2], L0_L0_cosThetaProdPlane_ME_eta_hist[etaBin1][etaBin2], 1, 1, "b");
      L0_L0_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      L0_L0_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->SetMarkerColor(kRed);
      L0_L0_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->SetLineColor(kRed);
      L0_L0_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0_L0_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->GetXaxis()->CenterTitle();
      L0_L0_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->GetYaxis()->SetTitle("Acceptance");
      L0_L0_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->GetYaxis()->CenterTitle();
      L0_L0_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->SetMinimum(0);
      L0_L0_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->Draw("p e");
      L0_L0_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->Write();


      L0_L0_text_eta->Draw("same");

      L0_L0_cosThetaProdPlane_ME_eta_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/efficiency/Lambda/correlations/L0_L0_cosThetaProdPlane_ME_eta1_%i_eta2_%i_eff.png", etaBin1, etaBin2));
      //_________________________________________________________________________________________________

      TCanvas *L0bar_L0bar_cosThetaProdPlane_ME_eta_corr_can = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_ME_eta_corr_can_eta1_%i_eta2_%i", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_ME_eta_corr_can_eta1_%i_eta2_%i", etaBin1, etaBin2), 1200, 1000);

      L0bar_L0bar_cosThetaProdPlane_ME_eta_hist[etaBin1][etaBin2] = (TH1D*)inFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));
      L0bar_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2] = (TH1D*)inFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_ME_cuts_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));

      L0bar_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->Sumw2();
      L0bar_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->SetNameTitle(Form("L0bar_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_eff", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_eff", etaBin1, etaBin2));
      L0bar_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->Divide(L0bar_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2], L0bar_L0bar_cosThetaProdPlane_ME_eta_hist[etaBin1][etaBin2], 1, 1, "b");
      L0bar_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      L0bar_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->SetMarkerColor(kRed);
      L0bar_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->SetLineColor(kRed);
      L0bar_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0bar_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->GetXaxis()->CenterTitle();
      L0bar_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->GetYaxis()->SetTitle("Acceptance");
      L0bar_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->GetYaxis()->CenterTitle();
      L0bar_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->SetMinimum(0);
      L0bar_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->Draw("p e");
      L0bar_L0bar_cosThetaProdPlane_ME_eta_cuts_hist[etaBin1][etaBin2]->Write();

      L0bar_L0bar_text_eta->Draw("same");

      L0bar_L0bar_cosThetaProdPlane_ME_eta_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/efficiency/Lambda/correlations/L0bar_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_eff.png", etaBin1, etaBin2));
      //_________________________________________________________________________________________________
    }
  }

  inFile->Close();

  outFile->Close();

  return;

}
