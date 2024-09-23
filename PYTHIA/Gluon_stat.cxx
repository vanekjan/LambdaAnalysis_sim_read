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
void Gluon_stat()
{
  TFile *inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/PYTHIA_tree/output_Lambda_pp_200_MB_10M_events_0.root", "READ");


  TH1F *has_g_to_s_sbar_hist = (TH1F*)inFile->Get("has_g_to_s_sbar_hist"); //number of g -> s-sbar in an event

  TH1F *has_L_Lbar_hist = (TH1F*)inFile->Get("has_L_Lbar_hist"); //number of events that have g -> s-sbar and also have L Lbar
  TH1F *has_L_Lbar_from_ssBar_hist = (TH1F*)inFile->Get("has_L_Lbar_from_ssBar_hist");
  TH1F *n_L_Lbar_tot_hist = (TH1F*)inFile->Get("n_L_Lbar_tot_hist"); //number of events with L-Lbar vs. without L-Lbar pari

  //n L pairs
  TH1F *n_L_hist = (TH1F*)inFile->Get("n_L_hist"); //number of L from g to s-sBar in an event
  TH1F *n_Lbar_hist = (TH1F*)inFile->Get("n_Lbar_hist"); //number of L from g to s-sBar in an event

  TH1F *n_L_from_ssbar_hist = (TH1F*)inFile->Get("n_L_from_ssbar_hist"); //number of L from g to s-sBar in an event
  TH1F *n_Lbar_from_ssbar_hist = (TH1F*)inFile->Get("n_Lbar_from_ssbar_hist"); //number of L from g to s-sBar in an event

  //s kine. vs. L kine.
  TH2F *s_vs_L_phi_hist = (TH2F*)inFile->Get("s_vs_L_phi_hist");
  TH2F *s_vs_L_eta_hist = (TH2F*)inFile->Get("s_vs_L_eta_hist");

  TH2F *sBar_vs_Lbar_phi_hist = (TH2F*)inFile->Get("sBar_vs_Lbar_phi_hist");
  TH2F *sBar_vs_Lbar_eta_hist = (TH2F*)inFile->Get("sBar_vs_Lbar_eta_hist");

  //--------------------------------------------------------------------------------

  TCanvas *nL_can = new TCanvas("nL_can", "nL_can", 2000, 2000);
  nL_can->Divide(2,2);

  nL_can->cd(1);

  gPad->SetLogy();

  n_L_hist->GetXaxis()->SetRangeUser(0,10);
  n_L_hist->Sumw2();
  n_L_hist->Draw("p e");


  nL_can->cd(2);

  gPad->SetLogy();

  n_L_from_ssbar_hist->GetXaxis()->SetRangeUser(0,10);
  n_L_from_ssbar_hist->Draw("p e");


  nL_can->cd(3);

  TH1F* n_L_ratio = (TH1F*)n_L_from_ssbar_hist->Clone("n_L_ratio");
  n_L_ratio->Sumw2();
  n_L_ratio->Divide(n_L_hist);
  n_L_ratio->GetXaxis()->SetRangeUser(0,10);
  n_L_ratio->Draw("p e");

  nL_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/PYTHIA_gluon/nL_can.png");

  //--------------------------------------------------------------------------------

  TCanvas *nLbar_can = new TCanvas("nLbar_can", "nLbar_can", 2000, 2000);
  nLbar_can->Divide(2,2);

  nLbar_can->cd(1);

  gPad->SetLogy();

  n_Lbar_hist->GetXaxis()->SetRangeUser(0,10);
  n_Lbar_hist->Sumw2();
  n_Lbar_hist->Draw("p e");


  nLbar_can->cd(2);

  gPad->SetLogy();

  n_Lbar_from_ssbar_hist->GetXaxis()->SetRangeUser(0,10);
  n_Lbar_from_ssbar_hist->Draw("p e");


  nLbar_can->cd(3);

  TH1F* n_Lbar_ratio = (TH1F*)n_Lbar_from_ssbar_hist->Clone("n_Lbar_ratio");
  n_Lbar_ratio->Sumw2();
  n_Lbar_ratio->Divide(n_Lbar_hist);
  n_Lbar_ratio->GetXaxis()->SetRangeUser(0,10);
  n_Lbar_ratio->Draw("p e");

  nLbar_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/PYTHIA_gluon/nLbar_can.png");

  //--------------------------------------------------------------------------------

  TCanvas *nKine_can = new TCanvas("nKine_can", "nKine_can", 2000, 2000);
  nKine_can->Divide(2,2);

  nKine_can->cd(1);

  s_vs_L_phi_hist->Draw("colz");


  nKine_can->cd(2);

  s_vs_L_eta_hist->Draw("colz");


  nKine_can->cd(3);

  sBar_vs_Lbar_phi_hist->Draw("colz");


  nKine_can->cd(4);

  sBar_vs_Lbar_eta_hist->Draw("colz");

  nKine_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/PYTHIA_gluon/s_vs_L_kine.png");


  cout<<"Number of g -> s-sbar: "<<has_g_to_s_sbar_hist->GetBinContent(2)<<endl;
  cout<<"Number of events with g -> s-sbar that also have L-Lbar: "<<has_L_Lbar_hist->GetBinContent(2)<<endl;
  cout<<"Number of events with L-Lbar that can be traced to g -> s-sbar: "<<has_L_Lbar_from_ssBar_hist->GetBinContent(2)<<endl;
  cout<<"Number of events with at least one L-Lbar pair: "<<n_L_Lbar_tot_hist->GetBinContent(2)<<endl;

  inFile->Close();

  return;
}
