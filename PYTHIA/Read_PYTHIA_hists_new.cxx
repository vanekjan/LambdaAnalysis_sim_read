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
void Read_PYTHIA_hists_new(const int energy = 510)
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
  else if(energy == 200)
  {
    //inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/output_Lambda_pp_200_MB_1B_events_new.root", "READ");

    inFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/input/Run12/Mother_ID/output_Lambda_pp_200_MB_1B_events_hists_mother_ID.root", "READ");
  }
  else
  {
    cout<<"Invalid energy!"<<endl;

    return;
  }

  const int nPtBins = 8;
  float const pT_bins[nPtBins+1] = { 0., 0.5, 1.,1.5, 2., 2.5, 3., 4., 5.};

  const int nEtaBins = 3;
  float const eta_bins[nEtaBins+1] = { -1, -0.4, 0.4, 1 };
  //______________________________________________________________________________________________

  //histograms
  //TH1D* L0_pt_hist = (TH1D*)inFile->Get("L0_pt_hist");
	//TH1D* L0_mass_hist = (TH1D*)inFile->Get("L0_mass_hist");

 // TH2D* L0_pT_vs_L0_y = (TH2D*)inFile->Get("L0_pT_vs_L0_y");

  TH1D *L0_thetaProdPlane_hist[nPtBins+1][nEtaBins+1];
  TH1D *L0_cosThetaProdPlane_hist[nPtBins+1][nEtaBins+1];

  TH2D *L0_y_vs_p_eta[nPtBins+1];
  TH2D *L0_y_vs_pi_eta[nPtBins+1];

  TH1D *L0_pz[nPtBins+1][nEtaBins+1];
  TH1D *L0_xF[nPtBins+1][nEtaBins+1];

  TH2D *L0_pT_vs_L0_pz[nEtaBins+1];

  TH2D *L0_p_eta_vs_pi_eta[nPtBins+1];



  //TH1D* L0bar_pt_hist = (TH1D*)inFile->Get("L0bar_pt_hist");
	//TH1D* L0bar_mass_hist = (TH1D*)inFile->Get("L0bar_mass_hist");

	//TH2D* L0bar_pT_vs_L0bar_y = (TH2D*)inFile->Get("L0bar_pT_vs_L0bar_y");

  TH1D *L0bar_thetaProdPlane_hist[nPtBins+1][nEtaBins+1];
  TH1D *L0bar_cosThetaProdPlane_hist[nPtBins+1][nEtaBins+1];

  TH2D *L0bar_y_vs_p_eta[nPtBins+1];
  TH2D *L0bar_y_vs_pi_eta[nPtBins+1];

  TH1D *L0bar_pz[nPtBins+1][nEtaBins+1];
  TH1D *L0bar_xF[nPtBins+1][nEtaBins+1];

  TH2D *L0bar_pT_vs_L0bar_pz[nEtaBins+1];

  TH2D *L0bar_p_eta_vs_pi_eta[nPtBins+1];


  TH1D *L0_L0bar_cosThetaProdPlane = (TH1D*)inFile->Get("L0_L0bar_cosThetaProdPlane");
  TH1D *L0_L0_cosThetaProdPlane = (TH1D*)inFile->Get("L0_L0_cosThetaProdPlane");
  TH1D *L0bar_L0bar_cosThetaProdPlane = (TH1D*)inFile->Get("L0bar_L0bar_cosThetaProdPlane");

  //L-L correlation QA histograms
  TH2D *L0_L0bar_delta_eta_vs_delta_phi_hist = (TH2D*)inFile->Get("L0_L0bar_delta_eta_vs_delta_phi_hist");

  TH2D *L0_L0bar_y1_vs_y2_hist = (TH2D*)inFile->Get("L0_L0bar_y1_vs_y2_hist");

  TH2D *L0_L0bar_pT1_vs_pT2_hist = (TH2D*)inFile->Get("L0_L0bar_pT1_vs_pT2_hist");

  TH2D *L0_L0bar_phi1_vs_phi2_hist = (TH2D*)inFile->Get("L0_L0bar_phi1_vs_phi2_hist");


  TH2D *L0_L0_delta_eta_vs_delta_phi_hist = (TH2D*)inFile->Get("L0_L0_delta_eta_vs_delta_phi_hist");

  TH2D *L0_L0_delta_eta_vs_delta_phi_zoom_hist = (TH2D*)inFile->Get("L0_L0_delta_eta_vs_delta_phi_zoom_hist");

  TH2D *L0_L0_y1_vs_y2_hist = (TH2D*)inFile->Get("L0_L0_y1_vs_y2_hist");

  TH2D *L0_p_L0_p_y1_vs_y2_hist = (TH2D*)inFile->Get("L0_p_L0_p_y1_vs_y2_hist");

  TH2D *L0_pi_L0_pi_y1_vs_y2_hist = (TH2D*)inFile->Get("L0_pi_L0_pi_y1_vs_y2_hist");

  TH2D *L0_L0_pT1_vs_pT2_hist = (TH2D*)inFile->Get("L0_L0_pT1_vs_pT2_hist");

  TH2D *L0_L0_phi1_vs_phi2_hist = (TH2D*)inFile->Get("L0_L0_phi1_vs_phi2_hist");

  TH2D *L0_p_L0_p_phi1_vs_phi2_hist = (TH2D*)inFile->Get("L0_p_L0_p_phi1_vs_phi2_hist");

  TH2D *L0_pi_L0_pi_phi1_vs_phi2_hist = (TH2D*)inFile->Get("L0_pi_L0_pi_phi1_vs_phi2_hist");


  TH2D *L0bar_L0bar_delta_eta_vs_delta_phi_hist = (TH2D*)inFile->Get("L0bar_L0bar_delta_eta_vs_delta_phi_hist");

  TH2D *L0bar_L0bar_delta_eta_vs_delta_phi_zoom_hist = (TH2D*)inFile->Get("L0bar_L0bar_delta_eta_vs_delta_phi_zoom_hist");

  TH2D *L0bar_L0bar_y1_vs_y2_hist = (TH2D*)inFile->Get("L0bar_L0bar_y1_vs_y2_hist");

  TH2D *L0bar_p_L0bar_p_y1_vs_y2_hist = (TH2D*)inFile->Get("L0bar_p_L0bar_p_y1_vs_y2_hist");

  TH2D *L0bar_pi_L0bar_pi_y1_vs_y2_hist = (TH2D*)inFile->Get("L0bar_pi_L0bar_pi_y1_vs_y2_hist");

  TH2D *L0bar_L0bar_pT1_vs_pT2_hist = (TH2D*)inFile->Get("L0bar_L0bar_pT1_vs_pT2_hist");

  TH2D *L0bar_L0bar_phi1_vs_phi2_hist = (TH2D*)inFile->Get("L0bar_L0bar_phi1_vs_phi2_hist");

  TH2D *L0bar_p_L0bar_p_phi1_vs_phi2_hist = (TH2D*)inFile->Get("L0bar_p_L0bar_p_phi1_vs_phi2_hist");

  TH2D *L0bar_pi_L0bar_pi_phi1_vs_phi2_hist = (TH2D*)inFile->Get("L0bar_pi_L0bar_pi_phi1_vs_phi2_hist");

  //L-l correlations after cuts
  TH1D *L0_L0bar_cosThetaProdPlane_cuts = (TH1D*)inFile->Get("L0_L0bar_cosThetaProdPlane_cuts");
  TH1D *L0_L0_cosThetaProdPlane_cuts = (TH1D*)inFile->Get("L0_L0_cosThetaProdPlane_cuts");
  TH1D *L0bar_L0bar_cosThetaProdPlane_cuts = (TH1D*)inFile->Get("L0bar_L0bar_cosThetaProdPlane_cuts");


  //L-L correlation QA histograms
  TH2D *L0_L0bar_delta_eta_vs_delta_phi_cuts_hist = (TH2D*)inFile->Get("L0_L0bar_delta_eta_vs_delta_phi_cuts_hist");

  TH2D *L0_L0bar_y1_vs_y2_cuts_hist = (TH2D*)inFile->Get("L0_L0bar_y1_vs_y2_cuts_hist");

  TH2D *L0_L0bar_pT1_vs_pT2_cuts_hist = (TH2D*)inFile->Get("L0_L0bar_pT1_vs_pT2_cuts_hist");

  TH2D *L0_L0bar_phi1_vs_phi2_cuts_hist = (TH2D*)inFile->Get("L0_L0bar_phi1_vs_phi2_cuts_hist");


  TH2D *L0_L0_delta_eta_vs_delta_phi_cuts_hist = (TH2D*)inFile->Get("L0_L0_delta_eta_vs_delta_phi_cuts_hist");

  TH2D *L0_L0_delta_eta_vs_delta_phi_zoom_cuts_hist = (TH2D*)inFile->Get("L0_L0_delta_eta_vs_delta_phi_zoom_cuts_hist");

  TH2D *L0_L0_y1_vs_y2_cuts_hist = (TH2D*)inFile->Get("L0_L0_y1_vs_y2_cuts_hist");

  TH2D *L0_p_L0_p_y1_vs_y2_cuts_hist = (TH2D*)inFile->Get("L0_p_L0_p_y1_vs_y2_cuts_hist");

  TH2D *L0_pi_L0_pi_y1_vs_y2_cuts_hist = (TH2D*)inFile->Get("L0_pi_L0_pi_y1_vs_y2_cuts_hist");

  TH2D *L0_L0_pT1_vs_pT2_cuts_hist = (TH2D*)inFile->Get("L0_L0_pT1_vs_pT2_cuts_hist");

  TH2D *L0_L0_phi1_vs_phi2_cuts_hist = (TH2D*)inFile->Get("L0_L0_phi1_vs_phi2_cuts_hist");

  TH2D *L0_p_L0_p_phi1_vs_phi2_cuts_hist = (TH2D*)inFile->Get("L0_p_L0_p_phi1_vs_phi2_cuts_hist");

  TH2D *L0_pi_L0_pi_phi1_vs_phi2_cuts_hist = (TH2D*)inFile->Get("L0_pi_L0_pi_phi1_vs_phi2_cuts_hist");


  TH2D *L0bar_L0bar_delta_eta_vs_delta_phi_cuts_hist = (TH2D*)inFile->Get("L0bar_L0bar_delta_eta_vs_delta_phi_cuts_hist");

  TH2D *L0bar_L0bar_delta_eta_vs_delta_phi_zoom_cuts_hist = (TH2D*)inFile->Get("L0bar_L0bar_delta_eta_vs_delta_phi_zoom_cuts_hist");

  TH2D *L0bar_L0bar_y1_vs_y2_cuts_hist = (TH2D*)inFile->Get("L0bar_L0bar_y1_vs_y2_cuts_hist");

  TH2D *L0bar_p_L0bar_p_y1_vs_y2_cuts_hist = (TH2D*)inFile->Get("L0bar_p_L0bar_p_y1_vs_y2_cuts_hist");

  TH2D *L0bar_pi_L0bar_pi_y1_vs_y2_cuts_hist = (TH2D*)inFile->Get("L0bar_pi_L0bar_pi_y1_vs_y2_cuts_hist");

  TH2D *L0bar_L0bar_pT1_vs_pT2_cuts_hist = (TH2D*)inFile->Get("L0bar_L0bar_pT1_vs_pT2_cuts_hist");

  TH2D *L0bar_L0bar_phi1_vs_phi2_cuts_hist = (TH2D*)inFile->Get("L0bar_L0bar_phi1_vs_phi2_cuts_hist");

  TH2D *L0bar_p_L0bar_p_phi1_vs_phi2_cuts_hist = (TH2D*)inFile->Get("L0bar_p_L0bar_p_phi1_vs_phi2_cuts_hist");

  TH2D *L0bar_pi_L0bar_pi_phi1_vs_phi2_cuts_hist = (TH2D*)inFile->Get("L0bar_pi_L0bar_pi_phi1_vs_phi2_cuts_hist");

  //________________________________________________________________________________________________


  //TCanvas *L0_pT_vs_L0_y_can = new TCanvas("L0_pT_vs_L0_y_can", "L0_pT_vs_L0_y_can", 1200, 1000);

  TCanvas *L0_thetaProdPlane_can[nPtBins+1][nEtaBins+1];
  TCanvas *L0_cosThetaProdPlane_can[nPtBins+1][nEtaBins+1];

  TCanvas *L0_y_vs_p_eta_can[nPtBins+1];
  TCanvas *L0_y_vs_pi_eta_can[nPtBins+1];

  TCanvas *L0_pz_can[nPtBins+1][nEtaBins+1];
  TCanvas *L0_xF_can[nPtBins+1][nEtaBins+1];

  TCanvas *L0_pT_vs_L0_pz_can[nEtaBins+1];


  TCanvas *L0_p_eta_vs_pi_eta_can[nPtBins+1];


  //TCanvas *L0bar_pT_vs_L0bar_y_can = new TCanvas("L0bar_pT_vs_L0bar_y_can", "L0bar_pT_vs_L0bar_y_can", 1200, 1000);

  TCanvas *L0bar_thetaProdPlane_can[nPtBins+1][nEtaBins+1];
  TCanvas *L0bar_cosThetaProdPlane_can[nPtBins+1][nEtaBins+1];

  TCanvas *L0bar_y_vs_p_eta_can[nPtBins+1];
  TCanvas *L0bar_y_vs_pi_eta_can[nPtBins+1];

  TCanvas *L0bar_pz_can[nPtBins+1][nEtaBins+1];
  TCanvas *L0bar_xF_can[nPtBins+1][nEtaBins+1];

  TCanvas *L0bar_pT_vs_L0bar_pz_can[nEtaBins+1];

  TCanvas *L0bar_p_eta_vs_pi_eta_can[nPtBins+1];



  TCanvas *L0_L0bar_cosThetaProdPlane_can = new TCanvas("L0_L0bar_cosThetaProdPlane_can", "L0_L0bar_cosThetaProdPlane_can", 1300, 1000);
  TCanvas *L0_L0_cosThetaProdPlane_can = new TCanvas("L0_L0_cosThetaProdPlane_can", "L0_L0_cosThetaProdPlane_can", 1300, 1000);
  TCanvas *L0bar_L0bar_cosThetaProdPlane_can = new TCanvas("L0bar_L0bar_cosThetaProdPlane_can", "L0bar_L0bar_cosThetaProdPlane_can", 1300, 1000);


  //L-L correlation QA histograms
  TCanvas *L0_L0bar_delta_eta_vs_delta_phi_can = new TCanvas("L0_L0bar_delta_eta_vs_delta_phi_can", "L0_L0bar_delta_eta_vs_delta_phi_can", 1200, 1000);

  TCanvas *L0_L0bar_y1_vs_y2_can = new TCanvas("L0_L0bar_y1_vs_y2_can", "L0_L0bar_y1_vs_y2_can", 1200, 1000);

  TCanvas *L0_L0bar_pT1_vs_pT2_can = new TCanvas("L0_L0bar_pT1_vs_pT2_can", "L0_L0bar_pT1_vs_pT2_can", 1200, 1000);

  TCanvas *L0_L0bar_phi1_vs_phi2_can = new TCanvas("L0_L0bar_phi1_vs_phi2_can", "L0_L0bar_phi1_vs_phi2_can", 1200, 1000);


  TCanvas *L0_L0_delta_eta_vs_delta_phi_can = new TCanvas("L0_L0_delta_eta_vs_delta_phi_can", "L0_L0_delta_eta_vs_delta_phi_can", 1200, 1000);

  TCanvas *L0_L0_delta_eta_vs_delta_phi_zoom_can = new TCanvas("L0_L0_delta_eta_vs_delta_phi_zoom_can", "L0_L0_delta_eta_vs_delta_phi_zoom_can", 1200, 1000);

  TCanvas *L0_L0_y1_vs_y2_can = new TCanvas("L0_L0_y1_vs_y2_can", "L0_L0_y1_vs_y2_can", 1200, 1000);

  TCanvas *L0_p_L0_p_y1_vs_y2_can = new TCanvas("L0_p_L0_p_y1_vs_y2_can", "L0_p_L0_p_y1_vs_y2_can", 1200, 1000);

  TCanvas *L0_pi_L0_pi_y1_vs_y2_can = new TCanvas("L0_pi_L0_pi_y1_vs_y2_can", "L0_pi_L0_pi_y1_vs_y2_can", 1200, 1000);

  TCanvas *L0_L0_pT1_vs_pT2_can = new TCanvas("L0_L0_pT1_vs_pT2_can", "L0_L0_pT1_vs_pT2_can", 1200, 1000);

  TCanvas *L0_L0_phi1_vs_phi2_can = new TCanvas("L0_L0_phi1_vs_phi2_can", "L0_L0_phi1_vs_phi2_can", 1200, 1000);

  TCanvas *L0_p_L0_p_phi1_vs_phi2_can = new TCanvas("L0_p_L0_p_phi1_vs_phi2_can", "L0_p_L0_p_phi1_vs_phi2_can", 1200, 1000);

  TCanvas *L0_pi_L0_pi_phi1_vs_phi2_can = new TCanvas("L0_pi_L0_pi_phi1_vs_phi2_can", "L0_pi_L0_pi_phi1_vs_phi2_can", 1200, 1000);


  TCanvas *L0bar_L0bar_delta_eta_vs_delta_phi_can = new TCanvas("L0bar_L0bar_delta_eta_vs_delta_phi_can", "L0bar_L0bar_delta_eta_vs_delta_phi_can", 1200, 1000);

  TCanvas *L0bar_L0bar_delta_eta_vs_delta_phi_zoom_can = new TCanvas("L0bar_L0bar_delta_eta_vs_delta_phi_zoom_can", "L0bar_L0bar_delta_eta_vs_delta_phi_zoom_can", 1200, 1000);

  TCanvas *L0bar_L0bar_y1_vs_y2_can = new TCanvas("L0bar_L0bar_y1_vs_y2_can", "L0bar_L0bar_y1_vs_y2_can", 1200, 1000);

  TCanvas *L0bar_p_L0bar_p_y1_vs_y2_can = new TCanvas("L0bar_p_L0bar_p_y1_vs_y2_can", "L0bar_p_L0bar_p_y1_vs_y2_can", 1200, 1000);

  TCanvas *L0bar_pi_L0bar_pi_y1_vs_y2_can = new TCanvas("L0bar_pi_L0bar_pi_y1_vs_y2_can", "L0bar_pi_L0bar_pi_y1_vs_y2_can", 1200, 1000);

  TCanvas *L0bar_L0bar_pT1_vs_pT2_can = new TCanvas("L0bar_L0bar_pT1_vs_pT2_can", "L0bar_L0bar_pT1_vs_pT2_can", 1200, 1000);

  TCanvas *L0bar_L0bar_phi1_vs_phi2_can = new TCanvas("L0bar_L0bar_phi1_vs_phi2_can", "L0bar_L0bar_phi1_vs_phi2_can", 1200, 1000);

  TCanvas *L0bar_p_L0bar_p_phi1_vs_phi2_can = new TCanvas("L0bar_p_L0bar_p_phi1_vs_phi2_can", "L0bar_p_L0bar_p_phi1_vs_phi2_can", 1200, 1000);

  TCanvas *L0bar_pi_L0bar_pi_phi1_vs_phi2_can = new TCanvas("L0bar_pi_L0bar_pi_phi1_vs_phi2_can", "L0bar_pi_L0bar_pi_phi1_vs_phi2_can", 1200, 1000);


  TCanvas *L0_L0bar_cosThetaProdPlane_cuts_can = new TCanvas("L0_L0bar_cosThetaProdPlane_cuts_can", "L0_L0bar_cosThetaProdPlane_cuts_can", 1200, 1000);
  TCanvas *L0_L0_cosThetaProdPlane_cuts_can = new TCanvas("L0_L0_cosThetaProdPlane_cuts_can", "L0_L0_cosThetaProdPlane_cuts_can", 1200, 1000);
  TCanvas *L0bar_L0bar_cosThetaProdPlane_cuts_can = new TCanvas("L0bar_L0bar_cosThetaProdPlane_cuts_can", "L0bar_L0bar_cosThetaProdPlane_cuts_can", 1200, 1000);

  //L-L correlation QA histograms
  TCanvas *L0_L0bar_delta_eta_vs_delta_phi_cuts_can = new TCanvas("L0_L0bar_delta_eta_vs_delta_phi_cuts_can", "L0_L0bar_delta_eta_vs_delta_phi_cuts_can", 1200, 1000);

  TCanvas *L0_L0bar_y1_vs_y2_cuts_can = new TCanvas("L0_L0bar_y1_vs_y2_cuts_can", "L0_L0bar_y1_vs_y2_cuts_can", 1200, 1000);

  TCanvas *L0_L0bar_pT1_vs_pT2_cuts_can = new TCanvas("L0_L0bar_pT1_vs_pT2_cuts_can", "L0_L0bar_pT1_vs_pT2_cuts_can", 1200, 1000);

  TCanvas *L0_L0bar_phi1_vs_phi2_cuts_can = new TCanvas("L0_L0bar_phi1_vs_phi2_cuts_can", "L0_L0bar_phi1_vs_phi2_cuts_can", 1200, 1000);


  TCanvas *L0_L0_delta_eta_vs_delta_phi_cuts_can = new TCanvas("L0_L0_delta_eta_vs_delta_phi_cuts_can", "L0_L0_delta_eta_vs_delta_phi_cuts_can", 1200, 1000);

  TCanvas *L0_L0_delta_eta_vs_delta_phi_zoom_cuts_can = new TCanvas("L0_L0_delta_eta_vs_delta_phi_zoom_cuts_can", "L0_L0_delta_eta_vs_delta_phi_zoom_cuts_can", 1200, 1000);

  TCanvas *L0_L0_y1_vs_y2_cuts_can = new TCanvas("L0_L0_y1_vs_y2_cuts_can", "L0_L0_y1_vs_y2_cuts_can", 1200, 1000);

  TCanvas *L0_p_L0_p_y1_vs_y2_cuts_can = new TCanvas("L0_p_L0_p_y1_vs_y2_cuts_can", "L0_p_L0_p_y1_vs_y2_cuts_can", 1200, 1000);

  TCanvas *L0_pi_L0_pi_y1_vs_y2_cuts_can = new TCanvas("L0_pi_L0_pi_y1_vs_y2_cuts_can", "L0_pi_L0_pi_y1_vs_y2_cuts_can", 1200, 1000);

  TCanvas *L0_L0_pT1_vs_pT2_cuts_can = new TCanvas("L0_L0_pT1_vs_pT2_cuts_can", "L0_L0_pT1_vs_pT2_cuts_can", 1200, 1000);

  TCanvas *L0_L0_phi1_vs_phi2_cuts_can = new TCanvas("L0_L0_phi1_vs_phi2_cuts_can", "L0_L0_phi1_vs_phi2_cuts_can", 1200, 1000);

  TCanvas *L0_p_L0_p_phi1_vs_phi2_cuts_can = new TCanvas("L0_p_L0_p_phi1_vs_phi2_cuts_can", "L0_p_L0_p_phi1_vs_phi2_cuts_can", 1200, 1000);

  TCanvas *L0_pi_L0_pi_phi1_vs_phi2_cuts_can = new TCanvas("L0_pi_L0_pi_phi1_vs_phi2_cuts_can", "L0_pi_L0_pi_phi1_vs_phi2_cuts_can", 1200, 1000);


  TCanvas *L0bar_L0bar_delta_eta_vs_delta_phi_cuts_can = new TCanvas("L0bar_L0bar_delta_eta_vs_delta_phi_cuts_can", "L0bar_L0bar_delta_eta_vs_delta_phi_cuts_can", 1200, 1000);

  TCanvas *L0bar_L0bar_delta_eta_vs_delta_phi_zoom_cuts_can = new TCanvas("L0bar_L0bar_delta_eta_vs_delta_phi_zoom_cuts_can", "L0bar_L0bar_delta_eta_vs_delta_phi_zoom_cuts_can", 1200, 1000);

  TCanvas *L0bar_L0bar_y1_vs_y2_cuts_can = new TCanvas("L0bar_L0bar_y1_vs_y2_cuts_can", "L0bar_L0bar_y1_vs_y2_cuts_can", 1200, 1000);

  TCanvas *L0bar_p_L0bar_p_y1_vs_y2_cuts_can = new TCanvas("L0bar_p_L0bar_p_y1_vs_y2_cuts_can", "L0bar_p_L0bar_p_y1_vs_y2_cuts_can", 1200, 1000);

  TCanvas *L0bar_pi_L0bar_pi_y1_vs_y2_cuts_can = new TCanvas("L0bar_pi_L0bar_pi_y1_vs_y2_cuts_can", "L0bar_pi_L0bar_pi_y1_vs_y2_cuts_can", 1200, 1000);

  TCanvas *L0bar_L0bar_pT1_vs_pT2_cuts_can = new TCanvas("L0bar_L0bar_pT1_vs_pT2_cuts_can", "L0bar_L0bar_pT1_vs_pT2_cuts_can", 1200, 1000);

  TCanvas *L0bar_L0bar_phi1_vs_phi2_cuts_can = new TCanvas("L0bar_L0bar_phi1_vs_phi2_cuts_can", "L0bar_L0bar_phi1_vs_phi2_cuts_can", 1200, 1000);

  TCanvas *L0bar_p_L0bar_p_phi1_vs_phi2_cuts_can = new TCanvas("L0bar_p_L0bar_p_phi1_vs_phi2_cuts_can", "L0bar_p_L0bar_p_phi1_vs_phi2_cuts_can", 1200, 1000);

  TCanvas *L0bar_pi_L0bar_pi_phi1_vs_phi2_cuts_can = new TCanvas("L0bar_pi_L0bar_pi_phi1_vs_phi2_cuts_can", "L0bar_pi_L0bar_pi_phi1_vs_phi2_cuts_can", 1200, 1000);
  //________________________________________________________________________________________________________

  TF1 * eta_one_line = new TF1("eta_one_line", "[0]", -1, 1);
  eta_one_line->SetParameter(0, 1);
  eta_one_line->SetLineColor(1);

  TF1 * eta_minus_one_line = new TF1("eta_minus_one_line", "[0]", -1, 1);
  eta_minus_one_line->SetParameter(0, -1);
  eta_minus_one_line->SetLineColor(1);



  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //___________________________________________________________________________________________________

  L0_L0bar_cosThetaProdPlane_can->cd();

  gPad->SetLeftMargin(0.1);

  L0_L0bar_cosThetaProdPlane->GetXaxis()->SetTitle("cos(#theta*_{12})");
  L0_L0bar_cosThetaProdPlane->GetXaxis()->CenterTitle();
  L0_L0bar_cosThetaProdPlane->GetXaxis()->SetTitleSize(0.04);
  L0_L0bar_cosThetaProdPlane->GetXaxis()->SetLabelSize(0.04);
  L0_L0bar_cosThetaProdPlane->GetYaxis()->SetTitle("1/#it{N}_{evt} d#it{N}/d cos(#theta*_{12})");
  L0_L0bar_cosThetaProdPlane->GetYaxis()->CenterTitle();
  L0_L0bar_cosThetaProdPlane->GetYaxis()->SetTitleSize(0.04);
  L0_L0bar_cosThetaProdPlane->GetYaxis()->SetLabelSize(0.04);
  L0_L0bar_cosThetaProdPlane->GetYaxis()->SetMaxDigits(3);
  L0_L0bar_cosThetaProdPlane->SetMarkerStyle(20);
  L0_L0bar_cosThetaProdPlane->SetMarkerSize(2);
  L0_L0bar_cosThetaProdPlane->SetMarkerColor(kRed);
  L0_L0bar_cosThetaProdPlane->SetLineColor(kRed);
  L0_L0bar_cosThetaProdPlane->SetLineWidth(2);
  L0_L0bar_cosThetaProdPlane->Scale(1./L0_L0bar_cosThetaProdPlane->GetXaxis()->GetBinWidth(1));
  L0_L0bar_cosThetaProdPlane->Scale(1./nEvents);
  L0_L0bar_cosThetaProdPlane->SetMinimum(0);
  //L0_L0bar_cosThetaProdPlane->SetMinimum(1.6e6); //for detail
  L0_L0bar_cosThetaProdPlane->Draw("p e");

  TF1 *fitL0_L0bar_pT_US_ThetaStar = new TF1("fitL0_L0bar_pT_US_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
  fitL0_L0bar_pT_US_ThetaStar->SetParameters(100, 0.5);

  //fit_res_gaus_wrong_sign = L_inv_mass_pT_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
  L0_L0bar_cosThetaProdPlane->Fit(fitL0_L0bar_pT_US_ThetaStar, "s i 0 r");

  float P_L0_L0bar_pT = fitL0_L0bar_pT_US_ThetaStar->GetParameter(1)/(L0_alpha*L0bar_alpha);
  float P_L0_L0bar_pT_err = fitL0_L0bar_pT_US_ThetaStar->GetParError(1)/(L0_alpha*L0bar_alpha);

  fitL0_L0bar_pT_US_ThetaStar->SetLineColor(1);
  fitL0_L0bar_pT_US_ThetaStar->Draw("same");

  TPaveText *L0_L0bar_text_pt = new TPaveText(0.5, 0.2, 0.8, 0.6, "NDC");
  L0_L0bar_text_pt->SetTextFont(42);
  L0_L0bar_text_pt->AddText("PYTHIA 8.3");
  L0_L0bar_text_pt->AddText(Form("p+p #sqrt{s} = %i GeV", energy));
  L0_L0bar_text_pt->AddText("Minimum bias");
  L0_L0bar_text_pt->AddText("#Lambda-#bar{#Lambda}");
  //L0_L0bar_text_pt->AddText("No cuts");
  L0_L0bar_text_pt->AddText("|#it{y}_{#Lambda}| < 1");
  L0_L0bar_text_pt->AddText(Form("P = %.3f #pm %.3f", P_L0_L0bar_pT, P_L0_L0bar_pT_err));
  //L0bar_text_pt->AddText(pT_range->Data());
  L0_L0bar_text_pt->SetFillColorAlpha(0, 0.01);
  L0_L0bar_text_pt->Draw("same");

  L0_L0bar_cosThetaProdPlane_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations/L0_L0bar_cosThetaProdPlane.png");


  L0_L0_cosThetaProdPlane_can->cd();

  gPad->SetLeftMargin(0.1);

  L0_L0_cosThetaProdPlane->GetXaxis()->SetTitle("cos(#theta*_{12})");
  L0_L0_cosThetaProdPlane->GetXaxis()->CenterTitle();
  L0_L0_cosThetaProdPlane->GetXaxis()->SetTitleSize(0.04);
  L0_L0_cosThetaProdPlane->GetXaxis()->SetLabelSize(0.04);
  L0_L0_cosThetaProdPlane->GetYaxis()->SetTitle("1/#it{N}_{evt} d#it{N}/d cos(#theta*_{12})");
  L0_L0_cosThetaProdPlane->GetYaxis()->CenterTitle();
  L0_L0_cosThetaProdPlane->GetYaxis()->SetTitleOffset(1.25);
  L0_L0_cosThetaProdPlane->GetYaxis()->SetTitleSize(0.04);
  L0_L0_cosThetaProdPlane->GetYaxis()->SetLabelSize(0.04);
  L0_L0_cosThetaProdPlane->GetYaxis()->SetMaxDigits(1);
  L0_L0_cosThetaProdPlane->SetMarkerStyle(20);
  L0_L0_cosThetaProdPlane->SetMarkerSize(2);
  L0_L0_cosThetaProdPlane->SetMarkerColor(kRed);
  L0_L0_cosThetaProdPlane->SetLineColor(kRed);
  L0_L0_cosThetaProdPlane->SetLineWidth(2);
  L0_L0_cosThetaProdPlane->Scale(1./L0_L0_cosThetaProdPlane->GetXaxis()->GetBinWidth(1));
  L0_L0_cosThetaProdPlane->Scale(1./nEvents);
  L0_L0_cosThetaProdPlane->SetMinimum(0);
  //L0_L0_cosThetaProdPlane->SetMinimum(0.2e6);
  L0_L0_cosThetaProdPlane->Draw("p e");

  TF1 *fitL0_L0_pT_US_ThetaStar = new TF1("fitL0_L0_pT_US_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
  fitL0_L0_pT_US_ThetaStar->SetParameters(100, 0.5);

  //fit_res_gaus_wrong_sign = L_inv_mass_pT_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
  L0_L0_cosThetaProdPlane->Fit(fitL0_L0_pT_US_ThetaStar, "s i 0 r");

  float P_L0_L0_pT = fitL0_L0_pT_US_ThetaStar->GetParameter(1)/(L0_alpha*L0_alpha);
  float P_L0_L0_pT_err = fitL0_L0_pT_US_ThetaStar->GetParError(1)/(L0_alpha*L0_alpha);

  fitL0_L0_pT_US_ThetaStar->SetLineColor(1);
  fitL0_L0_pT_US_ThetaStar->Draw("same");

  TPaveText *L0_L0_text_pt = new TPaveText(0.5, 0.2, 0.8, 0.6, "NDC");
  L0_L0_text_pt->SetTextFont(42);
  L0_L0_text_pt->AddText("PYTHIA 8.3");
  L0_L0_text_pt->AddText(Form("p+p #sqrt{s} = %i GeV", energy));
  L0_L0_text_pt->AddText("Minimum bias");
  L0_L0_text_pt->AddText("#Lambda-#Lambda");
  //L0_L0_text_pt->AddText("No cuts");
  L0_L0_text_pt->AddText("|#it{y}_{#Lambda}| < 1");
  L0_L0_text_pt->AddText(Form("P = %.3f #pm %.3f", P_L0_L0_pT, P_L0_L0_pT_err));
  //L0bar_text_pt->AddText(pT_range->Data());
  L0_L0_text_pt->SetFillColorAlpha(0, 0.01);
  L0_L0_text_pt->Draw("same");

  L0_L0_cosThetaProdPlane_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations/L0_L0_cosThetaProdPlane.png");


  L0bar_L0bar_cosThetaProdPlane_can->cd();

  gPad->SetLeftMargin(0.1);

  L0bar_L0bar_cosThetaProdPlane->GetXaxis()->SetTitle("cos(#theta*_{12})");
  L0bar_L0bar_cosThetaProdPlane->GetXaxis()->CenterTitle();
  L0bar_L0bar_cosThetaProdPlane->GetXaxis()->SetTitleSize(0.04);
  L0bar_L0bar_cosThetaProdPlane->GetXaxis()->SetLabelSize(0.04);
  L0bar_L0bar_cosThetaProdPlane->GetYaxis()->SetTitle("1/#it{N}_{evt} d#it{N}/d cos(#theta*_{12})");
  L0bar_L0bar_cosThetaProdPlane->GetYaxis()->CenterTitle();
  L0bar_L0bar_cosThetaProdPlane->GetYaxis()->SetTitleOffset(1.25);
  L0bar_L0bar_cosThetaProdPlane->GetYaxis()->SetTitleSize(0.04);
  L0bar_L0bar_cosThetaProdPlane->GetYaxis()->SetLabelSize(0.04);
  L0bar_L0bar_cosThetaProdPlane->GetYaxis()->SetMaxDigits(1);
  L0bar_L0bar_cosThetaProdPlane->SetMarkerStyle(20);
  L0bar_L0bar_cosThetaProdPlane->SetMarkerSize(2);
  L0bar_L0bar_cosThetaProdPlane->SetMarkerColor(kRed);
  L0bar_L0bar_cosThetaProdPlane->SetLineColor(kRed);
  L0bar_L0bar_cosThetaProdPlane->SetLineWidth(2);
  L0bar_L0bar_cosThetaProdPlane->Scale(1./L0bar_L0bar_cosThetaProdPlane->GetXaxis()->GetBinWidth(1));
  L0bar_L0bar_cosThetaProdPlane->Scale(1./nEvents);
  L0bar_L0bar_cosThetaProdPlane->SetMinimum(0);
  //L0bar_L0bar_cosThetaProdPlane->SetMinimum(0.2e6);
  L0bar_L0bar_cosThetaProdPlane->Draw("p e");

  TF1 *fitL0bar_L0bar_pT_US_ThetaStar = new TF1("fitL0bar_L0bar_pT_US_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
  fitL0bar_L0bar_pT_US_ThetaStar->SetParameters(100, 0.5);

  //fit_res_gaus_wrong_sign = L_inv_mass_pT_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
  L0bar_L0bar_cosThetaProdPlane->Fit(fitL0bar_L0bar_pT_US_ThetaStar, "s i 0 r");

  float P_L0bar_L0bar_pT = fitL0bar_L0bar_pT_US_ThetaStar->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
  float P_L0bar_L0bar_pT_err = fitL0bar_L0bar_pT_US_ThetaStar->GetParError(1)/(L0bar_alpha*L0bar_alpha);

  fitL0bar_L0bar_pT_US_ThetaStar->SetLineColor(1);
  fitL0bar_L0bar_pT_US_ThetaStar->Draw("same");


  TPaveText *L0bar_L0bar_text_pt = new TPaveText(0.5, 0.2, 0.8, 0.6, "NDC");
  L0bar_L0bar_text_pt->SetTextFont(42);
  L0bar_L0bar_text_pt->AddText("PYTHIA 8.3");
  L0bar_L0bar_text_pt->AddText(Form("p+p #sqrt{s} = %i GeV", energy));
  L0bar_L0bar_text_pt->AddText("Minimum bias");
  L0bar_L0bar_text_pt->AddText("#bar{#Lambda}-#bar{#Lambda}");
  //L0bar_L0bar_text_pt->AddText("No cuts");
  L0bar_L0bar_text_pt->AddText("|#it{y}_{#Lambda}| < 1");
  L0bar_L0bar_text_pt->AddText(Form("P = %.3f #pm %.3f", P_L0bar_L0bar_pT, P_L0bar_L0bar_pT_err));
  //L0bar_text_pt->AddText(pT_range->Data());
  L0bar_L0bar_text_pt->SetFillColorAlpha(0, 0.01);
  L0bar_L0bar_text_pt->Draw("same");

  L0bar_L0bar_cosThetaProdPlane_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations/L0bar_L0bar_cosThetaProdPlane.png");

  //_____________________________________________________________________

  L0_L0bar_cosThetaProdPlane_cuts_can->cd();

  L0_L0bar_cosThetaProdPlane_cuts->GetXaxis()->SetTitle("cos(#theta*)");
  L0_L0bar_cosThetaProdPlane_cuts->GetXaxis()->CenterTitle();
  L0_L0bar_cosThetaProdPlane_cuts->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0_L0bar_cosThetaProdPlane_cuts->GetYaxis()->CenterTitle();
  L0_L0bar_cosThetaProdPlane_cuts->SetMarkerStyle(20);
  L0_L0bar_cosThetaProdPlane_cuts->SetMarkerColor(kRed);
  L0_L0bar_cosThetaProdPlane_cuts->SetLineColor(kRed);
  L0_L0bar_cosThetaProdPlane_cuts->Scale(1./L0_L0bar_cosThetaProdPlane->GetXaxis()->GetBinWidth(1));
  L0_L0bar_cosThetaProdPlane_cuts->SetMinimum(0);
  L0_L0bar_cosThetaProdPlane_cuts->Draw("p e");


  TPaveText *L0_L0bar_text_pt_cuts = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
  L0_L0bar_text_pt_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  L0_L0bar_text_pt_cuts->AddText("Minimum bias");
  L0_L0bar_text_pt_cuts->AddText("#Lambda-#bar{#Lambda}");
  L0_L0bar_text_pt_cuts->AddText("Analysis cuts");
  //L0bar_text_pt_cuts->AddText(pT_range->Data());
  L0_L0bar_text_pt_cuts->SetFillColorAlpha(0, 0.01);
  L0_L0bar_text_pt_cuts->Draw("same");

  L0_L0bar_cosThetaProdPlane_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations/L0_L0bar_cosThetaProdPlane_cuts.png");


  L0_L0_cosThetaProdPlane_cuts_can->cd();

  L0_L0_cosThetaProdPlane_cuts->GetXaxis()->SetTitle("cos(#theta*)");
  L0_L0_cosThetaProdPlane_cuts->GetXaxis()->CenterTitle();
  L0_L0_cosThetaProdPlane_cuts->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0_L0_cosThetaProdPlane_cuts->GetYaxis()->CenterTitle();
  L0_L0_cosThetaProdPlane_cuts->SetMarkerStyle(20);
  L0_L0_cosThetaProdPlane_cuts->SetMarkerColor(kRed);
  L0_L0_cosThetaProdPlane_cuts->SetLineColor(kRed);
  L0_L0_cosThetaProdPlane_cuts->Scale(1./L0_L0_cosThetaProdPlane->GetXaxis()->GetBinWidth(1));
  L0_L0_cosThetaProdPlane_cuts->SetMinimum(0);
  L0_L0_cosThetaProdPlane_cuts->Draw("p e");


  TPaveText *L0_L0_text_pt_cuts = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
  L0_L0_text_pt_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  L0_L0_text_pt_cuts->AddText("Minimum bias");
  L0_L0_text_pt_cuts->AddText("#Lambda-#Lambda");
  L0_L0_text_pt_cuts->AddText("Analysis cuts");
  //L0bar_text_pt_cuts->AddText(pT_range->Data());
  L0_L0_text_pt_cuts->SetFillColorAlpha(0, 0.01);
  L0_L0_text_pt_cuts->Draw("same");

  L0_L0_cosThetaProdPlane_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations/L0_L0_cosThetaProdPlane_cuts.png");


  L0bar_L0bar_cosThetaProdPlane_cuts_can->cd();

  L0bar_L0bar_cosThetaProdPlane_cuts->GetXaxis()->SetTitle("cos(#theta*)");
  L0bar_L0bar_cosThetaProdPlane_cuts->GetXaxis()->CenterTitle();
  L0bar_L0bar_cosThetaProdPlane_cuts->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0bar_L0bar_cosThetaProdPlane_cuts->GetYaxis()->CenterTitle();
  L0bar_L0bar_cosThetaProdPlane_cuts->SetMarkerStyle(20);
  L0bar_L0bar_cosThetaProdPlane_cuts->SetMarkerColor(kRed);
  L0bar_L0bar_cosThetaProdPlane_cuts->SetLineColor(kRed);
  L0bar_L0bar_cosThetaProdPlane_cuts->Scale(1./L0bar_L0bar_cosThetaProdPlane->GetXaxis()->GetBinWidth(1));
  L0bar_L0bar_cosThetaProdPlane_cuts->SetMinimum(0);
  L0bar_L0bar_cosThetaProdPlane_cuts->Draw("p e");

  TPaveText *L0bar_L0bar_text_pt_cuts = new TPaveText(0.6, 0.2, 0.8, 0.4, "NDC");
  L0bar_L0bar_text_pt_cuts->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
  L0bar_L0bar_text_pt_cuts->AddText("Minimum bias");
  L0bar_L0bar_text_pt_cuts->AddText("#bar{#Lambda}-#bar{#Lambda}");
  L0bar_L0bar_text_pt_cuts->AddText("Analysis cuts");
  //L0bar_text_pt_cuts->AddText(pT_range->Data());
  L0bar_L0bar_text_pt_cuts->SetFillColorAlpha(0, 0.01);
  L0bar_L0bar_text_pt_cuts->Draw("same");

  L0bar_L0bar_cosThetaProdPlane_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations/L0bar_L0bar_cosThetaProdPlane_cuts.png");


  //L-L correlations QA histograms

  //L-Lbar
  L0_L0bar_delta_eta_vs_delta_phi_can->cd();

  L0_L0bar_delta_eta_vs_delta_phi_hist->GetXaxis()->SetTitle("#Delta#eta");
  L0_L0bar_delta_eta_vs_delta_phi_hist->GetXaxis()->CenterTitle();
  L0_L0bar_delta_eta_vs_delta_phi_hist->GetYaxis()->SetTitle("#Delta#phi");
  L0_L0bar_delta_eta_vs_delta_phi_hist->GetYaxis()->CenterTitle();
  L0_L0bar_delta_eta_vs_delta_phi_hist->Draw("colz");

  L0_L0bar_delta_eta_vs_delta_phi_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0_L0bar_delta_eta_vs_delta_phi.png");


  L0_L0bar_y1_vs_y2_can->cd();

  L0_L0bar_y1_vs_y2_hist->GetXaxis()->SetTitle("y_{1}");
  L0_L0bar_y1_vs_y2_hist->GetXaxis()->CenterTitle();
  L0_L0bar_y1_vs_y2_hist->GetYaxis()->SetTitle("y_{2}");
  L0_L0bar_y1_vs_y2_hist->GetYaxis()->CenterTitle();
  L0_L0bar_y1_vs_y2_hist->Draw("colz");

  L0_L0bar_y1_vs_y2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0_L0bar_y1_vs_y2.png");


  L0_L0bar_pT1_vs_pT2_can->cd();

  L0_L0bar_pT1_vs_pT2_hist->GetXaxis()->SetTitle("#it{p}_{T}^{1} (GeV)/#it{c}");
  L0_L0bar_pT1_vs_pT2_hist->GetXaxis()->CenterTitle();
  L0_L0bar_pT1_vs_pT2_hist->GetYaxis()->SetTitle("#it{p}_{T}^{2} (GeV)/#it{c}");
  L0_L0bar_pT1_vs_pT2_hist->GetYaxis()->CenterTitle();
  L0_L0bar_pT1_vs_pT2_hist->Draw("colz");

  L0_L0bar_pT1_vs_pT2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0_L0bar_pT1_vs_pT2.png");


  L0_L0bar_phi1_vs_phi2_can->cd();

  L0_L0bar_phi1_vs_phi2_hist->GetXaxis()->SetTitle("#phi_{1}");
  L0_L0bar_phi1_vs_phi2_hist->GetXaxis()->CenterTitle();
  L0_L0bar_phi1_vs_phi2_hist->GetYaxis()->SetTitle("#phi_{2}");
  L0_L0bar_phi1_vs_phi2_hist->GetYaxis()->CenterTitle();
  L0_L0bar_phi1_vs_phi2_hist->Draw("colz");

  L0_L0bar_phi1_vs_phi2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0_L0bar_phi1_vs_phi2.png");



  //L-L
  L0_L0_delta_eta_vs_delta_phi_can->cd();

  L0_L0_delta_eta_vs_delta_phi_hist->GetXaxis()->SetTitle("#Delta#eta");
  L0_L0_delta_eta_vs_delta_phi_hist->GetXaxis()->CenterTitle();
  L0_L0_delta_eta_vs_delta_phi_hist->GetYaxis()->SetTitle("#Delta#phi");
  L0_L0_delta_eta_vs_delta_phi_hist->GetYaxis()->CenterTitle();
  L0_L0_delta_eta_vs_delta_phi_hist->Draw("colz");

  L0_L0_delta_eta_vs_delta_phi_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0_L0_delta_eta_vs_delta_phi.png");


  L0_L0_y1_vs_y2_can->cd();

  L0_L0_y1_vs_y2_hist->GetXaxis()->SetTitle("y_{1}");
  L0_L0_y1_vs_y2_hist->GetXaxis()->CenterTitle();
  L0_L0_y1_vs_y2_hist->GetYaxis()->SetTitle("y_{2}");
  L0_L0_y1_vs_y2_hist->GetYaxis()->CenterTitle();
  L0_L0_y1_vs_y2_hist->Draw("colz");

  L0_L0_y1_vs_y2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0_L0_y1_vs_y2.png");


  L0_L0_pT1_vs_pT2_can->cd();

  L0_L0_pT1_vs_pT2_hist->GetXaxis()->SetTitle("#it{p}_{T}^{1} (GeV)/#it{c}");
  L0_L0_pT1_vs_pT2_hist->GetXaxis()->CenterTitle();
  L0_L0_pT1_vs_pT2_hist->GetYaxis()->SetTitle("#it{p}_{T}^{2} (GeV)/#it{c}");
  L0_L0_pT1_vs_pT2_hist->GetYaxis()->CenterTitle();
  L0_L0_pT1_vs_pT2_hist->Draw("colz");

  L0_L0_pT1_vs_pT2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0_L0_pT1_vs_pT2.png");


  L0_L0_phi1_vs_phi2_can->cd();

  L0_L0_phi1_vs_phi2_hist->GetXaxis()->SetTitle("#phi_{1}");
  L0_L0_phi1_vs_phi2_hist->GetXaxis()->CenterTitle();
  L0_L0_phi1_vs_phi2_hist->GetYaxis()->SetTitle("#phi_{2}");
  L0_L0_phi1_vs_phi2_hist->GetYaxis()->CenterTitle();
  L0_L0_phi1_vs_phi2_hist->Draw("colz");

  L0_L0_phi1_vs_phi2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0_L0_phi1_vs_phi2.png");


  L0_pi_L0_pi_y1_vs_y2_can->cd();

  L0_pi_L0_pi_y1_vs_y2_hist->GetXaxis()->SetTitle("y_{1}^{#pi}");
  L0_pi_L0_pi_y1_vs_y2_hist->GetXaxis()->CenterTitle();
  L0_pi_L0_pi_y1_vs_y2_hist->GetYaxis()->SetTitle("y_{2}^{#pi}");
  L0_pi_L0_pi_y1_vs_y2_hist->GetYaxis()->CenterTitle();
  L0_pi_L0_pi_y1_vs_y2_hist->Draw("colz");

  L0_pi_L0_pi_y1_vs_y2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0_pi_L0_pi_y1_vs_y2.png");


  L0_p_L0_p_y1_vs_y2_can->cd();

  L0_p_L0_p_y1_vs_y2_hist->GetXaxis()->SetTitle("y_{1}^{p}");
  L0_p_L0_p_y1_vs_y2_hist->GetXaxis()->CenterTitle();
  L0_p_L0_p_y1_vs_y2_hist->GetYaxis()->SetTitle("y_{2}^{p}");
  L0_p_L0_p_y1_vs_y2_hist->GetYaxis()->CenterTitle();
  L0_p_L0_p_y1_vs_y2_hist->Draw("colz");

  L0_p_L0_p_y1_vs_y2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0_p_L0_p_y1_vs_y2.png");


  L0_pi_L0_pi_phi1_vs_phi2_can->cd();

  L0_pi_L0_pi_phi1_vs_phi2_hist->GetXaxis()->SetTitle("#phi_{1}^{#pi}");
  L0_pi_L0_pi_phi1_vs_phi2_hist->GetXaxis()->CenterTitle();
  L0_pi_L0_pi_phi1_vs_phi2_hist->GetYaxis()->SetTitle("#phi_{2}^{#pi}");
  L0_pi_L0_pi_phi1_vs_phi2_hist->GetYaxis()->CenterTitle();
  L0_pi_L0_pi_phi1_vs_phi2_hist->Draw("colz");

  L0_pi_L0_pi_phi1_vs_phi2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0_pi_L0_pi_phi1_vs_phi2.png");


  L0_p_L0_p_phi1_vs_phi2_can->cd();

  L0_p_L0_p_phi1_vs_phi2_hist->GetXaxis()->SetTitle("#phi_{1}^{p}");
  L0_p_L0_p_phi1_vs_phi2_hist->GetXaxis()->CenterTitle();
  L0_p_L0_p_phi1_vs_phi2_hist->GetYaxis()->SetTitle("#phi_{2}^{p}");
  L0_p_L0_p_phi1_vs_phi2_hist->GetYaxis()->CenterTitle();
  L0_p_L0_p_phi1_vs_phi2_hist->Draw("colz");

  L0_p_L0_p_phi1_vs_phi2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0_p_L0_p_phi1_vs_phi2.png");


  //Lbar-Lbar
  L0bar_L0bar_delta_eta_vs_delta_phi_can->cd();

  L0bar_L0bar_delta_eta_vs_delta_phi_hist->GetXaxis()->SetTitle("#Delta#eta");
  L0bar_L0bar_delta_eta_vs_delta_phi_hist->GetXaxis()->CenterTitle();
  L0bar_L0bar_delta_eta_vs_delta_phi_hist->GetYaxis()->SetTitle("#Delta#phi");
  L0bar_L0bar_delta_eta_vs_delta_phi_hist->GetYaxis()->CenterTitle();
  L0bar_L0bar_delta_eta_vs_delta_phi_hist->Draw("colz");

  L0bar_L0bar_delta_eta_vs_delta_phi_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0bar_L0bar_delta_eta_vs_delta_phi.png");


  L0bar_L0bar_y1_vs_y2_can->cd();

  L0bar_L0bar_y1_vs_y2_hist->GetXaxis()->SetTitle("y_{1}");
  L0bar_L0bar_y1_vs_y2_hist->GetXaxis()->CenterTitle();
  L0bar_L0bar_y1_vs_y2_hist->GetYaxis()->SetTitle("y_{2}");
  L0bar_L0bar_y1_vs_y2_hist->GetYaxis()->CenterTitle();
  L0bar_L0bar_y1_vs_y2_hist->Draw("colz");

  L0bar_L0bar_y1_vs_y2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0bar_L0bar_y1_vs_y2.png");


  L0bar_L0bar_pT1_vs_pT2_can->cd();

  L0bar_L0bar_pT1_vs_pT2_hist->GetXaxis()->SetTitle("#it{p}_{T}^{1} (GeV)/#it{c}");
  L0bar_L0bar_pT1_vs_pT2_hist->GetXaxis()->CenterTitle();
  L0bar_L0bar_pT1_vs_pT2_hist->GetYaxis()->SetTitle("#it{p}_{T}^{2} (GeV)/#it{c}");
  L0bar_L0bar_pT1_vs_pT2_hist->GetYaxis()->CenterTitle();
  L0bar_L0bar_pT1_vs_pT2_hist->Draw("colz");

  L0bar_L0bar_pT1_vs_pT2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0bar_L0bar_pT1_vs_pT2.png");


  L0bar_L0bar_phi1_vs_phi2_can->cd();

  L0bar_L0bar_phi1_vs_phi2_hist->GetXaxis()->SetTitle("#phi_{1}");
  L0bar_L0bar_phi1_vs_phi2_hist->GetXaxis()->CenterTitle();
  L0bar_L0bar_phi1_vs_phi2_hist->GetYaxis()->SetTitle("#phi_{2}");
  L0bar_L0bar_phi1_vs_phi2_hist->GetYaxis()->CenterTitle();
  L0bar_L0bar_phi1_vs_phi2_hist->Draw("colz");

  L0bar_L0bar_phi1_vs_phi2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0bar_L0bar_phi1_vs_phi2.png");


  L0bar_pi_L0bar_pi_y1_vs_y2_can->cd();

  L0bar_pi_L0bar_pi_y1_vs_y2_hist->GetXaxis()->SetTitle("y_{1}^{#pi}");
  L0bar_pi_L0bar_pi_y1_vs_y2_hist->GetXaxis()->CenterTitle();
  L0bar_pi_L0bar_pi_y1_vs_y2_hist->GetYaxis()->SetTitle("y_{2}^{#pi}");
  L0bar_pi_L0bar_pi_y1_vs_y2_hist->GetYaxis()->CenterTitle();
  L0bar_pi_L0bar_pi_y1_vs_y2_hist->Draw("colz");

  L0bar_pi_L0bar_pi_y1_vs_y2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0bar_pi_L0bar_pi_y1_vs_y2.png");


  L0bar_p_L0bar_p_y1_vs_y2_can->cd();

  L0bar_p_L0bar_p_y1_vs_y2_hist->GetXaxis()->SetTitle("y_{1}^{p}");
  L0bar_p_L0bar_p_y1_vs_y2_hist->GetXaxis()->CenterTitle();
  L0bar_p_L0bar_p_y1_vs_y2_hist->GetYaxis()->SetTitle("y_{2}^{p}");
  L0bar_p_L0bar_p_y1_vs_y2_hist->GetYaxis()->CenterTitle();
  L0bar_p_L0bar_p_y1_vs_y2_hist->Draw("colz");

  L0bar_p_L0bar_p_y1_vs_y2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0bar_p_L0bar_p_y1_vs_y2.png");


  L0bar_pi_L0bar_pi_phi1_vs_phi2_can->cd();

  L0bar_pi_L0bar_pi_phi1_vs_phi2_hist->GetXaxis()->SetTitle("#phi_{1}^{#pi}");
  L0bar_pi_L0bar_pi_phi1_vs_phi2_hist->GetXaxis()->CenterTitle();
  L0bar_pi_L0bar_pi_phi1_vs_phi2_hist->GetYaxis()->SetTitle("#phi_{2}^{#pi}");
  L0bar_pi_L0bar_pi_phi1_vs_phi2_hist->GetYaxis()->CenterTitle();
  L0bar_pi_L0bar_pi_phi1_vs_phi2_hist->Draw("colz");

  L0bar_pi_L0bar_pi_phi1_vs_phi2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0bar_pi_L0bar_pi_phi1_vs_phi2.png");


  L0bar_p_L0bar_p_phi1_vs_phi2_can->cd();

  L0bar_p_L0bar_p_phi1_vs_phi2_hist->GetXaxis()->SetTitle("#phi_{1}^{p}");
  L0bar_p_L0bar_p_phi1_vs_phi2_hist->GetXaxis()->CenterTitle();
  L0bar_p_L0bar_p_phi1_vs_phi2_hist->GetYaxis()->SetTitle("#phi_{2}^{p}");
  L0bar_p_L0bar_p_phi1_vs_phi2_hist->GetYaxis()->CenterTitle();
  L0bar_p_L0bar_p_phi1_vs_phi2_hist->Draw("colz");

  L0bar_p_L0bar_p_phi1_vs_phi2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0bar_p_L0bar_p_phi1_vs_phi2.png");

  //____________________________________________________________________________________________________________________________________________________________________

  //L-L correlations QA histograms after cuts

  //L-Lbar
  L0_L0bar_delta_eta_vs_delta_phi_cuts_can->cd();

  L0_L0bar_delta_eta_vs_delta_phi_cuts_hist->GetXaxis()->SetTitle("#Delta#eta");
  L0_L0bar_delta_eta_vs_delta_phi_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0bar_delta_eta_vs_delta_phi_cuts_hist->GetYaxis()->SetTitle("#Delta#phi");
  L0_L0bar_delta_eta_vs_delta_phi_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0bar_delta_eta_vs_delta_phi_cuts_hist->Draw("colz");

  L0_L0bar_delta_eta_vs_delta_phi_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0_L0bar_delta_eta_vs_delta_phi_cuts.png");


  L0_L0bar_y1_vs_y2_cuts_can->cd();

  L0_L0bar_y1_vs_y2_cuts_hist->GetXaxis()->SetTitle("y_{1}");
  L0_L0bar_y1_vs_y2_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0bar_y1_vs_y2_cuts_hist->GetYaxis()->SetTitle("y_{2}");
  L0_L0bar_y1_vs_y2_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0bar_y1_vs_y2_cuts_hist->Draw("colz");

  L0_L0bar_y1_vs_y2_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0_L0bar_y1_vs_y2_cuts.png");


  L0_L0bar_pT1_vs_pT2_cuts_can->cd();

  L0_L0bar_pT1_vs_pT2_cuts_hist->GetXaxis()->SetTitle("#it{p}_{T}^{1} (GeV)/#it{c}");
  L0_L0bar_pT1_vs_pT2_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0bar_pT1_vs_pT2_cuts_hist->GetYaxis()->SetTitle("#it{p}_{T}^{2} (GeV)/#it{c}");
  L0_L0bar_pT1_vs_pT2_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0bar_pT1_vs_pT2_cuts_hist->Draw("colz");

  L0_L0bar_pT1_vs_pT2_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0_L0bar_pT1_vs_pT2_cuts.png");


  L0_L0bar_phi1_vs_phi2_cuts_can->cd();

  L0_L0bar_phi1_vs_phi2_cuts_hist->GetXaxis()->SetTitle("#phi_{1}");
  L0_L0bar_phi1_vs_phi2_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0bar_phi1_vs_phi2_cuts_hist->GetYaxis()->SetTitle("#phi_{2}");
  L0_L0bar_phi1_vs_phi2_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0bar_phi1_vs_phi2_cuts_hist->Draw("colz");

  L0_L0bar_phi1_vs_phi2_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0_L0bar_phi1_vs_phi2_cuts.png");



  //L-L
  L0_L0_delta_eta_vs_delta_phi_cuts_can->cd();

  L0_L0_delta_eta_vs_delta_phi_cuts_hist->GetXaxis()->SetTitle("#Delta#eta");
  L0_L0_delta_eta_vs_delta_phi_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0_delta_eta_vs_delta_phi_cuts_hist->GetYaxis()->SetTitle("#Delta#phi");
  L0_L0_delta_eta_vs_delta_phi_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0_delta_eta_vs_delta_phi_cuts_hist->Draw("colz");

  L0_L0_delta_eta_vs_delta_phi_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0_L0_delta_eta_vs_delta_phi_cuts.png");


  L0_L0_y1_vs_y2_cuts_can->cd();

  L0_L0_y1_vs_y2_cuts_hist->GetXaxis()->SetTitle("y_{1}");
  L0_L0_y1_vs_y2_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0_y1_vs_y2_cuts_hist->GetYaxis()->SetTitle("y_{2}");
  L0_L0_y1_vs_y2_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0_y1_vs_y2_cuts_hist->Draw("colz");

  L0_L0_y1_vs_y2_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0_L0_y1_vs_y2_cuts.png");


  L0_L0_pT1_vs_pT2_cuts_can->cd();

  L0_L0_pT1_vs_pT2_cuts_hist->GetXaxis()->SetTitle("#it{p}_{T}^{1} (GeV)/#it{c}");
  L0_L0_pT1_vs_pT2_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0_pT1_vs_pT2_cuts_hist->GetYaxis()->SetTitle("#it{p}_{T}^{2} (GeV)/#it{c}");
  L0_L0_pT1_vs_pT2_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0_pT1_vs_pT2_cuts_hist->Draw("colz");

  L0_L0_pT1_vs_pT2_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0_L0_pT1_vs_pT2_cuts.png");


  L0_L0_phi1_vs_phi2_cuts_can->cd();

  L0_L0_phi1_vs_phi2_cuts_hist->GetXaxis()->SetTitle("#phi_{1}");
  L0_L0_phi1_vs_phi2_cuts_hist->GetXaxis()->CenterTitle();
  L0_L0_phi1_vs_phi2_cuts_hist->GetYaxis()->SetTitle("#phi_{2}");
  L0_L0_phi1_vs_phi2_cuts_hist->GetYaxis()->CenterTitle();
  L0_L0_phi1_vs_phi2_cuts_hist->Draw("colz");

  L0_L0_phi1_vs_phi2_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0_L0_phi1_vs_phi2_cuts.png");


  L0_pi_L0_pi_y1_vs_y2_cuts_can->cd();

  L0_pi_L0_pi_y1_vs_y2_cuts_hist->GetXaxis()->SetTitle("y_{1}^{#pi}");
  L0_pi_L0_pi_y1_vs_y2_cuts_hist->GetXaxis()->CenterTitle();
  L0_pi_L0_pi_y1_vs_y2_cuts_hist->GetYaxis()->SetTitle("y_{2}^{#pi}");
  L0_pi_L0_pi_y1_vs_y2_cuts_hist->GetYaxis()->CenterTitle();
  L0_pi_L0_pi_y1_vs_y2_cuts_hist->Draw("colz");

  L0_pi_L0_pi_y1_vs_y2_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0_pi_L0_pi_y1_vs_y2_cuts.png");


  L0_p_L0_p_y1_vs_y2_cuts_can->cd();

  L0_p_L0_p_y1_vs_y2_cuts_hist->GetXaxis()->SetTitle("y_{1}^{p}");
  L0_p_L0_p_y1_vs_y2_cuts_hist->GetXaxis()->CenterTitle();
  L0_p_L0_p_y1_vs_y2_cuts_hist->GetYaxis()->SetTitle("y_{2}^{p}");
  L0_p_L0_p_y1_vs_y2_cuts_hist->GetYaxis()->CenterTitle();
  L0_p_L0_p_y1_vs_y2_cuts_hist->Draw("colz");

  L0_p_L0_p_y1_vs_y2_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0_p_L0_p_y1_vs_y2_cuts.png");


  L0_pi_L0_pi_phi1_vs_phi2_cuts_can->cd();

  L0_pi_L0_pi_phi1_vs_phi2_cuts_hist->GetXaxis()->SetTitle("#phi_{1}^{#pi}");
  L0_pi_L0_pi_phi1_vs_phi2_cuts_hist->GetXaxis()->CenterTitle();
  L0_pi_L0_pi_phi1_vs_phi2_cuts_hist->GetYaxis()->SetTitle("#phi_{2}^{#pi}");
  L0_pi_L0_pi_phi1_vs_phi2_cuts_hist->GetYaxis()->CenterTitle();
  L0_pi_L0_pi_phi1_vs_phi2_cuts_hist->Draw("colz");

  L0_pi_L0_pi_phi1_vs_phi2_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0_pi_L0_pi_phi1_vs_phi2_cuts.png");


  L0_p_L0_p_phi1_vs_phi2_cuts_can->cd();

  L0_p_L0_p_phi1_vs_phi2_cuts_hist->GetXaxis()->SetTitle("#phi_{1}^{p}");
  L0_p_L0_p_phi1_vs_phi2_cuts_hist->GetXaxis()->CenterTitle();
  L0_p_L0_p_phi1_vs_phi2_cuts_hist->GetYaxis()->SetTitle("#phi_{2}^{p}");
  L0_p_L0_p_phi1_vs_phi2_cuts_hist->GetYaxis()->CenterTitle();
  L0_p_L0_p_phi1_vs_phi2_cuts_hist->Draw("colz");

  L0_p_L0_p_phi1_vs_phi2_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0_p_L0_p_phi1_vs_phi2_cuts.png");


  //Lbar-Lbar
  L0bar_L0bar_delta_eta_vs_delta_phi_cuts_can->cd();

  L0bar_L0bar_delta_eta_vs_delta_phi_cuts_hist->GetXaxis()->SetTitle("#Delta#eta");
  L0bar_L0bar_delta_eta_vs_delta_phi_cuts_hist->GetXaxis()->CenterTitle();
  L0bar_L0bar_delta_eta_vs_delta_phi_cuts_hist->GetYaxis()->SetTitle("#Delta#phi");
  L0bar_L0bar_delta_eta_vs_delta_phi_cuts_hist->GetYaxis()->CenterTitle();
  L0bar_L0bar_delta_eta_vs_delta_phi_cuts_hist->Draw("colz");

  L0bar_L0bar_delta_eta_vs_delta_phi_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0bar_L0bar_delta_eta_vs_delta_phi_cuts.png");


  L0bar_L0bar_y1_vs_y2_cuts_can->cd();

  L0bar_L0bar_y1_vs_y2_cuts_hist->GetXaxis()->SetTitle("y_{1}");
  L0bar_L0bar_y1_vs_y2_cuts_hist->GetXaxis()->CenterTitle();
  L0bar_L0bar_y1_vs_y2_cuts_hist->GetYaxis()->SetTitle("y_{2}");
  L0bar_L0bar_y1_vs_y2_cuts_hist->GetYaxis()->CenterTitle();
  L0bar_L0bar_y1_vs_y2_cuts_hist->Draw("colz");

  L0bar_L0bar_y1_vs_y2_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0bar_L0bar_y1_vs_y2_cuts.png");


  L0bar_L0bar_pT1_vs_pT2_cuts_can->cd();

  L0bar_L0bar_pT1_vs_pT2_cuts_hist->GetXaxis()->SetTitle("#it{p}_{T}^{1} (GeV)/#it{c}");
  L0bar_L0bar_pT1_vs_pT2_cuts_hist->GetXaxis()->CenterTitle();
  L0bar_L0bar_pT1_vs_pT2_cuts_hist->GetYaxis()->SetTitle("#it{p}_{T}^{2} (GeV)/#it{c}");
  L0bar_L0bar_pT1_vs_pT2_cuts_hist->GetYaxis()->CenterTitle();
  L0bar_L0bar_pT1_vs_pT2_cuts_hist->Draw("colz");

  L0bar_L0bar_pT1_vs_pT2_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0bar_L0bar_pT1_vs_pT2_cuts.png");


  L0bar_L0bar_phi1_vs_phi2_cuts_can->cd();

  L0bar_L0bar_phi1_vs_phi2_cuts_hist->GetXaxis()->SetTitle("#phi_{1}");
  L0bar_L0bar_phi1_vs_phi2_cuts_hist->GetXaxis()->CenterTitle();
  L0bar_L0bar_phi1_vs_phi2_cuts_hist->GetYaxis()->SetTitle("#phi_{2}");
  L0bar_L0bar_phi1_vs_phi2_cuts_hist->GetYaxis()->CenterTitle();
  L0bar_L0bar_phi1_vs_phi2_cuts_hist->Draw("colz");

  L0bar_L0bar_phi1_vs_phi2_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0bar_L0bar_phi1_vs_phi2_cuts.png");


  L0bar_pi_L0bar_pi_y1_vs_y2_cuts_can->cd();

  L0bar_pi_L0bar_pi_y1_vs_y2_cuts_hist->GetXaxis()->SetTitle("y_{1}^{#pi}");
  L0bar_pi_L0bar_pi_y1_vs_y2_cuts_hist->GetXaxis()->CenterTitle();
  L0bar_pi_L0bar_pi_y1_vs_y2_cuts_hist->GetYaxis()->SetTitle("y_{2}^{#pi}");
  L0bar_pi_L0bar_pi_y1_vs_y2_cuts_hist->GetYaxis()->CenterTitle();
  L0bar_pi_L0bar_pi_y1_vs_y2_cuts_hist->Draw("colz");

  L0bar_pi_L0bar_pi_y1_vs_y2_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0bar_pi_L0bar_pi_y1_vs_y2_cuts.png");


  L0bar_p_L0bar_p_y1_vs_y2_cuts_can->cd();

  L0bar_p_L0bar_p_y1_vs_y2_cuts_hist->GetXaxis()->SetTitle("y_{1}^{p}");
  L0bar_p_L0bar_p_y1_vs_y2_cuts_hist->GetXaxis()->CenterTitle();
  L0bar_p_L0bar_p_y1_vs_y2_cuts_hist->GetYaxis()->SetTitle("y_{2}^{p}");
  L0bar_p_L0bar_p_y1_vs_y2_cuts_hist->GetYaxis()->CenterTitle();
  L0bar_p_L0bar_p_y1_vs_y2_cuts_hist->Draw("colz");

  L0bar_p_L0bar_p_y1_vs_y2_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0bar_p_L0bar_p_y1_vs_y2_cuts.png");


  L0bar_pi_L0bar_pi_phi1_vs_phi2_cuts_can->cd();

  L0bar_pi_L0bar_pi_phi1_vs_phi2_cuts_hist->GetXaxis()->SetTitle("#phi_{1}^{#pi}");
  L0bar_pi_L0bar_pi_phi1_vs_phi2_cuts_hist->GetXaxis()->CenterTitle();
  L0bar_pi_L0bar_pi_phi1_vs_phi2_cuts_hist->GetYaxis()->SetTitle("#phi_{2}^{#pi}");
  L0bar_pi_L0bar_pi_phi1_vs_phi2_cuts_hist->GetYaxis()->CenterTitle();
  L0bar_pi_L0bar_pi_phi1_vs_phi2_cuts_hist->Draw("colz");

  L0bar_pi_L0bar_pi_phi1_vs_phi2_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0bar_pi_L0bar_pi_phi1_vs_phi2_cuts.png");


  L0bar_p_L0bar_p_phi1_vs_phi2_cuts_can->cd();

  L0bar_p_L0bar_p_phi1_vs_phi2_cuts_hist->GetXaxis()->SetTitle("#phi_{1}^{p}");
  L0bar_p_L0bar_p_phi1_vs_phi2_cuts_hist->GetXaxis()->CenterTitle();
  L0bar_p_L0bar_p_phi1_vs_phi2_cuts_hist->GetYaxis()->SetTitle("#phi_{2}^{p}");
  L0bar_p_L0bar_p_phi1_vs_phi2_cuts_hist->GetYaxis()->CenterTitle();
  L0bar_p_L0bar_p_phi1_vs_phi2_cuts_hist->Draw("colz");

  L0bar_p_L0bar_p_phi1_vs_phi2_cuts_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/correlations_QA/L0bar_p_L0bar_p_phi1_vs_phi2_cuts.png");



  for(unsigned int pTbin = 0; pTbin < nPtBins+1; pTbin++)
  {

    TString *pT_range = new TString();
    if(pTbin < nPtBins) pT_range->Form("%0.1f < #it{p}_{T} < %0.1f GeV/#it{c}", pT_bins[pTbin], pT_bins[pTbin+1]);
    else pT_range->Form("p_{T} integrated");

    TPaveText *L0_text_y = new TPaveText(0.2, 0.7, 0.4, 0.85, "NDC");
    L0_text_y->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0_text_y->AddText("Minimum bias");
    L0_text_y->AddText("#Lambda^{0}");
    L0_text_y->AddText(pT_range->Data());
    L0_text_y->SetFillColorAlpha(0, 0.01);


    L0_y_vs_p_eta_can[pTbin] = new TCanvas(Form("L0_y_vs_p_eta_can_pT_%i", pTbin), Form("L0_y_vs_p_eta_can_pT_%i", pTbin), 1200, 1000);
    L0_y_vs_p_eta_can[pTbin]->cd();

    L0_y_vs_p_eta[pTbin] = (TH2D*)inFile->Get(Form("L0_y_vs_p_eta_pT_%i", pTbin));
    L0_y_vs_p_eta[pTbin]->GetYaxis()->SetTitle("#eta_{p}");
    L0_y_vs_p_eta[pTbin]->GetYaxis()->CenterTitle();
    L0_y_vs_p_eta[pTbin]->GetXaxis()->SetTitle("#it{y}_{#Lambda^{0}}");
    L0_y_vs_p_eta[pTbin]->GetXaxis()->CenterTitle();
    L0_y_vs_p_eta[pTbin]->Draw("colz");

    eta_one_line->Draw("same");
    eta_minus_one_line->Draw("same");

    L0_text_y->Draw("same");

    L0_y_vs_p_eta_can[pTbin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L0_y_vs_p_eta_pT_%i.png", pTbin));


    L0_y_vs_pi_eta_can[pTbin] = new TCanvas(Form("L0_y_vs_pi_eta_can_pT_%i", pTbin), Form("L0_y_vs_pi_eta_can_pT_%i", pTbin), 1200, 1000);
    L0_y_vs_pi_eta_can[pTbin]->cd();

    L0_y_vs_pi_eta[pTbin] = (TH2D*)inFile->Get(Form("L0_y_vs_pi_eta_pT_%i", pTbin));
    L0_y_vs_pi_eta[pTbin]->GetYaxis()->SetTitle("#eta_{#pi}");
    L0_y_vs_pi_eta[pTbin]->GetYaxis()->CenterTitle();
    L0_y_vs_pi_eta[pTbin]->GetXaxis()->SetTitle("#it{y}_{#Lambda^{0}}");
    L0_y_vs_pi_eta[pTbin]->GetXaxis()->CenterTitle();
    L0_y_vs_pi_eta[pTbin]->Draw("colz");

    L0_text_y->Draw("same");

    eta_one_line->Draw("same");
    eta_minus_one_line->Draw("same");

    L0_y_vs_pi_eta_can[pTbin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L0_y_vs_pi_eta_pT_%i.png", pTbin));



    TPaveText *L0bar_text_y = new TPaveText(0.2, 0.7, 0.4, 0.85, "NDC");
    L0bar_text_y->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
    L0bar_text_y->AddText("Minimum bias");
    L0bar_text_y->AddText("#bar{#Lambda^{0}}");
    L0bar_text_y->AddText(pT_range->Data());
    L0bar_text_y->SetFillColorAlpha(0, 0.01);

    L0bar_y_vs_p_eta_can[pTbin] = new TCanvas(Form("L0bar_y_vs_p_eta_can_pT_%i", pTbin), Form("L0bar_y_vs_p_eta_can_pT_%i", pTbin), 1200, 1000);
    L0bar_y_vs_p_eta_can[pTbin]->cd();

    L0bar_y_vs_p_eta[pTbin] = (TH2D*)inFile->Get(Form("L0bar_y_vs_p_eta_pT_%i", pTbin));
    L0bar_y_vs_p_eta[pTbin]->GetYaxis()->SetTitle("#eta_{p}");
    L0bar_y_vs_p_eta[pTbin]->GetYaxis()->CenterTitle();
    L0bar_y_vs_p_eta[pTbin]->GetXaxis()->SetTitle("#it{y}_{#bar{#Lambda^{0}}}");
    L0bar_y_vs_p_eta[pTbin]->GetXaxis()->CenterTitle();
    L0bar_y_vs_p_eta[pTbin]->Draw("colz");

    eta_one_line->Draw("same");
    eta_minus_one_line->Draw("same");

    L0bar_text_y->Draw("same");

    L0bar_y_vs_p_eta_can[pTbin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L0bar_y_vs_p_eta_pT_%i.png", pTbin));


    L0bar_y_vs_pi_eta_can[pTbin] = new TCanvas(Form("L0bar_y_vs_pi_eta_can_pT_%i", pTbin), Form("L0bar_y_vs_pi_eta_can_pT_%i", pTbin), 1200, 1000);

    L0bar_y_vs_pi_eta_can[pTbin]->cd();

    L0bar_y_vs_pi_eta[pTbin] = (TH2D*)inFile->Get(Form("L0bar_y_vs_pi_eta_pT_%i", pTbin));
    L0bar_y_vs_pi_eta[pTbin]->GetYaxis()->SetTitle("#eta_{#pi}");
    L0bar_y_vs_pi_eta[pTbin]->GetYaxis()->CenterTitle();
    L0bar_y_vs_pi_eta[pTbin]->GetXaxis()->SetTitle("#it{y}_{#bar{#Lambda^{0}}}");
    L0bar_y_vs_pi_eta[pTbin]->GetXaxis()->CenterTitle();
    L0bar_y_vs_pi_eta[pTbin]->Draw("colz");

    eta_one_line->Draw("same");
    eta_minus_one_line->Draw("same");

    L0bar_text_y->Draw("same");

    L0bar_y_vs_pi_eta_can[pTbin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/daughter_kine/L0bar_y_vs_pi_eta_pT_%i.png", pTbin));
    //______________________________________________________________________________________________________________________________________________________________________



    for(unsigned int etaBin = 0; etaBin < nEtaBins+1; etaBin++)
    {

      TString *eta_range = new TString();
      if(etaBin < nEtaBins) eta_range->Form("%0.1f < #eta < %0.1f", eta_bins[etaBin], eta_bins[etaBin+1]);
      else eta_range->Form("-1 < #eta < 1");


      //for histograms with bins in eta only
      //the outer loop is over pT
      if(pTbin == 0)
      {
        TPaveText *L0_text_pz = new TPaveText(0.2, 0.7, 0.4, 0.85, "NDC");
        L0_text_pz->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
        L0_text_pz->AddText("Minimum bias");
        L0_text_pz->AddText("#Lambda^{0}");
        L0_text_pz->AddText(eta_range->Data());
        L0_text_pz->SetFillColorAlpha(0, 0.01);


        L0_pT_vs_L0_pz_can[etaBin] = new TCanvas(Form("L0_pT_vs_L0_pz_can_eta_%i", etaBin), Form("L0_pT_vs_L0_pz_can_eta_%i", etaBin), 1200, 1000);
        L0_pT_vs_L0_pz_can[etaBin]->cd();

        L0_pT_vs_L0_pz[etaBin] = (TH2D*)inFile->Get(Form("L0_pT_vs_L0_pz_eta_%i", etaBin));
        L0_pT_vs_L0_pz[etaBin]->GetXaxis()->SetTitle("p_{T} (GeV)/#it{c}");
        L0_pT_vs_L0_pz[etaBin]->GetXaxis()->CenterTitle();
        L0_pT_vs_L0_pz[etaBin]->GetYaxis()->SetTitle("p_{z} (GeV)/#it{c}");
        L0_pT_vs_L0_pz[etaBin]->GetYaxis()->CenterTitle();
        L0_pT_vs_L0_pz[etaBin]->Draw("colz");

        L0_text_pz->Draw("same");

        L0_pT_vs_L0_pz_can[etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine/L0_pT_vs_L0_pz_eta_%i.png", etaBin));


        TPaveText *L0bar_text_pz = new TPaveText(0.2, 0.7, 0.4, 0.85, "NDC");
        L0bar_text_pz->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
        L0bar_text_pz->AddText("Minimum bias");
        L0bar_text_pz->AddText("#bar{#Lambda^{0}}");
        L0bar_text_pz->AddText(eta_range->Data());
        L0bar_text_pz->SetFillColorAlpha(0, 0.01);

        L0bar_pT_vs_L0bar_pz_can[etaBin] = new TCanvas(Form("L0bar_pT_vs_L0bar_pz_can_eta_%i", etaBin), Form("L0bar_pT_vs_L0bar_pz_can_eta_%i", etaBin), 1200, 1000);
        L0bar_pT_vs_L0bar_pz_can[etaBin]->cd();

        L0bar_pT_vs_L0bar_pz[etaBin] = (TH2D*)inFile->Get(Form("L0bar_pT_vs_L0bar_pz_eta_%i", etaBin));
        L0bar_pT_vs_L0bar_pz[etaBin]->GetXaxis()->SetTitle("p_{T} (GeV)/#it{c}");
        L0bar_pT_vs_L0bar_pz[etaBin]->GetXaxis()->CenterTitle();
        L0bar_pT_vs_L0bar_pz[etaBin]->GetYaxis()->SetTitle("p_{z} (GeV)/#it{c}");
        L0bar_pT_vs_L0bar_pz[etaBin]->GetYaxis()->CenterTitle();
        L0bar_pT_vs_L0bar_pz[etaBin]->Draw("colz");

        L0bar_text_pz->Draw("same");

        L0bar_pT_vs_L0bar_pz_can[etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine/L0bar_pT_vs_L0bar_pz_eta_%i.png", etaBin));
      }





      L0_thetaProdPlane_can[pTbin][etaBin] = new TCanvas(Form("L0_thetaProdPlane_can_pt_%i_eta_%i", pTbin, etaBin), Form("L0_thetaProdPlane_can_pt_%i_eta_%i", pTbin, etaBin), 1200, 1000);
      L0_cosThetaProdPlane_can[pTbin][etaBin] = new TCanvas(Form("L0_cosThetaProdPlane_can_pT_%i_eta_%i", pTbin, etaBin), Form("L0_cosThetaProdPlane_can_pT_%i_eta_%i", pTbin, etaBin), 1200, 1000);

      L0_pz_can[pTbin][etaBin] = new TCanvas(Form("L0_pz_can_pT_%i_eta_%i", pTbin, etaBin), Form("L0_pz_can_pT_%i_eta_%i", pTbin, etaBin), 1200, 1000);
      L0_xF_can[pTbin][etaBin] = new TCanvas(Form("L0_xF_can_pT_%i_eta_%i", pTbin, etaBin), Form("L0_xF_can_pT_%i_eta_%i", pTbin, etaBin), 1200, 1000);


      L0bar_thetaProdPlane_can[pTbin][etaBin] = new TCanvas(Form("L0bar_thetaProdPlane_can_pt_%i_eta_%i", pTbin, etaBin), Form("L0bar_thetaProdPlane_can_pt_%i_eta_%i", pTbin, etaBin), 1200, 1000);
      L0bar_cosThetaProdPlane_can[pTbin][etaBin] = new TCanvas(Form("L0bar_cosThetaProdPlane_can_pT_%i_eta_%i", pTbin, etaBin), Form("L0bar_cosThetaProdPlane_can_pT_%i_eta_%i", pTbin, etaBin), 1200, 1000);

      L0bar_pz_can[pTbin][etaBin] = new TCanvas(Form("L0bar_pz_can_pT_%i_eta_%i", pTbin, etaBin), Form("L0bar_pz_can_pT_%i_eta_%i", pTbin, etaBin), 1200, 1000);
      L0bar_xF_can[pTbin][etaBin] = new TCanvas(Form("L0bar_xF_can_pT_%i_eta_%i", pTbin, etaBin), Form("L0bar_xF_can_pT_%i_eta_%i", pTbin, etaBin), 1200, 1000);





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


      L0_thetaProdPlane_can[pTbin][etaBin]->cd();

      L0_thetaProdPlane_hist[pTbin][etaBin] = (TH1D*)inFile->Get(Form("L0_thetaProdPlane_pT_%i_eta_%i", pTbin, etaBin));
      L0_thetaProdPlane_hist[pTbin][etaBin]->SetMarkerStyle(20);
      L0_thetaProdPlane_hist[pTbin][etaBin]->SetMarkerColor(kRed);
      L0_thetaProdPlane_hist[pTbin][etaBin]->SetLineColor(kRed);
      L0_thetaProdPlane_hist[pTbin][etaBin]->GetXaxis()->SetTitle("#theta*");
      L0_thetaProdPlane_hist[pTbin][etaBin]->GetXaxis()->CenterTitle();
      L0_thetaProdPlane_hist[pTbin][etaBin]->GetYaxis()->SetTitle("d#{N}/d#theta*");
      L0_thetaProdPlane_hist[pTbin][etaBin]->GetYaxis()->CenterTitle();
      L0_thetaProdPlane_hist[pTbin][etaBin]->Sumw2();
      L0_thetaProdPlane_hist[pTbin][etaBin]->Scale(L0_thetaProdPlane_hist[pTbin][etaBin]->GetXaxis()->GetBinWidth(0)); //scale by bin width to obtain d N/d theta*
      L0_thetaProdPlane_hist[pTbin][etaBin]->SetMinimum(0);
      L0_thetaProdPlane_hist[pTbin][etaBin]->Draw("p e");

      TLegend *PYTHIA_legend = new TLegend(0.2, 0.4, 0.4, 0.5 );
      PYTHIA_legend->AddEntry(L0_thetaProdPlane_hist[pTbin][etaBin], "PYTHIA 8.162");
      PYTHIA_legend->SetBorderSize(0);
      PYTHIA_legend->SetFillColorAlpha(0, 0.01);
      PYTHIA_legend->Draw("same");

      cent_text_2->Draw("same");


      L0_cosThetaProdPlane_can[pTbin][etaBin]->cd();

      L0_cosThetaProdPlane_hist[pTbin][etaBin] = (TH1D*)inFile->Get(Form("L0_cosThetaProdPlane_pT_%i_eta_%i", pTbin, etaBin));
      L0_cosThetaProdPlane_hist[pTbin][etaBin]->SetMarkerStyle(20);
      L0_cosThetaProdPlane_hist[pTbin][etaBin]->SetMarkerColor(kRed);
      L0_cosThetaProdPlane_hist[pTbin][etaBin]->SetLineColor(kRed);
      L0_cosThetaProdPlane_hist[pTbin][etaBin]->GetXaxis()->SetTitle("cos(#theta*)");
      L0_cosThetaProdPlane_hist[pTbin][etaBin]->GetXaxis()->CenterTitle();
      L0_cosThetaProdPlane_hist[pTbin][etaBin]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      L0_cosThetaProdPlane_hist[pTbin][etaBin]->GetYaxis()->CenterTitle();
      L0_cosThetaProdPlane_hist[pTbin][etaBin]->Sumw2();
      L0_cosThetaProdPlane_hist[pTbin][etaBin]->Scale(L0_cosThetaProdPlane_hist[pTbin][etaBin]->GetXaxis()->GetBinWidth(0)); //scale by vin width to obtain d N/d cos(theta*)
      L0_cosThetaProdPlane_hist[pTbin][etaBin]->SetMinimum(0);
      L0_cosThetaProdPlane_hist[pTbin][etaBin]->Draw("p e");


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


      L0bar_thetaProdPlane_can[pTbin][etaBin]->cd();

      L0bar_thetaProdPlane_hist[pTbin][etaBin] = (TH1D*)inFile->Get(Form("L0bar_thetaProdPlane_pT_%i_eta_%i", pTbin, etaBin));
      L0bar_thetaProdPlane_hist[pTbin][etaBin]->SetMarkerStyle(20);
      L0bar_thetaProdPlane_hist[pTbin][etaBin]->SetMarkerColor(kRed);
      L0bar_thetaProdPlane_hist[pTbin][etaBin]->SetLineColor(kRed);
      L0bar_thetaProdPlane_hist[pTbin][etaBin]->GetXaxis()->SetTitle("#theta*");
      L0bar_thetaProdPlane_hist[pTbin][etaBin]->GetXaxis()->CenterTitle();
      L0bar_thetaProdPlane_hist[pTbin][etaBin]->GetYaxis()->SetTitle("d#{N}/d#theta*");
      L0bar_thetaProdPlane_hist[pTbin][etaBin]->GetYaxis()->CenterTitle();
      L0bar_thetaProdPlane_hist[pTbin][etaBin]->Sumw2();
      L0bar_thetaProdPlane_hist[pTbin][etaBin]->Scale(L0bar_thetaProdPlane_hist[pTbin][etaBin]->GetXaxis()->GetBinWidth(0)); //scale by bin width to obtain d N/d theta*
      L0bar_thetaProdPlane_hist[pTbin][etaBin]->SetMinimum(0);
      L0bar_thetaProdPlane_hist[pTbin][etaBin]->Draw("p e");

      PYTHIA_legend->Draw("same");

      cent_text_3->Draw("same");


      L0bar_cosThetaProdPlane_can[pTbin][etaBin]->cd();

      L0bar_cosThetaProdPlane_hist[pTbin][etaBin] = (TH1D*)inFile->Get(Form("L0bar_cosThetaProdPlane_pT_%i_eta_%i", pTbin, etaBin));
      L0bar_cosThetaProdPlane_hist[pTbin][etaBin]->SetMarkerStyle(20);
      L0bar_cosThetaProdPlane_hist[pTbin][etaBin]->SetMarkerColor(kRed);
      L0bar_cosThetaProdPlane_hist[pTbin][etaBin]->SetLineColor(kRed);
      L0bar_cosThetaProdPlane_hist[pTbin][etaBin]->GetXaxis()->SetTitle("cos(#theta*)");
      L0bar_cosThetaProdPlane_hist[pTbin][etaBin]->GetXaxis()->CenterTitle();
      L0bar_cosThetaProdPlane_hist[pTbin][etaBin]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      L0bar_cosThetaProdPlane_hist[pTbin][etaBin]->GetYaxis()->CenterTitle();
      L0bar_cosThetaProdPlane_hist[pTbin][etaBin]->Sumw2();
      L0bar_cosThetaProdPlane_hist[pTbin][etaBin]->Scale(L0bar_cosThetaProdPlane_hist[pTbin][etaBin]->GetXaxis()->GetBinWidth(0)); //scale by vin width to obtain d N/d cos(theta*)
      L0bar_cosThetaProdPlane_hist[pTbin][etaBin]->SetMinimum(0);
      L0bar_cosThetaProdPlane_hist[pTbin][etaBin]->Draw("p e");


      PYTHIA_legend->Draw("same");

      cent_text_3->Draw("same");


      //________________________________________________________________________________________________________________________________________________________________


      TLegend *PYTHIA_legend_2 = new TLegend(0.2, 0.8, 0.4, 0.89 );
      PYTHIA_legend_2->AddEntry(L0_thetaProdPlane_hist[pTbin][etaBin], "PYTHIA 8.162");
      PYTHIA_legend_2->SetBorderSize(0);
      PYTHIA_legend_2->SetFillColorAlpha(0, 0.01);
      PYTHIA_legend_2->Draw("same");

      TPaveText *cent_text_4 = new TPaveText(0.6, 0.7, 0.8, 0.85, "NDC");
      cent_text_4->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
      cent_text_4->AddText("Minimum bias");
      cent_text_4->AddText("#Lambda^{0}");
      cent_text_4->AddText(eta_range->Data());
      cent_text_4->AddText(pT_range->Data());
      cent_text_4->SetFillColorAlpha(0, 0.01);


      L0_pz_can[pTbin][etaBin]->cd();

      gPad->SetLogy();

      L0_pz[pTbin][etaBin] = (TH1D*)inFile->Get(Form("L0_pz_pT_%i_eta_%i", pTbin, etaBin));
      L0_pz[pTbin][etaBin]->GetXaxis()->SetTitle("p_{z} GeV/#it{c}");
      L0_pz[pTbin][etaBin]->GetXaxis()->CenterTitle();
      L0_pz[pTbin][etaBin]->GetXaxis()->SetRangeUser( -6, 6);
      L0_pz[pTbin][etaBin]->GetYaxis()->SetTitle("Counts");
      L0_pz[pTbin][etaBin]->GetYaxis()->CenterTitle();
      L0_pz[pTbin][etaBin]->SetMarkerStyle(20);
      L0_pz[pTbin][etaBin]->SetMarkerColor(kRed);
      L0_pz[pTbin][etaBin]->Draw("p e");

      PYTHIA_legend_2->Draw("same");

      cent_text_4->Draw("same");


      L0_xF_can[pTbin][etaBin]->cd();

      gPad->SetLogy();

      L0_xF[pTbin][etaBin] = (TH1D*)inFile->Get(Form("L0_xF_pT_%i_eta_%i", pTbin, etaBin));
      L0_xF[pTbin][etaBin]->GetXaxis()->SetTitle("x_{F}");
      L0_xF[pTbin][etaBin]->GetXaxis()->CenterTitle();
      L0_xF[pTbin][etaBin]->GetYaxis()->SetTitle("Counts");
      L0_xF[pTbin][etaBin]->GetYaxis()->CenterTitle();
      L0_xF[pTbin][etaBin]->SetMarkerStyle(20);
      L0_xF[pTbin][etaBin]->SetMarkerColor(kRed);
      L0_xF[pTbin][etaBin]->Draw("p e");

      PYTHIA_legend_2->Draw("same");

      cent_text_4->Draw("same");


      TPaveText *cent_text_5 = new TPaveText(0.6, 0.7, 0.8, 0.85, "NDC");
      //cent_text->AddText("STAR preliminary");
      //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
      cent_text_5->AddText(Form("MC p+p #sqrt{s} = %i GeV", energy));
      cent_text_5->AddText("Minimum bias");
      cent_text_5->AddText("#bar{#Lambda^{0}}");
      cent_text_5->AddText(eta_range->Data());
      cent_text_5->AddText(pT_range->Data());
      cent_text_5->SetFillColorAlpha(0, 0.01);

      L0bar_pz_can[pTbin][etaBin]->cd();

      gPad->SetLogy();

      L0bar_pz[pTbin][etaBin] = (TH1D*)inFile->Get(Form("L0bar_pz_pT_%i_eta_%i", pTbin, etaBin));
      L0bar_pz[pTbin][etaBin]->GetXaxis()->SetTitle("p_{z} GeV/#it{c}");
      L0bar_pz[pTbin][etaBin]->GetXaxis()->CenterTitle();
      L0bar_pz[pTbin][etaBin]->GetXaxis()->SetRangeUser( -6, 6);
      L0bar_pz[pTbin][etaBin]->GetYaxis()->SetTitle("Counts");
      L0bar_pz[pTbin][etaBin]->GetYaxis()->CenterTitle();
      L0bar_pz[pTbin][etaBin]->SetMarkerStyle(20);
      L0bar_pz[pTbin][etaBin]->SetMarkerColor(kRed);
      L0bar_pz[pTbin][etaBin]->Draw("p e");

      PYTHIA_legend_2->Draw("same");

      cent_text_5->Draw("same");


      L0bar_xF_can[pTbin][etaBin]->cd();

      gPad->SetLogy();

      L0bar_xF[pTbin][etaBin] = (TH1D*)inFile->Get(Form("L0bar_xF_pT_%i_eta_%i", pTbin, etaBin));
      L0bar_xF[pTbin][etaBin]->GetXaxis()->SetTitle("x_{F}");
      L0bar_xF[pTbin][etaBin]->GetXaxis()->CenterTitle();

      L0bar_xF[pTbin][etaBin]->GetYaxis()->SetTitle("Counts");
      L0bar_xF[pTbin][etaBin]->GetYaxis()->CenterTitle();
      L0bar_xF[pTbin][etaBin]->SetMarkerStyle(20);
      L0bar_xF[pTbin][etaBin]->SetMarkerColor(kRed);
      L0bar_xF[pTbin][etaBin]->Draw("p e");

      PYTHIA_legend_2->Draw("same");

      cent_text_5->Draw("same");
      //________________________________________________________________________________________________


      L0_cosThetaProdPlane_can[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/thetaProdlPlane/L0_cosThetaProdPlane_pT_%i_eta_%i.png", pTbin, etaBin));
      L0_thetaProdPlane_can[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/thetaProdlPlane/L0_thetaProdPlane_pT_%i_eta_%i.png", pTbin, etaBin));

      L0_pz_can[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine/L0_pz_pT_%i_eta_%i.png", pTbin, etaBin));
      L0_xF_can[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine/L0_xF_pT_%i_eta_%i.png", pTbin, etaBin));

      L0bar_cosThetaProdPlane_can[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/thetaProdlPlane/L0bar_cosThetaProdPlane_pT_%i_eta_%i.png", pTbin, etaBin));
      L0bar_thetaProdPlane_can[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/thetaProdlPlane/L0bar_thetaProdPlane_pT_%i_eta_%i.png", pTbin, etaBin));

      L0bar_pz_can[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine/L0bar_pz_pT_%i_eta_%i.png", pTbin, etaBin));
      L0bar_xF_can[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/L_kine/L0bar_xF_pT_%i_eta_%i.png", pTbin, etaBin));


      //save histograms to multi-page PDF
      if(pTbin == 0 && etaBin == 0)
      {
        L0_cosThetaProdPlane_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/thetaProdlPlane/L0_cosThetaProdPlane.pdf(", "pdf");
        L0_thetaProdPlane_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/thetaProdlPlane/L0_thetaProdPlane.pdf(", "pdf");

        L0bar_cosThetaProdPlane_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/thetaProdlPlane/L0bar_cosThetaProdPlane.pdf(", "pdf");
        L0bar_thetaProdPlane_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/thetaProdlPlane/L0bar_thetaProdPlane.pdf(", "pdf");


      }
      else if( pTbin == nPtBins && etaBin == nEtaBins)
      {
        L0_cosThetaProdPlane_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/thetaProdlPlane/L0_cosThetaProdPlane.pdf)", "pdf");
        L0_thetaProdPlane_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/thetaProdlPlane/L0_thetaProdPlane.pdf)", "pdf");

        L0bar_cosThetaProdPlane_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/thetaProdlPlane/L0bar_cosThetaProdPlane.pdf)", "pdf");
        L0bar_thetaProdPlane_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/thetaProdlPlane/L0bar_thetaProdPlane.pdf)", "pdf");

      }
      else
      {
        L0_cosThetaProdPlane_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/thetaProdlPlane/L0_cosThetaProdPlane.pdf", "pdf");
        L0_thetaProdPlane_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/thetaProdlPlane/L0_thetaProdPlane.pdf", "pdf");

        L0bar_cosThetaProdPlane_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/thetaProdlPlane/L0bar_cosThetaProdPlane.pdf", "pdf");
        L0bar_thetaProdPlane_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Simulation/PYTHIA/figs/thetaProdlPlane/L0bar_thetaProdPlane.pdf", "pdf");
      }


    }
  }


  cout<<"No. of L-Lbar pairs before cuts: "<<L0_L0bar_cosThetaProdPlane->Integral()*L0_L0bar_cosThetaProdPlane->GetXaxis()->GetBinWidth(1)<<endl;
  cout<<"No. of L-Lbar pairs after cuts: "<<L0_L0bar_cosThetaProdPlane_cuts->Integral()*L0_L0bar_cosThetaProdPlane_cuts->GetXaxis()->GetBinWidth(1)<<endl;
  cout<<endl;
  cout<<"No. of L-L pairs before cuts: "<<L0_L0_cosThetaProdPlane->Integral()*L0_L0_cosThetaProdPlane->GetXaxis()->GetBinWidth(1)<<endl;
  cout<<"No. of L-L pairs after cuts: "<<L0_L0_cosThetaProdPlane_cuts->Integral()*L0_L0_cosThetaProdPlane_cuts->GetXaxis()->GetBinWidth(1)<<endl;
  cout<<endl;
  cout<<"No. of Lbar-Lbar pairs before cuts: "<<L0bar_L0bar_cosThetaProdPlane->Integral()*L0bar_L0bar_cosThetaProdPlane->GetXaxis()->GetBinWidth(1)<<endl;
  cout<<"No. of Lbar-Lbar pairs after cuts: "<<L0bar_L0bar_cosThetaProdPlane_cuts->Integral()*L0bar_L0bar_cosThetaProdPlane_cuts->GetXaxis()->GetBinWidth(1)<<endl;



  inFile->Close();

  return;

}
