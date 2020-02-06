//
// Overwrite two histograms
//
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include <stdio.h>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

void nHits_2histo(std::string fname_w, std::string fname_wo){
  gStyle->SetOptStat(0);

  TCanvas* c = new TCanvas("canvas", "", 0, 0, 500, 400);

  TFile* f_w  = TFile::Open(fname_w.c_str(), "READ");
  TTree* t_w  = dynamic_cast<TTree*>(f_w->Get("HitTree"));
  int nevents_w  = t_w->GetEntries();

  TFile* f_wo = TFile::Open(fname_wo.c_str(), "READ");
  TTree* t_wo = dynamic_cast<TTree*>(f_wo->Get("HitTree"));
  int nevents_wo =  t_wo->GetEntries();

  int nhits_w, nhits_wo;
  t_w ->SetBranchAddress("nhits", &nhits_w);
  t_wo->SetBranchAddress("nhits", &nhits_wo);

  TH1F* h_nhits_w  = new TH1F("# of hits", ";# of hits;pixels", 1000, 0, 1000);
  TH1F* h_nhits_wo = new TH1F("# of hits", ";# of hits;pixels", 1000, 0, 1000);

  int ievent, Events_w, Events_wo = 0;
  for (ievent=0; ievent<nevents_w; ievent++) {
    t_w ->GetEntry(ievent);
    h_nhits_w ->Fill(nhits_w);
    if (nhits_w!=0) {
      Events_w += nhits_w;
    }
  }
  for (ievent=0; ievent<nevents_wo; ievent++) {
    t_wo->GetEntry(ievent);
    h_nhits_wo->Fill(nhits_wo);
  }

  gPad->SetTicks(1,1);
  gPad->SetLogy();
  gPad->SetLogx();
  h_nhits_w ->SetLineColor(kRed);
  h_nhits_wo->SetLineColor(kGray+3);
  h_nhits_w ->Draw("same");
  h_nhits_wo->Draw("same");

  TLegend* legend = new TLegend(0.62, 0.77, 0.8, 0.86);
  legend->AddEntry(h_nhits_w,"with source");
  legend->AddEntry(h_nhits_wo,"without source");
  legend->SetTextSize(0.03);
  legend->Draw();

  c->SaveAs("./figure/figure_nhits_2histo.png");
}
