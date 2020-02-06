#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <stdio.h>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

void Style(){
  int LineWidth = 2;
  gStyle->SetLineWidth(LineWidth);
  gStyle->SetGridWidth(LineWidth);

  int FontID = 62;
  gStyle->SetStatFont(FontID);
  gStyle->SetLabelFont(FontID,"XYZ");
  gStyle->SetLabelFont(FontID,"");
  gStyle->SetTitleFont(FontID,"XYZ");
  gStyle->SetTitleFont(FontID,"");
  gStyle->SetTextFont(FontID);
  gStyle->SetLegendFont(FontID);
}

void nHits(std::string ifname/*, std::string ofname*/){
  gStyle->SetOptStat(kFALSE);
  //gStyle->SetPalette(1);
  //gStyle->SetLabelOffset(1.2);
  //gStyle->SetLabelFont(72);

  TCanvas* c = new TCanvas("canvas", "", 0, 0, 500, 800);
  c->Divide(1, 2);

  TFile* f = TFile::Open(ifname.c_str(), "READ");
  TTree* t = dynamic_cast<TTree*>(f->Get("HitTree"));

  int nevents = t->GetEntries();
  int nhits, hit_col[1000], hit_row[1000];
  t->SetBranchAddress("nhits", &nhits);
  t->SetBranchAddress("hit_col", &hit_col);
  t->SetBranchAddress("hit_row", &hit_row);

  TH2* hitmap   = new TH2F("hit map", "hit map;col;row", 136, 264, 399, 192, 0, 192);
  TH1* h_nhits = new TH1F("# of hits", "# of hits;# of hits;pixels", 1000, 0, 1000);

  int ievent, ihit, nonZeroHits = 0;
  float trig = nevents/16;
  for (ievent=0; ievent<nevents; ievent++) {
    t->GetEntry(ievent);
    h_nhits->Fill(nhits);
    if (nhits!=0) {
    nonZeroHits += nhits;
    //cout << nonZeroHits << nhits <<endl;
    }
    for (ihit=0; ihit<nhits; ihit++){
      hitmap->Fill(hit_col[ihit], hit_row[ihit]);
    }
  }
  double trigRate = nonZeroHits/trig;

  cout << "number of trigger : " << trig << endl;
  cout << "number of hit     : " << nonZeroHits << endl;
  cout << "trigger rate      : " << trigRate << endl;

  c->cd(1);
  gPad->SetTicks(1,1);
  hitmap->Draw("colz");

  c->cd(2);
  gPad->SetTicks(1,1);
  gPad->SetLogy();
  gPad->SetLogx();
  h_nhits->Draw();

  c->SaveAs("./figure/figure_nhits.png");

}
