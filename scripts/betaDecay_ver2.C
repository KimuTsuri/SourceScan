#include <iostream>
#include <string>
#include <cmath>
#include <vector>

#include "TF1.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TRandom3.h"

#include "AtomList.h"
#include "RayList.h"
#include "BetheBloch.h"

Double_t func_decay(double x) {
  double pe = x; // kinetic energy [MeV]
  double E0 = 2.28; // beta decay energy [MeV]
  double N = 1.0 ; // coefficient
  double me = 0.511; // electron mass [MeV]

  double E = std::sqrt(pe*pe + me*me);
  double dE = E0 - E;
  double f = 0.0;
  if (pe < 0 || E > E0) {
    f = 0.0;
  } else {
    f = N*pe*E*dE*dE;
  }
  return f;
}

Int_t func_trans(double x) {
  BetheBloch BB;
  AtomList atom;

  double m = atom.emass; // electron mass [MeV]
  double pe = x;
  double E = std::sqrt(pe*pe + m*m);
  double beta = pe/E;

  double dE, rho, l;
  
  double dE_H = BB.func(beta, atom.H[0], atom.H[1], atom.emass);
  double dE_C = BB.func(beta, atom.C[0], atom.C[1], atom.emass);
  double dE_N = BB.func(beta, atom.N[0], atom.N[1], atom.emass);
  double dE_O = BB.func(beta, atom.O[0], atom.O[1], atom.emass);
  double dE_Al = BB.func(beta, atom.Al[0], atom.Al[1], atom.emass);
  double dE_Si = BB.func(beta, atom.Si[0], atom.Si[1], atom.emass);
  double dE_Cu = BB.func(beta, atom.Cu[0], atom.Cu[1], atom.emass);

  /*
  Single material
  */
  l = 0.1; // thickness [cm]
  rho = 8.96; // density [g/cm3]
  dE = dE_Cu * rho * l;

  /*
  Compound material
  */
  /*
  l = 0.05; // thickness [cm]
  rho = 1.40; // density [g/cm3]
  double dE_plastic = (dE_C*14 + dE_H*10 + dE_O*4) * rho * l; // Plastic C10H8O4
  l = 0.11; // thickness [cm]
  rho = 1.08; //density [g/cm3]
  double dE_ABS = (dE_C*15 + dE_H*17 + dE_N) * rho * l; // ABS resin
  //dE = dE_ABS + dE_plastic;
  */

  int f=0;
  if (dE < E) {
    f = 0;
  } else {
    f = 1;
  }

  return f;
}

void betaDecay_ver2() {
/*
  TCanvas* c = new TCanvas("c", "", 0, 0, 400, 1000);

  TPad* pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0);
  pad1->SetGridx();
  pad1->SetTicks(1,1);
  pad1->SetLogx();
  pad1->SetLogy();
  pad1->Draw();
  pad1->cd();
  TF1* f1 = new TF1("f_decay", "func_decay(x)", 0.0, 2.5);
  f1->SetNpx(100000);
  f1->SetLineColor(kGreen+1);
  //f1->SetMinimum(0.001);
  //f1->SetMaximum(2.5);
  //f1->GetYaxis()->SetLabelSize(0.);
  TGaxis* axis1 = new TGaxis( 0.0001, 0.2, 2.5, 1.6, 20,220,510,"");
  //axis1->SetLabelFont(43);
  axis1->Draw();

  c->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx();
  pad2->Draw();
  pad2->cd();
  TF1* f2 = new TF1("f_trans", "func_trans(x)", 0.0, 2.5);
  f2->SetNpx(100000);
  f2->SetLineColor(kGreen+2);
  f2->Draw();


  c->SaveAs("./figure/betaDecay_onScinti3.pdf");
  */
  double p = 2.3;
  double Decay_results = func_decay(p);
  double dE_results = func_trans(p);

  cout << Decay_results << endl;
  cout << dE_results << endl;

}
