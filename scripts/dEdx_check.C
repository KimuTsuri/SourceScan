#include <iostream>
#include <string>
#include <cmath>
#include <vector>

#include "TF1.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TLegend.h"

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
  return E;
}

Double_t func_trans(double x) {
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
  //l = 0.1; // thickness [cm]
  //rho = 8.96; // density [g/cm3]
  //dE = dE_Cu * rho * l;

  l = 0.05; // thickness [cm]
  rho = 1.40; // density [g/cm3]
  double dE_plastic = (dE_C*0.625017 + dE_H*0.625017 + dE_O*0.333025 ) * rho * l; // Plastic C10H8O4
  //l = 0.11; // thickness [cm]
  //rho = 1.08; //density [g/cm3]
  //double dE_ABS = (dE_C*15/atom.C[0] + dE_H*17/atom.H[0] + dE_N/atom.N[0]) * rho * l; // ABS resin
  dE = /*dE_ABS +*/ dE_plastic;

  return dE;
}

void dEdx_check() {

  TCanvas* c = new TCanvas("c", "", 0, 0, 600, 500);
  //c->SetLogy();
  //c->SetLogx();
  TF1* f_decay = new TF1("f_decay1", "func_decay(x)", 0.0, 2.5);
  f_decay->SetNpx(100000);
  f_decay->SetLineColor(kGreen+2);
  TF1* f = new TF1("f", "func_trans(x)", 0.0, 2.5);
  f->SetNpx(100000);
  f->SetLineColor(kRed);
  TH1* hf1 = 0;
  hf1 = c->DrawFrame(0.0001, 0.0, 2.5, 3);
  hf1->SetTitle(";p [MeV];Energy [MeV]");
  f_decay->Draw("same");
  f->Draw("same");

  TLegend *leg = new TLegend(0.2, 0.7, 0.5, 0.9, "");
  leg->AddEntry(f_decay, "Energy for momentum", "l");
  leg->AddEntry(f, "Energy calculated from dEdx", "l");
  leg->SetFillStyle(0);
  leg->Draw();


  c->SaveAs("./figure/betaDecay_Momentum-Energy_plastic.pdf");

}

