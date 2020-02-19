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
  double E = x; // energy [MeV]
  double E0 = 2.28; // beta decay energy [MeV]
  double N = 1.0 ; // coefficient
  double me = 0.511; // electron mass [MeV]

  //double E = std::sqrt(pe*pe + me*me);
  double pe = std::sqrt(E*E - me*me);
  double dE = E0 - E;
  double f = 0.0;
  if (pe < 0 || E > E0) {
    f = 0.0;
  } else {
    f = N*pe*E*dE*dE;
  }
  return f;
  cout <<f <<endl;
}

Int_t func_trans(double x) {
  BetheBloch BB;
  AtomList atom;

  double m = atom.emass; // electron mass [MeV]
  double E = x;
  //double E = std::sqrt(pe*pe + m*m);
  double pe = std::sqrt(E*E - m*m);
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

void betaDecay_ver3() {
  int nRow = 1;
  int nCal = 3;

  TCanvas* c = new TCanvas("c", "", 0, 0, 400, 1000);
  c->SetTopMargin(1.0);
  c->SetBottomMargin(1.0);
  c->Divide(nRow,nCal, 0, 0);


  c->cd(1);

  TF1* f_decay = new TF1("f_decay1", "func_decay(x)", 0.0, 2.5);
  f_decay->SetNpx(100000);
  gPad->SetTicks(1,1);
  //gPad->SetLogx();
  f_decay->SetLineColor(kGreen+2);
  TH1* hf1 = 0;
  hf1 = c->DrawFrame(0.0001, 0.0, 2.5, 1.8);
  f_decay->GetHistogram()->SetTitle(";;Intensity");
  f_decay->GetYaxis()->SetTitleSize(0.05);
  f_decay->Draw("same");
  double Events = f_decay->Integral(0.0,2.5);

  c->cd(2);

  TF1* f_trans = new TF1("f_trans", "func_trans(x)", 0.0, 2.5);
  f_trans->SetNpx(100000);
  gPad->SetTicks(1,1);
  //gPad->SetLogx();
  f_trans->SetLineColor(kGreen+2);
  TH1* hf2 = 0;
  hf2 = c->DrawFrame(0.0001, 0.0, 2.5, 1.8);
  f_trans->GetHistogram()->SetTitle(";;Transmission Probability");
  f_trans->Draw("same");

  c->cd(3);

  TF1* f = new TF1("f", "func_decay(x) * func_trans(x)", 0.0, 2.5);
  f->SetNpx(100000);
  gPad->SetTicks(1,1);
  //gPad->SetLogx();
  f->SetLineColor(kGreen+2);
  TH1* hf3 = 0;
  hf3 = c->DrawFrame(0.0001, 0.0, 2.5, 1.8);
  f->GetHistogram()->SetTitle(";Kinetic Energy [MeV];Intensity");
  f->Draw("same");
  double Transmission = f->Integral(0.0,2.5);

  std::cout << "|              | **data** |" << std::endl;
  std::cout << "| :----------- | :------: |" << std::endl;
  std::cout << "| All events   |  " << Events << " |" << std::endl;
  std::cout << "| Transmission |  " << Transmission << " |"<< std::endl;
  std::cout << "| Hit rate     |  " << Transmission/Events << " |" << std::endl;

  c->SaveAs("./figure/betaDecay_onCu.pdf");
}
