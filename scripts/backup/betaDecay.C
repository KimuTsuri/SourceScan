/*
  betaDecay.C
*/
#include <iostream>
#include <cmath>
#include "TF1.h"
#include "TH1.h"
#include "TCanvas.h"

Double_t func_decay(Double_t* x, Double_t *p) {
  // Unit: MeV
  double pe = x[0];
  double N = p[0];
  double E0 = p[1];
  double m = p[2];

  double E = std::sqrt(pe*pe + m*m);
  double de = E0 - E;
  double y = 0.0;
  if (pe < 0 || E > E0) {
    y = 0.0;
  } else {
    y = N*pe*E*de*de;
  }
  return y;
}

Double_t func_absoption(Double_t* x, DOuble_t* p) {
  double pe = x[0];
  double m = p[0];
  double x = p[1]; // rho*l

  double E = std::sqrt(pe*pe + m*m);
  double beta = pe/E;

  //double dE = func_BB(beta, x);
  double y=0;
  if (dE > E) {
    y = 1.0;
  } else {
    y = 0.0;
  }
  return y;
}

void betaDecay() {
  TF1* f_decay = new TF1("f_decay", &func_decay, 0.0, 5.0, 3);
  f_decay->SetParameter(0, 1.0);
  f_decay->SetParameter(1, 2.28);
  f_decay->SetParameter(2, 0.511);
  f_decay->SetNpx(10000);

  TCanvas* c = new TCanvas("c", "", 0, 0, 500, 500);
  TH1* hframe=0;
  hframe = c->DrawFrame(0.0, 0.0, 2.5, 1.8);
  f_decay->Draw("same");
}
