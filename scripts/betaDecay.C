/*
  betaDecay.C
*/
#include <iostream>
#include <cmath>
#include "TF1.h"
#include "TH1.h"
#include "TCanvas.h"

Double_t func_decay(Double_t* x, Double_t *y) {
  // Unit: MeV
  double pe = x[0];
  double E0 = y[0];
  double N = 1.0 ; //coefficient
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
/*
Double_t func_BB(Double_t* beta, Double_t* x, Double_t* a, Double_t* z) {
  double b = beta;
  double dx = x;
  double A = a;
  double Z = z;
  double K = 0.307075; // [MeVcm2/g]
  double me = 0.511; // electron mass [MeV]

  double I = 16 * std::pow(Z,0.9);
  double g = 1/std::sqrt(1 - b*b); // gamma
  double ln = std::log(2*me*b*g/I);
  double f = K*z*z* Z * ln / (A*b*b)
}

Double_t func_abs(Double_t* x, Double_t* y) {
  double pe = x[0];
  double m = y[0];
  double x = y[1]; // rho*l
  double A = y[2];
  double Z = y[3];

  double E = std::sqrt(pe*pe + m*m);
  double beta = pe/E;

  double dE = func_BB(beta, x, A, Z);

  double f=0.0;
  if (dE > E) {
    f = 1.0;
  } else {
    f = 0.0;
  }
  return f;
}
*/
void betaDecay() {
  TCanvas* c = new TCanvas("c", "", 0, 0, 500, 500);

  TF1* f_decay = new TF1("f_decay", &func_decay, 0.0, 5.0, 3);
  f_decay->SetParameter(0, 2.28);
  f_decay->SetNpx(10000);

  TH1* hframe=0;
  hframe = c->DrawFrame(0.0, 0.0, 2.5, 1.8);
  f_decay->Draw("same");





}
