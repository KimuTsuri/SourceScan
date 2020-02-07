#include <iostream>
#include <cmath>
#include "TF1.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TRandom3.h"

class Results {
  public:
    int nGenerate;
    int nTransmission;

  float Transmittance() const{
    if (nTransmission == 0) return -1.0;
    return static_cast<float> nTransmission/nGenerate;
  }
}

Double_t func_decay(Double_t* x, Double_t *c) {
  // Unit: MeV
  double pe = x[0];
  double E0 = c[0];
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

Double_t func_BB(Double_t b, Double_t A, Double_t Z) {
  //double b = beta;
  //double A = a;
  //double Z = z;
  double g = 1 / std::sqrt(1 - b*b);

  double K = 0.30701; // [MeVcm2/g]
  double me = 0.511; // electron mass [MeV]

  double I = 16 * std::pow(Z, 0.9);
  double ln = std::log( 2*me*b*b*g*g /I);

  double f = 0.0;
  if (b < 0) {
    f = 0.0;
  } else {
    f = - (K*Z / (A*b*b)) * (ln - b*b) ;
  }
  return f;
}

Double_t func_abs(Double_t* x, Double_t* c) {
  double pe = x[0];
  double m = c[0];
  double A = c[1];
  double Z = c[2];

  double E = std::sqrt(pe*pe + m*m);
  double beta = pe/E;

  double dE = func_BB(beta, A, Z);

  double f=0.0;
  if (dE > E) {
    f = 1.0;
  } else {
    f = 0.0;
  }
  return f;
}

Results Events() {
  Results results;

  double pe;
  pe = gRandom->Uniform(0.0, 2.5);


  bool Transmission = false;



  if (Transmission){
    results.Transmission ++;
  }

  return results;
}

void betaDecay() {
  TCanvas* c = new TCanvas("c", "", 0, 0, 600, 1000);
  c->Divide(1,2);

  c->cd(1);
  TF1* f_decay1 = new TF1("f_decay1;;Intensity", &func_decay, 0.0, 2.5, 1);
  f_decay1->SetParameter(0, 2.28);
  f_decay1->SetNpx(10000);
  TH1* hframe=0;
  hframe = c->DrawFrame(0.0, 0.0, 2.5, 1.8);
  f_decay1->Draw("same");

  c->cd(2);
  TF1* f_abs = new TF1(";Kinetic energy [MeV];Absorption rate", &func_abs, 0.0, 2.5, 4);
  f_abs->SetParameter(0, 0.511); //m
  f_abs->SetParameter(1, 28.0855); //A
  f_abs->SetParameter(2, 14); //Z
  f_abs->Draw();

}
