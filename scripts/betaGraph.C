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

class BetheBloch {
  public:
  static double func(double b, double A, double Z) {
    AtomList atom;
    
    double K = 0.30701; // [MeVcm2/g]
    double me = atom.emass; // electron mass [MeV]
    double g = 1 / std::sqrt(1 - b*b);
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
};

Double_t func_decay(Double_t* x, Double_t *c) {
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

Double_t func_abs(Double_t* x, Double_t* c) {
  BetheBloch BB;
  AtomList atom;

  double m = atom.emass; // electron mass
  double pe = x[0];
  double E = std::sqrt(pe*pe + m*m);
  double beta = pe/E;

  double dE = 0.0;
  
  double dE_H = BB.func(beta, atom.H[0], atom.H[1]);
  double dE_C = BB.func(beta, atom.C[0], atom.C[1]);
  double dE_N = BB.func(beta, atom.N[0], atom.N[1]);
  double dE_O = BB.func(beta, atom.O[0], atom.O[1]);
  double dE_Al = BB.func(beta, atom.Al[0], atom.Al[1]);
  double dE_Si = BB.func(beta, atom.Si[0], atom.Si[1]);
  double dE_Cu = BB.func(beta, atom.Cu[0], atom.Cu[1]);
  
  /*
  Single material
  */  
  //dE = dE_Si;
  
  /*
  Compound material
  */
  //dE = dE_C*14 + dE_H*10 + dE_O*4; // Plastic
  dE = dE_C*15 + dE_H*17 + dE_N;   // ABS resin


  double f=0.0;
  if (dE > E) {
    f = 1.0;
  } else {
    f = 0.0;
  }

  return f;
}

void betaGraph() {
  AtomList atom;
  RayList ray;
  
  TCanvas* c = new TCanvas("c", "", 0, 0, 600, 1000);
  c->Divide(1,2);

  c->cd(1);
  TF1* f_decay1 = new TF1("f_decay1;;Intensity", &func_decay, 0.0, 2.5, 1);
  f_decay1->SetParameter(0, 2.28);
  f_decay1->SetNpx(10000);
  TH1* hf1=0;
  hf1 = c->DrawFrame(0.0, 0.0, 2.5, 1.8);
  f_decay1->Draw();

  c->cd(2);
  TF1* f_abs = new TF1(";Kinetic energy [MeV];Absorption rate", &func_abs, 0.0, 2.5, 0);
  //f_abs->SetParameter(0, 1);
  //f_abs->SetParameter(1, atom.Si[0]); //m
  //f_abs->SetParameter(2, atom.Si[1]); //A
  //f_abs->SetParameter(3, atom.Si[2]); //Z
  f_abs->Draw();

}
