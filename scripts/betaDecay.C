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

Double_t func_decay(double x) {
  double pe = x; // kinetic energy [MeV]
  double E0 = 2.28; // beta decay energy [MeV]
  double N = 2.0 ; // coefficient
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

  double m = atom.emass; // electron mass
  double pe = x;
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
  //double l = 0.0001; // thickness [cm]
  //double rho = atom.Si[2] * 0.001 ; // density [g/cm3]
  //dE = dE_Si * rho * l;
  
  /*
  Compound material
  */
  // Plastic C10H8O4
  //double l = 0.05; // thickness [cm]
  //double rho = 1.40; // density [g/cm3]
  //dE = (dE_C*14 + dE_H*10 + dE_O*4) * rho * l;
  // ABS resin
  double l = 0.1; // thickness [cm]
  double rho = 1.08; //density [g/cm3]
  dE = (dE_C*15 + dE_H*17 + dE_N) * rho * l;

  int f=0;
  if (dE > E) {
    f = 0;
  } else {
    f = 1;
  }

  return f;
}

void betaDecay() {
  gStyle->SetOptStat(0);

  int nRow = 1;
  int nCal = 3;

  TCanvas* c = new TCanvas("c", "", 0, 0, 600, 1000);
  c->SetTopMargin(1.0);
  c->SetBottomMargin(1.0);
  c->Divide(nRow,nCal, 0, 0);

  c->cd(1);
  TF1* f_decay = new TF1("f_decay1", "func_decay(x)", 0.0, 2.5);
  f_decay->SetNpx(100000);
  f_decay->GetHistogram()->SetTitle(";;Intensity");
  f_decay->GetHistogram()->SetLineColor(kBlue);
  gPad->SetTicks(1,1);
  f_decay->Draw();
  double Events = f_decay->Integral(0.0,2.4);

  c->cd(2);
  TF1* f_trans = new TF1("f_trans", "func_trans(x)", 0.0, 2.5);
  f_trans->SetNpx(100000);
  f_trans->GetHistogram()->SetTitle(";;Transmission Probability");
  f_trans->GetHistogram()->SetLineColor(kGray);
  gPad->SetTicks(1,1);
  f_trans->Draw();

  c->cd(3);
  TF1* f = new TF1("f", "func_decay(x) * func_trans(x)", 0.0, 2.5);
  f->SetNpx(100000);
  f->GetHistogram()->SetTitle(";Kinetic Energy [MeV];Intensity");
  f->GetHistogram()->SetLineColor(kBlue);
  gPad->SetTicks(1,1);
  f->Draw();
  double Transmission = f->Integral(0.0,2.4);

  cout << "|              | **data** |" << endl;
  cout << "| :----------- | :------: |" << endl;
  cout << "| All events   |  " << Events << " |" << endl;
  cout << "| Transmission |  " << Transmission << " |"<< endl;
  cout << "| Hit rate     |  " << Transmission/Events << " |" << endl;

  c->SaveAs("./figure/betaDecay.pdf");
}
