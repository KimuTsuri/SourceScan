#include <iostream>
#include <cmath>
#include "TF1.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TGaxis.h"
#include "TMath.h"

#include "AtomList.h"

Double_t func_BB(Double_t* x, Double_t* c) {
  double p = x[0];
  double A = c[0];
  double Z = c[1];

  double K = 0.307075; // [MeVcm2/g]
  double me = 0.511; // electron mass [MeV]
  double M = c[2];

  double b2 = p*p / (me*me + p*p);
  double g = 1 / std::sqrt(1-b2);

  double I = 16 * std::pow(Z, 0.9) * std::pow(10, -6);
  double T = 2*me*b2*g*g / ( 1 + 2*g*me/M + (me/M)*(me/M) );
  double ln = std::log( 2.0*me*b2*g*g * T / (I*I));
  double f = (K*Z / (A*b2)) * (ln/2.0 - b2) ;

  return f;
}

void BetheBloch() {
  AtomList atom;

  double xmin = 0.1;
  double xmax = 1000.0;
  double ymin = 0.1;
  double ymax = 100.;

  TCanvas* c = new TCanvas("c", "", 0, 10, 600, 500);
  c->SetLogx();
  //c->SetLogy();
  c->SetTicks(1,1);

  TF1* f_BB = new TF1("Bethe Bloch", &func_BB, xmin, xmax, 3);
  f_BB->SetParameter(0, atom.Cu[0]); // A　[/mol]
  f_BB->SetParameter(1, atom.Cu[1]);  // Z [g/mol]
  f_BB->SetParameter(2, atom.emass/*atom.Cu[1] * std::pow(10,10)  /(1.78 * 6.02)*/); // mass [MeV]
  f_BB->SetNpx(100000);
  f_BB->SetLineColor(kBlack);
  f_BB->GetXaxis()->SetTitle("Momentum [MeV]");
  //f_BB->GetXaxis()->SetTitleSize(0.001);
  f_BB->GetXaxis()->CenterTitle();
  f_BB->GetYaxis()->SetTitle("Mass stopping power [MeV cm^{2}/g]");

  f_BB->Draw();

  c->SaveAs("./figure/BetheBloch_emass_Cu.pdf");
}
