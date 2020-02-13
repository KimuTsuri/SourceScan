#include <iostream>
#include <cmath>
#include "TF1.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TGaxis.h"

#include "AtomList.h"

Double_t func_BB(Double_t* x, Double_t* c) {
  double bg = x[0];
  double A = c[1];
  double Z = c[2];

  double K = 0.30701; // [MeVcm2/g]
  double me = 0.511; // electron mass [MeV]

  double b2 = bg*bg / (1 + bg*bg);
  double g = bg / std::sqrt(b2);

  double I = 16 * std::pow(Z, 0.9);
  double ln = std::log( 2*me*b2*g*g /I);

  double f = 0.0;
  if (bg < 0) {
    f = 0.0;
  } else {
    f = - (K*Z / (A*b2)) * (ln - b2) ;
  }
  return f;
}

void BetheBloch() {
  AtomList atom;

  double xmin = 0.01;
  double xmax = 10.0;
  double ymin = 0.1;
  double ymax = 10000.;

  TCanvas* c = new TCanvas("c", "", 0, 10, 600, 500);
  c->SetLogx();
  c->SetLogy();
  c->SetTicks(1,1);

  TF1* f_BB = new TF1("Bethe Bloch", &func_BB, xmin, xmax, 3);
  f_BB->SetParameter(1, atom.Cu[0]); // A　[/mol]
  f_BB->SetParameter(2, atom.Cu[1]);  // Z [g/mol]
  f_BB->SetNpx(100000);
  f_BB->GetXaxis()->SetTitle("#beta#gamma");
  f_BB->GetXaxis()->CenterTitle();
  f_BB->GetYaxis()->SetTitle("Mass stopping power [MeV cm^{2}/g]");

  f_BB->Draw();

  c->SaveAs("./figure/BetheBloch_Cu.pdf");
}
