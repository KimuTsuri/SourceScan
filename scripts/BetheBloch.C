#include <iostream>
#include <cmath>
#include "TF1.h"
#include "TH1.h"
#include "TCanvas.h"

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
  TCanvas* c = new TCanvas("c", "", 0, 10, 800, 600);

  TF1* f_BB = new TF1("Bethe Bloch;bg;Mass stopping power [MeV cm^2/g]", &func_BB, 0.001, 10.0, 3);
  f_BB->SetParameter(1, 28.0855); // A　[/mol]
  f_BB->SetParameter(2, 14);  // Z [g/mol]
  f_BB->SetNpx(100000);

  //TH1* hframe=0;
  //hframe = c->DrawFrame(0.0, 0.0, 2.5, 1.8);
  //f_BB->Draw("same");
  gPad->SetLogx();
  gPad->SetLogy();
  f_BB->Draw();

  c->SaveAs("./figure/BetheBloch_Si.png");
}
