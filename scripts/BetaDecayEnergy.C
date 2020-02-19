
#include "TF1.h"
#include "TH1.h"
#include "TCanvas.h"

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
  return f;
}

void BetaDecayEnergy() {
  TCanvas* c = new TCanvas("c", "", 0, 10, 600, 400);

  TF1* f_decay = new TF1("f_decay1", "func_decay(x)", 0.0, 2.5);
  f_decay->SetNpx(100000);
  gPad->SetTicks(1,1);
  f_decay->GetHistogram()->SetTitle("Beta Decay Energy;Kinetic Energy [MeV];Intensity");
  f_decay->SetLineColor(kBlack);


  f_decay->Draw();

  c->SaveAs("./figure/BetaDecayEnergy.pdf");
}
