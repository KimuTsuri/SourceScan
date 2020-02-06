/*
  ex_mc1.C
*/
#include <iostream>
#include <string>
#include <cmath>

#include "TH1.h"
#include "TFile.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TCanvas.h"

class McResults {
public:
  int nGenerated;
  int nOnTarget;

  McResults() : nGenerated(0), nOnTarget(0) {

  }
  float probability() const {
    if (nGenerated == 0) return -1.0;
    return static_cast<float>(nOnTarget)/nGenerated;
  }
};

McResults generateEvents(int nTotal, 
			 TH1* h_costheta0, TH1* h_phi0, 
			 TH1* h_costheta1, TH1* h_phi1) {
  McResults results;
  float sourcePosition[3] = { 0, 0, 20.0 };
  float targetCenter[2] = { 0, 0 };
  float targetHalfX = 20.0/2; /* rd53a x half size */
  float targetHalfY = 11.8/2; /* rd53a y half size */
  float costheta, phi;
  int ievent=0;
  float pi = TMath::Pi();
  float theta = 0.0;
  float v[3];
  float s, x, y;

  for (ievent=0; ievent<nTotal; ++ievent) {
    costheta = gRandom->Uniform(-1.0, 1.0);
    phi = gRandom->Uniform(-pi, pi);

    theta = std::acos(costheta);
    v[0] = std::sin(theta)*std::cos(phi);
    v[1] = std::sin(theta)*std::sin(phi);
    v[2] = std::cos(theta);
    s = -sourcePosition[2]/v[2];
    x = sourcePosition[0] + v[0]*s;
    y = sourcePosition[0] + v[1]*s;

    bool onTarget = false;
    if (s > 0 && std::fabs(x)<targetHalfX && std::fabs(y)<targetHalfY) {
      onTarget = true;
    }

    h_costheta0->Fill(costheta);
    h_phi0->Fill(phi);
    results.nGenerated ++;
    if (onTarget) {
      results.nOnTarget ++;
      h_costheta1->Fill(costheta);
      h_phi1->Fill(phi);
    }
  }

  return results;
}

void ex_mc1(const std::string& fname, int nToGenerate, ULong_t seed) {
  TRandom3 random;
  random.SetSeed(seed);

  TFile* fout=nullptr;
  if (fname != "") {
    fout = TFile::Open(fname.c_str(), "RECREATE");
  }
  const float pi = TMath::Pi();
  TH1* h_costheta0 = new TH1F("h_costheta0", "", 100, -1, 1);
  TH1* h_costheta1 = new TH1F("h_costheta1", "", 100, -1, 1);
  TH1* h_phi0 = new TH1F("h_phi0", "", 60, -pi, pi);
  TH1* h_phi1 = new TH1F("h_phi1", "", 60, -pi, pi);

  McResults results = generateEvents(nToGenerate, 
				     h_costheta0, h_phi0, 
				     h_costheta1, h_phi1);

  TCanvas* c = new TCanvas("c", "", 0, 0, 600, 600);
  TH1* hframe=0;
  float ymax=1.0;
  c->Divide(1, 2);

  c->cd(1);
  ymax = h_costheta0->GetMaximum();
  hframe = gPad->DrawFrame(-1, 0, 1, ymax*1.1);
  h_costheta1->SetLineColor(kRed);
  h_costheta0->Draw("same");
  h_costheta1->Draw("same");

  c->cd(2);
  ymax = h_phi0->GetMaximum();
  hframe = gPad->DrawFrame(-pi, 0, pi, ymax*1.1);
  h_phi1->SetLineColor(kRed);
  h_phi0->Draw("same");
  h_phi1->Draw("same");
  
  if (fout != nullptr) {
    fout->Write();
  }
  std::cout << "Probability for the particle to hit the target: " 
	    << results.probability() << std::endl;
}

void ex_mc1() {
  ex_mc1("mc1.root", 1000000, 20191101);
}
