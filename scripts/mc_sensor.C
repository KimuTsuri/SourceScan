/*
  Monte Carlo for source on scintillator
*/
#include <iostream>
#include <string>
#include <cmath>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TCanvas.h"

class McResults {
public:
  int nGenerate;
  int nOnScinti;
  int nOnSensor;

  McResults() : nGenerate(0), nOnScinti(0), nOnSensor(0) {

  }
  float probability() const {
    if (nOnScinti == 0) return -1.0;
    return static_cast<float>(nOnSensor)/nGenerate;
  }
};

McResults generateEvents(int nTotal, 
			 TH1* h_costheta0, TH1* h_phi0, 
			 TH1* h_costheta1, TH1* h_phi1,
			 TH1* h_costheta2, TH1* h_phi2,
			 TH2F* hit_map) {
  McResults results;

  float costheta, phi, psi, r;
  int ievent=0;
  float pi = TMath::Pi();
  float theta = 0.0;
  float v[3];
  float s, x, y;

  float sourceR = 25.0/2; /* source r half size */
  r = gRandom->Uniform(0, sourceR);
  psi = gRandom->Uniform(-pi, pi);
  float sourceX = std::sqrt(r)*std::cos(psi);
  float sourceY = std::sqrt(r)*std::sin(psi);
  float sourcePosition[3] = { sourceX, sourceY, 10.0 };

  float scintiCenter[3] = { 0, 0, 10.0 };
  float scintiHalfX = 28.0/2; /* scintillator x half size */
  float scintiHalfY = 28.0/2; /* scintillator x half size */

  float sensorCenter[3] = { 0, 0,  0};
  float sensorHalfX = 20.0/2; /* rd53a x half size */
  float sensorHalfY = 11.6/2; /* rd53a y half size */

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

    double hitx = 0.0;
    double hity = 0.0;

    bool onScinti = false;
    bool onSensor = false;
    if (s > scintiCenter[2] && std::fabs(x)<scintiHalfX && std::fabs(y)<scintiHalfY) {
      onScinti = true;
      if (s > sensorCenter[2] && std::fabs(x)<sensorHalfX && std::fabs(y)<sensorHalfY) {
        onSensor = true;
        hitx = x;
        hity = y;
      }
    }
    hit_map->Fill(hitx, hity);
    h_costheta0->Fill(costheta);
    h_phi0->Fill(phi);
    results.nGenerate ++;

    if (onScinti) {
      results.nOnScinti ++;
      h_costheta1->Fill(costheta);
      h_phi1->Fill(phi);
    }

    if (onSensor) {
      results.nOnSensor ++;
      h_costheta2->Fill(costheta);
      h_phi2->Fill(phi);
    }
  }

  return results;
}

void mc_sensor(const std::string& fname, int nToGenerate, ULong_t seed) {
  TRandom3 random;
  random.SetSeed(seed);

  TFile* fout=nullptr;
  if (fname != "") {
    fout = TFile::Open(fname.c_str(), "RECREATE");
  }
  const float pi = TMath::Pi();
  TH1* h_costheta0 = new TH1F("h_costheta0", "", 100, -1, 1);
  TH1* h_costheta1 = new TH1F("h_costheta1", "", 100, -1, 1);
  TH1* h_costheta2 = new TH1F("h_costheta2", "", 100, -1, 1);
  TH1* h_phi0 = new TH1F("h_phi0", "", 60, -pi, pi);
  TH1* h_phi1 = new TH1F("h_phi1", "", 60, -pi, pi);
  TH1* h_phi2 = new TH1F("h_phi2", "", 60, -pi, pi);
  TH2F* hit_map = new TH2F("hit_map", "Hit map by mc simulation;col;raw", 400, 20.0/2, 20.0/2,
                                                                          192, 11.6/2, 11.6/2);

  McResults results = generateEvents(nToGenerate, 
				     h_costheta0, h_phi0, 
				     h_costheta1, h_phi1,
				     h_costheta2, h_phi2,
				     hit_map);

  TCanvas* c1 = new TCanvas("canvas", "", 0, 0, 400, 300);
  c1->SetTicks(1,1);
  TH1* hframe=0;
  float ymax=1.0;

  ymax = h_costheta0->GetMaximum();
  hframe = gPad->DrawFrame(-1, 0, 1, ymax*1.1);
  hframe->GetXaxis()->SetTitle("cos #theta");
  hframe->GetXaxis()->SetTitleOffset(1.2);
  hframe->GetYaxis()->SetTitle("counts");
  hframe->GetYaxis()->SetTitleOffset(1.5);
  h_costheta0->SetLineColor(kGray+3);
  h_costheta0->Draw("same");
  //h_costheta1->SetLineColor(kBlue);
  //h_costheta1->Draw("same");
  h_costheta2->SetLineColor(kRed);
  h_costheta2->Draw("same");

  c1->SaveAs("./figure/mc_sensor.pdf");

  TCanvas* c2 = new TCanvas("canvas", "", 0, 0, 400, 380);
  c2->SetTicks(1,1);
  hit_map->SetMaximum(10);
  hit_map->SetStats(0);
  hit_map->Draw("colz");
  c2->SaveAs("./figure/mc_sensor_hitmap10_ver2.pdf");


  if (fout != nullptr) {
    fout->Write();
  }
  std::cout << "Probability for the particle to hit the target: " 
	    << results.probability() << std::endl;
}

void mc_sensor() {
  mc_sensor("./scripts/mc_sensor.root", 1000000, 20191101);
}
