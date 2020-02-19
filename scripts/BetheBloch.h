#ifndef BETHEBLOCH_H_
#define BETHEBLOCH_H_

#include <iostream>
#include <cmath>
#include "TF1.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TGaxis.h"
#include "TMath.h"

class BetheBloch {
  public:
    static double func(double b, double A, double Z, double M) {

      double K = 0.307075; // [MeVcm2/g]
      double me = 0.511; // electron mass [MeV]
      double g = 1 / std::sqrt(1 - b*b);
      double I = 16 * std::pow(Z, 0.9) * std::pow(10, -6);
      double T = 2*me*b*b*g*g / ( 1 + 2*g*me/M + (me/M)*(me/M) );
      double ln = std::log( 2.0*me*b*b*g*g * T / (I*I));
      double f = (K*Z / (A*b*b)) * (ln/2.0 - b*b);

      return f;
    }
};

#endif