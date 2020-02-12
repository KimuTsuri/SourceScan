#include <stdio.h>
#include <math.h>
#include <cmath>

#include "AtomList.h"
#include "RayList.h"

class Sigma {
  public:
    double E0, x, A, B, C, D;
  double PhotoElectric (double Z, double E) {
    if (E<50) {
      E0 = 30;
      A = (1.53 * std::exp(-0.0361*Z) + 0.510 * std::exp(-0.176*Z)) * std::pow(10,-4);
      B = 0.0515*Z - 0.18 * std::exp(-0.215*Z);
    } else {
      E0 = 100;
      A = (2.73 * std::exp(-0.02113*Z) + 0.860*std::exp(-0.128*Z)) / (1 - (2.63/Z + 0.473) * (E-100)*std::pow(10,-3)) * std::pow(10,-6);
      B = 0.008*Z - 0.18 * std::exp(-0.107*Z);
    }
    return A * std::pow(Z,5) * std::pow((E/E0), -3.32+B) * std::pow(10,-24);
  }
  double Compton (double Z, double E) {
    if (Z<4) {
      A = Z/5;
    } else {
      A = Z-3;
    }
    if (E<100) {
      B = 5.5/(Z+5) + 1.06;
    } else {
      B = 1.5;
    }
    double S = 1/ (1 + 5.184*std::pow(10,-3) * std::pow(A,0.8107) * std::pow(E/100,-B));
    double a = E/511;
    double sigmaKN = 0.49893 * ( (1+a)/(a*a) * ( 2*(1+a)/(1+2*a) - 1/a * std::log(1+2*a) + 1/(2*a) * std::log(1+2*a) - (1+3*a)/((1+2*a)*(1+2*a)) ) );
    //return ( Z * sigmaKN * S ) * std::pow(10,-24);
    return sigmaKN * std::pow(10,-24);
  }
  double Coherent (double Z, double E) {
    x = std::log(E/100);
    A = -6.6514 + 2.2802*std::log(Z) + 0.04174*std::log(Z)+std::log(Z);
    B = 3*std::pow(10,-4)*Z - 1.89;
    C = -0.03;
    D = -1.4 * std::exp((3*std::pow(10,-4)*Z - 0.08)*E);
    return std::exp( A + B*x + C*x*x + D ) * std::pow(10,-24);
  }
};


void gammaEstimate() {
  AtomList atom;
  RayList ray;
  Sigma sigma;

  double sigma_PE = sigma.PhotoElectric(atom.Al[1], ray.Am[1]);
  double sigma_Comp = sigma.Compton(atom.Al[1], ray.Am[1]);
  double sigma_Cohe = sigma.Coherent(atom.Al[1], ray.Am[1]);
  double sigma_tot = sigma_PE + sigma_Comp + sigma_Cohe;
  cout << "total sigma   : " << sigma_tot << endl;
  cout << "PhotoElectric : " << sigma_PE  << endl;
  cout << "Compton       : " << sigma_Comp << endl;
  cout << "Coherent      : " << sigma_Cohe << endl;
}
