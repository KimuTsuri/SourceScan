#include <stdio.h>
#include <math.h>
#include <cmath>

#include "AtomList.h"
#include "RayList.h"

class SigmaPhotoelectric {
  public:
  
  static double func(double Z, double E) {
    if (E<50) {
      E0 = 30;
      A = (1.53 * std::exp(-0.0361*Z) + 0.510 * std::exp(-0.176*Z)) * std::pow(10,-4);
      B = 0.0515*Z -0.18 *np.exp(-0.215*Z);
    } else {
      E0 = 100;
      A = (2.73 * std::exp(-0.02113*Z);
      B = 0.008*Z - 0.18 * std::exp(-0.107*Z);
    }
    return (A * std::pow(Z,5) * std::pow((E/E0), -3.32 + B)) * std::pow(10,-24);
  }
};

void gammaEstimate() {
  AtomList atom;
  RayList ray;
  SigmaPhotoelectric sigmaPE;

  sigmaPE.fumc(atom.Si[1], ray.Am[1]);
  cout << sigmaPE << endl;
}
