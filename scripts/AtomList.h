#ifndef ATOMLIST_H_
#define ATOMLIST_H_

class AtomList {
  // vector elements { A, Z, rho}
  public:
    double emass = 0.551; // [MeV]
    std::vector<double> H =  {  1.08,  1, 0.00000837};
    std::vector<double> C =  { 12.01,  6, 2.265};
    std::vector<double> N =  { 14.01,  7, 0.0001165};
    std::vector<double> O =  { 16.00,  8, 0.0001332};
    std::vector<double> Al = { 26.98, 13, 2.699};
    std::vector<double> Si = { 28.09, 14, 2.329};
    std::vector<double> Cu = { 63.55, 29, 8.960};
};

#endif
