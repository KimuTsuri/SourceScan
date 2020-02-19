#ifndef ATOMLIST_H_
#define ATOMLIST_H_

class AtomList {
  // vector elements { A, Z, rho}
  public:
    double emass = 0.551; // [MeV]
    std::vector<double> H =  { 1.008,  1, 0};
    std::vector<double> C =  { 12.01,  6, 0};
    std::vector<double> N =  { 14.01,  7, 0};
    std::vector<double> O =  { 16.00,  8, 0};
    std::vector<double> Al = { 26.98, 13, 0};
    std::vector<double> Si = { 28.09, 14, 0};
    std::vector<double> Cu = { 63.55, 29, 8.96};
};

#endif
