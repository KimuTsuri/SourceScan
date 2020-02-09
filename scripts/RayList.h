#ifndef RAYLIST_H_
#define RAYLIST_H_

class RayList {
  // vector elements { Bq, decay energy}
  public:
    // beta ray
    std::vector<double> Sr = { 0.0049, 2.28};
    // gamma ray
    std::vector<double> Cs = { 1.0, 0.546};
    std::vector<double> Am = { 1.0, 0.049};
};

#endif
