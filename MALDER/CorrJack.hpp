#ifndef CORRJACK_HPP
#define CORRJACK_HPP

#include <string>
#include <vector>
#include <utility>
namespace ALD {

  class Corr {
  public:
    double count, sum_x, sum_y, sum_xy, sum_x2, sum_y2;
    void add_term(double x, double y);
    void add_unbiased_sq_term(double sq_term); // augments both x2 and y2
  };

  class CorrJack {

    const int C; // number of chromosomes on which to jackknife
    const std::vector <std::string> &jack_ind_ids; // labels of left-out chroms for output purposes
    // central = true for usual corr, false for non-zeroed
    std::pair <double, double> jackknife_corr(bool central);
  
  public:
    std::vector <Corr> data;
  
    CorrJack(int _C, const std::vector <std::string> &_ids);

    std::pair <double, double> jackknife_corr(void);
    std::pair <double, double> jackknife_cos(void);
    std::pair <double, double> jackknife_cos_polyache_denom(const CorrJack &test_data,
							    const CorrJack &ref_data);
    std::pair <double, double> jackknife_x2_avg(void);
    double tot_count(void);
    double tot_sum_x2(void);
  };

}

#endif
