#ifndef ALDERPARAMS_HPP
#define ALDERPARAMS_HPP

#include <set>

namespace ALD {

  class AlderParams {

    void check_file_readable(const char *filename);
    void check_file_readable_if_specified(const char *filename);
    void check_file_writable(const char *filename);
    void check_pars(void);
    std::set <int> parse_to_set(char *str);
    void print_param_settings(void);

  public:
    static const double DEFAULT_BINSIZE;
    static const double MINDIS_NOT_SET;
    static const double DEFAULT_FIT_START;
    static const double DEFAULT_MAXDIS;
    static const int DEFAULT_MINCOUNT;

    char *genotypename, *snpname, *indivname, *badsnpname, *poplistname, *refpops, *admixpop,
      *admixlist, *raw_outname, *weightname, *chrom, *nochrom;
    int mincount;
    double mindis, maxdis, binsize; 
    int checkmap, verbose, num_threads;
    bool print_raw_jackknife, use_naive_algo, fast_snp_read, approx_ld_corr;
    std::set <int> chrom_set, nochrom_set;
    bool print_jackknife_fits;

    AlderParams(void);
    void readcommands(int argc, char **argv, const char *VERSION);
  };

}

#endif
