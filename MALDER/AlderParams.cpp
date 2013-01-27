#include <set>
#include <fstream>

#include "MiscUtils.hpp" // these have to be included before Nick's... some compiler error o/w
#include "AlderParams.hpp"

#include "nicklib.h"
#include "getpars.h"
#include "admutils.h"

namespace ALD {

  const double AlderParams::DEFAULT_BINSIZE = 0.0005;
  const double AlderParams::MINDIS_NOT_SET = -1;
  const double AlderParams::DEFAULT_FIT_START = 0.005;
  const double AlderParams::DEFAULT_MAXDIS = 0.5;
  const int AlderParams::DEFAULT_MINCOUNT = 4;

  void AlderParams::check_file_readable(const char *filename) {
    std::ifstream fin(filename);
    if (!fin.good())
      fatalx("input file does not exist or cannot be read: %s\n", filename);
    fin.close();
  }

  void AlderParams::check_file_readable_if_specified(const char *filename) {
    if (filename != NULL)
      check_file_readable(filename);
  }

  void AlderParams::check_file_writable(const char *filename) {
    std::ofstream fout(filename, std::ios::app);
    if (!fout.good())
      fatalx("output file cannot be written: %s\n", filename);
    fout.close();
  }

  void AlderParams::check_pars(void) {
    if (genotypename == NULL) fatalx("required parameter 'genotypename' not provided\n");
    if (snpname == NULL) fatalx("required parameter 'snpname' not provided\n");
    if (indivname == NULL) fatalx("required parameter 'indivname' not provided\n");
    check_file_readable(genotypename);
    check_file_readable(snpname);
    check_file_readable(indivname);
    check_file_readable_if_specified(badsnpname);
    check_file_readable_if_specified(admixlist);
    check_file_readable_if_specified(poplistname);
    check_file_readable_if_specified(weightname);
    if (raw_outname != NULL) check_file_writable(raw_outname);

    if ((weightname != NULL) + (poplistname != NULL) + (refpops != NULL) > 1)
      fatalx("cannot specify more than one of weight file, ref poplist file or refpops string\n") ;

    if (weightname == NULL && poplistname == NULL && refpops == NULL)
      fatalx("must specify weight file, ref poplist file or refpops string (semicolon-delim)\n") ;

    if (admixpop != NULL && admixlist != NULL)
      fatalx("cannot specify both admix pop and admixlist file\n") ;

    if (admixpop == NULL && admixlist == NULL)
      fatalx("must specify either admix pop or admixlist file\n") ;

    if (print_raw_jackknife && raw_outname == NULL)
      printf("warning: ignoring 'jackknife' parameter, which only applies to raw output\n");

    if (mincount < 2)
      fatalx("mincount param must be at least 2 to compute LD\n") ;

    if (fast_snp_read && badsnpname != NULL)
      fatalx("cannot do fast_snp_read if badsnpname specified\n") ;
    
    if (chrom != NULL && nochrom != NULL)
      fatalx("cannot specify both chrom list and nochrom list\n");
  }

  std::set <int> AlderParams::parse_to_set(char *str) {
    for (int i = 0; i < (int) strlen(str); i++)
      if (!(str[i] == ';' || isdigit(str[i])))
	fatalx("chrom and nochrom strings must be semicolon-delimited int lists\n");
    std::set <int> ret;
    for (int i = 0; i < (int) strlen(str); i++)
      if (i == 0 || str[i-1] == ';') {
	int c; sscanf(str+i, "%d", &c);
	ret.insert(c);
      }
    for (std::set <int>::iterator it = ret.begin(); it != ret.end(); it++)
      printf("%d ", *it);
    printf("\n");
    return ret;
  }

  void AlderParams::print_param_settings(void) {
    printf("---------- parameter settings used (with defaults for unspecified) ----------\n");
    printf("\nInput data files:\n");
    printf("%20s: %s\n", "genotypename", genotypename);
    printf("%20s: %s\n", "snpname", snpname);
    printf("%20s: %s\n", "indivname", indivname);
    if (badsnpname != NULL) printf("%20s: %s\n", "badsnpname", badsnpname);

    printf("\nAdmixed population:\n");
    printf("%20s: %s\n", "admixpop", admixpop);

    printf("\nReference populations/weights:\n");
    if (refpops != NULL) printf("%20s: %s\n", "refpops", refpops);
    if (poplistname != NULL) printf("%20s: %s\n", "poplistname", poplistname);
    if (weightname != NULL) printf("%20s: %s\n", "weightname", weightname);

    printf("\nRaw weighted LD curve output:\n");
    printf("%20s: %s\n", "raw_outname", raw_outname == NULL ? "(none)" : raw_outname);
    if (raw_outname != NULL) printf("%20s: %s\n", "jackknife", print_raw_jackknife ? "YES" : "NO");
    
    printf("\nData filtering:\n");
    printf("%20s: %d\n", "mincount", mincount);
    if (chrom != NULL) printf("%20s: %s\n", "chrom", chrom);
    if (nochrom != NULL) printf("%20s: %s\n", "nochrom", nochrom);
  
    printf("\nCurve fitting:\n");
    printf("%20s: %f\n", "binsize", binsize);
    printf("%20s: %f\n", "mindis", mindis);
    printf("%20s: %f\n", "maxdis", maxdis);
  
    printf("\nInput checks:\n");
    printf("%20s: %s\n", "fast_snp_read", fast_snp_read ? "YES" : "NO");
    printf("%20s: %s\n", "checkmap", checkmap ? "YES" : "NO");
    
    printf("\nComputational options:\n");
    printf("%20s: %d\n", "num_threads", num_threads);
    printf("%20s: %s\n", "approx_ld_corr", approx_ld_corr ? "YES" : "NO");
    printf("%20s: %s\n", "use_naive_algo", use_naive_algo ? "YES" : "NO");
    
    printf("\n");
  }
    
  // public functions

  AlderParams::AlderParams(void) {
    genotypename = NULL ;
    indivname = NULL ;
    badsnpname = NULL ;
    poplistname = NULL ;
    refpops = NULL ;
    admixpop = NULL ;
    admixlist = NULL ;
    raw_outname = NULL ;
    weightname = NULL ; 
    mincount = DEFAULT_MINCOUNT ;
    mindis = MINDIS_NOT_SET ;
    maxdis = DEFAULT_MAXDIS ;
    binsize = DEFAULT_BINSIZE ; 
    print_raw_jackknife = false ;
    checkmap = YES ;
    verbose = NO ;
    num_threads = 1 ;
    use_naive_algo = false ;
    fast_snp_read = false ;
    approx_ld_corr = true ;
    chrom = NULL ;
    nochrom = NULL ;
    print_jackknife_fits = false ;
  }

  void AlderParams::readcommands(int argc, char **argv, const char *VERSION) {
    char *parname = NULL ;
    phandle *ph ;

    int i;
    while ((i = getopt (argc, argv, "p:vV")) != -1) {
      switch (i)
	{
	case 'p':
	  parname = strdup(optarg) ;
	  break;
	case 'v':
	  printf("version: %s\n", VERSION) ; 
	  break; 
	case 'V':
	  verbose = YES ;
	  break; 
	case '?':
	  printf ("Usage: bad params.... \n") ;
	  fatalx("bad params\n") ;
	}
    }
         
    pcheck(parname,'p') ;
    ph = openpars(parname) ;
    dostrsub(ph) ;

    getstring(ph, "genotypename:", &genotypename) ;
    getstring(ph, "snpname:", &snpname) ;
    getstring(ph, "indivname:", &indivname) ;
    getstring(ph, "poplistname:", &poplistname) ;
    getstring(ph, "refpops:", &refpops) ;
    getstring(ph, "weightname:", &weightname) ;
    getstring(ph, "admixpop:", &admixpop) ;
    getstring(ph, "admixlist:", &admixlist) ; // deprecate?
    getstring(ph, "badsnpname:", &badsnpname) ;
    getstring(ph, "raw_outname:", &raw_outname) ;
    int jackknife_int = NO;
    getint(ph, "jackknife:", &jackknife_int) ; print_raw_jackknife = jackknife_int==YES;
    getdbl(ph, "mindis:", &mindis) ;
    getdbl(ph, "maxdis:", &maxdis) ;
    getdbl(ph, "binsize:", &binsize) ;
    getint(ph, "checkmap:", &checkmap) ;
    getint(ph, "mincount:", &mincount) ;
    getint(ph, "num_threads:", &num_threads) ;
    int use_naive_algo_int = NO, fast_snp_read_int = NO, approx_ld_corr_int = YES;
    getint(ph, "use_naive_algo:", &use_naive_algo_int) ; use_naive_algo = use_naive_algo_int==YES;
    getint(ph, "fast_snp_read:", &fast_snp_read_int) ; fast_snp_read = fast_snp_read_int==YES;
    getint(ph, "approx_ld_corr:", &approx_ld_corr_int) ; approx_ld_corr = approx_ld_corr_int==YES;
    getstring(ph, "chrom:", &chrom) ;
    getstring(ph, "nochrom:", &nochrom) ;
    int print_jackknife_fits_int = NO;
    getint(ph, "print_jackknife_fits:", &print_jackknife_fits_int) ; print_jackknife_fits = print_jackknife_fits_int==YES;
    

    check_pars();

    printf("---------- contents of parameter file: %s ----------\n", parname) ;
    writepars(ph) ;
    closepars(ph) ;
    
    print_param_settings();

    if (chrom != NULL) {
      printf("parsing chromosomes to use:");
      chrom_set = parse_to_set(chrom);
    }
    if (nochrom != NULL) {
      printf("parsing chromosomes to ignore:");
      nochrom_set = parse_to_set(nochrom);
    }
    printhline();
  }
}
