#include <string>
#include <vector>
#include <utility>
#include <set>

#include "mcio.h"

namespace ProcessInput {

  using std::string;
  using std::vector;
  using std::pair;
  using std::set;

  int cmap(SNP **snpmarkers, int numsnps);

  // returns locations (chrom, genpos) of valid (i.e., non-ignore) snps
  vector < pair <int, double> > process_snps(char *snpname, char *badsnpname, bool fast_snp_read,
					     SNP ***snpmarkers_ptr, int checkmap, int &numsnps,
					     const set <int> &chrom_set,
					     const set <int> &nochrom_set);

  // returns map from indices to groups (0, 1, 2, ... for refs, ADMIXPOP_ID for mixed, -1 o/w)
  vector <int> process_indivs(char *indivname, Indiv ***indivmarkers_ptr, char *admixlist,
			      char *admixpop, char *refpops, char *poplistname,
			      int &num_mixed_indivs, string &mixed_pop_name,
			      vector <int> &num_ref_indivs, vector <string> &ref_pop_names);
  
  // fills mixed_geno with valid_snps x num_mixed_indivs 0129-array
  // returns valid_snps x num_refs array of reference allele freqs
  // note: valid_snps includes those in chrom 1-22 and not badsnp file
  // it's imporant that mixed_geno and ref_genos are passed by value (copied) here!
  vector < vector <double> > process_geno(char *genotypename, const vector <int> &indiv_pop_inds,
					  char *mixed_geno, vector <char *> ref_genos,
					  SNP **snpmarkers, int numsnps);
  vector <double> process_weights(char *weightname, SNP **snpmarkers, int numsnps);

}
