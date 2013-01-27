#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <set>

#include "mcio.h"
#include "egsubs.h"

#include "ProcessInput.hpp"

namespace ProcessInput {

  using std::cout;
  using std::endl;
  using std::flush;
  using std::string;
  using std::vector;
  using std::pair;
  using std::make_pair;
  using std::set;
  using std::max;
  using std::count;

  const int ADMIXED_POP_IND = 9999; // for use in process_indivs, process_geno

  int cmap(SNP **snpmarkers, int numsnps) {
    int t, k ; 
    double y1, y2 ; 
    SNP *cupt ;  
    for (k=1 ; k<=10; ++k) { 
      t = ranmod(numsnps) ;
      cupt = snpmarkers[t] ; 
      y1 = cupt -> genpos ; 
      y2 = cupt -> physpos / 1.0e8  ; 
      if (fabs(y1-y2) > .001) return YES ;
    }
    return NO ;
  }

  // returns locations (chrom, genpos) of valid (i.e., non-ignore) snps
  vector < pair <int, double> > process_snps(char *snpname, char *badsnpname, bool fast_snp_read,
					     SNP ***snpmarkers_ptr, int checkmap, int &numsnps,
					     const set <int> &chrom_set,
					     const set <int> &nochrom_set) {
    int nignore = 0;
    if (fast_snp_read) {
      numsnps = 0;
      FILE *snp_file = fopen(snpname, "r");
      if (snp_file == NULL) fatalx("unable to open geno file\n");
      const int buf_size = 1024;
      char line[buf_size];
      while (fgets(line, buf_size, snp_file) != NULL)
	numsnps++;
      ZALLOC(*snpmarkers_ptr, numsnps, SNP *);
      ZALLOC((*snpmarkers_ptr)[0], numsnps, SNP);
      fseek(snp_file, 0, SEEK_SET);
      for (int s = 0; s < numsnps; s++) {
	SNP * &cupt = (*snpmarkers_ptr)[s]; // &cupt: need to allocate the actual pointer
	cupt = (*snpmarkers_ptr)[0] + s;
	fgets(line, buf_size, snp_file);
	sscanf(line, "%s%d%lf%lfd", cupt->ID, &cupt->chrom, &cupt->genpos, &cupt->physpos);
      }
      fclose(snp_file);
    }
    else
      numsnps = getsnps(snpname, snpmarkers_ptr, 0.0, badsnpname, &nignore, 0) ;

    if (checkmap && (cmap(*snpmarkers_ptr, numsnps) == NO)) { 
      printf("Error: genetic map from snp file looks fake.  "
	     "If you mean to do this set checkmap: NO\n") ;
      fatalx("no real map\n") ;
    }

    vector < pair <int, double> > snp_locs;
    for (int i=0; i<numsnps; ++i) {
      SNP *cupt = (*snpmarkers_ptr)[i] ;
      //cupt -> tagnumber = i ;
      int chrom = cupt -> chrom ; 
    
      if (chrom<1) cupt -> ignore = YES ;
      if (chrom>22) cupt -> ignore = YES ;
      //if ((xchrom>0) && (chrom != xchrom)) cupt -> ignore = YES ;
      //if (chrom == xnochrom) cupt -> ignore = YES ;
      if (!chrom_set.empty() && !chrom_set.count(chrom)) cupt -> ignore = YES;
      if (nochrom_set.count(chrom)) cupt -> ignore = YES;
      if (cupt -> ignore == NO) snp_locs.push_back(make_pair(chrom, cupt->genpos));
    }
    printf("using %d of %d snps in data set\n", (int) snp_locs.size(), numsnps);
    return snp_locs;
  }

  // returns map from indices to groups (0, 1, 2, ... for refs, ADMIXPOP_ID for mixed, -1 o/w)
  vector <int> process_indivs(char *indivname, Indiv ***indivmarkers_ptr, char *admixlist,
			      char *admixpop, char *refpops, char *poplistname,
			      int &num_mixed_indivs, string &mixed_pop_name,
			      vector <int> &num_ref_indivs, vector <string> &ref_pop_names) {
  
    const int MAXPOPS = 1000; // max pops for buffer allocation

    num_ref_indivs.clear(); ref_pop_names.clear();
    int numindivs = getindivs(indivname, indivmarkers_ptr) ;
    cout << "num indivs in full data set: " << numindivs << endl;

    // map indiv indices to groups

    // - get admixed population name, either from list or direct param

    char **admixpoplist;
    ZALLOC(admixpoplist, MAXPOPS, char *) ;  
    int nadmixpop = 0;
    if (admixlist != NULL) {
      printf("reading admixlist file: %s\n", admixlist);
      nadmixpop = loadlist(admixpoplist, admixlist) ;
    }
    if (nadmixpop == 0) {  
      if (admixpop == NULL) fatalx("no admixpop\n") ; 
      nadmixpop = 1 ;
      admixpoplist[0] = strdup(admixpop) ;
    }
    if (nadmixpop != 1) fatalx("only one admixpop allowed\n");
    mixed_pop_name = string(admixpoplist[0]);
  
    // - get reference pop names (if provided, as opposed to explicit weights)

    char **refpoplist;
    ZALLOC(refpoplist, MAXPOPS, char *) ; 
    int num_refs = 0 ;
    if (poplistname != NULL) { // reference pops provided in file
      printf("reading poplist file: %s\n", poplistname);
      num_refs = loadlist(refpoplist, poplistname) ;
    }
    else if (refpops != NULL) { // reference pops provided in list
      printf("parsing refpops string (semicolon-delimited; no spaces): %s\n", refpops);
      refpoplist[num_refs++] = &refpops[0];
      int refpops_strlen = strlen(refpops); // this will change!
      for (int i = 0; i < refpops_strlen; i++) {
	if (refpops[i] == ';') {
	  refpops[i] = '\0';
	  refpoplist[num_refs++] = &refpops[i+1];
	}
      }
    }

    // get rid of admixpop if it's in the list of refs
    int num_refs_save = num_refs; num_refs = 0;
    for (int r = 0; r < num_refs_save; r++) {
      if (strcmp(refpoplist[r], admixpoplist[0]) == 0)
	printf("WARNING: admixed pop is in list of reference pops; eliminating\n");
      else
	refpoplist[num_refs++] = refpoplist[r];
    }
  
    for (int r = 0; r < num_refs; r++)
      ref_pop_names.push_back(string(refpoplist[r]));

    vector <int> indiv_pop_inds(numindivs, -1);

    for (int i=0; i<numindivs; ++i) { 
      Indiv *indx = (*indivmarkers_ptr)[i] ;
      if (indx->ignore == YES) continue;
      int ind_reflist = indxindex(refpoplist, num_refs, indx -> egroup) ;
      int ind_admixlist = indxindex(admixpoplist, nadmixpop, indx -> egroup) ;

      if (ind_reflist >= 0) indiv_pop_inds[i] = ind_reflist;
      if (ind_admixlist >= 0) indiv_pop_inds[i] = ADMIXED_POP_IND;
      //indx -> affstatus = indiv_pop_inds[i] ;
      //if (indx -> affstatus < 0) indx -> ignore = YES ;
    }
  
    num_mixed_indivs = count(indiv_pop_inds.begin(), indiv_pop_inds.end(), ADMIXED_POP_IND);
    printf("  admixed pop: %s (%d indivs)\n", admixpoplist[0], num_mixed_indivs);
    if (num_mixed_indivs == 0)
      fatalx("population %s has no indivs... check spelling?\n", admixpoplist[0]);
    for (int k=0; k<num_refs; ++k) {
      num_ref_indivs.push_back(count(indiv_pop_inds.begin(), indiv_pop_inds.end(), k));
      printf("  ref pop %d: %s (%d indivs)\n", k, refpoplist[k], num_ref_indivs.back());
      if (num_ref_indivs.back() == 0)
	fatalx("population %s has no indivs... check spelling?\n", refpoplist[k]);
    }

    // todo: the admixpop and refpop lists really should be freed... memory leak

    return indiv_pop_inds;
  }

  // fills mixed_geno with valid_snps x num_mixed_indivs 0129-array
  // returns valid_snps x num_refs array of reference allele freqs
  // note: valid_snps includes those in chrom 1-22 and not badsnp file
  // it's imporant that mixed_geno and ref_genos are passed by value (copied) here!
  vector < vector <double> > process_geno(char *genotypename, const vector <int> &indiv_pop_inds,
					  char *mixed_geno, vector <char *> ref_genos,
					  SNP **snpmarkers, int numsnps) {
    int num_refs = 0, numindivs = indiv_pop_inds.size();
    vector <int> mixed_indivs;
    for (int i = 0; i < numindivs; i++) {
      if (indiv_pop_inds[i] == ADMIXED_POP_IND) mixed_indivs.push_back(i);
      else num_refs = max(num_refs, indiv_pop_inds[i]+1);
    }
    vector < vector <int> > ref_indivs(num_refs);
    for (int i = 0; i < numindivs; i++)
      if (0 <= indiv_pop_inds[i] && indiv_pop_inds[i] < num_refs)
	ref_indivs[indiv_pop_inds[i]].push_back(i);

    FILE *geno_file = fopen(genotypename, "r");
    if (geno_file == NULL) fatalx("unable to open geno file\n");

    // don't initialize; push_back because of ignored snps
    vector < vector <double> > ref_freqs(num_refs);
    cout << "reading genotype data" << flush;
    char line[numindivs+10];
    for (int s = 0; s < numsnps; s++) {
      if ((s & 0x3fff) == 0)
	cout << "." << flush;
      if (fgets(line, numindivs+10, geno_file) == NULL)
	fatalx("premature EOF (expected %d snps)\n", numsnps);
      if ((int) strlen(line) != numindivs+1)
	fatalx("geno file line has wrong length: expected %d, got %d\n",
	       numindivs, strlen(line)-1);

      if (snpmarkers[s]->ignore == NO) { // skip snps flagged to ignore in process_snps
	// admixed entries: add to mixed_geno buffer
	for (int j = 0; j < (int) mixed_indivs.size(); j++)
	  *mixed_geno++ = line[mixed_indivs[j]]-'0';
	// ref pops: compute freqs
	for (int r = 0; r < num_refs; r++) {
	  int gtype_tot = 0, gtype_ctr = 0;
	  for (int j = 0; j < (int) ref_indivs[r].size(); j++) {
	    int gtype = *ref_genos[r]++ = line[ref_indivs[r][j]]-'0';
	    if (gtype != 9) {
	      gtype_tot += gtype;
	      gtype_ctr++;
	    }
	  }
	  ref_freqs[r].push_back(0.5 * gtype_tot / gtype_ctr); // can be nan
	}
      }
    }
    cout << " done" << endl;

    if (fgetc(geno_file) != EOF)
      fatalx("expected EOF after %d snps, but file still has data\n", numsnps);

    fclose(geno_file);
    return ref_freqs;
  }

  vector <double> process_weights(char *weightname, SNP **snpmarkers, int numsnps) {
    vector <double> weights;
    printf("loading weights from file: %s\n", weightname) ; 
    printf("num weights set: %d\n", getweights_def(weightname, snpmarkers, numsnps, NAN)) ;
    // need to get rid of invalid snps
    for (int s = 0; s < numsnps; s++)
      if (snpmarkers[s]->ignore == NO)
	weights.push_back(snpmarkers[s]->weight);
    return weights;
  }

}
