
        |
        |      ALDER,   v1.0
     \..|./
    \ \  /       Admixture
     \ |/ /      Linkage
      \| /       Disequilibrium for
       |/        Evolutionary
       |         Relationships
       |

     Po-Ru Loh and Mark Lipson
          October 31, 2012



==== DOCUMENTATION FOR admixtools_src/ ====

This subdirectory (admixtools_src/) of the ALDER release contains
headers and source code implementing I/O routines provided in Nick
Patterson's ADMIXTOOLS software suite:

  http://genetics.med.harvard.edu/reich/Reich_Lab/Software.html

The source files provided here are almost identical to the
corresponding files in the src/ directory of the ADMIXTOOLS 1.0
release, with the following adaptations:

1. Enabling of C++ linkage by wrapping .h files with extern "C" { }.

2. Modification of the mcio.c and mcio.h files to add a function
   (getweights_def, which simply allows specifying a default weight).

If you wish, you may link to an existing installation of ADMIXTOOLS
instead of building the files here.  If you have a 1.0 install, you
can make it compatible (without breaking your ADMIXTOOLS build, of
course) by doing the following:

1. Run the script 'add_extern_C.sh' within the src/ directory of
   your ADMIXTOOLS install.

2. Overwrite mcio.c and mcio.h in src/ with the files provided here.

3. Rebuild ADMIXTOOLS.

4. Set the 'ADMIXDIR' variable in the ALDER base directory Makefile to
   the (absolute) path to your ADMIXTOOLS src/ directory.

5. Build ALDER.
