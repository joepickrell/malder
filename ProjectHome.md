MALDER. This is a version of ALDER (http://groups.csail.mit.edu/cb/alder/) that has been modified to allow multiple admixture events.

## Installation (unix only) ##

svn checkout http://malder.googlecode.com/svn/MALDER malder-read-only

cd malder-read-only/

make

## Usage ##

malder -p parfile

MALDER is a modified version of ALDER (http://groups.csail.mit.edu/cb/alder/), the parameter file should be identical to ALDER in its "multiple admixture test" mode

## Output ##

The output will be printed to the screen, with lines starting with RESULT\_X being the results from fitting X curves to the data. The program will always fit one more curve than is necessary, so if you get RESULT\_3, then the lines that start with RESULT\_2 are the ones you're interested in.

