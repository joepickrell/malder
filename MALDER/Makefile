CXX = g++
CXXOPT = -O2
CXXFLAGS = -fopenmp -Wall -I/opt/local/include -Wno-write-strings $(addprefix -I, ${IDIRS})
L = -L/opt/local/lib -lfftw3 -llapack -lgsl

LOCAL_ADMIXTOOLS_SRC = admixtools_src

# can replace ADMIXDIR with absolute path to the src/ directory of an ADMIXTOOLS install
ADMIXDIR = ${LOCAL_ADMIXTOOLS_SRC}
#ADMIXDIR = /home/pl88/admixtools/admixtools_2/src

IDIRS = ${ADMIXDIR} ${ADMIXDIR}/nicksrc
NLIB = ${ADMIXDIR}/libnick.a # note: Nick copies this to admixlib/nicklib.a
ADMIX_O = $(addprefix ${ADMIXDIR}/,  admutils.o  ldsubs.o  mcio.o  regsubs.o  egsubs.o)

T = malder
O = MalderMain.o Alder.o AlderParams.o CorrJack.o ExpFitALD.o ExpFit.o Jackknife.o MiscUtils.o ProcessInput.o Timer.o MultFitALD.o

.PHONY: libnick.a clean

%.o: %.cpp
	${CXX} ${CXXOPT} ${CXXFLAGS} -o $@ -c $<

$T: libnick.a $O $(ADMIX_O)
	${CXX} ${CXXOPT} ${CXXFLAGS} -o $T $O ${ADMIX_O} ${NLIB} $L

libnick.a:
ifeq ($(ADMIXDIR),$(LOCAL_ADMIXTOOLS_SRC))
	${MAKE} --directory=${ADMIXDIR} libnick.a # recurse and build local stripped-down ADMIXTOOLS library
endif

clean:
	rm -f $O
	rm -f $T
ifeq ($(ADMIXDIR),$(LOCAL_ADMIXTOOLS_SRC))
	${MAKE} --directory=${ADMIXDIR} clean # recurse and clean local stripped-down ADMIXTOOLS library
endif
