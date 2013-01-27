CC=gcc
HOMEL=.
DEBUG_OPTIONS= # -g -pg

#NLIB=$(HOMEL)/admixlib/nicklib.a
NLIB=$(HOMEL)/libnick.a

IDIRS=${HOMEL} ${HOMEL}/nicksrc
VPATH=.:nicksrc

CFLAGS= -c -O2  $(addprefix -I, ${IDIRS}) -Wimplicit ## -fprofile-arcs -ftest-coverage 
OBJ_LIBNICK=strsubs.o sortit.o vsubs.o statsubs.o linsubs.o getpars.o xsearch.o gauss.o	gds.o
OBJ_HOMEL=admutils.o  ldsubs.o  mcio.o  regsubs.o  egsubs.o
TWTAB=\"$(HOMEL)/twtable\"

libnick.a:	$(OBJ_LIBNICK) $(OBJ_HOMEL)
	ar -r libnick.a $(OBJ_LIBNICK)
#	mkdir -p $(HOMEL)/admixlib
#	cp libnick.a $(NLIB)

statsubs.o:     nicksrc/statsubs.c
	$(CC)  $(CFLAGS) -DTWTAB=$(TWTAB) -o statsubs.o nicksrc/statsubs.c

clean: 
	rm -f *.o 
	rm -f nicksrc/*.o
	rm -f libnick.a
#	rm -f ${NLIB}
