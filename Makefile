all: jewel-2.2.0-med jewel-2.2.0-simple jewel-2.2.0-vac jewel-2.2.0-reader

# path to LHAPDF library
LHAPDF_PATH := /sampa/fcanedo/lhapdf/lib
#LHAPDF_PATH :=/cvmfs/alice.cern.ch/el6-x86_64/Packages/lhapdf5/v5.9.1-alice1-2/lib
#LHAPDF_PATH := /sampa/alice/lhapdf/lib

FC := gfortran
FFLAGS := -g -static -Wno-extra

jewel-2.2.0-med: mymed.o medium-simple.o pythia6425mod.o meix.o
	$(FC) -o $@ -L$(LHAPDF_PATH) $^ -lLHAPDF

jewel-2.2.0-vac: jewel-2.2.0.o medium-vac.o pythia6425mod.o meix.o
	$(FC) -o $@ -L$(LHAPDF_PATH) $^ -lLHAPDF

jewel-2.2.0-simple: jewel-2.2.0.o medium-simple.o pythia6425mod.o meix.o
	$(FC) -o $@ -L$(LHAPDF_PATH) $^ -lLHAPDF

jewel-2.2.0-reader: jewel-2.2.0.o medium-reader.o reader.o pythia6425mod.o meix.o
	$(FC) -o $@ -L$(LHAPDF_PATH) $^ -lLHAPDF

clean:
	rm -f *.o 
	rm -f medium-*.o 
	rm -f mymed.o
	rm -f pythia6425mod.o meix.o
	rm -f *~

.PHONY: all
