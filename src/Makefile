# Program : SuSeFLAV: Supersymmetry seesaw and flavor violation calculator
# Version : 1.1.3
# Website : http://cts.iisc.ernet.in/Suseflav/main.html
# Authors : Debtosh Chowdhury  debtosh@cts.iisc.ernet.in 
#           Raghuveer Garani   rgarani@cts.iisc.ernet.in
#           Sudhir Vempati     vempati@cts.iisc.ernet.in
#
# Last Modified: Feb 29, 2012
# LAPACK dependency removed
# Makefile for the program SuSeFLAV


ifeq ($(FC),ifort)
  FFLAGS = -g -funroll-loops -check pointer -O3
  AR = xiar
else
  FFLAGS = -pg -funroll-loops -fstrength-reduce -O3 
  AR = ar
endif

FUNDOBJS = runrges.o  softspectrum.o  mueconvernew.o  SuSemain.o \
math.o  mssmrge.o  smrge.o  ewsbiterate.o  oneloopPV.o \
oneloopselfenergy.o  onelooptadpole.o  oneloophiggs.o \
oneloopsfermion.o  oneloopfermion.o  oneloopgauge.o \
oneloopneutralino.o  oneloopchargino.o  spectrumtl.o \
oneloopselfenergySM.o  oneloopfermionSM.o  oneloopgaugeSM.o \
rgeiterate.o  slha.o  slha2.o  sinsqthw.o  bsg.o  CEigensystem.o \
HEigensystem.o SVD.o twoloophiggs.o ProtondecaySU5goto.o   sinsqthwSM.o


.SUFFIXES : .o .f .F .a

#oneloopselfenergySM.o oneloopfermionSM.o oneloopgaugeSM.o sinsqthwSM.o \
# let make now that objs depend on inputfiles. make knows what to do then.
# the space before $(FC) MUST BE A TAB!!!

all: SuSeFLAV SuSeFLAV_slha SuSeFLAVscan suseflavlib.a 

SuSeFLAV: suseflavlib.a runonce.o 
	$(FC) -o suseflav  runonce.o ../lib/suseflavlib.a  
	mv -v suseflav ../bin/

SuSeFLAV_slha: suseflavlib.a runslha.o 
	$(FC) -o suseflavslha  runslha.o ../lib/suseflavlib.a
	mv -v suseflavslha ../bin/

SuSeFLAVscan: suseflavlib.a scanning.o 
	$(FC) ../lib/suseflavlib.a scanning.o -o suseflavscan
	mv -v suseflavscan ../bin/

suseflavlib.a:../lib/suseflavlib.a($(FUNDOBJS))

.f.a:
	$(FC) -c $(FFLAGS) $<
	$(AR) -ruc $@ $*.o
	rm -f $*.o

.F.a:
	$(FC) -c $(FFLAGS) $<
	$(AR) -ruc $@ $*.o
	rm -f $*.o

.f.o:
	$(FC) -c  $(FFLAGS) $< 

.F.o:
	$(FC) -c  $(FFLAGS) $< 

clean:
	rm -f  *.o *~ */*.o */*~ *.txt */*.txt *.out */*.out

cleanall: 
	rm -f  *.o *~ */*.o */*~ lib/*.a *.txt */*.txt *.out \
	*/*.out ../bin/suseflav ../bin/suseflavslha ../bin/suseflavscan 
