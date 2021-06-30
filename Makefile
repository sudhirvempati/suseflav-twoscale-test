# Program : SuSeFLAV: Supersymmetry seesaw and flavor violation calculator
# Version : 1.1.2
# Website : http://cts.iisc.ernet.in/Suseflav/main.html
# Authors : Debtosh Chowdhury  debtosh@cts.iisc.ernet.in 
#           Raghuveer Garani   rgarani@cts.iisc.ernet.in
#           Sudhir Vempati     vempati@cts.iisc.ernet.in
#
# Last Modified: Feb 25, 2012
# 
# Makefile for the program SuSeFLAV


#  FC = ifort
FC = gfortran


src/SuSeFLAV:
	$(MAKE) -C src FC=$(FC)

clean:
	rm -f  *.o *~ */*.o */*~ *.txt */*.txt *.out */*.out

cleanall: 
	rm -f  *.o *~ */*.o *.txt */*.txt *.out */*.out */*~ \
	lib/*.a  bin/suseflav bin/suseflavslha  bin/suseflavscan 
