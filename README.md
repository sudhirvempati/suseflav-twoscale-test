# SuSeFLAV

State of the art computational methods are essential to completely understand Supersymmetry. SuSeFLAV is one such numerical tool which is capable of investigating mSUGRA, GMSB, non-universal Higgs models and complete non-universal models. The program solves complete MSSM RGEs with complete 3 flavor mixing at 2-loop level and also adds one loop threshold corrections to all MSSM parameters by incorporating radiative electroweak symmetry breaking conditions, using standard model fermion masses and gauge couplings as inputs at the weak scale. The program has a provision to run three right handed neutrinos at user defined scales and mixing. Also, the program computes branching ratios and decay rates for various flavor violating processes such as μ → e γ, τ → e γ, τ → μ γ, μ → e e e, τ → μ μ μ, τ → e e e, b → s γ etc. and anomalous magnetic moment of muon.

# Compilation and Installation Instructions for SuSeFLAV-1.2.0

*NOTE: From this version, SUSEFLAV is independent of LAPACK*


Download the tar.gz file and unpack it in your `home` directory
using:

```tar -zxvf suseflav_1.x.x.tar.gz ```

This creates a directory in your `home` with name `suseflav_1.x.x`

To install SuSeFLAV you will need FORTRAN 90/95 compiler. SuSeFLAV has been
successfully compiled using GNU g95 and Intel ifort(12.x) on Linux and Unix
operating systems [1]. The compilation of the program is handled by the provided
Makefile. To make and install the distribution type the command 

[1]: From this version, SUSEFLAV is independent of LAPACK

` make `

This will compile SuSeFLAV with default FORTRAN compiler (gfortran) and
install the package in the `bin` sub-directory from the main directory
you compile it in. The `make` will build the library `libsuseflav.a`.

_Note! You might have problems in compiling if you are using some versions of g77_

---
##### **_Important Note_**

Some times, it might be required that you have to run make twice for the library to be created. Especially, if you have made some large changes in the original programs, please remove the library in `lib` directory and remake the library, by running make twice.

*Some compilers, especially gfortran also tend to give errors like*
`ld : symbol_xxx_ has length NNNN in file1.o whereas it has length NNNN in file2.o`

Please ignore all such warning signs, but they are compiler dependent.

----

To clean the object files and other output files type

`    make clean `

And finally, to uninstall SuSeFLAV (i.e. delete all compiled libraries and
files, but keep the sources), type

`  make cleanall `


## Running SuSeFLAV


SuSeFLAV package produces three executable files when compiled,
namely `suseflav`, `suseflavslha` and `suseflavscan`. To compute
the spectrum for a single point the usage of executable files
`suseflav` and `suseflavslha` is recommended. To scan the
parameter space the usage  of the executable `suseflavscan`
is recommended.

To run the program, go the `bin` directory and :

1. For a single point

i) `  ./suseflav  <filename `

example: filename = `sinputs.in` or `sinputs-gmsb.in` or `sinputs-nuhm.in`

The main file for the executable is `runonce.f`.

This executable takes the following input files in traditional
suseflav format,

`sinputs.in`       for msugra models
`sinputs-gmsb.in`  for gmsb models
`sinputs-nuhm.in`  for nuhm models


Example: To run mSUGRA

`./suseflav <sinputs.in`

The output is saved in `suseflav.out` as well as appearing in standard I/O.
The set of supersymmetric inputs and final observables are saved in `tmp/output.txt`.

The format for the output in the text file `output.txt` is :

### For mSUGRA

`tanbeta, m0, m12, a0, sgn(mu), mh(light neutral higgs mass), g_mu-2, Br(b->s+gamma),
Br(mu->e+gamma),Br(tau->mu+gamma),Br(tau ->e+gamma),Br(mu->3e), Br(tau->3mu),
Br(tau->3e),flags`

### For GMSB

`tanbeta,Lambda,Mmessenger,NMess,sgn(mu),mh(light neutral higgs mass), g_\mu-2,
Br(b->s+gamma),Br(mu->e+gamma) Br(tau->mu+gamma),Br(tau ->e+gamma),Br(mu->3e),
Br(tau->3mu),Br(tau->3e),flags`

### For NUHM

`tanbeta,m0, m12,a0,sgn(mu),mh10,mh20,mh(light neutral higgs mass),g_\mu-2,
Br(b->s+gamma),Br(mu->e+gamma), Br(tau->mu+gamma),Br(tau ->e+gamma),
Br(mu->3e), Br(tau->3mu),Br(tau->3e),flags`

ii) `./suseflavslha`

The main file for the executable is `runslha.f`. The input file for this
executable is slha.in, which explores mSUGRA. Abundant examples are provided
in the example folder. The user must rename the required slha file as
`slha.in` to use that particular file as input.
Example: To run mSUGRA

`./suseflavslha`

The output is saved in `slha.out`


2. Scanning parameter space

`./suseflavscan`

Main file for this executable is `scanning.f` and the corresponding input file
`sinputs_scan.in` (Inputs in traditional suseflav format). The main file
utilizes random number generator to assign values to input variables
and uses system call to run the executable `suseflav`.

The output of scanning run is saved in `scan.out`, where essentially the
main supersymmetry breaking input parameters and observables are listed.
The format is the same as in for single point run as mentioned above.
The user can change the parameters she/he wants to write in the
file `scanning.f`.

The example given has been run for 1000 points.
The user time is on an average 10m50secs. And the CPU time is
about 0m5secs on an average. So, approximately about 11-13 mins
on a Core 2 duo (*imagine a three year old laptop!!*) for about
1000 points. It is much faster on i5 processors with about
40 mins for 10000 points. On Intel sandy bridge 2.6 GHz, we
got speeds of less than 4 hours for 100000 points.
The speed also depends on the number of parameters one is varying.


## Precision of the spectrum.


In addition to the parameter spectrum tolerance in the input files (which
decides the accuracy at which the spectrum converges), precision of the
spectrum also crucially depends on the precision at which the RGE's are
evolved. This is set by three parameters `h1`, `hmin` and `eps` in the `runrges.f`
subroutine in `/src` folder. Users are advised to set
the precision according the problem at hand.

----------------------------------------------

EXAMPLES

Folder `examples` contains a number of example input files.

To run the example files, copy the file you want to run in to the
main directory and couple it with relevant executable as described
above.

----------------------------------------------
FILE STRUCTURE

``` Input files:
	1. slha.in
	2. sinputs.in
	3. sinputs-gmsb.in
	4. sinputs-nuhm.in
	5. sinputs-scan.in
	6. sinputs-cnum.in
	7. sinputs-nugm.in
	8. sinputs-rpar.in
	9. sinputs_scan.in
```

Several example inputs are contained in the folder `examples`.
Most of them reflect the benchmark points after the recent
LHC data.

Source files:  
       The directory `src` contains all the source files required
       by the program

Output files:  
       1. `suseflav.out` for the output in unformatted text.  
       2. `slha.out` for output in slha format.



## Documentation

If you use SuSeFLAV in your work please cite [D. Chowdhury et al., Comput. Phys. Commun. 184 (2013) 899](http://dx.doi.org/10.1016/j.cpc.2012.10.031), [[arXiv:1109.3551]](https://arxiv.org/abs/1109.3551). It will be regularly updated on arXiv and will serve as user manual.

## Web

<http://chep.iisc.ac.in/Suseflav/main.html>  
<https://suseflav.hepforge.org>

## Contact

Debtosh Chowdhury	&nbsp;		`Debtosh.Chowdhury AT polytechnique.edu`  
Raghuveer Garani	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;		`veergarani AT gmail.com`  
Sudhir K. Vempati	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;	`vempati AT iisc.ac.in`
