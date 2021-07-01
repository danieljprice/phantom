#
#--Makefile to remake the executable
#
FC      = gfortran
FFLAGS  = -O3 -fdefault-real-8  -fdefault-double-8  -fopenmp  -std=f2008 -fall-intrinsics
#FFLAGS  = -fdefault-real-8  -fdefault-double-8 -fopenmp  -std=f2008 -fall-intrinsics -g -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow,denormal -fcheck=all -finit-real=nan -finit-integer=nan -Wunused-parameter
#FFLAGS = 
#FC     = ifort
#FFLAGS = -O3 -inline-factor=500 -mcmodel=medium -shared-intel -warn uninitialized -warn unused -warn truncated_source -no-wrap-margin -r8 -openmp
.SUFFIXES: .F90 

%.o : %.F90
	$(FC) -c $(FFLAGS) $< -o $@

eta: nicil.o  nicil_ex_eta_sup.o nicil_ex_eta.o
	$(FC) $(FFLAGS) -o nicil_ex_eta nicil.o  nicil_ex_eta_sup.o nicil_ex_eta.o
	cp nicil_ex_eta ..

sph: nicil.o  nicil_ex_sph_sup.o nicil_ex_sph.o
	$(FC) $(FFLAGS) -o nicil_ex_sph nicil.o  nicil_ex_sph_sup.o nicil_ex_sph.o
	cp nicil_ex_sph ..

one: nicil.o  nicil_get_one_point.o
	$(FC) $(FFLAGS) -o nicil_get_one_point nicil.o  nicil_get_one_point.o
	cp nicil_get_one_point ..

disc: nicil.o  nicil_ex_eta_sup.o nicil_ex_disc.o
	$(FC) $(FFLAGS) -o nicil_ex_disc nicil.o  nicil_ex_eta_sup.o nicil_ex_disc.o
	cp nicil_ex_disc ..

clean:   
	rm -f *.o *.mod; rm -f plot_lib/*.pyc
	rm nicil_ex_eta; rm nicil_ex_sph; rm nicil_ex_disc; rm nicil_get_one_point; rm ../nicil_ex_eta; rm ../nicil_ex_sph; rm ../nicil_ex_disc; rm ../nicil_get_one_point

allclean:   
	rm -f *.o *.mod; rm -f plot_lib/*.pyc
	rm -f ../data/*.dat; 	rm -f ../data/*.log; rm -f ../disc_param_??.in
	rm -f ../Graphs/plot*.gp; rm -f ../Graphs/*eps; rm -f ../Graphs/*png
	rm nicil_ex_eta; rm nicil_ex_sph; rm nicil_ex_disc; rm nicil_get_one_point; rm ../nicil_ex_eta; rm ../nicil_ex_sph; rm ../nicil_ex_disc; rm ../nicil_get_one_point

