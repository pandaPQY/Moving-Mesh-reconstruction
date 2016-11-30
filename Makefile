

SOURCE= dimen.f90 gauss1.f90 main.f90 mg1024.f90 mg128.f90 mg16.f90 mg256.f90 mg32.f90 mg512.f90 mg64.f90 mg8.f90 mgutil.f90 multigrid.f90 relaxgl.f90 relaxing.f90 limiter.f90 mytools.f90

FPP= -DCOLD -DOMPTIME

OPT= -qopenmp 


FC=ifort -O3 -mcmodel=large

run:gauss1.o main.o mgutil.o multigrid.o recipes12.o gradient.o reconstruct.o relaxing.o limiter.o mytools.o
	$(FC) $^ -o $@ $(FPP) $(OPT)

gauss1.o:gauss1.f90 dimen.f90 relaxgl.f90
	$(FC)  -c -fpp $< $(FPP) $(OPT)

main.o:main.f90 dimen.f90
	$(FC) -c -fpp $< $(FPP) $(OPT)

mgutil.o:mgutil.f90
	$(FC) -c -fpp $< $(FPP) $(OPT)

relaxing.o:relaxing.f90 relaxgl.f90 cold.f90
	$(FC) -c -fpp $< $(FPP) $(OPT)

limiter.o:limiter.f90 relaxgl.f90 dimen.f90 globalpa.f90
	$(FC) -c -fpp $< $(FPP) $(OPT)

multigrid.o:multigrid.f90 dimen.f90 mg8.f90 mg16.f90 mg32.f90 mg64.f90 mg128.f90 mg256.f90 mg512.f90 mg1024.f90
	$(FC) -c -fpp $< $(FPP) $(OPT)

recipes12.o:recipes12.f90
	$(FC) -c -fpp $< $(FPP) $(OPT)

gradient.o:gradient.f90
	$(FC) -c -fpp $< $(FPP) $(OPT)

reconstruct.o:reconstruct.f90
	$(FC) -c -fpp $< $(FPP) $(OPT)

mytools.o:mytools.f90
	$(FC) -c -fpp $< $(FPP) $(OPT)

clean:
	rm -rf *.o run
