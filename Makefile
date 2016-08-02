CC=gfortran -g
OBJ= Set_Precision.o\
     Parameters.o\
     Variables.o\
     gauss_init.o\
     Initial_Conditions.o\
     Legendre_Functions.o\
     Quadrature.o\
     init.o\
     Boundary_Conditions.o\
     save_sol.o\
     Fluxes_2D.o\
     Conversion_Functions.o\
     residual.o\
     calc_time_step.o\
     update.o\
     output.o\
     main.o


%.o: %.f95
	$(CC) -c -o $@ $< 
2d_dg: $(OBJ)
	$(CC) -o $@ $^ 
.PHONY: clean

clean:
	rm -f *.o *.mod *~.core
.PHONY: cleanall

cleanall:
	rm -f *.o *.mod *~.core 1d_dg
