cmake_minimum_required (VERSION 2.8)
project(2dDG)

enable_language(Fortran)
add_executable(2dDG
     Set_Precision.f90 Parameters.f90 Variables.f90 gauss_init.f90 Initial_Conditions.f90 Legendre_Functions.f90
     Quadrature.f90 init.f90 Boundary_Conditions.f90 save_sol.f90 Fluxes_2D.f90 Conversion_Functions.f90
     residual.f90 calc_time_step.f90 update.f90 output.f90 main.f90)


