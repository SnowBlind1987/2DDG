program main
use Set_Precision
use Variables
use Initial_Conditions
implicit none
real(prec)::time,dt
integer::iter,rk,k,i

xmin=0.0_prec
xmax=10.0_prec
ymin=0.0_prec
ymax=10.0_prec
nCells_x=60
nCells_y=60

final_time=2.0_prec
cfl=0.5_prec
pord=3

!Initial Conditions Inputs

rho_inf=1.0_prec
u_inf=0.5_prec
v_inf=0.0_prec
p_inf=1.0_prec
vortex_strength=0.5
vortex_center(1)=5.0
vortex_center(2)=5.0
call gauss_init()
call init()
call output(0,0.0_prec)
cfl=cfl/(2.0_prec*pord)
time=0.0_prec
dt=0.0_prec
iter=0
do while(time<final_time)
   call save_sol()
   call calc_time_step(dt)

   if (time+dt>final_time)then
           dt=final_time-time
   endif

   do rk=1,3
      call Periodic_BC()
      call residual()
      call update(rk,dt)
   enddo
   time=time+dt
   iter=iter+1
   call output(iter,time)
   print *,time
enddo
end program main
