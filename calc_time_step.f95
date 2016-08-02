subroutine calc_time_step(dt)
use Variables
use Parameters
implicit none
real(prec)::dt
real(prec)::rho,vel_x,vel_y,p,spd,c,t
integer::i,j

dt=1.0e20_prec

do j=2,nCells_y+1
   do i=1,nCells_x+1
   rho=cell(i,j)%un(1,1)
   vel_x=cell(i,j)%un(2,1)/rho
   vel_y=cell(i,j)%un(3,1)/rho
   spd=sqrt(vel_x**2+vel_y**2)
   p=gm1*(cell(i,j)%un(4,1)-0.5_prec*rho*(spd**2))
   c=sqrt(gam*p/rho)
   t=cell(i,j)%hx/(spd+c)

   if (t<dt) then
           dt=t
   endif
   enddo
enddo
   dt=CFL*dt
end subroutine calc_time_step
