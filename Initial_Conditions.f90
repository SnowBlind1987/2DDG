module Initial_Conditions
use Set_Precision
use Parameters
!***********Common Variables****************
real(prec)::rho_inf,u_inf,v_inf,p_inf
real(prec)::vortex_strength,vortex_center(2)

contains

function Vortex_Convection(x,y) result(u)
implicit none
real(prec),intent(in)::x
real(prec),intent(in)::y
real(prec)::u(4)
real(prec)::rho, vel_x, vel_y, p
real(prec)::du,dv
real(prec)::xc,yc,r
real(prec)::b



u(:)=0.0_prec
b=vortex_strength
xc=vortex_center(1)
yc=vortex_center(2)

r=sqrt ( (x-xc)**2 + (y-yc)**2 )

rho=( 1.0_prec-( gm1*b**2/(8.0_prec*gam*pi**2) )* exp(1-r**2) )**(1.0_prec/gm1)
du=-b/(2.0_prec*pi)*exp((1-r**2)/2.0_prec)*(y-yc)
dv=b/(2.0_prec*pi)*exp((1-r**2)/2.0_prec)*(x-xc)
p=rho*gam

rho=rho_inf*rho
vel_x=u_inf+du
vel_y=v_inf+dv
p=p_inf*p

u(1)=rho
u(2)=rho*vel_x
u(3)=rho*vel_y
u(4)=p/gm1+0.5_prec*rho*( vel_x**2+vel_y**2 )

end function Vortex_Convection

end module Initial_Conditions
