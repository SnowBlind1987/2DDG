module Conversion_Functions
contains

function calc_uh(x,y,xc,yc,dx,dy,deg_fr,ord_deg,nfr) result(uh)

use Set_Precision
use Legendre_Functions
implicit none
real(prec)::x,y
real(prec)::xc,yc
real(prec)::dx,dy
real(prec)::deg_fr(4,nfr)
integer::ord_deg(2,6)
integer::nfr
real(prec)::uh(4)
integer::k,l

uh(:)=0.0_prec

do k=1,4
   do l=1,nfr
     uh(k)=uh(k)+deg_fr(k,l)*Legendre(x,xc,dx,ord_deg(1,l))*&
           Legendre(y,yc,dy,ord_deg(2,l))
    enddo
enddo

end function calc_uh

subroutine  Euler_Flux(U,Fx,Gy) 
use Set_Precision
use Parameters
implicit none
real(prec)::U(4)
real(prec)::Fx(4)
real(prec)::Gy(4)

real(prec)::rho,vel_x,vel_y,p,rhoH

rho=U(1)
vel_x=U(2)/rho
vel_y=U(3)/rho
p=gm1*(U(4)-0.5_prec*rho*(vel_x**2+vel_y**2))
rhoH=U(4)+p

Fx(1)=rho*vel_x
Fx(2)=rho*vel_x**2+p
Fx(3)=rho*vel_x*vel_y
Fx(4)=rhoH*vel_x

Gy(1)=rho*vel_y
Gy(2)=rho*vel_x*vel_y
Gy(3)=rho*vel_y**2+p
Gy(4)=rhoH*vel_y

end subroutine Euler_Flux

function cons2prim(U) result(w_prim)
use Set_Precision
use Parameters
implicit none
real(prec),intent(in)::U(4)
real(prec)::w_prim(4)
real(prec)::rho,vel_x,vel_y,p

rho=U(1)
vel_x=U(2)/rho
vel_y=U(3)/rho
p=gm1*(U(4)-0.5_prec*rho*(vel_x**2+vel_y**2))

w_prim(1)=rho
w_prim(2)=vel_x
w_prim(3)=vel_y
w_prim(4)=p

end function cons2prim

end module Conversion_Functions


