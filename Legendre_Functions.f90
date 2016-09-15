module Legendre_Functions
contains

function Legendre(x,xc,dx,nfunc) result(feval)
use Set_Precision
implicit none
real(prec),intent(in)::x
real(prec),intent(in)::xc
real(prec),intent(in)::dx
integer,intent(in)::nfunc
real(prec)::feval

real(prec)::ksi
feval=0.0_prec
ksi=2.0_prec*(x-xc)/dx

select case(nfunc)

case (1)
        feval=1
case (2)
        feval=ksi
case (3)
        feval=0.5_prec*( 3.0_prec*ksi**2-1.0_prec)
end select

end function Legendre

function dLegendre(x,xc,dx,nfunc) result (dfeval)
use Set_Precision
implicit none
real(prec),intent(in)::x
real(prec),intent(in)::xc
real(prec),intent(in)::dx
integer,intent(in)::nfunc
real(prec)::dfeval

real(prec)::ksi,dksi

dfeval=0.0_prec

ksi=2.0_prec*(x-xc)/dx
dksi=2.0_prec/dx

select case (nfunc)

case (1)
        dfeval=0.0_prec
case (2)
        dfeval=dksi
case (3)
        dfeval=3.0_prec*ksi*dksi
end select

end function dLegendre

end module Legendre_Functions
