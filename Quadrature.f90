module Quadrature
use Set_Precision
contains
function Integrate_2D(Mat_G,g_weight,jacob,ng) result(Integral)
implicit none
real(prec),intent(in)::Mat_G(ng,ng)
real(prec),intent(in)::g_weight(5,5)
real(prec),intent(in)::jacob
integer,intent(in)::ng
real(prec)::Integral
integer i,j,k
real(prec)::WJ(ng)
Integral=0.0_prec
do j=1,ng
   WJ(j)=g_weight(ng,j)
enddo

do j=1,ng
   do i=1,ng
      Integral=Integral+WJ(i)*Mat_G(i,j)*WJ(j)
   enddo
enddo

Integral=Integral*jacob
end function Integrate_2D

function Integrate_1D(Vect_G,g_weight,jacob,ng) result(Integral)
implicit none
real(prec),intent(in)::Vect_G(ng)
real(prec),intent(in)::g_weight(5,5)
real(prec),intent(in)::jacob
integer,intent(in)::ng
real(prec)::Integral
real(prec)::WJ(ng)
integer::i
Integral=0.0_prec
do i=1,ng
   Wj(i)=g_weight(ng,i)
enddo

do i=1,ng
   Integral=Integral+WJ(i)*Vect_G(i)
enddo

Integral=Integral*jacob

end function Integrate_1D

end module Quadrature


