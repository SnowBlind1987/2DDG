module fluxes
use Set_Precision
contains
!*****************************************************************************
!* A collection of two-dimensional Euler numerical fluxes, Version 2 (2010),
!*
!*        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
!*
!* the author of useful CFD books, "I do like CFD" (http://www.cfdbooks.com).
!*
!*------------------------
!* List of Flux Functions:
!
!   Roe
!   Rotated-RHLL
!
!*------------------------
!*
!* These F90 routines were written and made available for download
!* for an educational purpose. Detailed descripstion of each numerical
!* flux can be found in the original paper, or in popular textbooks, or
!* in the (will-be-available) second volume of "I do like CFD". 
!*
!* Note that all routines are not efficietly implemented for clarity; you can 
!* improve the efficiecy and also you can covert it to double precision 
!* version if you wish.
!*
!* NOTES: There were bugs in the Rotated-RHLL solver in version 1.
!*        It has been fixed and tested for a shock-diffraction problem.
!*
!* This file may be updated in future. (to add more flux functions.)
!*****************************************************************************

!*****************************************************************************
!* -- Roe's Flux Function ---
!*
!* P. L. Roe, Approximate Riemann Solvers, Parameter Vectors and Difference
!* Schemes, Journal of Computational Physics, 43, pp. 357-372.
!* 
!* Katate Masatsuka, February 2009. http://www.cfdbooks.com
!*****************************************************************************
 function Roe(uL, uR, nx, ny)
 real(prec) :: uL(4), uR(4) !  Input: conservative variables rho*[1, u, v, E]
 real(prec) :: nx, ny       !  Input: face normal vector, [nx, ny] (Left-to-Right)
 real(prec) :: Roe(4)       ! Output: Roe flux function (upwind)
!Local constants
 real(prec) :: gamma                          ! Ratio of specific heat.
 real(prec) :: zero, fifth, half, one, two    ! Numbers
!Local variables
 real(prec) :: tx, ty       ! Tangent vector (perpendicular to the face normal)
 real(prec) :: vxL, vxR, vyL, vyR             ! Velocity components.
 real(prec) :: rhoL, rhoR, pL, pR             ! Primitive variables.
 real(prec) :: vnL, vnR, vtL, vtR             ! Normal and tangent velocities
 real(prec) :: aL, aR, HL, HR                 ! Speeds of sound.
 real(prec) :: RT,rho,vx,vy,H,a,vn, vt        ! Roe-averages
 real(prec) :: drho,dvx,dvy,dvn,dvt,dp,dV(4)  ! Wave strenghs
 real(prec):: ws(4),dws(4), Rv(4,4)          ! Wave speeds and right-eigevectors
 real(prec) :: fL(4), fR(4), diss(4)          ! Fluxes ad dissipation term
 integer :: i, j

!Constants.
     gamma = 1.4
      zero = 0.0
     fifth = 0.2
      half = 0.5
       one = 1.0
       two = 2.0

!Tangent vector (Do you like it? Actually, Roe flux can be implemented 
! without any tangent vector. See "I do like CFD, VOL.1" for details.)
  tx = -ny
  ty = nx

!Primitive and other variables.
!  Left state
    rhoL = uL(1)
     vxL = uL(2)/uL(1)
     vyL = uL(3)/uL(1)
     vnL = vxL*nx+vyL*ny
     vtL = vxL*tx+vyL*ty
      pL = (gamma-one)*( uL(4) - half*rhoL*(vxL*vxL+vyL*vyL) )
      aL = sqrt(gamma*pL/rhoL)
      HL = ( uL(4) + pL ) / rhoL
!  Right state
    rhoR = uR(1)
     vxR = uR(2)/uR(1)
     vyR = uR(3)/uR(1)
     vnR = vxR*nx+vyR*ny
     vtR = vxR*tx+vyR*ty
      pR = (gamma-one)*( uR(4) - half*rhoR*(vxR*vxR+vyR*vyR) )
      aR = sqrt(gamma*pR/rhoR)
      HR = ( uR(4) + pR ) / rhoR

!First compute the Roe Averages
    RT = sqrt(rhoR/rhoL)
   rho = RT*rhoL
    vx = (vxL+RT*vxR)/(one+RT)
    vy = (vyL+RT*vyR)/(one+RT)
     H = ( HL+RT* HR)/(one+RT)
     a = sqrt( (gamma-one)*(H-half*(vx*vx+vy*vy)) )
    vn = vx*nx+vy*ny
    vt = vx*tx+vy*ty

!Wave Strengths
   drho = rhoR - rhoL 
     dp =   pR - pL
    dvn =  vnR - vnL
    dvt =  vtR - vtL

  dV(1) = (dp - rho*a*dvn )/(two*a*a)
  dV(2) = rho*dvt/a
  dV(3) =  drho - dp/(a*a)
  dV(4) = (dp + rho*a*dvn )/(two*a*a)

!Wave Speed
  ws(1) = abs(vn-a)
  ws(2) = abs(vn)
  ws(3) = abs(vn)
  ws(4) = abs(vn+a)

!Harten's Entropy Fix JCP(1983), 49, pp357-393:
! only for the nonlinear fields.
  dws(1) = fifth
   if ( ws(1) < dws(1) ) ws(1) = half * ( ws(1)*ws(1)/dws(1)+dws(1) )
  dws(4) = fifth
   if ( ws(4) < dws(4) ) ws(4) = half * ( ws(4)*ws(4)/dws(4)+dws(4) )

!Right Eigenvectors
  Rv(1,1) = one    
  Rv(2,1) = vx - a*nx
  Rv(3,1) = vy - a*ny
  Rv(4,1) =  H - vn*a

  Rv(1,2) = zero
  Rv(2,2) = a*tx
  Rv(3,2) = a*ty
  Rv(4,2) = vt*a

  Rv(1,3) = one
  Rv(2,3) = vx
  Rv(3,3) = vy 
  Rv(4,3) = half*(vx*vx+vy*vy)

  Rv(1,4) = one
  Rv(2,4) = vx + a*nx
  Rv(3,4) = vy + a*ny
  Rv(4,4) =  H + vn*a

!Dissipation Term
  diss = zero
  do i=1,4
   do j=1,4
    diss(i) = diss(i) + ws(j)*dV(j)*Rv(i,j)
   end do
  end do

!Compute the flux.
  fL(1) = rhoL*vnL
  fL(2) = rhoL*vnL * vxL + pL*nx
  fL(3) = rhoL*vnL * vyL + pL*ny
  fL(4) = rhoL*vnL *  HL

  fR(1) = rhoR*vnR
  fR(2) = rhoR*vnR * vxR + pR*nx
  fR(3) = rhoR*vnR * vyR + pR*ny
  fR(4) = rhoR*vnR *  HR

  Roe = half * (fL + fR - diss)

 end function Roe

end module fluxes
