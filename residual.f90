subroutine residual()
use Variables
use Fluxes
use Conversion_Functions
use Quadrature
use Legendre_Functions
implicit none
real(prec)::cons_coord
real(prec)::line_coord(pord)
real(prec)::xcL,ycL,xcR,ycR,hxL,hyL,hxR,hyR
real(prec)::nx,ny
real(prec)::UL(4),UR(4)
real(prec),allocatable::face_fluxes(:,:)
real(prec),allocatable::gauss_Fx(:,:,:),gauss_Gy(:,:,:)
real(prec)::FX(pord)
real(prec)::line_jacob,c_jacob
real(prec)::FXY(pord,pord)
real(prec)::dphix,dphiy
real(prec)::ctr_x,ctr_y,cdx,cdy
real(prec)::gauss_x(pord),gauss_y(pord)
real(prec)::wrk_uh(4)
real(prec)::Integral
integer::i,j,k,l,q,r

allocate(face_fluxes(1:pord,1:4))
allocate(gauss_Fx(1:pord,1:pord,1:4))
allocate(gauss_Gy(1:pord,1:pord,1:4))



UL(:)=0.0_prec
UR(:)=0.0_prec
wrk_uh(:)=0.0_prec
FX(:)=0.0_prec
FXY(:,:)=0.0_prec
face_fluxes(:,:)=0.0_prec
gauss_Fx(:,:,:)=0.0_prec
gauss_Gy(:,:,:)=0.0_prec
gauss_x(:)=0.0_prec
gauss_y(:)=0.0_prec
Integral=0.0_prec
do j=1,nCells_y+2
   do i=1,nCells_x+2
      do l=1,n_dg_fr
         do k=1,4
            cell(i,j)%res(k,l)=0.0_prec
         enddo
      enddo
   enddo
enddo


do j=2,nCells_y+1
   do i=2,nCells_x+1
   !**********Face One (Right Vertical)*****************
   cons_coord=cell(i,j)%face_const_coord(1)
   line_coord=cell(i,j)%gauss_pts_y(1:pord)
   xcL=cell(i,j)%xc
   ycL=cell(i,j)%yc
   
   xcR=cell(i+1,j)%xc
   ycR=cell(i+1,j)%yc
   
   hxL=cell(i,j)%hx
   hyL=cell(i,j)%hy

   hxR=cell(i+1,j)%hx
   hyR=cell(i+1,j)%hy

   nx=face_norm(1,1)
   ny=face_norm(1,2)


   line_jacob=0.5_prec*hyL
   do q=1,pord
      UL=calc_uh(cons_coord,line_coord(q),xcL,ycL,hxL,hyL,cell(i,j)%un(:,1:n_dg_fr),Aord,n_dg_fr)
      UR=calc_uh(cons_coord,line_coord(q),xcR,ycR,hxR,hyR,cell(i+1,j)%un(:,1:n_dg_fr),Aord,n_dg_fr)
      face_fluxes(q,:)=Roe(UL,UR,nx,ny)
   enddo
   
   do l=1,n_dg_fr
      do k=1,4
         do q=1,pord
            FX(q)=face_fluxes(q,k)*Legendre(cons_coord,xcL,hxL,Aord(1,l))*Legendre(line_coord(q),ycL,hyL,Aord(2,l))
         enddo
         Integral=Integrate_1D(FX,wg,line_jacob,pord)
         cell(i,j)%res(k,l)=cell(i,j)%res(k,l)+Word(l)/(hxL*hyL)*Integral
      enddo
    enddo
    !*************End Face One*******************************

    !*************Beging Face Two (Top Horizontal)***********
    cons_coord=cell(i,j)%face_const_coord(2)
    line_coord=cell(i,j)%gauss_pts_x(1:pord)
    xcL=cell(i,j)%xc
    ycL=cell(i,j)%yc
   
    xcR=cell(i,j+1)%xc
    ycR=cell(i,j+1)%yc
   
    hxL=cell(i,j)%hx
    hyL=cell(i,j)%hy

    hxR=cell(i,j+1)%hx
    hyR=cell(i,j+1)%hy

    nx=face_norm(2,1)
    ny=face_norm(2,2)

    line_jacob=0.5_prec*hxL
    do q=1,pord
       UL=calc_uh(line_coord(q),cons_coord,xcL,ycL,hxL,hyL,cell(i,j)%un(:,1:n_dg_fr),Aord,n_dg_fr)
       UR=calc_uh(line_coord(q),cons_coord,xcR,ycR,hxR,hyR,cell(i,j+1)%un(:,1:n_dg_fr),Aord,n_dg_fr)
       face_fluxes(q,:)=Roe(UL,UR,nx,ny)
    enddo
   
    do l=1,n_dg_fr
       do k=1,4
          do q=1,pord
             FX(q)=face_fluxes(q,k)*Legendre(line_coord(q),xcL,hxL,Aord(1,l))*Legendre(cons_coord,ycL,hyL,Aord(2,l))
          enddo
          Integral=Integrate_1D(FX,wg,line_jacob,pord)
          cell(i,j)%res(k,l)=cell(i,j)%res(k,l)+Word(l)/(hxL*hyL)*Integral
       enddo
     enddo
     !***************End Face Tw0****************************************
     !*******************************************************************
     !***************Beging Face Three (Left Vertical)*******************
    cons_coord=cell(i,j)%face_const_coord(3)
    line_coord=cell(i,j)%gauss_pts_y(1:pord)
    xcL=cell(i,j)%xc
    ycL=cell(i,j)%yc
   
    xcR=cell(i-1,j)%xc
    ycR=cell(i-1,j)%yc
   
    hxL=cell(i,j)%hx
    hyL=cell(i,j)%hy

    hxR=cell(i-1,j)%hx
    hyR=cell(i-1,j)%hy

    nx=face_norm(3,1)
    ny=face_norm(3,2)

    line_jacob=0.5_prec*hyL
    do q=1,pord
       UL=calc_uh(cons_coord,line_coord(q),xcL,ycL,hxL,hyL,cell(i,j)%un(:,1:n_dg_fr),Aord,n_dg_fr)
       UR=calc_uh(cons_coord,line_coord(q),xcR,ycR,hxR,hyR,cell(i-1,j)%un(:,1:n_dg_fr),Aord,n_dg_fr)
       face_fluxes(q,:)=Roe(UL,UR,nx,ny)
    enddo
   
    do l=1,n_dg_fr
       do k=1,4
          do q=1,pord
             FX(q)=face_fluxes(q,k)*Legendre(cons_coord,xcL,hxL,Aord(1,l))*Legendre(line_coord(q),ycL,hyL,Aord(2,l))
          enddo
          Integral=Integrate_1D(FX,wg,line_jacob,pord)
          cell(i,j)%res(k,l)=cell(i,j)%res(k,l)+Word(l)/(hxL*hyL)*Integral
       enddo
     enddo
     !***********************End Face Three***********************************
     !************************************************************************
     !********************Beging Face Four (Bottom Horizontal)****************
     cons_coord=cell(i,j)%face_const_coord(4)
     line_coord=cell(i,j)%gauss_pts_x(1:pord)
     xcL=cell(i,j)%xc
     ycL=cell(i,j)%yc
   
     xcR=cell(i,j-1)%xc
     ycR=cell(i,j-1)%yc
   
     hxL=cell(i,j)%hx
     hyL=cell(i,j)%hy

     hxR=cell(i,j-1)%hx
     hyR=cell(i,j-1)%hy

     nx=face_norm(4,1)
     ny=face_norm(4,2)

     line_jacob=0.5_prec*hxL
     do q=1,pord
        UL=calc_uh(line_coord(q),cons_coord,xcL,ycL,hxL,hyL,cell(i,j)%un(:,1:n_dg_fr),Aord,n_dg_fr)
        UR=calc_uh(line_coord(q),cons_coord,xcR,ycR,hxR,hyR,cell(i,j-1)%un(:,1:n_dg_fr),Aord,n_dg_fr)
        face_fluxes(q,:)=Roe(UL,UR,nx,ny)
     enddo
   
     do l=1,n_dg_fr
        do k=1,4
           do q=1,pord
              FX(q)=face_fluxes(q,k)*Legendre(line_coord(q),xcL,hxL,Aord(1,l))*Legendre(cons_coord,ycL,hyL,Aord(2,l))
           enddo
           Integral=Integrate_1D(FX,wg,line_jacob,pord)
           cell(i,j)%res(k,l)=cell(i,j)%res(k,l)+Word(l)/(hxL*hyL)*Integral
        enddo
      enddo
      
   enddo
enddo
!*******************************End Face Contribution**************

!*************Beging Integral Contribution*************************

do j=2,nCells_y+1
   do i=2,nCells_x+1
      ctr_x=cell(i,j)%xc
      ctr_y=cell(i,j)%yc
      cdx=cell(i,j)%hx
      cdy=cell(i,j)%hy
      gauss_x=cell(i,j)%gauss_pts_x(1:pord)
      gauss_y=cell(i,j)%gauss_pts_y(1:pord)
      c_jacob=cell(i,j)%jacob
      do r=1,pord
         do q=1,pord
            wrk_uh=calc_uh(gauss_x(q),gauss_y(r),ctr_x,ctr_y,cdx,&
                           cdy,cell(i,j)%un(:,1:n_dg_fr),Aord,n_dg_fr)
            call Euler_Flux(wrk_uh,gauss_Fx(q,r,:),gauss_Gy(q,r,:))
         enddo
      enddo

      do l=2,n_dg_fr
         do k=1,4
            do r=1,pord
               do q=1,pord
                  dphix=dLegendre(gauss_x(q),ctr_x,cdx,Bphix(l))*&
                        Legendre(gauss_y(r),ctr_y,cdy,Bphi_mix(l))
                  dphiy=dLegendre(gauss_y(r),ctr_y,cdy,Bphiy(l))*&
                        Legendre(gauss_x(q),ctr_x,cdx,Bphi_mix(l))
                  FXY(q,r)=gauss_Fx(q,r,k)*dphix+gauss_Gy(q,r,k)*dphiy
               enddo
            enddo
          Integral=Integrate_2D(FXY,wg,c_jacob,pord)
          cell(i,j)%res(k,l)=cell(i,j)%res(k,l)-Word(l)/(cdx*cdy)*&
                             Integral
         enddo
     enddo
    enddo
enddo
deallocate(face_fluxes)
deallocate(gauss_Fx)
deallocate(gauss_Gy)

end subroutine residual


