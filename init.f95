subroutine init()
use Variables 
use Initial_Conditions
use Legendre_Functions
use Quadrature
implicit none
real(prec)::ctr_x, ctr_y
real(prec)::cdx,cdy
integer:: i,j,k,l,m,q,r
real(prec)::gx,gy
real(prec),allocatable::gauss_u(:,:,:)
real(prec)::wrk_u(4)
real(prec)::FXY(pord,pord)
allocate(cell(1:nCells_x+2,1:nCells_y+2))

!****Legendre Function Order
Aord(1,1)=1
Aord(1,2)=2
Aord(1,3)=1
Aord(1,4)=2
Aord(1,5)=3
Aord(1,6)=1

Aord(2,1)=1
Aord(2,2)=1
Aord(2,3)=2
Aord(2,4)=2
Aord(2,5)=1
Aord(2,6)=3

Word(1)=1.0_prec
Word(2)=3.0_prec
Word(3)=3.0_prec
Word(4)=9.0_prec
Word(5)=5.0_prec
Word(6)=5.0_prec

Bphix(1)=0
Bphix(2)=2
Bphix(3)=1
Bphix(4)=2
Bphix(5)=3
Bphix(6)=1

Bphi_mix(1)=0
Bphi_mix(2)=1
Bphi_mix(3)=1
Bphi_mix(4)=2
Bphi_mix(5)=1
Bphi_mix(6)=1

Bphiy(1)=0
Bphiy(2)=1
Bphiy(3)=2
Bphiy(4)=2
Bphiy(5)=1
Bphiy(6)=3


!***************************************
!**************Face Norms***************
!***************************************
face_norm(1,1)=1.0_prec
face_norm(1,2)=0.0_prec
face_norm(2,1)=0.0_prec
face_norm(2,2)=1.0_prec
face_norm(3,1)=-1.0_prec
face_norm(3,2)=0.0_prec
face_norm(4,1)=0.0_prec
face_norm(4,2)=-1.0_prec
select case (pord)
case (1)
        n_dg_fr=1
case (2) 
        n_dg_fr=3
case (3) 
        n_dg_fr=6
end select
allocate(gauss_u(1:pord,1:pord,1:4))
gauss_u(:,:,:)=0.0_prec
wrk_u(:)=0.0_prec
!zero every variable 

do j=1,nCells_y+2
   do i=1,nCells_x+2
   cell(i,j)%xc=0.0_prec
   cell(i,j)%yc=0.0_prec
   cell(i,j)%jacob=0.0_prec
   do m=1,5
      cell(i,j)%gauss_pts_x(m)=0.0_prec
      cell(i,j)%gauss_pts_y(m)=0.0_prec
   enddo
   
   do m=1,4
      cell(i,j)%face_const_coord(m)=0.0_prec
   enddo

   do l=1,6
      do k=1,4
         cell(i,j)%un(k,l)=0.0_prec
         cell(i,j)%uo(k,l)=0.0_prec
         cell(i,j)%res(k,l)=0.0_prec
      enddo
   enddo

   cell(i,j)%hx=0.0_prec
   cell(i,j)%hy=0.0_prec

   enddo
enddo
dx=(xmax-xmin)/nCells_x
dy=(ymax-ymin)/nCells_y

!Building the mesh
do j=1,nCells_y+2
   do i=1,nCells_x+2
      cell(i,j)%xc=xmin+dx/2.0_prec+(i-2)*dx
      cell(i,j)%yc=ymin+dy/2.0_prec+(j-2)*dy
      cell(i,j)%hx=dx
      cell(i,j)%hy=dy
      cell(i,j)%face_const_coord(1)=cell(i,j)%xc+dx/2.0_prec
      cell(i,j)%face_const_coord(2)=cell(i,j)%yc+dy/2.0_prec
      cell(i,j)%face_const_coord(3)=cell(i,j)%xc-dx/2.0_prec
      cell(i,j)%face_const_coord(4)=cell(i,j)%yc-dy/2.0_prec

   enddo
enddo
!***********************************************************
do j=2,nCells_y+1
   do i=2,nCells_x+1
      ctr_x=cell(i,j)%xc
      ctr_y=cell(i,j)%yc
      cdx=cell(i,j)%hx
      cdy=cell(i,j)%hy
      !calculate the gauss points
      do m=1,pord
         cell(i,j)%gauss_pts_x(m)=ctr_x+0.5_prec*xg(pord,m)*cdx
         cell(i,j)%gauss_pts_y(m)=ctr_y+0.5_prec*xg(pord,m)*cdy
      enddo
      cell(i,j)%jacob=0.25_prec*cdx*cdy
    enddo
enddo

do j=2,nCells_y+1
   do i=2,nCells_x+1
   ctr_x=cell(i,j)%xc
   ctr_y=cell(i,j)%yc
   cdx=cell(i,j)%hx
   cdy=cell(i,j)%hy   
   do r=1,pord
      do q=1,pord
         gx=cell(i,j)%gauss_pts_x(q)
         gy=cell(i,j)%gauss_pts_y(r)
         gauss_u(q,r,:)=Vortex_Convection(gx,gy)
         !do k=1,4
           ! gauss_u(q,r,k)=wrk_u(k)
         !enddo
      enddo
   enddo
   
   do l=1,n_dg_fr
      do k=1,4
         do r=1,pord
            do q=1,pord
               gx=cell(i,j)%gauss_pts_x(q)
               gy=cell(i,j)%gauss_pts_y(r)
               FXY(q,r)=gauss_u(q,r,k)*Legendre(gx,ctr_x,cdx,Aord(1,l))&
                        *Legendre(gy,ctr_y,cdy,Aord(2,l))
            enddo
         enddo
        
       cell(i,j)%un(k,l)=Word(l)/(cdx*cdy)*Integrate_2D(FXY,wg,cell(i,j)%jacob,pord)

      enddo
    enddo

   enddo
enddo
deallocate(gauss_u)

!***************************************************************************************
!**************************************************************************************
!****************Initialize Runge-Kutta Variables******************************
ark(1)=0.0_prec
ark(2)=3.0_prec/4.0_prec
ark(3)=1.0_prec/3.0_prec

brk(1)=1.0_prec
brk(2)=1.0_prec/4.0_prec
brk(3)=2.0_prec/3.0_prec

end subroutine init
