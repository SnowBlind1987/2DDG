module Variables
use Set_Precision
implicit none
!****************Physical Variables*****************************
!***************************************************************
integer:: nCells_x,nCells_y
real(prec):: cfl,final_time 
integer:: pord !polynomial order
integer:: n_dg_fr
integer:: Mproj !Projection limiting constant
real(prec):: xmin,xmax
real(prec):: ymin,ymax
real(prec):: dx,dy
type cell_data
        real(prec)::xc !cell center x
        real(prec)::yc !cell center y
        real(prec)::gauss_pts_x(5)
        real(prec)::gauss_pts_y(5)
        real(prec)::face_const_coord(4) !the coordinate that stays constant
                                        !for each face (x for 1 and 3)
        real(prec)::jacob
        real(prec)::un(4,6) !conserved variables at time n+1
        real(prec)::uo(4,6) !conserved variables at time n
        real(prec)::res(4,6) !the resediual at each cell
        real(prec)::hx !length dx of each cell
        real(prec)::hy !length dy of each cell
end type cell_data

real(prec)::face_norm(4,2)

integer::Aord(2,6)
real(prec)::Word(6)
integer::Bphix(6), Bphiy(6),Bphi_mix(6)

type(cell_data),allocatable:: cell(:,:) !Array of cells
!*********************************************************
!*****************Quadrature Variables********************
!*********************************************************
real(kind=prec)::xg(5,5),wg(5,5)
!integer::ng
!*********************************************************
!********************Runge-Kutta Variables****************
real(prec)::ark(3),brk(3)


end module Variables
