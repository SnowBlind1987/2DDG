subroutine output(t,time)
use Variables
use Parameters
use Conversion_Functions
implicit none

integer::t
real(prec)::time
character(len=4)::sol_number
character(len=16)::file_name
character(len=4)::str_nc_x,str_nc_y
integer::i,j
real(prec)::w(4)
real(prec)::spd_snd, Mach
write(sol_number,'(i4)') t !convert t to string
write(str_nc_x,'(i4)') nCells_x
write(str_nc_y,'(i4)') nCells_y

sol_number=adjustl(sol_number)
str_nc_x=adjustl(str_nc_x)
str_nc_y=adjustl(str_nc_y)
if (t<10) then
        file_name='solution'//'000'//trim(sol_number)//'.dat'
elseif (t<100) then
        file_name='solution'//'00'//trim(sol_number)//'.dat'
elseif (t<1000) then
        file_name='solution'//'0'//trim(sol_number)//'.dat'
else
        file_name='solution'//trim(sol_number)//'.dat'
endif

w(:)=0.0_prec
open(21,file='solution/'//file_name)
write(21,'(A59)') 'VARIABLES = "X[m]","Y[m]","rho [kg/m^3]","U[m/s]","V[m/s]"'
write(21,'(A25)')'ZONE T="Rectangular zone"'
write(21,'(A6,A6,A17)') 'I='//trim(str_nc_x)//',','J='//trim(str_nc_y)//',','ZONETYPE=Ordered'
write(21,'(A18)') 'DATAPACKING=POINT'


do j=2,nCells_y+1
   do i=2,nCells_x+1
   w=cons2prim(cell(i,j)%un(:,1))
   spd_snd=sqrt(gam*w(4)/w(1))
   Mach=w(2)/spd_snd
   write(21,'(5F17.9)') cell(i,j)%xc,cell(i,j)%yc, w(1), w(2), w(3)

   enddo
enddo
close (21)

open(23,file='solution/time.dat',access='append')
write(23,'(F16.8)') time
close(23)
end subroutine output
