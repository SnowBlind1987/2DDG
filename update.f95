subroutine update(rk,dt)
use Variables
integer,intent(in)::rk
real(prec),intent(in)::dt
real(prec)::tau
integer::i,j,k,l
tau=dt
do j=2,nCells_y+1
   do i=2,nCells_x+1
      do k=1,4
         do l=1,n_dg_fr
            cell(i,j)%un(k,l)=ark(rk)*cell(i,j)%uo(k,l)+ &
                               brk(rk)*(cell(i,j)%un(k,l)-tau*cell(i,j)%Res(k,l))
         enddo
      enddo
   enddo
enddo

end subroutine update
