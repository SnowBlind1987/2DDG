subroutine save_sol()
use Variables
implicit none
integer::i,j,k,l

do j=1,nCells_y+2
   do i=1,nCells_x+2
      do l=1,n_dg_fr
         do k=1,4
            cell(i,j)%uo(k,l)=cell(i,j)%un(k,l)
         enddo
      enddo
   enddo
enddo

end subroutine save_sol

