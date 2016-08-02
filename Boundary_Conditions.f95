subroutine Periodic_BC()
use Variables
implicit none
integer::i,j,k,l

do i=1,nCells_x+2
   do l=1,n_dg_fr
      do k=1,4
         cell(i,1)%un(k,l)=cell(i,nCells_y+1)%un(k,l)
         cell(i,nCells_x+2)%un(k,l)=cell(i,2)%un(k,l)
      enddo
   enddo
enddo


do j=1,nCells_y+2 
   do l=1,n_dg_fr
      do k=1,4
         cell(1,j)%un(k,l)=cell(nCells_x+1,j)%un(k,l)
         cell(nCells_y+2,j)%un(k,l)=cell(2,j)%un(k,l)
      enddo
   enddo
enddo

end subroutine Periodic_BC
