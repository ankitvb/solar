SUBROUTINE BC_RHS(q, q_rhs)

 USE INIT

 implicit none

 real*8, dimension(mx,my,mr,5) :: q, q_rhs

 integer i,j,k,l

 if(nr.eq.1) return

! Zero Dirichlet boundary conditions
 if(rrank.eq.pr-1)then
  k=mr
  do l=1,5
  do i=1,mx
   do j=1,my
    q_rhs(i,j,k,l) = 0D0
   enddo
  enddo
  enddo
 endif

 if(rrank.eq.0)then
  k=1
  do l=1,5
  do i=1,mx
   do j=1,my
    q_rhs(i,j,k,l) = 0D0
   enddo
  enddo
  enddo
 endif

END SUBROUTINE BC_RHS
