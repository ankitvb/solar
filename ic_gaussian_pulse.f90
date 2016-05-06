!-----------------------------------------------------------------!
! Routines to generate initial conditions                         !
!-----------------------------------------------------------------!
SUBROUTINE GENERATE_IC(q)
 USE INIT

 implicit none
 
! include 'header'

 real*8, dimension(mx,my,mr,5) :: q

 integer i,j,k

 time = 0d0
 q    = 0d0

! Gaussian acoustic pulse 
 do k=1,mr
  do j=1,my
   do i=1,mx
    q(i,j,k,1) = exp(-(theta(i,j) - 0.5D0*pi)**2D0/0.1D0)
   enddo
  enddo
 enddo

 q(:,:,:,2) = 0D0

END SUBROUTINE GENERATE_IC

