!-----------------------------------------------------------------------
SUBROUTINE COMPUTE_TIMESTEP

 USE INIT
 USE AUXILIARY

 implicit none

 include 'mpif.h'

 integer i,j,k
 real*8, dimension(mx,my,mr,5) :: q

 COMMON /RKWORK/ q

 real*8, dimension(mx,my,mr) :: cs, dt_local

 if(nr.eq.1) then
  deltat = dxi/sqrt(2D0)
  return
 endif

 do k=1,mr
  do j=1,my
   do i=1,mx
    cs(i,j,k)        = sqrt(gamma(i,j,k)*(q(i,j,k,5)+p0(i,j,k))/(q(i,j,k,1)+rho0(i,j,k)))
    dt_local(i,j,k)  = 1D0 / (abs(q(i,j,k,1)/dxi) + abs(q(i,j,k,2)/deta) + abs(q(i,j,k,3)/dchi) &
                     + cs(i,j,k) * sqrt(1D0/dxi**2D0 + 1D0/deta**2D0 + 1D0/dchi**2D0))
   enddo
  enddo
 enddo

 call COMPUTE_MIN(dt_local, deltat)
 
 deltat = 5D0*dimc/diml
 
 return

END SUBROUTINE COMPUTE_TIMESTEP  
!------------------------------------------------------------------------
SUBROUTINE RK4

 USE INIT
 USE FILTER

 implicit none

 real*8, dimension(mx,my,mr,5) :: q, q_rhs, qw, qf

 COMMON /RKWORK/ q

 integer k,l

!  Each Runge-Kutta substep involves the following operations:
!  - update time variable
!  - compute dq/dt according to the governing equations
!  - time advance solutions
!  - filter fields if necessary
!  - set boundary conditions

 qw = 0D0

 do k=1,N_stage_rk4
  call compute_rhs(q, q_rhs)
  call bc_rhs(q, q_rhs)
  call sponge(q, q_rhs)

  qw = alpha_rk4(k) * qw + deltat * q_rhs
  q  = q + beta_rk4(k) * qw
 
 enddo

 if(mod(step,1).eq.0)then
   call FILTER_EXP(q,qf)
   q = qf
!  do l=1,5
!   call FILTER_R(q(1,1,1,l),qf(1,1,1,l))
!   q(:,:,:,l)=qf(:,:,:,l)! - sigma_d * qf
!  enddo
 endif

 time = time + deltat

 return

END SUBROUTINE RK4
!-----------------------------------------------------------------------
