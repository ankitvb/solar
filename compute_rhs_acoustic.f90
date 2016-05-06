SUBROUTINE COMPUTE_RHS(q, q_rhs)
 USE INIT
 USE GRAD 
 use auxiliary
 implicit none

 integer i,j,k,kk,m
 
 real*8, dimension(mx,my,mr,5) :: q, q_rhs 
 real*8, dimension(mx,my,mr)   :: div_v
 real*8, dimension(mx,my,mr,3) :: grad_p 
 real*8, dimension(mx,my,mr)   :: drho0dr
 real*8, dimension(mx,my)      :: forcing
 real*8  mag 

 COMMON /TIME_FORCING/ forcing

 q_rhs = 0D0                                                                          

! Primitive variables 
! q1 = rho 
! q2 = v_xi
! q3 = v_eta
! q4 = v_r
! q5 = p

! Computing divergences 
! div(v)
 call COMPUTE_DIVERGENCE(q(1,1,1,2), q(1,1,1,3), q(1,1,1,4), div_v)

! Computing pressure gradient
 call COMPUTE_GRADIENT(q(1,1,1,5), grad_p(1,1,1,1), grad_p(1,1,1,2), grad_p(1,1,1,3))

!  d rho            __
!  ----- = - rho0 ( \/ . v )
!   d t

 do k=1,mr
  do j=1,my
   do i=1,mx
    q_rhs(i,j,k,1) = -rho0(i,j,k) * div_v(i,j,k)
   enddo
  enddo
 enddo


!  d rho              d      
!  ----- = ... - v_r --- (rho0)
!   d t              d r

 call DDR(rho0, drho0dr)

 do k=1,mr
  do j=1,my
   do i=1,mx
    q_rhs(i,j,k,1) = q_rhs(i,j,k,1) - q(i,j,k,4) * drho0dr(i,j,k)
   enddo
  enddo
 enddo


! d v      __
! ---  = - \/ p / rho0
! d t   

 do k=1,mr
  do i=1,mx
   do j=1,my
    do m=1,3
     q_rhs(i,j,k,m+1) = -grad_p(i,j,k,m) / rho0(i,j,k)
    enddo
   enddo
  enddo
 enddo
 
! d v
! --- ... - rho g0 / rho0 e_r  
! d t

 do k=1,mr
  do j=1,my
   do i=1,mx
    q_rhs(i,j,k,4) = q_rhs(i,j,k,4) - q(i,j,k,1) * g0(i,j,k) / rho0(i,j,k)
   enddo
  enddo
 enddo

! d v
! --- ... + S(theta, phi, t) e_r 
! d t

  mag = exp(-(time*diml/dimc-650D0)**2D0/(2D0*200D0**2D0))*&
                    cos(4D0*pi*(time*diml/dimc-650D0)/700D0)
  do k=1,mr
   do j=1,my
    do i=1,mx
     kk = rrank*mr+k
     forcing(i,j) = mag*exp(-((theta(i,j)-pi/2D0)/(20D0*dxi))**2D0&
                            -((phi(i,j)-0D0)/(20D0*dxi))**2D0&
                            -((gr(kk)-gr(200))/(20D0*dchi))**2D0)
     q_rhs(i,j,k,4) = q_rhs(i,j,k,4) + forcing(i,j)
   enddo
  enddo
 enddo

! call LAGRANGE_INTERP

!  do k=1,mr
!   do j=1,my
!    do i=1,mx
!     kk = rrank*mr+k
!     if(kk.eq.e_rad) &
!      q_rhs(i,j,k,4) = q_rhs(i,j,k,4) + forcing(i,j)
!   enddo
!  enddo
! enddo

! mag = norm2(forcing)
 if(rank.eq.0)then
   open(unit=14, file='forcing_rms',status='unknown', form='formatted', action='write',&
           position='append')
   if(abs(mag).LT.1d-50) mag = 0D0 
   write(14,7) mag, time 
   close(14)
 endif

! d p               __
! --- = -c0^2 rho0 (\/ . v) 
! d t

 do k=1,mr
  do j=1,my
   do i=1,mx
    q_rhs(i,j,k,5) = -c2(i,j,k) * rho0(i,j,k) * div_v(i,j,k)
   enddo
  enddo
 enddo

! d p
! --- = ... + rho0 g0 v_r
! d t

 do k=1,mr
  do j=1,my
   do i=1,mx
    q_rhs(i,j,k,5) = q_rhs(i,j,k,5) + rho0(i,j,k) * g0(i,j,k) * q(i,j,k,4)
   enddo
  enddo
 enddo

! print *, minval(forcing(:,:)), maxval(forcing(:,:))
! print *, minval(q(:,:,:,4)), maxval(q(:,:,:,4)), 'vr'
! print *, minval(div_v(:,:,:)), maxval(div_v(:,:,:)), 'div_v'

 return

7 Format(1x, 2(ES13.5))

END SUBROUTINE COMPUTE_RHS
