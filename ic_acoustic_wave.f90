!-----------------------------------------------------------------!
! Routines to generate initial conditions                         !
!-----------------------------------------------------------------!
SUBROUTINE GENERATE_IC(q)
 USE INIT

 implicit none
 
 real*8, dimension(mx,my,mr,5) :: q

 integer i,j,k

 time = 0d0
 q    = 0d0

! Gaussian acoustic pulse 
 do k=1,mr
  do j=1,my
   do i=1,mx
    q(i,j,k,1) = 0D0     ! background density --add perturbation
    q(i,j,k,2) = 0D0     ! xi velocity
    q(i,j,k,3) = 0D0     ! eta velocity
    q(i,j,k,4) = 0D0     ! radial velocity
    q(i,j,k,5) = 0D0   !exp(-(theta(i,j)-0.5D0*pi)**2D0/0.01D0)!*(r(k) -gr(nr))*(r(k)-gr(1))!  * exp(-(r(k)-0.95D0)**2D0/0.0001D0)
   enddo
  enddo
 enddo

END SUBROUTINE GENERATE_IC
!---------------------------------------------------------------------------------------------
SUBROUTINE SOLAR_DATA
 USE INIT

 implicit none

 include 'mpif.h'

 integer i,j,k,q,kk
 real*8 data(nr,8)
 real*8 temp, temp2, temp_r

 integer ierr

! if(nr.eq.1) return

! Data in the file is arranged as 
! r/R, c_0, \rho_0,  p_0, T_0, g_0, \Gamma_1, random
! Solar radius,sound speed, density, pressure in cgs units. Should take care
! to convert them to SI. 10 dyne/cm^2 = 1 Pa, 1000 kg/m^3  = 1 g/cc

! Linear interpolation - could improve on this
! c rho p T g
! solar mass: 1.9891d33 g
! Rsun = 695.9894 *10.^8. cm
! mass at the middle of the convective zone (0.85R): 1.9672367d33g
! Everything is in cgs units. I have to multiply rsun by a 100 to produce
! cgs

 
! Root processor reads input file
 if(rank.eq.0)then
  open(7,file = 'solar_model',form = 'formatted',status = 'old')
  do k = 1,nr
    read(7,*) data(k,:)
  enddo
  close(7)
 endif


! Broadcasts it to all other processors
 call MPI_BCAST(data, nr*8, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

 q = 1
 dimc = data(q,2)
 dimrho = data(q,3)
 rsun = 695.9894D8

 diml = rsun
 gr   = data(:,1)

 if(nr.eq.1)then   ! 2d testing 
  do i=1,mx
   do j=1,my
     c_speed(i,j,1) = 1D0
     c2(i,j,1)      = 1D0 
     rho0(i,j,1)    = 1D0
     p0(i,j,1)      = 1D0
     g0(i,j,1)      = 0D0 
     gamma(i,j,1)   = 1.67D0 
     r(1)           = 1D0 
   enddo
  enddo
 else
  do k=1,mr
   kk = rrank*mr + k
   do i=1,mx
    do j=1,my
 
     c_speed(i,j,k) = data(kk,2)/dimc  
     c2(i,j,k)      = (data(kk,2)/dimc)**2.0 
     rho0(i,j,k)    = data(kk,3)/dimrho
     p0(i,j,k)      = data(kk,4)/(dimrho*dimc**2.0)
     g0(i,j,k)      = data(kk,6)*diml/dimc**2.0
     gamma(i,j,k)   = data(kk,7)
     r(k)           = data(kk,1)
    enddo
   enddo
  enddo
 endif

 excitdep = 1.0 - 100D5/rsun
 obsheight = 1.0 + 200D5/rsun 
 cadforcing = 60.0*dimc/diml

!===============================================================================================
 ! Determination of excitation location.
   temp = 1.0
   do k=1,nr
    if (abs(gr(k) - excitdep) .LE. temp) then
     temp = abs(gr(k) - excitdep)
     e_rad = k
    endif
   enddo

!===============================================================================================
   ! Determination of the output height
   temp = 1.0
   do k=1,nr
    if (abs(gr(k) - obsheight) .LE. temp) then
     temp = abs(gr(k) - obsheight)
     o_rad = k
    endif
   enddo

!  print *, o_rad

!===============================================================================================

END SUBROUTINE SOLAR_DATA
