MODULE INIT

 implicit none

 include 'header'


! This module contains declarations for all global arrays.
!-----------------------------------------------------------------------------------------------------------------------
! Cubed-sphere metric terms and coordinates
!-----------------------------------------------------------------------------------------------------------------------
 real*8, dimension(nx,ny,6) :: theta, phi
 real*8, dimension(nx)      :: xi, tanxi
 real*8, dimension(ny)      :: eta, taneta
 real*8, dimension(nr)      :: r, chi
 real*8, dimension(nx)      :: gxi, geta
 real*8, dimension(nr)      :: gr, gchi
 
 COMMON /COORDS/ theta, phi, xi, tanxi, eta, taneta, r, chi, gxi, geta, gr, gchi

 real*8, dimension(nx,ny)   :: metric_delta, XY, CD, sqrt_metric_delta
 real*8, dimension(nx)      :: metric_C
 real*8, dimension(ny)      :: metric_D
 real*8, dimension(nr)      :: drdchi
 real*8, dimension(nr)      :: gdrdchi

 COMMON /METRIC_TERMS/ metric_delta, XY, CD, sqrt_metric_delta, &
                 metric_C, metric_D, drdchi, gdrdchi
  
!-----------------------------------------------------------------------------------------------------------------------
! Interpolation
!-----------------------------------------------------------------------------------------------------------------------
 real*8, dimension(nx,noverlap)             :: interpoints
 real*8, dimension(-nord:nord,nx, noverlap) :: interpolant

 COMMON /INTERPOLATE/ interpoints, interpolant

! Mesh adjacency
 Integer interfaces(6,6)

 COMMON /INTERFACES/ interfaces

! Lagrange time interpolation stuff
 real*8  x0, x1, x2, x3, x4, x5, x6

 COMMON /TIME_INTERP/  x0, x1, x2, x3, x4, x5, x6
 
 real*8, dimension(nx,ny) :: LC0, LC1, LC2, LC3, LC4, LC5, LC6

 COMMON /LAGRANGE_TIME_INTERP/  LC0, LC1, LC2, LC3, LC4, LC5, LC6

 integer step_old

 COMMON /STEP_TRACK/ step_old

!-----------------------------------------------------------------------------------------------------------------------
! Solar parameters for radially stratified field
!-----------------------------------------------------------------------------------------------------------------------
 real*8, dimension(mx,my,mr) :: c_speed, c2, rho0, p0, g0, gamma

 COMMON /BACKGROUND_CONST/ c_speed, c2, rho0, p0, g0, gamma 

 real*8 dimc, dimrho, diml, rsun, excitdep, obsheight, cadforcing

 COMMON /SOLAR_PARAMS/ dimc, dimrho, diml, rsun, excitdep, obsheight, cadforcing

 integer e_rad, o_rad

 COMMON /EO_RADII/ e_rad, o_rad  

END MODULE INIT
