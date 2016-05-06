MODULE GRAD

  USE INIT
  USE DERIVATIVE
  USE INTERPOLATION
  USE TRANSFORM
  USE OUTPUT

CONTAINS

! ---------------------------------------------------------------------------------------------------------------------------

SUBROUTINE COMPUTE_GRADIENT(a, grad_xi, grad_eta, grad_r)

  implicit none

  include 'mpif.h'

  integer i,j,k,ii,jj,ierr
 
  real*8, dimension(mx,my,mr)     :: a, grad_xi, grad_eta, grad_r, deriv_xi, deriv_eta, deriv_r

  real*8, dimension(mx,noverlap,mr)  :: interpolated_south, interpolated_north
  real*8, dimension(noverlap,my,mr)  :: interpolated_south_transposed, interpolated_north_transposed
  real*8, dimension(noverlap, my,mr) :: interpolated_east, interpolated_west
  real*8, dimension(mx, noverlap,mr) :: interpolated_east_transposed, interpolated_west_transposed

  call INTERPOLATE_ALONG_X(a(:,2:noverlap+1,:), interpolated_south)
  call INTERPOLATE_ALONG_X(a(:,my-1:my-noverlap:-1,:), interpolated_north)
  call INTERPOLATE_ALONG_Y(a(2:noverlap+1,:,:), interpolated_west)
  call INTERPOLATE_ALONG_Y(a(mx-1:mx-noverlap:-1,:,:), interpolated_east)

! Transposing relevant arrays to make their shape conform to the one expected by derivative routines
  do k=1,mr
   do j=1,noverlap
    do i=1,mx
      interpolated_east_transposed(i,j,k) = interpolated_east(j,i,k)
      interpolated_west_transposed(i,j,k) = interpolated_west(j,i,k) 
      interpolated_north_transposed(j,i,k)  = interpolated_north(i,j,k)
      interpolated_south_transposed(j,i,k)  = interpolated_south(i,j,k)
    enddo
   enddo
  enddo

  call DDX(a, deriv_xi, interpolated_south_transposed, interpolated_north_transposed, interpolated_west, interpolated_east)
  call DDY(a, deriv_eta, interpolated_south, interpolated_north, interpolated_west_transposed, interpolated_east_transposed)
  call DDR(a, deriv_r)

  do k=1,mr
   do j=1,my
    do i=1,mx
     grad_xi(i,j,k)  = (1.0D0/r(k)) * (metric_D(j) * deriv_xi(i,j,k) + XY(i,j) * deriv_eta(i,j,k) /metric_D(j))
     grad_eta(i,j,k) = (1.0D0/r(k)) * (metric_C(i) * deriv_eta(i,j,k) + XY(i,j) * deriv_xi(i,j,k) /metric_C(i))
     grad_r(i,j,k)   = deriv_r(i,j,k)
    enddo
   enddo
  enddo

END SUBROUTINE COMPUTE_GRADIENT
! ---------------------------------------------------------------------------------------------------------------------------

SUBROUTINE COMPUTE_DIVERGENCE(a_xi, a_eta, a_r, divergence)
 
  implicit none

  include 'mpif.h'

  integer i,j,k, jj,ii, ierr
  real*8, dimension(mx,noverlap,mr) :: interpolated_north_theta, interpolated_north_phi
  real*8, dimension(mx,noverlap,mr) :: interpolated_south_theta, interpolated_south_phi
  real*8, dimension(noverlap,my,mr) :: interpolated_east_theta,  interpolated_east_phi
  real*8, dimension(noverlap,my,mr) :: interpolated_west_theta,  interpolated_west_phi

  COMMON /DIVERGENCE_POLAR/ interpolated_north_theta, interpolated_north_phi, &
                            interpolated_south_theta, interpolated_south_phi, &
                            interpolated_east_theta,  interpolated_east_phi, &
                            interpolated_west_theta,  interpolated_west_phi

  real*8, dimension(mx,noverlap,mr) :: interpolated_north_xi, interpolated_north_eta
  real*8, dimension(mx,noverlap,mr) :: interpolated_south_xi, interpolated_south_eta
  real*8, dimension(noverlap,my,mr) :: interpolated_east_xi,  interpolated_west_xi
  real*8, dimension(noverlap,my,mr) :: interpolated_east_eta, interpolated_west_eta

  COMMON /DIVERGENCE_CURVILINEAR/  interpolated_north_xi, interpolated_north_eta,&
                                   interpolated_south_xi, interpolated_south_eta,&
                                   interpolated_east_xi,  interpolated_west_xi,&
                                   interpolated_east_eta, interpolated_west_eta

  real*8, dimension(mx,noverlap,mr) :: interpolated_east_eta_transposed,  interpolated_west_eta_transposed
  real*8, dimension(noverlap,my,mr) :: interpolated_north_xi_transposed, interpolated_south_xi_transposed

  real*8, dimension(mx,my,mr) :: a_xi, a_eta, a_r, divergence, deriv_xi, deriv_eta, deriv_r
  real*8, dimension(mx,my,mr) :: a_theta, a_phi

  call CONVERT_XI_ETA_TO_THETA_PHI(a_xi, a_eta, a_theta, a_phi)

  call INTERPOLATE_ALONG_X(a_theta(:,2:noverlap+1,:), interpolated_south_theta)
  call INTERPOLATE_ALONG_X(a_theta(:,my-1:my-noverlap:-1,:), interpolated_north_theta)
  call INTERPOLATE_ALONG_Y(a_theta(2:noverlap+1,:,:), interpolated_west_theta)
  call INTERPOLATE_ALONG_Y(a_theta(mx-1:mx-noverlap:-1,:,:), interpolated_east_theta)
  
  call INTERPOLATE_ALONG_X(a_phi(:,2:noverlap+1,:), interpolated_south_phi)
  call INTERPOLATE_ALONG_X(a_phi(:,my-1:my-noverlap:-1,:), interpolated_north_phi)
  call INTERPOLATE_ALONG_Y(a_phi(mx-1:mx-noverlap:-1,:,:), interpolated_east_phi)
  call INTERPOLATE_ALONG_Y(a_phi(2:noverlap+1,:,:), interpolated_west_phi)

  call TRANSFORM_VECTOR_BETWEEN_MESHES_DIV

! Transposing relevant arrays to make their shape conform to that expected by derivative routines
  do k=1,mr
   do j=1,noverlap
    do i=1,mx
     interpolated_east_eta_transposed(i,j,k) = interpolated_east_eta(j,i,k)
     interpolated_west_eta_transposed(i,j,k) = interpolated_west_eta(j,i,k)
     interpolated_north_xi_transposed(j,i,k) = interpolated_north_xi(i,j,k)
     interpolated_south_xi_transposed(j,i,k) = interpolated_south_xi(i,j,k)
    enddo
   enddo
  enddo
 
  call DDX(a_xi, deriv_xi, interpolated_south_xi_transposed, interpolated_north_xi_transposed, interpolated_west_xi, interpolated_east_xi)
  call DDY(a_eta, deriv_eta, interpolated_south_eta, interpolated_north_eta, interpolated_west_eta_transposed, interpolated_east_eta_transposed)
  call DDR(a_r, deriv_r)

  do k=1,mr
   do j=1,my
    do i=1,mx
     divergence(i,j,k)  = (1.0D0/r(k))* (metric_delta(i,j)/CD(i,j) * ( deriv_xi(i,j,k)/metric_C(i) + &
                          deriv_eta(i,j,k)/metric_D(j) ) - tan(xi(i))/metric_D(j) * a_xi(i,j,k) - tan(eta(j))/metric_C(i) * a_eta(i,j,k) )+ &
                           deriv_r(i,j,k) + 2.0D0*a_r(i,j,k)/r(k)                 
    enddo
   enddo
  enddo


END SUBROUTINE COMPUTE_DIVERGENCE

! ---------------------------------------------------------------------------------------------------------------------------

SUBROUTINE COMPUTE_CURL(a_xi, a_eta, a_r, curl_xi, curl_eta, curl_r)

  implicit none

  include 'mpif.h'

  integer i,j,k, jj,ii, ierr

  real*8, dimension(mx,noverlap,mr)  :: interpolated_south_r, interpolated_north_r
  real*8, dimension(noverlap,my,mr)  :: interpolated_south_r_transposed, interpolated_north_r_transposed
  real*8, dimension(noverlap, my,mr) :: interpolated_east_r, interpolated_west_r
  real*8, dimension(mx, noverlap,mr) :: interpolated_east_r_transposed, interpolated_west_r_transposed

  real*8, dimension(mx,noverlap,mr) :: interpolated_north_theta, interpolated_north_phi
  real*8, dimension(mx,noverlap,mr) :: interpolated_south_theta, interpolated_south_phi
  real*8, dimension(noverlap,my,mr) :: interpolated_east_theta,  interpolated_east_phi
  real*8, dimension(noverlap,my,mr) :: interpolated_west_theta,  interpolated_west_phi

  COMMON /CURL_POLAR/ interpolated_north_theta, interpolated_north_phi, &
                      interpolated_south_theta, interpolated_south_phi, &
                      interpolated_east_theta,  interpolated_east_phi, &
                      interpolated_west_theta,  interpolated_west_phi

  real*8, dimension(mx,noverlap,mr) :: interpolated_north_xi, interpolated_north_eta
  real*8, dimension(mx,noverlap,mr) :: interpolated_south_xi, interpolated_south_eta
  real*8, dimension(noverlap,my,mr) :: interpolated_east_xi,  interpolated_west_xi
  real*8, dimension(noverlap,my,mr) :: interpolated_east_eta, interpolated_west_eta

  COMMON /CURL_CURVILINEAR/  interpolated_north_xi, interpolated_north_eta,&
                             interpolated_south_xi, interpolated_south_eta,&
                             interpolated_east_xi,  interpolated_west_xi,&
                             interpolated_east_eta, interpolated_west_eta

  real*8, dimension(mx,noverlap,mr) :: interpolated_east_eta_transposed,  interpolated_west_eta_transposed
  real*8, dimension(noverlap,my,mr) :: interpolated_north_xi_transposed, interpolated_south_xi_transposed

  real*8, dimension(mx,my,mr) :: a_xi, a_eta, a_r, curl_xi, curl_eta, curl_r 
  real*8, dimension(mx,my,mr) :: deriv_xi_r, deriv_eta_r
  real*8, dimension(mx,my,mr) :: deriv_r_xi, deriv_xi_xi, deriv_eta_xi
  real*8, dimension(mx,my,mr) :: deriv_r_eta, deriv_xi_eta, deriv_eta_eta
  
  real*8, dimension(mx,my,mr) :: a_theta, a_phi

  call INTERPOLATE_ALONG_X(a_r(:,2:noverlap+1,:), interpolated_south_r)
  call INTERPOLATE_ALONG_X(a_r(:,my-1:my-noverlap:-1,:), interpolated_north_r)
  call INTERPOLATE_ALONG_Y(a_r(2:noverlap+1,:,:), interpolated_west_r)
  call INTERPOLATE_ALONG_Y(a_r(mx-1:mx-noverlap:-1,:,:), interpolated_east_r)

! Transposing relevant arrays to make their shape conform to the one expected by derivative routines
  do k=1,mr
   do j=1,noverlap
    do i=1,mx
      interpolated_east_r_transposed(i,j,k) = interpolated_east_r(j,i,k)
      interpolated_west_r_transposed(i,j,k) = interpolated_west_r(j,i,k)
      interpolated_north_r_transposed(j,i,k)  = interpolated_north_r(i,j,k)
      interpolated_south_r_transposed(j,i,k)  = interpolated_south_r(i,j,k)
    enddo
   enddo
  enddo

  call CONVERT_XI_ETA_TO_THETA_PHI(a_xi, a_eta, a_theta, a_phi)

  call INTERPOLATE_ALONG_X(a_theta(:,2:noverlap+1,:), interpolated_south_theta)
  call INTERPOLATE_ALONG_X(a_theta(:,my-1:my-noverlap:-1,:), interpolated_north_theta)
  call INTERPOLATE_ALONG_Y(a_theta(2:noverlap+1,:,:), interpolated_west_theta)
  call INTERPOLATE_ALONG_Y(a_theta(mx-1:mx-noverlap:-1,:,:), interpolated_east_theta)

  call INTERPOLATE_ALONG_X(a_phi(:,2:noverlap+1,:), interpolated_south_phi)
  call INTERPOLATE_ALONG_X(a_phi(:,my-1:my-noverlap:-1,:), interpolated_north_phi)
  call INTERPOLATE_ALONG_Y(a_phi(mx-1:mx-noverlap:-1,:,:), interpolated_east_phi)
  call INTERPOLATE_ALONG_Y(a_phi(2:noverlap+1,:,:), interpolated_west_phi)

  call TRANSFORM_VECTOR_BETWEEN_MESHES_CURL

! Transposing relevant arrays to make their shape conform to that expected by derivative routines
  do k=1,mr
   do j=1,noverlap
    do i=1,mx
     interpolated_east_eta_transposed(i,j,k) = interpolated_east_eta(j,i,k)
     interpolated_west_eta_transposed(i,j,k) = interpolated_west_eta(j,i,k)
     interpolated_north_xi_transposed(j,i,k) = interpolated_north_xi(i,j,k)
     interpolated_south_xi_transposed(j,i,k) = interpolated_south_xi(i,j,k)
    enddo
   enddo
  enddo

  call DDX(a_r, deriv_r_xi, interpolated_south_r_transposed, interpolated_north_r_transposed, interpolated_west_r, interpolated_east_r)
  call DDY(a_r, deriv_r_eta, interpolated_south_r, interpolated_north_r, interpolated_west_r_transposed, interpolated_east_r_transposed)

  call DDX(a_xi, deriv_xi_xi, interpolated_south_xi_transposed, interpolated_north_xi_transposed, interpolated_west_xi, interpolated_east_xi)
  call DDY(a_eta, deriv_eta_eta, interpolated_south_eta, interpolated_north_eta, interpolated_west_eta_transposed, interpolated_east_eta_transposed)

  call DDX(a_eta, deriv_eta_xi, interpolated_south_eta, interpolated_north_eta, interpolated_west_eta_transposed, interpolated_east_eta_transposed)
  call DDY(a_xi, deriv_xi_eta, interpolated_south_xi_transposed, interpolated_north_xi_transposed, interpolated_west_xi, interpolated_east_xi)

  call DDR(a_xi, deriv_xi_r)
  call DDR(a_eta, deriv_eta_r)

  do k=1,mr
   do j=1,my
    do i=1,mx
     curl_xi(i,j,k)  =             XY(i,j)/sqrt(metric_delta(i,j))*deriv_xi_r(i,j,k) - CD(i,j)/sqrt(metric_delta(i,j))*deriv_eta_r(i,j,k) &
                     + (1D0/r(k))*(XY(i,j)/sqrt(metric_delta(i,j))*a_xi(i,j,k)       - CD(i,j)/sqrt(metric_delta(i,j))*a_eta(i,j,k) + &
                                   sqrt(metric_delta(i,j))/metric_D(j)*deriv_r_eta(i,j,k))
     curl_eta(i,j,k) =             CD(i,j)/sqrt(metric_delta(i,j))*deriv_xi_r(i,j,k) - XY(i,j)/sqrt(metric_delta(i,j))*deriv_eta_r(i,j,k) &
                     + (1D0/r(k))*(CD(i,j)/sqrt(metric_delta(i,j))*a_xi(i,j,k)       - XY(i,j)/sqrt(metric_delta(i,j))*a_eta(i,j,k) - &
                                   sqrt(metric_delta(i,j))/metric_C(i)*deriv_r_xi(i,j,k))
     curl_r(i,j,k)   = (sqrt(metric_delta(i,j))/r(k)) * ( XY(i,j)/CD(i,j) * ( deriv_eta_eta(i,j,k)/metric_D(j) - deriv_xi_xi(i,j,k)/metric_C(i) )&
                                                        - deriv_xi_eta(i,j,k)/metric_D(j) + deriv_eta_xi(i,j,k)/metric_C(i) )
    enddo
   enddo
  enddo
 


END SUBROUTINE COMPUTE_CURL

END MODULE GRAD


