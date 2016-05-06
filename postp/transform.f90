MODULE TRANSFORM

USE INIT

CONTAINS

!--------------------------------------------------------------------------------------------  
SUBROUTINE CONVERT_XI_ETA_TO_THETA_PHI(a_xi, a_eta, a_theta, a_phi)

 implicit none

 integer i,j,k
 real*8, dimension(nx, ny, 6) :: a_xi, a_eta, a_theta, a_phi
 real*8 factor


 do k=1,4
  do j=1,ny
   do i=1,nx
    a_theta(i,j,k) = XY(i,j)/CD(i,j) * a_xi(i,j,k) - a_eta(i,j,k)
    a_phi(i,j,k) = sqrt_metric_delta(i,j) / CD(i,j) * a_xi(i,j,k)
   enddo
  enddo
 enddo

 do j=1,ny
  do i=1,nx

   factor = 1.0/(metric_delta(i,j) - 1.0)**0.5
   a_theta(i,j,5) =  factor*( tanxi(i)/metric_D(j)* a_xi(i,j,5) + &
                taneta(j)/metric_C(i) * a_eta(i,j,5) )

   a_phi(i,j,5) = factor*( -sqrt_metric_delta(i,j) * taneta(j)/metric_D(j) * a_xi(i,j,5) + &
                sqrt_metric_delta(i,j) * tanxi(i)/metric_C(i) * a_eta(i,j,5))


   a_theta(i,j,6) =  -factor*( tanxi(i)/metric_D(j)* a_xi(i,j,6) + &
                taneta(j)/metric_C(i) * a_eta(i,j,6) )

   a_phi(i,j,6) = factor*( sqrt_metric_delta(i,j) * taneta(j)/metric_D(j) * a_xi(i,j,6) - &
                sqrt_metric_delta(i,j) * tanxi(i)/metric_C(i) * a_eta(i,j,6))

  enddo
 enddo


END SUBROUTINE CONVERT_XI_ETA_TO_THETA_PHI

! -------------------------------------------------------------------------------------------

SUBROUTINE CONVERT_THETA_PHI_TO_XI_ETA(a_theta, a_phi, a_xi, a_eta)

 implicit none

 integer i,j,k
 real*8, dimension(nx, ny, 6) :: a_xi, a_eta, a_theta, a_phi
 real*8 factor


 do k=1,4
  do j=1,ny
   do i=1,nx
    a_xi(i,j,k) = CD(i,j)/sqrt_metric_delta(i,j) * a_phi(i,j,k)
    a_eta(i,j,k) = - a_theta(i,j,k) +  XY(i,j)/sqrt_metric_delta(i,j) * a_phi(i,j,k)
   enddo
  enddo
 enddo

 do j=1,ny
  do i=1,nx
   factor = 1./(metric_delta(i,j) - 1.0)**0.5

   a_xi(i,j,5) = factor * (metric_D(j) * tanxi(i) * a_theta(i,j,5)  - &
                        metric_D(j) * taneta(j)/sqrt_metric_delta(i,j) * a_phi(i,j,5))

   a_eta(i,j,5) = factor * (metric_C(i) * taneta(j) * a_theta(i,j,5)  + &
                        metric_C(i) * tanxi(i)/sqrt_metric_delta(i,j) * a_phi(i,j,5))

   a_xi(i,j,6) = factor * (-metric_D(j) * tanxi(i) * a_theta(i,j,6)  + &
                        metric_D(j) * taneta(j)/sqrt_metric_delta(i,j) * a_phi(i,j,6))

   a_eta(i,j,6) = - factor * (metric_C(i) * taneta(j) * a_theta(i,j,6)  + &
                        metric_C(i) * tanxi(i)/sqrt_metric_delta(i,j) * a_phi(i,j,6))

  enddo
 enddo


END SUBROUTINE CONVERT_THETA_PHI_TO_XI_ETA

!------------------------------------------------------------------------------------------
SUBROUTINE TRANSFORM_VECTOR_BETWEEN_MESHES 
 implicit none

 real*8, dimension(nx,noverlap,6) :: interpolated_north_theta, interpolated_north_phi
 real*8, dimension(nx,noverlap,6) :: interpolated_south_theta, interpolated_south_phi
 real*8, dimension(noverlap,ny,6) :: interpolated_east_theta,  interpolated_east_phi
 real*8, dimension(noverlap,ny,6) :: interpolated_west_theta,  interpolated_west_phi

 COMMON /DIVERGENCE_POLAR/ interpolated_north_theta, interpolated_north_phi, &
                           interpolated_south_theta, interpolated_south_phi, &
                           interpolated_east_theta,  interpolated_east_phi, &
                           interpolated_west_theta,  interpolated_west_phi
 
  real*8, dimension(nx,noverlap,6) :: interpolated_north_xi, interpolated_north_eta
  real*8, dimension(nx,noverlap,6) :: interpolated_south_xi, interpolated_south_eta
  real*8, dimension(noverlap,ny,6) :: interpolated_east_xi,  interpolated_west_xi
  real*8, dimension(noverlap,ny,6) :: interpolated_east_eta, interpolated_west_eta

  COMMON /DIVERGENCE_CURVILINEAR/  interpolated_north_xi, interpolated_north_eta,&
                                   interpolated_south_xi, interpolated_south_eta,&
                                   interpolated_east_xi,  interpolated_west_xi,&
                                   interpolated_east_eta, interpolated_west_eta

 integer i,j,m,n
 real*8, dimension(-noverlap:noverlap) :: X, Y, C, D 

 do i=1,noverlap
   X(-i)     = tan(xi(1)-i*dxi)
   Y(-i)     = tan(eta(1)-i*deta)
   C(-i)     = sqrt(1D0 + X(-i)*X(-i))
   D(-i)     = sqrt(1D0 + Y(-i)*Y(-i))
 
   X(i)     = tan(xi(nx)+i*dxi)
   Y(i)     = tan(eta(ny)+i*deta)
   C(i)     = sqrt(1D0 + X(i)*X(i))
   D(i)     = sqrt(1D0 + Y(i)*Y(i))
 enddo

! vec_phi and vec_theta are the interpolated values from face m, i.e the face from which the vector is being 
! interpolated onto face n
! The vector transformation matrix is evaluated on n, i.e the face onto which the vector is being interpolated 

 do n=1,6
  do m = 1,6 
    if (abs(interfaces(n,m)) .eq. 1) then

      if(interfaces(m,n) .eq. +1)then  !(1,2), (2,3), (3,4), (4,1)
        do i=1,noverlap
          do j = 1,ny
            interpolated_east_xi(i,j,m) = C(-i)*metric_D(j)/sqrt(1D0+X(-i)*X(-i)+tan(eta(j))*tan(eta(j))) * interpolated_east_phi(i,j,m) 
          enddo
        enddo

      elseif(interfaces(m,n) .eq. -1)then !(2,1), (3,2), (4,3), (1,4)
        do i=1,noverlap
          do j = 1,ny
            interpolated_west_xi(i,j,m) = C(i)*metric_D(j)/sqrt(1D0+X(i)*X(i)+tan(eta(j))*tan(eta(j))) * interpolated_west_phi(i,j,m)
          enddo
        enddo
      endif
    endif
  enddo
 enddo

! There has got to be a more elegant way to do this....
! Interpolating from face 5 onto 1,2,3,4. Metrics evaluated on 1,2,3,4
! (5,1), (5,2), (5,3), (5,4)
  
 do j=1,noverlap
   do i=1,nx   
     interpolated_south_eta(i,j,5) = -interpolated_south_theta(i,j,5)&
                                   + tan(xi(i))*Y(j)/sqrt(1D0+tan(xi(i))*tan(xi(i))+Y(j)*Y(j)) * interpolated_south_phi(i,j,5)

     interpolated_east_eta(j,i,5) = -interpolated_east_theta(j,i,5)&
                                   + tan(xi(i))*Y(j)/sqrt(1D0+tan(xi(i))*tan(xi(i))+Y(j)*Y(j)) * interpolated_east_phi(j,i,5)
  
     interpolated_north_eta(i,j,5) = -interpolated_north_theta(i,j,5)&
                                   + tan(xi(i))*Y(j)/sqrt(1D0+tan(xi(i))*tan(xi(i))+Y(j)*Y(j)) * interpolated_north_phi(i,j,5)

     interpolated_west_eta(j,i,5)  = -interpolated_west_theta(j,i,5)&
                                   +  tan(xi(i))*Y(j)/sqrt(1D0+tan(xi(i))*tan(xi(i))+Y(j)*Y(j)) * interpolated_west_phi(j,i,5)
   enddo
 enddo
 
! Interpolating from face 6 onto 1,2,3,4. Metrics evaluated on 1,2,3,4
! (6,1), (6,2), (6,3), (6,4)
 do j=1,noverlap
   do i=1,nx
     interpolated_north_eta(i,j,6) = -interpolated_north_theta(i,j,6)&
                                   + tan(xi(i))*Y(-j)/sqrt(1D0+tan(xi(i))*tan(xi(i))+Y(-j)*Y(-j)) * interpolated_north_phi(i,j,6) 

     interpolated_east_eta(j,i,6)  = -interpolated_east_theta(j,i,6)&
                                   + tan(xi(i))*Y(-j)/sqrt(1D0+tan(xi(i))*tan(xi(i))+Y(-j)*Y(-j)) * interpolated_east_phi(j,i,6)
 
     interpolated_south_eta(i,j,6) = -interpolated_south_theta(i,j,6)&
                                   + tan(xi(i))*Y(-j)/sqrt(1D0+tan(xi(i))*tan(xi(i))+Y(-j)*Y(-j)) * interpolated_south_phi(i,j,6)

     interpolated_west_eta(j,i,6)  = -interpolated_west_theta(j,i,6)&
                                   + tan(xi(i))*Y(-j)/sqrt(1D0+tan(xi(i))*tan(xi(i))+Y(-j)*Y(-j)) * interpolated_west_phi(j,i,6)
   enddo
 enddo

! Faces 1,2,3,4 are now completely done for both xi and eta derivatives
! Interpolating from faces 1,2,3,4 onto face 5. Metric evaluated on face 5
! (1,5), (2,5), (3,5), (4,5)
  do j=1,noverlap
    do i=1,nx
      interpolated_north_eta(i,j,1) = metric_C(i)*Y(-j)/sqrt(tan(xi(i))*tan(xi(i))+Y(-j)*Y(-j)) * interpolated_north_theta(i,j,1)&
                                    + metric_C(i)/sqrt(tan(xi(i))*tan(xi(i))+Y(-j)*Y(-j)) * tan(xi(i))/sqrt(1D0+tan(xi(i))*tan(xi(i))+Y(-j)*Y(-j)) * &
                                      interpolated_north_phi(i,j,1)

      interpolated_north_xi(i,j,2)  = metric_D(i)*X(j)/sqrt(X(j)*X(j)+tan(eta(i))*tan(eta(i))) * interpolated_north_theta(i,j,2)&
                                    - metric_D(i)/sqrt(X(j)*X(j)+tan(eta(i))*tan(eta(i))) * tan(eta(i))/sqrt(1D0+X(j)*X(j)+tan(eta(i))*tan(eta(i))) * &
                                      interpolated_north_phi(i,j,2)

      interpolated_north_eta(i,j,3) = metric_C(i)*Y(j)/sqrt(tan(xi(i))*tan(xi(i))+Y(j)*Y(j)) * interpolated_north_theta(i,j,3)&
                                    + metric_C(i)/sqrt(tan(xi(i))*tan(xi(i))+Y(j)*Y(j)) * tan(xi(i))/sqrt(1D0+tan(xi(i))*tan(xi(i))+Y(j)*Y(j)) * &
                                      interpolated_north_phi(i,j,3)
 
      interpolated_north_xi(i,j,4)  = metric_D(i)*X(-j)/sqrt(X(-j)*X(-j)+tan(eta(i))*tan(eta(i))) * interpolated_north_theta(i,j,4)&
                                    - metric_D(i)/sqrt(X(-j)*X(-j)+tan(eta(i))*tan(eta(i))) * tan(eta(i))/sqrt(1D0+X(-j)*X(-j)+tan(eta(i))*tan(eta(i))) * &
                                      interpolated_north_phi(i,j,4)
    enddo
  enddo

! Interpolating from faces 1,2,3,4 onto face 6. Metric evaluated on face 6
! (1,6), (2,6), (3,6), (4,6)
  do j=1,noverlap
    do i=1,nx
      interpolated_south_eta(i,j,1) = -metric_C(i)*Y(j)/sqrt(tan(xi(i))*tan(xi(i))+Y(j)*Y(j)) * interpolated_south_theta(i,j,1)&
                                    -  metric_C(i)/sqrt(tan(xi(i))*tan(xi(i))+Y(j)*Y(j)) * tan(xi(i))/sqrt(1D0+tan(xi(i))*tan(xi(i))+Y(j)*Y(j)) * &
                                       interpolated_south_phi(i,j,1)

      interpolated_south_xi(i,j,2)  = -metric_D(i)*X(j)/sqrt(X(j)*X(j)+tan(eta(i))*tan(eta(i))) * interpolated_south_theta(i,j,2)&
                                    +  metric_D(i)/sqrt(X(j)*X(j)+tan(eta(i))*tan(eta(i))) * tan(eta(i))/sqrt(1D0+X(j)*X(j)+tan(eta(i))*tan(eta(i))) * &
                                       interpolated_south_phi(i,j,2)

      interpolated_south_eta(i,j,3) = -metric_C(i)*Y(-j)/sqrt(tan(xi(i))*tan(xi(i))+Y(-j)*Y(-j)) * interpolated_south_theta(i,j,3)&
                                    -  metric_C(i)/sqrt(tan(xi(i))*tan(xi(i))+Y(-j)*Y(-j)) * tan(xi(i))/sqrt(1D0+tan(xi(i))*tan(xi(i))+Y(-j)*Y(-j)) * &
                                       interpolated_south_phi(i,j,3)

      interpolated_south_xi(i,j,4)  = -metric_D(i)*X(-j)/sqrt(X(-j)*X(-j)+tan(eta(i))*tan(eta(i))) * interpolated_south_theta(i,j,4)&
                                    +  metric_D(i)/sqrt(X(-j)*X(-j)+tan(eta(i))*tan(eta(i))) * tan(eta(i))/sqrt(1D0+X(-j)*X(-j)+tan(eta(i))*tan(eta(i))) * &
                                       interpolated_south_phi(i,j,4)
    enddo
  enddo

END SUBROUTINE TRANSFORM_VECTOR_BETWEEN_MESHES
!------------------------------------------------------------------------------------------


END MODULE TRANSFORM
