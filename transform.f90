MODULE TRANSFORM

USE INIT

CONTAINS

!--------------------------------------------------------------------------------------------  
SUBROUTINE CONVERT_XI_ETA_TO_THETA_PHI(a_xi, a_eta, a_theta, a_phi)

 implicit none

 integer i,j,k
 real*8, dimension(mx,my,mr) :: a_xi, a_eta, a_theta, a_phi
 real*8 factor

 a_phi = 0.0d0
 a_theta = 0.0d0

 if((zoneid.ge.1).and.(zoneid.le.4))then
  do k=1,mr
   do j=1,my
    do i=1,mx
     a_theta(i,j,k) = XY(i,j)/CD(i,j) * a_xi(i,j,k) - a_eta(i,j,k)
     a_phi(i,j,k) = sqrt_metric_delta(i,j) / CD(i,j) * a_xi(i,j,k)
    enddo
   enddo
  enddo
 endif

 if(zoneid.eq.5)then
  do k=1,mr
   do j=1,my
    do i=1,mx
     factor = 1.0D0/(metric_delta(i,j) - 1.0D0)**0.5D0
     a_theta(i,j,k) =  factor*( tanxi(i)/metric_D(j)* a_xi(i,j,k) + &
                              taneta(j)/metric_C(i) * a_eta(i,j,k) )
     a_phi(i,j,k) = factor*( -sqrt_metric_delta(i,j) * taneta(j)/metric_D(j) * a_xi(i,j,k) + &
                            sqrt_metric_delta(i,j) * tanxi(i)/metric_C(i) * a_eta(i,j,k))
    enddo
   enddo
  enddo
 endif

 if(zoneid.eq.6)then
  do k=1,mr
   do j=1,my
    do i=1,mx
     factor = 1.0D0/(metric_delta(i,j) - 1.0D0)**0.5D0
     a_theta(i,j,k) =  -factor*( tanxi(i)/metric_D(j)* a_xi(i,j,k) + &
                                 taneta(j)/metric_C(i) * a_eta(i,j,k) )
     a_phi(i,j,k) = factor*( sqrt_metric_delta(i,j) * taneta(j)/metric_D(j) * a_xi(i,j,k) - &
                             sqrt_metric_delta(i,j) * tanxi(i)/metric_C(i) * a_eta(i,j,k))
    enddo
   enddo
  enddo
 endif


END SUBROUTINE CONVERT_XI_ETA_TO_THETA_PHI

! -------------------------------------------------------------------------------------------

SUBROUTINE CONVERT_THETA_PHI_TO_XI_ETA(a_theta, a_phi, a_xi, a_eta)

 implicit none

 integer i,j,k
 real*8, dimension(mx,my,mr) :: a_xi, a_eta, a_theta, a_phi
 real*8 factor

 if((zoneid.ge.1).and.(zoneid.le.4))then
  do k=1,mr
   do j=1,my
    do i=1,mx
     a_xi(i,j,k) = CD(i,j)/sqrt_metric_delta(i,j) * a_phi(i,j,k)
     a_eta(i,j,k) = - a_theta(i,j,k) +  XY(i,j)/sqrt_metric_delta(i,j) * a_phi(i,j,k)
    enddo
   enddo
  enddo
 endif

 if(zoneid.eq.5)then
 do k=1,mr
  do j=1,my
   do i=1,mx
    factor = 1./(metric_delta(i,j) - 1.0)**0.5
    a_xi(i,j,k) = factor * (metric_D(j) * tanxi(i) * a_theta(i,j,k)  - &
                          metric_D(j) * taneta(j)/sqrt_metric_delta(i,j) * a_phi(i,j,k))
    a_eta(i,j,k) = factor * (metric_C(i) * taneta(j) * a_theta(i,j,k)  + &
                           metric_C(i) * tanxi(i)/sqrt_metric_delta(i,j) * a_phi(i,j,k))
    enddo
   enddo
  enddo
 endif  

 if(zoneid.eq.6)then
  do k=1,mr
   do j=1,my
    do i=1,mx
     factor = 1./(metric_delta(i,j) - 1.0)**0.5
     a_xi(i,j,k) = factor * (-metric_D(j) * tanxi(i) * a_theta(i,j,k)  + &
                            metric_D(j) * taneta(j)/sqrt_metric_delta(i,j) * a_phi(i,j,k))
     a_eta(i,j,k) = - factor * (metric_C(i) * taneta(j) * a_theta(i,j,k)  + &
                              metric_C(i) * tanxi(i)/sqrt_metric_delta(i,j) * a_phi(i,j,k))
    enddo
   enddo
  enddo
 endif


END SUBROUTINE CONVERT_THETA_PHI_TO_XI_ETA

!------------------------------------------------------------------------------------------
SUBROUTINE TRANSFORM_VECTOR_BETWEEN_MESHES_DIV 
 implicit none

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

 integer i,j,k,m,n,ii,jj
 real*8, dimension(-noverlap:noverlap) :: X, Y, C, D 
 real*8, dimension(nx)                 :: global_C, global_D


 do i=1,noverlap
   X(-i)     = tan(gxi(1)-i*dxi)
   Y(-i)     = tan(geta(1)-i*deta)
   C(-i)     = sqrt(1D0 + X(-i)*X(-i))
   D(-i)     = sqrt(1D0 + Y(-i)*Y(-i))
 
   X(i)     = tan(gxi(nx)+i*dxi)
   Y(i)     = tan(geta(ny)+i*deta)
   C(i)     = sqrt(1D0 + X(i)*X(i))
   D(i)     = sqrt(1D0 + Y(i)*Y(i))
 enddo

 global_C = (1.0D0 + tan(gxi)**2.0D0)**0.5D0
 global_D = (1.0D0 + tan(geta)**2.0D0)**0.5D0

! vec_phi and vec_theta are the interpolated values from face m, i.e the face from which the vector is being 
! interpolated onto face n
! The vector transformation matrix is evaluated on n, i.e the face onto which the vector is being interpolated 

 if((zoneid.ge.1).and.(zoneid.le.4))then

  do k=1,mr
   do i=1,noverlap
    do j = 1,my
     jj = yrank*my+j
     interpolated_east_xi(i,j,k) = C(-i)*global_D(jj)/sqrt(1D0+X(-i)*X(-i)+tan(geta(jj))*tan(geta(jj))) * interpolated_east_phi(i,j,k) 
    enddo
   enddo

   do i=1,noverlap
    do j = 1,my
     jj = yrank*my+j
     interpolated_west_xi(i,j,k) = C(i)*global_D(jj)/sqrt(1D0+X(i)*X(i)+tan(geta(jj))*tan(geta(jj))) * interpolated_west_phi(i,j,k)
    enddo
   enddo
  enddo 
   
 endif

! Interpolating from face 5 onto 1,2,3,4. Metrics evaluated on 1,2,3,4
! (5,1), (5,2), (5,3), (5,4)
  
 if(zoneid.eq.5)then
  do k=1,mr
   do j=1,noverlap
    do i=1,mx
      ii = yrank*my+i   
      interpolated_south_eta(i,j,k) = -interpolated_south_theta(i,j,k)&
                                    + tan(xi(i))*Y(j)/sqrt(1D0+tan(xi(i))*tan(xi(i))+Y(j)*Y(j)) * interpolated_south_phi(i,j,k)
      interpolated_east_eta(j,i,k)  = -interpolated_east_theta(j,i,k)&
                                    + tan(gxi(ii))*Y(j)/sqrt(1D0+tan(gxi(ii))*tan(gxi(ii))+Y(j)*Y(j)) * interpolated_east_phi(j,i,k)
      interpolated_north_eta(i,j,k) = -interpolated_north_theta(i,j,k)&
                                    + tan(xi(i))*Y(j)/sqrt(1D0+tan(xi(i))*tan(xi(i))+Y(j)*Y(j)) * interpolated_north_phi(i,j,k)
      interpolated_west_eta(j,i,k)  = -interpolated_west_theta(j,i,k)&
                                    +  tan(gxi(ii))*Y(j)/sqrt(1D0+tan(gxi(ii))*tan(gxi(ii))+Y(j)*Y(j)) * interpolated_west_phi(j,i,k)
    enddo
   enddo
  enddo
 endif
 
! Interpolating from face 6 onto 1,2,3,4. Metrics evaluated on 1,2,3,4
! (6,1), (6,2), (6,3), (6,4)

 if(zoneid.eq.6)then
  do k=1,mr
   do j=1,noverlap
    do i=1,mx
      ii = yrank*my+i
      interpolated_north_eta(i,j,k) = -interpolated_north_theta(i,j,k)&
                                    + tan(xi(i))*Y(-j)/sqrt(1D0+tan(xi(i))*tan(xi(i))+Y(-j)*Y(-j)) * interpolated_north_phi(i,j,k) 
      interpolated_east_eta(j,i,k)  = -interpolated_east_theta(j,i,k)&
                                    + tan(gxi(ii))*Y(-j)/sqrt(1D0+tan(gxi(ii))*tan(gxi(ii))+Y(-j)*Y(-j)) * interpolated_east_phi(j,i,k)
      interpolated_south_eta(i,j,k) = -interpolated_south_theta(i,j,k)&
                                    + tan(xi(i))*Y(-j)/sqrt(1D0+tan(xi(i))*tan(xi(i))+Y(-j)*Y(-j)) * interpolated_south_phi(i,j,k)
      interpolated_west_eta(j,i,k)  = -interpolated_west_theta(j,i,k)&
                                    + tan(gxi(ii))*Y(-j)/sqrt(1D0+tan(gxi(ii))*tan(gxi(ii))+Y(-j)*Y(-j)) * interpolated_west_phi(j,i,k)
    enddo
   enddo
  enddo
 endif

! Faces 1,2,3,4 are now completely done for both xi and eta derivatives
! Interpolating from faces 1,2,3,4 onto face 5. Metric evaluated on face 5
! (1,5), (2,5), (3,5), (4,5)
  do k=1,mr
   do j=1,noverlap
    do i=1,mx
      ii = xrank*mx+i
      if(zoneid.eq.1)&
      interpolated_north_eta(i,j,k) = metric_C(i)*Y(-j)/sqrt(tan(xi(i))*tan(xi(i))+Y(-j)*Y(-j)) * interpolated_north_theta(i,j,k)&
                                  + metric_C(i)/sqrt(tan(xi(i))*tan(xi(i))+Y(-j)*Y(-j)) * tan(xi(i))/sqrt(1D0+tan(xi(i))*tan(xi(i))+Y(-j)*Y(-j)) * &
                                    interpolated_north_phi(i,j,k)
      if(zoneid.eq.2)&
      interpolated_north_xi(i,j,k)  = global_D(ii)*X(j)/sqrt(X(j)*X(j)+tan(geta(ii))*tan(geta(ii))) * interpolated_north_theta(i,j,k)&
                                  - global_D(ii)/sqrt(X(j)*X(j)+tan(geta(ii))*tan(geta(ii))) * tan(geta(ii))/sqrt(1D0+X(j)*X(j)+tan(geta(ii))*tan(geta(ii))) * &
                                    interpolated_north_phi(i,j,k)
      if(zoneid.eq.3)&
      interpolated_north_eta(i,j,k) = metric_C(i)*Y(j)/sqrt(tan(xi(i))*tan(xi(i))+Y(j)*Y(j)) * interpolated_north_theta(i,j,k)&
                                  + metric_C(i)/sqrt(tan(xi(i))*tan(xi(i))+Y(j)*Y(j)) * tan(xi(i))/sqrt(1D0+tan(xi(i))*tan(xi(i))+Y(j)*Y(j)) * &
                                    interpolated_north_phi(i,j,k)
      if(zoneid.eq.4)&
      interpolated_north_xi(i,j,k)  = global_D(ii)*X(-j)/sqrt(X(-j)*X(-j)+tan(geta(ii))*tan(geta(ii))) * interpolated_north_theta(i,j,k)&
                                  - global_D(ii)/sqrt(X(-j)*X(-j)+tan(geta(ii))*tan(geta(ii))) * tan(geta(ii))/sqrt(1D0+X(-j)*X(-j)+tan(geta(ii))*tan(geta(ii))) * &
                                    interpolated_north_phi(i,j,k)
    enddo
   enddo
  enddo

! Interpolating from faces 1,2,3,4 onto face 6. Metric evaluated on face 6
! (1,6), (2,6), (3,6), (4,6)
  do k=1,mr
   do j=1,noverlap
    do i=1,mx
      ii = xrank*mx+i
      if(zoneid.eq.1)&
      interpolated_south_eta(i,j,k) = -metric_C(i)*Y(j)/sqrt(tan(xi(i))*tan(xi(i))+Y(j)*Y(j)) * interpolated_south_theta(i,j,k)&
                                  -  metric_C(i)/sqrt(tan(xi(i))*tan(xi(i))+Y(j)*Y(j)) * tan(xi(i))/sqrt(1D0+tan(xi(i))*tan(xi(i))+Y(j)*Y(j)) * &
                                     interpolated_south_phi(i,j,k)
      if(zoneid.eq.2)&
      interpolated_south_xi(i,j,k)  = -global_D(ii)*X(j)/sqrt(X(j)*X(j)+tan(geta(ii))*tan(geta(ii))) * interpolated_south_theta(i,j,k)&
                                  +  global_D(ii)/sqrt(X(j)*X(j)+tan(geta(ii))*tan(geta(ii))) * tan(geta(ii))/sqrt(1D0+X(j)*X(j)+tan(geta(ii))*tan(geta(ii))) * &
                                     interpolated_south_phi(i,j,k)
      if(zoneid.eq.3)&
      interpolated_south_eta(i,j,k) = -metric_C(i)*Y(-j)/sqrt(tan(xi(i))*tan(xi(i))+Y(-j)*Y(-j)) * interpolated_south_theta(i,j,k)&
                                  -  metric_C(i)/sqrt(tan(xi(i))*tan(xi(i))+Y(-j)*Y(-j)) * tan(xi(i))/sqrt(1D0+tan(xi(i))*tan(xi(i))+Y(-j)*Y(-j)) * &
                                     interpolated_south_phi(i,j,k)
      if(zoneid.eq.4)&
      interpolated_south_xi(i,j,k)  = -global_D(ii)*X(-j)/sqrt(X(-j)*X(-j)+tan(geta(ii))*tan(geta(ii))) * interpolated_south_theta(i,j,k)&
                                  +  global_D(ii)/sqrt(X(-j)*X(-j)+tan(geta(ii))*tan(geta(ii))) * tan(geta(ii))/sqrt(1D0+X(-j)*X(-j)+tan(geta(ii))*tan(geta(ii))) * &
                                     interpolated_south_phi(i,j,k)
    enddo
   enddo
  enddo

END SUBROUTINE TRANSFORM_VECTOR_BETWEEN_MESHES_DIV
!------------------------------------------------------------------------------------------

SUBROUTINE TRANSFORM_VECTOR_BETWEEN_MESHES_CURL
 implicit none

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

 integer i,j,k,m,n,ii,jj
 real*8, dimension(-noverlap:noverlap) :: X, Y, C, D
 real*8, dimension(nx)                 :: global_C, global_D


 do i=1,noverlap
   X(-i)     = tan(gxi(1)-i*dxi)
   Y(-i)     = tan(geta(1)-i*deta)
   C(-i)     = sqrt(1D0 + X(-i)*X(-i))
   D(-i)     = sqrt(1D0 + Y(-i)*Y(-i))

   X(i)     = tan(gxi(nx)+i*dxi)
   Y(i)     = tan(geta(ny)+i*deta)
   C(i)     = sqrt(1D0 + X(i)*X(i))
   D(i)     = sqrt(1D0 + Y(i)*Y(i))
 enddo

 global_C = (1.0D0 + tan(gxi)**2.0D0)**0.5D0
 global_D = (1.0D0 + tan(geta)**2.0D0)**0.5D0

! vec_phi and vec_theta are the interpolated values from face m, i.e the face from which the vector is being
! interpolated onto face n
! The vector transformation matrix is evaluated on n, i.e the face onto which the vector is being interpolated

 if((zoneid.ge.1).and.(zoneid.le.4))then

  do k=1,mr
   do i=1,noverlap
    do j = 1,my
     jj = yrank*my+j
     interpolated_east_xi(i,j,k)  = C(-i)*global_D(jj)/sqrt(1D0+X(-i)*X(-i)+tan(geta(jj))*tan(geta(jj))) * interpolated_east_phi(i,j,k)
     interpolated_east_eta(i,j,k) = -interpolated_east_theta(i,j,k)&
                                  + X(-i)*tan(geta(jj))/sqrt(1D0+X(-i)*X(-i)+tan(geta(jj))*tan(geta(jj))) * interpolated_east_phi(i,j,k)
    enddo
   enddo

   do i=1,noverlap
    do j = 1,my
     jj = yrank*my+j
     interpolated_west_xi(i,j,k)  = C(i)*global_D(jj)/sqrt(1D0+X(i)*X(i)+tan(geta(jj))*tan(geta(jj))) * interpolated_west_phi(i,j,k)
     interpolated_west_eta(i,j,k) = -interpolated_west_theta(i,j,k)&
                                  + X(i)*tan(geta(jj))/sqrt(1D0+X(i)*X(i)+tan(geta(jj))*tan(geta(jj))) * interpolated_west_phi(i,j,k)
    enddo
   enddo
  enddo

 endif

! Interpolating from face 5 onto 1,2,3,4. Metrics evaluated on 1,2,3,4
! (5,1), (5,2), (5,3), (5,4)

 if(zoneid.eq.5)then
  do k=1,mr
   do j=1,noverlap
    do i=1,mx
      ii = yrank*my+i
      interpolated_south_eta(i,j,k) = -interpolated_south_theta(i,j,k)&
                                    + tan(xi(i))*Y(j)/sqrt(1D0+tan(xi(i))*tan(xi(i))+Y(j)*Y(j)) * interpolated_south_phi(i,j,k)
      interpolated_south_xi(i,j,k)  = metric_C(i)*D(j)/sqrt(1D0+tan(xi(i))*tan(xi(i))+Y(j)*Y(j)) * interpolated_south_phi(i,j,k)
  
      interpolated_east_eta(j,i,k)  = -interpolated_east_theta(j,i,k)&
                                    + tan(gxi(ii))*Y(j)/sqrt(1D0+tan(gxi(ii))*tan(gxi(ii))+Y(j)*Y(j)) * interpolated_east_phi(j,i,k)
      interpolated_east_xi(j,i,k)   = global_C(ii)*D(j)/sqrt(1D0+tan(gxi(ii))*tan(gxi(ii))+Y(j)*Y(j)) * interpolated_east_phi(j,i,k)

      interpolated_north_eta(i,j,k) = -interpolated_north_theta(i,j,k)&
                                    + tan(xi(i))*Y(j)/sqrt(1D0+tan(xi(i))*tan(xi(i))+Y(j)*Y(j)) * interpolated_north_phi(i,j,k)
      interpolated_north_xi(i,j,k)  = metric_C(i)*D(j)/sqrt(1D0+tan(xi(i))*tan(xi(i))+Y(j)*Y(j)) * interpolated_north_phi(i,j,k)

      interpolated_west_eta(j,i,k)  = -interpolated_west_theta(j,i,k)&
                                    + tan(gxi(ii))*Y(j)/sqrt(1D0+tan(gxi(ii))*tan(gxi(ii))+Y(j)*Y(j)) * interpolated_west_phi(j,i,k)
      interpolated_west_xi(j,i,k)   = global_C(ii)*D(j)/sqrt(1D0+tan(gxi(ii))*tan(gxi(ii))+Y(j)*Y(j)) * interpolated_west_phi(j,i,k)
    enddo
   enddo
  enddo
 endif

! Interpolating from face 6 onto 1,2,3,4. Metrics evaluated on 1,2,3,4
! (6,1), (6,2), (6,3), (6,4)

 if(zoneid.eq.6)then
  do k=1,mr
   do j=1,noverlap
    do i=1,mx
      ii = yrank*my+i
      interpolated_north_eta(i,j,k) = -interpolated_north_theta(i,j,k)&
                                    + tan(xi(i))*Y(-j)/sqrt(1D0+tan(xi(i))*tan(xi(i))+Y(-j)*Y(-j)) * interpolated_north_phi(i,j,k)
      interpolated_north_xi(i,j,k)  = metric_C(i)*D(-j)/sqrt(1D0+tan(xi(i))*tan(xi(i))+Y(-j)*Y(-j)) * interpolated_north_phi(i,j,k)

      interpolated_east_eta(j,i,k)  = -interpolated_east_theta(j,i,k)&
                                    + tan(gxi(ii))*Y(-j)/sqrt(1D0+tan(gxi(ii))*tan(gxi(ii))+Y(-j)*Y(-j)) * interpolated_east_phi(j,i,k)
      interpolated_east_xi(j,i,k)   = global_C(ii)*D(-j)/sqrt(1D0+tan(gxi(ii))*tan(gxi(ii))+Y(-j)*Y(-j)) * interpolated_east_phi(j,i,k)

      interpolated_south_eta(i,j,k) = -interpolated_south_theta(i,j,k)&
                                    + tan(xi(i))*Y(-j)/sqrt(1D0+tan(xi(i))*tan(xi(i))+Y(-j)*Y(-j)) * interpolated_south_phi(i,j,k)
      interpolated_south_xi(i,j,k)  = metric_C(i)*D(-j)/sqrt(1D0+tan(xi(i))*tan(xi(i))+Y(-j)*Y(-j)) * interpolated_south_phi(i,j,k)

      interpolated_west_eta(j,i,k)  = -interpolated_west_theta(j,i,k)&
                                    + tan(gxi(ii))*Y(-j)/sqrt(1D0+tan(gxi(ii))*tan(gxi(ii))+Y(-j)*Y(-j)) * interpolated_west_phi(j,i,k)
      interpolated_west_xi(j,i,k)   = global_C(ii)*D(-j)/sqrt(1D0+tan(gxi(ii))*tan(gxi(ii))+Y(-j)*Y(-j)) * interpolated_west_phi(j,i,k)
    enddo
   enddo
  enddo
 endif

! Faces 1,2,3,4 are now completely done for both xi and eta derivatives
! Interpolating from faces 1,2,3,4 onto face 5. Metric evaluated on face 5
! (1,5), (2,5), (3,5), (4,5)
  do k=1,mr
   do j=1,noverlap
    do i=1,mx
      ii = xrank*mx+i
      if(zoneid.eq.1)then
       interpolated_north_eta(i,j,k) = metric_C(i)*Y(-j)/sqrt(tan(xi(i))*tan(xi(i))+Y(-j)*Y(-j)) * interpolated_north_theta(i,j,k)&
                                     + metric_C(i)/sqrt(tan(xi(i))*tan(xi(i))+Y(-j)*Y(-j)) * tan(xi(i))/sqrt(1D0+tan(xi(i))*tan(xi(i))+Y(-j)*Y(-j)) * &
                                       interpolated_north_phi(i,j,k)
       interpolated_north_xi(i,j,k)  = D(-j)*tan(xi(i))/sqrt(tan(xi(i))*tan(xi(i))+Y(-j)*Y(-j)) * interpolated_north_theta(i,j,k)&
                                     - D(-j)/sqrt(tan(xi(i))*tan(xi(i))+Y(-j)*Y(-j)) * Y(-j)/sqrt(1D0+tan(xi(i))*tan(xi(i))+Y(-j)*Y(-j)) * &
                                       interpolated_north_phi(i,j,k)
      endif

      if(zoneid.eq.2)then
       interpolated_north_xi(i,j,k)  = global_D(ii)*X(j)/sqrt(X(j)*X(j)+tan(geta(ii))*tan(geta(ii))) * interpolated_north_theta(i,j,k)&
                                     - global_D(ii)/sqrt(X(j)*X(j)+tan(geta(ii))*tan(geta(ii))) * tan(geta(ii))/sqrt(1D0+X(j)*X(j)+tan(geta(ii))*tan(geta(ii))) * &
                                       interpolated_north_phi(i,j,k)
       interpolated_north_eta(i,j,k) = C(j)*tan(geta(ii))/sqrt(X(j)*X(j)+tan(geta(ii))*tan(geta(ii))) * interpolated_north_theta(i,j,k)&
                                     + C(j)/sqrt(X(j)*X(j)+tan(geta(ii))*tan(geta(ii))) * X(j)/sqrt(1D0+X(j)*X(j)+tan(geta(ii))*tan(geta(ii))) * &
                                       interpolated_north_phi(i,j,k)
      endif

      if(zoneid.eq.3)then
       interpolated_north_eta(i,j,k) = metric_C(i)*Y(j)/sqrt(tan(xi(i))*tan(xi(i))+Y(j)*Y(j)) * interpolated_north_theta(i,j,k)&
                                     + metric_C(i)/sqrt(tan(xi(i))*tan(xi(i))+Y(j)*Y(j)) * tan(xi(i))/sqrt(1D0+tan(xi(i))*tan(xi(i))+Y(j)*Y(j)) * &
                                       interpolated_north_phi(i,j,k)
       interpolated_north_xi(i,j,k)  = D(j)*tan(xi(i))/sqrt(tan(xi(i))*tan(xi(i))+Y(j)*Y(j)) * interpolated_north_theta(i,j,k)&
                                     - D(j)/sqrt(tan(xi(i))*tan(xi(i))+Y(j)*Y(j)) * Y(j)/sqrt(1D0+tan(xi(i))*tan(xi(i))+Y(j)*Y(j)) * &
                                       interpolated_north_phi(i,j,k)
      endif

      if(zoneid.eq.4)then
       interpolated_north_xi(i,j,k)  = global_D(ii)*X(-j)/sqrt(X(-j)*X(-j)+tan(geta(ii))*tan(geta(ii))) * interpolated_north_theta(i,j,k)&
                                     - global_D(ii)/sqrt(X(-j)*X(-j)+tan(geta(ii))*tan(geta(ii))) * tan(geta(ii))/sqrt(1D0+X(-j)*X(-j)+tan(geta(ii))*tan(geta(ii))) * &
                                       interpolated_north_phi(i,j,k)
       interpolated_north_eta(i,j,k) = C(-j)*tan(geta(ii))/sqrt(X(-j)*X(-j)+tan(geta(ii))*tan(geta(ii))) * interpolated_north_theta(i,j,k)&
                                     + C(-j)/sqrt(X(-j)*X(-j)+tan(geta(ii))*tan(geta(ii))) * X(-j)/sqrt(1D0+X(-j)*X(-j)+tan(geta(ii))*tan(geta(ii))) * &
                                       interpolated_north_phi(i,j,k)
      endif
    enddo
   enddo
  enddo

! Interpolating from faces 1,2,3,4 onto face 6. Metric evaluated on face 6
! (1,6), (2,6), (3,6), (4,6)
  do k=1,mr
   do j=1,noverlap
    do i=1,mx
      ii = xrank*mx+i
      if(zoneid.eq.1)then
       interpolated_south_eta(i,j,k) = -metric_C(i)*Y(j)/sqrt(tan(xi(i))*tan(xi(i))+Y(j)*Y(j)) * interpolated_south_theta(i,j,k)&
                                     -  metric_C(i)/sqrt(tan(xi(i))*tan(xi(i))+Y(j)*Y(j)) * tan(xi(i))/sqrt(1D0+tan(xi(i))*tan(xi(i))+Y(j)*Y(j)) * &
                                        interpolated_south_phi(i,j,k)
       interpolated_south_xi(i,j,k)  = -D(j)*tan(xi(i))/sqrt(tan(xi(i))*tan(xi(i))+Y(j)*Y(j)) * interpolated_south_theta(i,j,k)&
                                     +  D(j)/sqrt(tan(xi(i))*tan(xi(i))+Y(j)*Y(j)) * Y(j)/sqrt(1D0+tan(xi(i))*tan(xi(i))+Y(j)*Y(j)) * &
                                        interpolated_south_phi(i,j,k)
      endif

      if(zoneid.eq.2)then
       interpolated_south_xi(i,j,k)  = -global_D(ii)*X(j)/sqrt(X(j)*X(j)+tan(geta(ii))*tan(geta(ii))) * interpolated_south_theta(i,j,k)&
                                     +  global_D(ii)/sqrt(X(j)*X(j)+tan(geta(ii))*tan(geta(ii))) * tan(geta(ii))/sqrt(1D0+X(j)*X(j)+tan(geta(ii))*tan(geta(ii))) * &
                                        interpolated_south_phi(i,j,k)
       interpolated_south_eta(i,j,k) = -C(j)*tan(geta(ii))/sqrt(X(j)*X(j)+tan(geta(ii))*tan(geta(ii))) * interpolated_south_theta(i,j,k)&
                                     -  C(j)/sqrt(X(j)*X(j)+tan(geta(ii))*tan(geta(ii))) * X(j)/sqrt(1D0+X(j)*X(j)+tan(geta(ii))*tan(geta(ii))) * &
                                        interpolated_south_phi(i,j,k)
      endif

      if(zoneid.eq.3)then
       interpolated_south_eta(i,j,k) = -metric_C(i)*Y(-j)/sqrt(tan(xi(i))*tan(xi(i))+Y(-j)*Y(-j)) * interpolated_south_theta(i,j,k)&
                                     -  metric_C(i)/sqrt(tan(xi(i))*tan(xi(i))+Y(-j)*Y(-j)) * tan(xi(i))/sqrt(1D0+tan(xi(i))*tan(xi(i))+Y(-j)*Y(-j)) * &
                                        interpolated_south_phi(i,j,k)
       interpolated_south_xi(i,j,k)  = -D(-j)*tan(xi(i))/sqrt(tan(xi(i))*tan(xi(i))+Y(-j)*Y(-j)) * interpolated_south_theta(i,j,k)&
                                     +  D(-j)/sqrt(tan(xi(i))*tan(xi(i))+Y(-j)*Y(-j)) * Y(-j)/sqrt(1D0+tan(xi(i))*tan(xi(i))+Y(-j)*Y(-j)) * &
                                        interpolated_south_phi(i,j,k)
      endif

      if(zoneid.eq.4)then
       interpolated_south_xi(i,j,k)  = -global_D(ii)*X(-j)/sqrt(X(-j)*X(-j)+tan(geta(ii))*tan(geta(ii))) * interpolated_south_theta(i,j,k)&
                                     +  global_D(ii)/sqrt(X(-j)*X(-j)+tan(geta(ii))*tan(geta(ii))) * tan(geta(ii))/sqrt(1D0+X(-j)*X(-j)+tan(geta(ii))*tan(geta(ii))) * &
                                        interpolated_south_phi(i,j,k)
       interpolated_south_eta(i,j,k) = -C(-j)*tan(geta(ii))/sqrt(X(-j)*X(-j)+tan(geta(ii))*tan(geta(ii))) * interpolated_south_theta(i,j,k)&
                                     -  C(-j)/sqrt(X(-j)*X(-j)+tan(geta(ii))*tan(geta(ii))) * X(-j)/sqrt(1D0+X(-j)*X(-j)+tan(geta(ii))*tan(geta(ii))) * &
                                        interpolated_south_phi(i,j,k)
      endif
    enddo
   enddo
  enddo


END SUBROUTINE TRANSFORM_VECTOR_BETWEEN_MESHES_CURL

!-----------------------------------------------------------------------------------------

END MODULE TRANSFORM
