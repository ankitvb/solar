MODULE METRIC

USE INIT

CONTAINS

!------------------------------------------------------------------------------------------------
SUBROUTINE SETUP_COORDINATES_METRICS_ALL

  implicit none

  integer i,j

  do i=1,nx
   xi(i)  = -pi*0.25 + (DBLE(i)-1.0)* pi * 0.5/(DBLE(nx) - 1.0)
   eta(i) = -pi*0.25 + (DBLE(i)-1.0)* pi * 0.5/(DBLE(nx) - 1.0)
  enddo

  tanxi = tan(xi)
  taneta = tan(eta)

  dxi = xi(2) - xi(1)
  deta = dxi

  metric_C = (1.0 + tanxi**2.)**0.5
  metric_D = (1.0 + taneta**2.)**0.5

  do j=1,ny
   do i=1,nx
    metric_delta(i,j) = (1. + tanxi(i)**2. + taneta(j)**2.)
    sqrt_metric_delta(i,j) = (1. + tanxi(i)**2. + taneta(j)**2.)**(0.5)

    XY(i,j) = tanxi(i) * taneta(j)
    CD(i,j) = metric_C(i) * metric_D(j)
   enddo
  enddo

  ! INTERPOLATION POINTS:

  do i = 1,noverlap
   interpoints(:, i) = atan(taneta/tan(xi(nx) + i*dxi))
  enddo

  ! DETERMINING SPHERICAL COORDINATES

  !REGION 1

  do j=1,nx
   do i=1,nx
    phi(i,j,1) = xi(i) !+ pi
    theta(i,j,1) = atan(1./(taneta(j) * cos(phi(i,j,1))))
    if (theta(i,j,1) < 0) theta(i,j,1) = theta(i,j,1) + pi
!    if (phi(i,j,1) < 0) phi(i,j,1) = phi(i,j,1) + 2.*pi
   enddo
  enddo 

  !REGION 2

  ! phi(:,:,2) = phi(:,:,1) + pi*0.5
  ! theta(:,:,2) = theta(:,:,1)


  do j=1,nx
   do i=1,nx
    phi(i,j,2) = xi(i) + pi*0.5
    theta(i,j,2) = atan(1./(taneta(j) * sin(phi(i,j,2))))
    if (theta(i,j,2) < 0) theta(i,j,2) = theta(i,j,2) + pi
!    if (phi(i,j,2) < 0) phi(i,j,2) = phi(i,j,2) + 2.*pi

   enddo
  enddo 

  !REGION 3

  ! phi(:,:,3) = phi(:,:,2) + pi*0.5
  ! theta(:,:,3) = theta(:,:,2)

  do j=1,nx
   do i=1,nx
    phi(i,j,3) = xi(i) + pi
    theta(i,j,3) = atan(-1./(taneta(j) * cos(phi(i,j,3))))
    if (theta(i,j,3) < 0) theta(i,j,3) = theta(i,j,3) + pi
!    if (phi(i,j,3) < 0) phi(i,j,3) = phi(i,j,3) + 2.*pi

   enddo
  enddo 

  !REGION 4

  ! phi(:,:,4) = phi(:,:,3) + pi*0.5
  ! theta(:,:,4) = theta(:,:,3)

  do j=1,nx
   do i=1,nx
    phi(i,j,4) = xi(i) + pi*1.5
    theta(i,j,4) = atan(-1./(taneta(j) * sin(phi(i,j,4))))
    if (theta(i,j,4) < 0) theta(i,j,4) = theta(i,j,4) + pi
!    if (phi(i,j,4) < 0) phi(i,j,4) = phi(i,j,4) + 2.*pi
   enddo
  enddo 

  !REGION 5

  do j=1,nx
   do i=1,nx
    theta(i,j,5) = atan((tanxi(i)**2. + taneta(j)**2.)**0.5)
    phi(i,j,5) = atan2(tanxi(i),-taneta(j))! - pi*0.25 !+ pi
    if (phi(i,j,5) < -pi*0.25) phi(i,j,5) = phi(i,j,5) + 2.*pi
!    if (phi(i,j,5) > 1.75*pi) phi(i,j,5) = phi(i,j,5) - 2.*pi

   enddo
  enddo 


  !REGION 6

  do j=1,nx
   do i=1,nx
    theta(i,j,6) = pi - atan((tanxi(i)**2. + taneta(j)**2.)**0.5)
    phi(i,j,6) = atan2(tanxi(i),taneta(j)) !- pi*0.25!+ pi
    if (phi(i,j,6) < -pi*0.25) phi(i,j,6) = phi(i,j,6) + 2.*pi
!    if (phi(i,j,6) > 1.75*pi) phi(i,j,6) = phi(i,j,6) - 2.*pi

   enddo
  enddo 



END SUBROUTINE SETUP_COORDINATES_METRICS_ALL
!----------------------------------------------------------------------------------

END MODULE METRIC
