MODULE METRIC

USE INIT

CONTAINS

!------------------------------------------------------------------------------------------------
SUBROUTINE SETUP_COORDINATES_METRICS_ALL

  implicit none

  integer i,j,k,ii,jj,kk,iter

  do i=1,nx
   gxi(i)  = -pi*0.25D0 + (DBLE(i)-1.0D0)* pi * 0.5D0/(DBLE(nx)-1.0D0)
   geta(i) = -pi*0.25D0 + (DBLE(i)-1.0D0)* pi * 0.5D0/(DBLE(nx)-1.0D0)
  enddo
  do i=1,nr
!   gr(i)   = 0.0D0 + (DBLE(i)-1.0D0) * 2D0*pi/(DBLE(nr)-1.0D0)
   gchi(i) = 0.0D0 + (DBLE(i)-1.0D0) * 1D0/(DBLE(nr)-1.0D0)
  enddo

  do i=1,mx
   xi(i) = gxi(xrank*mx+i)
  enddo
  do i=1,my
   eta(i) = geta(yrank*my+i)
  enddo
  do i=1,mr
!   r(i)   = gr(rrank*mr+i)
   chi(i) = gchi(rrank*mr+i)
  enddo
 
! Nullifying the r effect for 2D case
  if(nr.eq.1)then
   r(1)  = 1.0D0
   gr(1) = 1.0D0
  endif
 
  tanxi = tan(xi)
  taneta = tan(eta)

  dxi  = xi(2) - xi(1)
  deta = dxi
  dchi = chi(2) - chi(1)

! Computing the metric for radial direction
 do kk=nst+1,nr-nst
  gdrdchi(kk) = 0.0
  do iter=-nst,nst
   gdrdchi(kk) = gdrdchi(kk) + (1D0/dchi) * coeffs1(iter) * gr(kk+iter)
  enddo
 enddo

 do kk=1,nst
  gdrdchi(kk) = 0.0
  do iter=1,nstencil
   gdrdchi(kk) = gdrdchi(kk) + (1D0/dchi) * bdy_coeffs1(kk,iter) * gr(iter)
  enddo
 enddo

 do kk=nr-nst+1,nr
  gdrdchi(kk) = 0.0
  do iter=1,nstencil
   gdrdchi(kk) = gdrdchi(kk) - (1D0/dchi) * bdy_coeffs1(nr-kk+1,nstencil-iter+1) * gr(nr-nstencil+iter)
  enddo
 enddo
 
 do k=1,mr
  drdchi(k) = gdrdchi(rrank*mr+k)
 enddo 

! Computing the metric for azimuthal directions
  metric_C = (1.0D0 + tanxi**2.0D0)**0.5D0
  metric_D = (1.0D0 + taneta**2.0D0)**0.5D0

  do j=1,my
   do i=1,mx
    metric_delta(i,j) = (1. + tanxi(i)**2. + taneta(j)**2.)
    sqrt_metric_delta(i,j) = (1. + tanxi(i)**2. + taneta(j)**2.)**(0.5)

    XY(i,j) = tanxi(i) * taneta(j)
    CD(i,j) = metric_C(i) * metric_D(j)
   enddo
  enddo

  ! INTERPOLATION POINTS:

  do i = 1,noverlap
   interpoints(:,i)  = atan(tan(geta)/tan(gxi(nx) + i*dxi))
  enddo

  ! DETERMINING SPHERICAL COORDINATES

  !REGION 1

  if(zoneid.eq.1)then
   do j=1,my
    do i=1,mx
     phi(i,j) = xi(i) !+ pi
     theta(i,j) = atan(1./(taneta(j) * cos(phi(i,j))))
     if (theta(i,j) < 0) theta(i,j) = theta(i,j) + pi
    enddo
   enddo 

  !REGION 2
  elseif(zoneid.eq.2)then
   do j=1,my
    do i=1,mx
     phi(i,j) = xi(i) + pi*0.5
     theta(i,j) = atan(1./(taneta(j) * sin(phi(i,j))))
     if (theta(i,j) < 0) theta(i,j) = theta(i,j) + pi
     if (phi(i,j) < 0) phi(i,j) = phi(i,j) + 2.*pi
    enddo
   enddo 

  !REGION 3
  elseif(zoneid.eq.3)then
   do j=1,my
    do i=1,mx
     phi(i,j) = xi(i) + pi
     theta(i,j) = atan(-1./(taneta(j) * cos(phi(i,j))))
     if (theta(i,j) < 0) theta(i,j) = theta(i,j) + pi
     if (phi(i,j) < 0) phi(i,j) = phi(i,j) + 2.*pi
    enddo
   enddo 

  !REGION 4
  elseif(zoneid.eq.4)then
   do j=1,my
    do i=1,mx
     phi(i,j) = xi(i) + pi*1.5
     theta(i,j) = atan(-1./(taneta(j) * sin(phi(i,j))))
     if (theta(i,j) < 0) theta(i,j) = theta(i,j) + pi
     if (phi(i,j) < 0) phi(i,j) = phi(i,j) + 2.*pi
    enddo
   enddo 

  !REGION 5
  elseif(zoneid.eq.5)then
   do j=1,my
    do i=1,mx
     theta(i,j) = atan((tanxi(i)**2. + taneta(j)**2.)**0.5)
     phi(i,j) = atan2(tanxi(i),-taneta(j)) 
     if (phi(i,j) < 0) phi(i,j) = phi(i,j) + 2.*pi
    enddo
   enddo 

  !REGION 6
  elseif(zoneid.eq.6)then
   do j=1,my
    do i=1,mx
     theta(i,j) = pi - atan((tanxi(i)**2. + taneta(j)**2.)**0.5)
     phi(i,j) = atan2(tanxi(i),taneta(j)) 
     if (phi(i,j) < 0) phi(i,j) = phi(i,j) + 2.*pi
    enddo
   enddo 
  endif

END SUBROUTINE SETUP_COORDINATES_METRICS_ALL
!----------------------------------------------------------------------------------

END MODULE METRIC
