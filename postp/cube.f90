PROGRAM CUBE

  USE INIT
  USE INTERPOLATION
  USE METRIC
  USE TRANSFORM
  USE OUTPUT
  USE BSPLINE
 
  implicit none

  integer i, j, k,mm, facenum,nx2,ny2
  real*8 xknot(nx+kord+2*noverlap),yknot(ny+kord+2*noverlap),thetaknot(ntheta+kord),phiknot(nphi+kord),temphi
  real*8 xietaface(nphi,ntheta,0:2),x,y, temtheta,phiknot2(21+kord),maxtheta,xi_ext(ny+2*noverlap),eta_ext(ny+2*noverlap)
! Globals
  real*8, allocatable, dimension(:,:,:) :: bb, a2, temp,b,btemp
  real*8, allocatable, dimension(:,:) :: bcoef, bcoef2pi,bcoefs

  real*8, dimension(:,:,:,:), allocatable :: a
  real*8 thta(ntheta), phiglobal(nphi), phi2(21)

  real*8, dimension(noverlap,ny,6) :: interpolated_east, interpolated_west, &
		interpolated_north_transposed, interpolated_south_transposed
  real*8, dimension(nx,noverlap,6) :: interpolated_north, interpolated_south, &
		interpolated_west_transposed, interpolated_east_transposed

  nprocs = 24
  nvis = 3 
  o_rad = 343

  call SETUP_COORDINATES_METRICS_ALL

  call SETUP_LAGRANGE_INTERPOLANT

  call SETUP_CONNECTIVITY

  call CONVERT_VIS_TO_FITS
  

  stop

  CALL GLATS(ntheta, thta)

  thta = pi/2. - thta
  do k=1,nphi
   phiglobal(k) = (k-1.0)/(nphi*1.0) * 2.0 * pi + minval(phi(:,:,1))
  enddo

  if (xieta_to_thetaphi) then

  allocate(bcoefs((nx+2*noverlap)*(ny+2*noverlap),6),a(nx,ny,6,nt), a2(nphi,ntheta,nt),b(nx+2*noverlap,ny+2*noverlap,6) )

  do j=1,ntheta
   do i=1,nphi
    do k=1,6
    
    if (k==1) then
!     x = tan(phiglobal(i))
     y = 1./(tan(thta(j))*cos(phiglobal(i)))

     if ((abs(phiglobal(i)) .le. 0.25*pi) .and. (abs(y) .le. 1)) then
       xietaface(i,j,0) = 1
       xietaface(i,j,1) = phiglobal(i)
       xietaface(i,j,2) = atan(y)
     endif

    elseif (k==2) then

!     x = -1./tan(phiglobal(i))
     y = 1./(tan(thta(j))*sin(phiglobal(i)))

     if ((abs(phiglobal(i)-pi*0.5) .le. 0.25*pi) .and. (abs(y) .le. 1)) then
       xietaface(i,j,0) = 2
       xietaface(i,j,1) = phiglobal(i) - pi*0.5
       xietaface(i,j,2) = atan(y)
     endif

    elseif (k==3) then

     !x = tan(phiglobal(i))
     y = -1./(tan(thta(j))*cos(phiglobal(i)))
 
     if ((abs(phiglobal(i)-pi) .le. 0.25*pi) .and. (abs(y) .le. 1)) then
       xietaface(i,j,0) = 3
       xietaface(i,j,1) = phiglobal(i) - pi
       xietaface(i,j,2) = atan(y)
     endif

    elseif (k==4) then

!     x = 1./tan(phiglobal(i))
     y = -1./(tan(thta(j))*sin(phiglobal(i)))

     if ((abs(phiglobal(i)-1.5*pi) .le. 0.25*pi) .and. (abs(y) .le. 1)) then
       xietaface(i,j,0) = 4 
       xietaface(i,j,1) = phiglobal(i) - pi*1.5
       xietaface(i,j,2) = atan(y)
     endif

    elseif (k==5) then

     x = tan(thta(j))*sin(phiglobal(i))
     y = -tan(thta(j))*cos(phiglobal(i))

     if ((abs(x) .le.1) .and. (abs(y) .le. 1) .and. (thta(j) .le. pi*0.35)) then
       xietaface(i,j,0) = 5 
       xietaface(i,j,1) = atan(x)
       xietaface(i,j,2) = atan(y)
     endif

    elseif (k==6) then
     x = -tan(thta(j))*sin(phiglobal(i))
     y = -tan(thta(j))*cos(phiglobal(i))

      if ((abs(x) .le.1) .and. (abs(y) .le. 1) .and. (thta(j) .ge. pi*0.65)) then
       xietaface(i,j,0) = 6
       xietaface(i,j,1) = atan(x)
       xietaface(i,j,2) = atan(y)
     endif
  
   endif
   enddo
  enddo
 enddo

 eta_ext(noverlap+1:nx+noverlap) = eta
 xi_ext(noverlap+1:nx+noverlap) = eta

 do i=1,noverlap
  eta_ext(i) = -pi*0.25 - dxi*(noverlap-i+1) 
  xi_ext(i) = -pi*0.25 - dxi*(noverlap-i+1) 
  eta_ext(i+nx+noverlap) = pi*0.25 + dxi*i 
  xi_ext(i+nx+noverlap) = pi*0.25 + dxi*i 
 enddo

  nx2 = nx + 2*noverlap
  ny2 = nx2 
  ! GOING FROM XI-ETA ---> XI-ETA
  call dbsnak(nx2,xi_ext,kord,xknot)
  call dbsnak(nx2,eta_ext,kord,yknot)

!  call readfits('cubesph.fits',a(:,:,:,1),nx,ny,6)
  
!  allocate(btemp(nphi,ntheta,1))
!  call readfits('sph.fits',btemp,nphi,ntheta,1)

  call readfits('Visual/q_4_165.fits',a(:,:,:,1),nx,ny,6)

!  do k=1,6
!   do j=1,ny
!    do i=1,nx
!     a(i,j,k,1) = sin(theta(i,j,k))**2. * cos(phi(i,j,k))
!    enddo
!   enddo
!  enddo

  do mm=1,nt

   call interpolate_all(a(:,:,:,mm),b)

   do k=1,6
    call dbs2in(nx2,xi_ext,ny2,eta_ext,b(:,:,k),nx2,kord,kord,    &
          xknot,yknot,bcoefs(:,k))
   enddo

   do j=ntheta,1,-1
    do i=1,nphi
     facenum = int(xietaface(i,j,0) )
     a2(i,j,mm) = dbs2vl(xietaface(i,j,1),xietaface(i,j,2),kord,kord,xknot,yknot,nx2,ny2,bcoefs(:,facenum))
     !if (abs(a2(i,j,mm)-btemp(i,j,mm)) > 0.01 .and. facenum==3) print *,abs(a2(i,j,mm)-btemp(i,j,mm)),xietaface(i,j,1),xietaface(i,j,2)
    enddo
   enddo

  enddo

 call writefits('outputsph.fits',a2,nphi,ntheta,nt)
 print *,'here'
!   do j=1,ntheta
!    do i=1,nphi
!     a2(i,j,1) = sin(thta(j))**2. * cos(phiglobal(i))
!    enddo
!   enddo

! call writefits('exactsph.fits',a2,nphi,ntheta,nt)

 endif

 if (.not. xieta_to_thetaphi) then

  allocate(bb(nphi,ntheta,nt),a(nx,ny,6,nt),bcoef(nphi,ntheta),bcoef2pi(21,ntheta))

  call readfits('/data/seismo1/hanasoge/forcing171.fits',bb,nphi,ntheta,nt)


!   do j=1,ntheta
!    do i=1,nphi
!     bb(i,j,1) = sin(thta(j))**2. * cos(phiglobal(i))
!    enddo
!   enddo

  call dbsnak(ntheta,thta,kord,thetaknot)
  call dbsnak(nphi,phiglobal,kord,phiknot)
 
  allocate(temp(21,ntheta,nt))
  do k=1,nt
   do j=1,ntheta
    do i=nphi-10,nphi
     temp(i-nphi+11,j,k) = bb(i,j,k)
     phi2(i-nphi+11) = phiglobal(i)
    enddo
    do i=1,10
     temp(i+11,j,k) = bb(i,j,k)
     phi2(i+11) = phiglobal(i)+2.*pi
    enddo

   enddo
  enddo
  
  call dbsnak(21,phi2,kord,phiknot2)
  maxtheta = maxval(abs(pi*0.5-thta))

  do mm=1,nt
   print *,mm
   call dbs2in(nphi,phiglobal,ntheta,thta,bb(:,:,mm),nphi,kord,kord,&
         phiknot,thetaknot,bcoef)

   call dbs2in(21,phi2,ntheta,thta,temp(:,:,mm),21,kord,kord,&
         phiknot2,thetaknot,bcoef2pi)

   do k=6,1,-1
    do j=1,ny
     do i=1,nx
      temphi = phi(i,j,k)

       if (maxtheta .ge. abs(pi*0.5- theta(i,j,k))) then

         if ((temphi < phi2(4)) .and. (temphi > (phi2(15)-2.*pi))) then
  	  a(i,j,k,mm) = dbs2vl(temphi,theta(i,j,k),kord,kord,phiknot,thetaknot,nphi,ntheta,bcoef)
         else
          if (temphi .lt. 0) temphi = temphi + 2.*pi
          a(i,j,k,mm) = dbs2vl(temphi,theta(i,j,k),kord,kord,phiknot2,thetaknot,21,ntheta,bcoef2pi)
         endif
       endif

     enddo
    enddo
   enddo
  enddo  
  !call writefits4d('/data/seismo1/hanasoge/forcing_cubesph.fits',a(1,1,1,1),nx,ny,6,nt)

  call write_random_field('/data/seismo1/hanasoge/forcing_cubesph',a,nx,ny,nt)


!  do k=1,6
!   do j=1,ny
!    do i=1,nx
!     a(i,j,k,1) = sin(theta(i,j,k))**2. * cos(phi(i,j,k))
!    enddo
!   enddo
!  enddo
!   call writefits('exactcubesph.fits',a(1,1,1,1),nx,ny,6)
 
 endif


  
  
END PROGRAM CUBE

!--------------------------

SUBROUTINE GLATS(NLAT, THTA)

      implicit none

      INTEGER NLAT
      REAL*8 EPS
      PARAMETER(EPS = 10.0**(-15.))
      REAL THTA(NLAT)
      REAL WTS(NLAT)
      REAL*8 FI, FI1, FN, DOT, DN, DN1, A, B, G, GM, GP, GT,&
          FTEMP, GTEMP,pi
      INTEGER IR, NL, IRP, IRM, NITER

      pi = acos(-1.0d0)

      IR  = NLAT                                                            
      FI  = REAL(IR)                                                           
      FI1 = FI+1.                                                               
      FN  = REAL(NLAT/2)                                                       

      DO 20 NL = 1, NLAT/2                                                      
         DOT   = REAL(NL-1)                                                    
         THTA(NL) = -PI*.5*(DOT+.5)/FN + PI*.5                                     
         THTA(NL) =  SIN(THTA(NL))                                                    
   20 CONTINUE                                                                  

      DN  = FI/SQRT(4.*FI*FI-1.)                                                
      DN1 = FI1/SQRT(4.*FI1*FI1-1.)                                             
      A   = DN1*FI                                                              
      B   = DN*FI1                                                              
      IRP = IR + 1                                                              
      IRM = IR - 1                                                              

      DO 50 NL = 1, NLAT/2                                                      
         NITER = 0                                                              
   30    NITER = NITER + 1                                                      

         IF(NITER .GE. 100) THEN                                                
            WRITE (0, 90) NL                                                    
   90       FORMAT(/,' STSWM: FATAL ERROR IN GLATS:',/,&
             ' NO CONVERGENCE IN NEWTON ITERATION FOR ROOT',/,& 
                  ' LATITUDE  NL = ', I3)                                  
            STOP                                                                

         ENDIF                                                                  

         CALL ORDLEG(THTA(NL), IR, G)                                            
         CALL ORDLEG(THTA(NL), IRM, GM)                                            
         CALL ORDLEG(THTA(NL), IRP, GP)                                            

         GT    = (A*GP-B*GM)/(THTA(NL)*THTA(NL)-1.)                                   
         FTEMP = THTA(NL) - G/GT                                                   
         GTEMP = THTA(NL) - FTEMP                                                  
         THTA(NL) = FTEMP                                                          

         IF (ABS(GTEMP) .GT. 2.0*EPS) GO TO 30  
   50 CONTINUE                                                                  

      DO 60 NL = 1, NLAT/2                                                      
         A      = 2.0*(1.0 - THTA(NL)**2)                                          
         CALL ORDLEG(THTA(NL), IRM, B)                                             
         B      = B*B*FI*FI                                                     
         WTS(NL) = A*(FI - 0.5)/B                                                
         THTA (NL) = PI/2.0-ACOS(THTA(NL)) 
   60 CONTINUE                                                                  

      DO 70 NL = 1,NLAT/2
         THTA (NLAT-NL+1) = - THTA(NL)
         WTS(NLAT-NL+1) = WTS(NL)
   70 CONTINUE

      RETURN                                                                    
      END SUBROUTINE GLATS                                                                       

!===============================================================================

SUBROUTINE ORDLEG(COA, IR, SX)

      implicit none

      REAL*8 COA
      INTEGER IR
      REAL*8 SX
      REAL*8 SQR2,DELTA,THETA,C1,FN,FN2,FN2SQ,ANG,S1,C4,A,B,FK
      INTEGER IRPP,IRPPM,N,N1,KK,K

      SQR2  = SQRT(2.)                                                          
      IRPP  = IR + 1                                                            
      IRPPM = IRPP - 1                                                          
      DELTA = ACOS(COA)                                                         
      THETA = DELTA                                                             
      C1    = SQR2                                                              

      DO 20 N = 1, IRPPM                                                        
         FN    = REAL(N)                                                       
         FN2   = 2.*FN                                                          
         FN2SQ = FN2*FN2                                                        
         C1    = C1* SQRT(1.0-1.0/FN2SQ)                                        
   20 CONTINUE                                                                  

      N   = IRPPM                                                               
      ANG = FN*THETA                                                            
      S1  = 0.0                                                                 
      C4  = 1.0                                                                 
      A   = -1.0                                                                
      B   = 0.0                                                                 
      N1  = N+1                                                                 

      DO 30 KK = 1, N1, 2                                                       
         K = KK-1                                                               
         IF (K .EQ. N) C4 = 0.5*C4   
         S1  = S1+C4* COS(ANG)                                                  
         A   = A+2.0                                                            
         B   = B+1.0                                                            
         FK  = REAL(K)                                                         
         ANG = THETA*(FN-FK-2.0)                                                
         C4  = (A*(FN-B+1.0)/(B*(FN2-A)))*C4                                    
   30 CONTINUE                                                                  

      SX = S1*C1                                                                
      RETURN                                                                    

      END SUBROUTINE ORDLEG                                                                       

!--------------------------




