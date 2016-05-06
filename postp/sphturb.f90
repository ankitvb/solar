PROGRAM SPHTURB

 use bspline
 implicit none

 integer i, j, k, nx, ny, nz, kord
 parameter(nx = 256, ny = nx, nz = nx, kord = 6)

 integer nphi, ntheta
 parameter(nphi = 128, ntheta = nphi/2)

 real*8  t1, t2, thta(ntheta), phi(nphi)

 real*8  radial_location, pi
 parameter(radial_location = 0.5)

 real*8, dimension(nphi, ntheta, 3) :: xnew
 real*8, dimension(nx) :: x
 real*8, dimension(ny) :: y 
 real*8, dimension(nz) :: z

 real*8, dimension(nx+kord) :: xknot, yknot, zknot 
 real*8, dimension(nx, nx, nz) :: bb, bcoef
 real*8, dimension(nphi,ntheta) :: a2

 pi = ACOS(-1.0D0)
 CALL GLATS(ntheta, thta)

 do i=1,nx
  x(i) = (i-1.0)/(nx*1.0) - 0.5
 enddo

 y = x
 z = x

 thta = pi/2. - thta
 do i=1,nphi
  phi(i) = (i-1.0)/(nphi*1.0) * 2.*pi
 enddo

 print *,thta
stop

 do j=1,ntheta
  do i=1,nphi
   xnew(i,j,1) = radial_location * cos(thta(j)) * cos(phi(i))
   xnew(i,j,2) = radial_location * cos(thta(j)) * sin(phi(i))
   xnew(i,j,3) = radial_location * sin(thta(j)) 
  enddo
 enddo

 call readfits('/tmp20/shravan/spag4/a2.fits',bb,nx,ny,nz)

 call dbsnak(nx,x,kord,xknot)
 yknot = xknot
 zknot = zknot

! call dbs3in(nx,x,ny,y,nz,z,a,nx,ny,kord,kord,kord,    &
!         & xknot,yknot,zknot,bcoef)


 call CPU_TIME(t1)

 call dbs3in(nx,x,ny,y,nz,z,bb,nx,ny,kord,kord,kord,&
          xknot,yknot,zknot,bcoef)

  do j=1,ntheta
   do i=1,nphi
    a2(i,j) = dbs3vl(xnew(i,j,1),xnew(i,j,2),xnew(i,j,3),kord,kord,kord,&
		xknot,yknot,zknot,nx,ny,nz,bcoef)
   enddo
  enddo

 call writefits('/tmp20/shravan/spag4/aspline.fits',a2,nphi,ntheta,1)

END PROGRAM SPHTURB

!================================================================================

   SUBROUTINE readfits(filename,readarr,nx,ny,dim3)

      implicit none
      integer status,unit,readwrite,blocksize,naxes(3), tag, k, dim3 
      integer group,firstpix
      integer disp, temptyp, sendtyp
      integer nelements,nx,ny
      real*8 nullval,readarr(nx,ny,dim3)
      real*8, dimension(:,:,:), allocatable :: temp
      logical anynull
      character*(*) filename

      
       allocate(temp(nx,ny,dim3))
       status=0
       call ftgiou(unit,status)
       readwrite=0
       print *,'Now reading the file: '//filename
       call ftopen(unit,filename,readwrite,blocksize,status)

       naxes(1) = nx
       naxes(2) = ny
       naxes(3) = dim3
       nelements=naxes(1)*naxes(2)*naxes(3)
       group=1
       firstpix=1
       nullval=-999

       call ftgpvd(unit,group,firstpix,nelements,nullval, &
       &            temp,anynull,status)


       call ftclos(unit, status)
       call ftfiou(unit, status)

!       if (status .gt. 0)call printerror(status)

       readarr = temp

    end SUBROUTINE readfits

!================================================================================


   SUBROUTINE writefits(filename,dump_array, nx, ny, dim3)

      implicit none
      integer blocksize,bitpix,naxes(3),unit1,dim3
      integer status1,group,fpixel, ierr
      integer temptype, sendtyp,nx,ny
      integer*8 nelements
      real*8 dump_array(nx,ny,dim3)
      character*(*) filename
      logical simple,extend

      status1 = 0
      call ftgiou(unit1,status1)
      blocksize=1
      call ftinit(unit1,filename,blocksize,status1)
      simple=.true.
      bitpix=-64
      naxes(1)=nx
      naxes(2)=ny
      naxes(3)=dim3
      nelements=naxes(1)*naxes(2)*naxes(3)
      extend=.false.
      group=1
      fpixel=1
                                                                                                                                                      
      call ftphpr(unit1,simple,bitpix,3,naxes,0,1,extend,status1)
      call ftpprd(unit1,group,fpixel,nelements,dump_array,status1)
      call ftclos(unit1, status1)
      call ftfiou(unit1, status1)

!      if (status1 .gt. 0) call printerror(status1)

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

     end SUBROUTINE writefits
!===============================================================================

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

