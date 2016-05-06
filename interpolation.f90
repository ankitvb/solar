MODULE INTERPOLATION

USE INIT

CONTAINS

!----------------------------------------------------------------------------------------

SUBROUTINE SETUP_LAGRANGE_INTERPOLANT

  implicit none

!  include 'header'

  integer i, j, k, ipoint, iter, lower, upper

!   g^{\rm interp}(j) = \sum_{k=-n_{ord}}^{n_{ord}} I(k) f(j+k)
!   I(k) is the interpolation function
!   f(k) is the set of points we interpolate from
!   g^{interp}(k) is the set of points we interpolate to

  interpolant = 0.0
  do k=1,noverlap ! NUMBER OF OVERLAP LAYERS

   interpolant(k,1,k) = 1.0
   interpolant(-k,nx,k) = 1.0

! INTERPOLATION POINTS ARE NOT UNIFORMLY SPACED. HOW TO ACCOUNT FOR ERRORS?

   do i= 2,nx-1  !INTERPOLATION POINT OF INTEREST

    lower = -nord
    upper = nord

! TO PREVENT INDICES FROM BECOMING NEGATIVE OR > nx
    if (i .lt. nord+1)  lower = -i + 1
    if (i .gt. nx-nord) upper = nx - i

    interpolant(lower:upper,i,k) = 1.0

!   LAGRANGE INTERPOLATION 
    do ipoint = lower,upper   !POLYNOMIAL COEFFICIENT OF point "i + ipoint"
     do iter= lower,upper

      if (iter .ne. ipoint) interpolant(ipoint, i,k) = (interpoints(i,k) - geta(i + iter)) * &
                interpolant(ipoint, i,k)/(geta(i + ipoint) - geta(i + iter))

     enddo
    enddo

   enddo
  enddo

END SUBROUTINE SETUP_LAGRANGE_INTERPOLANT

!---------------------------------------------------------------------------------------------

SUBROUTINE SETUP_CONNECTIVITY
 
 implicit none

! include 'header'

  ! CONNECTIVITY MATRIX

  ! eta - eta :: 
  ! eta - xi  :: EACH FOLLOW A SET OF TRANSFORMATION/ INTERPOLATION LAWS
  ! xi - xi   ::

  ! FOR EACH EDGE, ENCODE THE TYPE OF CONNECTIVITY
  ! KEEP THESE TRANSFORMATION LAWS FOR THE THREE
  ! AS PART OF THE INITIALIZATION, PRECOMPUTE THE INTERPOLATION POINTS 
  ! PRECOMPUTE METRIC TERMS 
  ! CUBED SPHERE DONE !

  ! WRITE SUBROUTINE THAT COMPUTE CURL, DIV, LAPLACIAN, AND OF COURSE
  ! INDIVIDUAL DERIVATIVES

  ! ALL SLICES HAVE THE SAME PDE. BOUNDARY TRANSFORMATION RULES ENSURE
  ! CONSISTENCY AT INTERFACES BETWEEN SLICES

  ! TYPES OF INTERFACES between the i-j faces: 1, 2
  ! + 1 is a xi coordinate on the +xi side of the face i
  ! - 1 is a xi coordinate on the -xi side of the face i

  ! + 2 is an eta coordinate on the +eta side of the face i
  ! - 2 is an eta coordinate on the -eta side of the face i

  ! Note for the xi-eta connections, if the product
  ! of the (i,j) and (j,i) elements is -2, then the
  ! +xi is connected to -eta (or vice versa)

  interfaces(1,1) = 0
  interfaces(1,2) = 1
  interfaces(1,3) = 0
  interfaces(1,4) = -1
  interfaces(1,5) = 2
  interfaces(1,6) = -2

  interfaces(2,1) = -1
  interfaces(2,2) = 0
  interfaces(2,3) = 1
  interfaces(2,4) = 0
  interfaces(2,5) = 2
  interfaces(2,6) = -2


  interfaces(3,1) = 0
  interfaces(3,2) = -1
  interfaces(3,3) = 0
  interfaces(3,4) = 1
  interfaces(3,5) = 2
  interfaces(3,6) = -2

  interfaces(4,1) = 1
  interfaces(4,2) = 0
  interfaces(4,3) = -1
  interfaces(4,4) = 0
  interfaces(4,5) = 2
  interfaces(4,6) = -2

  interfaces(5,1) = -2
  interfaces(5,2) = 1
  interfaces(5,3) = 2
  interfaces(5,4) = -1
  interfaces(5,5) = 0
  interfaces(5,6) = 0

  interfaces(6,1) = 2
  interfaces(6,2) = 1
  interfaces(6,3) = -2
  interfaces(6,4) = -1
  interfaces(6,5) = 0
  interfaces(6,6) = 0

END SUBROUTINE SETUP_CONNECTIVITY

!---------------------------------------------------------------------------------------

SUBROUTINE INTERPOLATE_ALONG_Y(input_func, func_int)

 implicit none

 include 'mpif.h'
 
 integer i, j, k, iter, lower, upper
 integer jj
 real*8 func_int(noverlap,my,mr)
 real*8 input_func(noverlap,my,mr)
  
 real*8 temp(noverlap,nord,mr,2), fp(noverlap,nord,mr), fn(noverlap,nord,mr)

 integer status_arr(MPI_STATUS_SIZE,4), ierr, req(4)
  
!  Putting first nord and last nord data planes into temp arrays for sending
 do k=1,mr
  do j=1,nord
   do i=1,noverlap
    temp(i,j,k,1) = input_func(i,j,k)
    temp(i,j,k,2) = input_func(i,my-nord+j,k)
   enddo
  enddo
 enddo 

! Transferring data from neighbouring processors
 if(yrank.gt.0)&
 call MPI_IRECV(fp, noverlap*nord*mr, MPI_REAL8,source_right_y,&
                source_right_y, MPI_XYR_COMM, req(1), ierr)               ! Recv last nord planes
 if(yrank.lt.py-1)&
 call MPI_IRECV(fn, noverlap*nord*mr, MPI_REAL8,source_left_y,&
                source_left_y, MPI_XYR_COMM, req(2), ierr)                ! Recv first nord planes
 if(yrank.lt.py-1)&
 call MPI_ISEND(temp(1,1,1,2), noverlap*nord*mr, MPI_REAL8,dest_right_y,&
                zonerank, MPI_XYR_COMM, req(3), ierr)                     ! Send last nord planes
 if(yrank.gt.0)&
 call MPI_ISEND(temp(1,1,1,1), noverlap*nord*mr, MPI_REAL8,dest_left_y,&
                zonerank, MPI_XYR_COMM, req(4), ierr)                     ! Send first nord planes

! Compute interior in the meantime
 do k=1,mr
  do i=1,noverlap

   if(yrank.eq.0) func_int(i,1,k) = input_func(i,i+1,k)
   if(yrank.eq.py-1) func_int(i,my,k) = input_func(i,my-i,k)

   do j=nord+1,my-nord
    jj = yrank*my+j
 
    lower = -nord
    upper = nord

    func_int(i,j,k) = 0.0
    do iter=lower,upper
     func_int(i,j,k) = func_int(i,j,k) + input_func(i,j+iter,k) * interpolant(iter,jj,i)
    enddo

   enddo
  enddo
 enddo

! Wait to recieve planes
 if(yrank.eq.0)then
  call MPI_WAITALL(2, req(2), status_arr(1,2), ierr)
 elseif(yrank.eq.py-1)then
  call MPI_WAIT(req(1), status_arr(1,1), ierr)
  call MPI_WAIT(req(4), status_arr(1,4), ierr)
 else
  call MPI_WAITALL(4, req, status_arr, ierr)
 endif
 
 lower=1
 upper=my
 if(yrank.eq.0)    lower=2
 if(yrank.eq.py-1) upper=my-1

! Compute at intra-zone boundaries   
! Left (lower) boundary
 do k=1,mr
  do i=1,noverlap
   do j=lower,nord
     jj = yrank*my+j
     func_int(i,j,k) = 0.0
 
     if(yrank.gt.0)then ! Smaller stencil from jj=1 to jj=nord-1
      do iter=-nord,-j
       func_int(i,j,k) = func_int(i,j,k) + fp(i,nord+iter+j,k) * interpolant(iter,jj,i)
      enddo
     endif
     do iter=-j+1,nord
      func_int(i,j,k) = func_int(i,j,k) + input_func(i,j+iter,k) * interpolant(iter,jj,i)
     enddo 

   enddo
  enddo
 enddo

! Right (upper) boundary
 do k=1,mr
  do i=1,noverlap
   do j=my-nord+1,upper
     jj = yrank*my+j
     func_int(i,j,k) = 0.0
 
     do iter=-nord,my-j
      func_int(i,j,k) = func_int(i,j,k) + input_func(i,j+iter,k) * interpolant(iter,jj,i)
     enddo
     if(yrank.lt.py-1)then ! Smaller stencil from jj=ny-nord+1 to jj=ny
      do iter=my-j+1,nord
       func_int(i,j,k) = func_int(i,j,k) + fn(i,iter+j-my,k) * interpolant(iter,jj,i) 
      enddo
     endif
 
   enddo
  enddo
 enddo

END SUBROUTINE INTERPOLATE_ALONG_Y

!------------------------------------------------------------------------------------

SUBROUTINE INTERPOLATE_ALONG_X(input_func, func_int)

 implicit none

 include 'mpif.h'

 integer i, j, k, iter, lower, upper
 integer ii
 real*8 func_int(mx,noverlap,mr)
 real*8 input_func(mx,noverlap,mr)

 real*8 temp(nord,noverlap,mr,2), fp(nord,noverlap,mr), fn(nord,noverlap,mr)

 integer status_arr(MPI_STATUS_SIZE,4), ierr, req(4)

!  Putting first nord and last nord data planes into temp arrays for sending
 do k=1,mr
  do j=1,noverlap
   do i=1,nord
    temp(i,j,k,1) = input_func(i,j,k)
    temp(i,j,k,2) = input_func(mx-nord+i,j,k)
   enddo
  enddo
 enddo

! Transferring data from neighbouring processors
 if(xrank.gt.0)&
 call MPI_IRECV(fp, noverlap*nord*mr, MPI_REAL8,source_right_x,&
                source_right_x, MPI_XYR_COMM, req(1), ierr)               ! Recv last nord planes 
 if(xrank.lt.px-1)&
 call MPI_IRECV(fn, noverlap*nord*mr, MPI_REAL8,source_left_x,&
                source_left_x, MPI_XYR_COMM, req(2), ierr)                ! Recv first nord planes
 if(xrank.lt.px-1)&
 call MPI_ISEND(temp(1,1,1,2), noverlap*nord*mr, MPI_REAL8,dest_right_x,&
                zonerank, MPI_XYR_COMM, req(3), ierr)                     ! Send last nord planes 
 if(xrank.gt.0)&
 call MPI_ISEND(temp(1,1,1,1),  noverlap*nord*mr, MPI_REAL8,dest_left_x,&
                zonerank, MPI_XYR_COMM, req(4), ierr)                     ! Send first nord planes

! Compute interior in the meantime
 do k=1,mr
  do j=1,noverlap

   if(xrank.eq.0) func_int(1,j,k) = input_func(j+1,j,k)
   if(xrank.eq.px-1) func_int(mx,j,k) = input_func(mx-j,j,k)

   do i=nord+1,mx-nord
    ii = xrank*mx+i

    lower = -nord
    upper = nord

    func_int(i,j,k) = 0.0
    do iter=lower,upper
     func_int(i,j,k) = func_int(i,j,k) + input_func(i+iter,j,k) * interpolant(iter,ii,j)
    enddo

   enddo
  enddo
 enddo

! Wait to recieve planes 
 if(xrank.eq.0)then
  call MPI_WAITALL(2, req(2), status_arr(1,2), ierr)
 elseif(xrank.eq.px-1)then
  call MPI_WAIT(req(1), status_arr(1,1), ierr)
  call MPI_WAIT(req(4), status_arr(1,4), ierr)
 else
  call MPI_WAITALL(4, req, status_arr, ierr)
 endif

 lower=1
 upper=mx
 if(xrank.eq.0)    lower=2
 if(xrank.eq.px-1) upper=mx-1

! Compute at intra-zone boundaries
! Left boundary
 do k=1,mr
  do j=1,noverlap
   do i=lower,nord
    ii = xrank*mx+i
    func_int(i,j,k) = 0.0
 
    if(xrank.gt.0)then ! Smaller stencil from ii=1 to ii=nord-1
     do iter=-nord,-i
      func_int(i,j,k) = func_int(i,j,k) + fp(nord+i+iter,j,k) * interpolant(iter,ii,j)
     enddo
    endif
    do iter=-i+1,nord
     func_int(i,j,k) = func_int(i,j,k) + input_func(iter+i,j,k) * interpolant(iter,ii,j)
    enddo

   enddo
  enddo 
 enddo

! Right boundary
 do k=1,mr
  do j=1,noverlap
   do i=mx-nord+1,upper
    ii=xrank*mx+i
    func_int(i,j,k) = 0.0
 
    do iter=-nord,mx-i
     func_int(i,j,k) = func_int(i,j,k) + input_func(i+iter,j,k) * interpolant(iter,ii,j)
    enddo
    if(xrank.lt.px-1)then ! Smaller stencil from ii=nx-nord+1 to ii=nx
     do iter=mx-i+1,nord
      func_int(i,j,k) = func_int(i,j,k) + fn(iter+i-mx,j,k) * interpolant(iter,ii,j)
     enddo
    endif

   enddo
  enddo
 enddo

END SUBROUTINE INTERPOLATE_ALONG_X

!---------------------------------------------------------------------------------------
SUBROUTINE LAGRANGE_INTERP

  USE INIT

  implicit none

  !integer step1
  real*8  xt
  real*8, dimension(mx,my) :: forcing

  COMMON /TIME_FORCING/ forcing

  real*8, dimension(mx,my,nt) :: S_t

  COMMON /FORCING/ S_t
 
  xt = DBLE(step)*deltat/cadforcing + 1.0d0

  if (step_old .LT. floor(xt)) then

   x0 = (xt-4.0)*cadforcing
   x1 = (xt-3.0)*cadforcing
   x2 = (xt-2.0)*cadforcing
   x3 = (xt-1.0)*cadforcing
   x4 = (xt)*cadforcing
   x5 = (xt+1.0)*cadforcing
   x6 = (xt+2.0)*cadforcing

   LC0 = 0.0
   LC1 = 0.0
   LC2 = 0.0
   LC3 = 0.0
   LC4 = 0.0
   LC5 = 0.0
   LC6 = 0.0

  if (x0 .GE. 0) &
    LC0 = S_t(:,:,(floor(xt)-3))/((x0-x1) * (x0-x2) * (x0-x3) * (x0-x4) * (x0-x5) *(x0-x6))

   if (x1 .GE. 0) &
    LC1 = S_t(:,:,(floor(xt)-2))/((x1-x0) * (x1-x2) * (x1-x3) * (x1-x4) * (x1-x5) *(x1-x6))

   if (x2 .GE. 0) &
    LC2 = S_t(:,:,(floor(xt)-1))/((x2-x0) * (x2-x1) * (x2-x3) * (x2-x4) * (x2-x5) *(x2-x6))

   if (x3 .GE. 0) &
    LC3 = S_t(:,:,(floor(xt)))/((x3-x0) * (x3-x1) * (x3-x2) * (x3-x4) * (x3-x5) *(x3-x6))

    LC4 = S_t(:,:,(floor(xt)+1))/((x4-x0) * (x4-x1) * (x4-x2) * (x4-x3) * (x4-x5) *(x4-x6))

   if (floor(xt) .LT. nstep) &
    LC5 = S_t(:,:,(floor(xt)+2))/((x5-x0) * (x5-x1) * (x5-x2) * (x5-x3) * (x5-x4) *(x5-x6))

   if (floor(xt) .LT. (nstep-1)) &
    LC6 = S_t(:,:,(floor(xt)+3))/((x6-x0) * (x6-x1) * (x6-x2) * (x6-x3) * (x6-x4) *(x6-x5))

   step_old = floor(xt) 
  endif

  xt = (DBLE(step)-1.0)*deltat

  forcing = LC0 * (xt-x1) * (xt-x2) * (xt-x3) * (xt-x4) * (xt-x5) *(xt-x6)+ &
            LC1 * (xt-x0) * (xt-x2) * (xt-x3) * (xt-x4) * (xt-x5) *(xt-x6)+ &
            LC2 * (xt-x0) * (xt-x1) * (xt-x3) * (xt-x4) * (xt-x5) *(xt-x6)+ &
            LC3 * (xt-x0) * (xt-x1) * (xt-x2) * (xt-x4) * (xt-x5) *(xt-x6)+ &
            LC4 * (xt-x0) * (xt-x1) * (xt-x2) * (xt-x3) * (xt-x5) *(xt-x6)+ &
            LC5 * (xt-x0) * (xt-x1) * (xt-x2) * (xt-x3) * (xt-x4) *(xt-x6)+ &
            LC6 * (xt-x0) * (xt-x1) * (xt-x2) * (xt-x3) * (xt-x4) *(xt-x5)


END SUBROUTINE LAGRANGE_INTERP
!================================================================================================ 

END MODULE INTERPOLATION
