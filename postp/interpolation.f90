MODULE INTERPOLATION

 USE INIT
 USE OUTPUT

CONTAINS

!----------------------------------------------------------------------------------------

SUBROUTINE SETUP_LAGRANGE_INTERPOLANT

  implicit none

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

      if (iter .ne. ipoint) interpolant(ipoint, i,k) = (interpoints(i,k) - eta(i + iter)) * &
                interpolant(ipoint, i,k)/(eta(i + ipoint) - eta(i + iter))

     enddo
    enddo

   enddo
  enddo

END SUBROUTINE SETUP_LAGRANGE_INTERPOLANT

!
!---------------------------------------------------------------------------------------------

SUBROUTINE SETUP_CONNECTIVITY
 implicit none

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

SUBROUTINE INTERPOLATE_ALONG_Y(input_func, func_int, interface_in, interface_out)

 implicit none

 integer i, j, k, lower, upper, interface_in, interface_out
 real*8 func_int(noverlap,ny)
 real*8 input_func(noverlap,ny)


 if (interface_in * interface_out .ne. -2) then

  do k=1, noverlap

   func_int(k,1) = input_func(k,k+1)
   func_int(k,nx) = input_func(k,nx-k)

   do j=2,nx-1
    lower = -nord
    upper = nord

    if (j .lt. nord+1) lower = -j + 1
    if (j .gt. nx-nord) upper = nx - j

    func_int(k,j) = 0.0
    do i=lower,upper
     func_int(k,j) = func_int(k,j) + input_func(k,j+i) * interpolant(i,j,k)
    enddo

   enddo
  enddo

 elseif (interface_in * interface_out .eq. -2) then

  do k=1, noverlap

   func_int(k,nx) = input_func(k,k+1)
   func_int(k,1) = input_func(k,nx-k)

   do j=2,nx-1
    lower = -nord
    upper = nord

    if (j .lt. nord+1) lower = -j + 1
    if (j .gt. nx-nord) upper = nx - j

    func_int(k,nx-j+1) = 0.0
    do i=lower,upper
     func_int(k,nx-j+1) = func_int(k,nx-j+1) + input_func(k,j+i) * interpolant(i,j,k)
    enddo

   enddo
  enddo

 endif

END SUBROUTINE INTERPOLATE_ALONG_Y

! ---------------------------------------------------------------------------------------------------------------------------

SUBROUTINE INTERPOLATE_ALL(a, aout)

  implicit none

  integer i,j,k
  real*8, dimension(nx,noverlap,6)  :: interpolated_south, interpolated_north
  real*8, dimension(noverlap,ny,6)  :: interpolated_south_transposed, interpolated_north_transposed
  real*8, dimension(noverlap, ny,6) :: interpolated_east, interpolated_west
  real*8, dimension(nx, noverlap,6) :: interpolated_east_transposed, interpolated_west_transposed
  real*8, dimension(nx,ny,6) :: a
  real*8, dimension(nx+2*noverlap,ny+2*noverlap,6) :: aout

  aout = 0.0
  aout(noverlap+1:nx+noverlap,noverlap+1:nx+noverlap,:) = a

  ! INTERPOLATION FROM  i --> j
  do j=1,6
   do i=1,6

    if (interfaces(i,j) .eq. +1) call INTERPOLATE_ALONG_Y(a(nx-1:nx-noverlap:-1,:,i), interpolated_east(:,:,i),&
                +1, interfaces(j,i))

    if (interfaces(i,j) .eq. -1) call INTERPOLATE_ALONG_Y(a(2:noverlap+1,:,i), interpolated_west(:,:,i),&
                -1, interfaces(j,i))

    if (interfaces(i,j) .eq. +2) call INTERPOLATE_ALONG_X(a(:,ny-1:ny-noverlap:-1,i), interpolated_north(:,:,i),&
                +2, interfaces(j,i))

    if (interfaces(i,j) .eq. -2) call INTERPOLATE_ALONG_X(a(:,2:noverlap+1,i), interpolated_south(:,:,i),&
                -2, interfaces(j,i))

   enddo
  enddo

! Transposing relevant arrays to make their shape conform to the one expected by derivative routines
  do j=1,noverlap
   do i=1,nx

     interpolated_east_transposed(i,j,5:6) = interpolated_east(j,i,5:6)
     interpolated_west_transposed(i,j,5:6) = interpolated_west(j,i,5:6) 
     interpolated_north_transposed(j,i,2)  = interpolated_north(i,j,2)
     interpolated_south_transposed(j,i,2)  = interpolated_south(i,j,2)
     interpolated_north_transposed(j,i,4)  = interpolated_north(i,j,4)
     interpolated_south_transposed(j,i,4)  = interpolated_south(i,j,4)

   enddo
  enddo
  

  !call DDX(a(:,:,1), deriv_xi(:,:,1), interpolated_east(:,:,4), interpolated_west(:,:,2))
  !call DDY(a(:,:,1), deriv_eta(:,:,1), interpolated_north(:,:,6), interpolated_south(:,:,5))

  aout(1:noverlap,noverlap+1:ny+noverlap,1) = interpolated_east(noverlap:1:-1,:,4)
  aout(nx+1+noverlap:nx+2*noverlap,noverlap+1:ny+noverlap,1) = interpolated_west(:,:,2)

  aout(noverlap+1:nx+noverlap,1:noverlap,1) = interpolated_north(:,noverlap:1:-1,6)
  aout(noverlap+1:nx+noverlap,ny+1+noverlap:ny+2*noverlap,1) =interpolated_south(:,:,5)


  !call DDX(a(:,:,2), deriv_xi(:,:,2), interpolated_east(:,:,1), interpolated_west(:,:,3))
  !call DDY(a(:,:,2), deriv_eta(:,:,2) ,interpolated_east_transposed(:,:,6), interpolated_east_transposed(:,:,5))

  aout(1:noverlap,noverlap+1:ny+noverlap,2) = interpolated_east(noverlap:1:-1,:,1)
  aout(nx+1+noverlap:nx+2*noverlap,noverlap+1:ny+noverlap,2) = interpolated_west(:,:,3)

  aout(noverlap+1:nx+noverlap,1:noverlap,2) = interpolated_east_transposed(:,noverlap:1:-1,6)
  aout(noverlap+1:nx+noverlap,ny+1+noverlap:ny+2*noverlap,2) = interpolated_east_transposed(:,:,5)

!  call DDX(a(:,:,3), deriv_xi(:,:,3), interpolated_east(:,:,2), interpolated_west(:,:,4))
!  call DDY(a(:,:,3), deriv_eta(:,:,3), interpolated_south(:,:,6), interpolated_north(:,:,5))

  aout(1:noverlap,noverlap+1:ny+noverlap,3) = interpolated_east(noverlap:1:-1,:,2)
  aout(nx+1+noverlap:nx+2*noverlap,noverlap+1:ny+noverlap,3) = interpolated_west(:,:,4)

  aout(noverlap+1:nx+noverlap,1:noverlap,3) = interpolated_south(:,noverlap:1:-1,6)
  aout(noverlap+1:nx+noverlap,ny+1+noverlap:ny+2*noverlap,3) = interpolated_north(:,:,5)

!  call DDX(a(:,:,4), deriv_xi(:,:,4), interpolated_east(:,:,3), interpolated_west(:,:,1))
!  call DDY(a(:,:,4), deriv_eta(:,:,4), interpolated_west_transposed(:,:,6), interpolated_west_transposed(:,:,5))

  aout(1:noverlap,noverlap+1:ny+noverlap,4) = interpolated_east(noverlap:1:-1,:,3)
  aout(nx+1+noverlap:nx+2*noverlap,noverlap+1:ny+noverlap,4) = interpolated_west(:,:,1)

  aout(noverlap+1:nx+noverlap,1:noverlap,4) = interpolated_west_transposed(:,noverlap:1:-1,6)
  aout(noverlap+1:nx+noverlap,ny+1+noverlap:ny+2*noverlap,4) = interpolated_west_transposed(:,:,5)

!  call DDX(a(:,:,5), deriv_xi(:,:,5), interpolated_north_transposed(:,:,4), interpolated_north_transposed(:,:,2))
!  call DDY(a(:,:,5), deriv_eta(:,:,5), interpolated_north(:,:,1), interpolated_north(:,:,3))

  aout(1:noverlap,noverlap+1:ny+noverlap,5) = interpolated_north_transposed(noverlap:1:-1,:,4)
  aout(nx+1+noverlap:nx+2*noverlap,noverlap+1:ny+noverlap,5) = interpolated_north_transposed(:,:,2)

  aout(noverlap+1:nx+noverlap,1:noverlap,5) = interpolated_north(:,noverlap:1:-1,1)
  aout(noverlap+1:nx+noverlap,ny+1+noverlap:ny+2*noverlap,5) = interpolated_north(:,:,3)

!  call DDX(a(:,:,6), deriv_xi(:,:,6), interpolated_south_transposed(:,:,4), interpolated_south_transposed(:,:,2))
!  call DDY(a(:,:,6), deriv_eta(:,:,6), interpolated_south(:,:,3), interpolated_south(:,:,1))

  aout(1:noverlap,noverlap+1:ny+noverlap,6) = interpolated_south_transposed(noverlap:1:-1,:,4)
  aout(nx+1+noverlap:nx+2*noverlap,noverlap+1:ny+noverlap,6) = interpolated_south_transposed(:,:,2)

  aout(noverlap+1:nx+noverlap,1:noverlap,6) = interpolated_south(:,noverlap:1:-1,3)
  aout(noverlap+1:nx+noverlap,ny+1+noverlap:ny+2*noverlap,6) = interpolated_south(:,:,1)

END SUBROUTINE INTERPOLATE_ALL

!------------------------------------------------------------------------------------

SUBROUTINE INTERPOLATE_ALONG_X(input_func, func_int, interface_in, interface_out)

 implicit none

 integer i, j, k, lower, upper, interface_in, interface_out
 real*8 func_int(nx,noverlap)
 real*8 input_func(nx, noverlap)

 if (interface_in * interface_out .ne. -2) then
  do k=1, noverlap

   func_int(1,k) = input_func(k+1,k)
   func_int(nx,k) = input_func(nx-k,k)

   do j=2,nx-1
    lower = -nord
    upper = nord

    if (j .lt. nord+1) lower = -j + 1
    if (j .gt. nx-nord) upper = nx - j

    func_int(j,k) = 0.0
    do i=lower,upper
     func_int(j,k) = func_int(j,k) + input_func(j+i,k) * interpolant(i,j,k)
    enddo

   enddo
  enddo

 elseif (interface_in * interface_out .eq. -2) then

  do k=1, noverlap

   func_int(nx,k) = input_func(k+1,k)
   func_int(1,k) = input_func(nx-k,k)

   do j=2,nx-1
    lower = -nord
    upper = nord

    if (j .lt. nord+1) lower = -j + 1
    if (j .gt. nx-nord) upper = nx - j

    func_int(nx-j+1,k) = 0.0
    do i=lower,upper
     func_int(nx-j+1,k) = func_int(nx-j+1,k) + input_func(j+i,k) * interpolant(i,j,k)
    enddo

   enddo
  enddo

 endif

END SUBROUTINE INTERPOLATE_ALONG_X

!---------------------------------------------------------------------------------------

END MODULE INTERPOLATION
