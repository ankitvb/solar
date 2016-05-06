MODULE DERIVATIVE

USE INIT

CONTAINS

!----------------------------------------------------------------------------------------
SUBROUTINE DDY(f, dy_f, interpol_low, interpol_high, interpol_left, interpol_right)

  implicit none

  include 'mpif.h'

  integer i, j, k, iter
  real*8 f(mx,my,mr), dy_f(mx,my,mr)
  real*8 interpol_low(mx,noverlap,mr), interpol_high(mx,noverlap,mr)
  real*8 interpol_left(mx,noverlap,mr), interpol_right(mx,noverlap,mr)
  real*8 fac

  real*8 temp(mx,nst,mr,2), fp(mx,nst,mr), fn(mx,nst,mr)  
  real*8 tmp1(mx,nst,mr), tmp2(mx,nst,mr)

  integer ierr
  integer status_arr1(MPI_STATUS_SIZE,4), status_arr2(MPI_STATUS_SIZE,6)
  integer req1(4), req2(6)

  fac = 1.0D0/deta  

! Putting first nst and last nst data planes into temp arrays for sending
 do k=1,mr
  do j=1,nst
   do i=1,mx
    temp(i,j,k,1) = f(i,j,k)                                                 
    temp(i,j,k,2) = f(i,my-j+1,k)
   enddo
  enddo
 enddo

! Intra-communication
! Transferring data from neighbouring processors
 if(yrank.gt.0)&
 call MPI_IRECV(fp, mx*nst*mr, MPI_REAL8, source_right_y,&
                source_right_y, MPI_XYR_COMM, req1(1), ierr)               ! Recv last nord planes
 if(yrank.lt.py-1)&
 call MPI_IRECV(fn, mx*nst*mr, MPI_REAL8, source_left_y,&
                source_left_y, MPI_XYR_COMM, req1(2), ierr)                ! Recv first nord planes
 if(yrank.lt.py-1)&
 call MPI_ISEND(temp(1,1,1,2), mx*nst*mr, MPI_REAL8, dest_right_y,&
                zonerank, MPI_XYR_COMM, req1(3), ierr)                     ! Send last nord planes
 if(yrank.gt.0)&
 call MPI_ISEND(temp(1,1,1,1), mx*nst*mr, MPI_REAL8, dest_left_y,&
                zonerank, MPI_XYR_COMM, req1(4), ierr)                     ! Send first nord planes

! Inter-communication
 if(zoneid.eq.1)then
  if(yrank.eq.0)&
  call MPI_IRECV(fp, mx*nst*mr, MPI_REAL8, icom_source_right_y,&
                 icom_source_right_y, MPI_INTERCOMM(1,6), req2(1), ierr)
  if(yrank.eq.py-1)&
  call MPI_IRECV(fn, mx*nst*mr, MPI_REAL8, icom_source_left_y,&
                 icom_source_left_y,  MPI_INTERCOMM(1,5), req2(2), ierr) 
  if(yrank.eq.0)&
  call MPI_ISEND(interpol_low, mx*nst*mr, MPI_REAL8, icom_dest_left_y,&
                 zonerank, MPI_INTERCOMM(1,6), req2(3), ierr)
  if(yrank.eq.py-1)&
  call MPI_ISEND(interpol_high, mx*nst*mr, MPI_REAL8, icom_dest_right_y,&
                 zonerank, MPI_INTERCOMM(1,5), req2(4), ierr)  
 elseif(zoneid.eq.2)then
  if(yrank.eq.0)&
  call MPI_IRECV(fp, mx*nst*mr, MPI_REAL8, icom_source_right_y,&
                 icom_source_right_y, MPI_INTERCOMM(2,6), req2(1), ierr)
  if(yrank.eq.py-1)&
  call MPI_IRECV(fn, mx*nst*mr, MPI_REAL8, icom_source_left_y,&
                 icom_source_left_y,  MPI_INTERCOMM(2,5), req2(2), ierr)
 elseif(zoneid.eq.3)then
  if(yrank.eq.0)&
  call MPI_IRECV(fp,  mx*nst*mr, MPI_REAL8, icom_source_right_y,&
                 icom_source_right_y, MPI_INTERCOMM(3,6), req2(1), ierr)
  if(yrank.eq.py-1)&
  call MPI_IRECV(fn, mx*nst*mr, MPI_REAL8, icom_source_left_y,&
                 icom_source_left_y,  MPI_INTERCOMM(3,5), req2(2), ierr)
  if(yrank.eq.0)&
  call MPI_ISEND(interpol_low, mx*nst*mr, MPI_REAL8, icom_dest_left_y,&
                 zonerank, MPI_INTERCOMM(3,6), req2(3), ierr)
  if(yrank.eq.py-1)&
  call MPI_ISEND(interpol_high, mx*nst*mr, MPI_REAL8, icom_dest_right_y,&
                 zonerank, MPI_INTERCOMM(3,5), req2(4), ierr)     
 elseif(zoneid.eq.4)then
  if(yrank.eq.0)&
  call MPI_IRECV(fp,  mx*nst*mr, MPI_REAL8, icom_source_right_y,&
                 icom_source_right_y, MPI_INTERCOMM(4,6), req2(1), ierr)
  if(yrank.eq.py-1)&
  call MPI_IRECV(fn, mx*nst*mr, MPI_REAL8, icom_source_left_y,&
                 icom_source_left_y,  MPI_INTERCOMM(4,5), req2(2), ierr)  
 elseif(zoneid.eq.5)then
  if(yrank.eq.0)&
  call MPI_IRECV(fp,  mx*nst*mr, MPI_REAL8, icom_source_right_y,&
                  icom_source_right_y, MPI_INTERCOMM(1,5), req2(1), ierr)
  if(yrank.eq.py-1)&
  call MPI_IRECV(fn,  mx*nst*mr, MPI_REAL8, icom_source_left_y,&
                  icom_source_left_y, MPI_INTERCOMM(3,5), req2(2), ierr)
  if(yrank.eq.0)&
  call MPI_ISEND(interpol_low, mx*nst*mr, MPI_REAL8, icom_dest_left_y,&
                  zonerank, MPI_INTERCOMM(1,5), req2(3), ierr)
  if(yrank.eq.py-1)&
  call MPI_ISEND(interpol_high, mx*nst*mr, MPI_REAL8, icom_dest_right_y,&
                  zonerank, MPI_INTERCOMM(3,5), req2(4), ierr)
  if(xrank.eq.px-1)&
  call MPI_ISEND(interpol_right, mx*nst*mr, MPI_REAL8, icom_dest_right_x,&
                  zonerank, MPI_INTERCOMM(2,5), req2(5), ierr)
  if(xrank.eq.0)&
  call MPI_ISEND(interpol_left, mx*nst*mr, MPI_REAL8, icom_dest_left_x,&
                  zonerank, MPI_INTERCOMM(4,5), req2(6), ierr)
 elseif(zoneid.eq.6)then
  if(yrank.eq.0)&
  call MPI_IRECV(fp, mx*nst*mr, MPI_REAL8, icom_source_right_y,&
                 icom_source_right_y, MPI_INTERCOMM(3,6), req2(1), ierr)
  if(yrank.eq.py-1)&
  call MPI_IRECV(fn, mx*nst*mr, MPI_REAL8, icom_source_left_y,&
                 icom_source_left_y, MPI_INTERCOMM(1,6), req2(2), ierr)
  if(yrank.eq.0)&
  call MPI_ISEND(interpol_low, mx*nst*mr, MPI_REAL8, icom_dest_left_y,&
                  zonerank, MPI_INTERCOMM(3,6), req2(3), ierr)
  if(yrank.eq.py-1)&
  call MPI_ISEND(interpol_high, mx*nst*mr, MPI_REAL8, icom_dest_right_y,&
                  zonerank, MPI_INTERCOMM(1,6), req2(4), ierr)
  if(xrank.eq.px-1)&
  call MPI_ISEND(interpol_right, mx*nst*mr, MPI_REAL8, icom_dest_right_x,&
                  zonerank, MPI_INTERCOMM(2,6), req2(5), ierr)
  if(xrank.eq.0)&
  call MPI_ISEND(interpol_left, mx*nst*mr, MPI_REAL8, icom_dest_left_x,&
                  zonerank, MPI_INTERCOMM(4,6), req2(6), ierr)
 endif

! Compute interior in the meantime                                                      
! ALL INTERIOR POINTS                                                                   
 do k=1,mr                                                                                       
  do j=nst+1,my-nst
   do i=1,mx
    dy_f(i,j,k) = 0.0                                                                     

    do iter=-nst,nst                                                                         
     dy_f(i,j,k) = dy_f(i,j,k) + fac * coeffs1(iter) * f(i,j+iter,k)
    enddo 

   enddo
  enddo
 enddo                                                                                       

! Waiting for intra communication to complete 
! Wait to recieve planes
 if(yrank.eq.0)then
  call MPI_WAITALL(2, req1(2), status_arr1(1,2), ierr)
 elseif(yrank.eq.py-1)then
  call MPI_WAIT(req1(1), status_arr1(1,1), ierr)
  call MPI_WAIT(req1(4), status_arr1(1,4), ierr)
 else
  call MPI_WAITALL(4, req1, status_arr1, ierr)
 endif

 if(zoneid.eq.1)then
  if(yrank.eq.0)&
   call MPI_WAIT(req2(1), status_arr2(1,1), ierr)
  if(yrank.eq.py-1)&
   call MPI_WAIT(req2(2), status_arr2(1,2), ierr)
  if(yrank.eq.0)&
   call MPI_WAIT(req2(3), status_arr2(1,3), ierr)
  if(yrank.eq.py-1)&
   call MPI_WAIT(req2(4), status_arr2(1,4), ierr)
 elseif(zoneid.eq.2)then
  if(yrank.eq.0)&
   call MPI_WAIT(req2(1), status_arr2(1,1), ierr)
  if(yrank.eq.py-1)&
   call MPI_WAIT(req2(2), status_arr2(1,2), ierr)
 elseif(zoneid.eq.3)then
  if(yrank.eq.0)&
   call MPI_WAIT(req2(1), status_arr2(1,1), ierr)
  if(yrank.eq.py-1)&
   call MPI_WAIT(req2(2), status_arr2(1,2), ierr)
  if(yrank.eq.0)&
   call MPI_WAIT(req2(3), status_arr2(1,3), ierr)
  if(yrank.eq.py-1)&
   call MPI_WAIT(req2(4), status_arr2(1,4), ierr)
 elseif(zoneid.eq.4)then
  if(yrank.eq.0)&
   call MPI_WAIT(req2(1), status_arr2(1,1), ierr)
  if(yrank.eq.py-1)&
   call MPI_WAIT(req2(2), status_arr2(1,2), ierr)
 elseif(zoneid.eq.5)then
  if(yrank.eq.0)&
   call MPI_WAIT(req2(1), status_arr2(1,1), ierr)
  if(yrank.eq.py-1)&
   call MPI_WAIT(req2(2), status_arr2(1,2), ierr)
  if(yrank.eq.0)&
   call MPI_WAIT(req2(3), status_arr2(1,3), ierr)
  if(yrank.eq.py-1)&
   call MPI_WAIT(req2(4), status_arr2(1,4), ierr)
  if(xrank.eq.px-1)&
   call MPI_WAIT(req2(5), status_arr2(1,5), ierr)
  if(xrank.eq.0)&
   call MPI_WAIT(req2(6), status_arr2(1,6), ierr)
 elseif(zoneid.eq.6)then
  if(yrank.eq.0)&
   call MPI_WAIT(req2(1), status_arr2(1,1), ierr)
  if(yrank.eq.py-1)&
   call MPI_WAIT(req2(2), status_arr2(1,2), ierr)
  if(yrank.eq.0)&
   call MPI_WAIT(req2(3), status_arr2(1,3), ierr)
  if(yrank.eq.py-1)&
   call MPI_WAIT(req2(4), status_arr2(1,4), ierr)
  if(xrank.eq.px-1)&
   call MPI_WAIT(req2(5), status_arr2(1,5), ierr)
  if(xrank.eq.0)&
   call MPI_WAIT(req2(6), status_arr2(1,6), ierr)
 endif

! Switching indices for reverse aligned arrays
 if(zoneid.eq.2)then
  if(yrank.eq.0)then
   tmp1=fp
   do k=1,mr
    do j=1,nst
     do i=1,mx
      fp(i,j,k) = tmp1(mx-i+1,j,k)
     enddo
    enddo
   enddo
  endif
 elseif(zoneid.eq.3)then
  tmp1=fp; tmp2=fn 
  if(yrank.eq.0)then
   do k=1,mr
    do j=1,nst
     do i=1,mx
      fp(i,j,k) = tmp1(mx-i+1,j,k)
     enddo
    enddo
   enddo
  endif
  if(yrank.eq.py-1)then
   do k=1,mr
    do j=1,nst
     do i=1,mx          
      fn(i,j,k) = tmp2(mx-i+1,j,k)
     enddo
    enddo 
   enddo
  endif
 elseif(zoneid.eq.4)then
  tmp2=fn
  if(yrank.eq.py-1)then
   do k=1,mr
    do j=1,nst
     do i=1,mx
      fn(i,j,k) = tmp2(mx-i+1,j,k)
     enddo
    enddo
   enddo
  endif
 elseif(zoneid.eq.5)then
  tmp2=fn
  if(yrank.eq.py-1)then
   do k=1,mr
    do j=1,nst
     do i=1,mx
      fn(i,j,k) = tmp2(mx-i+1,j,k)
     enddo
    enddo
   enddo
  endif
 elseif(zoneid.eq.6)then
  tmp1=fp
  if(yrank.eq.0)then
   do k=1,mr
    do j=1,nst
     do i=1,mx
      fp(i,j,k) = tmp1(mx-i+1,j,k)
     enddo
    enddo
   enddo
  endif 
 endif

! LOWER BOUNDARY POINTS

  do k=1,mr
   do j=1,nst
    do i=1,mx
     dy_f(i,j,k) = 0.0

     do iter=-nst,-j
      dy_f(i,j,k) = dy_f(i,j,k) + fac * coeffs1(iter) * fp(i,-iter-j+1,k)
     enddo

     do iter=-j+1,nst
      dy_f(i,j,k) = dy_f(i,j,k) + fac * coeffs1(iter) * f(i,iter+j,k)
     enddo

    enddo
   enddo
  enddo

! UPPER BOUNDARY POINTS

  do k=1,mr
   do j=my-nst+1,my
    do i=1,mx
     dy_f(i,j,k) = 0.0

     do iter=-nst,my-j
      dy_f(i,j,k) = dy_f(i,j,k) + fac * coeffs1(iter) * f(i,iter+j,k)
     enddo

     do iter=my-j+1,nst
      dy_f(i,j,k) = dy_f(i,j,k) + fac * coeffs1(iter) * fn(i,iter+j-my,k) 
     enddo

    enddo
   enddo
  enddo


END SUBROUTINE DDY

!-------------------------------------------------------------------------------

SUBROUTINE DDX(f, dx_f, interpol_low, interpol_high, interpol_left, interpol_right)

 implicit none

 include 'mpif.h'

 integer i, j, k, iter
 real*8 f(mx,my,mr), dx_f(mx,my,mr)
 real*8 interpol_left(noverlap,my,mr), interpol_right(noverlap,my,mr)
 real*8 interpol_low(noverlap,my,mr), interpol_high(noverlap,my,mr)
 real*8 fac

 real*8 temp(nst,my,mr,2), fp(nst,my,mr), fn(nst,my,mr)
 real*8 tmp1(nst,my,mr), tmp2(nst,my,mr)

 integer ierr
 integer status_arr1(MPI_STATUS_SIZE,4), status_arr2(MPI_STATUS_SIZE,6)
 integer req1(4), req2(6)

 fac = 1.0D0/dxi

! Putting first nst and last nst data planes into temp arrays for sending
 do k=1,mr
  do j=1,my
   do i=1,nst
    temp(i,j,k,1) = f(i,j,k)
    temp(i,j,k,2) = f(mx-i+1,j,k)
   enddo
  enddo
 enddo

! Intra-communication
! Transferring data from neighbouring processors
 if(xrank.gt.0)&
 call MPI_IRECV(fp, nst*my*mr, MPI_REAL8, source_right_x,&
                source_right_x, MPI_XYR_COMM, req1(1), ierr)           ! Recv last nst planes
 if(xrank.lt.px-1)&
 call MPI_IRECV(fn, nst*my*mr, MPI_REAL8, source_left_x,&
                source_left_x, MPI_XYR_COMM, req1(2), ierr)            ! Recv first nst planes
 if(xrank.lt.px-1)&
 call MPI_ISEND(temp(1,1,1,2), nst*my*mr, MPI_REAL8, dest_right_x,&
                zonerank, MPI_XYR_COMM, req1(3), ierr)                 ! Send last nst planes
 if(xrank.gt.0)&
 call MPI_ISEND(temp(1,1,1,1), nst*my*mr, MPI_REAL8, dest_left_x,&
                zonerank, MPI_XYR_COMM, req1(4), ierr)                 ! Send first nst planes 

! Inter-communication 
 if(zoneid.eq.1)then
  if(xrank.eq.0)&
  call MPI_IRECV(fp, nst*my*mr, MPI_REAL8, icom_source_right_x,&
                 icom_source_right_x, MPI_INTERCOMM(1,4), req2(1), ierr)
  if(xrank.eq.px-1)&
  call MPI_IRECV(fn, nst*my*mr, MPI_REAL8, icom_source_left_x,&
                 icom_source_left_x,  MPI_INTERCOMM(1,2), req2(2), ierr)
  if(xrank.eq.0)&
  call MPI_ISEND(interpol_left, nst*my*mr, MPI_REAL8, icom_dest_left_x,&
                 zonerank, MPI_INTERCOMM(1,4), req2(3), ierr)
  if(xrank.eq.px-1)&
  call MPI_ISEND(interpol_right, nst*my*mr, MPI_REAL8, icom_dest_right_x,&
                 zonerank, MPI_INTERCOMM(1,2), req2(4), ierr)
 elseif(zoneid.eq.2)then
  if(xrank.eq.0)&
  call MPI_IRECV(fp, nst*my*mr, MPI_REAL8, icom_source_right_x,&
                 icom_source_right_x, MPI_INTERCOMM(1,2), req2(1), ierr)
  if(xrank.eq.px-1)&
  call MPI_IRECV(fn, nst*my*mr, MPI_REAL8, icom_source_left_x,&
                 icom_source_left_x, MPI_INTERCOMM(2,3), req2(2), ierr)
  if(xrank.eq.0)&
  call MPI_ISEND(interpol_left, nst*my*mr, MPI_REAL8, icom_dest_left_x,&
                 zonerank, MPI_INTERCOMM(1,2), req2(3), ierr)
  if(xrank.eq.px-1)&
  call MPI_ISEND(interpol_right, nst*my*mr, MPI_REAL8, icom_dest_right_x,&
                 zonerank, MPI_INTERCOMM(2,3), req2(4), ierr)
  if(yrank.eq.0)&
  call MPI_ISEND(interpol_low, nst*my*mr, MPI_REAL8, icom_dest_left_y,&
                 zonerank, MPI_INTERCOMM(2,6), req2(5), ierr)
  if(yrank.eq.py-1)&
  call MPI_ISEND(interpol_high, nst*my*mr, MPI_REAL8, icom_dest_right_y,&
                 zonerank, MPI_INTERCOMM(2,5), req2(6), ierr)
 elseif(zoneid.eq.3)then
  if(xrank.eq.0)&
  call MPI_IRECV(fp, nst*my*mr, MPI_REAL8, icom_source_right_x,&
                 icom_source_right_x, MPI_INTERCOMM(2,3), req2(1), ierr)
  if(xrank.eq.px-1)&
  call MPI_IRECV(fn, nst*my*mr, MPI_REAL8, icom_source_left_x,&
                  icom_source_left_x,  MPI_INTERCOMM(3,4), req2(2), ierr)
  if(xrank.eq.0)&
  call MPI_ISEND(interpol_left, nst*my*mr, MPI_REAL8, icom_dest_left_x,&
                 zonerank, MPI_INTERCOMM(2,3), req2(3), ierr)
  if(xrank.eq.px-1)&
  call MPI_ISEND(interpol_right, nst*my*mr, MPI_REAL8, icom_dest_right_x,&
                 zonerank, MPI_INTERCOMM(3,4), req2(4), ierr)
 elseif(zoneid.eq.4)then
  if(xrank.eq.0)&
  call MPI_IRECV(fp, nst*my*mr, MPI_REAL8, icom_source_right_x,&
                 icom_source_right_x, MPI_INTERCOMM(3,4), req2(1), ierr)
  if(xrank.eq.px-1)&
  call MPI_IRECV(fn, nst*my*mr, MPI_REAL8, icom_source_left_x,&
                 icom_source_left_x, MPI_INTERCOMM(1,4), req2(2), ierr)
  if(xrank.eq.0)&
  call MPI_ISEND(interpol_left, nst*my*mr, MPI_REAL8, icom_dest_left_x,&
                 zonerank, MPI_INTERCOMM(3,4), req2(3), ierr)
  if(xrank.eq.px-1)&
  call MPI_ISEND(interpol_right, nst*my*mr, MPI_REAL8, icom_dest_right_x,&
                 zonerank, MPI_INTERCOMM(1,4), req2(4), ierr)
  if(yrank.eq.0)&
  call MPI_ISEND(interpol_low, nst*my*mr, MPI_REAL8, icom_dest_left_y,&
                 zonerank, MPI_INTERCOMM(4,6), req2(5), ierr)
  if(yrank.eq.py-1)&
  call MPI_ISEND(interpol_high, nst*my*mr, MPI_REAL8, icom_dest_right_y,&
                 zonerank, MPI_INTERCOMM(4,5), req2(6), ierr)
 elseif(zoneid.eq.5)then
  if(xrank.eq.0)&
  call MPI_IRECV(fp, nst*my*mr, MPI_REAL8, icom_source_right_x, &
                 icom_source_right_x, MPI_INTERCOMM(4,5), req2(1), ierr)
  if(xrank.eq.px-1)&
  call MPI_IRECV(fn, nst*my*mr, MPI_REAL8, icom_source_left_x, &
                 icom_source_left_x, MPI_INTERCOMM(2,5), req2(2), ierr)
 elseif(zoneid.eq.6)then
  if(xrank.eq.0)&
  call MPI_IRECV(fp, nst*my*mr, MPI_REAL8, icom_source_right_x,&
                 icom_source_right_x, MPI_INTERCOMM(4,6), req2(1), ierr)
  if(xrank.eq.px-1)&
  call MPI_IRECV(fn, nst*my*mr, MPI_REAL8, icom_source_left_x, &
                 icom_source_left_x, MPI_INTERCOMM(2,6), req2(2), ierr)
 endif

! Compute interior in the meantime
! ALL INTERIOR POINTS 
 do k=1,mr
  do j=1,my
   do i=nst+1,mx-nst
    dx_f(i,j,k) = 0.0

    do iter=-nst,nst
     dx_f(i,j,k) = dx_f(i,j,k) + fac * coeffs1(iter) * f(i+iter,j,k)
    enddo

   enddo
  enddo
 enddo

! Waiting for intra communication to complete
! Wait to recieve planes
 if(xrank.eq.0)then
  call MPI_WAITALL(2, req1(2), status_arr1(1,2), ierr)
 elseif(xrank.eq.px-1)then
  call MPI_WAIT(req1(1), status_arr1(1,1), ierr)
  call MPI_WAIT(req1(4), status_arr1(1,4), ierr)
 else
  call MPI_WAITALL(4, req1, status_arr1, ierr)
 endif

! Waiting for inter communication to complete
! Wait to recieve planes
 if(zoneid.eq.1)then
  if(xrank.eq.0)&
   call MPI_WAIT(req2(1), status_arr2(1,1), ierr)
  if(xrank.eq.px-1)&
   call MPI_WAIT(req2(2), status_arr2(1,2), ierr)
  if(xrank.eq.0)&
   call MPI_WAIT(req2(3), status_arr2(1,3), ierr)
  if(xrank.eq.px-1)&
   call MPI_WAIT(req2(4), status_arr2(1,4), ierr)
 elseif(zoneid.eq.2)then
  if(xrank.eq.0)&
   call MPI_WAIT(req2(1), status_arr2(1,1), ierr)
  if(xrank.eq.px-1)&
   call MPI_WAIT(req2(2), status_arr2(1,2), ierr)
  if(xrank.eq.0)&
   call MPI_WAIT(req2(3), status_arr2(1,3), ierr)
  if(xrank.eq.px-1)&
   call MPI_WAIT(req2(4), status_arr2(1,4), ierr)
  if(yrank.eq.0)&
   call MPI_WAIT(req2(5), status_arr2(1,5), ierr)
  if(yrank.eq.py-1)&
   call MPI_WAIT(req2(6), status_arr2(1,6), ierr)
 elseif(zoneid.eq.3)then
  if(xrank.eq.0)&
   call MPI_WAIT(req2(1), status_arr2(1,1), ierr)
  if(xrank.eq.px-1)&
   call MPI_WAIT(req2(2), status_arr2(1,2), ierr)
  if(xrank.eq.0)&
   call MPI_WAIT(req2(3), status_arr2(1,3), ierr)
  if(xrank.eq.px-1)&
   call MPI_WAIT(req2(4), status_arr2(1,4), ierr)
 elseif(zoneid.eq.4)then
  if(xrank.eq.0)&
   call MPI_WAIT(req2(1), status_arr2(1,1), ierr)
  if(xrank.eq.px-1)&
   call MPI_WAIT(req2(2), status_arr2(1,2), ierr)
  if(xrank.eq.0)&
   call MPI_WAIT(req2(3), status_arr2(1,3), ierr)
  if(xrank.eq.px-1)&
   call MPI_WAIT(req2(4), status_arr2(1,4), ierr)
  if(yrank.eq.0)&
   call MPI_WAIT(req2(5), status_arr2(1,5), ierr)
  if(yrank.eq.py-1)&
    call MPI_WAIT(req2(6), status_arr2(1,6), ierr)
 elseif(zoneid.eq.5)then
  if(xrank.eq.0)&
   call MPI_WAIT(req2(1), status_arr2(1,1), ierr)
  if(xrank.eq.px-1)&
   call MPI_WAIT(req2(2), status_arr2(1,2), ierr)
 elseif(zoneid.eq.6)then
  if(xrank.eq.0)&
   call MPI_WAIT(req2(1), status_arr2(1,1), ierr)
  if(xrank.eq.px-1)&
   call MPI_WAIT(req2(2), status_arr2(1,2), ierr)
 endif

! Switching indices for reverse aligned arrays
 if(zoneid.eq.5)then
  if(xrank.eq.0)then
   tmp1=fp
   do k=1,mr
    do j=1,my
     do i=1,nst
      fp(i,j,k) = tmp1(i,my-j+1,k)
     enddo
    enddo
   enddo
  endif
 elseif(zoneid.eq.6)then
  if(xrank.eq.px-1)then
   tmp2=fn
   do k=1,mr
    do j=1,my
     do i=1,nst
      fn(i,j,k) = tmp2(i,my-j+1,k)
     enddo
    enddo
   enddo
  endif
 endif

! LEFT BOUNDARY POINTS 

 do k=1,mr
  do j=1,my
   do i=1,nst
    dx_f(i,j,k) = 0.0

    do iter=-nst,-i
     dx_f(i,j,k) = dx_f(i,j,k) + fac * coeffs1(iter) * fp(-i-iter+1,j,k)
    enddo

    do iter=-i+1,nst
     dx_f(i,j,k) = dx_f(i,j,k) + fac * coeffs1(iter) * f(i+iter,j,k)
    enddo

   enddo 
  enddo
 enddo

! RIGHT BOUNDARY POINTS
 
 do k=1,mr 
  do j=1,my
   do i=mx-nst+1,mx
    dx_f(i,j,k) = 0.0

    do iter=-nst,mx-i
     dx_f(i,j,k) = dx_f(i,j,k) + fac * coeffs1(iter) * f(i+iter,j,k)
    enddo

    do iter=mx-i+1,nst
     dx_f(i,j,k) = dx_f(i,j,k) + fac * coeffs1(iter) * fn(i+iter-mx,j,k)
    enddo

   enddo
  enddo
 enddo

END SUBROUTINE DDX

!----------------------------------------------------------------------------------

SUBROUTINE DDR(f, dr_f)

 implicit none

 include 'mpif.h'

 integer i,j,k,iter
 real*8 f(mx,my,mr), dr_f(mx,my,mr)
 real*8 fac, fac2

 real*8 temp(mx,my,nst,2), fp(mx,my,nst), fn(mx,my,nst)

 integer status_arr(MPI_STATUS_SIZE,4), ierr, req(4)

! 2D case
 if(nr.eq.1)then 
  dr_f = 0.0
  return
 endif

 dr_f = 0.0
 fac  = 1.0D0/dchi

! Putting first nst and last nst data planes into temp arrays for sending 
 do k=1,nst
  do j=1,my
   do i=1,mx
    temp(i,j,k,1) = f(i,j,k)
    temp(i,j,k,2) = f(i,j,mr-k+1)
   enddo
  enddo
 enddo 

! Intra-communication
! Transferring data from neighbouring processors
 if(rrank.gt.0)&
 call MPI_IRECV(fp, mx*my*nst, MPI_REAL8, source_right_r,&
                source_right_r, MPI_XYR_COMM, req(1), ierr)               ! Recv last nst planes
 if(rrank.lt.pr-1)&
 call MPI_IRECV(fn, mx*my*nst, MPI_REAL8, source_left_r,&
                source_left_r, MPI_XYR_COMM, req(2), ierr)                ! Recv first nst planes
 if(rrank.lt.pr-1)&
 call MPI_ISEND(temp(1,1,1,2), mx*my*nst, MPI_REAL8, dest_right_r,&
                zonerank, MPI_XYR_COMM, req(3), ierr)                     ! Send last nst planes
 if(rrank.gt.0)&
 call MPI_ISEND(temp(1,1,1,1), mx*my*nst, MPI_REAL8, dest_left_r,&
                zonerank, MPI_XYR_COMM, req(4), ierr)                     ! Send first nst planes

! Compute interior in the meantime 
! ALL INTERIOR POINTS

 do k=nst+1,mr-nst
  do j=1,my
   do i=1,mx
    dr_f(i,j,k) = 0.0
    fac2 = fac/drdchi(k)

    do iter=-nst,nst
     dr_f(i,j,k) = dr_f(i,j,k) + fac2 * coeffs1(iter) * f(i,j,k+iter)
    enddo
 
  enddo
  enddo
 enddo

! Waiting for intra communication to complete before computing at boundaries
 if(rrank.eq.0)then
  call MPI_WAITALL(2, req(2), status_arr(1,2), ierr)
 elseif(rrank.eq.pr-1)then
  call MPI_WAIT(req(1), status_arr(1,1), ierr)
  call MPI_WAIT(req(4), status_arr(1,4), ierr)
 else
  call MPI_WAITALL(4, req, status_arr, ierr)
 endif

! LEFT BOUNDARY POINTS 

 if(rrank.gt.0)then

  do k=1,nst
   do j=1,my
    do i=1,mx
     dr_f(i,j,k) = 0.0
     fac2 = fac/drdchi(k) 

     do iter=-nst,-k
      dr_f(i,j,k) = dr_f(i,j,k) + fac2 * coeffs1(iter) * fp(i,j,-k-iter+1)
     enddo
 
     do iter=-k+1,nst
      dr_f(i,j,k) = dr_f(i,j,k) + fac2 * coeffs1(iter) * f(i,j,k+iter)
     enddo

    enddo
   enddo
  enddo

 else 

  do k=1,nst
   do j=1,my
    do i=1,mx
     dr_f(i,j,k) = 0.0
     fac2 = fac/drdchi(k)

     do iter=1,nstencil
      dr_f(i,j,k) = dr_f(i,j,k) + fac2 * bdy_coeffs1(k,iter) * f(i,j,iter)
     enddo

    enddo
   enddo
  enddo

 endif

 if(rrank.lt.pr-1)then
 
  do k=mr-nst+1,mr
   do j=1,my
    do i=1,mx
     dr_f(i,j,k) = 0.0
     fac2 = fac/drdchi(k) 

     do iter=-nst,mr-k
      dr_f(i,j,k) = dr_f(i,j,k) + fac2 * coeffs1(iter) * f(i,j,k+iter)
     enddo  

     do iter=mr-k+1,nst
      dr_f(i,j,k) = dr_f(i,j,k) + fac2 * coeffs1(iter) * fn(i,j,k+iter-mr)
     enddo
 
    enddo
   enddo
  enddo

 else

  do k=mr-nst+1,mr
   do j=1,my
    do i=1,mx
     dr_f(i,j,k) = 0.0
     fac2 = -fac/drdchi(k)

     do iter=1,nstencil
      dr_f(i,j,k) = dr_f(i,j,k) + fac2 * bdy_coeffs1(mr-k+1,nstencil-iter+1) * f(i,j,mr-nstencil+iter)
     enddo
 
    enddo
   enddo
  enddo

 endif

END SUBROUTINE DDR

!---------------------------------------------------------------------------------

END MODULE DERIVATIVE

