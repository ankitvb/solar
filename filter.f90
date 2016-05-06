MODULE FILTER

USE INIT
USE INTERPOLATION

CONTAINS

!--------------------------------------------------------------------------------------------------
SUBROUTINE FILTER_EXP(q, qf)

 implicit none

 real*8, dimension(mx,my,mr,5) :: q, qf, qf_temp1, qf_temp2

 integer i,j,k,m

 real*8, dimension(mx,noverlap,mr) :: interpolated_south, interpolated_north
 real*8, dimension(noverlap,my,mr) :: interpolated_south_transposed, interpolated_north_transposed
 real*8, dimension(noverlap,my,mr) :: interpolated_east, interpolated_west
 real*8, dimension(mx,noverlap,mr) :: interpolated_east_transposed, interpolated_west_transposed

! INTERPOLATION FROM  i --> j
 do m=1,5

  call INTERPOLATE_ALONG_X(q(:,2:noverlap+1,:,m), interpolated_south)
  call INTERPOLATE_ALONG_X(q(:,my-1:my-noverlap:-1,:,m), interpolated_north)
  call INTERPOLATE_ALONG_Y(q(2:noverlap+1,:,:,m), interpolated_west)
  call INTERPOLATE_ALONG_Y(q(mx-1:mx-noverlap:-1,:,:,m), interpolated_east)

! Transposing relevant arrays to make their shape conform to the one expected by the filter routines
 do k=1,mr
  do j=1,noverlap
   do i=1,mx
     interpolated_east_transposed(i,j,k)  = interpolated_east(j,i,k)
     interpolated_west_transposed(i,j,k)  = interpolated_west(j,i,k)
     interpolated_north_transposed(j,i,k) = interpolated_north(i,j,k)
     interpolated_south_transposed(j,i,k) = interpolated_south(i,j,k)
   enddo
  enddo
 enddo

  call FILTER_XI(q(:,:,:,m), qf_temp1(:,:,:,m), interpolated_south_transposed, interpolated_north_transposed, interpolated_west, interpolated_east)
  call FILTER_ETA(qf_temp1(:,:,:,m), qf_temp2(:,:,:,m), interpolated_south, interpolated_north, interpolated_west_transposed, interpolated_east_transposed)
  call FILTER_R(qf_temp2(:,:,:,m), qf(:,:,:,m))

 enddo
 
END SUBROUTINE FILTER_EXP
!--------------------------------------------------------------------------------------------------
SUBROUTINE FILTER_XI(g, gf, interpol_low, interpol_high, interpol_left, interpol_right)
 implicit none

 include 'mpif.h'

 integer i, j, k, iter
 real*8 g(mx,my,mr), gf(mx,my,mr)
 real*8 interpol_left(noverlap,my,mr), interpol_right(noverlap,my,mr)
 real*8 interpol_low(noverlap,my,mr), interpol_high(noverlap,my,mr)

 real*8 temp(nst,my,mr,2), gp(nst,my,mr), gn(nst,my,mr)
 real*8 tmp1(nst,my,mr), tmp2(nst,my,mr)

 integer ierr
 integer status_arr1(MPI_STATUS_SIZE,4), req1(4)
 integer status_arr2(MPI_STATUS_SIZE,6), req2(6)

! Putting first nst and last nst data planes into temp arrays for sending
 do k=1,mr
  do j=1,my
   do i=1,nst
    temp(i,j,k,1) = g(i,j,k)
    temp(i,j,k,2) = g(mx-i+1,j,k)
   enddo
  enddo
 enddo

! Intra-communication
! Transferring data from neighbouring processors
 if(xrank.gt.0)&
 call MPI_IRECV(gp, nst*my*mr, MPI_REAL8, source_right_x,&
                source_right_x, MPI_XYR_COMM, req1(1), ierr)           ! Recv last nst planes
 if(xrank.lt.px-1)&
 call MPI_IRECV(gn, nst*my*mr, MPI_REAL8, source_left_x,&
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
  call MPI_IRECV(gp, nst*my*mr, MPI_REAL8, icom_source_right_x,&
                 icom_source_right_x, MPI_INTERCOMM(1,4), req2(1), ierr)
  if(xrank.eq.px-1)&
  call MPI_IRECV(gn, nst*my*mr, MPI_REAL8, icom_source_left_x,&
                 icom_source_left_x,  MPI_INTERCOMM(1,2), req2(2), ierr)
  if(xrank.eq.0)&
  call MPI_ISEND(interpol_left, nst*my*mr, MPI_REAL8, icom_dest_left_x,&
                 zonerank, MPI_INTERCOMM(1,4), req2(3), ierr)
  if(xrank.eq.px-1)&
  call MPI_ISEND(interpol_right, nst*my*mr, MPI_REAL8, icom_dest_right_x,&
                 zonerank, MPI_INTERCOMM(1,2), req2(4), ierr)
 elseif(zoneid.eq.2)then
  if(xrank.eq.0)&
  call MPI_IRECV(gp, nst*my*mr, MPI_REAL8, icom_source_right_x,&
                 icom_source_right_x, MPI_INTERCOMM(1,2), req2(1), ierr)
  if(xrank.eq.px-1)&
  call MPI_IRECV(gn, nst*my*mr, MPI_REAL8, icom_source_left_x,&
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
  call MPI_IRECV(gp, nst*my*mr, MPI_REAL8, icom_source_right_x,&
                 icom_source_right_x, MPI_INTERCOMM(2,3), req2(1), ierr)
  if(xrank.eq.px-1)&
  call MPI_IRECV(gn, nst*my*mr, MPI_REAL8, icom_source_left_x,&
                  icom_source_left_x,  MPI_INTERCOMM(3,4), req2(2), ierr)
  if(xrank.eq.0)&
  call MPI_ISEND(interpol_left, nst*my*mr, MPI_REAL8, icom_dest_left_x,&
                 zonerank, MPI_INTERCOMM(2,3), req2(3), ierr)
  if(xrank.eq.px-1)&
  call MPI_ISEND(interpol_right, nst*my*mr, MPI_REAL8, icom_dest_right_x,&
                 zonerank, MPI_INTERCOMM(3,4), req2(4), ierr)
 elseif(zoneid.eq.4)then
  if(xrank.eq.0)&
  call MPI_IRECV(gp, nst*my*mr, MPI_REAL8, icom_source_right_x,&
                 icom_source_right_x, MPI_INTERCOMM(3,4), req2(1), ierr)
  if(xrank.eq.px-1)&
  call MPI_IRECV(gn, nst*my*mr, MPI_REAL8, icom_source_left_x,&
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
  call MPI_IRECV(gp, nst*my*mr, MPI_REAL8, icom_source_right_x, &
                 icom_source_right_x, MPI_INTERCOMM(4,5), req2(1), ierr)
  if(xrank.eq.px-1)&
  call MPI_IRECV(gn, nst*my*mr, MPI_REAL8, icom_source_left_x, &
                 icom_source_left_x, MPI_INTERCOMM(2,5), req2(2), ierr)
 elseif(zoneid.eq.6)then
  if(xrank.eq.0)&
  call MPI_IRECV(gp, nst*my*mr, MPI_REAL8, icom_source_right_x,&
                 icom_source_right_x, MPI_INTERCOMM(4,6), req2(1), ierr)
  if(xrank.eq.px-1)&
  call MPI_IRECV(gn, nst*my*mr, MPI_REAL8, icom_source_left_x, &
                 icom_source_left_x, MPI_INTERCOMM(2,6), req2(2), ierr)
 endif

! Compute interior in the meantime
! ALL INTERIOR POINTS 

 do k=1,mr
  do j=1,my
   do i=nst+1,mx-nst
    gf(i,j,k) = 0.0
    do iter=-nst,nst
      gf(i,j,k) = gf(i,j,k) + filt_coeffs(iter) * g(i+iter,j,k)
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
   tmp1=gp
   do k=1,mr
    do i=1,mx
     do j=1,my
      gp(i,j,k) = tmp1(i,my-j+1,k)
     enddo
    enddo
   enddo
  endif
 elseif(zoneid.eq.6)then
  if(xrank.eq.px-1)then
   tmp2=gn
   do k=1,mr
    do i=1,mx
     do j=1,my
      gn(i,j,k) = tmp2(i,my-j+1,k)
     enddo
    enddo
   enddo
  endif
 endif

! LEFT BOUNDARY POINTS 
 do k=1,mr
  do j=1,my
   do i=1,nst
    gf(i,j,k) = 0.0

    do iter=-nst,-i
     gf(i,j,k) = gf(i,j,k) + filt_coeffs(iter) * gp(-i-iter+1,j,k)
    enddo

    do iter=-i+1,nst
     gf(i,j,k) = gf(i,j,k) + filt_coeffs(iter) * g(i+iter,j,k)
    enddo

   enddo
  enddo
 enddo

! RIGHT BOUNDARY POINTS
 do k=1,mr 
  do j=1,my
   do i=mx-nst+1,mx
    gf(i,j,k) = 0.0
 
    do iter=-nst,mx-i
     gf(i,j,k) = gf(i,j,k) + filt_coeffs(iter) * g(i+iter,j,k)
    enddo

    do iter=mx-i+1,nst
     gf(i,j,k) = gf(i,j,k) + filt_coeffs(iter) * gn(i+iter-mx,j,k)
    enddo

   enddo
  enddo
 enddo


END SUBROUTINE FILTER_XI

!--------------------------------------------------------------------------------------------------
SUBROUTINE FILTER_ETA(g, gf, interpol_low, interpol_high, interpol_left, interpol_right)

 implicit none

 include 'mpif.h'

 integer i, j, k, iter
 real*8 g(mx,my,mr), gf(mx,my,mr)
 real*8 interpol_low(mx,noverlap,mr), interpol_high(mx,noverlap,mr)
 real*8 interpol_left(mx,noverlap,mr), interpol_right(mx,noverlap,mr)

 real*8 temp(mx,nst,mr,2), gp(mx,nst,mr), gn(mx,nst,mr)
 real*8 tmp1(mx,nst,mr), tmp2(mx,nst,mr)

 integer ierr
 integer status_arr1(MPI_STATUS_SIZE,4), req1(4)
 integer status_arr2(MPI_STATUS_SIZE,6), req2(6)

! Putting first nst and last nst data planes into temp arrays for sending
 do k=1,mr
  do j=1,nst
   do i=1,mx
    temp(i,j,k,1) = g(i,j,k)
    temp(i,j,k,2) = g(i,my-j+1,k)
   enddo
  enddo
 enddo

! Intra-communication
! Transferring data from neighbouring processors
 if(yrank.gt.0)&
 call MPI_IRECV(gp, mx*nst*mr, MPI_REAL8, source_right_y,&
                source_right_y, MPI_XYR_COMM, req1(1), ierr)               ! Recv last nst planes
 if(yrank.lt.py-1)&
 call MPI_IRECV(gn, mx*nst*mr, MPI_REAL8, source_left_y,&
                source_left_y, MPI_XYR_COMM, req1(2), ierr)                ! Recv first nst planes
 if(yrank.lt.py-1)&
 call MPI_ISEND(temp(1,1,1,2), mx*nst*mr, MPI_REAL8, dest_right_y,&
                zonerank, MPI_XYR_COMM, req1(3), ierr)                     ! Send last nst planes
 if(yrank.gt.0)&
 call MPI_ISEND(temp(1,1,1,1), mx*nst*mr, MPI_REAL8, dest_left_y,&
                zonerank, MPI_XYR_COMM, req1(4), ierr)                     ! Send first nst planes

! Inter-communication
 if(zoneid.eq.1)then
  if(yrank.eq.0)&
  call MPI_IRECV(gp, mx*nst*mr, MPI_REAL8, icom_source_right_y,&
                 icom_source_right_y, MPI_INTERCOMM(1,6), req2(1), ierr)
  if(yrank.eq.py-1)&
  call MPI_IRECV(gn, mx*nst*mr, MPI_REAL8, icom_source_left_y,&
                 icom_source_left_y,  MPI_INTERCOMM(1,5), req2(2), ierr)
  if(yrank.eq.0)&
  call MPI_ISEND(interpol_low, mx*nst*mr, MPI_REAL8, icom_dest_left_y,&
                 zonerank, MPI_INTERCOMM(1,6), req2(3), ierr)
  if(yrank.eq.py-1)&
  call MPI_ISEND(interpol_high, mx*nst*mr, MPI_REAL8, icom_dest_right_y,&
                 zonerank, MPI_INTERCOMM(1,5), req2(4), ierr)
 elseif(zoneid.eq.2)then
  if(yrank.eq.0)&
  call MPI_IRECV(gp, mx*nst*mr, MPI_REAL8, icom_source_right_y,&
                 icom_source_right_y, MPI_INTERCOMM(2,6), req2(1), ierr)
  if(yrank.eq.py-1)&
  call MPI_IRECV(gn, mx*nst*mr, MPI_REAL8, icom_source_left_y,&
                 icom_source_left_y,  MPI_INTERCOMM(2,5), req2(2), ierr)
 elseif(zoneid.eq.3)then
  if(yrank.eq.0)&
  call MPI_IRECV(gp,  mx*nst*mr, MPI_REAL8, icom_source_right_y,&
                 icom_source_right_y, MPI_INTERCOMM(3,6), req2(1), ierr)
  if(yrank.eq.py-1)&
  call MPI_IRECV(gn, mx*nst*mr, MPI_REAL8, icom_source_left_y,&
                 icom_source_left_y,  MPI_INTERCOMM(3,5), req2(2), ierr)
  if(yrank.eq.0)&
  call MPI_ISEND(interpol_low, mx*nst*mr, MPI_REAL8, icom_dest_left_y,&
                 zonerank, MPI_INTERCOMM(3,6), req2(3), ierr)
  if(yrank.eq.py-1)&
  call MPI_ISEND(interpol_high, mx*nst*mr, MPI_REAL8, icom_dest_right_y,&
                 zonerank, MPI_INTERCOMM(3,5), req2(4), ierr)
 elseif(zoneid.eq.4)then
  if(yrank.eq.0)&
  call MPI_IRECV(gp,  mx*nst*mr, MPI_REAL8, icom_source_right_y,&
                 icom_source_right_y, MPI_INTERCOMM(4,6), req2(1), ierr)
  if(yrank.eq.py-1)&
  call MPI_IRECV(gn, mx*nst*mr, MPI_REAL8, icom_source_left_y,&
                 icom_source_left_y,  MPI_INTERCOMM(4,5), req2(2), ierr)
 elseif(zoneid.eq.5)then
  if(yrank.eq.0)&
  call MPI_IRECV(gp,  mx*nst*mr, MPI_REAL8, icom_source_right_y,&
                  icom_source_right_y, MPI_INTERCOMM(1,5), req2(1), ierr)
  if(yrank.eq.py-1)&
  call MPI_IRECV(gn,  mx*nst*mr, MPI_REAL8, icom_source_left_y,&
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
  call MPI_IRECV(gp, mx*nst*mr, MPI_REAL8, icom_source_right_y,&
                 icom_source_right_y, MPI_INTERCOMM(3,6), req2(1), ierr)
  if(yrank.eq.py-1)&
  call MPI_IRECV(gn, mx*nst*mr, MPI_REAL8, icom_source_left_y,&
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
  do i=1,mx
   do j=nst+1,my-nst
    gf(i,j,k) = 0.0
    do iter=-nst,nst
     gf(i,j,k) = gf(i,j,k) + filt_coeffs(iter) * g(i,j+iter,k)
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

! Waiting for inter communication to complete
! Wait to recieve planes
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
   tmp1=gp
   do k=1,mr
    do j=1,my
     do i=1,mx
      gp(i,j,k) = tmp1(mx-i+1,j,k)
     enddo
    enddo
   enddo
  endif
 elseif(zoneid.eq.3)then
  tmp1=gp; tmp2=gn
  if(yrank.eq.0)then
   do k=1,mr
    do j=1,my
     do i=1,mx
      gp(i,j,k) = tmp1(mx-i+1,j,k)
     enddo
    enddo
   enddo
  endif
  if(yrank.eq.py-1)then
   do k=1,mr
    do j=1,my
     do i=1,mx
      gn(i,j,k) = tmp2(mx-i+1,j,k)
     enddo
    enddo
   enddo
  endif
 elseif(zoneid.eq.4)then
  tmp2=gn
  if(yrank.eq.py-1)then
   do k=1,mr
    do j=1,my
     do i=1,mx
      gn(i,j,k) = tmp2(mx-i+1,j,k)
     enddo
    enddo
   enddo
  endif
 elseif(zoneid.eq.5)then
  tmp2=gn
  if(yrank.eq.py-1)then
   do k=1,mr
    do j=1,my
     do i=1,mx
      gn(i,j,k) = tmp2(mx-i+1,j,k)
     enddo
    enddo
   enddo
  endif
 elseif(zoneid.eq.6)then
  tmp1=gp
  if(yrank.eq.0)then
   do k=1,mr
    do j=1,my
     do i=1,mx
      gp(i,j,k) = tmp1(mx-i+1,j,k)
     enddo
    enddo
   enddo
  endif
 endif

! LOWER BOUNDARY POINTS
  do k=1,mr
   do i=1,mx
    do j=1,nst
     gf(i,j,k) = 0.0

     do iter=-nst,-j
      gf(i,j,k) = gf(i,j,k) + filt_coeffs(iter) * gp(i,-iter-j+1,k)
     enddo

     do iter=-j+1,nst
      gf(i,j,k) = gf(i,j,k) + filt_coeffs(iter) * g(i,iter+j,k)
     enddo

    enddo
   enddo
  enddo

! UPPER BOUNDARY POINTS
  do k=1,mr
   do i=1,mx
    do j=my-nst+1,my
     gf(i,j,k) = 0.0
  
     do iter=-nst,my-j
      gf(i,j,k) = gf(i,j,k) + filt_coeffs(iter) * g(i,iter+j,k)
     enddo
  
     do iter=my-j+1,nst
      gf(i,j,k) = gf(i,j,k) + filt_coeffs(iter) * gn(i,iter+j-my,k)
     enddo
  
    enddo
   enddo
  enddo

END SUBROUTINE FILTER_ETA
!--------------------------------------------------------------------------------------------------
SUBROUTINE FILTER_R(g, gf)

 implicit none

 include 'mpif.h'

 integer i, j, k, iter
 real*8 g(mx,my,mr), gf(mx,my,mr)

 real*8 temp(mx,my,nstr,2), gp(mx,my,nstr), gn(mx,my,nstr)

 integer status_arr(MPI_STATUS_SIZE,4), ierr, req(4)

! 2D case
 if(nr.eq.1)then
  gf = g 
  return
 endif

 gf = 0.0

! Putting first nstr and last nstr data planes into temp arrays for sending 
 do k=1,nstr
  do j=1,my
   do i=1,mx
    temp(i,j,k,1) = g(i,j,k)
    temp(i,j,k,2) = g(i,j,mr-k+1)
   enddo
  enddo
 enddo

! Intra-communication
! Transferring data from neighbouring processors
 if(rrank.gt.0)&
 call MPI_IRECV(gp, mx*my*nstr, MPI_REAL8, source_right_r,&
                source_right_r, MPI_XYR_COMM, req(1), ierr)               ! Recv last nstr planes
 if(rrank.lt.pr-1)& 
 call MPI_IRECV(gn, mx*my*nstr, MPI_REAL8, source_left_r,&
                source_left_r, MPI_XYR_COMM, req(2), ierr)                ! Recv first nstr planes
 if(rrank.lt.pr-1)&
 call MPI_ISEND(temp(1,1,1,2), mx*my*nstr, MPI_REAL8, dest_right_r,&
                zonerank, MPI_XYR_COMM, req(3), ierr)                     ! Send last nstr planes
 if(rrank.gt.0)&
 call MPI_ISEND(temp(1,1,1,1), mx*my*nstr, MPI_REAL8, dest_left_r,&
                zonerank, MPI_XYR_COMM, req(4), ierr)                     ! Send first nstr planes

! Compute interior in the meantime 
! ALL INTERIOR POINTS
 
 do k=nstr+1,mr-nstr
  do j=1,my
   do i=1,mx
    gf(i,j,k) = 0.0

    do iter=-nstr,nstr
     gf(i,j,k) = gf(i,j,k) + filt_coeffs_r(iter) * g(i,j,k+iter)
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

  do k=1,nstr
   do j=1,my
    do i=1,mx
     gf(i,j,k) = 0.0

     do iter=-nstr,-k
      gf(i,j,k) = gf(i,j,k) + filt_coeffs_r(iter) * gp(i,j,-k-iter+1)
     enddo

     do iter=-k+1,nstr
      gf(i,j,k) = gf(i,j,k) + filt_coeffs_r(iter) * g(i,j,k+iter)
     enddo

    enddo
   enddo
  enddo

 else

  do k=1,nstr
   do j=1,my
    do i=1,mx
     gf(i,j,k) = 0.0

     do iter=1,4
!      gf(i,j,k) = gf(i,j,k) + bdy_filt_coeffs(k,iter) * g(i,j,iter)
     enddo

    enddo
   enddo
  enddo

 endif

 if(rrank.lt.pr-1)then

  do k=mr-nstr+1,mr
   do j=1,my
    do i=1,mx
     gf(i,j,k) = 0.0

     do iter=-nstr,mr-k
      gf(i,j,k) = gf(i,j,k) + filt_coeffs_r(iter) * g(i,j,k+iter)
     enddo

     do iter=mr-k+1,nstr
      gf(i,j,k) = gf(i,j,k) + filt_coeffs_r(iter) * gn(i,j,k+iter-mr)
     enddo

    enddo
   enddo
  enddo

 else

  do k=mr-nstr+1,mr
   do j=1,my
    do i=1,mx
     gf(i,j,k) = 0.0

     do iter=1,4
!      gf(i,j,k) = gf(i,j,k) + bdy_filt_coeffs(mr-k+1,4-iter+1) * g(i,j,mr-4+iter)
     enddo

    enddo
   enddo
  enddo

 endif

 return

END SUBROUTINE FILTER_R

!--------------------------------------------------------------------------------------------

END MODULE FILTER
