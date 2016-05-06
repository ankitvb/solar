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

 integer status_arr(MPI_STATUS_SIZE,6), ierr, req(6)

 fac = 1.0D0/dxi

! Putting first nst and last nst data planes into temp arrays for sending
 temp(:,:,:,1) = f(1:nst,:,:)
 temp(:,:,:,2) = f(mx:mx-nst+1:-1,:,:)

! Intra-communication
! Transferring data from neighbouring processors
 if(xrank.gt.0)&
 call MPI_IRECV(fp, nst*my*mr, MPI_REAL8, source_right_x,&
                source_right_x, MPI_XYR_COMM, req(1), ierr)           ! Recv last nst planes
 if(xrank.lt.px-1)&
 call MPI_IRECV(fn, nst*my*mr, MPI_REAL8, source_left_x,&
                source_left_x, MPI_XYR_COMM, req(2), ierr)            ! Recv first nst planes
 if(xrank.lt.px-1)&
 call MPI_ISEND(temp(1,1,1,2), nst*my*mr, MPI_REAL8, dest_right_x,&
                zonerank, MPI_XYR_COMM, req(3), ierr)                 ! Send last nst planes
 if(xrank.gt.0)&
 call MPI_ISEND(temp(1,1,1,1), nst*my*mr, MPI_REAL8, dest_left_x,&
                zonerank, MPI_XYR_COMM, req(4), ierr)                 ! Send first nst planes 

! Inter-communication 
 if(zoneid.eq.1)then
  if(xrank.eq.0)&
  call MPI_IRECV(fp, nst*my*mr, MPI_REAL8, icom_source_right_x,&
                 icom_source_right_x, MPI_INTERCOMM(1,4), req(1), ierr)
  if(xrank.eq.px-1)&
  call MPI_IRECV(fn, nst*my*mr, MPI_REAL8, icom_source_left_x,& 
                 icom_source_left_x,  MPI_INTERCOMM(1,2), req(2), ierr)
  if(xrank.eq.0)&
  call MPI_ISEND(interpol_left, nst*my*mr, MPI_REAL8, icom_dest_left_x,&
                 zonerank, MPI_INTERCOMM(1,4), req(3), ierr)
  if(xrank.eq.px-1)&
  call MPI_ISEND(interpol_right, nst*my*mr, MPI_REAL8, icom_dest_right_x,&
                 zonerank, MPI_INTERCOMM(1,2), req(4), ierr)
 elseif(zoneid.eq.2)then
  if(xrank.eq.0)&
  call MPI_IRECV(fp, nst*my*mr, MPI_REAL8, icom_source_right_x,&
                 icom_source_right_x, MPI_INTERCOMM(1,2), req(1), ierr)
  if(xrank.eq.px-1)&
  call MPI_IRECV(fn, nst*my*mr, MPI_REAL8, icom_source_left_x,&
                 icom_source_left_x, MPI_INTERCOMM(2,3), req(2), ierr)
  if(xrank.eq.0)&
  call MPI_ISEND(interpol_left, nst*my*mr, MPI_REAL8, icom_dest_left_x,&
                 zonerank, MPI_INTERCOMM(1,2), req(3), ierr)
  if(xrank.eq.px-1)&
  call MPI_ISEND(interpol_right, nst*my*mr, MPI_REAL8, icom_dest_right_x,&
                 zonerank, MPI_INTERCOMM(2,3), req(4), ierr)
  if(yrank.eq.0)&
  call MPI_ISEND(interpol_low, nst*my*mr, MPI_REAL8, icom_dest_left_y,&
                 zonerank, MPI_INTERCOMM(2,6), req(5), ierr)
  if(yrank.eq.py-1)&
  call MPI_ISEND(interpol_high, nst*my*mr, MPI_REAL8, icom_dest_right_y,&
                 zonerank, MPI_INTERCOMM(2,5), req(6), ierr)
 elseif(zoneid.eq.3)then
  if(xrank.eq.0)&
  call MPI_IRECV(fp, nst*my*mr, MPI_REAL8, icom_source_right_x,&
                 icom_source_right_x, MPI_INTERCOMM(2,3), req(1), ierr)
  if(xrank.eq.px-1)&
  call MPI_IRECV(fn, nst*my*mr, MPI_REAL8, icom_source_left_x,&
                  icom_source_left_x,  MPI_INTERCOMM(3,4), req(2), ierr)
  if(xrank.eq.0)&
  call MPI_ISEND(interpol_left, nst*my*mr, MPI_REAL8, icom_dest_left_x,&
                 zonerank, MPI_INTERCOMM(2,3), req(3), ierr)
  if(xrank.eq.px-1)&
  call MPI_ISEND(interpol_right, nst*my*mr, MPI_REAL8, icom_dest_right_x,&
                 zonerank, MPI_INTERCOMM(3,4), req(4), ierr)
 elseif(zoneid.eq.4)then
  if(xrank.eq.0)&
  call MPI_IRECV(fp, nst*my*mr, MPI_REAL8, icom_source_right_x,&
                 icom_source_right_x, MPI_INTERCOMM(3,4), req(1), ierr)
  if(xrank.eq.px-1)&
  call MPI_IRECV(fn, nst*my*mr, MPI_REAL8, icom_source_left_x,&
                 icom_source_left_x, MPI_INTERCOMM(1,4), req(2), ierr)
  if(xrank.eq.0)&
  call MPI_ISEND(interpol_left, nst*my*mr, MPI_REAL8, icom_dest_left_x,&
                 zonerank, MPI_INTERCOMM(3,4), req(3), ierr)
  if(xrank.eq.px-1)&
  call MPI_ISEND(interpol_right, nst*my*mr, MPI_REAL8, icom_dest_right_x,&
                 zonerank, MPI_INTERCOMM(1,4), req(4), ierr)
  if(yrank.eq.0)&
  call MPI_ISEND(interpol_low, nst*my*mr, MPI_REAL8, icom_dest_left_y,&
                 zonerank, MPI_INTERCOMM(4,6), req(5), ierr)
  if(yrank.eq.py-1)&
  call MPI_ISEND(interpol_high, nst*my*mr, MPI_REAL8, icom_dest_right_y,&
                 zonerank, MPI_INTERCOMM(4,5), req(6), ierr)
 elseif(zoneid.eq.5)then
  if(xrank.eq.0)&
  call MPI_IRECV(fp, nst*my*mr, MPI_REAL8, icom_source_right_x, &
                 icom_source_right_x, MPI_INTERCOMM(4,5), req(1), ierr)
  if(xrank.eq.px-1)&
  call MPI_IRECV(fn, nst*my*mr, MPI_REAL8, icom_source_left_x, &
                 icom_source_left_x, MPI_INTERCOMM(2,5), req(2), ierr)
 elseif(zoneid.eq.6)then
  if(xrank.eq.0)&
  call MPI_IRECV(fp, nst*my*mr, MPI_REAL8, icom_source_right_x,&
                 icom_source_right_x, MPI_INTERCOMM(4,6), req(1), ierr)
  if(xrank.eq.px-1)&
  call MPI_IRECV(fn, nst*my*mr, MPI_REAL8, icom_source_left_x, &
                 icom_source_left_x, MPI_INTERCOMM(2,6), req(2), ierr)
 endif

! Compute interior in the meantime
! ALL INTERIOR POINTS 

 do i=nst+1,mx-nst
  dx_f(i,:,:) = 0.0
  do j=-nst,nst
    dx_f(i,:,:) = dx_f(i,:,:) + fac * coeffs1(j) * f(i+j,:,:)
  enddo
 enddo

! Waiting for intra and inter communication to complete before computing at boundaries
 call MPI_WAITALL(4, req, status_arr, ierr)
 if((zoneid.eq.2).or.(zoneid.eq.4))then
  call MPI_WAITALL(2, req(5), status_arr(1,5), ierr)
 endif

! Switching indices for reverse aligned arrays
 if(zoneid.eq.5)then
  if(xrank.eq.0)then
   tmp1=fp
   do j=1,my
    fp(:,j,:) = tmp1(:,my-j+1,:)
   enddo
  endif
 elseif(zoneid.eq.6)then
  if(xrank.eq.px-1)then
   tmp2=fn
   do j=1,my
    fn(:,j,:) = tmp2(:,my-j+1,:)
   enddo
  endif
 endif

! LEFT BOUNDARY POINTS 
 
 do i=1,nst
  dx_f(i,:,:) = 0.0

  do j=-nst,-i
   dx_f(i,:,:) = dx_f(i,:,:) + fac * coeffs1(j) * fp(-i-j+1,:,:)
  enddo

  do j=-i+1,nst
   dx_f(i,:,:) = dx_f(i,:,:) + fac * coeffs1(j) * f(i+j,:,:)
  enddo

 enddo

! RIGHT BOUNDARY POINTS
 do i=mx-nst+1,mx
  dx_f(i,:,:) = 0.0

  do j=-nst,mx-i
   dx_f(i,:,:) = dx_f(i,:,:) + fac * coeffs1(j) * f(i+j,:,:)
  enddo

  do j=mx-i+1,nst
   dx_f(i,:,:) = dx_f(i,:,:) + fac * coeffs1(j) * fn(i+j-mx,:,:)
  enddo

 enddo

END SUBROUTINE DDX
