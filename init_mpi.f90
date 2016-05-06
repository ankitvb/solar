subroutine INITIALIZE_MPI
 
  implicit none

  include 'header'
  include 'mpif.h'

  integer :: i,j,p
  integer :: ierr

  !intarr - number of processors in each directions
  !coords - returned coordinates of current processor
  integer, dimension(3) :: intarr, coords
  logical :: reorder
  logical, dimension(3) :: tmplog
  integer :: color, key
  integer starts(6)

!  write(*,*)"MPI Setup"
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, ierr)
! Computing total number of processors
  nprocs=0
  do j=1,6
   nprocs = nprocs + px*py*pr
   np(j)  = px*py*pr
  enddo

! Check to make sure px,py,pz add up to correct number
  if (mpi_size.ne.nprocs) stop ' Number of processors allocated not consistent with header declaration'

! Assign zone ids to each process
  if((rank.ge.0).and.(rank.le.(np(1)-1))) zoneid=1
  do j=2,6
    if((rank.gt.((j-1)*np(j-1)-1)).and.(rank.le.(j*np(j)-1))) zoneid=j
  enddo
  
  do i=1,6
    starts(i) = px*py*pr*(i-1)
  enddo

!  print *, rank, zoneid 

! Split MPI_COMM_WORLD into 6 zones
  call MPI_COMM_SPLIT(MPI_COMM_WORLD, zoneid, 0, MPI_COMM_ZONE, ierr)
  
  call MPI_COMM_RANK(MPI_COMM_ZONE, zonerank, ierr)
!  print *, rank, zonerank, zoneid 

!Set the dimensions of the cartesian arrangement, made to be same order as breakup order
  intarr(1)=px
  intarr(2)=py
  intarr(3)=pr

!Set periodicity of points
  tmplog(1) = .false.
  tmplog(2) = .false.
  tmplog(3) = .false.
  reorder   = .false.

! Create a cartesian topology within each zone
  call MPI_CART_CREATE(MPI_COMM_ZONE, 3, intarr, tmplog, reorder, MPI_XYR_COMM,ierr)
  call MPI_CART_GET(MPI_XYR_COMM, 3, intarr, tmplog, coords, ierr)

! x-direction
  call MPI_CART_SHIFT(MPI_XYR_COMM, 0, -1, source_left_x, dest_left_x, ierr)
  call MPI_CART_SHIFT(MPI_XYR_COMM, 0, 1, source_right_x, dest_right_x, ierr)

! y-direction
  call MPI_CART_SHIFT(MPI_XYR_COMM, 1, -1, source_left_y, dest_left_y, ierr)
  call MPI_CART_SHIFT(MPI_XYR_COMM, 1, 1, source_right_y, dest_right_y, ierr)

! z-direction
  call MPI_CART_SHIFT(MPI_XYR_COMM, 2, -1, source_left_r, dest_left_r, ierr)
  call MPI_CART_SHIFT(MPI_XYR_COMM, 2, 1, source_right_r, dest_right_r, ierr)

  color = coords(1)*py*pr+coords(2)*pr
  key = coords(3)
  call MPI_COMM_SPLIT(MPI_XYR_COMM, color, key, MPI_R_COMM, ierr)
  call MPI_COMM_RANK(MPI_R_COMM, rrank, ierr)

  color = coords(1)*pr*py+coords(3)
  key = coords(2)
  call MPI_COMM_SPLIT(MPI_XYR_COMM, color, key, MPI_Y_COMM, ierr)
  call MPI_COMM_RANK(MPI_Y_COMM, yrank, ierr)
  
  color = coords(2)*pr+coords(3)
  key = coords(1)
  call MPI_COMM_SPLIT(MPI_XYR_COMM, color, key, MPI_X_COMM, ierr)
  call MPI_COMM_RANK(MPI_X_COMM, xrank, ierr)

!  write(*,"(i4 i4 i4 i4 i4 i4 i4)"),rank,coords(1),coords(2),coords(3),xrank,yrank,rrank

! Don't bother with intercomms if topology can't support it
! Should be used only for testing
  if(px.ne.py) return

! Creating intercommunicators

  if(zoneid.eq.1)then
   call MPI_INTERCOMM_CREATE(MPI_COMM_ZONE,0,MPI_COMM_WORLD,starts(2),1, &
                             MPI_INTERCOMM(1,2),ierr)
   call MPI_INTERCOMM_CREATE(MPI_COMM_ZONE,0,MPI_COMM_WORLD,starts(4),1, &
                             MPI_INTERCOMM(1,4),ierr)
   call MPI_INTERCOMM_CREATE(MPI_COMM_ZONE,0,MPI_COMM_WORLD,starts(5),1, &
                             MPI_INTERCOMM(1,5),ierr)
   call MPI_INTERCOMM_CREATE(MPI_COMM_ZONE,0,MPI_COMM_WORLD,starts(6),1, &
                             MPI_INTERCOMM(1,6),ierr)  
  endif
  if(zoneid.eq.2)then
   call MPI_INTERCOMM_CREATE(MPI_COMM_ZONE,0,MPI_COMM_WORLD,starts(1),1, &
                             MPI_INTERCOMM(1,2),ierr)
   call MPI_INTERCOMM_CREATE(MPI_COMM_ZONE,0,MPI_COMM_WORLD,starts(3),1, &
                             MPI_INTERCOMM(2,3),ierr) 
   call MPI_INTERCOMM_CREATE(MPI_COMM_ZONE,0,MPI_COMM_WORLD,starts(5),1, &
                             MPI_INTERCOMM(2,5),ierr)
   call MPI_INTERCOMM_CREATE(MPI_COMM_ZONE,0,MPI_COMM_WORLD,starts(6),1, &
                             MPI_INTERCOMM(2,6),ierr) 
  endif
  if(zoneid.eq.3)then
   call MPI_INTERCOMM_CREATE(MPI_COMM_ZONE,0,MPI_COMM_WORLD,starts(2),1, &
                             MPI_INTERCOMM(2,3),ierr)
   call MPI_INTERCOMM_CREATE(MPI_COMM_ZONE,0,MPI_COMM_WORLD,starts(4),1, &
                             MPI_INTERCOMM(3,4),ierr)
   call MPI_INTERCOMM_CREATE(MPI_COMM_ZONE,0,MPI_COMM_WORLD,starts(5),1, &
                             MPI_INTERCOMM(3,5),ierr)
   call MPI_INTERCOMM_CREATE(MPI_COMM_ZONE,0,MPI_COMM_WORLD,starts(6),1, &
                             MPI_INTERCOMM(3,6),ierr)
  endif
  if(zoneid.eq.4)then
   call MPI_INTERCOMM_CREATE(MPI_COMM_ZONE,0,MPI_COMM_WORLD,starts(1),1, &
                             MPI_INTERCOMM(1,4),ierr)
   call MPI_INTERCOMM_CREATE(MPI_COMM_ZONE,0,MPI_COMM_WORLD,starts(3),1, &
                             MPI_INTERCOMM(3,4),ierr)
   call MPI_INTERCOMM_CREATE(MPI_COMM_ZONE,0,MPI_COMM_WORLD,starts(5),1, &
                             MPI_INTERCOMM(4,5),ierr)
   call MPI_INTERCOMM_CREATE(MPI_COMM_ZONE,0,MPI_COMM_WORLD,starts(6),1, &
                             MPI_INTERCOMM(4,6),ierr)
  endif
  if(zoneid.eq.5)then
   call MPI_INTERCOMM_CREATE(MPI_COMM_ZONE,0,MPI_COMM_WORLD,starts(1),1, &
                             MPI_INTERCOMM(1,5),ierr)
   call MPI_INTERCOMM_CREATE(MPI_COMM_ZONE,0,MPI_COMM_WORLD,starts(2),1, &
                             MPI_INTERCOMM(2,5),ierr)
   call MPI_INTERCOMM_CREATE(MPI_COMM_ZONE,0,MPI_COMM_WORLD,starts(3),1, &
                             MPI_INTERCOMM(3,5),ierr)
   call MPI_INTERCOMM_CREATE(MPI_COMM_ZONE,0,MPI_COMM_WORLD,starts(4),1, &
                             MPI_INTERCOMM(4,5),ierr)
  endif
  if(zoneid.eq.6)then
   call MPI_INTERCOMM_CREATE(MPI_COMM_ZONE,0,MPI_COMM_WORLD,starts(1),1, &
                             MPI_INTERCOMM(1,6),ierr)
   call MPI_INTERCOMM_CREATE(MPI_COMM_ZONE,0,MPI_COMM_WORLD,starts(2),1, &
                             MPI_INTERCOMM(2,6),ierr)
   call MPI_INTERCOMM_CREATE(MPI_COMM_ZONE,0,MPI_COMM_WORLD,starts(3),1, &
                             MPI_INTERCOMM(3,6),ierr)
   call MPI_INTERCOMM_CREATE(MPI_COMM_ZONE,0,MPI_COMM_WORLD,starts(4),1, &
                             MPI_INTERCOMM(4,6),ierr)
  endif

  
! Setting communication partners for boundary nodes
  if(zoneid.eq.1)then
   if(xrank.eq.0)then
    coords(1)=px-1; coords(2)=yrank; coords(3)=rrank
    call MPI_CART_RANK(MPI_XYR_COMM, coords, icom_source_right_x, ierr)
    icom_dest_left_x = icom_source_right_x
   endif
   if(xrank.eq.(px-1))then
    coords(1)=0; coords(2)=yrank; coords(3)=rrank
    call MPI_CART_RANK(MPI_XYR_COMM, coords, icom_source_left_x, ierr)
    icom_dest_right_x = icom_source_left_x
   endif
   if(yrank.eq.0)then
    coords(1)=xrank; coords(2)=py-1; coords(3)=rrank
    call MPI_CART_RANK(MPI_XYR_COMM, coords, icom_source_right_y, ierr)
    icom_dest_left_y = icom_source_right_y
   endif
   if(yrank.eq.(py-1))then
    coords(1)=xrank; coords(2)=0; coords(3)=rrank
    call MPI_CART_RANK(MPI_XYR_COMM, coords, icom_source_left_y, ierr)
    icom_dest_right_y = icom_source_left_y
   endif
  endif

  if(zoneid.eq.2)then
   if(xrank.eq.0)then
    coords(1)=px-1; coords(2)=yrank; coords(3)=rrank
    call MPI_CART_RANK(MPI_XYR_COMM, coords, icom_source_right_x, ierr)
    icom_dest_left_x = icom_source_right_x
   endif
   if(xrank.eq.(px-1))then
    coords(1)=0; coords(2)=yrank; coords(3)=rrank
    call MPI_CART_RANK(MPI_XYR_COMM, coords, icom_source_left_x, ierr)
    icom_dest_right_x = icom_source_left_x
   endif
   if(yrank.eq.0)then
   coords(1)=px-1; coords(2)=px-1-xrank; coords(3)=rrank
    call MPI_CART_RANK(MPI_XYR_COMM, coords, icom_source_right_y, ierr)
    icom_dest_left_y = icom_source_right_y
   endif
   if(yrank.eq.(py-1))then
    coords(1)=px-1; coords(2)=xrank; coords(3)=rrank
    call MPI_CART_RANK(MPI_XYR_COMM, coords, icom_source_left_y, ierr)
    icom_dest_right_y = icom_source_left_y
   endif
  endif
 
  if(zoneid.eq.3)then
   if(xrank.eq.0)then
    coords(1)=px-1; coords(2)=yrank; coords(3)=rrank
    call MPI_CART_RANK(MPI_XYR_COMM, coords, icom_source_right_x, ierr)
    icom_dest_left_x = icom_source_right_x
   endif
   if(xrank.eq.(px-1))then
    coords(1)=0; coords(2)=yrank; coords(3)=rrank
    call MPI_CART_RANK(MPI_XYR_COMM, coords, icom_source_left_x, ierr)
    icom_dest_right_x = icom_source_left_x
   endif
   if(yrank.eq.0)then
    coords(1)=px-1-xrank; coords(2)=0; coords(3)=rrank
    call MPI_CART_RANK(MPI_XYR_COMM, coords, icom_source_right_y, ierr)
    icom_dest_left_y = icom_source_right_y
   endif
   if(yrank.eq.(py-1))then
    coords(1)=px-1-xrank; coords(2)=py-1; coords(3)=rrank
    call MPI_CART_RANK(MPI_XYR_COMM, coords, icom_source_left_y, ierr)
    icom_dest_right_y = icom_source_left_y
   endif 
  endif 
 
  if(zoneid.eq.4)then
   if(xrank.eq.0)then
    coords(1)=px-1; coords(2)=yrank; coords(3)=rrank
    call MPI_CART_RANK(MPI_XYR_COMM, coords, icom_source_right_x, ierr)
    icom_dest_left_x = icom_source_right_x
   endif
   if(xrank.eq.(px-1))then
    coords(1)=0; coords(2)=yrank; coords(3)=rrank
    call MPI_CART_RANK(MPI_XYR_COMM, coords, icom_source_left_x, ierr)
    icom_dest_right_x = icom_source_left_x
   endif
   if(yrank.eq.0)then
    coords(1)=0; coords(2)=xrank; coords(3)=rrank
    call MPI_CART_RANK(MPI_XYR_COMM, coords, icom_source_right_y, ierr)
    icom_dest_left_y = icom_source_right_y
   endif
   if(yrank.eq.(py-1))then
    coords(1)=0; coords(2)=px-1-xrank; coords(3)=rrank
    call MPI_CART_RANK(MPI_XYR_COMM, coords, icom_source_left_y, ierr)
    icom_dest_right_y = icom_source_left_y
   endif
  endif

  if(zoneid.eq.5)then
   if(xrank.eq.0)then
    coords(1)=py-1-yrank; coords(2)=py-1; coords(3)=rrank
    call MPI_CART_RANK(MPI_XYR_COMM, coords, icom_source_right_x, ierr)
    icom_dest_left_x = icom_source_right_x
   endif
   if(xrank.eq.(px-1))then
    coords(1)=yrank; coords(2)=py-1; coords(3)=rrank
    call MPI_CART_RANK(MPI_XYR_COMM, coords, icom_source_left_x, ierr)
    icom_dest_right_x = icom_source_left_x
   endif
   if(yrank.eq.0)then
    coords(1)=xrank; coords(2)=py-1; coords(3)=rrank
    call MPI_CART_RANK(MPI_XYR_COMM, coords, icom_source_right_y, ierr)
    icom_dest_left_y = icom_source_right_y
   endif
   if(yrank.eq.(py-1))then
    coords(1)=px-1-xrank; coords(2)=py-1; coords(3)=rrank
    call MPI_CART_RANK(MPI_XYR_COMM, coords, icom_source_left_y, ierr)
    icom_dest_right_y = icom_source_left_y
   endif
  endif
  
  if(zoneid.eq.6)then
   if(xrank.eq.0)then
    coords(1)=yrank; coords(2)=0; coords(3)=rrank
    call MPI_CART_RANK(MPI_XYR_COMM, coords, icom_source_right_x, ierr)
    icom_dest_left_x = icom_source_right_x
!    print *, xrank, yrank, icom_source_right_x,'rechts'
   endif
   if(xrank.eq.(px-1))then
    coords(1)=py-1-yrank; coords(2)=0; coords(3)=rrank
    call MPI_CART_RANK(MPI_XYR_COMM, coords, icom_source_left_x, ierr)
    icom_dest_right_x = icom_source_left_x
!    print *, xrank, yrank, icom_source_left_x,'links'
   endif
   if(yrank.eq.0)then
    coords(1)=px-1-xrank; coords(2)=0; coords(3)=rrank
    call MPI_CART_RANK(MPI_XYR_COMM, coords, icom_source_right_y, ierr)
    icom_dest_left_y = icom_source_right_y
!    print *, xrank, yrank, icom_source_right_y,'rechts'
   endif
   if(yrank.eq.(py-1))then
    coords(1)=xrank; coords(2)=0; coords(3)=rrank
    call MPI_CART_RANK(MPI_XYR_COMM, coords, icom_source_left_y, ierr)
    icom_dest_right_y = icom_source_left_y
!    print *, xrank, yrank, icom_source_left_y,'links'
   endif
  endif

  return
 
end subroutine Initialize_MPI
