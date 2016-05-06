MODULE OUTPUT

USE INIT

CONTAINS

!-----------------------------------------------------------------------

SUBROUTINE writefits(filename, dump_array, nx, ny, nz)

      implicit none
      integer blocksize,bitpix,naxes(3),unit1,nx,ny,nz
      integer status1,group,fpixel, ierr
      integer temptype, sendtyp
      integer*8 nelements
      real*8 dump_array(nx,ny,nz)
      character*(*) filename
      logical simple,extend

      status1 = 0
      call ftgiou_(unit1,status1)
      blocksize=1
      call ftinit_(unit1,filename,blocksize,status1)
      simple=.true.
      bitpix=-64
      naxes(1)=nx
      naxes(2)=ny
      naxes(3)=nz
      nelements=naxes(1)*naxes(2)*naxes(3)
      extend=.false.
      group=1
      fpixel=1
                                                                                                                                                      
      call ftphpr_(unit1,simple,bitpix,3,naxes,0,1,extend,status1)
      call ftpprd_(unit1,group,fpixel,nelements,dump_array,status1)
      call ftclos_(unit1, status1)
      call ftfiou_(unit1, status1)

      if (status1 .gt. 0) call printerror(status1)


     end SUBROUTINE writefits

!================================================================================


   SUBROUTINE readfits(filename,read,dim1, dim2, dim3)

      implicit none
      integer status,unit,readwrite,blocksize,naxes(3), dim1, dim2, dim3 
      integer group,firstpix
      integer disp, temptyp, sendtyp
      integer nelements
      real*8 nullval,read(dim1, dim2,dim3)
      logical anynull
      character*(*) filename

      
       status=0
       call ftgiou_(unit,status)
       readwrite=0
       print *,'Now reading the file: '//filename
       call ftopen_(unit,filename,readwrite,blocksize,status)

       naxes(1) = dim1
       naxes(2) = dim2
       naxes(3) = dim3
       nelements=naxes(1)*naxes(2)*naxes(3)
       group=1
       firstpix=1
       nullval=-999

       call ftgpvd_(unit,group,firstpix,nelements,nullval, &
       &            read,anynull,status)


       call ftclos_(unit, status)
       call ftfiou_(unit, status)

       if (status .gt. 0)call printerror(status)


    end SUBROUTINE readfits

!================================================================================
    SUBROUTINE WRITE_RESTART

    implicit none

    include 'mpif.h'

!  Globals
    real*8, dimension(mx,my,mr,5) :: q

    COMMON  /RKWORK/ q

!  Locals
    character*128 index, zoneid_str
    character*128 filename, restart_path, dirname
    character*128 index_str(5)
    Integer       I, J, K, L, M, N
    Integer       nx_read, ny_read, nr_read, px_read
    Integer       status(MPI_STATUS_SIZE), ierr
    Integer       gsizes(4), start_indices(4) 
    Integer       lsizes(4), ndims
    Integer       count1, count2, count3
    Integer       fh, filetype, local_array_size
    Integer(kind=MPI_OFFSET_KIND) disp
    Integer intSize, realSize

    gsizes(1) = nx 
    gsizes(2) = ny
    gsizes(3) = nr
    gsizes(4) = 5

    count1 = 4
    count2 = 1
    count3 = 1

    restart_path = 'Restart'
    WRITE(index_str(1), '(I128)') FLOOR(REAL(cue_restart/100))
    WRITE(index_str(2), '(I128)') FLOOR(REAL(MOD(cue_restart,100)/10))
    WRITE(index_str(3), '(I128)') MOD(cue_restart,10)
    dirname = TRIM(ADJUSTL(index_str(1))) //&
              TRIM(ADJUSTL(index_str(2))) //&
              TRIM(ADJUSTL(index_str(3)))

! Generating file name
    WRITE(zoneid_str, '(I128)') zoneid
    filename = TRIM(ADJUSTL(restart_path)) // '/' // &
               'restart.'// TRIM(ADJUSTL(dirname)) &
               //'.'//TRIM(ADJUSTL(zoneid_str))

    Call MPI_FILE_OPEN(MPI_COMM_ZONE, filename,&
                       MPI_MODE_CREATE + MPI_MODE_WRONLY,&
                       MPI_INFO_NULL, fh, ierr)
    if(ierr.NE.MPI_SUCCESS) Stop 'Error opening file for write'

    Call MPI_TYPE_SIZE(MPI_INTEGER, intSize, ierr)
    Call MPI_TYPE_SIZE(MPI_REAL8, realSize, ierr)

! Only processor 0 writes the header
    if(rank .EQ. 0)then
      Call MPI_FILE_WRITE(fh, gsizes, count1, MPI_INTEGER,&
                          MPI_STATUS_IGNORE, ierr)
      Call MPI_FILE_WRITE(fh, time, count2, MPI_REAL8,&
                          MPI_STATUS_IGNORE, ierr)
      Call MPI_FILE_WRITE(fh, step, count3, MPI_INTEGER,&
                          MPI_STATUS_IGNORE, ierr)
    endif

    ndims = 4

    lsizes(1) = mx
    lsizes(2) = my
    lsizes(3) = mr
    lsizes(4) = 5

    start_indices(1) = xrank * mx
    start_indices(2) = yrank * my
    start_indices(3) = rrank * mr
    start_indices(4) = 0

    local_array_size = mx * my * mr * 5

    disp = 5*intSize + realSize

    Call MPI_TYPE_CREATE_SUBARRAY(ndims,gsizes,lsizes,start_indices,&
                                  MPI_ORDER_FORTRAN, MPI_REAL8,&
                                  filetype, ierr)

    Call MPI_TYPE_COMMIT(filetype, ierr)

    Call MPI_BARRIER(MPI_COMM_ZONE, ierr)

    Call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL8, filetype,&
                          "native", MPI_INFO_NULL, ierr)
    If(ierr.NE.MPI_SUCCESS)Stop 'Error setting file view'

    Call MPI_FILE_WRITE_ALL(fh, Q, local_array_size,&
                            MPI_REAL8, status, ierr)
    if(ierr.NE.MPI_SUCCESS)&
    stop 'Error during collective write of solution file.'

    Call MPI_FILE_CLOSE(fh, ierr)
    if(ierr.NE.MPI_SUCCESS)Stop 'Error closing file'

    Return

    END SUBROUTINE WRITE_RESTART 
!================================================================================
    SUBROUTINE WRITE_VISUAL

    character*128 filename, dirname, vis_path, rank_str, cue_str
    character*128 index_str(5)
    character*20 vis_path_C

    real*8, dimension(mx,my,mr,5) :: q
    real*8, dimension(mx,my,5)    :: q_out

    COMMON /RKWORK/ q

    integer i,j,k,l

    vis_path = 'Visual'
    WRITE(rank_str, '(I128)') rank
    WRITE(cue_str, '(I128)') cue_vis

    vis_path_C = TRIM(ADJUSTL(vis_path)) // CHAR(0)
    call Make_Dir(cue_vis, vis_path_C)

    WRITE(index_str(1), '(I128)') FLOOR(REAL(cue_vis/100))
    WRITE(index_str(2), '(I128)') FLOOR(REAL(MOD(cue_vis,100)/10))
    WRITE(index_str(3), '(I128)') MOD(cue_vis,10)

    dirname = TRIM(ADJUSTL(index_str(1))) // &
              TRIM(ADJUSTL(index_str(2))) // &
              TRIM(ADJUSTL(index_str(3)))

! Generating file name
    filename = TRIM(ADJUSTL(vis_path)) // '/' // &
               TRIM(ADJUSTL(dirname))  // '/' // &
               'vis_'// TRIM(ADJUSTL(rank_str))

    do k=1,mr
     if((rrank*mr+k).eq.o_rad)then
      do i=1,mx
       do j=1,my
        do l=1,5
         q_out(i,j,l) = q(i,j,k,l)
        enddo
       enddo
      enddo
      open(unit=14, file=filename, status='unknown', form='unformatted', action='write')
      write(14) (((q_out(i,j,l),i=1,mx),j=1,my),l=1,5)
      close(14)
     endif
    enddo

    END SUBROUTINE WRITE_VISUAL
!================================================================================ 
    SUBROUTINE WRITE_GRID

    implicit none

    character*128 filename, vis_path, rank_str
    integer start_xi, end_xi, start_eta, end_eta, start_r, end_r

    vis_path = 'Visual/grid'
    WRITE(rank_str, '(I128)') rank

    filename = TRIM(ADJUSTL(vis_path)) // '/' // 'grid_' // TRIM(ADJUSTL(rank_str))

    start_xi  = xrank*mx
    end_xi    = start_xi + mx - 1
    start_eta = yrank*my
    end_eta   = start_eta + my - 1
    start_r   = rrank*mr
    end_r     = start_r + mr - 1

    open(unit=12, file=filename, status='unknown', form='formatted', action='write')     
    write(12,6) start_xi, end_xi, start_eta, end_eta, start_r, end_r, zoneid 
    close(12)

6   FORMAT(1x, i5,i5,i5,i5,i5,i5,i5)

   END SUBROUTINE WRITE_GRID

!================================================================================

    SUBROUTINE READ_FORCING_FUNC
 
    implicit none
 
    real*8, dimension(:,:,:), allocatable :: f
    real*8, dimension(mx,my,nt) :: S_t

    COMMON /FORCING/ S_t

    character*128 filename, zoneid_str
    integer i,j,n,ii,jj
    integer err
 
    WRITE(zoneid_str, '(I128)') zoneid
   
    filename='forcing/forcing_cubesph' // '.' // TRIM(ADJUSTL(zoneid_str))
   
    allocate(f(nx,ny,nt),STAT=err)
    if(err.ne.0) Stop '*** Array allocation failed!'

! Each processor reads in the input file for its zone
    open(unit=10, file=filename, status='unknown', form='unformatted', action='read')
    read(unit=10)(((f(i,j,n),n=1,nt),j=1,ny),i=1,nx)
    close(10)
 
! Each processor allocates its slice
    do i=1,mx
     do j=1,my
      do n=1,nt
       ii=xrank*mx+i
       jj=yrank*my+j
       S_t(i,j,n) = f(ii,jj,n)
      enddo
    enddo
   enddo
 
   deallocate(f)

   END SUBROUTINE READ_FORCING_FUNC
!================================================================================

  SUBROUTINE WRITE_RADIAL_TRACE

  implicit none

  include 'mpif.h'

  character*128 filename, zone_str
  real*8, dimension(mx,my,mr,5) :: q
  integer i, k, l
  integer ierr

  COMMON /RKWORK/ q

  WRITE(zone_str, '(I128)') zoneid

  filename = 'trace_' // TRIM(ADJUSTL(zone_str))

  if((xrank.eq.px/2).and.(yrank.eq.py/2))then
   do i=0,pr-1
    if(rrank.eq.i)then
      open(unit=zoneid+10, file=filename, status='unknown', form='formatted', action='write',&
           position='append')
      write(zoneid+10,7) ((q(1,1,k,l), k=1,mr), l=1,5)
      close(zoneid+10)
     endif
     call MPI_BARRIER(MPI_R_COMM, ierr)
   enddo
  endif

7 Format(1x, 5(ES13.5))

  END SUBROUTINE WRITE_RADIAL_TRACE

!================================================================================


!================================================================================
    SUBROUTINE printerror(status)

    ! See cookbook.f on FITSIO for details

    integer status
    character errtext*30,errmessage*80

    if (status .le. 0)return

    call ftgerr_(status,errtext)
    print *,'FITSIO Error Status =',status,': ',errtext

    call ftgmsg_(errmessage)
    do while (errmessage .ne. ' ')
        print *,errmessage
        call ftgmsg_(errmessage)
    end do

    end SUBROUTINE printerror

!================================================================================
      SUBROUTINE deletefile(filename,status)

      integer status,unit,blocksize
      character*(*) filename

      if (status .gt. 0)return
      call ftgiou_(unit,status)
      call ftopen_(unit,filename,1,blocksize,status)
      

      if (status .eq. 0)then
          call ftdelt_(unit,status)
      endif
      if (status .eq. 103)then
          status=0
          call ftcmsg_
      endif 
      if ((status .NE. 0) .and. (status .ne. 103)) then 

          status=0
          call ftcmsg_
          call ftdelt_(unit,status)
      end if

      call ftfiou_(unit, status)

      end SUBROUTINE deletefile

!================================================================================



 SUBROUTINE convert_to_string(numbe,sting,length_string)

  implicit none
  integer i,length_string,numbe,n(1:length_string),number_temp
  character*(length_string) sting
  character*1 charc(10)

  charc(1)  = '0'
  charc(2)  = '1'
  charc(3)  = '2'
  charc(4)  = '3'
  charc(5)  = '4'
  charc(6)  = '5'
  charc(7)  = '6'
  charc(8)  = '7'
  charc(9)  = '8'
  charc(10) = '9'


  number_temp = numbe
  do i=length_string,1,-1

    n(length_string-i+1) = floor(number_temp*10.0**(-(i-1.0)))
    number_temp = number_temp - n(length_string-i+1)*10**(i-1)
    sting(length_string-i+1:length_string-i+1) = charc(n(length_string-i+1)+1)

  enddo

  end SUBROUTINE convert_to_string

!================================================================================

END MODULE OUTPUT 
