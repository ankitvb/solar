MODULE OUTPUT

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
      call ftgiou(unit1,status1)
      blocksize=1
      call ftinit(unit1,filename,blocksize,status1)
      simple=.true.
      bitpix=-64
      naxes(1)=nx
      naxes(2)=ny
      naxes(3)=nz
      nelements=naxes(1)*naxes(2)*naxes(3)
      extend=.false.
      group=1
      fpixel=1
                                                                                                                                                      
      call ftphpr(unit1,simple,bitpix,3,naxes,0,1,extend,status1)
      call ftpprd(unit1,group,fpixel,nelements,dump_array,status1)
      call ftclos(unit1, status1)
      call ftfiou(unit1, status1)

      if (status1 .gt. 0) call printerror(status1)


     end SUBROUTINE writefits

!================================================================================


SUBROUTINE writefits4d(filename, dump_array, nx, ny, nz,nt)

      implicit none
      integer blocksize,bitpix,naxes(4),unit1,nx,ny,nz,nt
      integer status1,group,fpixel, ierr
      integer temptype, sendtyp
      integer*8 nelements
      real*8 dump_array(nx,ny,nz,nt)
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
      naxes(3)=nz
      naxes(4) = nt
      nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)
      extend=.false.
      group=1
      fpixel=1
                                                                                                                                                      
      call ftphpr(unit1,simple,bitpix,4,naxes,0,1,extend,status1)
      call ftpprd(unit1,group,fpixel,nelements,dump_array,status1)
      call ftclos(unit1, status1)
      call ftfiou(unit1, status1)

      if (status1 .gt. 0) call printerror(status1)


     end SUBROUTINE writefits4d

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
       call ftgiou(unit,status)
       readwrite=0
       print *,'Now reading the file: '//filename
       call ftopen(unit,filename,readwrite,blocksize,status)

       naxes(1) = dim1
       naxes(2) = dim2
       naxes(3) = dim3
       nelements=naxes(1)*naxes(2)*naxes(3)
       group=1
       firstpix=1
       nullval=-999

       call ftgpvd(unit,group,firstpix,nelements,nullval, &
       &            read,anynull,status)


       call ftclos(unit, status)
       call ftfiou(unit, status)

       if (status .gt. 0)call printerror(status)


    end SUBROUTINE readfits

!================================================================================
    SUBROUTINE CONVERT_VIS_TO_FITS

    USE INIT

    character*128 gridfname, visfname, fitsfname, dirname
    character*128 vis_path, rank_str, cue_str, l_str
    character*128 index_str(5)

    integer start_xi, end_xi, start_eta, end_eta, start_r, end_r, zoneid 

    real*8, dimension(mx,my,5)   :: q_in
    real*8, dimension(nx,ny,6,5) :: q_out

    integer i,j,k,l,n,ii,jj

    vis_path = 'Visual'

    do cue_vis=1,800

     print *, 'Writing fits file for cue_vis #',cue_vis

     WRITE(index_str(1), '(I128)') FLOOR(REAL(cue_vis/100))
     WRITE(index_str(2), '(I128)') FLOOR(REAL(MOD(cue_vis,100)/10))
     WRITE(index_str(3), '(I128)') MOD(cue_vis,10)

     dirname =  TRIM(ADJUSTL(index_str(1))) // &
                TRIM(ADJUSTL(index_str(2))) // &
                TRIM(ADJUSTL(index_str(3)))

     do n=0,nprocs-1
      WRITE(rank_str, '(I128)') n 
      visfname  = TRIM(ADJUSTL(vis_path)) // '/' // &
                  TRIM(ADJUSTL(dirname))  // '/' // &
                  'vis_'// TRIM(ADJUSTL(rank_str)) 
      gridfname = TRIM(ADJUSTL(vis_path)) // '/' //&
                  'grid/' // 'grid_' // TRIM(ADJUSTL(rank_str)) 

      open(unit=8, file=gridfname, status='unknown', form='formatted', action='read')
      read(8,6) start_xi, end_xi, start_eta, end_eta, start_r, end_r, zoneid
      close(8)

      if((o_rad.ge.start_r).and.(o_rad.le.end_r))then
       open(unit=10, file=visfname, status='unknown', form='unformatted', action='read')
       read(10) (((q_in(i,j,l),i=1,mx),j=1,my),l=1,5)
       close(10)
       do i=1,mx
        do j=1,my
         do l=1,5
          ii = start_xi + i 
          jj = start_eta + j 
          q_out(ii,jj,zoneid,l) = q_in(i,j,l)
         enddo
        enddo
       enddo

      endif
     enddo

     do l=4,4
      WRITE(l_str, '(I128)') l
       fitsfname = TRIM(ADJUSTL(vis_path)) // '/' // &
                  'q' // TRIM(ADJUSTL(l_str)) // '_' // &
                   TRIM(ADJUSTL(dirname)) // '.fits'
      call writefits(fitsfname, q_out(1,1,1,l), nx,ny,6)
     enddo
 
    enddo

6   FORMAT(1x, i5,i5,i5,i5,i5,i5,i5)

    END SUBROUTINE CONVERT_VIS_TO_FITS

!================================================================================
    SUBROUTINE CONVERT_VIS_TO_BINARY

    USE INIT

    character*128 gridfname, visfname, binfname, dirname
    character*128 vis_path, rank_str, cue_str, l_str
    character*128 index_str(5)
 
    integer start_xi, end_xi, start_eta, end_eta, start_r, end_r, zoneid 

    real*8, dimension(mx,my,5)   :: q_in
    real*8, dimension(nx,ny,6,5) :: q_out

    integer i,j,k,l,n,ii,jj

    vis_path = 'Visual'

    do cue_vis=1,1000

     print *, 'Writing binary file for cue_vis #',cue_vis

     WRITE(index_str(1), '(I128)') FLOOR(REAL(cue_vis/100))
     WRITE(index_str(2), '(I128)') FLOOR(REAL(MOD(cue_vis,100)/10))
     WRITE(index_str(3), '(I128)') MOD(cue_vis,10)

     dirname =  TRIM(ADJUSTL(index_str(1))) // &
                TRIM(ADJUSTL(index_str(2))) // &
                TRIM(ADJUSTL(index_str(3)))

     do n=0,nprocs-1
      WRITE(rank_str, '(I128)') n 
      visfname  = TRIM(ADJUSTL(vis_path)) // '/' // &
                  TRIM(ADJUSTL(dirname))  // '/' // &
                  'vis_'// TRIM(ADJUSTL(rank_str)) 
      gridfname = TRIM(ADJUSTL(vis_path)) // '/' //&
                  'grid/' //'grid_' // TRIM(ADJUSTL(rank_str)) 

      open(unit=8, file=gridfname, status='unknown', form='formatted', action='read')
      read(8,6) start_xi, end_xi, start_eta, end_eta, start_r, end_r, zoneid
      close(8)

      if((o_rad.ge.start_r).and.(o_rad.le.end_r))then
       open(unit=10, file=visfname, status='unknown', form='unformatted', action='read')
       read(10) (((q_in(i,j,l),i=1,mx),j=1,my),l=1,5)
       close(10)
       do i=1,mx
        do j=1,my
         do l=1,5
          ii = start_xi + i 
          jj = start_eta + j 
          q_out(ii,jj,zoneid,l) = q_in(i,j,l)
         enddo
        enddo
       enddo

      endif
     enddo

    binfname = TRIM(ADJUSTL(vis_path)) // '/' // &
               'q_' // TRIM(ADJUSTL(dirname)) // '_2d.dat'
   open(unit=10, file=binfname, status='unknown', form='unformatted', action='write')
   write(10)((((q_out(i,j,k,l),i=1,nx),j=1,ny),k=1,6),l=1,5) 
   close(10)
 
   enddo

6   FORMAT(1x, i5,i5,i5,i5,i5,i5,i5)

    END SUBROUTINE CONVERT_VIS_TO_BINARY

!================================================================================


      SUBROUTINE printerror(status)

      ! See cookbook.f on FITSIO for details

      integer status
      character errtext*30,errmessage*80

      if (status .le. 0)return

      call ftgerr(status,errtext)
      print *,'FITSIO Error Status =',status,': ',errtext

      call ftgmsg(errmessage)
      do while (errmessage .ne. ' ')
          print *,errmessage
          call ftgmsg(errmessage)
      end do

      end SUBROUTINE printerror

!================================================================================
      SUBROUTINE deletefile(filename,status)

      integer status,unit,blocksize
      character*(*) filename

      if (status .gt. 0)return
      call ftgiou(unit,status)
      call ftopen(unit,filename,1,blocksize,status)
      

      if (status .eq. 0)then
          call ftdelt(unit,status)
      endif
      if (status .eq. 103)then
          status=0
          call ftcmsg
      endif 
      if ((status .NE. 0) .and. (status .ne. 103)) then 

          status=0
          call ftcmsg
          call ftdelt(unit,status)
      end if

      call ftfiou(unit, status)

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
    SUBROUTINE WRITE_RANDOM_FIELD(filename,f,nx,ny,nt)
    implicit none
    integer nx,ny,nt
    real*8, dimension(nx,ny,6,nt) :: f

    character*(*) filename
    integer i,j,l,m

    open(unit=10, file=filename//'.1', status='unknown', form='unformatted')    
    write(unit=10)(((f(i,j,1,m),m=1,nt),j=1,ny),i=1,nx) 
    close(10)

    open(unit=10, file=filename//'.2', status='unknown', form='unformatted')    
    write(unit=10)(((f(i,j,2,m),m=1,nt),j=1,ny),i=1,nx) 
    close(10)

    open(unit=10, file=filename//'.3', status='unknown', form='unformatted')    
    write(unit=10)(((f(i,j,3,m),m=1,nt),j=1,ny),i=1,nx) 
    close(10)

    open(unit=10, file=filename//'.4', status='unknown', form='unformatted')    
    write(unit=10)(((f(i,j,4,m),m=1,nt),j=1,ny),i=1,nx) 
    close(10)

    open(unit=10, file=filename//'.5', status='unknown', form='unformatted')    
    write(unit=10)(((f(i,j,5,m),m=1,nt),j=1,ny),i=1,nx) 
    close(10)

    open(unit=10, file=filename//'.6', status='unknown', form='unformatted')    
    write(unit=10)(((f(i,j,6,m),m=1,nt),j=1,ny),i=1,nx) 
    close(10)

    END SUBROUTINE WRITE_RANDOM_FIELD
!================================================================================
    SUBROUTINE READ_RANDOM_FIELD(filename,f,nx,ny,nt)
    implicit none
 
    integer nx,ny,nt
    real*8, dimension(nx,ny,6,nt) :: f

    character*(*) filename
    integer i,j,l,m


    open(unit=10, file=filename, status='unknown', form='unformatted')
    read(unit=10)((((f(i,j,l,m),m=1,nt),l=1,6),j=1,ny),i=1,nx)
    close(10)

6   Format(1x, 6(ES13.5))

    END SUBROUTINE READ_RANDOM_FIELD


!================================================================================


END MODULE OUTPUT 
