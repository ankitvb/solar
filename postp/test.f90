PROGRAM TEST

 use bspline
 implicit none

 integer i, j, k, nx, ny, nz, kord,mm,nx2,ny2,nt
 parameter(nx = 768, ny = nx, nz = nx, kord = 6)
 parameter(nx2 = 800, ny2 = nx2, nt = 643)
 real*8  t1, t2, xx, yy, zz
 real*8, dimension(nx2) :: xnew
 real*8, dimension(nx) :: x, y, z, bcoef1, ydata
 real*8, dimension(nx, nx) :: a, bcoef, c, val
 real*8, dimension(nx+kord) :: xknot, yknot, zknot 
 real*8, dimension(nx, nx, nt) :: bb
 real*8, dimension(nx2, nx2, nt) :: a2


 

 call dbsnak(nx,x,kord,xknot)
 call dbsnak(nx,x(nx:-1:1),kord,yknot)

! yknot = xknot

! call dbs3in(nx,x,ny,y,nz,z,a,nx,ny,kord,kord,kord,    &
!         & xknot,yknot,zknot,bcoef)

 call dbsint(nx,x,ydata,kord,xknot,bcoef1)

 call CPU_TIME(t1)
 do mm=1,nt
! do k=1,nz-15
 call dbs2in(nx,x,ny,y,bb(:,:,mm),nx,kord,kord,    &
         & xknot,yknot,bcoef)
 print *,mm

  do j=1,nx2!205,602
   do i=1,nx2!205,602
!    xx = x(i) + 0.012
 !   yy = y(j) + 0.0033
!    zz = z(k) + 0.0009312

!    val(1,1) = dbsval(xx,kord,xknot,nx,bcoef1)
!    print *,i,j,xnew(i),xnew(j)

    a2(i,j,mm) = dbs2vl(xnew(i),xnew(j),kord,kord,xknot,yknot,nx,ny,bcoef)
   enddo
  enddo
! enddo
 enddo
 call CPU_TIME(t2)
 print*,t2 - t1
! print *,val,0.6334**2.0 + 0.3**2.0 + 0.7**2.0

 call writefits('/tmp20/shravan/spag4/aspline.fits',a2,nx2,ny2,nt)
END PROGRAM TEST

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

