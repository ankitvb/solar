Module all_modules

! --------------------------------------------------------------------------
! MPI Version of the Spherical Acoustic Sun Simulator.
! Copyright 2006, Shravan Hanasoge
                                                                                                                                                        
! Hansen Experimental Physics Laboratory
! 455 Via Palou way, Stanford
! CA 94305, USA
! Email: shravan@stanford.edu
! --------------------------------------------------------------------------
!
! More Subroutines.
!

use initialize
implicit none

Contains
!==================================================================================

function norm2(matrix)
  
   implicit none
   integer i,j,ierr
   real*8 matrix(nx,dim2(rank))
   real*8 norm2, sum

   sum = 0.0  
   norm2  = 0.0
   do j =1,dim2(rank)
    do i =1,nx

      sum = sum + matrix(i,j)**2.0
    end do     
   end do     

   call MPI_REDUCE(sum, norm2, 1, MPI_DOUBLE_PRECISION, &
				MPI_SUM, 0, MPI_COMM_WORLD, ierr)
   

   norm2 = (norm2/(DBLE(nx)*DBLE(ny)))**0.5

end function norm2

!================================================================================
      SUBROUTINE writefits()
      
      ! Generates 2D images of velocities (vr, vtheta and vphi) at the surface 
      ! Also computes doppler maps of the front and back. It also includes important
      ! SOI keywords.

      integer blocksize,bitpix,naxes(2),unit1,unit2,unit3,unit4,unit5
      integer i,status1,status2,status3,status4,status5,k,ierr,disp,tag
      integer j,group,fpixel,nelements, stat(MPI_STATUS_SIZE),naxes_dopp(2)
      real*8 prod1, prod2
      character*80 filename1, filename2, filename3, filename4, filename5
      character T_LAST*27, numer*5
      logical simple,extend
      real*8, dimension(:,:,:), allocatable :: temp, temp_dopp
 
    if (rank ==0) then

      T_LAST = date_time()

      call convert_to_string(time,numer,5)
      filename2 = directory//'vx_at_'//numer//'_deltat.fits'
      filename3 = directory//'vy_at_'//numer//'_deltat.fits'
      filename4 = directory//'vz_at_'//numer//'_deltat.fits'

      blocksize=1


      status2=0
      status3=0
      status4=0

      call ftgiou(unit2,status2)
      call ftgiou(unit3,status3)
      call ftgiou(unit4,status4)

      blocksize=1

      call ftinit(unit2,filename2,blocksize,status2)
      call ftinit(unit3,filename3,blocksize,status3)
      call ftinit(unit4,filename4,blocksize,status4)

      simple=.true.
      bitpix=-32
      naxes(1)=nx
      naxes(2)=ny
!      naxes_dopp(1) = nx/2
 !     naxes_dopp(2) = ny
      extend=.false.
                                                                                                                                 
      call ftphpr(unit2,simple,bitpix,2,naxes,0,1,extend,status2)
      call ftphpr(unit3,simple,bitpix,2,naxes,0,1,extend,status3)
      call ftphpr(unit4,simple,bitpix,2,naxes,0,1,extend,status4)



      ! END OF THE KEYWORD SECTION  

         allocate(temp(nx,ny,3),temp_dopp(nx,ny,2))


         temp(:, 1:dim2(0), 1) = a(:,:,o_Rad,2)
         temp(:, 1:dim2(0), 2) = a(:,:,o_Rad,3)
         temp(:, 1:dim2(0), 3) = a(:,:,o_Rad,4)
         disp = dim2(0) + 1

 !      if (ALSO_THERMO_VARIABLES) then
        filename1 = directory//'bz_at_'//numer//'_deltat.fits'
        filename5 = directory//'p_at_'//numer//'_deltat.fits'
        status1 = 0
        status5 = 0
        call ftgiou(unit1,status1)
        call ftgiou(unit5,status5)
        call ftinit(unit1,filename1,blocksize,status1)
        call ftinit(unit5,filename5,blocksize,status5)
        call ftphpr(unit1,simple,bitpix,2,naxes,0,1,extend,status1)
        call ftphpr(unit5,simple,bitpix,2,naxes,0,1,extend,status5)
        temp_dopp(:, 1:dim2(0), 1) = a(:,:,o_Rad,8)
        temp_dopp(:, 1:dim2(0), 2) = a(:,:,o_Rad,5)
!       endif        


      endif

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      do k=1,numtasks-1
          if (rank == k) then
  	   nelements = nx*dim2(rank)

!           if (ALSO_THERMO_VARIABLES) then
 	    tag = 4
            call MPI_SEND(a(:,:,o_Rad,8), nelements, MPI_DOUBLE_PRECISION, &
			0, tag, MPI_COMM_WORLD, ierr)

	    tag = 5
            call MPI_SEND(a(:,:,o_Rad,5), nelements, MPI_DOUBLE_PRECISION, &
			0, tag, MPI_COMM_WORLD, ierr)
       !    endif

	   tag = 0
           call MPI_SEND(a(:,:,o_Rad,2), nelements, MPI_DOUBLE_PRECISION, &
			0, tag, MPI_COMM_WORLD, ierr)

	   tag = 1
           call MPI_SEND(a(:,:,o_Rad,3), nelements, MPI_DOUBLE_PRECISION, &
			0, tag, MPI_COMM_WORLD, ierr)

	   tag = 3
           call MPI_SEND(a(:,:,o_Rad,4), nelements, MPI_DOUBLE_PRECISION, &
			0, tag, MPI_COMM_WORLD, ierr)
          endif
	  if (rank ==0) then
            
           nelements = nx*dim2(k)
!           if (ALSO_THERMO_VARIABLES) then
  	    tag = 4
	    call MPI_RECV(temp_dopp(1,disp,1), nelements, MPI_DOUBLE_PRECISION,&
		  k, tag, MPI_COMM_WORLD, stat, ierr)

	    tag = 5
  	    call MPI_RECV(temp_dopp(1,disp,2), nelements, MPI_DOUBLE_PRECISION,&
			k, tag, MPI_COMM_WORLD, stat, ierr)
 !          endif

 	   tag = 0
	   call MPI_RECV(temp(1,disp,1), nelements, MPI_DOUBLE_PRECISION,&
		  k, tag, MPI_COMM_WORLD, stat, ierr)

	   tag = 1
  	   call MPI_RECV(temp(1,disp,2), nelements, MPI_DOUBLE_PRECISION,&
			k, tag, MPI_COMM_WORLD, stat, ierr)

	   tag = 3
 	   call MPI_RECV(temp(1,disp,3), nelements, MPI_DOUBLE_PRECISION,&
			k, tag, MPI_COMM_WORLD, stat, ierr)

 	   disp = disp + dim2(k)

	 endif
	 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     enddo


     if (rank == 0) then

      temp = temp * dimc

      group=1
      fpixel=1
      nelements = nx*ny
      call ftpprd(unit2,group,fpixel,nelements,temp(:,:,1),status2)
      call ftpprd(unit3,group,fpixel,nelements,temp(:,:,2),status3)
      call ftpprd(unit4,group,fpixel,nelements,temp(:,:,3),status4)

      call ftclos(unit2, status2)
      call ftfiou(unit2, status2)

      call ftclos(unit3, status3)
      call ftfiou(unit3, status3)

      call ftclos(unit4, status4)
      call ftfiou(unit4, status4)

      if (status2 .gt. 0)call printerror(status2)

      if (status3 .gt. 0)call printerror(status3)

      if (status4 .gt. 0)call printerror(status4)
    
!      if (ALSO_THERMO_VARIABLES) then
        
       temp_dopp(:,:,1)  = temp_dopp(:,:,1)*dimrho
       temp_dopp(:,:,2)  = temp_dopp(:,:,2)*dimrho*dimc**2.0

       call ftpprd(unit1,group,fpixel,nelements,temp_dopp(:,:,1),status3)
       call ftpprd(unit5,group,fpixel,nelements,temp_dopp(:,:,2),status4)

       call ftclos(unit1, status1)
       call ftfiou(unit1, status1)

       call ftclos(unit5, status5)
       call ftfiou(unit5, status5)

       if (status1 .gt. 0)call printerror(status1)

       if (status5 .gt. 0)call printerror(status5)

 !     endif
 
      deallocate(temp, temp_dopp)

     endif

     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   end SUBROUTINE writefits
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

   SUBROUTINE readfits(filename,read,dim3)

      implicit none
      integer status,unit,readwrite,blocksize,naxes(3), tag, k, dim3 
      integer group,firstpix, ierr, stat(MPI_STATUS_SIZE)
      integer disp, temptyp, sendtyp
      integer nelements
      real*8 nullval,read(nx,dim2(rank),dim3)
      real*8, dimension(:,:,:), allocatable :: temp
      logical anynull
      character*(*) filename
      integer (KIND=MPI_ADDRESS_KIND) ext, spacing

      
      if (rank == 0) then
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

       if (status .gt. 0)call printerror(status)

       read = temp(:,1:dim2(0),:)  

       call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION, ext, ierr)
       spacing = nx*ny*ext
       disp = dim2(0) + 1
      endif

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)       
      do k=1,numtasks-1
                     
        if (rank ==0) then
         print *,k
 	 call MPI_TYPE_VECTOR(nx*dim2(k),1,1,MPI_DOUBLE_PRECISION,temptyp,ierr)
	 call MPI_TYPE_CREATE_HVECTOR(dim3,1,spacing,temptyp,sendtyp,ierr)
	 call MPI_TYPE_COMMIT(sendtyp, ierr)
	 tag = 0
	 call MPI_SEND(temp(1,disp,1), 1, sendtyp, &
			k, tag, MPI_COMM_WORLD, ierr)

	 disp = disp + dim2(k)
        endif
        if (rank  == k) then
 	 print*,k
  	 nelements = nx * dim2(k) * dim3
	 tag = 0
	 call MPI_RECV(read, nelements, MPI_DOUBLE_PRECISION,&
			0, tag, MPI_COMM_WORLD, stat, ierr)
        endif

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
 
      enddo

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if (rank == 0) deallocate(temp)

    end SUBROUTINE readfits

!================================================================================


function norm(matrix)

   implicit none
   integer i, j, k, ierr
   real*8 matrix(nx,dim2(rank),nz)
   real*8 norm, sum
  
   norm  = 0.0
   sum = 0.0

   do k = 1,nz
    do j =1,dim2(rank)
     do i =1,nx
       sum = sum + matrix(i,j,k)**2.0
     end do     
    end do     
   end do 
 
   call MPI_REDUCE(sum, norm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr) 

   norm = (norm/(DBLE(nx)*DBLE(ny)*DBLE(nz)))**0.5d0

   call MPI_BARRIER(MPI_COMM_WORLD, ierr)
end function norm

!================================================================================

  SUBROUTINE convert_to_string(number,string,length_string)

  implicit none
  integer i,length_string,number,n(1:length_string),number_temp
  character*length_string string
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


  number_temp = number
  do i=length_string,1,-1
  
    n(length_string-i+1) = floor(number_temp*10.0**(-(i-1.0)))
    number_temp = number_temp - n(length_string-i+1)*10**(i-1)
    string(length_string-i+1:length_string-i+1) = charc(n(length_string-i+1)+1)

  enddo

  end SUBROUTINE convert_to_string

!================================================================================


function norm_1D(matrix)

   use initialize   
   implicit none
   integer i, j, k, ierr
   real*8 matrix(nz)
   real*8 norm_1d, sum

  
   norm_1d  = 0.0
   sum = 0.0
   do k = 1,nz
     sum = sum + matrix(k)**2.0
   end do 

   norm_1d = (sum/DBLE(nz))**0.5


end function norm_1d

!================================================================================

     function date_time()
  
      implicit none 
      character*27 date_time
      integer i, j, one_day, day, month, hour, minutes
      parameter(one_day = 1440)
      character*2 monthchar, daychar, hourchar, minutechar, secondchar
      character*4 yearchar

      ! Assuming the timestep (in seconds) divides 60 seconds.

      yearchar = '2007'
      monthchar = '06'
      secondchar = '00'
      day = time/(steps*one_day) ! Number of days of computation
      hour = (time/steps - day*1440)/60 ! Remainder hours
      minutes =  (time/steps - day*1440 - hour*60) ! Remainder number of minutes
      day = day + 1 ! Starting day is first day of the month

      call convert_to_string(day, daychar, 2)
      call convert_to_string(hour, hourchar, 2)
      call convert_to_string(minutes, minutechar, 2)

!      date_time ='2015.05.11_'//t1//':'//t2//':'01':00.0_UT'
    
      date_time = yearchar//'.'//monthchar//'.'//daychar//'_'//hourchar//':'//minutechar//':'//secondchar//':00.0_UT'
     end function date_time

!================================================================================
     SUBROUTINE writefits_local(filename, dump_array, dim3)

      implicit none
      integer blocksize,bitpix,naxes(3),unit1,dim3
      integer status1,group,fpixel, ierr
      integer temptype, sendtyp
      integer*8 nelements
      real*8 dump_array(nx,dim2(rank),dim3)
      character*(*) filename
      logical simple,extend

      status1 = 0
      call ftgiou(unit1,status1)
      blocksize=1
      call ftinit(unit1,filename,blocksize,status1)
      simple=.true.
      bitpix=-64
      naxes(1)=nx
      naxes(2)=dim2(rank)
      naxes(3)=dim3
      nelements=naxes(1)*naxes(2)*naxes(3)
      extend=.false.
      group=1
      fpixel=1
                                                                                                                                                      
      call ftphpr(unit1,simple,bitpix,3,naxes,0,1,extend,status1)
      call ftpprd(unit1,group,fpixel,nelements,dump_array,status1)
      call ftclos(unit1, status1)
      call ftfiou(unit1, status1)

      if (status1 .gt. 0) call printerror(status1)

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

     end SUBROUTINE writefits_local
!===============================================================================

     SUBROUTINE writefits_3D(filename,dump_array,dim3)

      implicit none 
      integer blocksize,bitpix,naxes(3),unit1,k,disp, tag, dim3
      integer status1,group,fpixel,ierr,stat(MPI_STATUS_SIZE)
      integer temptype, sendtyp, nelements, ext
      real*8 dump_array(nx,dim2(rank),dim3)
      real*8, allocatable, dimension(:,:,:) :: temp
      character*(*) filename
      logical simple,extend
      integer (KIND = MPI_ADDRESS_KIND) space
 
      if (rank == 0) then
	 allocate(temp(nx,ny,dim3))
         temp(:, 1:dim2(0), :) = dump_array
         disp = dim2(0) + 1
      endif
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      do k=1,numtasks-1

       if (rank == k) then
	 tag = 0
         nelements = nx*dim3*dim2(rank)
         call MPI_SEND(dump_array, nelements, MPI_DOUBLE_PRECISION, &
			0, tag, MPI_COMM_WORLD, ierr)
       endif
       if (rank == 0) then

	  call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION,ext,ierr)
	  space = nx*ny*ext

	  tag = 0
	  call MPI_TYPE_VECTOR(nx*dim2(k),1,1,MPI_DOUBLE_PRECISION,temptype,ierr)
	  call MPI_TYPE_CREATE_HVECTOR(dim3,1,space,temptype,sendtyp,ierr)
	  call MPI_TYPE_COMMIT(sendtyp, ierr)
	  call MPI_RECV(temp(1,disp,1), 1, sendtyp, &
			k, tag, MPI_COMM_WORLD, stat, ierr)

	  disp = disp + dim2(k)
        endif
 	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       enddo
      
       if (rank == 0) then
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
         call ftpkyj(unit1,'TIME',time,'Also, init = time+1', status1)
         call ftpprd(unit1,group,fpixel,nelements,temp,status1)
         call ftclos(unit1, status1)
         call ftfiou(unit1, status1)

         deallocate(temp)
         if (status1 .gt. 0) call printerror(status1)

       endif

       call MPI_BARRIER(MPI_COMM_WORLD,ierr)

     end SUBROUTINE writefits_3D
!================================================================================
function sum2(matrix, dima, dimb, option)
  
   implicit none
   integer i,j,ierr, option, dima, dimb
   real*8 matrix(dima, dimb)
   real*8 sum2, summer

   summer = 0.0  
   sum2  = 0.0
   if (option == 0) then
    do j =1,dimb
     do i =1,dima
       summer = summer + matrix(i,j)
     end do     
    end do     
   else
    do j =1,dimb
     do i =1,dima
       summer = summer + abs(matrix(i,j))
     end do     
    end do     
   endif

   call MPI_REDUCE(summer, sum2, 1, MPI_DOUBLE_PRECISION, &
				MPI_SUM, 0, MPI_COMM_WORLD, ierr)
   

   sum2 = (sum2/(DBLE(dima)*DBLE(dimb)))

end function sum2

!================================================================================
function sum2_normal(matrix, dima, dimb, option)
  
   implicit none
   integer i,j,ierr, option, dima, dimb
   real*8 matrix(dima, dimb)
   real*8 sum2_normal, summer

   summer = 0.0  
   sum2_normal  = 0.0
   if (option == 0) then
    do j =1,dimb
     do i =1,dima
       summer = summer + matrix(i,j)
     end do     
    end do     
   else
    do j =1,dimb
     do i =1,dima
       summer = summer + abs(matrix(i,j))
     end do     
    end do     
   endif

   sum2_normal = (summer/(DBLE(dima)*DBLE(dimb)))

end function sum2_normal


!===============================================================================

  SUBROUTINE SETUP_COORDINATES_METRICS_ALL(nx, ny, noverlap, xi, eta, theta, phi, interpoints,&
                     metric_C, metric_D, metric_delta, sqrt_metric_delta, XY, CD)
                     
  implicit none
  integer i,j,nx, ny, noverlap
  real*8 dxi, deta
  real*8, dimension(nx, noverlap) :: interpoints
  real*8 theta(nx,nx,6), phi(nx,nx,6), tanxi(nx), taneta(nx), pi, xi(nx), eta(nx)
  real*8, dimension(nx) :: metric_C, metric_D
  real*8, dimension(nx,nx) :: metric_delta, sqrt_metric_delta, XY, CD
 

  pi = acos(-1.0d0)

  do i=1,nx
   xi(i) = -pi*0.25 + (i-1.0)* pi * 0.5/(nx - 1.0)
   eta(i) = -pi*0.25 + (i-1.0)* pi * 0.5/(nx-1.0)
  enddo
  
  tanxi = tan(xi)
  taneta = tan(eta) 

  dxi = xi(2) - xi(1)
  deta = dxi
 
  metric_C = (1.0 + tanxi**2.)**0.5
  metric_D = (1.0 + taneta**2.)**0.5

  do j=1,ny
   do i=1,nx
    metric_delta(i,j) = (1. + tanxi(i)**2. + taneta(j)**2.)
    sqrt_metric_delta(i,j) = (1. + tanxi(i)**2. + taneta(j)**2.)**(0.5)

    XY(i,j) = tanxi(i) * taneta(j)
    CD(i,j) = metric_C(i) * metric_D(j)
   enddo
  enddo

  ! INTERPOLATION POINTS:

  do i = 1,noverlap
   interpoints(:, i) = atan(taneta/tan(xi(nx) + i*dxi)) 
  enddo

  ! DETERMINING SPHERICAL COORDINATES

  !REGION 1

  do j=1,nx
   do i=1,nx
    phi(i,j,1) = xi(i) !+ pi
    theta(i,j,1) = atan(1./(taneta(j) * cos(phi(i,j,1))))
    if (theta(i,j,1) < 0) theta(i,j,1) = theta(i,j,1) + pi
!    if (phi(i,j,1) < 0) phi(i,j,1) = phi(i,j,1) + 2.*pi
   enddo
  enddo	

  !REGION 2

  ! phi(:,:,2) = phi(:,:,1) + pi*0.5
  ! theta(:,:,2) = theta(:,:,1)


  do j=1,nx
   do i=1,nx
    phi(i,j,2) = xi(i) + pi*0.5
    theta(i,j,2) = atan(1./(taneta(j) * sin(phi(i,j,2))))
    if (theta(i,j,2) < 0) theta(i,j,2) = theta(i,j,2) + pi
    if (phi(i,j,2) < 0) phi(i,j,2) = phi(i,j,2) + 2.*pi

   enddo
  enddo	

  !REGION 3

  ! phi(:,:,3) = phi(:,:,2) + pi*0.5
  ! theta(:,:,3) = theta(:,:,2)

  do j=1,nx
   do i=1,nx
    phi(i,j,3) = xi(i) + pi
    theta(i,j,3) = atan(-1./(taneta(j) * cos(phi(i,j,3))))
    if (theta(i,j,3) < 0) theta(i,j,3) = theta(i,j,3) + pi
    if (phi(i,j,3) < 0) phi(i,j,3) = phi(i,j,3) + 2.*pi

   enddo
  enddo	

  !REGION 4

  ! phi(:,:,4) = phi(:,:,3) + pi*0.5
  ! theta(:,:,4) = theta(:,:,3)

  do j=1,nx
   do i=1,nx
    phi(i,j,4) = xi(i) + pi*1.5
    theta(i,j,4) = atan(-1./(taneta(j) * sin(phi(i,j,4))))
    if (theta(i,j,4) < 0) theta(i,j,4) = theta(i,j,4) + pi
    if (phi(i,j,4) < 0) phi(i,j,4) = phi(i,j,4) + 2.*pi
   enddo
  enddo	

  !REGION 5

  do j=1,nx
   do i=1,nx
    theta(i,j,5) = atan((tanxi(i)**2. + taneta(j)**2.)**0.5)
    phi(i,j,5) = atan2(tanxi(i),-taneta(j)) !+ pi
    if (phi(i,j,5) < 0) phi(i,j,5) = phi(i,j,5) + 2.*pi
   enddo
  enddo	

  !REGION 6

  do j=1,nx
   do i=1,nx
    theta(i,j,6) = pi - atan((tanxi(i)**2. + taneta(j)**2.)**0.5)
    phi(i,j,6) = atan2(tanxi(i),taneta(j)) !+ pi
    if (phi(i,j,6) < 0) phi(i,j,6) = phi(i,j,6) + 2.*pi
   enddo
  enddo	

END SUBROUTINE SETUP_COORDINATES_METRICS_ALL

!--------------------------

SUBROUTINE SETUP_CONNECTIVITY(interfaces)
 implicit none
 integer interfaces(6,6)

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

!--------------------------

SUBROUTINE INTERPOLATE_ALONG_Y(input_func, func_int, interface_in, interface_out, &
			interpolant, eta, nx, nord, noverlap)

 implicit none
 integer nx, nord, i, j, k, noverlap, lower, upper, interface_in, interface_out
 real*8 func_int(noverlap,nx), interpolant(-nord:nord, nx, noverlap)
 real*8 input_func(noverlap,nx), eta(nx)


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

!--------------------------

SUBROUTINE INTERPOLATE_ALONG_X(input_func, func_int, interface_in, interface_out, &
			interpolant, eta, nx, nord, noverlap)

 implicit none
 integer nx, nord, i, j, k, noverlap, lower, upper, interface_in, interface_out
 real*8 func_int(nx,noverlap), interpolant(-nord:nord, nx, noverlap)
 real*8 input_func(nx, noverlap), eta(nx)

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

!--------------------------

SUBROUTINE SETUP_LAGRANGE_INTERPOLANT(interpolant, interpoints, eta, nx, noverlap, nord)
 
  implicit none
  integer i, j, k, ipoint, nord, noverlap, nx, iter, lower, upper
  real*8 interpolant(-nord:nord, nx, noverlap), interpoints(nx, noverlap), eta(nx)

!   g^{\rm interp}(j) = \sum_{k=-n_{ord}}^{n_{ord}} I(k) f(j+k)
!   I(k) is the interpolation function
!   f(k) is the set of points we interpolate from
!   g^{interp}(k) is the set of points we interpolate to

  interpolant = 1.0
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

! ---------------------------


end module all_modules

