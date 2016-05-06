  real*8, dimension(mx,nst,mr)  :: int_north, int_south, int_east_t, int_west_t
  real*8, dimension(nst,my,mr)  :: int_east, int_west, int_north_t, int_south_t
  real*8, dimension(mx,my,mr)   :: test_func, dr_func, filt_func
  real*8, dimension(mx,my,mr)   :: b, b_theta, b_phi, b_xi, b_eta, b_r, div


!------------------------------------------------------------------------------------
! Checking interpolation indices
  if(zoneid.eq.1)then
   if((yrank.eq.0).and.(xrank.eq.0))then
!   do k=1,noverlap                                                                    
    do j=1,nord
      jj = yrank*my+j
      do i=-nord,-j
!       print *, jj, j, i, nord+i+j                                                    
      enddo
      do i=-j+1,nord
        print *, jj, i, interpolant(i,jj,1)!j, i+j                         
      enddo                                                                            
    enddo
!   enddo
    endif        
  endif
                                                                                        
! Right (upper) boundary
 if(zoneid.eq.1)then                                                                   
  if((yrank.eq.1).and.(xrank.eq.0))then                                                
!  do k=1,noverlap
   do j=my-nord+1,my
     jj = yrank*my+j                                                                   
     do i=-nord,my-j
      print *, jj, i, interpolant(i,jj,1)!j, j+i                                       
     enddo                                                                             
     do i=my-j+1,nord                                                                  
!      print *, jj, j, i, i+j-my                                                       
     enddo                                                                             
   enddo
!  enddo
  endif                                                                                
 endif

!------------------------------------------------------------------------------------
! Initial field for interpolation  and derivative tests
  do i=1,mx                                                                             
   do j=1,my
    do k=1,mr                                                                           
     jj = yrank*my+j
     test_func(i,j,k) = sin(theta(i,j))                                                 
    enddo
   enddo
  enddo

! Y-interpolation
  do i=1,mx
   do j=1,my         
    do k=1,mr
     jj = yrank*my+j
     test_func(i,j,k) = sin(theta(i,j)) 
    enddo
   enddo
  enddo                  
   
  call INTERPOLATE_ALONG_Y(test_func(2:noverlap+1,:,:), func_int)

  if(zoneid.eq.1)then
   if(xrank.eq.0)then
    do i=0,py-1
     if(yrank.eq.i)then
      do j=1,my
       jj = yrank*my+j
       print *, jj, test_func(5,j,1), func_int(4,j,1)
      enddo
     endif
     call MPI_BARRIER(MPI_Y_COMM,ierr)
    enddo    
   endif
  endif

  call MPI_FINALIZE(ierr)

! X-interpolation
  call INTERPOLATE_ALONG_X(test_func(:,2:noverlap+1,:), func_int)                       

  if(zoneid.eq.1)then                                                                   
   if(yrank.eq.0)then                                                                   
    do j=0,px-1
     if(xrank.eq.j)then  
      do i=1,mx                                                                         
       ii = xrank*mx+i                                                                  
       print *, ii, test_func(i,4,1), func_int(i,3,1)                                   
      enddo                                                                             
     endif 
     call MPI_BARRIER(MPI_X_COMM,ierr)                                                  
    enddo                                                                               
   endif
  endif

!------------------------------------------------------------------------------------
! Y-derivative
  call INTERPOLATE_ALONG_X(test_func(:,2:noverlap+1,:), int_south)                      
  call INTERPOLATE_ALONG_X(test_func(:,my-1:my-noverlap:-1,:), int_north)               
  call INTERPOLATE_ALONG_Y(test_func(2:noverlap+1,:,:), int_west)                       
  call INTERPOLATE_ALONG_Y(test_func(mx-1:mx-noverlap:-1,:,:), int_east)                
   
  do i=1,nst
   do j=1,my
    int_east_t(j,i,:) = int_east(i,j,:)                                                 
    int_west_t(j,i,:) = int_west(i,j,:)                                                 
   enddo                                                                                
  enddo
                                                                                        
  call DDY(test_func, dy_func, int_south, int_north, int_west_t, int_east_t)            
   
  if(zoneid.eq.6)then                                                                   
   if(yrank.eq.0)then
    do j=0,px-1                                                                         
     if(xrank.eq.j)then                                                                 
      do i=1,mx
       ii = xrank*mx+i                                                                  
!       print *, ii, dy_func(i,my-2,1), -metric_C(i)*tan(eta(my-2))*(metric_D(my-2)**2D0)&
!                *(metric_delta(i,my-2)**(-1.5D0))   ! Zones 1-4
        print *, ii, dy_func(i,my,1), tan(eta(my))*(metric_D(my)**2D0)&
                *( metric_delta(i,my)**(-0.5)*(metric_delta(i,my)-1.0)**(-0.5) -&       
                   metric_delta(i,my)**(-1.5)*(metric_delta(i,my)-1.0)**(0.5) ) ! Zones 5-6         
      enddo                                                                             
     endif
     call MPI_BARRIER(MPI_X_COMM,ierr)                                                  
    enddo                                                                               
   endif                                                                                
  endif                                             

! X-derivative
  call INTERPOLATE_ALONG_X(test_func(:,2:noverlap+1,:), int_south)                      
  call INTERPOLATE_ALONG_X(test_func(:,my-1:my-noverlap:-1,:), int_north)               
  call INTERPOLATE_ALONG_Y(test_func(2:noverlap+1,:,:), int_west)                       
  call INTERPOLATE_ALONG_Y(test_func(mx-1:mx-noverlap:-1,:,:), int_east)                
   
  do i=1,nst
   do j=1,my
    int_north_t(i,j,:) = int_north(j,i,:)                                               
    int_south_t(i,j,:) = int_south(j,i,:)                                               
   enddo                                                                                
  enddo
                                                                                        
  call DDX(test_func, dx_func, int_south_t, int_north_t, int_west, int_east)            
   
  if(zoneid.eq.6)then                                                                   
   if(xrank.eq.0)then
    do i=0,px-1                                                                         
     if(yrank.eq.i)then                                                                 
      do j=1,my
       jj = yrank*my+j                                                                  
!       print *,ii,dx_func(i,my/2,1),metric_C(i)*tan(xi(i))*(metric_delta(i,my/2)**(-0.5))&                                                                                      !                *(1.0 - (metric_C(i)**2.0)/metric_delta(i,my/2))
       print *, jj, dx_func(1,j,1), tan(xi(1))*(metric_C(1)**2D0)&
               *( metric_delta(1,j)**(-0.5)*(metric_delta(1,j)-1.0)**(-0.5) -&          
                  metric_delta(1,j)**(-1.5)*(metric_delta(1,j)-1.0)**(0.5) )            
      enddo                                                                             
     endif
     call MPI_BARRIER(MPI_Y_COMM,ierr)                                                  
    enddo                                                                               
   endif                                                                                
  endif                                        

!------------------------------------------------------------------------------------
! Divergence
  do i=1,mx
   do j=1,my
    do k=1,mr
     b_theta(i,j,k) = sin(theta(i,j))
    enddo
   enddo
  enddo
  b_phi = 0.0; b_r = 0.0

  call CONVERT_THETA_PHI_TO_XI_ETA(b_theta, b_phi, b_xi, b_eta)

  call COMPUTE_DIVERGENCE(b_xi, b_eta, b_r, div) 

  if(zoneid.eq.3)then
   if(yrank.eq.0)then
    do j=0,px-1
     if(xrank.eq.j)then
      do i=1,mx
       ii = xrank*mx+i
!       print *, ii, div(i,1,1)
      enddo
     endif
     call MPI_BARRIER(MPI_X_COMM,ierr)
    enddo
   endif
  endif

  if(zoneid.eq.5)then
   if(xrank.eq.px-1)then
    do i=0,py-1
     if(yrank.eq.i)then
      do j=1,my
       jj = yrank*my+j
!       print *, jj, div(mx,j,1)
      enddo
     endif
     call MPI_BARRIER(MPI_Y_COMM,ierr)
    enddo
   endif
  endif


  call MPI_FINALIZE(ierr)


!------------------------------------------------------------------------------------
! Filters

! Initial field for filter test
  do i=1,mx
   do j=1,my
    do k=1,mr
     jj = yrank*my+j
     test_func(i,j,k) = sin(theta(i,j))
    enddo
   enddo
  enddo

  call INTERPOLATE_ALONG_X(test_func(:,2:noverlap+1,:), int_south)
  call INTERPOLATE_ALONG_X(test_func(:,my-1:my-noverlap:-1,:), int_north)
  call INTERPOLATE_ALONG_Y(test_func(2:noverlap+1,:,:), int_west)
  call INTERPOLATE_ALONG_Y(test_func(mx-1:mx-noverlap:-1,:,:), int_east)

  do j=1,noverlap
   do i=1,mx
    int_east_t(i,j,:)  = int_east(j,i,:)
    int_west_t(i,j,:)  = int_west(j,i,:)
    int_north_t(j,i,:) = int_north(i,j,:)
    int_south_t(j,i,:) = int_south(i,j,:)
   enddo
  enddo

!  call FILTER_XI(test_func, filt_func, int_south_t, int_north_t, int_west, int_east)
  call FILTER_ETA(test_func, filt_func, int_south, int_north, int_west_t, int_east_t)

  if(zoneid.eq.6)then
   if(xrank.eq.0)then
    do i=0,px-1
     if(xrank.eq.i)then
      do j=1,my
       jj = yrank*my+j
       print *, jj, test_func(1,j,1), filt_func(1,j,1)
      enddo
     endif
     call MPI_BARRIER(MPI_Y_COMM,ierr)
    enddo
   endif
  endif

  call MPI_FINALIZE(ierr)

!------------------------------------------------------------------------------------
! radial derivative
  do i=1,mx
   do j=1,my
    do k=1,mr
     test_func(i,j,k) = sin(r(k))
    enddo 
   enddo
  enddo

  call DDR(test_func, dr_func)

  if(zoneid.eq.1)then
   if((xrank.eq.0).and.(yrank.eq.0))then
    do i=0,pr-1
     if(rrank.eq.i)then
      do k=1,mr
       print *, mr*rrank+k, dr_func(1,1,k), cos(r(k))
      enddo
     endif
     call MPI_BARRIER(MPI_R_COMM, ierr)
    enddo
   endif
  endif
  
  call MPI_FINALIZE(ierr)

!------------------------------------------------------------------------------------
! base parameters
  if(zoneid.eq.1)then
   if((xrank.eq.0).and.(yrank.eq.0))then
    do i=0,pr-1
     if(rrank.eq.i)then
      do k=1,mr
       print *, c_speed(1,1,k)*dimc, g0(1,1,k)*dimc**2.0/diml, gamma(1,1,k)
      enddo
     endif
     call MPI_BARRIER(MPI_R_COMM, ierr)
    enddo
   endif
  endif

  call MPI_FINALIZE(ierr)
!------------------------------------------------------------------------------------
! radial metric

  if(zoneid.eq.1) print *, k, gr(k), drdchi(k)
