PROGRAM CUBE

 USE INIT
 USE DERIVATIVE
 USE INTERPOLATION
 USE METRIC
 USE TRANSFORM
 USE GRAD
 USE OUTPUT
 USE FILTER
 USE AUXILIARY

 implicit none

 include 'mpif.h'

! Globals
 real*8, dimension(mx,my,mr,5) :: q

 COMMON /RKWORK/ q

 real*8, dimension(mx,my,nt) :: S_t

 COMMON /FORCING/ S_t

 real*8 q1_min, q1_max, tt, rms(nt)
 
! Locals
 integer i, j, k, n, ii, jj
 integer ierr

! MPI timing
 real*8 start_time, step_time

 call INITIALIZE_MPI

 call SOLAR_DATA 

 call SETUP_COORDINATES_METRICS_ALL

 call SETUP_LAGRANGE_INTERPOLANT

 call SETUP_CONNECTIVITY

 call COMPUTE_TIMESTEP   

 call GENERATE_IC(q)

 call READ_FORCING_FUNC

 call WRITE_RADIAL_TRACE

 nstep = 6000  ! to be read from an input file
 if(rank.eq.0) print *,'nstep = ',nstep

!  Start timing
 start_time = MPI_WTIME()

 cue_restart=0
 cue_vis=0
 step=0
 step_old=0 
 
! Write initial restart file
! call WRITE_RESTART 
! cue_restart = cue_restart + 1

 call WRITE_GRID
 
! Evolve in time
 do step=0,nstep
   start_time = MPI_WTIME()

   call RK4 

   if(rank.eq.0) print *, step, deltat, time
   tt = norm(q(:,:,:,4))
   if (rank==0) print *,tt

!  Determine step timing
   step_time = MPI_WTIME()

   if(rank.eq.0) print *, 'mpitime = ',step_time-start_time  

   if (MOD(step,6).eq.0) then
    call WRITE_VISUAL
    call WRITE_RADIAL_TRACE
    cue_vis = cue_vis + 1
!    call WRITE_RESTART
!    cue_restart = cue_restart + 1
   endif
 enddo 

 call MPI_FINALIZE(ierr)

 stop

  
END PROGRAM CUBE

!--------------------------

