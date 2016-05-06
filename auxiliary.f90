MODULE AUXILIARY

USE INIT

CONTAINS

!***********************************************************************
! Subroutine: Compute_Max
!
! Computes the maximum value of the elements in F, which can be
! distributed over several nodes.
!
! In: F: array with dimensions ny, nz, mx
!     mx, ny, nz: array dimensions
!
! Out: F_max: max value
!***********************************************************************
SUBROUTINE COMPUTE_MAX(f, f_max)

 Include 'mpif.h'

!  Procedure parameters
 Real*8  f(mx,my,mr), f_max

!  Locals
 Integer ierr
 Real*8  f_localmax

!  Compute the maximum value on this processor
 f_localmax = maxval(f)

!  Find the maximum value over all processors
 Call MPI_ALLREDUCE(f_localmax, f_max, 1, MPI_REAL8, MPI_MAX,&
                    MPI_COMM_WORLD, ierr)


END SUBROUTINE COMPUTE_MAX

!***********************************************************************
! Subroutine: Compute_Min
!
! Computes the minimum value of the elements in F, which can be
! distributed over several nodes.
!
! In: F: array with dimensions ny, nz, mx
!     mx, ny, nz: array dimensions
!
! Out: F_min: min value
!***********************************************************************
SUBROUTINE COMPUTE_MIN(f, f_min)

 Include 'mpif.h'

! Procedure parameters
 Real*8  f(mx,my,mr), f_min

! Locals
 Integer ierr
 Real*8  f_localmin

! Compute the minimum value on this processor
 f_localmin = minval(f)

! Find the minimum value over all processors
 Call MPI_ALLREDUCE(f_localmin, f_min, 1, MPI_REAL8, MPI_MIN,&
                    MPI_COMM_WORLD, ierr)


END SUBROUTINE COMPUTE_MIN
!************************************************************************


function norm2(matrix)
   implicit none
   include "mpif.h"   

   integer i,j,ierr
   real*8 matrix(mx,my)
   real*8 norm2, summ

   summ = 0.0  
   norm2  = 0.0
   do j =1,my
    do i =1,mx
      summ = summ + matrix(i,j)**2.0
    end do     
   end do     

   call MPI_REDUCE(summ, norm2, 1, MPI_REAL8, &
				MPI_SUM, 0, MPI_COMM_WORLD, ierr)
   

   norm2 = (norm2/(DBLE(nx)*DBLE(ny)*6.0))**0.5

end function norm2
!==================================================================================



function norm(matrix)

   implicit none
   include "mpif.h"   

   integer i, j, k, ierr
   real*8 matrix(mx,my,mr)
   real*8 norm, summ
  
   norm  = 0.0
   summ = 0.0

   do k = 1,mr
    do j =1,my
     do i =1,mx
       summ = summ + matrix(i,j,k)**2.0
     end do     
    end do     
   end do 
 
   call MPI_REDUCE(summ, norm, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr) 

   norm = (norm/(DBLE(nx)*DBLE(ny)*DBLE(nr))*6.0)**0.5d0

   call MPI_BARRIER(MPI_COMM_WORLD, ierr)
end function norm

!================================================================================

END MODULE AUXILIARY



