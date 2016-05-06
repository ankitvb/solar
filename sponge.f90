SUBROUTINE SPONGE(q, q_rhs)

 USE INIT

 implicit none

 real*8, dimension(mx,my,mr,5) :: q, q_rhs, q_ref 
 real*8, dimension(mr)         :: sigma_r
 
! Sponge parameters
 real*8 a_sponge, s_exp
 real*8 vel_ref, rho_ref, p_ref
  
 integer k_sponge_length, k_sponge_start_lower, k_sponge_start_upper
 integer i,j,k
 integer, dimension(mr) :: gk

 if(nr.eq.1) return

 a_sponge      = 0.1D0                        ! Sponge coefficient
 s_exp         = 3.0D0                        ! Power of sponge polynomial

 vel_ref = 0D0
 rho_ref = 0D0 
 p_ref   = 0D0

 k_sponge_length = 30                         ! Number of boundary points to be covered with sponge 
 k_sponge_start_lower = k_sponge_length 
 k_sponge_start_upper = nr-k_sponge_length+1
 
 do i=1,mx
  do j=1,my
   do k=1,mr
    gk(k) = rrank*mr+k

    q_ref(i,j,k,1) = rho_ref
    q_ref(i,j,k,4) = vel_ref
    q_ref(i,j,k,5) = p_ref

    if(gk(k).ge.k_sponge_start_upper)then
     sigma_r(k) = a_sponge*( REAL(gk(k)-k_sponge_start_upper)/REAL(k_sponge_length) )**s_exp
    elseif(gk(k).le.k_sponge_start_lower)then
     sigma_r(k) = a_sponge*( REAL(k_sponge_start_lower-gk(k))/REAL(k_sponge_length) )**s_exp
    else
     sigma_r(k) = 0D0
    endif

    q_rhs(i,j,k,1) = q_rhs(i,j,k,1) - sigma_r(k) * (q(i,j,k,1)-q_ref(i,j,k,1))
    q_rhs(i,j,k,4) = q_rhs(i,j,k,4) - sigma_r(k) * (q(i,j,k,4)-q_ref(i,j,k,4))
    q_rhs(i,j,k,5) = q_rhs(i,j,k,5) - sigma_r(k) * (q(i,j,k,5)-q_ref(i,j,k,5))

   enddo
  enddo
 enddo

 return

END SUBROUTINE SPONGE
