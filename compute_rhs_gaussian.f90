SUBROUTINE COMPUTE_RHS(q, q_rhs)

 USE INIT
 USE GRAD 

 implicit none

 real*8, dimension(mx,my,mr,5) :: q, q_rhs 

 real*8, dimension(mx,my,mr,3) :: temp 
 
 q_rhs = 0D0                                                                          
                                                                                        
 call COMPUTE_GRADIENT(q(1,1,1,1), temp(1,1,1,1), temp(1,1,1,2), temp(1,1,1,3))                      
 call COMPUTE_DIVERGENCE(temp(1,1,1,1), temp(1,1,1,2), temp(1,1,1,3), q_rhs(1,1,1,2))                
                                                                                        
 q_rhs(:,:,:,1) = q(:,:,:,2)                         

END SUBROUTINE COMPUTE_RHS
