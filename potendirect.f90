module potendirect
 use octreetype
 implicit none 
 contains 
 subroutine get_accel(x,a,m,np)
  integer, intent(in) :: np
  real, intent(in) :: x(3,np)
  real, intent(out) :: a(3,np)
  real, intent(in) :: m(np)
  !real :: massratio
  integer :: i, j
  real :: dx(3), r, r2
  real :: h 

  h = 50.0

  a = 0.0
 	do i=1, np
 		do j=1, np
 			if (j/= i) then ! Don't want to calculate the force from a particle on itself
             dx = x(:,i) - x(:,j)
             r2 = dot_product(dx,dx)
             r  = sqrt(r2)
             a(:,i) = a(:,i) - m(j)*dx*(1/(r2 + h**2)**1.5)
 			endif 
 		enddo
 	enddo

end subroutine get_accel

subroutine get_poten (x,poten,m,np,particlesindex1,particlesindex2)
  ! the index's of particles (BODIES)
  integer, intent(in) :: particlesindex1(10), particlesindex2(10), np
  real, intent(in) :: x(3,np)
  ! technially this is potential 
  real, intent(inout) :: poten(np)
  real, intent(in) :: m(np)
  real :: r,r2,dx(3), h
  integer :: i, j
  integer :: indexi,indexj

  h = 50.0
  ! in this case we don't want to reset a
  do i=1, 10
   indexi = particlesindex1(i)
   do j=1, 10
    indexj = particlesindex2(j)
    if (j/=i .AND. i /= 0 .AND. j /= 0) then
     dx = x(:,indexi) - x(:,indexj)
     r2 = dot_product(dx,dx)
     r  = sqrt(r2)
     poten(indexi) = poten(indexi) - m(indexj)*(1/(sqrt(r + h**2)))
    endif
  enddo
 enddo 
end subroutine get_poten

subroutine get_poten_node(x,poten,m,np,particlesindex,node)
 integer, intent(in) :: particlesindex(10),np
 real, intent(in) :: x(3,np)
 real, intent(inout) :: poten(np)   
 real, intent(in) :: m(np)
 type(octreenode) :: node
 real :: nodecm(3), nodemass
 real :: dx,r2,r
 integer :: i,indexi 
 real :: h 
 h = 50.0

 nodecm = node % centerofmass
 nodemass = node % totalmass

 do i=1,10
  indexi = particlesindex(i)


enddo 

end subroutine get_poten_node

end module potendirect
