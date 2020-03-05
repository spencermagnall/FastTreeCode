module potendirect
 implicit none 
 contains 
 subroutine get_accel(x,a,m,np)
 	integer, intent(in) :: np
 	real, intent(in) :: x(3,np)
 	real, intent(out) :: a(3,np)
 	real, intent(in) :: m(np)
 	real :: massratio
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

end module potendirect
