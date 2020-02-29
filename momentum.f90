module momentum
 implicit none 
	contains 
 subroutine get_momentum(x,v,m,p,angm,np)
  real, intent(in) :: x(3,np)
  real, intent(in) :: v(3,np)
  real, intent(in) :: m(np)
  !real :: p(3)
  real, intent(out) :: p(3)
  real, intent(out) :: angm(3)
  real :: p1(3), p2(3)
  integer, intent(in) :: np
  integer :: i
  integer :: start, end
  real :: angi(3)

  p = 0.0
  p1= 0.
  p2 =0.
  angm = 0.0
  angi = 0.0
  end = np/2
  do i=1, end
  	! should be a m in here 
  	p1= p1 + m(i)*v(:,i)
  enddo

  start = end+1
  end = np
  do i=start, end
  	p2 = p2 + m(i)*v(:,i)
  !pmag = p(1)*p(1) + p(2)*p(2) + p(3)*p(3)
  !pmag = sqrt(pmag)
  enddo

  do i=1, np
     call cross_product(x(:,i),v(:,i),angi)
     angm = angm +m(i)*angi
  enddo 
  write(*,*) p1
  write(*,*) p2
  p = p1 + p2 
  write(*,*) p
 end subroutine get_momentum

 subroutine cross_product(a,b,c)
   real, intent(in) :: a(3), b(3)
   real, intent(out) :: c(3)

   c(1) = a(2)*b(3) - a(3)*b(2)
   c(2) = a(3)*b(1) - a(1)*b(3)
   c(3) = a(1)*b(2) - a(2)*b(1)

 end subroutine cross_product
end module momentum

