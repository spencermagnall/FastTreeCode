module testdirect
 use potendirect

 implicit none 

 contains 
 subroutine test_direct(x,m,np)
  integer, intent(in) :: np 
  real, intent(in) :: x(3,np), m(np)
  real :: aintract(3,np), adirect(3,np)
  integer :: particleindex1(10), particleindex2(10),i


  aintract = 0.
  adirect = 0. 

  do i=1,10
   particleindex1(i) = i
   particleindex2(i) = i + 10
  enddo

  call get_accel_test(x,adirect,m,np)

  print*, adirect

  call get_accel(x,aintract,m,np,particleindex1)

  print*, aintract

  print*, "Delta a: ", adirect - aintract

end subroutine test_direct


end module testdirect
