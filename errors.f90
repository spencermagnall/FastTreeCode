module errors
implicit none 
contains
subroutine get_rms(a,aexact,np,rms)
  integer, intent(in) :: np
  real, intent(in) :: a(3,np), aexact(3,np)
  real, intent(out) :: rms
  integer :: i


  rms = 0.
  do i=1,np
     rms = (norm2(aexact(:,i)-a(:,i)))**2
  enddo 

  rms = sqrt(rms/np)

end subroutine get_rms 

end module errors 