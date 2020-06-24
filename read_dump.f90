module read_dump
 implicit none 
 contains 
 subroutine read_ascii(x,v,m,np,filename)
  character (len = *), intent(in) :: filename 
  integer, intent(in) :: np
  real, intent(inout) :: x(3,np), v(3,np), m(np)
  integer :: i
  real :: a(3,np)

  open(unit=11,file = filename,status='OLD',action='read')
  !print*
  !read(11,*)
  !read(11,*)
  print*, filename
  !STOP
  do i=1,np
   read(11,*) x(1,i), x(2,i), x(3,i), v(1,i), v(2,i), v(3,i) !, a(1,i), a(2,i), a(3,i), m(i)
   write(*,*) x(1,i), x(2,i), x(3,i), v(1,i), v(2,i), v(3,i) !, a(1,i), a(2,i), a(3,i), m(i)
  enddo 

  print*, "pos:"
  print*,x
  print*,"Velocity: "
  print*, v
  close(11)
  m = 1.0/np 
  !STOP


 end subroutine read_ascii
end module read_dump 