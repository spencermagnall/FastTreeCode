module output 
 implicit none
 integer, private :: nfile = 1


contains
 subroutine write_output(x,v,m,np,time)
  integer, intent(in) :: np 
  real, intent(in) :: x(3,np), v(3,np)
  real, intent(in) :: m(np)
  real, intent(in) :: time
  integer :: i

  character(len=100) :: filename
  
  !write to a sequences of files called snap_00000 , snap_000001 etc
  write(filename, "(a,i5.5)") 'snap_',nfile
  nfile = nfile + 1

  print "(a,f8.3)", ' writing '//trim(filename)//' t = ',time
  open(unit=67, file=filename,status='replace')
  write(67,*) time
  do i=1, np
     write(67,*) x(:,i), v(:,i), m(i)
  enddo
  close(unit=67)

 end subroutine write_output
end module output