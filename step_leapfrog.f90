module step_leapfrog
 use poten
 implicit none
 contains
  subroutine step(x,v,a,m,dt,np,nodes)

   integer, intent(in) :: np
   real, intent(inout) :: x(3,np), v(3,np), a(3, np)
   real, intent(in) :: m(np)
   real, intent(in) :: dt
   type(octreenode), optional, intent(in) :: nodes(np)
   integer :: i 

    ! LEAPFROG 
	do i=1, np
       v(:,i) = v(:,i) + 0.5*dt*a(:,i)
       x(:,i) = x(:,i) + dt*v(:,i)
	enddo
    
    if (present(nodes)) then
        call get_accel(x,a,m,np,nodes)
    endif 

    do i=1, np
       v(:,i) = v(:,i) + 0.5*dt*a(:,i)
    enddo

  end subroutine step

end module step_leapfrog
