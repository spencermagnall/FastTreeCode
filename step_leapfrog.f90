module step_leapfrog
 use potendirect
 use octree 
 use computemass 
 implicit none
 contains
  subroutine step(x,v,a,m,dt,np,nodes)

   integer, intent(in) :: np
   real, intent(inout) :: x(3,np), v(3,np), a(3, np)
   real, intent(in) :: m(np)
   real, intent(in) :: dt
   type(octreenode), allocatable, optional, intent(out) :: nodes(:)
   integer :: i, rootnode
   real :: sumMass, cm(3)

   rootnode = 1
   sumMass = 0.0
   cm = 0.0
    ! LEAPFROG 
	do i=1, np
       v(:,i) = v(:,i) + 0.5*dt*a(:,i)
       x(:,i) = x(:,i) + dt*v(:,i)
	enddo
    
    if (present(nodes)) then
        !call cleartree(nodes)
        call maketree(nodes,x,v,a,np)
        call get_com(x,v,m,np,nodes,rootnode,sumMass,cm)
        call get_accel(x,a,m,np)
    endif 

    do i=1, np
       v(:,i) = v(:,i) + 0.5*dt*a(:,i)
    enddo

  end subroutine step

end module step_leapfrog
