module step_leapfrog
 use potendirect
 use octree 
 use computemass
 use opening_criterion
 use evaluate 
 use interaction  
 implicit none
 contains
  subroutine step(x,v,a,m,dt,np,nodes)

   integer, intent(in) :: np
   real, intent(inout) :: x(3,np), v(3,np), a(3, np)
   real, intent(in) :: m(np)
   real, intent(in) :: dt
   type(octreenode), allocatable, optional, intent(out) :: nodes(:)
   real :: c0,c1(3),c2(3,3),c3(3,3,3)
   integer :: i, rootnode
   real :: sumMass, cm(3),rmax

   rootnode = 1
   sumMass = 0.0
   cm = 0.0
   rmax = 0.0 
   c0 = 0.
   c1 = 0.
   c2 = 0.
   c3 = 0.

   

    ! LEAPFROG 
	do i=1, np
       v(:,i) = v(:,i) + 0.5*dt*a(:,i)
       x(:,i) = x(:,i) + dt*v(:,i)
	enddo
    a = 0.
    if (present(nodes)) then
        !call cleartree(nodes)
        call maketree(nodes,x,v,a,np)
        call get_com(x,v,m,np,nodes,rootnode,sumMass,cm)
        rootnode = 1
        rmax = 0.
        cm = 0.
        call find_rmax(x,nodes,rootNode,rmax)
        call interact(nodes(1),nodes(1),nodes,x,m,a,np)
        print*, "Accel:"
        print*, a
        !STOP
        c0 = 0.
        c1 = 0.
        c2 = 0.
        c3 = 0.
        call evaluate_gravity(nodes(1),nodes,cm,c0,c1,c2,c3,x,a)
        print*, "Accel:"
        print*, a
        STOP
    endif 

    do i=1, np
       v(:,i) = v(:,i) + 0.5*dt*a(:,i)
    enddo

  end subroutine step

end module step_leapfrog
