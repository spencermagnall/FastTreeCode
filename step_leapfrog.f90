module step_leapfrog
 use potendirect
 use octree 
 use computemass
 use opening_criterion
 use evaluate 
 use interaction  
 use delta_t
 implicit none
 contains
  subroutine step(x,v,a,m,dt,np,nodes)

   integer, intent(in) :: np
   real, intent(inout) :: x(3,np), v(3,np), a(3,np)
   real, intent(inout) :: m(np)
   real, intent(inout) :: dt
   type(octreenode), allocatable, optional, intent(out) :: nodes(:)
   integer :: endnode
   real :: c0,c1(3),c2(3,3),c3(3,3,3), atest(3,np)
   integer :: i, rootnode
   real :: sumMass, cm(3),rmax,asum(3),boxsize,origin(3)

   rootnode = 1
   sumMass = 0.0
   cm = 0.0
   rmax = 0.0 
   c0 = 0.
   c1 = 0.
   c2 = 0.
   c3 = 0.
   atest = 0.
   asum = 0.
   endnode = 0
   i = 0.
   boxsize = 0.
   origin = 0.

   

    ! LEAPFROG 
	do i=1, np
       v(:,i) = v(:,i) + 0.5*dt*a(:,i)
       x(:,i) = x(:,i) + dt*v(:,i)
	enddo
    a = 0.
    if (present(nodes)) then
        !call cleartree(nodes)
        ! Find the boxsize and origin for the tree 
        call get_optimal_boxsize(x,m,np,boxsize,origin)
        print*, boxsize
        !stop
        call maketree(nodes,x,v,m,a,np,endnode)
        print*, "Endnode: ", endnode
        !if (endnode /= 65) stop
        !stop
        call get_com(x,v,m,np,nodes,rootnode,sumMass,cm)
        !call compute_quads(x,m,np,nodes,endnode)
        !call get_accel_test(x,atest,m,np)
        rootnode = 1
        rmax = 0.
        cm = 0.
        call find_rmax(x,nodes,rootNode,rmax)
        rootnode = 1
        call interact(rootNode,rootNode,nodes,x,m,a,np)
        !stop
        !print*, "Accel:"
        !print*, a
        !STOP
        c0 = 0.
        c1 = 0.
        c2 = 0.
        c3 = 0.
        cm = 0.
        rootnode = 1
        !call evaluate_gravity(rootnode,nodes,cm,c0,c1,c2,c3,x,a,asum)
        call evaluate_gravity_stack(nodes,x,a)
        !stop
        print*, "Accel:"
        do i=1, np
        print*, a(:,i)
        enddo 
        print*, "Accel dir:"
        do i=1,np
        print*, atest(:,i)
        enddo 
        print*, "delta a: "
        do i=1,np
        print*, atest(:,i)-a(:,i)
        enddo 
        print*, "Net accel from tree: ", asum 
        !STOP
    endif 


    !call get_dtnew(0.1,a,dt,np)

    do i=1, np
       v(:,i) = v(:,i) + 0.5*dt*a(:,i)
    enddo

  end subroutine step

end module step_leapfrog
