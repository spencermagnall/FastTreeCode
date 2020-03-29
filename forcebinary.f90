module poten
     use octreetype
     implicit none
     contains 
     subroutine get_accel(x,a,m,np,nodes)
      integer, intent(in) :: np
      real, intent(in) :: x(3,np)
      real, intent(out) :: a(3,np)
      real, intent(in) :: m(np)
      type(octreenode), intent(in) :: nodes(:)
      real :: mass1, mass2
      real :: cm1(3), cm2(3)
      integer :: leafnodesfound, nodeindex,i,j,childindex
      type(octreenode) :: childnodes(8)
      real :: dx(3), r2, r, a1(3), a2(3), anodes(3,8)
      
      leafnodesfound = 0
      nodeindex = 1
      !childnodes = 0.0
      do i=1,8
       call null_node(childnodes(i))
  	  enddo

      ! START FROM THE ROOT AND PROCESS THE TWO NODES
      ! SHOULD ONLY BE TWO LEAF NODES

      !do i=1, 8
       !childindex = nodes(1) % children(i)
       !if (nodes(childindex)% isLeaf .EQV. .true. .AND. leafnodesfound .EQ. 0 .AND. nodes(childindex) %data(1) .NE. 0) then
        !write(*, *) "Child index 1: ", childindex
        !child1 = nodes(childindex)
        !cm1(:) = nodes(childindex) % centerofmass
        !mass1 = nodes(childindex) % totalmass
        !leafnodesfound = leafnodesfound + 1
       !else if (nodes(childindex) % isLeaf .EQV. .true. .AND. leafnodesfound .EQ. 1 .AND. nodes(childindex) %data(1) .NE. 0 ) then
       	!write(*, *) "Child index 2: ", childindex
        !child2 = nodes(childindex)
        !cm2(:) = nodes(childindex) % centerofmass
        !mass2 = nodes(childindex) % totalmass
        !leafnodesfound = leafnodesfound + 1 
       !endif 
      !enddo

      ! Get all of the child nodes

      do i=1, 8
       ! Start from the root node
       childindex = nodes(1) % children(i)
       ! If  child node exists and has data 
       if (nodes(childindex) %isLeaf .EQV. .true. .AND. nodes(childindex) % data(1) .NE. 0) then
          childnodes(i) = nodes(childindex)
       endif
      enddo 

      ! FORCE CALCULATION

      !a = 0.0
      !a1 = 0.0
      !a2 =0.0 

      !dx = cm1-cm2
      !r2 = dot_product(dx,dx)
      !r = sqrt(r2)

      !a1 = -mass2*dx*(1/(r2*r+(50)**1.5))
      !dx = cm2-cm1
      !r2 = dot_product(dx,dx)
      !r = sqrt(r2)
      !a2 = -mass1*dx*(1/(r2*r+(50)**1.5))
      
      anodes = 0.0
      ! for all nodes
      do i=1, 8
       do j=1, 8
        if (childnodes(i)%data(1) .NE. 0 .AND. childnodes(j)%data(1) .NE. 0) then
         cm1 = childnodes(i) % centerofmass
         cm2 = childnodes(j) % centerofmass
         mass1 = childnodes(i) % totalmass
         mass2 = childnodes(j) % totalmass
         dx = cm1-cm2
         r2 = dot_product(dx,dx)
         r = sqrt(r2)
         write(*,*) "r is: ", r

         ! If nodes are not the same and are well seperated 
         if (i /= j .AND. r > 1000.0) then
          anodes(:,i) = anodes(:,i) - mass2*(1/(r2*r))*dx

         endif
        endif  
       enddo
      enddo  

      do i=1,8
       do j=1,10
        ! Node and Particle exist; propigate accel to particles
        if (childnodes(i)%data(1) .NE. 0 .AND. childnodes(i) % data(j) .NE. 0) then
          nodeindex = childnodes(i) % data(j)
          a(:,nodeindex) = anodes(:,i)
        endif  
       enddo
      enddo
      ! propigate this to all particles 
      !do i=1, 10
      !	write(*,*) "I: ",i
       !nodeindex = child1%data(i)
       !write(*,*) "NodeIndex: ",nodeindex
       !a(:,nodeindex) = a1
      !enddo

      !do i=1, 10
      !	write(*,*) "I: ",i
       !nodeindex = child2%data(i)
       !write(*,*) "NodeIndex: ",nodeindex
       !a(:,nodeindex) = a2

      !enddo 







      end subroutine get_accel

      

end module poten 
