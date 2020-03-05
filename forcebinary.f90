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
      integer :: leafnodesfound, nodeindex,i,childindex
      type(octreenode) :: child1
      type(octreenode) :: child2
      real :: dx(3), r2, r, a1(3), a2(3)
      leafnodesfound = 0
      nodeindex = 1


      ! START FROM THE ROOT AND PROCESS THE TWO NODES
      ! SHOULD ONLY BE TWO LEAF NODES

      do i=1, 8
       childindex = nodes(1) % children(i)
       if (nodes(childindex)% isLeaf .EQV. .true. .AND. leafnodesfound .EQ. 0) then
       	write(*, *) "Child index 1: ", childindex
        child1 = nodes(childindex)
        cm1(:) = nodes(childindex) % centerofmass
        mass1 = nodes(childindex) % totalmass
        leafnodesfound = leafnodesfound + 1
       else if (nodes(childindex) % isLeaf .EQV. .true. .AND. leafnodesfound .EQ. 1) then
       	write(*, *) "Child index 2: ", childindex
        child2 = nodes(childindex)
        cm2(:) = nodes(childindex) % centerofmass
        mass2 = nodes(childindex) % totalmass
        leafnodesfound = leafnodesfound + 1 
       endif 
      enddo

      ! FORCE CALCULATION

      a = 0.0

      dx = cm1-cm2
      r2 = dot_product(dx,dx)
      r = sqrt(r2)

      a1 = -mass2*dx*(1/(r2*r))
      a2 = -a1


      ! propigate this to all particles 
      do i=1, 10
      	write(*,*) "I: ",i
       nodeindex = child1%data(i)
       write(*,*) "NodeIndex: ",nodeindex
       a(:,nodeindex) = a1
      enddo

      do i=1, 10
       nodeindex = child2%data(i)
       a(:,nodeindex) = a2
      enddo 







      end subroutine get_accel

      

end module poten 
