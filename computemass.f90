module computemass
use octreetype
implicit none 
contains 
RECURSIVE subroutine get_com(x,v,m,np,nodes,currentnode,sumMass,cm)
 integer, intent(in) :: np
 real, intent(in) :: x(3,np)
 real, intent(in) :: v(3,np)
 real, intent(in) :: m(np)
 type(octreenode), intent(inout) :: nodes(:)
 integer, intent(inout) :: currentnode
 integer :: i, dataIndex, childindex
 real, intent(out) :: sumMass, cm(3)
 real :: massChild(8), cmChild(3,8)
 
 !print*, "Mass: ", m
 write(*, *) "CurrentNode: ", currentnode
 print*, "Summass"
 print*, "where are we "
 ! BASE CASE 
 ! If we have a leaf node 


 massChild = 0.0
 cmChild = 0.0
 if (nodes(currentnode) % isLeaf .EQV. .TRUE.) then
    sumMass = 0.0
    cm = 0.0

   ! Total Mass is the sum of all masses at leaf node
   do i=1, 10
    ! if there is a data point
    if (nodes(currentnode)%data(i) .NE. 0) then
     ! Get the index of the datapoint and then sum the Mass
     dataIndex = nodes(currentnode)%data(i)
     sumMass = sumMass + m(dataIndex)
     !write(*,*) "Summass"
    endif 
   enddo

   ! set the Total Mass for the node
   nodes(currentnode)%totalmass = sumMass 
   print*, "Total mass: ", sumMass

  ! Center of Mass is the center of Mass of all the data points 
  ! Weighted sum

  do i=1,10
   ! if there is a data point
    if (nodes(currentnode)%data(i) .NE. 0) then
     ! Get the index of the datapoint and then sum the Mass
     dataIndex = nodes(currentnode)%data(i)
     cm = cm + m(dataIndex)*x(:,dataIndex)
    endif 

  enddo
  ! Store Center of Mass
  ! This is a singularity if the totalmass
  ! of the node is 0
  if (nodes(currentnode)%totalmass > 0.) then
    cm = cm/(nodes(currentnode)%totalmass)
  else 
    cm = 0.
  endif 
  nodes(currentnode) % centerofmass = cm 

  return  ! Return sumMass and cm 


 else
  print*,"Else condition"
  ! For all children

   do i=1, 8
    ! GET COM OF CHILDREN AND TOTAL MASS
    childindex = nodes(currentnode) % children(i)
    if (childindex .NE. 0) then 
      call get_com(x,v,m,np,nodes,childindex,sumMass,cm)

      massChild(i) = sumMass
      cmChild(:,i) = cm
    
    endif 
   enddo 
 
   ! TOTAL MASS OF NODE IS SUM OF CHILDREN
   sumMass = SUM(massChild)

   cm = 0.0

   do i=1, 8
    cm = cm + massChild(i)*cmChild(:,i)
   enddo
   cm = cm / sumMass

   ! Store these values 

   nodes(currentnode) % totalmass = sumMass
   nodes(currentnode) % centerofmass = cm 

   print*, "Center of mass: ",cm
   print*, "Total mass: ", sumMass

   return   


 endif 



end subroutine get_com

subroutine compute_quads(x,m,np,nodes,endnode)
  integer, intent(in) :: np,endnode
  real, intent(in) :: x(3,np), m(np)
  type(octreenode), intent(inout) :: nodes(endnode)
  real :: quads(6), dx(3)
  integer :: i,j, particleindex 

  quads = 0.

  do i=1, endnode
    quads = 0.
    ! For all bodychildren of the node 
    do j=1, nodes(i) % bodychildpont
       particleindex = nodes(i) % bodychildren(j)
       print*, "particleindex: ",particleindex
       ! Particle distance from center of mass 
       dx = x(:, particleindex) - nodes(i) % centerofmass
       
       ! UNROLLED QUADS, DONT NEED ALL OF THESE BECAUSE OF SYM
       quads(1) = quads(1) + m(particleindex)*(3.*dx(1)*dx(1) - dot_product(dx,dx))
       quads(2) = quads(2) + m(particleindex)*(3.*dx(1)*dx(2))
       quads(3) = quads(3) + m(particleindex)*(3.*dx(1)*dx(3))
       !quads(2,1) = quads(2,1) + m(j)*(3.*dx(2)*dx(1))
       quads(4) = quads(4) + m(particleindex)*(3.*dx(2)*dx(2) - dot_product(dx,dx))
       quads(5) = quads(5) + m(particleindex)*(3.*dx(2)*dx(3))
       quads(6) = quads(6) + m(particleindex)*(3.*dx(3)*dx(3)- dot_product(dx,dx))
       


    enddo 

    nodes(i) % quads = quads 
    print*, "Quads are: ", quads 

  enddo 

  !STOP


  

end subroutine compute_quads
 
end module computemass
