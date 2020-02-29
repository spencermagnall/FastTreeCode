module interaction
 use octree

 contains 

 RECURSIVE subroutine interact(node1,node2,nodes,nopart)
  integer :: nopart 
  type(octreenode), intent(inout) :: node1, node2, nodes(nopart),newnode1,newnode2
  integer :: i, j 
  ! BODY SELF INTERACTION IS IGNORED

  if (node1 == node2 ) then
   if (node1 % isBody .EQ. true) then  
    return
   ! CELL SELF INTERACTION IS SPLIT INTO MI BETWEEN SUBNODES
   else
    ! AT MOST 36 Interactions 
    do i=1,8
      do j=1,8
        newnode1 = node1
      enddo 
    enddo


  end if    

  

  ! WELL SEPERATED NODE MI IS CALCULATED: TAYLOR COEFFs computed and added to 
  ! node data fields

  ! THE NODE WITH THE LARGER RMAX IS SPLIT; up to 8 new MI are created and processed 

 end subroutine interact
 
end module interaction 