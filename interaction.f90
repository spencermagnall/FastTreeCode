module interaction
 use octree

 contains 

 RECURSIVE subroutine interact(node1,node2,nodes,nopart)
  integer :: nopart 
  type(octreenode), intent(inout) :: node1, node2, nodes(nopart),newnode1,newnode2
  integer :: i, j, nodeindex1, nodeindex2 
  ! BODY SELF INTERACTION IS IGNORED

  if (node1 == node2 ) then
   if (node1 % isBody .EQ. true) then  
    return
   ! CELL SELF INTERACTION IS SPLIT INTO MI BETWEEN SUBNODES
   else
    ! AT MOST 36 Interactions 
    do i=1,8
      do j=1,8
        ! Get the index of the sub-cells

        ! If children exist
        if (node1 % children(i) .NE. 0 .AND. node2 % children(j) .NE. 0) then 
         nodeindex1 = node1 % children(i)
         nodeindex2 = node2 % children(j)

          ! Get the sub-cells
          newnode1 = nodes(nodeindex1)
          newnode2 = nodes(nodeindex2)

          ! call the interact for each of the sub-cells
          call interact(newnode1, newnode2, nodes, nopart)

        endif 

      enddo 
    enddo


  else     

  ! WELL SEPERATED NODE MI IS CALCULATED: TAYLOR COEFFs computed and added to 
  ! node data fields

  ! THE NODE WITH THE LARGER RMAX IS SPLIT; up to 8 new MI are created and processed 

 endif 

 end subroutine interact
 
end module interaction 