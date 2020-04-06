module interaction
 use octree
 use opening_criterion
 use taylor_expansions
 use potendirect  

 contains 

 RECURSIVE subroutine interact(node1,node2,nodes,x,m,poten,nopart)
  integer, intent(in) :: nopart 
  type(octreenode), intent(inout) :: node1, node2, nodes(:)
  real, intent(in) :: x(3,nopart), m(3,nopart)
  real, intent(inout) :: poten(nopart)
  type(octreenode) :: newnode1,newnode2, splitnode, regnode
  integer :: i, j, nodeindex1, nodeindex2,counter 
  real :: rmax1, rmax2,cm1(3),cm2(3)
  real :: fnode(20),quads(6)
  real :: dr,dx,dy,dz,totmass
  integer :: particleindex(10),particleindex2(10) 

  logical :: nodesAreEqual
  
  rmax1 = node1 % rmax
  rmax2 = node2 % rmax 

  quads(:) = 0.0
  counter = 0

  fnode(:) = 0.0
  particleindex = 0.0
  particleindex2 = 0.0

  ! STILL NEED TO COMPUTE BODY-BODY, BODY-NODE, NODE_BODY 
  ! BODY SELF INTERACTION IS IGNORED
  nodesAreEqual = nodes_equal(node1,node2)
  if (nodesAreEqual .EQV. .true. ) then
    ! If we are doing leafnode leafnode 
   if (node1 % isLeaf .EQV. .true.) then 
     ! Get index of all bodies 
     do i=1, 10
       if (node1 % data(i) /= 0) then
         particleindex(i) = node1 % data(i)
       endif  
     enddo 
    ! CALL DIRECT SUM 
    call get_poten(x,poten,m,np,particleindex,particleindex)
   
   ! CELL SELF INTERACTION IS SPLIT INTO MI BETWEEN SUBNODES
   else
    ! AT MOST 36 Interactions 
    ! 64 in this case but some will be duplicates ???
    do i=1,8
      do j=1,8

        ! Get the index of the sub-cells

        ! If children exist
        if (node1 % children(i) .NE. 0 .AND. node2 % children(j) .NE. 0) then 
         nodeindex1 = node1 % children(i)
         nodeindex2 = node2 % children(j)
         print*, nodeindex1
         print*, nodeindex2
         counter = counter + 1
         print*,"Counter: ", counter

          ! Get the sub-cells
          newnode1 = nodes(nodeindex1)
          newnode2 = nodes(nodeindex2)

          ! call the interact for each of the sub-cells
          call interact(newnode1, newnode2, nodes,x,m,poten,nopart)

        endif 

      enddo 
    enddo
    endif 




  elseif (well_seperated(node1,node2)) then 


  

  ! WELL SEPERATED NODE MI IS CALCULATED: TAYLOR COEFFs computed and added to 
  ! node data fields

 
    ! Call taylor COEFFs
    print*, "Calling taylor coeff: "

    ! PUT MULTIPOLE STUFF HERE 
    ! ------------------------
    ! ------------------------
    ! ------------------------
    ! ------------------------
    cm1 = node1 % centerofmass
    cm2 = node2 % centerofmass
    dx = cm1(1) - cm2(1)
    print*, "Dx: ", dx
    dy = cm1(2) - cm2(2)
    print*, "Dy: ", dy
    dz = cm1(3) - cm2(3)
    print*, "Dz: ", dz
    dr = 1./(norm2(cm1-cm2))
    print*, "Dr: ", dr 

    fnode = 0.0

    totmass = node2 % totalmass
    !totmass = 1.0

    call compute_fnode(dx,dy,dz,dr,totmass,quads,fnode)

    ! store coeff for walk phase 
    node1 % fnode = node1 % fnode + fnode

    ! print poten
    print*, "Poten is: ", fnode(20)


  ! THE NODE WITH THE LARGER RMAX IS SPLIT; up to 8 new MI are created and processed 
  else
    print*, "Split nodes"

    if (rmax1 > rmax2) then
     splitnode = node1
     regnode = node2
    else 
     splitnode = node2
     regnode = node1
    endif 

    ! process MI's on node and splitnode children
    if (.not. splitnode % isLeaf ) then 
    do i=1,8
      nodeindex1 = splitnode % children(i)
      newnode1 = nodes(nodeindex1)
      call interact(regnode,newnode1,nodes,x,m,poten,nopart)
    enddo

    ! LEAF-NODE NODE 
    elseif(splitnode % isLeaf .AND. .NOT. regnode % isLeaf) then

      ! Call poten function
      ! Need a way to get all of the children of a node for direct sum

    ! LEAF-NODE LEAF-NODE
    else 

      ! This is direct sum as before 

       !Get index of all bodies 
     do i=1, 10
       if (node1 % data(i) /= 0) then
         particleindex(i) = node1 % data(i)
       endif
       if (node2 % data(i) /= 0) then
         particleindex2(i) = node2 % data(i)
       endif

     enddo

     ! CALL DIRECT SUM 
    call get_poten(x,poten,m,np,particleindex,particleindex2)


    endif 


  endif 


 end subroutine interact

 LOGICAL function nodes_equal(node1, node2) result(bool)
  type(octreenode), intent(in) :: node1, node2
  real :: cm1(3), cm2(3)

  ! Compare the cm of nodes, if equal nodes are equal
  cm1 = node1 % centerofmass
  cm2 = node2 % centerofmass

  if (cm1(1) .eq. cm2(1) .AND. cm1(2) .eq. cm2(2) .AND. cm1(3) .eq. cm2(3)) then
    bool = .true. 
  else 
    bool = .false.
  endif 

end function nodes_equal
 
end module interaction 