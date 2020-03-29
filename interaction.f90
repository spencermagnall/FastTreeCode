module interaction
 use octree
 use opening_criterion
 use taylor_expansions 

 contains 

 RECURSIVE subroutine interact(node1,node2,nodes,nopart)
  integer :: nopart 
  type(octreenode), intent(inout) :: node1, node2, nodes(:)
  type(octreenode) :: newnode1,newnode2, splitnode, regnode
  integer :: i, j, nodeindex1, nodeindex2,counter 
  real :: rmax1, rmax2,cm1(3),cm2(3)
  real :: fnode(20),quads(6)
  real :: dr,dx,dy,dz,totmass
  logical :: nodesAreEqual
  
  rmax1 = node1 % rmax
  rmax2 = node2 % rmax 

  quads(:) = 0.0
  counter = 0

  fnode(:) = 0.0

  ! BODY SELF INTERACTION IS IGNORED
  nodesAreEqual = nodes_equal(node1,node2)
  if (nodesAreEqual .EQV. .true. ) then
   if (node1 % isBody .EQV. .true.) then  
    return
   ! CELL SELF INTERACTION IS SPLIT INTO MI BETWEEN SUBNODES
   else
    ! AT MOST 36 Interactions 
    ! 64 in this case but some will be duplicates ???
    do i=1,8
      do j=1,8

        ! Get the index of the sub-cells

        ! If children exist
        if (i /= j .AND. node1 % children(i) .NE. 0 .AND. node2 % children(j) .NE. 0) then 
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
          call interact(newnode1, newnode2, nodes, nopart)

        endif 

      enddo 
    enddo
    endif 




  else 

  

  ! WELL SEPERATED NODE MI IS CALCULATED: TAYLOR COEFFs computed and added to 
  ! node data fields

  if (well_seperated(node1,node2)) then 

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
    do i=1,8
      nodeindex1 = splitnode % children(i)
      newnode1 = nodes(nodeindex1)
      call interact(regnode,newnode1,nodes,nopart)
    enddo


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