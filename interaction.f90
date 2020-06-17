module interaction
 use octree
 use opening_criterion
 use taylor_expansions
 use potendirect 
 ! This is just for testing remove later
 use evaluate  

 contains 

 RECURSIVE subroutine interact(nodeindex1,nodeindex2,nodes,x,m,a,nopart,interactionlist)
  integer, intent(in) :: nopart 
  type(octreenode), intent(inout) :: nodes(:)
  real, intent(in) :: x(3,nopart), m(3,nopart)
  real, intent(inout) :: a(3,nopart)
  real, optional , intent(inout) :: interactionlist(:,:)
  integer, intent(inout) :: nodeindex1,nodeindex2
  type(octreenode) :: newnode1,newnode2, splitnode, regnode
  integer :: i, j, counter 
  real :: rmax1, rmax2,cm1(3),cm2(3)
  real :: fnode(20),quads(6)
  real :: dr,dx(3),dy,dz,totmass,r,r2
  integer :: particleindex(10),particleindex2(10) 
  real :: c0,c1(3),c2(3,3),c3(3,3,3)
  logical :: nodesAreEqual, flag
  integer :: regnodechild(128), splitnodeindex,regnodeindex,newnodeindex1,newnodeindex2,start
  integer, allocatable :: particlesforsum(:)
  logical :: node1empty, node2empty
  
  rmax1 = nodes(nodeindex1) % rmax
  rmax2 = nodes(nodeindex2) % rmax 

  quads(:) = 0.0
  counter = 0

  fnode(:) = 0.0
  particleindex = 0.0
  particleindex2 = 0.0

  regnodechild(:) = 0.
  c0 = 0.
  c1 = 0.
  c2 = 0.
  c3 = 0.

  print*, "Node indexes: ", nodeindex1, nodeindex2

  ! STILL NEED TO COMPUTE BODY-BODY, BODY-NODE, NODE_BODY 
  ! BODY SELF INTERACTION IS IGNORED
  nodesAreEqual = nodes_equal(nodes(nodeindex1),nodes(nodeindex2))
  if (nodesAreEqual) then
    ! If we are doing leafnode leafnode 
   if (nodes(nodeindex1) % isLeaf) then 
     ! Get index of all bodies 
     do i=1, 10
       if (nodes(nodeindex1) % data(i) /= 0) then
         particleindex(i) = nodes(nodeindex1) % data(i)
       endif  
     enddo 
    ! CALL DIRECT SUM 
    print*,"Direct sum"
    !print*, x
    call get_accel(x,a,m,np,particleindex,particleindex2)
    !print*,a
  !stop
    !return 
   
   ! CELL SELF INTERACTION IS SPLIT INTO MI BETWEEN SUBNODES
   else
    ! AT MOST 36 Interactions 
    ! 64 in this case but some will be duplicates ???
    start = 1
    do i=start,8
      do j=start,8

        ! Get the index of the sub-cells
         print*, "i, j: "
         print*, i
         print*, j
         !print*, nodes(nodeindex1) % children(i)
         !print*, nodes(nodeindex2) % children(j)
        ! If children exist
        if (nodes(nodeindex1) % children(i) .NE. 0 .AND. nodes(nodeindex2) % children(j) .NE. 0) then 
         !nodeindex1 = nodes(nodeindex1) % children(i)
         !nodeindex2 = nodes(nodeindex2) % children(j)
         !print*, nodeindex1
         !print*, nodeindex2
         counter = counter + 1
         print*,"Counter: ", counter

          ! Get the sub-cells
          !newnode1 = nodes(nodeindex1)
          !newnode2 = nodes(nodeindex2)
          newnodeindex1 = nodes(nodeindex1) % children(i)
          newnodeindex2 = nodes(nodeindex2) % children(j)

          !print*, "Interact called"
          ! call the interact for each of the sub-cells
          call interact(newnodeindex1, newnodeindex2, nodes,x,m,a,nopart)

        endif 

      enddo 
      start = start + 1
    enddo
    endif 




 !elseif (well_separated(nodes(nodeindex1),nodes(nodeindex2))) then 
  !elseif (.false.) then 
  elseif (.true.) then
    print*, "Well separated!"
   

  

  ! WELL SEPARATED NODE MI IS CALCULATED: TAYLOR COEFFs computed and added to 
  ! node data fields

 
    ! Call taylor COEFFs
    !print*, "Calling taylor coeff: "

    ! PUT MULTIPOLE STUFF HERE 
    ! ------------------------
    ! ------------------------
    ! ------------------------
    ! ------------------------

    if (nodes(nodeindex1) % totalmass /= 0. .AND. nodes(nodeindex2) %totalmass /= 0.) then
    cm1 = nodes(nodeindex1) % centerofmass
    cm2 = nodes(nodeindex2) % centerofmass

    print*, "particles in node 1"
    print*, nodes(nodeindex1) % data
    print*, "particles in node 2"
    print*, nodes(nodeindex2) % data


    call get_dx_dr(cm1,cm2,dx,dr)
    print*, "Cm of node1: ",cm1
    print*, "Cm of node2: ",cm2
    !dx = cm1(1) - cm2(1)
    print*, "Dx: ", dx
    !dy = cm1(2) - cm2(2)
    !print*, "Dy: ", dy
    !dz = cm1(3) - cm2(3)
    !print*, "Dz: ", dz
    !dr = 1./(sqrt(dot_product(cm1-cm2,cm1-cm2)))
    !print*, "Dr: ", dr 

    fnode = 0.0

    totmass = nodes(nodeindex2) % totalmass
    print*, "Totalmass: ",totmass
    !totmass = 1.0

    ! calculate real accel here 
    print*, "Real accel: "
     r2 = dot_product(dx,dx)
     r  = sqrt(r2)
     print*, -totmass*(1/((r2)**1.5))*dx 
     print*, "Real force: "
     print*, -totmass*(1/((r2)**1.5))*dx * nodes(nodeindex1) %totalmass

    c0 = 0.
    c1 = 0.
    c2 = 0.
    c3 = 0.
    call compute_fnode(dx(1),dx(2),dx(3),dr,totmass,quads,fnode)
    call compute_coeff(dx(1),dx(2),dx(3),dr,totmass,quads,c0,c1,c2,c3)

    ! store coeff for walk phase 
    !node1 % fnode = node1 % fnode + fnode
    print*, "Stored acccel: "
    nodes(nodeindex1) % c0 =  nodes(nodeindex1)%c0 + c0
    print*, nodes(nodeindex1) % c1
    print*, "calc accel: "
    print*,c1
    nodes(nodeindex1) % c1 = nodes(nodeindex1)%c1 + c1
    nodes(nodeindex1) % c2 = nodes(nodeindex1) % c2 + c2
    nodes(nodeindex1) % c3 = nodes(nodeindex1) % c3 + c3

    ! print poten
    !print*, "Poten is: ", fnode(20)
    !call poten_at_bodypos(cm1,cm2,c0,c1,c2,c3,poten(20))
    !print*, "Poten is: ", poten(20)

    !return 

    open(unit=77,file="wellseperated.txt")
    write(77,*), "Node 1:", nodeindex1, "Node 2: ", nodeindex2

    cm2 = nodes(nodeindex1) % centerofmass
    cm1 = nodes(nodeindex2) % centerofmass

    print*, "particles in node 1"
    print*, nodes(nodeindex2) % data
    print*, "particles in node 2"
    print*, nodes(nodeindex1) % data


    call get_dx_dr(cm1,cm2,dx,dr)
    print*, "Cm of node1: ",cm1
    print*, "Cm of node2: ",cm2
    !dx = cm1(1) - cm2(1)
    print*, "Dx: ", dx
    !dy = cm1(2) - cm2(2)
    !print*, "Dy: ", dy
    !dz = cm1(3) - cm2(3)
    !print*, "Dz: ", dz
    !dr = 1./(sqrt(dot_product(cm1-cm2,cm1-cm2)))
    !print*, "Dr: ", dr 

    fnode = 0.0

    totmass = nodes(nodeindex1) % totalmass
    print*, "Totalmass: ",totmass
    !totmass = 1.0

    ! calculate real accel here 
    print*, "Real accel: "
     r2 = dot_product(dx,dx)
     r  = sqrt(r2)
     print*, -totmass*(1/((r2)**1.5))*dx 
     print*, "Real force: "
     print*, -totmass*(1/((r2)**1.5))*dx * nodes(nodeindex2) %totalmass

    c0 = 0.
    c1 = 0.
    c2 = 0.
    c3 = 0.
    call compute_fnode(dx(1),dx(2),dx(3),dr,totmass,quads,fnode)
    call compute_coeff(dx(1),dx(2),dx(3),dr,totmass,quads,c0,c1,c2,c3)

    ! store coeff for walk phase 
    !node1 % fnode = node1 % fnode + fnode
    print*, "Stored acccel: "
    nodes(nodeindex2) % c0 =  nodes(nodeindex2)%c0 + c0
    print*, nodes(nodeindex2) % c1
    print*, "calc accel: "
    print*,c1
    nodes(nodeindex2) % c1 = nodes(nodeindex2)%c1 + c1
    nodes(nodeindex2) % c2 = nodes(nodeindex2) % c2 + c2
    nodes(nodeindex2) % c3 = nodes(nodeindex2) % c3 + c3

    ! print poten
    !print*, "Poten is: ", fnode(20)
    !call poten_at_bodypos(cm1,cm2,c0,c1,c2,c3,poten(20))
    !print*, "Poten is: ", poten(20)

    !return 
    !STOP

    open(unit=77,file="wellseperated.txt")
    write(77,*), "Node 1:", nodeindex2, "Node 2: ", nodeindex1

  endif 

  ! THE NODE WITH THE LARGER RMAX IS SPLIT; up to 8 new MI are created and processed 
  else
    print*, "Split nodes"

    if (rmax1 > rmax2) then
     splitnode = nodes(nodeindex1)
     splitnodeindex = nodeindex1
     regnode = nodes(nodeindex2)
     regnodeindex = nodeindex2
    else 
     splitnode = nodes(nodeindex2)
     splitnodeindex = nodeindex2
     regnode = nodes(nodeindex1)
     regnodeindex = nodeindex1
    endif 

    ! process MI's on node and splitnode children
    if (.not. splitnode % isLeaf ) then 
    do i=1,8
      splitnodeindex = splitnode % children(i)
      !newnode1 = nodes(nodeindex1)
      print*, "Internal splitnode case"
      call interact(regnodeindex,splitnodeindex,nodes,x,m,a,nopart)
      !call interact(splitnodeindex,regnodeindex,nodes,x,m,a,nopart)
    enddo

    ! LEAF-NODE NODE 
    elseif(splitnode % isLeaf .AND. .NOT. regnode % isLeaf) then
      !print*, "FIX THIS"

      ! Call poten function
      ! Need a way to get all of the children of a node for direct sum

      ! Get children of regular node 
      node1empty = .true.
      nodeindex2 = .true.

      do i=1,128
        if (regnode % bodychildren(i) /= 0) then
           node1empty = .false.
          !print*, i
          regnodechild(i) = regnode % bodychildren(i)
        endif 
      enddo 

      do i=1,10
        if (splitnode % data(i) /= 0) then
          node2empty = .false.
          particleindex(i) = splitnode % data(i)
        endif 
      enddo 

      ! CALL DIRECT SUM

      print*, node1empty, node2empty


      if (.NOT. node1empty .and. .not. node2empty) then
        !call get_accel(x,a,m,np,particleindex,regnodechild)
        call get_accel_leafnode(x,a,m,np,particleindex,regnodechild)
        call get_accel_leafnode(x,a,m,np,regnodechild,particleindex)
      endif 
     !call get_accel_leafnode(x,a,m,np,particleindex,regnodechild)

     !return

      ! Get children of a leaf node 

    ! LEAF-NODE LEAF-NODE
    else 
      print*, "Leafnode-Leafnode"

      node1empty = .true.
      node2empty = .true.

      ! This is direct sum as before 

       !Get index of all bodies 
     do i=1, 10
      !print*, i
       if (nodes(nodeindex1) % data(i) /= 0) then
          node1empty = .false.
         !print*, "Crash 1"
         particleindex(i) = nodes(nodeindex1) % data(i)
         !print*, "particleindex: ", particleindex(i)
       endif
       if (nodes(nodeindex2) % data(i) /= 0) then
          node2empty = .false.
         !print*, "Crash2"
         particleindex2(i) = nodes(nodeindex2) % data(i)
         !print*, "particleindex2: ", particleindex2(i)
       endif

     enddo

     print*, "Finished getting bodies"
     print*, node1empty
     print*, node2empty
     ! CALL DIRECT SUM 
    if (.NOT. node1empty .and. .not. node2empty) then
      !call get_accel(x,a,m,np,particleindex,particleindex2)
      call get_accel_leafnode(x,a,m,np,particleindex,particleindex2)
      call get_accel_leafnode(x,a,m,np,particleindex2,particleindex)
    endif 
    !call get_accel_leafnode(x,a,m,np,particleindex,particleindex2)

    !return


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

function get_direct_sum_particles(particleindex1,particleindex2) result(particles)
  integer,intent(in) :: particleindex1(:),particleindex2(:)
  integer :: index1, index2, listsize,i
  integer, allocatable :: particles(:)

  index1 = size(particleindex1)
  index2 = size(particleindex2)
  listsize = size(particleindex1) + size(particleindex2)
  allocate(particles(listsize))
  particles = 0.
  
  do i=1, index1
    particles(i) = particleindex1(i)
  enddo 

  i = index1 + 1

  do i=1, index2
    if (.NOT. in_list(particles,particleindex2(i))) then
      particles(i+index1) = particleindex2(i)
    endif 
  enddo 




end function get_direct_sum_particles

LOGICAL function in_list(list,index)
 integer, intent(in) :: list(:),index
 integer :: listsize,i

 listsize = size(list)

 ! This will be n^2 but n should be small so its fine 
 in_list = .FALSE.
 do i=1,listsize
  if (list(i) == index) then
    in_list = .TRUE.
  endif 
 enddo 

end function in_list

subroutine get_dx_dr(x1,x2,dx,dr)
  real, intent(in) :: x1(3),x2(3)
  real, intent(out) :: dx(3),dr
  real :: h


  h = 0.
  !h = 0.
  dx = x1 - x2
  dr = 1./(sqrt(dot_product(dx,dx) + h**2))

 end subroutine get_dx_dr

 
end module interaction 