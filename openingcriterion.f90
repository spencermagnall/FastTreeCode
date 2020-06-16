module opening_criterion
 use octreetype
 implicit none
 contains

 logical function well_separated(node1,node2) result(bool)
  type(octreenode), intent(in) :: node1, node2
  real :: cm1(3), cm2(3), dx(3), zmag
  real :: theta,rmax1,rmax2
  !logical :: bool

  !bool =.true.
  !return 
  theta = 0.3

  cm1 = node1 % centerofmass
  cm2 = node2 % centerofmass

  ! Get the magnitude of the difference between CoM's
  dx = cm1 - cm2
  print*,dx
  zmag = norm2(dx)

  ! get the rmax for nodes 1 and 2
  rmax1 = node1 % rmax
  rmax2 = node2 % rmax 


  print*,"Zmag: ", zmag
  print*, "rmax1 + rmax2/theta: ", (rmax1 + rmax2)/theta
  if (zmag > (rmax1+rmax2)/theta) then
    bool = .TRUE.
  else 
    bool = .FALSE.
  endif  
 end function well_separated

 subroutine distance_to_corner(node,rcorn)
  type(octreenode), intent(in) :: node
  real, intent(out) :: rcorn
  ! Stores the corners of the Node
  real :: corns(3,8)
  real :: nodeorigin(3), nodecm(3), nodesize,rmax,dx(3), normx
  integer i

  corns = 0.0
  nodeorigin = node % origin
  nodesize = node % size
  nodecm = node % centerofmass
  rmax = 0.0
  dx = 0.0
  normx = 0.0

  ! first find out what the corners are
  ! + + + 
  corns(:,1) = nodeorigin + nodesize

  ! - + +
  corns(:,2) = nodeorigin + nodesize
  corns(1,2) = nodeorigin(1) - nodesize

  ! - - +
  corns(:,3) = nodeorigin - nodesize
  corns(3,3) = nodeorigin(3) + nodesize

  ! + - + 
  corns(:,4) = nodeorigin + nodesize
  corns(2,4) = nodeorigin(2) - nodesize

  ! + + -
  corns(:,5) = nodeorigin + nodesize
  corns(3,5) = nodeorigin(3) - nodesize

  ! - + -
  corns(:,6) = nodeorigin - nodesize
  corns(2,6) = nodeorigin(2) + nodesize

  ! - - -
  corns(:,7) = nodeorigin - nodesize

  ! + - -
  corns(:,8) = nodeorigin - nodesize
  corns(1,8) = nodeorigin(1) + nodesize 

  ! find the largest corner 
  do i=1,8

   ! find magnitude of vector 
   dx = nodecm-corns(:,i)
   normx = norm2(dx)

   if (normx > rmax) then
    rmax = normx
   endif 
  enddo

  ! Assign rcorn the distance to the furthest corner from CoM
  rcorn = rmax 


 end subroutine distance_to_corner

 ! Make this recursive 
 RECURSIVE subroutine find_rmax(x,nodes,currentnode,rmax)
 ! This should be run after tree construction
 type(octreenode), intent(out) :: nodes(:)
 integer, intent(out) :: currentnode
 real, intent(out) :: rmax
 real, intent(in) :: x(:,:)
 real :: rcorn, rmaxsubnode,rmaxcurrent
 ! Shouldn't be a fixed value 
 real :: children(8),nodecm(3),dx(3),dxnorm,dnode
 !real :: size
 integer i,childindex

 rcorn = 0.0
 
 rmaxsubnode = 0.0
 children = 0.0
 dx = 0.0
 dxnorm = 0.0 
 childindex = 0
 dnode = 0.0
 rmaxcurrent = HUGE(0.0)
 nodecm = nodes(currentnode) % centerofmass

 write(*,*) "CurrentNode is: ", currentnode

 ! Base case leaf node 
 if (nodes(currentnode) %isLeaf .EQV. .TRUE.) then
      rmax = 0.0
     ! Get the distance to corner
     call distance_to_corner(nodes(currentnode),rcorn)
     
     do i=1,10
       ! probaly need to check for null nodes
       if (nodes(currentnode)%data(i) .NE. 0) then
       ! find the distance between particle and node cm 
          dx = x(:,nodes(currentnode)%data(i)) - nodecm
          dxnorm = norm2(dx)
          if (dxnorm > rmaxsubnode) then 
              rmaxsubnode = dxnorm
          endif
       endif 
     enddo 

     ! rmax is the smaller of the two 
     if (rmaxsubnode < rcorn) then
         rmax = rmaxsubnode
     else 
         rmax = rcorn
     endif 

    ! Now store this value at the node 
    nodes(currentnode) % rmax = rmax 

    ! return up the tree
    return 
! We are at interior node so collect rmax from children 
else 
    ! find the distance_to_corner 
    call distance_to_corner(nodes(currentnode),rcorn)

    ! now we are finding rmax of all the children
    do  i=1,8
      print*,i
      ! if node is not null 
      if (nodes(currentnode)%children(i) .NE. 0) then
         ! get node CoM
         childindex = nodes(currentnode)%children(i)
         nodecm = nodes(childindex) % centerofmass
         dnode = norm2(nodecm - nodes(currentnode)%centerofmass)
         call find_rmax(x,nodes,nodes(currentnode)%children(i),rmaxsubnode)
         ! EQ [3] from Benz et al. 1990
         if (rmaxsubnode + dnode .LT. rmaxcurrent) then 
            rmaxcurrent = rmaxsubnode + dnode
            write(*,*) "rmaxcurrent: ",rmaxcurrent
         endif
      endif 
    enddo 

    ! is corner sphere or particle sphere smaller 
    if (rmaxcurrent < rcorn) then
        rmax = rmaxcurrent
    else 
        rmax  = rcorn 
    endif 

    ! Now store the value at this Node
    nodes(currentnode) % rmax = rmax

    ! Now return up the tree
    return
endif     



end subroutine find_rmax


end module opening_criterion