module octree
 use octreetype
 implicit none 

 ! the opening criterion 
 real, public :: theta = 0.5

 integer:: maxnodes = 10000
 integer :: totalnodes  = 0

 public :: maketree

 contains 

  subroutine maketree(nodes,x,v,a,np)
   integer, intent(in) :: np
   type(octreenode), allocatable, intent(out) :: nodes(:)
   real, intent(in) :: x(:,:), v(:,:), a(:,:)
   integer :: i
   integer :: currentnode
   integer :: endnode

   ! allocate nodes for tree
   allocate(nodes(maxnodes))

   ! Data and children for each node = 0 
   ! i.e dosen't exist 
   do i=1, maxnodes
    call null_node(nodes(i))
  enddo



   ! setup the root node
   call new_node(nodes(1),100000.,(/0.0,0.0,0.0/))
   

   ! endnode is now root
   endnode = 1
   ! currentnode is also root
   currentnode = 1

   ! iterate through all particles 
   do i =1, np
    write(*,*) "Particle: ",i
    write(*,*) "Endnode: " , endnode
    write(*,*) "currentnode: ", currentnode
    ! if tree is full expand it 
    !if (endnode + 10 >= maxnodes) then
    !   oldnode = endnode
    !   write(*,*) "Working before resize"
     !  call resize_nodes(nodes)
     !  write(*,*) "resize works"
     !   do j=oldnode, maxnodes
            ! this is wrong and breaking code
            ! go from oldsize to maxnodes
      !      nodes(j) % children(:) = 0
      !      nodes(j) % data = 0
       !     nodes(j) % isLeaf = .FALSE. 
       ! enddo
    !endif

    !write(*,*) "Resize finished"

    ! insert particle in tree
    currentnode = 1
    call insert_particle(nodes,x,v,a,currentnode,i,endnode)
    write(*,*) " Insert finished"
    !call print_tree(nodes,x,0,1)
    write(*,*) endnode

   enddo

   write(*,*) "Tree Built!"     

  end subroutine maketree

  subroutine sort_by_octant(x, np, origin, x1, x2, x3, x4, x5, x6, x7, x8)
    ! work out what octant each particle is in
    ! so we can paralize  
    real, intent(in) :: x(:,:)
    real, intent(in) :: origin(3)
    integer, intent(in) :: np
    real, allocatable, intent(out) :: x1(:,:), x2(:,:),x3(:,:),x4(:,:),x5(:,:),x6(:,:),x7(:,:),x8(:,:)
    integer :: octant, i, length
    integer :: dim

    dim = 3
    length = 1

    allocate(x1(dim,np),x2(dim,np),x3(dim,np),x4(dim,np),x5(dim,np),x6(dim,np),x7(dim,np),x8(dim,np))
    ! init all particles to 0 = null
    x1 = 0 
    x2 = 0 
    x3 = 0 
    x4 = 0 
    x5 = 0
    x6 = 0
    x7 = 0
    x8 = 0

    

    do i=1, np
      call get_containingbox(x,i,octant,origin)
      if (octant .EQ. 1) then
         x1(:,i) = x(:,i)
        !x x(:,i)
      else if (octant .EQ. 2) then
        x2(:,i) = x(:,i)
      else if (octant .EQ. 3) then
        x3(:,i) = x(:,i)
      else if (octant .EQ. 4) then 
        x4(:,i) = x(:,i)
      else if (octant .EQ. 5) then
        x5(:,i) = x(:,i)
      else if (octant .EQ. 6) then
        x6(:,i) = x(:,i)
      else if (octant .EQ. 7) then 
        x7(:,i) = x(:,i)
      else
        x8(:,i) = x(:,i)
      end if 
    end do

    ! now put these in fixed length array 


  end subroutine sort_by_octant

  subroutine resize_nodes(node)
   type(octreenode), allocatable, intent(out) :: node(:)
   type(octreenode) :: temp(size(node))
   integer :: newsize,oldsize, i


   ! supress warnings 
   oldsize = 1
   ! increase the number of possible nodes by factor of 2
   maxnodes = oldsize
   newsize = maxnodes*2
   maxnodes = newsize

   
   ! copy old data into temp array
   do i=1, oldsize
    temp(i) = node(i)
   enddo 
   write(*,*) "Copy crashing"

   ! not sure if this is needed but just to be safe
   !deallocate(node)

   ! allocate new array 
   allocate(node(newsize))

   do i=i, oldsize
    node(i) = temp(i)
   enddo 

  end subroutine resize_nodes

  RECURSIVE subroutine insert_particle(nodes,x,v,a,currentnode,currentparticle,endnode)
   ! ENDNODE AND CURRENT NODE ARE NOT THE SAME THING 
   type(octreenode), allocatable, intent(inout) :: nodes(:)
   real, intent(in) :: x(:,:), v(:,:), a(:,:)
   integer, intent(in) ::  currentparticle
   ! currentnode is the node which we are trying to insert
   ! endnode is the current latest node 
   integer, intent(out) :: currentnode, endnode
   ! 10 particles allowed in bucket 
   integer :: olddata(10)
   integer :: i, octant
   ! index of the child node 
   integer :: child, bodychildindex
   real :: origin(3), size 

   origin =0.0
   size = 0.0
   olddata = 0
   octant =0
   child = 0
   bodychildindex = 0

   write(*,*) "CurrentNode is: ", currentnode
    write(*,*) "endnode is: ", endnode
    write(*,*) "Currentparticle is: ", currentparticle
    write(*,*) "Part pos is: ", x(:,currentparticle)
    write(*,*) "Size is: ", nodes(currentnode) % size

    if (nodes(currentnode) % size  == 0 ) then
      write(*,*) "Cannot insert particle, underflow error!"
      write(*,*) "Is your box big enough?"
      stop
    endif 

   ! if we have a leaf node 
   ! and it is empty
   if (nodes(currentnode)%isLeaf .EQV. .TRUE. ) then
    !write(*,*) "First if"
   !	thedata = nodes(currentnode) % data
    !if (thedata == 0) then
    !write(*,*) nodes(currentnode) % data
    if (node_full(nodes(currentnode)) .EQV. .FALSE. ) then
     write(*,*) "Leaf node Met"
     ! set data as particle data
     ! nodes(currentnode)%data = currentparticle
     call insert_data(nodes(currentnode), currentparticle)
     return 
   ! node is a leaf but has a particle 
    
    else
      write(*,*) "Splitting part 1 "
      ! save old data so it can be reinserted later
      olddata(:) =  nodes(currentnode) % data(:)
      nodes(currentnode)% data = 0
      ! no longer a leaf node 
      nodes(currentnode)% isLeaf = .FALSE.

      ! Add particle to body children
      call insert_bodychild(nodes(currentnode),currentparticle)
      
      do i = 1, 8
        ! gen new octants 
        origin = nodes(currentnode) % origin
        size = nodes(currentnode) % size
        write(*,*) "Origin: ", origin
        write(*, *) "Size: ", size 

        ! find new origin and size for octants
        call gen_octant(origin,size,i)
        write(*,*) "Origin: ", origin
        write(*, *) "Size: ", size 
        ! store this new information
        !nodes(endnode+i) % origin = origin
        !nodes(endnode+i) % size =  size
        !write(*,*) size
        !nodes(endnode+i) % isLeaf = .TRUE.
        !nodes(endnode+i) % data = 0
        ! let the parent node know its children
        !nodes(currentnode) % children(i) = endnode + i
        !print*, nodes(currentnode) % children(i)
        call new_node(nodes(endnode+i), size, origin)
        nodes(currentnode) % children(i) = endnode + i 
      enddo
      endnode = endnode + 8 

      !write(*, *) "Octants created"

      ! work out what octant the current olddata and new data 
      ! should be in
      !
      !do i=1,8
        !write(*,*) nodes(currentnode) % children(i)
      !enddo
      ! old data 

      ! I THINK I NEED TO USE CHILD ORIGIN RATHER THAN CURRENTNODE

      ! INSERT ALL PARTICLES IN THE OLD BUCKET

      do i=1,10
        origin = nodes(currentnode) % origin
        write(*,*) "Origin is: ", origin
        call get_containingbox(x,olddata(i),octant,origin)
        write(*,*) x(:,olddata(i))
        write(*,*) x(:, currentparticle)
        write(*,*) "Delta is: ", x(:,olddata(i))-x(:,currentparticle)
        write(*,*) "Octant is: ", octant


      ! index of the child node that is being inserted
        child = nodes(currentnode) % children(octant)
        write(*,*) "Child is: ", child
       

        call insert_particle(nodes,x,v,a,child,olddata(i),endnode)
        !write(*,*) nodes(child) % data
        !write(*,*) nodes(child) % isLeaf
        !write(*,*) "Old particle finished"

      enddo
        
      ! new data
      call get_containingbox(x,currentparticle,octant,origin)
      write(*,*) "Octant is: ",  octant
      write(*,*) "Third condition reached"
      !write(*,*) "Octant is: ", octant
      ! index of the child node that is being inserted
      child = nodes(currentnode) % children(octant)
      write(*,*) "Child is 0?? ", child
      call insert_particle(nodes,x,v,a,child,currentparticle,endnode)
      write(*,*) x(:, currentparticle)
      

      write(*,*) "Insertion finished"


      
      end if
    else
      ! we are at an interior node 
      ! insert recursively until we reach leaf
      ! Should be able to insert from interior nodes 

      ! HERE IS THE PROBLEM
      origin = nodes(currentnode) % origin
      write(*,*) "Interior"
      call get_containingbox(x,currentparticle,octant,origin)
      write(*,*) "Octant:", octant
      child = nodes(currentnode) % children(octant)
      write(*,*) "CHild: ", child

      ! Add particle to body children
      call insert_bodychild(nodes(currentnode),currentparticle)
      
      call insert_particle(nodes,x,v,a,child,currentparticle,endnode) 
    end if    

  end subroutine insert_particle

  subroutine get_containingbox(x,currentparticle,octant,origin)
    real, intent(in) :: x(:,:), origin(3)
    integer, intent(in) :: currentparticle
    integer, intent(out) :: octant

    octant = 0

    ! if z is smaller than origin
    ! split into two quadtrees based on z value 
  	if (x(3,currentparticle) < origin(3)) then
     octant = 4
  	end if 

   ! This isn't needed 
  	if (x(1,currentparticle) >= origin(1) .AND. x(2,currentparticle) >= origin(2)) then 
      octant = octant + 0
  	else if (x(1,currentparticle) < origin(1) .AND. x(2,currentparticle) >= origin(2)) then 
      octant = octant  + 1
  	else if (x(1,currentparticle) < origin(1) .AND. x(2,currentparticle) < origin(2)) then 
      octant = octant + 2
  	else if (x(1,currentparticle) >= origin(1) .AND. x(2,currentparticle) < origin(2)) then 
      octant = octant + 3
  	end if

    octant = octant + 1

  end subroutine get_containingbox  



 subroutine gen_octant(origin,size,octant)
   real, intent(out) :: origin(3)
   real, intent(out) :: size
   integer, intent(in) :: octant
   real :: halfsize

   ! THIS CODE IS CORRECT 
   halfsize = size / 2.0 
   size = halfsize
   !write(*,*) size

   ! This is very ugly but I can't think of a better way
   if (octant .EQ. 1) then
    origin(1) = origin(1) + halfsize*0.5
    origin(2) = origin(2) + halfsize*0.5
    origin(3) = origin(3) + halfsize*0.5


   else if (octant .EQ. 2) then
    origin(1) = origin(1) - halfsize * 0.5
    origin(2) = origin(2) + halfsize * 0.5
    origin(3) = origin(3) + halfsize * 0.5
   else if (octant .EQ. 3) then
    origin(1) = origin(1) - halfsize * 0.5
    origin(2) = origin(2) - halfsize * 0.5
    origin(3) = origin(3) + halfsize * 0.5
   else if (octant .EQ. 4) then
    origin(1) = origin(1) + halfsize * 0.5
    origin(2) = origin(2) - halfsize * 0.5
    origin(3) = origin(3) + halfsize * 0.5
   else if (octant .EQ. 5) then
    origin(1) = origin(1) +halfsize * 0.5
    origin(2) = origin(2) + halfsize * 0.5
    origin(3) = origin(3) - halfsize * 0.5
   else if (octant .EQ. 6) then
    origin(1) = origin(1) - halfsize * 0.5
    origin(2) = origin(2) + halfsize * 0.5
    origin(3) = origin(3) - halfsize * 0.5
   else if (octant .EQ. 7) then
    origin(1) = origin(1) - halfsize * 0.5
    origin(2) = origin(2) - halfsize * 0.5
    origin(3) = origin(3) - halfsize * 0.5
   else if (octant .EQ. 8) then
    origin(1) = origin(1) + halfsize * 0.5
    origin(2) = origin(2) - halfsize * 0.5
    origin(3) = origin(3) - halfsize * 0.5
   end if
    
 end subroutine gen_octant

 recursive subroutine print_tree(nodes, x,depth,currentnode)
  type(octreenode), allocatable, intent(in) :: nodes(:)
  real, intent(in) :: x(:,:)
  integer, intent(in) :: depth, currentnode
  integer :: i, newdepth,j

  write(*,*) "Depth: ", depth
  write(*,*) "Isleaf: ",nodes(currentnode) % isLeaf 
  write(*,*) "Box Size: ", nodes(currentnode) % size
  write(*,*) "Total Mass: ", nodes(currentnode) % totalMass
  write(*,*) "Center of Mass: ", nodes(currentnode) % centerofmass
  write(*,*) "Rmax: ", nodes(currentnode) % rmax

  print*, "Current node index is: ", currentnode

  do i=1, 10
    if (nodes(currentnode)% data(i) .NE. 0 ) then 
      write(*,*) x(:,nodes(currentnode) % data(i))
    endif 
  enddo 

  print*, "Body Children: "
  do j=1, 2000
        if (nodes(currentnode) % bodychildren(j) /= 0) then

          write(*,*) nodes(currentnode) % bodychildren(j)
        endif 
      enddo 
  do i=1, 8
    if (nodes(currentnode) % children(i) .NE. 0) then
      newdepth = depth + 1
      
      call print_tree(nodes,x,newdepth,nodes(currentnode) % children(i))
    endif
  enddo 

end subroutine print_tree

subroutine cleartree(nodes)
  type(octreenode), allocatable, intent(out) :: nodes(:)
  integer :: i,sizeof

  sizeof = size(nodes)

  do i=1,sizeof
    call null_node(nodes(i))
  enddo 

end subroutine cleartree

subroutine insert_bodychild(currentnode,currentparticle)
  type(octreenode), intent(inout) :: currentnode
  integer, intent(in) :: currentparticle
  integer  :: bodychildindex

  ! INSERT 
  currentnode % bodychildpont = currentnode % bodychildpont + 1
  bodychildindex = currentnode % bodychildpont 
  currentnode % bodychildren(bodychildindex) = currentparticle 



end subroutine insert_bodychild

end module octree
  
