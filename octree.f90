module octree
 use octreetype, only:octreenode
 implicit none 

 ! the opening criterion 
 real, public :: theta = 0.5

 integer:: maxnodes = 1000
 integer :: totalnodes  = 0

 public :: maketree

 contains 

  subroutine maketree(nodes,x,v,a,np)
   integer, intent(in) :: np
   type(octreenode), allocatable, intent(out) :: nodes(:)
   real, intent(in) :: x(:,:), v(:,:), a(:,:)
   integer :: i, root,j
   integer :: currentnode
   integer :: endnode

   ! allocate nodes for tree
   allocate(nodes(maxnodes))

   ! Data and children for each node = 0 
   ! i.e dosen't exist 
   do i=1, maxnodes
    nodes(i) % children(:) = 0
    nodes(i) % data = 0 
  enddo

   ! setup the root node
   nodes(1) % parent = 0
   nodes(1) % origin = (/0.0, 0.0, 0.0/)
   nodes(1) % size = 100
   nodes(1) % isLeaf = .TRUE.
   nodes(1) % data = 0 
   root = 1

   ! endnode is now root
   endnode = 1

   ! iterate through all particles 
   do i =1, np
   	write(*,*) i
    ! if tree is full expand it 
    if (endnode + 10 >= maxnodes) then
       call resize_nodes(nodes)
        do j=1, maxnodes
            ! this is wrong and breaking code
            ! go from oldsize to maxnodes
            nodes(j) % children(:) = 0
            nodes(j) % data = 0 
        enddo
    endif

    write(*,*) "Resize finished"

    ! insert particle in tree
    call insert_particle(nodes,x,v,a,root,i,endnode)
    write(*,*) " Insert finished"

   enddo     

  end subroutine maketree

  subroutine resize_nodes(node)
   type(octreenode), allocatable, intent(out) :: node(:)
   type(octreenode) :: temp(size(node))
   integer :: newsize 

   ! increase the number of possible nodes by factor of 2
   newsize = maxnodes*2
   maxnodes = newsize

   ! copy old data into temp array
   temp(:) = node(:)

   ! not sure if this is needed but just to be safe
   deallocate(node)

   ! allocate new array 
   allocate(node(newsize))

  end subroutine resize_nodes

  RECURSIVE subroutine insert_particle(nodes,x,v,a,currentnode,currentparticle,endnode)
   ! ENDNODE AND CURRENT NODE ARE NOT THE SAME THING 
   type(octreenode), allocatable, intent(inout) :: nodes(:)
   real, intent(in) :: x(:,:), v(:,:), a(:,:)
   integer, intent(in) ::  currentparticle
   ! currentnode is the node which we are trying to insert
   ! endnode is the current latest node 
   integer, intent(out) :: currentnode, endnode
   integer :: olddata
   integer :: i, octant
   ! index of the child node 
   integer :: child
   real :: origin(3), size 

    write(*,*) "CurrentNode is: ", currentnode
    write(*,*) "endnode is: ", endnode
    write(*,*) "Currentparticle is: ", currentparticle

   ! if we have a leaf node 
   ! and it is empty
   if (nodes(currentnode)%isLeaf .EQV. .TRUE. ) then
   	!write(*,*) "First if"
   !	thedata = nodes(currentnode) % data
   	!if (thedata == 0) then
    write(*,*) nodes(currentnode) % data
   if (nodes(currentnode)%data == 0) then
   	write(*,*) "Leaf node Met"
    ! set data as particle data
     nodes(currentnode)%data = currentparticle
   ! node is a leaf but has a particle 
    
   	else
   		write(*,*) "Splitting part 1 "
      ! save old data so it can be reinserted later
      olddata =  nodes(currentnode) % data
      nodes(currentnode)% data = 0
      ! no longer a leaf node 
      nodes(currentnode)% isLeaf = .FALSE.

      
      do i = 1, 8
        ! gen new octants 
        origin = nodes(currentnode) % origin
        size = nodes(currentnode) % size

      	! find new origin and size for octants
      	call gen_octant(origin,size,i)
      	! store this new information
      	nodes(endnode+i) % origin = origin
      	nodes(endnode+i) % size =  size
        write(*,*) size
      	nodes(endnode+i) % isLeaf = .TRUE.
        nodes(endnode+i) % data = 0
      	! let the parent node know its children
      	nodes(currentnode) % children(i) = endnode + i
      	!print*, nodes(currentnode) % children(i)
      enddo
      endnode = endnode + 8 

      write(*, *) "Octants created"

      ! work out what octant the current olddata and new data 
      ! should be in
      !
      do i=1,8
        !write(*,*) nodes(currentnode) % children(i)
      enddo
      ! old data 
      origin = nodes(currentnode) % origin
      call get_containingbox(x,olddata,octant,origin)
      write(*,*) x(:,olddata)
      write(*,*) octant

      ! index of the child node that is being inserted
      child =  nodes(currentnode) % children(octant)
      write(*,*) child 

      call insert_particle(nodes,x,v,a,child,olddata,endnode)
      
      write(*,*) "Old particle finished"

      ! new data
      call get_containingbox(x,currentparticle,octant,origin)
      write(*,*) "Octant is: ", octant
      ! index of the child node that is being inserted
      child = nodes(currentnode) % children(octant)
      call insert_particle(nodes,x,v,a,child,currentparticle,endnode)

      write(*,*) "Insertion finished"


      
      end if
    else
    	! we are at an interior node 
    	! insert recursively until we reach leaf
    	call get_containingbox(x,currentparticle,octant,origin)
    	child = nodes(currentnode) % children(octant)
    	call insert_particle(nodes,x,v,a,child,currentparticle,endnode) 
    end if    

  end subroutine insert_particle

  subroutine get_containingbox(x,currentparticle,octant,origin)
  	real, intent(in) :: x(:,:), origin(3)
  	integer, intent(in) :: currentparticle
  	integer, intent(out) :: octant

    octant = 1

  	! if z is smaller than origin
  	! split into two quadtrees based on z value 
  	if (x(3,currentparticle) < origin(3)) then
  		octant = 5
  	end if 

  	! This isn't needed 
  	if (x(1,currentparticle) >= origin(1) .AND. x(2,currentparticle) >= origin(2)) then 
  		octant = octant + 0
  	else if (x(1,currentparticle) < origin(1) .AND. x(2,currentparticle) >= origin(2)) then 
  		octant = octant  + 1
  	else if (x(1,currentparticle) < origin(1) .AND. x(2,currentparticle) < origin(2)) then 
  		octant = octant + 2
  	else 
  		octant = octant + 3
  	end if

  end subroutine get_containingbox  



 subroutine gen_octant(origin,size,octant)
   real, intent(out) :: origin(3)
   real, intent(out) :: size
   integer, intent(in) :: octant
   real :: halfsize

   halfsize = size / 2.0 
   size = halfsize

   ! This is very ugly but I can't think of a better way
   if (octant .EQ. 1) then
   	origin(:) = origin(:) + halfsize*0.5

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
  integer :: i, newdepth

  write(*,*) "Depth: ", depth

  if (nodes(currentnode)% data .NE. 0 ) then 
    write(*,*) x(:,nodes(currentnode) % data)
  endif 

  do i=1, 8
    if (nodes(currentnode) % children(i) .NE. 0) then
      newdepth = depth + 1
      call print_tree(nodes,x,newdepth,nodes(currentnode) % children(i))
    endif
  enddo 

end subroutine print_tree

end module octree
  
