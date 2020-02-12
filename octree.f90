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
   integer :: i
   integer :: currentnode

   ! allocate nodes for tree
   allocate(nodes(maxnodes))

   ! iterate through all particles 
   do i =1, np
    ! if tree is full expand it 
    if (totalnodes >= maxnodes) then
       call resize_nodes(nodes)
    endif

    ! insert particle in tree
    call insert_particle(nodes,x,v,a,currentnode,i)

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

  subroutine insert_particle(nodes,x,v,a,currentnode,currentparticle)
   type(octreenode), allocatable, intent(inout) :: nodes(:)
   real, intent(in) :: x(:,:), v(:,:), a(:,:)
   integer :: currentnode, currentparticle
   integer :: olddata

   ! if we have a leaf node 
   ! and it is empty
   if (nodes(currentnode)%isLeaf .EQV. .TRUE.) then
   !	thedata = nodes(currentnode) % data
   	!if (thedata == 0) then
   if (nodes(currentnode)%data == 0) then
    ! set data as particle data
     nodes(currentnode)%data = currentparticle
   ! node is a leaf but has a particle 
   	else
      ! save old data so it can be reinserted later
      olddata =  nodes(currentnode) % data
      nodes(currentnode)% data = 0
      ! work out what octant the current olddata and new data 
      ! should be in
    end if
    end if    

  end subroutine insert_particle

 subroutine gen_octant(origin,size)
   real, intent(in) :: origin(3)
   real, intent(in) :: size
   real :: halfsize

   halfsize = size / 2.0 

    
 end subroutine gen_octant

end module octree
  
