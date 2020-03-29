module contrivedtree
 
 use octreetype
 use octree, only : gen_octant, print_tree 
 implicit none 

 integer :: maxnodes = 3
 
 contains 
  subroutine maketreecontrived(nodes,x,v,a,np)
   integer, intent(in) :: np
   type(octreenode), allocatable, intent(out) :: nodes(:)
   real, intent(in) :: x(:,:), v(:,:), a(:,:)
   integer :: i, root,j
   integer :: currentnode
   integer :: endnode, oldnode
   real :: origin(3), size

   ! allocate nodes for tree
   allocate(nodes(maxnodes))

    ! Data and children for each node = 0 
   ! i.e dosen't exist 
   do i=1, maxnodes
    call null_node(nodes(i))
  enddo

  ! endnode is now root
   endnode = 1
   ! currentnode is also root
   currentnode = 1

   ! setup the root node
   call new_node(nodes(1),1.e6,(/0.0,0.0,0.0/))

   ! Setup left and right children of root 
   
   ! right child 
   size = 1.e6
   origin = (/0.0,0.0,0.0/)
   call gen_octant(origin,size,1)
   call new_node(nodes(2),size,origin)
   
   ! left child 
   size = 1.e6 
   origin = (/0.0,0.0,0.0/)
   call gen_octant(origin,size,2)
   call new_node(nodes(3),size,origin)

   ! Now connect nodes to root node
   nodes(1) % children(1) = 2
   nodes(1) % children(2) = 3
   nodes(1) % isLeaf = .FALSE.



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


  end subroutine maketreecontrived


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
   integer :: child
   real :: origin(3), size 

   origin = 0.0
   size = 0.0
   olddata = 0
   octant =0
   child = 0

   ! Divide particles into +x and -x 

   if (x(1,currentparticle) .GE. 0) then
    ! insert into first child of the root 
     call insert_data(nodes(2), currentparticle)
   else 
    ! insert into second child of the root
     call insert_data(nodes(3), currentparticle)
   endif

  end subroutine insert_particle 

end module contrivedtree