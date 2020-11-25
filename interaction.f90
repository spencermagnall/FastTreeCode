module interaction
 use octree
 use opening_criterion
 use taylor_expansions
 use potendirect 
 ! This is just for testing remove later
 use evaluate  
 use timing

 implicit none  

 type interact_global_stack
  integer :: nodeindex1
  integer :: interactions(8) = 0 

 end type 

 type interact_stack_data

  integer :: nodeindex1,nodeindex2 

 end type 

 contains 


 subroutine interact_nolock(nodes,x,m,a,nopart)
  use omp_lib
  !$ use omputils
  integer, intent(in) :: nopart
  type(octreenode), intent(inout) :: nodes(:)
  real, intent(in) :: x(3,nopart),m(nopart)
  real, intent(inout) :: a(3,nopart)
  integer :: stacksize, top, iter,nodeindex1,nodeindex2
  type(interact_stack_data) :: stacklocal(1000)
  type(interact_global_stack) :: stack(1000) 
  type(octreenode) :: newnode1,newnode2, splitnode, regnode
  integer :: i, j, counter,k
  real :: rmax1, rmax2,cm1(3),cm2(3)
  real :: fnode(20),quads(6)
  real :: dr,dx(3),dy,dz,totmass,r,r2
  integer :: particleindex(10),particleindex2(10) 
  real :: c0,c1(3),c2(3,3),c3(3,3,3), c0new,c1new(3),c2new(3,3),c3new(3,3,3)
  logical :: nodesAreEqual, flag
  integer :: regnodechild(2000), splitnodeindex,regnodeindex,newnodeindex1,newnodeindex2,start,numthreads
  integer, allocatable :: particlesforsum(:),istacklocal(:)
  logical :: node1empty, node2empty, updatestackpont
  logical, allocatable :: threadworking(:)
  integer, allocatable :: nodeforthread(:)



  quads(:) = 0.0
  counter = 0

  fnode(:) = 0.0
  particleindex = 0.0
  particleindex2 = 0.0

  node1empty = .true.
  node2empty = .true.

  regnodechild(:) = 0.
  c0 = 0.
  c1 = 0.
  c2 = 0.
  c3 = 0.
  c0new = 0.
  c1new = 0.
  c2new = 0.
  c3new = 0.

  ! Push the root node to the top of the stack 
  top = 1
  stack(top) % nodeindex1 = 1
  stack(top) % interactions(1) = 1

   print*, "top is: ", top 
   print*, "Before openmp loop"
   numthreads = 1
   ! get number of OpenMPthreads
   !$omp parallel default(none) shared(numthreads)
     numthreads = omp_get_num_threads()
   !$omp end parallel
   allocate(threadworking(numthreads))
   allocate(istacklocal(numthreads))
   allocate(nodeforthread(numthreads))
   threadworking = .true.

   print*, "Number of threads is: ", numthreads
   istacklocal = 0
   !stacklocal = 0
   nodeforthread = 0

   ! Initialise omp locks 
   !$ call init_omp()

   print*,"Wall time start: ", wallclock()

   !$omp parallel default(none) &
   !$omp shared(stack,top,x,a,m,nopart,nodes,numthreads,threadworking,istacklocal) &
   !$omp shared(ipart_omp_lock) &
   !$omp private(c0,c1,c2,c3,c0new,c1new,c2new,c3new) &
   !$omp private(splitnode,dx,k,nodeindex1,nodeindex2) &
   !$omp private(node1empty,node2empty,particleindex,start) &
   !$omp private(counter,newnodeindex1,newnodeindex2,cm1,cm2,dr) &
   !$omp private(fnode,totmass,r2,r,quads,rmax1,rmax2,splitnodeindex) &
   !$omp private(regnode,regnodeindex,regnodechild,particleindex2,stacklocal,nodeforthread,updatestackpont)

   do while (any(istacklocal > 0) .or. top > 0)

    !$ k=omp_get_thread_num() + 1

    !print*,"Stack local",istacklocal
    !print*, "K: ",k

    if (istacklocal(k) > 0 ) then ! pop off local stack
        !!$omp critical 
        !print*,"istacklocal", istacklocal(k)
        nodeindex1 = stacklocal(istacklocal(k)) % nodeindex1
        nodeindex2 = stacklocal(istacklocal(k)) % nodeindex2
        !print*, "Node indexes : ", nodeindex1, nodeindex2

        ! If the nodes are free then we can do work with them 
        !if (nodes(nodeindex1) % nodefree .and. nodes(nodeindex2) % nodefree) then 
          istacklocal(k) = istacklocal(k) - 1
          threadworking(k) = .true.
          nodes(nodeindex1) % nodefree  = .false.
          nodes(nodeindex2) % nodefree  = .false.
        !else 
        !  threadworking(k) = .false.
        !endif 
        !!$omp critical  

    else 
      !$omp critical (stack)
      if (top > 0) then 
        !!$OMP TASK 
       ! nodeindex1 = stack(top) % nodeindex1
        !nodeindex2 = stack(top) % nodeindex2
          ! If the nodes are free then we can do work with them 
        !if (nodes(nodeindex1) % nodefree .and. nodes(nodeindex2) % nodefree) then 
          !top = top - 1
          call global_to_local(stack,stacklocal,k,istacklocal,top)
          ! pop some work of the local stack 
          nodeindex1 = stacklocal(istacklocal(k)) % nodeindex1
          nodeindex2 = stacklocal(istacklocal(k)) % nodeindex2
          istacklocal(k) = istacklocal(k) - 1
          threadworking(k) = .true.
          !nodes(nodeindex1) % nodefree  = .false.
          !nodes(nodeindex2) % nodefree  = .false.
      else 
        threadworking(k) = .false.
      endif
      !$omp end critical (stack)
    endif 

    if (threadworking(k)) then 
    rmax1 = nodes(nodeindex1) % rmax
    rmax2 = nodes(nodeindex2) % rmax 
    quads(:) = 0.0
    counter = 0

    fnode(:) = 0.0
    particleindex = 0.0
    particleindex2 = 0.0
  !particleindex1 = 0.

    node1empty = .true.
    node2empty = .true.

    regnodechild(:) = 0.
    c0 = 0.
    c1 = 0.
    c2 = 0.
    c3 = 0.
    c0new = 0.
    c1new = 0.
    c2new = 0.
    c3new = 0.
    updatestackpont = .true. 
    !print*, "nodeindex1: ", nodeindex1
    !print*, "nodeindex2: ", nodeindex2

   

    
    ! START NODES ARE EQUAL SECTION 

    if (nodeindex1 == nodeindex2) then
      ! IF NODES ARE LEAFS PERFORM DIRECT SUM ON NODE 
      if (nodes(nodeindex1) % isLeaf) then 
         do i=1, 10
          if (nodes(nodeindex1) % data(i) /= 0) then
              node1empty = .false.
              particleindex(i) = nodes(nodeindex1) % data(i)
          endif  
          enddo 
          ! CALL DIRECT SUM 
          !print*,"Direct sum"
          !print*, x
          if (.not. node1empty) then 
            !!$omp critical 
            call get_accel(x,a,m,nopart,particleindex)
             !print*, "nodeindex1: ", nodeindex1
             !print*, "nodeindex2: ", nodeindex2
            !!$omp end critical 
          endif 

      ! OTHERWISE SPLIT NODES 
      else
        !print*, "Splitting"
        do i=1,8
          ! THREAD HAS TO TAKE ALL INTERACTIONS FOR THAT NODE 
          newnodeindex1 = nodes(nodeindex1) % children(i)
          !$omp critical (stack)
          top = top + 1
          !!$omp end critical (stack)

          do j=1,8
            if (nodes(nodeindex1) % children(i) .NE. 0 .AND. nodes(nodeindex2) % children(j) .NE. 0) then 
              
              ! Get the sub-cells
              !newnodeindex1 = nodes(nodeindex1) % children(i)
              newnodeindex2 = nodes(nodeindex2) % children(j)

              !if (istacklocal(k) < 3) then 
              !if (all(threadworking)) then 
              !istacklocal(k) = istacklocal(k) + 1
              !stacklocal(istacklocal(k)) % nodeindex1 = newnodeindex1
              !stacklocal(istacklocal(k)) % nodeindex2 = newnodeindex2
              
              !else
              ! Push new data to stack 
                !print*, "top", top 
                stack(top) % nodeindex1 = newnodeindex1
                stack(top) % interactions(j) = newnodeindex2
            endif
          enddo 
          !$omp end critical (stack)
        enddo 

        !print*, "stack: ", stack(1:top)
        !stop 

      endif

    ! END NODES ARE EQUAL SECTION 


    ! ARE NODES WELL SEPARATED? 
    elseif (well_separated(nodes(nodeindex1),nodes(nodeindex2))) then
      ! SYMMETRY CHECK 
     if ( .not. well_separated(nodes(nodeindex2),nodes(nodeindex1))) then 
        print*, "MAC BROKEN"
        stop
     endif  
   

  

      ! WELL SEPARATED NODE MI IS CALCULATED: TAYLOR COEFFs computed and added to 
      ! node data fields

    if (nodes(nodeindex1) % totalmass /= 0. .AND. nodes(nodeindex2) %totalmass /= 0.) then !&
     !.and. (.NOT. nodes(nodeindex1) % isLeaf) .and. (.not. nodes(nodeindex2) % isLeaf)) then
      cm1 = nodes(nodeindex1) % centerofmass
      cm2 = nodes(nodeindex2) % centerofmass


      call get_dx_dr(cm1,cm2,dx,dr)

      fnode = 0.0

      totmass = nodes(nodeindex2) % totalmass

      c0 = 0.
      c1 = 0.
      c2 = 0.
      c3 = 0.
      quads = 0.
      !quads = nodes(nodeindex2) % quads 
      call compute_coeff(dx(1),dx(2),dx(3),dr,totmass,quads,c0,c1,c2,c3)
      !nodes(nodeindex1) % fnode = fnode 

      !print*, "Crashing on lock section "
      !print*, "Size: ", ipart_omp_lock(nodeindex1)
      !$ call omp_set_lock(ipart_omp_lock(nodeindex1))
      !print*, "Entered locking"
      nodes(nodeindex1) % c0 =  nodes(nodeindex1)%c0 + c0 
      nodes(nodeindex1) % c1 = nodes(nodeindex1)%c1 + c1 
      nodes(nodeindex1) % c2 = nodes(nodeindex1) % c2 + c2 
      nodes(nodeindex1) % c3 = nodes(nodeindex1) % c3 + c3 
      !$ call omp_unset_lock(ipart_omp_lock(nodeindex1))


    endif 

    ! END NODES WELL SEPARATED SECTION


    ! NODES ARE BEING SPLIT 
    else 
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
    !!$omp critical (stack)
    !print*, "Splitting 2 "
    do i=1,8
      !print*, i
      !print*, splitnode % children(i)
      !!$omp critical (stack)
      if (splitnode % children(i) /= 0 .and. nodes(splitnode % children(i)) % totalmass /= 0.) then 
        splitnodeindex = splitnode % children(i)
        !print*, "Splitnode index: ", splitnodeindex
        !print*,"regnodeindex: ", regnodeindex
        ! CHECK FOR WHICH NODE IS NODE 1
        if (regnodeindex == nodeindex1) then 

            ! NODE 1 IS NOT SPLIT SO PUSH TO GLOBAL STACK 
            istacklocal(k) = istacklocal(k) + 1
            stacklocal(istacklocal(k)) % nodeindex1 = regnodeindex
            stacklocal(istacklocal(k)) % nodeindex2 = splitnodeindex
        else 
            !$omp critical (stack)
            !if (updatestackpont) then
              top = top + 1
              updatestackpont = .false.
              stack(top) % interactions = 0
            !endif
            !print*, "top: ", top 
            !print*, "i: ",i  
            stack(top) % nodeindex1 = splitnodeindex
            stack(top) % interactions(i) = regnodeindex 
            !print*, "These values pushed: ", stack(top) % nodeindex1
            !print*, stack(top) % interactions
            !print*, "top is: ", top 
            !$omp end critical (stack) 
            !print*,"stack: ", stack(1:top)
            !stop 
        endif 
      endif  
    enddo
    !!$omp end critical (stack)

    ! LEAF-NODE NODE 
    elseif(splitnode % isLeaf .AND. .NOT. regnode % isLeaf) then
      !print*, "FIX THIS"

      ! Call poten function
      ! Need a way to get all of the children of a node for direct sum

      ! Get children of regular node 
      node1empty = .true.
      node2empty = .true.

      do i=1,regnode % bodychildpont
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

      !print*, node1empty, node2empty


      if (.NOT. node1empty .and. .not. node2empty) then
        if (splitnodeindex == nodeindex1) then 
          call get_accel_leafnode(x,a,m,nopart,particleindex,regnodechild)
        else
          call get_accel_leafnode(x,a,m,nopart,regnodechild,particleindex)
        endif 
      endif 
     !call get_accel_leafnode(x,a,m,np,particleindex,regnodechild)

      ! Get children of a leaf node 

    ! LEAF-NODE LEAF-NODE
    else 
      !print*, "Leafnode-Leafnode"

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

     !print*, "Finished getting bodies"
     !print*, node1empty
     !print*, node2empty
     ! CALL DIRECT SUM 
    if (.NOT. node1empty .and. .not. node2empty) then
      call get_accel_leafnode(x,a,m,nopart,particleindex,particleindex2)
    endif 


    endif 
    
    endif  

    endif 

   enddo 

   !$omp end parallel 

   print*,"Wall time : ", wallclock()




 end subroutine interact_nolock

 subroutine interact_yokota(nodes,x,m,a,nopart)
  use omp_lib
  integer, intent(in) :: nopart
  type(octreenode), intent(inout) :: nodes(:)
  real, intent(in) :: x(3,nopart),m(nopart)
  real, intent(inout) :: a(3,nopart)
  integer :: stacksize, top, iter,nodeindex1,nodeindex2
  type(interact_stack_data) :: stack(100000), stacklocal(100000)
  type(octreenode) :: newnode1,newnode2, splitnode, regnode
  integer :: i, j, counter,k
  real :: rmax1, rmax2,cm1(3),cm2(3)
  real :: fnode(20),quads(6)
  real :: dr,dx(3),dy,dz,totmass,r,r2
  integer :: particleindex(10),particleindex2(10) 
  real :: c0,c1(3),c2(3,3),c3(3,3,3), c0new,c1new(3),c2new(3,3),c3new(3,3,3)
  logical :: nodesAreEqual, flag
  integer :: regnodechild(100000), splitnodeindex,regnodeindex,newnodeindex1,newnodeindex2,start,numthreads
  integer, allocatable :: particlesforsum(:),istacklocal(:)
  logical :: node1empty, node2empty
  logical, allocatable :: threadworking(:)



  quads(:) = 0.0
  counter = 0

  fnode(:) = 0.0
  particleindex = 0.0
  particleindex2 = 0.0

  node1empty = .true.
  node2empty = .true.

  regnodechild(:) = 0.
  c0 = 0.
  c1 = 0.
  c2 = 0.
  c3 = 0.
  c0new = 0.
  c1new = 0.
  c2new = 0.
  c3new = 0.

  ! Push the root node to the top of the stack 
  stack(top)% nodeindex1 = 1
  stack(top)% nodeindex2 = 1

   print*, "top is: ", top 
   print*, "Before openmp loop"
   numthreads = 1
   ! get number of OpenMPthreads
   !$omp parallel default(none) shared(numthreads)
     numthreads = omp_get_num_threads()
   !$omp end parallel
   allocate(threadworking(numthreads))
   allocate(istacklocal(numthreads))
   threadworking = .true.

   print*, "Number of threads is: ", numthreads
   istacklocal = 0

   !$omp parallel default(shared) &
   !$omp shared(top,nodes,x,a,m)
   !$omp single 
    do while (top > 0)
      ! Migrate to lock free stack at some point 
      !$omp critical 
        nodeindex1 = stack(top) % nodeindex1
        nodeindex2 = stack(top) % nodeindex2
        top = top - 1
      !$omp end critical 
      ! NODES ARE EQUAL but not leaf nodes then split into MI between children (i.e root)
      if (nodeindex1 == nodeindex2 .and. nodes(nodeindex1) % isLeaf .eqv. .false.) then 
        ! AT MOST 36 Interactions 
      ! 64 in this case but some will be duplicates ???
        start = 1
        do i=start,8
          do j=start,8

        ! Get the index of the sub-cells
         !print*, "i, j: "
         !print*, i
         !print*, j
         !print*, nodes(nodeindex1) % children(i)
         !print*, nodes(nodeindex2) % children(j)
        ! If children exist
        if (nodes(nodeindex1) % children(i) .NE. 0 .AND. nodes(nodeindex2) % children(j) .NE. 0) then 
         !nodeindex1 = nodes(nodeindex1) % children(i)
         !nodeindex2 = nodes(nodeindex2) % children(j)
         !print*, nodeindex1
         !print*, nodeindex2
         !counter = counter + 1
         !print*,"Counter: ", counter

          ! Get the sub-cells
          !newnode1 = nodes(nodeindex1)
          !newnode2 = nodes(nodeindex2)
          newnodeindex1 = nodes(nodeindex1) % children(i)
          newnodeindex2 = nodes(nodeindex2) % children(j)

          !print*, "Nodeindex1: ", newnodeindex1, "Nodeindex2: ", newnodeindex2

          !print*, "Interact called"
          ! call the interact for each of the sub-cells
          !call interact(newnodeindex1, newnodeindex2, nodes,x,m,a,nopart)
          ! If local stack is sufficently full push to global stack 
          !$omp critical 
          top = top + 1
          stack(top) % nodeindex1 = newnodeindex1
          stack(top) % nodeindex2 = newnodeindex2 
          !$omp end critical 

        endif 

          enddo 
          start = start + 1
        enddo

      endif 

      enddo 
   !$omp end single 
   !$omp end parallel 


end subroutine interact_yokota

 subroutine interact_stack(nodes,x,m,a,nopart)
  use omp_lib
  !$ use omputils
  integer, intent(in) :: nopart
  type(octreenode), intent(inout) :: nodes(:)
  real, intent(in) :: x(3,nopart),m(nopart)
  real, intent(inout) :: a(3,nopart)
  integer :: stacksize, top, iter,nodeindex1,nodeindex2
  type(interact_stack_data) :: stack(10000), stacklocal(10000)
  type(octreenode) :: newnode1,newnode2, splitnode, regnode
  integer :: i, j, counter,k
  real :: rmax1, rmax2,cm1(3),cm2(3)
  real :: fnode(20),quads(6)
  real :: dr,dx(3),dy,dz,totmass,r,r2
  integer :: particleindex(10),particleindex2(10) 
  real :: c0,c1(3),c2(3,3),c3(3,3,3), c0new,c1new(3),c2new(3,3),c3new(3,3,3)
  logical :: nodesAreEqual, flag
  integer :: regnodechild(100000), splitnodeindex,regnodeindex,newnodeindex1,newnodeindex2,start,numthreads
  integer, allocatable :: particlesforsum(:),istacklocal(:)
  logical :: node1empty, node2empty
  logical, allocatable :: threadworking(:)
  real :: threadcounts, idlecounts 
  

  top = 1
  iter = 1

  ! Push the root node to the top of the stack 
  stack(top)% nodeindex1 = 1
  stack(top)% nodeindex2 = 1

   print*, "top is: ", top 
   print*, "Before openmp loop"
   numthreads = 1
   ! get number of OpenMPthreads
   !$omp parallel default(none) shared(numthreads)
     numthreads = omp_get_num_threads()
   !$omp end parallel
   allocate(threadworking(numthreads))
   allocate(istacklocal(numthreads))
   threadworking = .true.

   print*, "Number of threads is: ", numthreads
   istacklocal = 0

    ! Initialise omp locks 
   !$ call init_omp()

   ! Want to start timing from the parallel loop 

  print*,"Wall time start: ", wallclock()
  !$omp parallel default(none) &
  !!$omp threadprivate()
  !$omp shared(ipart_omp_lock) &
  !$omp shared(stack,top,x,a,m,nopart,nodes,numthreads,threadworking,istacklocal) &
  !$omp private(c0,c1,c2,c3,c0new,c1new,c2new,c3new) &
  !$omp private(splitnode,dx,k,nodeindex1,nodeindex2) &
  !$omp private(node1empty,node2empty,particleindex,start) &
  !$omp private(counter,newnodeindex1,newnodeindex2,cm1,cm2,dr) &
  !$omp private(fnode,totmass,r2,r,quads,rmax1,rmax2,splitnodeindex) &
  !$omp private(regnode,regnodeindex,regnodechild,particleindex2,stacklocal) &
  !$omp private(threadcounts,idlecounts)

  !!$omp single 

  !quads(:) = 0.0
  !counter = 0

  !fnode(:) = 0.0
  !particleindex = 0.0
  !particleindex2 = 0.0

  !node1empty = .true.
  !node2empty = .true.

  !regnodechild(:) = 0.
  !c0 = 0.
  !c1 = 0.
  !c2 = 0.
  !c3 = 0.
  !c0new = 0.
  !c1new = 0.
  !c2new = 0.
  !c3new = 0.
  !rmax1 = 0.
  !rmax2 = 0.
  !nodeindex1 = 0
  !nodeindex2 = 0

  threadcounts = 0
  idlecounts = 0

  do while (any(istacklocal > 0) .or. top > 0)
    threadcounts = threadcounts + 1 
  !do while (top > 0)
    !print*, "Current itteration is: ",iter
    !iter = iter + 1

    !$ k=omp_get_thread_num() + 1
    !print*,"Top",top
    !print*,"Stack top:", istacklocal(k)


    ! POP ITEM FROM STACK
    !!$omp critical 
   if (istacklocal(k) > 0 ) then ! pop off local stack
        !!$omp critical 
        nodeindex1 = stacklocal(istacklocal(k)) % nodeindex1
        nodeindex2 = stacklocal(istacklocal(k)) % nodeindex2

        ! If the nodes are free then we can do work with them 
        !if (nodes(nodeindex1) % nodefree .and. nodes(nodeindex2) % nodefree) then 
          istacklocal(k) = istacklocal(k) - 1
          threadworking(k) = .true.
          nodes(nodeindex1) % nodefree  = .false.
          nodes(nodeindex2) % nodefree  = .false.
        !else 
        !  threadworking(k) = .false.
        !endif 
        !!$omp critical  

    else 
      !$omp critical (stack)
      if (top > 0) then 
        !!$OMP TASK 
        nodeindex1 = stack(top) % nodeindex1
        nodeindex2 = stack(top) % nodeindex2
          ! If the nodes are free then we can do work with them 
        !if (nodes(nodeindex1) % nodefree .and. nodes(nodeindex2) % nodefree) then 
          top = top - 1
          threadworking(k) = .true.
          nodes(nodeindex1) % nodefree  = .false.
          nodes(nodeindex2) % nodefree  = .false.
          
        !endif 
        
        !!$OMP END TASK 
        !else
        !  threadworking(k) = .false.
        !endif
      else 
        threadworking(k) = .false.
        idlecounts = idlecounts + 1
      endif
      !$omp end critical (stack)
    endif 
    !!$omp end critical 

     !print*, "Thread number", k
     !print*, "Working? :", threadworking(k)
     

    ! START OF INTERACTION 
    !if (nodeindex1 /= 0 .and. nodeindex2 /= 0) then
    !print*, "Node index 1: ", nodeindex1
    !print*,"Node index 2: ", nodeindex2 
    !print*,"Working (1) (2): ", nodes(nodeindex1) % nodefree, nodes(nodeindex2) % nodefree 
    !endif

    if (threadworking(k)) then 
     rmax1 = nodes(nodeindex1) % rmax
     rmax2 = nodes(nodeindex2) % rmax 
     quads(:) = nodes(nodeindex2) % quads
     counter = 0

     fnode(:) = 0.0
     particleindex = 0.0
     particleindex2 = 0.0
  !particleindex1 = 0.

     node1empty = .true.
     node2empty = .true.

     regnodechild(:) = 0.
     c0 = 0.
     c1 = 0.
     c2 = 0.
     c3 = 0.
     c0new = 0.
     c1new = 0.
     c2new = 0.
     c3new = 0.
    ! if nodes are equal 
    if (nodeindex1 == nodeindex2) then 
      ! if nodes are leaf node 
      if (nodes(nodeindex1)%isLeaf) then 
         do i=1, 10
          if (nodes(nodeindex1) % data(i) /= 0) then
              node1empty = .false.
              particleindex(i) = nodes(nodeindex1) % data(i)
          endif  
          enddo 
        ! CALL DIRECT SUM 
        !print*,"Direct sum"
        !print*, x
        if (.not. node1empty) then 
          !!$omp critical 
          call get_accel(x,a,m,nopart,particleindex)
          !!$omp end critical 
          !UNLOCK NODES 
          !print*, "Nodes unlocked"
          !nodes(nodeindex1) % nodefree = .true.
          !nodes(nodeindex2) % nodefree = .true. 

          !!$omp end critical 
        endif 
        !nodes(nodeindex1) % nodefree = .true.


        ! CELL SELF INTERACTION IS SPLIT INTO MI BETWEEN SUBNODES
      else
      ! AT MOST 36 Interactions 
      ! 64 in this case but some will be duplicates ???
        start = 1
        do i=start,8
          do j=start,8

        ! Get the index of the sub-cells
         !print*, "i, j: "
         !print*, i
         !print*, j
         !print*, nodes(nodeindex1) % children(i)
         !print*, nodes(nodeindex2) % children(j)
        ! If children exist
        if (nodes(nodeindex1) % children(i) .NE. 0 .AND. nodes(nodeindex2) % children(j) .NE. 0) then 
         !nodeindex1 = nodes(nodeindex1) % children(i)
         !nodeindex2 = nodes(nodeindex2) % children(j)
         !print*, nodeindex1
         !print*, nodeindex2
         counter = counter + 1
         !print*,"Counter: ", counter

          ! Get the sub-cells
          !newnode1 = nodes(nodeindex1)
          !newnode2 = nodes(nodeindex2)
          newnodeindex1 = nodes(nodeindex1) % children(i)
          newnodeindex2 = nodes(nodeindex2) % children(j)

          !print*, "Nodeindex1: ", newnodeindex1, "Nodeindex2: ", newnodeindex2

          !print*, "Interact called"
          ! call the interact for each of the sub-cells
          !call interact(newnodeindex1, newnodeindex2, nodes,x,m,a,nopart)


          ! PUSH NEW INTERACTIONS TO STACK 
          ! Push to local stack 
          ! Evenly divide new work between threads 
          if (istacklocal(k) < int(36/numthreads)) then 
          !if (all(threadworking)) then 
            istacklocal(k) = istacklocal(k) + 1
            stacklocal(istacklocal(k)) % nodeindex1 = newnodeindex1
            stacklocal(istacklocal(k)) % nodeindex2 = newnodeindex2 
          else
          ! If local stack is sufficently full push to global stack 
          !$omp critical (stack)
          top = top + 1
          stack(top) % nodeindex1 = newnodeindex1
          stack(top) % nodeindex2 = newnodeindex2 
          !$omp end critical (stack)
          endif 

        endif 

          enddo 
          start = start + 1
        enddo
        ! UNLOCK THE NODES 
        !print*,"Nodes unlocked"
        !nodes(nodeindex1) % nodefree  = .true.
        !nodes(nodeindex2) % nodefree  = .true.

    endif

   

  elseif (well_separated(nodes(nodeindex1),nodes(nodeindex2))) then 
  !elseif (.false.) then 
  !elseif (.true.) then
    !print*, "Well separated!"

    ! SYMMETRY CHECK 
     if ( .not. well_separated(nodes(nodeindex2),nodes(nodeindex1))) then 
        print*, "MAC BROKEN"
        stop
     endif  

  

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

    !print*, "particles in node 1"
    !print*, nodes(nodeindex1) % data
    !print*, "particles in node 2"
    !print*, nodes(nodeindex2) % data
    print*, "nodeindex1: ", nodeindex1
    print*, "nodeindex2: ", nodeindex2

    call get_dx_dr(cm1,cm2,dx,dr)
    !print*, "Cm of node1: ",cm1
    !print*, "Cm of node2: ",cm2
    !dx = cm1(1) - cm2(1)
    !print*, "Dx: ", dx
    !dy = cm1(2) - cm2(2)
    !print*, "Dy: ", dy
    !dz = cm1(3) - cm2(3)
    !print*, "Dz: ", dz
    !dr = 1./(sqrt(dot_product(cm1-cm2,cm1-cm2)))
    !print*, "Dr: ", dr 

    fnode = 0.0

    totmass = nodes(nodeindex2) % totalmass
    !print*, "Totalmass: ",totmass
    !totmass = 1.0
    print*, "dx: ", dx
    print*, "dr: ", dr 
    ! calculate real accel here 
    !print*, "Real accel: "
     ! r2 = dot_product(dx,dx)
     ! r  = 1./sqrt(r2)
     !print*, -totmass*(1/((r2)**1.5))*dx 
     !print*, "Real force: "
     !print*, -totmass*(1/((r2)**1.5))*dx * nodes(nodeindex1) %totalmass

    c0 = 0.
    c1 = 0.
    c2 = 0.
    c3 = 0.
    quads = 0.
    quads = nodes(nodeindex2) % quads 
    call compute_fnode(dx(1),dx(2),dx(3),dr,totmass,quads,fnode)
    !call compute_coeff(dx(1),dx(2),dx(3),dr,totmass,quads,c0,c1,c2,c3)

    !!$omp critical (node)
    

    ! store coeff for walk phase 
    !node1 % fnode = node1 % fnode + fnode
    !print*, "Stored acccel: "
    !$ call omp_set_lock(ipart_omp_lock(nodeindex1))
    !nodes(nodeindex1) % c0 =  nodes(nodeindex1)%c0 + c0 
    nodes(nodeindex1) % fnode =nodes(nodeindex1) % fnode + fnode 
    print*, "fnode(1:20) 1: ", fnode(1:20)* nodes(nodeindex1) % totalmass
    ! print*, "fnode: ", nodes(nodeindex1) % fnode
    !print*, nodes(nodeindex1) % c1
    !print*, "calc accel: "
    !print*,c1
    !nodes(nodeindex1) % c1 = nodes(nodeindex1)%c1 + c1 
    !nodes(nodeindex1) % c2 = nodes(nodeindex1) % c2 + c2 
    !nodes(nodeindex1) % c3 = nodes(nodeindex1) % c3 + c3 
    !$ call omp_unset_lock(ipart_omp_lock(nodeindex1))
    !!$omp end critical (node)


    call get_dx_dr(cm2,cm1,dx,dr)

    print*, "dx: ", dx
    print*, "dr: ", dr 

    fnode = 0.0

    totmass = nodes(nodeindex1) % totalmass
    !totmass = 1
    quads = nodes(nodeindex1) % quads
    !quads = 0.  
    call compute_fnode(dx(1),dx(2),dx(3),dr,totmass,quads,fnode)
    !$ call omp_set_lock(ipart_omp_lock(nodeindex1))
    print*, "fnode(1:20) 2: ", fnode(1:20)*nodes(nodeindex2) % totalmass
    nodes(nodeindex2) % fnode = nodes(nodeindex2) % fnode + fnode 
    
    !$ call omp_unset_lock(ipart_omp_lock(nodeindex1))

    ! ! Calculate "force"
    ! c0 = c0 * nodes(nodeindex1) % totalmass
    ! c1 = c1 * nodes(nodeindex1) % totalmass 
    ! c2 = c2 * nodes(nodeindex1) % totalmass
    ! c3 = c3 * nodes(nodeindex1) % totalmass


    ! ! force is equal and opposite 
    ! c0new = -c0
    ! c1new = -c1
    ! c2new = -c2 
    ! c3new = -c3  

    ! ! divide by node mass to get 

    ! c0 = c0new / (nodes(nodeindex2) % totalmass)
    ! c1 = c1new / (nodes(nodeindex2) % totalmass)
    ! c2 = c2new / (nodes(nodeindex2) % totalmass)
    ! c3 = c3new / (nodes(nodeindex2) % totalmass)

    ! print poten
    !print*, "Poten is: ", fnode(20)
    !call poten_at_bodypos(cm1,cm2,c0,c1,c2,c3,poten(20))
    !print*, "Poten is: ", poten(20)

    !return 

    !open(unit=77,file="wellseperated.txt")
    !write(77,*), "Node 1:", nodeindex1, "Node 2: ", nodeindex2

    !cm2 = nodes(nodeindex1) % centerofmass
    !cm1 = nodes(nodeindex2) % centerofmass

    !print*, "particles in node 1"
    !print*, nodes(nodeindex2) % data
    !print*, "particles in node 2"
    !print*, nodes(nodeindex1) % data


    !call get_dx_dr(cm1,cm2,dx,dr)
    !print*, "Cm of node1: ",cm1
    !print*, "Cm of node2: ",cm2
    !dx = cm1(1) - cm2(1)
    !print*, "Dx: ", dx
    !dy = cm1(2) - cm2(2)
    !print*, "Dy: ", dy
    !dz = cm1(3) - cm2(3)
    !print*, "Dz: ", dz
    !dr = 1./(sqrt(dot_product(cm1-cm2,cm1-cm2)))
    !print*, "Dr: ", dr 

    !fnode = 0.0

    !totmass = nodes(nodeindex1) % totalmass
    !print*, "Totalmass: ",totmass
    !totmass = 1.0

    ! calculate real accel here 
    !print*, "Real accel: "
     !r2 = dot_product(dx,dx)
     !r  = sqrt(r2)
     !print*, -totmass*(1/((r2)**1.5))*dx 
     !print*, "Real force: "
     !print*, -totmass*(1/((r2)**1.5))*dx * nodes(nodeindex2) %totalmass

    !c0new = c0
    !c1new = c1
    !c2new = c2 
    !c3new = c3 
    !c0 = 0.
    !c1 = 0.
    !c2 = 0.
    !c3 = 0.
    !quads = 0.
    !quads = nodes(nodeindex1) % quads 


    !call compute_fnode(dx(1),dx(2),dx(3),dr,totmass,quads,fnode)
    !call compute_coeff(dx(1),dx(2),dx(3),dr,totmass,quads,c0,c1,c2,c3)

    !!$omp critical (node)

    ! store coeff for walk phase 
    ! !$ call omp_set_lock(ipart_omp_lock(nodeindex2))
    ! !nodes(nodeindex2) % fnode = nodes(nodeindex2) % fnode + fnode
    ! !print*, "Stored acccel: "
    ! nodes(nodeindex2) % c0 =  nodes(nodeindex2)%c0 + c0
    ! !print*, nodes(nodeindex2) % c1
    ! !print*, "calc accel: "
    ! !print*,c1
    ! nodes(nodeindex2) % c1 = nodes(nodeindex2)%c1 + c1
    ! nodes(nodeindex2) % c2 = nodes(nodeindex2) % c2 + c2
    ! nodes(nodeindex2) % c3 = nodes(nodeindex2) % c3 + c3
    ! !$ call omp_unset_lock(ipart_omp_lock(nodeindex2))
    !!$omp end critical (node)

    ! print poten
    !print*, "Poten is: ", fnode(20)
    !call poten_at_bodypos(cm1,cm2,c0,c1,c2,c3,poten(20))
    !print*, "Poten is: ", poten(20)

    !return 
    !STOP

    !open(unit=77,file="wellseperated.txt")
    !write(77,*), "Node 1:", nodeindex2, "Node 2: ", nodeindex1

    !print*,"Values for interaciton"
    !print*, c1 
    !print*,abs(c0 * nodes(nodeindex2) % totalmass) -  abs(c0new * totmass)
    !print*, abs(c1 *nodes(nodeindex2) % totalmass) - abs(c1new * totmass)
    !print*,abs(c2 * nodes(nodeindex2) % totalmass) - abs(c2new * totmass)
    !print*, abs(c3 * nodes(nodeindex2) % totalmass) - abs(c3new * totmass)

    !STOP

    ! UNLOCK THE NODES 
    !print*,"Nodes unlocked"
    

  endif 
  ! UNLOCK THE NODES 
  !print*,"Nodes unlocked"
  !nodes(nodeindex1) % nodefree  = .true.
  !nodes(nodeindex2) % nodefree  = .true.

  ! THE NODE WITH THE LARGER RMAX IS SPLIT; up to 8 new MI are created and processed 
  else
    !print*, "Split nodes"

    !print*,"rmax1", rmax1
    !print*, "rmax2",rmax2

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
      !print*, i
      !print*, splitnode % children(i)
      if (splitnode % children(i) /= 0 .and. nodes(splitnode % children(i)) % totalmass /= 0.) then 
        splitnodeindex = splitnode % children(i)
        !newnode1 = nodes(nodeindex1)
        !print*, "Internal splitnode case"
        !print*, "regnode index: ", regnodeindex
        !print*, "Split node index: ",splitnodeindex

        if (istacklocal(k) < int(8/numthreads)) then 
         !if (all(threadworking)) then 
            istacklocal(k) = istacklocal(k) + 1
            stacklocal(istacklocal(k)) % nodeindex1 = regnodeindex
            stacklocal(istacklocal(k)) % nodeindex2 = splitnodeindex
         else
          ! If local stack is sufficently full push to global stack 
          !$omp critical (stack)
          top = top + 1
          stack(top) % nodeindex1 = regnodeindex
          stack(top) % nodeindex2 = splitnodeindex 
          !$omp end critical (stack)
        endif 
        !call interact(regnodeindex,splitnodeindex,nodes,x,m,a,nopart)
        !call interact(splitnodeindex,regnodeindex,nodes,x,m,a,nopart)
      endif 
    enddo

     ! UNLOCK THE NODES 
    !print*,"Nodes unlocked"
    !nodes(nodeindex1) % nodefree  = .true.
    !nodes(nodeindex2) % nodefree  = .true.

    ! LEAF-NODE NODE 
    elseif(splitnode % isLeaf .AND. .NOT. regnode % isLeaf) then
      !print*, "FIX THIS"

      ! Call poten function
      ! Need a way to get all of the children of a node for direct sum

      ! Get children of regular node 
      node1empty = .true.
      node2empty = .true.

      do i=1,regnode % bodychildpont
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

      !print*, node1empty, node2empty


      if (.NOT. node1empty .and. .not. node2empty) then
        !call get_accel(x,a,m,np,particleindex,regnodechild)
        !print*, "Crashing here"
        !print*, "Particleindex: ", particleindex
        !print*, "Regnode children: ", regnodechild
        !!$omp critical 
        call get_accel_leafnode(x,a,m,nopart,particleindex,regnodechild)
        !!$omp end critical 
        !print*, "Crashing on second loop"
        !!$omp critical 
        call get_accel_leafnode(x,a,m,nopart,regnodechild,particleindex)
        !!$omp end critical 
        !print*, "Works fine "
        !!$omp end critical 

        
      endif 
     !call get_accel_leafnode(x,a,m,np,particleindex,regnodechild)

     !return
     ! UNLOCK THE NODES 
     !nodes(nodeindex1) % nodefree  = .true.
     !nodes(nodeindex2) % nodefree  = .true.

      ! Get children of a leaf node 

    ! LEAF-NODE LEAF-NODE
    else 
      !print*, "Leafnode-Leafnode"

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

     !print*, "Finished getting bodies"
     !print*, node1empty
     !print*, node2empty
     ! CALL DIRECT SUM 
    if (.NOT. node1empty .and. .not. node2empty) then
      !call get_accel(x,a,m,np,particleindex,particleindex2)
      !!$omp critical 
      call get_accel_leafnode(x,a,m,nopart,particleindex,particleindex2)
      !!$omp end critical 
      !!$omp critical 
      call get_accel_leafnode(x,a,m,nopart,particleindex2,particleindex)
      !!$omp end critical 
      ! UNLOCK THE NODES 
      !print*,"Nodes unlocked"
      !nodes(nodeindex1) % nodefree  = .true.
      !nodes(nodeindex2) % nodefree  = .true.
    endif 
    !call get_accel_leafnode(x,a,m,np,particleindex,particleindex2)

    !return
    !nodes(nodeindex1) % nodefree  = .true.
    !nodes(nodeindex2) % nodefree  = .true.

    endif 


  endif
  !nodes(nodeindex1) % nodefree  = .true.
  !nodes(nodeindex2) % nodefree  = .true.  

  endif 
  
  enddo  
   print*, "Threadcounts: ", threadcounts
   print*, "Idlecounts: ", idlecounts
   !if (idlecounts /= 0 .and. threadcounts /= 0) then 
    print*,"Fraction of time spent idle: "
    print*, real(idlecounts/threadcounts)
   !endif 

   
  !!$OMP end single 
  


  !$omp end parallel 
  print*,"Wall time total: ", wallclock()

 end subroutine interact_stack

 RECURSIVE subroutine interact(nodeindex1,nodeindex2,nodes,x,m,a,nopart,interactionlist)
  integer, intent(in) :: nopart 
  type(octreenode), intent(inout) :: nodes(:)
  real, intent(inout) :: x(3,nopart), m(nopart)
  real, intent(inout) :: a(3,nopart)
  real, optional , intent(inout) :: interactionlist(:,:)
  integer, intent(inout) :: nodeindex1,nodeindex2
  type(octreenode) :: newnode1,newnode2, splitnode, regnode
  integer :: i, j, counter 
  real :: rmax1, rmax2,cm1(3),cm2(3)
  real :: fnode(20),quads(6)
  real :: dr,dx(3),dy,dz,totmass,r,r2
  integer :: particleindex(10),particleindex2(10) 
  real :: c0,c1(3),c2(3,3),c3(3,3,3), c0new,c1new(3),c2new(3,3),c3new(3,3,3)
  logical :: nodesAreEqual, flag
  integer :: regnodechild(10000), splitnodeindex,regnodeindex,newnodeindex1,newnodeindex2,start
  integer, allocatable :: particlesforsum(:)
  logical :: node1empty, node2empty
  
  rmax1 = nodes(nodeindex1) % rmax
  rmax2 = nodes(nodeindex2) % rmax 

  quads(:) = 0.0
  counter = 0

  fnode(:) = 0.0
  particleindex = 0.0
  particleindex2 = 0.0
  !particleindex1 = 0.

  node1empty = .true.
  node2empty = .true.

  regnodechild(:) = 0.
  c0 = 0.
  c1 = 0.
  c2 = 0.
  c3 = 0.
  c0new = 0.
  c1new = 0.
  c2new = 0.
  c3new = 0.



  !print*, "Node indexes: ", nodeindex1, nodeindex2

  ! STILL NEED TO COMPUTE BODY-BODY, BODY-NODE, NODE_BODY 
  ! BODY SELF INTERACTION IS IGNORED
  nodesAreEqual = nodes_equal(nodes(nodeindex1),nodes(nodeindex2))
  if (nodeindex1 == nodeindex2) then
    ! If we are doing leafnode leafnode 
   if (nodes(nodeindex1) % isLeaf) then 
     ! Get index of all bodies 
     do i=1, 10
       if (nodes(nodeindex1) % data(i) /= 0) then
         node1empty = .false.
         particleindex(i) = nodes(nodeindex1) % data(i)
       endif  
     enddo 
    ! CALL DIRECT SUM 
    !print*,"Direct sum"
    !print*, x
    if (.not. node1empty) then
      call get_accel(x,a,m,nopart,particleindex)
    endif 
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
         !print*, "i, j: "
         !print*, i
         !print*, j
         !print*, nodes(nodeindex1) % children(i)
         !print*, nodes(nodeindex2) % children(j)
        ! If children exist
        if (nodes(nodeindex1) % children(i) .NE. 0 .AND. nodes(nodeindex2) % children(j) .NE. 0) then 
         !nodeindex1 = nodes(nodeindex1) % children(i)
         !nodeindex2 = nodes(nodeindex2) % children(j)
         !print*, nodeindex1
         !print*, nodeindex2
         counter = counter + 1
         !print*,"Counter: ", counter

          ! Get the sub-cells
          !newnode1 = nodes(nodeindex1)
          !newnode2 = nodes(nodeindex2)
          newnodeindex1 = nodes(nodeindex1) % children(i)
          newnodeindex2 = nodes(nodeindex2) % children(j)

          !print*, "Nodeindex1: ", newnodeindex1, "Nodeindex2: ", newnodeindex2

          !print*, "Interact called"
          ! call the interact for each of the sub-cells
          call interact(newnodeindex1, newnodeindex2, nodes,x,m,a,nopart)

        endif 

      enddo 
      start = start + 1
    enddo
    endif 




 elseif (well_separated(nodes(nodeindex1),nodes(nodeindex2))) then 
  !elseif (.false.) then 
  !elseif (.true.) then
    !print*, "Well separated!"


    ! SYMMETRY CHECK 
     if ( .not. well_separated(nodes(nodeindex2),nodes(nodeindex1))) then 
        print*, "MAC BROKEN"
        stop
     endif  
   

  

  ! WELL SEPARATED NODE MI IS CALCULATED: TAYLOR COEFFs computed and added to 
  ! node data fields

 
    ! Call taylor COEFFs
    !print*, "Calling taylor coeff: "

    ! PUT MULTIPOLE STUFF HERE 
    ! ------------------------
    ! ------------------------
    ! ------------------------
    ! ------------------------

    if (nodes(nodeindex1) % totalmass /= 0. .AND. nodes(nodeindex2) %totalmass /= 0.) then !&
     !.and. (.NOT. nodes(nodeindex1) % isLeaf) .and. (.not. nodes(nodeindex2) % isLeaf)) then
    cm1 = nodes(nodeindex1) % centerofmass
    cm2 = nodes(nodeindex2) % centerofmass

    !print*, "particles in node 1"
    !print*, nodes(nodeindex1) % data
    !print*, "particles in node 2"
    !print*, nodes(nodeindex2) % data


    call get_dx_dr(cm1,cm2,dx,dr)
    !print*, "Cm of node1: ",cm1
    !print*, "Cm of node2: ",cm2
    !dx = cm1(1) - cm2(1)
    !print*, "Dx: ", dx
    !dy = cm1(2) - cm2(2)
    !print*, "Dy: ", dy
    !dz = cm1(3) - cm2(3)
    !print*, "Dz: ", dz
    !dr = 1./(sqrt(dot_product(cm1-cm2,cm1-cm2)))
    !print*, "Dr: ", dr 

    fnode = 0.0

    totmass = nodes(nodeindex2) % totalmass
    !print*, "Totalmass: ",totmass
    !totmass = 1.0

    ! calculate real accel here 
    !print*, "Real accel: "
     !r2 = dot_product(dx,dx)
     !r  = sqrt(r2)
     !print*, -totmass*(1/((r2)**1.5))*dx 
     !print*, "Real force: "
     !print*, -totmass*(1/((r2)**1.5))*dx * nodes(nodeindex1) %totalmass

    c0 = 0.
    c1 = 0.
    c2 = 0.
    c3 = 0.
    quads = 0.
    !quads = nodes(nodeindex2) % quads 
    !call compute_fnode(dx(1),dx(2),dx(3),dr,totmass,quads,fnode)
    call compute_coeff(dx(1),dx(2),dx(3),dr,totmass,quads,c0,c1,c2,c3)
    nodes(nodeindex1) % fnode = fnode 

    !print*,"C1: ", c1 * nodes(nodeindex1) % totalmass

    ! store coeff for walk phase 
    !node1 % fnode = node1 % fnode + fnode
    !print*, "Stored acccel: "
    nodes(nodeindex1) % c0 =  nodes(nodeindex1)%c0 + c0 
    !print*, nodes(nodeindex1) % c1
    !print*, "calc accel: "
    !print*,c1
    nodes(nodeindex1) % c1 = nodes(nodeindex1)%c1 + c1 
    nodes(nodeindex1) % c2 = nodes(nodeindex1) % c2 + c2 
    nodes(nodeindex1) % c3 = nodes(nodeindex1) % c3 + c3 


    ! Find the "force" from the coefficents to exploit symmetry 
    c0 = c0 * nodes(nodeindex1) % totalmass
    c1 = c1 * nodes(nodeindex1) % totalmass
    c2 = c2 * nodes(nodeindex1) % totalmass
    c3 = c3 * nodes(nodeindex1) % totalmass 

    ! "Force" is equal and opposite 
    c0new = -c0
    c1new = -c1
    c2new = -c2 
    c3new = -c3 

    ! print poten
    !print*, "Poten is: ", fnode(20)
    !call poten_at_bodypos(cm1,cm2,c0,c1,c2,c3,poten(20))
    !print*, "Poten is: ", poten(20)

    !return 

    !open(unit=77,file="wellseperated.txt")
    !write(77,*), "Node 1:", nodeindex1, "Node 2: ", nodeindex2

    !cm2 = nodes(nodeindex1) % centerofmass
    !cm1 = nodes(nodeindex2) % centerofmass

    !print*, "particles in node 1"
    !print*, nodes(nodeindex2) % data
    !print*, "particles in node 2"
    !print*, nodes(nodeindex1) % data


    !call get_dx_dr(cm1,cm2,dx,dr)
    !print*, "Cm of node1: ",cm1
    !print*, "Cm of node2: ",cm2
    !dx = cm1(1) - cm2(1)
    !print*, "Dx: ", dx
    !dy = cm1(2) - cm2(2)
    !print*, "Dy: ", dy
    !dz = cm1(3) - cm2(3)
    !print*, "Dz: ", dz
    !dr = 1./(sqrt(dot_product(cm1-cm2,cm1-cm2)))
    !print*, "Dr: ", dr 

    !fnode = 0.0

    !totmass = nodes(nodeindex1) % totalmass
    !print*, "Totalmass: ",totmass
    !totmass = 1.0

    ! calculate real accel here 
    !print*, "Real accel: "
     !r2 = dot_product(dx,dx)
     !r  = sqrt(r2)
     !print*, -totmass*(1/((r2)**1.5))*dx 
     !print*, "Real force: "
     !print*, -totmass*(1/((r2)**1.5))*dx * nodes(nodeindex2) %totalmass

    !c0new = c0
    !c1new = c1
    !c2new = c2 
    !c3new = c3 
    !c0 = 0.
    !c1 = 0.
    !c2 = 0.
    !c3 = 0.
    !quads = 0.
    !quads = nodes(nodeindex1) % quads 
    !call compute_fnode(dx(1),dx(2),dx(3),dr,totmass,quads,fnode)
    !call compute_coeff(dx(1),dx(2),dx(3),dr,totmass,quads,c0,c1,c2,c3)
    !print*, "c1: ",c1* nodes(nodeindex2) % totalmass

    ! store coeff for walk phase 
    !nodes(nodeindex2) % fnode = nodes(nodeindex2) % fnode + fnode
    !print*, "Stored acccel: "

    ! Divide by mass to get "accel"
    c0 = c0new / (nodes(nodeindex2) % totalmass)
    c1 = c1new / (nodes(nodeindex2) % totalmass) 
    c2 = c2new / (nodes(nodeindex2) % totalmass)
    c3 = c3new / (nodes(nodeindex2) % totalmass)

    nodes(nodeindex2) % c0 =  nodes(nodeindex2)%c0 + c0
    !print*, nodes(nodeindex2) % c1
    !print*, "calc accel: "
    !print*,c1
    nodes(nodeindex2) % c1 = nodes(nodeindex2)%c1 + c1
    nodes(nodeindex2) % c2 = nodes(nodeindex2) % c2 + c2
    nodes(nodeindex2) % c3 = nodes(nodeindex2) % c3 + c3

    ! print poten
    !print*, "Poten is: ", fnode(20)
    !call poten_at_bodypos(cm1,cm2,c0,c1,c2,c3,poten(20))
    !print*, "Poten is: ", poten(20)

    !return 
    !STOP

   !open(unit=77,file="wellseperated.txt")
    !write(77,*), "Node 1:", nodeindex2, "Node 2: ", nodeindex1

   !print*,"Values for interaciton"
    !print*, c1 
    !print*,abs(c0 * nodes(nodeindex2) % totalmass) -  abs(c0new * totmass)
    !print*, abs(c1 *nodes(nodeindex2) % totalmass) - abs(c1new * totmass)
    !print*,abs(c2 * nodes(nodeindex2) % totalmass) - abs(c2new * totmass)
    !print*, abs(c3 * nodes(nodeindex2) % totalmass) - abs(c3new * totmass)

    !print*, "Are values symmetric ? ", c1, c1new

    !if ( abs(c1(1) * nodes(nodeindex2) % totalmass) - abs(c1new(1) * totmass) > 1.e-18) stop

    !STOP

  endif 

  ! THE NODE WITH THE LARGER RMAX IS SPLIT; up to 8 new MI are created and processed 
  else
    !print*, "Split nodes"

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
      !print*, i
      !print*, splitnode % children(i)
      if (splitnode % children(i) /= 0 .and. nodes(splitnode % children(i)) % totalmass /= 0.) then 
        splitnodeindex = splitnode % children(i)
        !newnode1 = nodes(nodeindex1)
        !print*, "Internal splitnode case"
        !print*, "regnode index: ", regnodeindex
        !print*, "Split node index: ",splitnodeindex
        call interact(regnodeindex,splitnodeindex,nodes,x,m,a,nopart)
        !call interact(splitnodeindex,regnodeindex,nodes,x,m,a,nopart)
      endif 
    enddo

    ! LEAF-NODE NODE 
    elseif(splitnode % isLeaf .AND. .NOT. regnode % isLeaf) then
      !print*, "FIX THIS"

      ! Call poten function
      ! Need a way to get all of the children of a node for direct sum

      ! Get children of regular node 
      node1empty = .true.
      node2empty = .true.

      do i=1,regnode % bodychildpont
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

      !print*, node1empty, node2empty


      if (.NOT. node1empty .and. .not. node2empty) then
        !call get_accel(x,a,m,np,particleindex,regnodechild)
        !print*, "Crashing here"
        !print*, "Particleindex: ", particleindex
        !print*, "Regnode children: ", regnodechild
        call get_accel_leafnode(x,a,m,nopart,particleindex,regnodechild)
        !print*, "Crashing on second loop"
        call get_accel_leafnode(x,a,m,nopart,regnodechild,particleindex)
        !print*, "Works fine "
      endif 
     !call get_accel_leafnode(x,a,m,np,particleindex,regnodechild)

     return

      ! Get children of a leaf node 

    ! LEAF-NODE LEAF-NODE
    else 
      !print*, "Leafnode-Leafnode"

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

     !print*, "Finished getting bodies"
     !print*, node1empty
     !print*, node2empty
     ! CALL DIRECT SUM 
    if (.NOT. node1empty .and. .not. node2empty) then
      !call get_accel(x,a,m,np,particleindex,particleindex2)
      call get_accel_leafnode(x,a,m,nopart,particleindex,particleindex2)
      call get_accel_leafnode(x,a,m,nopart,particleindex2,particleindex)
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


  !h = 0.
  h = 0.1
  dx = x1 - x2
  dr = 1./(sqrt(dot_product(dx,dx)))

 end subroutine get_dx_dr

 !function thread_free(threadworking)
 ! logical :: threadworking(:)
 ! integer :: thesize,i

  !thesize = size(threadworking)

  !do i=1,thesize 
  !  if (threadworking(i) .eqv. .false.) then 
  !    return i 
  !  endif 
  !enddo 

!end function thread_free


subroutine global_to_local(stack,stacklocal,threadnumber,istacklocal,top)
  type(interact_global_stack), intent(inout) :: stack(:)
  type(interact_stack_data), intent(inout) :: stacklocal(:)
  integer, intent(in) :: threadnumber
  integer, intent(inout) :: istacklocal(:),top 
  integer :: theinteractions(8), thenode,i 


  ! POP FROM GLOBAL AND PUSH TO LOCAL STACK 

  thenode = 0 
  theinteractions = 0

  ! POP FROM GLOBAL
  !!$omp critical (stack)
  thenode = stack(top) % nodeindex1
  theinteractions =  stack(top) % interactions
  ! cleanup interactions
  stack(top) % interactions = 0
  top = top - 1 
  !!$omp end critical (stack)


  ! PUSH TO LOCAL 
  do i=1,8
    if (theinteractions(i) /= 0) then 
      !print*,"i is: ", i
      istacklocal(threadnumber) = istacklocal(threadnumber) + 1
      stacklocal(istacklocal(threadnumber)) % nodeindex1 = thenode
      stacklocal(istacklocal(threadnumber)) % nodeindex2 = theinteractions(i)
      !print*, "the node", thenode
      !print*, "the interactions ", theinteractions(i)
    endif 
  enddo 

  !print*, "stacklocal", stacklocal(1:istacklocal(threadnumber))
  !stop 

end subroutine global_to_local



 
end module interaction 