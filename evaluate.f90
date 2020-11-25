module evaluate
 use octreetype
 use omp_lib 
 implicit none 

type evaluate_stack_data

 type(octreenode) :: node 
 real :: z0(3)
 ! real :: c0,c1(3),c2(3,3),c3(3,3,3)
 real :: fnode(20)

end type 

contains 

 subroutine evaluate_gravity_parallel(nodes,x,accel)
  type(octreenode), intent(inout) :: nodes(:)
   real, intent(in) :: x(:,:)
   real, intent(inout) :: accel(:,:)
   integer :: stacksize, top,numthreads,iter,nodeindex
   ! Openmp disables heap allocation so this should be allocatable 
   type(evaluate_stack_data) :: stack(1000)
   real :: c0,c1(3),c2(3,3),c3(3,3,3),z0(3)


   top = 1

   ! Push the root node 
   !stack(1) % node = nodes(1)
   z0 = 0.
   c0 = 0.
   c1 = 0.
   c2 = 0.
   c3 = 0.
   
   print*, "Before openmp loop"
   numthreads = 1
   ! get number of OpenMPthreads
   !$omp parallel default(none) shared(numthreads)
     numthreads = omp_get_num_threads()
   !$omp end parallel

   print*, "Number of threads is: ", numthreads

   !$omp parallel default(none) &
   !$omp shared(nodes,x,accel) &
   !$omp private(nodeindex,z0,c0,c1,c2,c3)
   !$omp single 
   !!$omp task shared(nodeindex,nodes,x,accel) &
   
   print*,"nodeindex: ",nodeindex
   nodeindex = 1 
   !!$omp task 
    call eval(nodeindex,nodes,z0,c0,c1,c2,c3,x,accel)
        !if (top == 0) stop 
   !!$omp end task 
      !!$omp taskwait 

   !$omp end single 
   !$omp end parallel 


 end subroutine evaluate_gravity_parallel


 subroutine eval(nodeindex,nodes,z0,c0,c1,c2,c3,x,accel)
  type(octreenode), intent(inout) :: nodes(:)
  real,intent(inout) :: c0,c1(3),c2(3,3),c3(3,3,3),z0(3),accel(:,:)
  real, intent(in) :: x(:,:)!,c0total(:),c1(3,:),c2(3,3,:),c3(3,3,3,:)
  !real, intent(inout) :: asum(3) 
  integer, intent(inout) :: nodeindex
  integer :: childnode
  real :: z1(3),xbod(3),bodpot,bodaccel(3),accelbef(3)
  integer :: i,nochild,bodyindex
  real :: c0new, c1new(3), c2new(3,3),c3new(3,3,3),z0new(3),fnode(20),dx(3),m,nodemass
  real :: c0copy,c1copy(3),c2copy(3,3),c3copy(3,3,3)
  integer  :: nthreads
  real :: nodeaccelsum(3)
 

  bodpot = 0.
  bodaccel = 0.
  bodyindex = 0.

  c0new = 0.
  c1new = 0.
  c2new = 0.
  c3new = 0.
  nodeaccelsum = 0.
  childnode = 0
  z1 = 0.
  xbod = 0.
  bodpot = 0.
  bodaccel = 0.
  accelbef = 0.
  i = 0
  nochild = 0
  bodyindex = 0
  z0new = 0.
  fnode = 0.
  dx = 0. 
  nodemass = 0.
  m = 1./(size(accel)/3.)

   print*, "Node index is: ", nodeindex

   print*, "Center of mass (old): ", z0
  ! Get CoM of current node
   z1 = nodes(nodeindex) % centerofmass
   print*, "Center of Mass new: "
   print*, z1

   nodemass = nodes(nodeindex) % totalmass
   print*, "Node mass: ", nodes(nodeindex) % totalmass
   print*, "c1 is: ", c1

  ! TRANSLATE TAYLOR SERIES T0 TO CENTER OF MASS OF A
  
  ! ! COPY VALUES SO THEY ARENT CHANGED BY TRANSLATION
  ! c0new = c0
  ! c1new = c1
  ! c2new = c2
  ! c3new = c3 
  dx(1) = z1(1)-z0(1)
  dx(2) = z1(2)-z0(2)
  dx(3) = z1(3)-z0(3)
  !call translate_expansion_center(z0,z1,c0new,c1new,c2new,c3new)
  call translate_fgrav_in_taylor_series(fnode,dx(1),dx(2),dx(3))
  !if (c0 /= 0. .and. nodes(nodeindex) % totalmass /= 0.) call translate_expansion_center(z0,z1,c0,c1,c2,c3)

  !dr = z0-z1

  !call expand_fgrav_in_taylor_series()
  
  ! TA += T0
  ! Sum up the field tensors in a taylor series to get poten
  ! Accumulate field tensors 

  !print*, "Node taylor coeffs: "
  !print*, "c0 node: "
  !print*, node % c0 
  ! print*, "c1 node: "
  ! print*, nodes(nodeindex) %c1 * nodes(nodeindex) % totalmass 
  ! !print*, "c2 node: "
  ! !print*, nodes(nodeindex) % c2
  ! !print*, "c3 node: "
  ! !print*, nodes(nodeindex) % c3


  ! !if (node % c0 /= 0.) then
  ! c0new = (nodes(nodeindex) % c0) + c0new
  ! !c0new = c0
  ! !print*, "C0: "
  ! !print*,c0new
  ! c1new = (nodes(nodeindex) % c1) + c1new
  ! !c1new = c1 
  ! !print*, "C1: "
  ! !print*,c1new
  ! !c2new = c2
  ! c2new = (nodes(nodeindex) % c2) + c2new
  ! !print*, "C2: "
  ! !print*,c2new
  ! c3new = (nodes(nodeindex) % c3) + c3new
  !c3new = c3 
  !print*, "C3: "
  !print*, c3new
  ! FOR BODY CHILDREN OF A 

  !else
  !  c0new = c0
  !  c1new = c1
  !  c2new = c2
  !  c3new = c3 
  !endif 
 !STOP

 nodes(nodeindex) % fnode = nodes(nodeindex) % fnode + fnode 

  !nochild = node % bodychildpont
  if (nodes(nodeindex) % isLeaf) then
  do i=1, 10
    
    !if (node%data(i)==0) then
    !  EXIT
    !endif 
    ! Get the index of the body 
    !bodyindex = node % bodychildren(i)

    if (nodes(nodeindex) % data(i) /= 0 .and. c0new /= 0.) then
    !print*, nodes(nodeindex) % data(i) 
    bodyindex = nodes(nodeindex) % data(i)
    !print*, "Bodyindex: ",bodyindex
    xbod = x(:,bodyindex)



   ! Evaluate TA at body's position
   !call poten_at_bodypos(xbod,z1,c0,c1,c2,c3,bodpot)
   bodaccel = 0.
   print*, "c0new: ", c0new
   dx = xbod - z1 
   bodaccel = 0.
   fnode = nodes(nodeindex) % fnode
   !call accel_at_bodypos(xbod,z1,c0,c1,c2,c3,bodaccel)
   call expand_fgrav_in_taylor_series(fnode,dx(1),dx(2),dx(3),bodaccel(1),bodaccel(2),bodaccel(3),bodpot)
   !accel(:,bodyindex) = accel(:,bodyindex) + bodaccel
   !call accel_at_bodypos(xbod,z1,c0new,c1new,c2new,c3new,bodaccel)

   ! add to body's potential and acceleration

   !poten(bodyindex) = poten(bodyindex) + bodpot
   ! print*, "accel bef:",accel(:,bodyindex)
   ! !accelbef = accel(:,bodyindex)
   ! !accel(:,bodyindex) = accel(:,bodyindex) + c1 
   ! accel(:,bodyindex) = accel(:,bodyindex) + bodaccel
   ! print*, "Accel now: ", accel(:,bodyindex)
   ! !asum = asum + bodaccel*m
   ! nodeaccelsum = nodeaccelsum + bodaccel*m
   !print*, "Asum: ", asum
   !print*, "bodaccel: ", bodaccel
   !print*, "delta accel:"
   !print*, accel(:,bodyindex) - accelbef

    


   endif 


  enddo 
  print*, "Part accelsum: ",nodeaccelsum
  print*, "c1: ", c1new * nodes(nodeindex) % totalmass   
  !write(88,*) "Node index: ", nodeindex, " Force sum: ", c1new * nodes(nodeindex) % totalmass, "Part accelsum: ",nodeaccelsum, &
  ! "Center of mass: ", nodes(nodeindex) %centerofmass, "Node mass: ", nodes(nodeindex) % totalmass
 else 

  ! FOR CHILDREN OF A 
  ! change center of mass
  z0new = z1
  !!$OMP DO 
  ! print*,"Reached this point"
  ! if (nodeindex == 1) then 
  !   open(88,file="forcesums.txt",position="append")
  ! endif 
  !write(88,*) "Node index: ", nodeindex, " Force sum: ", c1new, &
   !"Center of mass: ", nodes(nodeindex) %centerofmass, "Node mass: ", nodes(nodeindex) % totalmass, &
   !"c2: ", nodes(nodeindex) % c2 
  
  !close(88)
  !!$omp parallel default(none) &
  !!$omp shared(nodes,x,accel,nodeindex) &
  !!$omp private(childnode,z0new,c0new,c1new,c2new,c3new,i) 
  !!$omp single 
  do i=1, 8
    !print*, "Childnode index: ", nodes(nodeindex) % children(i)
    !print*, "Has c1 changed: ",c1new
    !print*, "Has z0new changed: ",z0new
  if (nodes(nodeindex) % children(i) /= 0 .and. norm2(nodes(nodes(nodeindex) %children(i)) % centerofmass) /= 0. ) then
      
    
    !!$omp task 
    print*, "Childnode index: ", nodes(nodeindex) % children(i)
    childnode = nodes(nodeindex) % children(i)
    call eval(childnode,nodes,z0new,c0new,c1new,c2new,c3new,x,accel)
    !!$omp end task 
    !print*, nodes(nodeindex) % children(i)
    !c1new = c1copy 
   endif 

   !if (nodeindex == 1) stop 

  enddo
  !!$omp end single 
  !!$omp end parallel 
  !stop 
  
  !!$OMP ENDDO   

  endif 



 end subroutine eval

  

 subroutine evaluate_gravity_stack(nodes,x,accel)
   type(octreenode), intent(inout) :: nodes(:)
   real, intent(in) :: x(:,:)
   real, intent(inout) :: accel(:,:)
   integer :: stacksize, top,toplocal 
   ! Openmp disables heap allocation so this should be allocatable 
   type(evaluate_stack_data) :: stack(1000)
   type(octreenode) :: currentnode,childnode
   real :: c0,c1(3),c2(3,3),c3(3,3,3),z0(3),z1(3)
   real :: bodaccel(3),xbod(3),fnode(20),bodpot,dx(3)
   integer :: bodyindex,i,j,iter,numthreads,toptemp,k
   logical, allocatable :: threadworking(:)
   integer, allocatable :: istacklocal(:)

   top = 1
   stacksize = 1
   iter = 1

   ! Push the root node 
   stack(1) % node = nodes(1)
   stack(1) % z0 = 0.
   stack(1) % fnode = 0.
   ! stack(1) % c0 = 0.
   ! stack(1) % c1 = 0.
   ! stack(1) % c2 = 0.
   ! stack(1) % c3 = 0.
   currentnode = nodes(1)
   print*, "Root node pushed"

    
  !print*, "top is: ", top 
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

  !$omp parallel default(none) &
  !$omp shared(stack,top,x,accel,nodes,numthreads,threadworking,istacklocal) &
  !$omp private(z0,z1,c0,c1,c2,c3) &
  !$omp private(currentnode,bodyindex,xbod,bodaccel,childnode,dx,k,fnode,bodpot)
  
   !$ k=omp_get_thread_num() + 1
   ! local stack is currently empty 
   istacklocal(k) = 0 
   ! Just for compiler errors 
   currentnode = nodes(1)
   !print*,"Thread is: ", k
   ! here top is the top of the global stack 
  over_stack: do while (any(threadworking) .or. top > 0)  
    print*, top 
   !do j=1, top
    !print*, "Thread number is: ", j
    !print*, "The current iteration is: ",iter
    !print*, "Top is: ", top
    !print*, "local stack top: ", istacklocal(k)
    !iter = iter + 1 
    ! POP ITEM FROM STACK 
    !print*,"istacklocal: ",istacklocal(k)
    ! pop of local stack 
    !$omp critical(globalstack) 
        if (top > 0) then 
        currentnode = stack(top) % node
        z0 = stack(top) % z0
        fnode = stack(top) % fnode
        top = top - 1 
        threadworking(k) = .true. 
        !print*,"working on thread",k
      else
       threadworking(k) = .false.
        
      endif 
    !$omp end critical(globalstack)
    !if top
  ! If thread has work to do, do it 
  if (threadworking(k)) then 
   

    ! Perform evaluate on popped item 

   z1 = currentnode % centerofmass
   print*, "z1:", z1 
   !stop 

   dx(1) = z1(1)-z0(1)
   dx(2) = z1(2)-z0(2)
   dx(3) = z1(3)-z0(3)
  

  ! TRANSLATE TAYLOR SERIES T0 TO CENTER OF MASS OF A
  !print*, "Translating expansion: "
  ! call translate_expansion_center(z0,z1,c0,c1,c2,c3)

  ! c0 = currentnode % c0 + c0
  ! !c0new = c0
  ! !print*, "C0: ",c0
  ! !print*,c0new
  ! c1 = currentnode % c1 + c1
  ! !c1new = c1 
  ! !print*, "C1: ",c1
  ! !print*,c1new
  ! !c2new = c2
  ! c2 = currentnode % c2 + c2
  ! !print*, "C2: "
  ! !print*,c2new
  ! c3 = currentnode % c3 + c3

  print*, "fnodebefore: ", fnode
  print*, "currentnode fnode", currentnode % fnode
  call translate_fgrav_in_taylor_series(fnode,dx(1),dx(2),dx(3))
  print*, "fnode after: ", fnode
  currentnode % fnode = currentnode % fnode + fnode 
  fnode = currentnode % fnode
  !print*, "currentnode fnode: ", currentnode % fnode

  !nochild = node % bodychildpont
  if (currentnode % isLeaf) then
  do i=1, 10

  if (currentnode % data(i) /= 0) then
   bodyindex = currentnode % data(i)
   xbod = x(:,bodyindex)
   dx = xbod - z1 
   bodaccel = 0.
   !fnode = currentnode % fnode
   !print*,"calculated accel: "
   !call accel_at_bodypos(xbod,z1,c0,c1,c2,c3,bodaccel)
   !print*,bodaccel
   call expand_fgrav_in_taylor_series(fnode,dx(1),dx(2),dx(3),bodaccel(1),bodaccel(2),bodaccel(3),bodpot)
   !$omp critical (accel)
   print*, "bodaccel: ", bodaccel
   accel(:,bodyindex) = accel(:,bodyindex) + bodaccel
   print*, "accel: ",accel(:,bodyindex)
   !stop 
   !$omp end critical(accel)
  endif 
  enddo 
 else 

  ! FOR CHILDREN OF A 
  ! change center of mass
  z0 = z1
  
  do i=1, 8



   if (currentnode % children(i) /= 0) then  
    !print*, "Childnode index: ", currentnode % children(i)
    childnode = nodes(currentnode % children(i))
    !nodeindex = node % children(i)
    !print*, "Node index: ",nodeindex
    !call evaluate_gravity(childnode,nodes,z0new,c0new,c1new,c2new,c3new,x,accel) 

    ! if threads are waiting push to global stack 
    !if (any(threadworking) .eqv. .false.) then 
      !$omp critical(globalstack)
      top = top + 1 
      ! Push Children onto stack 
      stack(top) % node = childnode
      stack(top) % z0 = z0
      stack(top) % fnode = fnode 
      !$omp end critical(globalstack)


   endif 

  enddo

 

  endif 

endif 

    
  !print*, "top is: ", top 
  enddo over_stack
  !$OMP end parallel

  !stop 


 end subroutine evaluate_gravity_stack


 ! BE VERY CAREFUL WITH YOUR RECUSRIVE VARIABLES 

 recursive subroutine evaluate_gravity(nodeindex,nodes,z0,c0,c1,c2,c3,x,accel,asum)
  type(octreenode), intent(inout) :: nodes(:)
  real,intent(inout) :: c0,c1(3),c2(3,3),c3(3,3,3),z0(3),accel(:,:),x(:,:)!,c0total(:),c1(3,:),c2(3,3,:),c3(3,3,3,:)
  real, intent(inout) :: asum(3) 
  integer, intent(inout) :: nodeindex
  integer :: childnode
  real :: z1(3),xbod(3),bodpot,bodaccel(3),accelbef(3)
  integer :: i,nochild,bodyindex
  real :: c0new, c1new(3), c2new(3,3),c3new(3,3,3),z0new(3),fnode(20),dx(3),m,nodemass
  real :: c0copy,c1copy(3),c2copy(3,3),c3copy(3,3,3)
  integer  :: nthreads
  real :: nodeaccelsum(3)
  !real :: c0old,c1old(3),c2old(3,3),c3old(3,3,3)

  bodpot = 0.
  bodaccel = 0.
  bodyindex = 0.

  c0new = 0.
  c1new = 0.
  c2new = 0.
  c3new = 0.
  nodeaccelsum = 0.
  childnode = 0
  z1 = 0.
  xbod = 0.
  bodpot = 0.
  bodaccel = 0.
  accelbef = 0.
  i = 0
  nochild = 0
  bodyindex = 0
  z0new = 0.
  fnode = 0.
  dx = 0. 
  nodemass = 0.
  m = 1./(size(accel)/3.)

   print*, "Center of mass (old): ", z0
  ! Get CoM of current node
   z1 = nodes(nodeindex) % centerofmass
   print*, "Center of Mass new: "
   print*, z1

   nodemass = nodes(nodeindex) % totalmass
   print*, "Node mass: ", nodes(nodeindex) % totalmass
   print*, "c1 is: ", c1

  ! TRANSLATE TAYLOR SERIES T0 TO CENTER OF MASS OF A
  
  ! COPY VALUES SO THEY ARENT CHANGED BY TRANSLATION
  c0new = c0
  c1new = c1
  c2new = c2
  c3new = c3 
  call translate_expansion_center(z0,z1,c0new,c1new,c2new,c3new)
  !if (c0 /= 0. .and. nodes(nodeindex) % totalmass /= 0.) call translate_expansion_center(z0,z1,c0,c1,c2,c3)

  !dr = z0-z1

  !call expand_fgrav_in_taylor_series()
  
  ! TA += T0
  ! Sum up the field tensors in a taylor series to get poten
  ! Accumulate field tensors 

  !print*, "Node taylor coeffs: "
  !print*, "c0 node: "
  !print*, node % c0 
  print*, "c1 node: "
  print*, nodes(nodeindex) %c1 * nodes(nodeindex) % totalmass 
  !print*, "c2 node: "
  !print*, nodes(nodeindex) % c2
  !print*, "c3 node: "
  !print*, nodes(nodeindex) % c3


  !if (node % c0 /= 0.) then
  c0new = (nodes(nodeindex) % c0) + c0new
  !c0new = c0
  !print*, "C0: "
  !print*,c0new
  c1new = (nodes(nodeindex) % c1) + c1new
  !c1new = c1 
  !print*, "C1: "
  !print*,c1new
  !c2new = c2
  c2new = (nodes(nodeindex) % c2) + c2new
  !print*, "C2: "
  !print*,c2new
  c3new = (nodes(nodeindex) % c3) + c3new
  !c3new = c3 
  !print*, "C3: "
  !print*, c3new
  ! FOR BODY CHILDREN OF A 

  !else
  !  c0new = c0
  !  c1new = c1
  !  c2new = c2
  !  c3new = c3 
  !endif 
 !STOP

  !nochild = node % bodychildpont
  if (nodes(nodeindex) % isLeaf) then
  do i=1, 10
    
    !if (node%data(i)==0) then
    !  EXIT
    !endif 
    ! Get the index of the body 
    !bodyindex = node % bodychildren(i)

    if (nodes(nodeindex) % data(i) /= 0 .and. c0new /= 0.) then
    !print*, nodes(nodeindex) % data(i) 
    bodyindex = nodes(nodeindex) % data(i)
    !print*, "Bodyindex: ",bodyindex
    xbod = x(:,bodyindex)



   ! Evaluate TA at body's position
   !call poten_at_bodypos(xbod,z1,c0,c1,c2,c3,bodpot)
   bodaccel = 0.
   print*, "c0new: ", c0new
   dx = xbod - z1 
   bodaccel = 0.
   !fnode = node % fnode
   !call accel_at_bodypos(xbod,z1,c0,c1,c2,c3,bodaccel)
   !call expand_fgrav_in_taylor_series(fnode,dx(1),dx(2),dx(3),bodaccel(1),bodaccel(2),bodaccel(3),bodpot)
   !accel(:,bodyindex) = accel(:,bodyindex) + bodaccel
   call accel_at_bodypos(xbod,z1,c0new,c1new,c2new,c3new,bodaccel)

   ! add to body's potential and acceleration

   !poten(bodyindex) = poten(bodyindex) + bodpot
   print*, "accel bef:",accel(:,bodyindex)
   !accelbef = accel(:,bodyindex)
   !accel(:,bodyindex) = accel(:,bodyindex) + c1 
   accel(:,bodyindex) = accel(:,bodyindex) + bodaccel
   print*, "Accel now: ", accel(:,bodyindex)
   !asum = asum + bodaccel*m
   nodeaccelsum = nodeaccelsum + bodaccel*m
   !print*, "Asum: ", asum
   !print*, "bodaccel: ", bodaccel
   !print*, "delta accel:"
   !print*, accel(:,bodyindex) - accelbef

    


   endif 


  enddo 
  asum = asum + nodeaccelsum
  print*,"Force sum: ", asum
  print*, "Part accelsum: ",nodeaccelsum
  print*, "c1: ", c1new * nodes(nodeindex) % totalmass   
  write(88,*) "Node index: ", nodeindex, " Force sum: ", c1new * nodes(nodeindex) % totalmass, "Part accelsum: ",nodeaccelsum, &
   "Center of mass: ", nodes(nodeindex) %centerofmass, "Node mass: ", nodes(nodeindex) % totalmass
 else 

  ! FOR CHILDREN OF A 
  ! change center of mass
  z0new = z1
  !!$OMP DO 
  print*,"Reached this point"
  if (nodeindex == 1) then 
    open(88,file="forcesums.txt",position="append")
  endif 
  write(88,*) "Node index: ", nodeindex, " Force sum: ", c1new, &
   "Center of mass: ", nodes(nodeindex) %centerofmass, "Node mass: ", nodes(nodeindex) % totalmass, &
   "c2: ", nodes(nodeindex) % c2 
  
  !close(88)
  do i=1, 8
    print*, "Childnode index: ", nodes(nodeindex) % children(i)
    print*, "Has c1 changed: ",c1new
    print*, "Has z0new changed: ",z0new
  if (nodes(nodeindex) % children(i) /= 0 .and. norm2(nodes(nodes(nodeindex) %children(i)) % centerofmass) /= 0. ) then
      
    print*, "Childnode index: ", nodes(nodeindex) % children(i)
    childnode = nodes(nodeindex) % children(i)
    c1copy = c1new 
    !nodeindex = node % children(i)
    !print*, "Node index: ",nodeindex
    call evaluate_gravity(childnode,nodes,z0new,c0new,c1new,c2new,c3new,x,accel,asum)
    print*, nodes(nodeindex) % children(i)
    c1new = c1copy 
   endif 

   !if (nodeindex == 1) stop 

  enddo
  !stop 
  write(88,*) "Asum after summing: ",asum
  asum = 0.
  !!$OMP ENDDO   

  endif 

 end subroutine evaluate_gravity


pure subroutine translate_fgrav_in_taylor_series(fnode,dx,dy,dz)
 real, intent(inout)  :: fnode(20)
 real, intent(in)  :: dx,dy,dz
 real :: fnodeold(20)
 real :: dfxx,dfxy,dfxz,dfyy,dfyz,dfzz
 real :: d2fxxx,d2fxxy,d2fxxz,d2fxyy,d2fxyz,d2fxzz,d2fyyy,d2fyyz,d2fyzz,d2fzzz

 ! fxi = fnode(1)
 ! fyi = fnode(2)
 ! fzi = fnode(3)
 dfxx = fnode(4)
 dfxy = fnode(5)
 dfxz = fnode(6)
 dfyy = fnode(7)
 dfyz = fnode(8)
 dfzz = fnode(9)

 

 d2fxxx = fnode(10)
 d2fxxy = fnode(11)
 d2fxxz = fnode(12)
 d2fxyy = fnode(13)
 d2fxyz = fnode(14)
 d2fxzz = fnode(15)
 d2fyyy = fnode(16)
 d2fyyz = fnode(17)
 d2fyzz = fnode(18)
 d2fzzz = fnode(19)
 ! poti = fnode(20)

 fnodeold = fnode 

 fnode(1) = fnodeold(1) + dx*(dfxx  + 0.5*(dx*d2fxxx + dy*d2fxxy + dz*d2fxxz)) &
           + dy*(dfxy  + 0.5*(dx*d2fxxy + dy*d2fxyy + dz*d2fxyz)) &
           + dz*(dfxz    + 0.5*(dx*d2fxxz + dy*d2fxyz + dz*d2fxzz))
 fnode(2) = fnodeold(2) + dx*(dfxy   + 0.5*(dx*d2fxxy + dy*d2fxyy + dz*d2fxyz)) &
           + dy*(dfyy   + 0.5*(dx*d2fxyy + dy*d2fyyy + dz*d2fyyz)) &
           + dz*(dfyz   + 0.5*(dx*d2fxyz + dy*d2fyyz + dz*d2fyzz))
 fnode(3) = fnodeold(3) + dx*(dfxz   + 0.5*(dx*d2fxxz + dy*d2fxyz + dz*d2fxzz)) &
           + dy*(dfyz   + 0.5*(dx*d2fxyz + dy*d2fyyz + dz*d2fyzz)) &
           + dz*(dfzz   + 0.5*(dx*d2fxzz + dy*d2fyzz + dz*d2fzzz))

 ! ! The c2 field tensor to be translated 

  fnode(4) = fnodeold(4) + dx*(d2fxxx + d2fxxy + d2fxxz)
  fnode(5) = fnodeold(5) + dy*(d2fxxy + d2fxyy + d2fxyz)
  fnode(6) = fnodeold(6) + dz*(d2fxxz + d2fxyz + d2fxzz)
  fnode(7) = fnodeold(7) + dy*(d2fxyy + d2fyyy + d2fyyz)
  fnode(8) = fnodeold(8) + dz*(d2fxyz + d2fyyz + d2fyzz)
  fnode(9) = fnodeold(9) + dz*(d2fxzz + d2fyzz + d2fzzz)
 ! endif 
   fnode(10:20) = fnodeold(10:20)
  !poti = poti - (dx*fxi + dy*fyi + dz*fzi)

 return
end subroutine translate_fgrav_in_taylor_series

 subroutine translate_expansion_center(z0,z1,c0,c1,c2,c3)
  real, intent(in) :: z0(3), z1(3)
  real, intent(inout) :: c0, c1(3), c2(3,3), c3(3,3,3)
  real :: c0old, c1old(3), c2old(3,3), c3old(3,3,3)
  real :: sep1(3), sep2(3,3), sep3(3,3,3)
  real :: sep1c2comp(3), sep2c3comp(3), sep1c3comp(3,3)
  
  c0old = c0
  c1old = c1
  c2old = c2
  c3old = c3
  !print*, "c0"
  !print*, c0old
  !print*, "c1"
  print*, c1old
  !print*, "c2"
  print*, c2old
  !print*, "c3"
  !print*, c3old(2,1,1)
  sep1 = 0.
  sep1 = z1 - z0
  !print*, "sep1: ", sep1
  call outer_product1(sep1,sep1,sep2)
  call outer_product2(sep2,sep1,sep3)

  !print*, c2 

  sep1c2comp(1) = sep1(1)*c2(1,1) + sep1(2)*c2(1,2) + sep1(3)*c2(1,3)
  sep1c2comp(2) = sep1(1)*c2(2,1) + sep1(2)*c2(2,2) + sep1(3)*c2(2,3)
  sep1c2comp(3) = sep1(1)*c2(3,1) + sep1(2)*c2(3,2) + sep1(3)*c2(3,3)


  sep2c3comp(1) = sep1(1)*(sep1(1)*c3(1,1,1) + sep1(2)*c3(1,1,2) + sep1(3)*c3(1,1,3)) &
                + sep1(2)*(sep1(1)*c3(1,2,1)+ sep1(2)*c3(1,2,2) + sep1(3)*c3(1,2,3)) &
                + sep1(3)*(sep1(1)*c3(1,3,1)+ sep1(2)*c3(1,3,2) + sep1(3)*c3(1,3,3))

  sep2c3comp(2) = sep1(1)*(sep1(1)*c3(2,1,1) + sep1(2)*c3(2,1,2) + sep1(3)*c3(2,1,3)) &
                + sep1(2)*(sep1(1)*c3(2,2,1)+ sep1(2)*c3(2,2,2) + sep1(3)*c3(2,2,3)) &
                + sep1(3)*(sep1(1)*c3(2,3,1)+ sep1(2)*c3(2,3,2) + sep1(3)*c3(2,3,3))

  sep2c3comp(3) = sep1(1)*(sep1(1)*c3(3,1,1) + sep1(2)*c3(3,1,2) + sep1(3)*c3(3,1,3)) &
                + sep1(2)*(sep1(1)*c3(3,2,1)+ sep1(2)*c3(3,2,2) + sep1(3)*c3(3,2,3)) &
                + sep1(3)*(sep1(1)*c3(3,3,1)+ sep1(2)*c3(3,3,2) + sep1(3)*c3(3,3,3))


  sep1c3comp(1,1) = sep1(1)*c3(1,1,1) + sep1(2)*c3(1,1,2) + sep1(3)*c3(1,1,3)
  sep1c3comp(1,2) = sep1(1)*c3(1,2,1) + sep1(2)*c3(1,2,2) + sep1(3)*c3(1,2,3)
  sep1c3comp(1,3) = sep1(1)*c3(1,3,1) + sep1(2)*c3(1,3,2) + sep1(3)*c3(1,3,3)
  sep1c3comp(2,1) = sep1(1)*c3(2,1,1) + sep1(2)*c3(2,1,2) + sep1(3)*c3(2,1,3)
  sep1c3comp(2,2) = sep1(1)*c3(2,2,1) + sep1(2)*c3(2,2,2) + sep1(3)*c3(2,2,3)
  sep1c3comp(2,3) = sep1(1)*c3(2,3,1) + sep1(2)*c3(2,3,2) + sep1(3)*c3(2,3,3)
  sep1c3comp(3,1) = sep1(1)*c3(3,1,1) + sep1(2)*c3(3,1,2) + sep1(3)*c3(3,1,3)
  sep1c3comp(3,2) = sep1(1)*c3(3,2,1) + sep1(2)*c3(3,2,2) + sep1(3)*c3(3,2,3)
  sep1c3comp(3,3) = sep1(1)*c3(3,3,1) + sep1(2)*c3(3,3,2) + sep1(3)*c3(3,3,3)

  !print*, "Second term value: "
  !print*, 0.5*inner_product2(sep2,c2old)
  !print*, 0.5*dot_product(sep2,c2old)
  !print*, "C2 translated"
  !print*, sep1c2comp

  ! The components of these sums should all have the save order as the coefficent i.e  c0 = scalar, c1 = vector
  c0 = c0old + dot_product(c1old,sep1) + 0.5*inner_product2(sep2,c2old) + 1./6.*inner_product3(sep3,c3)
  c1 = c1old + sep1c2comp + 0.5*sep2c3comp
  !print*, "c1 is: ",c1
  c2 = c2old + sep1c3comp !+ 0.5*inner_product31_to_2(sep1,c3old)
  c3 = c3old

 end subroutine translate_expansion_center

 subroutine accel_at_bodypos(x,com,c0,c1,c2,c3,accel)
  real, intent(inout) :: accel(3)
  real, intent(in) :: x(3),com(3)
  real, intent(in) :: c0,c1(3),c2(3,3),c3(3,3,3)
  real :: sep1(3), sep2(3,3) !, sep3(3,3,3)
  real :: sep1c2comp(3),sep2c3comp(3)

  ! The N-fold outer products of the seperation
  sep1 = x - com 

  ! Replace outer product by matrix mul
  !sep2 = matmul(RESHAPE(sep1,(/3,1/)), RESHAPE(sep1,(/1,3/)))
  !call outer_product1(sep1,sep1,sep2)

  sep1c2comp(1) = sep1(1)*c2(1,1) + sep1(2)*c2(1,2) + sep1(3)*c2(1,3)
  sep1c2comp(2) = sep1(1)*c2(2,1) + sep1(2)*c2(2,2) + sep1(3)*c2(2,3)
  sep1c2comp(3) = sep1(1)*c2(3,1) + sep1(3)*c2(3,2) + sep1(3)*c2(3,3)

  sep2c3comp(1) = sep1(1)*(sep1(1)*c3(1,1,1) + sep1(2)*c3(1,1,2) + sep1(3)*c3(1,1,3)) &
                + sep1(2)*(sep1(1)*c3(1,2,1)+ sep1(2)*c3(1,2,2) + sep1(3)*c3(1,2,3)) &
                + sep1(3)*(sep1(1)*c3(1,3,1)+ sep1(2)*c3(1,3,2) + sep1(3)*c3(1,3,3))

  sep2c3comp(2) = sep1(1)*(sep1(1)*c3(2,1,1) + sep1(2)*c3(2,1,2) + sep1(3)*c3(2,1,3)) &
                + sep1(2)*(sep1(1)*c3(2,2,1)+ sep1(2)*c3(2,2,2) + sep1(3)*c3(2,2,3)) &
                + sep1(3)*(sep1(1)*c3(2,3,1)+ sep1(2)*c3(2,3,2) + sep1(3)*c3(2,3,3))

  sep2c3comp(3) = sep1(1)*(sep1(1)*c3(3,1,1) + sep1(2)*c3(3,1,2) + sep1(3)*c3(3,1,3)) &
                + sep1(2)*(sep1(1)*c3(3,2,1)+ sep1(2)*c3(3,2,2) + sep1(3)*c3(3,2,3)) &
                + sep1(3)*(sep1(1)*c3(3,3,1)+ sep1(2)*c3(3,3,2) + sep1(3)*c3(3,3,3))

  !accel = accel  !+ (c1 &
  !+ inner_product2_to_vec(sep1,c2)) &
 !+ 0.5*inner_product23_to_vec(sep2,c3)
 accel = c1 + sep1c2comp  + 0.5*sep2c3comp !+ 0.5*inner_product23_to_vec(sep2,c3)

 end subroutine accel_at_bodypos
 subroutine poten_at_bodypos(x,com,c0,c1,c2,c3,poten)
  real, intent(inout) :: poten
  real, intent(in) :: x(3),com(3)
  real, intent(in) :: c0,c1(3),c2(3,3),c3(3,3,3)
  real :: sep1(3), sep2(3,3), sep3(3,3,3)

  ! The N-fold outer products of the seperation
  sep1 = x - com 

  ! Replace outer product by matrix mul
  !sep2 = matmul(RESHAPE(sep1,(/3,1/)), RESHAPE(sep1,(/1,3/)))
  call outer_product1(sep1,sep1,sep2)
  call outer_product2(sep2,sep1,sep3)
  !print(matmul(RESHAPE(sep1,(/3,1/)),sep2))

  poten = poten -(c0 + dot_product((x-com),c1) + 0.5*inner_product2(sep2,c2)) !+ 1./6.*inner_product3(sep3,c3))
 end subroutine poten_at_bodypos

 subroutine outer_product1(tens1,tens2,theproduct)
  real, intent(in) :: tens1(3),tens2(3)
  real,intent(out) :: theproduct(3,3)
  integer :: i,j

  do j=1,3
    do i=1,3
      theproduct(i,j) = tens1(i)*tens2(j)
    enddo 
  enddo 

 end subroutine outer_product1

 subroutine outer_product2(tens1,tens2,theproduct)
  real, intent(in) :: tens1(3,3), tens2(3)
  real,intent(out) :: theproduct(3,3,3)
  integer :: i,j,k

  do k=1,3
    do j=1,3
      do i=1,3
        theproduct(i,j,k) = tens1(i,j)*tens2(k)
      enddo 
    enddo 
  enddo 

 end subroutine outer_product2

 real function inner_product2(x2,c2) result(scalar)
 ! Returns the inner product of  a rank 2 tensor 
  real, intent(in) :: x2(3,3), c2(3,3)
  integer :: i, j

  scalar = 0.

  do j=1,3
    do i=1,3
      !print*, i, j 
      scalar = scalar + x2(i,j)*c2(i,j)
    enddo 
  enddo 

 end function inner_product2

 real function inner_product3(x3,c3) result(scalar)
 ! Returns the inner product of a rank 3 tensor 
  real, intent(in) :: x3(3,3,3), c3(3,3,3)
  integer :: i, j, k

  scalar = 0.

  do k=1,3
    do j=1,3
      do i=1,3
        scalar = scalar + x3(i,j,k)*c3(i,j,k)
      enddo
    enddo 
  enddo 

 end function inner_product3

 function inner_product2_to_vec(vec,tens) result(thevector)
  real, intent(in) :: vec(3),tens(3,3)
  real, dimension(3) :: thevector
  integer :: i,j

  thevector = 0.

  do j=1,3
    do i=1,3
      thevector(i) = thevector(i) + vec(j)*tens(j,i)
    enddo 
  enddo 
  !print*,"The vector: "
  !print*, thevector


 end function inner_product2_to_vec

 function inner_product23_to_vec(tens1,tens2) result(thevector)
  real, intent(in) :: tens1(3,3),tens2(3,3,3)
  real, dimension(3) :: thevector
  integer :: i,j,k

  thevector = 0.

  do k=1,3
    do j=1,3
      do i=1,3 
        thevector(i) = thevector(i) + tens1(j,k)*tens2(j,k,i)
      enddo
    enddo 
  enddo 


 end function inner_product23_to_vec

 function inner_product31_to_2(tens1,tens2) result(thetens)
  real, intent(in) :: tens1(3), tens2(3,3,3)
  real, dimension(3,3) :: thetens
  integer :: i,j,k

  thetens = 0.

  do k=1,3
    do j=1,3
      do i=1,3
        thetens(i,j) = thetens(i,j) + tens1(k)* tens2(i,j,k)

      enddo 
    enddo 
  enddo 

 end function inner_product31_to_2
 !----------------------------------------------------------------
!+
!  Internal subroutine to compute the Taylor-series expansion
!  of the gravitational force, given the force acting on the
!  centre of the node and its derivatives
!
! INPUT:
!   fnode: array containing force on node due to distant nodes
!          and first derivatives of f (i.e. Jacobian matrix)
!          and second derivatives of f (i.e. Hessian matrix)
!   dx,dy,dz: offset of the particle from the node centre of mass
!
! OUTPUT:
!   fxi,fyi,fzi : gravitational force at the new position
!+
!----------------------------------------------------------------
subroutine expand_fgrav_in_taylor_series(fnode,dx,dy,dz,fxi,fyi,fzi,poti)
 real, intent(in)  :: fnode(20)
 real, intent(in)  :: dx,dy,dz
 real, intent(out) :: fxi,fyi,fzi,poti
 real :: dfxx,dfxy,dfxz,dfyy,dfyz,dfzz
 real :: d2fxxx,d2fxxy,d2fxxz,d2fxyy,d2fxyz,d2fxzz,d2fyyy,d2fyyz,d2fyzz,d2fzzz

 fxi = fnode(1)
 fyi = fnode(2)
 fzi = fnode(3)
 !print*, "Correct components: "
 !print*, fnode(4),fnode(5),fnode(6),fnode(7),fnode(8),fnode(9)
 dfxx = fnode(4)
 dfxy = fnode(5)
 dfxz = fnode(6)
 dfyy = fnode(7)
 dfyz = fnode(8)
 dfzz = fnode(9)
 d2fxxx = fnode(10)
 d2fxxy = fnode(11)
 d2fxxz = fnode(12)
 d2fxyy = fnode(13)
 d2fxyz = fnode(14)
 d2fxzz = fnode(15)
 d2fyyy = fnode(16)
 d2fyyz = fnode(17)
 d2fyzz = fnode(18)
 d2fzzz = fnode(19)
 poti = fnode(20)

 fxi = fxi + dx*(dfxx + 0.5*(dx*d2fxxx + dy*d2fxxy + dz*d2fxxz)) &
            + dy*(dfxy + 0.5*(dx*d2fxxy + dy*d2fxyy + dz*d2fxyz)) &
            + dz*(dfxz + 0.5*(dx*d2fxxz + dy*d2fxyz + dz*d2fxzz))
  fyi = fyi + dx*(dfxy + 0.5*(dx*d2fxxy + dy*d2fxyy + dz*d2fxyz)) &
            + dy*(dfyy + 0.5*(dx*d2fxyy + dy*d2fyyy + dz*d2fyyz)) &
            + dz*(dfyz + 0.5*(dx*d2fxyz + dy*d2fyyz + dz*d2fyzz))
  fzi = fzi + dx*(dfxz + 0.5*(dx*d2fxxz + dy*d2fxyz + dz*d2fxzz)) &
            + dy*(dfyz + 0.5*(dx*d2fxyz + dy*d2fyyz + dz*d2fyzz)) &
            + dz*(dfzz + 0.5*(dx*d2fxzz + dy*d2fyzz + dz*d2fzzz))
 poti = poti - (dx*fxi + dy*fyi + dz*fzi)

  !print*, "The vector 2"
  !print*, fxi,fyi,fzi
 return
end subroutine expand_fgrav_in_taylor_series
 
end module evaluate