module evaluate
 use octreetype
 implicit none 

type evaluate_stack_data

 type(octreenode) :: node 
 real :: z0(3)
 real :: c0,c1(3),c2(3,3),c3(3,3,3)

end type 

contains 

 subroutine evaluate_gravity_stack(nodes,x,accel)
   type(octreenode), intent(inout) :: nodes(:)
   real, intent(in) :: x(:,:)
   real, intent(inout) :: accel(:,:)
   integer :: stacksize, top 
   type(evaluate_stack_data) :: stack(100000)
   type(octreenode) :: currentnode,childnode
   real :: c0,c1(3),c2(3,3),c3(3,3,3),z0(3),z1(3)
   real :: bodaccel(3),xbod(3)
   integer :: bodyindex,i,iter

   top = 1
   stacksize = 1
   iter = 1

   ! Push the root node 
   stack(1) % node = nodes(1)
   stack(1) % z0 = 0.
   stack(1) % c0 = 0.
   stack(1) % c1 = 0.
   stack(1) % c2 = 0.
   stack(1) % c3 = 0.

   do while (top /= 0.)

    print*, "The current iteration is: ",iter
    ! POP ITEM FROM STACK 
    currentnode = stack(top) % node
    z0 = stack(top) % z0
    c0 = stack(top) % c0
    c1 = stack(top) % c1
    c2 = stack(top) % c2 
    c3 = stack(top) % c3 

    top = top - 1 

    ! Perform evaluate on popped item 

   z1 = currentnode % centerofmass
  

  ! TRANSLATE TAYLOR SERIES T0 TO CENTER OF MASS OF A
  
  call translate_expansion_center(z0,z1,c0,c1,c2,c3)

  c0 = currentnode % c0 + c0
  !c0new = c0
  !print*, "C0: "
  !print*,c0new
  c1 = currentnode % c1 + c1
  !c1new = c1 
  !print*, "C1: "
  !print*,c1new
  !c2new = c2
  c2 = currentnode % c2 + c2
  !print*, "C2: "
  !print*,c2new
  c3 = currentnode % c3 + c3
  

  !nochild = node % bodychildpont
  if (currentnode % isLeaf) then
  do i=1, 10
    bodyindex = currentnode % data(i)
    xbod = x(:,bodyindex)
   bodaccel = 0.
   call accel_at_bodypos(xbod,z1,c0,c1,c2,c3,bodaccel)
   accel(:,bodyindex) = accel(:,bodyindex) + bodaccel
  enddo 
 else 

  ! FOR CHILDREN OF A 
  ! change center of mass
  z0 = z1
  
  do i=1, 8

   if (currentnode % children(i) /= 0 .and. c0 /= 0.) then  
    print*, "Childnode index: ", currentnode % children(i)
    childnode = nodes(currentnode % children(i))
    !nodeindex = node % children(i)
    !print*, "Node index: ",nodeindex
    !call evaluate_gravity(childnode,nodes,z0new,c0new,c1new,c2new,c3new,x,accel)
    top = top + 1 

    ! Push Children onto stack 
    stack(top) % node = childnode
    stack(top) % z0 = z0
    stack(top) % c0 = c0
    stack(top) % c1 = c1 
    stack(top) % c2 = c2
    stack(top) % c3  = c3 
   endif 

  enddo
 

  endif 

    


  enddo 


 end subroutine evaluate_gravity_stack

 recursive subroutine evaluate_gravity(node,nodes,z0,c0,c1,c2,c3,x,accel)
  type(octreenode), intent(inout) :: node, nodes(:)
  real,intent(inout) :: c0,c1(3),c2(3,3),c3(3,3,3),z0(3),accel(:,:),x(:,:)!,c0total(:),c1(3,:),c2(3,3,:),c3(3,3,3,:)
  type(octreenode) :: childnode
  real :: z1(3),xbod(3),bodpot,bodaccel(3),accelbef(3)
  integer :: i,nochild,bodyindex
  real :: c0new, c1new(3), c2new(3,3),c3new(3,3,3),z0new(3)
  integer  :: nodeindex,nthreads
  !real :: c0old,c1old(3),c2old(3,3),c3old(3,3,3)

  bodpot = 0.
  bodaccel = 0.
  bodyindex = 0.

  c0new = 0.
  c1new = 0.
  c2new = 0.
  c3new = 0.

   print*, "Center of mass (old): ", z0
  ! Get CoM of current node
   z1 = node % centerofmass
   print*, "Center of Mass new: "
   print*, z1

   print*, "Node mass: ", node % totalmass
   print*, "c0 is: ", c0 

  ! TRANSLATE TAYLOR SERIES T0 TO CENTER OF MASS OF A
  
  if (c0 /= 0. .and. node % totalmass /= 0.) call translate_expansion_center(z0,z1,c0,c1,c2,c3)

  !dr = z0-z1

  !call expand_fgrav_in_taylor_series()
  
  ! TA += T0
  ! Sum up the field tensors in a taylor series to get poten
  ! Accumulate field tensors 

  !print*, "Node taylor coeffs: "
  !print*, "c0 node: "
  !print*, node % c0 
  print*, "c1 : "
  print*, c1
  !print*, "c2 node: "
  !print*, node % c2
  !print*, "c3 node: "
  !print*, node % c3


  if (node % c0 /= 0.) then
  c0new = node % c0 + c0
  !c0new = c0
  !print*, "C0: "
  !print*,c0new
  c1new = node % c1 + c1
  !c1new = c1 
  !print*, "C1: "
  !print*,c1new
  !c2new = c2
  c2new = node % c2 + c2
  !print*, "C2: "
  !print*,c2new
  c3new = node % c3 + c3
  !c3new = c3 
  !print*, "C3: "
  !print*, c3new
  ! FOR BODY CHILDREN OF A 

  else
    c0new = c0
    c1new = c1
    c2new = c2
    c3new = c3 
  endif 
 !STOP

  !nochild = node % bodychildpont
  if (node % isLeaf) then
  do i=1, 10
    
    !if (node%data(i)==0) then
    !  EXIT
    !endif 
    ! Get the index of the body 
    !bodyindex = node % bodychildren(i)

    if (node % data(i) /= 0 .and. c0new /= 0.) then
    print*, node % data(i) 
    bodyindex = node % data(i)
    xbod = x(:,bodyindex)



   ! Evaluate TA at body's position
   !call poten_at_bodypos(xbod,z1,c0,c1,c2,c3,bodpot)
   bodaccel = 0.
   print*, "c0new: ", c0new
   call accel_at_bodypos(xbod,z1,c0new,c1new,c2new,c3new,bodaccel)

   ! add to body's potential and acceleration

   !poten(bodyindex) = poten(bodyindex) + bodpot
   print*, "accel bef:",accel(:,bodyindex)
   !accelbef = accel(:,bodyindex)
   accel(:,bodyindex) = accel(:,bodyindex) + bodaccel
   print*, "Accel now: ", accel(:,bodyindex)
   !print*, "bodaccel: ", bodaccel
   !print*, "delta accel:"
   !print*, accel(:,bodyindex) - accelbef


   endif 


  enddo 
 else 

  ! FOR CHILDREN OF A 
  ! change center of mass
  z0new = z1
  !!$OMP DO 
  do i=1, 8

  if (node % children(i) /= 0 .and. norm2(nodes(node %children(i)) % centerofmass) /= 0. ) then  
    print*, "Childnode index: ", node % children(i)
    childnode = nodes(node % children(i))
    !nodeindex = node % children(i)
    !print*, "Node index: ",nodeindex
    call evaluate_gravity(childnode,nodes,z0new,c0new,c1new,c2new,c3new,x,accel)
   endif 

  enddo
  !!$OMP ENDDO   

  endif 

 end subroutine evaluate_gravity




 subroutine translate_expansion_center(z0,z1,c0,c1,c2,c3)
  real, intent(in) :: z0(3), z1(3)
  real, intent(inout) :: c0, c1(3), c2(3,3), c3(3,3,3)
  real :: c0old, c1old(3), c2old(3,3), c3old(3,3,3)
  real :: sep1(3), sep2(3,3), sep3(3,3,3)
  real :: sep1c2comp(3)
  
  c0old = c0
  c1old = c1
  c2old = c2
  c3old = c3
  !print*, "c0"
  !print*, c0old
  !print*, "c1"
  !print*, c1old
  !print*, "c2"
  !print*, c2old
  !print*, "c3"
  !print*, c3old(2,1,1)

  sep1 = z1-z0
  call outer_product1(sep1,sep1,sep2)
  call outer_product2(sep2,sep1,sep3)

  print*, c2 

  sep1c2comp(1) = sep1(1)*c2(1,1) + sep1(2)*c2(1,2) + sep1(3)*c2(1,3)
  sep1c2comp(2) = sep1(1)*c2(2,1) + sep1(2)*c2(2,2) + sep1(3)*c2(2,3)
  sep1c2comp(3) = sep1(1)*c2(3,1) + sep1(2)*c2(3,2) + sep1(3)*c2(3,3)

  !print*, "Second term value: "
  !print*, 0.5*inner_product2(sep2,c2old)
  !print*, 0.5*dot_product(sep2,c2old)
  print*, "C2 translated"
  print*, sep1c2comp

  ! The components of these sums should all have the save order as the coefficent i.e  c0 = scalar, c1 = vector
  c0 = c0old + dot_product(c1old,sep1) + 0.5*inner_product2(sep2,c2old) !+ 1./6.*inner_product3(sep3,c3)
  c1 = c1old + sep1c2comp!+ inner_product2_to_vec(sep1,c2) + 0.5*inner_product23_to_vec(sep2,c3old)
  c2 = c2old !+ inner_product31_to_2(sep1,c3old)
  c3 = c3old

 end subroutine translate_expansion_center

 subroutine accel_at_bodypos(x,com,c0,c1,c2,c3,accel)
  real, intent(inout) :: accel(3)
  real, intent(in) :: x(3),com(3)
  real, intent(in) :: c0,c1(3),c2(3,3),c3(3,3,3)
  real :: sep1(3), sep2(3,3) !, sep3(3,3,3)
  real :: sep1c2comp(3),sep2c3comp

  ! The N-fold outer products of the seperation
  sep1 = x-com 

  ! Replace outer product by matrix mul
  !sep2 = matmul(RESHAPE(sep1,(/3,1/)), RESHAPE(sep1,(/1,3/)))
  call outer_product1(sep1,sep1,sep2)

  sep1c2comp(1) = sep1(1)*c2(1,1) + sep1(2)*c2(1,2) + sep1(3)*c2(1,3)
  sep1c2comp(2) = sep1(1)*c2(2,1) + sep1(2)*c2(2,2) + sep1(3)*c2(2,3)
  sep1c2comp(3) = sep1(1)*c2(3,1) + sep1(2)*c2(3,2) + sep1(3)*c2(3,3)

  !accel = accel  !+ (c1 &
  !+ inner_product2_to_vec(sep1,c2)) &
 !+ 0.5*inner_product23_to_vec(sep2,c3)
 accel = accel + c1 + sep1c2comp !+ 0.5*inner_product23_to_vec(sep2,c3)

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
        thetens(i,j) = thetens(i,j) + tens2(i,j,k)*tens1(k)

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
 print*, "Correct components: "
 print*, fnode(4),fnode(5),fnode(6),fnode(7),fnode(8),fnode(9)
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

 fxi = fxi + dx*dfxx & 
           + dy*dfxy & 
           + dz*dfxz  
 fyi = fyi + dx*dfxy & 
           + dy*dfyy & 
           + dz*dfyz 
 fzi = fzi + dx*dfxz & 
           + dy*dfyz & 
           + dz*dfzz  
 poti = poti - (dx*fxi + dy*fyi + dz*fzi)

  !print*, "The vector 2"
  !print*, fxi,fyi,fzi
 return
end subroutine expand_fgrav_in_taylor_series
 
end module evaluate