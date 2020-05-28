module evaluate
 use octreetype
 implicit none 
 contains

 recursive subroutine evaluate_gravity(node,nodes,z0,c0,c1,c2,c3,x,poten)
  type(octreenode), intent(inout) :: node, nodes(:)
  real,intent(inout) :: c0,c1(3),c2(3,3),c3(3,3,3),z0(3),x(:,:),poten(:,:)
  type(octreenode) :: childnode
  real :: z1(3),xbod(3),bodpot
  integer :: i,nochild,bodyindex
  !real :: c0old,c1old(3),c2old(3,3),c3old(3,3,3)

  ! TAYLOR SERIES OF CELL A
  ! TA 

  ! Get CoM of current node
   z1 = node % centerofmass

  ! TRANSLATE TAYLOR SERIES T0 TO CENTER OF MASS OF A
  call translate_expansion_center(z0,z1,c0,c1,c2,c3)
  
  ! TA += T0
  ! Sum up the field tensors in a taylor series to get poten
  ! Accumulate field tensors 
  !node % c0 = node % c0 + c0
  !node % c1 = node % c1 + c1
  !node % c2 = node % c2 + c2
  !node % c3 = node % c3 + c3
  ! FOR BODY CHILDREN OF A 

  nochild = size(node % bodychildren)
  do i=1, nochild
    ! Get the index of the body 
    bodyindex = node % bodychildren(i)
    xbod = x(:,bodyindex)



   ! Evaluate TA at body's position
   call poten_at_bodypos(xbod,z1,c0,c1,c2,c3,bodpot)

   ! add to body's potential and acceleration

   poten(:,bodyindex) = poten(:,bodyindex) + bodpot

  enddo 

  ! FOR CHILDREN OF A 

  do i=1, 8

   if (node % children(i) /= 0) then  
    childnode = nodes(node % children(i))
    call evaluate_gravity(childnode,nodes,z0,c0,c1,c2,c3,x,poten)
   endif 

  enddo  

 end subroutine evaluate_gravity


 subroutine translate_expansion_center(z0,z1,c0,c1,c2,c3)
  real, intent(in) :: z0(3), z1(3)
  real, intent(inout) :: c0, c1(3), c2(3,3), c3(3,3,3)
  real :: c0old, c1old(3), c2old(3,3), c3old(3,3,3)
  real :: sep1(3), sep2(3,3), sep3(3,3,3)
  
  c0old = c0
  c1old = c1
  c2old = c2
  c3old = c3

  sep1 = z0-z1
  call outer_product1(sep1,sep1,sep2)
  call outer_product2(sep2,sep1,sep3)

  ! The components of these sums should all have the save order as the coefficent i.e  c0 = scalar, c1 = vector
  c0 = c0old + dot_product(c1old,sep1) + 0.5*inner_product2(sep2,c2old) + 1./6.*inner_product3(sep3,c3)
  c1 = c1old + inner_product2_to_vec(sep1,c2old) + 0.5*inner_product23_to_vec(sep2,c3old)
  c2 = c2old + inner_product31_to_2(sep1,c3old)
  c3 = c3old

 end subroutine translate_expansion_center

 subroutine poten_at_bodypos(x,com,c0,c1,c2,c3,poten)
  real, intent(inout) :: poten
  real, intent(in) :: x(3),com(3)
  real, intent(in) :: c0,c1(3),c2(3,3),c3(3,3,3)
  real :: sep1(3), sep2(3,3), sep3(3,3,3)

  ! The N-fold outer products of the seperation
  sep1 = x-com

  ! Replace outer product by matrix mul
  sep2 = matmul(RESHAPE(sep1,(/3,1/)), RESHAPE(sep1,(/1,3/)))
  !call outer_product1(sep1,sep1,sep2)
  call outer_product2(sep2,sep1,sep3)
  !print(matmul(RESHAPE(sep1,(/3,1/)),sep2))

  poten = poten -(c0 + dot_product((x-com),c1) + 0.5*inner_product2(sep2,c2) + 1./6.*inner_product3(sep3,c3))
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
      thevector(i) = thevector(i) + vec(i)*tens(i,j)
    enddo 
  enddo 


 end function inner_product2_to_vec

 function inner_product23_to_vec(tens1,tens2) result(thevector)
  real, intent(in) :: tens1(3,3),tens2(3,3,3)
  real, dimension(3) :: thevector
  integer :: i,j,k

  thevector = 0.

  do k=1,3
    do j=1,3
      do i=1,3 
        thevector(i) = thevector(i) + tens1(i,j)*tens2(i,j,k)
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
pure subroutine expand_fgrav_in_taylor_series(fnode,dx,dy,dz,fxi,fyi,fzi,poti)
 real, intent(in)  :: fnode(20)
 real, intent(in)  :: dx,dy,dz
 real, intent(out) :: fxi,fyi,fzi,poti
 real :: dfxx,dfxy,dfxz,dfyy,dfyz,dfzz
 real :: d2fxxx,d2fxxy,d2fxxz,d2fxyy,d2fxyz,d2fxzz,d2fyyy,d2fyyz,d2fyzz,d2fzzz

 fxi = fnode(1)
 fyi = fnode(2)
 fzi = fnode(3)
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

 return
end subroutine expand_fgrav_in_taylor_series
 
end module evaluate