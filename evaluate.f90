module evaluate
 use octreetype
 implicit none 
 contains

 recursive subroutine evaluate_gravity(node,nodes,fnode)
  type(octreenode), intent(inout) :: node, nodes(:)
  real,intent(inout) :: fnode(20)
  type(octreenode) :: childnode
  integer :: i,nochild

  ! TAYLOR SERIES OF CELL A
  ! TA 

  ! TRANSLATE TAYLOR SERIES T0 TO CENTER OF MASS OF A
  
  ! TA += T0

  ! FOR BODY CHILDREN OF A 

  nochild = size(node % bodychildren)
  do i=1, nochild
   ! Evaluate TA at body's position

   ! add to body's potential and acceleration

  enddo 

  ! FOR CHILDREN OF A 

  do i=1, 8

   if (node % children(i) /= 0) then  
    childnode = nodes(node % children(i))
    call evaluate_gravity(childnode,nodes,fnode)
   endif 

  enddo  

 end subroutine evaluate_gravity

end module evaluate