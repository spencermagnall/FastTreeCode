module evaluate
 use octreetype
 implicit none 
 contains

 recursive subroutine evaluate_gravity(node,nodes,fnode)
  type(octreenode), intent(inout) :: node, nodes(:)
  real,intent(inout) :: fnode(20)
  type(octreenode) :: childnode 
  ! TAYLOR SERIES OF CELL A
  ! TA 

  ! TRANSLATE TAYLOR SERIES T0 TO CENTER OF MASS OF A
  
  ! TA += T0

  ! FOR BODY CHILDREN OF A 

  ! FOR CHILDREN OF A 

 end subroutine evaluate_gravity

end module evaluate