
module octreetype
    type octreenode

    ! the size of the node
    real :: size

    ! total mass of the node 
    real :: totalmass

    ! the index of the parent node
    ! root node points to 0
    integer :: parent

    ! the index children of the node  
    integer :: children(8)

    ! the center of the node
    real :: origin(3)

    ! the center of mass of the node
    real :: centerofmass(3) 

    ! the index of the data
    ! 0 = no data 
    integer :: data

    ! if leaf node
    ! 0 false, 1 true 
    logical :: isleaf

    end type
    contains
    subroutine new_node(this,size,origin)
     type(octreenode), intent(out) :: this
     real, intent(in) :: size, origin(3)
     this % isleaf = .TRUE.
     this % children(:) = 0
     this % size = size
     this % origin = origin
    end subroutine new_node

    subroutine null_node(this)
    type(octreenode), intent(out) :: this
    ! subroutine for uninitialized node
    this % isleaf = .FALSE. 
    this % children(:) = 0
    this % data = 0
    end subroutine null_node
end module octreetype