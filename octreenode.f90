
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
    ! set bucket size to 10 
    integer :: data(10)

    ! if leaf node
    ! 0 false, 1 true 
    logical :: isleaf

    ! if body
    logical :: isBody

    ! rmax for MAC
    real :: rmax

    ! taylor series coeff for the node 
    real :: fnode(20)

    end type
    contains
    subroutine new_node(this,size,origin)
     type(octreenode), intent(out) :: this
     real, intent(in) :: size, origin(3)
     this % isleaf = .TRUE.
     this % children(:) = 0
     this % size = size
     write(*,*) "Size set as: ", size
     this % origin = origin
    end subroutine new_node

    subroutine null_node(this)
    type(octreenode), intent(out) :: this
    ! subroutine for uninitialized node
    this % isleaf = .FALSE. 
    this % children(:) = 0
    this % data(:) = 0
    this % size = 0.0
    this % totalmass = 0.0
    end subroutine null_node

    logical function node_full(this) result(flag)
     type(octreenode), intent(in) :: this
      flag = this % data(10) .NE. 0
    end function node_full

    subroutine insert_data(this,data)
     type(octreenode), intent(out) :: this
     integer, intent(in) :: data
     integer :: i

     do i=1, 10
      if (this % data(i) .EQ. 0) then 
       this % data(i) = data 
       EXIT
      end if 
     enddo 
    
    end subroutine insert_data   


end module octreetype