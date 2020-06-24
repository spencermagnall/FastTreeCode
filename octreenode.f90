
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

    ! Coeff 0, i.e potential
    real :: c0

    ! Coeff 1, i.e accel, a vector 
    real :: c1(3)

    ! Coeff 2, rank2 tensor
    real :: c2(3,3)

    ! Coeff 3, rank3 tensor
    real :: c3(3,3,3)

    ! The quadrupole moment of this node 
    real :: quads(3,3)


    ! Body children for the node
    ! I've just set this to a large value for now 
    ! But may need to do some dynamic array allocation
    integer,  dimension(:),allocatable :: bodychildren
    ! points to the index of the last inserted body child
    integer :: bodychildpont

    end type
    contains
    subroutine new_node(this,size,origin)
     type(octreenode), intent(inout) :: this
     real, intent(in) :: size, origin(3)
     this % isleaf = .TRUE.
     this % children(:) = 0
     if (.not. allocated(this % bodychildren)) allocate ( this % bodychildren(100) )
     this % bodychildren(:) = 0
     this % bodychildpont = 0
     this % size = size
     this % c0 = 0
     this % c1 = 0
     this % c2 = 0 
     this % c3 = 0
     write(*,*) "Size set as: ", size
     this % origin = origin
     this % quads = 0.
    end subroutine new_node

    subroutine null_node(this)
    type(octreenode), intent(inout) :: this
    ! subroutine for uninitialized node
    this % isleaf = .FALSE. 
    this % children(:) = 0
    this % data(:) = 0
    !deallocate(this % bodychildren)
    if (.not. allocated(this % bodychildren)) allocate ( this % bodychildren(100) )
    this % bodychildren(:) = 0
    this % bodychildpont = 0
    this % size = 0.0
    this % totalmass = 0.0
    this % c0 = 0
    this % c1 = 0
    this % c2 = 0 
    this % c3 = 0
    this % quads = 0. 
    end subroutine null_node

    logical function node_full(this) result(flag)
     type(octreenode), intent(in) :: this
      flag = this % data(10) .NE. 0
    end function node_full

    subroutine insert_data(this,data)
     type(octreenode), intent(inout) :: this
     integer, intent(in) :: data
     integer :: i,bodyindex

     ! insert into body children as well
     ! This simplfiies the process of spliting a leaf node
     print*, "bodychilpont: ", this % bodychildpont + 1
     this % bodychildpont = this % bodychildpont + 1
     print*,"Size: ", size(this % bodychildren)

     bodyindex = this % bodychildpont
     print*, "Bodychildpont: ", this % bodychildpont
     print*, "bodychildren: ", (this % bodychildren)

    this % bodychildren(bodyindex) = data
     print*, "Crashing here!! "

     ! Find an empty space for data
     do i=1, 10
        print*, i
      if (this % data(i) .EQ. 0) then 
       this % data(i) = data 
       EXIT
      end if 
     enddo 

     print*, "Insert finished"
    
    end subroutine insert_data   

    subroutine resize_bodychildren(this)
        type(octreenode), intent(inout) :: this
        integer, allocatable :: oldchildcopy(:)
        integer :: sizeofarray,i

        ! get the size of the current array 
        sizeofarray = size(this % bodychildren)
        ! allocate a copy 
        allocate(oldchildcopy(sizeofarray))
        ! copy array 
        oldchildcopy = this % bodychildren
        ! deallocate old array 
        deallocate(this % bodychildren)
        ! Double the size of the array
        sizeofarray = sizeofarray * 2
        ! Allocate new larger array 
        allocate(this % bodychildren(sizeofarray)) 
        this % bodychildren = 0.
        do i=1, this % bodychildpont
            if (oldchildcopy(i) /= 0 ) then
                this % bodychildren(i) = oldchildcopy(i)
            endif 
        enddo 
        !this % bodychildren = this % bodychildren + oldchildcopy


        !if this 
    end subroutine resize_bodychildren


end module octreetype