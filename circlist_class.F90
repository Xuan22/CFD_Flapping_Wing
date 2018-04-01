module circlist_class

    implicit none

    type node
        private
        type(node), pointer :: next
        real :: val
    end type node

    type circlist
        private
        integer :: n
        type(node), pointer :: curr
        type(node), pointer :: first
    end type circlist

    interface new
        module procedure circlist_new
    end interface new

    interface del
        module procedure circlist_del
    end interface del

    interface size
        module procedure circlist_size
    end interface size

    interface append
        module procedure circlist_append
    end interface append

    interface first
        module procedure circlist_first
    end interface first

    interface next
        module procedure circlist_next
    end interface next

    interface reset
        module procedure circlist_reset
    end interface reset

    interface valadvance
        module procedure circlist_valadvance
    end interface valadvance


contains

subroutine circlist_new(self, valarray)
    !Constructor for circlist class

    type(circlist), intent(out) :: self
    real, dimension(:), intent(in), optional :: valarray

    ! local var
    integer :: j

    self%n = 0
    self%curr => null()
    self%first => null()

    if(present(valarray))then
        do j = 1, size(valarray)
            call circlist_append(self, valarray(j))
        end do
    end if

end subroutine circlist_new

subroutine circlist_del(self)
! Circlist destructor

    type(circlist), intent(inout) :: self

    ! Break loop...
    self%curr => self%first%next
    self%first%next => null()

    ! Deallocate nodes...
    do while(associated(self%curr%next))
        self%first => self%curr%next
        deallocate(self%curr)
        self%curr => self%first
    end do

    deallocate(self%curr)
    self%curr => null()
    self%first => null()

    ! Update node counter...
    self%n = 0

end subroutine circlist_del

pure function circlist_size(self)
! Returns the number of nodes in the list
    type(circlist), intent(in) :: self
    integer circlist_size
    circlist_size = self%n
end function circlist_size

subroutine circlist_append(self, val)
! Appends a new node into the circularly linked list. The new value
! appended between the current element pointed to by 'curr' and the
! next element pointed to by 'curr%next'. That is,
!
!              (new node)
!                [node] 
!                  | 
!     [ node ] ----x---- [ node ]
!      (curr)           (curr%next)
!

    type(circlist), intent(inout) :: self
    real, intent(in) :: val

    !local var
    type(node), pointer :: ptr

    ! Create new node...
    allocate(ptr)
    ptr%val = val

    if(associated(self%curr))then

        ! list is non-empty
        
        ! insert new node between curr and curr%next ...
        ptr%next => self%curr%next
        self%curr%next => ptr

    else
        ! list is emtpy... alloc first element

        ptr%next => ptr
        self%curr => ptr
        self%first => ptr

    end if

    ! Advance curr to next...
    self%curr => self%curr%next

    ! update node count...
    self%n = self%n + 1

    ! nullify temporary pointer (technically this shouldn't be needed)
    ptr => null()

end subroutine circlist_append

function circlist_first(self)
    ! Sets the curr to the first allocated node and returns its value

    type(circlist), intent(inout) :: self
    real :: circlist_first

    self%curr=>self%first
    circlist_first = self%curr%val

end function circlist_first

function circlist_next(self)
    ! Sets the curr to the next node and returns its value

    type(circlist), intent(inout) :: self
    real :: circlist_next

    self%curr=>self%curr%next
    circlist_next = self%curr%val

end function circlist_next

subroutine circlist_reset(self)
! Sets curr to first allocated node in list
    type(circlist), intent(inout) :: self
    self%curr => self%first
end subroutine circlist_reset

function circlist_valadvance(self)
    ! Returns the value of the current node, then advances curr to the next
    ! node.

    type(circlist), intent(inout) :: self
    real :: circlist_valadvance
    
    circlist_valadvance = self%curr%val
    self%curr => self%curr%next

end function circlist_valadvance

end module circlist_class

#ifdef TEST_CIRCLIST

!--------------------------------------------------------------------------
!                             Test program
!--------------------------------------------------------------------------

program test_circlist

    use circlist_class

    implicit none

    type(circlist) :: c

    integer, parameter :: n = 10
    integer :: i
    real :: val
    real, dimension(n) :: valarray

    ! Create circularly linked list...
    call init(c)

    ! populate list...
    print*, "Populating list..."
    do i = 1, n
        val = real(i)
        print*, "i, val", i, val
        call append(c, val)
    end do
    print*, "...done."

    ! inquire size of list...
    print*, "size of list = ", size(c)

    ! print out list twice...
    print*, "Printing contents of list twice..."
    val = first(c)
    do i = 1, 2*n
        print*, "i, val", i, val
        val = next(c)
    end do
    print*, "...done."

    ! print out list twist using reset and valadvance...
    print*, "Printing contents of list twist with reset and val advance..."
    call reset(c)
    do i = 1,n
        print*, "i, val", i, valadvance(c)
    end do
    do i = 1, n/2
        print*, "i, val", i, valadvance(c)
    end do
    print*, "call reset(c)"
    call reset(c)
    do i = n/2+1, n
        print*, "i, val", i, valadvance(c)
    end do
    print*, "...done."


    ! deallocate list..
    print*, "Deallocating list..."
    call del(c)
    print*, "...done."

    ! Now try initializing from an array...

    do i = 1, n
        valarray(i) = real(i)
    end do

    print*, "Re-allocating list using valarray..."
    call init(c, valarray)
    print*, "...done."

    print*, "Printing contents of list..."
    val = first(c)
    do i = 1, n
        print*, "i, valarray(i), val", i, valarray(i),  val
        val = next(c)
    end do
    print*, "...done."

    print*, "Deallocating list..."
    call del(c)
    print*, "...done."

end program test_circlist

#endif

