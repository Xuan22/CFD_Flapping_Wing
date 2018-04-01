module bcparam
  
  implicit none

  !Maximum number of BCs per grid
  integer, parameter, private :: MAX_BC=25

  type bc_t
     integer :: nbc ! Total number of BCs for a mesh
    
     !Type of BC and direction
     integer, dimension(:), allocatable :: ibtyp, ibdir

     !Start and end points in all 3 directions 
     integer, dimension(:), allocatable :: jbcs, jbce, kbcs, kbce, &
          & lbcs, lbce
     
     !Processor flag for exchanging data across processors
     integer, dimension(:), allocatable :: ibproc
  end type bc_t

contains 
  subroutine init_bc(this,size)
    !! Constructor routine to initialize the arrays in the grid bc
    !! object
    type(bc_t), intent(inout) :: this
    integer, intent(in) :: size

    this%nbc=size
    
    allocate(this%ibtyp(size),this%ibdir(size))
    allocate(this%jbcs(size),this%jbce(size))
    allocate(this%kbcs(size),this%kbce(size))
    allocate(this%lbcs(size),this%lbce(size))
    allocate(this%ibproc(size))
  end subroutine init_bc
  
  subroutine del_bc(this)
    !! Destructor routine to deallocate arrays in grid bc object
    
    type(bc_t), intent(inout) :: this
    
    this%nbc=-1 
    deallocate(this%ibtyp,this%ibdir)
    deallocate(this%jbcs,this%jbce)
    deallocate(this%kbcs,this%kbce)
    deallocate(this%lbcs,this%lbce)
    deallocate(this%ibproc)
  end subroutine del_bc
  
  subroutine read_bc(this,fid,flnam)
    !! Read BC information for a single mesh from input file
    !!
    !! Inputs: 
    !!    this - Grid boundary conditions object
    !!    fid - File unit to read BCs from
    !!
 
    type(bc_t), intent(inout) :: this
    character(len=*), intent(in) :: flnam
    integer, intent(in) :: fid

    integer :: nbc, i
    integer, dimension(MAX_BC) :: jbcs, jbce, kbcs, kbce, lbcs, lbce
    integer, dimension(MAX_BC) :: ibtyp, ibdir, ibproc

    namelist/bcinp/ nbc,ibtyp,ibdir,jbcs,jbce,kbcs,kbce,lbcs,lbce,&
         &ibproc

    read(fid,bcinp)

    if (nbc .gt. MAX_BC) then 
       write(*,*) "Error reading BC in file: "//flnam
       stop
    end if
    
    call init_bc(this,nbc)
    
    do i=1,nbc
       this%ibtyp(i)=ibtyp(i)
       this%ibdir(i)=ibdir(i)
       this%jbcs(i)=jbcs(i)
       this%jbce(i)=jbce(i)
       this%kbcs(i)=kbcs(i)
       this%kbce(i)=kbce(i)
       this%lbcs(i)=lbcs(i)
       this%lbce(i)=lbce(i)
       this%ibproc(i)=ibproc(i)
    end do
  end subroutine read_bc
     
end module bcparam
