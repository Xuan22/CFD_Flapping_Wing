module mesh_types
!!! This module contains the data types used to represent meshes,
!!! partitioned mesh blocks and their BCs and their relationships in a
!!! computational domain. These data types provide the global picture
!!! so that the flow computational routines can be concerned with just
!!! a specific computational block. The data structures are then used
!!! by connectivity routines to determine the flow of information
!!! between the mesh blocks.
   
   use kind_param
   use ioutils
   use string_utils
   use TURNS_options

   implicit none

   !! Length of names used in data types
   integer, parameter, private :: MAX_NAME_LEN=32

   type simple_bc
      !! Datatype to represent a simple BC condition, i.e., one that
      !! does not depend on data exchange with another face.

      ! bc_name - Name of the BC 
      character(len=MAX_NAME_LEN) :: bc_name

      ! ibtyp - Type code for the BC (see TURNS_options)
      ! ibdir - Direction of the BC 
      integer :: ibtyp, ibdir

      ! Jbcs, Kbcs, Lbcs - Starting index for each dimension
      ! Jbce, Kbce, Lbce - Ending index for each dimension
      integer :: jbcs, jbce, kbcs, kbce, lbcs, lbce
   end type simple_bc

   type interface_bc
      !! Datatype to represent an interface (1-to-1 connectivity)
      !! information. Used to represent surfaces generated as part of
      !! mesh splitting as well as mesh faces in the wake region, or
      !! C-O mesh interfaces.

      ! bc_name - Name of the BC 
      character (len=MAX_NAME_LEN) :: bc_name

      ! ibtyp - Type code of the BC
      ! ibdir - Direction in which BC is applied
      integer :: ibtyp, ibdir

      ! Jbcs, Kbcs, Lbcs - Starting index for each dimension
      ! Jbce, Kbce, Lbce - Ending index for each dimension
      integer :: jbcs, jbce, kbcs, kbce, lbcs, lbce

      ! donor_mesh - Name of the donor mesh
      ! donorid - Donor mesh's block ID.
      ! didir - Donor's direction
      character(len=MAX_NAME_LEN) :: donor_mesh
      integer :: donorid, didir

      ! dJs, dKs, dLs - Donor surface starting indices
      ! dJe, dKe, dLe - Donor surface end index
      integer :: djs, dje, dks, dke, dls, dle
   end type interface_bc

   type mesh_info
      ! Default mesh information type to store physical domain information

      ! Name of the mesh (used in I/O) 
      character(len=MAX_NAME_LEN) :: mesh_name

      ! mesh_id - Unique internal ID for tracking meshes
      ! dom_type - Physical domain type represented by the mesh
      ! jm, km, lm - Max indices of the mesh in 3 directions
      ! nbc_sim - No. of simple BC conditions
      ! nbc_int - No. of interface BC conditions
      ! ncells - Number of cells
      ! nfaces - Number of faces
      integer :: mesh_id, dom_type
      integer :: jm, km, lm
      integer :: nbc_sim, nbc_int
      integer :: ncells, nfaces

      ! bcinfo - Array of simple BC objects
      ! ifaceinfo - Array of interface BC objects
      type(simple_bc), dimension(:), pointer :: bcinfo
      type(interface_bc), dimension(:), pointer :: ifaceinfo

      ! visc_dir - Track if any direction has a viscous BC
      logical, dimension(3) :: visc_dir

      ! grid_file - File to read the grid from
      character(len=FLNLEN) :: grid_file

      ! nblocks - No. of computational blocks corresponding to this
      !           mesh
      ! start(3,nblocks) - Starting indices for the computational blocks.
      ! end(3,nblocks) - Ending indices for the computational blocks. 
      integer :: nblocks
      integer, dimension(:), pointer :: block_id
      integer, dimension(:,:), pointer :: start, end

      ! master_id - ID of the processor that will act as the MASTER
      ! for this mesh
      integer :: master_id
   end type mesh_info

   type mesh_block
      ! Object representing a computational block

      ! block_name - Name of the block (used in I/O)
      ! mesh_id - The global ID of the physical mesh to which this
      !           block belongs to
      ! dom_type - Physical object that this mesh block represents
      ! block_id - Unique ID representing this mesh block
      ! jm, km, lm - Max indices of this block

      ! gjs, gks, gls - Starting indices in global mesh where this
      !                 block is located.
      ! gje, gke, gle - Ending indices in global mesh
      character(len=MAX_NAME_LEN) :: block_name
      integer :: mesh_id, dom_type, block_id
      integer :: jm, km, lm
      integer :: gjs, gje, gks, gke, gls, gle
	  
      ! nbc_sim - Number of simple BCs in this mesh block
      ! nbc_int - Number of interface BCs in this mesh block
      ! bcinfo - Simple BC object array (size=nbc_sim)
      ! ifaceinfo - Interface BC object array (size=nbc_int)
      integer :: nbc_sim, nbc_int
      type(simple_bc), dimension(:), pointer :: bcinfo
      type(interface_bc), dimension(:), pointer :: ifaceinfo

      ! is_master - Is this block the master proc for the mesh?
      logical :: is_master
   end type mesh_block

   type overset_info
      ! nmesh - Number of meshes overlapping each other
      ! mesh_ids - Global mesh IDs of these meshes
      ! master_id - ID of the processor that will act as master for
      !             this set of overlapping meshes in connectivity
      !             etc.
      integer :: nmesh
      integer, dimension(:), pointer :: mesh_ids
      integer :: master_id
   end type overset_info

   private :: set_mesh_domain_type, read_simple_bc, set_simple_bc_type, &
         read_interface_bc, set_iface_bc_type, set_iface_params
         
contains
   function get_global_mesh_id(self,mname) result (mid)
      !! Given an array of global meshes and a name, return the
      !! corresponding mesh ID
      type(mesh_info), dimension(:), intent(in) :: self
      character(len=*), intent(in) :: mname
      integer :: mid

      character(len=MAX_NAME_LEN) :: mshnam
      integer :: i

      mshnam=upcase(trim(adjustl(mname)))

      mid=0
      do i=1,size(self,1)
         if (self(i)%mesh_name == mshnam) then
            mid=i
            exit
         end if
      end do

      if (mid == 0) call stop_execution('get_global_mesh_id', &
            'Cannot find mesh : '//mshnam)
   end function get_global_mesh_id

   subroutine set_mesh_domain_type(self,mshtyp)
      !! Determine the physical domain ID represented by the mesh 
      type(mesh_info), intent(inout) :: self
      character(len=*), intent(in) :: mshtyp

      select case (upcase(trim(adjustl(mshtyp))))
      case ("BLADE")
         self%dom_type=DOM_BLADE
      case ("WAKE")
         self%dom_type=DOM_WAKE
      case ("GROUNDWAKE")
         self%dom_type=DOM_GROUNDWAKE
      case ("SLAT")
         self%dom_type=DOM_SLAT
      ! Math: add IBC
      case ("IBC")
         self%dom_type=DOM_IBC
      ! Math: adding DOM_FUS
      case ("FUS")
         self%dom_type=DOM_FUS
      case default
         call stop_execution('set_mesh_domain_type', &
               'Unknown mesh domain type: '//mshtyp)
      end select
   end subroutine set_mesh_domain_type

   subroutine read_mesh_info(self,flnam)
      !! Read mesh information from input file and initialize the
      !! global mesh array object.
      type(mesh_info), dimension(:), allocatable, intent(inout) :: self
      character(len=*), optional :: flnam

      integer :: un, nmesh, i
      character(len=RECLEN) :: mshtyp
      
      if (.not. present(flnam)) then
         un=open_file("mesh.inp")
      else
         un=open_file(flnam)
      end if

      ! Number of meshes
      read(un,*) nmesh
      allocate(self(nmesh))

      do i=1,nmesh
         ! Mesh ID for internal tracking
         self(i)%mesh_id=i
         ! Preset visc_dir to false
         self(i)%visc_dir=.false.
         
         ! Name of mesh, mesh type, grid input file
         read(un,*) self(i)%mesh_name, mshtyp, self(i)%grid_file
         call set_mesh_domain_type(self(i),mshtyp)
         
         ! Grid dimensions
         read(un,*) self(i)%jm, self(i)%km, self(i)%lm

         ! BC information
         read(un,*) self(i)%nbc_sim, self(i)%nbc_int
         call read_simple_bc(self(i),un)
         call read_interface_bc(self(i),un)

         ! Calculate inferred parameters
         self(i)%ncells=(self(i)%jm-1)*(self(i)%km-1)*(self(i)%lm-1)
         self(i)%nfaces=self(i)%ncells+ &
                        (self(i)%km-1)*(self(i)%lm-1) + &
                        (self(i)%lm-1)*(self(i)%jm-1) + &
                        (self(i)%jm-1)*(self(i)%km-1)

      end do
      close(un)
		
      ! Set up donor IDs for interface BCs and fix negative indices. 
      call set_iface_params(self,nmesh)
      call check_simple_bc(self,nmesh)
      call check_interface_bc(self,nmesh)
   end subroutine read_mesh_info

   subroutine read_simple_bc(self,un)
      !! For a given mesh, parse the simple BC information
      !! corresponding to that mesh.
      type(mesh_info), intent(inout) :: self
      integer, intent(in) :: un

      integer :: i
      type(simple_bc), pointer :: sbc

      if (self%nbc_sim <= 0) return

      nullify(self%bcinfo)
      allocate(self%bcinfo(self%nbc_sim))
      do i=1,self%nbc_sim
         sbc=>self%bcinfo(i)
         read(un,*) sbc%bc_name, sbc%ibdir, &
               sbc%jbcs, sbc%jbce, &
               sbc%kbcs, sbc%kbce, &
               sbc%lbcs, sbc%lbce

         ! Set the BC ID
         call set_simple_bc_type(sbc)

         ! Fix negative indices
         if (sbc%jbcs < 0) &
               sbc%jbcs=self%jm+sbc%jbcs+1
         if (sbc%jbce < 0) &
               sbc%jbce=self%jm+sbc%jbce+1
         if (sbc%kbcs < 0) &
               sbc%kbcs=self%km+sbc%kbcs+1
         if (sbc%kbce < 0) &
               sbc%kbce=self%km+sbc%kbce+1
         if (sbc%lbcs < 0) &
               sbc%lbcs=self%lm+sbc%lbcs+1
         if (sbc%lbce < 0) &
               sbc%lbce=self%lm+sbc%lbce+1

         if (sbc%ibtyp == BC_WALL) &
            self%visc_dir(abs(sbc%ibdir))=.true.
      end do

   end subroutine read_simple_bc

   subroutine set_simple_bc_type(self)
      !! Determine the ID of the BC type read from the input files. 
      type(simple_bc), intent(inout) :: self

      select case (upcase(trim(adjustl(self%bc_name))))
      case ("INVISCWALL")
         self%ibtyp=BC_INVISCWALL
      case ("WALL")
         self%ibtyp=BC_WALL
      case ("GROUND")
         self%ibtyp=BC_GROUND
      case ("EXTRAPOLATE")
         self%ibtyp=BC_EXTRAPOLATE
      case ("SYMMETRIC")
         self%ibtyp=BC_SYMMETRIC
      !Math: periodic BC
      case ("PERIODIC")
         self%ibtyp=BC_PERIODIC
      case ("AVERAGE")
         self%ibtyp=BC_AVERAGE
      case ("AXISYM")
         self%ibtyp=BC_AXISYM
      case ("FARFIELD")
         self%ibtyp=BC_FARFIELD
      case ("HOVER")
         self%ibtyp=BC_HOVER
      case default
         call stop_execution('set_simple_bc_type', &
               "Unknown BC type: "//self%bc_name)
      end select
   end subroutine set_simple_bc_type

   subroutine read_interface_bc(self,un)
      !! For a given mesh read the interface BC information 
      type(mesh_info), intent(inout) :: self
      integer :: un

      integer :: i
      type(interface_bc), pointer :: ifbc
      
      if (self%nbc_int <= 0 ) return

      nullify(self%ifaceinfo)
      allocate(self%ifaceinfo(self%nbc_int))
      do i=1,self%nbc_int
         ifbc=>self%ifaceinfo(i)

         read(un,*)  ifbc%bc_name, ifbc%ibdir, &
               ifbc%jbcs, ifbc%jbce, &
               ifbc%kbcs, ifbc%kbce, &
               ifbc%lbcs, ifbc%lbce

         ! Determine BC type ID and process it further
         call set_iface_bc_type(ifbc,un,self%mesh_name, &
               self%dom_type)
         
      end do
   end subroutine read_interface_bc

   subroutine set_iface_bc_type(self, un, mesh_name, dom_type)
      !! Set up the BC Type for this given interface BC. This
      !! subroutine does a little more than what the Simple BC
      !! equivalent does. It tries to read/set up donor information
      !! which is not present in the Simple BC types.
      !!
      !! Inputs:
      !!   un - File unit to read BC donor info (if applicable)
      !!   mesh_name - Mesh name for Wake Averaging BC
      !!   dom_type - Physical domain type ID for error checks
      !!
      !! Modifies:
      !!   self - Interface BC object
      !!
      type(interface_bc), intent(inout) :: self
      integer, intent(in) :: un, dom_type
      character(len=MAX_NAME_LEN), intent(in) :: mesh_name

      select case (upcase(trim(adjustl(self%bc_name))))
      case ("INTERFACE")
         ! A general interface so we need to read the donor information
         self%ibtyp=BC_INTERFACE
         read(un,*) self%donor_mesh, self%didir, &
               self%djs, self%dje, &
               self%dks, self%dke, &
               self%dls, self%dle
      case ("WAKEAVG")
         ! Wake Averaging (so we know this is a C-H/C-O style mesh, we
         ! can set up the values ourselves

         ! Sanity check, we cannot guess for other domain types
         if ((dom_type /= DOM_BLADE) .and. &
             (dom_type /= DOM_SLAT) .and. &
             ! Math: adding DOM_FUS (wake average bc ok for fuselage mesh)
             (dom_type /= DOM_FUS)) then
            call stop_execution('set_iface_bc_type', &
                  "WAKE Averaging BC not defined for "// &
                  trim(adjustl(mesh_name)))
         end if
         
         self%ibtyp=BC_WAKEAVG
         self%donor_mesh=mesh_name

         ! Set the donor information
         self%didir=self%ibdir
         self%djs=-self%jbcs
         self%dje=-self%jbce
         self%dks=self%kbcs
         self%dke=self%kbce
         self%dls=self%lbcs
         self%dle=self%lbce
      case default
         call stop_execution('set_iface_bc_type', &
               "Unknown BC type: "//self%bc_name)
      end select
   end subroutine set_iface_bc_type

   subroutine set_iface_params(self,nmesh)
      !! Once all the meshes are read, go through the interface BC for
      !! each mesh and set up inter-mesh BC donor IDs as well as fix
      !! negative indices in self or donor BC limits.
      type(mesh_info), dimension(:), intent(inout) :: self
      integer, intent(in) :: nmesh

      integer :: i, n, did
      type(interface_bc), pointer :: ifbc

      do i=1,nmesh
         do n=1,self(i)%nbc_int
            ifbc=>self(i)%ifaceinfo(n)
            ifbc%donorid=get_global_mesh_id(self, &
                  ifbc%donor_mesh)
            did=ifbc%donorid

            ! Fix negative indices on self
            if (ifbc%jbcs < 0) ifbc%jbcs=self(i)%jm+ifbc%jbcs+1
            if (ifbc%jbce < 0) ifbc%jbce=self(i)%jm+ifbc%jbce+1
            if (ifbc%kbcs < 0) ifbc%kbcs=self(i)%km+ifbc%kbcs+1
            if (ifbc%kbce < 0) ifbc%kbce=self(i)%km+ifbc%kbce+1
            if (ifbc%lbcs < 0) ifbc%lbcs=self(i)%lm+ifbc%lbcs+1
            if (ifbc%lbce < 0) ifbc%lbce=self(i)%lm+ifbc%lbce+1

            ! Fix negative indices in donor
            if (ifbc%djs < 0) ifbc%djs=self(did)%jm+ifbc%djs+1
            if (ifbc%dje < 0) ifbc%dje=self(did)%jm+ifbc%dje+1
            if (ifbc%dks < 0) ifbc%dks=self(did)%km+ifbc%dks+1
            if (ifbc%dke < 0) ifbc%dke=self(did)%km+ifbc%dke+1
            if (ifbc%dls < 0) ifbc%dls=self(did)%lm+ifbc%dls+1
            if (ifbc%dle < 0) ifbc%dle=self(did)%lm+ifbc%dle+1
            
         end do
      end do

   end subroutine set_iface_params

   subroutine check_simple_bc(self,nmesh)
      !! Perform some sanity checks on the BC supplied. Only
      !! rudimentary checks are implemented here.
      type(mesh_info), dimension(:), intent(in) :: self
      integer, intent(in) :: nmesh
      
      type(simple_bc), pointer :: sbc
      integer :: i, n, adir, idir
      integer, dimension(3) :: dims, start, end

      do i=1,nmesh
         dims=(/ self(i)%jm, self(i)%km, self(i)%lm /)
         do n=1,self(i)%nbc_sim
            sbc=>self(i)%bcinfo(n)
            start=(/ sbc%jbcs, sbc%kbcs, sbc%lbcs /)
            end=(/ sbc%jbce, sbc%kbce, sbc%lbce /)
            adir=abs(sbc%ibdir)
            idir=sbc%ibdir

            !if (start(adir) /= end(adir)) then
            !   write(STDERR,'(/,A)') 'Mesh: &
            !         &'//trim(adjustl(self(i)%mesh_name))//' BC not &
            !         &on a plane'
            !   call stop_execution('check_simple_bc', &
            !         'Incorrect BC information provided for mesh')
            !end if

            select case (idir)
            case (1:3)
               if (start(adir) /= 1) then
                  write(STDERR,'(A)') "BC not on the MIN plane"
               end if
            case (-3:-1)
               if (end(adir) /= dims(adir)) then
                  write(STDERR,'(A)') "BC not on the MAX plane"
               end if
            case default
               call stop_execution('check_simple_bc', &
                     'Incorrect BC information provided for mesh'&
                     //trim(adjustl(self(i)%mesh_name)))
            end select
         end do
      end do
   end subroutine check_simple_bc

   subroutine check_interface_bc(self,nmesh)
      !! Perform simple checks on the Interface BC conditions
      !! given. The checks are quite rudimentary and only handles
      !! interfaces along the same direction. A more generic version
      !! of this needs to be implemented. However, for meshes used
      !! presently this tests will work fine.
      type(mesh_info), dimension(:), intent(in) :: self
      integer, intent(in) :: nmesh
      
      type(interface_bc), pointer :: sbc
      integer :: i, n, adir, idir, did
      integer, dimension(3) :: dims1, start1, end1
      integer, dimension(3) :: dims2, start2, end2
      
      do i=1,nmesh
         do n=1,self(i)%nbc_int
            sbc=>self(i)%ifaceinfo(n)

            ! First check indices on self
            dims1=(/ self(i)%jm, self(i)%km, self(i)%lm /)
            start1=(/ sbc%jbcs, sbc%kbcs, sbc%lbcs /)
            end1=(/ sbc%jbce, sbc%kbce, sbc%lbce /)
            adir=abs(sbc%ibdir)
            idir=sbc%ibdir
            call bc_limits_check(dims1,start1,end1)

            ! Now check the donor face information
            did=sbc%donorid
            dims2=(/ self(did)%jm, self(did)%km, self(did)%lm /)
            start2=(/ sbc%djs, sbc%dks, sbc%dls /)
            end2=(/ sbc%dje, sbc%dke, sbc%dle /)
            adir=abs(sbc%didir)
            idir=sbc%didir
            call bc_limits_check(dims2,start2,end2)

            call iface_match_check
         end do
      end do

   contains
      subroutine iface_match_check
         !! Check to make sure that the interfaces have the right
         !! number of 1-to-1 grid point overlap. This check is not
         !! 100% correct. We will need to modify this later to provide
         !! a more generic 1-to-1 interface.
         integer :: nn
         integer :: nrecv, ndonor
         
         do nn=1,3
            if (start1(nn) > end1(nn)) then
               nrecv=(start1(nn)-end1(nn)+1)
            else
               nrecv=(end1(nn)-start1(nn)+1)
            end if
            if (start2(nn) > end2(nn)) then
               ndonor=(start2(nn)-end2(nn)+1)
            else
               ndonor=(end2(nn)-start2(nn)+1)
            end if
            if (nrecv /= ndonor) then
               call stop_execution('iface_match_check', &
                     "1-to-1 donor/receiver faces do not match")
            end if
         end do
      end subroutine iface_match_check
      
      subroutine bc_limits_check(dims,start,end)
         !! Helper routine to check for BC on a plane on both receiver
         !! and donor planes.
         integer, dimension(3), intent(in) :: dims, start, end
         
         if (abs(start(adir)-end(adir)) > 2) then
            write(STDERR,'(/,A)') 'Mesh: &
                  &'//trim(adjustl(self(i)%mesh_name))//' BC not &
                  &on a plane'
            call stop_execution('check_interface_bc', &
                  'Incorrect BC information provided for mesh')
         end if
         
         select case (idir)
         case (1:3)
            if (start(adir) /= 1) then
               write(STDERR,'(A)') "BC not on the MIN plane"
            end if
         case (-3:-1)
            if (end(adir) /= dims(adir)) then
               write(STDERR,'(A)') "BC not on the MAX plane"
            end if
         case default
            call stop_execution('check_interface_bc', &
                  'Incorrect BC information provided for mesh'&
                  //trim(adjustl(self(i)%mesh_name)))
         end select
      end subroutine bc_limits_check
   end subroutine check_interface_bc

   subroutine copy_mesh_to_block(mesh,mblock)
      type(mesh_info), intent(in) :: mesh
      type(mesh_block), intent(inout) :: mblock

      mblock%mesh_id=mesh%mesh_id
      mblock%dom_type=mesh%dom_type
      mblock%jm=mesh%jm
      mblock%km=mesh%km
      mblock%lm=mesh%lm

      mblock%gjs=1
      mblock%gks=1
      mblock%gls=1
      mblock%gje=mesh%jm
      mblock%gke=mesh%km
      mblock%gle=mesh%lm

      mblock%nbc_sim=0
      mblock%nbc_int=0
   end subroutine copy_mesh_to_block

   subroutine read_overset_info(self,meshes,flnam)
      type(overset_info), dimension(:), allocatable, intent(inout) :: self
      type(mesh_info), dimension(:), intent(in) :: meshes
      character(len=*), optional :: flnam

      integer :: un, nolap, i, nm, j
      character(len=MAX_NAME_LEN), dimension(:), allocatable :: gnames
      
      if (.not. present(flnam)) then
         un=open_file("overset.inp")
      else
         un=open_file(flnam)
      end if

      read(un,*) nolap
      allocate(self(nolap))

      do i=1,nolap
         read(un,*) nm
         self(i)%nmesh=nm
         allocate(self(i)%mesh_ids(nm))
         allocate(gnames(nm))
         read(un,*) gnames
         do j=1,nm
            self(i)%mesh_ids(j)=get_global_mesh_id(meshes,gnames(j))
         end do
         deallocate(gnames)
         self(i)%master_id=meshes(self(i)%mesh_ids(1))%master_id
      end do
      close(un)
   end subroutine read_overset_info
   
   subroutine print_mesh_info(mesh,un)
      type(mesh_info), intent(in) :: mesh
      integer, intent(in), optional :: un
      
      integer :: i, un1

      if (present(un)) then
         un1=un
      else
         un1=STDOUT
      end if
      
      write(un1,*) 
      write(un1,'(A)') mesh%mesh_name
      write(un1,'(3I5)') mesh%mesh_id, mesh%dom_type, mesh%master_id
      write(un1,'(3I5)') mesh%jm, mesh%km, mesh%lm
      ! write(un1,'(2I10)') mesh%ncells, mesh%nfaces
      write(un1,*) mesh%visc_dir
      write(un1,'(A,2I5)') "BC information: ", mesh%nbc_sim,mesh%nbc_int
      do i=1,mesh%nbc_sim
         call print_simple_bc_info(mesh%bcinfo(i),un1)
         write(un1,'(A)') '------------------------------------------------------------------------'
      end do
      do i=1,mesh%nbc_int
         call print_interface_bc_info(mesh%ifaceinfo(i),un1)
         write(un1,'(A)') '------------------------------------------------------------------------'
      end do
      write(un1,'(/,A,I5)') "Subblock Information: nblocks=",&
            mesh%nblocks
      do i=1,mesh%nblocks
         write(un1,'(I3,":",6I8)') i, mesh%start(JDIR,i), mesh%end(JDIR,i), &
               mesh%start(KDIR,i), mesh%end(KDIR,i), &
               mesh%start(LDIR,i), mesh%end(LDIR,i)
      end do
      write(un1,'(A)') '========================================================================'
      
   end subroutine print_mesh_info

   subroutine print_simple_bc_info(sbc, un)
      type(simple_bc), intent(in) :: sbc
      integer, intent(in), optional :: un
      
      integer :: un1

      if (present(un)) then
         un1=un
      else
         un1=STDOUT
      end if

      write(un1,*)
      write(un1,'(A,/,I4,I4)') trim(sbc%bc_name),sbc%ibtyp,sbc%ibdir
      write(un1,'(6I6)') sbc%jbcs, sbc%jbce, &
            sbc%kbcs, sbc%kbce, sbc%lbcs, sbc%lbce
   end subroutine print_simple_bc_info

   subroutine print_interface_bc_info(sbc, un)
      type(interface_bc), intent(in) :: sbc
      integer, intent(in), optional :: un
      
      integer :: un1

      if (present(un)) then
         un1=un
      else
         un1=STDOUT
      end if
      
      write(un1,*)
      write(un1,'(A,/,I4,I4)') trim(sbc%bc_name),sbc%ibtyp,sbc%ibdir
      write(un1,'(6I6)') sbc%jbcs, sbc%jbce, &
            sbc%kbcs, sbc%kbce, sbc%lbcs, sbc%lbce
      write(un1,'(A,I8,I4)') trim(sbc%donor_mesh),sbc%donorid,sbc%didir
      write(un1,'(6I6)') sbc%djs, sbc%dje, &
            sbc%dks, sbc%dke, sbc%dls, sbc%dle
   end subroutine print_interface_bc_info

   subroutine print_mesh_block_info(self,un)
      type(mesh_block), intent(in) :: self
      integer, intent(in), optional :: un

      integer :: i, un1

      if (present(un)) then
         un1=un
      else
         un1=STDOUT
      end if

      write(un1,*)
      write(un1,'(A)') self%block_name
      write(un1,'(3I4)') self%block_id,self%mesh_id, self%dom_type
      write(un1,'(3I8)') self%jm, self%km, self%lm
      write(un1,'(6I8)') self%gjs, self%gje, self%gks, self%gke, &
            self%gls, self%gle
      write(un1,*) self%is_master
      write(un1,'(A,2I5)') "BC information: ", self%nbc_sim,self%nbc_int
      do i=1,self%nbc_sim
         call print_simple_bc_info(self%bcinfo(i),un1)
         write(un1,'(A)') '------------------------------------------------------------------------'
      end do
      do i=1,self%nbc_int
         call print_interface_bc_info(self%ifaceinfo(i),un1)
         write(un1,'(A)') '------------------------------------------------------------------------'
      end do
      write(un1,'(A)') '========================================================================'
    end subroutine print_mesh_block_info
end module mesh_types

