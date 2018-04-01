module domain_connectivity
   !! A generic domain connectivity module.
   !!
   !! This module defines the datatypes that will be used to store
   !! information about connectivity between overset meshes, and
   !! contains the subrotuines that can perform the connectivity.
   !!
   !! 1. There are three main functions to the connectivity module:
   !!    setup, connectivity information update, interpolation of the
   !!    solution variables.
   !! 2. The connectivity happens in 2 major loops: inter-mesh
   !!    communication, and intra-mesh (or inter-block)
   !!    communication.
   !! 3. In the first (outer) loop, only the master processors of each
   !!    base mesh participate. Amongst these, one master processor is
   !!    chosen as the master for the overset group. Only the overset
   !!    master stores data related to all the meshes in the
   !!    group. The other meshes merely send their mesh coordinates to
   !!    this master and wait for it to perform connectivity
   !!    information.
   !! 4. In the second (inner) loop, the master processor for each
   !!    base mesh, relays the overset information to the
   !!    sub-blocks. The sub-blocks then interpolate the variables and
   !!    relay it back to the master processor. The master processor
   !!    then aggregates the information and relays it back to the
   !!    overset master.
   !! 5. The overset master then (in the outer loop) sends the
   !!    interpolated data to the other master processors and the
   !!    information trickles down to the sub-blocks where the fringe
   !!    points are updated.
   !!
   !! Additional comments:
   !!
   !! 0. The data structures will be used up in only a few
   !!    processors. At most, nmeshes processors.
   !!
   !! 1. The master ID of the first mesh in the overset group is
   !!    chosen as the overset master by default. So you can
   !!    manipulate which processor becomes the master by ordering the
   !!    meshes appropriately.
   !! 2. While the datatypes have been setup to handle the possibility
   !!    of one mesh being part of multiple overset groups, the
   !!    functionality has not been tested rigorously, and it is very
   !!    likely that the code will hang with MPI deadlocks. While some
   !!    clever ordering of the overset groups input in the input file
   !!    might remove the issue, the fix is temporary and a more
   !!    permanent solution should be implemented.
   !!
   !! 3. The data structures in this module are nested, and involve a
   !!    lot of indices to track different information. To perserve
   !!    sanity, the following convention is used to specify indices
   !!    in the subroutines while working with the overset group data
   !!    structures.
   !!
   !!   ii - Index for the overset groups that this mesh belongs to.
   !!   jj - Index for the global overset data if this proc is a
   !!        master, else this is 0
   !!   kk - Index of the overset metadata
   !!   mm - Index of the current mesh in the current overset group(ii)
   !!   mid - Index of the mesh in the mesh metadata array
   !!
   !! 3a. ii and kk are different the mesh might not be part of all
   !!     the overset groups. jj is 0 if the processor is not the
   !!     master of the current overset group; its value can vary
   !!     depending on how many overset groups have this processor as
   !!     their overset master. 
   use domain_info
   implicit none

   integer, parameter :: NOGMAX=2
   
   type ovset_data
      !! The per-mesh data for a given overset mesh. 
      
      ! nfringe - No. of fringe cells where information is
      !           interpolated from donor meshes.
      ! ndonor - No. of donor cells from which information is
      !          interpolated to send to the fringe layer mesh.
      ! imesh - Indices of the fringe cells for each mesh. 
      ! idonor - Indices of the donor cells for each mesh. 
      ! frac - fraction info in the three directions. 
      ! imap - Mapping of the fringe point with the receiver point
      real(kind=rdp), dimension(:,:), pointer :: frac
      integer, dimension(:,:), pointer :: imesh, idonor
      integer, dimension(:,:), pointer :: imap
      integer :: nfringe, ndonor
   end type ovset_data

   type ovset_blk_data
      !! Connectivity data per mesh block, this aggregates all the
      !! group information into one per-block connectivity dataset.

      ! nfringe - No. of fringe points in this mesh block.
      ! ndonor - No. of donor cells in this mesh block.
      ! imst, imend - Pointers to start/end locations for each overset group.
      ! dost, doend - Pointers to donor start/end indices per ovset group.
      ! imesh - j,k,l indices where data is to be inserted.
      ! idonor, frac - j,k,l indices of cells where data is interpolated from.
      ! fridx, doidx - Pointers to the global mesh indices list. 
      integer :: nfringe, ndonor
      integer, dimension(:), pointer :: imst, imend, dost, doend
      integer, dimension(:), pointer :: fridx, doidx
      integer, dimension(:), pointer :: imesh, idonor
      real(kind=rdp), dimension(:), pointer :: frac
   end type ovset_blk_data
   
   type ovset_group
      !! The data for an overset group.

      ! Per-mesh data for all the participating meshes. 
      type(ovset_data), dimension(:), pointer :: doninfo

      ! xyz - coordinates of the grids (3,JM1,KM1,LM1,nm)
      ! iblank - iblank info (JM1,KM1,LM1,nm)
      real(kind=rdp), dimension(:,:,:,:,:), pointer :: xyz
      integer, dimension(:,:,:,:), pointer :: iblank

      ! ibc - BC/mesh type for IHC routines
      ! jhimmune - Immunization for IHC routines
      ! jte - JTAIL1 for IHC routines
      integer, dimension(:), pointer :: ibc, jhimmune, jte
      
      ! nm - No. of meshes in this overset group
      ! jmx, kmx, lmx - Size(nm). Dimensions of the participating
      !                 grids. Note: JM1=maxval(jmx), etc. 
      integer, dimension(:), pointer :: jmx, kmx, lmx
   end type ovset_group

   ! ovset_grids - the overset groups for which this proc is a master
   ! o_ngroups - No. of overset groups in which this mesh is participating
   ! o_nrecvs - No. of overset groups for which this proc is a master
   ! this_ovset - The overset groups in which this mesh is participating
   ! o_grpptr - Pointer to the grid data (if 0, this is not a master)
   ! o_mmax - Maximum number of participating meshes amongst overset groups
   type(ovset_group), dimension(:), pointer :: ovset_grids
   integer, dimension(:), allocatable :: this_ovset
   integer, dimension(:), allocatable :: o_grpptr
   integer :: o_ngroups, o_nrecvs, o_mmax

   ! m_conn_info - Local mesh's connectivity data.
   ! b_conn_info - Connectivity data sorted per mesh block
   ! b_nfringe - No. of fringe points in the local mesh block
   ! b_ndonor - no. of donor points in the local mesh block
   ! b_imesh, b_idon, b_frac - Fringe and donor indices, and fraction
   !                           information for the local mesh block.
   type(ovset_data), dimension(:), pointer :: m_conn_info
   type(ovset_blk_data), dimension(:), pointer :: b_conn_info
   integer :: b_nfringe, b_ndonor
   integer, dimension(:,:), pointer :: b_imesh, b_idon
   real(kind=rdp), dimension(:,:), pointer :: b_frac

   ! Additional notes:
   !    Not all the above declared variables exist or make sense in
   !    all processors. Only some of those arrays are allocated and
   !    are active in all processors. Some are active only on the mesh
   !    masters, and some are active only on overset masters.
   ! 1. Data on overset masters only (all ovset masters are mesh masters):
   !    ovset_grids, o_nrecvs
   ! 2. Data on mesh masters only:
   !    m_conn_info, o_ngroups, this_ovset, o_grpptr,
   ! 3. Data on block processors:
   !    b_conn_info, b_nfringe, b_ndonor, b_imesh, b_idon, b_frac
   !
   ! Every overset master is also a mesh master, every mesh master is
   ! also a block processor, so (1) has (2) & (3), and (2) has (3).
   
   private :: init_ovset_data, init_ovset_group, copy_xyz_tmp_arrays, &
         reset_iblanks, copy_ovset_data, init_ovset_blk_data

   private :: assemble_overset_grids, &
         setup_connectivity_jaina, do_connectivity_jaina, &
         relay_conn_info_to_master, relay_iblanks_to_master, &
         update_block_iblanks, determine_block_conn_info, &
         relay_conn_info_to_blocks
   
contains
   subroutine init_ovset_data(self)
      !! Helper functions to allocate all array members in the
      !! ovset_data datatype.
      type(ovset_data), intent(inout) :: self

      self%nfringe=0
      self%ndonor=0
      allocate(self%imesh(3,dsize3),self%idonor(3,dsize3))
      allocate(self%frac(3,dsize3),self%imap(2,dsize3))
   end subroutine init_ovset_data

   subroutine copy_ovset_data(lhs,rhs)
      !! Helper routine to copy internal data from one structure to
      !! another.
      type(ovset_data), intent(inout) :: lhs
      type(ovset_data), intent(in) :: rhs

      integer :: i, j
      
      lhs%nfringe=rhs%nfringe
      lhs%ndonor=rhs%ndonor

      do i=1,lhs%nfringe
         do j=1,3
            lhs%imesh(j,i)=rhs%imesh(j,i)
         end do
         lhs%imap(1,j)=rhs%imap(1,j)
         lhs%imap(2,j)=rhs%imap(2,j)
      end do

      do i=1,lhs%ndonor
         do j=1,3
            lhs%idonor(j,i)=rhs%idonor(j,i)
            lhs%frac(j,i)=rhs%frac(j,i)
         end do
      end do
   end subroutine copy_ovset_data

   subroutine init_ovset_blk_data(self,ngrp)
      !! Internal helper to initialize arrays in the overset block
      !! data structures.
      type(ovset_blk_data), intent(inout) :: self
      integer, intent(in) :: ngrp

      self%nfringe=0
      self%ndonor=0
      allocate(self%imst(ngrp),self%imend(ngrp))
      allocate(self%dost(ngrp),self%doend(ngrp))
      allocate(self%fridx(dsize3*ngrp),self%doidx(dsize3*ngrp))
      allocate(self%imesh(dsize*ngrp),self%idonor(dsize*ngrp))
      allocate(self%frac(dsize*ngrp))
   end subroutine init_ovset_blk_data
      
   subroutine init_ovset_group(self,nov)
      !! Helper function to allocate all array members. Note: This
      !! uses internal data from the domain_info module and cannot be
      !! a standalone constructor. Hence, delcaring it as private.
      type(ovset_group), intent(inout) :: self
      integer, intent(in) :: nov

      integer :: i
      allocate(self%xyz(3,m_jm,m_km,m_lm,nov))
      allocate(self%iblank(m_jm,m_km,m_lm,nov))
      allocate(self%jmx(nov),self%kmx(nov),self%lmx(nov))
      allocate(self%doninfo(nov))
      allocate(self%ibc(nov),self%jhimmune(nov),self%jte(nov))

      do i=1,nov
         call init_ovset_data(self%doninfo(i))
      end do
      self%iblank=1
   end subroutine init_ovset_group
   
   subroutine setup_connectivity_data
      !! The connectivity module initialization routine. Sets up all
      !! the data structures necessary for performing connectivity and
      !! determines the overset groups that a particular mesh belongs
      !! to.
      integer, dimension(:), allocatable :: sends, recvs
      integer :: nsends, nrecvs
      integer :: n, tag

      tag=101
      ! 1. We need to first determine how many of the overset groups
      ! this mesh is involved in. Then for each of the groups it is
      ! involved in, is it the overset master or if it just a slave
      ! that will be sending data. The connectivity will be performed
      ! only on the master processor and the connectivity data will be
      ! relayed to the slave grids. Once inside the slave grid master
      ! processors they'll be further split to be relayed to the mesh
      ! blocks.
      o_ngroups=100
      nsends=0
      nrecvs=0
      allocate(sends(novset),recvs(novset))
      sends=0
      recvs=0
      ! Only masters need to participate in the overset groups
      if (this_isMaster) then 
         call determine_master_slave_relations

         if ((nsends > 1) .or. (nrecvs > 1) .or. &
             ((nsends > 0) .and. (nrecvs > 0))) then
            if (conn_type == CONN_JAINA) then
               call stop_execution('setup_local_overset_data',&
                     "Only one group per mesh is possible with Jaina's &
                     &connectivity routine.")
            else
               write(STDERR,'(A)') "WARNING: One mesh involved in more than&
                     &two overset groups. This can possibly lead to MPI &
                     &deadlocks and has not been debugged fully!"
            end if
         end if

         call setup_local_overset_data

         do n=1,this_mesh%nblocks
            if (this_mesh%block_id(n) /= (PID+1)) then
               call mpi_bsend(o_ngroups,1,INT_TYPE,(this_mesh%block_id(n)-1),&
                     tag,DEFAULT_COMM,ierr)
            end if
         end do
      else
         call mpi_recv(o_ngroups,1,INT_TYPE,(this_mesh%master_id-1),&
               tag,DEFAULT_COMM,stats_mpi,ierr)
         o_nrecvs=0
      end if
      deallocate(sends,recvs)

      allocate(b_imesh(3,dsize*NOGMAX),b_idon(3,dsize*NOGMAX), &
            b_frac(3,dsize*NOGMAX))
      b_nfringe=0
      b_ndonor=0
   contains
      subroutine determine_master_slave_relations
         !! Wrapper to loop through overset metadata and determine the
         !! groups to which the current mesh master belongs to. Then
         !! decide whether it should be the overset master or just a
         !! slave.
         integer :: i, j, mid, mmid

         o_mmax=0
         do i=1,novset
            mmid=ovset(i)%master_id
            if (ovset(i)%nmesh > o_mmax) o_mmax=ovset(i)%nmesh
            do j=1,ovset(i)%nmesh
               mid=ovset(i)%mesh_ids(j)
               ! Am I involved in this overset group?
               if (this_mesh%mesh_id == mid) then 
                  ! Am I the master for this group?
                  if (this_mesh%master_id == mmid) then 
                     nrecvs=nrecvs+1
                     recvs(nrecvs)=i
                  else
                     nsends=nsends+1
                     sends(nsends)=i
                  end if
               end if
            end do
         end do
      end subroutine determine_master_slave_relations

      subroutine setup_local_overset_data
         !! The groups have been identified and an overset master for
         !! each group has been nominated. Now we loop through this
         !! information and allocate data for both the master and its
         !! slaves. The master gets additional data structures
         !! initialized to hold data for all the meshes in the overset
         !! groups.
         integer :: ii, r, s, mid, jj, mm, kk

         o_nrecvs=nrecvs
         o_ngroups=(nrecvs+nsends)

         ! If I am not part of any overset group, there is nothing for
         ! me to do here.
         if (o_ngroups == 0) return

         ! I am part of one or more overset groups. I need to set up
         ! data.
         nullify(m_conn_info,b_conn_info)
         allocate(m_conn_info(o_ngroups),b_conn_info(this_mesh%nblocks))
         do ii=1,o_ngroups
            call init_ovset_data(m_conn_info(ii))
         end do
         do ii=1,this_mesh%nblocks
            call init_ovset_blk_data(b_conn_info(ii),o_ngroups)
         end do
         allocate(this_ovset(o_ngroups),o_grpptr(o_ngroups))
         if (nrecvs > 0) then 
            allocate(ovset_grids(nrecvs))
            do ii=1,nrecvs
               call init_ovset_group(ovset_grids(ii),ovset(recvs(ii))%nmesh)
            end do
         end if
         
         r=1; s=1
         ngroups: do ii=1,o_ngroups
            if ((r <= nrecvs) .and. (s <= nsends)) then
               send_recv: if (sends(s) > recvs(r)) then 
                  this_ovset(ii)=recvs(r)
                  o_grpptr(ii)=r
                  r=r+1
               else if (sends(s) < recvs(r)) then
                  this_ovset(ii)=sends(s)
                  o_grpptr(ii)=0
                  s=s+1
               else
                  call stop_execution('setup_local_overset_data',&
                        "Cannot have same processing doing send and recv.")
               end if send_recv
            else if (s > nsends) then 
               this_ovset(ii)=recvs(r)
               o_grpptr(ii)=r
               r=r+1
            else if (r > nrecvs) then 
               this_ovset(ii)=sends(s)
               o_grpptr(ii)=0
               s=s+1
            else
               call stop_execution('setup_local_overset_data', &
                     "Failed logic error when determining connectivity.")
            end if
         end do ngroups

         ! Set the J,K,L dimensions for each of the meshes in the
         ! overset group...

         !-> There are several indices that come into play with this
         !-> data structure. To preserve sanity, we will follow the
         !-> same indices in all the routines inside this module with
         !-> the following meanings:
         ! ii - Counter for number of groups this mesh belongs to.
         ! jj - Index for the master data if this proc is a master.
         ! kk - Index of the overset metadata array.
         ! mm - Index of the current mesh in the overset group. 
         ! mid - Index of the mesh in the mesh metadata array.
         do ii=1,o_ngroups
            jj=o_grpptr(ii)
            ! If 0 then this mesh is a slave, so nothing to do. 
            if (jj == 0) cycle

            kk=this_ovset(ii)
            do mm=1,ovset(kk)%nmesh
               mid=ovset(kk)%mesh_ids(mm)
               ovset_grids(jj)%jmx(mm)=meshes(mid)%jm
               ovset_grids(jj)%kmx(mm)=meshes(mid)%km
               ovset_grids(jj)%lmx(mm)=meshes(mid)%lm
            end do
         end do
      end subroutine setup_local_overset_data
   end subroutine setup_connectivity_data

   subroutine assemble_overset_grids
      !! For each overset group, assemble the member grids in the
      !! overset master processor. All the data necessary is already
      !! available in the domain_info and domain_connectivity modules,
      !! so nothing needs to be sent in as arguments.
      
      !-> There are several indices that come into play with this
      !-> data structure. To preserve sanity, we will follow the
      !-> same indices in all the routines inside this module with
      !-> the following meanings:
      ! ii - Counter for number of groups this mesh belongs to.
      ! jj - Index for the master data if this proc is a master.
      ! kk - Index of the overset metadata array.
      ! mm - Index of the current mesh in the overset group. 
      ! mid - Proc ID of the mesh master
      ! o_master - Proc ID of the overset master
      integer :: ii, jj, kk, mm, mid, o_master

      integer :: bufsize, i, j, k, l, tag
      real(kind=rdp), dimension(:), allocatable :: xbuffer

      tag=201
      bufsize=m_jm*m_km*m_lm*3
      allocate(xbuffer(bufsize))

      call barrier_mpi
      if (.not. this_isMaster) go to 100
      !call print_message("Assembling overset meshes")
      ovset_groups: do ii=1,o_ngroups
         jj=o_grpptr(ii)
         kk=this_ovset(ii)

         master_proc: if (jj == 0) then
            ! I am a slave, I should send data
            o_master=ovset(kk)%master_id-1
            mid=this_mesh%master_id-1
            call copy_mesh_to_buffer
            bufsize=this_mesh%jm*this_mesh%km*this_mesh%lm*3
            call mpi_bsend(xbuffer,bufsize,REAL_TYPE,o_master,tag,&
                  DEFAULT_COMM,ierr)
         else master_proc
            group_meshes: do mm=1,ovset(kk)%nmesh
               mid=ovset(kk)%mesh_ids(mm)
               o_master=meshes(mid)%master_id-1
               remote_mesh: if (mid == this_mesh%mesh_id) then
                  call copy_local_mesh
               else remote_mesh
                  bufsize=meshes(mid)%jm*meshes(mid)%km*meshes(mid)%lm*3
                  call mpi_recv(xbuffer,bufsize,REAL_TYPE,o_master,tag, &
                        DEFAULT_COMM,stats_mpi,ierr)
                  call copy_buffer_to_overset_group
               end if remote_mesh
            end do group_meshes
         end if master_proc
      end do ovset_groups

100   call barrier_mpi
      deallocate(xbuffer)
   contains
      !! Wrapper routines so that we can see the if-then-else and
      !! do-loop logic in the main routine.
      subroutine copy_mesh_to_buffer
         !! Copy the local mesh into the MPI send buffer
         i=0
         do l=1,this_mesh%lm
            do k=1,this_mesh%km
               do j=1,this_mesh%jm
                  xbuffer(i+1)=m_grid%x(j,k,l)
                  xbuffer(i+2)=m_grid%y(j,k,l)
                  xbuffer(i+3)=m_grid%z(j,k,l)
                  i=i+3
               end do
            end do
         end do
      end subroutine copy_mesh_to_buffer

      subroutine copy_local_mesh
         !! Copy the local mesh into the overset group data set
         do l=1,this_mesh%lm
            do k=1,this_mesh%km
               do j=1,this_mesh%jm
                  ovset_grids(jj)%xyz(1,j,k,l,mm)=m_grid%x(j,k,l)
                  ovset_grids(jj)%xyz(2,j,k,l,mm)=m_grid%y(j,k,l)
                  ovset_grids(jj)%xyz(3,j,k,l,mm)=m_grid%z(j,k,l)
               end do
            end do
         end do
      end subroutine copy_local_mesh

      subroutine copy_buffer_to_overset_group
         !! Copy grid from remote processor (MPI recv buffer) into the
         !! appropriate overset data set
         i=0
         ldims: do l=1,ovset_grids(jj)%lmx(mm)
            kdims: do k=1,ovset_grids(jj)%kmx(mm)
               jdims: do j=1,ovset_grids(jj)%jmx(mm)
                  ovset_grids(jj)%xyz(1,j,k,l,mm)=xbuffer(i+1)
                  ovset_grids(jj)%xyz(2,j,k,l,mm)=xbuffer(i+2)
                  ovset_grids(jj)%xyz(3,j,k,l,mm)=xbuffer(i+3)
                  i=i+3
               end do jdims
            end do kdims
         end do ldims
      end subroutine copy_buffer_to_overset_group
   end subroutine assemble_overset_grids

   subroutine setup_connectivity_jaina
      !! Jaina's connectivity routines can only handle certain
      !! cases. So impose several checks to ensure that connectivity
      !! makes sense, and set up the octree data for the background
      !! mesh.
      integer :: ii, jj, kk, mm

      type(ovset_group), pointer :: cur_group
      real, dimension(:,:,:), allocatable :: xtmp, ytmp, ztmp

      ! Am I a mesh master? If not, I've nothing to do here
      if ((.not. this_isMaster) .or. (o_ngroups == 0)) return

      ! Perform sanity checks to make sure we can do connectivity.
      if (o_ngroups > 1) call stop_execution('setup_connectivity_jaina',&
            "One mesh can belong to only one group for overset connectivity")

      ii=1
      jj=o_grpptr(ii)
      kk=this_ovset(ii)
      mm=2
      ! Am I the overset master? If not, I don't have to bother with
      ! connectivity routines. I'll just wait for the data from the
      ! master.
      if (ovset(kk)%master_id /= this_mesh%master_id) return

      if (ovset(kk)%nmesh /= 2) then
         call stop_execution('setup_jaina_connectivity', &
               "Only two mesh groups can be used with this connectivity")
      end if
      ! Math: This IHC routine wasn't modified to deal with DOM_IBC or DOM_FUS
      if ((meshes(ovset(kk)%mesh_ids(1))%dom_type == DOM_IBC) .or. &
          (meshes(ovset(kk)%mesh_ids(2))%dom_type == DOM_FUS)) then
         call stop_execution('setup_connectivity_jaina', &
         "Jaina's IHC routine wasn't modified to deal &
          &with DOM_IBC or DOM_FUS")
      end if
      if ((meshes(ovset(kk)%mesh_ids(1))%dom_type /= DOM_BLADE) .or. &
          (meshes(ovset(kk)%mesh_ids(2))%dom_type /= DOM_WAKE) .or. &
          (meshes(ovset(kk)%mesh_ids(2))%dom_type /= DOM_GROUNDWAKE)) then
         call stop_execution('setup_connectivity_jaina', &
               "Overset master has to be a blade mesh")
      end if

      !call print_message("Building octree for Jaina's connectivity")
      cur_group=>ovset_grids(jj)
      call copy_xyz_tmp_arrays(cur_group,mm,xtmp,ytmp,ztmp)
      call make_octree(xtmp,ytmp,ztmp,cur_group%jmx(mm), &
            cur_group%kmx(mm),cur_group%lmx(mm))

      deallocate(xtmp,ytmp,ztmp)
   end subroutine setup_connectivity_jaina

   subroutine copy_xyz_tmp_arrays(cur_group,mm,xtmp,ytmp,ztmp)
      !! Wrapper routine to reshape arrays into format that Jaina's
      !! connectivity routine takes.
      type(ovset_group), pointer, intent(in) :: cur_group
      integer, intent(in) :: mm
      real, dimension(:,:,:), allocatable, intent(inout) :: xtmp
      real, dimension(:,:,:), allocatable, intent(inout) :: ytmp
      real, dimension(:,:,:), allocatable, intent(inout) :: ztmp

      integer :: j, k, l
      
      allocate(xtmp(cur_group%jmx(mm),cur_group%kmx(mm),cur_group%lmx(mm)))
      allocate(ytmp(cur_group%jmx(mm),cur_group%kmx(mm),cur_group%lmx(mm)))
      allocate(ztmp(cur_group%jmx(mm),cur_group%kmx(mm),cur_group%lmx(mm)))
      do l=1,cur_group%lmx(mm)
         do k=1,cur_group%kmx(mm)
            do j=1,cur_group%jmx(mm)
               xtmp(j,k,l)=cur_group%xyz(1,j,k,l,mm)
               ytmp(j,k,l)=cur_group%xyz(2,j,k,l,mm)
               ztmp(j,k,l)=cur_group%xyz(3,j,k,l,mm)
            end do
         end do
      end do
   end subroutine copy_xyz_tmp_arrays

   subroutine perform_domain_connectivity(x,y,z,iblank,is_init)
      !! Assemble the blocks in the mesh master, relay the mesh
      !! information to the overset master, pre-process connectivity
      !! information, perform connectivity and determine
      !! donors/receivers, perform post-processing and relay data all
      !! the way down to the mesh blocks.
      !!
      !! Inputs:
      !!    x,y,z - Coordinates of the local mesh block
      !!    iblank - I-blanking for the local mesh block
      !!    is_init - Are we in the initialization phase? If true then
      !!              the data setup routines are called, if false
      !!              then it is assumed that the arrays have been all
      !!              allocated properly.

      real, dimension(jmax,kmax,lmax), intent(in) :: x,y,z
      integer, dimension(jmax,kmax,lmax), intent(inout) :: iblank
      logical, intent(in) :: is_init

      call assemble_mesh(x,y,z)
      call assemble_overset_grids
      call reset_iblanks

      select case(conn_type)
      case (CONN_JAINA)
         if (is_init) call setup_connectivity_jaina
         call do_connectivity_jaina
      case (CONN_IHC)
         if (is_init) call setup_connectivity_ihc

         ! Math: send IBC iblank info from IBC blocks to IBC mesh master
         !       and then to overset master for IHC (create iblk_ibc)
         ! If no IBC mesh, it will just return
         ! Update: not necessary because IBC init was done on IBC master 
         !         and iblank info was sent to overset master then
         !if (is_init) call assemble_iblk_ibc(iblank)

         call do_connectivity_ihc
      case default
         call stop_execution('perform_domain_connectivity', &
               "Unknown connectivity type specified.")
      end select

      call relay_iblanks_to_master
      call update_block_iblanks(iblank)
      call relay_conn_info_to_master
      call relay_conn_info_to_blocks
      
   end subroutine perform_domain_connectivity

   subroutine do_connectivity_jaina
      !! Performs Jaina's connectivity algorithm. We make several
      !! assumptions here: 1. each mesh belongs to only one overset
      !! group, 2. each overset group contains exactly two meshes (a
      !! blade and a background mesh).
      use rootVariables
      
      integer :: ii, jj, kk
      type(ovset_group), pointer :: og
      type(ovset_data), pointer :: bld,wke
      integer :: i, j, k, jjm, kkm, llm
      real, dimension(:,:,:), allocatable :: x1b, y1b, z1b
      real, dimension(:,:,:), allocatable :: x2w, y2w, z2w

      nfringe=0
      ndonor=0

      if ((.not. this_isMaster) .or. (o_ngroups == 0)) go to 100
      ii=1
      jj=o_grpptr(ii)
      kk=this_ovset(ii)

      if (ovset(kk)%master_id /= this_mesh%master_id) go to 100
      og=>ovset_grids(jj)
      call copy_xyz_tmp_arrays(og,1,x1b,y1b,z1b)
      call copy_xyz_tmp_arrays(og,2,x2w,y2w,z2w)
      jjm=og%jmx(2)
      kkm=og%kmx(2)
      llm=og%lmx(2)

      ! In this routine nfringe and ndonor take on different meanings
      ! depending on whether we are talking about the blade or the
      ! background mesh. Another quirk is that the routine uses
      ! imesh(:,3) as the array to store indices. This needs to be
      ! fixed too.
      call do_connect(x1b,y1b,z1b,og%jmx(1),og%kmx(1),og%lmx(1), &
            x2w,y2w,z2w, &
            og%iblank(1:jjm,1:kkm,1:llm,2),og%jmx(2),og%kmx(2),og%lmx(2),&
            imesh,idonor,frac,imesh1,idonor1,frac1,ndonor,nfringe,dsize3, &
            psi,m_root_co)      

      ! Duality of nfringe/ndonor is resolved here. 
      bld=>og%doninfo(1)
      wke=>og%doninfo(2)
      bld%nfringe=ndonor
      bld%ndonor=nfringe
      wke%nfringe=nfringe
      wke%ndonor=ndonor

      ! Copy wake->blade interpolation information
      do i=1,ndonor
         do j=1,3
            k=i+(j-1)*dsize3
            bld%imesh(j,i)=imesh(k)
            wke%idonor(j,i)=idonor(k)
            wke%frac(j,i)=frac(k)
         end do
         bld%imap(1,i)=i
         bld%imap(2,i)=2
      end do

      ! Copy blade->wake interpolation information
      do i=1,nfringe
         do j=1,3
            k=i+(j-1)*dsize3
            wke%imesh(j,i)=imesh1(k)
            bld%idonor(j,i)=idonor1(k)
            bld%frac(j,i)=frac1(k)
         end do
         wke%imap(1,i)=i
         wke%imap(2,i)=1
      end do

      deallocate(x1b,y1b,z1b,x2w,y2w,z2w)
100   call barrier_mpi
   end subroutine do_connectivity_jaina

   subroutine reset_iblanks
      !! Simple wrapper to reset iblank information for the global
      !! connectivity data in the overset master processor.
      integer :: ii, jj, kk, mm
      integer :: j, k, l
      type(ovset_group), pointer :: og

      if (.not. this_isMaster) return
      ovset_groups: do ii=1,o_ngroups
         jj=o_grpptr(ii)
         kk=this_ovset(ii)

         if (jj == 0) cycle

         og=>ovset_grids(jj)
         og%iblank=1
         !-??? Should we just do og%iblank=1? Explicitly stating the
         !     loop bounds will be useful for meshes with skewed ratios of
         !     cells :-/
         ! grp_meshes: do mm=1,ovset(kk)%nmesh
         !    do l=1,og%lmx(mm)
         !       do k=1,og%kmx(mm)
         !          do j=1,og%jmx(mm)
         !             og%iblank(j,k,l,mm)=1
         !          end do
         !       end do
         !    end do
         ! end do grp_meshes
      end do ovset_groups
   end subroutine reset_iblanks

   subroutine relay_conn_info_to_master
      !! Relay per-mesh connectivity information from the overset
      !! master to the mesh masters.
      integer :: ii, jj, kk, mm, mid, o_master

      integer :: bufsize, i, j
      integer, dimension(:), allocatable :: imbuf, donbuf
      real(kind=rdp), dimension(:), allocatable :: fbuf
      type(ovset_data), pointer :: condat, ogdon

      bufsize=dsize
      allocate(imbuf(dsize),donbuf(dsize),fbuf(dsize))

      call barrier_mpi
      if (.not. this_isMaster) go to 100
      !call print_message("Relay connectivity info to mesh masters")
      ovset_groups: do ii=1,o_ngroups
         jj=o_grpptr(ii)
         kk=this_ovset(ii)
         condat=>m_conn_info(ii)

         master_proc: if (jj == 0) then
            o_master=ovset(kk)%master_id-1
            mid=this_mesh%master_id-1
            call mpi_recv(condat%nfringe,1,INT_TYPE,o_master,501, &
                  DEFAULT_COMM,stats_mpi,ierr)
            call mpi_recv(condat%ndonor,1,INT_TYPE,o_master,502, &
                  DEFAULT_COMM,stats_mpi,ierr)
            call mpi_recv(imbuf,condat%nfringe*3,INT_TYPE,o_master,503, &
                  DEFAULT_COMM,stats_mpi,ierr)
            call mpi_recv(donbuf,condat%ndonor*3,INT_TYPE,o_master,504, &
                  DEFAULT_COMM,stats_mpi,ierr)
            call mpi_recv(fbuf,condat%ndonor*3,REAL_TYPE,o_master,505, &
                  DEFAULT_COMM,stats_mpi,ierr)
            call copy_buffer_to_mdata
         else master_proc
            group_meshes: do mm=1,ovset(kk)%nmesh
               mid=ovset(kk)%mesh_ids(mm)
               o_master=meshes(mid)%master_id-1
               ogdon=>ovset_grids(jj)%doninfo(mm)
               remote_mesh: if (mid == this_mesh%mesh_id) then
                  call copy_ovset_data(condat,ogdon)
               else
                  call copy_data_to_buffer
                  call mpi_bsend(ogdon%nfringe,1,INT_TYPE,o_master,501, &
                        DEFAULT_COMM,ierr)
                  call mpi_bsend(ogdon%ndonor,1,INT_TYPE,o_master,502, &
                        DEFAULT_COMM,ierr)
                  call mpi_bsend(imbuf,ogdon%nfringe*3,INT_TYPE,o_master,503, &
                        DEFAULT_COMM,ierr)
                  call mpi_bsend(donbuf,ogdon%ndonor*3,INT_TYPE,o_master,504, &
                        DEFAULT_COMM,ierr)
                  call mpi_bsend(fbuf,ogdon%ndonor*3,REAL_TYPE,o_master,505, &
                        DEFAULT_COMM,ierr)
               end if remote_mesh
            end do group_meshes
         end if master_proc
      end do ovset_groups
      
      deallocate(imbuf,donbuf,fbuf)
100   call barrier_mpi
   contains
      subroutine copy_data_to_buffer
         ! Copy the connectivity data to temporary buffers used for
         ! MPI send/recv operations.
         i=0
         do j=1,ogdon%nfringe
            imbuf(i+1)=ogdon%imesh(1,j)
            imbuf(i+2)=ogdon%imesh(2,j)
            imbuf(i+3)=ogdon%imesh(3,j)
            i=i+3
         end do
         i=0
         do j=1,ogdon%ndonor
            donbuf(i+1)=ogdon%idonor(1,j)
            donbuf(i+2)=ogdon%idonor(2,j)
            donbuf(i+3)=ogdon%idonor(3,j)
            fbuf(i+1)=ogdon%frac(1,j)
            fbuf(i+2)=ogdon%frac(2,j)
            fbuf(i+3)=ogdon%frac(3,j)
            i=i+3
         end do
      end subroutine copy_data_to_buffer

      subroutine copy_buffer_to_mdata
         ! Helper routine to copy data from temporary MPI recv buffers
         ! to the master connectivity data structures.
         i=0
         do j=1,condat%nfringe
            condat%imesh(1,j)=imbuf(i+1)
            condat%imesh(2,j)=imbuf(i+2)
            condat%imesh(3,j)=imbuf(i+3)
            i=i+3
         end do
         i=0
         do j=1,condat%ndonor
            condat%idonor(1,j)=donbuf(i+1)
            condat%idonor(2,j)=donbuf(i+2)
            condat%idonor(3,j)=donbuf(i+3)
            condat%frac(1,j)=fbuf(i+1)
            condat%frac(2,j)=fbuf(i+2)
            condat%frac(3,j)=fbuf(i+3)
            i=i+3
         end do
      end subroutine copy_buffer_to_mdata
   end subroutine relay_conn_info_to_master

   subroutine relay_iblanks_to_master
      !! Relay per-mesh iblank information to the mesh masters. Note
      !! that this routine needs to account for the possibility that
      !! one mesh might be part of more than one overset
      !! group. Therefore, we cannot just assign the connectivity
      !! iblank information without checks. While the grid point might
      !! be a field point w.r.t. one group, it might be a hole/fringe
      !! point w.r.t. another overset group. Currently this routine
      !! only handles hole/field point distinction. We need to
      !! implement the hole-fringe-field point distinctions here.
      integer :: ii, jj, kk, mm, mid, o_master

      integer :: bufsize, i, j, k, l, tag
      integer, dimension(:), allocatable :: iblbuf

      tag=201
      bufsize=m_jm*m_km*m_lm
      allocate(iblbuf(bufsize))

      call barrier_mpi
      if (.not. this_isMaster) go to 100

      do l=1,this_mesh%lm
         do k=1,this_mesh%km
            do j=1,this_mesh%jm
               m_iblank(j,k,l)=1
            end do
         end do
      end do
      
      !call print_message("Relay iblanks to masters")
      ovset_groups: do ii=1,o_ngroups
         jj=o_grpptr(ii)
         kk=this_ovset(ii)

         master_proc: if (jj == 0) then
            o_master=ovset(kk)%master_id-1
            mid=this_mesh%master_id-1
            call mpi_recv(iblbuf,bufsize,INT_TYPE,o_master,tag, &
                  DEFAULT_COMM,stats_mpi,ierr)
            
            i=0
            do l=1,this_mesh%lm
               do k=1,this_mesh%km
                  do j=1,this_mesh%jm
                     i=i+1
                     m_iblank(j,k,l)=min(iblbuf(i),m_iblank(j,k,l))
                  end do
               end do
            end do
         else master_proc
            group_meshes: do mm=1,ovset(kk)%nmesh
               mid=ovset(kk)%mesh_ids(mm)
               o_master=meshes(mid)%master_id-1
               remote_mesh: if (mid == this_mesh%mesh_id) then
                  do l=1,this_mesh%lm
                     do k=1,this_mesh%km
                        do j=1,this_mesh%jm
                           m_iblank(j,k,l)=min(ovset_grids(jj)%iblank(j,k,l,mm)&
                                 ,m_iblank(j,k,l))
                        end do
                     end do
                  end do
               else remote_mesh
                  i=0
                  ldims: do l=1,ovset_grids(jj)%lmx(mm)
                     kdims: do k=1,ovset_grids(jj)%kmx(mm)
                        jdims: do j=1,ovset_grids(jj)%jmx(mm)
                           i=i+1
                           iblbuf(i)=ovset_grids(jj)%iblank(j,k,l,mm)
                        end do jdims
                     end do kdims
                  end do ldims
                  call mpi_bsend(iblbuf,bufsize,INT_TYPE,o_master,tag, &
                        DEFAULT_COMM,ierr)
               end if remote_mesh
            end do group_meshes
         end if master_proc
      end do ovset_groups

      deallocate(iblbuf)
100   call barrier_mpi
   end subroutine relay_iblanks_to_master

   subroutine update_block_iblanks(iblank)
      !! Relay per-block iblank information to the individual blocks
      !! from the mesh masters. Here we also send the iblank
      !! information for the ghost points as well.
      integer, dimension(jmax,kmax,lmax), intent(inout) :: iblank

      integer :: j, k, l, bufsize, i, tag, n
      integer :: j1, k1, l1
      integer, dimension(:), allocatable :: iblbuf
      integer, dimension(:), pointer :: bids
      integer, dimension(:,:), pointer :: bst, bend
      integer, dimension(6) :: blk_lim

      tag=101
      bufsize=b_jm*b_km*b_lm
      allocate(iblbuf(bufsize))

      !call print_message("Updating block iblanks")
      call barrier_mpi
      master_check: if (this_isMaster) then
         bids=>this_mesh%block_id
         bst=>this_mesh%start
         bend=>this_mesh%end
         mesh_blocks: do n=1,this_mesh%nblocks
            if (bids(n) == (PID+1)) then
               do l=jkl_lims(LDIR),jkl_lims(LDIR+3)
                  l1=l-jkl_lims(LDIR)+1
                  do k=jkl_lims(KDIR),jkl_lims(KDIR+3)
                     k1=k-jkl_lims(KDIR)+1
                     do j=jkl_lims(JDIR),jkl_lims(JDIR+3)
                        j1=j-jkl_lims(JDIR)+1
                        iblank(j1,k1,l1)=m_iblank(j,k,l)
                     end do
                  end do
               end do
            else
               call determine_send_limits
               i=0
               do l=blk_lim(LDIR),blk_lim(LDIR+3)
                  do k=blk_lim(KDIR),blk_lim(KDIR+3)
                     do j=blk_lim(JDIR),blk_lim(JDIR+3)
                        i=i+1
                        iblbuf(i)=m_iblank(j,k,l)
                     end do
                  end do
               end do
               call mpi_bsend(iblbuf,bufsize,INT_TYPE,bids(n)-1,tag, &
                     DEFAULT_COMM,ierr)
            end if
         end do mesh_blocks
      else master_check
         call mpi_recv(iblbuf,bufsize,INT_TYPE,(this_mesh%master_id-1),tag,&
               DEFAULT_COMM,stats_mpi,ierr)
         i=0
         do l=1,lmax
            do k=1,kmax
               do j=1,jmax
                  i=i+1
                  iblank(j,k,l)=iblbuf(i)
               end do
            end do
         end do
      end if master_check
      call barrier_mpi
      deallocate(iblbuf)
   contains
      subroutine determine_send_limits
         !! Determine the lower and upper j,k,l limits of each block
         !! taking ghost points into account. This information isn't
         !! stored in the master processor, so we need to calculate
         !! them.
         integer, dimension(6) :: extents
         integer, dimension(3) :: maxes

         extents=(/ bst(JDIR,n), bst(KDIR,n), bst(LDIR,n), &
               bend(JDIR,n), bend(KDIR,n), bend(LDIR,n) /)
         maxes=(/ this_mesh%jm, this_mesh%km, this_mesh%lm /)

         do i=1,3
            if (extents(i) == 1) then 
               blk_lim(i) = extents(i)
            else
               blk_lim(i) = extents(i)-NGHOST
            end if
         end do
         
         ! The max faces
         do i=4,6
            if (extents(i) == maxes(i-3)) then 
               blk_lim(i) = extents(i)
            else
               blk_lim(i) = extents(i)+NGHOST
            end if
         end do
      end subroutine determine_send_limits
   end subroutine update_block_iblanks

   subroutine determine_block_conn_info
      !! For each overset group's mesh parse through the
      !! receiver/donor point information and sort them out into bins
      !! for each mesh block. At the end of the routine, the block
      !! connectivity data will contain the receiver/donor point
      !! information for that mesh block from all the overset groups
      !! that the block's mesh is a member of.
      integer :: ii, n, i, iidx, j
      type(ovset_data), pointer :: condat
      type(ovset_blk_data), pointer :: blk

      ! Limits for donors
      integer, dimension(:,:), pointer :: bst, bend
      ! Limits for fringe points
      integer, dimension(:,:), allocatable :: bstf, bendf
      integer, dimension(3) :: maxes
      integer :: idir
      
      maxes=(/ this_mesh%jm, this_mesh%km, this_mesh%lm /)
      allocate(bstf(3,this_mesh%nblocks),bendf(3,this_mesh%nblocks))
      nullify(bst,bend)
      bst=>this_mesh%start
      bend=>this_mesh%end
      bstf=bst
      bendf=bend
      ! Set up preliminary data 
      do n=1,this_mesh%nblocks
         !if (n > 1) bstf(KDIR,n)=bstf(KDIR,n)-NGHOST
         !if (n < this_mesh%nblocks) bendf(KDIR,n)=bendf(KDIR,n)+NGHOST

         ! Update fringe points to include ghost point info
         do idir=1,3
            if (bstf(idir,n) > 1) bstf(idir,n)=bstf(idir,n)-NGHOST
            if (bendf(idir,n) < maxes(idir)) &
                  bendf(idir,n)=bendf(idir,n)+NGHOST
         end do
         ! Do some housekeeping, reset all indices. 
         b_conn_info(n)%nfringe=0
         b_conn_info(n)%ndonor=0
         do ii=1,o_ngroups
            b_conn_info(n)%imst(ii)=0
            b_conn_info(n)%imend(ii)=-1
            b_conn_info(n)%dost(ii)=0
            b_conn_info(n)%doend(ii)=-1
         end do
      end do

      ! Go through each overset group and determine the recv/send
      ! points that belong to the individual sub-blocks.
      do ii=1,o_ngroups
         condat=>m_conn_info(ii)

         ! Pass through the fringe points
         do i=1,condat%nfringe
            do n=1,this_mesh%nblocks
               blk=>b_conn_info(n)
               
               if ((bstf(JDIR,n) <= condat%imesh(JDIR,i)) .and. &
                   (condat%imesh(JDIR,i) <= bendf(JDIR,n)) .and. &
                   (bstf(KDIR,n) <= condat%imesh(KDIR,i)) .and. &
                   (condat%imesh(KDIR,i) <= bendf(KDIR,n)) .and. &
                   (bstf(LDIR,n) <= condat%imesh(LDIR,i)) .and. &
                   (condat%imesh(LDIR,i) <= bendf(LDIR,n))) then
                  iidx=blk%nfringe*3
                  blk%nfringe=blk%nfringe+1
                  do j=1,3
                     blk%imesh(iidx+j)=condat%imesh(j,i)
                  end do
                  blk%fridx(blk%nfringe)=i
                  !exit
               end if
            end do
         end do
         ! Pass through donor points
         do i=1,condat%ndonor
            do n=1,this_mesh%nblocks
               blk=>b_conn_info(n)
               
               if ((bst(JDIR,n) <= condat%idonor(JDIR,i)) .and. &
                   (condat%idonor(JDIR,i) <= bend(JDIR,n)) .and. &
                   (bst(KDIR,n) <= condat%idonor(KDIR,i)) .and. &
                   (condat%idonor(KDIR,i) <= bend(KDIR,n)) .and. &
                   (bst(LDIR,n) <= condat%idonor(LDIR,i)) .and. &
                   (condat%idonor(LDIR,i) <= bend(LDIR,n))) then
                  iidx=blk%ndonor*3
                  blk%ndonor=blk%ndonor+1
                  do j=1,3
                     blk%idonor(iidx+j)=condat%idonor(j,i)
                     blk%frac(iidx+j)=condat%frac(j,i)
                  end do
                  blk%doidx(blk%ndonor)=i
                  exit
               end if
            end do
         end do

         ! Now pass through the mesh blocks and determine if a
         ! particular mesh block contributes/needs data from a
         ! particular overset group.
         do n=1,this_mesh%nblocks
            blk=>b_conn_info(n)

            if (ii == 1) then
               if (blk%nfringe > 0) then
                  blk%imst(ii)=1
                  blk%imend(ii)=blk%nfringe
               end if
               if (blk%ndonor > 0) then
                  blk%dost(ii)=1
                  blk%doend(ii)=blk%ndonor
               end if
            else
               if (blk%nfringe > blk%imend(ii-1)) then
                  blk%imst(ii)=blk%imend(ii-1)+1
                  blk%imend(ii)=blk%nfringe
               end if
               if (blk%ndonor > blk%doend(ii-1)) then
                  blk%dost=blk%doend(ii-1)+1
                  blk%doend=blk%ndonor
               end if
            end if
         end do
      end do

      deallocate(bstf,bendf)
   end subroutine determine_block_conn_info

   subroutine relay_conn_info_to_blocks
      !! Relay the connectivity information to individual blocks from
      !! mesh masters.
      integer :: n, i, j, i1
      integer, dimension(:), allocatable :: imbuf, idbuf
      real(kind=rdp), dimension(:), allocatable :: fbuf
      type(ovset_blk_data), pointer :: blk

      allocate(imbuf(3*dsize*NOGMAX),idbuf(3*dsize*NOGMAX),&
            fbuf(3*dsize*NOGMAX))
      call barrier_mpi
      if (o_ngroups == 0) goto 100
      !call print_message("Relay connectivity info to blocks")
      master_check: if (this_isMaster) then
         call determine_block_conn_info
         mesh_blocks: do n=1,this_mesh%nblocks
            blk=>b_conn_info(n)
            if (this_mesh%block_id(n) == (PID+1)) then
               b_nfringe=blk%nfringe
               b_ndonor=blk%ndonor
               call insert_conn_info(blk%nfringe,blk%ndonor, &
                     blk%imesh,blk%idonor,blk%frac)
            else
               call mpi_bsend(blk%nfringe,1,INT_TYPE, &
                     this_mesh%block_id(n)-1,101,DEFAULT_COMM,ierr)
               call mpi_bsend(blk%ndonor,1,INT_TYPE, &
                     this_mesh%block_id(n)-1,102,DEFAULT_COMM,ierr)
               call mpi_bsend(blk%imesh,blk%nfringe*3,INT_TYPE, &
                     this_mesh%block_id(n)-1,103,DEFAULT_COMM,ierr)
               call mpi_bsend(blk%idonor,blk%ndonor*3,INT_TYPE, &
                     this_mesh%block_id(n)-1,104,DEFAULT_COMM,ierr)
               call mpi_bsend(blk%frac,blk%ndonor*3,REAL_TYPE, &
                     this_mesh%block_id(n)-1,105,DEFAULT_COMM,ierr)
            end if
         end do mesh_blocks
      else
         call mpi_recv(b_nfringe,1,INT_TYPE,(this_mesh%master_id-1), &
               101,DEFAULT_COMM,stats_mpi,ierr)
         call mpi_recv(b_ndonor,1,INT_TYPE,(this_mesh%master_id-1), &
               102,DEFAULT_COMM,stats_mpi,ierr)
         call mpi_recv(imbuf,b_nfringe*3,INT_TYPE,(this_mesh%master_id-1), &
               103,DEFAULT_COMM,stats_mpi,ierr)
         call mpi_recv(idbuf,b_ndonor*3,INT_TYPE,(this_mesh%master_id-1), &
               104,DEFAULT_COMM,stats_mpi,ierr)
         call mpi_recv(fbuf,b_ndonor*3,REAL_TYPE,(this_mesh%master_id-1), &
               105,DEFAULT_COMM,stats_mpi,ierr)
         call insert_conn_info(b_nfringe,b_ndonor,imbuf,idbuf,fbuf)
      end if master_check

100   call barrier_mpi
      deallocate(imbuf,idbuf,fbuf)
   contains
      subroutine insert_conn_info(nfr,ndr,imt,idt,frt)
         !! Setup connectivity information on local arrays (in each
         !! processor) and transform global indices to local mesh
         !! indices.
         integer :: nfr, ndr
         integer, dimension(:), intent(in) :: imt, idt
         real(kind=rdp), dimension(:), intent(in) :: frt
         
         do i=1,nfr
            i1=3*(i-1)
            do j=1,3
               b_imesh(j,i)=imt(i1+j)-jkl_lims(j)+1
            end do
         end do

         do i=1,ndr
            i1=3*(i-1)
            do j=1,3
               b_idon(j,i)=idt(i1+j)-jkl_lims(j)+1
               b_frac(j,i)=frt(i1+j)
            end do
         end do
      end subroutine insert_conn_info
   end subroutine relay_conn_info_to_blocks

   subroutine perform_overset_interpolations(q,vnu)
      real, dimension(jmax,kmax,lmax,6) :: q
      real, dimension(jmax,kmax,lmax) :: vnu

      ! Interpolated data from the mesh blocks
      real(kind=rdp), dimension(:), allocatable :: qblk, vnublk
      real(kind=rdp), dimension(:,:), allocatable :: qmsh, vnumsh
      real(kind=rdp), dimension(:,:,:), allocatable :: qog, vnuog

      !if (o_ngroups == 0) goto 100
      
      allocate(qblk(dsize3*5),vnublk(dsize3))
      if (this_isMaster) then
         allocate(qmsh(dsize3*5,o_ngroups),vnumsh(dsize3,o_ngroups))
      end if
      if (o_nrecvs > 0) then
         allocate(qog(dsize3*5,o_mmax,o_nrecvs),&
               vnuog(dsize3,o_mmax,o_nrecvs))
      end if
      

      call donor_block_interpolations
      call relay_interps_to_master
      call assemble_interpolated_variables
      call relay_fringe_vars_to_master
      call relay_fringe_vars_to_blocks
      call update_fringe_point_solution

      if (allocated(qblk)) deallocate(qblk,vnublk)
      if (allocated(qmsh)) deallocate(qmsh,vnumsh)
      if (allocated(qog)) deallocate(qog,vnuog)

100   call barrier_mpi
   contains
      subroutine donor_block_interpolations
         integer :: j, k, l, jp, kp, lp, i, i1, n
         real(kind=rdp) :: dj0, dk0, dl0, dj1, dk1, dl1
         real(kind=rdp) :: w1, w2, w3, w4, w5, w6, w7, w8

         do i=1,b_ndonor
            i1=(i-1)*5
            j=b_idon(1,i); jp=min(j+1,jmax)
            k=b_idon(2,i); kp=min(k+1,kmax)
            l=b_idon(3,i); lp=min(l+1,lmax)

            dj1=b_frac(1,i); dj0=one-dj1
            dk1=b_frac(2,i); dk0=one-dk1
            dl1=b_frac(3,i); dl0=one-dl1

            w1=(dj0*dk0)*dl0
            w2=(dj1*dk0)*dl0
            w3=(dj0*dk1)*dl0
            w4=(dj1*dk1)*dl0
            w5=(dj0*dk0)*dl1
            w6=(dj1*dk0)*dl1
            w7=(dj0*dk1)*dl1
            w8=(dj1*dk1)*dl1
 
            do n=1,5
               qblk(i1+n)=w1*q(j,k,l,n)*q(j,k,l,6) &
                     +w2*q(jp,k,l,n)*q(jp,k,l,6) &
                     +w3*q(j,kp,l,n)*q(j,kp,l,6) &
                     +w4*q(jp,kp,l,n)*q(jp,kp,l,6) &
                     +w5*q(j,k,lp,n)*q(j,k,lp,6) &
                     +w6*q(jp,k,lp,n)*q(jp,k,lp,6) &
                     +w7*q(j,kp,lp,n)*q(j,kp,lp,6) &
                     +w8*q(jp,kp,lp,n)*q(jp,kp,lp,6)
            end do
            vnublk(i)=w1*vnu(j,k,l) &
                  +w2*vnu(jp,k,l) &
                  +w3*vnu(j,kp,l) &
                  +w4*vnu(jp,kp,l) &
                  +w5*vnu(j,k,lp) &
                  +w6*vnu(jp,k,lp) &
                  +w7*vnu(j,kp,lp) &
                  +w8*vnu(jp,kp,lp)
         end do
      end subroutine donor_block_interpolations

      subroutine relay_interps_to_master
         real(kind=rdp), dimension(:), allocatable :: qtmp, vnutmp
         integer, dimension(:), pointer :: bids
         type(ovset_blk_data), pointer :: blk
         integer :: n, i, i1, j, i2, i0, g

         allocate(qtmp(dsize3*5),vnutmp(dsize3))

         call barrier_mpi
         !!call print_message("Relay interps to master")
         if (o_ngroups == 0) goto 100
         master_check: if (this_isMaster) then
            bids=>this_mesh%block_id
            mesh_blocks: do n=1,this_mesh%nblocks
               blk=>b_conn_info(n)
               local_block: if (bids(n) == (PID+1)) then
                  do g=1,o_ngroups
                     do i=blk%dost(g),blk%doend(g)
                        i0=(i-1)*5
                        i1=blk%doidx(i)
                        i2=(i1-1)*5
                        
                        do j=1,5
                           qmsh(i2+j,g)=qblk(i0+j)
                        end do
                        vnumsh(i1,g)=vnublk(i)
                     end do !i=blk%dost(g),blk%doend(g)
                  end do ! g=1,o_ngroups
               else local_block
                  call mpi_recv(qtmp,blk%ndonor*5,REAL_TYPE,bids(n)-1,101, &
                        DEFAULT_COMM,stats_mpi,ierr)
                  call mpi_recv(vnutmp,blk%ndonor,REAL_TYPE,bids(n)-1,102, &
                        DEFAULT_COMM,stats_mpi,ierr)
                  do g=1,o_ngroups
                     do i=blk%dost(g),blk%doend(g)
                        i0=(i-1)*5
                        i1=blk%doidx(i)
                        i2=(i1-1)*5
                        
                        do j=1,5
                           qmsh(i2+j,g)=qtmp(i0+j)
                        end do
                        vnumsh(i1,g)=vnutmp(i)
                     end do !i=blk%dost(g),blk%doend(g)
                  end do ! g=1,o_ngroups
               end if local_block
            end do mesh_blocks
         else master_check
            call mpi_bsend(qblk,b_ndonor*5,REAL_TYPE,this_mesh%master_id-1, &
                  101,DEFAULT_COMM,ierr)
            call mpi_bsend(vnublk,b_ndonor,REAL_TYPE,this_mesh%master_id-1, &
                  102,DEFAULT_COMM,ierr)
         end if master_check
100      call barrier_mpi
         deallocate(qtmp,vnutmp)
      end subroutine relay_interps_to_master

      subroutine assemble_interpolated_variables
         integer :: ii, jj, kk, mm, mid, o_master
         integer :: i, i1, n, j, j1
         real(kind=rdp), dimension(:,:), allocatable :: qtmp, vnutmp
         type(ovset_data), pointer :: dinfo

         if (o_nrecvs > 0) then
            !I am an overset master
            allocate(qtmp(dsize3*5,o_mmax), &
                  vnutmp(dsize3,o_mmax))
         end if
         
         call barrier_mpi
         if (.not. this_isMaster) go to 100
         ovset_groups: do ii=1,o_ngroups
            jj=o_grpptr(ii)
            kk=this_ovset(ii)

            master_proc: if (jj == 0) then
               dinfo=>m_conn_info(ii)
               o_master=ovset(kk)%master_id-1
               mid=this_mesh%master_id-1
               call mpi_bsend(qmsh(:,ii),dinfo%ndonor*5,REAL_TYPE,o_master,501,&
                     DEFAULT_COMM,ierr)
               call mpi_bsend(vnumsh(:,ii),dinfo%ndonor,REAL_TYPE,o_master,502,&
                     DEFAULT_COMM,ierr)
            else
               ! 1. Assemble interpolated data in the overset master
               group_meshes: do mm=1,ovset(kk)%nmesh
                  mid=ovset(kk)%mesh_ids(mm)
                  o_master=meshes(mid)%master_id-1
                  dinfo=>ovset_grids(jj)%doninfo(mm)

                  remote_mesh: if (mid == this_mesh%mesh_id) then
                     do i=1,dinfo%ndonor
                        i1=(i-1)*5
                        do n=1,5
                           qtmp(i1+n,mm)=qmsh(i1+n,ii)
                        end do
                        vnutmp(i,mm)=vnumsh(i,ii)
                     end do
                  else remote_mesh
                     call mpi_recv(qtmp(1,mm),dinfo%ndonor*5,REAL_TYPE, &
                           o_master,501,DEFAULT_COMM,stats_mpi,ierr)
                     call mpi_recv(vnutmp(1,mm),dinfo%ndonor,REAL_TYPE, &
                           o_master,502,DEFAULT_COMM,stats_mpi,ierr)
                  end if remote_mesh
               end do group_meshes

               ! 2. Redistribute interpolated data to reciever points
               do mm=1,ovset(kk)%nmesh
                  dinfo=>ovset_grids(jj)%doninfo(mm)
                  do i=1,dinfo%nfringe
                     i1=(i-1)*5
                     j=dinfo%imap(1,i)
                     j1=(j-1)*5
                     mid=dinfo%imap(2,i)
                     do n=1,5
                        qog(i1+n,mm,jj)=qtmp(j1+n,mid)
                     end do
                     vnuog(i,mm,jj)=vnutmp(j,mid)
                  end do
               end do
            end if master_proc
         end do ovset_groups

100      call barrier_mpi
      end subroutine assemble_interpolated_variables

      subroutine relay_fringe_vars_to_master
         integer :: ii, jj, kk, mm, mid, o_master
         integer :: i, i1, n
         type(ovset_data), pointer :: dinfo

         call barrier_mpi
         if (.not. this_isMaster) go to 100
         ovset_groups: do ii=1,o_ngroups
            jj=o_grpptr(ii)
            kk=this_ovset(ii)

            master_proc: if (jj == 0) then
               dinfo=>m_conn_info(ii)
               o_master=ovset(kk)%master_id-1
               mid=this_mesh%master_id-1
               call mpi_recv(qmsh(:,ii),dinfo%nfringe*5,REAL_TYPE,o_master, &
                     601,DEFAULT_COMM,stats_mpi,ierr)
               call mpi_recv(vnumsh(:,ii),dinfo%nfringe,REAL_TYPE,o_master, &
                     602,DEFAULT_COMM,stats_mpi,ierr)
            else master_proc
               group_meshes: do mm=1,ovset(kk)%nmesh
                  mid=ovset(kk)%mesh_ids(mm)
                  o_master=meshes(mid)%master_id-1
                  dinfo=>ovset_grids(jj)%doninfo(mm)

                  remote_mesh: if (mid == this_mesh%mesh_id) then
                     do i=1,dinfo%nfringe
                        i1=(i-1)*5
                        do n=1,5
                           qmsh(i1+n,ii)=qog(i1+n,mm,ii)
                        end do
                        vnumsh(i,ii)=vnuog(i,mm,ii)
                     end do
                  else remote_mesh
                     call mpi_bsend(qog(:,mm,ii),dinfo%nfringe*5,REAL_TYPE, &
                           o_master,601,DEFAULT_COMM,ierr)
                     call mpi_bsend(vnuog(:,mm,ii),dinfo%nfringe,REAL_TYPE, &
                           o_master,602,DEFAULT_COMM,ierr)
                  end if remote_mesh
               end do group_meshes
            end if master_proc
         end do ovset_groups
100      call barrier_mpi
      end subroutine relay_fringe_vars_to_master

      subroutine relay_fringe_vars_to_blocks
         real(kind=rdp), dimension(:), allocatable :: qtmp, vnutmp
         integer, dimension(:), pointer :: bids
         type(ovset_blk_data), pointer :: blk
         integer :: n, g, i, i0, i1, i2, j

         allocate(qtmp(dsize3*5),vnutmp(dsize3))

         call barrier_mpi
         if (o_ngroups == 0) goto 100
         !call print_message("Relay fringe vars to blocks")
         master_check: if (this_isMaster) then
            bids=>this_mesh%block_id
            mesh_blocks: do n=1,this_mesh%nblocks
               blk=>b_conn_info(n)
               local_block: if (bids(n) == (PID+1)) then
                  do g=1,o_ngroups
                     do i=blk%imst(g),blk%imend(g)
                        i0=(i-1)*5
                        i1=blk%fridx(i)
                        i2=(i1-1)*5
                        do j=1,5
                           qblk(i0+j)=qmsh(i2+j,g)
                        end do
                        vnublk(i)=vnumsh(i1,g)
                     end do
                  end do
               else local_block
                  do g=1,o_ngroups
                     do i=blk%imst(g),blk%imend(g)
                        i0=(i-1)*5
                        i1=blk%fridx(i)
                        i2=(i1-1)*5
                        do j=1,5
                           qtmp(i0+j)=qmsh(i2+j,g)
                        end do
                        vnutmp(i)=vnumsh(i1,g)
                     end do
                  end do
                  call mpi_bsend(qtmp,blk%nfringe*5,REAL_TYPE,bids(n)-1, &
                        701,DEFAULT_COMM,ierr)
                  call mpi_bsend(vnutmp,blk%nfringe,REAL_TYPE,bids(n)-1, &
                        702,DEFAULT_COMM,ierr)
               end if local_block
            end do mesh_blocks
         else master_check
            call mpi_recv(qblk,b_nfringe*5,REAL_TYPE,this_mesh%master_id-1, &
                  701,DEFAULT_COMM,stats_mpi,ierr)
            call mpi_recv(vnublk,b_nfringe,REAL_TYPE,this_mesh%master_id-1, &
                  702,DEFAULT_COMM,stats_mpi,ierr)
         end if master_check
100      call barrier_mpi
         deallocate(qtmp,vnutmp)
      end subroutine relay_fringe_vars_to_blocks

      subroutine update_fringe_point_solution
         integer :: i, j, k, l, n, i1

         do i=1,b_nfringe
            i1=(i-1)*5

            j=b_imesh(1,i)
            k=b_imesh(2,i)
            l=b_imesh(3,i)

            do n=1,5
               q(j,k,l,n)=qblk(i1+n)/q(j,k,l,6)
            end do
            vnu(j,k,l)=vnublk(i)
         end do
      end subroutine update_fringe_point_solution
   end subroutine perform_overset_interpolations

   subroutine setup_connectivity_ihc
      !! Determine some more data necessary for IHC subroutines.

      integer :: ii, jj, kk, mm, mid, nn
      type(mesh_info), pointer :: cur_mesh
      type(simple_bc), pointer :: sbc

      ! tmp_immune - Immunization limits in L-direction
      ! tmp_te     - Immunization limits JTAIL
      integer, dimension(:), allocatable :: tmp_immune, tmp_te

      nullify(cur_mesh,sbc)

      if ((.not. this_isMaster) .or. (o_ngroups == 0)) return
      allocate(tmp_immune(nmeshes),tmp_te(nmeshes))
      call read_ihc_inputs
      
      ogroups: do ii=1,o_ngroups
         jj=o_grpptr(ii)
         ! I am a slave, I don't need to generate this data
         if (jj == 0) cycle

         kk=this_ovset(ii)
         do mm=1,ovset(kk)%nmesh
            mid=ovset(kk)%mesh_ids(mm)
            cur_mesh=>meshes(mid)
            select case (cur_mesh%dom_type)
            ! Math: add DOM_IBC
            case (DOM_WAKE,DOM_GROUNDWAKE,DOM_IBC)
               ovset_grids(jj)%ibc(mm)=IHC_FARFIELD
               if (tmp_immune(mid) > 0) then
                  ovset_grids(jj)%jhimmune(mm)=tmp_immune(mid)
                  ovset_grids(jj)%jte(mm)=tmp_te(mid)
               else
                  ovset_grids(jj)%jhimmune(mm)=1
                  ovset_grids(jj)%jte(mm)=1
               endif
            ! Math: adding DOM_FUS, not fuse if necessary here...
            case (DOM_BLADE,DOM_SLAT,DOM_FUS)
               ! By default let's assume we get a C-O mesh
               ovset_grids(jj)%ibc(mm)=IHC_CO
               if (tmp_immune(mid) > 0) then
                  ovset_grids(jj)%jhimmune(mm)=tmp_immune(mid)
                  ovset_grids(jj)%jte(mm)=tmp_te(mid)
               else
                  ovset_grids(jj)%jhimmune(mm)=45
                  ovset_grids(jj)%jte(mm)=5
               endif
               do nn=1,cur_mesh%nbc_sim
                  sbc=>cur_mesh%bcinfo(nn)
                  if (sbc%ibtyp == BC_WALL) then
                     !Now data is read from user input file
                     !ovset_grids(jj)%jte(mm)=5!sbc%jbcs

                     if ((sbc%kbcs > 1) .or. &
                         (sbc%kbce < cur_mesh%km)) &
                       ovset_grids(jj)%ibc(mm)=IHC_CH
                     exit
                  end if
               end do
            case default
               call stop_execution('setup_connectivity_ihc', &
                     'Unknown mesh type for IHC connectivity')
            end select
         end do
      end do ogroups

      deallocate(tmp_immune,tmp_te)

   contains
      subroutine read_ihc_inputs
         !! Helper routine to read IHC input parameters 
         integer :: un, n, ntot, tmp1, tmp2
         logical :: has_file
         character(len=128) :: mname
         
         ! Initialize immunization limits to bogus values, so that we
         ! can differentiate whether user has specified these limits
         ! or if we need to initialize them to default values.
         tmp_immune=-1
         tmp_te=-1

         inquire(file=ihc_inp,exist=has_file)

         if (.not. has_file) then
            write(STDOUT,'(A)') "WARNING!! IHC input file not &
                  &found. Using default values for body meshes."
            return
         end if
         un=open_file(ihc_inp)
         read(un,*) ntot
         do n=1,ntot
            read(un,*) mname, tmp1, tmp2
            
            mid=get_global_mesh_id(meshes,mname)
            tmp_immune(mid)=tmp1
            tmp_te(mid)=tmp2
         end do
         close(un)
      end subroutine read_ihc_inputs
   end subroutine setup_connectivity_ihc

   subroutine do_connectivity_ihc
      !! Prepare data structures and call Jaina's connectivity
      !! routines. Process the data returned by IHC routines and
      !! insert them in appropriate mesh metadata structures.

      use ihc
      integer :: ii, jj, kk, mm, nmb, i, j
      type(ovset_group), pointer :: og
      type(ovset_data), pointer :: dinfo
      integer, dimension(:), allocatable :: nfringe, ndonor
      integer, dimension(:,:,:), allocatable :: imesh, idonor
      integer, dimension(:,:,:), allocatable :: imap
      real(kind=rdp), dimension(:,:,:), allocatable :: frac

      ! Math: add IBC
      integer :: mid, mibc

      if ((.not. this_isMaster) .or. (o_ngroups == 0)) goto 100

      ovset_groups: do ii=1,o_ngroups
         jj=o_grpptr(ii)
         kk=this_ovset(ii)

         if (jj == 0) cycle
         og=>ovset_grids(jj)
         nmb=ovset(kk)%nmesh
         allocate(nfringe(nmb),ndonor(nmb))
         allocate(imesh(3,dsize3,nmb),idonor(3,dsize3,nmb),frac(3,dsize3,nmb))
         allocate(imap(2,dsize3,nmb))
         nfringe=0
         ndonor=0

         ! Math: get id of ibc grid
         mibc = -1
         do mm=1,ovset(kk)%nmesh
            mid=ovset(kk)%mesh_ids(mm)
            if (meshes(mid)%dom_type==DOM_IBC) mibc=mm
         end do

!         call do_connectihc(og%xyz,og%iblank,og%jmx,og%kmx,og%lmx,imesh,&
!               idonor,frac,imap,og%ibc,og%jhimmune,og%jte,ndonor,nfringe,&
!               dsize3,m_jm,m_km,m_lm,nmb)
         call do_connectihc(og%xyz,og%iblank,og%jmx,og%kmx,og%lmx,imesh,&
               idonor,frac,imap,ndonor,nfringe,dsize3,m_jm,m_km,m_lm,nmb,&
               mibc) ! Math: add IBC

         group_meshes: do mm=1,nmb
            dinfo=>ovset_grids(jj)%doninfo(mm)
            dinfo%nfringe=nfringe(mm)
            dinfo%ndonor=ndonor(mm)

            do i=1,dinfo%nfringe
               do j=1,3
                  dinfo%imesh(j,i)=imesh(j,i,mm)
               end do
               dinfo%imap(1,i)=imap(1,i,mm)
               dinfo%imap(2,i)=imap(2,i,mm)
            end do

            do i=1,dinfo%ndonor
               do j=1,3
                  dinfo%idonor(j,i)=idonor(j,i,mm)
                  dinfo%frac(j,i)=frac(j,i,mm)
               end do
            end do
         end do group_meshes
         deallocate(nfringe,ndonor,imesh,idonor,frac,imap)
      end do ovset_groups

100   call barrier_mpi
   end subroutine do_connectivity_ihc

!   ! Math: assemble iblanks of IBC blocks on IBC mesh master
!   ! Note: this assumes only one IBC mesh per overset group
!   subroutine assemble_iblk_ibc(iblank)
!      !use immersedBoundVars
!
!      integer, dimension(jmax,kmax,lmax), intent(in) :: iblank
!
!      INTEGER(KIND=idp),ALLOCATABLE,DIMENSION(:,:,:) :: iblk
!      integer :: j, k, l, bufsize, i, tag, n
!      integer :: jj, kk, ll
!      real(kind=rdp), dimension(:), allocatable :: xbuffer
!      integer, dimension(:), pointer :: bids
!      integer, dimension(:,:), pointer :: bst, bend
!   
!      call barrier_mpi
!      if (this_mesh%dom_type.ne.DOM_IBC) goto 100 ! if not IBC block, skip
!      print*,'!!!! There must only be ONE IBC MESH per overset group (make this more generic) !!!'
!   
!      nullify(bids,bst,bend)
!      tag=101
!      bufsize=b_jm*b_km*b_lm
!      allocate(xbuffer(bufsize))
!   
!      master_check: if (this_isMaster) then ! IBC mesh master
!         !allocate(iblk(jmax,kmax,lmax))
!         allocate(iblk(this_mesh%jm,this_mesh%km,this_mesh%lm))
!         bids=>this_mesh%block_id
!         bst=>this_mesh%start
!         bend=>this_mesh%end 
!         mesh_blocks: do n=1,this_mesh%nblocks
!            if (bids(n) == (PID+1)) then
!               !Local so copy data
!                        print*,'1setting on master block',PID,bids(n)-1
!               do l=bst(LDIR,n),bend(LDIR,n)
!                  ll=send_lims(LDIR)+(l-bst(LDIR,n))
!                  do k=bst(KDIR,n),bend(KDIR,n)
!                     kk=send_lims(KDIR)+(k-bst(KDIR,n))
!                     do j=bst(JDIR,n),bend(JDIR,n)
!                        jj=send_lims(JDIR)+(j-bst(JDIR,n))
!                        !iblk(j,k,l)=iblk_ibc(jj,kk,ll)
!                        iblk(j,k,l)=iblank(jj,kk,ll)
!                     end do
!                  end do
!               end do
!                        print*,'set on master block',PID,bids(n)-1
!            else ! Receive data from remote processor
!               print*,'receiving on / from',PID,bids(n)-1
!               call mpi_recv(xbuffer,bufsize,INT_TYPE,bids(n)-1,tag,&
!                     DEFAULT_COMM,stats_mpi,ierr)
!               print*,'received on / from',PID,bids(n)-1
!               i=0
!               do l=bst(LDIR,n),bend(LDIR,n)
!                  do k=bst(KDIR,n),bend(KDIR,n)
!                     do j=bst(JDIR,n),bend(JDIR,n)
!                        iblk(j,k,l)=xbuffer(i+1)
!                        i=i+1
!                     end do
!                  end do
!              end do
!            end if
!         end do mesh_blocks
!      else  master_check ! We will send information
!         i=0
!         print*,'preparing to send on',PID
!         do l=send_lims(LDIR),send_lims(LDIR+3)
!            do k=send_lims(KDIR),send_lims(KDIR+3)
!               do j=send_lims(JDIR),send_lims(JDIR+3)
!                  !xbuffer(i+1)=iblk_ibc(j,k,l)
!                  xbuffer(i+1)=iblank(j,k,l)
!                  i=i+1
!               end do
!            end do
!         end do
!         print*,'sending from / to',PID,this_mesh%master_id-1
!         call mpi_bsend(xbuffer,bufsize,INT_TYPE,(this_mesh%master_id-1),&
!               tag,DEFAULT_COMM,ierr)
!         print*,'sent from / to',PID,this_mesh%master_id-1
!      end if master_check
!
!      deallocate(xbuffer)
!
!      print*,'assemble OK'
!100   call barrier_mpi
!
!      call send_iblk_ibc(iblk)
!      if(allocated(iblk)) deallocate(iblk)
!
!   end subroutine assemble_iblk_ibc
!
!   ! Math: add IBC
!   subroutine send_iblk_ibc(iblk)
!      use immersedBoundVars
!
!      INTEGER(KIND=idp),DIMENSION(:,:,:), intent(in) :: iblk
!      integer :: ii, jj, kk, mm, mid, o_master
!
!      integer :: bufsize, i, j, k, l, tag
!      real(kind=rdp), dimension(:), allocatable :: xbuffer
!
!      integer :: has_ibc
!
!      call barrier_mpi
!      if (.not. this_isMaster) go to 100
!
!      tag=201
!      bufsize=m_jm*m_km*m_lm
!      allocate(xbuffer(bufsize))
!
!      ovset_groups: do ii=1,o_ngroups
!         jj=o_grpptr(ii)
!         kk=this_ovset(ii)
!         o_master=ovset(kk)%master_id-1 ! overset master
!
!         ! Check if overset groups has IBC mesh
!         has_ibc=0
!         do mm=1,ovset(kk)%nmesh
!            mid=ovset(kk)%mesh_ids(mm)
!            if (meshes(mid)%dom_type==DOM_IBC) has_ibc=1
!         enddo
!         if (has_ibc.eq.0) go to 100
!
!         master_proc: if (this_mesh%dom_type==DOM_IBC) then
!            ! I am the IBC mesh master, I should send iblk to overset master (if I'm not it)
!            ovset_master: if (PID.ne.o_master) then
!               !! Copy iblk into the MPI send buffer
!               i=0
!               do l=1,this_mesh%lm
!                  do k=1,this_mesh%km
!                     do j=1,this_mesh%jm
!                        xbuffer(i+1)=iblk(j,k,l)
!                        i=i+1
!                     end do
!                  end do
!               end do
!               bufsize=this_mesh%jm*this_mesh%km*this_mesh%lm
!               call mpi_bsend(xbuffer,bufsize,INT_TYPE,o_master,tag,&
!                     DEFAULT_COMM,ierr)
!            else ovset_master
!               print*,'setting iblk on overset master (==IBC master)'
!               if(allocated(iblk_ibc)) deallocate(iblk_ibc)
!               allocate(iblk_ibc(this_mesh%jm,this_mesh%km,this_mesh%lm))
!               iblk_ibc=iblk
!               print*,'set iblk on overset master (==IBC master)'
!            end if ovset_master
!         else master_proc
!            ovset_master2: if (PID==o_master) then ! I'm the overset master but not IBC mesh master, I have to receive iblk of IBC mesh
!               group_meshes: do mm=1,ovset(kk)%nmesh
!                  mid=ovset(kk)%mesh_ids(mm)
!                  o_master=meshes(mid)%master_id-1
!                  ibc_mesh: if (meshes(mid)%dom_type==DOM_IBC) then ! this is the IBC mesh master, receive iblk from it
!                     bufsize=meshes(mid)%jm*meshes(mid)%km*meshes(mid)%lm
!                     call mpi_recv(xbuffer,bufsize,INT_TYPE,o_master,tag, &
!                           DEFAULT_COMM,stats_mpi,ierr)
!                     !! Copy iblk from the MPI send buffer
!                     if(allocated(iblk_ibc)) deallocate(iblk_ibc)
!                     allocate(iblk_ibc(meshes(mid)%jm,meshes(mid)%km,meshes(mid)%lm))
!                     i=0
!                     do l=1,meshes(mid)%lm
!                        do k=1,meshes(mid)%km
!                           do j=1,meshes(mid)%jm
!                              iblk_ibc(j,k,l)=xbuffer(i+1)
!                              i=i+1
!                           end do
!                        end do
!                     end do
!                  end if ibc_mesh
!               end do group_meshes
!            end if ovset_master2
!         end if master_proc
!      end do ovset_groups
!
!      deallocate(xbuffer)
!100   call barrier_mpi
!   end subroutine send_iblk_ibc

end module domain_connectivity
