module airloads
!!! This module provides an airloads exchange interface between the
!!! driver code and the parallelized overset code. The subroutines in
!!! this module parse through the mesh metadata and identify the
!!! meshes which have airloads computations performed, and then at the
!!! end of the computations aggregate all the mesh airloads into the
!!! root processor. This allows the driver code to query airloads data
!!! for any specific mesh via proper subroutines.
   
   use ioutils
   use mpi_wrapper
   use domain_info
   use params_global
   use pyMOD
   use domain_connectivity

   implicit none

   ! Dimension of the airloads array 
   integer, parameter, private :: NFXMOM=7

   ! nmloads - Number of meshes that have airloads associated with
   !           them (for efficiency reasons)
   ! mindex - Array (size: nmesh) storing index in the airloads array
   ! nkdims - No. of spanwise stations on the blade/slat for each mesh
   ! spdist - Radial positions for each mesh
   ! for_mom - forces and moments for all the meshes 
   integer, private :: nmloads
   integer, dimension(:), allocatable, private :: mindex
   integer, dimension(:), allocatable, private :: nkdims
   real(kind=rdp), dimension(:,:), allocatable, private :: spdist
   real(kind=rdp), dimension(:,:,:,:), allocatable, private :: for_mom
   
contains
   subroutine initialize_airloads_data
      !! Loop through the meshes and determine which meshes have
      !! airloads associated with them (currently BLADE/SLAT
      !! meshes). Based on these meshes, initialize the working arrays
      !! that will be used to store data in the zeroth processor for
      !! interfacing with Python etc.
      integer :: n, ii

      ! We only do this on PID==0
      if (.not. isMaster) goto 100
      allocate(mindex(nmeshes))
      mindex=0
      nmloads=0
      ii=1
      do n=1,nmeshes
         if (check_wing_domain_type(meshes(n)%dom_type)) then
            mindex(n)=ii
            ii=ii+1
            nmloads=nmloads+1
         endif
      end do
      allocate(nkdims(nmloads),spdist(m_km,nmloads))
100   call barrier_mpi
   end subroutine initialize_airloads_data
   
   subroutine assemble_surface_forces(mforce,nsteps)
      !! For each participating mesh, assemble the full spanwise
      !! loading from its sub-blocks.
      real(kind=rdp), dimension(:,:,:), intent(inout) :: mforce
      integer, intent(in) :: nsteps
      
      integer :: i, j, k,l, n, bufsize
      integer :: jj, kk, ll, kr1, kt1
      real(kind=rdp), dimension(:), allocatable :: fbuf
      integer, dimension(:), pointer :: bids
      integer, dimension(:,:), pointer :: bst, bend

      call barrier_mpi
      if (.not. is_wing) goto 100

      nullify(bids,bst,bend)
      bufsize=b_km*nsteps*NFXMOM
      allocate(fbuf(bufsize))

      master_check: if (this_isMaster) then
         bids=>this_mesh%block_id
         bst=>this_mesh%start
         bend=>this_mesh%end
         mesh_blocks: do n=1,this_mesh%nblocks
            kr1=max(m_kroot,bst(KDIR,n))
            kt1=min(m_ktip,bend(KDIR,n))
            if (bids(n) == (PID+1)) then
               do l=1,nsteps
                  do k=kr1,kt1
                     kk=send_lims(KDIR)+(k-bst(KDIR,n))
                     do j=1,NFXMOM
                        mforce(j,k,l)=loads(j,kk,l)
                     end do
                  end do
               end do
            else ! Receive data from remote processor
               call mpi_recv(fbuf,bufsize,REAL_TYPE,bids(n)-1,701, &
                     DEFAULT_COMM,stats_mpi,ierr)
               i=1
               do l=1,nsteps
                  do k=kr1,kt1
                     do j=1,NFXMOM
                        mforce(j,k,l)=fbuf(i)
                        i=i+1
                     end do
                  end do
               end do
            end if
         end do mesh_blocks
      else master_check
         kr1=max(kroot,send_lims(KDIR))
         kt1=min(ktip,send_lims(KDIR+3))
         i=1
         do l=1,nsteps
            do k=kr1,kt1
               do j=1,NFXMOM
                  fbuf(i)=loads(j,k,l)
                  i=i+1
               end do
            end do
         end do
         call mpi_bsend(fbuf,bufsize,REAL_TYPE,(this_mesh%master_id-1), &
               701,DEFAULT_COMM,ierr)
      end if master_check

      deallocate(fbuf)
100   call barrier_mpi
   end subroutine assemble_surface_forces

   subroutine collect_mesh_forces(nsteps)
      !! Collect the spanwise airloads information from all meshes in
      !! the root processor.
      integer, intent(in) :: nsteps

      integer :: i, j, k, l, n, bufsize, idx, mid
      real, dimension(:,:,:), allocatable :: mforce
      real, dimension(:), allocatable :: mfbuf

      bufsize=NFXMOM*m_km*nsteps
      allocate(mforce(NFXMOM,m_km,nsteps))
      call assemble_surface_forces(mforce,nsteps)

      if ((.not. is_wing) .or. &
          (.not. this_isMaster)) goto 100

      allocate(mfbuf(bufsize))
      if (isMaster) then ! PID == 0
         if (allocated(for_mom)) deallocate(for_mom)
         allocate(for_mom(NFXMOM,m_km,nsteps,nmloads))
         mesh_loop: do n=1,nmeshes
            if (mindex(n) == 0) cycle

            idx=mindex(n)
            mid=meshes(n)%master_id-1
            if (mid == MASTER) then
               nkdims(idx)=m_nsa
               do k=1,m_nsa
                  spdist(k,idx)=m_span(k+m_kroot-1)
               end do
               do l=1,nsteps
                  do k=1,m_nsa
                     do j=1,NFXMOM
                        for_mom(j,k,l,idx)=mforce(j,k,l)
                     end do
                  end do
               end do
            else
               call mpi_recv(nkdims(idx),1,INT_TYPE,mid,701, &
                     DEFAULT_COMM,stats_mpi,ierr)
               call mpi_recv(spdist(:,idx),nkdims(idx),REAL_TYPE,mid,702, &
                     DEFAULT_COMM,stats_mpi,ierr)
               call mpi_recv(mfbuf,bufsize,REAL_TYPE,mid,703, &
                     DEFAULT_COMM,stats_mpi,ierr)
               i=1
               do l=1,nsteps
                  do k=1,nkdims(idx)
                     do j=1,NFXMOM
                        for_mom(j,k,l,idx)=mfbuf(i)
                        i=i+1
                     end do
                  end do
               end do
            end if
         end do mesh_loop
      else ! Not zeroth processor so we send data
         i=1
         do l=1,nsteps
            do k=1,m_nsa
               do j=1,NFXMOM
                  mfbuf(i)=mforce(j,k,l)
                  i=i+1
               end do
            end do
         end do
         call mpi_bsend(m_nsa,1,INT_TYPE,MASTER,701,DEFAULT_COMM,ierr)
         call mpi_bsend(m_span,m_nsa,REAL_TYPE,MASTER,702,DEFAULT_COMM,ierr)
         call mpi_bsend(mfbuf,bufsize,REAL_TYPE,MASTER,703,DEFAULT_COMM,ierr)
      end if

      deallocate(mfbuf)
100   call barrier_mpi
      deallocate(mforce)
   end subroutine collect_mesh_forces

   pure function get_num_span_loc(msh_id) result(nsa)
      !! Return the number of spanwise segments for a given mesh
      !!
      !! Inputs:
      !!   msh_id - Mesh ID
      !! Outputs:
      !!   nsa - Number of spanwise sections (KTIP-KROOT+1)
      integer, intent(in) :: msh_id
      integer :: nsa

      nsa=nkdims(mindex(msh_id))
   end function get_num_span_loc

   subroutine send_airloads(mid,nstps,rdist,fmom)
      !! Return the radial positions and the corresponding airloads
      !! information for a given mesh ID.
      !!
      !! Inputs:
      !!    mid - Mesh ID
      !!    nstps - Number of timesteps
      !!
      !! Outputs:
      !!    rdist - Radial positions of the segments
      !!    fmom - Force/Moment array (7,nsa,nstps)
      integer, intent(in) :: mid,nstps
      real(kind=rdp), dimension(:), intent(out) :: rdist
      real(kind=rdp), dimension(:,:,:), intent(out) :: fmom

      integer :: j, k, l, idx

      idx=mindex(mid)
      do l=1,nstps
         do k=1,nkdims(idx)
            rdist(k)=spdist(k,idx)
            do j=1,NFXMOM
               fmom(j,k,l)=for_mom(j,k,l,idx)
            end do
         end do
      end do
   end subroutine send_airloads
end module airloads

