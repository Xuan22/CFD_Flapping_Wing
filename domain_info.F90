module domain_info
!!! This module contains the data and the routines that manipulate the
!!! data for the global mesh information. It contains the metadata of
!!! all the meshes participating in a particular problem, the
!!! partitioning information. In addition, this module also contains
!!! additional data about the base global mesh whose partitioned chunk
!!! this processor is solving. The computational block information in
!!! params_global is initialized using the data in this module.

   use constants
   use mpi_wrapper
   use ioutils
   use mesh_types
   use domain_partitioning
   use plot3d
   use params_global
   use io_filenames
   use bcparam

   implicit none

   ! NGHOST - Number of overlapping points on one face of the mesh
   ! NOVERLAP - Total no. of overlapping points in a given direction
   integer, parameter :: NGHOST=3
   integer, parameter :: NOVERLAP=2*NGHOST+1

   ! meshes - Metadata for all the physical meshes involved in the
   !          computation
   ! mblks - Metadata concerning all the computational blocks after
   !         partitioning
   ! ovset - Metadata for overset connectivity between meshes
   ! nmeshes - No. of meshes
   ! nmblks - Total no. of computational blocks
   ! novset - Total no. of overset connectivities
   type(mesh_info), dimension(:), target, allocatable :: meshes
   type(mesh_block), dimension(:), target, allocatable :: mblks
   type(overset_info), dimension(:), target, allocatable :: ovset
   integer :: nmeshes, nmblks, novset

   ! this_* - pointers to metadata relevant to this processor
   ! this_mesh - Metadata for the mesh whose computational chunk this
   !             processor is solving.
   ! this_block - Metadata for the computational block on this
   !              processor
   type(mesh_info), pointer :: this_mesh
   type(mesh_block), pointer :: this_block
   logical :: this_ismaster

   ! m_* - Values of the global mesh, the meanings are the same as
   ! that for hte params_global variables without the leading m_
   integer :: m_kroot, m_ktip, m_nsa, m_jtail1, m_jle, m_jtail2, m_root_co
   real, dimension(:), allocatable :: m_span, m_twist
#ifdef USE_SINGLE_PRECISION
   type(p3dsingle), target :: m_grid
#else
   type(p3d), target :: m_grid
#endif
   integer, dimension(:,:,:), allocatable :: m_iblank
   real(kind=rdp), dimension(:,:,:,:), allocatable :: m_q
   real(kind=rdp), dimension(:,:,:), allocatable :: m_vnu

   ! b_jm, b_km, b_lm - block buffer dimensions
   ! m_jm, m_km, m_lm - mesh buffer dimensions
   integer :: b_jm, b_km, b_lm
   integer :: m_jm, m_km, m_lm

   ! jkl_lims - {J,k,L}MIN/MAX info including ghost points
   ! send_lims - limits for sending data to master processor
   integer, dimension(6) :: jkl_lims, send_lims

   ! psi_offset - Azimuthal position of the grids
   real(kind=rdp) :: psi_offset

contains
   subroutine setup_flow_domain_info
      !! The basic entry point for this module. This routine reads in
      !! the global mesh information, calls the domain partitioning
      !! routine to determine the partitioning information, then sets
      !! up pointers to the mesh and the block relevant to this
      !! processor. It then proceeds to compute array dimensioning
      !! variables for this processor.
      integer :: un, i, j
      integer, dimension(:,:), allocatable :: b_dims, m_dims

      ! Read in meshes and partition them
      !call set_basedir("grids/")
      call read_mesh_info(meshes,mesh_inp)
      select case (split_type)
      case (SPLIT_AUTO)
         call determine_proc_alloc(meshes,NPROCS,mblks,palloc_file)
      case (SPLIT_MANUAL_3D)
         call manual_split_3d(meshes,NPROCS,mblks,palloc_file)
      end select
      call read_overset_info(ovset,meshes,overset_inp)
      !call set_basedir(' ')

      nmeshes=size(meshes,1)
      nmblks=size(mblks,1)
      novset=size(ovset,1)

      ! Set up pointers to metadata in this processor
      nullify(this_mesh,this_block)
      this_block=>mblks(PID+1)
      this_mesh=>meshes(this_block%mesh_id)
      write(STDOUT,2001) trim(this_block%block_name),trim(this_mesh%mesh_name)
2001  format(/,"Solving mesh block ",A," of grid ",A)
      ! Print out partitioning information for debugging
      if (isMASTER) then
         un=open_file("meshes.out")
         do i=1,nmeshes
            call print_mesh_info(meshes(i),un)
         end do
         close(un)
         
         un=open_file("blocks.out")
         do i=1,nmblks
            call print_mesh_block_info(mblks(i),un)
         end do
         close(un)
      end if
      ! Set up domain type of this block
      dom_type=this_block%dom_type

      ! Check if this is a wing type
      is_wing=check_wing_domain_type(dom_type)

      ! Check if this is a ground wake type
      is_ground=check_groundwake_domain_type(dom_type)

      ! Check if this is a body type
      is_body=check_body_domain_type(dom_type)
      
      ! Setup dimensions pertinent to this processor
      call determine_jkl_limits
      call determine_span_chord_limits
      jmax=jkl_lims(JDIR+3)-jkl_lims(JDIR)+1
      kmax=jkl_lims(KDIR+3)-jkl_lims(KDIR)+1
      lmax=jkl_lims(LDIR+3)-jkl_lims(LDIR)+1

      ! This will be reset later when checking BC types
      jtail1=1
      jtail2=jmax
      jle=0
      kroot=1
      ktip=kmax

      ! Variables necessary for array allocations
      gsize=jmax*kmax*lmax
      qsize=gsize*6
      ssize=gsize*5
      size2d=jmax*kmax
      size2dp=lmax*kmax*2
      size1d=kmax
      mdim=max(jmax,kmax,lmax)

      ! Read in azimuthal position of the meshes 
      call read_mesh_azimuth_positions

      ! If this is the master ID for this mesh, setup global grid
      ! arrays
      if ((PID+1) == this_mesh%master_id) then 
         this_ismaster=.true.
         m_grid%jm=this_mesh%jm
         m_grid%km=this_mesh%km
         m_grid%lm=this_mesh%lm
         call allocate_grid_arrays(m_grid)
         allocate(m_iblank(this_mesh%jm,this_mesh%km,this_mesh%lm))
         allocate(m_q(this_mesh%jm,this_mesh%km,this_mesh%lm,nv))
         allocate(m_vnu(this_mesh%jm,this_mesh%km,this_mesh%lm))
         m_iblank=1
      else
         this_ismaster=.false.
      end if
      
      ! Determine the buffer size we will need to assemble chunks
      allocate(b_dims(this_mesh%nblocks,3),m_dims(nmeshes,3))
      do i=1,this_mesh%nblocks
         do j=1,3
            b_dims(i,j)=this_mesh%end(j,i)-this_mesh%start(j,i)+1
         end do
      end do
      b_jm=maxval(b_dims(:,JDIR))+NOVERLAP
      b_km=maxval(b_dims(:,KDIR))+NOVERLAP
      b_lm=maxval(b_dims(:,LDIR))+NOVERLAP

      do i=1,nmeshes
         m_dims(i,JDIR)=meshes(i)%jm
         m_dims(i,KDIR)=meshes(i)%km
         m_dims(i,LDIR)=meshes(i)%lm
      end do
      m_jm=maxval(m_dims(:,JDIR))
      m_km=maxval(m_dims(:,KDIR))
      m_lm=maxval(m_dims(:,LDIR))

      deallocate(b_dims,m_dims)
   end subroutine setup_flow_domain_info

   subroutine determine_jkl_limits
      !! Determine the {j,k,l}min/max information for this
      !! computational block w.r.t. the base mesh. This information is
      !! already available in the blocks metadata, however, we need to
      !! include the ghost cell information to determine the overlap
      !! we need.

      ! Array storing the block limits for easier manipulation
      integer, dimension(6) :: extents
      ! The J,K,L MAX values for the base mesh
      integer, dimension(3) :: maxes
      integer :: i
      
      jkl_lims=0
      extents=(/ this_block%gjs, this_block%gks, this_block%gls, &
            this_block%gje, this_block%gke, this_block%gle /)

      maxes=(/ this_mesh%jm, this_mesh%km, this_mesh%lm /)

      ! The min faces
      do i=1,3
         if (extents(i) == 1) then 
            jkl_lims(i) = extents(i)
         else
            jkl_lims(i) = extents(i)-NGHOST
         end if
         send_lims(i)=1+(extents(i)-jkl_lims(i))
      end do
      
      ! The max faces
      do i=4,6
         if (extents(i) == maxes(i-3)) then 
            jkl_lims(i) = extents(i)
         else
            jkl_lims(i) = extents(i)+NGHOST
         end if
         send_lims(i)=extents(i)-jkl_lims(i-3)+1
      end do

   end subroutine determine_jkl_limits

   subroutine setup_local_grid(x,y,z,iblank,grbc,twist,rdist)
      !! Read in the grid coordinates from the plot3d file and
      !! initialize the computational block coordinates used by the
      !! solver subroutines. Also initializes the iblank and the BC
      !! data used by the solver.

      use pymod, only: nsa
      
      ! Grid coordinates used by the solver
      real, dimension(jmax,kmax,lmax), intent(inout) :: x, y, z
      ! I-blank information used by the solver
      integer, dimension(jmax,kmax,lmax), intent(inout) :: iblank
      ! BC data used by the BC subroutine
      type(bc_t), intent(inout) :: grbc
      ! Twist for the computational chunk
      real, dimension(kmax), intent(inout) :: twist, rdist

      ! Plot-3D data 
#ifdef USE_SINGLE_PRECISION
      type(p3dsingle) :: griddata
#else
      type(p3d) :: griddata
#endif
      integer :: j, k, l, j1, k1, l1, i, ibc
      type(simple_bc), pointer :: sbc
      type(interface_bc), pointer :: ifbc

      nullify(sbc,ifbc)

      ! Read in grid data from the input file
      !call set_basedir("grids/")
      call read_p3d(griddata,this_mesh%grid_file)
      !call set_basedir(' ')

      ! Make sure that the dimensions are correct
      if (.not. check_mesh_dimensions(griddata)) &
         call stop_execution('setup_local_grid', &
            "Dimensions don't match inputs for "//trim(adjustl(this_mesh%grid_file)))
      
      ! Extract the block pertinent to this processor
      do l=jkl_lims(LDIR),jkl_lims(LDIR+3)
         l1=l-jkl_lims(LDIR)+1
         do k=jkl_lims(KDIR),jkl_lims(KDIR+3)
            k1=k-jkl_lims(KDIR)+1
            do j=jkl_lims(JDIR),jkl_lims(JDIR+3)
               j1=j-jkl_lims(JDIR)+1
               x(j1,k1,l1)=griddata%x(j,k,l)
               y(j1,k1,l1)=griddata%y(j,k,l)
               z(j1,k1,l1)=griddata%z(j,k,l)
               iblank(j1,k1,l1)=1
            end do
         end do
      end do

      ! Initialize the BC data for this processor
      call init_bc(grbc,this_block%nbc_sim+this_block%nbc_int)
      ibc=0

      ! Copy over the simple BC information first
      do i=1,this_block%nbc_sim
         ibc=ibc+1
         sbc=>this_block%bcinfo(i)
         grbc%ibtyp(ibc)=sbc%ibtyp
         grbc%ibdir(ibc)=sbc%ibdir
         grbc%jbcs(ibc)=sbc%jbcs+send_lims(JDIR)-1
         grbc%kbcs(ibc)=sbc%kbcs+send_lims(KDIR)-1
         grbc%lbcs(ibc)=sbc%lbcs+send_lims(LDIR)-1
         grbc%jbce(ibc)=sbc%jbce+send_lims(JDIR)-1
         grbc%kbce(ibc)=sbc%kbce+send_lims(KDIR)-1
         grbc%lbce(ibc)=sbc%lbce+send_lims(LDIR)-1
         grbc%ibproc(ibc)=0

         if (sbc%ibtyp == BC_WALL) then 
            ! THIS IS A HACK... BUT SO IS JTAIL1
            jtail1=grbc%jbcs(ibc)
            jtail2=grbc%jbce(ibc)
            jle=(jtail2+jtail1)/2
         end if
      end do

      ! Now process the interface information
      do i=1,this_block%nbc_int
         ibc=ibc+1
         ifbc=>this_block%ifaceinfo(i)
         select case (ifbc%ibtyp)
         case (BC_INTERNAL)
            ! This needs to be tweaked a bit over what is already done
            ! by the domain partitioning code, first we make it
            ! INTERFACE BC and then adjust the end points to include
            ! the ghost points
            grbc%ibtyp(ibc)=BC_INTERFACE
            grbc%ibdir(ibc)=ifbc%ibdir
            grbc%ibproc(ibc)=ifbc%donorid-1
            grbc%jbcs(ibc)=ifbc%jbcs
            grbc%kbcs(ibc)=ifbc%kbcs
            grbc%lbcs(ibc)=ifbc%lbcs
            grbc%jbce(ibc)=ifbc%jbce
            grbc%kbce(ibc)=ifbc%kbce
            grbc%lbce(ibc)=ifbc%lbce
            
            select case (ifbc%ibdir)
            case (JDIR)
               grbc%jbcs(ibc)=1
               grbc%jbce(ibc)=NGHOST
            case (KDIR)
               grbc%kbcs(ibc)=1
               grbc%kbce(ibc)=NGHOST
            case (LDIR)
               grbc%lbcs(ibc)=1
               grbc%lbce(ibc)=NGHOST
            case (-JDIR)
               grbc%jbcs(ibc)=ifbc%jbcs+1+send_lims(JDIR)-1
               grbc%jbce(ibc)=ifbc%jbce+NGHOST+send_lims(JDIR)-1
            case (-KDIR)
               grbc%kbcs(ibc)=ifbc%kbcs+1+send_lims(KDIR)-1
               grbc%kbce(ibc)=ifbc%kbce+NGHOST+send_lims(KDIR)-1
            case (-LDIR)
               grbc%lbcs(ibc)=ifbc%lbcs+1+send_lims(LDIR)-1
               grbc%lbce(ibc)=ifbc%lbce+NGHOST+send_lims(LDIR)-1
            end select

         case default
            ! For all other cases we can just copy over the information.
            grbc%ibtyp(ibc)=ifbc%ibtyp
            grbc%ibdir(ibc)=ifbc%ibdir
            grbc%ibproc(ibc)=ifbc%donorid-1
            grbc%jbcs(ibc)=ifbc%jbcs+send_lims(JDIR)-1
            grbc%kbcs(ibc)=ifbc%kbcs+send_lims(KDIR)-1
            grbc%lbcs(ibc)=ifbc%lbcs+send_lims(LDIR)-1
            grbc%jbce(ibc)=ifbc%jbce+send_lims(JDIR)-1
            grbc%kbce(ibc)=ifbc%kbce+send_lims(KDIR)-1
            grbc%lbce(ibc)=ifbc%lbce+send_lims(LDIR)-1
         end select
      end do

      do i=1,ibc
         if (grbc%jbcs(i) == 1+NGHOST) grbc%jbcs(i)=1
         if (grbc%jbce(i) == jmax-NGHOST) grbc%jbce(i)=jmax
         if (grbc%kbcs(i) == 1+NGHOST) grbc%kbcs(i)=1
         if (grbc%kbce(i) == kmax-NGHOST) grbc%kbce(i)=kmax
         if (grbc%lbcs(i) == 1+NGHOST) grbc%lbcs(i)=1
         if (grbc%lbce(i) == lmax-NGHOST) grbc%lbce(i)=lmax

         if (grbc%ibtyp(i) == BC_WALL) then
            kroot=grbc%kbcs(i)
            ktip=grbc%kbce(i)
            nsa=ktip-kroot+1
         end if
      end do

      write(STDOUT,'(A)') "Boundary Conditions for block: "//trim(this_block%block_name)
      write(STDOUT,'(10(I6,1x))') grbc%ibtyp
      write(STDOUT,'(10(I6,1x))') grbc%ibdir
      write(STDOUT,'(10(I6,1x))') grbc%jbcs
      write(STDOUT,'(10(I6,1x))') grbc%jbce
      write(STDOUT,'(10(I6,1x))') grbc%kbcs
      write(STDOUT,'(10(I6,1x))') grbc%kbce
      write(STDOUT,'(10(I6,1x))') grbc%lbcs
      write(STDOUT,'(10(I6,1x))') grbc%lbce
      write(STDOUT,'(10(I6,1x),/)') grbc%ibproc
      
      ! Setup the global twist and span data
      call setup_span_twist_dist(griddata)
      ! Now initialize the arrays used by the solver
      do k=jkl_lims(KDIR),jkl_lims(KDIR+3)
         k1=k-jkl_lims(KDIR)+1
         twist(k1)=m_twist(k)
         rdist(k1)=m_span(k)
      end do
      call free_grid_arrays(griddata)
   end subroutine setup_local_grid

   subroutine read_mesh_azimuth_positions
      integer :: un, i, mid, rottmp, ilimtmp
      real(kind=rdp) :: psitmp
      character(len=32) :: gnam

      !call set_basedir("grids/")
      ! Math: read ilim for each mesh in mesh_azi.inp
      if (PID==0) print*,'reading mesh_azi.inp: name psi_offset irot_dir ilim'
      un=open_file(mesh_positions,form='formatted',status='old')
      do i=1,nmeshes
         read(un,*) gnam, psitmp, rottmp, ilimtmp
         mid=get_global_mesh_id(meshes,gnam)
         
         if (mid == this_mesh%mesh_id) then
            psi_offset=psitmp*d2r
            irot_dir=rottmp
            rf=rf*irot_dir
            ilim=ilimtmp
            !ideform=ideformtmp
            exit
         end if
      end do
      close(un)
      !call set_basedir(' ')
	
   end subroutine read_mesh_azimuth_positions

   subroutine determine_span_chord_limits
      integer :: i
      type(simple_bc), pointer :: sbc

      root_co=1
      m_root_co=0
      if (.not. is_wing) then 
         m_kroot=1
         m_ktip=this_mesh%km
         m_jtail1=1
         m_jle=1
      else
         do i=1,this_mesh%nbc_sim
            sbc=>this_mesh%bcinfo(i)
            if (sbc%ibtyp == BC_WALL) then
               m_kroot=sbc%kbcs
               m_ktip=sbc%kbce
               m_jtail1=sbc%jbcs
               m_jtail2=sbc%jbce
               m_jle=(m_jtail1+m_jtail2)/2
               if (m_kroot == 1) m_root_co=1
               exit
            end if
         end do
      end if
      m_nsa=m_ktip-m_kroot+1
   end subroutine determine_span_chord_limits

   subroutine setup_span_twist_dist(bgrid)
#ifdef USE_SINGLE_PRECISION
      type(p3dsingle), intent(in) :: bgrid
#else
      type(p3d), intent(in) :: bgrid
#endif
      
      integer :: k
      real(kind=rdp) :: dx, dy, dz, dnorm, denom

      allocate(m_span(bgrid%km),m_twist(bgrid%km))
      m_twist=zero

      ! This kludge is implemented so that C-O meshes have r/R~1.0
      ! when K=KTIP.
      if (.not. is_wing) then 
         m_span=zero
         return
      else if (dom_type == DOM_BLADE) then
         denom=bgrid%y(m_jle,m_ktip,1)
      else if (dom_type == DOM_SLAT) then 
         denom=rartio
      end if

      do k=1,bgrid%km
         dx=bgrid%x(m_jle,k,1)-bgrid%x(m_jtail1,k,1)
         dy=bgrid%y(m_jle,k,1)-bgrid%y(m_jtail1,k,1)
         dz=bgrid%z(m_jle,k,1)-bgrid%z(m_jtail1,k,1)   
         
         dnorm=sqrt(dx**2+dy**2+dz**2)
         m_twist(k)=asin(dz/dnorm)
         m_span(k)=bgrid%y(m_jle,k,1)/denom!bgrid%y(m_jle,m_ktip,1)
      end do
   end subroutine setup_span_twist_dist

   subroutine assemble_mesh(x,y,z)
      real(kind=rdp), dimension(jmax,kmax,lmax) :: x, y, z

      integer :: j, k, l, bufsize, i, tag, n
      integer :: jj, kk, ll
      real(kind=rdp), dimension(:), allocatable :: xbuffer
      integer, dimension(:), pointer :: bids
      integer, dimension(:,:), pointer :: bst, bend

      nullify(bids,bst,bend)
      tag=101
      bufsize=b_jm*b_km*b_lm*3
      allocate(xbuffer(bufsize))

      call barrier_mpi
      master_check: if (this_isMaster) then 
         !call print_message("Assembling mesh: "//this_mesh%mesh_name)
         bids=>this_mesh%block_id
         bst=>this_mesh%start
         bend=>this_mesh%end
         mesh_blocks: do n=1,this_mesh%nblocks
            if (bids(n) == (PID+1)) then 
               !Local so copy data
               do l=bst(LDIR,n),bend(LDIR,n)
                  ll=send_lims(LDIR)+(l-bst(LDIR,n))
                  do k=bst(KDIR,n),bend(KDIR,n)
                     kk=send_lims(KDIR)+(k-bst(KDIR,n))
                     do j=bst(JDIR,n),bend(JDIR,n)
                        jj=send_lims(JDIR)+(j-bst(JDIR,n))
                        m_grid%x(j,k,l)=x(jj,kk,ll)
                        m_grid%y(j,k,l)=y(jj,kk,ll)
                        m_grid%z(j,k,l)=z(jj,kk,ll)
                     end do
                  end do
               end do
            else ! Receive data from remote processor
               call mpi_recv(xbuffer,bufsize,REAL_TYPE,bids(n)-1,tag,&
                     DEFAULT_COMM,stats_mpi,ierr)
               i=0
               do l=bst(LDIR,n),bend(LDIR,n)
                  do k=bst(KDIR,n),bend(KDIR,n)
                     do j=bst(JDIR,n),bend(JDIR,n)
                        m_grid%x(j,k,l)=xbuffer(i+1)
                        m_grid%y(j,k,l)=xbuffer(i+2)
                        m_grid%z(j,k,l)=xbuffer(i+3)
                        i=i+3
                     end do
                  end do
              end do
            end if
         end do mesh_blocks
      else  master_check ! We will send information
         i=0
         do l=send_lims(LDIR),send_lims(LDIR+3)
            do k=send_lims(KDIR),send_lims(KDIR+3)
               do j=send_lims(JDIR),send_lims(JDIR+3)
                  xbuffer(i+1)=x(j,k,l)
                  xbuffer(i+2)=y(j,k,l)
                  xbuffer(i+3)=z(j,k,l)
                  i=i+3
               end do
            end do
         end do
         call mpi_bsend(xbuffer,bufsize,REAL_TYPE,(this_mesh%master_id-1),&
               tag,DEFAULT_COMM,ierr)
      end if master_check

      call barrier_mpi
      deallocate(xbuffer)
!       if (isMaster) then
!          call write_p3d(m_grid,"testgrid1.p3d")
!       endif
   end subroutine assemble_mesh

   subroutine assemble_soln(q,vnu,iblank)
      real(kind=rdp), dimension(jmax,kmax,lmax) :: vnu
      real(kind=rdp), dimension(jmax,kmax,lmax,nd) :: q
      integer, dimension(jmax,kmax,lmax) :: iblank

      integer :: j, k, l, bufsize, i, ii, in, n
      integer :: jj, kk, ll
      real(kind=rdp), dimension(:), allocatable :: qbuf, vbuf
      integer, dimension(:), allocatable :: ibbuf
      integer, dimension(:), pointer :: bids
      integer, dimension(:,:), pointer :: bst, bend

      nullify(bids,bst,bend)
      bufsize=b_jm*b_km*b_lm
      allocate(qbuf(bufsize*nv),vbuf(bufsize),ibbuf(bufsize))

      call barrier_mpi
      master_check: if (this_isMaster) then
         bids=>this_mesh%block_id
         bst=>this_mesh%start
         bend=>this_mesh%end
         mesh_blocks: do n=1,this_mesh%nblocks
            if (bids(n) == (PID+1)) then 
               !Local so copy data
               do l=bst(LDIR,n),bend(LDIR,n)
                  ll=send_lims(LDIR)+(l-bst(LDIR,n))
                  do k=bst(KDIR,n),bend(KDIR,n)
                     kk=send_lims(KDIR)+(k-bst(KDIR,n))
                     do j=bst(JDIR,n),bend(JDIR,n)
                        jj=send_lims(JDIR)+(j-bst(JDIR,n))
                        do in=1,nv
                           m_q(j,k,l,in)=q(jj,kk,ll,in)*q(jj,kk,ll,6)
                        end do
                        m_vnu(j,k,l)=vnu(jj,kk,ll)
                        m_iblank(j,k,l)=iblank(jj,kk,ll)
                     end do
                  end do
               end do
            else ! Receive data from remote processor
               call mpi_recv(qbuf,bufsize*nv,REAL_TYPE,bids(n)-1,100,&
                     DEFAULT_COMM,stats_mpi,ierr)
               call mpi_recv(vbuf,bufsize,REAL_TYPE,bids(n)-1,101,&
                     DEFAULT_COMM,stats_mpi,ierr)
               call mpi_recv(ibbuf,bufsize,INT_TYPE,bids(n)-1,102,&
                     DEFAULT_COMM,stats_mpi,ierr)
               i=1
               ii=0
               do l=bst(LDIR,n),bend(LDIR,n)
                  do k=bst(KDIR,n),bend(KDIR,n)
                     do j=bst(JDIR,n),bend(JDIR,n)
                        do in=1,nv
                           m_q(j,k,l,in)=qbuf(ii+in)
                        end do
                        m_vnu(j,k,l)=vbuf(i)
                        m_iblank(j,k,l)=ibbuf(i)
                        i=i+1
                        ii=ii+nv
                     end do
                  end do
              end do
            end if
         end do mesh_blocks
      else  master_check ! We will send information
         i=1
         ii=0
         do l=send_lims(LDIR),send_lims(LDIR+3)
            do k=send_lims(KDIR),send_lims(KDIR+3)
               do j=send_lims(JDIR),send_lims(JDIR+3)
                  do in=1,nv
                     qbuf(ii+in)=q(j,k,l,in)*q(j,k,l,6)
                  end do
                  vbuf(i)=vnu(j,k,l)
                  ibbuf(i)=iblank(j,k,l)
                  i=i+1
                  ii=ii+nv
               end do
            end do
         end do
         call mpi_bsend(qbuf,bufsize*nv,REAL_TYPE,(this_mesh%master_id-1),&
               100,DEFAULT_COMM,ierr)
         call mpi_bsend(vbuf,bufsize,REAL_TYPE,(this_mesh%master_id-1),&
               101,DEFAULT_COMM,ierr)
         call mpi_bsend(ibbuf,bufsize,INT_TYPE,(this_mesh%master_id-1),&
               102,DEFAULT_COMM,ierr)
      end if master_check

      call barrier_mpi
      deallocate(qbuf,vbuf,ibbuf)
   end subroutine assemble_soln
   
   subroutine store_solutions_blocks(x,y,z,iblank,q,vnu)
      !! Store per-block solutions and the grid for restart purposes
      !!
      real(kind=rdp), dimension(jmax,kmax,lmax) :: x,y,z,vnu
      real(kind=rdp), dimension(jmax,kmax,lmax,nd) :: q
      integer, dimension(jmax,kmax,lmax) :: iblank

      integer :: j, k, l, n
      integer :: logq, logg

      !Remove Jacobian scaling
      do l=1,lmax
         do k=1,kmax
            do j=1,jmax
               do n=1,5
                  q(j,k,l,n)=q(j,k,l,n)*q(j,k,l,6)
               end do
            end do
         end do
      end do

      call set_basedir(soln_dir)
      logg=open_file("g."//trim(adjustl(this_block%block_name)),&
            form='unformatted')
      logq=open_file("q."//trim(adjustl(this_block%block_name)),&
            form='unformatted')

      write(logq) jmax,kmax,lmax
      write(logq) fsmach,alf,rey*(fsmach+fmtip),totime
      write(logq) ((((q(j,k,l,n),j=1,jmax),k=1,kmax),l=1,lmax),n=1,nv)
      write(logq) (((vnu(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax)
      write(logq) istep0

      write(logg) jmax,kmax,lmax
      write(logg) (((x(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax), &
            (((y(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax), &
            (((z(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax), &
            (((iblank(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax)

      close(logg)
      close(logq)
      call set_basedir(' ')

      !Restore Jacobian scaling
      do l=1,lmax
         do k=1,kmax
            do j=1,jmax
               do n=1,5
                  q(j,k,l,n)=q(j,k,l,n)/q(j,k,l,6)
               end do
            end do
         end do
      end do

   end subroutine store_solutions_blocks

   subroutine store_solutions
      !! Store solutions and the grid for restart purposes
      !!

      integer :: j, k, l, n
      integer :: logq, logg
      character(len=32) :: int_string
      ! Math: add store solution with different names
      ! will work for IBC_x, otherwise change A5 in format string
      character(len=1024) :: gfilename,qfilename
      character(len=1024) :: format_string

      if (istep0 < 10) then
          format_string = "(A2,A6,A1,I1)"
      elseif (istep0.ge.10.and.istep0.lt.100) then
          format_string = "(A2,A6,A1,I2)"
      elseif (istep0.ge.100.and.istep0.lt.1000) then
          format_string = "(A2,A6,A1,I3)"
      elseif (istep0.ge.1000.and.istep0.lt.10000) then
          format_string = "(A2,A6,A1,I4)"
      elseif (istep0.ge.10000.and.istep0.lt.100000) then
          format_string = "(A2,A6,A1,I5)"
      endif

      write (gfilename,format_string) "g.",this_mesh%mesh_name,"_",mod(istep0,int(360/(rf*dt)))
      write (qfilename,format_string) "q.",this_mesh%mesh_name,"_",mod(istep0,int(360/(rf*dt)))

      if (restart_counter.ne.0) then
         write(int_string,*) restart_counter
         int_string = '_'//trim(adjustl(int_string))
      else
         int_string = ''
      endif

      ! Note that grid/solution must be previously assembled for this
      ! subroutine to work. The logic is handled in
      ! src/overturns/turns_api.F90:process_outputs.
      if (.not. this_isMaster) goto 100
      call set_basedir(soln_dir)
      !logg=open_file("g."//trim(adjustl(this_mesh%mesh_name))&
      !      //int_string, form='unformatted')
      logg=open_file(trim(adjustl(gfilename)),form='unformatted')
      ! Math: change q filename
      !logq=open_file("q."//trim(adjustl(this_mesh%mesh_name))&
      !      //int_string, form='unformatted')
      logq=open_file(trim(adjustl(qfilename)),form='unformatted')

      write(logq) this_mesh%jm,this_mesh%km,this_mesh%lm
      write(logq) fsmach,alf,rey*(fsmach+fmtip),totime
      write(logq) ((((m_q(j,k,l,n),j=1,this_mesh%jm),k=1,this_mesh%km),&
            l=1,this_mesh%lm),n=1,nv)
      write(logq) (((m_vnu(j,k,l),j=1,this_mesh%jm),k=1,this_mesh%km),&
            l=1,this_mesh%lm)
      write(logq) istep0

      write(logg) this_mesh%jm,this_mesh%km,this_mesh%lm
      write(logg) (((m_grid%x(j,k,l),j=1,this_mesh%jm),k=1,this_mesh%km),&
            l=1,this_mesh%lm), &
            (((m_grid%y(j,k,l),j=1,this_mesh%jm),k=1,this_mesh%km),&
            l=1,this_mesh%lm), &
            (((m_grid%z(j,k,l),j=1,this_mesh%jm),k=1,this_mesh%km),&
            l=1,this_mesh%lm), &
            (((m_iblank(j,k,l),j=1,this_mesh%jm),k=1,this_mesh%km),&
            l=1,this_mesh%lm)

      close(logg)
      close(logq)
      call set_basedir(' ')

100   call barrier_mpi
   end subroutine store_solutions

   subroutine restart_solutions_blocks(q,vnu)
      real(kind=rdp), dimension(jmax,kmax,lmax) :: vnu
      real(kind=rdp), dimension(jmax,kmax,lmax,nd) :: q

      integer :: j, k, l, n
      integer :: logg, logq
      integer :: jtmp, ktmp, ltmp
      real(kind=rdp) :: fs1,al1,re1

      call set_basedir(rest_dir)
      logq=open_file("q."//trim(adjustl(this_block%block_name)),&
            form='unformatted')
      
      read(logq) jtmp,ktmp,ltmp

      if ((jmax /= jtmp) .or. &
          (kmax /= ktmp) .or. &
          (lmax /= ltmp)) then
         call stop_execution('restart_solutions', &
               'Restart grid dimensions not equal to original dimensions.')
      end if
      read(logq) fs1,al1,re1,totime
      read(logq) ((((q(j,k,l,n),j=1,jmax),k=1,kmax),l=1,lmax),n=1,nv)
      read(logq) (((vnu(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax)
      read(logq) istep0

      close(logq)
      call set_basedir(' ')
   end subroutine restart_solutions_blocks
   
   subroutine restart_solutions(q,vnu)
      real(kind=rdp), dimension(jmax,kmax,lmax) :: vnu
      real(kind=rdp), dimension(jmax,kmax,lmax,nd) :: q

      integer :: j, k, l, n
      integer :: j1, k1, l1
      integer :: logg, logq
      integer :: jtmp, ktmp, ltmp
      real(kind=rdp) :: fs1,al1,re1
      real(kind=rdp), dimension(:,:,:,:), allocatable :: qtmp
      real(kind=rdp), dimension(:,:,:), allocatable :: vnutmp

      allocate(qtmp(this_mesh%jm,this_mesh%km,this_mesh%lm,nv))
      allocate(vnutmp(this_mesh%jm,this_mesh%km,this_mesh%lm))
      
      call set_basedir(rest_dir)
      logq=open_file("q."//trim(adjustl(this_mesh%mesh_name)),&
            form='unformatted')
      
      read(logq) jtmp,ktmp,ltmp

      if ((this_mesh%jm /= jtmp) .or. &
          (this_mesh%km /= ktmp) .or. &
          (this_mesh%lm /= ltmp)) then
         call stop_execution('restart_solutions', &
               'Restart grid dimensions not equal to original dimensions.')
      end if
      read(logq) fs1,al1,re1,totime
      read(logq) ((((qtmp(j,k,l,n),j=1,jtmp),k=1,ktmp),l=1,ltmp),n=1,nv)
      read(logq) (((vnutmp(j,k,l),j=1,jtmp),k=1,ktmp),l=1,ltmp)
      read(logq) istep0
      close(logq)
      call set_basedir(' ')

      ! Extract the block pertinent to this processor
      do l=jkl_lims(LDIR),jkl_lims(LDIR+3)
         l1=l-jkl_lims(LDIR)+1
         do k=jkl_lims(KDIR),jkl_lims(KDIR+3)
            k1=k-jkl_lims(KDIR)+1
            do j=jkl_lims(JDIR),jkl_lims(JDIR+3)
               j1=j-jkl_lims(JDIR)+1
               do n=1,nv
                  q(j1,k1,l1,n)=qtmp(j,k,l,n)
               end do
               vnu(j1,k1,l1)=vnutmp(j,k,l)
            end do
         end do
      end do
      
      deallocate(qtmp,vnutmp)
   end subroutine restart_solutions

   pure function check_mesh_dimensions(p3dgrid) result(matched)
      !! Wrapper function to check if mesh dimensions mentioned in the
      !! input file are the same as the ones inferred from the P3D
      !! grid file.
#ifdef USE_SINGLE_PRECISION
      type(p3dsingle), intent(in) :: p3dgrid
#else
      type(p3d), intent(in) :: p3dgrid
#endif
      logical :: matched

      matched=.true.
      if ((this_mesh%jm /= p3dgrid%jm) .or. &
          (this_mesh%km /= p3dgrid%km) .or. &
          (this_mesh%lm /= p3dgrid%lm)) then
         matched=.false.
      endif
   end function check_mesh_dimensions

   function check_wing_domain_type(dtype) result(yesno)
      !! Simple function to determine if the domain type is a blade or
      !! a slat type.
      integer, intent(in) :: dtype
      logical :: yesno

      select case (dtype)
      case (DOM_BLADE, DOM_SLAT)
         yesno=.true.
      case default
         yesno=.false.
      end select
   end function check_wing_domain_type
   function check_groundwake_domain_type(dtype) result(yesno)
      !! Simple function to determine if the mesh simulates ground plane
      integer, intent(in) :: dtype
      logical :: yesno

      select case (dtype)
      case (DOM_GROUNDWAKE)
         yesno=.true.
      case default
         yesno=.false.
      end select
   end function check_groundwake_domain_type

! function added by Camli.b 10/2014, to check the FUS grid type
   function check_body_domain_type(dtype) result(yesno)
      !! Simple function to determine if the domain type is a FUS (INSECT BODY) 
      integer, intent(in) :: dtype
      logical :: yesno

      select case (dtype)
      case (DOM_FUS)
         yesno=.true.
      case default
         yesno=.false.
      end select
   end function check_body_domain_type
end module domain_info
