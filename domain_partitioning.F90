module domain_partitioning
!!! This module implements the domain partitioning code. There are two
!!! strategies implemented: 1. Automated mesh splitting, and 2. User
!!! defined splitting strategy.
!!!
!!! The automated mesh splitting involves the following steps:
!!!
!!! 1. Given N_m meshes, it first determines the relative density of
!!!    each meshes (based on the ratio of number of cells N_c in each
!!!    mesh), and uses this information to calculate an effective
!!!    number of computational blocks N_l (nload). N_l is such that
!!!    each block obtained after domain partitioning will have
!!!    approximately the same number of cells. This means that some
!!!    physical domain meshes will be split more than the other
!!!    meshes. In addition, the code also computes an array of load
!!!    ratios LR whose elements are:
!!!           LR(m)=N_c(m)/min(N_c[:]); m=[1..N_m]
!!!
!!! 2. Based on the N_l, the code then tries to determine how to
!!!    allocate the given number of processors N_p amongst the
!!!    meshes. To simplify the domain splitting strategy, N_p is
!!!    restricted such that N_p=N_l*2^n, where n is a whole number. If
!!!    such a combination is not possible, then it tries to see if
!!!    N_p=N_m*2^n. However, this is not probably the best strategy.
!!!
!!! 3. If the load-balanced split strategy is possible, then N_p is
!!!    split into fractions determined by N_ppm(m)=N_p*LR(m)/N_l;
!!!    m=[1..N_m]. This determines the total number of processors
!!!    allocated to each of the participating meshes.
!!!
!!! 4. The mesh partition strategy involves two steps: (a)
!!!    pre-partitioning, and (b) an octree-like recursive bisection.
!!!
!!! 4a. Pre-partitioning :
!!!       For any given mesh 0<m<N_m, if LR(m)>1 then the mesh is
!!!       first split into LR(m) sections in the direction with the
!!!       maximum number of grid points. At the end of this step, any
!!!       mesh block of mesh m has approximately the same number of
!!!       cells as the coarsest participating mesh.
!!!
!!! 4b. Recursive bisection:
!!!       If LR(m)=1, or after the pre-bisection split, the mesh
!!!       blocks are then recursively bisected along the optimum
!!!       direction upto n levels, where n=log(N_p/N_m)/log(2). The
!!!       optimum direction is chosen such that the number of new
!!!       faces (normal to that direction) introduced by the split is
!!!       a minimum. Since we primarily use rotor blade meshes, the
!!!       selection of optimal direction is biased to choose the
!!!       spanwise (K-direction) over the other two directions if the
!!!       difference is less than 10%. At each level the number of
!!!       mesh blocks are twice the number at the parent level with
!!!       roughly half the number of cells, such that at the n-th
!!!       level there are 2^n mesh blocks than we initially started
!!!       with. Therefore, a physical mesh m will have LR(m)*2^n
!!!       computational blocks.
!!!
!!! 4c. Viscous direction preservation
!!!       In both steps 4a and 4b, there are additional checks
!!!       introduced to choose the next optimal direction if the
!!!       current optimum is a viscous wall normal direction. This
!!!       check helps to preserve implicitness in the viscous
!!!       direction.
!!!
!!!
!!! The manual mesh splitting requires an input file, and is
!!! implemented in the subroutine manual_split_3d. See this subroutine
!!! for more details.
!!!
!!! Additional comments:
!!!    1. The module is private. The entry point to all the
!!!        subroutines/data in this module is the subroutine
!!!        determine_proc_alloc or manual_split_3d.
!!!    2. Need to revisit the pre-bisection split logic. The
!!!       determination of wake-averaging BC might fail if the
!!!       partitioning happens along JDIR in that subroutine. This
!!!       might happen if KDIR is flagged as a viscous normal
!!!       direction.
!!!    3. The aforementioned description is the ideal case. However,
!!!       when the required number of processors are unavailable, and
!!!       still some sort of partitioning is desired, the user can
!!!       manually specify how many processors each mesh must be split
!!!       into. The split_user_defined subroutine will then just
!!!       partition each mesh along the K-direction.
!!!
   
   use constants
   use mesh_types

   implicit none

   ! Minium number of points along any given direction. 
   integer, parameter :: split_tol=16

   ! Array (size=nmesh) to store the relative ratio of number of cells
   ! in each mesh.
   integer, dimension(:), allocatable :: load_ratio

   ! nmesh - Number of participating  meshes
   ! nload - No. of near load-balanced blocks
   ! nlevels - Number of hierarchial divisions possible
   integer :: nmesh, nload, nlevels

   ! Really bad design decision! FIX ME!!!!
   ! See determine_opt_direction, pre_bisection_split,
   ! split_load_balanced.
   logical, dimension(3) :: visc_dir
   
   private
   public :: determine_proc_alloc, get_load_info, manual_split_3d, &
         get_nprocs_user_defined
   
contains
   subroutine print_split_strategy_help(nprocs)
      !! Print out an informational splitting strategy help message. 
      integer, intent(in) :: nprocs

      write(STDOUT,999) nprocs
      write(STDOUT,1000) nload
      write(STDOUT,1001) !nmesh

999   format(/,'Cannot partition domain over ',I3, ' processor(s)')       
1000  format('For load balanced domain partitioning assign',/'nprocs = ',&
            I3,'* 2**n, n=1,2,3,...')
1001  format('For per mesh splitting create a file with processor allocation')
   end subroutine print_split_strategy_help
      
   subroutine can_split_binary(nprocs,nm,splitok)
      !! Determine if the meshes can be split amongst N processors
      !! (nprocs) using a given strategy.
      !!
      !! Inputs:
      !!    nprocs - Number of processors
      !!    nm - number of meshes or near load-balanced chunks
      !!
      !! Outputs:
      !!    splitok - Returns true if splitting is possible
      !!
      !! Modifies:
      !!    nlevels - Binary split levels
      !!
      integer, intent(in) :: nprocs, nm
      logical, intent(out) :: splitok

      integer :: nfrac
      
      nlevels=-1
      nfrac=int(real(nprocs,rdp)/real(nm,rdp))
      if (abs(nprocs-nfrac*nm) > 0) then
         splitok=.false.
         return
      end if
      
      nlevels=nint(log(real(nfrac,rdp))/log(two))
      if (abs(nfrac-ishft(1,nlevels)) > 0) then
         splitok=.false.
         nlevels=-1
      else
         splitok=.true.
      end if
   end subroutine can_split_binary

   subroutine determine_mesh_load(meshes)
      !! Determine the ratio of mesh cells to determine optimal split
      !! ratio.
      !!
      !! Inputs:
      !!    meshes - array mesh_info objects
      !!
      !! Modifies:
      !!    load_ratio - Integer ratio of cells per mesh
      type(mesh_info), dimension(:), intent(in) :: meshes

      integer :: i
      integer :: mincells
      real(kind=rdp) :: mcells, tmpval

      mincells=minval(meshes(:)%ncells)
      mcells=real(mincells,rdp)

      write(STDOUT,'(/,A)') "Computing mesh cell ratios"
      do i=1,nmesh
         tmpval=(real(meshes(i)%ncells,rdp)/mcells)
         load_ratio(i)=nint(tmpval)
         if (abs(tmpval-real(load_ratio(i),rdp)) > 0) then
            write(STDOUT,1000) i, abs(tmpval-real(load_ratio(i),rdp))*100_rdp
         end if
      end do
1000  format('Mesh ID: ',I4,'; Load imbalance: ',F8.3,'%')
   end subroutine determine_mesh_load

   subroutine determine_opt_direction(nfx,nfy,nfz,opt_dir)
      !! Determine the optimum split direction based on the number of
      !! faces normal to each direction. By default we want splitting
      !! along K direction (which is spanwise for most meshes), so
      !! check if that is possible first. We retain the second and
      !! third optimum directions in the event that the optimum
      !! direction has a viscous BC, so that we can swap to the next
      !! preferred direction.
      !!
      !! Inputs:
      !!   nfx, nfy, nfz - No. of faces in J,K,L directions
      !!
      !! Outputs:
      !!   opt_dir - Array of optimum directions
      !!
      integer, intent(in) ::  nfx, nfy, nfz
      integer, dimension(3), intent(inout) :: opt_dir

      real(kind=rdp), PARAMETER :: tol_eps=0.3_rdp
      integer :: i

      ! We prefer to split along K direction, so check on that
      ! first. Even if nfx is less than nfy, prefer K-direction if the
      ! difference is less than 10%. This is necessary when JMAX and
      ! KMAX differ by 1 or so.
      if ((nfy <= nfx) .or. &
          (abs((real(nfx,rdp)/real(nfy,rdp))-one) < tol_eps)) then
         if ((nfy <= nfz) .or. &
             (abs((real(nfz,rdp)/real(nfy,rdp))-one) < tol_eps))  then
            opt_dir(1)=KDIR
            if (nfx < nfz) then
               opt_dir(2)=JDIR
               opt_dir(3)=LDIR
            else
               opt_dir(2)=LDIR
               opt_dir(3)=JDIR
            end if
         else
            opt_dir(1)=LDIR
            opt_dir(2)=KDIR
            opt_dir(3)=JDIR
         end if
      else if ((nfx <= nfz) .or. &
               (abs((real(nfz,rdp)/real(nfx,rdp))-one) < tol_eps)) then
         opt_dir(1)=JDIR
         if ((nfy <= nfz) .or. &
             (abs((real(nfz,rdp)/real(nfy,rdp))-one) < tol_eps))  then
            opt_dir(2)=KDIR
            opt_dir(3)=LDIR
         else
            opt_dir(2)=LDIR
            opt_dir(3)=KDIR
         end if
      else
         opt_dir(1)=LDIR
         opt_dir(2)=JDIR
         opt_dir(3)=KDIR
      end if

      ! Check to make sure that the optimum direction selected is not
      ! a viscous normal direction. If it is, then switch to the
      ! second best.
      do i=1,3
         if (visc_dir(opt_dir(1))) then
            opt_dir=cshift(opt_dir,1,1)
         else
            exit
         end if
      end do
      
   end subroutine determine_opt_direction
   
   subroutine bisect_block(blk1,blk2)
      !! Bisect a block (blk1) into two blocks (blk1, blk2) along the
      !! optimum direction.
      !!
      !! Inputs:
      !!    blk1 - Parent block to be split
      !!
      !! Outputs:
      !!    blk1, blk2 - Child blocks as a result of split
      !!
      type(mesh_block), intent(inout) :: blk1, blk2

      integer :: jm, km, lm, nfx, nfy, nfz, prefdir, tmp
      integer, dimension(3) :: opt_dir, dims1, dims2
      integer, dimension(3) :: start1, start2, end1, end2

      ! Set preferred direction, will be overwritten later
      opt_dir=(/ KDIR, JDIR, LDIR /)

      ! Store parent mesh dimensions for easier access
      jm=blk1%jm
      km=blk1%km
      lm=blk1%lm

      ! Set up dimensions and extents of child meshes to parent mesh
      ! values. Remember these meshes are split recursively, so the
      ! parent is not always the base mesh. This will be overwritten
      ! for direction where the mesh is split.
      dims1=(/ jm, km, lm /)
      dims2=(/ jm, km, lm /)
      start1=(/ blk1%gjs, blk1%gks, blk1%gls /)
      start2=(/ blk1%gjs, blk1%gks, blk1%gls /)
      end1=(/ blk1%gje, blk1%gke, blk1%gle /)
      end2=(/ blk1%gje, blk1%gke, blk1%gle /)

      ! We need to determine the optimum direction to split. For this
      ! we check the number of faces in each direction and choose the
      ! direction that will introduce the 1-to-1 partition with the
      ! least number of faces.
      ! Faces along J, K, and L directions
      nfx=(km-1)*(lm-1)   
      nfy=(lm-1)*(jm-1)
      nfz=(jm-1)*(km-1)

      ! Determine optimum direction based on number of faces
      ! Math: comment this for now, to prescrib splitting
      if (blk1%dom_type .eq. DOM_WAKE .or. &
          blk1%dom_type .eq. DOM_GROUNDWAKE) &
            call determine_opt_direction(nfx,nfy,nfz,opt_dir)

      ! Now that we've a preferred direction, bisect along this
      ! direction and set the new dimensions for each block. Note that
      ! here we do not yet set the overlap and ghost cell
      ! information. We will do that later.
      prefdir=opt_dir(1)
      tmp=dims1(prefdir)
      dims1(prefdir)=tmp/2+mod(tmp,2)
#ifdef TURNS_EVEN_OVERLAP      
      dims2(prefdir)=tmp/2
      end1(prefdir)=end1(prefdir)-dims2(prefdir)
      start2(prefdir)=end1(prefdir)+1
#else
      dims2(prefdir)=tmp/2+1
      end1(prefdir)=end1(prefdir)-dims2(prefdir)+1
      start2(prefdir)=end1(prefdir)
#endif
      
      if (dims1(prefdir) < split_tol) then
         call stop_execution('bisect_block', &
               "Bisect split violates minium points per block. Try &
               &with half the number of currently &
               &available processors")
      end if
      
      ! Update the mesh information for the new split blocks
      blk1%jm=dims1(1); blk1%km=dims1(2); blk1%lm=dims1(3)
      blk2%jm=dims2(1); blk2%km=dims2(2); blk2%lm=dims2(3)
      
      ! Update relative position information w.r.t. global mesh
      blk1%gjs=start1(1); blk1%gks=start1(2); blk1%gls=start1(3)
      blk2%gjs=start2(1); blk2%gks=start2(2); blk2%gls=start2(3)
      blk1%gje=end1(1); blk1%gke=end1(2); blk1%gle=end1(3)
      blk2%gje=end2(1); blk2%gke=end2(2); blk2%gle=end2(3)

      ! Set global mesh information to this mesh block
      blk2%mesh_id=blk1%mesh_id
      blk2%dom_type=blk1%dom_type
      blk2%nbc_sim=blk1%nbc_sim
      blk2%nbc_int=blk1%nbc_int
      
   end subroutine bisect_block

   subroutine recursive_bisection(mblks,nprocs,noffset)
      !! Given a mesh block perform recursive bisection up to nlevels
      !! (calculated previously) such that nprocs=2**nlevels
      !!
      !! Inputs:
      !!
      !!   nprocs - Number of processors over which binary split is
      !!            performed. Note nprocs=2**nlevels in this routine!
      !!   noffset - The offset in array mblks where the split data
      !!             must be stored
      !!
      !! Inputs & Outputs
      !!   mblks - mblks(noffset+1) contains the parent mesh that is
      !!           recursively split. After splitting,
      !!           mblks(noffset+1:noffset+nprocs) contains the split
      !!           mesh information.
      !!
      !! The logic: We have to split to 2**nlevels blocks. At a given
      !!     level n=0:(nlevels-1) the block information of the array
      !!     mblks at 1:2**nlevels:2**(nlevels-n) entries are bisected
      !!     and the entries at
      !!     (1+2**(nlevels-n+1)):2**nlevels:2**(nlevels-n+1) are
      !!     updated with the second half of the partitioned mesh. The
      !!     original values at the initial entries are overwritten by
      !!     the first half of the partitioned mesh. Remember that
      !!     this is for a single isolated block, the actual entries
      !!     need to be offset by noffset which is a function of
      !!     number of meshes and the pre-bisection partitioning.
      type(mesh_block), dimension(:), intent(inout) :: mblks
      integer, intent(in) :: nprocs, noffset

      integer :: i, j, nm, nskip

      nskip=nprocs
            print*,'nprocs = ',nprocs
      do i=0,(nlevels-1)
         nm=nskip
         nskip=nprocs/ishft(1,i+1)
         do j=1+noffset,nprocs+noffset,nm
            mblks(j)%block_id=j
            mblks(j+nskip)%block_id=j+nskip
            call bisect_block(mblks(j),mblks(j+nskip))
         end do
      end do
   end subroutine recursive_bisection

   subroutine pre_bisection_split(mblks,np,np1,noffset,jkl_dir)
      !! Perform a naive uni-directional split along the optimum
      !! direction. This is a pre-processing step before recursive
      !! bisection. This step is necessary when there are different
      !! mesh types which a difference in the number of cells per
      !! mesh. In this situation, the mesh is first split in chunks
      !! which have approximately the same number of cells as the
      !! smallest mesh such that there are 'nload' blocks with roughly
      !! equal number of computational cells. These blocks can then be
      !! recursively bisected.
      !!
      !! The preferred direction of splitting can be enforced using
      !! the optional parameter jkl_dir. If this is not provided, the
      !! optimal direction is calculated internally.
      !!
      !! Inputs:
      !!   np - number of blocks to split the mesh into
      !!   np1 - Offset between each block. This is the total number
      !!         of bisections that will be performed on this block
      !!         (2**nlevels)
      !!   noffset - Offset of the mesh block in the mblks array
      !!   jkl_dir - (Optional) Preferred direction of splitting
      !!
      !! Inputs & outputs:
      !!   mblks - Working mesh_blocks object array
      type(mesh_block), dimension(:), intent(inout) :: mblks
      integer, intent(in) :: np, np1, noffset
      integer, intent(in), optional :: jkl_dir

      integer :: jm, km, lm, opt_dir, kavg, extra, offset, n, ke, nn
      integer ::  kst, kend, ktot, i
      integer, dimension(3) :: dims1, start1, end1, dims2, start2, end2

      ! Set up dimensions and extents of the base mesh. Here we don't
      ! split recursively, so the base mesh is always the global mesh.
      jm=mblks(noffset)%jm
      km=mblks(noffset)%km
      lm=mblks(noffset)%lm
      dims1=(/ jm, km, lm /)

      if (present(jkl_dir)) then
         opt_dir=jkl_dir
         start1=(/ mblks(noffset)%gjs, mblks(noffset)%gks, mblks(noffset)%gls /)
         end1=(/ mblks(noffset)%gje, mblks(noffset)%gke, mblks(noffset)%gle /)
      else
         start1=(/ 1, 1, 1 /)
         end1=(/ jm, km, lm /)
         ! Determine the optimum direction for splitting
         opt_dir=maxloc(dims1,1)
         ! For blade-type meshes we prefer splitting in K-direction
         if ((opt_dir == JDIR) .and. &
               (jm <= 4*km) .and. &
               (mblks(noffset)%dom_type /= (DOM_WAKE .or. DOM_GROUNDWAKE))) then
            opt_dir=KDIR
         end if

         ! Swap with next optimal direction if the best direction is a
         ! viscous normal.
         do i=1,3
            if (visc_dir(opt_dir)) then
               opt_dir=maxloc(dims1,1,dims1<dims1(opt_dir))
            else
               exit
            end if
         end do
      end if

      ! Average number of points per split block
      kavg=dims1(opt_dir)/np
      if (kavg < split_tol) then
         call stop_execution('pre_bisection_split', &
               'Pre-bisection split failed tolerance test')
      end if
      extra=mod(dims1(opt_dir),np)
      offset=0

      ! Distribute points across the split blocks
      do n=1,np
         nn=noffset+np1*(n-1)
         if (n > 1) then
#ifdef TURNS_EVEN_OVERLAP
            kst=kend+1
#else            
            kst=kend!+1
#endif            
         else
            kst=1
         end if
         ke=offset+kavg
         if (n <= extra) ke=ke+1
         kend=ke
         ktot=ke-kst+1
         offset=ke

         dims2=dims1
         start2=start1
         end2=end1
         dims2(opt_dir)=ktot
         start2(opt_dir)=kst
         end2(opt_dir)=kend

         ! Update the mesh information for the new split blocks
         mblks(nn)%jm=dims2(1); mblks(nn)%km=dims2(2); mblks(nn)%lm=dims2(3)
      
         ! Update relative position information w.r.t. global mesh
         mblks(nn)%gjs=start2(1)
         mblks(nn)%gks=start2(2)
         mblks(nn)%gls=start2(3)
         mblks(nn)%gje=end2(1); mblks(nn)%gke=end2(2); mblks(nn)%gle=end2(3)

         mblks(nn)%mesh_id=mblks(noffset)%mesh_id
         mblks(nn)%dom_type=mblks(noffset)%dom_type
         mblks(nn)%nbc_sim=mblks(noffset)%nbc_sim
         mblks(nn)%nbc_int=mblks(noffset)%nbc_int
         mblks(nn)%block_id=nn
      end do
   end subroutine pre_bisection_split

   subroutine set_block_info_global(nmesh,meshes,mblks)
      integer, intent(in) :: nmesh
      type(mesh_info), dimension(:), intent(inout) :: meshes
      type(mesh_block), dimension(:), allocatable, intent(inout) :: mblks

      integer :: m, b, bid, nblocks, ii
      integer, dimension(:), allocatable :: blkidx
      character(len=10) :: intstr

      nblocks=size(mblks,1)
      ! We need to track the last index per mesh that was updated, so
      ! have an array of indices that will be updated to the next
      ! unfilled index in the mesh_block.
      allocate(blkidx(nmesh))
      blkidx=1 ! Initialize all to 1 at start
      
      do m=1,nmesh
         nullify(meshes(m)%block_id,meshes(m)%start,meshes(m)%end)
         allocate(meshes(m)%block_id(meshes(m)%nblocks))
         allocate(meshes(m)%start(3,meshes(m)%nblocks))
         allocate(meshes(m)%end(3,meshes(m)%nblocks))
      end do

      do b=1,nblocks
         m=mblks(b)%mesh_id
         bid=mblks(b)%block_id
         ii=blkidx(m)

         meshes(m)%block_id(ii)=bid
         meshes(m)%start(JDIR,ii)=mblks(b)%gjs
         meshes(m)%start(KDIR,ii)=mblks(b)%gks
         meshes(m)%start(LDIR,ii)=mblks(b)%gls
         meshes(m)%end(JDIR,ii)=mblks(b)%gje
         meshes(m)%end(KDIR,ii)=mblks(b)%gke
         meshes(m)%end(LDIR,ii)=mblks(b)%gle

         if (blkidx(m) == 1) then
            mblks(b)%is_master=.true.
            meshes(m)%master_id=bid
         else
            mblks(b)%is_master=.false.
         end if
         ! Set up the block name
         write(intstr,'(I10)') ii
         mblks(b)%block_name=trim(adjustl(meshes(m)%mesh_name))//"_"//trim(adjustl(intstr))
         
         blkidx(m)=blkidx(m)+1
      end do

      deallocate(blkidx)
   end subroutine set_block_info_global
   
   subroutine split_load_balanced(nmesh,meshes,nprocs,mblks)
      !! Recursively bisect the meshes across the given number of
      !! processors. Also perform pre-splitting if necessary to
      !! balance loads between disparate meshes.
      !!
      !! Inputs:
      !!    nmesh - number of meshes
      !!    nprocs - Total number of processors available
      !!
      !! Inputs & outputs:
      !!    meshes - The global mesh object array, this gets updated
      !!             with the number of subblocks it has been split
      !!             into.
      !!    mblks - The computational block information. Almost all
      !!            information in this object is set by this
      !!            subroutine (and the ones it calls).
      !!
      integer, intent(in) :: nmesh, nprocs
      type(mesh_info), dimension(:), intent(inout) :: meshes
      type(mesh_block), dimension(:), allocatable, intent(inout) :: mblks

      integer :: m, np1, n, npused, nppm

      character(len=32) :: nam

      allocate(mblks(nprocs))

      ! Index of the last processor that has a block assigned to it. 
      npused=0

      ! Split over meshes
      do n=1,nmesh
         ! No. of procs per mesh
         nppm=(load_ratio(n)*nprocs)/nload
         meshes(n)%nblocks=nppm

         ! No. of procs per load balanced block of this mesh
         np1=nppm/load_ratio(n)

         ! Mesh name for printouts
         nam=(meshes(n)%mesh_name)

         write(STDOUT,1000) trim(nam), nppm

         ! Set viscous direction information
         visc_dir=meshes(n)%visc_dir
         
         ! Initialize first block with the global mesh information
         call copy_mesh_to_block(meshes(n),mblks(npused+1))
         ! Set starting block ID for this partitioning
         mblks(npused+1)%block_id=npused+1

         ! Perform pre-bisection split if necessary
         if (load_ratio(n) > 1) then
            write(STDOUT,1001) trim(nam), load_ratio(n)
            call pre_bisection_split(mblks,load_ratio(n),np1,npused+1)
         end if

         ! Perform recursive bisection of meshes
         if (nlevels > 0) then
            write(STDOUT,1002) trim(nam), nlevels
         end if
         do m=1,load_ratio(n) 
            call recursive_bisection(mblks,np1,npused+np1*(m-1))
         end do
         npused=npused+nppm
         write(STDOUT,1003) trim(nam)
      end do

      ! Update global mesh object with subblock information
      call set_block_info_global(nmesh,meshes,mblks)
      
      
1000  format(/,'Attempting partition of mesh ',A, ' into ',I3,' parts')
1001  format('|-> Pre-bisection split of mesh ',A,' into ', I3,' parts')      
1002  format('|-> Recursive bisection of mesh ',A,' upto ', I3,' levels')      
1003  format('|-> Mesh ',A,' partitioning completed successfully')
   end subroutine split_load_balanced

   subroutine process_block_bcs(meshes,mblks)
      type(mesh_info), dimension(:), intent(in) :: meshes
      type(mesh_block), dimension(:), intent(inout) :: mblks

      integer :: m, n, b, f

      ! Total number of meshes and subblocks 
      integer :: nmesh, nblocks
      ! Maximum value of nbc_sim and nbc_int for the global meshes.
      integer :: nsbc, nifbc

      ! Track number of simple BCs and interface BCs per block
      integer :: nsimbcs, nifbcs
      ! Track corresponding BC IDs for the blocks
      integer, dimension(:), allocatable :: sbcid, ifbcid

      ! To process the BC planes in a generic way, a normal direction
      ! index Zdir is always set to the index normal to the plane. X-
      ! and Y-dir are set to the other two indices. These indices can
      ! be then used to correctly choose the right values from the
      ! start, end, bcst, and bcend arrays. See function block_has_bc
      ! for actual implementation.
      integer :: xdir, ydir, zdir

      ! start, end - Array storing the J, K, and L indices of the
      !              starting index (in the global mesh) and the
      !              ending index respectively.
      ! bcst, bcend - Array storing the J, K, and L limits for the BC
      !               in the global mesh plane.
      ! gldims - Max dimensions of the global base mesh
      integer, dimension(3) :: start, end, bcst, bcend, gldims

      ! mid - Global mesh id
      ! blkid - ID of the current block
      ! bcdir - Direction of the BC under consideration
      ! bcid - ID of the BC under consideration
      integer :: mid, blkid, bcdir, bcid

      ! Temporary data to identify new interfaces created by domain
      ! partitioning.
      ! nnbc - Total no. of new BC faces created per subblock
      ! donid - Index of the donor block
      ! nnbdir - Direction index of this BC patch
      integer :: nnbc, donid
      integer, dimension(6) :: nnbdir

      type(simple_bc), pointer :: sbc
      type(interface_bc), pointer :: ifbc

      nmesh=size(meshes,1)
      nblocks=size(mblks,1)
      nullify(sbc,ifbc)

      ! Calculate maximum number of BCs that might be allocated
      call calc_max_bc_nums
      allocate(sbcid(nsbc),ifbcid(nifbc))

      ! Process each block and determine if they need to have base
      ! mesh's BC information.
      do n=1,nblocks
         nsimbcs=0
         nifbcs=0
         sbcid=0
         ifbcid=0
         nnbc=0
         donid=0
         nnbdir=0
         
         mid=mblks(n)%mesh_id
         gldims=(/ meshes(mid)%jm, meshes(mid)%km, meshes(mid)%lm /)
         blkid=mblks(n)%block_id
         start=(/ mblks(n)%gjs, mblks(n)%gks, mblks(n)%gls /)
         end=(/ mblks(n)%gje, mblks(n)%gke, mblks(n)%gle /)

         ! Check and process simple BC information
         call preprocess_simple_bcs
         call process_simple_bcs
         
         ! Check and process interface BC information
         call preprocess_iface_bcs
         ! Determine the internal faces created because of domain
         ! partitioning.
         call identify_newface_bcs
         call process_iface_bcs

      end do

      ! do n=1,nblocks
      !    print '(4I4)', n, mblks(n)%jm, mblks(n)%km, mblks(n)%lm
      !    do b=1,mblks(n)%nbc_sim
      !       sbc=>mblks(n)%bcinfo(b)
      !       print '(6I4)', sbc%jbcs, &
      !             sbc%jbce, &
      !             sbc%kbcs, &
      !             sbc%kbce, &
      !             sbc%lbcs, &
      !             sbc%lbce
                  
      !    end do
      !    print *
      ! end do

      ! do n=1,nblocks
      !    print '(4I4)', n, mblks(n)%jm, mblks(n)%km, mblks(n)%lm
      !    do b=1,mblks(n)%nbc_int
      !       ifbc=>mblks(n)%ifaceinfo(b)
      !       print *, ifbc%bc_name, ifbc%donorid, ifbc%ibdir
      !       print *, ifbc%jbcs, ifbc%jbce, &
      !             ifbc%kbcs, ifbc%kbce, ifbc%lbcs, ifbc%lbce
      !       print *, ifbc%djs, ifbc%dje, &
      !             ifbc%dks, ifbc%dke, ifbc%dls, ifbc%dle
      !    end do
      !    print *
      ! end do
      
      deallocate(sbcid,ifbcid)
   contains
      subroutine calc_max_bc_nums
         !! Helper routine to calculate the max simple BC and
         !! interface BC amongst the given meshes.

         nsbc=0
         nifbc=0
         do m=1,nmesh
            if (meshes(m)%nbc_sim > nsbc) nsbc=meshes(m)%nbc_sim
            if (meshes(m)%nbc_int > nifbc) nifbc=meshes(m)%nbc_int
         end do
      end subroutine calc_max_bc_nums

      function block_has_bc result(bcok)
         !! Determine if the block face is coincident with the base
         !! mesh's face where the BC was defined originally. If yes,
         !! then check if there is an overlap of the planes. Return
         !! true, if both conditions hold true, or false otherwise.
         logical :: bcok
         integer :: fid1, fid2

         bcok=.false.

         ! Set up indices according to the direction of BC
         select case (bcdir)
         case (1,-1)
            zdir=JDIR; xdir=KDIR; ydir=LDIR
         case (2,-2)
            zdir=KDIR; xdir=KDIR; ydir=JDIR
         case (3,-3)
            zdir=LDIR; xdir=JDIR; ydir=KDIR
         case default
            call stop_execution('process_block_bcs', &
                  "Invalid BC direction supplied")
         end select

         ! Set up index of the face plane
         if (bcdir > 0) then
            fid1=start(zdir)
            fid2=bcst(zdir)
         else
            fid1=end(zdir)
            fid2=bcend(zdir)
         end if

         if (fid1 == fid2) then
            ! We are on the same face, but do the planes overlap?
            if ((start(xdir) <= bcend(xdir)) .and. &
                (bcst(xdir)  <= end(xdir)) .and. &
                (start(ydir) <= bcend(ydir)) .and. &
                (bcst(ydir)  <= end(ydir))) then
               bcok=.true.
            end if
         end if
      end function block_has_bc

      subroutine preprocess_simple_bcs
         !! Simple helper routine to loop through all the simple BCs
         !! of the base mesh and check if they need to be transferred
         !! to the subblocks. If yes, then the number of subblock BCs
         !! is updated, and the BC ID from the base mesh is stored for
         !! future processing. At the end of execution, we know the
         !! total number of simple BCs that need to be transferred to
         !! this subblock.
         !!
         do b=1,meshes(mid)%nbc_sim
            sbc=>meshes(mid)%bcinfo(b)
            bcdir=sbc%ibdir
            bcst=(/ sbc%jbcs, sbc%kbcs, sbc%lbcs /)
            bcend=(/sbc%jbce, sbc%kbce, sbc%lbce /)
            if (block_has_bc()) then
               nsimbcs=nsimbcs+1
               sbcid(nsimbcs)=b
            end if
         end do
      end subroutine preprocess_simple_bcs

      subroutine process_simple_bcs
         !! Setup simple BC information for mesh blocks.
         !!
         !! The preprocess step has determined the number of simple
         !! BCs in a sub-block and their corresponding IDs in the base
         !! mesh. So make another pass through those BCs and set those
         !! up for the mesh blocks. Remember, we need to reset the
         !! limits such that they correspond to the local mesh
         !! indices.

         mblks(n)%nbc_sim=nsimbcs
         nullify(mblks(n)%bcinfo)
         allocate(mblks(n)%bcinfo(nsimbcs))
         do b=1,nsimbcs
            sbc=>mblks(n)%bcinfo(b)
            bcid=sbcid(b)
            ! Basic BC information
            sbc%bc_name=meshes(mid)%bcinfo(bcid)%bc_name
            sbc%ibtyp=meshes(mid)%bcinfo(bcid)%ibtyp
            sbc%ibdir=meshes(mid)%bcinfo(bcid)%ibdir

            ! Shorter version for easier writing
            bcst=(/ meshes(mid)%bcinfo(bcid)%jbcs, &
                  meshes(mid)%bcinfo(bcid)%kbcs, &
                  meshes(mid)%bcinfo(bcid)%lbcs /)
            bcend= (/ meshes(mid)%bcinfo(bcid)%jbce, &
                  meshes(mid)%bcinfo(bcid)%kbce, &
                  meshes(mid)%bcinfo(bcid)%lbce /)
            ! Set up limits 
            sbc%jbcs=max(start(JDIR),bcst(JDIR)) -start(JDIR)+1
            sbc%jbce=min(end(JDIR),bcend(JDIR)) -start(JDIR)+1
            sbc%kbcs=max(start(KDIR),bcst(KDIR)) -start(KDIR)+1
            sbc%kbce=min(end(KDIR),bcend(KDIR))-start(KDIR)+1
            sbc%lbcs=max(start(LDIR),bcst(LDIR))-start(LDIR)+1
            sbc%lbce=min(end(LDIR),bcend(LDIR))-start(LDIR)+1
         end do
      end subroutine process_simple_bcs

      subroutine preprocess_iface_bcs
         !! Helper routine for interface BCs, similar to
         !! preprocess_simple_bcs. However, this one does some extra work
         !! to treat wake averaging BCs correctly.
         !!
         integer :: adir, tmp
         
         do b=1,meshes(mid)%nbc_int
            ifbc=>meshes(mid)%ifaceinfo(b)
            bcdir=ifbc%ibdir
            bcst=(/ ifbc%jbcs, ifbc%kbcs, ifbc%lbcs /)
            bcend=(/ifbc%jbce, ifbc%kbce, ifbc%lbce /)
            if (block_has_bc()) then
               nifbcs=nifbcs+1
               ifbcid(nifbcs)=b

            ! Even if the check failed, we still need to check if this
            ! block can participate in wake averaging BC. We didn't
            ! write them out as two interfaces, so need to introduce
            ! this hack.
            else if (ifbc%ibtyp == BC_WAKEAVG) then
               ! This time check with the donor face information that
               ! was automatically generated earlier. NOTE: needs a
               ! fix for JLE point.
               bcdir=ifbc%didir
               bcst=(/ ifbc%djs, ifbc%dks, ifbc%dls /)
               bcend= (/ifbc%dje, ifbc%dke, ifbc%dle /)
               ! Swap out JBCS and JBCE because they run -ve JDIR
               tmp=bcst(JDIR)
               bcst(JDIR)=bcend(JDIR)
               bcend(JDIR)=tmp
               ! Is the donor face part of this block
               if (block_has_bc()) then
                  nifbcs=nifbcs+1
                  ! Set negative index to remind us to use donor
                  ! blocks.
                  ifbcid(nifbcs)=-b 
               end if
            end if
         end do
      end subroutine preprocess_iface_bcs

      subroutine identify_newface_bcs
         logical :: endface
         integer :: f1

         ! Process each face for a particular subblock
         do f=1,6
            ! Check whether it is coincident with edge faces of the
            ! global mesh.
            endface=.false.
            select case (f)
            case (1:3)
               ! Process the left faces 
               if (start(f) == 1) endface=.true.
               f1=f
            case (4:6)
               ! Process the right faces
               f1=3-f
               if (end(f-3) == gldims(f-3)) endface=.true.
            end select

            if (.not. endface) then
               nnbc=nnbc+1
               nnbdir(nnbc)=f1
            end if
         end do
      end subroutine identify_newface_bcs

      subroutine process_iface_bcs
         integer :: ii, dir, nfcur, mdon_id
         type(interface_bc), pointer :: ifbc0
         integer, dimension(3) :: bs1,be1,bs2,be2
         
         ! For a given block, the no. of interface BCs is the sum of
         ! edge BCs inherited from 0-level parent, as well as the
         ! internal faces created by mesh splitting.
         mblks(n)%nbc_int=nnbc+nifbcs
         nullify(mblks(n)%ifaceinfo)
         allocate(mblks(n)%ifaceinfo(mblks(n)%nbc_int))
         nfcur=0

         end_faces: do b=1,nifbcs
            nfcur=nfcur+1
            bcid=ifbcid(b)
            ifbc=>mblks(n)%ifaceinfo(nfcur)
            ifbc0=>meshes(mid)%ifaceinfo(abs(bcid))
            ifbc%bc_name=ifbc0%bc_name
            ifbc%ibtyp=ifbc0%ibtyp
            
            select case (ifbc0%ibtyp)
            case (BC_WAKEAVG) ! This  needs special treatment
               ! NOTE: The implementation makes a VERY SIMPLIFYING
               ! ASSUMPTION that the mesh is split up only along the
               ! K-direction. This is how it will be for the
               ! blade/slat meshes, but must be aware of this.
               ifbc%ibdir=ifbc0%ibdir
               ifbc%jbcs=ifbc0%jbcs
               ifbc%jbce=ifbc0%jbce
               ifbc%kbcs=max(start(KDIR),ifbc0%kbcs)-start(KDIR)+1
               ifbc%kbce=min(end(KDIR),ifbc0%kbce)-start(KDIR)+1
               ifbc%lbcs=ifbc0%lbcs
               ifbc%lbce=ifbc0%lbce
               ifbc%donor_mesh=mblks(n)%block_name
               ifbc%donorid=n
               ifbc%didir=ifbc0%didir
               ifbc%djs=ifbc0%djs
               ifbc%dje=ifbc0%dje
               ifbc%dks=max(start(KDIR),ifbc0%dks)-start(KDIR)+1
               ifbc%dke=min(ifbc0%dke,end(KDIR))-start(KDIR)+1
               ifbc%dls=ifbc0%dls
               ifbc%dle=ifbc0%dle
            case default
               donid=0
               bcdir=ifbc0%ibdir
               zdir=abs(bcdir)
               mdon_id=ifbc0%donorid

               select case (zdir)
               case (JDIR)
                  xdir=KDIR; ydir=LDIR
               case (KDIR)
                  xdir=LDIR; ydir=JDIR
               case (LDIR)
                  xdir=JDIR; ydir=KDIR
               end select
               ! Shorter version for easier writing
               bcst=(/ ifbc0%jbcs, ifbc0%kbcs, ifbc0%lbcs /)
               bcend=(/ ifbc0%jbce, ifbc0%kbce, ifbc0%lbce /)
               ! Set up limits 
               bs1(JDIR)=max(start(JDIR),bcst(JDIR))
               be1(JDIR)=min(end(JDIR),bcend(JDIR))
               bs1(KDIR)=max(start(KDIR),bcst(KDIR))
               be1(KDIR)=min(end(KDIR),bcend(KDIR))
               bs1(LDIR)=max(start(LDIR),bcst(LDIR))
               be1(LDIR)=min(end(LDIR),bcend(LDIR))
               bs2(1)=ifbc0%djs+sign(1,ifbc0%dje-ifbc0%djs)* &
                     (bs1(1)-ifbc0%jbcs)
               bs2(2)=ifbc0%dks+sign(1,ifbc0%dke-ifbc0%dks)* &
                     (bs1(2)-ifbc0%kbcs)
               bs2(3)=ifbc0%dls+sign(1,ifbc0%dle-ifbc0%dls)* &
                     (bs1(3)-ifbc0%lbcs)
               be2(1)=ifbc0%djs+sign(1,ifbc0%dje-ifbc0%djs)* &
                     (be1(1)-ifbc0%jbcs)
               be2(2)=ifbc0%dks+sign(1,ifbc0%dke-ifbc0%dks)*&
                     (be1(2)-ifbc0%kbcs)
               be2(3)=ifbc0%dls+sign(1,ifbc0%dle-ifbc0%dls)*&
                        (be1(3)-ifbc0%lbcs)

               do ii=1,meshes(mdon_id)%nblocks
                  select case (ifbc0%didir)
                  case (1:3)
                     if ((bs2(zdir) == meshes(mdon_id)%start(zdir,ii)) .and. &
                         (bs2(xdir) == meshes(mdon_id)%start(xdir,ii)) .and. &
                         (be2(xdir) == meshes(mdon_id)%end(xdir,ii)) .and. &
                         (bs2(ydir) == meshes(mdon_id)%start(ydir,ii)) .and. &
                         (be2(ydir) == meshes(mdon_id)%end(ydir,ii))) then
                        donid=ii
                     end if
                  case (-3:-1)
                     if ((be2(zdir) == meshes(mdon_id)%end(zdir,ii)) .and. &
                         (bs2(xdir) == meshes(mdon_id)%start(xdir,ii)) .and. &
                         (be2(xdir) == meshes(mdon_id)%end(xdir,ii)) .and. &
                         (bs2(ydir) == meshes(mdon_id)%start(ydir,ii)) .and. &
                         (be2(ydir) == meshes(mdon_id)%end(ydir,ii))) then
                        donid=ii
                     end if
                     if (donid > 0) exit
                  end select
               end do
               if (donid <= 0) then
                  call stop_execution('process_iface_bcs', &
                        "Cannot find donor for interface face, possible error&
                        &in domain partitioning logic")
               end if

               ifbc%ibdir=bcdir
               ifbc%jbcs=bs1(JDIR)-start(JDIR)+1
               ifbc%kbcs=bs1(KDIR)-start(KDIR)+1
               ifbc%lbcs=bs1(LDIR)-start(LDIR)+1
               ifbc%jbce=be1(JDIR)-start(JDIR)+1
               ifbc%kbce=be1(KDIR)-start(KDIR)+1
               ifbc%lbce=be1(LDIR)-start(LDIR)+1
               ifbc%donor_mesh=mblks(meshes(mdon_id)%block_id(donid))%block_name
               ifbc%donorid=meshes(mdon_id)%block_id(donid)
               ifbc%didir=ifbc0%didir
               ifbc%djs=bs2(JDIR)-meshes(mdon_id)%start(JDIR,donid)+1
               ifbc%dje=be2(JDIR)-meshes(mdon_id)%start(JDIR,donid)+1
               ifbc%dks=bs2(KDIR)-meshes(mdon_id)%start(KDIR,donid)+1
               ifbc%dke=be2(KDIR)-meshes(mdon_id)%start(KDIR,donid)+1
               ifbc%dls=bs2(LDIR)-meshes(mdon_id)%start(LDIR,donid)+1
               ifbc%dle=be2(LDIR)-meshes(mdon_id)%start(LDIR,donid)+1
            end select
         end do end_faces
         
         ! Now process the new faces that were created by partitioning
         new_faces: do b=1,nnbc
            donid=0
            dir=nnbdir(b)
            zdir=abs(dir)

            select case (zdir)
            case (JDIR)
               xdir=KDIR; ydir=LDIR
            case (KDIR)
               xdir=LDIR; ydir=JDIR
            case (LDIR)
               xdir=JDIR; ydir=KDIR
            end select
            ! Set limits of this BC plane
            bcst(xdir)=1
            bcst(ydir)=1
            bcend(xdir)=end(xdir)-start(xdir)+1
            bcend(ydir)=end(ydir)-start(ydir)+1
            if (dir > 0) then
               bcst(zdir)=1
               bcend(zdir)=1
            else
               bcst(zdir)=end(zdir)-start(zdir)+1
               bcend(zdir)=bcst(zdir)
            end if

            ! Find the matching 1-to-1 face amongst the sibling blocks
            siblings: do ii=1,meshes(mid)%nblocks
               ! Am I my own sibling?
               if (meshes(mid)%block_id(ii) == n) cycle
               facematch: select case (dir)
               case (1:3)
#ifdef TURNS_EVEN_OVERLAP
                  if (((start(zdir)-1) == meshes(mid)%end(zdir,ii)) .and. &
#else
                  if (((start(zdir)) == meshes(mid)%end(zdir,ii)) .and. &
#endif
                      (start(xdir) == meshes(mid)%start(xdir,ii)) .and. &
                      (end(xdir) == meshes(mid)%end(xdir,ii)) .and. &
                      (start(ydir) == meshes(mid)%start(ydir,ii)) .and. &
                      (end(ydir) == meshes(mid)%end(ydir,ii))) then
                     donid=ii
                  end if
               case (-3:-1)
#ifdef TURNS_EVEN_OVERLAP
                  if (((end(zdir+1)) == meshes(mid)%start(zdir,ii)) .and. &
#else
                  if (((end(zdir)) == meshes(mid)%start(zdir,ii)) .and. &
#endif
                      (start(xdir) == meshes(mid)%start(xdir,ii)) .and. &
                      (end(xdir) == meshes(mid)%end(xdir,ii)) .and. &
                      (start(ydir) == meshes(mid)%start(ydir,ii)) .and. &
                      (end(ydir) == meshes(mid)%end(ydir,ii))) then
                     donid=ii
                  end if
               end select facematch
               if (donid > 0) exit
            end do siblings

            if (donid <= 0) then
               call stop_execution('process_iface_bcs', &
                     "Cannot find donor for internal face, possible error&
                     &in domain partitioning logic")
            end if

            ! Set up BC info in the blocks object
            nfcur=nfcur+1
            ifbc=>mblks(n)%ifaceinfo(nfcur)
            ifbc%bc_name="INTERNAL"
            ifbc%ibtyp=BC_INTERNAL
            ifbc%ibdir=dir
            ifbc%jbcs=bcst(JDIR)
            ifbc%kbcs=bcst(KDIR)
            ifbc%lbcs=bcst(LDIR)
            ifbc%jbce=bcend(JDIR)
            ifbc%kbce=bcend(KDIR)
            ifbc%lbce=bcend(LDIR)
            ifbc%donor_mesh=mblks(meshes(mid)%block_id(donid))%block_name
            ifbc%donorid=meshes(mid)%block_id(donid)
            ifbc%didir=-dir
            ! Fix Donor plane information
            if (dir > 0) then
               bcst(zdir)=meshes(mid)%end(zdir,donid) &
                     -meshes(mid)%start(zdir,donid)+1
               bcend(zdir)=bcst(zdir)
            else
               bcst(zdir)=1
               bcend(zdir)=1
            end if
            ifbc%djs=bcst(JDIR)
            ifbc%dks=bcst(KDIR)
            ifbc%dls=bcst(LDIR)
            ifbc%dje=bcend(JDIR)
            ifbc%dke=bcend(KDIR)
            ifbc%dle=bcend(LDIR)
         end do new_faces

      end subroutine process_iface_bcs

   end subroutine process_block_bcs

   subroutine read_proc_alloc_info(meshes,nprocs,proc_alloc)
      !! If load balanced split is not possible, see if we can load up
      !! the splitting logic from user inputs.
      !!
      !! Inputs:
      !!    nprocs - Number of processors available
      !!    meshes - Information about the computational domain
      !!    proc_alloc - Filename containing the processor allocation
      !!
      integer, intent(in) :: nprocs
      type(mesh_info), dimension(:), intent(in) :: meshes
      character(len=*) :: proc_alloc

      integer :: un, i, nmesh, npa, mid
      character(len=32) :: gnam

      ! Reset load_ratio array which was filled with load balanced
      ! ratio information
      load_ratio=0
      nmesh=size(meshes,1)
      un=open_file(proc_alloc,form='formatted',status='old')
      do i=1,nmesh
         read(un,*) gnam, npa
         mid=get_global_mesh_id(meshes,gnam)
         load_ratio(mid)=npa
      end do
      close(un)
      nload=sum(load_ratio)

      print*,'nload, nprocs = ',nload,nprocs
      if (nload /= nprocs) call stop_execution('read_proc_alloc_info',&
            "Total number of procs in "//trim(adjustl(proc_alloc))//" is &
            &not equal to number of processors available")
   end subroutine read_proc_alloc_info

   subroutine split_user_defined(nmesh,meshes,nprocs,mblks)
      !! Split meshes based on user defined partitioning scheme. This
      !! doesn't guarantee near load-balanced execution, but is
      !! provided so that the user can run the code under situations
      !! when the optimal number of processors are unavailable. This
      !! subroutine is invoked only if the load balanced split is not
      !! possible with the given number of processors.
      !!
      !! Inputs:
      !!    nmesh - number of meshes
      !!    nprocs - Total number of processors available
      !!
      !! Inputs & outputs:
      !!    meshes - The global mesh object array, this gets updated
      !!             with the number of subblocks it has been split
      !!             into.
      !!    mblks - The computational block information. Almost all
      !!            information in this object is set by this
      !!            subroutine (and the ones it calls).
      !!
      integer, intent(in) :: nmesh, nprocs
      type(mesh_info), dimension(:), intent(inout) :: meshes
      type(mesh_block), dimension(:), allocatable, intent(inout) :: mblks

      integer :: n, npused, nppm
      character(len=32) :: nam

      allocate(mblks(nprocs))
      ! Index of the last processor that has a block assigned to it.
      npused=0

      do n=1,nmesh
         nppm=load_ratio(n)
         nam=(meshes(n)%mesh_name)

         if (nppm == 0) call stop_execution('split_user_defined',&
               trim(nam)//" has 0 processors allocated to it.")
         meshes(n)%nblocks=nppm

         write(STDOUT,1000) trim(nam), nppm

         ! Set viscous direction information
         visc_dir=meshes(n)%visc_dir

         ! Initialize first block with the global mesh information
         call copy_mesh_to_block(meshes(n),mblks(npused+1))
         ! Set starting block ID for this partitioning
         mblks(npused+1)%block_id=npused+1

         ! User defined splitting into n parts
         if (nppm > 1) call pre_bisection_split(mblks,nppm,1,npused+1)
         npused=npused+nppm
         write(STDOUT,1003) trim(nam)
      end do

      ! Update global mesh object with subblock information
      call set_block_info_global(nmesh,meshes,mblks)
      
1000  format(/,'Attempting partition of mesh ',A, ' into ',I3,' parts')
1003  format('|-> Mesh ',A,' partitioning completed successfully')
    end subroutine split_user_defined
   
   subroutine determine_proc_alloc(meshes,nprocs,mblks,palloc_file)
      !! Determine the domain partitioning strategy based on the mesh
      !! information and the number of processors available.
      !!
      !! Inputs:
      !!    nprocs - Number of processors available
      !!
      !! Modifies:
      !!    nload - Total number of near load-balanced chunks
      !!
      !! Inputs & Outputs:
      !!    meshes - gets updated with the subblock information
      !!
      !! Outputs:
      !!    mblks - gets populated with the computational block info.
      !!
      type(mesh_info), dimension(:), intent(inout) :: meshes
      integer, intent(in) :: nprocs
      type(mesh_block), dimension(:), allocatable, intent(inout) :: mblks
      character(len=*), intent(in) :: palloc_file

      logical :: splitok

      nmesh=size(meshes,1)
      allocate(load_ratio(nmesh))

      ! Determine the ratio of cells in participating meshes
      call determine_mesh_load(meshes)
      nload=sum(load_ratio)

      ! Can we achieve a near load-balanced split? 
      call can_split_binary(nprocs,nload,splitok)
      if (splitok) then
         write(STDOUT,'(/,A)') "Near load balanced split possible; &
               &attempting partition"
         call split_load_balanced(nmesh,meshes,nprocs,mblks)
      else
         inquire(file=palloc_file,exist=splitok)
         if (splitok) then
            write(STDOUT,'(/,A,/,A)') "Load balanced split not possible...",&
                  "Resorting to per mesh splitting"
            call read_proc_alloc_info(meshes,nprocs,palloc_file)
            call split_user_defined(nmesh,meshes,nprocs,mblks)
         else
            ! None of the strategies work, print suggestions for
            ! possible choices of number of processors and quit.
            call print_split_strategy_help(nprocs)
            call stop_execution('determine_proc_alloc', &
                  "Failed domain partitioning. Run ot_numprocs to &
                  &determine preferred number of processors.")
         end if
      end if

      call process_block_bcs(meshes,mblks)
   end subroutine determine_proc_alloc

   function get_load_info(meshes) result(nload)
      type(mesh_info), dimension(:), intent(in) :: meshes
      integer :: nload

      nmesh=size(meshes,1)
      allocate(load_ratio(nmesh))

      ! Determine the ratio of cells in participating meshes
      call determine_mesh_load(meshes)
      nload=sum(load_ratio)
   end function get_load_info

   subroutine manual_split_3d(meshes,nprocs,mblks,palloc_file)
      !! Split the participating meshes based on user defined
      !! strategy.
      !!
      !! Unlike the determine_proc_alloc strategies, this subroutine
      !! does not try to guess the splitting strategy; instead it
      !! requires input from the user regarding how each mesh must be
      !! split along the three directions. This subroutine is provided
      !! for meshes with disparate sizes, or special splitting
      !! requirements to address certain constraints in mesh topology.
      !!
      !! Inputs:
      !!    nprocs - Number of processors available
      !!    palloc_file - File which contains the user defined
      !!                  splitting strategy
      !!
      !! Outputs:
      !!    mblks - gets populated with the compuational block info.
      !!
      !! Inputs & Outputs:
      !!    meshes - gets updated with the sub-block information
      !!
      type(mesh_info), dimension(:), intent(inout) :: meshes
      integer, intent(in) :: nprocs
      type(mesh_block), dimension(:), allocatable, intent(inout) :: mblks
      character(len=*), intent(in) :: palloc_file

      ! Number of participating meshes
      integer :: nmesh
      ! Splitting information for each direction for each mesh
      integer, dimension(:,:), allocatable :: split_stats

      nmesh=size(meshes,1)
      allocate(split_stats(3,nmesh))
      call parse_proc_allocations
      call perform_manual_split
      deallocate(split_stats)
   contains
      subroutine parse_proc_allocations
         !! Helper subroutine to read the user defined splitting input
         !! file.
         integer :: un, i, mid
         integer :: sj, sk, sl, counter
         character(len=32) :: gnam

         ! Initialize the running counter for allocated processors and
         ! the split statistics to allocate 1 processor per mesh
         counter=0
         split_stats=1
         
         un=open_file(palloc_file,form='formatted',status='old')
         do i=1,nmesh
            read(un,*) gnam, sj, sk, sl
            mid=get_global_mesh_id(meshes,gnam)
            split_stats(JDIR,mid)=sj
            split_stats(KDIR,mid)=sk
            split_stats(LDIR,mid)=sl

            ! Determine the number of processors allocated to this
            ! mesh, and add it to a running counter.
            counter=counter+sj*sk*sl
         end do
         close(un)

         if (counter /= nprocs) call stop_execution('manual_split_3d', &
               "Total number of processors in "//trim(adjustl(palloc_file))//&
               " is not equal to the number of processors available")
      end subroutine parse_proc_allocations

      subroutine perform_manual_split
         !! Split meshes based on the user-defined strategy.
         !!
         !! The splitting strategy in this subroutine is implemented
         !! such that all the blocks that are split along the J
         !! direction will be in contiguous processors. The J-split
         !! blocks are ordered based on the number of K-dir
         !! split. Finally, these blocks are arranged along the
         !! L-dir. This is achieved by splitting first in the
         !! L-direction, and inserting them the mesh blocks array with
         !! pre-caclulated offsets (NJDIR*NKDIR). Then each of these
         !! blocks are split along K-dir and inserted with
         !! offset=NJDIR. Finally the J-dir split is performed.
         !!
         integer :: n, npused, nppm
         integer :: sj, sk, sl
         integer :: j, k, l, np1, np2
         character(len=32) :: nam

         allocate(mblks(nprocs))
         ! Indes of the last processor that has a block assigned to it.
         npused=0

         do n=1,nmesh
            nam=meshes(n)%mesh_name
            sj=split_stats(JDIR,n)
            sk=split_stats(KDIR,n)
            sl=split_stats(LDIR,n)
            nppm=sj*sk*sl
            meshes(n)%nblocks=nppm

            if (nppm == 0) call stop_execution('manual_split_3d', &
                  trim(nam)//" has 0 processors allocated to it.")

            write(STDOUT,1000) trim(nam), nppm
            ! Set viscous direction information
            visc_dir=meshes(n)%visc_dir
            ! Initialize first block with the global mesh information
            call copy_mesh_to_block(meshes(n),mblks(npused+1))
            ! Set starting block ID for this partitioning
            mblks(npused+1)%block_id=npused+1

            if (sl > 1) &
                  call pre_bisection_split(mblks,sl,sj*sk,npused+1,LDIR)
            do l=1,sl
               np1=sj*sk*(l-1)
               if (sk > 1) &
                     call pre_bisection_split(mblks,sk,sj,npused+1+np1,KDIR)
               do k=1,sk
                  np2=np1+sj*(k-1)
                  if (sj > 1) &
                        call pre_bisection_split(mblks,sj,1,npused+1+np2,JDIR)
               end do
            end do

            npused=npused+nppm
            write(STDOUT,1003) trim(nam)
         end do
         
         ! Update global mesh object with subblock information
         call set_block_info_global(nmesh,meshes,mblks)
         call process_block_bcs(meshes,mblks)
             
1000     format(/,'Attempting partition of mesh ',A, ' into ',I3,' parts')
1003     format('|-> Mesh ',A,' partitioning completed successfully')
      end subroutine perform_manual_split
   end subroutine manual_split_3d

   function get_nprocs_user_defined(meshes,palloc_file) result(nprocs)
      !! Parse the user defined splitting strategy input file and
      !! determine the total number of processors that will be
      !! necessary to perform the run.
      !!
      type(mesh_info), dimension(:), intent(in) :: meshes
      character(len=*), intent(in) :: palloc_file
      integer :: nprocs

      integer :: un, i, mid
      integer :: sj, sk, sl
      character(len=32) :: gnam
      ! Number of participating meshes
      integer :: nmesh
      ! Splitting information for each direction for each mesh
      integer, dimension(:,:), allocatable :: split_stats

      nmesh=size(meshes,1)
      allocate(split_stats(3,nmesh))
      ! Initialize the running nprocs for allocated processors and
      ! the split statistics to allocate 1 processor per mesh
      nprocs=0
      split_stats=1
      
      un=open_file(palloc_file,form='formatted',status='old')
      do i=1,nmesh
         read(un,*) gnam, sj, sk, sl
         mid=get_global_mesh_id(meshes,gnam)
         split_stats(JDIR,mid)=sj
         split_stats(KDIR,mid)=sk
         split_stats(LDIR,mid)=sl
         ! Determine the number of processors allocated to this
         ! mesh, and add it to a running nprocs.
         nprocs=nprocs+sj*sk*sl
      end do
      close(un)
   end function get_nprocs_user_defined
end module domain_partitioning

