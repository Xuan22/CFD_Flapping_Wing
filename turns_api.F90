module turns_api
!!! The public API for Fortran and Python drivers. The purpose of this
!!! module is to aggregate all major functionality under a few
!!! subroutines which can be called without needing to pass a lot of
!!! arguments along with it. The API should be fairly constant (and
!!! stable) and deal with modifications of internals gracefully.

   use ioutils
   use mpi_wrapper
   use domain_info
   use params_global
   use rootVariables
   use work
   use tightcouple
   use dualtime
   use deflections
   use arf_mod
   use pyMOD
   use domain_connectivity
   use airloads
   use movie_utils
! Math: add IBC, copy files: module/ibc_parameters.f90
!                            overturns/set_ibc.f90
!                            overturns/ibc_module.f90
   use SET_IMMERSED_BOUNDARY
   use IBC_MODULE

   implicit none

   private :: redirect_stdout_mpi, allocate_array_memory

contains 
   subroutine init(fname,comm)
      !! The entry point for any driver routine. This subroutine
      !! performs the following tasks:
      !!
      !! 0. Initializes MPI 
      !! 1. Reads in the input file and sets up the flow parameters
      !! 2. Reads the computational domain information and sets up the
      !!    global metadata for the meshes, and overset groups.
      !! 3. Performs domain partitioning and sets up the per-proc
      !!    computational block information
      !! 4. Reads in the global meshes, computes span/twist data.
      !! 5. Sets up local variables used by the solver.
      !!
      !! Inputs:
      !!    fname - Input file name
      !!    comm - MPI Communicator if running with several COMMs. 
      character(len=*), intent(in) :: fname
      integer, intent(in), optional :: comm

      integer :: un

      ! MPI boilerplate & and per-proc output redirection
      if (present(comm)) then 
         call initialize_mpi(comm)
      else
         call initialize_mpi_standalone
      end if
#ifdef HAVE_MPIF_H
      call redirect_stdout_mpi
#endif
      call overturns_intro
      
      ! Read inputs file & setup flight/rotor parameters
      un=open_file(fname,status='old')
      write(STDOUT,'(a)') "Reading inputs from... "//trim(adjustl(fname))
      call read_inputs(un)
      close(un)
      call ot_print_inputs
      call setup

      ! Read in information for all the meshes and its overset
      ! connectivity information
      ! also modifies ilim value from mesh_azi.inp, eg to use WENO on specific grid
      call setup_flow_domain_info 
      !print*,'     WARNING: Optimum Domain Partitioning Disabled For Wake Grids!!!!!!'
      call allocate_array_memory
      call setup_local_grid(x,y,z,iblank,gridbc,twist,spandist)

      ! Initialize some more variables
      call init_stage1(s,q,x,y,z,ug,vg,wg,kkr,kkp,xg,yg,zg,turmu,&
            vnu,twist,this_block%block_name)

      if (idual == 1) call dualtime_init
      if (imovie == 1) call load_movie_info

      if (novset > 0) call setup_connectivity_data
      write(STDOUT,*) 
      call print_message("Initialized OVERTURNS successfully")
      call barrier_mpi
   end subroutine init

   subroutine run(nsteps)
      !! Take a CFD timestep. Actually take 'nsteps' timesteps. This
      !! assumes that the structure of the participating meshes
      !! remains unchanged during the Newton sub-iterations.
      !!
      !! Inputs:
      !!    nsteps - Number of time-steps to take
      integer, intent(in) :: nsteps

      integer :: ntac_loc, nstart

      if (turnsInit == 0) call init_final_phase
      call barrier_mpi

      nstart=1
      nrc2=nsteps
      srot=rf*dt ! srot now represents (\Delta \psi)
      itn=1 ! This is needed in deform routine
      ! Iterations begin here
      time_steps: do istep=nstart,nrc2
               ! Camli _init
      	 if (isolno==1)then
             call process_outputs(dotest=.false.)
             stop
      	 end if

         istep0=istep0+1
         ntac_loc=1
         if ((ntac == 2) .and. (istep > 1)) ntac_loc=2

         ! Modify space and time metrics for unsteady flight
         call update_metrics_unsteady(ntac_loc)
         ! Store previous solution before re-scaling by Jacobian
         if (itnmax > 1) call stqol(q,qnewt,qtn,vnu,vnu0)

         write(STDOUT,101) istep0,rf*totime*r2d,totime
         if (idual == 1) call dualtime_reset_sequence()

         ! Math: Perform IBC interpolation of q variables
         call barrier_mpi
         if (this_mesh%dom_type.eq.DOM_IBC) call set_ibc(q)

         call barrier_mpi

         ! Perform Newton sub-iterations
         newton_iter: do itn=1,itnmax
            call perform_newton_subiter
         end do newton_iter

         call compute_forces(istep)
         call process_outputs(dotest=.true.)
      end do time_steps

      call barrier_mpi
      ! Store deflections for the next coupling cycle. This is needed
      ! for smoothed deflection transition
      if (is_wing .and. &
          (ideform == DEFORMING)) then
         defl_dat_prev=defl_dat
         def_first=.false.
      end if

      call process_outputs(dotest=.false.)
      if (pythonMode == 1) call collect_mesh_forces(nsteps)
101   format(/,' istep,azimuth,time =',1x,i5,2(1x,f10.5))

      call barrier_mpi
   
   end subroutine run

   subroutine redirect_stdout_mpi
      !! In MPI mode redirect output from each processor into a
      !! separate file: proc.out_PID. 
      character(len=10) :: intstr
      character(len=FLNLEN) :: outfile

      write(intstr,'(I10)') PID
      outfile='proc.out_'//trim(adjustl(intstr))
      open(unit=6,file=outfile,status='unknown')
   end subroutine redirect_stdout_mpi

   subroutine allocate_array_memory
      !! Allocate all the arrays declared in rootVariables and
      !! tightcouple
      integer :: bufsize, n

      call print_message("Allocating arrays")
      allocate(s(ssize),q(qsize))
      allocate(qnewt(qsize),qtn(qsize))
      allocate(x(gsize),y(gsize),z(gsize),bt(gsize))
      allocate(xg(gsize),yg(gsize),zg(gsize))
      allocate(xt2(gsize),yt2(gsize),zt2(gsize))

      allocate(xt2ref(gsize), yt2ref(gsize), zt2ref(gsize))
      allocate(svjref(gsize), svkref(gsize), svlref(gsize))
      
      allocate(iblank(gsize))

      allocate(xx(gsize),xy(gsize),xz(gsize))
      allocate(yx(gsize),yy(gsize),yz(gsize))
      allocate(zx(gsize),zy(gsize),zz(gsize))
      
      allocate(ug(gsize),vg(gsize),wg(gsize))
      allocate(turmu(gsize),vnu(gsize),vnu0(gsize))

      allocate(zx0(size2d),zy0(size2d),zz0(size2d))
      allocate(zt0(size2d))

      allocate(wvec(mdim,25))
      allocate(wgrid(gsize))

      allocate(z1(size2dp),z2(size2dp))
      allocate(loc1(size2dp),loc2(size2dp))

      allocate(kkr(size1d),kkp(size1d),twist(size1d),spandist(size1d))
      allocate(uge(gsize),vge(gsize),wge(gsize),ppp(gsize),d(gsize))

      ! Allocate GCL/FD data
      allocate(svj(gsize), svk(gsize), svl(gsize))
      allocate(vnaj(gsize), vnak(gsize), vnal(gsize))
      allocate(xaj(gsize), yaj(gsize), zaj(gsize))
      allocate(xak(gsize), yak(gsize), zak(gsize))
      allocate(xal(gsize), yal(gsize), zal(gsize))

      ! attach buffer for MPI_Bsend
      bufsize=0
      dsize3=0
      n=0
      do n=1,nmeshes
         bufsize=max(bufsize,meshes(n)%nfaces)
         dsize3=max(dsize3,meshes(n)%ncells)
      end do
      bufsize=bufsize*40
      dsize=dsize3*3
      !bufsize=m_jm*m_km*m_lm*10*8
      allocate(imesh(dsize),idonor(dsize),frac(dsize))
      allocate(imesh1(dsize),idonor1(dsize),frac1(dsize))
      allocate(buffer(bufsize))
      call mpi_buffer_attach(buffer,bufsize,ierr)

   end subroutine allocate_array_memory

   subroutine setup_blade_defl_radial(deftype,rstruct,bldefl,nr,ploc)
      !! Set up the blade deformation information at K-planes based on
      !! the data obtained from structural model. The data is
      !! interpolated onto the base mesh and then the relevant data is
      !! inserted into the computational block arrays to be used by
      !! the deform routines. This subroutine only updates for one
      !! azimuthal position. This allows it to be used both in
      !! periodic and time-accurate runs.

      ! deftype - Type of deflection (see TURNS_options DEFL_*)
      ! nr - number of radial stations where structural data is available
      ! ploc - Azimuth location where structural data is given
      integer, intent(in) :: deftype, nr, ploc

      ! rstruct - {x,y,z} coords where structural data is given
      ! bldefl - size: (6,nr) blade deformations from structures
      real(kind=rdp), dimension(:,:), intent(inout) :: rstruct, bldefl

      ! KTFIX, KRFIX - fudge factors to prevent bad metrics at
      ! tip/root caps.
      integer, parameter :: KTFIX=7, KRFIX=5
      integer :: k, i, km, k1, kr, kt
      real(kind=rdp), dimension(:), allocatable :: rdef,rdef1
      real(kind=rdp), dimension(:), allocatable :: rst

      ! Interpolate data to the base mesh locations first
      km=this_mesh%km
      allocate(rdef1(km),rdef(nr),rst(nr))

      ! For C-O type meshes, we shouldn't interpolate for the curved
      ! section of the tip/root caps. We force them to deform the same
      ! as the last almost-normal plane. THIS IS A FIX!!!
      if (m_ktip == km) then 
         kt=m_ktip-KTFIX
      else
         kt=m_ktip
      end if
      if (m_kroot == 1) then 
         kr=KRFIX
      else
         kr=m_kroot
      end if

      ! Fix data if they are of the RCAS/Dymore style
      if (deftype /= DEFL_TDU) then 
         do k=1,nr
            do i=4,6
               bldefl(i,k)=bldefl(i,k)*d2r
            end do
         end do
      end if

      ! Interpolating deformations
      rst=rstruct(1,:)
      do i=1,6
         rdef=bldefl(i,:)

         call interp1(rst,rdef,m_span,rdef1,nr,km)
         !call spline(rst,rdef,m_span,rdef1,nr,km)

         ! Insert relevant section into deform routine arrays
         do k=jkl_lims(KDIR),jkl_lims(KDIR+3)
            k1=k-jkl_lims(KDIR)+1
            if ((k-kr)*(k-kt) <= 0) then 
               defl_dat(i,k1,ploc)=rdef1(k)
            else if (k <= kr) then 
               defl_dat(i,k1,ploc)=rdef1(kr)
            else if (k >= kt) then 
               defl_dat(i,k1,ploc)=rdef1(kt)
            end if
         end do
      end do

      ! Setting the reference span position
      do k=jkl_lims(KDIR),jkl_lims(KDIR+3)
         k1=k-jkl_lims(KDIR)+1
         if ((k-kr)*(k-kt) <= 0) then 
            xyzref(1,k1)=m_span(k)*rartio
         else if (k <= kr) then 
            xyzref(1,k1)=m_span(kr)*rartio
         else if (k >= kt) then 
            xyzref(1,k1)=m_span(kt)*rartio
         end if
      end do
      
      ! Interpolating y, z coordinates
      do i=2,3
         rdef=rstruct(i,:)
         
         call interp1(rst,rdef,m_span,rdef1,nr,km)
         do k=jkl_lims(KDIR),jkl_lims(KDIR+3)
            k1=k-jkl_lims(KDIR)+1
            if ((k-kr)*(k-kt) <= 0) then 
               xyzref(i,k1)=rdef1(k)
            else if (k <= kr) then 
               xyzref(i,k1)=rdef1(kr)
            else if (k >= kt) then 
               xyzref(i,k1)=rdef1(kt)
            end if
         end do
      end do

      deallocate(rst,rdef,rdef1)
   end subroutine setup_blade_defl_radial

   subroutine setup_blade_defl_periodic(deftype,rstruct,bldefl,nr,np,&
         iazimuth)
      !! Set up blade deformation information at all radial stations
      !! and azimuthal positions based on the structural data provided
      !! by the inputs.

      ! deftype - deflection type (TDU,EULER,RCAS)
      ! nr - number of radial positions where deflections are given
      ! np - number of azimuthal positions where deflections are given
      ! iazimuth - number of azimuthal positions (CFD) where data is
      !            to be interpolated.
      ! rstruct - {x,y,z} positions where deflections are given
      ! bldefl - All the blade deflections (nd=6)
      integer, intent(in) :: deftype, nr, np, iazimuth
      real(kind=rdp), dimension(:,:), intent(inout) :: rstruct
      real(kind=rdp), dimension(:,:,:), intent(inout) :: bldefl

      integer :: p, i, k
      real(kind=rdp), dimension(:,:,:), allocatable :: def_tmp
      real(kind=rdp), dimension(:), allocatable :: defl, defl1
      
      allocate(defl(np),defl1(iazimuth),def_tmp(6,nr,iazimuth))

      ! Azimuthal interpolation to CFD stations
      do i=1,6
         do k=1,nr
            do p=1,np
               defl(p)=bldefl(i,k,p)
            end do
            call fourier_extend(defl,defl1,np,iazimuth)
            do p=1,iazimuth
               def_tmp(i,k,p)=defl1(p)
            end do
         end do
      end do

      ! At each azimuth, interpolate to CFD radial stations
      do p=1,iazimuth
         call setup_blade_defl_radial(deftype,rstruct,def_tmp(:,:,p), &
               nr,p)
      end do
      deallocate(defl,defl1,def_tmp)
   end subroutine setup_blade_defl_periodic

   subroutine update_metrics_rigid_rotation(ntac_loc)
      !! Perform rigid blade rotation of the participating meshes and
      !! update the spatial and time metrics.
      !!
      !! Inputs:
      !!    ntac_loc - Current time accuracy.
      integer, intent(in) :: ntac_loc

      tef: if ((iteflap == 1) .and. (dom_type == DOM_BLADE)) then
         call rotate_rigid_flap(srot,x,y,z,xx,xy,xz,ug,yx,yy,yz,vg, &
               zx,zy,zz,wg,zx0,zy0,zz0,zt0,0)

         call qmulj(q)
         usegcl1: if (usegcl) then
            call metgcl(ntac_loc,kkp,kkr,x,y,z,xt2,yt2,zt2,svj,svk,svl,&
                  xx,xy,xz,yx,yy,yz,zx,zy,zz,xaj,yaj,zaj,xak,yak,zak, &
                  xal,yal,zal,q,wgrid,vnaj,vnak,vnal,0)
         else
            call metfv(q,x,y,z,xx,xy,xz,yx,yy,yz,zx,zy,zz,wgrid,kkp,kkr,0)
         end if usegcl1
         call qdivj(q)
      else tef
         usegcl2: if (usegcl) then
            call rotate_grid(srot,x,y,z,xt2,yt2,zt2)
            call metgcl(ntac_loc,kkp,kkr,x,y,z,xt2,yt2,zt2,svj,svk,svl,&
                  xx,xy,xz,yx,yy,yz,zx,zy,zz,xaj,yaj,zaj,xak,yak,zak, &
                  xal,yal,zal,q,wgrid,vnaj,vnak,vnal,0)
         else usegcl2
            !if (irot_dir.eq.0) then
            !   srot=0.0
            !endif
            !print*,'PID irot_dir ilim rf srot',PID,irot_dir,ilim,rf,srot
            call rotate(srot,x,y,z,xx,xy,xz,ug,yx,yy,yz,vg, &
                  zx,zy,zz,wg,zx0,zy0,zz0,zt0)
         end if usegcl2
         call mett(x,y,z,xx,xy,xz,yx,yy,yz,zx,zy,zz,ug,vg,wg)
      end if tef
   end subroutine update_metrics_rigid_rotation

   ! Math: add kinematics
   subroutine update_metrics_flap_pitch_plunge(ntac_loc)
      !! Perform rigid blade rotation of the participating meshes and
      !! update the spatial and time metrics.
      !!
      !! Inputs:
      !!    ntac_loc - Current time accuracy.
      integer, intent(in) :: ntac_loc

	  
      if (dom_type == DOM_BLADE) then 
         if (iunst == FL_PITCH) call pitch_plunge_wing(x,y,z,xg,yg,zg,ug,vg,wg,zx,zy,zz,zx0,zy0,zz0,zt0,0)
         if (iunst == FL_FLAP)  call flap_pitch_wing(x,y,z,xg,yg,zg,ug,vg,wg,zx,zy,zz,zx0,zy0,zz0,zt0,0)

         call qmulj(q)
         usegcl1: if (usegcl) then
            call metgcl(ntac_loc,kkp,kkr,x,y,z,xt2,yt2,zt2,svj,svk,svl,&
                  xx,xy,xz,yx,yy,yz,zx,zy,zz,xaj,yaj,zaj,xak,yak,zak, &
                  xal,yal,zal,q,wgrid,vnaj,vnak,vnal,0)
         else
            call metfv(q,x,y,z,xx,xy,xz,yx,yy,yz,zx,zy,zz,wgrid,kkp,kkr,0)
         end if usegcl1
         call qdivj(q)
      endif

   end subroutine update_metrics_flap_pitch_plunge

   subroutine update_metrics_deforming(ntac_loc)
      !! Rotate and deform the participating meshes, and update the
      !! spatial and time metrics.
      !!
      !! Inputs:
      !!    ntac_loc - Current time accuracy
      integer, intent(in) :: ntac_loc

      tef: if (iteflap == 1) then
         tef_def: if (flpdef == 0) then
            tef_ovh: if (flpovh > 0.0) then
               call rigid_flapoh(xg,yg,zg,psi_rot,psi_rot-srot,.false.)
            else
               call rigid_flap(xg,yg,zg,psi_rot,psi_rot-srot,.false.)
            end if tef_ovh
         else
            call deflect_flap(xg,yg,zg,psi_rot,psi_rot-srot,.false.)
         end if tef_def
      end if tef

      call deform(srot,psi_rot,x,y,z,xx,xy,xz,ug,yx,yy,yz,vg,zx,zy,zz, &
            wg,zx0,zy0,zz0,zt0,xt2,yt2,zt2,xg,yg,zg)
      
      call qmulj(q)
      usegcl1: if (usegcl) then
         call metgcl(ntac_loc,kkp,kkr,x,y,z,xt2,yt2,zt2,svj,svk,svl,&
               xx,xy,xz,yx,yy,yz,zx,zy,zz,xaj,yaj,zaj,xak,yak,zak, &
               xal,yal,zal,q,wgrid,vnaj,vnak,vnal,0)
      else
         call metfv(q,x,y,z,xx,xy,xz,yx,yy,yz,zx,zy,zz,wgrid,kkp,kkr,0)
      end if usegcl1
      call qdivj(q)
   end subroutine update_metrics_deforming

   subroutine update_metrics_unsteady(ntac_loc)
      !! Update the time, mesh coordinates, spatial and time metrics
      !! in unsteady flight mode.
      integer, intent(in) :: ntac_loc

      unsteady: if ((iunst == FL_UNSTEADY) .or. &
                    (iunst == FL_HOVER)) then
         totime=totime+dt
         psi=rf*(totime-totime0)
         psi_rot=rf*totime+psi_offset
         
         if (ideform /= DEFORMING) then
            call update_metrics_rigid_rotation(ntac_loc)
         else
            call update_metrics_deforming(ntac_loc)
         end if
         
         if (arf_opt == 2) call add_arf_gridvel(x,y,z,ug,vg,wg)
         if (ipert == -3) call pert3(q,x,y,z,ug,vg,wg,zt0,zx,zy,zz,psi_rot)

         if ((novset > 0) .and. (conn_mode == CONN_DYNAMIC)) then
            call perform_domain_connectivity(x,y,z,iblank,.false.)
            call fix_iblanks(iblank,gridbc)
         end if

      ! Math: add kinematics
      else if (iunst == FL_PITCH .or. iunst == FL_FLAP) then unsteady
         totime=totime+dt
         call update_metrics_flap_pitch_plunge(ntac_loc)

         if ((novset > 0) .and. (conn_mode == CONN_DYNAMIC)) then
            call perform_domain_connectivity(x,y,z,iblank,.false.)
            call fix_iblanks(iblank,gridbc)
         end if

      end if unsteady

      if (iunst == FL_HOVER) call rotation_term(srot,q,qtn)
      
      ! if ((iunst == FL_UNSTEADY) .and. (novset > 0)) &
      !       call perform_domain_connectivity(x,y,z,iblank,.false.)
   end subroutine update_metrics_unsteady

   subroutine perform_newton_subiter
      !! Perform one Newton sub-iteration
      
      call step(x,y,z,iblank,q,qnewt,s,turmu,vnu,vnu0,xx,xy,xz, &
            ug,yx,yy,yz,vg,zx,zy,zz,wg,zx0,zy0,zz0,zt0,bt,&
            xaj,yaj,zaj,xak,yak,zak,xal,yal,zal,vnaj,vnak,vnal,&
            kkp,kkr,wgrid,wvec,resmax,resrho,rsum,gridbc)
      
      print '(I5,A,I8,I5)', PID, ": Step completed ",istep0,itn
      call monitor(mstop,q)
      if (mstop > 0) then
         call stop_execution('perform_newton_subiter',&
               'Negative speed of sound')
      endif

      call bc(x,y,z,q,vnu,xx,xy,xz,ug,yx,yy,yz,vg, &
            zx,zy,zz,wg,bt,zx0,zy0,zz0,zt0,kkr,kkp, &
            gridbc)
      call barrier_mpi
      
      if (novset > 0) call perform_overset_interpolations(q,vnu)
      call barrier_mpi
      if (idual == 1) call dualtime_next_dtpseudo()
   end subroutine perform_newton_subiter

   subroutine compute_forces(fmidx)
      !! Compute forces for the current mesh block.
      !!
      !! Inputs:
      !!   fmidx - Index where the airloads should be inserted in the
      !!           local loads array.
      integer, intent(in) :: fmidx

      logical, save :: fopened=.false.
      integer, save :: ifile, ifile1, ifile2, ctf
      !!!ifile2 added by Camli B. 9/2014

      ! Compute forces makes sense only if we're a blade/slat mesh
      !!! added by camli.b calculate forces for the body
      if (.not. is_wing .AND. .not. is_body) return 

      if (.not. fopened) then
         ifile=open_file(trim(adjustl(this_block%block_name))//'.128')
         ifile1=open_file(trim(adjustl(this_block%block_name))//'.138')
         ifile2=open_file(trim(adjustl(this_block%block_name))//'.10')
         ctf=open_file(trim(adjustl(this_block%block_name))//'.11')
         fopened=.true.
      end if

      call force2d(1,zero,zero,zero,ct,cq,figofmerit,&
            x,y,z,q,xx,xy,xz,zx,zy,zz,twist,fmidx,ifile,ifile1,ifile2)
      ! Math: add kinematics, copy force3.f file and add to Makefile
      call formom(1,ct,cq,figofmerit,&
            x,y,z,q,xx,xy,xz,yx,yy,yz,zx,zy,zz,twist,fmidx,ifile,ifile1,ifile2)

      write(ctf,100) istep0,totime,ct,cq,figofmerit
100   format(I5,4(E12.5,x))
   end subroutine compute_forces

   subroutine process_outputs(dotest)
      !! Wrapper function around movie/solution output subroutines.
      !!
      !! Inputs:
      !!    dotest - Perform tests to determine if we want to output
      !!    data. If .false. then we force output of solutions
      logical, intent(in) :: dotest

      logical :: do_mov, do_sol

      do_sol=.false.
      do_mov=.false.
      if (dotest) then
         if (mod(istep0,nrest) == 0) do_sol=.true.

       ! Added by James Lankford: Additional if statement
         ! to save restart files for a given interval (2/18/2014)
         if ((mod(istep0,ninterval) == 0) .and. &
             (istep0.ge.nstart) .and. &
             (istep0.le.nstop)) do_sol=.true.

         if ((imovie == 1) .and. &
             (mod(istep0,nmovie) == 0)) do_mov=.true.
         restart_counter = restart_counter + 1
         if (restart_counter.eq.2) restart_counter = 0
      else
         do_sol=.true.
         restart_counter = 0
      end if

      if (do_sol .or. do_mov) then
         call assemble_mesh(x,y,z)
         call assemble_soln(q,vnu,iblank)
      end if
      if (do_sol) call store_solutions
      if (do_mov) call store_movie

      call barrier_mpi
   end subroutine process_outputs

   subroutine init_final_phase
      ! Deform and pre-rotate the grids 
      call init_stage2(q,x,y,z,xx,xy,xz,yx,yy,yz,zx,zy,zz, &
            ug,vg,wg,xaj,yaj,zaj,xak,yak,zak,xal,yal,zal, &
            vnaj,vnak,vnal,svj,svk,svl,wgrid,kkr,kkp, &
            zx0,zy0,zz0,zt0,xg,yg,zg,xt2,yt2,zt2)
      
      ! Set up connectivity information
      psi=zero
      if (novset > 0) then
         call perform_domain_connectivity(x,y,z,iblank,.true.)
         call fix_iblanks(iblank,gridbc)
      end if
      call initialize_airloads_data
      
      ! Math: initialize IBC (cut hole and find interpolation coefs)
      ! Need overset data so has to do connectivity first
      ! It's ok if iblk_ibc not set for first connectivity call
      !if (this_mesh%dom_type.eq.DOM_IBC) call init_ibc(x,y,z,iblank)
      call init_ibc(x,y,z,iblank) ! ok even if no IBC grid

      psi=rf*totime
      istep=0
      turnsInit=1
   end subroutine init_final_phase

   subroutine ta_new_iter(dpsi)
      !! This subroutine is invoked at the start of every timestep. It
      !! updates the global time indexes, the total time variables,
      !! and other variables. For undeforming meshes, it rotates and
      !! updates the grid metrics. Finally, processes movie outputs. 
      real(kind=rdp), intent(in) :: dpsi

      integer :: ntac_loc

      dt=dpsi*rf*d2r
      if (itnmax == 1) itnmax=2
      ntac_loc=1
      if ((ntac ==2) .and. (istep > 1)) ntac_loc=2
      itn=0

      if (turnsInit == 0) call init_final_phase

      istep=istep+1
      istep0=istep0+1
      srot=rf*dt
      totime=totime+dt
      psi=rf*(totime-totime0)
      psi_rot=rf*totime+psi_offset

      ! Store q variables for the next timestep.
      call stqol(q,qnewt,qtn,vnu,vnu0)
      ! Update wake meshes only
      if (.not. is_wing) then
         call update_metrics_rigid_rotation(ntac_loc)
         if (arf_opt == 2) call add_arf_gridvel(x,y,z,ug,vg,wg)
      else
         ! For other meshes store their values for metrics calculation.
         call store_gridconfig(xt2,yt2,zt2,xt2ref,yt2ref,zt2ref, &
               svj,svk,svl,svjref,svkref,svlref)
      end if

      if (idual == 1) call dualtime_reset_sequence()
      write(STDOUT,101) istep0,rf*totime*r2d,totime

      if ((imovie == 1) .and. (mod(istep0,nmovie) == 0)) then
         call store_movie
      endif
101   format(/,' istep,azimuth,time =',1x,i5,2(1x,f10.5))
   end subroutine ta_new_iter

   subroutine ta_exec_subiter
      !! Execute one sub-iteration in time-accurate coupling mode. In
      !! this case, the deforming meshes are redeformed based on the
      !! latest estimate of the blade structural deformations, then
      !! the domain connectivity information is updated. The step, BC,
      !! and overset interpolation routines are called next. Finally,
      !! a new estimate of aerodynamic forces are obtained.
      integer :: ntac_loc

      ntac_loc=1
      if ((ntac ==2) .and. (istep > 1)) ntac_loc=2
      itn=itn+1

      ! Mesh deformation for blade/slat meshes
      if (is_wing) then
         if (itn > 1) call restore_gridconfig(x,y,z,xt2,yt2,zt2,&
               xt2ref,yt2ref,zt2ref,svj,svk,svl,svjref,svkref,svlref)
         call update_metrics_deforming(ntac_loc)
         if (arf_opt == 2) call add_arf_gridvel(x,y,z,ug,vg,wg)
      end if

      if (novset > 0) then
         call perform_domain_connectivity(x,y,z,iblank,.false.)
         call fix_iblanks(iblank,gridbc)
      end if
      call perform_newton_subiter
      call compute_forces(1)
      if (idual == 1) call dualtime_next_dtpseudo
   end subroutine ta_exec_subiter

   subroutine overturns_intro
      character(len=16) :: date, time, zone
      character(len=6) :: cprocs
      integer, dimension(8) :: values

      write(STDOUT,'(/,A)') 'OVERset Transonic Unsteady Rotor Navier Stokes'
      write(STDOUT,'(A)') 'Version: '//PACKAGE_STRING
      write(STDOUT,'(A)') 'Compiled: '//__DATE__//' '//__TIME__
      write(STDOUT,'(A)') 'Send bug reports to: '//PACKAGE_BUGREPORT
      write(STDOUT,'(A,/)') 'See AUTHORS file for developer information'

      call date_and_time(date,time,zone,values)
      write(STDOUT,'(A)') 'Run started: '//date(1:4)//'-'//date(5:6)//'-'//date(7:8)//" at "//time(1:2)//':'//time(3:4)//':'//time(5:6)

#ifdef HAVE_MPIF_H
      write(cprocs,'(I6)') NPROCS
      write(STDOUT,'(A)') 'Parallel execution running on '//trim(adjustl(cprocs))//' processors'
#endif
      write(STDOUT,'(A,/)') '==================================================================='
   end subroutine overturns_intro

   subroutine ot_print_inputs
      integer :: istop

      istop=0
      write(STDOUT,'(A)') 'Summary of inputs: '
      select case (IREAD)
      case (0)
         write(STDOUT,1000) "Initial run without restart"
      case (1)
         write(STDOUT,1000) "Restart from previous solution"
         write(STDOUT,1002) "Solutions to be read from"//trim(adjustl(rest_dir))
      case default
         istop=istop+1
         write(STDOUT,1001) "IREAD must be 0 or 1"
      end select
      write(STDOUT,'(4x,A,I8,A)') "... Write restart files every ",nrest," timesteps" 
      if (imovie == 1) then
         write(STDOUT,'(4x,A,I8,A)') "... Output movie files every ",nmovie," timesteps"
      end if

      write(STDOUT,1000) 'Solver options'
      select case (iunst)
      case (FL_STEADY)
         write(STDOUT,1002) 'Steady flow conditions'
      case (FL_UNSTEADY)
         write(STDOUT,1002) 'Unsteady flow conditions'
         if (timeac /= one) then
            istop=istop+1
            write(STDOUT,1001) 'Invalid TIMEAC for unsteady flow'
         end if
      case (FL_HOVER)
         write(STDOUT,1002) 'Hovering rotor (using source terms)'
      case (FL_QUASI)
         write(STDOUT,1002) 'Rotor in quasi-steady flight'
      ! Math: add kinematics
      case (FL_PITCH)
         write(STDOUT,1002) 'Pitching and plunging wing'
      case (FL_FLAP)
         write(STDOUT,1002) 'Flaping and pitching wing'
      case default
         istop=istop+1
         write(STDOUT,1001) 'IUNST must be between 0 and 6' ! Math: FL_PITCH=5, FL_FLAP=6
      end select

      if (invisc) then
         write(STDOUT,1002) 'Inviscid flow conditions'
      else
         write(STDOUT,1002) 'Viscous flow'
         if (ithin) then
            write(STDOUT,1003) 'Thin boundary layer assumptions'
         else
            write(STDOUT,1003) 'Fully viscous flow'
         end if
         if (lamin) then
            write(STDOUT,1003) 'Laminar flow'
         else
            select case (iturb)
            case (BALDWIN_LOMAX)
               write(STDOUT,1003) 'Turbulent flow: Baldwin Lomax model'
            case (SPALART_ALLMARAS)
               write(STDOUT,1003) 'Turbulent flow: Spalart Allmaras model'
            case default
               write(STDOUT,1001) 'Invalid ITURB: either 0 or 1'
            end select
         end if
      end if

      select case (irhsy)
      case (-1)
         write(STDOUT,1002) 'First order in space'
      case (-2)
         write(STDOUT,1002) 'Second order in space'
      case (-3)
         write(STDOUT,1002) 'Third order in space'
      case default
         istop=istop+1
         write(STDOUT,1001) 'Invalid IRHSY: -3 <= IRHSY -1'
      end select

      select case (ntac)
      case (1)
         write(STDOUT,1002) 'First order in time'
      case (2)
         write(STDOUT,1002) 'Second order in time'
      case default
         istop=istop+1
         write(STDOUT,1001) 'NTAC must be either 1 or 2'
      end select

      select case (itnmax)
      ! Math: add kinematics
      case (0)
         write(STDOUT,'(8x,a,I3)') "|-> Newton subiterations: ",itnmax
         write(STDOUT,'(8x,a)') "|--> Used only for testing grid motion" 
      case (1:)
         write(STDOUT,'(8x,a,I3)') "|-> Newton subiterations: ",itnmax
      case default
         istop=istop+1
         ! Math: itnmax=0 for testing grid motion
         write(STDOUT,1001) 'Invalid ITNMAX: must be greater than or equal to zero'
      end select

      select case (idual)
      case (0)
         write(STDOUT,1002) 'No dual time sequencing'
      case (1)
         write(STDOUT,1002) 'Using dual time sequencing from dtpseduo.inp'
      case default
         istop=istop+1
         write(STDOUT,1001) 'Invalid IDUAL; must be 0 or 1'
      end select

      select case (ilim)
      case (LIM_ALL)
         write(STDOUT,1002) 'Limiting in all 3 directions'
      case (LIM_JDIR)
         write(STDOUT,1002) 'Limiting in J-direction only'
      case (LIM_JTIP)
         write(STDOUT,1002) 'Limiting in tip region for J-direction only'
      case (LIM_NONE)
         write(STDOUT,1002) 'No limiting'
      case default
         istop=istop+1
         write(STDOUT,1001) 'Invalid ILIM: -1 <= ILIM < 2'
      end select

      if (usegcl) then
         write(STDOUT,1002) 'GCL based mesh metrics calculation'
      else
         write(STDOUT,1002) 'Default mesh metrics calculation procedure'
      end if

      select case (conn_type)
      case (CONN_JAINA)
         write(STDOUT,1002) 'Using traditional hole-cutting for overset runs'
      case (CONN_IHC)
         write(STDOUT,1002) 'Using IHC for overset runs'
      case default
         istop=istop+1
         write(STDOUT,1001) 'Invalid CONN_TYPE: either 0 or 1'
      end select

      if (pythonmode == 0) then 

         write(STDOUT,1000) "Flight conditions"
         write(STDOUT,"(8x,A,I7)") "|-> Total number of timesteps: ",nsteps
         if (dt < 0) then
            write(STDOUT,'(8x,A,F8.4,A)') "|-> Timestep size: ",abs(dt)," degrees"
         else
            write(STDOUT,'(8x,A,F12.6,A)') "|-> Timestep size: ",dt," seconds"
         end if
         write(STDOUT,"(8x,A,F12.2)") "|-> Flow Reynolds number: ", rey
         write(STDOUT,1004) "Aspect Ratio: ", rartio
         write(STDOUT,1004) "Freestream Mach number: ", fsmach
         if (fmtip > zero) then 
            write(STDOUT,1004) "Tip Mach number: ", fmtip
            write(STDOUT,1004) "Shaft tilt angle: ",alfa
         else
            write(STDOUT,1004) "Angle of attack: ", alfa
         end if
         select case (ideform)
         case (NO_DEFORM)
            write(STDOUT,1002) "No deformations: Rigid blade rotation"
         case (DEFORMING)
            write(STDOUT,1002) "Flexible blade dynamics with..."
            select case (idefl)
            case (DEFL_TDU)
               write(STDOUT,1003) "UMARC deformations"
            case (DEFL_EULER)
               write(STDOUT,1003) "Dymore deformations"
            case (DEFL_RCAS)
               write(STDOUT,1003) "RCAS/CAMARAD deformations"
            case default
               write(STDOUT,1001) "Invalid IDEFL; 0 <= IDEFL <= 2"
            end select
         case default
            istop=istop+1
            write(STDOUT,1001) "Invalid IDEFORM; must be 0 or 1"
         end select
      else
         write(STDOUT,1000) "Flow parameters controlled via Python interface"
      end if

      write(STDOUT,*)
      write(STDOUT,'(A,/)') '==================================================================='
      if (istop > 0) then
         write(STDOUT,'(I3,a)') istop, " error(s) in input file. Exiting!"
         write(STDERR,'(I3,a)') istop, " error(s) in input file. Exiting!"
         stop
      end if

1000  format(4x,'... ',A)
1001  format('***ERROR: ',A,' ***')
1002  format(8x,'|-> ',A)
1003  format(8x,'|   |-> ',A)
1004  format(8x,'|-> ',A,F12.6)
   end subroutine ot_print_inputs
end module turns_api
