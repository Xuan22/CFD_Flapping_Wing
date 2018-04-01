!
!           MAIN PYTHON-OVERTURNS_ARF INTERFACE PROCEDURES
!
! DEVELOPERS:
!   Jay Sitaraman: original work, 02/08/05
!   Benjamin Silbaugh: modified for tight-coupling and Accelerating
!     Reference Frame (ARF)
!
! LAST MODIFIED:
!   7/3/07 (Modified execute and set_blade_def)
!==========================================================================


      subroutine time_step(dpsi)
      !DESC: Executes a full time step with the number of prescribed
      !      newton subs. This should should be used in place of 
      !      new_iteration() and exec_subitertion() if grid deformations
      !      and airloads are NOT to be exechanged every sub-iteration.

      use params_global, only: dt, pi, rf
      use pyMOD
      use turns_api, only : run
      implicit none

      real dpsi
      integer :: n

      !Compute time step
      dt = dpsi*pi/rf/180.

      ! Execute one time step
      n = 1
      call run(n)
      
      end subroutine time_step

      subroutine new_iteration(dpsi)
      !DESC: Call this at the begining of a new sub-iteration
      !      sequence. This (along with exec_subiteration) should
      !      be used in place of time_step() if airloads and grid
      !      deformations are to be exchanged every sub-iteration.

      use turns_api, only: ta_new_iter
      implicit none

      real dpsi

      call ta_new_iter(dpsi)

      end subroutine new_iteration

      subroutine exec_subiteration()
      !DESC: Exectutes a single sub-iteration
      use turns_api, only: ta_exec_subiter
      
        call ta_exec_subiter()
      end subroutine exec_subiteration

      
      subroutine startup_arf(inp_file, arf_option)
      ! Essentially the same start-up procedure as in pyinterface.f
      ! except setting arf params to default values and setting arf =
      ! .TRUE.
      !
      ! arf_option = 1 => use source terms to account for vehicle motion
      ! arf_option = 2 => use grid velocities to account for vehicle motion

      use turns_api, only: init
      use pyMOD
      use arf_mod

      implicit none

      character*128 inp_file
      integer :: arf_option

      
cf2py intent(in) :: inp_file

      pythonMode=1

      arf_opt = arf_option

      call init(inp_file)

      call init_arf_mod()

      if(allocated(loads)) deallocate(loads)
      allocate(loads(7,nsa, 1))
      write(6,*) 'Allocated array for collecting airloads'

      end subroutine startup_arf

      function get_azimuth()
      ! An inquire function that returns the current azimuth position
      ! of "this" blade. On output the azimuth, psi, is expresed in
      ! radians. This is function is useful when performing restarts.

      use rootVariables, only: psi_rot
    
      implicit none

      real get_azimuth

      get_azimuth = psi_rot

      end function get_azimuth

      function get_procid()
      ! Returns process id number.

      use mpi_wrapper, only: PID
      implicit none

      integer get_procid

      get_procid = PID

      end function get_procid

      subroutine init_blade_motion(mname,ideff, xyzref1, defl_dat1, nrhat)
      !DESC: Allocates deformation arrays and stores initial
      !deformations.

      use params_global
      use deflections
      use rootVariables
      use pyMOD
      use mesh_types, only: get_global_mesh_id
      use domain_info, only: meshes
      use turns_api

      implicit none

      character*128 mname
      integer ideff, nrhat
      real xyzref1(3,nrhat), defl_dat1(6,nrhat)

! local var
      integer i, k
      integer mid

      ! If this is wake mesh, then we quit without doing anything
      if (.not. is_wing) then
         ideform=0
         return
      endif

      ! Check if the deflections are meant for this mesh block
      mid=get_global_mesh_id(meshes,mname)
      if (mid /= this_mesh%mesh_id) return
      
c..   allocate a few things
      
      write(6,*) 'def_alloc_flag=',def_alloc_flag
      write(6,*) 'nsa=',nsa
      write(6,*) 'kmax=',kmax

      if(allocated(defl_dat)) deallocate(defl_dat)
      if(allocated(xyzref)) deallocate(xyzref)
      if(allocated(defl_dat_prev)) deallocate(defl_dat_prev)
      allocate(defl_dat(6,kmax,2), xyzref(3,kmax))
      allocate(defl_dat_prev(6,kmax,2))
      def_alloc_flag = 1

      call setup_blade_defl_radial(ideff,xyzref1,defl_dat1,nrhat,2)

      do k = 1, kmax
        do i = 1, 6
          defl_dat(i,k,1) = defl_dat(i,k,2)
        end do
      end do
      defl_dat_prev=0.0D0

      ideform=1
      idefl=ideff
      end subroutine init_blade_motion


c********************************************************************
      subroutine set_blade_motion(mname,ideff,xyzref1,defl_dat1,nrhat)
      !DESC:
      ! Interpolates blade deformation from CSD spanwise stations
      ! to CFD spanwise stations. Spanwise interpolation is via
      ! cubic splines. Deformations are assumed to be in the rotating
      ! frame. (Grid deformation scheme will map deformations to
      ! the TURNS non-rotating frame.)
      !
      !INPUT:
      !  mname....... String identifier for the blade mesh
      !  ideff....... Format used to define beam displacement field
      !  xyzref1..... Undeformed reference points for beam sections
      !               shape = (3,nrhat)
      !  defl_dat1... Beam displacement and rotation parameters
      !               shape = (6,nrhat)

      use params_global
      use deflections
      use rootVariables
      use pyMOD
      use mesh_types, only: get_global_mesh_id
      use domain_info, only: meshes
      use turns_api

c*****************  list of global variables used *******************
c     
c     jmax,kmax,lmax,iazimuth,pi,rf,dt,ktip,theta0,theta1c,theta1s
c     idefl,jle,rartio -> (params_global)
c
c     xyzref,defl_dat,iazimuth -> (deflections)
c
c********************************************************************

      implicit none

      character*128 mname
      integer ideff,nrhat
      real xyzref1(3,nrhat),defl_dat1(6,nrhat)

cf2py intent(in) :: nrhat,xyzref1,defl_dat1

      
c..   local variables


!      real, allocatable:: temp(:)
!      real, allocatable :: defl(:),defl1(:)
!      real, allocatable :: rad(:),rad1(:)
!      real, allocatable :: rdef(:),rdef1(:)
!      
!      real thet,pif
!      integer i,j,k,l,kt,kr
!
!      real xpsi, dpsi !Added local var (bS)
      integer mid
      
c**** first executable statement
      
      ! If this is wake mesh, then we quit without doing anything
      if (.not. is_wing) then
         ideform=0
         return
      endif

      ! Check if the deflections are meant for this mesh block
      mid=get_global_mesh_id(meshes,mname)
      if (mid /= this_mesh%mesh_id) return

      call setup_blade_defl_radial(ideff,xyzref1,defl_dat1,nrhat,2)
c$$$      
c$$$      pif = pi/180.
c$$$
c$$$      call set_mesh_to_global(jgmx,kgmx,lgmx,
c$$$     &     jtail1mx,jtail2mx,jbotmx,jlemx,ktipmx,krootmx,nmesh,1)
c$$$
c$$$
c$$$      allocate(rad(nrhat),rad1(kmax))
c$$$      allocate(rdef(nrhat),rdef1(ktip))
c$$$
c$$$c..   modify blade motions appropriately
c$$$
c$$$!NOTE: For UMARC Motions, control pitch was originally superimposed
c$$$!      on phi here. Will now assume control motions have already been
c$$$!      incorporated as part of the blade deformation (as they should
c$$$!      be).
c$$$
c$$$         ! Convert rotation angles from deg to rad for euler angle
c$$$         ! inputs (e.g. RCAS motions).
c$$$
c$$$         if (idefl .ne. 0)then
c$$$	   do j=1,nrhat
c$$$	     defl_dat1(4,j)=defl_dat1(4,j)*pif
c$$$	     defl_dat1(5,j)=defl_dat1(5,j)*pif
c$$$	     defl_dat1(6,j)=defl_dat1(6,j)*pif
c$$$	   enddo	
c$$$	endif
c$$$
c$$$c..   now interpolate radially
c$$$
c$$$      ! set spanwise coords
c$$$
c$$$      do k=1,kmax
c$$$         rad1(k)=spandist(k)
c$$$      enddo
c$$$      
c$$$      do k=1,nrhat
c$$$         rad(k)=xyzref1(1,k)
c$$$      enddo
c$$$
c$$$      ! small fix to prevent the bad metrics at the tip (<= Why? BS)
c$$$
c$$$      if (ktip.eq.kmax) then
c$$$         !kt=kmax-7     ! original value by Jaina (works for some meshes)
c$$$         kt = kmax - 20 ! changed by BS (may wish to define as input grid param)
c$$$      else
c$$$         kt=ktip
c$$$      endif
c$$$
c$$$      if (root_co.eq.1) then
c$$$	kr=5
c$$$      else
c$$$	kr=kroot
c$$$      endif
c$$$	
c$$$
c$$$      do l=1,6
c$$$            do k=1,nrhat
c$$$               rdef(k)=defl_dat1(l,k)
c$$$            enddo
c$$$            
c$$$            !call interp1(rad,rdef,rad1,rdef1,nrhat,kt)
c$$$            
c$$$            call spline(rad,rdef,rad1,rdef1,nrhat,kt)
c$$$            
c$$$            do k=1,kt
c$$$               defl_dat(l,k,2)=rdef1(k)
c$$$            enddo
c$$$
c$$$c..   setting all outside tip to be equal to the tip
c$$$
c$$$            do k=kt+1,kmax
c$$$               defl_dat(l,k,2)=defl_dat(l,kt,2)
c$$$            enddo
c$$$	 
c$$$            do k=1,kr-1
c$$$	      defl_dat(l,k,2)=defl_dat(l,kr,2)
c$$$	    enddo
c$$$
c$$$      enddo
c$$$
c$$$c..   interpolate the point of rotation, i.e where the deflections
c$$$c..   were supplied also radially
c$$$
c$$$      ! x-coord same as grid
c$$$
c$$$      do k=1,kmax
c$$$         xyzref(1,k)=rad1(k)*rartio
c$$$      enddo
c$$$
c$$$      ! y,z coords need interpolation
c$$$
c$$$      do i=2,3
c$$$         do k=1,nrhat
c$$$            rdef(k)=xyzref1(i,k)
c$$$         enddo
c$$$         call spline(rad,rdef,rad1,rdef1,nrhat,kt)
c$$$
c$$$         do k=1,kt
c$$$            xyzref(i,k)=rdef1(k)*rartio
c$$$         enddo
c$$$      enddo
c$$$      
c$$$      do k=kt+1,kmax
c$$$         xyzref(1,k)=xyzref(1,kt)
c$$$         xyzref(2,k)=xyzref(2,kt)
c$$$         xyzref(3,k)=xyzref(3,kt)
c$$$      enddo
c$$$
c$$$      do k=1,kr-1
c$$$         xyzref(1,k)=xyzref(1,kr)
c$$$         xyzref(2,k)=xyzref(2,kr)
c$$$         xyzref(3,k)=xyzref(3,kr)
c$$$      enddo
c$$$
c$$$190   format(6(1X,E14.8))
c$$$
c$$$      deallocate(rad,rad1,rdef,rdef1)
	
      !write(6,*) 'set blade motions in turns'
      ideform=1
      idefl=ideff
      return
      end subroutine set_blade_motion


      subroutine set_frame_motion(uvw, uvw_dot, pqr) !dpsi,fsm,fmtt,alp)

      use params_global
      use arf_mod

      implicit none

      real uvw(3), uvw_dot(3), pqr(3)

!      real dpsi,fsm,alp,fmtt

!** cf2py intent (in) :: dpsi,fsm,alp,fmtt

!      fmtip=fmtt
!      fsmach=fsm
!      alf=alp
      
!      cs=cos(alf*pi/180.)
!      ss=sin(alf*pi/180.)

!      if(iunst.eq.0.or.fmtip.eq.0) then
!        uinf = fsmach*cs
!        vinf = uinf*tan(beta1*pi/180.)
!      else
!        uinf  = 0.0
!        vinf  = fsmach*cs
!      endif

      ! Scale input frame motion to be consistent with TURNS non-dim
      ! scheme and NR-frame:

      frame_vel(1) = -uvw(2)*fmtip
      frame_vel(2) =  uvw(1)*fmtip
      frame_vel(3) =  uvw(3)*fmtip

      frame_acc(1) = -uvw_dot(2)*fmtip**2
      frame_acc(2) =  uvw_dot(1)*fmtip**2
      frame_acc(3) =  uvw_dot(3)*fmtip**2

      frame_avel(1) = -pqr(2)*rf
      frame_avel(2) =  pqr(1)*rf
      frame_avel(3) =  pqr(3)*rf

      ! Set free-stream:

      if(arf_opt .eq. 1)then

        ! Using source terms, set to relative free-stream

        uinf = -frame_vel(1)
        vinf = -frame_vel(2)
        winf = -frame_vel(3)

      elseif(arf_opt .eq. 2)then

        ! Using grid velocities for all motion, free-stream is zero.

        ! Will treat unsteady vehicle motions as perturbation from
        ! steady. The effective free-steam remains constant.

        ! REMEMBER: (free-stream) = -(frame velocity)

        frame_vel(1) = frame_vel(1) + uinf
        frame_vel(2) = frame_vel(2) + vinf
        frame_vel(3) = frame_vel(3) + winf

      else
        print*, "ERROR: arf_opt not set!"
        stop
      end if

      end subroutine set_frame_motion

      subroutine write_restartfiles()

      use turns_api

      call process_outputs(dotest=.false.)
      
      end subroutine write_restartfiles

      subroutine ta_get_mesh_span_info(mname,mid,nsa)

      use domain_info, only: meshes,get_global_mesh_id
      use airloads, only: get_num_span_loc
      
      character*128 :: mname
      integer :: mid, nsa

cf2py intent(in) :: mname
cf2py intent(out) :: mid, nsa      
      
      mid=get_global_mesh_id(meshes,mname)
      nsa=get_num_span_loc(mid)

      end

      subroutine ta_get_airloads(mid,nstps,nsa,rdist,fmom)

      use airloads, only: send_airloads, collect_mesh_forces
      
      integer, intent(in) :: mid,nstps,nsa
      real, dimension(nsa), intent(out) :: rdist
      real, dimension(7,nsa,nstps), intent(out) :: fmom

      call collect_mesh_forces(nstps)
      call send_airloads(mid,nstps,rdist,fmom)

      end
      
