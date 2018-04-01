      subroutine set_defaults

      use params_global
      use pyMOD
      use rootVariables
      use deflections
      use freewake

      implicit none

      nmesh=1
      dx1=1.
      dy1=1.
      dz1=1.

      gamma=1.4
      pr=0.72
      rmue=1.0
      tinf=540.

      ksym=0.
      kroot=1.
      ibcwp=1.
      cnbr=1.
      iread=1.
      jmax=109
      kmax=36
      lmax=31
      jtail1=19
      ktip=27
      half=0
      nsteps=500
      nrest=1000
      npnorm=25
      nmovie=4000
      fsmach=0.16
      alfa=0.
      rey=3.9e6
      invisc=.FALSE.
      LAMIN=.FALSE.
      ITHIN=.FALSE.
      ITURB=1
! Math: add kinematics, Added by Ria

       isroot = 0
       root = 0
       flap_amp = 0
       pitch_angle = 20.0
       phase = 90.0
       pitch_knot = 0.0
       plunge_amp = 0.5
       omega = 0.6
       xac = 0.25
       zac = 0.0

      kin_rf = 0.1
      kin_freq = 20.0
      kin_phiMax = 60.0
      kin_phiOff = 0.0
      kin_phi_phase = 0.0
      kin_alphaMax = 45.0
      kin_alphaOff = 0.0
      kin_alpha_phase = 0.0
      kin_beta = 0.0 
      xoff = 0.0
      yoff = 0.0
      zoff = 0.0
      kin_phiInit = 0.0
      kin_alphaInit = 0.0
      kin_beta_new = 0.0
      kin_phiOff_new = 0.0
      kin_phiMax_new = 0.0
      isolno = 0
      ismth  = 0
      Ncyc = 0

      xcg = 0.0
      ycg = 0.0
      zcg = 0.0 

      iunst=3
      ntac=1
      itnmax=1
      dt=1.
      timeac=1.
      indicial=1.
      
      epse=0.01
      irhsy=-3
      ilim=0
      totime=0.0
      angmax=0.0

      fmtip=0.8
      rartio=8.0
      alpha_dot=0.
      beta1=0.

      iwake=0.
      ctinp=0.00666
      nblade=1

      nacou=100
      nacoup=20
      nakp1=1
      nakp2=1
      nakp3=1
      nakp4=1
      nakp5=1
      nww=1
      nwwq=1000

      jint=1
      kint=1
      lint=1

      ipert=0
      initpe=1
      xvor=-5
      zvor=-0.26
      rcore=0.05
      vorgam=0.2

      imth=1
      ivec=0
      ntot=72
      nrev=3
      
      theta0=0.
      theta1c=0.
      theta1s=0.
      xlam=0.

      irotate=0
      totime0=0.
      ideform=0
      iconfine=0
      irefine=0
      idefl=0
      idecay=0

      beta0=0
      beta1c=0
      beta1s=0
      eflap=0
      
      root_co=0
      
      pxc=0.8
      flap_ys=1.44
      flap_ye=2.31
      dely=0.2
      dela=30
      rf_tef=1.0
      angmax_tef=5.0
      theta_flap=0.
      theta_f0=0.
      iteflap=0
      ampFlap=0
      phiFlap=0
      flpdef=0

      pythonMode=0
      turnsInit=0
      dzeta=5
      icom=1
      def_alloc_flag=0

      ! Added by James Lankford for file reading routines (03/02/18)
      read_flap_kine = 0
      read_pitch_kine = 0

      return
      end

c*********************************************************************
      subroutine check_inputs

c***  Prologue : ***
c
c     Write summary of inputs (and check validity!)
c     
c     last updated by jaina 07/13/04
c***  end prologue ***************************************************

      use params_global
      use pyMOD, only: pythonMode

c*******************  list of global variables used *******************
c
c     iread,jmax,kmax,lmax,iwake,half,ktip,jtail1,jtail2,jle,nsteps
c     fsmach,alpha,invisc,rey,lamin,iunst,ntac,itnmax,dt,timeac
c     epse,irhsy,ilim,totime,fmtip,rartio,rf,nblade,ctinp,nacoup,nacou
c     jint,kint,lint
c     
c     basically almost all of params_global !!
c
c*********************************************************************
      implicit none

      integer istop,jtest

c***  first executable statment

      ISTOP = 0
      WRITE(STDOUT,*) ' '
      WRITE(STDOUT,*) 'HERE IS YOUR INPUT:'
     
      WRITE(STDOUT,*) 'Restart'
      IF(IREAD.EQ.0) THEN
        WRITE(STDOUT,*) '  Initial run w/o restart for this case '
      ELSEIF(IREAD.EQ.1) THEN
        WRITE(STDOUT,*) '  Restart run for this case '
      ELSE
        WRITE(STDOUT,*) '**Invalid IREAD, must be 0 or 1**'
        ISTOP = ISTOP+1
      ENDIF
     
c      WRITE(STDOUT,*) 'Grid Sizes'
c      JTEST = JMAX
c      IF(IWAKE.EQ.2) JTEST=JMAX+(2-HALF)
c
c      IF(KMAX.LT.3) THEN
c        WRITE(STDOUT,*) '**KMAX can not be smaller than 3**'
c        ISTOP = ISTOP+1
c      ELSEIF(KMAX.EQ.3) THEN
c        WRITE(STDOUT,*) '  2-D flow over an airfoil (with 3 planes)'
c        WRITE(STDOUT,*) '    ',JMAX,' points in wrap-around direction',
c     <            ' and ',LMAX,' points in the normal direction'
c        KTIP = 3
c      ELSE
c        WRITE(STDOUT,*) '  3-D flow over a wing or rotor'
c        WRITE(STDOUT,*) '    ',JMAX,' points in wrap-around direction',
c     <                ' ',KMAX,' points in span-wise direction'
c        WRITE(STDOUT,*) '    and ',LMAX,' points in the normal direction'
c        WRITE(STDOUT,*) '    Tip located at k=',KTIP
c      ENDIF
c
c      WRITE(STDOUT,*) '  Airfoil cross section from j=',JTAIL1,' to j=',
c     <     JTAIL2,' with jle=',JLE
c      IF(HALF.EQ.1) WRITE(STDOUT,*) '    Thus, symmetrical airfoil assumed'
c
      if (pythonMode /= 1) then 
      WRITE(STDOUT,*) '  This run will last until ',NSTEPS,' time steps'
     
      WRITE(STDOUT,*) 'Flow Parameters'
      WRITE(STDOUT,*) '  Free Stream Mach Number is ',FSMACH
      IF(ALFA.NE.0.) WRITE(STDOUT,*) '  Flow at ',
     <               ALFA,' degrees angle of attack'
      endif
      IF(INVISC) THEN
        WRITE(STDOUT,*) '  Inviscid Flow'
      ELSE
        WRITE(STDOUT,*) '  Viscous Flow'
        WRITE(STDOUT,*) '    Reynolds Number = ',REY
        IF(LAMIN) THEN
          WRITE(STDOUT,*) '      Laminar Flow'
        ELSE
          WRITE(STDOUT,*) '      Turbulent Flow'
        ENDIF
      ENDIF

      WRITE(STDOUT,*) 'Time info'
      
      IF(IUNST.EQ. FL_STEADY) THEN
        WRITE(STDOUT,*) '  Steady flow '
      ELSEIF(IUNST.EQ. FL_UNSTEADY) THEN
        WRITE(STDOUT,*) '  Rotor in forward flight '
      ELSEIF(IUNST.EQ. FL_HOVER) THEN
        WRITE(STDOUT,*) '  Rotor in hover (using source term) '
      ELSEIF(IUNST.EQ. FL_QUASI) THEN
        WRITE(STDOUT,*) '  Rotor in quasi-steady flight '
      ELSE
        WRITE(STDOUT,*) '**Invalid IUNST, must be between 0 and 3**'
        ISTOP = ISTOP+1
      ENDIF
      
      IF(NTAC.EQ.1) THEN
        WRITE(STDOUT,*) '  1st order in time '
      ELSEIF(NTAC.EQ.2) THEN
        WRITE(STDOUT,*) '  2nd order in time '
      ELSE
        WRITE(STDOUT,*) '**Invalid NTAC, must be 1 or 2**'
        ISTOP = ISTOP+1
      ENDIF
      IF(ITNMAX.EQ.1) THEN
        WRITE(STDOUT,*) '  No use of Newton Iterations '
      ELSEIF(ITNMAX.GT.1) THEN
        WRITE(STDOUT,*) '  Using Newton Iterations with ',
     <        ITNMAX,' iterations'
      ELSE
        WRITE(STDOUT,*) '**Invalid ITNMAX, must be 0 or greater**'
        ISTOP = ISTOP+1
      ENDIF
      IF(IUNST.EQ. FL_UNSTEADY .AND. DT.GT.0.5) THEN
        WRITE(STDOUT,*) '**Invalid DT for unsteady flow**'
        ISTOP=ISTOP+1
      ELSEIF(DT.GT.0.) THEN
        WRITE(STDOUT,*) '  DT is ',DT
      ENDIF
      IF(TIMEAC.NE.1.0 .AND. IUNST.EQ. FL_UNSTEADY) THEN
        WRITE(STDOUT,*) '**Invalid TIMEAC for unsteady flow**'
        ISTOP=ISTOP+1
      ELSEIF(TIMEAC.LE.1.0 .AND. TIMEAC.GE.0.0) THEN
        WRITE(STDOUT,*) '  TIMEAC for Jacobian Scaling is ',TIMEAC
      ELSE
        WRITE(STDOUT,*) '**Invalid TIMEAC, must be between 0 and 1**'
        ISTOP=ISTOP+1
      ENDIF
     
      WRITE(STDOUT,*) 'Algorithm stuff'
      IF(EPSE.GT. 0.0 .AND. EPSE.LT. 1.0) THEN
        WRITE(STDOUT,*) '  Dissipation for ILU3D is ',EPSE
      ELSE
        WRITE(STDOUT,*) '**Invalid EPSE, must be between 0 and 1**'
        ISTOP=ISTOP+1
      ENDIF
      IF(IRHSY.EQ.-1) THEN
        WRITE(STDOUT,*) '  1st order in space '
      ELSEIF(IRHSY.EQ.-2) THEN
        WRITE(STDOUT,*) '  2nd order in space '
      ELSEIF(IRHSY.EQ.-3) THEN
        WRITE(STDOUT,*) '  3rd order in space '
      ELSE
        WRITE(STDOUT,*) '**Invalid IRHSY, must be between -1 and -3**'
        ISTOP=ISTOP+1
      ENDIF
      IF(ILIM.EQ. LIM_ALL) THEN
        WRITE(STDOUT,*) '  Limiting in all three directions'
      ELSEIF(ILIM.EQ. LIM_JDIR) THEN
        WRITE(STDOUT,*) '  Limiting in only the J-direction'
      ELSEIF(ILIM.EQ. LIM_JTIP) THEN
        WRITE(STDOUT,*) '  Limiting in tip region in only the J-direction'
      ELSE
        WRITE(STDOUT,*) '  No Limiting'
      ENDIF
      IF(IREAD.EQ.0) THEN
        WRITE(STDOUT,*) '  TOTIME = ',TOTIME
      ELSE
        WRITE(STDOUT,*) '  TOTIME = ',TOTIME,' will be reset from Q-file'
      ENDIF

      IF(IUNST .ne. FL_STEADY) THEN
        WRITE(STDOUT,*) 'Rotor Flow Parameters'
        WRITE(STDOUT,*) '  Rotational tip Mach number is ',FMTIP
        WRITE(STDOUT,*) '  Aspect Ratio is ',RARTIO
        WRITE(STDOUT,*) '    Thus, reduced frequency is ',RF
      ENDIF
     
c      IF(IUNST.GE.2) WRITE(STDOUT,*) 'Wake Parameters'
c      IF(IWAKE.NE.0 .AND. IUNST.LT.2) THEN
c        WRITE(STDOUT,*) '**Invalid IWAKE, can only simulate wake for hover**'
c        ISTOP=ISTOP+1
c      ELSEIF(IWAKE.EQ.0) THEN
c        WRITE(STDOUT,*) '  No Wake Capturing '
c      ELSEIF(IWAKE.EQ.1) THEN
c        WRITE(STDOUT,*) '  Wake Capturing for ',NBLADE,' blades '
c        WRITE(STDOUT,*) '             CTINP = ',CTINP
c      ELSEIF(IWAKE.EQ.2) THEN
c        WRITE(STDOUT,*) '  Wake Capturing for ',NBLADE,' blades, add pts'
c        WRITE(STDOUT,*) '             CTINP = ',CTINP
c      ELSE
c        WRITE(STDOUT,*) '**Invalid IWAKE, must be between 0 and 2**'
c        ISTOP=ISTOP+1
c      ENDIF
c     
c      WRITE(STDOUT,*) 'Acoustic Parameters'
c      IF(IUNST.EQ. FL_UNSTEADY) THEN
c        WRITE(STDOUT,*) '  Write p-file info every ',NACOUP,' iterations'
c        WRITE(STDOUT,*) '  Write p-file for flow every ',NACOU,' iterations'
c      ENDIF
c     
      IF(JINT+KINT+LINT.GT.3) THEN
        WRITE(STDOUT,*) 'Coarsening of grid for this run, use with caution'
        IF(JINT.NE.1)
     <     WRITE(STDOUT,*) ' Use every ',JINT,' points in j-direction'
        IF(KINT.NE.1)
     <     WRITE(STDOUT,*) ' Use every ',KINT,' points in k-direction'
        IF(LINT.NE.1)
     <     WRITE(STDOUT,*) ' Use every ',LINT,' points in l-direction'
      ENDIF

      IF(ISTOP.GT.0) WRITE(STDOUT,*) 'NOTE: ',ISTOP,' errors in input******'
      WRITE(STDOUT,*) ' '

C*************************************************************C
C..    End write summary of inputs (and check validity!)
C*************************************************************C
      return
      end

c********************************************************************
      subroutine read_inputs(ifile)

c**   Prologue : **
c     
c     read inputs using namelists from piped input files
c
c     last updated 07/13/04 by jaina
c
c***  end prologue **************************************************

      use ioutils, only: stop_execution
      use params_global
      use freewake
      use deflections
      use io_filenames

c*****************  list of global variables used *******************
c      
c     almost all of params_global !!
c
c********************************************************************
      implicit none

      integer ifile
      integer ihar

c***  split input lists here ?
c***********************************
      namelist/inputs/ iread,nmesh,
     & jmax,kmax,lmax,jtail1,ktip,half,nsteps,nrest,npnorm,nmovie,
     & fsmach,alfa,rey,invisc,lamin,ithin,iturb,iunst,ntac,itnmax,dt,
     & timeac,ilhs,idual,iprecon,MP,epse,ibcwp,irhsy,ilim,
     & usegcl,testgcl,totime,angmax,alpha_dot,indicial,beta1,     
     & fmtip,rartio,zground,iwake,ctinp,nblade,
     & nacou,nacoup,nakp1,nakp2,nakp3,nakp4,nakp5,nacouq,
     & jint,kint,lint,
     & ipert,initpe,xvor,zvor,rcore,vorgam,
     & ntot,nrev,
     & imth,dzeta,ivec,iconfine,
     & theta0,theta1c,theta1s,xlam,ideform,irotate,totime0,
     & beta0,beta1c,beta1s,eflap,idefl,idecay,
     & nww,nwwq,def_file,icom,root_co,iteflap,flpdef,imovie,conn_type,
     & split_type, conn_mode, isroot, root, flap_amp, pitch_angle,
     & phase, pitch_knot,plunge_amp, omega,
     & xac, zac,
     & nstart, nstop, ninterval, 
     & kin_rf, kin_freq, kin_phiMax, kin_phiOff, kin_phi_phase, 
     & kin_alphaMax, kin_alphaOff, kin_alpha_phase,
     & kin_beta, kin_phiInit, kin_alphaInit,
     & kin_beta_new, kin_phiOff_new , kin_phiMax_new,
     & xoff, yoff, zoff,
     & isolno, ismth, Ncyc,
     & read_flap_kine, read_pitch_kine,  ! Added by James Lankford
     & xcg, ycg, zcg
c Math: add kinematics

      namelist/flapinp/pxc,flap_ys,flap_ye,dely,dela,angmax_tef,rf_tef,
     &                 nharmflap,ampFlap,phiFlap,theta_f0

      namelist/inpfiles/mesh_inp,overset_inp,mesh_positions,def_file,
     $     tef_inp,movie_inp,soln_dir,rest_dir,palloc_file,
     &     flap_kine_file,pitch_kine_file
c***********************************

c***  first executable statement

      !call set_defaults

      pi = 4.0*atan(1.)

      rewind(ifile)
      read (ifile, inputs)
      rewind(ifile)
      read (ifile, inpfiles)
      if (iteflap.eq.1) then
        rewind(ifile)
        read(ifile,flapinp)
        dela=dela*pi/180
        flap_ys=flap_ys*rartio
        flap_ye=flap_ye*rartio
        dely=dely*rartio
        do ihar=1,nharmflap
           ampFlap(ihar)=ampFlap(ihar)*pi/180.
           phiFlap(ihar)=phiFlap(ihar)*pi/180.
        enddo
        
        !HQP TEF deformation implementation
        if (flpdef .eq. 1) then
	    	call init_flap_params
	    	dely=dely/rartio
		endif
      endif
	
      nd=6
      nv=5

      !Fix for the hover BC nblade
      nblade_tot=nblade
      return
      end

c*********************************************************************
      subroutine init_freewake1

c***  Prologue : **
c     
c     Initialize the freewake geometry by appropriately interpolating
c     the data from 'FWGEOM.dat' which is the output of Maryland Free Wake
c     the subroutine also needs information of the core radius variation
c     (caused by diffusion) with wake age. This information is read in
c     from a file called 'rcore' which is also an output of MFW.
c
c     last updated  07/13/04 by jaina
c
c***  End prologue ***************************************************
      
      use params_global
      use freewake

c*****************  list of global variables used *******************
c
c     pi,rf,dt,alf,rartio,fmtip -> params_global
c     pcx,pcy,pcz,circ,np,nz,ft -> freewake
c
c*********************************************************************

      use ioutils
      implicit none

c..   local variables

      real, allocatable :: pcxt(:,:,:),pcyt(:,:,:),pczt(:,:,:),gv(:,:)
      real, allocatable :: ps(:),ps1(:),rc(:)
      integer i,r,w,p,z
      real xj1,xj2,xjunk,junk,rden,dz,ca,sa
      real tempx,tempy,tempz
      integer un

c***  First executable statement

      un=open_file('FWGEOM.dat',status='old')
      read(un,*) np,nz,nw,nr,ft

      allocate(pcxt(np,nz,nw),pcyt(np,nz,nw),pczt(np,nz,nw),gv(np,nw))
      allocate(ps(nz),rc(nz))
      
      do r=1,nr
         do w=1,nw
            do p=1,np
               do z=1,nz
                  read(un,10) xj1,
     &                 xj2,
     &                 pcxt(p,z,w),
     &                 pcyt(p,z,w),
     &                 pczt(p,z,w),
     &                 xjunk
                  
                  if (z.eq.1) gv(p,w)=xjunk
               end do
               read(un,*)
            end do
            read(un,*)
         end do
         read(un,*)
      enddo
      close(un)
      
 10   format(6(f15.8,1x))

      un=open_file('rcore',status='unknown')
      do i=1,nz
	 read(un,*) junk,rc(i)
      enddo
      close(un)
      
      rden=rc(1)
      do i=1,nz
	 rc(i)=rc(i)*0.05/rden
      enddo


      dz=360./((nz-1)/ft)
      iadim=nint(2*pi/(rf*dt))
      izdim=nint(360./dzeta)*ft+1
      
      allocate(pcx(iadim,izdim,nw),pcy(iadim,izdim,nw),
     &     pcz(iadim,izdim,nw),circ(iadim,nw),rrc(nw,izdim),ps1(izdim))
      
      do z=1,izdim
         ps1(z)=(z-1)*dzeta
      enddo

      do z=1,nz
         ps(z)=(z-1)*dz
      enddo

      do w=1,nw
         call twodinterp(pcxt(1,1,w),pcx(1,1,w),ps,ps1,
     &        nz,np,izdim,iadim)
         call twodinterp(pcyt(1,1,w),pcy(1,1,w),ps,ps1,
     &        nz,np,izdim,iadim)
         call twodinterp(pczt(1,1,w),pcz(1,1,w),ps,ps1,
     &        nz,np,izdim,iadim)
         call fourier_extend(gv(1,w),circ(1,w),np,iadim)
      enddo
      
      call spline(ps,rc,ps1,rrc,nz,izdim)

      np=iadim
      nz=izdim

c..   rotate wake to rotor plane

      ca=cos(-alf*pi/180.)
      sa=sin(-alf*pi/180.)

      print *,ca,sa
      
      do w=1,nw
         do p=1,np
            do z=1,nz
               tempx=-pcy(p,z,w)*rartio
               tempy=(pcx(p,z,w)*ca+pcz(p,z,w)*sa)*rartio
               tempz=(-pcx(p,z,w)*sa+pcz(p,z,w)*ca)*rartio
               
               pcx(p,z,w)=tempx
               pcy(p,z,w)=tempy
               pcz(p,z,w)=tempz
            enddo
         enddo


         do p=1,np   
            circ(p,w)=circ(p,w)*fmtip/(2.*PI)
         enddo
      enddo

      print *,'initialized freewake'
      
      return
      end
 
c..   write out induced velocity at one span station

      subroutine writewg(psi,x,y,z,ug,vg,wg)

      use params_global
      implicit none

      real psi
      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax),vg(jmax,kmax,lmax),wg(jmax,kmax,lmax)

      integer k

 1220 format(26(1X,E10.4))
 1221 format(10(1X,E10.4))

      k=8
      write(101,1221) psi,x(jle,k,1),y(jle,k,1),z(jle,k,1),
     &     x(jtail1,k,1),y(jtail1,k,1),z(jtail1,k,1),ug(jle,k,1),
     &     vg(jle,k,1),wg(jle,k,1)        

      return
      end

c***************************************************************
      subroutine read_kine_def()

c  Subroutine to read in the kinematic deflections
c  Added by James Lankford (03-02-18)

      use ioutils
      use params_global
      use io_filenames

      implicit none

c.. local variables
      integer un,i,io,nlines

c.. first executable statement
      if (read_flap_kine.eq.1) then
      ! read number of lines in data file
      nlines = 0
      un = open_file(flap_kine_file,status='old')
      rewind(un)

      do
        read(un,*,iostat=io)
        if (io.ne.0) exit
        nlines = nlines + 1       
      enddo

      rewind(un)
      allocate(flap_kine_data(nlines))

c.. read in kinematics data
      print *, 'Reading in flap kinematics from file'
	
      do i = 1,nlines
        read(un,*) flap_kine_data(i)
      enddo
	
      close(un)
      endif
	
      if (read_pitch_kine.eq.1) then
      ! read number of lines in data file
      nlines = 0
      un = open_file(pitch_kine_file,status='old')
      rewind(un)
	
      do
        read(un,*,iostat=io)
        if (io.ne.0) exit
        nlines = nlines + 1
      enddo

      rewind(un)
      allocate(pitch_kine_data(nlines))

c.. read in kinematics data
      print *, 'Reading in pitch kinematics from file'

      do i = 1,nlines
        read(un,*) pitch_kine_data(i)
      enddo

      close(un)
      endif

      return
      end

c********************************************************************
      subroutine read_blade_defl(x,y,z)
c***  Prologue : ****
c     
c     Read the deformation data and interpolate it azimuthally and
c     radially for the current resolution. The azimuthal interpolation
c     is performed spectrally (using fourier transforms), while the
c     the radial interpolation is performed using cubic splines.
c
c     the data is read from the 'defs.mod' file
c     the data can be of format in terms of euler parameters (dymore),
c     rotation angles (rcas and camrad) or radial derivatives (umarc)
c
c***  End prologue ****
      use ioutils
      use params_global
      use deflections
      use io_filenames
      use turns_api, only: setup_blade_defl_periodic

c*****************  list of global variables used *******************
c     
c     jmax,kmax,lmax,iazimuth,pi,rf,dt,ktip,theta0,theta1c,theta1s
c     idefl,jle,rartio -> (params_global)
c
c     xyzref,defl_dat,iazimuth -> (deflections)
c
c********************************************************************
      implicit none
      real  x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)
      
c..   local variables
      real, allocatable :: defl_dat1(:,:,:),xyzref1(:,:)
      real, allocatable :: defl_dat2(:,:,:),xyzref2(:,:)
      real, allocatable :: temp(:)
      real, allocatable :: defl(:),defl1(:)
      real, allocatable :: rad(:),rad1(:)
      real, allocatable :: rdef(:),rdef1(:)
      
      real thet,psi,pif
      integer nrhat,npsi
      integer i,j,k,l,kt,kr
      integer un
      
c**** first executable statement
      pif = pi/180.

c..   read in the deflect.dat file
      write(STDOUT,*) 'reading deflections set from ', def_file
      un=open_file(def_file,status='old')

c..   skip comment lines
      if (icom.eq.1) read(un,*)
      if (icom.eq.1) read(un,*)	
      if (icom.eq.1) read(un,*)
      read(un,*) nrhat,npsi
      if (icom.eq.1) read(un,*)

c..   allocate a few things
      iazimuth=nint(2*pi/(rf*dt))

      allocate(defl_dat(6,kmax,iazimuth),xyzref(3,kmax))
      allocate(defl_dat_prev(6,kmax,iazimuth))
      allocate(defl_dat1(6,nrhat,iazimuth))
      allocate(defl_dat2(6,nrhat,npsi),xyzref2(3,nrhat))
      allocate(rad(nrhat),rad1(kmax))
      allocate(defl(npsi),defl1(iazimuth))
      allocate(rdef(nrhat),rdef1(ktip))

c..   begin read
      defl_dat_prev=0.0D0
      read(un,*) xyzref2
      do i=1,npsi
         psi=(i-1)*2*pi/npsi
         thet=theta0+theta1c*cos(psi)+theta1s*sin(psi)

c..	 skip comment lines if necessary
         if (icom.eq.1) then
            read(un,*)
            read(un,*)
            read(un,*)
         endif

         read(un,*) defl_dat2(:,:,i)

         if (idefl.eq.0) then  ! deflections of form u,v,w,vp,wp,phi
            do j=1,nrhat
               defl_dat2(6,j,i)=defl_dat2(6,j,i)+thet
            enddo
         else
		    do j=1,nrhat
			   defl_dat2(4,j,i)=defl_dat2(4,j,i)*pif
			   defl_dat2(5,j,i)=defl_dat2(5,j,i)*pif
			   defl_dat2(6,j,i)=defl_dat2(6,j,i)*pif
		    enddo	
		 endif
      enddo

      close(un)
c..   end read

      call setup_blade_defl_periodic(idefl,xyzref2,defl_dat2,nrhat,
     $     npsi,iazimuth)

Cc..   first interpolate azimuthally
C      
Cc..   now interpolate radially
C
Cc..   interpolate the point of rotation, i.e where the deflections
Cc..   were supplied also radially
C
      return
      end

      subroutine print_credits

      WRITE(STDOUT,*) ' '
      WRITE(STDOUT,*) ' WELCOME TO UM_TURNS_F90  version 1.1'
      WRITE(STDOUT,*) '  This research code should help you with some'
      WRITE(STDOUT,*) '  of your helicopter problems (hopefully :)) '
      WRITE(STDOUT,*)
      WRITE(STDOUT,*) ' This code is the cleaned up fortran90 version with' 
      WRITE(STDOUT,*) ' many additional capabilites developed from '
      WRITE(STDOUT,*) ' TURNS v1.50 which was '
      WRITE(STDOUT,*) ' developed by Srini Srinivasan and cleaned up by'
      WRITE(STDOUT,*) ' J.D. Beader in 1993'
      WRITE(STDOUT,*)
      WRITE(STDOUT,*) '  This code can simulate' 
      WRITE(STDOUT,*) '   1. Steady or Forward flight'
      WRITE(STDOUT,*) '      (deforming meshes, field velocity approach for'
      WRITE(STDOUT,*) '       farwake)'
      WRITE(STDOUT,*)
      WRITE(STDOUT,*) '   2. Wing or 2-D airfoil in steady flight'
      WRITE(STDOUT,*) 
      WRITE(STDOUT,*) ' Supports C-H or C-O mesh topologies'
      WRITE(STDOUT,*) ' single block meshes only supported in this release'
      WRITE(STDOUT,*) ' Overset mesh code is under beta testing'
      WRITE(STDOUT,*)
      WRITE(STDOUT,*) 'Primary Developers:'
      WRITE(STDOUT,*) 'Jay Sitaraman          : jaina@glue.umd.edu'
      WRITE(STDOUT,*) 'Karthikeyan Duraisamy  : dkarthik@glue.umd.edu'
      WRITE(STDOUT,*) 'James Baeder           : baeder@eng.umd.edu'
      WRITE(STDOUT,*) 
      WRITE(STDOUT,*) 'Last updated 10/07/2004'

      end
      
