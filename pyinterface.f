c..   pyinterface routines
c...  copyright jay sitaraman
c..   last modified 02/08/05


      subroutine execute(nsteps)

      use TURNS_options
      use params_global, only: is_wing
      use pyMOD
      use turns_api, only: run
      implicit none
      
      integer:: nsteps

cf2py intent(in) :: nsteps

      if (is_wing) then 
         if (allocated(loads)) deallocate(loads)
         allocate(loads(7,nsa,nsteps))
         write(6,*) 'Allocated array for collecting airloads'
      endif
      
      call run(nsteps)
      
      return
      end

      
      subroutine startup(inp_file,procid)

      use pyMOD
      use mpi_wrapper, only: PID
      use turns_api, only: init
      implicit none

      character*128 inp_file
      integer procid
      
cf2py intent(in) :: inp_file
cf2py intent(out) :: procid
      
      pythonMode=1
      call init(inp_file)
      procid=PID

      return
      end


c********************************************************************
      subroutine set_blade_defl(ideff,xyzref2,defl_dat2,nrhat,npsi)

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

      use params_global
      use deflections
      use pyMOD
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
      
      integer ideff,nrhat,npsi
      real xyzref2(3,nrhat),defl_dat2(6,nrhat,npsi)

cf2py intent(in) :: nrhat,npsi,xyzref2,defl_dat2

      
c..   local variables


      real, allocatable:: defl_dat1(:,:,:)
      real, allocatable :: defl(:),defl1(:)
      real, allocatable :: rad(:),rad1(:)
      real, allocatable :: rdef(:),rdef1(:)
      
      real thet,pif
      integer i,j,k,l,kt,kr
      real psi
      
c**** first executable statement

      pif = pi/180.

      if (.not. is_wing) then
         ideform=0
         return
      endif

      !call set_mesh_to_global

c..   allocate a few things
      
      iazimuth=nint(2*pi/(rf*dt))
      write(6,*) pi,rf,dt,pif
      write(6,*) 'dpsi=',rf*dt/pif
      write(6,*) 'iazimuth=',iazimuth
      write(6,*) 'def_alloc_flag=',def_alloc_flag
      write(6,*) 'nsa=',nsa
      write(6,*) 'kmax=',kmax
      if (def_alloc_flag.eq.0) then
       write(6,*) 'allocating arrays for deflections'
       allocate(defl_dat(6,kmax,iazimuth),xyzref(3,kmax))
       allocate(defl_dat_prev(6,kmax,iazimuth))
       defl_dat_prev=0.0D0
       def_alloc_flag=1
      endif

      call setup_blade_defl_periodic(ideff,xyzref2,defl_dat2,nrhat,
     $     npsi,iazimuth)
      ideform=1
      idefl=ideff
      return

      allocate(defl_dat1(6,nrhat,iazimuth))
      allocate(rad(nrhat),rad1(kmax))
      allocate(defl(npsi),defl1(iazimuth))
      allocate(rdef(nrhat),rdef1(ktip))

c..   modify blade motions appropriately

      do i=1,npsi
         psi=(i-1)*2*pi/npsi
         thet=theta0+theta1c*cos(psi)+theta1s*sin(psi)

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


c..   first interpolate azimuthally

      do l=1,6
         do k=1,nrhat

            do i=1,npsi
               defl(i)=defl_dat2(l,k,i)
            enddo
            
c..   calling fourier extension, spectral interpolation

            call fourier_extend(defl,defl1,npsi,iazimuth)
            
            do i=1,iazimuth
               defl_dat1(l,k,i)=defl1(i)
            enddo
            
         enddo
      enddo


c..   now interpolate radially

      ! set spanwise coords

      do k=1,kmax
         rad1(k)=spandist(k)
      enddo
      
      do k=1,nrhat
         rad(k)=xyzref2(1,k)
      enddo

      ! small fix to prevent the bad metrics at the tip

      if (ktip.eq.kmax) then
         kt=kmax-7
      else
         kt=ktip
      endif

      if (root_co.eq.1) then
	kr=5
      else
	kr=kroot
      endif
	

      do l=1,6
         do i=1,iazimuth

            do k=1,nrhat
               rdef(k)=defl_dat1(l,k,i)
            enddo
            
            call interp1(rad,rdef,rad1,rdef1,nrhat,kt)
            
            !call spline(rad,rdef,rad1,rdef1,nrhat,kt)
            
            do k=1,kt
               defl_dat(l,k,i)=rdef1(k)
            enddo

c..   setting all outside tip to be equal to the tip

            do k=kt+1,kmax
               defl_dat(l,k,i)=defl_dat(l,kt,i)
            enddo
	 
            do k=1,kr-1
	      defl_dat(l,k,i)=defl_dat(l,kr,i)
	    enddo

         enddo
      enddo

c..   interpolate the point of rotation, i.e where the deflections
c..   were supplied also radially

      ! x-coord same as grid

      do k=1,kmax
         xyzref(1,k)=rad1(k)*rartio
      enddo

      ! y,z coords need interpolation

      do i=2,3
         do k=1,nrhat
            rdef(k)=xyzref2(i,k)
         enddo
c         call interp1(rad,rdef,rad1,rdef1,nrhat,kt)
         call spline(rad,rdef,rad1,rdef1,nrhat,kt)

         do k=1,kt
            xyzref(i,k)=rdef1(k)*rartio
         enddo
      enddo
      
      do k=kt+1,kmax
         xyzref(1,k)=xyzref(1,kt)
         xyzref(2,k)=xyzref(2,kt)
         xyzref(3,k)=xyzref(3,kt)
      enddo

      do k=1,kr-1
         xyzref(1,k)=xyzref(1,kr)
         xyzref(2,k)=xyzref(2,kr)
         xyzref(3,k)=xyzref(3,kr)
      enddo

      
!      write(20,*) kmax,iazimuth
!      do k=1,kmax
!      write(20,*) xyzref(1,k),xyzref(2,k),xyzref(3,k)
!      enddo
!      do i=1,iazimuth
!	do k=1,kmax
!	 write(20,190) (defl_dat(l,k,i),l=1,6)
!	enddo	
!      enddo

190   format(6(1X,E14.8))

      deallocate(defl_dat1,defl,defl1,rad,rad1,rdef,rdef1)
      write(6,*) 'set blade motions in turns'
      ideform=1
      idefl=ideff
      return
      end


      subroutine setFrstrm(dpsi,fsm,fmtt,alp)

      use params_global

      implicit none

      real dpsi,fsm,alp,fmtt, cs, ss

cf2py intent (in) :: dpsi,fsm,alp,fmtt

      fmtip=fmtt
      fsmach=fsm
      alf=alp
      
      cs=cos(alf*pi/180.)
      ss=sin(alf*pi/180.)

      if(iunst.eq. FL_STEADY .or.fmtip.eq.0) then
        uinf = fsmach*cs
        vinf = uinf*tan(beta1*pi/180.)
      else
        uinf  = 0.0
        vinf  = fsmach*cs
      endif

      winf  = fsmach*ss
      
      if (fmtip.gt.0) then
         rf = fmtip/rartio
      else
         rf=2*fsmach*rf
      endif

      dt = dpsi*pi/rf/180.

      return
      end

      subroutine get_mesh_span_info(mname,mid,nsa)

      use domain_info, only: meshes,get_global_mesh_id
      use airloads, only: get_num_span_loc
      
      character*128 :: mname
      integer :: mid, nsa

cf2py intent(in) :: mname
cf2py intent(out) :: mid, nsa      
      
      mid=get_global_mesh_id(meshes,mname)
      nsa=get_num_span_loc(mid)

      end

      subroutine get_airloads(mid,nstps,nsa,rdist,fmom)

      use airloads, only: send_airloads
      
      integer, intent(in) :: mid,nstps,nsa
      real, dimension(nsa), intent(out) :: rdist
      real, dimension(7,nsa,nstps), intent(out) :: fmom

      call send_airloads(mid,nstps,rdist,fmom)

      end

      subroutine setWake(pcxt,pcyt,pczt,gv,rc,np1,nz1,nw1)

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

      implicit none

c..   local variables
      integer np1,nz1,nw1

      real :: pcxt(np1,nz1,nw1),pcyt(np1,nz1,nw1),pczt(np1,nz1,nw1)
      real :: gv(np1,nw1)
      real :: rc(nw1,nz1)

      real, allocatable :: ps(:),ps1(:)
      integer i,r,w,p,z
      real xj1,xj2,xjunk,junk,rden,dz,ca,sa,dp
      real tempx,tempy,tempz

c***  First executable statement

      if (.not. is_wing) then
         ipert=0
         return
      endif
      
      np=np1
      nz=nz1
      nw=nw1
      
      ft=(nz-1)/np
      dz=360./((nz-1)/ft)
      iadim=nint(2*pi/(rf*dt))
      dzeta=dz
      izdim=nint(360./dzeta)*ft+1

      if (fw_alloc.eq.0) then
         allocate(pcx(iadim,izdim,nw),pcy(iadim,izdim,nw),
     &        pcz(iadim,izdim,nw),circ(iadim,nw),rrc(nw,izdim))
         fw_alloc=1
      endif

      allocate(ps(nz),ps1(izdim))
      
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
     
      do w=1,nw 
       call spline(ps,rc(w,:),ps1,rrc(w,:),nz,izdim)
      enddo
   
      np=iadim
      nz=izdim

c..   rotate wake to rotor plane

c      ca=cos(-alf*pi/180.)
c      sa=sin(-alf*pi/180.)

      ca=1.0
      sa=0.

!      print *,ca,sa
      
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

      !call print_message('initialized freewake')

c      if (procid .eq. MASTER) then 
c      open(unit=86,file='FWrefine.dat',form='formatted')
c      write(86,*) np,nz,nw,1,ft
c
c      dp=360./np
c
c      do r=1,1
c         do w=1,nw
c            do p=1,np
c               do z=1,nz
c                  write(86,10) dp,
c     &                   dp,
c     &                   pcx(p,z,w),
c     &                   pcy(p,z,w),
c     &                   pcz(p,z,w),
c     &                   circ(p,w)
c               end do
c               write(86,*)
c            end do
c            write(86,*)
c         end do
c         write(86,*)
c      end do
c
c 10          format(6(f15.8,1x))
c      close(86)
c      endif
c      
      deallocate(ps)
      deallocate(ps1)

      ipert=-3
      
      return
      end
      
      subroutine set_tef_angle(dlt0,dlt1c,dlt1s)

      use params_global

      implicit none

      real :: dlt0, dlt1c, dlt1s

cf2py intent(in) :: dlt0, dlt1c, dlt1s

      theta_f0=dlt0

      end 
