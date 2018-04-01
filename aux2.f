      subroutine init_stage1(s,q,x,y,z,ug,vg,wg,kkr,kkp,xg,yg,zg,turmu
     $     ,vnu,twist,rest_file)
c
c  initialize the variables to default values, read inputs, check them
c  and set up the working files.
c
c***********************************************************************

      use ioutils
      use params_global
      use refmesh
      use bcparam
      use domain_info, only: restart_solutions, restart_solutions_blocks

      implicit none

      real s(jmax,kmax,lmax,nv)
      real q(jmax,kmax,lmax,nd)
      
      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax),vg(jmax,kmax,lmax),wg(jmax,kmax,lmax)
      real xg(jmax,kmax,lmax),yg(jmax,kmax,lmax),zg(jmax,kmax,lmax)
      real turmu(jmax,kmax,lmax)

      real vnu(jmax,kmax,lmax)
      real twist(kmax)

      integer kkp(kmax),kkr(kmax)

      character(len=FLNLEN) rest_file

c..   local variables

      integer j,k,l,n
      integer jd,kd,ld,md

c..   end local variables

c***  first executable statement

      jd = jmax
      kd = kmax
      ld = lmax
      md = mdim

      jle    = (jmax+1)/2
      jtail2 = jmax - jtail1 + 1
      jbot   = jle - lmax + 1      

      do l=1,lmax
         do k=1,kmax
            do j=1,jmax
               turmu(j,k,l)=0.
               ug(j,k,l)=0.
               vg(j,k,l)=0.
               wg(j,k,l)=0.

               do n=1,nv
                  s(j,k,l,n)=0.
               enddo

            enddo
         enddo
      enddo

      jm = jmax - 1
      km = kmax - 1
      lm = lmax - 1

      do 1 k=1,kmax
        kkp(k) = k + 1
        kkr(k) = k - 1
    1 continue
      kkp(kmax) = kmax +ksym*(1-kmax)
      kkr(1) = 1 + ksym*(kmax-1)
     
c..   setup q and x
      istep0 = 0

c..   copy undeformed grid to base grid
      do l=1,lmax
         do k=1,kmax
            do j=1,jmax
               xg(j,k,l)=x(j,k,l)
               yg(j,k,l)=y(j,k,l)
               zg(j,k,l)=z(j,k,l)
            enddo
         enddo
      enddo

      call qzero( q,vnu)

c..   restart if necessary (get flow information from restart file)
      if( iread .gt. 0) then

         call restart_solutions(q,vnu)
         write(STDOUT,*) 'totime= ', totime
         write(STDOUT,*) 'Restarting at timestep: ', istep0
c$$$        call restr2(q,vnu,rest_file)
c$$$
c$$$	write(13, *) (nsteps-istep0)
c$$$        write(133,*) (nsteps-istep0)
c$$$	!totime=(istep0-360)*rf*dt
c$$$        print *,'totime=',totime
c$$$	print *,'restarting at time step :',istep0
	
      endif

c..   option to coarsen the grid

      if(jint.ne.1 .or. kint.ne.1 .or. lint.ne.1) then
      	call coarse(x,y,z,q,kkp,kkr)
      endif

      return
      end subroutine init_stage1

c***********************************************************************
      subroutine init_stage2(q,x,y,z,xx,xy,xz,yx,yy,yz,zx,zy,zz,
     &     ug,vg,wg,xaj,yaj,zaj,xak,yak,zak,xal,yal,zal,vnaj,vnak,vnal,
     &     svj,svk,svl,wgrid,kkr,kkp,
     &     zx0,zy0,zz0,zt0,xg,yg,zg,xt2,yt2,zt2)
c
c  initialize the variables to default values, read inputs, check them
c  and set up the working files.
c
c***********************************************************************
      use mpi_wrapper
      use params_global
      use refmesh
      use arf_mod
      use pyMOD, only: pythonMode
      use domain_info, only: psi_offset, reshst, this_block
	  
      implicit none
	  
      real q(jmax,kmax,lmax,nd)
      
      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)
      real xx(jmax,kmax,lmax),xy(jmax,kmax,lmax),xz(jmax,kmax,lmax)
      real yx(jmax,kmax,lmax),yy(jmax,kmax,lmax),yz(jmax,kmax,lmax)
      real zx(jmax,kmax,lmax),zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax),vg(jmax,kmax,lmax),wg(jmax,kmax,lmax)
      real xg(jmax,kmax,lmax),yg(jmax,kmax,lmax),zg(jmax,kmax,lmax)
      real xt2(jmax,kmax,lmax),yt2(jmax,kmax,lmax),zt2(jmax,kmax,lmax)
      real wgrid(jmax,kmax,lmax)
	  
      real svj(jmax,kmax,lmax), svk(jmax,kmax,lmax), svl(jmax,kmax,lmax)
      real vnaj(jmax,kmax,lmax), vnak(jmax,kmax,lmax), vnal(jmax,kmax,lmax)
      real xaj(jmax,kmax,lmax), yaj(jmax,kmax,lmax), zaj(jmax,kmax,lmax)
      real xak(jmax,kmax,lmax), yak(jmax,kmax,lmax), zak(jmax,kmax,lmax)
      real xal(jmax,kmax,lmax), yal(jmax,kmax,lmax), zal(jmax,kmax,lmax)
	  
      real zx0(jmax,kmax),zy0(jmax,kmax),zz0(jmax,kmax)
      real zt0(jmax,kmax)
      
      integer kkp(kmax),kkr(kmax)
	  
c..   local variables

      integer jd,kd,ld, imesh
      real srot,psi

      imesh=-100
      srot=rf*totime+psi_offset
      
      psi = srot
      totime0=totime
      
! Read blade deflections only if running standalone mode
      if ((pythonMode .eq. 0) .and.
     $     (ideform .eq. 1)) then
         if (is_wing) then
            call read_blade_defl(x,y,z)
         else
            ! We don't want wake meshes deforming (yet!)
            ideform=NO_DEFORM
         endif
      endif
      
      deform: if (ideform.eq. DEFORMING) then
         if (iteflap.eq.1) then
            if (flpdef .eq. 0) then
               if (flpovh .gt. 0.0) then
                  call rigid_flapoh(xg,yg,zg,srot,0.0,.true.)
               else
                  call rigid_flap(xg,yg,zg,srot,0.0,.true.)
               endif
            else
               call deflect_flap(xg,yg,zg,srot,0.0,.true.)
            endif
         endif

         call init_deform(srot,x,y,z,xx,xy,xz,ug,yx,yy,yz,vg,zx,zy,zz
     $        ,wg,zx0,zy0,zz0,zt0,xg,yg,zg,xt2,yt2,zt2)
         
c..if things are fine lets do the new gcl consistent way
c..compute the metrics and the jacobians using finite volume formula

        if(usegcl)then
            call metgcl(1,kkp,kkr,x,y,z,xt2,yt2,zt2,svj,svk,svl,
     <                xx,xy,xz,yx,yy,yz,zx,zy,zz,
     <                xaj,yaj,zaj,xak,yak,zak,xal,yal,zal,
     <                q,wgrid,vnaj,vnak,vnal,imesh)
        else
            call metfv(q,x,y,z,xx,xy,xz,yx,yy,yz,zx,zy,zz,wgrid,
     &               kkp,kkr,imesh)
        end if

c.. Grid velocities are not computed in init_deform. So let's intialize
c.. using rotational velocities only.
        call mett( x,y,z,xx,xy,xz,yx,yy,yz,zx,zy,zz, ug,vg,wg)
        
      endif deform

      if ((ipert.eq.-3) .and. pythonMode == 0) then
         if (imth.eq.2.or.imth.eq.3) call init_freewake1
      endif

c..    rotate the grid to starting position, if not already done by deform

      ! Math: add kinematics
      if (ideform.ne. DEFORMING .and.
     &   (iunst == FL_UNSTEADY) .or. (iunst == FL_HOVER)) then
         ! We might still want to deflect the flap
         if ((iteflap.eq.1) .and. is_wing) then
            if (flpdef .eq. 0) then
               if (flpovh .gt. 0.0) then
                  call rigid_flapoh(xg,yg,zg,srot,0.0,.true.)
               else
                  call rigid_flap(xg,yg,zg,srot,0.0,.true.)
               endif
            else
               call deflect_flap(xg,yg,zg,srot,0.0,.true.)
            endif
         endif
		
         if(usegcl)then
            call rotate_grid(srot-rf*dt,x,y,z,xt2,yt2,zt2)
            call rotate_grid(rf*dt,x,y,z,xt2,yt2,zt2)
            call metgcl(1,kkp,kkr,x,y,z,xt2,yt2,zt2,svj,svk,svl,
     <                xx,xy,xz,yx,yy,yz,zx,zy,zz,
     <                xaj,yaj,zaj,xak,yak,zak,xal,yal,zal,
     <                q,wgrid,vnaj,vnak,vnal,imesh)
         else
            ! Initialize space metrics...
            call metfv( q,x,y,z,xx,xy,xz,yx,yy,yz,zx,zy,zz,wgrid,
     &           kkp,kkr,imesh)
            ! Rotate to starting position...
            call rotate(srot,x,y,z,xx,xy,xz,ug,yx,yy,yz,vg,zx,zy,zz,wg,
     &           zx0,zy0,zz0,zt0)
         end if
		
         ! Compute grid velocities...
         call mett( x,y,z,xx,xy,xz,yx,yy,yz,zx,zy,zz,ug,vg,wg)
		
      ! Math: add kinematics
      elseif (iunst == FL_PITCH .or. iunst == FL_FLAP) then
		 
		! if statement added to read in kinematics by James 03/02/2018
         if (read_flap_kine.eq.1 .or. read_pitch_kine.eq.1) then
           call read_kine_def(read_flap_kine,read_pitch_kine)
         endif
		
         if (is_wing) then
           if (iunst == FL_PITCH) call pitch_plunge_wing(x,y,z,xg,yg,zg,
     &                            ug,vg,wg,zx,zy,zz,zx0,zy0,zz0,zt0,1)
           if (iunst == FL_FLAP)  call flap_pitch_wing(x,y,z,xg,yg,zg,
     &                            ug,vg,wg,zx,zy,zz,zx0,zy0,zz0,zt0,1)
         endif
         call metfv( q,x,y,z,xx,xy,xz,yx,yy,yz,zx,zy,zz,wgrid,
     &               kkp,kkr,imesh)
		
      endif

c..   time
      if(iunst.eq. FL_UNSTEADY)then
          call tim1(x, y, z, jd, kd, ld)
      end if

c..   divide q by the jacobian, q/q6 
      call qdivj( q)

c..    Add grid velocities due to frame motion, if maneuvering flight
c..    and not using source terms.

      if(arf_opt .eq. 2)then
        call add_arf_gridvel(x,y,z,ug,vg,wg)
      end if

c..field velocities with wake if any

      if ((ipert.eq.-3) .and. is_wing) then
         call pert3(q,x,y,z,ug,vg,wg,zt0,zx,zy,zz,psi)
      endif

c..   calculate dt if cnbr given or else calculate cnbrmax
      call eigen( q,xx,xy,xz,yx,yy,yz,zx,zy,zz,ug,vg,wg)

c..   Open some files for output
      reshst=open_file(trim(adjustl(this_block%block_name))//'.7')
      
      return
      end subroutine init_stage2

c**********************************************************************
      subroutine rotation_term(srot,q,qtn)
c
c**********************************************************************
      use params_global
      implicit none

      real srot,q(jmax,kmax,lmax,nd),qtn(jmax,kmax,lmax,nd)
      integer j,k,l
      real qtmp2,qtmp3,cs,ss

      cs = cos(srot)
      ss = sin(srot)

      do j=1,jmax
      do k=1,kmax
      do l=1,lmax
        qtmp2 = q(j,k,l,2)*cs - q(j,k,l,3)*ss
        qtmp3 = q(j,k,l,3)*cs + q(j,k,l,2)*ss
        q(j,k,l,2) = qtmp2
        q(j,k,l,3) = qtmp3
        qtmp2 = qtn(j,k,l,2)*cs - qtn(j,k,l,3)*ss
        qtmp3 = qtn(j,k,l,3)*cs + qtn(j,k,l,2)*ss
        qtn(j,k,l,2) = qtmp2
        qtn(j,k,l,3) = qtmp3
      enddo
      enddo
      enddo

      end subroutine rotation_term

c**********************************************************************

      subroutine fix_iblanks(iblank,grbc)

      use params_global
      use bcparam
      implicit none

      integer iblank(jmax,kmax,lmax)
      type(bc_t) :: grbc

      integer js,je,ks,ke,ls,le,ib,idir,iproc
      integer j, k, l

      do ib=1,grbc%nbc
        js = grbc%jbcs(ib)
        je = grbc%jbce(ib)
        ks = grbc%kbcs(ib)
        ke = grbc%kbce(ib)
        ls = grbc%lbcs(ib)
        le = grbc%lbce(ib)
        if(js.lt.0) js = jmax+js+1
        if(ks.lt.0) ks = kmax+ks+1
        if(ls.lt.0) ls = lmax+ls+1
        if(je.lt.0) je = jmax+je+1
        if(ke.lt.0) ke = kmax+ke+1
        if(le.lt.0) le = lmax+le+1
        do j = js,je
        do k = ks,ke
        do l = ls,le
           if (iblank(j,k,l) .eq. 1) iblank(j,k,l) = -1
        enddo
        enddo
        enddo
      enddo

      return
      end
