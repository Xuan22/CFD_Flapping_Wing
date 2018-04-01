!
!                 OVERTURNS TIGHT-COUPLING INTERFACE
!
! DESCRIPTION:
! Alternative OVERTURNS interface file for time-accurate tight-coupling
! between flow solver and comprehensive rotor model. Primarily this 
! interface provides a means of exchanging blade deformations and 
! airloads every sub-iteration of the flow solver.
!
! DEVELOPER(S):
!   Benjamin Silbaugh
!
! LAST MODIFIED:
!   Aug 8, 2007 (file created)
!
!==========================================================================

c$$$      subroutine new_iter(dpsi)
c$$$      ! DESC: Call at the begining of a new sub-iteration sequence
c$$$      ! INPUT:
c$$$      !   dpsi... azimuthal increment (deg) for sub iteration sequence
c$$$      !           to follow
c$$$      
c$$$      use params_global
c$$$      use work
c$$$      use rootVariables
c$$$      use tightcouple
c$$$      use arf_mod
c$$$      use dualtime
c$$$
c$$$      implicit none
c$$$
c$$$      integer ntac_loc
c$$$      real dpsi
c$$$
c$$$      ! Be sure that step() routine calls newton()
c$$$      if(itnmax .eq. 1) itnmax = 2
c$$$
c$$$      ntac_loc = 1
c$$$      if((ntac.eq.2).and.(istep.gt.1)) ntac_loc = 2
c$$$
c$$$      ! Update time step
c$$$      dt = dpsi*pi/rf/180.
c$$$
c$$$      ! Reset sub-iteration index
c$$$      itn = 0
c$$$
c$$$      ! *** Execute final initialization procedure if not already done so
c$$$      !
c$$$      ! This includes, (see intia_py2, in aux1.f):
c$$$      !   0) srot, psi and totime0 are set to initial/restart values
c$$$      !   1) setting initial grid deformations (be sure initial grid 
c$$$      !      defs have been set if this is a restart)
c$$$      !   2) Computation of space metrics using finite volume routine
c$$$      !   3) Computation of time metrics due to grid rotation only
c$$$      !   4) Freewake initialization if performing wake-coupling
c$$$      !   5) Tranformation of convservative var to computational var
c$$$      !      (i.e. division by jacobian)
c$$$      !   6) Computation of cnbr
c$$$      !   7) Write out of data for acoustics via call to tim1()
c$$$      !   8) Rotation of background grids to starting position
c$$$      !   9) If restart, airloads are computed
c$$$      !  10) If wake coupling, freewake v-field is set
c$$$      !  11) Initial solution is stored in restart files
c$$$      !  12) Grid connectivity data structures are created
c$$$      !      and initialized
c$$$      !
c$$$      ! Hence, its presence here is to ensure that the initial grid
c$$$      ! configuration (and initialization procedures in which it is a
c$$$      ! dependent) are properly treated. In principle, it should be
c$$$      ! possible to move this to the init_blade_motion() 
c$$$      ! procedure; however, will keep here for now since this is 
c$$$      ! consistent with with the implimentation of run().
c$$$
c$$$      if(turnsInit.eq.0)then
c$$$
c$$$        do im=1,nmesh
c$$$          
c$$$            call set_pointers(ipointer,filenames,
c$$$     &           ig,ig1,ig2,ig21,igq,igs,grid_file,
c$$$     &           soln_file,rest_file,nmesh,im,myid)
c$$$
c$$$            call set_mesh_to_global(jgmx,kgmx,lgmx,
c$$$     &           jtail1mx,jtail2mx,jbotmx,jlemx,ktipmx,krootmx,nmesh,im)
c$$$
c$$$             call initia_py2(s(igs),q(igq),x(ig),y(ig),z(ig),iblank(ig),
c$$$     &        xx(ig),xy(ig),
c$$$     &        xz(ig),yx(ig),yy(ig),yz(ig),zx(ig),zy(ig),zz(ig),
c$$$     &        ug(ig),vg(ig),wg(ig),
c$$$     &        xaj(ig),yaj(ig),zaj(ig),
c$$$     &        xak(ig),yak(ig),zak(ig),
c$$$     &        xal(ig),yal(ig),zal(ig),
c$$$     &        vnaj(ig),vnak(ig),vnal(ig),
c$$$     &        svj(ig),svk(ig),svl(ig),
c$$$     &        z1(ig21),z2(ig21),
c$$$     &        loc1(ig21),loc2(ig21),
c$$$     &        wgrid(ig),kkr(ig1),kkp(ig1),
c$$$     &        zx0(ig2),zy0(ig2),zz0(ig2),zt0(ig2),
c$$$     &        xg(ig),yg(ig),zg(ig),xt2(ig),yt2(ig),zt2(ig),
c$$$     &        turmu(ig),vnu(ig),twist(ig1),
c$$$     &        grid_file,rest_file,im)
c$$$
c$$$            call store(q(igq),vnu(ig),x(ig),y(ig),z(ig),iblank(ig),
c$$$     &        xx(ig),xy(ig),xz(ig),ug(ig),yx(ig),yy(ig),
c$$$     &        yz(ig),vg(ig),zx(ig),zy(ig),zz(ig),wg(ig),
c$$$     &        zx0(ig2),zy0(ig2),zz0(ig2),zt0(ig2),
c$$$     &        kkp(ig1),kkr(ig1),grid_file,soln_file)
c$$$
c$$$
c$$$        end do
c$$$
c$$$         write(6,*) "Initialized q-variables in processor", myid
c$$$
c$$$c..  construct background mesh octree
c$$$
c$$$         if (nmesh.gt.1) then
c$$$            igl=1
c$$$            
c$$$            call set_iblanks_to_one(iblank(igl),jgmx(nmesh),
c$$$     >           kgmx(nmesh),lgmx(nmesh))
c$$$
c$$$c..   make octree for the background mesh
c$$$
c$$$            call make_octree(x(igl),y(igl),z(igl),jgmx(nmesh),
c$$$     >           kgmx(nmesh),lgmx(nmesh))
c$$$
c$$$         end if ! nmesh .gt. 1
c$$$
c$$$c..  do the connectivity once before starting
c$$$
c$$$         psi=rf*totime
c$$$         istep=0
c$$$         
c$$$         if (nmesh.gt.1) then
c$$$            igl=1
c$$$            do im=1,nblade
c$$$               ig=1
c$$$               igd=1
c$$$               
c$$$               call do_connect(x(ig),y(ig),z(ig),
c$$$     <              jgmx(im),kgmx(im),lgmx(im),
c$$$     <              x(igl),y(igl),z(igl),iblank(igl),
c$$$     <              jgmx(nmesh),kgmx(nmesh),lgmx(nmesh),
c$$$     <              imesh(igd),idonor(igd),frac(igd),
c$$$     <              imesh1(igd),idonor1(igd),frac1(igd),ndonor,
c$$$     <              nfringe,dsize3,0.,root_co)
c$$$               
c$$$            enddo
c$$$         end if ! nmesh .gt. 1
c$$$
c$$$         turnsInit=1
c$$$
c$$$      end if ! turnsInit .eq. 0
c$$$
c$$$      ! *** START OF TIME STEP INITIALIZATION ***
c$$$
c$$$      ! Setting params_global data in a maner consistent with run()
c$$$
c$$$
c$$$
c$$$      istep = istep + 1   ! Increment as if in do... loop
c$$$      nrc2 = istep        ! Since executing one step
c$$$      istep0 = istep0 + 1
c$$$
c$$$      srot = rf*dt
c$$$
c$$$      totime = totime + dt
c$$$      psi=rf*(totime-totime0)
c$$$
c$$$      ! Move the background meshes
c$$$
c$$$
c$$$      do im = 1, nmesh
c$$$
c$$$           call set_pointers(ipointer,filenames,
c$$$     &          ig,ig1,ig2,ig21,igq,igs,grid_file,
c$$$     &          soln_file,rest_file,nmesh,im,myid)
c$$$
c$$$           call set_mesh_to_global(jgmx,kgmx,lgmx,
c$$$     &          jtail1mx,jtail2mx,jbotmx,jlemx,ktipmx,krootmx,
c$$$     &          nmesh,im)           
c$$$
c$$$           ! Store previous solution prior to execute of sub-iterations
c$$$           call stqol(q(igq),qnewt(igq),qtn(igq),vnu(ig),vnu0(ig))
c$$$
c$$$c..   rotate
c$$$
c$$$              if (ideform.ne.1.or.im.gt.nblade) then
c$$$                if (iteflap.eq.1.and.im.le.nblade) then
c$$$                  call rotate_rigid_flap(srot,x(ig),y(ig),z(ig),
c$$$     <                xx(ig),xy(ig),xz(ig),ug(ig),yx(ig),
c$$$     <                yy(ig),yz(ig),vg(ig),
c$$$     <                zx(ig),zy(ig),zz(ig),wg(ig),
c$$$     <                zx0(ig2),zy0(ig2),zz0(ig2),zt0(ig2),0)
c$$$
c$$$                  call qmulj(q)
c$$$                  if(usegcl)then
c$$$                    call metgcl(ntac_loc,kkp(ig1),kkr(ig1),
c$$$     <                   x(ig),y(ig),z(ig),xt2(ig),yt2(ig),zt2(ig),
c$$$     <                   svj(ig),svk(ig),svl(ig),
c$$$     <                   xx(ig),xy(ig),xz(ig),yx(ig),yy(ig),yz(ig),
c$$$     <                   zx(ig),zy(ig),zz(ig),xaj(ig),yaj(ig),zaj(ig),
c$$$     <                   xak(ig),yak(ig),zak(ig),xal(ig),yal(ig),zal(ig),
c$$$     <                   q(igq),wgrid(ig),vnaj(ig),vnak(ig),vnal(ig),imesh)
c$$$                  else
c$$$		            call metfv(q(igq),x(ig),y(ig),z(ig),
c$$$     <                   xx(ig),xy(ig),xz(ig),yx(ig),yy(ig),
c$$$     <                   yz(ig),zx(ig),zy(ig),zz(ig),wgrid(ig),
c$$$     <                   kkp(ig1),kkr(ig1),im)
c$$$                  end if
c$$$                  call qdivj(q)
c$$$
c$$$                else
c$$$
c$$$                   if(usegcl)then
c$$$                        call rotate_grid(srot,x(ig),y(ig),z(ig),
c$$$     <                       xt2(ig),yt2(ig),zt2(ig))
c$$$                        call metgcl(ntac_loc,kkp(ig1),kkr(ig1),
c$$$     <                   x(ig),y(ig),z(ig),xt2(ig),yt2(ig),zt2(ig),
c$$$     <                   svj(ig),svk(ig),svl(ig),
c$$$     <                   xx(ig),xy(ig),xz(ig),yx(ig),yy(ig),yz(ig),
c$$$     <                   zx(ig),zy(ig),zz(ig),xaj(ig),yaj(ig),zaj(ig),
c$$$     <                   xak(ig),yak(ig),zak(ig),xal(ig),yal(ig),zal(ig),
c$$$     <                   q(igq),wgrid(ig),vnaj(ig),vnak(ig),vnal(ig),imesh)
c$$$                   else
c$$$                        call rotate(srot,x(ig),y(ig),z(ig),
c$$$     <                       xx(ig),xy(ig),xz(ig),ug(ig),yx(ig),
c$$$     <                       yy(ig),yz(ig),vg(ig),
c$$$     <                       zx(ig),zy(ig),zz(ig),wg(ig),
c$$$     <                       zx0(ig2),zy0(ig2),zz0(ig2),zt0(ig2))
c$$$                   endif
c$$$
c$$$                   call mett(x(ig),y(ig),z(ig),
c$$$     &                  xx(ig),xy(ig),xz(ig),yx(ig),yy(ig),
c$$$     &                  yz(ig),zx(ig),zy(ig),zz(ig),
c$$$     &                  ug(ig),vg(ig),wg(ig))
c$$$
c$$$                endif
c$$$
c$$$                 if(arf_opt .eq. 2)then
c$$$                   ! Using grid velocities insteady of source terms
c$$$                   call add_arf_gridvel( x(ig),  y(ig),  z(ig), 
c$$$     &                                  ug(ig), vg(ig), wg(ig) )
c$$$                 end if
c$$$
c$$$         else
c$$$
c$$$           ! Store grid configuration at time level n-1; needed
c$$$           ! for subsequent sub-iterations with deformations.
c$$$           ! (only need to do this for the blade meshes)
c$$$
c$$$            call store_gridconfig(xt2(ig), yt2(ig), zt2(ig),
c$$$     <            xt2ref(ig), yt2ref(ig), zt2ref(ig),
c$$$     <            svj(ig), svk(ig), svl(ig), 
c$$$     <            svjref(ig), svkref(ig), svlref(ig))
c$$$
c$$$         end if
c$$$
c$$$
c$$$c..   write some outputs
c$$$c-----------------------------------------------------------------------
c$$$        if (mod(istep0,nmovie).eq.0.and.myid.eq.0) then
c$$$           call movie(q(igq),x(ig),y(ig),z(ig),im)
c$$$        endif
c$$$
c$$$        
c$$$!        if( mod(istep0,nrest).eq.0 .or.istep0.eq.2 )
c$$$!     &       call store(q(igq),vnu(ig),x(ig),
c$$$!     &       y(ig),z(ig),iblank(ig),xx(ig),xy(ig),xz(ig),ug(ig),
c$$$!     &       yx(ig),yy(ig),yz(ig),vg(ig),zx(ig),zy(ig),zz(ig),wg(ig),
c$$$!     &       zx0(ig2),zy0(ig2),zz0(ig2),zt0(ig2),
c$$$!     &       kkp(ig1),kkr(ig1),grid_file,soln_file)
c$$$!
c$$$        if( mod(istep,npnorm).eq.0 )
c$$$     &        write(6,101) istep0,rf*totime*180./pi,totime
c$$$
c$$$c-----------------------------------------------------------------------
c$$$c..   end write outputs
c$$$
c$$$        end do !do im=1,nmesh
c$$$
c$$$           write(6,700) ct,cq,figofmerit
c$$$           write(11,710) istep0,totime,ct,cq,figofmerit
c$$$
c$$$ 101  format(/,' istep,azimuth,time =',1x,i5,2(1x,f10.5))
c$$$ 700  format ('  ct,cq,figofmerit=   ',/, 3g12.5)
c$$$ 710  format (i5,4g13.5)
c$$$        
c$$$        call dualtime_reset_sequence()
c$$$
c$$$      end subroutine new_iter
c$$$
c$$$      subroutine exec_subiter()
c$$$!DESC: Executes a single sub-iteration. The sub-iteration sequence
c$$$!      begins with deformation of the body-fitted grid using the
c$$$!      data currently contained in defl_dat array. Then the flow
c$$$!      solver step() procedure is invoked. Lastly, airloads are
c$$$!      computed thus completing the sub-iteration.
c$$$
c$$$      use params_global
c$$$      use work
c$$$      use rootVariables
c$$$      use tightcouple
c$$$      use arf_mod
c$$$      use dualtime
c$$$
c$$$      implicit none
c$$$
c$$$      integer :: ntac_loc
c$$$
c$$$      character*3 :: idstr
c$$$      character*50 :: conlog
c$$$
c$$$      write(idstr, "(I3.3)") myid
c$$$      conlog = "ConvHist_"//idstr//".dat"
c$$$
c$$$      ntac_loc = 1
c$$$      if((ntac.eq.2).and.(istep.gt.1)) ntac_loc = 2
c$$$
c$$$      itn = itn + 1
c$$$
c$$$! *** DEFORM BODY FITTED GRIDS ***
c$$$
c$$$
c$$$
c$$$      do im = 1, nblade
c$$$
c$$$           psi_rot = rf*totime+myid*dblad
c$$$
c$$$           call set_pointers(ipointer,filenames,
c$$$     &          ig,ig1,ig2,ig21,igq,igs,grid_file,
c$$$     &          soln_file,rest_file,nmesh,im,myid)
c$$$
c$$$           call set_mesh_to_global(jgmx,kgmx,lgmx,
c$$$     &          jtail1mx,jtail2mx,jbotmx,jlemx,ktipmx,krootmx,
c$$$     &          nmesh,im)           
c$$$
c$$$
c$$$           if(itn > 1)then
c$$$             ! Re-set reference grids to configuration at prior time 
c$$$             ! level t_n. This is needed for correct time-metrics
c$$$             call restore_gridconfig(x(ig),y(ig),z(ig),xt2(ig),yt2(ig),
c$$$     <             zt2(ig), xt2ref(ig), yt2ref(ig), zt2ref(ig),
c$$$     <             svj(ig), svk(ig), svl(ig), 
c$$$     <             svjref(ig), svkref(ig), svlref(ig))
c$$$           end if
c$$$
c$$$! .. now (re)deform grids
c$$$
c$$$        if (iteflap.eq.1)
c$$$     <     call rigid_flap(xg(ig),yg(ig),zg(ig), psi_rot,psi_rot-srot)
c$$$
c$$$           call deform(srot,psi_rot,x(ig),y(ig),z(ig),
c$$$     <          xx(ig),xy(ig),xz(ig),ug(ig),yx(ig),yy(ig),
c$$$     <          yz(ig),vg(ig),zx(ig),zy(ig),zz(ig),wg(ig),
c$$$     <          zx0(ig2),zy0(ig2),zz0(ig2),zt0(ig2),
c$$$     <          xt2(ig),yt2(ig),zt2(ig),
c$$$     <          xg(ig),yg(ig),zg(ig))
c$$$
c$$$           call qmulj(q)
c$$$           if(usegcl)then
c$$$                call metgcl(ntac_loc,kkp(ig1),kkr(ig1),
c$$$     <               x(ig),y(ig),z(ig),xt2(ig),yt2(ig),zt2(ig),
c$$$     <               svj(ig),svk(ig),svl(ig),
c$$$     <               xx(ig),xy(ig),xz(ig),yx(ig),yy(ig),yz(ig),
c$$$     <               zx(ig),zy(ig),zz(ig),xaj(ig),yaj(ig),zaj(ig),
c$$$     <               xak(ig),yak(ig),zak(ig),xal(ig),yal(ig),zal(ig),
c$$$     <               q(igq),wgrid(ig),vnaj(ig),vnak(ig),vnal(ig),imesh)
c$$$           else
c$$$		        call metfv(q(igq),x(ig),y(ig),z(ig),
c$$$     <               xx(ig),xy(ig),xz(ig),yx(ig),yy(ig),
c$$$     <               yz(ig),zx(ig),zy(ig),zz(ig),wgrid(ig),
c$$$     <               kkp(ig1),kkr(ig1),im)
c$$$           end if
c$$$           call qdivj(q)
c$$$                 
c$$$        if (arf_opt .eq. 2) then
c$$$          ! Using grid velocities instead of source terms
c$$$          call add_arf_gridvel(  x(ig),  y(ig),  z(ig), 
c$$$     &                          ug(ig), vg(ig), wg(ig) )
c$$$        end if
c$$$
c$$$      end do ! im=1,nblade
c$$$
c$$$!  *** UPDATE GRID CONNECTIVITY ***
c$$$
c$$$      call cpu_time(t1)
c$$$      if (nmesh.gt.1) then
c$$$         igl=1
c$$$         call set_iblanks_to_one(iblank(igl),jgmx(nmesh),kgmx(nmesh),
c$$$     <        lgmx(nmesh))
c$$$         do im=1,nblade
c$$$            ig=1
c$$$            igd=1
c$$$            
c$$$            call do_connect(x(ig),y(ig),z(ig),
c$$$     <           jgmx(im),kgmx(im),lgmx(im),
c$$$     <           x(igl),y(igl),z(igl),iblank(igl),
c$$$     <           jgmx(nmesh),kgmx(nmesh),lgmx(nmesh),
c$$$     <           imesh(igd),idonor(igd),frac(igd),
c$$$     <           imesh1(igd),idonor1(igd),frac1(igd),ndonor,
c$$$     <           nfringe,dsize3,psi,root_co)
c$$$            
c$$$         enddo
c$$$      endif
c$$$      call cpu_time(t2)
c$$$
c$$$      print *, 'connectivity time', t2-t1
c$$$
c$$$! *** EXECUTE FLOW SOLVER SUB-ITERATION ***
c$$$
c$$$      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c$$$
c$$$      do im=1,nmesh
c$$$
c$$$              call set_pointers(ipointer,filenames,
c$$$     &             ig,ig1,ig2,ig21,igq,igs,grid_file,
c$$$     &             soln_file,rest_file,nmesh,im,myid)
c$$$
c$$$              call set_mesh_to_global(jgmx,kgmx,lgmx,
c$$$     &             jtail1mx,jtail2mx,jbotmx,jlemx,ktipmx,krootmx,
c$$$     &             nmesh,im)
c$$$
c$$$              call step(x(ig),y(ig),z(ig),iblank(ig),q(igq),qnewt(igq),
c$$$     <          s(igs),turmu(ig),vnu(ig),vnu0(ig),xx(ig),xy(ig),xz(ig),
c$$$     <          ug(ig),yx(ig),yy(ig),yz(ig),vg(ig),
c$$$     <          zx(ig),zy(ig),zz(ig),wg(ig),
c$$$     <          zx0(ig2),zy0(ig2),zz0(ig2),zt0(ig2),
c$$$     <          xaj(ig),yaj(ig),zaj(ig),xak(ig),yak(ig),zak(ig),
c$$$     <          xal(ig),yal(ig),zal(ig),vnaj(ig),vnak(ig),vnal(ig),
c$$$     <          kkp(ig1),kkr(ig1),wgrid(ig),wvec,
c$$$     <          resmax,resrho,rsum,im)
c$$$
c$$$              open(unit=200, file=conlog, position='append')
c$$$              write(200, *) im, totime, itn, rsum
c$$$              close(200)
c$$$
c$$$              call bc(x(ig),y(ig),z(ig),q(igq),
c$$$     <             xx(ig),xy(ig),xz(ig),
c$$$     <             ug(ig),yx(ig),yy(ig),yz(ig),vg(ig),
c$$$     <             zx(ig),zy(ig),zz(ig),wg(ig),
c$$$     <             zx0(ig2),zy0(ig2),zz0(ig2),zt0(ig2),
c$$$     <             kkr(ig1),kkp(ig1),ireq,mpistatus,im,0)
c$$$
c$$$      end do ! nmesh = 1, nmesh
c$$$
c$$$!  *** PERFORM CHIMERA INTERPOLATIONS ***
c$$$
c$$$      print*, "Performing chimera interpolations, proc = ", myid
c$$$
c$$$      if (nmesh .gt. 1) then
c$$$        igl=1
c$$$        igql=1
c$$$
c$$$        do im = 1, nblade
c$$$
c$$$           ig=1
c$$$           igq=1
c$$$           igd=1
c$$$
c$$$           ! background mesh to blade mesh
c$$$
c$$$           call do_interpolations(q(igq),vnu(ig),
c$$$     <          jgmx(im),kgmx(im),lgmx(im),
c$$$     <          q(igql),vnu(igl),
c$$$     <          jgmx(nmesh),kgmx(nmesh),lgmx(nmesh),
c$$$     <          imesh(igd),idonor(igd),frac(igd),
c$$$     <          ndonor,dsize3)
c$$$
c$$$           ! blade mesh to background mesh
c$$$                 
c$$$            call do_interpolations(q(igql),vnu(igl),
c$$$     <           jgmx(nmesh),kgmx(nmesh),lgmx(nmesh),
c$$$     <           q(igq),vnu(ig),jgmx(im),kgmx(im),lgmx(im),
c$$$     <           imesh1(igd),idonor1(igd),frac1(igd),
c$$$     <           nfringe,dsize3)
c$$$
c$$$        end do
c$$$
c$$$      end if
c$$$
c$$$      ! Call boundary conditions again
c$$$
c$$$      mstop=0
c$$$
c$$$      do im = 1, nmesh
c$$$
c$$$              call set_pointers(ipointer,filenames,
c$$$     &             ig,ig1,ig2,ig21,igq,igs,grid_file,
c$$$     &             soln_file,rest_file,nmesh,im,myid)
c$$$
c$$$              call set_mesh_to_global(jgmx,kgmx,lgmx,
c$$$     &             jtail1mx,jtail2mx,jbotmx,jlemx,ktipmx,krootmx,
c$$$     &             nmesh,im)
c$$$
c$$$              call bc(x(ig),y(ig),z(ig),q(igq),
c$$$     <             xx(ig),xy(ig),xz(ig),
c$$$     <             ug(ig),yx(ig),yy(ig),yz(ig),vg(ig),
c$$$     <             zx(ig),zy(ig),zz(ig),wg(ig),
c$$$     <             zx0(ig2),zy0(ig2),zz0(ig2),zt0(ig2),
c$$$     <             kkr(ig1),kkp(ig1),ireq,mpistatus,im,1)
c$$$
c$$$              call monitor(mstopTmp,q(igq))
c$$$
c$$$              mstopNode=mstopNode+mstopTmp
c$$$
c$$$      end do ! im = 1, nmesh
c$$$
c$$$           call MPI_REDUCE(mstopNode,mstop,1,MPI_INTEGER,MPI_SUM,0,
c$$$     $          MPI_COMM_WORLD,ierr)
c$$$
c$$$           if(mstop .gt. 0 ) stop
c$$$
c$$$	   call flush(6)
c$$$
c$$$! *** COMPUTE AIRLOADS ***
c$$$
c$$$      do im = 1, nblade
c$$$
c$$$           call set_pointers(ipointer,filenames,
c$$$     &          ig,ig1,ig2,ig21,igq,igs,grid_file,
c$$$     &          soln_file,rest_file,nmesh,im,myid)
c$$$
c$$$           call set_mesh_to_global(jgmx,kgmx,lgmx,
c$$$     &          jtail1mx,jtail2mx,jbotmx,jlemx,
c$$$     &          ktipmx,krootmx,nmesh,im)
c$$$
c$$$           ! These are artifacts from older version?
c$$$           ! Will keep for now... (BS)
c$$$           x0 = 0.
c$$$           y0 = 0.
c$$$           z0 = 0.
c$$$
c$$$           call force2d(1,x0,y0,z0,
c$$$     <          ct,cq,figofmerit,
c$$$     <          x(ig),y(ig),z(ig),q(igq),xx(ig),xy(ig),xz(ig),
c$$$     <          zx(ig),zy(ig),zz(ig),twist,128+myid,2)
c$$$
c$$$
c$$$      end do
c$$$
c$$$      call dualtime_next_dtpseudo()
c$$$
c$$$      call mpi_barrier(mpi_comm_world,ierr)
c$$$
c$$$      end subroutine exec_subiter
c$$$
c$$$      subroutine write_restart()
c$$$!DESC: Writes restart files
c$$$
c$$$      use params_global
c$$$      use refmesh
c$$$      use work
c$$$      use rootVariables
c$$$
c$$$      implicit none
c$$$
c$$$      print*, "WRITTING RESTART FILES, proc = ", myid
c$$$
c$$$      do im=1,nmesh
c$$$         
c$$$         call set_pointers(ipointer,filenames,
c$$$     &        ig,ig1,ig2,ig21,igq,igs,grid_file,
c$$$     &        soln_file,rest_file,nmesh,im,myid)
c$$$         
c$$$         call set_mesh_to_global(jgmx,kgmx,lgmx,
c$$$     &        jtail1mx,jtail2mx,jbotmx,jlemx,
c$$$     &        ktipmx,krootmx,nmesh,im)
c$$$         
c$$$         call store(q(igq),vnu(ig),x(ig),y(ig),z(ig),iblank(ig),
c$$$     &        xx(ig),xy(ig),xz(ig),ug(ig),yx(ig),yy(ig),
c$$$     &        yz(ig),vg(ig),zx(ig),zy(ig),zz(ig),wg(ig),
c$$$     &        zx0(ig2),zy0(ig2),zz0(ig2),zt0(ig2),
c$$$     &        kkp(ig1),kkr(ig1),grid_file,soln_file)
c$$$
c$$$
c$$$      enddo
c$$$
c$$$      end subroutine write_restart

      subroutine store_gridconfig(xt2, yt2, zt2, xt2ref, yt2ref, zt2ref,
     <           svj, svk, svl, svjref, svkref, svlref)

      use params_global

      real, dimension(jmax, kmax, lmax) ::  xt2, yt2, zt2
      real, dimension(jmax, kmax, lmax) ::  xt2ref, yt2ref, zt2ref
      real, dimension(jmax, kmax, lmax) ::  svj, svk, svl
      real, dimension(jmax, kmax, lmax) ::  svjref, svkref, svlref

      !local var
      integer j, k, l

      do l = 1, lmax
        do k = 1, kmax
          do j = 1, jmax

            xt2ref(j,k,l) = xt2(j,k,l)
            yt2ref(j,k,l) = yt2(j,k,l)
            zt2ref(j,k,l) = zt2(j,k,l)

            svjref(j,k,l) = svj(j,k,l)
            svkref(j,k,l) = svk(j,k,l)
            svlref(j,k,l) = svl(j,k,l)

          end do
        end do
      end do

      end subroutine store_gridconfig

      subroutine restore_gridconfig(x, y, z, xt2, yt2, zt2, xt2ref,
     < yt2ref, zt2ref, svj, svk, svl, svjref, svkref, svlref)

      use params_global

      implicit none

      real, dimension(jmax, kmax, lmax) ::  x, y, z
      real, dimension(jmax, kmax, lmax) ::  xt2, yt2, zt2
      real, dimension(jmax, kmax, lmax) ::  xt2ref, yt2ref, zt2ref
      real, dimension(jmax, kmax, lmax) ::  svj, svk, svl
      real, dimension(jmax, kmax, lmax) ::  svjref, svkref, svlref

      !local var
      integer j, k, l

      do l = 1, lmax
        do k = 1, kmax
          do j = 1, jmax

             x(j,k,l) = xt2(j,k,l)
             y(j,k,l) = yt2(j,k,l)
             z(j,k,l) = zt2(j,k,l)

             xt2(j,k,l) = xt2ref(j,k,l)
             yt2(j,k,l) = yt2ref(j,k,l)
             zt2(j,k,l) = zt2ref(j,k,l)

             svj(j,k,l) = svjref(j,k,l)
             svk(j,k,l) = svkref(j,k,l)
             svl(j,k,l) = svlref(j,k,l)

          end do
        end do
      end do

      end subroutine restore_gridconfig
