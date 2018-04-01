
c***********************************************************************
      subroutine bc(x,y,z,q,vnut,xx,xy,xz,ug,yx,yy,yz,vg,zx,zy,zz,wg,
     <     bt,zx0,zy0,zz0,zt0,kkr,kkp,grbc)

c** prologue : ****
c  control function to update the mesh boundaries
c  calls wall bc, characteristic extrapolation, averaging at wake-cuts
c  for c-o, c-h mesh types
c
c  last update 07/13/04 by jaina
c***********************************************************************
      
      use params_global
      use bcparam

c*******************  list of global variables used ********************
c
c    jmax, kmax, lmax, nd, ktip
c    jtail1,jtail2,kroot,half,iunst,iwake
c
c***********************************************************************

      implicit none

      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)
      real q(jmax,kmax,lmax,nd),vnut(jmax,kmax,lmax)
      real xx(jmax,kmax,lmax),xy(jmax,kmax,lmax),xz(jmax,kmax,lmax),
     <     yx(jmax,kmax,lmax),yy(jmax,kmax,lmax),yz(jmax,kmax,lmax),
     <     zx(jmax,kmax,lmax),zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax),vg(jmax,kmax,lmax),wg(jmax,kmax,lmax)
      real bt(jmax,kmax,lmax)
      integer kkr(kmax),kkp(kmax)
      real zx0(jmax,kmax),zy0(jmax,kmax),zz0(jmax,kmax),zt0(jmax,kmax)

      integer j,k,l,n,ierr,ireq(4)
      integer js,je,ks,ke,ls,le,ib,idir,iproc
      character*128 outfile,integer_string

      type(bc_t) :: grbc

c**** first executable statement

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
         idir = grbc%ibdir(ib)
         iproc = grbc%ibproc(ib)

c.. inviscid windtunnel wall bc
         if (grbc%ibtyp(ib).eq.4) then
         call bcwall_internal(q,xx,xy,xz,ug,yx,yy,yz,vg,zx,zy,zz,wg,
     <        zx0,zy0,zz0,zt0,kkr,kkp,js,je,ks,ke,ls,le,idir)

c.. wall bc at l = 1 (only interior portion of wall)
         elseif (grbc%ibtyp(ib).eq.5) then
         call bcwall(q,xx,xy,xz,ug,yx,yy,yz,vg,zx,zy,zz,wg,
     <        zx0,zy0,zz0,zt0,kkr,kkp,js,je,ks,ke,ls,le,idir)

c.. wall bc at l = 1 ground mesh (only interior portion of wall)
         elseif (grbc%ibtyp(ib).eq.6) then
         call bcwall(q,xx,xy,xz,ug-ug,yx,yy,yz,vg-vg,zx,zy,zz,wg-wg,
     <        zx0,zy0,zz0,zt0,kkr,kkp,js,je,ks,ke,ls,le,idir)

c.. extrapolate bc at k = 1 
         elseif (grbc%ibtyp(ib).eq.10) then
         call bcextp(q,x,y,z,yx,yy,yz,ug,vg,wg,js,je,ks,ke,ls,le,idir)

c.. symmetric bc for fixed wing at k = 1
         elseif (grbc%ibtyp(ib).eq.11) then
         call bcsym(q,js,je,ks,ke,ls,le,idir)

c.. Math: periodic bc
		 elseif (grbc%ibtyp(ib).eq.12) then
		 call bcperiodic(q,js,je,ks,ke,ls,le,idir)

c.. averaging bc at k = 1
         elseif (grbc%ibtyp(ib).eq.14) then
         call bcav(q,x,y,z,js,je,ks,ke,ls,le,idir)

c.. symmetric bc j (3 planes)
         elseif (grbc%ibtyp(ib).eq.22) then
         call bcaxisym(q,js,je,ks,ke,ls,le,idir)

c.. freesream bc
         elseif (grbc%ibtyp(ib).eq.47) then
          if (.not.iprecon) then
             call bcout(q,x,y,z,xx,xy,xz,yx,yy,yz,zx,zy,zz,
     >            ug,vg,wg,js,je,ks,ke,ls,le,idir)
          else
             call bcouttrkl(q,x,y,z,xx,xy,xz,yx,yy,yz,zx,zy,zz,
     >            ug,vg,wg,bt,js,je,ks,ke,ls,le,idir)
          endif

c.. freesream bc for hover
         elseif (grbc%ibtyp(ib).eq.48) then
          if(.not.iprecon) then
           call bcout_hover(q,x,y,z,xx,xy,xz,yx,yy,yz,zx,zy,zz,
     >             ug,vg,wg,js,je,ks,ke,ls,le,idir)
          else
           call bcouttrkl_hover(q,x,y,z,xx,xy,xz,yx,yy,yz,zx,zy,zz,
     >             ug,vg,wg,bt,js,je,ks,ke,ls,le,idir)
          endif

c.. averaging bc for wake
         elseif (grbc%ibtyp(ib).eq.51) then
         call bcwake(q,x,y,z,js,je,ks,ke,ls,le,idir)
     
c.. bc for parallel runs
         elseif (grbc%ibtyp(ib).eq.87) then
         call bcparallel(q,vnut,js,je,ks,ke,ls,le,idir,iproc)

c.. bc for parallel runs, with velocities rotation
         elseif (grbc%ibtyp(ib).eq.88) then
         call bcparallel_bg(q,vnut,js,je,ks,ke,ls,le,idir,iproc)

c.. cylindrical mesh bc for coaxial rotor
         elseif (grbc%ibtyp(ib).eq.90) then
         call bccoaxialcub(q,vnut,x,y,z,js,je,ks,ke,ls,le,idir,iproc)

         endif
      enddo
      return
      end

c***********************************************************************
      subroutine bcextp(q,x,y,z,yx,yy,yz,ug,vg,wg,js,je,ks,ke,ls,le,idir)

c**** Prologue : ******
c
c  boundary condition at k = 1 (root station)
c  first-order extrapolation
c  the call to this subroutine is reduntant with the use of characteristic 
c  extrapolation at bcoutf_wing
c  subroutine still included for completeness
c
c  last updated 07/13/04 by jaina
c***********************************************************************

      use params_global

c*******************  list of global variables used ********************
c
c     jmax,kmax,lmax,nd,jm,gm1
c***********************************************************************

      implicit none

      real q(jmax,kmax,lmax,nd)
      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)
      real yx(jmax,kmax,lmax),yy(jmax,kmax,lmax),yz(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax),vg(jmax,kmax,lmax),wg(jmax,kmax,lmax)
      integer js,je,ks,ke,ls,le,idir

c..   local variables

      integer j,k,l
      integer j1,j2,k1,k2,iadd,iadir
      real rk,rkp1,rkp2,h1,h2,foso
      real rhoext,uext,vext,wext,eext,pext
      real press

c**** first executable statement
      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then
        print*,'idir = ',idir,' is not implemented in bcextp'

      elseif(iadir.eq.2) then
       k  = ks
       k1 = k + iadd
       k2 = k1 + iadd
c     
c..first-order extrapolation of all quantities?
c
       do l = ls,le
         do j = js,je
            rk = 1./q(j,k,l,6)
            rkp1  = 1./q(j,k1,l,1)
            rkp2   = 1./q(j,k2,l,1)
            h1 = (x(j,k1,l)-x(j,k,l))**2
     <           +(y(j,k1,l)-y(j,k,l))**2
     <           +(z(j,k1,l)-z(j,k,l))**2
            h2 = (x(j,k2,l)-x(j,k1,l))**2
     <           +(y(j,k2,l)-y(j,k1,l))**2
     <           +(z(j,k2,l)-z(j,k1,l))**2
            foso = 0.5*sqrt(h1/h2)
c     foso = 0.99*foso
cjaina why not higher order?
            foso = 0.00*foso
c     
            rhoext = (1.+foso)*q(j,k1,l,1)*q(j,k1,l,6)-
     <           foso*q(j,k2,l,1)*q(j,k2,l,6)
            uext   = (1.+foso)*q(j,k1,l,2)*rkp1-foso*q(j,k2,l,2)*rkp2
            vext   = (1.+foso)*q(j,k1,l,3)*rkp1-foso*q(j,k2,l,3)*rkp2
            wext   = (1.+foso)*q(j,k1,l,4)*rkp1-foso*q(j,k2,l,4)*rkp2
            eext   = (1.+foso)*q(j,k1,l,5)*q(j,k1,l,6)-
     <           foso*q(j,k2,l,5)*q(j,k2,l,6)
            pext  = gm1*(eext - 0.5*rhoext*(uext**2+vext**2+wext**2))
c     
            q(j,k,l,1) = rhoext*rk
            q(j,k,l,2) = rhoext*uext*rk
            q(j,k,l,3) = rhoext*vext*rk
            q(j,k,l,4) = rhoext*wext*rk
            press = pext
            q(j,k,l,5) = press/gm1*rk + 0.5*(q(j,k,l,2)**2 +
     <           q(j,k,l,3)**2+q(j,k,l,4)**2)/q(j,k,l,1)
         enddo
       enddo

       elseif(iadir.eq.3) then 
         print*,'idir = ',idir,' is not implemented in bcextp'

       endif
c
      return
      end

c**********************************************************************
      subroutine bcav(q,x,y,z,js,je,ks,ke,ls,le,idir)

c**** Prologue : ******
c
c  boundary condition at k = 1 (root station)
c  averaging the values of k=1 when it collapses to a point
c***********************************************************************

      use params_global

c*******************  list of global variables used ********************
c
c     jmax,kmax,lmax,nd,jm,gm1
c***********************************************************************

      implicit none

      real q(jmax,kmax,lmax,nd)
      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)
      integer js,je,ks,ke,ls,le,idir

c..   local variables

      integer j,k,l
      integer k1,k2,iadd,iadir 
      real rk,rkp1,rkp2,h1,h2,foso
      real rhoext,uext,vext,wext,eext,pext
      real press
      real sumq1,sumq2,sumq3,sumq4,sumq5

c**** first executable statement

      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then
        print*,'idir = ',idir,' is not implemented in bcav'

      elseif(iadir.eq.2) then
       k  = ks
       k1 = k + iadd
       k2 = k1 + iadd
c     
c..first-order extrapolation of all quantities?
c
       do l = ls,le
         do j = js,je
            rk = 1./q(j,k,l,6)
            rkp1  = 1./q(j,k1,l,1)
            rkp2   = 1./q(j,k2,l,1)
            h1 = (x(j,k1,l)-x(j,k,l))**2
     <           +(y(j,k1,l)-y(j,k,l))**2
     <           +(z(j,k1,l)-z(j,k,l))**2
            h2 = (x(j,k2,l)-x(j,k1,l))**2
     <           +(y(j,k2,l)-y(j,k1,l))**2
     <           +(z(j,k2,l)-z(j,k1,l))**2
            foso = 0.5*sqrt(h1/h2)
c     foso = 0.99*foso
cjaina why not higher order?
            foso = 0.00*foso
c     
            rhoext = (1.+foso)*q(j,k1,l,1)*q(j,k1,l,6)-
     <           foso*q(j,k2,l,1)*q(j,k2,l,6)
            uext   = (1.+foso)*q(j,k1,l,2)*rkp1-foso*q(j,k2,l,2)*rkp2
            vext   = (1.+foso)*q(j,k1,l,3)*rkp1-foso*q(j,k2,l,3)*rkp2
            wext   = (1.+foso)*q(j,k1,l,4)*rkp1-foso*q(j,k2,l,4)*rkp2
            eext   = (1.+foso)*q(j,k1,l,5)*q(j,k1,l,6)-
     <           foso*q(j,k2,l,5)*q(j,k2,l,6)
            pext  = gm1*(eext - 0.5*rhoext*(uext**2+vext**2+wext**2))
c     
            q(j,k,l,1) = rhoext*rk
            q(j,k,l,2) = rhoext*uext*rk
            q(j,k,l,3) = rhoext*vext*rk
            q(j,k,l,4) = rhoext*wext*rk
            press = pext
            q(j,k,l,5) = press/gm1*rk + 0.5*(q(j,k,l,2)**2 +
     <           q(j,k,l,3)**2+q(j,k,l,4)**2)/q(j,k,l,1)
         enddo
       enddo
     
c..averaging the quantities

      do l = ls,le
         sumq1 = 0.0
         sumq2 = 0.0
         sumq3 = 0.0
         sumq4 = 0.0
         sumq5 = 0.0
         do j = js,je
            sumq1 = sumq1 + q(j,k,l,1)
            sumq2 = sumq2 + q(j,k,l,2)
            sumq3 = sumq3 + q(j,k,l,3)
            sumq4 = sumq4 + q(j,k,l,4)
            sumq5 = sumq5 + q(j,k,l,5)
         enddo
         sumq1 = sumq1/(je-js+1)
         sumq2 = sumq2/(je-js+1)
         sumq3 = sumq3/(je-js+1)
         sumq4 = sumq4/(je-js+1)
         sumq5 = sumq5/(je-js+1)
         do j = js,je
            q(j,k,l,1) = sumq1
            q(j,k,l,2) = sumq2
            q(j,k,l,3) = sumq3
            q(j,k,l,4) = sumq4
            q(j,k,l,5) = sumq5
         enddo
      enddo

       elseif(iadir.eq.3) then 
         print*,'idir = ',idir,' is not implemented in bcav'

       endif
         
      return
      end

c***********************************************************************
      subroutine bcout_hover(q,x,y,z,xx,xy,xz,yx,yy,yz,zx,zy,zz,
     >     ug,vg,wg,js,je,ks,ke,ls,le,idir)
c*** Prologue : ***********
c     
c  subroutine was written by dkarthik (2003)
c  have to add sink bc

c**********************************************************************

      use params_global

c*******************  list of global variables used ********************
c
c     jmax,kmax,lmax,nd,gm1,gamma,uinf,vinf,winf,rinf,pinf,ktip
c***********************************************************************

      implicit none

      real q(jmax,kmax,lmax,nd)
      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)
      real xx(jmax,kmax,lmax),xy(jmax,kmax,lmax),xz(jmax,kmax,lmax),
     <     yx(jmax,kmax,lmax),yy(jmax,kmax,lmax),yz(jmax,kmax,lmax),
     <     zx(jmax,kmax,lmax),zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax),vg(jmax,kmax,lmax),wg(jmax,kmax,lmax)
      integer js,je,ks,ke,ls,le,idir

c..local variables

      integer j,k,l,k1,l1,n,iadd,iadir
      real gm1i,gi,rl,rlm1,rhoext,uext,vext,wext,eext,pext
      real uind,vind,wind,mind2,rind,pind,snorm
      real zxn,zyn,zzn
      real xxn,xyn,xzn
      real yxn,yyn,yzn
      real uinn,r1,uexn,r2,qn,cspe,c2
      real vel2,vel2ext,velcheck
      real stj,entro,u,v,w,press,rj,rjp1,rjm1
      real ct,fnu,fnub,radb,radp,radd,radpo,theta,phi
      real radvel,uuvel,vvvel,wwvel
      real tcos,tsin

c**** first executable statement

      gm1i = 1./gm1
      gi   = 1./gamma
c
      ct = ctinp
c..induced velocity at the plane of the rotor based on thrust
      fnu = -0.5*winf+sqrt(0.25*winf**2+(fmtip*sqrt(ct/2.0))**2)
      if(fnu.eq.0) then
        radb = (1./sqrt(2.)+4.65*ct)*rartio
      else
        radb = (sqrt( (1-winf/fnu)/(2-winf/fnu) ) +4.65*ct)*rartio
      endif
c..correct, assuming that winf is zero!!!!!
      fnub = (0.5*rartio**2)/(radb**2)*fnu

      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then
         print*,'idir = ',idir,' is not implemented in bcout_hover'
      elseif(idir.eq.2) then
         print*,'idir = ',idir,' is not implemented in bcout_hover'
      elseif(idir.eq.-2) then
      k  = ks
      k1 = k + iadd
c     
      do  l = ls,le
      do  j = js,je
        rl    = 1./q(j,k,l,6)
        rlm1  = 1./q(j,k1,l,1)
               
c..   zeroeth-order extrapolation

        rhoext= q(j,k1,l,1)*q(j,k1,l,6)
        uext  = q(j,k1,l,2)*rlm1
        vext  = q(j,k1,l,3)*rlm1
        wext  = q(j,k1,l,4)*rlm1
        eext  = q(j,k1,l,5)*q(j,k1,l,6)
        pext  = gm1*(eext - 0.5*rhoext*(uext**2+vext**2+wext**2))
C..Calculate coordinates for sink flow
        radp  = sqrt(x(j,k,l)**2 + y(j,k,l)**2)
        radd  = sqrt(x(j,k,l)**2 + y(j,k,l)**2 + z(j,k,l)**2)
        theta = atan2(z(j,k,l),radp)
        phi   = atan2(x(j,k,l),y(j,k,l))
c.free stream
        radvel = 0.25*fnu*(rartio/radd)**2
        wwvel = radvel*sin(theta)
        vvvel = radvel*cos(theta)*cos(phi)
        uuvel = radvel*cos(theta)*sin(phi)
        uind = uinf - uuvel
        vind = vinf - vvvel
        wind = winf - wwvel
        mind2 = uind**2+vind**2+wind**2
        rind = rinf*(1.-.5*rinf/pinf*gm1*gi*(mind2-fsmach**2))**gm1i
        pind = pinf*(rind/rinf)**gamma
        
        snorm = -1.*iadd/sqrt(yx(j,k,l)**2+yy(j,k,l)**2+yz(j,k,l)**2)
        
        yxn = yx(j,k,l)*snorm
        yyn = yy(j,k,l)*snorm
        yzn = yz(j,k,l)*snorm
        
c..calculate riemann invariants

        uinn = (uind-ug(j,k,l))*yxn + (vind-vg(j,k,l))*yyn
     <       + (wind-wg(j,k,l))*yzn
        r1 = uinn -2.*sqrt(gamma*pind/rind)*gm1i
        uexn = (uext-ug(j,k,l))*yxn + (vext-vg(j,k,l))*yyn
     <       + (wext-wg(j,k,l))*yzn
        r2 = uexn +2.*sqrt(gamma*pext/rhoext)*gm1i

c..calculate normal velocity and speed of sound based on riemann

        qn = 0.5*(r1+r2)
        cspe = (r2-r1)*gm1*0.25
        c2 = cspe**2
               
c..is flow relatively subsonic or supersonic?
               
        vel2 = (uind-ug(j,k,l))**2+(vind-vg(j,k,l))**2
     <       +(wind-wg(j,k,l))**2
        vel2ext = (uext-ug(j,k,l))**2+(vext-vg(j,k,l))**2
     <       +(wext-wg(j,k,l))**2
        velcheck = abs(qn) !0.5*(vel2+vel2ext)
               
c..calculate contributions from interior and exterior

        if(qn .lt. 0) then
c..inflow boundary
           if(velcheck .lt. 1.0) then
c..   fix four and extrapolate one (pressure)
              stj = qn - uinn
              entro = rind**gamma/pind
              u = uind + stj*yxn
              v = vind + stj*yyn
              w = wind + stj*yzn
              q(j,k,l,1) = (c2*entro*gi)**gm1i
              press = c2*q(j,k,l,1)*gi
              q(j,k,l,1) = q(j,k,l,1)*rl
              q(j,k,l,2) = q(j,k,l,1)*u
              q(j,k,l,3) = q(j,k,l,1)*v
              q(j,k,l,4) = q(j,k,l,1)*w
              q(j,k,l,5) = press/gm1*rl + 0.5*(q(j,k,l,2)**2 +
     <             q(j,k,l,3)**2+q(j,k,l,4)**2)/q(j,k,l,1)
           else
c..fix five
              q(j,k,l,1) = rind*rl
              q(j,k,l,2) = rind*uind*rl
              q(j,k,l,3) = rind*vind*rl
              q(j,k,l,4) = rind*wind*rl
              press = pind
              q(j,k,l,5) = press/gm1*rl + 0.5*(q(j,k,l,2)**2 +
     <             q(j,k,l,3)**2+q(j,k,l,4)**2)/q(j,k,l,1)
           endif
        else
                  
c..outflow boundary
          
           if(velcheck .lt. 1.0) then
                     
c..prescribe p and extrapolate rho,u,v,and w
                     
              stj = qn - uexn
              entro = rhoext**gamma/pext
              u = uext + stj*yxn
              v = vext + stj*yyn
              w = wext + stj*yzn
              q(j,k,l,1) = (c2*entro*gi)**gm1i
              press = c2*q(j,k,l,1)*gi
              q(j,k,l,1) = q(j,k,l,1)*rl
              q(j,k,l,2) = q(j,k,l,1)*u
              q(j,k,l,3) = q(j,k,l,1)*v
              q(j,k,l,4) = q(j,k,l,1)*w
              q(j,k,l,5) = press/gm1*rl + 0.5*(q(j,k,l,2)**2 +
     <             q(j,k,l,3)**2+q(j,k,l,4)**2)/q(j,k,l,1)
           else
            
c..extrapolate five

              q(j,k,l,1) = rhoext*rl
              q(j,k,l,2) = rhoext*uext*rl
              q(j,k,l,3) = rhoext*vext*rl
              q(j,k,l,4) = rhoext*wext*rl
              press = pext
              q(j,k,l,5) = press/gm1*rl + 0.5*(q(j,k,l,2)**2 +
     <             q(j,k,l,3)**2+q(j,k,l,4)**2)/q(j,k,l,1)
           endif
        endif
               
      enddo
      enddo

      elseif(idir.eq.3) then

      l  = ls
      l1 = l + iadd
c
      do j = js,je
        radp = 0.d0
      do k = ks,ke
        radpo = radp
        rl    = 1./q(j,k,l,6)
        rlm1  = 1./q(j,k,l1,1)
c..zeroeth-order extrapolation
        rhoext= q(j,k,l1,1)*q(j,k,l1,6)
        uext  = q(j,k,l1,2)*rlm1
        vext  = q(j,k,l1,3)*rlm1
        wext  = q(j,k,l1,4)*rlm1
        eext  = q(j,k,l1,5)*q(j,k,l1,6)
        pext  = gm1*(eext - 0.5*rhoext*(uext**2+vext**2+wext**2))
c..Calculate coordinates for sink flow
        radp  = sqrt(x(j,k,l)**2 + y(j,k,l)**2)
        radd  = sqrt(x(j,k,l)**2 + y(j,k,l)**2 + z(j,k,l)**2)
        theta = atan2(z(j,k,l),radp)
        phi   = atan2(x(j,k,l),y(j,k,l))
c.free stream
        radvel = 0.25*fnu*(rartio/radd)**2
        wwvel = radvel*sin(theta)
        vvvel = radvel*cos(theta)*cos(phi)
        uuvel = radvel*cos(theta)*sin(phi)
	uind = uinf - uuvel 
        vind = vinf - vvvel
        wind = winf - wwvel
        mind2 = uind**2+vind**2+wind**2
        rind = rinf*(1.-.5*rinf/pinf*gm1*gi*(mind2-fsmach**2))**gm1i
        pind = pinf*(rind/rinf)**gamma
c..Calculate velocities due to actuator disc at bottom
        if(radp .le. radb.and.z(j,k,l).lt.0.) wind = wind - 2.*fnub
        if(radp .gt. radb. and. radpo .le. radb .and. z(j,k,l).lt.0.)
     <          wind = wind - fnub

        snorm = -1.*iadd/sqrt(zx(j,k,l)**2+zy(j,k,l)**2+zz(j,k,l)**2)
        zxn = zx(j,k,l)*snorm
        zyn = zy(j,k,l)*snorm
        zzn = zz(j,k,l)*snorm
c..calculate riemann invariants
        uinn = (uind-ug(j,k,l))*zxn + (vind-vg(j,k,l))*zyn
     <       + (wind-wg(j,k,l))*zzn
        r1 = uinn -2.*sqrt(gamma*pind/rind)*gm1i
        uexn = (uext-ug(j,k,l))*zxn + (vext-vg(j,k,l))*zyn
     <       + (wext-wg(j,k,l))*zzn
        r2 = uexn +2.*sqrt(gamma*pext/rhoext)*gm1i
c..calculate normal velocity and speed of sound based on riemann
        qn = 0.5*(r1+r2)
        cspe = (r2-r1)*gm1*0.25
        c2 = cspe**2
c..is flow relatively subsonic or supersonic?
        vel2 = (uind-ug(j,k,l))**2+(vind-vg(j,k,l))**2
     <                            +(wind-wg(j,k,l))**2
        vel2ext = (uext-ug(j,k,l))**2+(vext-vg(j,k,l))**2
     <                            +(wext-wg(j,k,l))**2
        velcheck = abs(qn) !0.5*(vel2+vel2ext)
c..calculate contributions from interior and exterior
        if(qn .lt. 0) then
c..inflow boundary
          if(velcheck .lt. 1.0) then
c..fix four and extrapolate one (pressure)
            stj = qn - uinn
            entro = rind**gamma/pind
            u = uind + stj*zxn
            v = vind + stj*zyn
            w = wind + stj*zzn
            q(j,k,l,1) = (c2*entro*gi)**gm1i
            press = c2*q(j,k,l,1)*gi
            q(j,k,l,1) = q(j,k,l,1)*rl
            q(j,k,l,2) = q(j,k,l,1)*u
            q(j,k,l,3) = q(j,k,l,1)*v
            q(j,k,l,4) = q(j,k,l,1)*w
            q(j,k,l,5) = press/gm1*rl + 0.5*(q(j,k,l,2)**2 +
     <                   q(j,k,l,3)**2+q(j,k,l,4)**2)/q(j,k,l,1)
          else
c..fix five
            q(j,k,l,1) = rind*rl
            q(j,k,l,2) = rind*uind*rl
            q(j,k,l,3) = rind*vind*rl
            q(j,k,l,4) = rind*wind*rl
            press = pind
            q(j,k,l,5) = press/gm1*rl + 0.5*(q(j,k,l,2)**2 +
     <                   q(j,k,l,3)**2+q(j,k,l,4)**2)/q(j,k,l,1)
          endif
        else
c..outflow boundary
          if(velcheck .lt. 1.0) then
c..prescribe p and extrapolate rho,u,v,and w
            stj = qn - uexn
            entro = rhoext**gamma/pext
            u = uext + stj*zxn
            v = vext + stj*zyn
            w = wext + stj*zzn
            q(j,k,l,1) = (c2*entro*gi)**gm1i
            press = c2*q(j,k,l,1)*gi
            q(j,k,l,1) = q(j,k,l,1)*rl
            q(j,k,l,2) = q(j,k,l,1)*u
            q(j,k,l,3) = q(j,k,l,1)*v
            q(j,k,l,4) = q(j,k,l,1)*w
            q(j,k,l,5) = press/gm1*rl + 0.5*(q(j,k,l,2)**2 +
     <                   q(j,k,l,3)**2+q(j,k,l,4)**2)/q(j,k,l,1)
          else
c..extrapolate five
            q(j,k,l,1) = rhoext*rl
            q(j,k,l,2) = rhoext*uext*rl
            q(j,k,l,3) = rhoext*vext*rl
            q(j,k,l,4) = rhoext*wext*rl
            press = pext
            q(j,k,l,5) = press/gm1*rl + 0.5*(q(j,k,l,2)**2 +
     <                   q(j,k,l,3)**2+q(j,k,l,4)**2)/q(j,k,l,1)
          endif
        endif
c
      enddo
      enddo

      elseif (idir.eq.-3) then

      l  = ls
      l1 = l + iadd

c..remove l=lmax condition to check

      do j = js,je
      do k = ks,ke
        rl    = 1./q(j,k,l,6)
        rlm1  = 1./q(j,k,l1,1)
c..zeroeth-order extrapolation
        rhoext= q(j,k,l1,1)*q(j,k,l1,6)
        uext  = q(j,k,l1,2)*rlm1
        vext  = q(j,k,l1,3)*rlm1
        wext  = q(j,k,l1,4)*rlm1
        eext  = q(j,k,l1,5)*q(j,k,l1,6)
        pext  = gm1*(eext - 0.5*rhoext*(uext**2+vext**2+wext**2))
C..Calculate coordinates for sink flow
        radp  = sqrt(x(j,k,l)**2 + y(j,k,l)**2)
        radd  = sqrt(x(j,k,l)**2 + y(j,k,l)**2 + z(j,k,l)**2)
        theta = atan2(z(j,k,l),radp)
        phi   = atan2(x(j,k,l),y(j,k,l))
c.free stream
        radvel = 0.25*fnu*(rartio/radd)**2
        wwvel = radvel*sin(theta)
        vvvel = radvel*cos(theta)*cos(phi)
        uuvel = radvel*cos(theta)*sin(phi)
        uind = uinf - uuvel
        vind = vinf - vvvel
        wind = winf - wwvel
        mind2 = uind**2+vind**2+wind**2
        rind = rinf*(1.-.5*rinf/pinf*gm1*gi*(mind2-fsmach**2))**gm1i
        pind = pinf*(rind/rinf)**gamma

        snorm = -1.*iadd/sqrt(zx(j,k,l)**2+zy(j,k,l)**2+zz(j,k,l)**2)
        zxn = zx(j,k,l)*snorm
        zyn = zy(j,k,l)*snorm
        zzn = zz(j,k,l)*snorm
c..calculate riemann invariants
        uinn = (uind-ug(j,k,l))*zxn + (vind-vg(j,k,l))*zyn
     <       + (wind-wg(j,k,l))*zzn
        r1 = uinn -2.*sqrt(gamma*pind/rind)*gm1i
        uexn = (uext-ug(j,k,l))*zxn + (vext-vg(j,k,l))*zyn
     <       + (wext-wg(j,k,l))*zzn
        r2 = uexn +2.*sqrt(gamma*pext/rhoext)*gm1i
c..calculate normal velocity and speed of sound based on riemann
        qn = 0.5*(r1+r2)
        cspe = (r2-r1)*gm1*0.25
        c2 = cspe**2
c..is flow relatively subsonic or supersonic?
        vel2 = (uind-ug(j,k,l))**2+(vind-vg(j,k,l))**2
     <                            +(wind-wg(j,k,l))**2
        vel2ext = (uext-ug(j,k,l))**2+(vext-vg(j,k,l))**2
     <                            +(wext-wg(j,k,l))**2
        velcheck = abs(qn) !0.5*(vel2+vel2ext)
c..calculate contributions from interior and exterior
        if(qn .lt. 0) then
c..inflow boundary
          if(velcheck .lt. 1.0) then
c..fix four and extrapolate one (pressure)
            stj = qn - uinn
            entro = rind**gamma/pind
            u = uind + stj*zxn
            v = vind + stj*zyn
            w = wind + stj*zzn
            q(j,k,l,1) = (c2*entro*gi)**gm1i
            press = c2*q(j,k,l,1)*gi
            q(j,k,l,1) = q(j,k,l,1)*rl
            q(j,k,l,2) = q(j,k,l,1)*u
            q(j,k,l,3) = q(j,k,l,1)*v
            q(j,k,l,4) = q(j,k,l,1)*w
            q(j,k,l,5) = press/gm1*rl + 0.5*(q(j,k,l,2)**2 +
     <                   q(j,k,l,3)**2+q(j,k,l,4)**2)/q(j,k,l,1)
          else
c..fix five
            q(j,k,l,1) = rind*rl
            q(j,k,l,2) = rind*uind*rl
            q(j,k,l,3) = rind*vind*rl
            q(j,k,l,4) = rind*wind*rl
            press = pind
            q(j,k,l,5) = press/gm1*rl + 0.5*(q(j,k,l,2)**2 +
     <                   q(j,k,l,3)**2+q(j,k,l,4)**2)/q(j,k,l,1)
          endif
        else
c..outflow boundary
          if(velcheck .lt. 1.0) then
c..prescribe p and extrapolate rho,u,v,and w
            stj = qn - uexn
            entro = rhoext**gamma/pext
            u = uext + stj*zxn
            v = vext + stj*zyn
            w = wext + stj*zzn
            q(j,k,l,1) = (c2*entro*gi)**gm1i
            press = c2*q(j,k,l,1)*gi
            q(j,k,l,1) = q(j,k,l,1)*rl
            q(j,k,l,2) = q(j,k,l,1)*u
            q(j,k,l,3) = q(j,k,l,1)*v
            q(j,k,l,4) = q(j,k,l,1)*w
            q(j,k,l,5) = press/gm1*rl + 0.5*(q(j,k,l,2)**2 +
     <                   q(j,k,l,3)**2+q(j,k,l,4)**2)/q(j,k,l,1)
          else
c..extrapolate five
            q(j,k,l,1) = rhoext*rl
            q(j,k,l,2) = rhoext*uext*rl
            q(j,k,l,3) = rhoext*vext*rl
            q(j,k,l,4) = rhoext*wext*rl
            press = pext
            q(j,k,l,5) = press/gm1*rl + 0.5*(q(j,k,l,2)**2 +
     <                   q(j,k,l,3)**2+q(j,k,l,4)**2)/q(j,k,l,1)
          endif
        endif
c
      enddo
      enddo

      endif

      return
      end

c**********************************************************************
      subroutine bcouttrkl_hover(q,x,y,z,xx,xy,xz,yx,yy,yz,zx,zy,zz,
     >     ug,vg,wg,bt,js,je,ks,ke,ls,le,idir)
c*** Prologue : ***********
c     
c  Characteristically extrapolate the flow variables at the outer
c  boundaries using the Riemann invariants
c
c**********************************************************************

      use params_global

c*******************  list of global variables used ********************
c
c     jmax,kmax,lmax,nd,gm1,gamma,uinf,vinf,winf,rinf,pinf,ktip
c***********************************************************************

      implicit none

      real q(jmax,kmax,lmax,nd)
      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)
      real xx(jmax,kmax,lmax),xy(jmax,kmax,lmax),xz(jmax,kmax,lmax),
     <     yx(jmax,kmax,lmax),yy(jmax,kmax,lmax),yz(jmax,kmax,lmax),
     <     zx(jmax,kmax,lmax),zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax),vg(jmax,kmax,lmax),wg(jmax,kmax,lmax)
      real bt(jmax,kmax,lmax)
      integer js,je,ks,ke,ls,le,idir

c..local variables
      real,allocatable :: qext(:),qbar(:),delq(:),delw(:)
      integer j,k,l,k1,l1,n,iadd,iadir
      real rhoext,uext,vext,wext,eext,pext
      real uind,vind,wind,mind2,rind,pind,snorm
      real zxn,zyn,zzn
      real xxn,xyn,xzn
      real yxn,yyn,yzn
      real rho,rrho,uu,vv,ww,e,uvw,cjkl,c2i,ge,qq,gi,gm1i
      real ct,fnu,fnub,radb,radp,radd,radpo,theta,phi
      real radvel,uuvel,vvvel,wwvel
      real tcos,tsin
      real Z1,R,S,Zi,bSq

      allocate(qext(5),qbar(5),delq(5),delw(5))

c**** first executable statement

      ct = ctinp
c..induced velocity at the plane of the rotor based on thrust
      fnu = -0.5*winf+sqrt(0.25*winf**2+(fmtip*sqrt(ct/2.0))**2)
      if(fnu.eq.0) then
        radb = (1./sqrt(2.)+4.65*ct)*rartio
      else
        radb = (sqrt( (1-winf/fnu)/(2-winf/fnu) ) +4.65*ct)*rartio
      endif
c..correct, assuming that winf is zero!!!!!
      fnub = (0.5*rartio**2)/(radb**2)*fnu

      iadd = sign(1,idir)
      iadir  = abs(idir)
 
      if(iadir.eq.1) then
         print*,'idir = ',idir,' is not implemented in bcouttrkl_hover'
      elseif(idir.eq.2) then
         print*,'idir = ',idir,' is not implemented in bcouttrkl_hover'
      elseif(idir.eq.-2) then

      k  = ks
      k1 = k + iadd
c     
      do  l = ls,le
      do  j = js,je

C..Calculate coordinates for sink flow
        radp  = sqrt(x(j,k,l)**2 + y(j,k,l)**2)
        radd  = sqrt(x(j,k,l)**2 + y(j,k,l)**2 + z(j,k,l)**2)
        theta = atan2(z(j,k,l),radp)
        phi   = atan2(x(j,k,l),y(j,k,l))

c..free stream
        radvel = 0.25*fnu*(rartio/radd)**2
        wwvel = radvel*sin(theta)
        vvvel = radvel*cos(theta)*cos(phi)
        uuvel = radvel*cos(theta)*sin(phi)
        uind = uinf - uuvel
        vind = vinf - vvvel
        wind = winf - wwvel
        mind2 = uind**2+vind**2+wind**2
        rind = rinf*(1.-.5*rinf/pinf*gm1*gi*(mind2-fsmach**2))**gm1i
        pind = pinf*(rind/rinf)**gamma

        snorm = -1.*iadd/sqrt(yx(j,k,l)**2+yy(j,k,l)**2+yz(j,k,l)**2)
        yxn = yx(j,k,l)*snorm
        yyn = yy(j,k,l)*snorm
        yzn = yz(j,k,l)*snorm

c..zeroeth-order extrapolation
        qext(1) = q(j,k1,l,1)*q(j,k1,l,6)
        qext(2) = q(j,k1,l,2)*q(j,k1,l,6)
        qext(3) = q(j,k1,l,3)*q(j,k1,l,6)
        qext(4) = q(j,k1,l,4)*q(j,k1,l,6)
        qext(5) = q(j,k1,l,5)*q(j,k1,l,6)

        qbar(1) = 0.5*(qext(1)+rinf)
        qbar(2) = 0.5*(qext(2)+rinf*uind)
        qbar(3) = 0.5*(qext(3)+rinf*vind)
        qbar(4) = 0.5*(qext(4)+rinf*wind)
        qbar(5) = 0.5*(qext(5)+einf)

        delq(1) = 0.5*(rinf-qext(1))
        delq(2) = 0.5*(rinf*uind-qext(2))
        delq(3) = 0.5*(rinf*vind-qext(3))
        delq(4) = 0.5*(rinf*wind-qext(4))
        delq(5) = 0.5*(einf-qext(5))

        rho = qbar(1)
        rrho = 1./rho
        uu = qbar(2)*rrho
        vv = qbar(3)*rrho
        ww = qbar(4)*rrho
        e  = qbar(5)*rrho
        uvw = 0.5*(uu*uu+vv*vv+ww*ww)
        cjkl = sqrt(ggm1*(e-uvw))
        c2i = 1./(cjkl*cjkl)
        ge = gamma*e - gm1*uvw
        qq = (uu-ug(j,k,l))*yxn+(vv-vg(j,k,l))*yyn+(ww-wg(j,k,l))*yzn

        bSq = Mp**2/(bt(j,k,l)-Mp**2*(bt(j,k,l)-1))
        
        Z1 = sqrt((1.-bSq)*qq*(1.-bSq)*qq+4.*bSq*cjkl*cjkl)
        Zi = 1./Z1
        R = 0.5*((1.-bSq)*qq+Z1)
        S = 0.5*((1.-bSq)*qq-Z1)

c..multiplying by sign(lamda_pa)*inv(X_pa)

        delw(1) = (yxn*(1.-uvw*gm1*c2i)-(vv*yzn-ww*yyn)*rrho)*
     c                                    sign(1.,qq)*delq(1)
        delw(1) = delw(1) +(uu*yxn*c2i*gm1)*sign(1.,qq)*delq(2)
        delw(1) = delw(1) +(vv*yxn*c2i*gm1+yzn*rrho)*sign(1.,qq)*delq(3)
        delw(1) = delw(1) +(ww*yxn*c2i*gm1-yyn*rrho)*sign(1.,qq)*delq(4)
        delw(1) = delw(1) - yxn*gm1*c2i*sign(1.,qq)*delq(5)
       
        delw(2) = (yyn*(1.-uvw*gm1*c2i)-(ww*yxn-uu*yzn)*rrho)*
     c                                    sign(1.,qq)*delq(1)
        delw(2) = delw(2) +(uu*yyn*c2i*gm1-yzn*rrho)*sign(1.,qq)*delq(2)
        delw(2) = delw(2) +(vv*yyn*c2i*gm1)*sign(1.,qq)*delq(3)
        delw(2) = delw(2) +(ww*yyn*c2i*gm1+yxn*rrho)*sign(1.,qq)*delq(4)
        delw(2) = delw(2) - yyn*gm1*c2i*sign(1.,qq)*delq(5)

        delw(3) = (yzn*(1.-uvw*gm1*c2i)-(uu*yyn-vv*yxn)*rrho)*
     c                                    sign(1.,qq)*delq(1)
        delw(3) = delw(3) +(uu*yzn*c2i*gm1+yyn*rrho)*sign(1.,qq)*delq(2)
        delw(3) = delw(3) +(vv*yzn*c2i*gm1-yxn*rrho)*sign(1.,qq)*delq(3)
        delw(3) = delw(3) +(ww*yzn*c2i*gm1)*sign(1.,qq)*delq(4)
        delw(3) = delw(3) - yzn*gm1*c2i*sign(1.,qq)*delq(5)
       
        delw(4) = (-S*uvw*gm1*c2i-qq*bSq)*Zi*sign(1.,R+qq*bSq)*delq(1)
        delw(4) = delw(4) + (yxn*bSq + S*uu*gm1*c2i)*Zi*
     c                      sign(1.,R+qq*bSq)*delq(2)
        delw(4) = delw(4) + (yyn*bSq + S*vv*gm1*c2i)*Zi*
     c                      sign(1.,R+qq*bSq)*delq(3)
        delw(4) = delw(4) + (yzn*bSq + S*ww*gm1*c2i)*Zi*
     c                      sign(1.,R+qq*bSq)*delq(4)
        delw(4) = delw(4) - S*gm1*c2i*Zi*sign(1.,R+qq*bSq)*delq(5)

        delw(5) = (R*uvw*gm1*c2i+qq*bSq)*Zi*sign(1.,S+qq*bSq)*delq(1)
        delw(5) = delw(5) - (yxn*bSq + R*uu*gm1*c2i)*Zi*
     c                      sign(1.,S+qq*bSq)*delq(2)
        delw(5) = delw(5) - (yyn*bSq + R*vv*gm1*c2i)*Zi*
     c                      sign(1.,S+qq*bSq)*delq(3)
        delw(5) = delw(5) - (yzn*bSq + R*ww*gm1*c2i)*Zi*
     c                      sign(1.,S+qq*bSq)*delq(4)
        delw(5) = delw(5) +  R*gm1*c2i*Zi*sign(1.,S+qq*bSq)*delq(5)

c.. multiplying by (X_pa)

        delq(1) = yxn*delw(1)+yyn*delw(2)+yzn*delw(3)+delw(4)+delw(5)

        delq(2) = uu*yxn*delw(1)+(uu*yyn-yzn*rho)*delw(2)
        delq(2) = delq(2) + (uu*yzn+yyn*rho)*delw(3)
        delq(2) = delq(2) + (uu+R*yxn/bSq)*delw(4)
        delq(2) = delq(2) + (uu+S*yxn/bSq)*delw(5)

        delq(3) = (vv*yxn+yzn*rho)*delw(1)+vv*yyn*delw(2)
        delq(3) = delq(3) + (vv*yzn-yxn*rho)*delw(3)
        delq(3) = delq(3) + (vv+R*yyn/bSq)*delw(4)
        delq(3) = delq(3) + (vv+S*yyn/bSq)*delw(5)

        delq(4) = (ww*yxn-yyn*rho)*delw(1)+(ww*yyn+yxn*rho)*delw(2)
        delq(4) = delq(4) + ww*yzn*delw(3)
        delq(4) = delq(4) + (ww+R*yzn/bSq)*delw(4)
        delq(4) = delq(4) + (ww+S*yzn/bSq)*delw(5)

        delq(5) = (uvw*yxn+rho*(vv*yzn-ww*yyn))*delw(1)
        delq(5) = delq(5) + (uvw*yyn+rho*(ww*yxn-uu*yzn))*delw(2)
        delq(5) = delq(5) + (uvw*yzn+rho*(uu*yyn-vv*yxn))*delw(3)
        delq(5) = delq(5) + (ge+R*qq/bSq)*delw(4)+(ge+S*qq/bSq)*delw(5)

c.. qbar - delq
        q(j,k,l,1) = (qbar(1)-delq(1))/q(j,k,l,6)
        q(j,k,l,2) = (qbar(2)-delq(2))/q(j,k,l,6)
        q(j,k,l,3) = (qbar(3)-delq(3))/q(j,k,l,6)
        q(j,k,l,4) = (qbar(4)-delq(4))/q(j,k,l,6)
        q(j,k,l,5) = (qbar(5)-delq(5))/q(j,k,l,6)

      enddo
      enddo

      elseif(idir.eq.3) then

      l  = ls
      l1 = l + iadd
c
      do j = js,je
        radp = 0.d0
      do k = ks,ke
        radpo = radp
c..Calculate coordinates for sink flow
        radp  = sqrt(x(j,k,l)**2 + y(j,k,l)**2)
        radd  = sqrt(x(j,k,l)**2 + y(j,k,l)**2 + z(j,k,l)**2)
        theta = atan2(z(j,k,l),radp)
        phi   = atan2(x(j,k,l),y(j,k,l))
c.free stream
        radvel = 0.25*fnu*(rartio/radd)**2
        wwvel = radvel*sin(theta)
        vvvel = radvel*cos(theta)*cos(phi)
        uuvel = radvel*cos(theta)*sin(phi)
	uind = uinf - uuvel 
        vind = vinf - vvvel
        wind = winf - wwvel
        mind2 = uind**2+vind**2+wind**2
        rind = rinf*(1.-.5*rinf/pinf*gm1*gi*(mind2-fsmach**2))**gm1i
        pind = pinf*(rind/rinf)**gamma
c..Calculate velocities due to actuator disc at bottom
        if(radp .le. radb.and.z(j,k,l).lt.0.) wind = wind - 2.*fnub
        if(radp .gt. radb. and. radpo .le. radb .and. z(j,k,l).lt.0.)
     <          wind = wind - fnub

        snorm = -1.*iadd/sqrt(zx(j,k,l)**2+zy(j,k,l)**2+zz(j,k,l)**2)
        zxn = zx(j,k,l)*snorm
        zyn = zy(j,k,l)*snorm
        zzn = zz(j,k,l)*snorm

c..zeroeth-order extrapolation
        qext(1) = q(j,k,l1,1)*q(j,k,l1,6)
        qext(2) = q(j,k,l1,2)*q(j,k,l1,6)
        qext(3) = q(j,k,l1,3)*q(j,k,l1,6)
        qext(4) = q(j,k,l1,4)*q(j,k,l1,6)
        qext(5) = q(j,k,l1,5)*q(j,k,l1,6)

        qbar(1) = 0.5*(qext(1)+rinf)
        qbar(2) = 0.5*(qext(2)+rinf*uind)
        qbar(3) = 0.5*(qext(3)+rinf*vind)
        qbar(4) = 0.5*(qext(4)+rinf*wind)
        qbar(5) = 0.5*(qext(5)+einf)

        delq(1) = 0.5*(rinf-qext(1))
        delq(2) = 0.5*(rinf*uind-qext(2))
        delq(3) = 0.5*(rinf*vind-qext(3))
        delq(4) = 0.5*(rinf*wind-qext(4))
        delq(5) = 0.5*(einf-qext(5))

        rho = qbar(1)
        rrho = 1./rho
        uu = qbar(2)*rrho
        vv = qbar(3)*rrho
        ww = qbar(4)*rrho
        e  = qbar(5)*rrho
        uvw = 0.5*(uu*uu+vv*vv+ww*ww)
        cjkl = sqrt(ggm1*(e-uvw))
        c2i = 1./(cjkl*cjkl)
        ge = gamma*e - gm1*uvw
        qq = (uu-ug(j,k,l))*zxn+(vv-vg(j,k,l))*zyn+(ww-wg(j,k,l))*zzn
        
        bSq = Mp**2/(bt(j,k,l)-Mp**2*(bt(j,k,l)-1))
        
        Z1 = sqrt((1.-bSq)*qq*(1.-bSq)*qq+4.*bSq*cjkl*cjkl)
        Zi = 1./Z1
        R = 0.5*((1.-bSq)*qq+Z1)
        S = 0.5*((1.-bSq)*qq-Z1)

c..multiplying by sign(lamda_pa)*inv(X_pa)
        delw(1) = (zxn*(1.-uvw*gm1*c2i)-(vv*zzn-ww*zyn)*rrho)*
     c                                    sign(1.,qq)*delq(1)
        delw(1) = delw(1) +(uu*zxn*c2i*gm1)*sign(1.,qq)*delq(2)
        delw(1) = delw(1) +(vv*zxn*c2i*gm1+zzn*rrho)*sign(1.,qq)*delq(3)
        delw(1) = delw(1) +(ww*zxn*c2i*gm1-zyn*rrho)*sign(1.,qq)*delq(4)
        delw(1) = delw(1) - zxn*gm1*c2i*sign(1.,qq)*delq(5)
       
        delw(2) = (zyn*(1.-uvw*gm1*c2i)-(ww*zxn-uu*zzn)*rrho)*
     c                                    sign(1.,qq)*delq(1)
        delw(2) = delw(2) +(uu*zyn*c2i*gm1-zzn*rrho)*sign(1.,qq)*delq(2)
        delw(2) = delw(2) +(vv*zyn*c2i*gm1)*sign(1.,qq)*delq(3)
        delw(2) = delw(2) +(ww*zyn*c2i*gm1+zxn*rrho)*sign(1.,qq)*delq(4)
        delw(2) = delw(2) - zyn*gm1*c2i*sign(1.,qq)*delq(5)

        delw(3) = (zzn*(1.-uvw*gm1*c2i)-(uu*zyn-vv*zxn)*rrho)*
     c                                    sign(1.,qq)*delq(1)
        delw(3) = delw(3) +(uu*zzn*c2i*gm1+zyn*rrho)*sign(1.,qq)*delq(2)
        delw(3) = delw(3) +(vv*zzn*c2i*gm1-zxn*rrho)*sign(1.,qq)*delq(3)
        delw(3) = delw(3) +(ww*zzn*c2i*gm1)*sign(1.,qq)*delq(4)
        delw(3) = delw(3) - zzn*gm1*c2i*sign(1.,qq)*delq(5)
       
        delw(4) = (-S*uvw*gm1*c2i-qq*bSq)*Zi*sign(1.,R+qq*bSq)*delq(1)
        delw(4) = delw(4) + (zxn*bSq + S*uu*gm1*c2i)*Zi*
     c                      sign(1.,R+qq*bSq)*delq(2)
        delw(4) = delw(4) + (zyn*bSq + S*vv*gm1*c2i)*Zi*
     c                      sign(1.,R+qq*bSq)*delq(3)
        delw(4) = delw(4) + (zzn*bSq + S*ww*gm1*c2i)*Zi*
     c                      sign(1.,R+qq*bSq)*delq(4)
        delw(4) = delw(4) - S*gm1*c2i*Zi*sign(1.,R+qq*bSq)*delq(5)

        delw(5) = (R*uvw*gm1*c2i+qq*bSq)*Zi*sign(1.,S+qq*bSq)*delq(1)
        delw(5) = delw(5) - (zxn*bSq + R*uu*gm1*c2i)*Zi*
     c                      sign(1.,S+qq*bSq)*delq(2)
        delw(5) = delw(5) - (zyn*bSq + R*vv*gm1*c2i)*Zi*
     c                      sign(1.,S+qq*bSq)*delq(3)
        delw(5) = delw(5) - (zzn*bSq + R*ww*gm1*c2i)*Zi*
     c                      sign(1.,S+qq*bSq)*delq(4)
        delw(5) = delw(5) +  R*gm1*c2i*Zi*sign(1.,S+qq*bSq)*delq(5)

c.. multiplying by (X_pa)

        delq(1) = zxn*delw(1)+zyn*delw(2)+zzn*delw(3)+delw(4)+delw(5)

        delq(2) = uu*zxn*delw(1)+(uu*zyn-zzn*rho)*delw(2)
        delq(2) = delq(2) + (uu*zzn+zyn*rho)*delw(3)
        delq(2) = delq(2) + (uu+R*zxn/bSq)*delw(4)
        delq(2) = delq(2) + (uu+S*zxn/bSq)*delw(5)

        delq(3) = (vv*zxn+zzn*rho)*delw(1)+vv*zyn*delw(2)
        delq(3) = delq(3) + (vv*zzn-zxn*rho)*delw(3)
        delq(3) = delq(3) + (vv+R*zyn/bSq)*delw(4)
        delq(3) = delq(3) + (vv+S*zyn/bSq)*delw(5)

        delq(4) = (ww*zxn-zyn*rho)*delw(1)+(ww*zyn+zxn*rho)*delw(2)
        delq(4) = delq(4) + ww*zzn*delw(3)
        delq(4) = delq(4) + (ww+R*zzn/bSq)*delw(4)
        delq(4) = delq(4) + (ww+S*zzn/bSq)*delw(5)

        delq(5) = (uvw*zxn+rho*(vv*zzn-ww*zyn))*delw(1)
        delq(5) = delq(5) + (uvw*zyn+rho*(ww*zxn-uu*zzn))*delw(2)
        delq(5) = delq(5) + (uvw*zzn+rho*(uu*zyn-vv*zxn))*delw(3)
        delq(5) = delq(5) + (ge+R*qq/bSq)*delw(4)+(ge+S*qq/bSq)*delw(5)

c.. qbar - delq
        q(j,k,l,1) = (qbar(1)-delq(1))/q(j,k,l,6)
        q(j,k,l,2) = (qbar(2)-delq(2))/q(j,k,l,6)
        q(j,k,l,3) = (qbar(3)-delq(3))/q(j,k,l,6)
        q(j,k,l,4) = (qbar(4)-delq(4))/q(j,k,l,6)
        q(j,k,l,5) = (qbar(5)-delq(5))/q(j,k,l,6)

      enddo
      enddo

      elseif (idir.eq.-3) then

      l  = ls
      l1 = l + iadd

      do j = js,je
      do k = ks,ke

C..Calculate coordinates for sink flow
        radp  = sqrt(x(j,k,l)**2 + y(j,k,l)**2)
        radd  = sqrt(x(j,k,l)**2 + y(j,k,l)**2 + z(j,k,l)**2)
        theta = atan2(z(j,k,l),radp)
        phi   = atan2(x(j,k,l),y(j,k,l))

c.free stream
        radvel = 0.25*fnu*(rartio/radd)**2
        wwvel = radvel*sin(theta)
        vvvel = radvel*cos(theta)*cos(phi)
        uuvel = radvel*cos(theta)*sin(phi)
        uind = uinf - uuvel
        vind = vinf - vvvel
        wind = winf - wwvel
        mind2 = uind**2+vind**2+wind**2
        rind = rinf*(1.-.5*rinf/pinf*gm1*gi*(mind2-fsmach**2))**gm1i
        pind = pinf*(rind/rinf)**gamma

        snorm = -1.*iadd/sqrt(zx(j,k,l)**2+zy(j,k,l)**2+zz(j,k,l)**2)
        zxn = zx(j,k,l)*snorm
        zyn = zy(j,k,l)*snorm
        zzn = zz(j,k,l)*snorm

c..zeroeth-order extrapolation
        qext(1) = q(j,k,l1,1)*q(j,k,l1,6)
        qext(2) = q(j,k,l1,2)*q(j,k,l1,6)
        qext(3) = q(j,k,l1,3)*q(j,k,l1,6)
        qext(4) = q(j,k,l1,4)*q(j,k,l1,6)
        qext(5) = q(j,k,l1,5)*q(j,k,l1,6)

        qbar(1) = 0.5*(qext(1)+rinf)
        qbar(2) = 0.5*(qext(2)+rinf*uind)
        qbar(3) = 0.5*(qext(3)+rinf*vind)
        qbar(4) = 0.5*(qext(4)+rinf*wind)
        qbar(5) = 0.5*(qext(5)+einf)

        delq(1) = 0.5*(rinf-qext(1))
        delq(2) = 0.5*(rinf*uind-qext(2))
        delq(3) = 0.5*(rinf*vind-qext(3))
        delq(4) = 0.5*(rinf*wind-qext(4))
        delq(5) = 0.5*(einf-qext(5))

        rho = qbar(1)
        rrho = 1./rho
        uu = qbar(2)*rrho
        vv = qbar(3)*rrho
        ww = qbar(4)*rrho
        e  = qbar(5)*rrho
        uvw = 0.5*(uu*uu+vv*vv+ww*ww)
        cjkl = sqrt(ggm1*(e-uvw))
        c2i = 1./(cjkl*cjkl)
        ge = gamma*e - gm1*uvw
        qq = (uu-ug(j,k,l))*zxn+(vv-vg(j,k,l))*zyn+(ww-wg(j,k,l))*zzn
        
        bSq = Mp**2/(bt(j,k,l)-Mp**2*(bt(j,k,l)-1))
        
        Z1 = sqrt((1.-bSq)*qq*(1.-bSq)*qq+4.*bSq*cjkl*cjkl)
        Zi = 1./Z1
        R = 0.5*((1.-bSq)*qq+Z1)
        S = 0.5*((1.-bSq)*qq-Z1)

c..multiplying by sign(lamda_pa)*inv(X_pa)
        delw(1) = (zxn*(1.-uvw*gm1*c2i)-(vv*zzn-ww*zyn)*rrho)*
     c                                    sign(1.,qq)*delq(1)
        delw(1) = delw(1) +(uu*zxn*c2i*gm1)*sign(1.,qq)*delq(2)
        delw(1) = delw(1) +(vv*zxn*c2i*gm1+zzn*rrho)*sign(1.,qq)*delq(3)
        delw(1) = delw(1) +(ww*zxn*c2i*gm1-zyn*rrho)*sign(1.,qq)*delq(4)
        delw(1) = delw(1) - zxn*gm1*c2i*sign(1.,qq)*delq(5)
       
        delw(2) = (zyn*(1.-uvw*gm1*c2i)-(ww*zxn-uu*zzn)*rrho)*
     c                                    sign(1.,qq)*delq(1)
        delw(2) = delw(2) +(uu*zyn*c2i*gm1-zzn*rrho)*sign(1.,qq)*delq(2)
        delw(2) = delw(2) +(vv*zyn*c2i*gm1)*sign(1.,qq)*delq(3)
        delw(2) = delw(2) +(ww*zyn*c2i*gm1+zxn*rrho)*sign(1.,qq)*delq(4)
        delw(2) = delw(2) - zyn*gm1*c2i*sign(1.,qq)*delq(5)

        delw(3) = (zzn*(1.-uvw*gm1*c2i)-(uu*zyn-vv*zxn)*rrho)*
     c                                    sign(1.,qq)*delq(1)
        delw(3) = delw(3) +(uu*zzn*c2i*gm1+zyn*rrho)*sign(1.,qq)*delq(2)
        delw(3) = delw(3) +(vv*zzn*c2i*gm1-zxn*rrho)*sign(1.,qq)*delq(3)
        delw(3) = delw(3) +(ww*zzn*c2i*gm1)*sign(1.,qq)*delq(4)
        delw(3) = delw(3) - zzn*gm1*c2i*sign(1.,qq)*delq(5)
       
        delw(4) = (-S*uvw*gm1*c2i-qq*bSq)*Zi*sign(1.,R+qq*bSq)*delq(1)
        delw(4) = delw(4) + (zxn*bSq + S*uu*gm1*c2i)*Zi*
     c                      sign(1.,R+qq*bSq)*delq(2)
        delw(4) = delw(4) + (zyn*bSq + S*vv*gm1*c2i)*Zi*
     c                      sign(1.,R+qq*bSq)*delq(3)
        delw(4) = delw(4) + (zzn*bSq + S*ww*gm1*c2i)*Zi*
     c                      sign(1.,R+qq*bSq)*delq(4)
        delw(4) = delw(4) - S*gm1*c2i*Zi*sign(1.,R+qq*bSq)*delq(5)

        delw(5) = (R*uvw*gm1*c2i+qq*bSq)*Zi*sign(1.,S+qq*bSq)*delq(1)
        delw(5) = delw(5) - (zxn*bSq + R*uu*gm1*c2i)*Zi*
     c                      sign(1.,S+qq*bSq)*delq(2)
        delw(5) = delw(5) - (zyn*bSq + R*vv*gm1*c2i)*Zi*
     c                      sign(1.,S+qq*bSq)*delq(3)
        delw(5) = delw(5) - (zzn*bSq + R*ww*gm1*c2i)*Zi*
     c                      sign(1.,S+qq*bSq)*delq(4)
        delw(5) = delw(5) +  R*gm1*c2i*Zi*sign(1.,S+qq*bSq)*delq(5)

c.. multiplying by (X_pa)

        delq(1) = zxn*delw(1)+zyn*delw(2)+zzn*delw(3)+delw(4)+delw(5)

        delq(2) = uu*zxn*delw(1)+(uu*zyn-zzn*rho)*delw(2)
        delq(2) = delq(2) + (uu*zzn+zyn*rho)*delw(3)
        delq(2) = delq(2) + (uu+R*zxn/bSq)*delw(4)
        delq(2) = delq(2) + (uu+S*zxn/bSq)*delw(5)

        delq(3) = (vv*zxn+zzn*rho)*delw(1)+vv*zyn*delw(2)
        delq(3) = delq(3) + (vv*zzn-zxn*rho)*delw(3)
        delq(3) = delq(3) + (vv+R*zyn/bSq)*delw(4)
        delq(3) = delq(3) + (vv+S*zyn/bSq)*delw(5)

        delq(4) = (ww*zxn-zyn*rho)*delw(1)+(ww*zyn+zxn*rho)*delw(2)
        delq(4) = delq(4) + ww*zzn*delw(3)
        delq(4) = delq(4) + (ww+R*zzn/bSq)*delw(4)
        delq(4) = delq(4) + (ww+S*zzn/bSq)*delw(5)

        delq(5) = (uvw*zxn+rho*(vv*zzn-ww*zyn))*delw(1)
        delq(5) = delq(5) + (uvw*zyn+rho*(ww*zxn-uu*zzn))*delw(2)
        delq(5) = delq(5) + (uvw*zzn+rho*(uu*zyn-vv*zxn))*delw(3)
        delq(5) = delq(5) + (ge+R*qq/bSq)*delw(4)+(ge+S*qq/bSq)*delw(5)

c.. qbar - delq
        q(j,k,l,1) = (qbar(1)-delq(1))/q(j,k,l,6)
        q(j,k,l,2) = (qbar(2)-delq(2))/q(j,k,l,6)
        q(j,k,l,3) = (qbar(3)-delq(3))/q(j,k,l,6)
        q(j,k,l,4) = (qbar(4)-delq(4))/q(j,k,l,6)
        q(j,k,l,5) = (qbar(5)-delq(5))/q(j,k,l,6)

      enddo
      enddo

      endif

      return
      end

! Math: added periodic BC (eg for O mesh overlap)
c***********************************************************************
	subroutine bcperiodic(q,js,je,ks,ke,ls,le,idir)
c***********************************************************************
      use params_global
      implicit none

      real		:: q(jmax,kmax,lmax,nd)
      integer	:: js,je,ks,ke,ls,le,idir

c.. local variables

      integer	:: j,k,l,j1,k1,l1,n
      integer	:: iadd,iadir,n_periodic!,add
      real		:: scale,press

c**** first executable statement

      iadd	= sign(1,idir)
      iadir = abs(idir)
	  
      if (iadir.eq.1) then
		  
		 n_periodic = je - js + 1
 
		 if (iadd .eq. 1) then
			 ! Minimum Boundary
			 !add	= 0
			 do j = js, je 
				 n = j - js
				 j1 = jmax - 2*n_periodic + n
				 !j1 = jmax - n_periodic -1 + add
				 do l = ls, le
					 do k = ks, ke
						 scale		= q(j1,k,l,6)/q(j,k,l,6)
						 q(j,k,l,1) = q(j1,k,l,1)*scale
						 q(j,k,l,2) = q(j1,k,l,2)*scale
						 q(j,k,l,3) = q(j1,k,l,3)*scale
						 q(j,k,l,4) = q(j1,k,l,4)*scale
						 q(j,k,l,5) = q(j1,k,l,5)*scale
					 enddo
				  enddo
				  !add = add + 1
			 enddo
		  else
			  ! Maximum Boundary
			  !add	= 0
			  do j = js, je
				  n = je - j
				  j1 = 1 + 2*n_periodic - n 
				  !j1 = n_periodic + 1 + add
				  do l = ls, le
					  do k = ks, ke
						 scale		= q(j1,k,l,6)/q(j,k,l,6)
						 q(j,k,l,1) = q(j1,k,l,1)*scale
						 q(j,k,l,2) = q(j1,k,l,2)*scale
						 q(j,k,l,3) = q(j1,k,l,3)*scale
						 q(j,k,l,4) = q(j1,k,l,4)*scale
						 q(j,k,l,5) = q(j1,k,l,5)*scale
					  enddo
				   enddo
				   !add = add + 1
			  enddo
		  endif		
		
      elseif (iadir.eq.2) then
	  
		 n_periodic = ke - ks + 1
		 if (iadd .eq. 1) then
			 ! Minimum Boundary
			 !add	= 0
			 do k = ks, ke 
				 n = k - ks
				 k1 = kmax - 2*n_periodic + n
				 !k1 = kmax - n_periodic -1 + add
				 do l = ls, le
					 do j = js, je
						 scale		= q(j,k1,l,6)/q(j,k,l,6)
						 q(j,k,l,1) = q(j,k1,l,1)*scale
						 q(j,k,l,2) = q(j,k1,l,2)*scale
						 q(j,k,l,3) = q(j,k1,l,3)*scale
						 q(j,k,l,4) = q(j,k1,l,4)*scale
						 q(j,k,l,5) = q(j,k1,l,5)*scale
					 enddo
				  enddo
				  !add = add + 1
			 enddo
		  else
			  ! Maximum Boundary
			  !add	= 0
			  do k = ks, ke
				  n = ke - k
				  k1 = 1 + 2*n_periodic - n 
				  !k1 = n_periodic + 1 + add
				  do l = ls, le
					  do j = js, je
						 scale		= q(j,k1,l,6)/q(j,k,l,6)
						 q(j,k,l,1) = q(j,k1,l,1)*scale
						 q(j,k,l,2) = q(j,k1,l,2)*scale
						 q(j,k,l,3) = q(j,k1,l,3)*scale
						 q(j,k,l,4) = q(j,k1,l,4)*scale
						 q(j,k,l,5) = q(j,k1,l,5)*scale
					  enddo
				   enddo
				   !add = add + 1
			  enddo
		  endif		
      elseif (iadir.eq.3) then

		 n_periodic = le - ls + 1
		 if (iadd .eq. 1) then
			 ! Minimum Boundary
			 !add	= 0
			 do l = ls, le 
				 n = l - ls
				 l1 = lmax - 2*n_periodic + n
				 !l1 = lmax - n_periodic -1 + add
				 do j = js, je
					 do k = ks, ke
						 scale		= q(j,k,l1,6)/q(j,k,l,6)
						 q(j,k,l,1) = q(j,k,l1,1)*scale
						 q(j,k,l,2) = q(j,k,l1,2)*scale
						 q(j,k,l,3) = q(j,k,l1,3)*scale
						 q(j,k,l,4) = q(j,k,l1,4)*scale
						 q(j,k,l,5) = q(j,k,l1,5)*scale
					 enddo
				  enddo
				  !add = add + 1
			 enddo
		  else
			  ! Maximum Boundary
			  !add	= 0
			  do l = ls, le
				  n = le - l
				  l1 = 1 + 2*n_periodic - n 
				  !l1 = n_periodic + 1 + add
				  do j = js, je
					  do k = ks, ke
						 scale		= q(j,k,l1,6)/q(j,k,l,6)
						 q(j,k,l,1) = q(j,k,l1,1)*scale
						 q(j,k,l,2) = q(j,k,l1,2)*scale
						 q(j,k,l,3) = q(j,k,l1,3)*scale
						 q(j,k,l,4) = q(j,k,l1,4)*scale
						 q(j,k,l,5) = q(j,k,l1,5)*scale
					  enddo
				   enddo
				   !add = add + 1
			  enddo
		  endif		

      endif
c
      return
	  
      end
c**********************************************************************
      subroutine bcaxisym(q,js,je,ks,ke,ls,le,idir)
c
c..periodic boundary for overlapping planes
c
c***********************************************************************

      use params_global

c*******************  list of global variables used ******************
c
c     jmax,kmax,lmax,nd,ksym,gm1,invisc,half
c
c***********************************************************************

      implicit none

      real q(jmax,kmax,lmax,nd)
      integer js,je,ks,ke,ls,le,idir

c.. local variables

c      integer j,k,l,j1,j2,j3,jc1,jc2,jc3
      integer j,k,l,jc,jj,jj1
      integer iadd
      real tcos,tsin
      
      iadd = sign(1,idir)

      pi = 4*atan(1.0)
      blang = 2*pi/float(nblade_tot)
      tcos  = cos(blang)
      tsin  = sin(blang)
      
      jj  = je - js + 1
      if(idir.eq.1) then
      do j = js,je
         jj1 = j - js
         jc = jmax - 2*jj + jj1
 
         do k = ks,ke
         do l = ls,le
           q(j,k,l,1) = q(jc,k,l,1)
           q(j,k,l,2) = q(jc,k,l,2)*tcos - q(jc,k,l,3)*tsin
           q(j,k,l,3) = q(jc,k,l,2)*tsin + q(jc,k,l,3)*tcos
           q(j,k,l,4) = q(jc,k,l,4)
           q(j,k,l,5) = q(jc,k,l,5)
         enddo
         enddo
      enddo

      elseif(idir.eq.-1) then
      do j = js,je
         jj1 = je - j
         jc = 1 + 2*jj - jj1 
 
         do k = ks,ke
         do l = ls,le
           q(j,k,l,1) = q(jc,k,l,1)
           q(j,k,l,2) = q(jc,k,l,2)*tcos + q(jc,k,l,3)*tsin
           q(j,k,l,3) =-q(jc,k,l,2)*tsin + q(jc,k,l,3)*tcos
           q(j,k,l,4) = q(jc,k,l,4)
           q(j,k,l,5) = q(jc,k,l,5)
         enddo
         enddo
      enddo
      else 
         print*,'idir = ',idir,' is not implemented in bcaxisym'
      endif

c      if(idir.eq.1) then
c      j1  = js
c      j2  = j1 + iadd
c      j3  = j2 + iadd
c      jc1 = jmax-6*iadd 
c      jc2 = jmax-5*iadd 
c      jc3 = jmax-4*iadd 
c 
c      do k = ks,ke
c      do l = ls,le
c        q(j1,k,l,1) = q(jc1,k,l,1)
c        q(j1,k,l,2) = q(jc1,k,l,2)*tcos - q(jc1,k,l,3)*tsin
c        q(j1,k,l,3) = q(jc1,k,l,2)*tsin + q(jc1,k,l,3)*tcos
c        q(j1,k,l,4) = q(jc1,k,l,4)
c        q(j1,k,l,5) = q(jc1,k,l,5)
c        q(j2,k,l,1) = q(jc2,k,l,1)
c        q(j2,k,l,2) = q(jc2,k,l,2)*tcos - q(jc2,k,l,3)*tsin
c        q(j2,k,l,3) = q(jc2,k,l,2)*tsin + q(jc2,k,l,3)*tcos
c        q(j2,k,l,4) = q(jc2,k,l,4)
c        q(j2,k,l,5) = q(jc2,k,l,5)
c        q(j3,k,l,1) = q(jc3,k,l,1)
c        q(j3,k,l,2) = q(jc3,k,l,2)*tcos - q(jc3,k,l,3)*tsin
c        q(j3,k,l,3) = q(jc3,k,l,2)*tsin + q(jc3,k,l,3)*tcos
c        q(j3,k,l,4) = q(jc3,k,l,4)
c        q(j3,k,l,5) = q(jc3,k,l,5)
c      enddo
c      enddo
c
c      elseif (idir.eq.-1) then 
c      j1  = js
c      j2  = j1 + iadd
c      j3  = j2 + iadd
c      jc1 = 7 
c      jc2 = 6 
c      jc3 = 5 
c 
c      do k = ks,ke
c      do l = ls,le
c        q(j1,k,l,1) = q(jc1,k,l,1)
c        q(j1,k,l,2) = q(jc1,k,l,2)*tcos - q(jc1,k,l,3)*tsin
c        q(j1,k,l,3) = q(jc1,k,l,2)*tsin + q(jc1,k,l,3)*tcos
c        q(j1,k,l,4) = q(jc1,k,l,4)
c        q(j1,k,l,5) = q(jc1,k,l,5)
c        q(j2,k,l,1) = q(jc2,k,l,1)
c        q(j2,k,l,2) = q(jc2,k,l,2)*tcos - q(jc2,k,l,3)*tsin
c        q(j2,k,l,3) = q(jc2,k,l,2)*tsin + q(jc2,k,l,3)*tcos
c        q(j2,k,l,4) = q(jc2,k,l,4)
c        q(j2,k,l,5) = q(jc2,k,l,5)
c        q(j3,k,l,1) = q(jc3,k,l,1)
c        q(j3,k,l,2) = q(jc3,k,l,2)*tcos - q(jc3,k,l,3)*tsin
c        q(j3,k,l,3) = q(jc3,k,l,2)*tsin + q(jc3,k,l,3)*tcos
c        q(j3,k,l,4) = q(jc3,k,l,4)
c        q(j3,k,l,5) = q(jc3,k,l,5)
c      enddo
c      enddo
c      endif

      return
      end

c***********************************************************************
      subroutine bck2d(q)
c
c*** Prologue : ************
c  2-d airfoil bc (kmax=3)
c
c  as is from turnsv1.50 except for f90 changes
c  last updated jaina 07/13/04
c***********************************************************************

      use params_global

c*******************  list of global variables used ********************
c
c     jmax,kmax,lmax,nd
c***********************************************************************

      implicit none

      real q(jmax,kmax,lmax,nd)
      
c..   local variables
      
      integer j,l

c**** first executable statement

      do l = 1,lmax
         do j = 1,jmax
            q(j,3,l,1) = q(j,2,l,1)
            q(j,3,l,2) = q(j,2,l,2)
            q(j,3,l,3) = 0.0
            q(j,3,l,4) = q(j,2,l,4)
            q(j,3,l,5) = q(j,2,l,5)
            q(j,2,l,3) = 0.0
            q(j,1,l,1) = q(j,2,l,1)
            q(j,1,l,2) = q(j,2,l,2)
            q(j,1,l,3) = 0.0
            q(j,1,l,4) = q(j,2,l,4)
            q(j,1,l,5) = q(j,2,l,5)
         enddo
      enddo

      return
      end


c*********************************************************************
      subroutine bckhal(q,x,y,z)

c**** Prologue : *************
c     if (half.eq.0), i.e. for symmetric airfoil simulations
c     obsolete now, havent tested a case in few years
c     for half of cgrid   
c     treatment for beyond the tip and the trailing wake at l = 1
c     
c  as is from turnsv1.50 except for f90 changes
c     last updated 07/13/04 by jaina
c*********************************************************************

      use params_global

c*******************  list of global variables used ********************
c
c     jmax,kmax,lmax,nd,ktip,jtail1,kroot
c***********************************************************************

      implicit none

      real q(jmax,kmax,lmax,nd)
      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)
      
c..   local variables

      integer j,k,l
      integer l1,jj,j1,j2
      real scale1,scale2,h1,h2,f1,f2

c**** first executable statement

      l  = 1
      l1 = l + 1

c..   set w-velocity equal to zero
c..   beyond the tip zeroeth order extrapolation of other variables
      
      if(kmax.gt.3) then

        k=ktip
        do 10 j = jtail1, jmax
          q(j,k,l,4) = 0.0
 10    continue

       do 25 k = ktip, kmax
          do 20 j = 1, jmax
          scale1 = q(j,k,l1,6)/q(j,k,l,6)
          q(j,k,l,1) = q(j,k,l1,1)*scale1
          q(j,k,l,2) = q(j,k,l1,2)*scale1
          q(j,k,l,3) = q(j,k,l1,3)*scale1
          q(j,k,l,4) = 0.0
          q(j,k,l,5) = q(j,k,l1,5)*scale1
   20   continue

c..   average at leading edge

        j = jmax-1
        jj = j-1
        scale1= q(jj,k,l,6)/q(j,k,l,6)
        h1 = sqrt( (x(j,k,l1)-x(j,k,l))**2 +
     <             (y(j,k,l1)-y(j,k,l))**2 +
     <             (z(j,k,l1)-z(j,k,l))**2 )
        h2 = sqrt( (x(j,k,l)-x(jj,k,l))**2 +
     <             (y(j,k,l)-y(jj,k,l))**2 +
     <             (z(j,k,l)-z(jj,k,l))**2 )
        f1 = h2/(h1+h2)
        f2 = 1.-f1
        q(j,k,l,1) = f1*q(j,k,l,1)+f2*q(jj,k,l,1)*scale1
        q(j,k,l,2) = f1*q(j,k,l,2)+f2*q(jj,k,l,2)*scale1
        q(j,k,l,3) = f1*q(j,k,l,3)+f2*q(jj,k,l,3)*scale1
        q(j,k,l,4) = f1*q(j,k,l,4)+f2*q(jj,k,l,4)*scale1
        q(j,k,l,5) = f1*q(j,k,l,5)+f2*q(jj,k,l,5)*scale1
   25   continue

      endif
  
c..   set w-velocity equal to zero
c..   in the wake zeroeth order extrapolation of other variables
      
      do k = kroot, ktip
         j = jtail1
         q(j,k,l,4) = 0.0
         do j = 1, jtail1-1
            scale1 = q(j,k,l1,6)/q(j,k,l,6)
            q(j,k,l,1) = q(j,k,l1,1)*scale1
            q(j,k,l,2) = q(j,k,l1,2)*scale1
            q(j,k,l,3) = q(j,k,l1,3)*scale1
            q(j,k,l,4) = 0.0
            q(j,k,l,5) = q(j,k,l1,5)*scale1
         enddo
      enddo
  
c..   in the front part
         
         j = jmax
         j1 = j-1
         j2 = j-2
         do  k = kroot, kmax
            do l = 1, lmax
               q(j1,k,l,4) = 0.0
               scale2 = q(j2,k,l,6)/q(j,k,l,6)
               q(j,k,l,1) = q(j2,k,l,1)*scale2
               q(j,k,l,2) = q(j2,k,l,2)*scale2
               q(j,k,l,3) = q(j2,k,l,3)*scale2
               q(j,k,l,4) = -q(j2,k,l,4)*scale2
               q(j,k,l,5) = q(j2,k,l,5)*scale2
            enddo
         enddo
         
            
      return
      end
      
c***********************************************************************
      subroutine bctany( q,xx,xy,xz,yx,yy,yz,zx,zy,zz,
     <                   ug,vg,wg,js,je,ks,ke,ls,le,idir)

c*** Prologue : *******
c
c  inviscid : tangency b.c. at a zeta equal constant surface
c  viscous  : no slip condition
c  presently, set up for surface at l = 1.
c
c  as is from turnsv1.50 except for f90 changes
c  last updated 07/13/04 by jaina
c***********************************************************************
      
      use params_global

c*******************  list of global variables used ******************
c
c     jmax,kmax,lmax,nd,mdim,invisc
c
c***********************************************************************
      implicit none

      real q(jmax,kmax,lmax,nd)
      real xx(jmax,kmax,lmax),xy(jmax,kmax,lmax),xz(jmax,kmax,lmax),
     <     yx(jmax,kmax,lmax),yy(jmax,kmax,lmax),yz(jmax,kmax,lmax),
     <     zx(jmax,kmax,lmax),zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax),vg(jmax,kmax,lmax),wg(jmax,kmax,lmax)
      real u(mdim),v(mdim),w(mdim)
      integer js,je,ks,ke,ls,le,idir,iadir

c..   local variables

      real u1,u2,u3,v1,v2,v3,w1,w2,w3,foso
      integer j,k,l,iadd
      integer j1,j2,j3,k1,k2,k3,l1,l2,l3,ihigh

c**** first executable statement      

      iadd = sign(1,idir)
      iadir = abs(idir)
      foso = 0.

      if(iadir.eq.1) then
        j  = js
        j1 = j + iadd
        j2 = j1 + iadd
        j3 = j2 + iadd
        ihigh = 0

        if( invisc ) then

c..for inviscid:
c..extrapolate contravariant velocities to the surface
c..solve for the surface cartesian velocity components
cjdb..note: really we are extrapolating the physical values!

          do l = ls,le
            do k = ks,ke
              v1 = (q(j1,k,l,2)*yx(j,k,l)+q(j1,k,l,3)*yy(j,k,l)
     $             +q(j1,k,l,4)*yz(j,k,l))/q(j1,k,l,1)
              v2 = (q(j2,k,l,2)*yx(j,k,l)+q(j2,k,l,3)*yy(j,k,l)
     $             +q(j2,k,l,4)*yz(j,k,l))/q(j2,k,l,1)
              v3 = (q(j3,k,l,2)*yx(j,k,l)+q(j3,k,l,3)*yy(j,k,l)
     $             +q(j3,k,l,4)*yz(j,k,l))/q(j3,k,l,1)
              w1 = (q(j1,k,l,2)*zx(j,k,l)+q(j1,k,l,3)*zy(j,k,l)
     $             +q(j1,k,l,4)*zz(j,k,l))/q(j1,k,l,1)
              w2 = (q(j2,k,l,2)*zx(j,k,l)+q(j2,k,l,3)*zy(j,k,l)
     $             +q(j2,k,l,4)*zz(j,k,l))/q(j2,k,l,1)
              w3 = (q(j3,k,l,2)*zx(j,k,l)+q(j3,k,l,3)*zy(j,k,l)
     $             +q(j3,k,l,4)*zz(j,k,l))/q(j3,k,l,1)
              v(k)  = (1.+foso)*v1 - foso*v2
              w(k)  = (1.+foso)*w1 - foso*w2
c...  higher order bc?
              if(ihigh.eq.1) v(k)  = 3.*v1-3.*v2+v3
              if(ihigh.eq.1) w(k)  = 3.*w1-3.*w2+w3
              u(k)  = ug(j,k,l)*xx(j,k,l)+vg(j,k,l)*xy(j,k,l)
     <             +wg(j,k,l)*xz(j,k,l)
            enddo
            call uvw(j,j,ks,ke,l,l,q,xx,xy,xz,yx,yy,yz,zx,zy,zz,
     $           u,v,w,2)
          enddo
c     
        else
c
c..for viscous flows the no slip condition at the surface requires
c  setting the contravariant velocities u = 0., v = 0., & w = 0. for
c  no suction or injection at the surface. this satisfies tangency
c  at the walls  
c..grs
c  solve for the surface cartesian velocity components
c
          do l = ls,le
            do k = ks,ke
              u(k) = ug(j,k,l)*xx(j,k,l) + vg(j,k,l)*xy(j,k,l)
     <             + wg(j,k,l)*xz(j,k,l)
              v(k) = ug(j,k,l)*yx(j,k,l) + vg(j,k,l)*yy(j,k,l)
     <             + wg(j,k,l)*yz(j,k,l)
              w(k) = ug(j,k,l)*zx(j,k,l) + vg(j,k,l)*zy(j,k,l)
     <             + wg(j,k,l)*zz(j,k,l)
            enddo
            call uvw(j,j,ks,ke,l,l,q,xx,xy,xz,yx,yy,yz,zx,zy,zz,
     $               u,v,w,2)
          enddo
        endif

      elseif(iadir.eq.2) then
        k  = ks
        k1 = k + iadd
        k2 = k1 + iadd
        k3 = k2 + iadd
        ihigh = 0

        if( invisc ) then

c..for inviscid:
c..extrapolate contravariant velocities to the surface
c..solve for the surface cartesian velocity components
cjdb..note: really we are extrapolating the physical values!

          do j = js,je
            do l = ls,le
              u1 = (q(j,k1,l,2)*xx(j,k,l)+q(j,k1,l,3)*xy(j,k,l)
     $             +q(j,k1,l,4)*xz(j,k,l))/q(j,k1,l,1)
              u2 = (q(j,k2,l,2)*xx(j,k,l)+q(j,k2,l,3)*xy(j,k,l)
     $             +q(j,k2,l,4)*xz(j,k,l))/q(j,k2,l,1)
              u3 = (q(j,k3,l,2)*xx(j,k,l)+q(j,k3,l,3)*xy(j,k,l)
     $             +q(j,k3,l,4)*xz(j,k,l))/q(j,k3,l,1)
              w1 = (q(j,k1,l,2)*zx(j,k,l)+q(j,k1,l,3)*zy(j,k,l)
     $             +q(j,k1,l,4)*zz(j,k,l))/q(j,k1,l,1)
              w2 = (q(j,k2,l,2)*zx(j,k,l)+q(j,k2,l,3)*zy(j,k,l)
     $             +q(j,k2,l,4)*zz(j,k,l))/q(j,k2,l,1)
              w3 = (q(j,k3,l,2)*zx(j,k,l)+q(j,k3,l,3)*zy(j,k,l)
     $             +q(j,k3,l,4)*zz(j,k,l))/q(j,k3,l,1)
              u(l)  = (1.+foso)*u1 - foso*u2
              w(l)  = (1.+foso)*w1 - foso*w2
c...  higher order bc?
              if(ihigh.eq.1) u(l)  = 3.*u1-3.*u2+u3
              if(ihigh.eq.1) w(l)  = 3.*w1-3.*v2+v3
              v(l)  = ug(j,k,l)*yx(j,k,l)+vg(j,k,l)*yy(j,k,l)
     <             +wg(j,k,l)*yz(j,k,l)
            enddo
            call uvw(j,j,k,k,ls,le,q,xx,xy,xz,yx,yy,yz,zx,zy,zz,
     $           u,v,w,3)
          enddo
c     
        else
c
c..for viscous flows the no slip condition at the surface requires
c  setting the contravariant velocities u = 0., v = 0., & w = 0. for
c  no suction or injection at the surface. this satisfies tangency
c  at the walls  
c..grs
c  solve for the surface cartesian velocity components
c
          do j = js,je
            do l = ls,le
              u(l) = ug(j,k,l)*xx(j,k,l) + vg(j,k,l)*xy(j,k,l)
     <             + wg(j,k,l)*xz(j,k,l)
              v(l) = ug(j,k,l)*yx(j,k,l) + vg(j,k,l)*yy(j,k,l)
     <             + wg(j,k,l)*yz(j,k,l)
              w(l) = ug(j,k,l)*zx(j,k,l) + vg(j,k,l)*zy(j,k,l)
     <             + wg(j,k,l)*zz(j,k,l)
            enddo
            call uvw(j,j,k,k,ls,le,q,xx,xy,xz,yx,yy,yz,zx,zy,zz,
     $               u,v,w,3)
          enddo
        endif

      elseif(iadir.eq.3) then
        l  = ls
        l1 = l + iadd
        l2 = l1 + iadd
        l3 = l2 + iadd
        ihigh = 0

        if( invisc ) then

c..for inviscid:
c..extrapolate contravariant velocities to the surface
c..solve for the surface cartesian velocity components
cjdb..note: really we are extrapolating the physical values!

          do k = ks,ke
            do j = js,je
               u1 = (q(j,k,l1,2)*xx(j,k,l)+q(j,k,l1,3)*xy(j,k,l)
     $              +q(j,k,l1,4)*xz(j,k,l))/q(j,k,l1,1)
               u2 = (q(j,k,l2,2)*xx(j,k,l)+q(j,k,l2,3)*xy(j,k,l)
     $              +q(j,k,l2,4)*xz(j,k,l))/q(j,k,l2,1)
               u3 = (q(j,k,l3,2)*xx(j,k,l)+q(j,k,l3,3)*xy(j,k,l)
     $              +q(j,k,l3,4)*xz(j,k,l))/q(j,k,l3,1)
               v1 = (q(j,k,l1,2)*yx(j,k,l)+q(j,k,l1,3)*yy(j,k,l)
     $              +q(j,k,l1,4)*yz(j,k,l))/q(j,k,l1,1)
               v2 = (q(j,k,l2,2)*yx(j,k,l)+q(j,k,l2,3)*yy(j,k,l)
     $              +q(j,k,l2,4)*yz(j,k,l))/q(j,k,l2,1)
               v3 = (q(j,k,l3,2)*yx(j,k,l)+q(j,k,l3,3)*yy(j,k,l)
     $              +q(j,k,l3,4)*yz(j,k,l))/q(j,k,l3,1)
               u(j)  = 2.0*u1 - u2
               v(j)  = 2.0*v1 - v2
               u(j)  = u1 
               v(j)  = v1
c...  higher order bc?
               if(ihigh.eq.1) u(j)  = 3.*u1-3.*u2+u3
               if(ihigh.eq.1) v(j)  = 3.*v1-3.*v2+v3
               w(j)  = ug(j,k,l)*zx(j,k,l)+vg(j,k,l)*zy(j,k,l)
     <              +wg(j,k,l)*zz(j,k,l)
            enddo
            call uvw(js,je,k,k,l,l,q,xx,xy,xz,yx,yy,yz,zx,zy,zz,
     $           u,v,w,1)
          enddo
c     
        else
c
c..for viscous flows the no slip condition at the surface requires
c  setting the contravariant velocities u = 0., v = 0., & w = 0. for
c  no suction or injection at the surface. this satisfies tangency
c  at the walls  
c..grs
c  solve for the surface cartesian velocity components
c
          do k = ks, ke
            do j = js, je
              u(j) = ug(j,k,l)*xx(j,k,l) + vg(j,k,l)*xy(j,k,l)
     <             + wg(j,k,l)*xz(j,k,l)
              v(j) = ug(j,k,l)*yx(j,k,l) + vg(j,k,l)*yy(j,k,l)
     <             + wg(j,k,l)*yz(j,k,l)
              w(j) = ug(j,k,l)*zx(j,k,l) + vg(j,k,l)*zy(j,k,l)
     <             + wg(j,k,l)*zz(j,k,l)
            enddo
            call uvw(js,je,k,k,l,l,q,xx,xy,xz,yx,yy,yz,zx,zy,zz,
     $               u,v,w,1)
          enddo
        endif

      endif

      return
      end

c***********************************************************************
      subroutine uvw(js,je,ks,ke,ls,le,q,xx,xy,xz,yx,yy,yz,zx,
     $               zy,zz,u,v,w,idir)

c*** Prologue : ******
c
c  solve for the cartesian momemtum components (ru,rv,rw) from
c  the contravariant velocity components (u,v,w) along lines of j.
c     
c  note: consistent with metric calculations in metfv
c
c  as is from turnsv1.50 except for f90 changes
c  last updated 07/13/04 by jaina
c***********************************************************************

      use params_global

c*******************  list of global variables used ******************
c
c     jmax,kmax,lmax,nd,mdim
c
c***********************************************************************

      implicit none
            
      real q(jmax,kmax,lmax,nd)
      real xx(jmax,kmax,lmax),xy(jmax,kmax,lmax),xz(jmax,kmax,lmax),
     <     yx(jmax,kmax,lmax),yy(jmax,kmax,lmax),yz(jmax,kmax,lmax),
     <     zx(jmax,kmax,lmax),zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax),vg(jmax,kmax,lmax),wg(jmax,kmax,lmax)
      real u(mdim),v(mdim),w(mdim)

      integer js,je,ks,ke,ls,le,idir

c..   local variables

      integer j,k,l
      real t11,t12,t13,t21,t22,t23,t31,t32,t33,vol
      real uw,vw,ww,uwall,vwall,wwall
      real dj,djinv,epsilon

c**** first executable statement

      epsilon=1.d-30

      do j = js,je
      do k = ks,ke
      do l = ls,le
c
        t11   = (yy(j,k,l)*zz(j,k,l) -yz(j,k,l)*zy(j,k,l))
        t12   =-(xy(j,k,l)*zz(j,k,l) -xz(j,k,l)*zy(j,k,l))
        t13   = (xy(j,k,l)*yz(j,k,l) -xz(j,k,l)*yy(j,k,l))
c
        t21   =-(yx(j,k,l)*zz(j,k,l) -yz(j,k,l)*zx(j,k,l))
        t22   = (xx(j,k,l)*zz(j,k,l) -xz(j,k,l)*zx(j,k,l))
        t23   =-(xx(j,k,l)*yz(j,k,l) -xz(j,k,l)*yx(j,k,l))
c
        t31   = (yx(j,k,l)*zy(j,k,l) -yy(j,k,l)*zx(j,k,l))
        t32   =-(xx(j,k,l)*zy(j,k,l) -xy(j,k,l)*zx(j,k,l))
        t33   = (xx(j,k,l)*yy(j,k,l) -xy(j,k,l)*yx(j,k,l))
c
c..   according to karthik

        dj=xx(j,k,l)*t11+xy(j,k,l)*t21+xz(j,k,l)*t31
        djinv=sign(1.,dj)/(abs(dj) + epsilon)
        vol=q(j,k,l,1)*djinv
c        vol=q(j,k,l,1)/(xx(j,k,l)*t11+xy(j,k,l)*t21+xz(j,k,l)*t31)
c
c..the following are the physical plane velocities u, v,& w
c
        if (idir.eq.1) then
          uw = u(j)
          vw = v(j)
          ww = w(j)
        elseif (idir.eq.2) then
          uw = u(k)
          vw = v(k)
          ww = w(k)
        elseif (idir.eq.3) then
          uw = u(l)
          vw = v(l)
          ww = w(l)
        endif

        uwall = (t11*uw +t12*vw +t13*ww)*vol
        vwall = (t21*uw +t22*vw +t23*ww)*vol
        wwall = (t31*uw +t32*vw +t33*ww)*vol
c
c..physical plane velocities multiplied by the density to from q's
c
        q(j,k,l,2) = uwall
        q(j,k,l,3) = vwall
        q(j,k,l,4) = wwall

      enddo
      enddo
      enddo
c
      return
      end

c***********************************************************************
      subroutine bcsym(q,js,je,ks,ke,ls,le,idir)

c*** Prologue : ***********
c
c  3-d wing bc at symmetry plane (k=1)
c
c  as is from turnsv1.50 except for f90 changes
c  last updated 07/13/04 by jaina
c***********************************************************************

      use params_global

c*******************  list of global variables used ******************
c
c     jmax,kmax,lmax,nd,ksym,gm1,invisc,half
c
c***********************************************************************

      implicit none

      real q(jmax,kmax,lmax,nd)
      integer js,je,ks,ke,ls,le,idir

c.. local variables

      integer j,k,l,j1,k1,l1
      integer iadd,iadir
      real scale,press

c**** first executable statement

      iadd = sign(1,idir)
      iadir = abs(idir)
      if (iadir.eq.1) then
         j = js
         j1 = j + iadd
         do k = ks,ke
         do l = ls,le
           scale = q(j1,k,l,6)/q(j,k,l,6)
           q(j,k,l,1) = q(j1,k,l,1)*scale
           q(j,k,l,2) = 0.0 
           q(j,k,l,3) = q(j1,k,l,3)*scale
           q(j,k,l,4) = q(j1,k,l,4)*scale
           q(j,k,l,5) = q(j1,k,l,5)*scale
         enddo
         enddo
      elseif (iadir.eq.2) then
         k = ks
         k1 = k + iadd
         do l = ls,le
         do j = js,je
           scale = q(j,k1,l,6)/q(j,k,l,6)
           q(j,k,l,1) = q(j,k1,l,1)*scale
           q(j,k,l,2) = q(j,k1,l,2)*scale
           q(j,k,l,3) = 0.0
           q(j,k,l,4) = q(j,k1,l,4)*scale
           q(j,k,l,5) = q(j,k1,l,5)*scale
         enddo
         enddo
      elseif (iadir.eq.3) then
         l = ls
         l1 = l + iadd
         do j = js,je
         do k = ks,ke
           scale = q(j,k,l1,6)/q(j,k,l,6)
           q(j,k,l,1) = q(j,k,l1,1)*scale
           q(j,k,l,2) = q(j,k,l1,2)*scale
           q(j,k,l,3) = q(j,k,l1,3)*scale
           q(j,k,l,4) = 0.0     
           q(j,k,l,5) = q(j,k,l1,5)*scale
         enddo
         enddo
      endif
c
      return
      end

c***********************************************************************
      subroutine bcwing_old(q)

c*** Prologue : ***********
c
c  3-d wing bc at symmetry plane (k=1)
c
c  as is from turnsv1.50 except for f90 changes
c  last updated 07/13/04 by jaina
c***********************************************************************

      use params_global

c*******************  list of global variables used ******************
c
c     jmax,kmax,lmax,nd,ksym,gm1,invisc,half
c
c***********************************************************************

      implicit none

      real q(jmax,kmax,lmax,nd)

c.. local variables

      integer j,l
      real scale,press

c**** first executable statement
c
c..set v-velocity equal to zero
c..at k=1 zeroeth order extrapolation of other variables
c
      do l = 1,lmax
      do j = 1,jmax
        scale = q(j,2,l,6)/q(j,1,l,6)
        q(j,1,l,1) = q(j,2,l,1)*scale
        q(j,1,l,2) = q(j,2,l,2)*scale
        q(j,1,l,3) = 0.0
        q(j,1,l,4) = q(j,2,l,4)*scale
        press = gm1*(q(j,2,l,5) -0.5*(q(j,2,l,2)**2+q(j,2,l,3)**2
     $             + q(j,2,l,4)**2)/q(j,2,l,1))*q(j,2,l,6)
        q(j,1,l,5) = press/(gm1*q(j,1,l,6)) +0.5*(((q(j,1,l,2)**2)
     $           +q(j,1,l,3)**2) +q(j,1,l,4)**2)/q(j,1,l,1)
      enddo
      enddo
c
      return
      end



c***********************************************************************
      subroutine bcwp0( p,q,kkr,kkp,js,je,ks,ke,ls,le,idir )

c*** Prologue : **********
c
c  compute the total energy by using p1 = p2 in lieu of
c  the normal momentum equation for the pressure at the wall. 
c  the pressure and and density are used to compute the total energy. 
c  the present implementation assumes the wall is at l = 1.
c  subroutine by jdb
c
c  as is from turnsv1.50 except for f90 changes
c  last updated jaina 07/13/04
c***********************************************************************'

      use params_global

c*******************  list of global variables used ******************
c
c     jmax,kmax,lmax,nd,ksym,gm1,invisc,half
c
c***********************************************************************'

      implicit none

      real p(jmax,kmax,3), q(jmax,kmax,lmax,nd)
      integer kkr(kmax),kkp(kmax)
      integer js,je,ks,ke,ls,le,idir,iadir

c..   local variables

      integer j,k,l,jj,ll
      integer l1,l2,ka,kb,iadd
      real eps

c**** first executable statement

      iadd = sign(1,idir)
      iadir = abs(idir)

      eps=1.0e-20

      if(iadir.eq.1) then
         print*,'idir = ',idir,' is not implemented in bcwp0'
      elseif(iadir.eq.2) then
         print*,'idir = ',idir,' is not implemented in bcwp0'
      elseif(iadir.eq.3) then
        ka = ks+1 - ksym
        kb = ke-1 + ksym

c..   compute the pressure field at l = 1,3
      
        do l = ls,ls+2*iadd
           do k = ks,ke
              do j = js,je
                 p(j,k,l) = gm1*(q(j,k,l,5) 
     $                -0.5*(q(j,k,l,2)**2+q(j,k,l,3)**2
     $                + q(j,k,l,4)**2)/q(j,k,l,1))*q(j,k,l,6)
              enddo
           enddo
        enddo

        l  = ls
        l1 = l  + iadd
        l2 = l1 + iadd

c..   extrapolate pressure to the surface
     
        if( invisc ) then
           do k = ks,ke
              do j = js,je
                 p(j,k,l) = 4./3.*p(j,k,l1)-1./3.*p(j,k,l2)
              enddo
           enddo
        else
           do k = ks,ke
              do j = js,je
                 p(j,k,l) = p(j,k,l1)
              enddo
           enddo
        endif
           
c..now bcs

        if (half.ne.1) then
          do 101 k = ka,kb
            p(js,k,l) = 0.5*(p(js,k,l1)+p(je,k,l1))
            p(je,k,l) = 0.5*(p(js,k,l1)+p(je,k,l1))
 101      continue
        else
          do 102 k = ka,kb
            p(js,k,l) = p(js,k,l1)
            p(je,k,l) = p(je,k,l1)
 102      continue
        endif
        if(ka.eq.kb) then
          do 111 j=js,je
            p(j,ks,l) = p(j,ks+1,l)
            p(j,ke,l) = p(j,ke-1,l)
 111      continue
        else
          if(half.ne.1) then
            do 112 j=js,je
              jj = je + js - j
              p(j,ks,l) = p(j,ks+1,l)
              p(j,ke,l) = 0.5*(p(j,ke,l1)+p(jj,ke,l1))
 112        continue
          else
            do 113 j=js,je
              p(j,ks,l) = p(j,ks+1,l)
              p(j,ke,l) = p(j,ke,l1)
 113        continue
          endif
        endif

c..update the total energy

        do k = ks,ke
           do j = js,je
              q(j,k,l,5) = p(j,k,l)/(gm1*q(j,k,l,6)) +0.5*(q(j,k,l,2)**2
     $             +q(j,k,l,3)**2 +q(j,k,l,4)**2)/q(j,k,l,1)
           enddo
        enddo
      
      endif
 
      return
      end

c***********************************************************************
      subroutine bcwp(a,b,rhs,p,q,xx,xy,xz,yx,yy,yz,zx,zy,zz,
     $     zx0,zy0,zz0,zt0,kkr,kkp,ug,vg,wg,js,je,ks,ke,ls,le,idir)

c*** Prologue **************
c
c  compute the total energy by solving the normal momentum
c  equation for the pressure at the wall. the pressure and
c  and density are used to compute the total energy. the wall
c  is at l = 1.
c
c  note that the wall pressure is computed only in the interior
c  of the l = constant surface.
c
c  as is from turnsv1.50 except for f90 changes
c  last updated jaina 07/13/04
c***********************************************************************

      use params_global

c*******************  list of global variables used ******************
c
c     jmax,kmax,lmax,nd,ksym,dx1,dy1,dz1,gm1,dt,iunst,rf,half
c
c***********************************************************************
      implicit none

      real p(jmax,kmax,3),q(jmax,kmax,lmax,nd)
      
      real xx(jmax,kmax,lmax),xy(jmax,kmax,lmax),xz(jmax,kmax,lmax),
     <     yx(jmax,kmax,lmax),yy(jmax,kmax,lmax),yz(jmax,kmax,lmax),
     <     zx(jmax,kmax,lmax),zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax),vg(jmax,kmax,lmax),wg(jmax,kmax,lmax)

      real a(jmax,kmax),b(jmax,kmax),rhs(jmax,kmax)
      real zx0(jmax,kmax),zy0(jmax,kmax),zz0(jmax,kmax),zt0(jmax,kmax)
      integer kkr(kmax),kkp(kmax)
      integer js,je,ks,ke,ls,le,idir

c..   local variables

      real eps
      integer ka,kb,ja,jb,jam1,jbp1,kam1,kbp1,kp,kr
      integer l,l1,l2,i,j,j1,jj,ll,lp,k,k1,kk,iadd,iadir
      
      real dx2,dy2,dz2,foso,cross,isign
      real drj1,drj2,drk1,drk2,ru,rv
      real dux,dvx,dwx,duy,dvy,dwy,dtinv
      real dzdx,dzdy,dzdz
      real rzt1,rzt2,rzt3,rzt4,rztt,beta

c**** first executable statement
      iadd = sign(1,idir)
      iadir = abs(idir)
      eps=1.0e-20

      if(iadir.eq.1) then
         print*,'idir = ',idir,' is not implemented in bcwp'
      elseif(iadir.eq.2) then
         print*,'idir = ',idir,' is not implemented in bcwp'
      elseif(iadir.eq.3) then
        l = ls
        l1 = l + iadd
        l2 = l1+ iadd

        ka = ks+1 - ksym
        kb = ke-1 + ksym
        ja = js
        jb = je

        jam1 = ja - 1
        jbp1 = jb + 1
        kam1 = ka - 1
        kbp1 = kb + 1

        dx2 = 2.0*dx1
        dy2 = 2.0*dy1
        dz2 = 2.0*dz1

c..compute pressure field

        lp = 0
        do ll = ls,ls+2*iadd
        lp = lp + 1
        do k = kam1,kbp1
        do  j = jam1,jbp1
           p(j,k,lp)  = gm1*(q(j,k,ll,5) -0.5*(q(j,k,ll,2)**2
     &          +q(j,k,ll,3)**2+q(j,k,ll,4)**2)/q(j,k,ll,1))*q(j,k,ll,6)
        enddo
        enddo
        enddo

c..set up rhs in p(j,k,1) and compute difference coefficients a and b

        foso = 0.
        do 24 k = ka,kb
           
        kp = kkp(k)
        kr = kkr(k)
        
        do 21 j = ja,jb
          drj1 = 1./q(j+1,k,l,1)
          drj2 = 1./q(j-1,k,l,1)
          drk1 = 1./q(j,kp,l,1)
          drk2 = 1./q(j,kr,l,1)
          ru   = (((q(j,k,l,2)*xx(j,k,l)) +q(j,k,l,3)*xy(j,k,l)) +
     $              q(j,k,l,4)*xz(j,k,l))*q(j,k,l,6)
          rv   = (((q(j,k,l,2)*yx(j,k,l)) +q(j,k,l,3)*yy(j,k,l)) +
     $         q(j,k,l,4)*yz(j,k,l))*q(j,k,l,6)

          dux  = q(j+1,k,l,2)*drj1 -q(j-1,k,l,2)*drj2
          dvx  = q(j+1,k,l,3)*drj1 -q(j-1,k,l,3)*drj2
          dwx  = q(j+1,k,l,4)*drj1 -q(j-1,k,l,4)*drj2
          
          duy  = q(j,kp,l,2)*drk1 -q(j,kr,l,2)*drk2
          dvy  = q(j,kp,l,3)*drk1 -q(j,kr,l,3)*drk2
          dwy  = q(j,kp,l,4)*drk1 -q(j,kr,l,4)*drk2
        
c..   we have to account for the unsteady motion of the blade
c..   the following is the first term on the right hand side of
c..   eq. (11) of pulliam and stegers aiaa j. paper, feb. 1980.
c..   for quasi-steady rztt = 0.
c..   for unsteady flow rztt is non-zero

          dtinv = 1.0/dt

          if(iunst.eq. FL_HOVER) then
            rzt1 = 0.
            rzt2 = -rf*zy(j,k,l)*q(j,k,l,2)
            rzt3 = +rf*zx(j,k,l)*q(j,k,l,3)
            rzt4 = 0.
          elseif(iunst.eq. FL_QUASI) then
            rzt1 = 0.
            rzt2 = 0.
            rzt3 = 0.
            rzt4 = 0.
c..        10/18/2013  added by camli acceleration term=0 
c..          elseif(iunst.eq. FL_FLAP) then
c..            rzt1 = 0.
c..            rzt2 = ((zx(j,k,l) - zx0(j,k))/dt) * q(j,k,l,2)
c..            rzt3 = ((zy(j,k,l) - zy0(j,k))/dt) * q(j,k,l,3)
c..            rzt4 = ((zz(j,k,l) - zz0(j,k))/dt) * q(j,k,l,4)
c..  done editing        
          else
            rzt1 = ((-ug(j,k,l)*zx(j,k,l)-vg(j,k,l)*zy(j,k,l)
     <               -wg(j,k,l)*zz(j,k,l) - zt0(j,k))/dt) * q(j,k,l,1)
            rzt2 = ((zx(j,k,l) - zx0(j,k))/dt) * q(j,k,l,2)
            rzt3 = ((zy(j,k,l) - zy0(j,k))/dt) * q(j,k,l,3)
            rzt4 = ((zz(j,k,l) - zz0(j,k))/dt) * q(j,k,l,4)
          endif

          rztt = rzt1 + rzt2 + rzt3 + rzt4
          rztt = rztt * q(j,k,l,6)

          rhs(j,k)   = ( -zx(j,k,l)*(ru*dux +rv*duy)
     $                   -zy(j,k,l)*(ru*dvx +rv*dvy)
     $                   -zz(j,k,l)*(ru*dwx +rv*dwy) )*.5 + rztt

   21   continue

        rhs(jam1,k) = 0.0
        rhs(jbp1,k) = 0.0
        
        do 22 j = jam1,jbp1

           dzdx       = (((zx(j,k,l)*xx(j,k,l)) +zy(j,k,l)*xy(j,k,l))
     $          +zz(j,k,l)*xz(j,k,l)) +eps
           dzdy       = (((zx(j,k,l)*yx(j,k,l)) +zy(j,k,l)*yy(j,k,l))
     $          +zz(j,k,l)*yz(j,k,l)) +eps
           dzdz       = (((zx(j,k,l)*zx(j,k,l)) +zy(j,k,l)*zy(j,k,l))
     $                  +zz(j,k,l)*zz(j,k,l)) +eps
           rhs(j,k)   = -((dz2*rhs(j,k)/dzdz) +( -4.0*p(j,k,2)
     $               +p(j,k,3)) )/3.0
           a(j,k)     = dz2*dzdx/(3.0*dzdz*dx2)
           b(j,k)     = dz2*dzdy/(3.0*dzdz*dy2)

 22     continue

        rhs(jam1,k)   = 0.0
        rhs(jbp1,k) = 0.0

        do 23 j = ja,jb
           cross      = a(j,k)*b(j,k)*( p(j+1,kp,1) -p(j+1,kr,1)
     $          -p(j-1,kp,1) +p(j-1,kr,1) )
           rhs(j,k)   = rhs(j,k) +cross
 23     continue
 24     continue
      
c..solve for pstar
c..set boundary conditions

        do 32 i = 1,2

        if( i.eq.1 ) then
           j          = jam1
           j1         = j+1
           isign       = 1.0
        else
           j          = jbp1
           j1          = j-1
           isign       = -1.0
        end if

        do 31 k = ka,kb
           kp = kkp(k)
           kr = kkr(k)
           rhs(j1,k)  = -isign*a(j1,k)*( p(j,k,1) -b(j,k)*(p(j,kp,1) -
     $          p(j,kr,1)) ) +rhs(j1,k)
 31     continue
 32     continue
        
c..     invert xi operator
        
        do 34 j = ja+1,jb
           do 34 k = ka,kb
              beta     = (1.0) +a(j,k)*a(j-1,k)
              rhs(j,k) = (rhs(j,k) -a(j,k)*rhs(j-1,k))/beta
              a(j,k)   = a(j,k)/beta
 34        continue
           do 35 jj  = ja+1,jb
              j = ja + jb -jj
              do 35 k = ka,kb
                 rhs(j,k)   = rhs(j,k) +a(j,k)*rhs(j+1,k)
 35           continue
            
c..solve for p on surface
c..non-periodic adi solver 

        do 42 i = 1,2
           if( i.eq.1) then
            k          = kam1
            k1         = k+1
            isign       = 1.0
        else
            k          = kbp1
            k1         = k-1
            isign       = -1.0
        end if

        do 41 j = ja,jb
        rhs(j,k1)  = -isign*b(j,k1)*p(j,k,1) +rhs(j,k1)
   41   continue
   42   continue

c..invert eta operator
     
        if(kmax.ge.4) then
        do 44 k = ka+1,kb
        do 44 j = ja,jb
          beta       = (1.0) +b(j,k)*b(j,k-1)
          rhs(j,k)   = (rhs(j,k) -b(j,k)*rhs(j,k-1))/beta
          b(j,k)     = b(j,k)/beta
   44   continue

        do 45 kk = ka+1,kb
          k = ka + kb -kk
          do 45 j = ja,jb
            rhs(j,k)   = rhs(j,k) +b(j,k)*rhs(j,k+1)
   45   continue
        endif

        do 46 k = ka,kb
        do 46 j = ja,jb
          p(j,k,1)   = rhs(j,k)
   46   continue
     
c..now bcs
     
        if(half.ne.1) then
          do 101 k = ka,kb
            p(js,k,1) = 0.5*(p(js,k,2)+p(je,k,2))
            p(je,k,1) = 0.5*(p(js,k,2)+p(je,k,2))
 101      continue
        else
          do 102 k = ka,kb
            p(js,k,1) = p(js,k,2)
            p(je,k,1) = p(je,k,2)
 102      continue
        endif
        if(ka.eq.kb) then
          do 111 j=js,je
            p(j,ks,1) = p(j,ks+1,1)
            p(j,ke,1) = p(j,ke-1,1)
 111      continue
        else
          if(half.ne.1) then
            do 112 j=js,je
              jj = je + js - j
              p(j,ks,1) = p(j,ks+1,1)
              p(j,ke,1) = 0.5*(p(j,ke,2)+p(jj,ke,2))
 112        continue
          else
            do 113 j=js,je
              p(j,ks,1) = p(j,ks+1,1)
              p(j,ke,1) = p(j,ke,2)
 113        continue
          endif
        endif

c..update the total energy
     
        do 51 k = ks,ke
        do 51 j = js,je
          q(j,k,l,5) = p(j,k,1)/(gm1*q(j,k,l,6)) +0.5*(q(j,k,l,2)**2
     $                 +q(j,k,l,3)**2 +q(j,k,l,4)**2)/q(j,k,l,1)
   51   continue

      endif

      return
      end

c***********************************************************************
      subroutine bcwall_internal(q,xx,xy,xz,ug,yx,yy,yz,vg,zx,zy,zz,wg,
     <     zx0,zy0,zz0,zt0,kkr,kkp,js,je,ks,ke,ls,le,idir)

c*** Prologue : ************
c
c  solid wall boundary condition for a l = 1 solid surface. 
c  calls the tangential bc subroutine bctany and computes wall pressure
c  need to still validate the higher order interpolation (ihigh=1)
c
c  last updated 07/13/04 by jaina
c***********************************************************************
      
      use params_global

c*******************  list of global variables used ******************
c
c     jmax,kmax,lmax,nd,half,invisc,ibcwp,gm1,istep0
c     rinf,uinf,vinf,winf,einf
c
c***********************************************************************

      implicit none

      real q(jmax,kmax,lmax,nd)

      real zx0(jmax,kmax),zy0(jmax,kmax),zz0(jmax,kmax),zt0(jmax,kmax)
      integer kkr(kmax),kkp(kmax)

      integer js,je,ks,ke,ls,le,idir

      real xx(jmax,kmax,lmax),xy(jmax,kmax,lmax),xz(jmax,kmax,lmax),
     <     yx(jmax,kmax,lmax),yy(jmax,kmax,lmax),yz(jmax,kmax,lmax),
     <     zx(jmax,kmax,lmax),zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax),vg(jmax,kmax,lmax),wg(jmax,kmax,lmax)

c..   local variables

      integer j,k,l
      integer iadd,iadir,j1,j2,k1,k2,l1,l2
      real foso 
      real rinv,u,v,w,p1,p2,press
      logical tmp

c**** first executable statement      

c...setting inviscid to true for wind tunnel walls
      tmp = invisc
      invisc = .true.

      iadd = sign(1,idir)
      iadir = abs(idir)

      foso = 0.

      if(iadir.eq.1) then
        j  = js
        j1 = j  + iadd
        j2 = j1 + iadd

c..compute surface velocities

        call bctany( q,xx,xy,xz,yx,yy,yz,zx,zy,zz,ug,vg,wg,
     <               js,je,ks,ke,ls,le,idir)
c
c..extrapolate pressure and density to the surface
c     
        do l = ls,le
        do k = ks,ke
          rinv = 1./q(j,k,l,1)
          u = q(j,k,l,2)*rinv
          v = q(j,k,l,3)*rinv
          w = q(j,k,l,4)*rinv
          q(j,k,l,1) = (1.+foso)*q(j1,k,l,1)*q(j1,k,l,6)/q(j,k,l,6)
     <                  - foso*q(j2,k,l,1)*q(j2,k,l,6)/q(j,k,l,6)
          q(j,k,l,1) = q(j1,k,l,1)*q(j1,k,l,6)/q(j,k,l,6)
          q(j,k,l,2) = u*q(j,k,l,1)
          q(j,k,l,3) = v*q(j,k,l,1)
          q(j,k,l,4) = w*q(j,k,l,1)
          p1 = gm1*(q(j1,k,l,5) 
     $         -0.5*(q(j1,k,l,2)**2+q(j1,k,l,3)**2+q(j1,k,l,4)**2)
     $             /q(j1,k,l,1))*q(j1,k,l,6)
          p2 = gm1*(q(j2,k,l,5)
     $         -0.5*(q(j2,k,l,2)**2+q(j2,k,l,3)**2+q(j2,k,l,4)**2)
     $             /q(j2,k,l,1))*q(j2,k,l,6)

          press = (1.+foso)*p1 - foso*p2

          q(j,k,l,5) = press/(gm1*q(j,k,l,6))+0.5*(q(j,k,l,2)**2
     $                +q(j,k,l,3)**2+q(j,k,l,4)**2)/q(j,k,l,1)
        enddo
        enddo

      elseif(iadir.eq.2) then
        k  = ks
        k1 = k  + iadd
        k2 = k1 + iadd

c..compute surface velocities

        call bctany( q,xx,xy,xz,yx,yy,yz,zx,zy,zz,ug,vg,wg,
     <               js,je,ks,ke,ls,le,idir)
c
c..extrapolate pressure and density to the surface
c     
        do j = js,je
        do l = ls,le
          rinv = 1./q(j,k,l,1)
          u = q(j,k,l,2)*rinv
          v = q(j,k,l,3)*rinv
          w = q(j,k,l,4)*rinv
          q(j,k,l,1) = (1.+foso)*q(j,k1,l,1)*q(j,k1,l,6)/q(j,k,l,6)
     <                  - foso*q(j,k2,l,1)*q(j,k2,l,6)/q(j,k,l,6)
          q(j,k,l,1) = q(j,k1,l,1)*q(j,k1,l,6)/q(j,k,l,6)
          q(j,k,l,2) = u*q(j,k,l,1)
          q(j,k,l,3) = v*q(j,k,l,1)
          q(j,k,l,4) = w*q(j,k,l,1)
          p1 = gm1*(q(j,k1,l,5) 
     $         -0.5*(q(j,k1,l,2)**2+q(j,k1,l,3)**2+q(j,k1,l,4)**2)
     $             /q(j,k1,l,1))*q(j,k1,l,6)
          p2 = gm1*(q(j,k2,l,5)
     $         -0.5*(q(j,k2,l,2)**2+q(j,k2,l,3)**2+q(j,k2,l,4)**2)
     $             /q(j,k2,l,1))*q(j,k2,l,6)

          press = (1.+foso)*p1 - foso*p2

          q(j,k,l,5) = press/(gm1*q(j,k,l,6))+0.5*(q(j,k,l,2)**2
     $                +q(j,k,l,3)**2+q(j,k,l,4)**2)/q(j,k,l,1)
        enddo
        enddo

      elseif(iadir.eq.3) then
        l  = ls
        l1 = l  + iadd
        l2 = l1 + iadd

c..compute surface velocities

        call bctany( q,xx,xy,xz,yx,yy,yz,zx,zy,zz,ug,vg,wg,
     <               js,je,ks,ke,ls,le,idir)
c
c..extrapolate pressure and density to the surface
c     
        do k = ks,ke
        do j = js,je
          rinv = 1./q(j,k,l,1)
          u = q(j,k,l,2)*rinv
          v = q(j,k,l,3)*rinv
          w = q(j,k,l,4)*rinv
          q(j,k,l,1) = (1.+foso)*q(j,k,l1,1)*q(j,k,l1,6)/q(j,k,l,6)
     <                  - foso*q(j,k,l2,1)*q(j,k,l2,6)/q(j,k,l,6)
          q(j,k,l,1) = q(j,k,l1,1)*q(j,k,l1,6)/q(j,k,l,6)
          q(j,k,l,2) = u*q(j,k,l,1)
          q(j,k,l,3) = v*q(j,k,l,1)
          q(j,k,l,4) = w*q(j,k,l,1)
          p1 = gm1*(q(j,k,l1,5) 
     $         -0.5*(q(j,k,l1,2)**2+q(j,k,l1,3)**2+q(j,k,l1,4)**2)
     $             /q(j,k,l1,1))*q(j,k,l1,6)
          p2 = gm1*(q(j,k,l2,5)
     $         -0.5*(q(j,k,l2,2)**2+q(j,k,l2,3)**2+q(j,k,l2,4)**2)
     $             /q(j,k,l2,1))*q(j,k,l2,6)

          press = (1.+foso)*p1 - foso*p2

          q(j,k,l,5) = press/(gm1*q(j,k,l,6))+0.5*(q(j,k,l,2)**2
     $                +q(j,k,l,3)**2+q(j,k,l,4)**2)/q(j,k,l,1)
        enddo
        enddo
 
      endif

c..resetting the inviscid value to original
      invisc = tmp

      return
      end

c***********************************************************************
      subroutine bcwall(q,xx,xy,xz,ug,yx,yy,yz,vg,zx,zy,zz,wg,
     <     zx0,zy0,zz0,zt0,kkr,kkp,js,je,ks,ke,ls,le,idir)

c*** Prologue : ************
c
c  solid wall boundary condition for a l = 1 solid surface. 
c  calls the tangential bc subroutine bctany and computes wall pressure
c  need to still validate the higher order interpolation (ihigh=1)
c
c  last updated 07/13/04 by jaina
c***********************************************************************
      
      use params_global

c*******************  list of global variables used ******************
c
c     jmax,kmax,lmax,nd,half,invisc,ibcwp,gm1,istep0
c     rinf,uinf,vinf,winf,einf
c
c***********************************************************************

      implicit none

      real q(jmax,kmax,lmax,nd)

      real zx0(jmax,kmax),zy0(jmax,kmax),zz0(jmax,kmax),zt0(jmax,kmax)
      integer kkr(kmax),kkp(kmax)

      integer js,je,ks,ke,ls,le,idir

      real xx(jmax,kmax,lmax),xy(jmax,kmax,lmax),xz(jmax,kmax,lmax),
     <     yx(jmax,kmax,lmax),yy(jmax,kmax,lmax),yz(jmax,kmax,lmax),
     <     zx(jmax,kmax,lmax),zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax),vg(jmax,kmax,lmax),wg(jmax,kmax,lmax)

c..   local variables

      real,allocatable :: p(:,:,:),a(:),b(:),c(:)
      integer j,k,l
      integer iadd,iadir,l1,l2,l3,ihigh
      real use,vse,wse
      real rinv,u,v,w,rj,us,vs,ws,ue,ve,we,t
      real scal,rscal,ajacinv

      allocate(p(jmax,kmax,3),a(jmax*kmax),b(jmax*kmax),c(jmax*kmax))

c**** first executable statement      
      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then
         print*,'idir = ',idir,' is not implemented in bcwall'
      elseif(iadir.eq.2) then
         print*,'idir = ',idir,' is not implemented in bcwall'
      elseif(iadir.eq.3) then
        l  = ls
        l1 = l  + iadd
        l2 = l1 + iadd
        l3 = l2 + iadd

        ihigh = 0

c       if (invisc) ihigh = 1

c..compute surface velocities

        call bctany( q,xx,xy,xz,yx,yy,yz,zx,zy,zz,ug,vg,wg,
     <               js,je,ks,ke,ls,le,idir)

c..at the trailing edge (don't do if half-plane!)

        if(half.ne.1 .and. invisc) then

         do 211 k = ks,ke

c...  higher order?

         if(ihigh.eq.1) then
          use = (-q(je,k,l2,2)/q(je,k,l2,1)+4.*q(je,k,l1,2)/q(je,k,l1,1)
     <        +4.*q(js,k,l1,2)/q(js,k,l1,1)-q(js,k,l2,2)/q(js,k,l2,1))/6.
          vse = (-q(je,k,l2,3)/q(je,k,l2,1)+4.*q(je,k,l1,3)/q(je,k,l1,1)
     <        +4.*q(js,k,l1,3)/q(js,k,l1,1)-q(js,k,l2,3)/q(js,k,l2,1))/6.
          wse = (-q(je,k,l2,4)/q(je,k,l2,1)+4.*q(je,k,l1,4)/q(je,k,l1,1)
     <        +4.*q(js,k,l1,4)/q(js,k,l1,1)-q(js,k,l2,4)/q(js,k,l2,1))/6.
         else
            use = 0.5*(q(je,k,l1,2)/q(je,k,l1,1)+q(js,k,l1,2)
     <           /q(js,k,l1,1))
            vse = 0.5*(q(je,k,l1,3)/q(je,k,l1,1)+q(js,k,l1,3)
     <           /q(js,k,l1,1))
            wse = 0.5*(q(je,k,l1,4)/q(je,k,l1,1)+q(js,k,l1,4)
     <           /q(js,k,l1,1))
         endif

         q(js,k,l,2) = use*q(js,k,l,1)
         q(je,k,l,2) = use*q(je,k,l,1)
         q(js,k,l,3) = vse*q(js,k,l,1)
         q(je,k,l,3) = vse*q(je,k,l,1)
         q(js,k,l,4) = wse*q(js,k,l,1)
         q(je,k,l,4) = wse*q(je,k,l,1)
  211    continue
        endif

        if(half.eq.1) then
          do 221 k = ks,ke
            q(js,k,l,4) = 0.0
  221     continue
        endif
c
c..compute total energy and static wall pressure
c
        if(ibcwp.eq.0) then
          call bcwp0(p,q,kkr,kkp,js,je,ks,ke,ls,le,idir)
        else
          call bcwp(a,b,c,p,q,xx,xy,xz,yx,yy,yz,zx,zy,zz,
     $       zx0,zy0,zz0,zt0,kkr,kkp,ug,vg,wg,js,je,ks,ke,ls,le,idir)

        endif
c
c..extrapolate density to the surface
c     
        if( invisc ) then
         do 11 k = ks,ke
         do 11 j = js,je
           rinv = 1./q(j,k,l,1)
           u = q(j,k,l,2)*rinv
           v = q(j,k,l,3)*rinv
           w = q(j,k,l,4)*rinv
           q(j,k,l,1) = 2.*q(j,k,l1,1)*q(j,k,l1,6)/q(j,k,l,6)
     <                   - q(j,k,l2,1)*q(j,k,l2,6)/q(j,k,l,6)
           q(j,k,l,1) = q(j,k,l1,1)*q(j,k,l1,6)/q(j,k,l,6)
c...higher order?
           if(ihigh.eq.1)
     <     q(j,k,l,1) = 3.*q(j,k,l1,1)*q(j,k,l1,6)/q(j,k,l,6)
     <                - 3.*q(j,k,l2,1)*q(j,k,l2,6)/q(j,k,l,6)
     <                    +q(j,k,l3,1)*q(j,k,l3,6)/q(j,k,l,6)
           q(j,k,l,2) = u*q(j,k,l,1)
           q(j,k,l,3) = v*q(j,k,l,1)
           q(j,k,l,4) = w*q(j,k,l,1)
           q(j,k,l,5) = p(j,k,1)/(gm1*q(j,k,l,6))+0.5*(q(j,k,l,2)**2
     $                 +q(j,k,l,3)**2+q(j,k,l,4)**2)/q(j,k,l,1)
 11     continue
        else
         do 12 k = ks,ke
         do 12 j = js,je
           rinv = 1./q(j,k,l,1)
           u = q(j,k,l,2)*rinv
           v = q(j,k,l,3)*rinv
           w = q(j,k,l,4)*rinv
           q(j,k,l,1) = q(j,k,l1,1)*q(j,k,l1,6)/q(j,k,l,6)
           q(j,k,l,2) = u*q(j,k,l,1)
           q(j,k,l,3) = v*q(j,k,l,1)
           q(j,k,l,4) = w*q(j,k,l,1)
           q(j,k,l,5) = p(j,k,1)/(gm1*q(j,k,l,6))+0.5*(q(j,k,l,2)**2
     $                 +q(j,k,l,3)**2+q(j,k,l,4)**2)/q(j,k,l,1)
 12     continue
        endif
c
c..at the trailing edge (don't do if half-plane!)
c
        if(half.ne.1) then
          do 111 k = ks,ke
          rj = q(js,k,l,6)/q(je,k,l,6)
          rinv = 1./q(js,k,l,1)
          us = q(js,k,l,2)*rinv
          vs = q(js,k,l,3)*rinv
          ws = q(js,k,l,4)*rinv
          rinv = 1./q(je,k,l,1)
          ue = q(je,k,l,2)*rinv
          ve = q(je,k,l,3)*rinv
          we = q(je,k,l,4)*rinv
          q(js,k,l,1) = 0.5*( q(js,k,l,1)+q(je,k,l,1)/rj )
          q(js,k,l,2) = us*q(js,k,l,1)
          q(js,k,l,3) = vs*q(js,k,l,1)
          q(js,k,l,4) = ws*q(js,k,l,1)
          q(js,k,l,5) = p(js,k,1)/(gm1*q(js,k,l,6))+0.5*(q(js,k,l,2)**2
     $                 +q(js,k,l,3)**2+q(js,k,l,4)**2)/q(js,k,l,1)
          q(je,k,l,1) = q(js,k,l,1)*rj
          q(je,k,l,2) = ue*q(je,k,l,1)
          q(je,k,l,3) = ve*q(je,k,l,1)
          q(je,k,l,4) = we*q(je,k,l,1)
          q(je,k,l,5) = p(je,k,1)/(gm1*q(je,k,l,6))+0.5*(q(je,k,l,2)**2
     $                 +q(je,k,l,3)**2+q(je,k,l,4)**2)/q(je,k,l,1)
  111     continue
        endif
c
c..slowly turn on the wall boundary condition
c  over 20 time steps
c
        t = (istep0-1.)/30.
        if(t.gt.1..or.iread.gt.0) t=1.
        !t = 1
        scal = (10.-15.*t+6.*t*t)*t**3
        rscal = 1.-scal
c
        do k=ks,ke
        do j=js,je
          ajacinv = 1./q(j,k,l,6)
          q(j,k,l,1) = q(j,k,l,1)*scal+rscal*rinf*ajacinv
          q(j,k,l,2) = q(j,k,l,2)*scal+rscal*rinf*uinf*ajacinv
          q(j,k,l,3) = q(j,k,l,3)*scal+rscal*rinf*vinf*ajacinv
          q(j,k,l,4) = q(j,k,l,4)*scal+rscal*rinf*winf*ajacinv
          q(j,k,l,5) = q(j,k,l,5)*scal+rscal*einf*ajacinv
        enddo
        enddo

      endif
c
      return
      end


c***********************************************************************
        subroutine bcwake(q,x,y,z,js,je,ks,ke,ls,le,idir)
c     
c*** Prologue ********
c
c  for complete c-o grid   
c  treatment for plane k =kroot=1
c  dkarthik (2003)
c
c  last updated by jaina 07/13/04
c*********************************************************************

      use params_global

c*******************  list of global variables used ******************
c
c     jmax,kmax,lmax,nd,ktip,kroot,jtail1
c
c*********************************************************************

      implicit none
      
      real q(jmax,kmax,lmax,nd)
      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)
      integer js,je,ks,ke,ls,le,idir
      
c..   local variables
      
      integer j,k,l,jj,k1,l1,l2,iadd,iadir,ihigh
      real sc1,sc2,scale1,scale2,sc,qav1,qav2,qav3,qav4,qav5
      real h1,h2,f1,f2

c..   in the wake, average points above and below
       
      iadd = sign(1,idir)
      iadir = abs(idir)
      if(iadir.eq.1) then
         print*,'idir = ',idir,' is not implemented in bcwake'
        
      elseif(iadir.eq.2) then
         k=ks
         k1=k+iadd

         do l=ls,le
            do j=js,je
               
               jj=jmax-j+1

               sc1 = q(j,k1,l,6)/q(j,k,l,6)
               scale1= q(jj,k1,l,6)/q(j,k,l,6)
               sc= q(j,k,l,6)/q(jj,k,l,6)
               qav1 = 0.5*(q(j,k1,l,1)*sc1+q(jj,k1,l,1)*scale1)
               qav2 = 0.5*(q(j,k1,l,2)*sc1+q(jj,k1,l,2)*scale1)
               qav3 = 0.5*(q(j,k1,l,3)*sc1+q(jj,k1,l,3)*scale1)
               qav4 = 0.5*(q(j,k1,l,4)*sc1+q(jj,k1,l,4)*scale1)
               qav5 = 0.5*(q(j,k1,l,5)*sc1+q(jj,k1,l,5)*scale1)
               q(j,k,l,1)= qav1
               q(j,k,l,2)= qav2
               q(j,k,l,3)= qav3
               q(j,k,l,4)= qav4
               q(j,k,l,5)= qav5
               q(jj,k,l,1)= qav1*sc
               q(jj,k,l,2)= qav2*sc
               q(jj,k,l,3)= qav3*sc
               q(jj,k,l,4)= qav4*sc
               q(jj,k,l,5)= qav5*sc
               
            enddo
            
c..ave   rage at leading edge
c..now    that q(j,k,l)=q(j,k1,l)*sc

            j = je
            jj = j-1
            scale1= q(jj,k,l,6)/q(j,k,l,6)
            h1 = sqrt( (x(j,k1,l)-x(j,k,l))**2 +
     <           (y(j,k1,l)-y(j,k,l))**2 +
     <           (z(j,k1,l)-z(j,k,l))**2 )
            h2 = sqrt( (x(j,k,l)-x(jj,k,l))**2 +
     <           (y(j,k,l)-y(jj,k,l))**2 +
     <           (z(j,k,l)-z(jj,k,l))**2 )
            f1 = h2/(h1+h2)
            f2 = 1.-f1

            q(j,k,l,1) = f1*q(j,k,l,1)+f2*q(jj,k,l,1)*scale1
            q(j,k,l,2) = f1*q(j,k,l,2)+f2*q(jj,k,l,2)*scale1
            q(j,k,l,3) = f1*q(j,k,l,3)+f2*q(jj,k,l,3)*scale1
            q(j,k,l,4) = f1*q(j,k,l,4)+f2*q(jj,k,l,4)*scale1
            q(j,k,l,5) = f1*q(j,k,l,5)+f2*q(jj,k,l,5)*scale1
            
         enddo

      elseif(iadir.eq.3) then
         ihigh = 0
         l=ls
         l1=l+iadd
         l2=l1+iadd
         do k = ks,ke
            do j = js,je
               jj = jmax - j + 1
               sc1 = q(j,k,l1,6)/q(j,k,l,6)
               sc2 = q(j,k,l2,6)/q(j,k,l,6)
               scale1= q(jj,k,l1,6)/q(j,k,l,6)
               scale2= q(jj,k,l2,6)/q(j,k,l,6)
               sc= q(j,k,l,6)/q(jj,k,l,6)
               qav1 = 0.5*(q(j,k,l1,1)*sc1+q(jj,k,l1,1)*scale1)
               qav2 = 0.5*(q(j,k,l1,2)*sc1+q(jj,k,l1,2)*scale1)
               qav3 = 0.5*(q(j,k,l1,3)*sc1+q(jj,k,l1,3)*scale1)
               qav4 = 0.5*(q(j,k,l1,4)*sc1+q(jj,k,l1,4)*scale1)
               qav5 = 0.5*(q(j,k,l1,5)*sc1+q(jj,k,l1,5)*scale1)
               if (ihigh.eq.1) then
                 qav1 = (-q(j,k,l2,1)*sc2 + 4*q(j,k,l1,1)*sc1 +
     c                  4*q(jj,k,l1,1)*scale1 - q(jj,k,l2,1)*scale2)/6.
                 qav2 = (-q(j,k,l2,2)*sc2 + 4*q(j,k,l1,2)*sc1 +
     c                  4*q(jj,k,l1,2)*scale1 - q(jj,k,l2,2)*scale2)/6.
                 qav3 = (-q(j,k,l2,3)*sc2 + 4*q(j,k,l1,3)*sc1 +
     c                  4*q(jj,k,l1,3)*scale1 - q(jj,k,l2,3)*scale2)/6.
                 qav4 = (-q(j,k,l2,4)*sc2 + 4*q(j,k,l1,4)*sc1 +
     c                  4*q(jj,k,l1,4)*scale1 - q(jj,k,l2,4)*scale2)/6.
                 qav5 = (-q(j,k,l2,5)*sc2 + 4*q(j,k,l1,5)*sc1 +
     c                  4*q(jj,k,l1,5)*scale1 - q(jj,k,l2,5)*scale2)/6.
               endif
               q(j,k,l,1)= qav1
               q(j,k,l,2)= qav2
               q(j,k,l,3)= qav3
               q(j,k,l,4)= qav4
               q(j,k,l,5)= qav5
               q(jj,k,l,1)= qav1*sc
               q(jj,k,l,2)= qav2*sc
               q(jj,k,l,3)= qav3*sc
               q(jj,k,l,4)= qav4*sc
               q(jj,k,l,5)= qav5*sc
            enddo
         enddo
      endif

      return
      end

c**********************************************************************
      subroutine bcout(q,x,y,z,xx,xy,xz,yx,yy,yz,zx,zy,zz,
     >     ug,vg,wg,js,je,ks,ke,ls,le,idir)
c*** Prologue : ***********
c     
c  Characteristically extrapolate the flow variables at the outer
c  boundaries using the Riemann invariants
c
c  subroutine was written by dkarthik (2003)
c  last update 07/13/04 by jaina with an if loop to support c-o mesh
c  also updates for l=lmax and k=kmax boundaries for c-h mesh
c**********************************************************************

      use params_global

c*******************  list of global variables used ********************
c
c     jmax,kmax,lmax,nd,gm1,gamma,uinf,vinf,winf,rinf,pinf,ktip
c***********************************************************************

      implicit none

      real q(jmax,kmax,lmax,nd)
      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)
      real xx(jmax,kmax,lmax),xy(jmax,kmax,lmax),xz(jmax,kmax,lmax),
     <     yx(jmax,kmax,lmax),yy(jmax,kmax,lmax),yz(jmax,kmax,lmax),
     <     zx(jmax,kmax,lmax),zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax),vg(jmax,kmax,lmax),wg(jmax,kmax,lmax)
      integer js,je,ks,ke,ls,le,idir

c..local variables

      integer j,k,l,iadd,iadir,j1,j2,k1,k2,l1,l2
      real gm1i,gi,rl,rlm1,rhoext,uext,vext,wext,eext,pext
      real uind,vind,wind,mind2,rind,pind,snorm
      real zxn,zyn,zzn
      real xxn,xyn,xzn
      real yxn,yyn,yzn
      real uinn,r1,uexn,r2,qn,cspe,c2
      real vel2,vel2ext,velcheck
      real stj,entro,u,v,w,press

c**** first executable statement

      gm1i = 1./gm1
      gi   = 1./gamma

      iadd = sign(1,idir)
      iadir = abs(idir)
      
      if (iadir.eq.1) then
      j  = js
      j1 = j + iadd

c note that rl, rlm1 etc.. are nothing but rj,rjp1 etc..
c for ease of copy and paste
c dkarthik

      do  l = ls,le
      do  k = ks,ke
        rl    = 1./q(j,k,l,6)
        rlm1  = 1./q(j1,k,l,1)
c..zeroeth-order extrapolation
        rhoext= q(j1,k,l,1)*q(j1,k,l,6)
        uext  = q(j1,k,l,2)*rlm1
        vext  = q(j1,k,l,3)*rlm1
        wext  = q(j1,k,l,4)*rlm1
        eext  = q(j1,k,l,5)*q(j1,k,l,6)
        pext  = gm1*(eext - 0.5*rhoext*(uext**2+vext**2+wext**2))
c..free stream
	uind = uinf 
        vind = vinf 
        wind = winf 
        mind2 = uind**2+vind**2+wind**2
c        rind = rinf*(1.-.5*rinf/pinf*gm1*gi*(mind2-fsmach**2))**gm1i
c        pind = pinf*(rind/rinf)**gamma
        rind = rinf
        pind = pinf

        snorm = -1.*iadd/sqrt(xx(j,k,l)**2+xy(j,k,l)**2+xz(j,k,l)**2)
        xxn = xx(j,k,l)*snorm
        xyn = xy(j,k,l)*snorm
        xzn = xz(j,k,l)*snorm
c..calculate riemann invariants
        uinn = (uind-ug(j,k,l))*xxn + (vind-vg(j,k,l))*xyn
     <       + (wind-wg(j,k,l))*xzn
        r1 = uinn -2.*sqrt(gamma*pind/rind)*gm1i
        uexn = (uext-ug(j,k,l))*xxn + (vext-vg(j,k,l))*xyn
     <       + (wext-wg(j,k,l))*xzn
        r2 = uexn +2.*sqrt(gamma*pext/rhoext)*gm1i
c..calculate normal velocity and speed of sound based on riemann
        qn = 0.5*(r1+r2)
        cspe = (r2-r1)*gm1*0.25
        c2 = cspe**2
c..is flow relatively subsonic or supersonic?
        velcheck = abs(qn)
c..calculate contributions from interior and exterior
        if (qn .lt. 0) then
c..inflow boundary
          if(velcheck .lt. 1.0) then
c..fix four and extrapolate one (pressure)
            stj = qn - uinn
            entro = rind**gamma/pind
            u = uind + stj*xxn
            v = vind + stj*xyn
            w = wind + stj*xzn
            q(j,k,l,1) = (c2*entro*gi)**gm1i
            press = c2*q(j,k,l,1)*gi
            q(j,k,l,1) = q(j,k,l,1)*rl
            q(j,k,l,2) = q(j,k,l,1)*u
            q(j,k,l,3) = q(j,k,l,1)*v
            q(j,k,l,4) = q(j,k,l,1)*w
            q(j,k,l,5) = press/gm1*rl + 0.5*(q(j,k,l,2)**2 +
     <                   q(j,k,l,3)**2+q(j,k,l,4)**2)/q(j,k,l,1)
          else
c..fix five
            q(j,k,l,1) = rind*rl
            q(j,k,l,2) = rind*uind*rl
            q(j,k,l,3) = rind*vind*rl
            q(j,k,l,4) = rind*wind*rl
            press = pind
            q(j,k,l,5) = press/gm1*rl + 0.5*(q(j,k,l,2)**2 +
     <                   q(j,k,l,3)**2+q(j,k,l,4)**2)/q(j,k,l,1)
          endif
        else
c..outflow boundary
          if(velcheck .lt. 1.0) then
c..prescribe p and extrapolate rho,u,v,and w
            stj = qn - uexn
            entro = rhoext**gamma/pext
            u = uext + stj*xxn
            v = vext + stj*xyn
            w = wext + stj*xzn
            q(j,k,l,1) = (c2*entro*gi)**gm1i
            press = c2*q(j,k,l,1)*gi
            q(j,k,l,1) = q(j,k,l,1)*rl
            q(j,k,l,2) = q(j,k,l,1)*u
            q(j,k,l,3) = q(j,k,l,1)*v
            q(j,k,l,4) = q(j,k,l,1)*w
            q(j,k,l,5) = press/gm1*rl + 0.5*(q(j,k,l,2)**2 +
     <                   q(j,k,l,3)**2+q(j,k,l,4)**2)/q(j,k,l,1)
          else
c..extrapolate five
            q(j,k,l,1) = rhoext*rl
            q(j,k,l,2) = rhoext*uext*rl
            q(j,k,l,3) = rhoext*vext*rl
            q(j,k,l,4) = rhoext*wext*rl
            press = pext
            q(j,k,l,5) = press/gm1*rl + 0.5*(q(j,k,l,2)**2 +
     <                   q(j,k,l,3)**2+q(j,k,l,4)**2)/q(j,k,l,1)
          endif
        endif
c
	enddo
	enddo

      elseif (iadir.eq.2) then
       k  = ks
       k1 = k + iadd
c     
       do  l = ls,le
       do  j = js,je
          rl    = 1./q(j,k,l,6)
          rlm1  = 1./q(j,k1,l,1)
               
c..   zeroeth-order extrapolation

          rhoext= q(j,k1,l,1)*q(j,k1,l,6)
          uext  = q(j,k1,l,2)*rlm1
          vext  = q(j,k1,l,3)*rlm1
          wext  = q(j,k1,l,4)*rlm1
          eext  = q(j,k1,l,5)*q(j,k1,l,6)
          pext  = gm1*(eext - 0.5*rhoext*(uext**2+vext**2+wext**2))
               
c..   free stream
               
          uind = uinf 
          vind = vinf 
          wind = winf 
          mind2 = uind**2+vind**2+wind**2
          rind = rinf
          pind = pinf
             
          snorm = -1.*iadd/sqrt(yx(j,k,l)**2+yy(j,k,l)**2+yz(j,k,l)**2)
          
          yxn = yx(j,k,l)*snorm
          yyn = yy(j,k,l)*snorm
          yzn = yz(j,k,l)*snorm
        
c..calculate riemann invariants

          uinn = (uind-ug(j,k,l))*yxn + (vind-vg(j,k,l))*yyn
     <         + (wind-wg(j,k,l))*yzn
          r1 = uinn -2.*sqrt(gamma*pind/rind)*gm1i
          uexn = (uext-ug(j,k,l))*yxn + (vext-vg(j,k,l))*yyn
     <         + (wext-wg(j,k,l))*yzn
          r2 = uexn +2.*sqrt(gamma*pext/rhoext)*gm1i

c..calculate normal velocity and speed of sound based on riemann

          qn = 0.5*(r1+r2)
          cspe = (r2-r1)*gm1*0.25
          c2 = cspe**2
            
c..is flow relatively subsonic or supersonic?
               
          velcheck = abs(qn)
            
c..calculate contributions from interior and exterior

          if (qn .lt. 0) then
c..inflow boundary
             if(velcheck .lt. 1.0) then
c..   fix four and extrapolate one (pressure)
                stj = qn - uinn
                entro = rind**gamma/pind
                u = uind + stj*yxn
                v = vind + stj*yyn
                w = wind + stj*yzn
                q(j,k,l,1) = (c2*entro*gi)**gm1i
                press = c2*q(j,k,l,1)*gi
                q(j,k,l,1) = q(j,k,l,1)*rl
                q(j,k,l,2) = q(j,k,l,1)*u
                q(j,k,l,3) = q(j,k,l,1)*v
                q(j,k,l,4) = q(j,k,l,1)*w
                q(j,k,l,5) = press/gm1*rl + 0.5*(q(j,k,l,2)**2 +
     <               q(j,k,l,3)**2+q(j,k,l,4)**2)/q(j,k,l,1)
             else
c..fix five
                q(j,k,l,1) = rind*rl
                q(j,k,l,2) = rind*uind*rl
                q(j,k,l,3) = rind*vind*rl
                q(j,k,l,4) = rind*wind*rl
                press = pind
                q(j,k,l,5) = press/gm1*rl + 0.5*(q(j,k,l,2)**2 +
     <               q(j,k,l,3)**2+q(j,k,l,4)**2)/q(j,k,l,1)
             endif
          else
                  
c..outflow boundary
        
             if(velcheck .lt. 1.0) then
                     
c..prescribe p and extrapolate rho,u,v,and w
                     
                stj = qn - uexn
                entro = rhoext**gamma/pext
                u = uext + stj*yxn
                v = vext + stj*yyn
                w = wext + stj*yzn
                q(j,k,l,1) = (c2*entro*gi)**gm1i
                press = c2*q(j,k,l,1)*gi
                q(j,k,l,1) = q(j,k,l,1)*rl
                q(j,k,l,2) = q(j,k,l,1)*u
                q(j,k,l,3) = q(j,k,l,1)*v
                q(j,k,l,4) = q(j,k,l,1)*w
                q(j,k,l,5) = press/gm1*rl + 0.5*(q(j,k,l,2)**2 +
     <               q(j,k,l,3)**2+q(j,k,l,4)**2)/q(j,k,l,1)
             else
            
c..extrapolate five

                q(j,k,l,1) = rhoext*rl
                q(j,k,l,2) = rhoext*uext*rl
                q(j,k,l,3) = rhoext*vext*rl
                q(j,k,l,4) = rhoext*wext*rl
                press = pext
                q(j,k,l,5) = press/gm1*rl + 0.5*(q(j,k,l,2)**2 +
     <               q(j,k,l,3)**2+q(j,k,l,4)**2)/q(j,k,l,1)
             endif
          endif
               
       enddo
       enddo

      elseif (iadir.eq.3) then
      l  = ls
      l1 = l + iadd

c..remove l=lmax condition to check

      do j = js,je
      do k = ks,ke
        rl    = 1./q(j,k,l,6)
        rlm1  = 1./q(j,k,l1,1)
c..zeroeth-order extrapolation
        rhoext= q(j,k,l1,1)*q(j,k,l1,6)
        uext  = q(j,k,l1,2)*rlm1
        vext  = q(j,k,l1,3)*rlm1
        wext  = q(j,k,l1,4)*rlm1
        eext  = q(j,k,l1,5)*q(j,k,l1,6)
        pext  = gm1*(eext - 0.5*rhoext*(uext**2+vext**2+wext**2))
c.free stream
	uind = uinf 
        vind = vinf 
        wind = winf 
        mind2 = uind**2+vind**2+wind**2
        rind = rinf
        pind = pinf

        snorm = -1.*iadd/sqrt(zx(j,k,l)**2+zy(j,k,l)**2+zz(j,k,l)**2)
        zxn = zx(j,k,l)*snorm
        zyn = zy(j,k,l)*snorm
        zzn = zz(j,k,l)*snorm
c..calculate riemann invariants
        uinn = (uind-ug(j,k,l))*zxn + (vind-vg(j,k,l))*zyn
     <       + (wind-wg(j,k,l))*zzn
        r1 = uinn -2.*sqrt(gamma*pind/rind)*gm1i
        uexn = (uext-ug(j,k,l))*zxn + (vext-vg(j,k,l))*zyn
     <       + (wext-wg(j,k,l))*zzn
        r2 = uexn +2.*sqrt(gamma*pext/rhoext)*gm1i
c..calculate normal velocity and speed of sound based on riemann
        qn = 0.5*(r1+r2)
        cspe = (r2-r1)*gm1*0.25
        c2 = cspe**2
c..is flow relatively subsonic or supersonic?
        velcheck = abs(qn)
c..calculate contributions from interior and exterior
        if (qn .lt. 0) then
c..inflow boundary
          if(velcheck .lt. 1.0) then
c..fix four and extrapolate one (pressure)
            stj = qn - uinn
            entro = rind**gamma/pind
            u = uind + stj*zxn
            v = vind + stj*zyn
            w = wind + stj*zzn
            q(j,k,l,1) = (c2*entro*gi)**gm1i
            press = c2*q(j,k,l,1)*gi
            q(j,k,l,1) = q(j,k,l,1)*rl
            q(j,k,l,2) = q(j,k,l,1)*u
            q(j,k,l,3) = q(j,k,l,1)*v
            q(j,k,l,4) = q(j,k,l,1)*w
            q(j,k,l,5) = press/gm1*rl + 0.5*(q(j,k,l,2)**2 +
     <                   q(j,k,l,3)**2+q(j,k,l,4)**2)/q(j,k,l,1)
          else
c..fix five
            q(j,k,l,1) = rind*rl
            q(j,k,l,2) = rind*uind*rl
            q(j,k,l,3) = rind*vind*rl
            q(j,k,l,4) = rind*wind*rl
            press = pind
            q(j,k,l,5) = press/gm1*rl + 0.5*(q(j,k,l,2)**2 +
     <                   q(j,k,l,3)**2+q(j,k,l,4)**2)/q(j,k,l,1)
          endif
        else
c..outflow boundary
          if(velcheck .lt. 1.0) then
c..prescribe p and extrapolate rho,u,v,and w
            stj = qn - uexn
            entro = rhoext**gamma/pext
            u = uext + stj*zxn
            v = vext + stj*zyn
            w = wext + stj*zzn
            q(j,k,l,1) = (c2*entro*gi)**gm1i
            press = c2*q(j,k,l,1)*gi
            q(j,k,l,1) = q(j,k,l,1)*rl
            q(j,k,l,2) = q(j,k,l,1)*u
            q(j,k,l,3) = q(j,k,l,1)*v
            q(j,k,l,4) = q(j,k,l,1)*w
            q(j,k,l,5) = press/gm1*rl + 0.5*(q(j,k,l,2)**2 +
     <                   q(j,k,l,3)**2+q(j,k,l,4)**2)/q(j,k,l,1)
          else
c..extrapolate five
            q(j,k,l,1) = rhoext*rl
            q(j,k,l,2) = rhoext*uext*rl
            q(j,k,l,3) = rhoext*vext*rl
            q(j,k,l,4) = rhoext*wext*rl
            press = pext
            q(j,k,l,5) = press/gm1*rl + 0.5*(q(j,k,l,2)**2 +
     <                   q(j,k,l,3)**2+q(j,k,l,4)**2)/q(j,k,l,1)
          endif
        endif
c
      enddo
      enddo

      endif

      return
      end

c**********************************************************************
      subroutine bcouttrkl(q,x,y,z,xx,xy,xz,yx,yy,yz,zx,zy,zz,
     >     ug,vg,wg,bt,js,je,ks,ke,ls,le,idir)
c*** Prologue : ***********
c     
c  Characteristically extrapolate the flow variables at the outer
c  boundaries using the Riemann invariants
c
c**********************************************************************

      use params_global

c*******************  list of global variables used ********************
c
c     jmax,kmax,lmax,nd,gm1,gamma,uinf,vinf,winf,rinf,pinf,ktip
c***********************************************************************

      implicit none

      real q(jmax,kmax,lmax,nd)
      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)
      real xx(jmax,kmax,lmax),xy(jmax,kmax,lmax),xz(jmax,kmax,lmax),
     <     yx(jmax,kmax,lmax),yy(jmax,kmax,lmax),yz(jmax,kmax,lmax),
     <     zx(jmax,kmax,lmax),zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax),vg(jmax,kmax,lmax),wg(jmax,kmax,lmax)
      real bt(jmax,kmax,lmax)
      integer js,je,ks,ke,ls,le,idir

c..local variables
      real,allocatable :: qext(:),qbar(:),delq(:),delw(:)
      integer j,k,l,j1,j2,k1,k2,l1,l2,iadd,iadir
      real rhoext,uext,vext,wext,eext,pext
      real uind,vind,wind,snorm
      real zxn,zyn,zzn
      real xxn,xyn,xzn
      real yxn,yyn,yzn
      real rho,rrho,uu,vv,ww,e,uvw,cjkl,c2i,ge,qq
      real Z1,R,S,Zi,bSq

      allocate(qext(5),qbar(5),delq(5),delw(5))

c**** first executable statement

      iadd = sign(1,idir)
      iadir = abs(idir)
      

      if (iadir.eq.1) then

      j  = js
      j1 = j + iadd

      do  l = ls,le
      do  k = ks,ke
c..free stream
        uind = uinf 
        vind = vinf 
        wind = winf 

        snorm = -1.*iadd/sqrt(xx(j,k,l)**2+xy(j,k,l)**2+xz(j,k,l)**2)
        xxn = xx(j,k,l)*snorm
        xyn = xy(j,k,l)*snorm
        xzn = xz(j,k,l)*snorm

c..zeroeth-order extrapolation
        qext(1) = q(j1,k,l,1)*q(j1,k,l,6)
        qext(2) = q(j1,k,l,2)*q(j1,k,l,6)
        qext(3) = q(j1,k,l,3)*q(j1,k,l,6)
        qext(4) = q(j1,k,l,4)*q(j1,k,l,6)
        qext(5) = q(j1,k,l,5)*q(j1,k,l,6)

        qbar(1) = 0.5*(qext(1)+rinf)
        qbar(2) = 0.5*(qext(2)+rinf*uind)
        qbar(3) = 0.5*(qext(3)+rinf*vind)
        qbar(4) = 0.5*(qext(4)+rinf*wind)
        qbar(5) = 0.5*(qext(5)+einf)

        delq(1) = 0.5*(rinf-qext(1))
        delq(2) = 0.5*(rinf*uind-qext(2))
        delq(3) = 0.5*(rinf*vind-qext(3))
        delq(4) = 0.5*(rinf*wind-qext(4))
        delq(5) = 0.5*(einf-qext(5))

        rho = qbar(1)
        rrho = 1./rho
        uu = qbar(2)*rrho
        vv = qbar(3)*rrho
        ww = qbar(4)*rrho
        e  = qbar(5)*rrho
        uvw = 0.5*(uu*uu+vv*vv+ww*ww)
        cjkl = sqrt(ggm1*(e-uvw))
        c2i = 1./(cjkl*cjkl)
        ge = gamma*e - gm1*uvw
        qq = (uu-ug(j,k,l))*xxn+(vv-vg(j,k,l))*xyn+(ww-wg(j,k,l))*xzn
        
        bSq = Mp**2/(bt(j,k,l)-Mp**2*(bt(j,k,l)-1))
        
        Z1 = sqrt((1.-bSq)*qq*(1.-bSq)*qq+4.*bSq*cjkl*cjkl)
        Zi = 1./Z1
        R = 0.5*((1.-bSq)*qq+Z1)
        S = 0.5*((1.-bSq)*qq-Z1)

c..multiplying by sign(lamda_pa)*inv(X_pa)

        delw(1) = (xxn*(1.-uvw*gm1*c2i)-(vv*xzn-ww*xyn)*rrho)*
     c                                    sign(1.,qq)*delq(1)
        delw(1) = delw(1) +(uu*xxn*c2i*gm1)*sign(1.,qq)*delq(2)
        delw(1) = delw(1) +(vv*xxn*c2i*gm1+xzn*rrho)*sign(1.,qq)*delq(3)
        delw(1) = delw(1) +(ww*xxn*c2i*gm1-xyn*rrho)*sign(1.,qq)*delq(4)
        delw(1) = delw(1) - xxn*gm1*c2i*sign(1.,qq)*delq(5)
       
        delw(2) = (xyn*(1.-uvw*gm1*c2i)-(ww*xxn-uu*xzn)*rrho)*
     c                                    sign(1.,qq)*delq(1)
        delw(2) = delw(2) +(uu*xyn*c2i*gm1-xzn*rrho)*sign(1.,qq)*delq(2)
        delw(2) = delw(2) +(vv*xyn*c2i*gm1)*sign(1.,qq)*delq(3)
        delw(2) = delw(2) +(ww*xyn*c2i*gm1+xxn*rrho)*sign(1.,qq)*delq(4)
        delw(2) = delw(2) - xyn*gm1*c2i*sign(1.,qq)*delq(5)

        delw(3) = (xzn*(1.-uvw*gm1*c2i)-(uu*xyn-vv*xxn)*rrho)*
     c                                    sign(1.,qq)*delq(1)
        delw(3) = delw(3) +(uu*xzn*c2i*gm1+xyn*rrho)*sign(1.,qq)*delq(2)
        delw(3) = delw(3) +(vv*xzn*c2i*gm1-xxn*rrho)*sign(1.,qq)*delq(3)
        delw(3) = delw(3) +(ww*xzn*c2i*gm1)*sign(1.,qq)*delq(4)
        delw(3) = delw(3) - xzn*gm1*c2i*sign(1.,qq)*delq(5)
       
        delw(4) = (-S*uvw*gm1*c2i-qq*bSq)*Zi*sign(1.,R+qq*bSq)*delq(1)
        delw(4) = delw(4) + (xxn*bSq + S*uu*gm1*c2i)*Zi*
     c                      sign(1.,R+qq*bSq)*delq(2)
        delw(4) = delw(4) + (xyn*bSq + S*vv*gm1*c2i)*Zi*
     c                      sign(1.,R+qq*bSq)*delq(3)
        delw(4) = delw(4) + (xzn*bSq + S*ww*gm1*c2i)*Zi*
     c                      sign(1.,R+qq*bSq)*delq(4)
        delw(4) = delw(4) - S*gm1*c2i*Zi*sign(1.,R+qq*bSq)*delq(5)

        delw(5) = (R*uvw*gm1*c2i+qq*bSq)*Zi*sign(1.,S+qq*bSq)*delq(1)
        delw(5) = delw(5) - (xxn*bSq + R*uu*gm1*c2i)*Zi*
     c                      sign(1.,S+qq*bSq)*delq(2)
        delw(5) = delw(5) - (xyn*bSq + R*vv*gm1*c2i)*Zi*
     c                      sign(1.,S+qq*bSq)*delq(3)
        delw(5) = delw(5) - (xzn*bSq + R*ww*gm1*c2i)*Zi*
     c                      sign(1.,S+qq*bSq)*delq(4)
        delw(5) = delw(5) +  R*gm1*c2i*Zi*sign(1.,S+qq*bSq)*delq(5)

c.. multiplying by (X_pa)

        delq(1) = xxn*delw(1)+xyn*delw(2)+xzn*delw(3)+delw(4)+delw(5)

        delq(2) = uu*xxn*delw(1)+(uu*xyn-xzn*rho)*delw(2)
        delq(2) = delq(2) + (uu*xzn+xyn*rho)*delw(3)
        delq(2) = delq(2) + (uu+R*xxn/bSq)*delw(4)
        delq(2) = delq(2) + (uu+S*xxn/bSq)*delw(5)

        delq(3) = (vv*xxn+xzn*rho)*delw(1)+vv*xyn*delw(2)
        delq(3) = delq(3) + (vv*xzn-xxn*rho)*delw(3)
        delq(3) = delq(3) + (vv+R*xyn/bSq)*delw(4)
        delq(3) = delq(3) + (vv+S*xyn/bSq)*delw(5)

        delq(4) = (ww*xxn-xyn*rho)*delw(1)+(ww*xyn+xxn*rho)*delw(2)
        delq(4) = delq(4) + ww*xzn*delw(3)
        delq(4) = delq(4) + (ww+R*xzn/bSq)*delw(4)
        delq(4) = delq(4) + (ww+S*xzn/bSq)*delw(5)

        delq(5) = (uvw*xxn+rho*(vv*xzn-ww*xyn))*delw(1)
        delq(5) = delq(5) + (uvw*xyn+rho*(ww*xxn-uu*xzn))*delw(2)
        delq(5) = delq(5) + (uvw*xzn+rho*(uu*xyn-vv*xxn))*delw(3)
        delq(5) = delq(5) + (ge+R*qq/bSq)*delw(4)+(ge+S*qq/bSq)*delw(5)

c.. qbar - delq
        q(j,k,l,1) = (qbar(1)-delq(1))/q(j,k,l,6)
        q(j,k,l,2) = (qbar(2)-delq(2))/q(j,k,l,6)
        q(j,k,l,3) = (qbar(3)-delq(3))/q(j,k,l,6)
        q(j,k,l,4) = (qbar(4)-delq(4))/q(j,k,l,6)
        q(j,k,l,5) = (qbar(5)-delq(5))/q(j,k,l,6)

      enddo
      enddo

      elseif (iadir.eq.2) then

      k  = ks
      k1 = k + iadd
c     
      do  l = ls,le
      do  j = js,je
c..   free stream
        uind = uinf 
        vind = vinf 
        wind = winf 
        
        snorm = -1.*iadd/sqrt(yx(j,k,l)**2+yy(j,k,l)**2+yz(j,k,l)**2)
        yxn = yx(j,k,l)*snorm
        yyn = yy(j,k,l)*snorm
        yzn = yz(j,k,l)*snorm
        
c..zeroeth-order extrapolation
        qext(1) = q(j,k1,l,1)*q(j,k1,l,6)
        qext(2) = q(j,k1,l,2)*q(j,k1,l,6)
        qext(3) = q(j,k1,l,3)*q(j,k1,l,6)
        qext(4) = q(j,k1,l,4)*q(j,k1,l,6)
        qext(5) = q(j,k1,l,5)*q(j,k1,l,6)

        qbar(1) = 0.5*(qext(1)+rinf)
        qbar(2) = 0.5*(qext(2)+rinf*uind)
        qbar(3) = 0.5*(qext(3)+rinf*vind)
        qbar(4) = 0.5*(qext(4)+rinf*wind)
        qbar(5) = 0.5*(qext(5)+einf)

        delq(1) = 0.5*(rinf-qext(1))
        delq(2) = 0.5*(rinf*uind-qext(2))
        delq(3) = 0.5*(rinf*vind-qext(3))
        delq(4) = 0.5*(rinf*wind-qext(4))
        delq(5) = 0.5*(einf-qext(5))

        rho = qbar(1)
        rrho = 1./rho
        uu = qbar(2)*rrho
        vv = qbar(3)*rrho
        ww = qbar(4)*rrho
        e  = qbar(5)*rrho
        uvw = 0.5*(uu*uu+vv*vv+ww*ww)
        cjkl = sqrt(ggm1*(e-uvw))
        c2i = 1./(cjkl*cjkl)
        ge = gamma*e - gm1*uvw
        qq = (uu-ug(j,k,l))*yxn+(vv-vg(j,k,l))*yyn+(ww-wg(j,k,l))*yzn
        
        bSq = Mp**2/(bt(j,k,l)-Mp**2*(bt(j,k,l)-1))
        
        Z1 = sqrt((1.-bSq)*qq*(1.-bSq)*qq+4.*bSq*cjkl*cjkl)
        Zi = 1./Z1
        R = 0.5*((1.-bSq)*qq+Z1)
        S = 0.5*((1.-bSq)*qq-Z1)

c..multiplying by sign(lamda_pa)*inv(X_pa)

        delw(1) = (yxn*(1.-uvw*gm1*c2i)-(vv*yzn-ww*yyn)*rrho)*
     c                                    sign(1.,qq)*delq(1)
        delw(1) = delw(1) +(uu*yxn*c2i*gm1)*sign(1.,qq)*delq(2)
        delw(1) = delw(1) +(vv*yxn*c2i*gm1+yzn*rrho)*sign(1.,qq)*delq(3)
        delw(1) = delw(1) +(ww*yxn*c2i*gm1-yyn*rrho)*sign(1.,qq)*delq(4)
        delw(1) = delw(1) - yxn*gm1*c2i*sign(1.,qq)*delq(5)
       
        delw(2) = (yyn*(1.-uvw*gm1*c2i)-(ww*yxn-uu*yzn)*rrho)*
     c                                    sign(1.,qq)*delq(1)
        delw(2) = delw(2) +(uu*yyn*c2i*gm1-yzn*rrho)*sign(1.,qq)*delq(2)
        delw(2) = delw(2) +(vv*yyn*c2i*gm1)*sign(1.,qq)*delq(3)
        delw(2) = delw(2) +(ww*yyn*c2i*gm1+yxn*rrho)*sign(1.,qq)*delq(4)
        delw(2) = delw(2) - yyn*gm1*c2i*sign(1.,qq)*delq(5)

        delw(3) = (yzn*(1.-uvw*gm1*c2i)-(uu*yyn-vv*yxn)*rrho)*
     c                                    sign(1.,qq)*delq(1)
        delw(3) = delw(3) +(uu*yzn*c2i*gm1+yyn*rrho)*sign(1.,qq)*delq(2)
        delw(3) = delw(3) +(vv*yzn*c2i*gm1-yxn*rrho)*sign(1.,qq)*delq(3)
        delw(3) = delw(3) +(ww*yzn*c2i*gm1)*sign(1.,qq)*delq(4)
        delw(3) = delw(3) - yzn*gm1*c2i*sign(1.,qq)*delq(5)
       
        delw(4) = (-S*uvw*gm1*c2i-qq*bSq)*Zi*sign(1.,R+qq*bSq)*delq(1)
        delw(4) = delw(4) + (yxn*bSq + S*uu*gm1*c2i)*Zi*
     c                      sign(1.,R+qq*bSq)*delq(2)
        delw(4) = delw(4) + (yyn*bSq + S*vv*gm1*c2i)*Zi*
     c                      sign(1.,R+qq*bSq)*delq(3)
        delw(4) = delw(4) + (yzn*bSq + S*ww*gm1*c2i)*Zi*
     c                      sign(1.,R+qq*bSq)*delq(4)
        delw(4) = delw(4) - S*gm1*c2i*Zi*sign(1.,R+qq*bSq)*delq(5)

        delw(5) = (R*uvw*gm1*c2i+qq*bSq)*Zi*sign(1.,S+qq*bSq)*delq(1)
        delw(5) = delw(5) - (yxn*bSq + R*uu*gm1*c2i)*Zi*
     c                      sign(1.,S+qq*bSq)*delq(2)
        delw(5) = delw(5) - (yyn*bSq + R*vv*gm1*c2i)*Zi*
     c                      sign(1.,S+qq*bSq)*delq(3)
        delw(5) = delw(5) - (yzn*bSq + R*ww*gm1*c2i)*Zi*
     c                      sign(1.,S+qq*bSq)*delq(4)
        delw(5) = delw(5) +  R*gm1*c2i*Zi*sign(1.,S+qq*bSq)*delq(5)

c.. multiplying by (X_pa)

        delq(1) = yxn*delw(1)+yyn*delw(2)+yzn*delw(3)+delw(4)+delw(5)

        delq(2) = uu*yxn*delw(1)+(uu*yyn-yzn*rho)*delw(2)
        delq(2) = delq(2) + (uu*yzn+yyn*rho)*delw(3)
        delq(2) = delq(2) + (uu+R*yxn/bSq)*delw(4)
        delq(2) = delq(2) + (uu+S*yxn/bSq)*delw(5)

        delq(3) = (vv*yxn+yzn*rho)*delw(1)+vv*yyn*delw(2)
        delq(3) = delq(3) + (vv*yzn-yxn*rho)*delw(3)
        delq(3) = delq(3) + (vv+R*yyn/bSq)*delw(4)
        delq(3) = delq(3) + (vv+S*yyn/bSq)*delw(5)

        delq(4) = (ww*yxn-yyn*rho)*delw(1)+(ww*yyn+yxn*rho)*delw(2)
        delq(4) = delq(4) + ww*yzn*delw(3)
        delq(4) = delq(4) + (ww+R*yzn/bSq)*delw(4)
        delq(4) = delq(4) + (ww+S*yzn/bSq)*delw(5)

        delq(5) = (uvw*yxn+rho*(vv*yzn-ww*yyn))*delw(1)
        delq(5) = delq(5) + (uvw*yyn+rho*(ww*yxn-uu*yzn))*delw(2)
        delq(5) = delq(5) + (uvw*yzn+rho*(uu*yyn-vv*yxn))*delw(3)
        delq(5) = delq(5) + (ge+R*qq/bSq)*delw(4)+(ge+S*qq/bSq)*delw(5)

c.. qbar - delq
        q(j,k,l,1) = (qbar(1)-delq(1))/q(j,k,l,6)
        q(j,k,l,2) = (qbar(2)-delq(2))/q(j,k,l,6)
        q(j,k,l,3) = (qbar(3)-delq(3))/q(j,k,l,6)
        q(j,k,l,4) = (qbar(4)-delq(4))/q(j,k,l,6)
        q(j,k,l,5) = (qbar(5)-delq(5))/q(j,k,l,6)

      enddo
      enddo

      elseif(iadir.eq.3) then
      l  = ls
      l1 = l + iadd

c..remove l=lmax condition to check

      do j = js,je
      do k = ks,ke
c.free stream
        uind = uinf 
        vind = vinf 
        wind = winf 

        snorm = -1.*iadd/sqrt(zx(j,k,l)**2+zy(j,k,l)**2+zz(j,k,l)**2)
        zxn = zx(j,k,l)*snorm
        zyn = zy(j,k,l)*snorm
        zzn = zz(j,k,l)*snorm

c..zeroeth-order extrapolation
        qext(1) = q(j,k,l1,1)*q(j,k,l1,6)
        qext(2) = q(j,k,l1,2)*q(j,k,l1,6)
        qext(3) = q(j,k,l1,3)*q(j,k,l1,6)
        qext(4) = q(j,k,l1,4)*q(j,k,l1,6)
        qext(5) = q(j,k,l1,5)*q(j,k,l1,6)

        qbar(1) = 0.5*(qext(1)+rinf)
        qbar(2) = 0.5*(qext(2)+rinf*uind)
        qbar(3) = 0.5*(qext(3)+rinf*vind)
        qbar(4) = 0.5*(qext(4)+rinf*wind)
        qbar(5) = 0.5*(qext(5)+einf)

        delq(1) = 0.5*(rinf-qext(1))
        delq(2) = 0.5*(rinf*uind-qext(2))
        delq(3) = 0.5*(rinf*vind-qext(3))
        delq(4) = 0.5*(rinf*wind-qext(4))
        delq(5) = 0.5*(einf-qext(5))

        rho = qbar(1)
        rrho = 1./rho
        uu = qbar(2)*rrho
        vv = qbar(3)*rrho
        ww = qbar(4)*rrho
        e  = qbar(5)*rrho
        uvw = 0.5*(uu*uu+vv*vv+ww*ww)
        cjkl = sqrt(ggm1*(e-uvw))
        c2i = 1./(cjkl*cjkl)
        ge = gamma*e - gm1*uvw
        qq = (uu-ug(j,k,l))*zxn+(vv-vg(j,k,l))*zyn+(ww-wg(j,k,l))*zzn
        
        bSq = Mp**2/(bt(j,k,l)-Mp**2*(bt(j,k,l)-1))
        
        Z1 = sqrt((1.-bSq)*qq*(1.-bSq)*qq+4.*bSq*cjkl*cjkl)
        Zi = 1./Z1
        R = 0.5*((1.-bSq)*qq+Z1)
        S = 0.5*((1.-bSq)*qq-Z1)

c..multiplying by sign(lamda_pa)*inv(X_pa)
        delw(1) = (zxn*(1.-uvw*gm1*c2i)-(vv*zzn-ww*zyn)*rrho)*
     c                                    sign(1.,qq)*delq(1)
        delw(1) = delw(1) +(uu*zxn*c2i*gm1)*sign(1.,qq)*delq(2)
        delw(1) = delw(1) +(vv*zxn*c2i*gm1+zzn*rrho)*sign(1.,qq)*delq(3)
        delw(1) = delw(1) +(ww*zxn*c2i*gm1-zyn*rrho)*sign(1.,qq)*delq(4)
        delw(1) = delw(1) - zxn*gm1*c2i*sign(1.,qq)*delq(5)
       
        delw(2) = (zyn*(1.-uvw*gm1*c2i)-(ww*zxn-uu*zzn)*rrho)*
     c                                    sign(1.,qq)*delq(1)
        delw(2) = delw(2) +(uu*zyn*c2i*gm1-zzn*rrho)*sign(1.,qq)*delq(2)
        delw(2) = delw(2) +(vv*zyn*c2i*gm1)*sign(1.,qq)*delq(3)
        delw(2) = delw(2) +(ww*zyn*c2i*gm1+zxn*rrho)*sign(1.,qq)*delq(4)
        delw(2) = delw(2) - zyn*gm1*c2i*sign(1.,qq)*delq(5)

        delw(3) = (zzn*(1.-uvw*gm1*c2i)-(uu*zyn-vv*zxn)*rrho)*
     c                                    sign(1.,qq)*delq(1)
        delw(3) = delw(3) +(uu*zzn*c2i*gm1+zyn*rrho)*sign(1.,qq)*delq(2)
        delw(3) = delw(3) +(vv*zzn*c2i*gm1-zxn*rrho)*sign(1.,qq)*delq(3)
        delw(3) = delw(3) +(ww*zzn*c2i*gm1)*sign(1.,qq)*delq(4)
        delw(3) = delw(3) - zzn*gm1*c2i*sign(1.,qq)*delq(5)
       
        delw(4) = (-S*uvw*gm1*c2i-qq*bSq)*Zi*sign(1.,R+qq*bSq)*delq(1)
        delw(4) = delw(4) + (zxn*bSq + S*uu*gm1*c2i)*Zi*
     c                      sign(1.,R+qq*bSq)*delq(2)
        delw(4) = delw(4) + (zyn*bSq + S*vv*gm1*c2i)*Zi*
     c                      sign(1.,R+qq*bSq)*delq(3)
        delw(4) = delw(4) + (zzn*bSq + S*ww*gm1*c2i)*Zi*
     c                      sign(1.,R+qq*bSq)*delq(4)
        delw(4) = delw(4) - S*gm1*c2i*Zi*sign(1.,R+qq*bSq)*delq(5)

        delw(5) = (R*uvw*gm1*c2i+qq*bSq)*Zi*sign(1.,S+qq*bSq)*delq(1)
        delw(5) = delw(5) - (zxn*bSq + R*uu*gm1*c2i)*Zi*
     c                      sign(1.,S+qq*bSq)*delq(2)
        delw(5) = delw(5) - (zyn*bSq + R*vv*gm1*c2i)*Zi*
     c                      sign(1.,S+qq*bSq)*delq(3)
        delw(5) = delw(5) - (zzn*bSq + R*ww*gm1*c2i)*Zi*
     c                      sign(1.,S+qq*bSq)*delq(4)
        delw(5) = delw(5) +  R*gm1*c2i*Zi*sign(1.,S+qq*bSq)*delq(5)

c.. multiplying by (X_pa)

        delq(1) = zxn*delw(1)+zyn*delw(2)+zzn*delw(3)+delw(4)+delw(5)

        delq(2) = uu*zxn*delw(1)+(uu*zyn-zzn*rho)*delw(2)
        delq(2) = delq(2) + (uu*zzn+zyn*rho)*delw(3)
        delq(2) = delq(2) + (uu+R*zxn/bSq)*delw(4)
        delq(2) = delq(2) + (uu+S*zxn/bSq)*delw(5)

        delq(3) = (vv*zxn+zzn*rho)*delw(1)+vv*zyn*delw(2)
        delq(3) = delq(3) + (vv*zzn-zxn*rho)*delw(3)
        delq(3) = delq(3) + (vv+R*zyn/bSq)*delw(4)
        delq(3) = delq(3) + (vv+S*zyn/bSq)*delw(5)

        delq(4) = (ww*zxn-zyn*rho)*delw(1)+(ww*zyn+zxn*rho)*delw(2)
        delq(4) = delq(4) + ww*zzn*delw(3)
        delq(4) = delq(4) + (ww+R*zzn/bSq)*delw(4)
        delq(4) = delq(4) + (ww+S*zzn/bSq)*delw(5)

        delq(5) = (uvw*zxn+rho*(vv*zzn-ww*zyn))*delw(1)
        delq(5) = delq(5) + (uvw*zyn+rho*(ww*zxn-uu*zzn))*delw(2)
        delq(5) = delq(5) + (uvw*zzn+rho*(uu*zyn-vv*zxn))*delw(3)
        delq(5) = delq(5) + (ge+R*qq/bSq)*delw(4)+(ge+S*qq/bSq)*delw(5)

c.. qbar - delq
        q(j,k,l,1) = (qbar(1)-delq(1))/q(j,k,l,6)
        q(j,k,l,2) = (qbar(2)-delq(2))/q(j,k,l,6)
        q(j,k,l,3) = (qbar(3)-delq(3))/q(j,k,l,6)
        q(j,k,l,4) = (qbar(4)-delq(4))/q(j,k,l,6)
        q(j,k,l,5) = (qbar(5)-delq(5))/q(j,k,l,6)

      enddo
      enddo
      endif

      return
      end

c***********************************************************************
      subroutine bcparallel(q,vnut,js,je,ks,ke,ls,le,idir,iproc)
c     
c  for parallel runs   
c*********************************************************************
      use mpi_wrapper
      use params_global
c*********************************************************************

      implicit none
      
      real q(jmax,kmax,lmax,nd),vnut(jmax,kmax,lmax)
      integer js,je,ks,ke,ls,le,idir,iproc
      
c..   local variables
      
      integer j,k,l,n,iadir,iadd
      integer j1,k1,l1,j2,k2,l2,jmax1,kmax1,lmax1
      integer tag
      real scale1
      real,allocatable :: q1(:,:,:,:),vnut1(:,:,:)

      jmax1 = je-js+1
      kmax1 = ke-ks+1
      lmax1 = le-ls+1

      allocate(q1(jmax1,kmax1,lmax1,nd),vnut1(jmax1,kmax1,lmax1))
      
      iadd = sign(1,idir)
      iadir = abs(idir)

      tag = 0
      if (iadir.eq.1) then
         do j = js,je
         do k = ks,ke
         do l = ls,le
            j1 = j-js+1
            k1 = k-ks+1
            l1 = l-ls+1
            j2 = j+iadd*(jmax1+1)
            do n=1,6
               q1(j1,k1,l1,n) = q(j2,k,l,n)
            enddo
            vnut1(j1,k1,l1) = vnut(j2,k,l)
         enddo
         enddo
         enddo
      elseif (iadir.eq.2) then
         do j = js,je
         do k = ks,ke
         do l = ls,le
            j1 = j-js+1
            k1 = k-ks+1
            l1 = l-ls+1
            k2 = k+iadd*(kmax1+1)
            do n=1,6
               q1(j1,k1,l1,n) = q(j,k2,l,n)
            enddo
            vnut1(j1,k1,l1) = vnut(j,k2,l)
         enddo
         enddo
         enddo
      elseif (iadir.eq.3) then
         do j = js,je
         do k = ks,ke
         do l = ls,le
            j1 = j-js+1
            k1 = k-ks+1
            l1 = l-ls+1
            l2 = l+iadd*(lmax1+1)
            do n=1,6
               q1(j1,k1,l1,n) = q(j,k,l2,n)
            enddo
            vnut1(j1,k1,l1) = vnut(j,k,l2)
         enddo
         enddo
         enddo
      endif

      call mpi_bsend(q1,jmax1*kmax1*lmax1*nd,REAL_TYPE,iproc,tag,
     c               DEFAULT_COMM,ierr)
      call mpi_bsend(vnut1,jmax1*kmax1*lmax1,REAL_TYPE,iproc,tag,
     c               DEFAULT_COMM,ierr)
      call mpi_recv(q1,jmax1*kmax1*lmax1*nd,REAL_TYPE,iproc,tag,
     c              DEFAULT_COMM,stats_mpi,ierr)
      call mpi_recv(vnut1,jmax1*kmax1*lmax1,REAL_TYPE,iproc,tag,
     c              DEFAULT_COMM,stats_mpi,ierr)
         
      do j = js,je
      do k = ks,ke
      do l = ls,le
         j1 = j-js+1
         k1 = k-ks+1
         l1 = l-ls+1
         scale1 = q1(j1,k1,l1,6)/q(j,k,l,6)
         do n=1,5
            q(j,k,l,n) = q1(j1,k1,l1,n)*scale1
         enddo
         vnut(j,k,l) = vnut1(j1,k1,l1)
      enddo
      enddo
      enddo

      return
      end

c*********************************************************************
      subroutine bcparallel_bg(q,vnut,js,je,ks,ke,ls,le,idir,iproc)
c     
c  for parallel runs   
c*********************************************************************
      use mpi_wrapper
      use params_global
c*********************************************************************

      implicit none
      
      real q(jmax,kmax,lmax,nd),vnut(jmax,kmax,lmax)
      integer js,je,ks,ke,ls,le,idir,iproc
      
c..   local variables
      
      integer j,k,l,n,iadir,iadd
      integer j1,k1,l1,j2,k2,l2,jmax1,kmax1,lmax1
      integer tag
      real scale1
      real tcos,tsin
      real,allocatable :: q1(:,:,:,:),vnut1(:,:,:)

      jmax1 = je-js+1
      kmax1 = ke-ks+1
      lmax1 = le-ls+1

      allocate(q1(jmax1,kmax1,lmax1,nd),vnut1(jmax1,kmax1,lmax1))
      
      iadd = sign(1,idir)
      iadir = abs(idir)

      ! rotation terms
      pi = 4*atan(1.0)
      blang = 2*pi/float(nblade)
      tcos  = cos(blang)
      tsin  = sin(blang)

      tag = 0
      if (iadir.eq.1) then
         do j = js,je
         do k = ks,ke
         do l = ls,le
            j1 = j-js+1
            k1 = k-ks+1
            l1 = l-ls+1
            j2 = j+iadd*(jmax1+1)
            q1(j1,k1,l1,1) = q(j2,k,l,1)
            q1(j1,k1,l1,2) = q(j2,k,l,2)*tcos - q(j2,k,l,3)*tsin
            q1(j1,k1,l1,3) = q(j2,k,l,2)*tsin + q(j2,k,l,3)*tcos
            q1(j1,k1,l1,4) = q(j2,k,l,4)
            q1(j1,k1,l1,5) = q(j2,k,l,5)
            q1(j1,k1,l1,6) = q(j2,k,l,6)
            vnut1(j1,k1,l1) = vnut(j2,k,l)
         enddo
         enddo
         enddo
      elseif (iadir.eq.2) then
         do j = js,je
         do k = ks,ke
         do l = ls,le
            j1 = j-js+1
            k1 = k-ks+1
            l1 = l-ls+1
            k2 = k+iadd*(kmax1+1)
            q1(j1,k1,l1,1) = q(j,k2,l,1)
            q1(j1,k1,l1,2) = q(j,k2,l,2)*tcos - q(j,k2,l,3)*tsin
            q1(j1,k1,l1,3) = q(j,k2,l,2)*tsin + q(j,k2,l,3)*tcos
            q1(j1,k1,l1,4) = q(j,k2,l,4)
            q1(j1,k1,l1,5) = q(j,k2,l,5)
            q1(j1,k1,l1,6) = q(j,k2,l,6)
            vnut1(j1,k1,l1) = vnut(j,k2,l)
         enddo
         enddo
         enddo
      elseif (iadir.eq.3) then
         do j = js,je
         do k = ks,ke
         do l = ls,le
            j1 = j-js+1
            k1 = k-ks+1
            l1 = l-ls+1
            l2 = l+iadd*(lmax1+1)
            q1(j1,k1,l1,1) = q(j,k,l2,1)
            q1(j1,k1,l1,2) = q(j,k,l2,2)*tcos - q(j,k,l2,3)*tsin
            q1(j1,k1,l1,3) = q(j,k,l2,2)*tsin + q(j,k,l2,3)*tcos
            q1(j1,k1,l1,4) = q(j,k,l2,4)
            q1(j1,k1,l1,5) = q(j,k,l2,5)
            q1(j1,k1,l1,6) = q(j,k,l2,6)
            vnut1(j1,k1,l1) = vnut(j,k,l2)
         enddo
         enddo
         enddo
      endif

      call mpi_bsend(q1,jmax1*kmax1*lmax1*nd,REAL_TYPE,iproc,tag,
     c               DEFAULT_COMM,ierr)
      call mpi_bsend(vnut1,jmax1*kmax1*lmax1,REAL_TYPE,iproc,tag,
     c               DEFAULT_COMM,ierr)
      call mpi_recv(q1,jmax1*kmax1*lmax1*nd,REAL_TYPE,iproc,tag,
     c              DEFAULT_COMM,stats_mpi,ierr)
      call mpi_recv(vnut1,jmax1*kmax1*lmax1,REAL_TYPE,iproc,tag,
     c              DEFAULT_COMM,stats_mpi,ierr)
         
      do j = js,je
      do k = ks,ke
      do l = ls,le
         j1 = j-js+1
         k1 = k-ks+1
         l1 = l-ls+1
         scale1 = q1(j1,k1,l1,6)/q(j,k,l,6)
         do n=1,5
            q(j,k,l,n) = q1(j1,k1,l1,n)*scale1
         enddo
         vnut(j,k,l) = vnut1(j1,k1,l1)
      enddo
      enddo
      enddo

      return
      end

c**********************************************************************
      subroutine bccoaxial(q,vnut,x,y,z,js,je,ks,ke,ls,le,idir,iproc)
c     
c  for parallel runs   
c*********************************************************************
      use mpi_wrapper
      use params_global
c*********************************************************************

      implicit none
      
      real q(jmax,kmax,lmax,nd),vnut(jmax,kmax,lmax)
      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)
      integer js,je,ks,ke,ls,le,idir,iproc
      
c..   local variables
      
      integer j,k,l,n,iadir,iadd
      integer j1,k1,l1,j2,k2,l2,jmax1,kmax1,lmax1
      integer jmax2,kmax2,lmax2,ll
      integer tag
      integer j2bndry,nrot_spl
      real psitmp,psitmp1,psitmp2,dpsi,dpsi1,dpsi2,tcos,tsin,u,v
      real,allocatable :: q1(:,:,:,:),vnut1(:,:,:),z1(:)
      real,allocatable :: psigrid(:),psigrid1(:),fac1(:),fac2(:)
      integer,allocatable :: nrot1(:),nrot2(:),ji(:)

      jmax1 = je-js+1
      kmax1 = ke-ks+1
      lmax1 = le-ls+1

      allocate(q1(jmax1,kmax1,lmax1,nd),vnut1(jmax1,kmax1,lmax1))
      
      iadd = sign(1,idir)
      iadir = abs(idir)

      tag = 0
      call mpi_bsend(jmax,1,mpi_integer,iproc,tag,
     c              DEFAULT_COMM,ierr)
      call mpi_bsend(kmax,1,mpi_integer,iproc,tag,
     c              DEFAULT_COMM,ierr)
      call mpi_bsend(lmax,1,mpi_integer,iproc,tag,
     c              DEFAULT_COMM,ierr)

      call mpi_recv(jmax2,1,mpi_integer,iproc,tag,
     c              DEFAULT_COMM,stats_mpi,ierr)
      call mpi_recv(kmax2,1,mpi_integer,iproc,tag,
     c              DEFAULT_COMM,stats_mpi,ierr)
      call mpi_recv(lmax2,1,mpi_integer,iproc,tag,
     c              DEFAULT_COMM,stats_mpi,ierr)

      if (iadir.eq.1.or.iadir.eq.2) then
        print*,'idir = ',idir,' is not implemented in bccoaxial'

      elseif (iadir.eq.3) then
         allocate(z1(lmax))
         do l=1,lmax
            z1(l) = z(js,ks,l)
         enddo
         call mpi_bsend(z1,lmax,REAL_TYPE,iproc,tag,
     c                 DEFAULT_COMM,ierr)
         deallocate(z1)
         allocate(z1(lmax2))
         call mpi_recv(z1,lmax2,REAL_TYPE,iproc,tag,
     c                 DEFAULT_COMM,stats_mpi,ierr)

         do l=lmax2,1,-1
           if(z(js,ks,ls).eq.z1(l)) then
              ll = l
              exit
           endif
         enddo
         call mpi_bsend(ll,1,mpi_integer,iproc,tag,
     c                 DEFAULT_COMM,ierr)
         call mpi_recv(ll,1,mpi_integer,iproc,tag,
     c                 DEFAULT_COMM,stats_mpi,ierr)

         do j = js,je
         do k = ks,ke
         do l = ls,le
            j1 = j-js+1
            k1 = k-ks+1
            l1 = l-ls+1
            l2 = ll+l1-1
            do n=1,6
               q1(j1,k1,l1,n) = q(j,k,l2,n)
            enddo
            vnut1(j1,k1,l1) = vnut(j,k,l2)
         enddo
         enddo
         enddo

         call mpi_bsend(q1,jmax1*kmax1*lmax1*nd,REAL_TYPE,iproc,tag,
     c                 DEFAULT_COMM,ierr)
         call mpi_bsend(vnut1,jmax1*kmax1*lmax1,REAL_TYPE,iproc,tag,
     c                 DEFAULT_COMM,ierr)
         call mpi_recv(q1,jmax1*kmax1*lmax1*nd,REAL_TYPE,iproc,tag,
     c                 DEFAULT_COMM,stats_mpi,ierr)
         call mpi_recv(vnut1,jmax1*kmax1*lmax1,REAL_TYPE,iproc,tag,
     c                 DEFAULT_COMM,stats_mpi,ierr)

         allocate(psigrid(jmax1),psigrid1(jmax1))
         allocate(fac1(jmax1),fac2(jmax1))
         allocate(nrot1(jmax1),nrot2(jmax1))
         allocate(ji(jmax1))
            
c...finding the azimuth location of each grid point in cylindrical mesh
c   along the with azimuth location of cylindrical mesh of other blade

         nrot1 = 0
         do j = js,je
            j1 = j-js+1
            psigrid(j1) = atan(y(j,ke,ls)/x(j,ke,ls))
            if(x(j,ke,ls).lt.0) then
               psigrid(j1) = psigrid(j1) + pi
            endif
            do while (psigrid(j1).gt.2*pi/Nblade)
               psigrid(j1) = psigrid(j1)-2*pi/Nblade
               nrot1(j1) = nrot1(j1)+1
            enddo
            do while (psigrid(j1).lt.0)
               psigrid(j1) = psigrid(j1)+2*pi/Nblade
               nrot1(j1) = nrot1(j1)-1
            enddo
         enddo

         dpsi = -2*rf*totime
         nrot2 = -nrot1
         do j = js,je
            j1 = j-js+1
            psigrid1(j1) = psigrid(j1)+dpsi
            do while (psigrid1(j1).gt.2*pi/Nblade)
               psigrid1(j1) = psigrid1(j1)-2*pi/Nblade
               nrot2(j1) = nrot2(j1)-1
            enddo
            do while (psigrid1(j1).lt.0)
               psigrid1(j1) = psigrid1(j1)+2*pi/Nblade
               nrot2(j1) = nrot2(j1)+1
            enddo
         enddo

c...this part assumes that psi(j) is greater than psi(j+1) for cylindrical mesh
c  and finds the interpolation factor

         j2bndry = 0
         do j1=1,jmax1
         nrot_spl = 0
         do j2=1,jmax1-1
            psitmp = psigrid(j1)
            psitmp1 = psigrid1(j2)
            psitmp2 = psigrid1(j2+1)
            if(psitmp2.gt.pi/Nblade.and.psitmp1.lt.pi/Nblade) then
               j2bndry = j2+1
               psitmp2 = psitmp2 - 2*pi/Nblade
               if(psitmp.gt.pi/Nblade) then
                  psitmp = psitmp - 2*pi/Nblade
                  nrot_spl = 1
               endif
            endif
            if(psitmp.le.psitmp1.and.psitmp.ge.psitmp2) then
               dpsi1 = psitmp1 - psitmp 
               dpsi2 = psitmp - psitmp2
               fac1(j1) = dpsi1/(dpsi1+dpsi2)
               fac2(j1) = dpsi2/(dpsi1+dpsi2)
               ji(j1) = j2
               nrot1(j1) = nrot1(j1)+nrot_spl
               go to 110
            endif
            nrot_spl = 0
         enddo
  110    continue
         enddo
         
         do j = js,je
         do k = ks,ke
         do l = ls,le
            j1 = j-js+1
            k1 = k-ks+1
            l1 = l-ls+1
            if(nrot2(j1).ne.0) then
            n = nrot2(j1)
            tcos  = cos(n*blang)
            tsin  = sin(n*blang)
            u = q1(j1,k1,l1,2)
            v = q1(j1,k1,l1,3)
            q1(j1,k1,l1,2) = u*tcos - v*tsin
            q1(j1,k1,l1,3) = u*tsin + v*tcos
            endif
         enddo 
         enddo 
         enddo 

c...computing the interpolated solution

         do j = js,je
         do k = ks,ke
         do l = ls,le

         j1 = j-js+1
         k1 = k-ks+1
         l1 = l-ls+1

         if(ji(j1).eq.j2bndry-1) then
            tcos  = cos(blang)
            tsin  = sin(blang)
            u = q1(j2bndry,k1,l1,2)
            v = q1(j2bndry,k1,l1,3)
            q1(j2bndry,k1,l1,2) = u*tcos + v*tsin
            q1(j2bndry,k1,l1,3) =-u*tsin + v*tcos
         endif

         do n = 1,5
           q(j,k,l,n)=fac1(j1)*q1(ji(j1)+1,k1,l1,n)*q1(ji(j1)+1,k1,l1,6)+
     c                fac2(j1)*q1(ji(j1),k1,l1,n)*q1(ji(j1),k1,l1,6)
           q(j,k,l,n) = q(j,k,l,n)/q(j,k,l,6)
         enddo

         if(ji(j1).eq.j2bndry-1) then
            tcos  = cos(blang)
            tsin  = sin(blang)
            u = q1(j2bndry,k1,l1,2)
            v = q1(j2bndry,k1,l1,3)
            q1(j2bndry,k1,l1,2) = u*tcos - v*tsin
            q1(j2bndry,k1,l1,3) = u*tsin + v*tcos
         endif

         if(nrot1(j1).ne.0) then
           n = nrot1(j1)
           tcos  = cos(n*blang)
           tsin  = sin(n*blang)
           u = q(j,k,l,2)
           v = q(j,k,l,3)
           q(j,k,l,2) = u*tcos - v*tsin
           q(j,k,l,3) = u*tsin + v*tcos
         endif

         vnut(j,k,l)=fac1(j1)*vnut1(ji(j1)+1,k1,l1)+
     c               fac2(j1)*vnut1(ji(j1),k1,l1)
         enddo
         enddo
         enddo

      endif

      return
      end

c**********************************************************************
      subroutine bccoaxialcub(q,vnut,x,y,z,js,je,ks,ke,ls,le,idir,iproc)
c     
c  for parallel runs   
c*********************************************************************
      use mpi_wrapper
      use params_global
c*********************************************************************

      implicit none
      
      real q(jmax,kmax,lmax,nd),vnut(jmax,kmax,lmax)
      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)
      integer js,je,ks,ke,ls,le,idir,iproc
      
c..   local variables
      
      integer j,k,l,n,iadir,iadd
      integer j1,k1,l1,j2,k2,l2,jmax1,kmax1,lmax1,jp1,jp2
      integer jmax2,kmax2,lmax2,ll
      integer tag
      integer ibndry,nrot_spl
      real psitmp,psitmp1,psitmp2,dpsi,ds,dpsi1,dpsi2,tcos,tsin,u,v
      real,allocatable :: q1(:,:,:,:),vnut1(:,:,:),z1(:)
      real,allocatable :: qtmp(:,:,:,:)
      real,allocatable :: psigrid(:),psigrid1(:)
      integer,allocatable :: nrot1(:),nrot2(:),ji(:)

      jmax1 = je-js+1
      kmax1 = ke-ks+1
      lmax1 = le-ls+1

      allocate(q1(jmax1,kmax1,lmax1,nd),vnut1(jmax1,kmax1,lmax1))
      
      iadd = sign(1,idir)
      iadir = abs(idir)

      tag = 0
      call mpi_bsend(jmax,1,mpi_integer,iproc,tag,
     c              DEFAULT_COMM,ierr)
      call mpi_bsend(kmax,1,mpi_integer,iproc,tag,
     c              DEFAULT_COMM,ierr)
      call mpi_bsend(lmax,1,mpi_integer,iproc,tag,
     c              DEFAULT_COMM,ierr)

      call mpi_recv(jmax2,1,mpi_integer,iproc,tag,
     c              DEFAULT_COMM,stats_mpi,ierr)
      call mpi_recv(kmax2,1,mpi_integer,iproc,tag,
     c              DEFAULT_COMM,stats_mpi,ierr)
      call mpi_recv(lmax2,1,mpi_integer,iproc,tag,
     c              DEFAULT_COMM,stats_mpi,ierr)

      if (iadir.eq.1.or.iadir.eq.2) then
        print*,'idir = ',idir,' is not implemented in bccoaxialcub'

      elseif (iadir.eq.3) then
        allocate(z1(lmax))
        do l=1,lmax
           z1(l) = z(js,ks,l)
        enddo
        call mpi_bsend(z1,lmax,REAL_TYPE,iproc,tag,
     c                DEFAULT_COMM,ierr)
        deallocate(z1)
        allocate(z1(lmax2))
        call mpi_recv(z1,lmax2,REAL_TYPE,iproc,tag,
     c                DEFAULT_COMM,stats_mpi,ierr)

        do l=lmax2,1,-1
          if(z(js,ks,ls).eq.z1(l)) then
             ll = l
             exit
          endif
        enddo
        call mpi_bsend(ll,1,mpi_integer,iproc,tag,
     c                DEFAULT_COMM,ierr)
        call mpi_recv(ll,1,mpi_integer,iproc,tag,
     c                DEFAULT_COMM,stats_mpi,ierr)

        do j = js,je
        do k = ks,ke
        do l = ls,le
          j1 = j-js+1
          k1 = k-ks+1
          l1 = l-ls+1
          l2 = ll+l1-1
          do n=1,6
             q1(j1,k1,l1,n) = q(j,k,l2,n)
          enddo
          vnut1(j1,k1,l1) = vnut(j,k,l2)
        enddo
        enddo
        enddo

        call mpi_bsend(q1,jmax1*kmax1*lmax1*nd,REAL_TYPE,iproc,tag,
     c                DEFAULT_COMM,ierr)
        call mpi_bsend(vnut1,jmax1*kmax1*lmax1,REAL_TYPE,iproc,tag,
     c                DEFAULT_COMM,ierr)
        call mpi_recv(q1,jmax1*kmax1*lmax1*nd,REAL_TYPE,iproc,tag,
     c                DEFAULT_COMM,stats_mpi,ierr)
        call mpi_recv(vnut1,jmax1*kmax1*lmax1,REAL_TYPE,iproc,tag,
     c                DEFAULT_COMM,stats_mpi,ierr)

        allocate(qtmp(jmax1+4,kmax1,lmax1,6))
        allocate(psigrid(jmax1),psigrid1(jmax1))
        allocate(nrot1(jmax1),nrot2(jmax1))
        allocate(ji(jmax1))
            
c...finding the azimuth location of each grid point in cylindrical mesh
c   along the with azimuth location of cylindrical mesh of other blade

        nrot1 = 0
        do j = js,je
          j1 = j-js+1
          psigrid(j1) = atan(y(j,ke,ls)/x(j,ke,ls))
          if(x(j,ke,ls).lt.0) then
             psigrid(j1) = psigrid(j1) + pi
          endif
          do while (psigrid(j1).gt.2*pi/Nblade)
             psigrid(j1) = psigrid(j1)-2*pi/Nblade
             nrot1(j1) = nrot1(j1)+1
          enddo
          do while (psigrid(j1).lt.0)
             psigrid(j1) = psigrid(j1)+2*pi/Nblade
             nrot1(j1) = nrot1(j1)-1
          enddo
        enddo

        dpsi = -2*rf*totime
        nrot2 = -nrot1
        do j = js,je
           j1 = j-js+1
           psigrid1(j1) = psigrid(j1)+dpsi
           do while (psigrid1(j1).gt.2*pi/Nblade)
              psigrid1(j1) = psigrid1(j1)-2*pi/Nblade
              nrot2(j1) = nrot2(j1)-1
           enddo
           do while (psigrid1(j1).lt.0)
              psigrid1(j1) = psigrid1(j1)+2*pi/Nblade
              nrot2(j1) = nrot2(j1)+1
           enddo
        enddo

c...finding the angle between two adjacent grid lines
c...Here we assume that the gridlines are equispaced.

        dpsi = psigrid(1) - psigrid(2)
        if(dpsi.lt.0.) dpsi = 2*pi/Nblade + dpsi

c...finding ds for cubic interpolation

        ds = mod(2*pi+mod(2*rf*totime,2*pi),dpsi)
        if(ds.gt.1e-6) ds = dpsi - ds
        if(ds.le.1e-6) ds = 0.
        ds = ds/dpsi
       
c...Assigning extra gridlines using periodicity to qtmp 
c...to have 3rd interpolation even at boundaries.

        do j = 1,jmax1-1
          if(abs(psigrid1(jmax1)-psigrid1(j)).le.1e-6.or.
     c       abs(psigrid1(jmax1)-psigrid1(j)+2*pi/Nblade).le.1e-6) then
            jp1 = j
            go to 108
          endif
        enddo
108     continue

        do j = jmax1,2,-1
          if(abs(psigrid1(1)-psigrid1(j)).le.1e-6.or.
     c       abs(psigrid1(1)-psigrid1(j)-2*pi/Nblade).le.1e-6) then
            jp2 = j
            go to 109
          endif
        enddo
109     continue

        tcos  = cos(blang)
        tsin  = sin(blang)
        do k = ks,ke
        do l = ls,le
          k1 = k-ks+1
          l1 = l-ls+1

          do n = 1,5
            qtmp(1,k1,l1,n) = q1(jp2-2,k1,l1,n)
            qtmp(2,k1,l1,n) = q1(jp2-1,k1,l1,n)
          enddo
          u = q1(jp2-2,k1,l1,2)
          v = q1(jp2-2,k1,l1,3)
          qtmp(1,k1,l1,2) = u*tcos - v*tsin
          qtmp(1,k1,l1,3) = u*tsin + v*tcos
          u = q1(jp2-1,k1,l1,2)
          v = q1(jp2-1,k1,l1,3)
          qtmp(2,k1,l1,2) = u*tcos - v*tsin
          qtmp(2,k1,l1,3) = u*tsin + v*tcos
          qtmp(1,k1,l1,6) = vnut1(jp2-2,k1,l1)
          qtmp(2,k1,l1,6) = vnut1(jp2-1,k1,l1)

          do j = js,je
            j1 = j-js+1
            do n = 1,5
              qtmp(j1+2,k1,l1,n) = q1(j1,k1,l1,n)
            enddo
            qtmp(j1+2,k1,l1,6) = vnut1(j1,k1,l1)
          enddo

          do n = 1,5
            qtmp(jmax1+3,k1,l1,n) = q1(jp1+1,k1,l1,n)
            qtmp(jmax1+4,k1,l1,n) = q1(jp1+2,k1,l1,n)
          enddo
          u = q1(jp1+1,k1,l1,2)
          v = q1(jp1+1,k1,l1,3)
          qtmp(jmax1+3,k1,l1,2) = u*tcos + v*tsin
          qtmp(jmax1+3,k1,l1,3) = -u*tsin + v*tcos
          u = q1(jp1+2,k1,l1,2)
          v = q1(jp1+2,k1,l1,3)
          qtmp(jmax1+4,k1,l1,2) = u*tcos + v*tsin
          qtmp(jmax1+4,k1,l1,3) = -u*tsin + v*tcos
          qtmp(jmax1+3,k1,l1,6) = vnut1(jp1+1,k1,l1)
          qtmp(jmax1+4,k1,l1,6) = vnut1(jp1+2,k1,l1)
        enddo
        enddo

c...doing cubic interpolation

        do k = ks,ke
        do l = ls,le
          k1 = k-ks+1
          l1 = l-ls+1
          do n = 1,5
            call cubic_interpolate(qtmp(:,k1,l1,n),q1(:,k1,l1,n),ds,jmax1+4,
     c                        jmax1,3,jmax1+2)
          enddo
          call cubic_interpolate(qtmp(:,k1,l1,6),vnut1(:,k1,l1),ds,jmax1+4,
     c                      jmax1,3,jmax1+2)
        enddo
        enddo

c...this part assumes that psi(j) is greater than psi(j+1) for cylindrical mesh
c...and finds the donor grid point. 

        do j1=1,jmax1
        nrot_spl = 0
        do j2=1,jmax1
          psitmp = psigrid(j1)
          psitmp1 = psigrid1(j2)
          psitmp2 = psigrid1(j2+1)

          ibndry = 0
          if(psitmp2.gt.pi/Nblade.and.psitmp1.lt.pi/Nblade) then
            psitmp2 = psitmp2 - 2*pi/Nblade
            ibndry = 1
            if(psitmp.gt.pi/Nblade) then
              psitmp = psitmp - 2*pi/Nblade
              nrot_spl = 1
            endif
          endif

          if(psitmp.le.psitmp1.and.psitmp.ge.psitmp2) then
            dpsi1 = psitmp1 - psitmp 
            dpsi2 = psitmp - psitmp2
            ji(j1) = j2
            if(dpsi2/(dpsi1+dpsi2).le.1e-6) then
              ji(j1) = j2+1
              if(ibndry.eq.1) nrot_spl = nrot_spl-1
            endif
            nrot1(j1) = nrot1(j1)+nrot_spl
            go to 110
          endif
          nrot_spl = 0
        enddo
  110   continue
        enddo
        
        do j = js,je
        do k = ks,ke
        do l = ls,le
          j1 = j-js+1
          k1 = k-ks+1
          l1 = l-ls+1
          if(nrot2(j1).ne.0) then
            n = nrot2(j1)
            tcos  = cos(n*blang)
            tsin  = sin(n*blang)
            u = q1(j1,k1,l1,2)
            v = q1(j1,k1,l1,3)
            q1(j1,k1,l1,2) = u*tcos - v*tsin
            q1(j1,k1,l1,3) = u*tsin + v*tcos
          endif
        enddo 
        enddo 
        enddo 

c...assigning the interpolated solution

        do j = js,je
        do k = ks,ke
        do l = ls,le

        j1 = j-js+1
        k1 = k-ks+1
        l1 = l-ls+1

        do n = 1,5
          q(j,k,l,n)=q1(ji(j1),k1,l1,n)*q1(ji(j1),k1,l1,6)
          q(j,k,l,n)=q(j,k,l,n)/q(j,k,l,6)
        enddo

        if(nrot1(j1).ne.0) then
          n = nrot1(j1)
          tcos  = cos(n*blang)
          tsin  = sin(n*blang)
          u = q(j,k,l,2)
          v = q(j,k,l,3)
          q(j,k,l,2) = u*tcos - v*tsin
          q(j,k,l,3) = u*tsin + v*tcos
        endif

        vnut(j,k,l)=vnut1(ji(j1),k1,l1)
        enddo
        enddo
        enddo

      endif

      return
      end

c**********************************************************************
      subroutine cubic_interpolate(f,fint,ds,mdim,mdim1,is,ie)
c
c m3-quartic scheme of hyunh et al to determine the slope
c                                               and midpoints
c*************************************************************************
      implicit none

      integer :: mdim,mdim1,is,ie
      real    :: ds

      real f(mdim),fint(mdim1)

      ! local variables

      real ammd,fl,fr
      integer i,i1
      real at,ati1,s1,t1,sl,tmax
      real,allocatable :: f0(:),f1(:),f2(:),f2m(:),quar(:),fdot(:)
c*************************************************************************
      ammd(fl,fr) = 0.5*(sign(1.,fl)+sign(1.,fr))*amin1(abs(fl),abs(fr))
c*************************************************************************

      allocate(f0(mdim),f1(mdim),f2(mdim),f2m(mdim),quar(mdim),fdot(mdim))
c***  first executable statement

c..the m3-quartic interpolation scheme of hyunh follows


c..let's load up a few difference arrays
c
        do i=1,mdim
          f0(i) = f(i)
        enddo
c..1st difference at i+1/2
        do i=1,mdim-1
          f1(i) = f0(i+1) - f0(i)
        enddo
c..2nd difference at i
        do i=2,mdim-1
          f2(i)  = f1(i) - f1(i-1)
        enddo
c..limit 2nd difference to i+1/2
        do i = 2,mdim-2
          f2m(i) = ammd(f2(i),f2(i+1))
        enddo
c..quartic slope
      do i=is,ie
        quar(i) = (+f(i-2)-8.*f(i-1)+8.*f(i+1)-f(i+2))/12.
      enddo
c
c..now combine everything to get the slope
c
      do i=is,ie
c..include limited curvatures to calculate new slopes
        at    = f1(i)   - 0.5*f2m(i)
        ati1  = f1(i-1) + 0.5*f2m(i-1)
c..limit slopes at i
        s1    = ammd(f1(i),f1(i-1))
        t1    = ammd(at,ati1)
c..now find appropriate slope
c       sl = 0.5*(at+ati1)
        sl = quar(i)
        sl = sl+ammd(at-sl,ati1-sl)
        tmax  = sign(1.,t1)*amax1(3.*abs(s1),1.5*abs(t1))
        sl = ammd(sl,tmax)
c..use slope to calculate ql and qr
        fdot(i) = sl
      enddo
c
      do i=is,ie-1
        i1 = i-is+1
        fint(i1) = f(i)+fdot(i)*ds+(3.*f1(i)-2.*fdot(i)-fdot(i+1))*ds**2
     <           +(fdot(i)+fdot(i+1)-2.*f1(i))*ds**3
      enddo
      i = ie
      i1 = i-is+1
      fint(i1) = f(i)+fdot(i)*ds
c
      return
      end

c***********************************************************************
