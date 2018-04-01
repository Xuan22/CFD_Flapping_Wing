c**********************************************************************
      subroutine formom (lsurf,ct,cq,figofmerit,x,y,z,q,
     &        xx,xy,xz,yx,yy,yz,zx,zy,zz,twist,ldidx,ifile,ifile1,ifile2)
c
c  this subroutine evaluates the force and moment coefficients.
c  the nondimensional coefficient variables are defined below:
c
c                ****  cx  = force in x-direction ****
c                ****  cy  = force in y-direction ****
c                ****  cz  = force in z-direction ****
c                ****  cl  = lift                 ****
c                ****  cd  = drag                 ****
c
c    a constant laminar viscosity = 1.0 is assumed.  the moments are
c  first taken about cartesian axies centered at the origin and then
c  shifted to an axis system whose origin is at (x0,y0,z0).  the thin
c  layer viscous stress tensor is also valid for a general newtonian
c  fluid where the no slip condition on a surface applies.  the
c  reynolds number used in this subroutine is based on the free
c  stream speed of sound.
c                                                        (12/09/85)
c                                                       --- nmc ---
c j1,j2,k1,k2,lsurf ... j,k,l indicies defining surface
c x0,y0,z0 ... origin of moment reference system
c cx,cy,cz,cl,cd ... force and moment coefficients
c x,y,z ... grid
c q(j,k,l,6) ... conservative variables scaled by jacobian and jacobian
c re, fsmach, sref, alp, gamma, gami ... reynolds number (based on free
c     stream speed of sound), free stream mach number, reference area,
c     angle of attack, ratio of specific heats, gamma - 1
c ksym = 1 if periodicity in k, =0 otherwise
c invisc = true if inviscid, = false if viscous 
c
c**********************************************************************
      use params_global
      use deflections
c      use pyMOD
      use arf_mod
      use rootVariables, only: psi_rot
      use domain_info, only: send_lims

      implicit none
!
!     subroutine arguments
!
      integer lsurf,ifile,ifile1,ifile2,ldidx
      real cx,cy,cz
!   add it by camli 11/15/2013
      real cp
!  add power coeff.
      real ct,cq,figofmerit
      real q(jmax,kmax,lmax,6),x(jmax,kmax,lmax)
      real y(jmax,kmax,lmax),z(jmax,kmax,lmax)
      real xx(jmax,kmax,lmax),xy(jmax,kmax,lmax),xz(jmax,kmax,lmax)
      real yx(jmax,kmax,lmax),yy(jmax,kmax,lmax),yz(jmax,kmax,lmax)
      real zx(jmax,kmax,lmax),zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)
      real twist(kmax)
!     
!     local variables
!

!     allocatable arrays

      real, allocatable :: sxx(:,:),syy(:,:),szz(:,:)
      real, allocatable :: sxy(:,:),sxz(:,:),syz(:,:)
      real, allocatable :: vmu(:,:)
      real, allocatable :: cxsec(:),cmsec(:),czsec(:),cysec(:)

!     floats

      real rt(3,3)
      real c2b,c2bp,gigm1,gami,sw,sf,ca,sa,cost,sint
      real area,xxi,yxi,zxi,xet,yet,zet,ax,ay,az,sref
      real dj,rhoji,rho,u,v,w,e,qs,p,ei,ttinf,rei
!  added by camlib 11/15/2013 defining velx,vely,velyz & rhij2
      real velx,vely,velz,rhoij2
! done adding 
      real q1i1,q1i2,q1i3
      real usi,vsi,wsi
      real ueta,veta,weta
      real uzet,vzet,wzet,sfdiv
      real ux,uy,uz,vx,vy,vz,wx,wy,wz
      real xjle,yjle,zjle,xjtp,yjtp,zjtp
      real xptu,yptu,zptu,xptl,yptl,zptl
      real dx,dy,dz,dnorm
      real cchord,cvecx,cvecy,cvecz
      real pnx,pny,pnz,pnorm
      real fnx,fny,fnz,fnorm
      real xqc,yqc,zqc,roverr
      real sxxave,syyave,szzave
      real sxyave,sxzave,syzave
      real xave,yave,zave
      real cxs,cys,czs
      real rx,ry,rz
      real mx,my,mz
      real fac
!     added by camli 11/28/2013
      real fac1
      
!
      real dyinv
      real fmac,xpsi,srot,e1,e2,e3
      real cfx,cfy,cfz,cmpx
      real f1,f2
      real chorde,chordt,chordp,sigma,sigmat,sigmap
      real ctbysigma,cqbysigma

      ! Tip Mach scaling for Fortran coupling routine
      real fmac1

!     integers

      integer j1,j2,k1,k2,j2r,k2r,is
      integer l2,l3,kp,km1,jp,jm1
      integer j,k,l
      integer ii
!
!     first executable statement
!
      allocate(sxx(jmax,kmax),syy(jmax,kmax),szz(jmax,kmax))
      allocate(sxy(jmax,kmax),sxz(jmax,kmax),syz(jmax,kmax))
      allocate(vmu(jmax,kmax))
      allocate(cxsec(kmax),cmsec(kmax),czsec(kmax),cysec(kmax))

c..initialize constants and variables

      gami = gamma -1.
      j1 = jtail1
      j2 = jtail2
      k1 = kroot
      k2 = ktip
      j2r = j2 - 1
      k2r = k2 - 1
      is  = 1
      l   = lsurf
      l2  = l  + is
      l3  = l2 + is
      sw  = .5 * float( is )
      sf = 2.0 / 3.0

      cx  = 0.0
      cy  = 0.0
      cz  = 0.0
! added by camli badrya 11/15
      cp  = 0.0
! done   
c ... intilize mx,my,mz
c..... added by Camli B 9/2014
      mx  = 0.0
      my  = 0.0 
      mz  = 0.0 
       
      cost = cos(psi_rot)
      sint = sin(psi_rot)

      ca  = cos( pi * alf / 180.0 )
      sa  = sin( pi * alf / 180.0 )
cv
c      write(ifile,*) !k2r-k1+1
c     write(ifile,*) !totime,istep

c     write(ifile1,*) !k2r-k1+1
c     write(ifile1,*) !totime,istep

c     
c..evaluate area vector of the near & far triangles times two
c..then compute total area vector
c
      area = 0.0
      do 8 k=k1,k2r
        kp = k + 1
        if( k .eq. k2) kp = k1
        j = j1
        jp = jle
        xxi = x(jp,k,l) - x(j,k,l)
        yxi = y(jp,k,l) - y(j,k,l)
        xet = x(j,kp,l) - x(j,k,l)
        yet = y(j,kp,l) - y(j,k,l)
        az = xxi * yet - yxi * xet
        xxi = x(jp,kp,l) - x(j,kp,l)
        yxi = y(jp,kp,l) - y(j,kp,l)
        xet = x(jp,kp,l) - x(jp,k,l)
        yet = y(jp,kp,l) - y(jp,k,l)
        az = .5 * ( az + xxi * yet - yxi * xet )
        area = area + az
 8    continue
     
      if (fmtip.eq.0) then
         fac  = 2.0/((fsmach**2))
         fac1 = 2.0/((fsmach**3))  
      else
         fac  = 2.0/((fmtip**2))
!       added by camlib 11/28/20113
         fac1 = 2.0/((fmtip**3))
! 
      endif

c..set up the stress tensor s = - p d + t, where p is the pressure,
c..d is the kronecker delta tensor and t is the viscous stress tensor
c
c..pressure contribution
      c2b =199./tinf
      c2bp=c2b+1.
      gm1      = gamma - 1.
      gigm1    = gamma*gm1

      do 20 k=k1,k2
        do 10 j=j1,j2
          dj       = q(j,k,l,6)
          rhoji    = 1.0 / q(j,k,l,1)
          rho      = q(j,k,l,1) * dj
          u        = q(j,k,l,2) * rhoji
          v        = q(j,k,l,3) * rhoji
          w        = q(j,k,l,4) * rhoji
          e        = q(j,k,l,5) * dj
          qs       = .5 * rho * ( u*u + v*v + w*w )
          p        = gami * ( e - qs )
          ei       = (e-qs)/rho
          ttinf    = gigm1*ei
          vmu(j,k) = c2bp*(ttinf*sqrt(ttinf))/(c2b + ttinf)
          sxx(j,k) = - p + pinf
          syy(j,k) = - p + pinf
          szz(j,k) = - p + pinf
          sxy(j,k) = 0.0
          sxz(j,k) = 0.0
          syz(j,k) = 0.0
   10   continue
   20 continue
c
      if( .not. invisc) then

c..jaina
c..viscous contribution
c..this is assuming thin layer derivatives are only significant
c..have to be changed for full viscous terms

c.. trace of stress tensor tau_ii= mu*(2*du_i/dx_i-2/3*div(V))
c.. traceless part tau_ij = mu*(du_i/dx_j+du_j/dx_i)

        rei  = 1.0 / rey

        do 60 k=k1,k2
          kp = k + 1
          km1 = k - 1
          if (k.eq.k1) km1 = kp
          if (k.eq.k2) kp = km1
          do 50 j=j1,j2
            jp = j + 1
            jm1 = j - 1
            q1i1 = sw / q(j,k,l ,1)
            q1i2 = sw / q(j,k,l2,1)
            q1i3 = sw / q(j,k,l3,1)

c.. 2nd order in space

            usi  = 0.5*(q(jp,k,l,2)/q(jp,k,l,1)-q(jm1,k,l,2)/q(jm1,k,l,1))
            vsi  = 0.5*(q(jp,k,l,3)/q(jp,k,l,1)-q(jm1,k,l,3)/q(jm1,k,l,1))
            wsi  = 0.5*(q(jp,k,l,4)/q(jp,k,l,1)-q(jm1,k,l,4)/q(jm1,k,l,1))

            ueta = 0.5*(q(j,kp,l,2)/q(j,kp,l,1)-q(j,km1,l,2)/q(j,km1,l,1))
            veta = 0.5*(q(j,kp,l,3)/q(j,kp,l,1)-q(j,km1,l,3)/q(j,km1,l,1))
            weta = 0.5*(q(j,kp,l,4)/q(j,kp,l,1)-q(j,km1,l,4)/q(j,km1,l,1))

            uzet = - 3.0 * q1i1 * q(j,k,l,2) + 4.0 * q1i2 * q(j,k,l2,2)
     ,                                       -       q1i3 * q(j,k,l3,2)
            vzet = - 3.0 * q1i1 * q(j,k,l,3) + 4.0 * q1i2 * q(j,k,l2,3)
     ,                                       -       q1i3 * q(j,k,l3,3)
            wzet = - 3.0 * q1i1 * q(j,k,l,4) + 4.0 * q1i2 * q(j,k,l2,4)
     ,                                       -       q1i3 * q(j,k,l3,4)
c
            ux = xx(j,k,l)*usi + yx(j,k,l)*ueta + zx(j,k,l)*uzet
            uy = xy(j,k,l)*usi + yy(j,k,l)*ueta + zy(j,k,l)*uzet
            uz = xz(j,k,l)*usi + yz(j,k,l)*ueta + zz(j,k,l)*uzet

            vx = xx(j,k,l)*vsi + yx(j,k,l)*veta + zx(j,k,l)*vzet
            vy = xy(j,k,l)*vsi + yy(j,k,l)*veta + zy(j,k,l)*vzet
            vz = xz(j,k,l)*vsi + yz(j,k,l)*veta + zz(j,k,l)*vzet

            wx = xx(j,k,l)*wsi + yx(j,k,l)*weta + zx(j,k,l)*wzet
            wy = xy(j,k,l)*wsi + yy(j,k,l)*weta + zy(j,k,l)*wzet
            wz = xz(j,k,l)*wsi + yz(j,k,l)*weta + zz(j,k,l)*wzet

            sfdiv = sf*(ux + vy + wz)
c
            sxx(j,k) = sxx(j,k) + vmu(j,k) * rei * ( 2.0 * ux - sfdiv )
            syy(j,k) = syy(j,k) + vmu(j,k) * rei * ( 2.0 * vy - sfdiv )
            szz(j,k) = szz(j,k) + vmu(j,k) * rei * ( 2.0 * wz - sfdiv )
            sxy(j,k) = vmu(j,k) * rei * ( uy + vx )
            sxz(j,k) = vmu(j,k) * rei * ( uz + wx )
            syz(j,k) = vmu(j,k) * rei * ( vz + wy )
  210       format(2g13.5)
   50     continue
   60   continue
      endif
c
c..evaluate scaled force and moment coefficients
c..evaluate area vector of the near & far triangles times two,
c..then compute total area vector
c
      do 80 k=k1,k2r

        kp = k + 1

        cmsec(k)=0.0 
        cxsec(k)=0.0
        cysec(k)=0.0
        czsec(k)=0.0
        area=0.0

        xjle=(x(jle,k,l)+x(jle,kp,l))*0.5
        yjle=(y(jle,k,l)+y(jle,kp,l))*0.5
        zjle=(z(jle,k,l)+z(jle,kp,l))*0.5

        xjtp=(x(j2,k,l)+x(j2,kp,l))*0.5
        yjtp=(y(j2,k,l)+y(j2,kp,l))*0.5
        zjtp=(z(j2,k,l)+z(j2,kp,l))*0.5

        xptu=(x(j2-1,k,l)+x(j2-1,kp,l))*0.5
        yptu=(y(j2-1,k,l)+y(j2-1,kp,l))*0.5
        zptu=(z(j2-1,k,l)+z(j2-1,kp,l))*0.5

        xptl=(x(j1+1,k,l)+x(j1+1,kp,l))*0.5
        yptl=(y(j1+1,k,l)+y(j1+1,kp,l))*0.5
        zptl=(z(j1+1,k,l)+z(j1+1,kp,l))*0.5

c..     chord line vector

        cchord=sqrt((xjle-xjtp)**2+(zjle-zjtp)**2+(yjle-yjtp)**2)

        cvecx=(xjle-xjtp)/cchord
        cvecy=(yjle-yjtp)/cchord
        cvecz=(zjle-zjtp)/cchord

c..     normal to plane of airfoil, in the spanwise direction

        pnx=-(yptu-yjtp)*(zptl-zjtp)+(zptu-zjtp)*(yptl-yjtp)
        pny=-(zptu-zjtp)*(xptl-xjtp)+(xptu-xjtp)*(zptl-zjtp)
        pnz=-(xptu-xjtp)*(yptl-yjtp)+(yptu-yjtp)*(xptl-xjtp)
        
        pnorm=sqrt(pnx**2+pny**2+pnz**2)
        
        pnx=pnx/pnorm
        pny=pny/pnorm
        pnz=pnz/pnorm

c..     cross chord line vector with this to get normal force direction 
        
        fnx=-cvecy*pnz+cvecz*pny
        fny=-cvecz*pnx+cvecx*pnz
        fnz=-cvecx*pny+cvecy*pnx

        fnorm=sqrt(fnx**2+fny**2+fnz**2)
        
        fnx=fnx/fnorm
        fny=fny/fnorm
        fnz=fnz/fnorm

cvvc..     twist distribution 
cvv
cvv        dx=xjle-xjtp
cvv        dy=yjle-yjtp
cvv        dz=zjle-zjtp
cvv        
cvv        dnorm=sqrt(dx**2+dy**2+dz**2)
cvv        twist(k)=asin(dz/dnorm)
        
c..     quarter chord coordinates

        xqc=(xjle*0.75+xjtp*0.25)
        yqc=(yjle*0.75+yjtp*0.25)
        zqc=(zjle*0.75+zjtp*0.25)

        roverr=(yqc*cost-xqc*sint)/rartio

        do 70 j=j1,j2r

          jp = j + 1

c..   find the area vector

          xxi = x(jp,kp,l) - x(j,k,l)
          yxi = y(jp,kp,l) - y(j,k,l)
          zxi = z(jp,kp,l) - z(j,k,l)
          xet = x(j,kp,l) - x(jp,k,l)
          yet = y(j,kp,l) - y(jp,k,l)
          zet = z(j,kp,l) - z(jp,k,l)
          ax = 0.5*(yxi * zet - zxi * yet)
          ay = 0.5*(zxi * xet - xxi * zet)
          az = 0.5*(xxi * yet - yxi * xet)

!..original area calculation written by Jaina
!
!          xxi = x(jp,k,l) - x(j,k,l)
!          yxi = y(jp,k,l) - y(j,k,l)
!          zxi = z(jp,k,l) - z(j,k,l)
!          xet = x(j,kp,l) - x(j,k,l)
!          yet = y(j,kp,l) - y(j,k,l)
!          zet = z(j,kp,l) - z(j,k,l)
!          ax = yxi * zet - zxi * yet
!          ay = zxi * xet - xxi * zet
!          az = xxi * yet - yxi * xet
!          xxi = x(jp,kp,l) - x(j,kp,l)
!          yxi = y(jp,kp,l) - y(j,kp,l)
!          zxi = z(jp,kp,l) - z(j,kp,l)
!          xet = x(jp,kp,l) - x(jp,k,l)
!          yet = y(jp,kp,l) - y(jp,k,l)
!          zet = z(jp,kp,l) - z(jp,k,l)
!          ax = .5 * ( ax + yxi * zet - zxi * yet )
!          ay = .5 * ( ay + zxi * xet - xxi * zet )
!          az = .5 * ( az + xxi * yet - yxi * xet )
          
c
c..evaluate cell center values of the stress tensor and coordinates
c
          sxxave = .25 * ( sxx(j,k)+sxx(jp,k)+sxx(jp,kp)+sxx(j,kp) )
          syyave = .25 * ( syy(j,k)+syy(jp,k)+syy(jp,kp)+syy(j,kp) )
          szzave = .25 * ( szz(j,k)+szz(jp,k)+szz(jp,kp)+szz(j,kp) )
          sxyave = .25 * ( sxy(j,k)+sxy(jp,k)+sxy(jp,kp)+sxy(j,kp) )
          sxzave = .25 * ( sxz(j,k)+sxz(jp,k)+sxz(jp,kp)+sxz(j,kp) )
          syzave = .25 * ( syz(j,k)+syz(jp,k)+syz(jp,kp)+syz(j,kp) )
          xave = .25 * ( x( j,k,l)+x( jp,k,l)+x( jp,kp,l)+x( j,kp,l) )
          yave = .25 * ( y( j,k,l)+y( jp,k,l)+y( jp,kp,l)+y( j,kp,l) )
          zave = .25 * ( z( j,k,l)+z( jp,k,l)+z( jp,kp,l)+z( j,kp,l) )
c
          cxs = sxxave * ax + sxyave * ay + sxzave * az
          cys = sxyave * ax + syyave * ay + syzave * az
          czs = sxzave * ax + syzave * ay + szzave * az

          area=area+sqrt(ax*ax+ay*ay+az*az)
c

         if (fmtip.eq.0) then
            if (k.ge.send_lims(KDIR).and.k.le.send_lims(KDIR+3)-1) then
              cx  = cx + fac*cxs
              cy  = cy + fac*cys
              cz  = cz + fac*czs

!  added by camli 8/20/2013 Calculating power
              rhoij2 = 1.0/q(j,k,l,1)
              velx = q(j,k,l,2)*rhoij2
              vely = q(j,k,l,3)*rhoij2
              velz = q(j,k,l,4)*rhoij2
              cp = cp +fac1*(cxs*velx + cys*vely + czs*velz)

            endif
          else
            if (k.ge.send_lims(KDIR).and.k.le.send_lims(KDIR+3)-1) then
              cx  = cx + fac*cxs!*roverr
              cy  = cy + fac*cys!*roverr
              cz  = cz + fac*czs

!  added by camli 11/15/2013 Calculating power
              rhoij2 = 1.0/q(j,k,l,1)
              velx = q(j,k,l,2)*rhoij2
              vely = q(j,k,l,3)*rhoij2
              velz = q(j,k,l,4)*rhoij2
              cp = cp +fac1*(cxs*velx + cys*vely + czs*velz)
! done calculating power coeff.     
 
            endif
          endif


c .... Moments calculate...........................................................
c .... Sep. 2014, Camli.B ..........................................................
c .... M=rxf vector calculation ....................................................
c .... r=[rx,ry,rz]... f=[cxs,cys,czs]... mx,my,mz -total moments in x,y,z directions

            rx = xave-xcg
            ry = yave-ycg
            rz = zave-zcg

            mx = mx + fac*(ry*czs-rz*cys)
            my = my + fac*(rz*cxs-rx*czs)
            mz = mz + fac*(rx*cys-ry*cxs)

!            write(*,*) mx,my,mz
c .....done calculating moments ....................................................





cvin          cxsec(k) = cxsec(k) + (cxs*pnx+cys*pny+czs*pnz)
cvin          cysec(k) = cysec(k) + (cxs*cvecx+cys*cvecy+czs*cvecz)
cvin          czsec(k) = czsec(k) + (cxs*fnx+cys*fny+czs*fnz)
cvin
cvin          rx=xave-xqc
cvin          ry=yave-yqc
cvin          rz=zave-zqc
cvin
cvin          mx=ry*czs-rz*cys
cvin          my=rz*cxs-rx*czs
cvin          mz=rx*cys-ry*cxs
cvin
cvin          cmsec(k)=cmsec(k)+ (pnx*mx+pny*my+pnz*mz)

   70   continue
        
cvin        dyinv=2.0/(0.5*area*rinf)
cvin
cvin        cxsec(k)=cxsec(k)*dyinv
cvin        cysec(k)=cysec(k)*dyinv
cvin        czsec(k)=czsec(k)*dyinv
cvin        cmsec(k)=cmsec(k)*dyinv
cvin
cvinc..   these are in true deformed frame
cvin       if (fmtip.ne.0) then
cvin
cvin         !Use fmac1 for legacy fortran coupling (Anubhav-Jaina setup)
cvinc         fmac1=(fmm/fmtip)**2
cvin         !Use fmac for python-coupling (Jaina's setup)
cvinc         fmac=fmm**2
cvinc         write(ifile,2001) roverr,cxsec(k)*fmac1,cysec(k)*fmac1,czsec(k)*fmac1,cmsec(k)*fmac1
cvin       else
cvinc         write(ifile,2001) roverr,cxsec(k),cysec(k),czsec(k),cmsec(k)
cvin         fmac=1
cvin       endif
cvin
cvinc..   lets find it in true hub frame also now
cvin
cvin        xpsi=psi_rot
cvin        srot=rf*dt
cvin
cvin        if (idefl.eq.1.or.idefl.eq.2) then
cvin          e1=twist(k)
cvin          e2=0.
cvin          e3=0.
cvin        else
cvin         e1=0.
cvin         e2=0.
cvin         e3=twist(k)
cvin        endif
cvin        
cvinc        ipyloads = ldidx ! Default value; may be overwritten below
cvinc        if (ideform.eq.1) then
cvinc           if (srot.ne.0.0) then
cvinc
cvinc              if(arf_opt .eq. 0)then
cvinc                 ii=mod(nint(xpsi/srot),iazimuth)+1
cvinc              else
cvinc                 ii = 2
cvinc                 ipyloads = 1
cvinc              end if
cvinc                 
cvinc              e1=e1+defl_dat(4,k,ii)
cvinc              e2=e2+defl_dat(5,k,ii)
cvinc              e3=e3+defl_dat(6,k,ii)
cvinc
cvinc           endif
cvinc        endif
cvin
cvin        if (idefl.eq.1) then
cvin           call rmat_extract(e1,e2,e3,rt)
cvin        elseif (idefl.eq.2) then
cvin           call rmat_extract_rcas(e1,e2,e3,rt)
cvin        else
cvin           call rmat_extract_tdu(e1,e2,e3,rt)
cvin        endif
cvin        
cvin        cfx=rt(1,1)*cxsec(k)+rt(1,2)*cysec(k)+rt(1,3)*czsec(k)
cvin        cfy=rt(2,1)*cxsec(k)+rt(2,2)*cysec(k)+rt(2,3)*czsec(k)
cvin        cfz=rt(3,1)*cxsec(k)+rt(3,2)*cysec(k)+rt(3,3)*czsec(k)
cvin        cmpx=rt(1,1)*cmsec(k)
cvin
cvinc        ! for python
cvinc        if (pythonMode.eq.1) then
cvinc           loads(1,k-k1+1,ipyloads)=cfx*fmac
cvinc           loads(2,k-k1+1,ipyloads)=cfy*fmac
cvinc           loads(3,k-k1+1,ipyloads)=cfz*fmac
cvinc           loads(4,k-k1+1,ipyloads)=cmpx*fmac
cvinc	   loads(5,k-k1+1,ipyloads)=cd2dt2*fmac
cvinc	   loads(6,k-k1+1,ipyloads)=cl2dt2*fmac
cvinc	   loads(7,k-k1+1,ipyloads)=cmqc2*fmac
cvinc        endif
cvinc        ! end for python
cvin
cvinc        write(ifile1,2001) roverr,cfx*fmac,cfy*fmac,cfz*fmac,cmpx*fmac
cvin
   80 continue
c.... added by Camli.B adding ifile2 that plot the moments mx,my,mz

      write(ifile1,2001) cx,cy,cz,cp
      write(ifile2,2002) mx,my,mz

      if (fmtip.ne.0) then
         ct=float(nblade_tot)*cz/(2*pi*rartio**2)
         cq=float(nblade_tot)*(cx*cost+cy*sint)/(2*pi*rartio**2)
         sigma=float(nblade_tot)/pi/rartio
         ctbysigma=ct/sigma
         cqbysigma=cq/sigma
         figofmerit=sqrt(abs(ct)**3)/sqrt(2.0)/cq
      else
         f1 = cz/rartio
         f2 = (cx*cost+cy*sint)/rartio
         ct = f1*ca - f2*sa
         cq = f1*sa + f2*ca
         figofmerit = ct/cq
      endif

      deallocate(sxx,syy,szz)
      deallocate(sxy,sxz,syz)
      deallocate(cxsec,cysec,czsec,cmsec)

!2001   format(5(E12.5,x))
2001   format(4(E12.5,x))
2002   format(3(E12.5,x))

      return
      end

