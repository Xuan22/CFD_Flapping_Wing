c**********************************************************************
      subroutine force2d(lsurf,x0,y0,z0,
     &     ct,cq,figofmerit,x,y,z,q,x1x,x1y,x1z,z1x,z1y,z1z,twist,
     $     ldidx,ifile,ifile1)

c..   compute forces and moments both in deformed and undeformed frames
c**********************************************************************

      use params_global
      use deflections
      use pyMOD
      use arf_mod
      use rootVariables, only: psi_rot
      use domain_info, only: send_lims

      implicit none

      integer lsurf
      real x0,y0,z0
      
      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)
      real q(jmax,kmax,lmax,nd)
      real x1x(jmax,kmax,lmax),x1y(jmax,kmax,lmax),x1z(jmax,kmax,lmax)
      real z1x(jmax,kmax,lmax),z1y(jmax,kmax,lmax),z1z(jmax,kmax,lmax)
      real twist(kmax)

      real ct,cq,figofmerit
      integer ifile, ldidx, ifile1

c.. local variables
      
      real,allocatable :: rt(:,:)
      real,allocatable :: zz(:)
      real,allocatable :: df1(:),df2(:),df3(:),df4(:),df5(:),df6(:)
      real,allocatable :: df7(:),df8(:),r(:)
		
      real cl2d,cd2d,cm2d,clvv2d,cdvv2d,cmvv2d,cmqc,cmqcv
      real ca,sa,cost,sint,roverr,fmm,roverrr,rr,rho,u,v,w,e,vsq,pp,pp1
      real cp,cn,cc,cmle,cnv,ccv,cmlev,cn1,cc1,cmle1,cnv1,ccv1,cmlev1
      real xjle,yjle,zjle,xjtp,yjtp,zjtp,xpt,zpt,ypt,cchord,chord
      real cvecx,cvecy,cvecz,pnx,pny,pnz,pnorm,fnx,fny,fnz,fnorm
      real xqc,yqc,zqc
      real cfav,cpav,xj,yj,zj,xjm1,yjm1,zjm1
      real cnx,cny,cnz
      real rx,ry,rz,cmx,cmy,cmz
      real alngth,amu,uinf2,uxi,vxi,wxi,ueta,veta,weta
      real xix,xiy,xiz,etax,etay,etaz,tauw
      real fmac,xpsi,srot,e1,e2,e3,cfx,cfy,cfz,cmpx,fmtt
      real cl2dt,cd2dt,cl2dt2,cd2dt2,cmqc2
      real ct2d,cq2d,gamr
      real f1,f2,f3,f4,f5,f6,f7,f8
      real chorde,chordt,chordp,sigma,sigmat,sigmap
      real ctbysigma,cqbysigma

      ! Tip Mach scaling for Fortran coupling routine
      real fmac1

      integer ii,j,k,l,l2,l3,j1,j2,jtp,k1,k2,ibin,jm1,ja,jb,jp1
      integer ipyloads ! Added 10/25/07 BS

      allocate(rt(3,3))
      allocate(zz(jmax))
      allocate(df1(kmax),df2(kmax),df3(kmax),df4(kmax),df5(kmax),df6(kmax))
      allocate(df7(kmax),df8(kmax),r(kmax))

c***  first executable statement

c.. initialize constants and variables

      if(istep.eq.0)then
c         write(ifile,*)jtail1,jle,jtail2,kroot,ktip
c         write(ifile1,*)jtail1,jle,jtail2,kroot,ktip
      endif

      pi  = 4.0 * atan( 1.0 )
      ca  = cos( pi * alf / 180.0 )
      sa  = sin( pi * alf / 180.0 )

      cl2d = 0.
      cd2d = 0.
      cm2d = 0.
      clvv2d = 0.
      cdvv2d = 0.
      cmvv2d = 0.

c..set limits of integration
     
      l = lsurf
      l2 = l+1
      l3 = l2+1
      j1 = jtail1
      j2 = jtail2
      jtp = jtail1+1
      k1 = kroot
      k2 = ktip
      ibin = 1
      
c      write(ifile,*) !k2-k1+1
c      write(ifile,*) !totime,istep
      
c      write(ifile1,*) !k2-k1+1
c      write(ifile1,*) !totime,istep

c..compute cp at grid points and store in an array 

      cost = cos(psi_rot)
      sint = sin(psi_rot)
      do 1000 k=k1,k2

        roverr=(y(jle,k,l)*cost-x(jle,k,l)*sint)/rartio

c..   for a fixed wing simulation

        if (fmtip.eq.0) then 
           fmm=fsmach
        else
!******* ORIGINAL: LOOKS LIKE A BUG? ************************
!           fmm = sqrt( (fmtip*roverr+vinf*sin(rf*totime))**2
!     <          + winf**2 + (vinf*cos(rf*totime))**2 )
!************************************************************
            fmm = sqrt( (fmtip*roverr + vinf*sint + uinf*cost)**2
     <           + winf**2 + (vinf*cost - uinf*sint)**2)
        endif

          do j=j1,j2
            rr=1./q(j,k,l,1)
            rho=q(j,k,l,1)*q(j,k,l,6)
            u=q(j,k,l,2)*rr
            v=q(j,k,l,3)*rr
            w=q(j,k,l,4)*rr
            e=q(j,k,l,5)*q(j,k,l,6)

            vsq=u*u+v*v+w*w
            pp=gm1*(e-rho*vsq/2.)
            pp1=pp/pinf
            cp=2.*(pp1-1.)/(gamma*fmm**2)
            zz(j) = cp
          enddo

c..compute normal force coefficient and chord directed force coeff
c..chord taken as one in all cases, modified later

        cn = 0.
        cc = 0.
        cmle = 0.

        cn1=0.
        cc1=0.
        cmle1=0.

c jaina
c for xjle need not be zero for blades with elastic lag or sweep back
c like uh-60
c.. modifying for true deformed frame of reference
c.. 05/19/03

        xjle=sint*y(jle,k,l)+cost*x(jle,k,l)
        yjle=cost*y(jle,k,l)-sint*x(jle,k,l)
        zjle=z(jle,k,l)

        xjtp=sint*y(j2,k,l)+cost*x(j2,k,l)
        yjtp=cost*y(j2,k,l)-sint*x(j2,k,l)
        zjtp=z(j2,k,l)
        
        xpt=sint*y(j2-1,k,l)+cost*x(j2-1,k,l)
        ypt=cost*y(j2-1,k,l)-sint*x(j2-1,k,l)
        zpt=z(j2-1,k,l)
        
        cchord=sqrt((xjle-xjtp)**2+(zjle-zjtp)**2+(yjle-yjtp)**2)

c..     chord line vector
        
        cvecx=(xjle-xjtp)/cchord
        cvecy=(yjle-yjtp)/cchord
        cvecz=(zjle-zjtp)/cchord

c..     normal to plane of airfoil, in the spanwise direction

        pnx=-(ypt-yjtp)*cvecz+(zpt-zjtp)*cvecy
        pny=-(zpt-zjtp)*cvecx+(xpt-xjtp)*cvecz
        pnz=-(xpt-xjtp)*cvecy+(ypt-yjtp)*cvecx
        
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


c..     quarter chord coordinates

        xqc=(xjle*0.75+xjtp*0.25)
        yqc=(yjle*0.75+yjtp*0.25)
        zqc=(zjle*0.75+zjtp*0.25)


        do 11 j=jtp,j2

          jm1 = j-1
          cpav = (zz(j) + zz(jm1))*.5

          xj = sint*y(j,k,l)+cost*x(j,k,l)
          yj = cost*y(j,k,l)-sint*x(j,k,l)
          zj = z(j,k,l)

          xjm1 = sint*y(jm1,k,l)+cost*x(jm1,k,l)
          yjm1 = cost*y(jm1,k,l)-sint*x(jm1,k,l)
          zjm1 = z(jm1,k,l)

          cnx=(zj-zjm1)*pny-(yj-yjm1)*pnz
          cny=(xj-xjm1)*pnz-(zj-zjm1)*pnx
          cnz=(yj-yjm1)*pnx-(xj-xjm1)*pny

c..   in the deformed frame

          cn1=cn1+cpav*(cnx*fnx+cny*fny+cnz*fnz)
          cc1=cc1-cpav*(cnx*cvecx+cny*cvecy+cnz*cvecz)

c..   normal force and chord force in the classical 2-d sense

          cn = cn - cpav*( xj - xjm1 )
          cc = cc + cpav*( z(j,k,l) - z(jm1,k,l))

          cmle = cmle + cpav*
     &       (  (xj+xjm1-2*xjle)*.5*(xj-xjm1) +
     &          ( z(j,k,l)+z(jm1,k,l)-2*zjle)*.5*(z(j,k,l)-z(jm1,k,l)) )

          rx=(xj+xjm1)*0.5-xqc
          ry=(yj+yjm1)*0.5-yqc
          rz=(zj+zjm1)*0.5-zqc

          cmx=cnz*ry-rz*cny
          cmy=cnx*rz-rx*cnz
          cmz=cny*rx-ry*cnx

c..   pitching moment in the direction of true deformed axis

          cmle1=cmle1+cpav*(cmx*pnx+cmy*pny+cmz*pnz)

   11   continue

        cn = cn - half*cn
        cmle = cmle - half*cmle
        cc = cc + half*cc
        cl2d = cn
        cd2d = cc
        cmqc = cmle +.25*cn
        cm2d = cmqc

c..   viscous coefficent of friction calculation
c..   re already has fsmach scaling

        cnv = 0.
        ccv = 0.
        cmlev = 0.

        cnv1=0.
        ccv1=0.
        cmlev1=0.

        if(.not.invisc) then

          alngth= 1.0
          amu   = rinf*alngth/rey
          uinf2 = fmm**2
          !uinf2 = 1.0
          j = j1
          jp1 = j+1
          uxi = q(jp1,k,l,2)/q(jp1,k,l,1)-q(j,k,l,2)/q(j,k,l,1)
          vxi = q(jp1,k,l,3)/q(jp1,k,l,1)-q(j,k,l,3)/q(j,k,l,1)
          wxi = q(jp1,k,l,4)/q(jp1,k,l,1)-q(j,k,l,4)/q(j,k,l,1)
          ueta= -1.5*q(j,k,l,2)/q(j,k,l,1)+2.*q(j,k,l2,2)/q(j,k,l2,1)
     *         -.5*q(j,k,l3,2)/q(j,k,l3,1)
          veta= -1.5*q(j,k,l,3)/q(j,k,l,1)+2.*q(j,k,l2,3)/q(j,k,l2,1)
     *         -.5*q(j,k,l3,3)/q(j,k,l3,1)
          weta= -1.5*q(j,k,l,4)/q(j,k,l,1)+2.*q(j,k,l2,4)/q(j,k,l2,1)
     *         -.5*q(j,k,l3,4)/q(j,k,l3,1)
          xix = x1x(j,k,l)
          xiy = x1y(j,k,l)
          xiz = x1z(j,k,l)
          etax = z1x(j,k,l)
          etay = z1y(j,k,l)
          etaz = z1z(j,k,l)
          tauw= amu*((uxi*xiz+ueta*etaz)-(wxi*xix+weta*etax))*cost +
     *          amu*((vxi*xiz+veta*etaz)-(wxi*xiy+weta*etay))*sint
          zz(j)= tauw/(.5*rinf*uinf2)

          j = j2
          jm1 = j-1
          uxi = q(j,k,l,2)/q(j,k,l,1)-q(jm1,k,l,2)/q(jm1,k,l,1)
          vxi = q(j,k,l,3)/q(j,k,l,1)-q(jm1,k,l,3)/q(jm1,k,l,1)
          wxi = q(j,k,l,4)/q(j,k,l,1)-q(jm1,k,l,4)/q(jm1,k,l,1)
          ueta= -1.5*q(j,k,l,2)/q(j,k,l,1)+2.*q(j,k,l2,2)/q(j,k,l2,1)
     *         -.5*q(j,k,l3,2)/q(j,k,l3,1)
          veta= -1.5*q(j,k,l,3)/q(j,k,l,1)+2.*q(j,k,l2,3)/q(j,k,l2,1)
     *         -.5*q(j,k,l3,3)/q(j,k,l3,1)
          weta= -1.5*q(j,k,l,4)/q(j,k,l,1)+2.*q(j,k,l2,4)/q(j,k,l2,1)
     *         -.5*q(j,k,l3,4)/q(j,k,l3,1)
          xix = x1x(j,k,l)
          xiy = x1y(j,k,l)
          xiz = x1z(j,k,l)
          etax = z1x(j,k,l)
          etay = z1y(j,k,l)
          etaz = z1z(j,k,l)
          tauw= amu*((uxi*xiz+ueta*etaz)-(wxi*xix+weta*etax))*cost +
     *          amu*((vxi*xiz+veta*etaz)-(wxi*xiy+weta*etay))*sint
          zz(j)= tauw/(.5*rinf*uinf2)

c..set limits

          ja = j1+1
          jb = j2-1
          do 110 j = ja,jb
            jp1 = j+1
            jm1 = j-1
            uxi=.5*(q(jp1,k,l,2)/q(jp1,k,l,1)-q(jm1,k,l,2)/q(jm1,k,l,1))
            vxi=.5*(q(jp1,k,l,3)/q(jp1,k,l,1)-q(jm1,k,l,3)/q(jm1,k,l,1))
            wxi=.5*(q(jp1,k,l,4)/q(jp1,k,l,1)-q(jm1,k,l,4)/q(jm1,k,l,1))
            ueta= -1.5*q(j,k,l,2)/q(j,k,l,1)+2.*q(j,k,l2,2)/q(j,k,l2,1)
     *                 -.5*q(j,k,l3,2)/q(j,k,l3,1)
            veta= -1.5*q(j,k,l,3)/q(j,k,l,1)+2.*q(j,k,l2,3)/q(j,k,l2,1)
     *                 -.5*q(j,k,l3,3)/q(j,k,l3,1)
            weta= -1.5*q(j,k,l,4)/q(j,k,l,1)+2.*q(j,k,l2,4)/q(j,k,l2,1)
     *                 -.5*q(j,k,l3,4)/q(j,k,l3,1)
            xix = x1x(j,k,l)
            xiy = x1y(j,k,l)
            xiz = x1z(j,k,l)
            etax = z1x(j,k,l)
            etay = z1y(j,k,l)
            etaz = z1z(j,k,l)
            tauw= amu*((uxi*xiz+ueta*etaz)-(wxi*xix+weta*etax))*cost +
     *            amu*((vxi*xiz+veta*etaz)-(wxi*xiy+weta*etay))*sint
            zz(j)= tauw/(.5*rinf*uinf2)

110       continue

c..compute viscous normal and axial forces

          do 111 j=jtp,j2
            jm1 = j-1
            cfav = ( zz(j) + zz(jm1))*.5

            xj = sint*y(j,k,l)+cost*x(j,k,l)
            yj = cost*y(j,k,l)-sint*x(j,k,l)
            zj = z(j,k,l)

            xjm1 = sint*y(jm1,k,l)+cost*x(jm1,k,l)
            yjm1 = cost*y(jm1,k,l)-sint*x(jm1,k,l)
            zjm1 = z(jm1,k,l)

            cnx=xj-xjm1
            cny=yj-yjm1
            cnz=zj-zjm1

c..   in true deformed frame

            ccv1=ccv1-(cnx*cvecx+cny*cvecy+cnz*cvecz)*cfav
            cnv1=cnv1+(cnx*fnx+cny*fny+cnz*fnz)*cfav

            rx=(xj+xjm1)*0.5-xqc
            ry=(yj+yjm1)*0.5-yqc
            rz=(zj+zjm1)*0.5-zqc
            
            cmx=cnz*ry-rz*cny
            cmy=cnx*rz-rx*cnz
            cmz=cny*rx-ry*cnx
            
            cmlev1=cmlev1-cfav*(cmx*pnx+cmy*pny+cmz*pnz)

c..   in classical 2-d sense

            ccv = ccv + cfav*(( xj-xjm1))
            cnv = cnv + cfav*( z(j,k,l) - z(jm1,k,l))

            cmlev = cmlev - cfav*
     *       (  (xj+xjm1-2*xjle)*.5*(z(j,k,l) -z(jm1,k,l)) -
     *          (z(j,k,l)+z(jm1,k,l)-2*zjle)*.5*(xj -xjm1)   )

111       continue

          cnv = cnv - half*cnv
          cmlev = cmlev - half*cmlev
          ccv = ccv + half*ccv
          clvv2d = cnv
          cdvv2d = ccv
          cmqcv = cmlev +.25*cnv
          cmvv2d = cmqcv

        endif

c..modified for chord being non-unity

        chord = sqrt(((x(jtail1,k,l)-x(jle,k,l))**2
     &        +      (y(jtail1,k,l)-y(jle,k,l))**2
     &        +      (z(jtail1,k,l)-z(jle,k,l))**2))

        cl2d=cl2d/chord
        clvv2d=clvv2d/chord
        cd2d=cd2d/chord
        cdvv2d=cdvv2d/chord

        cl2dt = cl2d+clvv2d
        cd2dt = cd2d+cdvv2d
        
        cl2dt2=(cn1+cnv1)/chord
        cd2dt2=(cc1+ccv1)/chord
        cmqc2=(cmle1+cmlev1)/chord

c..   if 2-d type forces required comment these

c        write(ifile,*) roverr,cl2dt,cd2dt,
c     &       (cmqc+cmqcv)/chord

c..   these are in true deformed frame
       if (fmtip.ne.0) then

         !Use fmac1 for legacy fortran coupling (Anubhav-Jaina setup)
         fmac1=(fmm/fmtip)**2
         !Use fmac for python-coupling (Jaina's setup)
         fmac=fmm**2
         write(ifile,'(4(E15.6,x))') roverr,cl2dt2*fmac1,cd2dt2*fmac1,cmqc2*fmac1
       else
         write(ifile,*) roverr,cl2dt2,cd2dt2,cmqc2
         fmac=1
       endif

c..   lets find it in true hub frame also now

       xpsi=psi_rot
       srot=rf*dt

       if (idefl.eq.1.or.idefl.eq.2) then
          e1=twist(k)
          e2=0.
          e3=0.
        else
         e1=0.
         e2=0.
         e3=twist(k)
	 !e3=0.
        endif
        
        ipyloads = ldidx ! Default value; may be overwritten below
        if (ideform.eq.1) then
           if (srot.ne.0.0) then

              if(arf_opt .eq. 0)then
                 ii=mod(nint(xpsi/srot),iazimuth)+1
              else
                 ii = 2
                 ipyloads = 1
              end if
                 
              e1=e1+defl_dat(4,k,ii)
              e2=e2+defl_dat(5,k,ii)
              e3=e3+defl_dat(6,k,ii)

           endif
        endif
        
        if (idefl.eq.1) then
           call rmat_extract(e1,e2,e3,rt)
        elseif (idefl.eq.2) then
           call rmat_extract_rcas(e1,e2,e3,rt)
        else
           call rmat_extract_tdu(e1,e2,e3,rt)
        endif
        
        !NOTE-TO-SELF: Here, the x-coor in the TURNS NR-frame
        ! is switched to 
        cfx=-rt(1,2)*cd2dt2+rt(1,3)*cl2dt2
        cfy=-rt(2,2)*cd2dt2+rt(2,3)*cl2dt2
        cfz=-rt(3,2)*cd2dt2+rt(3,3)*cl2dt2
        cmpx=cmqc2*rt(1,1)

        ! for python
        if (pythonMode.eq.1) then
           loads(1,k-k1+1,ipyloads)=cfx*fmac
           loads(2,k-k1+1,ipyloads)=cfy*fmac
           loads(3,k-k1+1,ipyloads)=cfz*fmac
           loads(4,k-k1+1,ipyloads)=cmpx*fmac
	   loads(5,k-k1+1,ipyloads)=cd2dt2*fmac
	   loads(6,k-k1+1,ipyloads)=cl2dt2*fmac
	   loads(7,k-k1+1,ipyloads)=cmqc2*fmac
        endif
        ! end for python

        !write(ifile1,2001) roverr,cfx*fmac,cfy*fmac,cfz*fmac,cmpx*fmac
c        write(ifile1,2001) roverr,cfx,cfy,cfz,cmpx
 2001   format(5(E12.5,x))

        if (fmtip.ne.0) then

           fmtt=fmtip
           df1(k) = cfz*chord*(fmm/fmtt)**2
           df2(k) = -cfy*chord*(fmm/fmtt)**3
           df3(k) = chord
           df4(k) = 1.
           df5(k) = chord*(fmm/fmtt)**2
           df6(k) = (fmm/fmtt)**2
           df7(k) = chord*(fmm/fmtt)**3
           df8(k) = (fmm/fmtt)**3
           r(k) = rartio*roverr
           ct2d = cfz*(fmm/fmtt)**2
           cq2d = -cfy*(fmm/fmtt)**3
           gamr = chord*cfz*r(k)/2./rartio**2

        else
           df1(k) = cfz*chord
           df2(k) = -cfy*chord
           r(k) = roverr
        endif

 1000 continue

c..trapezoidal integration
c..if computing rotor flows

      f1 = 0.0
      f2 = 0.0
      f3 = 0.0
      f4 = 0.0
      f5 = 0.0
      f6 = 0.0
      f7 = 0.0
      f8 = 0.0

      if (fmtip.ne.0) then
         do 30 k=send_lims(KDIR),send_lims(KDIR+3)-1
            f1=f1+0.5*(df1(k)+df1(k+1))*(r(k+1)-r(k))
            f2=f2+0.5*(df2(k)+df2(k+1))*(r(k+1)-r(k))
            f3=f3+0.5*(df3(k)+df3(k+1))*(r(k+1)-r(k))
            f4=f4+0.5*(df4(k)+df4(k+1))*(r(k+1)-r(k))
            f5=f5+0.5*(df5(k)+df5(k+1))*(r(k+1)-r(k))
            f6=f6+0.5*(df6(k)+df6(k+1))*(r(k+1)-r(k))
            f7=f7+0.5*(df7(k)+df7(k+1))*(r(k+1)-r(k))
            f8=f8+0.5*(df8(k)+df8(k+1))*(r(k+1)-r(k))
 30      continue
     
         ct=f1*float(nblade_tot)/pi/rartio**2/2.0
         cq=f2*float(nblade_tot)/pi/rartio**2/2.0
         chorde=f3/f4
         chordt=f5/f6
         chordp=f7/f8
         sigma=float(nblade_tot)*chorde/pi/rartio
         sigmat=float(nblade_tot)*chordt/pi/rartio
         sigmap=float(nblade_tot)*chordp/pi/rartio

         ctbysigma=ct/sigma
         cqbysigma=cq/sigma
         figofmerit=sqrt(abs(ct)**3)/sqrt(2.0)/cq
      else
         do k=send_lims(KDIR),send_lims(KDIR+3)-1
            f1=f1+0.5*(df1(k)+df1(k+1))*(r(k+1)-r(k))
            f2=f2+0.5*(df2(k)+df2(k+1))*(r(k+1)-r(k))
         enddo
         ct = f1*ca - f2*sa
         cq = f1*sa + f2*ca
         figofmerit = ct/cq
      endif

      return
      end
