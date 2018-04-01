c***********************************************************************
      subroutine pert3(q,x,y,z,ug,vg,wg,zt0,zx,zy,zz,psi)

c**** Prologue : *****
c
c     subroutine which controls the wake inclusion using the
c     field velocity approach. 
c
c     imth= 1 (Kitaplioglu and Caradonna rotor: BVI at 180 deg)
c     imth= 2 (brute force induced velocity calculation using freewake)
c     imth= 3 (fast induced velocity calculation)
c     
c     last updated by jaina 07/13/04
c
c**** end prologue ****************************************************

      use mpi_wrapper
      use params_global
      use refmesh
      use freewake

c*******************  list of global variables used *******************
c
c     jmax,kmax,lmax,nd,totime,rf,imth,irefine  -> (params_global)
c     nw -> (freewake)
c     ug1,vg1,wg1 -> (refmesh)
c
c**********************************************************************

      implicit none

      real q(jmax,kmax,lmax,nd)
      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax),vg(jmax,kmax,lmax),wg(jmax,kmax,lmax)
      real zt0(jmax,kmax)
      real zx(jmax,kmax,lmax),zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)
      real psi

c..   local variables

      integer j,k,l,j21,k21,l21,w,npts
      real rc,xa,za,xv,zv,ra
      real, dimension(:,:,:), allocatable :: uge,vge,wge
      
c****first executable statement

      allocate(uge(jmax,kmax,lmax),vge(jmax,kmax,lmax),
     $     wge(jmax,kmax,lmax))
      npts=jmax*kmax*lmax

c..   setting stuff to zero

      do l=1,lmax
         do k=1,kmax
            do j=1,jmax
               uge(j,k,l)=0.
               vge(j,k,l)=0.
               wge(j,k,l)=0.
            enddo
         enddo
      enddo

      ! if using refined mesh

      if (irefine.eq.1) then
         do l=1,2*lmax-1
            do k=1,2*kmax-1
               do j=1,2*jmax-1
                  ugp(j,k,l)=ug1(j,k,l)
                  vgp(j,k,l)=vg1(j,k,l)
                  wgp(j,k,l)=wg1(j,k,l)
               enddo
            enddo
         enddo
      endif
      
      if (imth.eq.1) then  ! kit and caradonna test

       ! See Caradonna and Tung test report
       ! NASA/TM-1999-208790
 
       vorgam=3*0.374*fsmach*fsmach/fmtip/(2*pi)

        rc=0.162
        zv=0.25
        xv=0.0

        do l=1,lmax
          do k=1,kmax
            do j=1,jmax
             xa=x(j,k,l)-xv
             za=z(j,k,l)-zv
             ra=xa**2+za**2+rc**2
             uge(j,k,l)=-vorgam*za/ra
             vge(j,k,l)=0.
             wge(j,k,l)=vorgam*xa/ra
            enddo
          enddo
        enddo

      else if (imth.eq.FV_BIOTSAVART) then  ! brute force induced velocity evaluation
         do w=1,nw
            call fwind(w,x,y,z,uge,vge,wge,psi)
         enddo
      elseif (imth.eq.FV_FAST) then ! fast velocity evaluation
         do w=1,nw
!            call fwind_fast(w,x,y,z,uge,vge,wge,psi)
            call fwind_fast3(w,x,y,z,uge,vge,wge,npts,psi)
         enddo
      else if (imth.eq.FV_USER) then
	call user_specified_field(x,y,z,uge,vge,wge,npts)
      endif

c..   if using refined mesh for metrics and rhsup

      if (irefine.eq.1) then

c..   interpolating to find the grid velocities in refined mesh
c..   1-d interpolation is only used

         do l=1,lmax
            l21=2*l-1
            do k=1,kmax
               k21=2*k-1
               do j=1,jmax
                  ug1(2*j-1,k21,l21)=uge(j,k,l)
                  vg1(2*j-1,k21,l21)=vge(j,k,l)
                  wg1(2*j-1,k21,l21)=wge(j,k,l)
               enddo
            enddo
         enddo
      
         do k=1,kmax
            k21=2*k-1
            
            do l=1,lmax
               l21=2*l-1
               do j=1,jmax-1
                  ug1(2*j,k21,l21)=(uge(j,k,l)+uge(j+1,k,l))*0.5
                  vg1(2*j,k21,l21)=(vge(j,k,l)+vge(j+1,k,l))*0.5
                  wg1(2*j,k21,l21)=(wge(j,k,l)+wge(j+1,k,l))*0.5
               enddo
            enddo
            
            do j=1,2*jmax-1
               do l=1,lmax-1
                  ug1(j,k21,2*l)=(ug1(j,k21,2*l-1)+ug1(j,k21,2*l+1))*0.5
                  vg1(j,k21,2*l)=(vg1(j,k21,2*l-1)+vg1(j,k21,2*l+1))*0.5
                  wg1(j,k21,2*l)=(wg1(j,k21,2*l-1)+wg1(j,k21,2*l+1))*0.5
               enddo
            enddo
         enddo

      endif
         
c..   accumulating the total grid time metrics in ug,vg,wg
c..   to be used in the lhs and rhs (if not using refined meshes)

      do j=1,jmax
         do k=1,kmax
            do l=1,lmax
               ug(j,k,l)=ug(j,k,l)+uge(j,k,l)
               vg(j,k,l)=vg(j,k,l)+vge(j,k,l)
               wg(j,k,l)=wg(j,k,l)+wge(j,k,l)
            enddo
         enddo
      enddo

      deallocate(uge,vge,wge)
      return
      end


c************************************************************************
      subroutine fwind(w,x,y,z,ug,vg,wg,psi)

c**** Prologue : *****
c
c     Brute force calculation of induced velocities
c     use the wake geometry to caluclate field velocities at all
c     grid points                                                    
c
c     last updated 07/13/04 by jaina
c**** end prologue ****************************************************


      use params_global
      use freewake

c*******************  list of global variables used **********************
c
c     jmax,kmax,lmax,nblade,izdim,pi -> (params_global)
c     pcx,pcy,pcz,np,ft,nz -> (freewake)
c
c*************************************************************************

      implicit none

      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax),vg(jmax,kmax,lmax),wg(jmax,kmax,lmax)
      real psi

c..   local variables

      integer windex(nblade),ll,idz
      integer inwg(2*izdim),istart(nblade)
      integer w,i,j,k,l,n,p3,p4

      real wgtemp(2*izdim)
      real wgh(2*izdim)
      real dblade
      real ugt,vgt,wgt,gam1,gam2,rc1,hh
      real ugx,ugy,ugz
      real x1,y1,z1,x2,y2,z2

c***  first executable statement

      dblade=2*pi/nblade
      !idz=nint(1.0*np*ft/nz)
      idz=1

      do i=1,nblade
         ll=nint((psi+dblade*(i-1))*np/(2*pi))
         windex(i)=mod(ll,np)+1
      enddo

      istart(1)=nint(nz/ft/12)
      do i=2,nblade
	istart(i)=1
      enddo

      do l=1,lmax
         do k=1,kmax
            do j=1,jmax

	       !write(23,*) x(j,k,l),y(j,k,l),z(j,k,l)

               ugt=0.
               vgt=0.
               wgt=0.

               do i=1,nblade
                  do n=istart(i),nz-1

                     x1=pcx(windex(i),n,w)
                     y1=pcy(windex(i),n,w)
                     z1=pcz(windex(i),n,w)
                     
                     x2=pcx(windex(i),n+1,w)
                     y2=pcy(windex(i),n+1,w)
                     z2=pcz(windex(i),n+1,w)
                     
                     p3=mod(windex(i)-(n-1)*idz-1,np)+1
                     p4=mod(windex(i)-n*idz-1,np)+1
                     
                     if (p3.le.0) p3=np+p3
                     if (p4.le.0) p4=np+p4
                     
                     gam1=circ(p3,w)
                     gam2=circ(p4,w)
                     rc1=rrc(w,n)
		     !write(24,1014) x1,y1,z1,gam1,gam2,rc1                    
                     
                     call influence(x1,y1,z1,x2,y2,z2,gam1,gam2,rc1,
     &                    x(j,k,l),y(j,k,l),z(j,k,l),ugx,ugy,ugz,hh)
                     
                     ugt=ugt+ugx
                     vgt=vgt+ugy
                     wgt=wgt+ugz
                     !wgtemp((i-1)*(nz-1)+n)=ugz
                     !inwg((i-1)*(nz-1)+n)=(i-1)*(nz-1)+n
                     !wgh((i-1)*(nz-1)+n)=hh

                  enddo
               enddo
               
               ug(j,k,l)=ug(j,k,l)+ugt
               vg(j,k,l)=vg(j,k,l)+vgt
               wg(j,k,l)=wg(j,k,l)+wgt
               

            enddo
         enddo
      enddo
      
1014  format(10(1X,E14.8))
      
      return
      end

c***********************************************************************
      subroutine influence(x1,y1,z1,x2,y2,z2,gam1,gam2,rc,xp,yp,zp,
     &     ugx,ugy,ugz,h)

c**** Prologue : *****
c
c     Velocity induced by a vortex filament with end point coordinates
c     (x1,y1,z1) and (x2,y2,z2) and strength varying linearly from gam1
c     to gam2 at point (xp,yp,zp)
c
c     last updated 07/13/04 by jaina
c**** end prologue ****************************************************

      implicit none

      real x1,y1,z1,x2,y2,z2,xp,yp,zp,gam1,gam2,rc,ugx,ugy,ugz
      real l,nx,ny,nz,cz,hx,hy,hz,xa,xb,vvp
      real h,vortex
      
      cz=(x1-xp)*(y2-yp)-(x2-xp)*(y1-yp)

      if (cz.eq.0) then
         ugx=0.
         ugy=0.
         ugz=0.
      else
         l=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1))
         nx=(x2-x1)/l
         ny=(y2-y1)/l
         nz=(z2-z1)/l

         hx=(yp-y1)*nz-ny*(zp-z1)
         hy=(zp-z1)*nx-nz*(xp-x1)
         hz=(xp-x1)*ny-nx*(yp-y1)

         h=sqrt(hx*hx+hy*hy+hz*hz)
         xa=(x1-xp)*nx+(y1-yp)*ny+(z1-zp)*nz
         xb=(x2-xp)*nx+(y2-yp)*ny+(z2-zp)*nz

         vvp=vortex(xa,xb,h,rc)
         vvp=vvp*(gam1+gam2)*0.5

         ugx=vvp*hx/h
         ugy=vvp*hy/h
         ugz=vvp*hz/h
      endif

      return
      end


c***********************************************************************
      function vortex(xa,xb,h,rcore)

c**** Prologue : *****
c     
c     Scully's model for tangential velocity induced by a vortex in 2-D
c     
c     last updated by jaina 07/13/04
c**** end prologue ****************************************************

      implicit none
      
      real xa,xb,h,rcore
      real vortex
      real cos1,cos2
      
      cos1=xa/sqrt(xa*xa+h*h)
      cos2=xb/sqrt(xb*xb+h*h)

      vortex=0.5*(cos2-cos1)*h/(h*h+rcore*rcore)

      return
      end


***********************************************************************
      subroutine fwind_fast3(w,x,y,z,ug,vg,wg,npts,psi)

c***** Prologue : *********
c
c     subroutine for fast evaluation of induced velocities. It was originally
c     insired from the fast multipole method due to Greengard and Rokhlin.
c     The present implementation however does not use multipoles! it uses
c     near field and far field separation, exact evaluation of near field
c     and evaluation of far field effects using trilinear interpolation 
c     from values from quadrature points. The induced velocities at the
c     quadrature points were evaluated exactly.
c
c     This approach was found to give the best performance than many
c     variants of Barnes-Hut, Adaptive FMM etc which used multipoles.
c     Although there is no theoretical bound for the error estimates
c     which makes the accuracy of the approach a suspect.
c     But extensive computations have shown that except for the highly
c     BVI dominated cases the implementation is as accurate as the 
c     methods using multipoles. For the BVI dominated case small oscillations
c     at the inboard stations were observed and were suspected to 
c     originate because of the inaccuracy of this subroutine
c
c     last updated 07/13/04 by jaina
c**** end prologue ****************************************************

      use freewake
      use params_global

c*******************  list of global variables used ********************
c
c     jmax,kmax,lmax,nblade,pi -> (params_global)
c      np,ft,nz,pcx,pcy,pcz -> (free wake)
c
c***********************************************************************


      implicit none

      integer w,npts
      real x(npts),y(npts),z(npts)
      real ug(npts),vg(npts),wg(npts)
      real psi

c..   local variables

      integer jjdim,kkdim,lldim,nndim

      integer, allocatable :: indx(:,:,:,:),numpts(:,:,:)
      real, allocatable :: xq(:),yq(:),zq(:)
      real, allocatable :: uug(:,:,:),vvg(:,:,:),wwg(:,:,:)

      integer ins(10*izdim)

      integer jvert(8),kvert(8),lvert(8),win
      integer i,j,k,l,jj,kk,ll,nn,jq,kq,lq,n,pp,iter,ipts,isum,is1
      integer ji,ki,li,llow,lhig,klow,khig,jlow,jhig,p3,p4
      integer windex(10)
      integer istart(10)
      integer inwg(2*izdim),idz

      real ugp(8),vgp(8),wgp(8),sh(8)
      real xxmax,yymax,zzmax,xxmin,yymin,zzmin
      real dxmax,dymax,dzmax,cs,ss
      real gam1,gam2,rc1,hh
      real dblade
      real xq1,yq1,zq1,ugt,vgt,wgt
      real x1,y1,z1,x2,y2,z2,ugx,ugy,ugz,xp,yp,zp
      real xlen,ylen,zlen,zeta,eta,nu,ug1,vg1,wg1
      real sbox,xp1,yp1,zp1
      real xmax,xmin
      real wgtemp(2*izdim)
      real wgh(2*izdim)

      real xbar,ybar,zbar
      real aa(3,3),eigenv(3),sume,trace
      integer ier
      real core_fac
      

c**** first executable statement


      dblade=2*pi/nblade
      idz=nint(1.0*np*ft/nz)

      cs=cos(psi)
      ss=sin(psi)

      do i=1,nblade
         ll=nint((psi+dblade*(i-1))*np/(2*pi))
         windex(i)=mod(ll,np)+1
      enddo

c..   use only after 30 degrees of wake-age for the first blade
c..   because CFD gives effect of near wake anyway (hopefully ;-)) )
c..   factor 12, because 30=360/12
      
      istart(1)=nint(nz/ft/12.0)
      do i=2,nblade
         istart(i)=1
      enddo

c..   find mean

      xbar=0.
      ybar=0.
      zbar=0.

      do ipts=1,npts
         xbar=xbar+x(ipts)
         ybar=ybar+y(ipts)
         zbar=zbar+z(ipts)
      enddo

      xbar=xbar/npts
      ybar=ybar/npts
      zbar=zbar/npts
      
      aa=0.

      do ipts=1,npts
         aa(1,1)=aa(1,1)+(x(ipts)-xbar)*(x(ipts)-xbar)
         aa(1,2)=aa(1,2)+(x(ipts)-xbar)*(y(ipts)-ybar)
         aa(1,3)=aa(1,3)+(x(ipts)-xbar)*(z(ipts)-zbar)
         
         aa(2,1)=aa(1,2)+(y(ipts)-ybar)*(x(ipts)-xbar)
         aa(2,2)=aa(2,2)+(y(ipts)-ybar)*(y(ipts)-ybar)
         aa(2,3)=aa(2,3)+(y(ipts)-ybar)*(z(ipts)-zbar)
         
         aa(3,1)=aa(3,1)+(z(ipts)-zbar)*(x(ipts)-xbar)
         aa(3,2)=aa(3,2)+(z(ipts)-zbar)*(y(ipts)-ybar)
         aa(3,3)=aa(3,3)+(z(ipts)-zbar)*(z(ipts)-zbar)
      enddo

c..  kaisers method for finding eigen values of real,symmetric matrices
c..  remember eigen vectors are principal axes (inertial bisection)      
 
      call kaiser(aa,3,3,eigenv,trace,sume,ier)

c..  now find the largest distance in the principal directions

      xxmax=x(1)*aa(1,1)+y(1)*aa(2,1)+z(1)*aa(3,1)
      xxmin=xxmax
      yymax=x(1)*aa(1,2)+y(1)*aa(2,2)+z(1)*aa(3,2)
      yymin=yymax
      zzmax=x(1)*aa(1,3)+y(1)*aa(2,3)+z(1)*aa(3,3)
      zzmin=zzmax
      
      do ipts=1,npts

         xq1=x(ipts)*aa(1,1)+y(ipts)*aa(2,1)+z(ipts)*aa(3,1)
         yq1=x(ipts)*aa(1,2)+y(ipts)*aa(2,2)+z(ipts)*aa(3,2)
         zq1=x(ipts)*aa(1,3)+y(ipts)*aa(2,3)+z(ipts)*aa(3,3)
         
         if (xxmax.lt.xq1)  xxmax=xq1
         if (xxmin.gt.xq1)  xxmin=xq1
         
         if (yymax.lt.yq1)  yymax=yq1
         if (yymin.gt.yq1)  yymin=yq1
         
         if (zzmax.lt.zq1)  zzmax=zq1
         if (zzmin.gt.zq1)  zzmin=zq1
         
      enddo
      

c..   heuristic non-dimensional box size (found to be optimal)

c      sbox=(yymax-yymin)/15.

      sbox=2. ! two chords (seems to be ideal)

      xxmax=xxmax+sbox
      xxmin=xxmin-sbox
      yymax=yymax+sbox
      yymin=yymin-sbox
      zzmax=zzmax+sbox
      zzmin=zzmin-sbox

      jq=aint((xxmax-xxmin)/sbox)+1
      kq=aint((yymax-yymin)/sbox)+1
      lq=aint((zzmax-zzmin)/sbox)+1

      allocate(indx(jq,kq,lq,2),numpts(jq,kq,lq))
      allocate(uug(jq,kq,lq),vvg(jq,kq,lq),wwg(jq,kq,lq))
      allocate(xq(jq),yq(kq),zq(lq))

      dxmax=(xxmax-xxmin)/(jq-1)
      dymax=(yymax-yymin)/(kq-1)
      dzmax=(zzmax-zzmin)/(lq-1)

      ! find box vertices

      do j=1,jq
         xq(j)=xxmin+dxmax*(j-1)
      enddo
      
      do k=1,kq
         yq(k)=yymin+dymax*(k-1)
      enddo

      do l=1,lq
         zq(l)=zzmin+dzmax*(l-1)
      enddo

      ! do exact evaluation at the box vertices

      do l=1,lq
         do k=1,kq
            do j=1,jq

               ugt=0.
               vgt=0.
               wgt=0.

               xq1=xq(j)*aa(1,1)+yq(k)*aa(1,2)+zq(l)*aa(1,3)
               yq1=xq(j)*aa(2,1)+yq(k)*aa(2,2)+zq(l)*aa(2,3)
               zq1=xq(j)*aa(3,1)+yq(k)*aa(3,2)+zq(l)*aa(3,3)

               do i=1,nblade
                  do n=istart(i),nz-1
                     
                     x1=pcx(windex(i),n,w)
                     y1=pcy(windex(i),n,w)
                     z1=pcz(windex(i),n,w)
                     
                     x2=pcx(windex(i),n+1,w)
                     y2=pcy(windex(i),n+1,w)
                     z2=pcz(windex(i),n+1,w)
                     
                     p3=mod(windex(i)-(n-1)*idz-1,np)+1
                     p4=mod(windex(i)-n*idz-1,np)+1

                     
                     if (p3.le.0) p3=np+p3
                     if (p4.le.0) p4=np+p4
                     
                     gam1=circ(p3,w)
                     gam2=circ(p4,w)
                     rc1=rrc(w,n)

                     call influence(x1,y1,z1,x2,y2,z2,gam1,gam2,rc1,
     &                    xq1,yq1,zq1,ugx,ugy,ugz,hh)
                     
                     ugt=ugt+ugx
                     vgt=vgt+ugy
                     wgt=wgt+ugz

                  enddo
               enddo
               
               uug(j,k,l)=ugt
               vvg(j,k,l)=vgt
               wwg(j,k,l)=wgt
               
            enddo
         enddo
      enddo

      do j=1,jq
         do k=1,kq
            do l=1,lq
               indx(j,k,l,1)=0
            enddo
         enddo
      enddo

      do i=1,nblade
         do n=istart(i),nz-1
            
            xp=pcx(windex(i),n,w)*aa(1,1)
     &         +pcy(windex(i),n,w)*aa(2,1)
     &         +pcz(windex(i),n,w)*aa(3,1)

            yp=pcx(windex(i),n,w)*aa(1,2)
     &         +pcy(windex(i),n,w)*aa(2,2)
     &         +pcz(windex(i),n,w)*aa(3,2)

            zp=pcx(windex(i),n,w)*aa(1,3)
     &         +pcy(windex(i),n,w)*aa(2,3)
     &         +pcz(windex(i),n,w)*aa(3,3)
 
            ji=aint((xp-xxmin)/sbox)+1
            ki=aint((yp-yymin)/sbox)+1
            li=aint((zp-zzmin)/sbox)+1
           
            if (ji.ge.1.and.ji.lt.jq.and.ki.ge.1.and.ki.lt.kq.and.
     &           li.ge.1.and.li.lt.lq) then
               
                 indx(ji,ki,li,1)=indx(ji,ki,li,1)+1
            endif
         enddo
      enddo

      isum=1
      do j=1,jq
         do k=1,kq
            do l=1,lq
               indx(j,k,l,2)=isum
               isum=isum+indx(j,k,l,1)
               numpts(j,k,l)=indx(j,k,l,1)
            enddo
         enddo
      enddo

      do i=1,nblade
         do n=istart(i),nz-1
            
            xp=pcx(windex(i),n,w)*aa(1,1)
     &         +pcy(windex(i),n,w)*aa(2,1)
     &         +pcz(windex(i),n,w)*aa(3,1)

            yp=pcx(windex(i),n,w)*aa(1,2)
     &         +pcy(windex(i),n,w)*aa(2,2)
     &         +pcz(windex(i),n,w)*aa(3,2)

            zp=pcx(windex(i),n,w)*aa(1,3)
     &         +pcy(windex(i),n,w)*aa(2,3)
     &         +pcz(windex(i),n,w)*aa(3,3)

            ji=aint((xp-xxmin)/sbox)+1
            ki=aint((yp-yymin)/sbox)+1
            li=aint((zp-zzmin)/sbox)+1
     
           if (ji.ge.1.and.ji.lt.jq.and.ki.ge.1.and.ki.lt.kq.and.
     &           li.ge.1.and.li.lt.lq) then
               
               is1=indx(ji,ki,li,2)+indx(ji,ki,li,1)-numpts(ji,ki,li)
               numpts(ji,ki,li)=numpts(ji,ki,li)-1
               ins(is1)=(i-1)*(nz-1)+n
            endif
         enddo
      enddo

c..  end reordering
      
      do ipts=1,npts

         xp=x(ipts)*aa(1,1)+y(ipts)*aa(2,1)+z(ipts)*aa(3,1)
         yp=x(ipts)*aa(1,2)+y(ipts)*aa(2,2)+z(ipts)*aa(3,2)
         zp=x(ipts)*aa(1,3)+y(ipts)*aa(2,3)+z(ipts)*aa(3,3)

         ji=aint((xp-xxmin)/sbox)+1
         ki=aint((yp-yymin)/sbox)+1
         li=aint((zp-zzmin)/sbox)+1
         
         if (ji.eq.jq) ji=jq-1
         if (ki.eq.kq) ki=kq-1
         if (li.eq.lq) li=lq-1
         
         ugt=0.
         vgt=0.
         wgt=0.
         
         call setvertices(jvert,kvert,lvert,ji,ki,li,jq,kq,lq)
         
         do pp=1,8
            ugp(pp)=0.
            vgp(pp)=0.
            wgp(pp)=0.
         enddo
         
         llow=li-1
         lhig=li+1
         
         klow=ki-1
         khig=ki+1
         
         jlow=ji-1
         jhig=ji+1
         
         if (llow.lt.1) llow=1
         if (lhig.gt.lq-1) lhig=lq-1
         
         if (klow.lt.1) klow=1
         if (khig.gt.kq-1) khig=kq-1
         
         if (jlow.lt.1) jlow=1
         if (jhig.gt.jq-1) jhig=jq-1
         
         iter=0
         
         do ll=llow,lhig
            do kk=klow,khig
               do jj=jlow,jhig
                  
                  do nn=1,indx(jj,kk,ll,1)
                     
                     is1=indx(jj,kk,ll,2)+nn-1

                     win=(ins(is1)-1)/(nz-1)+1
                     n=mod(ins(is1)-1,nz-1)+1
                     
                     x1=pcx(windex(win),n,w)
                     y1=pcy(windex(win),n,w)
                     z1=pcz(windex(win),n,w)
                     
                     x2=pcx(windex(win),n+1,w)
                     y2=pcy(windex(win),n+1,w)
                     z2=pcz(windex(win),n+1,w)
                     
                     p3=mod(windex(win)-(n-1)*idz-1,np)+1
                     p4=mod(windex(win)-n*idz-1,np)+1
                     
                     if (p3.le.0) p3=np+p3
                     if (p4.le.0) p4=np+p4
                     
                     gam1=circ(p3,w)
                     gam2=circ(p4,w)
                     rc1=rrc(w,n)
                     
                     xq1=x(ipts)
                     yq1=y(ipts)
                     zq1=z(ipts)
                     
                     call influence(x1,y1,z1,x2,y2,z2,gam1,gam2,
     &                    rc1,xq1,yq1,zq1,ugx,ugy,ugz,hh)
                     
                     ugt=ugt+ugx
                     vgt=vgt+ugy
                     wgt=wgt+ugz
                     
                     iter=iter+1
                     
                     wgtemp(iter)=ugz
                     inwg(iter)=ins(is1)
                     wgh(iter)=hh

                     do pp=1,8
                        
                         
                        xq1=xq(jvert(pp))*aa(1,1)
     &                       +yq(kvert(pp))*aa(1,2)
     &                       +zq(lvert(pp))*aa(1,3)
                            
                        yq1=xq(jvert(pp))*aa(2,1)
     &                      +yq(kvert(pp))*aa(2,2)
     &                      +zq(lvert(pp))*aa(2,3)

                        zq1=xq(jvert(pp))*aa(3,1)
     &                      +yq(kvert(pp))*aa(3,2)
     &                      +zq(lvert(pp))*aa(3,3)
                        
                        call influence(x1,y1,z1,x2,y2,z2,gam1,
     $                       gam2,
     &                       rc1,xq1,yq1,zq1,ugx,ugy,ugz,hh)
                        
                        ugp(pp)=ugp(pp)+ugx
                        vgp(pp)=vgp(pp)+ugy
                        wgp(pp)=wgp(pp)+ugz
                        
                     enddo
                     
                  enddo
               enddo
            enddo
         enddo
         
         xlen=abs(xq(jvert(1))-xq(jvert(2)))
         ylen=abs(yq(kvert(1))-yq(kvert(5)))
         zlen=abs(zq(lvert(1))-zq(lvert(3)))

         xp1=x(ipts)*aa(1,1)+y(ipts)*aa(2,1)+z(ipts)*aa(3,1)
         yp1=x(ipts)*aa(1,2)+y(ipts)*aa(2,2)+z(ipts)*aa(3,2)
         zp1=x(ipts)*aa(1,3)+y(ipts)*aa(2,3)+z(ipts)*aa(3,3)
         
         zeta=(xp1-xq(jvert(1)))/xlen
         eta=(yp1-yq(kvert(1)))/ylen
         nu=(zp1-zq(lvert(1)))/zlen
         
         do pp=1,8
            ugp(pp)=uug(jvert(pp),kvert(pp),lvert(pp))-ugp(pp)
            vgp(pp)=vvg(jvert(pp),kvert(pp),lvert(pp))-vgp(pp)
            wgp(pp)=wwg(jvert(pp),kvert(pp),lvert(pp))-wgp(pp)
         enddo
	       
         sh(1)=(1-zeta)*(1-eta)*(1-nu)
         sh(2)=zeta*(1-eta)*(1-nu)
         sh(3)=zeta*(1-eta)*nu
         sh(4)=(1-zeta)*(1-eta)*nu
         
         sh(5)=(1-zeta)*eta*(1-nu)
         sh(6)=zeta*eta*(1-nu)
         sh(7)=zeta*eta*nu
         sh(8)=(1-zeta)*eta*nu
         
         ug1=0.
         vg1=0.
         wg1=0.
         
         do pp=1,8
            ug1=ug1+ugp(pp)*sh(pp)
            vg1=vg1+vgp(pp)*sh(pp)
            wg1=wg1+wgp(pp)*sh(pp)
         enddo
         
         ug(ipts)=ug(ipts)+ug1+ugt
         vg(ipts)=vg(ipts)+vg1+vgt
         wg(ipts)=wg(ipts)+wg1+wgt
         
      enddo

      deallocate(indx,numpts)
      deallocate(xq,yq,zq)
      deallocate(uug,vvg,wwg)

      return
      end

c***********************************************************************
      subroutine fwind_fast(w,x,y,z,ug,vg,wg,psi)

c***** Prologue : *********
c
c     subroutine for fast evaluation of induced velocities. It was originally
c     insired from the fast multipole method due to Greengard and Rokhlin.
c     The present implementation however does not use multipoles! it uses
c     near field and far field separation, exact evaluation of near field
c     and evaluation of far field effects using trilinear interpolation 
c     from values from quadrature points. The induced velocities at the
c     quadrature points were evaluated exactly.
c
c     This approach was found to give the best performance than many
c     variants of Barnes-Hut, Adaptive FMM etc which used multipoles.
c     Although there is no theoretical bound for the error estimates
c     which makes the accuracy of the approach a suspect.
c     But extensive computations have shown that except for the highly
c     BVI dominated cases the implementation is as accurate as the
c     methods using multipoles. For the BVI dominated case small oscillations
c     at the inboard stations were observed and were suspected to
c     originate because of the inaccuracy of this subroutine
c
c     last updated 07/13/04 by jaina
c**** end prologue ****************************************************
      use params_global
      use freewake
	
c*******************  list of global variables used ********************
c
c     jmax,kmax,lmax,nblade,pi -> (params_global)
c      np,ft,nz,pcx,pcy,pcz -> (free wake)
c
c***********************************************************************

      implicit none

      integer w
      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax),vg(jmax,kmax,lmax),wg(jmax,kmax,lmax)
      real psi

c..   local variables

      integer jjdim,kkdim,lldim,nndim
      parameter(jjdim=60,kkdim=60,lldim=60,nndim=60)

      integer ins(jjdim,kkdim,lldim,nndim)
      integer jvert(8),kvert(8),lvert(8),win
      integer i,j,k,l,jj,kk,ll,nn,jq,kq,lq,n,pp,iter
      integer ji,ki,li,llow,lhig,klow,khig,jlow,jhig,p3,p4
      integer windex(nblade)
      integer istart(nblade)
      integer inwg(2*izdim),idz

      real xq(jjdim),yq(kkdim),zq(lldim)
      real uug(jjdim,kkdim,lldim),vvg(jjdim,kkdim,lldim)
      real wwg(jjdim,kkdim,lldim)
      real ugp(8),vgp(8),wgp(8),sh(8)
      real xxmax,yymax,zzmax,xxmin,yymin,zzmin
      real dxmax,dymax,dzmax,cs,ss
      real gam1,gam2,rc1,hh
      real dblade
      real xq1,yq1,zq1,ugt,vgt,wgt
      real x1,y1,z1,x2,y2,z2,ugx,ugy,ugz,xp,yp,zp
      real xlen,ylen,zlen,zeta,eta,nu,ug1,vg1,wg1
      real sbox,xp1,yp1,zp1
      real xmax,xmin
      real wgtemp(2*izdim)
      real wgh(2*izdim)


c**** first executable statement

      dblade=2*pi/nblade
      idz=nint(1.0*np*ft/nz)

      cs=cos(psi)
      ss=sin(psi)

      do i=1,nblade
         ll=nint((psi+dblade*(i-1))*np/(2*pi))
         windex(i)=mod(ll,np)+1
      enddo

c..   use only after 30 degrees of wake-age for the first blade
c..   because CFD gives effect of near wake anyway (hopefully ;-)) )
c..   factor 12, because 30=360/12
      
      istart(1)=nint(nz/ft/12.0)
      do i=2,nblade
         istart(i)=1
      enddo

      xxmax=x(1,1,1)*cs+y(1,1,1)*ss
      xxmin=xxmax
      yymax=y(1,1,1)*cs-x(1,1,1)*ss
      yymin=yymax
      zzmax=z(1,1,1)
      zzmin=zzmax

      do l=1,lmax
         do k=1,kmax
            do j=1,jmax
               
               xq1=x(j,k,l)*cs+y(j,k,l)*ss
               yq1=y(j,k,l)*cs-x(j,k,l)*ss
               zq1=z(j,k,l)

               if (xxmax.lt.xq1)  xxmax=xq1
               if (xxmin.gt.xq1)  xxmin=xq1
               
               if (yymax.lt.yq1)  yymax=yq1
               if (yymin.gt.yq1)  yymin=yq1
               
               if (zzmax.lt.zq1)  zzmax=zq1
               if (zzmin.gt.zq1)  zzmin=zq1
            enddo
         enddo
      enddo

c..   heuristic non-dimensional box size (found to be optimal)

      sbox=2.

      jq=aint((xxmax-xxmin)/sbox)+1
      kq=aint((yymax-yymin)/sbox)+1
      lq=aint((zzmax-zzmin)/sbox)+1

      dxmax=(xxmax-xxmin)/(jq-1)
      dymax=(yymax-yymin)/(kq-1)
      dzmax=(zzmax-zzmin)/(lq-1)

      do j=1,jq
         xq(j)=xxmin+dxmax*(j-1)
      enddo
      
      do k=1,kq
         yq(k)=yymin+dymax*(k-1)
      enddo

      do l=1,lq
         zq(l)=zzmin+dzmax*(l-1)
      enddo

      do l=1,lq
         do k=1,kq
            do j=1,jq

               ugt=0.
               vgt=0.
               wgt=0.

               xq1=xq(j)*cs-yq(k)*ss
               yq1=yq(k)*cs+xq(j)*ss
               zq1=zq(l)
               
               do i=1,nblade
                  do n=istart(i),nz-1
                     
                     x1=pcx(windex(i),n,w)
                     y1=pcy(windex(i),n,w)
                     z1=pcz(windex(i),n,w)
                     
                     x2=pcx(windex(i),n+1,w)
                     y2=pcy(windex(i),n+1,w)
                     z2=pcz(windex(i),n+1,w)
                     
c                     p3=mod(windex(i)-n,np)+1
c                     p4=mod(windex(i)-n-1,np)+1

                     p3=mod(windex(i)-(n-1)*idz-1,np)+1
                     p4=mod(windex(i)-n*idz-1,np)+1
                     
                     if (p3.le.0) p3=np+p3
                     if (p4.le.0) p4=np+p4
                     
                     gam1=circ(p3,w)
                     gam2=circ(p4,w)
                     rc1=rrc(w,n)

                     call influence(x1,y1,z1,x2,y2,z2,gam1,gam2,rc1,
     &                    xq1,yq1,zq1,ugx,ugy,ugz,hh)
                     
                     ugt=ugt+ugx
                     vgt=vgt+ugy
                     wgt=wgt+ugz

                  enddo
               enddo
               
               uug(j,k,l)=ugt
               vvg(j,k,l)=vgt
               wwg(j,k,l)=wgt
               
            enddo
         enddo
      enddo

      do j=1,jq
         do k=1,kq
            do l=1,lq
               ins(j,k,l,1)=0
            enddo
         enddo
      enddo

      do i=1,nblade
         do n=istart(i),nz-1
            
            xp=pcx(windex(i),n,w)*cs+pcy(windex(i),n,w)*ss
            yp=pcy(windex(i),n,w)*cs-pcx(windex(i),n,w)*ss
            zp=pcz(windex(i),n,w)

            ji=aint((xp-xxmin)/sbox)+1
            ki=aint((yp-yymin)/sbox)+1
            li=aint((zp-zzmin)/sbox)+1
            
            if (ji.ge.1.and.ji.lt.jq.and.ki.ge.1.and.ki.lt.kq.and.
     &           li.ge.1.and.li.lt.lq) then
               
               ins(ji,ki,li,1)=ins(ji,ki,li,1)+1
               ins(ji,ki,li,ins(ji,ki,li,1)+1)=(i-1)*(nz-1)+n
            endif
         enddo
      enddo

      
      do l=1,lmax
         do k=1,kmax
            do j=1,jmax

               ji=aint((x(j,k,l)*cs+y(j,k,l)*ss-xxmin)/sbox)+1
               ki=aint((y(j,k,l)*cs-x(j,k,l)*ss-yymin)/sbox)+1
               li=aint((z(j,k,l)-zzmin)/sbox)+1

               if (ji.eq.jq) ji=jq-1
               if (ki.eq.kq) ki=kq-1
               if (li.eq.lq) li=lq-1

               ugt=0.
               vgt=0.
               wgt=0.

               call setvertices(jvert,kvert,lvert,ji,ki,li,jq,kq,lq)

               do pp=1,8
                  ugp(pp)=0.
                  vgp(pp)=0.
                  wgp(pp)=0.
               enddo
               
               llow=li-1
               lhig=li+1
               
               klow=ki-1
               khig=ki+1
               
               jlow=ji-1
               jhig=ji+1

               if (llow.lt.1) llow=1
               if (lhig.gt.lq-1) lhig=lq-1

               if (klow.lt.1) klow=1
               if (khig.gt.kq-1) khig=kq-1
               
               if (jlow.lt.1) jlow=1
               if (jhig.gt.jq-1) jhig=jq-1

               iter=0

               do ll=llow,lhig
                  do kk=klow,khig
                     do jj=jlow,jhig

                        do nn=2,ins(jj,kk,ll,1)+1
                           
			   win=(ins(jj,kk,ll,nn)-1)/(nz-1)+1
			   n=mod(ins(jj,kk,ll,nn)-1,nz-1)+1

			   x1=pcx(windex(win),n,w)
			   y1=pcy(windex(win),n,w)
			   z1=pcz(windex(win),n,w)
                           
			   x2=pcx(windex(win),n+1,w)
			   y2=pcy(windex(win),n+1,w)
			   z2=pcz(windex(win),n+1,w)
                           
			   p3=mod(windex(win)-(n-1)*idz-1,np)+1
			   p4=mod(windex(win)-n*idz-1,np)+1
                           
			   if (p3.le.0) p3=np+p3
			   if (p4.le.0) p4=np+p4
                           
			   gam1=circ(p3,w)
			   gam2=circ(p4,w)
			   rc1=rrc(w,n)

			   xq1=x(j,k,l)
			   yq1=y(j,k,l)
			   zq1=z(j,k,l)

			   call influence(x1,y1,z1,x2,y2,z2,gam1,gam2,
     &                          rc1,xq1,yq1,zq1,ugx,ugy,ugz,hh)

			   ugt=ugt+ugx
			   vgt=vgt+ugy
			   wgt=wgt+ugz

			   iter=iter+1
			   
			   wgtemp(iter)=ugz
			   inwg(iter)=ins(jj,kk,ll,nn)
			   wgh(iter)=hh

			   do pp=1,8

			      xq1=xq(jvert(pp))*cs-yq(kvert(pp))*ss
			      yq1=yq(kvert(pp))*cs+xq(jvert(pp))*ss
			      zq1=zq(lvert(pp))
			      
			      call influence(x1,y1,z1,x2,y2,z2,gam1,
     $                             gam2,
     &                             rc1,xq1,yq1,zq1,ugx,ugy,ugz,hh)
			      
			      ugp(pp)=ugp(pp)+ugx
			      vgp(pp)=vgp(pp)+ugy
			      wgp(pp)=wgp(pp)+ugz

			   enddo
                           
			enddo
		     enddo
		  enddo
	       enddo

	       xlen=abs(xq(jvert(1))-xq(jvert(2)))
	       ylen=abs(yq(kvert(1))-yq(kvert(5)))
	       zlen=abs(zq(lvert(1))-zq(lvert(3)))
               
	       xp1=x(j,k,l)*cs+y(j,k,l)*ss
	       yp1=y(j,k,l)*cs-x(j,k,l)*ss
	       zp1=z(j,k,l)

	       zeta=(xp1-xq(jvert(1)))/xlen
	       eta=(yp1-yq(kvert(1)))/ylen
	       nu=(zp1-zq(lvert(1)))/zlen
	       
	       do pp=1,8
		  ugp(pp)=uug(jvert(pp),kvert(pp),lvert(pp))-ugp(pp)
		  vgp(pp)=vvg(jvert(pp),kvert(pp),lvert(pp))-vgp(pp)
		  wgp(pp)=wwg(jvert(pp),kvert(pp),lvert(pp))-wgp(pp)
	       enddo
	       
	       sh(1)=(1-zeta)*(1-eta)*(1-nu)
	       sh(2)=zeta*(1-eta)*(1-nu)
	       sh(3)=zeta*(1-eta)*nu
	       sh(4)=(1-zeta)*(1-eta)*nu
	       
	       sh(5)=(1-zeta)*eta*(1-nu)
	       sh(6)=zeta*eta*(1-nu)
	       sh(7)=zeta*eta*nu
	       sh(8)=(1-zeta)*eta*nu
	       
	       ug1=0.
	       vg1=0.
	       wg1=0.
	       
	       do pp=1,8
		  ug1=ug1+ugp(pp)*sh(pp)
		  vg1=vg1+vgp(pp)*sh(pp)
		  wg1=wg1+wgp(pp)*sh(pp)
	       enddo
	       

	       ug(j,k,l)=ug(j,k,l)+ug1+ugt
	       vg(j,k,l)=vg(j,k,l)+vg1+vgt
	       wg(j,k,l)=wg(j,k,l)+wg1+wgt
	       
	    enddo
	 enddo
      enddo

      return
      end

c***********************************************************************      
      subroutine setvertices(jvert,kvert,lvert,ji,ki,li,jq,kq,lq)

c***** Prologue : *********
c
c     set vertices of imaginary cube with body diagonal verices
c     (ji,ki,li) and (ji+1,ki+1,li+1)
c     stuff done in a facewise counter-clockwise sense
c
c     last updated 07/13/04 by jaina
c**** end prologue *****************************************************

      implicit none

      integer jvert(8),kvert(8),lvert(8)
      integer ji,ki,li,jq,kq,lq,pp

      jvert(1)=ji
      jvert(4)=ji
      jvert(5)=ji
      jvert(8)=ji
      
      jvert(2)=ji+1
      jvert(3)=ji+1
      jvert(6)=ji+1
      jvert(7)=ji+1
      
      kvert(1)=ki
      kvert(2)=ki
      kvert(3)=ki
      kvert(4)=ki
      
      kvert(5)=ki+1
      kvert(6)=ki+1
      kvert(7)=ki+1
      kvert(8)=ki+1
      
      lvert(1)=li
      lvert(2)=li
      lvert(5)=li
      lvert(6)=li
      
      lvert(3)=li+1
      lvert(4)=li+1
      lvert(7)=li+1
      lvert(8)=li+1
      
      do pp=1,8
         if (lvert(pp).gt.lq) lvert(pp)=lq
         if (kvert(pp).gt.kq) kvert(pp)=kq
         if (jvert(pp).gt.jq) jvert(pp)=jq
         
         if (lvert(pp).lt.1) lvert(pp)=1
         if (kvert(pp).lt.1) kvert(pp)=1
         if (jvert(pp).lt.1) jvert(pp)=1
      enddo
      
      return
      end
      

      subroutine user_specified_field(x,y,z,ug,vg,wg,npts)

      use params_global

      implicit none
!
!     subroutine arguments
!
      integer npts
      real x(npts),y(npts),z(npts)
      real ug(npts),vg(npts),wg(npts)
      
!
!     local variables
!
      integer i
      real xx,yy,zz,uug,vvg,wwg

      do i=1,npts
	
	xx=y(i)/rartio
	yy=-x(i)/rartio
	zz=z(i)/rartio
	
	!!call gust(xx,yy,zz,uug,vvg,wwg)

	ug(i)=vvg*fmtip
	vg(i)=-uug*fmtip
	wg(i)=wwg*fmtip
      enddo	

      return
      end
