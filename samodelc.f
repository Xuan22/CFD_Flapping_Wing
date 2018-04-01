c***********************************************************************
      subroutine vmu_sa(x,y,z,q,turmu,xx,xy,xz,yx,yy,yz,zx,zy,zz,
     >			ug,vg,wg,vnu,vnu0,vnul,tscale,iblank,grbc)
c  turbulent eddy viscosity. model is one equation spalart-
c  allmaras. ref{ aiaa 92-0439}.
c
c..notes: can manage memory much better than present
c	     - remove aj...cl and u,v,w
c***********************************************************************

      use params_global
      use bcparam
      implicit none

      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)
      real turmu(jmax,kmax,lmax)
      real q(jmax,kmax,lmax,nd)
      real xx(jmax,kmax,lmax), xy(jmax,kmax,lmax),xz(jmax,kmax,lmax),
     $     yx(jmax,kmax,lmax), yy(jmax,kmax,lmax),yz(jmax,kmax,lmax),
     $     zx(jmax,kmax,lmax), zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax),vg(jmax,kmax,lmax),wg(jmax,kmax,lmax)
      real vnu(jmax,kmax,lmax),vnul(jmax,kmax,lmax)
      real vnu0(jmax,kmax,lmax),tscale(jmax,kmax,lmax)
      integer iblank(jmax,kmax,lmax)
      type(bc_t) :: grbc


c***********************************

c..   local variables

      real, allocatable :: sn(:,:,:),ts(:,:,:),u(:,:,:),v(:,:,:),
     >     w(:,:,:),vort(:,:,:),fv1(:,:,:),trip(:,:,:),
     >     sjmat(:,:,:),amat(:,:),bmat(:,:),cmat(:,:),tscal(:,:,:)

      real relfac,vnulim,dl2norm,resl2,dtpseudo_turb

      integer nturiter,niter
      integer j,k,l

      allocate(sn(jmax,kmax,lmax),ts(jmax,kmax,lmax),u(jmax,kmax,lmax),
     <     v(jmax,kmax,lmax),w(jmax,kmax,lmax),vort(jmax,kmax,lmax),
     <     fv1(jmax,kmax,lmax),
     <     trip(jmax,kmax,lmax),sjmat(jmax,kmax,lmax),amat(mdim,mdim),
     <     bmat(mdim,mdim),cmat(mdim,mdim),tscal(jmax,kmax,lmax))


      relfac=1.0
      vnulim=1.0e-10
      !nturiter=itnmax
      nturiter=1
c      cfltur=1.0
c      idual=1
c      dtpseudo=1.0	
	
c..   laminar co-efficient of viscosity calculation
      
      call lamvis(q,vnul)

c..   for compatibility with turbulence model

      do  l = 1,lmax
         do  k = 1,kmax
            do  j = 1,jmax
               u(j,k,l)  = q(j,k,l,2)/q(j,k,l,1)
               v(j,k,l)  = q(j,k,l,3)/q(j,k,l,1)
               w(j,k,l)  = q(j,k,l,4)/q(j,k,l,1)
               ts(j,k,l) = 0.0
               vort(j,k,l)  = 0.0
               sn(j,k,l)  = 1.e10
               trip(j,k,l)= 0.0
               tscal(j,k,l)=tscale(j,k,l)
            enddo
         enddo
      enddo

      call turbc(vnu,u,v,w,ug,vg,wg,xx,xy,xz,yx,yy,yz,zx,zy,zz,grbc)

c..   compute vorticity & distance function

      call vortic(vort,u,v,w,xx,xy,xz,yx,yy,yz,zx,zy,zz)
      
      if (is_wing) call dist(sn,x,y,z)
      if (is_ground) call dist_ground(sn,x,y,z)
      call c_fv1(vnu,q,vnul,fv1)

c..   modify pseudo time-step size for time accurate runs for stability
      
      if (timeac.eq.1.and.1.eq.0) then
        do  l = 1,lmax
          do  k = 1,kmax
            do  j = 1,jmax
              dtpseudo_turb = 1.0
              tscal(j,k,l) = max(iblank(j,k,l),0)* 
     &              ( 1.0 + 0.002*sqrt(q(j,k,l,6)))/(1.+sqrt(q(j,k,l,6)))
              tscal(j,k,l) = tscal(j,k,l)*dtpseudo_turb
              tscal(j,k,l) = tscal(j,k,l)/(1.+tscal(j,k,l)/h)
            enddo
          enddo
        enddo
      endif
      
c..   start sub-iterations (note that rhs*dt is done in rhs)

      do 30 niter=1,nturiter

        call rhstur(q,fv1,sn,vort,u,v,w,ug,vg,wg,vnu,vnu0,vnul,trip,
     &	ts,xx,xy,xz,yx,yy,yz,zx,zy,zz,x,y,z,sjmat,dl2norm,
     &  tscal)
	  
        call lhstur(q,fv1,sn,vort,u,v,w,ug,vg,wg,vnu,vnul,trip,ts,
     &	xx,xy,xz,yx,yy,yz,zx,zy,zz,x,y,z,sjmat,amat,bmat,cmat,
     &	tscal)
	
c..   Update and find norms

        resl2=0.0	
        do l = 2,lmax-1
           do k = 2,kmax-1
              do j = 2,jmax-1
                 ts(j,k,l) = relfac*ts(j,k,l)
                 vnu(j,k,l) = max((vnu(j,k,l) + ts(j,k,l)),vnulim)
                 resl2=resl2+ts(j,k,l)**2
              enddo
           enddo
        enddo

	resl2  =sqrt(resl2/jmax/kmax/lmax)
	dl2norm=sqrt(dl2norm/jmax/kmax/lmax)

        call turbc(vnu,u,v,w,ug,vg,wg,xx,xy,xz,yx,yy,yz,zx,zy,zz,grbc)

30	continue
        
        call c_fv1(vnu,q,vnul,fv1)

        call c_turm(fv1,turmu,vnu,q)
        
      ! Let's clean up memory before exiting
      deallocate(sn,ts,u,v,w,vort,fv1,trip,sjmat,amat,bmat,cmat)
      
      return
      end


c*************************************************************
	subroutine lamvis(q,vnul)


        use params_global
        implicit none
        
        real q(jmax,kmax,lmax,nd),vnul(jmax,kmax,lmax)

c..   local variables

        real gkpr,prtr,dre,c2b,c2bp
        real ra,uvel,vvel,wvel,ei,tt,vmul
        
        integer j,k,l
        
c..   first executable statement

c..   sutherland's formulae

	gkpr = gamma/pr
        prtr = pr/0.9
        dre  = .5/rey
        c2b  =198.6/tinf
        c2bp = c2b +1.
	do 10 l = 1,lmax
        do 10 k = 1,kmax
        do 10 j = 1,jmax
           ra    = 1./q(j,k,l,1)
           uvel  = q(j,k,l,2)*ra
           vvel  = q(j,k,l,3)*ra
           wvel  = q(j,k,l,4)*ra
           ei    = q(j,k,l,5)*ra-0.5*(uvel**2+vvel**2+wvel**2)
           tt    = ggm1*ei
           vmul  = c2bp*tt*sqrt(tt)/(c2b + tt)
	   vnul(j,k,l)=vmul*ra/q(j,k,l,6)
10      continue

	return
        end

c*************************************************************
      subroutine dist(sn,x,y,z)

      use params_global
      implicit none
      
      real sn(jmax,kmax,lmax)
      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)

c..   local variables
      
      integer j,k,l

c..   first exectubable statement
c..   calculation very humble for now. more sophisitication later.

      do l=2,lmax
         do k=kroot,ktip
            do j=jtail1,jtail2
               sn(j,k,l)=sqrt((x(j,k,l)-x(j,k,1))**2
     >              +(y(j,k,l)-y(j,k,1))**2
     >              +(z(j,k,l)-z(j,k,1))**2)
            enddo
         enddo
      enddo
      
      return
      end


c*************************************************************
      subroutine dist_ground(sn,x,y,z)

      use params_global
      implicit none
      
      real sn(jmax,kmax,lmax)
      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)

c..   local variables
      
      integer j,k,l

c..   first exectubable statement
c..   calculation very humble for now. more sophisitication later.

      do l=2,lmax
         do k=kroot,ktip
            do j=jtail1,jtail2
               sn(j,k,l)=abs(z(j,k,l)+zground*rartio)
            enddo
         enddo
      enddo
      
      return
      end


c*************************************************************
      subroutine vortic(vort,u,v,w,xx,xy,xz,yx,yy,yz,zx,zy,zz)
   
      use params_global
      implicit none
      
c*************************************************************
      
      real vort(jmax,kmax,lmax)
      real u(jmax,kmax,lmax),v(jmax,kmax,lmax),w(jmax,kmax,lmax)
      real xx(jmax,kmax,lmax), xy(jmax,kmax,lmax),xz(jmax,kmax,lmax),
     $     yx(jmax,kmax,lmax), yy(jmax,kmax,lmax),yz(jmax,kmax,lmax),
     $     zx(jmax,kmax,lmax), zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)
      
      
c..   local variables
      
      integer j,k,l
      real dx2,dy2,dz2,third
      real ux,vx,wx,uy,vy,wy,uz,vz,wz,tx,ty,tz
      real vort1,dudx,dudy,dudz
      real dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
      real fdiv,stra1
      real s11,s22,s33,s12,s13,s23
      
c..   first exectubable statement
      
      dx2     = 0.5
      dy2     = 0.5
      dz2     = 0.5
      third   = 1./3.
      
      do  l = 2,lmax-1
         do k   = 2,kmax-1
            do  j   = 2,jmax-1
               
               ux  = (u(j+1,k,l) - u(j-1,k,l))*dx2
               vx  = (v(j+1,k,l) - v(j-1,k,l))*dx2
               wx  = (w(j+1,k,l) - w(j-1,k,l))*dx2
               
               uy  = (u(j,k+1,l) - u(j,k-1,l))*dy2
               vy  = (v(j,k+1,l) - v(j,k-1,l))*dy2
               wy  = (w(j,k+1,l) - w(j,k-1,l))*dy2
               
               uz  = (u(j,k,l+1) - u(j,k,l-1))*dz2
               vz  = (v(j,k,l+1) - v(j,k,l-1))*dz2
               wz  = (w(j,k,l+1) - w(j,k,l-1))*dz2

               tx  =  xy(j,k,l)*wx + yy(j,k,l)*wy +zy(j,k,l)*wz
     *               -xz(j,k,l)*vx - yz(j,k,l)*vy -zz(j,k,l)*vz
               ty  =  xz(j,k,l)*ux + yz(j,k,l)*uy +zz(j,k,l)*uz
     *               -xx(j,k,l)*wx - yx(j,k,l)*wy -zx(j,k,l)*wz
               tz  =  xx(j,k,l)*vx + yx(j,k,l)*vy +zx(j,k,l)*vz
     *               -xy(j,k,l)*ux - yy(j,k,l)*uy -zy(j,k,l)*uz
          
               vort1       = sqrt(tx**2 +ty**2 +tz**2)

               dudx        = ux*xx(j,k,l) + uy*yx(j,k,l)
     &                                    + uz*zx(j,k,l)
               dvdy        = vx*xy(j,k,l) + vy*yy(j,k,l)
     &                                    + vz*zy(j,k,l)
               dwdz        = wx*xz(j,k,l) + wy*yz(j,k,l)
     &                                    + wz*zz(j,k,l)

               fdiv        = (dudx + dvdy + dwdz)*third

               dudy        = ux*xy(j,k,l) + uy*yy(j,k,l)
     &                                    + uz*zy(j,k,l)
               dudz        = ux*xz(j,k,l) + uy*yz(j,k,l)
     &                                    + uz*zz(j,k,l)
               dvdx        = vx*xx(j,k,l) + vy*yx(j,k,l)
     &                                    + vz*zx(j,k,l)
               dvdz        = vx*xz(j,k,l) + vy*yz(j,k,l)
     &                                    + vz*zz(j,k,l)
               dwdx        = wx*xx(j,k,l) + wy*yx(j,k,l)
     &                                    + wz*zx(j,k,l)
               dwdy        = wx*xy(j,k,l) + wy*yy(j,k,l)
     &                                    + wz*zy(j,k,l)
               s11         = dudx - fdiv
               s22         = dvdy - fdiv
               s33         = dwdz - fdiv
               s12         = dudy + dvdx
               s13         = dudz + dwdx
               s23         = dwdy + dvdz
               stra1 = sqrt(( 2.0*(s11**2 + s22**2 + s33**2)
     &                           + s12**2 + s13**2 + s23**2 ))
               
              vort(j,k,l)= (vort1+2.*min(0.,stra1-vort1))
c 	      vort(j,k,l)=vort1
            enddo
         enddo
      enddo

c	print*,'NEW PRODUCTION TERM'
      
      return
      end
      
      

c*************************************************************
      subroutine c_fv1(vnu,q,vnul,fv1)
      
      use params_global
      implicit none

      real vnu(jmax,kmax,lmax),q(jmax,kmax,lmax,nd)
      real vnul(jmax,kmax,lmax),fv1(jmax,kmax,lmax)
      
c..   local variables

      real cb1,sigma,cb2,akt,cw1,cw2,cw3,cv1,ct1,ct2,ct3,ct4
      data cb1,sigma,cb2,akt,cw1,cw2,cw3,cv1,ct1,ct2,ct3,ct4/0.1355,
     &  0.6666667,0.622,0.41,3.2390678,0.3,2.,7.1,1.,2.,1.2,0.5/

      real chi
      integer j,k,l


      do  l   = 1,lmax
         do  k   = 1,kmax
            do  j   = 1,jmax
               chi	  =vnu(j,k,l)/vnul(j,k,l)	
               fv1(j,k,l)=chi**3/(chi**3+cv1**3)
            enddo
         enddo
      enddo
            
            
      return
      end


c*************************************************************
      subroutine c_turm(fv1,turmu,vnu,q)
      
      use params_global
      implicit none
      
      real vnu(jmax,kmax,lmax),q(jmax,kmax,lmax,nd)
      real turmu(jmax,kmax,lmax),fv1(jmax,kmax,lmax)
      
      integer j,k,l
c*************************************************************
      do  l   = 1,lmax
         do  k   = 1,kmax
            do   j   = 1,jmax
               turmu(j,k,l)=q(j,k,l,1)*q(j,k,l,6)*vnu(j,k,l)*fv1(j,k,l)
            enddo
         enddo
      enddo
            
      return
      end

c***********************************************************************
      subroutine rhstur(q,fv1,sn,vort,u,v,w,ug,vg,wg,vnu,vnu0,vnul,
     & trip,s,xx,xy,xz,yx,yy,yz,zx,zy,zz,x,y,z,sjmat,dl2norm,
     & tscale)
c***********************************************************************

      use params_global
      implicit none


      real q(jmax,kmax,lmax,nd),tscale(jmax,kmax,lmax)
      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)
      real xx(jmax,kmax,lmax), xy(jmax,kmax,lmax),xz(jmax,kmax,lmax),
     $     yx(jmax,kmax,lmax), yy(jmax,kmax,lmax),yz(jmax,kmax,lmax),
     $     zx(jmax,kmax,lmax), zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax), vg(jmax,kmax,lmax),wg(jmax,kmax,lmax)
      real u(jmax,kmax,lmax),v(jmax,kmax,lmax),w(jmax,kmax,lmax)
      real vnul(jmax,kmax,lmax),vort(jmax,kmax,lmax),sn(jmax,kmax,lmax)
      real s(jmax,kmax,lmax),vnu(jmax,kmax,lmax),vnu0(jmax,kmax,lmax)
      real fv1(jmax,kmax,lmax),trip(jmax,kmax,lmax)
      real sjmat(jmax,kmax,lmax)
      real dl2norm
            
c..   local variables

      integer l,k,j
	
      do l = 1,lmax
         do k = 1,kmax
            do j = 1,jmax
               s(j,k,l)=0.0
               sjmat(j,k,l)=0.0
            enddo
         enddo
      enddo
      

      call convec(q,vnu,u,v,w,ug,vg,wg,s,
     &     xx,xy,xz,yx,yy,yz,zx,zy,zz,sjmat)

      
      call diffus(q,vnu,u,v,w,s,vnul,
     &     xx,xy,xz,yx,yy,yz,zx,zy,zz,sjmat)
	
      
      call sourc(q,vnu,u,v,w,s,vnul,vort,sn,
     &     xx,xy,xz,yx,yy,yz,zx,zy,zz,sjmat)

      if(itnmax.gt.1) then
       do l = 1,lmax
         do k = 1,kmax
            do j = 1,jmax
               s(j,k,l)=s(j,k,l)+(vnu0(j,k,l)-vnu(j,k,l))/h
            enddo
         enddo
       enddo
      endif
     

      dl2norm=0.0
      do l = 1,lmax
         do k = 1,kmax
            do j = 1,jmax
               dl2norm    =dl2norm+s(j,k,l)*s(j,k,l)
               s(j,k,l)  = s(j,k,l)*tscale(j,k,l)
               sjmat(j,k,l)= 1.+sjmat(j,k,l)*tscale(j,k,l)
            enddo
         enddo
      enddo
      
      return
      end


c***********************************************************************
      subroutine lhstur(q,fv1,sn,vort,u,v,w,ug,vg,wg,vnu,vnul,trip,s,
     &     xx,xy,xz,yx,yy,yz,zx,zy,zz,x,y,z,sjmat,amat,bmat,cmat,tscale)
c***********************************************************************
      use params_global
      implicit none

      real q(jmax,kmax,lmax,nd),tscale(jmax,kmax,lmax)
      real fv1(jmax,kmax,lmax),sn(jmax,kmax,lmax),vort(jmax,kmax,lmax)
      real u(jmax,kmax,lmax),v(jmax,kmax,lmax),w(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax),vg(jmax,kmax,lmax),wg(jmax,kmax,lmax)
      real vnu(jmax,kmax,lmax),vnul(jmax,kmax,lmax),trip(jmax,kmax,lmax)
      real s(jmax,kmax,lmax)
      real xx(jmax,kmax,lmax), xy(jmax,kmax,lmax),xz(jmax,kmax,lmax),
     $     yx(jmax,kmax,lmax), yy(jmax,kmax,lmax),yz(jmax,kmax,lmax),
     $     zx(jmax,kmax,lmax), zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)
      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)

      real sjmat(jmax,kmax,lmax),amat(mdim,mdim),bmat(mdim,mdim)
      real cmat(mdim,mdim)

c..   local variables

      real, allocatable :: up(:),um(:),xmh(:),ymh(:),zmh(:),chp(:)
      real fwd,bck
      real cb1,sigma,cb2,akt,cw1,cw2,cw3,cv1,ct1,ct2,ct3,ct4,eps
      data cb1,sigma,cb2,akt,cw1,cw2,cw3,cv1,ct1,ct2,ct3,ct4,eps/0.1355,
     &  0.6666667,0.622,0.41,3.2390678,0.3,2.,7.1,1.,2.,1.2,0.5,0.0/

      integer j,k,l,jp,kp,lp,jm1,km1,lm1
      
      real sigreyinv,vnulh,vnuh,dmetp,dmetm
      real c2,dcp,dcm,uu
     
      allocate(up(mdim),um(mdim),xmh(mdim),ymh(mdim),
     <     zmh(mdim),chp(mdim))
 
      sigreyinv=1./sigma/REY


	eps=0.
	do 30 l=2,lmax-1


c*******************************
c	   j-direction
c*******************************
	

	   do k=1,kmax
	   do j=1,jmax
	      amat(j,k)=0.0	
	      bmat(j,k)=0.0	
	      cmat(j,k)=0.0	
	   enddo
	   enddo

	   do k=2,kmax-1
	     do j=1,jmax-1
		uu=xx(j,k,l)*(u(j,k,l)-ug(j,k,l))+xy(j,k,l)*(v(j,k,l)-vg(j,k,l))
     >	  			                 +xz(j,k,l)*(w(j,k,l)-wg(j,k,l))
		up(j)=0.5*(uu+abs(uu))
		um(j)=0.5*(uu-abs(uu))
		xmh(j) =0.5*(  xx(j,k,l)+  xx(j+1,k,l))
		ymh(j) =0.5*(  xy(j,k,l)+  xy(j+1,k,l))
		zmh(j) =0.5*(  xz(j,k,l)+  xz(j+1,k,l))
		vnulh  =0.5*(vnul(j,k,l)+vnul(j+1,k,l))
		vnuh   =0.5*( vnu(j,k,l)+ vnu(j+1,k,l))
		chp(j) =(1.+cb2)*sigreyinv*(vnulh+vnuh)	
	     enddo
		j=jmax
		uu=xx(j,k,l)*(u(j,k,l)-ug(j,k,l))+xy(j,k,l)*(v(j,k,l)-vg(j,k,l))
     >	  			                 +xz(j,k,l)*(w(j,k,l)-wg(j,k,l))
		up(j)=0.5*(uu+abs(uu))
		um(j)=0.5*(uu-abs(uu))

	     do j=2,jmax-1
	        jp=j+1
	        jm=j-1
c	        convection contribution
		if(up(jp).gt.eps) then
	 	  fwd=1.
 		else
	          fwd=0.
		endif
		if(um(jm).lt.-eps) then
	 	  bck=1.
 		else
	          bck=0.
		endif
		amat(jp,k)=amat(jp,k)-fwd*(up(j)+um(j))
c		bmat(j ,k)=bmat(j ,k)+(up(j)-um(j))
		cmat(jm,k)=cmat(jm,k)+bck*(up(j)+um(j))


	        dmetp=xmh(j )*xx(j,k,l)+ymh(j )*xy(j,k,l)+zmh(j )*xz(j,k,l)
	        dmetm=xmh(jm)*xx(j,k,l)+ymh(jm)*xy(j,k,l)+zmh(jm)*xz(j,k,l)

	        c2  =cb2*sigreyinv
	        c2  =c2*(vnul(j,k,l)+vnu(j,k,l))
c	 limit
	        dcp = max(dmetp*(chp(j )-c2),0.)
	        dcm = max(dmetm*(chp(jm)-c2),0.)
c		diffusion contribution
		amat(j,k) = amat(j,k)-dcm
c		bmat(j,k) = bmat(j,k)+dcm+dcp
		cmat(j,k) = cmat(j,k)-dcp
	     enddo
	   enddo

	   do k = 1,kmax
             do j = 1,jmax
               amat(j,k) =    amat(j,k)*tscale(j,k,l)
               bmat(j,k) =    sjmat(j,k,l)
               cmat(j,k) =    cmat(j,k)*tscale(j,k,l)
             enddo
           enddo

c	   call lsolvej(amat,bmat,cmat,s,jmax-1,kmax-1,l,mdim)
	   call lsolvej(amat,bmat,cmat,s,jmax-1,kmax-1,l,jmax,kmax,lmax,mdim)

c***********************************************
c	   k-direction 
c***********************************************

	   do k=1,kmax
	   do j=1,jmax
	      amat(j,k)=0.0	
	      bmat(j,k)=0.0	
	      cmat(j,k)=0.0	
	      s(j,k,l) =sjmat(j,k,l)*s(j,k,l)
	   enddo
	   enddo

	   do j=2,jmax-1
	     do k=1,kmax-1
		uu=yx(j,k,l)*(u(j,k,l)-ug(j,k,l))+yy(j,k,l)*(v(j,k,l)-vg(j,k,l))
     >	  			                 +yz(j,k,l)*(w(j,k,l)-wg(j,k,l))
		up(k)=0.5*(uu+abs(uu))
		um(k)=0.5*(uu-abs(uu))
		xmh(k) =0.5*(  yx(j,k,l)+  yx(j,k+1,l))
		ymh(k) =0.5*(  yy(j,k,l)+  yy(j,k+1,l))
		zmh(k) =0.5*(  yz(j,k,l)+  yz(j,k+1,l))
		vnulh  =0.5*(vnul(j,k,l)+vnul(j,k+1,l))
		vnuh   =0.5*( vnu(j,k,l)+ vnu(j,k+1,l))
		chp(k) =(1.+cb2)*sigreyinv*(vnulh+vnuh)	
	     enddo
		k=kmax
		uu=yx(j,k,l)*(u(j,k,l)-ug(j,k,l))+yy(j,k,l)*(v(j,k,l)-vg(j,k,l))
     >	  			                 +yz(j,k,l)*(w(j,k,l)-wg(j,k,l))
		up(k)=0.5*(uu+abs(uu))
		um(k)=0.5*(uu-abs(uu))

	     do k=2,kmax-1
	        kp=k+1
	        km=k-1
c	        convection contribution
		if(up(kp).gt.eps) then
	 	  fwd=1.
 		else
	          fwd=0.
		endif
		if(um(km).lt.-eps) then
	 	  bck=1.
 		else
	          bck=0.
		endif

		amat(j,kp)=amat(j,kp)-fwd*(up(k)+um(k))
c		bmat(j,k )=bmat(j,k )+up(k)-um(k)
		cmat(j,km)=cmat(j,km)+bck*(up(k)+um(k))


	        dmetp=xmh(k )*yx(j,k,l)+ymh(k )*yy(j,k,l)+zmh(k )*yz(j,k,l)
	        dmetm=xmh(km)*yx(j,k,l)+ymh(km)*yy(j,k,l)+zmh(km)*yz(j,k,l)

	        c2  =cb2*sigreyinv
	        c2  =c2*(vnul(j,k,l)+vnu(j,k,l))
c	 limit
	        dcp = max(dmetp*(chp(k )-c2),0.)
	        dcm = max(dmetm*(chp(km)-c2),0.)
c		diffusion contribution
		amat(j,k) = amat(j,k)-dcm
c		bmat(j,k) = bmat(j,k)+dcm+dcp
		cmat(j,k) = cmat(j,k)-dcp
	     enddo
	   enddo

	   do k = 1,kmax
             do j = 1,jmax
               amat(j,k) =    amat(j,k)*tscale(j,k,l)
               bmat(j,k) =    sjmat(j,k,l)
               cmat(j,k) =    cmat(j,k)*tscale(j,k,l)
             enddo
           enddo

c	  call lsolvek(amat,bmat,cmat,s,jmax-1,kmax-1,l,mdim)
          call lsolvek(amat,bmat,cmat,s,jmax-1,kmax-1,l,
     <         jmax,kmax,lmax,mdim)

	   do k = 1,kmax
             do j = 1,jmax
	      s(j,k,l) =sjmat(j,k,l)*s(j,k,l)
	     enddo
	   enddo
30	continue





c*******************************
c	   l-direction
c*******************************

	do 20 k=2,kmax-1

	   do l=1,lmax
	   do j=1,jmax
	      amat(j,l)=0.0	
	      bmat(j,l)=0.0	
	      cmat(j,l)=0.0	
	   enddo
	   enddo

	   do j=2,jmax-1
	     do l=1,lmax-1
		uu=zx(j,k,l)*(u(j,k,l)-ug(j,k,l))+zy(j,k,l)*(v(j,k,l)-vg(j,k,l))
     >	  			                 +zz(j,k,l)*(w(j,k,l)-wg(j,k,l))
		up(l)=0.5*(uu+abs(uu))
		um(l)=0.5*(uu-abs(uu))
		xmh(l) =0.5*(  zx(j,k,l)+  zx(j,k,l+1))
		ymh(l) =0.5*(  zy(j,k,l)+  zy(j,k,l+1))
		zmh(l) =0.5*(  zz(j,k,l)+  zz(j,k,l+1))
		vnulh  =0.5*(vnul(j,k,l)+vnul(j,k,l+1))
		vnuh   =0.5*( vnu(j,k,l)+ vnu(j,k,l+1))
		chp(l) =(1.+cb2)*sigreyinv*(vnulh+vnuh)	
	     enddo
		l=lmax
		uu=zx(j,k,l)*(u(j,k,l)-ug(j,k,l))+zy(j,k,l)*(v(j,k,l)-vg(j,k,l))
     >	  			                 +zz(j,k,l)*(w(j,k,l)-wg(j,k,l))
		up(l)=0.5*(uu+abs(uu))
		um(l)=0.5*(uu-abs(uu))

	     do l=2,lmax-1
	        lp=l+1
	        lm=l-1
c		need to make this better later.
		if(up(lp).gt.eps) then
	 	  fwd=1.
 		else
	          fwd=0.
		endif
		if(um(lm).lt.-eps) then
	 	  bck=1.
 		else
	          bck=0.
		endif
		amat(j,lp)=amat(j,lp)-fwd*(up(l)+um(l))
c		bmat(j,l )=bmat(j,l )+up(l)-um(l)
		cmat(j,lm)=cmat(j,lm)+bck*(up(l)+um(l))

	        dmetp=xmh(l )*zx(j,k,l)+ymh(l )*zy(j,k,l)+zmh(l )*zz(j,k,l)
	        dmetm=xmh(lm)*zx(j,k,l)+ymh(lm)*zy(j,k,l)+zmh(lm)*zz(j,k,l)

	        c2  =cb2*sigreyinv
	        c2  =c2*(vnul(j,k,l)+vnu(j,k,l))
c	 limit
	        dcp = max(dmetp*(chp(l )-c2),0.)
	        dcm = max(dmetm*(chp(lm)-c2),0.)
c		diffusion contribution
		amat(j,l) = amat(j,l)-dcm
c		bmat(j,l) = bmat(j,l)+dcm+dcp
		cmat(j,l) = cmat(j,l)-dcp
	     enddo
	   enddo

	   do l = 1,lmax
             do j = 1,jmax
               amat(j,l) =    amat(j,l)*tscale(j,k,l)
               bmat(j,l) =    sjmat(j,k,l)
               cmat(j,l) =    cmat(j,l)*tscale(j,k,l)
             enddo
           enddo

c	   call lsolvel(amat,bmat,cmat,s,jmax-1,k,lmax-1,mdim)
           call lsolvel(amat,bmat,cmat,s,jmax-1,k,lmax-1,
     <          jmax,kmax,lmax,mdim)
                        
20	continue

      lm=lmax-1
      km=kmax-1
      jm=jmax-1
      return
      end

c***********************************************************************
      subroutine lhstur1(q,fv1,sn,vort,u,v,w,ug,vg,wg,vnu,vnul,trip,s,
     &     xx,xy,xz,yx,yy,yz,zx,zy,zz,x,y,z,sjmat,amat,bmat,cmat,tscale)
c***********************************************************************
      use params_global
      implicit none

      real q(jmax,kmax,lmax,nd),tscale(jmax,kmax,lmax)
      real fv1(jmax,kmax,lmax),sn(jmax,kmax,lmax),vort(jmax,kmax,lmax)
      real u(jmax,kmax,lmax),v(jmax,kmax,lmax),w(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax),vg(jmax,kmax,lmax),wg(jmax,kmax,lmax)
      real vnu(jmax,kmax,lmax),vnul(jmax,kmax,lmax),trip(jmax,kmax,lmax)
      real s(jmax,kmax,lmax)
      real xx(jmax,kmax,lmax), xy(jmax,kmax,lmax),xz(jmax,kmax,lmax),
     $     yx(jmax,kmax,lmax), yy(jmax,kmax,lmax),yz(jmax,kmax,lmax),
     $     zx(jmax,kmax,lmax), zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)
      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)

      real sjmat(jmax,kmax,lmax),amat(mdim,mdim),bmat(mdim,mdim)
      real cmat(mdim,mdim)

c..   local variables

      real, allocatable :: up(:),um(:),xmh(:),ymh(:),zmh(:),chp(:)
      real fwd,bck
      real cb1,sigma,cb2,akt,cw1,cw2,cw3,cv1,ct1,ct2,ct3,ct4,eps
      data cb1,sigma,cb2,akt,cw1,cw2,cw3,cv1,ct1,ct2,ct3,ct4,eps/0.1355,
     &  0.6666667,0.622,0.41,3.2390678,0.3,2.,7.1,1.,2.,1.2,0.5,0.0/

      integer j,k,l,jp,kp,lp,jm1,km1,lm1
      
      real sigreyinv,vnulh,vnuh,metp,metm
      real c2,dcp,dcm,uu
     
      allocate(up(mdim),um(mdim),xmh(mdim),ymh(mdim),
     <     zmh(mdim),chp(mdim))
 
      sigreyinv=1./sigma/REY

c*******************************
c..   j-direction
c*******************************

      do 10 l=2,lmax-1
         
         do k=1,kmax
            do j=1,jmax
               amat(j,k)=0.0	
               cmat(j,k)=0.0	
            enddo
         enddo
         
         do k=2,kmax-1
            do j=1,jmax-1
               uu  =xx(j,k,l)*(u(j,k,l)-ug(j,k,l))
     >              +xy(j,k,l)*(v(j,k,l)-vg(j,k,l))
     >              +xz(j,k,l)*(w(j,k,l)-wg(j,k,l))

               up(j)=0.5*(uu+abs(uu))
               um(j)=0.5*(uu-abs(uu))
               xmh(j) =0.5*(  xx(j,k,l)+  xx(j+1,k,l))
               ymh(j) =0.5*(  xy(j,k,l)+  xy(j+1,k,l))
               zmh(j) =0.5*(  xz(j,k,l)+  xz(j+1,k,l))
               vnulh  =0.5*(vnul(j,k,l)+vnul(j+1,k,l))
               vnuh   =0.5*( vnu(j,k,l)+ vnu(j+1,k,l))
               chp(j) =(1.+cb2)*sigreyinv*(vnulh+vnuh)	
            enddo
            j=jmax
            uu  =xx(j,k,l)*(u(j,k,l)-ug(j,k,l))
     >           +xy(j,k,l)*(v(j,k,l)-vg(j,k,l))
     >           +xz(j,k,l)*(w(j,k,l)-wg(j,k,l))
               up(j)=0.5*(uu+abs(uu))
               um(j)=0.5*(uu-abs(uu))
            
            do j=2,jmax-1
               jp=j+1
               jm1=j-1
c..   convection contribution
		if(up(jp).gt.EPS) then
                  fwd=1.
                else
                  fwd=0.
                endif
                if(um(jm1).lt.-EPS) then
                  bck=1.
                else
                  bck=0.
                endif
		amat(jp,k)=amat(jp,k)-fwd*up(j)
		cmat(jm1,k)=cmat(jm1,k)+bck*um(j)

	        metp=xmh(j )*xx(j,k,l)+ymh(j )*xy(j,k,l)+zmh(j )*xz(j,k,l)
	        metm=xmh(jm1)*xx(j,k,l)+ymh(jm1)*xy(j,k,l)+zmh(jm1)*xz(j,k,l)

	        c2  =cb2*sigreyinv
	        c2  =c2*(vnul(j,k,l)+vnu(j,k,l))
c..   limit
	        dcp = max(metp*(chp(j )-c2),0.)
	        dcm = max(metm*(chp(jm1)-c2),0.)
c..   diffusion contribution
		amat(j,k) = amat(j,k)-dcm
		cmat(j,k) = cmat(j,k)-dcp
	     enddo
          enddo
          
	   do k = 1,kmax
              do j = 1,jmax
                 amat(j,k) =    amat(j,k)*tscale(j,k,l)
                 bmat(j,k) = 	sjmat(j,k,l)
                 cmat(j,k) =    cmat(j,k)*tscale(j,k,l)
              enddo
           enddo
           
	   call lsolvej(amat,bmat,cmat,s,jmax-1,kmax-1,l,jmax,kmax,lmax,mdim)


             
             do k=1,kmax
                do j=1,jmax
                   amat(j,k)=0.0	
                   bmat(j,k)=0.0	
                   cmat(j,k)=0.0	
		   s(j,k,l) =sjmat(j,k,l)*s(j,k,l)
                enddo
             enddo
             
             do j=2,jmax-1
                do k=1,kmax-1
                   uu=yx(j,k,l)*(u(j,k,l)-ug(j,k,l))
     >                  +yy(j,k,l)*(v(j,k,l)-vg(j,k,l))
     >                  +yz(j,k,l)*(w(j,k,l)-wg(j,k,l))
               up(k)=0.5*(uu+abs(uu))
               um(k)=0.5*(uu-abs(uu))
                   xmh(k) =0.5*(  yx(j,k,l)+  yx(j,k+1,l))
                   ymh(k) =0.5*(  yy(j,k,l)+  yy(j,k+1,l))
                   zmh(k) =0.5*(  yz(j,k,l)+  yz(j,k+1,l))
                   vnulh  =0.5*(vnul(j,k,l)+vnul(j,k+1,l))
                   vnuh   =0.5*( vnu(j,k,l)+ vnu(j,k+1,l))
                   chp(k) =(1.+cb2)*sigreyinv*(vnulh+vnuh)	
                enddo
		k=kmax
		uu=yx(j,k,l)*(u(j,k,l)-ug(j,k,l))
     >               +yy(j,k,l)*(v(j,k,l)-vg(j,k,l))
     >               +yz(j,k,l)*(w(j,k,l)-wg(j,k,l))
               up(k)=0.5*(uu+abs(uu))
               um(k)=0.5*(uu-abs(uu))
                
                do k=2,kmax-1
                   kp=k+1
                   km1=k-1
c..   convection contribution
		if(up(kp).gt.EPS) then
                  fwd=1.
                else
                  fwd=0.
                endif
                if(um(km1).lt.-EPS) then
                  bck=1.
                else
                  bck=0.
                endif
                   amat(j,kp)=amat(j,kp)-fwd*up(k)
                   cmat(j,km1)=cmat(j,km1)+bck*um(k)
                   
                   metp=xmh(k )*yx(j,k,l)+ymh(k )*yy(j,k,l)
     >                  +zmh(k )*yz(j,k,l)
                   metm=xmh(km1)*yx(j,k,l)+ymh(km1)*yy(j,k,l)
     >                  +zmh(km1)*yz(j,k,l)
                   
                   c2  =cb2*sigreyinv
                   c2  =c2*(vnul(j,k,l)+vnu(j,k,l))
c..   limit
                   dcp = max(metp*(chp(k )-c2),0.)
                   dcm = max(metm*(chp(km1)-c2),0.)
c..   diffusion contribution
                   amat(j,k) = amat(j,k)-dcm
                   cmat(j,k) = cmat(j,k)-dcp
	     enddo
          enddo

          do k = 1,kmax
             do j = 1,jmax
                amat(j,k) =    amat(j,k)*tscale(j,k,l)
                bmat(j,k) =    sjmat(j,k,l)
                cmat(j,k) =    cmat(j,k)*tscale(j,k,l)
             enddo
          enddo
           
          call lsolvek(amat,bmat,cmat,s,jmax-1,kmax-1,l,
     <         jmax,kmax,lmax,mdim)

	    do k = 1,kmax
             do j = 1,jmax
              s(j,k,l) =sjmat(j,k,l)*s(j,k,l)
             enddo
           enddo

10	continue




c*******************************
c..   l-direction
c*******************************

	do 20 k=2,kmax-1
           
	   do l=1,lmax
              do j=1,jmax
                 amat(j,l)=0.0	
                 bmat(j,l)=0.0
                 cmat(j,l)=0.0	
              enddo
	   enddo

	   do j=2,jmax-1
              do l=1,lmax-1
                 uu=zx(j,k,l)*(u(j,k,l)-ug(j,k,l))
     >                +zy(j,k,l)*(v(j,k,l)-vg(j,k,l))
     >                +zz(j,k,l)*(w(j,k,l)-wg(j,k,l))
               up(l)=0.5*(uu+abs(uu))
               um(l)=0.5*(uu-abs(uu))
		xmh(l) =0.5*(  zx(j,k,l)+  zx(j,k,l+1))
		ymh(l) =0.5*(  zy(j,k,l)+  zy(j,k,l+1))
		zmh(l) =0.5*(  zz(j,k,l)+  zz(j,k,l+1))
		vnulh  =0.5*(vnul(j,k,l)+vnul(j,k,l+1))
		vnuh   =0.5*( vnu(j,k,l)+ vnu(j,k,l+1))
		chp(l) =(1.+cb2)*sigreyinv*(vnulh+vnuh)	
	     enddo
		l=lmax
		uu=zx(j,k,l)*(u(j,k,l)-ug(j,k,l))
     >               +zy(j,k,l)*(v(j,k,l)-vg(j,k,l))
     >               +zz(j,k,l)*(w(j,k,l)-wg(j,k,l))

               up(l)=0.5*(uu+abs(uu))
               um(l)=0.5*(uu-abs(uu))

                do l=2,lmax-1
                   lp=l+1
                   lm1=l-1
c..   need to make this better later.
	        if(up(lp).gt.EPS) then
                  fwd=1.
                else
                  fwd=0.
                endif
                if(um(lm1).lt.-EPS) then
                  bck=1.
                else
                  bck=0.
                endif

                   amat(j,lp)=amat(j,lp)-fwd*up(l)
                   cmat(j,lm1)=cmat(j,lm1)+bck*um(l)
                   
                   metp=xmh(l )*zx(j,k,l)+ymh(l )*zy(j,k,l)
     >                  +zmh(l )*zz(j,k,l)
                   metm=xmh(lm1)*zx(j,k,l)+ymh(lm1)*zy(j,k,l)
     >                  +zmh(lm1)*zz(j,k,l)
                   
                   c2  =cb2*sigreyinv
                   c2  =c2*(vnul(j,k,l)+vnu(j,k,l))
c..   limit
                   dcp = max(metp*(chp(l )-c2),0.)
                   dcm = max(metm*(chp(lm1)-c2),0.)
c..   diffusion contribution
                   amat(j,l) = amat(j,l)-dcm
                   cmat(j,l) = cmat(j,l)-dcp
                enddo
             enddo
             
             do l = 1,lmax
                do j = 1,jmax
                   amat(j,l) =    amat(j,l)*tscale(j,k,l)
   		   bmat(j,l) =    sjmat(j,k,l)
                   cmat(j,l) =    cmat(j,l)*tscale(j,k,l)
                enddo
             enddo
             
             call lsolvel(amat,bmat,cmat,s,jmax-1,k,lmax-1,
     <            jmax,kmax,lmax,mdim)

 20       continue
          
       return
       end



c***********************************************************************
      subroutine convec(q,vnu,u,v,w,ug,vg,wg,s,
     &  xx,xy,xz,yx,yy,yz,zx,zy,zz,sjmat)
c***********************************************************************
      use params_global
      implicit none
      
      real q(jmax,kmax,lmax,nd),sjmat(jmax,kmax,lmax)
      real xx(jmax,kmax,lmax), xy(jmax,kmax,lmax),xz(jmax,kmax,lmax),
     $     yx(jmax,kmax,lmax), yy(jmax,kmax,lmax),yz(jmax,kmax,lmax),
     $     zx(jmax,kmax,lmax), zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax), vg(jmax,kmax,lmax),wg(jmax,kmax,lmax)
      real u(jmax,kmax,lmax),v(jmax,kmax,lmax),w(jmax,kmax,lmax)
      real s(jmax,kmax,lmax),vnu(jmax,kmax,lmax)

c..   local variables

      integer j,k,l
      real uu,vv,ww
      real up,um,vp,vm,wp,wm


      do 10 l=2,lmax-1
      do 10 k=2,kmax-1
      do 10 j=2,jmax-1
        uu=xx(j,k,l)*(u(j,k,l)-ug(j,k,l))+xy(j,k,l)*(v(j,k,l)-vg(j,k,l))
     >	  			         +xz(j,k,l)*(w(j,k,l)-wg(j,k,l))
	vv=yx(j,k,l)*(u(j,k,l)-ug(j,k,l))+yy(j,k,l)*(v(j,k,l)-vg(j,k,l))
     >	  			         +yz(j,k,l)*(w(j,k,l)-wg(j,k,l))
	ww=zx(j,k,l)*(u(j,k,l)-ug(j,k,l))+zy(j,k,l)*(v(j,k,l)-vg(j,k,l))
     >	  			         +zz(j,k,l)*(w(j,k,l)-wg(j,k,l))
	up=0.5*(uu+abs(uu))
	um=0.5*(uu-abs(uu))
	vp=0.5*(vv+abs(vv))
	vm=0.5*(vv-abs(vv))
	wp=0.5*(ww+abs(ww))
	wm=0.5*(ww-abs(ww))

	s(j,k,l)=s(j,k,l)-up*(vnu(j  ,k  ,l  )-vnu(j-1,k  ,l  ))
     <		         -um*(vnu(j+1,k  ,l  )-vnu(j  ,k  ,l  ))
	s(j,k,l)=s(j,k,l)-vp*(vnu(j  ,k  ,l  )-vnu(j  ,k-1,l  ))
     <		         -vm*(vnu(j  ,k+1,l  )-vnu(j  ,k  ,l  ))
	s(j,k,l)=s(j,k,l)-wp*(vnu(j  ,k  ,l  )-vnu(j  ,k  ,l-1))
     <		         -wm*(vnu(j  ,k  ,l+1)-vnu(j  ,k  ,l  ))
	SJMAT(J,K,L)=SJMAT(J,K,L)+UP-UM
        SJMAT(J,K,L)=SJMAT(J,K,L)+VP-VM
        SJMAT(J,K,L)=SJMAT(J,K,L)+WP-WM
10	continue

        return
        end


c***********************************************************************
      subroutine diffus(q,vnu,u,v,w,s,vnul,
     &     xx,xy,xz,yx,yy,yz,zx,zy,zz,sjmat)
c***********************************************************************
      use params_global
      implicit none

      real q(jmax,kmax,lmax,nd)
      real xx(jmax,kmax,lmax), xy(jmax,kmax,lmax),xz(jmax,kmax,lmax),
     $     yx(jmax,kmax,lmax), yy(jmax,kmax,lmax),yz(jmax,kmax,lmax),
     $     zx(jmax,kmax,lmax), zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)
      real u(jmax,kmax,lmax),v(jmax,kmax,lmax),w(jmax,kmax,lmax)
      real s(jmax,kmax,lmax),sjmat(jmax,kmax,lmax)
      
      real vnu(jmax,kmax,lmax),vnul(jmax,kmax,lmax)


c..   local variables

      real, allocatable :: xmh(:),ymh(:),zmh(:),chp(:),dnuhp(:)
      real cb1,sigma,cb2,akt,cw1,cw2,cw3,cv1,ct1,ct2,ct3,ct4
      data cb1,sigma,cb2,akt,cw1,cw2,cw3,cv1,ct1,ct2,ct3,ct4/0.1355,
     &  0.6666667,0.622,0.41,3.2390678,0.3,2.,7.1,1.,2.,1.2,0.5/
      integer j,k,l
      real metp,metm,c2,dcp,dcm,sigreyinv,vnuh,vnulh

c..   first executable statement

      allocate(xmh(mdim),ymh(mdim),zmh(mdim),chp(mdim),dnuhp(mdim))
      sigreyinv=1./sigma/REY
      
c.. j fluxes

      do 10 l=2,lmax-1
         do 10 k=2,kmax-1
            
            do j=1,jmax-1
		xmh(j) =0.5*(  xx(j,k,l)+  xx(j+1,k,l))
		ymh(j) =0.5*(  xy(j,k,l)+  xy(j+1,k,l))
		zmh(j) =0.5*(  xz(j,k,l)+  xz(j+1,k,l))
		vnulh  =0.5*(vnul(j,k,l)+vnul(j+1,k,l))
		vnuh   =0.5*( vnu(j,k,l)+ vnu(j+1,k,l))
		chp(j) =(1.+cb2)*sigreyinv*(vnulh+vnuh)	
	      dnuhp(j) =vnu(j+1,k,l)-vnu(j,k,l)
           enddo

           do j=2,jmax-1
              
           metp=xmh(j  )*xx(j,k,l)+ymh(j  )*xy(j,k,l)+zmh(j  )*xz(j,k,l)
           metm=xmh(j-1)*xx(j,k,l)+ymh(j-1)*xy(j,k,l)+zmh(j-1)*xz(j,k,l)
              
              c2  =cb2*sigreyinv
              c2  =c2*(vnul(j,k,l)+vnu(j,k,l))
c..  limit
              dcp = max(metp*(chp(j)  -c2),0.)
              dcm = max(metm*(chp(j-1)-c2),0.)
              
              s(j,k,l)=s(j,k,l)+dcp*dnuhp(j)-dcm*dnuhp(j-1)
          sjmat(j,k,l)=sjmat(j,k,l)+dcp+dcm
           enddo

10	continue

c..  k fluxes
	do 20 l=2,lmax-1
           do 20 j=2,jmax-1

              do k=1,kmax-1
                 xmh(k) =0.5*(  yx(j,k,l)+  yx(j,k+1,l))
                 ymh(k) =0.5*(  yy(j,k,l)+  yy(j,k+1,l))
                 zmh(k) =0.5*(  yz(j,k,l)+  yz(j,k+1,l))
                 vnulh  =0.5*(vnul(j,k,l)+vnul(j,k+1,l))
                 vnuh   =0.5*( vnu(j,k,l)+ vnu(j,k+1,l))
                 chp(k) =(1.+cb2)*sigreyinv*(vnulh+vnuh)	
                 dnuhp(k) =vnu(j,k+1,l)-vnu(j,k,l)
              enddo

              do k=2,kmax-1
                 
                 metp=xmh(k  )*yx(j,k,l)+ymh(k  )*yy(j,k,l)
     <                                  +zmh(k  )*yz(j,k,l)
                 metm=xmh(k-1)*yx(j,k,l)+ymh(k-1)*yy(j,k,l)
     <                                  +zmh(k-1)*yz(j,k,l)

                 c2  =cb2*sigreyinv
                 c2  =c2*(vnul(j,k,l)+vnu(j,k,l))
c..  limit
                 dcp = max(metp*(chp(k)  -c2),0.)
                 dcm = max(metm*(chp(k-1)-c2),0.)

                 s(j,k,l)=s(j,k,l)+dcp*dnuhp(k)-dcm*dnuhp(k-1)
                 sjmat(j,k,l)=sjmat(j,k,l)+dcp+dcm
              enddo
              
 20        continue

c.. l fluxes
           do 30 k=2,kmax-1
              do 30 j=2,jmax-1
                 
                 do l=1,lmax-1
                    xmh(l) =0.5*(  zx(j,k,l)+  zx(j,k,l+1))
                    ymh(l) =0.5*(  zy(j,k,l)+  zy(j,k,l+1))
                    zmh(l) =0.5*(  zz(j,k,l)+  zz(j,k,l+1))
                    vnulh  =0.5*(vnul(j,k,l)+vnul(j,k,l+1))
                    vnuh   =0.5*( vnu(j,k,l)+ vnu(j,k,l+1))
                    chp(l) =(1.+cb2)*sigreyinv*(vnulh+vnuh)	
                    dnuhp(l) =vnu(j,k,l+1)-vnu(j,k,l)
                 enddo

                 do l=2,lmax-1
                    
                    metp=xmh(l  )*zx(j,k,l)+ymh(l  )*zy(j,k,l)
     <                                     +zmh(l  )*zz(j,k,l)
                    metm=xmh(l-1)*zx(j,k,l)+ymh(l-1)*zy(j,k,l)
     <                                     +zmh(l-1)*zz(j,k,l)

                    c2  =cb2*sigreyinv
                    c2  =c2*(vnul(j,k,l)+vnu(j,k,l))
c..  limit
                    dcp = max(metp*(chp(l)  -c2),0.)
                    dcm = max(metm*(chp(l-1)-c2),0.)
                    
                    s(j,k,l)=s(j,k,l)+dcp*dnuhp(l)-dcm*dnuhp(l-1)
                    sjmat(j,k,l)=sjmat(j,k,l)+dcp+dcm
                 enddo
                 
 30   continue
              
      return	
      end


c***********************************************************************
      subroutine sourc(q,vnu,u,v,w,s,vnul,vort,sn,
     &     xx,xy,xz,yx,yy,yz,zx,zy,zz,sjmat)
c***********************************************************************
      use params_global
      implicit none


      real q(jmax,kmax,lmax,nd)
      real u(jmax,kmax,lmax),v(jmax,kmax,lmax),w(jmax,kmax,lmax)
      real s(jmax,kmax,lmax)
      real vnu(jmax,kmax,lmax),vnul(jmax,kmax,lmax),sn(jmax,kmax,lmax)
      real vort(jmax,kmax,lmax),sjmat(jmax,kmax,lmax)
      real xx(jmax,kmax,lmax), xy(jmax,kmax,lmax),xz(jmax,kmax,lmax),
     $     yx(jmax,kmax,lmax), yy(jmax,kmax,lmax),yz(jmax,kmax,lmax),
     $     zx(jmax,kmax,lmax), zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)


c..   local variables
      
      real cb1,sigma,cb2,akt,cw1,cw2,cw3,cv1,ct1,ct2,ct3,ct4
      data cb1,sigma,cb2,akt,cw1,cw2,cw3,cv1,ct1,ct2,ct3,ct4/0.1355,
     &     0.6666667,0.622,0.41,3.2390678,0.3,2.,7.1,1.,2.,1.2,0.5/
      
      integer j,k,l
      real rcv2,d2min,stilim,rmax,chilim,fturf,cappa2
      real chi,d,fv1,fv2,fv3,ft2,dchi,dfv1,dfv2,dfv3
      real dft2,d2,stilda,r,g,fw,dstild,dr,dg,dfw
      real pro,des,prod,dest,dpro,ddes,tk1,tk2,r5,g6

      rcv2  =1./5.
      d2min =1.e-12
      stilim=1.e-10
      rmax  =10.0
      chilim=1.e-12
      fturf=0.0
      cappa2=akt*akt
      
      do 10 l=2,lmax-1
         do 10 k=2,kmax-1
            do 10 j=2,jmax-1
               chi  = vnu(j,k,l)/vnul(j,k,l)
               d    = sn(j,k,l)
               chi  = max(chi,chilim)
               fv1  = (chi**3)/(chi**3+cv1**3)
               fv2  = 1./(1.+chi*rcv2)
               fv2  = fv2**3
               fv3  = (1.+chi*fv1)*(1.-fv2)/chi
               ft2  = fturf*ct3*exp(-1.*ct4*chi*chi)
               
               dchi = 1./vnul(j,k,l)
               dfv1 = (3.*cv1**3)*(chi**2)*dchi*(1./(chi**3+cv1**3))**2
               dfv2 = (-3.*rcv2)*fv2*dchi/(1.+chi*rcv2)
               dfv3 = ((chi*dfv1 + dchi*fv1)*(1.-fv2) 
     &               - (1.+chi*fv1)*dfv2- dchi*fv3)/chi 
               dft2 = (-2.*ct4)*chi*dchi*ft2
               d2   = max(d**2,d2min)
           
c	for new definition of s_{tilda}, refer aiaa-95-0312

               stilda  = vort(j,k,l)*fv3 
     &                + vnul(j,k,l)/(d2*cappa2*rey)*chi*fv2
               stilda  = max(stilda,stilim)

               r   = vnu(j,k,l)/(d2*cappa2*rey*stilda)	
               r   = min(r,rmax)

               if(r.gt.1.e-8) then
                  r5	= r**5
                  g	= r*(1. + cw2*(r5 - 1.))
                  g6      = g**6
               else
                  r5	= 0.0
                  g	= r*(1. + cw2*(r5 - 1.))
                  g6	= 0.0
               endif

               fw  = (1.+cw3**6)/(g6+cw3**6)
               fw  = g*(fw**(1./6.))

               dstild = vort(j,k,l)*dfv3+vnul(j,k,l)*(dchi*fv2+chi*dfv2)
     &              /(d2*cappa2*rey)
               dr     = vnul(j,k,l)*(dchi-chi*dstild/stilda)
     &              /stilda/(d2*cappa2*rey)
               dg     = dr*(1.+cw2*(6.*r5-1.))
               dfw    = ((1.+cw3**6)/(g6+cw3**6))**(1./6.)
               dfw    =  dfw*dg*(1.- g6/(g6+cw3**6))
               
               pro   = cb1*stilda*(1.-ft2)
               des   = (cw1*fw-cb1/cappa2*ft2)/d2/rey

               prod   = cb1*stilda*(1.-ft2)*vnu(j,k,l)
               dest   = (cw1*fw-cb1/cappa2*ft2)*vnu(j,k,l)
     &              *vnu(j,k,l)/d2/rey

               dpro = pro*dstild/stilda - cb1*stilda*dft2
  
               ddes = (cw1*dfw-cb1/cappa2*dft2)/d2/rey*vnul(j,k,l)*chi
               ddes = ddes + des
               ddes = ddes*vnu(j,k,l)
               dpro = dpro*vnu(j,k,l)
               
               s(j,k,l) = s(j,k,l) + prod-dest
               tk1=max(des*vnu(j,k,l)-pro,0.)
               tk2=max(ddes-dpro,0.)
               sjmat(j,k,l) = sjmat(j,k,l)+tk1+tk2

10	continue


        return
        end

c*************************************************************
      subroutine turbc(vnu,u,v,w,ug,vg,wg,xx,xy,xz,yx,yy,yz,zx,zy,zz,grbc)

      use params_global
      use bcparam
      implicit none

      real vnu(jmax,kmax,lmax)
      real u(jmax,kmax,lmax),v(jmax,kmax,lmax),w(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax),vg(jmax,kmax,lmax),wg(jmax,kmax,lmax) 
      real xx(jmax,kmax,lmax),xy(jmax,kmax,lmax),xz(jmax,kmax,lmax),
     &     yx(jmax,kmax,lmax),yy(jmax,kmax,lmax),yz(jmax,kmax,lmax),
     &     zx(jmax,kmax,lmax),zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)
      type(bc_t) :: grbc

c..   local variables
      integer js,je,ks,ke,ls,le,idir,iproc
      integer j,k,l,ib

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

c.. wall bc at l = 1 (only interior portion of wall)
         if (grbc%ibtyp(ib).eq.5.or.grbc%ibtyp(ib).eq.6) then
         call turbc_wall(vnu,js,je,ks,ke,ls,le,idir)

c.. extrapolate bc at k = 1 
         elseif (grbc%ibtyp(ib).eq.4.or.grbc%ibtyp(ib).eq.10) then
         call turbc_extp(vnu,js,je,ks,ke,ls,le,idir)

c.. symmetric bc for fixed wing at k = 1
         elseif (grbc%ibtyp(ib).eq.11) then
         call turbc_sym(vnu,js,je,ks,ke,ls,le,idir)

c.. averaging bc at k = 1
         elseif (grbc%ibtyp(ib).eq.14) then
         call turbc_av(vnu,js,je,ks,ke,ls,le,idir)

c.. symmetric bc j (3 planes)
         elseif (grbc%ibtyp(ib).eq.22) then
         call turbc_axisym(vnu,js,je,ks,ke,ls,le,idir)

c.. freesream bc
         elseif (grbc%ibtyp(ib).eq.47.or.grbc%ibtyp(ib).eq.48) then
         call turbc_out(vnu,u,v,w,ug,vg,wg,xx,xy,xz,
     &                  yx,yy,yz,zx,zy,zz,js,je,ks,ke,ls,le,idir)

c.. averaging bc for wake
         elseif (grbc%ibtyp(ib).eq.51) then
         call turbc_wake(vnu,js,je,ks,ke,ls,le,idir)
     
c.. bc for parallel runs
         elseif (grbc%ibtyp(ib).eq.87) then
         call turbc_parallel(vnu,js,je,ks,ke,ls,le,idir,iproc)

         elseif (grbc%ibtyp(ib).eq.90) then
         call turbc_coaxial(vnu,js,je,ks,ke,ls,le,idir,iproc)

         endif
      enddo

      return
      end

c*************************************************************
      subroutine turbc_out(vnu,u,v,w,ug,vg,wg,xx,xy,xz,
     &                      yx,yy,yz,zx,zy,zz,js,je,ks,ke,ls,le,idir)


      use params_global
      implicit none


      real vnu(jmax,kmax,lmax)
      real u(jmax,kmax,lmax),v(jmax,kmax,lmax),w(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax),vg(jmax,kmax,lmax),wg(jmax,kmax,lmax) 
      real xx(jmax,kmax,lmax),xy(jmax,kmax,lmax),xz(jmax,kmax,lmax),
     &     yx(jmax,kmax,lmax),yy(jmax,kmax,lmax),yz(jmax,kmax,lmax),
     &     zx(jmax,kmax,lmax),zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)
      integer js,je,ks,ke,ls,le,idir

c..   local variables
      integer j,k,l,j1,k1,l1,iadd,iadir
      real vnuav,uu

      iadd = sign(1,idir)
      iadir = abs(idir)
 
      if(iadir.eq.1) then

         j  = js
         j1 = j + iadd 
         do k=ks,ke
         do l=ls,le
          vnu(j,k,l)=0.1
          uu=   (u(j,k,l)-ug(j,k,l)+u(j1,k,l)-ug(j1,k,l))*xx(j,k,l)
          uu=uu+(v(j,k,l)-vg(j,k,l)+v(j1,k,l)-vg(j1,k,l))*xy(j,k,l)
          uu=uu+(w(j,k,l)-wg(j,k,l)+w(j1,k,l)-wg(j1,k,l))*xz(j,k,l)
          uu=uu*iadd
          if(uu.lt.0.) vnu(j,k,l)=vnu(j1,k,l)
         enddo
         enddo

      elseif(iadir.eq.2) then

         k  = ks
         k1 = k + iadd 
         do l=ls,le
         do j=js,je
          vnu(j,k,l)=0.1
          uu=   (u(j,k,l)-ug(j,k,l)+u(j,k1,l)-ug(j,k1,l))*yx(j,k,l)
          uu=uu+(v(j,k,l)-vg(j,k,l)+v(j,k1,l)-vg(j,k1,l))*yy(j,k,l)
          uu=uu+(w(j,k,l)-wg(j,k,l)+w(j,k1,l)-wg(j,k1,l))*yz(j,k,l)
          uu=uu*iadd
          if(uu.lt.0.) vnu(j,k,l)=vnu(j,k1,l)
         enddo
         enddo

      elseif(iadir.eq.3) then

         l  = ls
         l1 = l + iadd
         do k=ks,ke
         do j=js,je
          vnu(j,k,l)=0.1
          uu=   (u(j,k,l)-ug(j,k,l)+u(j,k,l1)-ug(j,k,l1))*zx(j,k,l)
          uu=uu+(v(j,k,l)-vg(j,k,l)+v(j,k,l1)-vg(j,k,l1))*zy(j,k,l)
          uu=uu+(w(j,k,l)-wg(j,k,l)+w(j,k,l1)-wg(j,k,l1))*zz(j,k,l)
          uu=uu*iadd
          if(uu.lt.0.) vnu(j,k,l)=vnu(j,k,l1)
         enddo
         enddo
 
       endif

      return
      end

c*************************************************************
      subroutine turbc_sym(vnu,js,je,ks,ke,ls,le,idir)

      use params_global
      implicit none

      real vnu(jmax,kmax,lmax)
      integer js,je,ks,ke,ls,le,idir


c..   local variables

      integer j,k,l,j1,k1,l1,iadd,iadir
      real vnuav

      iadd = sign(1,idir)
      iadir = abs(idir)
 
      if(iadir.eq.1) then
         j = js
         j1= j + iadd
         do k=ks,ke
         do l=ls,le
          vnu(j,k,l)=vnu(j1,k,l)
         enddo
         enddo
      elseif(iadir.eq.2) then
         k = ks
         k1= k + iadd
         do j=js,je
         do l=ls,le
          vnu(j,k,l)=vnu(j,k1,l)
         enddo
         enddo
      elseif(iadir.eq.3) then
         l = ls
         l1= l + iadd
         do j=js,je
         do k=ks,ke
          vnu(j,k,l)=vnu(j,k,l1)
         enddo
         enddo
      endif

      return
      end

c*************************************************************
      subroutine turbc_extp(vnu,js,je,ks,ke,ls,le,idir)

      use params_global
      implicit none

      real vnu(jmax,kmax,lmax)
      integer js,je,ks,ke,ls,le,idir


c..   local variables

      integer j,k,l,j1,k1,l1,iadd,iadir
      real vnuav

      iadd = sign(1,idir)
      iadir = abs(idir)
 
      if(iadir.eq.1) then
         j = js
         j1= j + iadd
         do k=ks,ke
         do l=ls,le
          vnu(j,k,l)=vnu(j1,k,l)
         enddo
         enddo
      elseif(iadir.eq.2) then
         k = ks
         k1= k + iadd
         do l=ls,le
         do j=js,je
          vnu(j,k,l)=vnu(j,k1,l)
         enddo
         enddo
      elseif(iadir.eq.3) then
         l = ls
         l1= l + iadd
         do j=js,je
         do k=ks,ke
          vnu(j,k,l)=vnu(j,k,l1)
         enddo
         enddo
      endif

      return
      end

c*************************************************************
      subroutine turbc_av(vnu,js,je,ks,ke,ls,le,idir)

      use params_global
      implicit none

      real vnu(jmax,kmax,lmax)
      integer js,je,ks,ke,ls,le,idir


c..   local variables

      integer j,k,l,j1,k1,l1,iadd,iadir
      real sumv

      iadd = sign(1,idir)
      iadir = abs(idir)
 
      if(iadir.eq.1) then
         j = js
         j1= j + iadd
         do k=ks,ke
         do l=ls,le
          vnu(j,k,l)=vnu(j1,k,l)
         enddo
         enddo

         do k = ks,ke
            sumv = 0.0
            do l = ls,le
               sumv = sumv + vnu(j,k,l)
            enddo
            sumv = sumv/(le-ls+1)
            do l = ls,le
               vnu(j,k,l) = sumv
            enddo
         enddo

      elseif(iadir.eq.2) then
         k = ks
         k1= k + iadd
         do j=js,je
         do l=ls,le
          vnu(j,k,l)=vnu(j,k1,l)
         enddo
         enddo

         do l = ls,le
            sumv = 0.0
            do j = js,je
               sumv = sumv + vnu(j,k,l)
            enddo
            sumv = sumv/(je-js+1)
            do j = js,je
               vnu(j,k,l) = sumv
            enddo
         enddo

      elseif(iadir.eq.3) then
         l = ls
         l1= l + iadd
         do j=js,je
         do k=ks,ke
          vnu(j,k,l)=vnu(j,k,l1)
         enddo
         enddo

         do j = js,je
            sumv = 0.0
            do k = ks,ke
               sumv = sumv + vnu(j,k,l)
            enddo
            sumv = sumv/(ke-ks+1)
            do k = ks,ke
               vnu(j,k,l) = sumv
            enddo
         enddo

      endif

      return
      end

c*************************************************************
      subroutine turbc_vortex(vnu,vnul)

      use params_global
      implicit none

      real vnul(jmax,kmax,lmax),vnu(jmax,kmax,lmax)

c..   local variables
       
      integer j,k,l,j1


c      print*,'vortex turbc'

c     extrapolation for j=jmax. actually should
c     use characteristic info, but using bare
c     minimum now.
       
      j=jmax
      j1=jmax-1
      do l=1,lmax
         do k=1,kmax
            vnu(j,k,l)=vnu(j1,k,l)
         enddo
      enddo
       
       
      return
      end

c*************************************************************
      subroutine turbc_wake(vnu,js,je,ks,ke,ls,le,idir)

      use params_global
      implicit none

      real vnu(jmax,kmax,lmax)
      integer js,je,ks,ke,ls,le,idir


c..   local variables

      integer j,k,l,j1,k1,l1,jj,iadd,iadir
      real vnuav

      iadd = sign(1,idir)
      iadir = abs(idir)
 
      if(iadir.eq.1) then
         print*,'idir = ',idir,' is not implemented in turbc_wake'
      elseif(iadir.eq.2) then
         k  = ks
         k1 = k + iadd
         do l=ls,le
            do j=js,je
               jj = jmax - j + 1
               vnuav = 0.5*(vnu(j,k1,l)+vnu(jj,k1,l))
               vnu(j,k,l)  = vnuav
               vnu(jj,k,l) = vnuav
            enddo
         enddo
      elseif(iadir.eq.3) then
         l  = ls
         l1 = l + iadd
         do k = ks,ke
            do j = js, je
               jj = jmax - j + 1
               vnuav = 0.5*(vnu(j,k,l1)+vnu(jj,k,l1))
               vnu(j,k,l)  = vnuav
               vnu(jj,k,l) = vnuav
            enddo
         enddo
      endif
      
      return
      end

c*************************************************************
      subroutine turbc_wall(vnu,js,je,ks,ke,ls,le,idir)

      use params_global
      implicit none

      real vnu(jmax,kmax,lmax)
      integer js,je,ks,ke,ls,le,idir


c..   local variables

      integer j,k,l,iadd,iadir
      
      iadd = sign(1,idir)
      iadir = abs(idir)
 
      if(iadir.eq.1) then
         j = js
         do k=ks,ke
            do l=ls,le
               vnu(j,k,l)=0.0
            enddo
         enddo
      elseif(iadir.eq.2) then
         k = ks
         do l=ls,le
            do j=js,je
               vnu(j,k,l)=0.0
            enddo
         enddo
      elseif(iadir.eq.3) then
         l = ls
         do k=ks,ke
            do j=js,je
               vnu(j,k,l)=0.0
            enddo
         enddo
      endif
      
      return
      end

c*************************************************************
      subroutine turbc_axisym(vnu,js,je,ks,ke,ls,le,idir)

      use params_global
      implicit none

      real vnu(jmax,kmax,lmax)
      integer js,je,ks,ke,ls,le,idir


c..   local variables

      integer j,k,l,jc,jj,jj1,iadd,iadir
c      integer j1,j2,j3,jc1,jc2,jc3,k1,k2
      
      iadd = sign(1,idir)
      iadir = abs(idir)
 
      jj  = je - js + 1
      if(idir.eq.1) then
      do j = js,je
         jj1 = j - js
         jc = jmax - 2*jj + jj1
 
         do k = ks,ke
         do l = ls,le
            vnu(j,k,l) = vnu(jc,k,l)
         enddo
         enddo
      enddo

      elseif(idir.eq.-1) then
      do j = js,je
         jj1 = je - j
         jc = 1 + 2*jj - jj1 
 
         do k = ks,ke
         do l = ls,le
            vnu(j,k,l) = vnu(jc,k,l)
         enddo
         enddo
      enddo

      else 
        print*,'idir = ',idir,' is not implemented in turbc_axisym'
      endif


c      if(idir.eq.1) then
c      j1  = 1
c      j2  = 2
c      j3  = 3
c      jc1 = jmax-6 
c      jc2 = jmax-5 
c      jc3 = jmax-4 
c 
c      k1 = 2
c      k2 = km
c
c      do k = k1,k2
c      do l = 1,lmax
c        vnu(j1,k,l) = vnu(jc1,k,l)
c        vnu(j2,k,l) = vnu(jc2,k,l)
c        vnu(j3,k,l) = vnu(jc3,k,l)
c      enddo
c      enddo
c       
c      elseif(idir.eq.-1) then
c
c      j1  = jmax
c      j2  = jmax-1
c      j3  = jmax-2
c      jc1 = 7 
c      jc2 = 6 
c      jc3 = 5 
c 
c      k1 = 2
c      k2 = km
c
c      do k = k1,k2
c      do l = 1,lmax
c        vnu(j1,k,l) = vnu(jc1,k,l)
c        vnu(j2,k,l) = vnu(jc2,k,l)
c        vnu(j3,k,l) = vnu(jc3,k,l)
c      enddo
c      enddo
c
c      endif
      
c      j  = 1
c      jp = j + 1
c      jpp= jp + 1
c      jj = jmax
c      jm1 = jj - 1
c      jmm= jm1 - 1
c      sixth=1./6.
c      
c      do l = 1, lmax
c         do k = 1, kmax
c            vnu(j,k,l) = (-vnu(jpp,k,l)+4.*vnu(jp,k,l)
c     >           +4.*vnu(jm1,k,l) -vnu(jmm,k,l))*sixth
c            vnu(jj,k,l)= vnu(j,k,l)
c         enddo
c      enddo

      return
      end

c*************************************************************
      subroutine turbc_parallel(vnut,js,je,ks,ke,ls,le,idir,iproc)

c,,   THIS IS IMPLEMENTED IN THE BC SUBROUTINE TO SAVE COMPUTATIONAL TIME 

c*************************************************************
      return
      end

c*************************************************************
      subroutine turbc_coaxial(vnut,js,je,ks,ke,ls,le,idir,iproc)

c,,   THIS IS IMPLEMENTED IN THE BC SUBROUTINE TO SAVE COMPUTATIONAL TIME 

c*************************************************************
      return
      end

c*************************************************************
      subroutine lsolvej(a,b,c,s,jm,km,l,jmax,kmax,lmax,mdim)
c*************************************************************
      
      implicit none
      
      integer jm,km,l,jmax,kmax,lmax,mdim
      real a(mdim,mdim),b(mdim,mdim),c(mdim,mdim)
      real s(jmax,kmax,lmax)
      
c..   local variables

      integer k,j
      real bb

c..   first executable statement

      do k = 2,km
         bb        = 1./b(2,k)
         s(2,k,l)  = s(2,k,l)*bb
         c(2,k)    = c(2,k)*bb
      enddo
 
      do  j = 3,jm
         do  k = 2,km
            bb      = 1./(b(j,k) - a(j,k)*c(j-1,k))
            s(j,k,l)= (s(j,k,l)  - a(j,k)*s(j-1,k,l))*bb
            c(j,k)  = c(j,k)*bb
         enddo
      enddo
      
      do j = jm-1,2,-1
         do k = 2,km
            s(j,k,l) = s(j,k,l) - c(j,k)*s(j+1,k,l)
         enddo
      enddo
      
      
      return
      end
      
c*************************************************************
      subroutine lsolvel(a,b,c,s,jm,k,lm,jmax,kmax,lmax,mdim)
c*************************************************************

      implicit none

      integer jm,k,lm,jmax,kmax,lmax,mdim
      real a(mdim,mdim),b(mdim,mdim),c(mdim,mdim)
      real s(jmax,kmax,lmax)
      
c..   local variables
      
      integer j,l
      real bb
	

      
      do j = 2,jm
         bb        = 1./b(j,2)
         s(j,k,2)  = s(j,k,2)*bb
         c(j,2)    = c(j,2)*bb
      enddo
      
      do  l = 3,lm
         do  j = 2,jm
            bb      = 1./(b(j,l) - a(j,l)*c(j,l-1))
            s(j,k,l)= (s(j,k,l) - a(j,l)*s(j,k,l-1))*bb
            c(j,l)  = c(j,l)*bb
         enddo
      enddo
      
      do l = lm-1,2,-1
         do j = 2,jm
            s(j,k,l) = s(j,k,l) - c(j,l)*s(j,k,l+1)
         enddo
      enddo
      
	
      return
      end

c*************************************************************
      subroutine lsolvek(a,b,c,s,jm,km,l,jmax,kmax,lmax,mdim)
c*************************************************************

      implicit none

      integer jm,km,l,jmax,kmax,lmax,mdim
      
      real a(mdim,mdim),b(mdim,mdim),c(mdim,mdim)
      real s(jmax,kmax,lmax)
      
c..   local variables

      integer j,k
      real bb
     
      do j = 2,jm
         bb        = 1./b(j,2)
         s(j,2,l)  = s(j,2,l)*bb
         c(j,2)    = c(j,2)*bb
      enddo
      
      do  k = 3,km
         do  j = 2,jm
            bb      = 1./(b(j,k) - a(j,k)*c(j,k-1))
            s(j,k,l)= (s(j,k,l) - a(j,k)*s(j,k-1,l))*bb
            c(j,k)  = c(j,k)*bb
         enddo
      enddo
      
      do k = km-1,2,-1
         do j = 2,jm
            s(j,k,l) = s(j,k,l) - c(j,k)*s(j,k+1,l)
         enddo
      enddo
      
	
      return
      end


