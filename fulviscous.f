C***********************************************************************
      subroutine fulvisrhs(vnu,vnu0,vmut,vmul,q,s,xx,xy,
     $		xz,yx,yy,yz,zx,zy,
     $          zz,x,y,z,ug,vg,wg,tscale,kkr,kkp,iblank,grbc)
C
C  Compute the Full viscous RHS 
C  
C  Code by dkarthik
C  cleaned up a little bit by jaina
C
C***********************************************************************
      use params_global
      use bcparam
      implicit none

      real q(jmax,kmax,lmax,nd),s(jmax,kmax,lmax,nv)
      real xx(jmax,kmax,lmax), xy(jmax,kmax,lmax),xz(jmax,kmax,lmax),
     $     yx(jmax,kmax,lmax), yy(jmax,kmax,lmax),yz(jmax,kmax,lmax),
     $     zx(jmax,kmax,lmax), zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)

      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax),vg(jmax,kmax,lmax),wg(jmax,kmax,lmax)
      real vnu(jmax,kmax,lmax),vmut(jmax,kmax,lmax),vmul(jmax,kmax,lmax)
      
      real vnu0(jmax,kmax,lmax),tscale(jmax,kmax,lmax)

      integer kkr(kmax),kkp(kmax),iblank(jmax,kmax,lmax)
      type(bc_t) :: grbc

c..   local variables
      real, allocatable :: f1(:,:,:),f2(:,:,:)

      integer icross,iviscj,ivisck,iviscl
      integer js,je,ks,ke,ls,le
      integer j,k,l,n
      integer jm1,km1,lm1
      real st1,st,st2,st3
      real gkpr,prt,dre,c2b,c2bp
      real ra,uvel,vvel,wvel,ei,tt

      allocate(f1(mdim,mdim,5),f2(mdim,mdim,5))

c**** first exectutable statement
      
      icross=1
      iviscj=1
      ivisck=1
      iviscl=1
      
      js=2
      je=jmax-1
      ks=2
      ke=kmax-1
      ls=2
      le=lmax-1

      if(lamin) then
         do l = 1,lmax
            do k = 1,kmax
               do j = 1,jmax
                  vmut(j,k,l)=0.0
               end do
            end do
         end do
      else
         if (iturb.eq.1) then
            call vmutur( x,y,z,q,s,vmut,xx,xy,xz,yx,yy,yz,
     <           zx,zy,zz,ug,vg,wg,kkr,kkp)
         else
c	    call cpu_time(st)
            call vmu_sa(x,y,z,q,vmut,xx,xy,xz,yx,yy,yz,zx,zy,zz,
     >           ug,vg,wg,vnu,vnu0,vmul,tscale,iblank,grbc)

	    !insurance
c	    call cpu_time(st1)
c	    print *,'cpu time for SA model=', st1-st
         endif
      endif
      jm=jmax-1
      km=kmax-1
      lm=lmax-1

c..   calculate laminar viscosity using Sutherland's formula

      gkpr = gamma/pr
      prt = 0.9
      dre  = .5/rey
      c2b  =198.6/tinf
      c2bp = c2b +1.
      
      do l = 1,lmax
         do k = 1,kmax
            do j = 1,jmax
               q(j,k,l,1)=q(j,k,l,1)*q(j,k,l,6)	
               q(j,k,l,2)=q(j,k,l,2)*q(j,k,l,6)	
               q(j,k,l,3)=q(j,k,l,3)*q(j,k,l,6)	
               q(j,k,l,4)=q(j,k,l,4)*q(j,k,l,6)	
               q(j,k,l,5)=q(j,k,l,5)*q(j,k,l,6)	
               ra   = 1./q(j,k,l,1)
               uvel = q(j,k,l,2)*ra
               vvel = q(j,k,l,3)*ra
               wvel = q(j,k,l,4)*ra
               ei   = q(j,k,l,5)*ra-.5*(uvel**2+vvel**2+wvel**2)
               tt    = ggm1*ei
               vmul(j,k,l)  = c2bp*tt*sqrt(tt)/(c2b + tt)
            end do
         end do
      end do


c..   viscous terms in j-direction (wrap around)

      if(iviscj.eq.1) then
         do l=ls,le

            call vflj(l,rey,pr,prt,q,gamma,vmul,vmut,xx,xy,xz,f1,
     >           jmax,kmax,lmax,nd,mdim)

!            call vflj(l,rey,pr,prt,q,gamma,vmul,vmut,xx,xy,xz,f1,
!     >  tmp2(:,:,1),tmp2(:,:,2),tmp2(:,:,3),tmp2(:,:,4),tmp2(:,:,5),
!     >	tmp2(:,:,6),tmp2(:,:,7),tmp2(:,:,8),tmp2(:,:,9),tmp2(:,:,10),
!     >  tmp2(:,:,11),jmax,kmax,lmax,nd,mdim)

            do k = ks,ke
               do j = js,je
                  jm1         = j-1
                  s(j,k,l,2) = s(j,k,l,2) + (f1(j,k,2) - f1(jm1,k,2))
                  s(j,k,l,3) = s(j,k,l,3) + (f1(j,k,3) - f1(jm1,k,3))
                  s(j,k,l,4) = s(j,k,l,4) + (f1(j,k,4) - f1(jm1,k,4))
                  s(j,k,l,5) = s(j,k,l,5) + (f1(j,k,5) - f1(jm1,k,5))
               end do	
            end do	
         end do	
      endif

c..   viscouse terms in k-direction (spanwise)

      if(ivisck.eq.1) then

         do l=ls,le
            call vflk(l,rey,pr,prt,q,gamma,vmul,vmut,yx,yy,yz,f1,
     >           jmax,kmax,lmax,nd,mdim)

            do k = ks,ke
               do j = js,je
                  km1         = k-1
                  s(j,k,l,2) = s(j,k,l,2) + (f1(j,k,2) - f1(j,km1,2))
                  s(j,k,l,3) = s(j,k,l,3) + (f1(j,k,3) - f1(j,km1,3))
                  s(j,k,l,4) = s(j,k,l,4) + (f1(j,k,4) - f1(j,km1,4))
                  s(j,k,l,5) = s(j,k,l,5) + (f1(j,k,5) - f1(j,km1,5))
               end do	
            end do	
         end do	

      endif

c..   viscous terms in l-direction (normal)

      if(iviscl.eq.1) then
         do k=ks,ke
            call vfll(k,rey,pr,prt,q,gamma,vmul,vmut,zx,zy,zz,f1,
     >           jmax,kmax,lmax,nd,mdim)
            
            do l = ls,le
               do j = js,je
                  lm1         = l-1
                  s(j,k,l,2) = s(j,k,l,2) + (f1(j,l,2) - f1(j,lm1,2))
                  s(j,k,l,3) = s(j,k,l,3) + (f1(j,l,3) - f1(j,lm1,3))
                  s(j,k,l,4) = s(j,k,l,4) + (f1(j,l,4) - f1(j,lm1,4))
                  s(j,k,l,5) = s(j,k,l,5) + (f1(j,l,5) - f1(j,lm1,5))
               end do	
            end do	
         end do	
      endif

c..   cross derivative terms (refer hoffman and chang)
      if(icross.eq.1) then

c..   j-k terms

      do l=ls,le
         call vfljk(l,rey,pr,prt,q,gamma,vmul,vmut,xx,xy,xz,yx,yy,yz,
     >        f1,f2,jmax,kmax,lmax,nd,mdim)


         do k = ks,ke
            km1         = k-1
            do j = js,je
               jm1         = j-1
               s(j,k,l,2) = s(j,k,l,2) + (f1(j,k,2) - f1(jm1,k,2))
     &                                 + (f2(j,k,2) - f2(j,km1,2))
               s(j,k,l,3) = s(j,k,l,3) + (f1(j,k,3) - f1(jm1,k,3))
     &                                 + (f2(j,k,3) - f2(j,km1,3))
               s(j,k,l,4) = s(j,k,l,4) + (f1(j,k,4) - f1(jm1,k,4))
     &                                 + (f2(j,k,4) - f2(j,km1,4))
               s(j,k,l,5) = s(j,k,l,5) + (f1(j,k,5) - f1(jm1,k,5))
     &                                 + (f2(j,k,5) - f2(j,km1,5))
            end do
         end do
         
      end do	
	

c..   j-l terms

      do k=ks,ke
         call vfljl(k,rey,pr,prt,q,gamma,vmul,vmut,xx,xy,xz,zx,zy,zz,
     >        f1,f2,jmax,kmax,lmax,nd,mdim)

         
         do l = ls,le
            lm1         = l-1
            do j = js,je
               jm1         = j-1
               s(j,k,l,2) = s(j,k,l,2) + (f1(j,l,2) - f1(jm1,l,2))
     &                                 + (f2(j,l,2) - f2(j,lm1,2))
               s(j,k,l,3) = s(j,k,l,3) + (f1(j,l,3) - f1(jm1,l,3))
     &                                 + (f2(j,l,3) - f2(j,lm1,3))
               s(j,k,l,4) = s(j,k,l,4) + (f1(j,l,4) - f1(jm1,l,4))
     &                                 + (f2(j,l,4) - f2(j,lm1,4))
               s(j,k,l,5) = s(j,k,l,5) + (f1(j,l,5) - f1(jm1,l,5))
     &                                 + (f2(j,l,5) - f2(j,lm1,5))
            enddo
         enddo
      
      end do	


c..   k-l terms

      do l=ls,le
	call vflkl(l,rey,pr,prt,q,gamma,vmul,vmut,yx,yy,yz,zx,zy,zz,
     >        f1,f2,jmax,kmax,lmax,nd,mdim)
      
	
        do k = ks,ke
           km1         = k-1
           do j = js,je
              s(j,k,l,2) = s(j,k,l,2) + (f1(j,k,2) - f1(j,km1,2))
     &                                +  f2(j,k,2)
              s(j,k,l,3) = s(j,k,l,3) + (f1(j,k,3) - f1(j,km1,3))
     &                                +  f2(j,k,3)
              s(j,k,l,4) = s(j,k,l,4) + (f1(j,k,4) - f1(j,km1,4))
     &                                +  f2(j,k,4)
              s(j,k,l,5) = s(j,k,l,5) + (f1(j,k,5) - f1(j,km1,5))
     &                                +  f2(j,k,5)
           end do
        end do
      end do

      endif

c..   rescale by jacobian

      do l = 1,lmax
         do k = 1,kmax
            do j = 1,jmax
               
               q(j,k,l,1)=q(j,k,l,1)/q(j,k,l,6)	
               q(j,k,l,2)=q(j,k,l,2)/q(j,k,l,6)	
               q(j,k,l,3)=q(j,k,l,3)/q(j,k,l,6)	
               q(j,k,l,4)=q(j,k,l,4)/q(j,k,l,6)	
               q(j,k,l,5)=q(j,k,l,5)/q(j,k,l,6)	
               
            end do
         end do
      end do

      return
      end
