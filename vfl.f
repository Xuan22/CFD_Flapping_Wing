      subroutine vflj  ( l,
     &                   rey,pr,prt,
     &                   q,vgamma,vmul,vmut,
     &                   xx,xy,xz,
     &                   f,jmax,kmax,lmax,nd,mdim )
c
c   compute a jk-plane of j-direction viscous fluxes.
c   follows the arc2d approach for getting values at half grid points.
c
c*******************************************************************

      implicit none

      integer l,jmax,kmax,lmax,nd,mdim
      real rey,pr,prt,vgamma
      
      real q(jmax,kmax,lmax,6)
      real vmul(jmax,kmax,lmax),vmut(jmax,kmax,lmax)
      real xx(jmax,kmax,lmax),xy(jmax,kmax,lmax),xz(jmax,kmax,lmax)
      real f(mdim,mdim,5)

      real, allocatable :: s0(:,:),s1(:,:),s2(:,:)
      real, allocatable :: s3(:,:),s4(:,:),s5(:,:),s6(:,:)
      real, allocatable :: u(:,:),v(:,:),w(:,:),ei(:,:)
!      real s0(mdim,mdim),s1(mdim,mdim),s2(mdim,mdim)
!      real s3(mdim,mdim),s4(mdim,mdim),s5(mdim,mdim),s6(mdim,mdim)
!      real u(mdim,mdim),v(mdim,mdim),w(mdim,mdim),ei(mdim,mdim)

c..   local variables

      integer j,k,jm,jp
      real r3,prtr,dre,dprre,drho,xs2
      real e0,v2,xx2,xy2,xz2,rv,xs2rv,r3rv
      real gamav2,vmula2,vmuta2,vmud2,gamkap,uav2,vav2,wav2
      real s0av2,s1av2,s2av2,s3av2,s4av2,s5av2,s6av2
      real du,dv,dw,dei

      allocate(s0(mdim,mdim),s1(mdim,mdim),s2(mdim,mdim))
      allocate(s3(mdim,mdim),s4(mdim,mdim),s5(mdim,mdim),s6(mdim,mdim))
      allocate(u(mdim,mdim),v(mdim,mdim),w(mdim,mdim),ei(mdim,mdim))

      r3=1.d0/3.d0
      prtr  = pr/prt
      dre   = 0.25/rey
      dprre = 0.25/(pr*rey)
c
      do 10 k = 1,kmax
      do 10 j = 1,jmax
         drho       = 1./q(j,k,l,1)
         u(j,k)     = q(j,k,l,2)*drho
         v(j,k)     = q(j,k,l,3)*drho
         w(j,k)     = q(j,k,l,4)*drho
         e0         = q(j,k,l,5)*drho
         v2         = 0.5*(u(j,k)**2 + v(j,k)**2 + w(j,k)**2)
         ei(j,k)    = e0 - v2
         xx2        = xx(j,k,l)**2
         xy2        = xy(j,k,l)**2
         xz2        = xz(j,k,l)**2
         xs2        = xx2 + xy2 + xz2
         rv         = 1./q(j,k,l,6)
         xs2rv      = xs2*rv
         r3rv       = r3*rv
         s0(j,k)    = xs2rv
         s1(j,k)    = xs2rv + r3rv*xx2
         s2(j,k)    = xs2rv + r3rv*xy2
         s3(j,k)    = xs2rv + r3rv*xz2
         s4(j,k)    = (r3rv*xx(j,k,l))*xy(j,k,l)
         s5(j,k)    = (r3rv*xx(j,k,l))*xz(j,k,l)
         s6(j,k)    = (r3rv*xy(j,k,l))*xz(j,k,l)
   10    continue
c
      do 20 k = 1,kmax
      do 20 j = 1,jmax-1
         jp         = j+1
         gamav2     = vgamma+vgamma
         vmula2     = vmul(jp,k,l) + vmul(j,k,l)
         vmuta2     = vmut(jp,k,l) + vmut(j,k,l)
         vmud2      = (vmula2 + vmuta2)*dre
         gamkap     = gamav2*(vmula2 + prtr*vmuta2)*dprre
         uav2       = u(jp,k) + u(j,k)
         vav2       = v(jp,k) + v(j,k)
         wav2       = w(jp,k) + w(j,k)
         s0av2      = s0(jp,k) + s0(j,k)
         s1av2      = s1(jp,k) + s1(j,k)
         s2av2      = s2(jp,k) + s2(j,k)
         s3av2      = s3(jp,k) + s3(j,k)
         s4av2      = s4(jp,k) + s4(j,k)
         s5av2      = s5(jp,k) + s5(j,k)
         s6av2      = s6(jp,k) + s6(j,k)
         du         = u(jp,k)  - u(j,k)
         dv         = v(jp,k)  - v(j,k)
         dw         = w(jp,k)  - w(j,k)
         dei        = ei(jp,k) - ei(j,k)
c        f(j,k,1)   = 0.
         f(j,k,2)   = (s1av2*du + s4av2*dv + s5av2*dw)*vmud2
         f(j,k,3)   = (s4av2*du + s2av2*dv + s6av2*dw)*vmud2
         f(j,k,4)   = (s5av2*du + s6av2*dv + s3av2*dw)*vmud2
         f(j,k,5)   = 0.5*(f(j,k,2)*uav2 + f(j,k,3)*vav2 + f(j,k,4)*wav2
     &                                   + s0av2*dei*gamkap)
   20    continue
c
c
      deallocate(s0,s1,s2,s3,s4,s5,s6,u,v,w,ei)

      return
      end


      subroutine vflk  ( l,
     &                   rey,pr,prt,
     &                   q,vgamma,vmul,vmut,
     &                   yx,yy,yz,
     &                   f,jmax,kmax,lmax,nd,mdim)
c
c   compute a jk-plane of k-direction viscous fluxes.
c
c*******************************************************************

      integer l,jmax,kmax,lmax,nd,mdim
      real rey,pr,prt,vgamma
      
      real q(jmax,kmax,lmax,nd)
      real vmul(jmax,kmax,lmax),vmut(jmax,kmax,lmax)
      real yx(jmax,kmax,lmax),yy(jmax,kmax,lmax),yz(jmax,kmax,lmax)
      
      real f(mdim,mdim,5)
!      real s0(mdim,mdim),s1(mdim,mdim),s2(mdim,mdim),s3(mdim,mdim),
!     &     s4(mdim,mdim),s5(mdim,mdim),s6(mdim,mdim)
!      real u(mdim,mdim),v(mdim,mdim),w(mdim,mdim),ei(mdim,mdim)

      real, allocatable :: s0(:,:),s1(:,:),s2(:,:),s3(:,:)
      real, allocatable :: s4(:,:),s5(:,:),s6(:,:)
      real, allocatable :: u(:,:),v(:,:),w(:,:),ei(:,:)

c..   local variables

      integer k,j
      real r3,prtr,dre,dprre
      real drho,e0,v2,yx2,yy2,yz2,ys2,rv,ys2rv,r3rv
      real gamav2,vmula2,vmuta2,vmud2,gamkap,uav2
      real vav2,wav2,s0av2,s1av2,s2av2,s3av3,s4av2
      real s5av2,s6av2,du,dv,dw,dei

      allocate(s0(mdim,mdim),s1(mdim,mdim),s2(mdim,mdim))
      allocate(s3(mdim,mdim),s4(mdim,mdim),s5(mdim,mdim),s6(mdim,mdim))
      allocate(u(mdim,mdim),v(mdim,mdim),w(mdim,mdim),ei(mdim,mdim))


      r3=1.d0/3.d0
      prtr  = pr/prt
      dre   = 0.25/rey
      dprre = 0.25/(pr*rey)
c
      do 10 k = 1,kmax
      do 10 j = 1,jmax
         drho       = 1./q(j,k,l,1)
         u(j,k)     = q(j,k,l,2)*drho
         v(j,k)     = q(j,k,l,3)*drho
         w(j,k)     = q(j,k,l,4)*drho
         e0         = q(j,k,l,5)*drho
         v2         = 0.5*(u(j,k)**2 + v(j,k)**2 + w(j,k)**2)
         ei(j,k)    = e0 - v2
         yx2        = yx(j,k,l)**2
         yy2        = yy(j,k,l)**2
         yz2        = yz(j,k,l)**2
         ys2        = yx2 + yy2 + yz2
         rv         = 1./q(j,k,l,6)
         ys2rv      = ys2*rv
         r3rv       = r3*rv
         s0(j,k)    = ys2rv
         s1(j,k)    = ys2rv + r3rv*yx2
         s2(j,k)    = ys2rv + r3rv*yy2
         s3(j,k)    = ys2rv + r3rv*yz2
         s4(j,k)    = (r3rv*yx(j,k,l))*yy(j,k,l)
         s5(j,k)    = (r3rv*yx(j,k,l))*yz(j,k,l)
         s6(j,k)    = (r3rv*yy(j,k,l))*yz(j,k,l)
   10    continue
c
      do 20 k = 1,kmax-1
      do 20 j = 1,jmax
         kp         = k+1
         gamav2     = vgamma + vgamma
         vmula2     = vmul(j,kp,l) + vmul(j,k,l)
         vmuta2     = vmut(j,kp,l) + vmut(j,k,l)
         vmud2      = (vmula2 + vmuta2)*dre
         gamkap     = gamav2*(vmula2 + prtr*vmuta2)*dprre
         uav2       = u(j,kp) + u(j,k)
         vav2       = v(j,kp) + v(j,k)
         wav2       = w(j,kp) + w(j,k)
         s0av2      = s0(j,kp) + s0(j,k)
         s1av2      = s1(j,kp) + s1(j,k)
         s2av2      = s2(j,kp) + s2(j,k)
         s3av2      = s3(j,kp) + s3(j,k)
         s4av2      = s4(j,kp) + s4(j,k)
         s5av2      = s5(j,kp) + s5(j,k)
         s6av2      = s6(j,kp) + s6(j,k)
         du         = u(j,kp)  - u(j,k)
         dv         = v(j,kp)  - v(j,k)
         dw         = w(j,kp)  - w(j,k)
         dei        = ei(j,kp) - ei(j,k)
c        f(j,k,1)   = 0.
         f(j,k,2)   = (s1av2*du + s4av2*dv + s5av2*dw)*vmud2
         f(j,k,3)   = (s4av2*du + s2av2*dv + s6av2*dw)*vmud2
         f(j,k,4)   = (s5av2*du + s6av2*dv + s3av2*dw)*vmud2
         f(j,k,5)   = 0.5*(f(j,k,2)*uav2 + f(j,k,3)*vav2 + f(j,k,4)*wav2
     &                                   + s0av2*dei*gamkap)
   20    continue
c
c
         deallocate(s0,s1,s2,s3,s4,s5,s6,u,v,w,ei)
      return
      end


      subroutine vfll  ( k,
     &                   rey,pr,prt,
     &                   q,vgamma,vmul,vmut,
     &                   zx,zy,zz,
     &                   f,jmax,kmax,lmax,nd,mdim )
c
c   compute a jl-plane of l-direction viscous fluxes.
c   follows the arc2d approach for getting values at half grid points.
c
c*******************************************************************

      implicit none
      
      integer k,jmax,kmax,lmax,nd,mdim

      real rey,pr,prt,vgamma
      real q(jmax,kmax,lmax,nd)
      real vmul(jmax,kmax,lmax),vmut(jmax,kmax,lmax)
      real zx(jmax,kmax,lmax),zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)
      
      real f(mdim,mdim,5)
!      real s0(mdim,mdim),s1(mdim,mdim),s2(mdim,mdim),s3(mdim,mdim)
!      real s4(mdim,mdim),s5(mdim,mdim),s6(mdim,mdim)

!      real u(mdim,mdim),v(mdim,mdim),w(mdim,mdim),ei(mdim,mdim)

      real, allocatable :: s0(:,:),s1(:,:),s2(:,:)
      real, allocatable :: s3(:,:),s4(:,:),s5(:,:),s6(:,:)
      real, allocatable :: u(:,:),v(:,:),w(:,:),ei(:,:)


c..   local variables

      integer l,j,lp
      real r3,prtr,dre,dprre,drho
      real e0,v2,zx2,zy2,zz2,zs2,rv,zs2rv,r3rv
      real gamav2,vmula2,vmuta2,vmud2,gamkap
      real uav2,vav2,wav2
      real s0av2,s1av2,s2av2,s3av2,s4av2,s5av2,s6av2
      real du,dv,dw,dei

      allocate(s0(mdim,mdim),s1(mdim,mdim),s2(mdim,mdim))
      allocate(s3(mdim,mdim),s4(mdim,mdim),s5(mdim,mdim),s6(mdim,mdim))
      allocate(u(mdim,mdim),v(mdim,mdim),w(mdim,mdim),ei(mdim,mdim))

      r3=1.d0/3.d0
      prtr  = pr/prt
      dre   = 0.25/rey
      dprre = 0.25/(pr*rey)
c
      do 10 l = 1,lmax
      do 10 j = 1,jmax
         drho       = 1./q(j,k,l,1)
         u(j,l)     = q(j,k,l,2)*drho
         v(j,l)     = q(j,k,l,3)*drho
         w(j,l)     = q(j,k,l,4)*drho
         e0         = q(j,k,l,5)*drho
         v2         = 0.5*(u(j,l)**2 + v(j,l)**2 + w(j,l)**2)
         ei(j,l)    = e0 - v2
         zx2        = zx(j,k,l)**2
         zy2        = zy(j,k,l)**2
         zz2        = zz(j,k,l)**2
         zs2        = zx2 + zy2 + zz2
         rv         = 1./q(j,k,l,6)
         zs2rv      = zs2*rv
         r3rv       = r3*rv
         s0(j,l)    = zs2rv
         s1(j,l)    = zs2rv + r3rv*zx2
         s2(j,l)    = zs2rv + r3rv*zy2
         s3(j,l)    = zs2rv + r3rv*zz2
         s4(j,l)    = (r3rv*zx(j,k,l))*zy(j,k,l)
         s5(j,l)    = (r3rv*zx(j,k,l))*zz(j,k,l)
         s6(j,l)    = (r3rv*zy(j,k,l))*zz(j,k,l)
   10    continue
c
      do 20 l = 1,lmax-1
      do 20 j = 1,jmax
         lp         = l+1
         gamav2     = vgamma+vgamma
         vmula2     = vmul(j,k,lp) + vmul(j,k,l)
         vmuta2     = vmut(j,k,lp) + vmut(j,k,l)
         vmud2      = (vmula2 + vmuta2)*dre
         gamkap     = gamav2*(vmula2 + prtr*vmuta2)*dprre
         uav2       = u(j,lp) + u(j,l)
         vav2       = v(j,lp) + v(j,l)
         wav2       = w(j,lp) + w(j,l)
         s0av2      = s0(j,lp) + s0(j,l)
         s1av2      = s1(j,lp) + s1(j,l)
         s2av2      = s2(j,lp) + s2(j,l)
         s3av2      = s3(j,lp) + s3(j,l)
         s4av2      = s4(j,lp) + s4(j,l)
         s5av2      = s5(j,lp) + s5(j,l)
         s6av2      = s6(j,lp) + s6(j,l)
         du         = u(j,lp)  - u(j,l)
         dv         = v(j,lp)  - v(j,l)
         dw         = w(j,lp)  - w(j,l)
         dei        = ei(j,lp) - ei(j,l)
c        f(j,l,1)   = 0.
         f(j,l,2)   = (s1av2*du + s4av2*dv + s5av2*dw)*vmud2
         f(j,l,3)   = (s4av2*du + s2av2*dv + s6av2*dw)*vmud2
         f(j,l,4)   = (s5av2*du + s6av2*dv + s3av2*dw)*vmud2
         f(j,l,5)   = 0.5*(f(j,l,2)*uav2 + f(j,l,3)*vav2 + f(j,l,4)*wav2
     &                                   + s0av2*dei*gamkap)
   20    continue
c
c
         deallocate(s0,s1,s2,s3,s4,s5,s6,u,v,w,ei)
      return
      end


      subroutine vfljk ( l,
     &                   rey,pr,prt,
     &                   q,vgamma,vmul,vmut,
     &                   xx,xy,xz,yx,yy,yz,
     &                   fj,fk,jmax,kmax,lmax,nd,mdim )
c
c   compute a jk-plane of jk (mixed derivative) viscous fluxes.
c   get fluxes at grid points j+1/2 in fj, k+1/2 in fk.
c
c   called by: rvjk2
c
c   calls: none
c
c   181 ops/pt.
c
c*******************************************************************

      integer l,jmax,kmax,lmax,nd,mdim
      
      real rey,pr,prt,vgamma
      
      real q(jmax,kmax,lmax,nd),vmul(jmax,kmax,lmax)
      real vmut(jmax,kmax,lmax)
      real xx(jmax,kmax,lmax),xy(jmax,kmax,lmax),xz(jmax,kmax,lmax),
     &     yx(jmax,kmax,lmax),yy(jmax,kmax,lmax),yz(jmax,kmax,lmax)
      real fj(mdim,mdim,5),fk(mdim,mdim,5)

!      real s0(mdim,mdim),s1(mdim,mdim),s2(mdim,mdim),
!     &     s3(mdim,mdim),s4(mdim,mdim),s5(mdim,mdim),
!     &     s6(mdim,mdim),s7(mdim,mdim),s8(mdim,mdim),
!     &     s9(mdim,mdim)
!      real u(mdim,mdim),v(mdim,mdim),w(mdim,mdim),ei(mdim,mdim)

      real, allocatable :: s0(:,:),s1(:,:),s2(:,:)
      real, allocatable :: s3(:,:),s4(:,:),s5(:,:),s6(:,:),s7(:,:)
      real, allocatable :: s8(:,:),s9(:,:)
      real, allocatable :: u(:,:),v(:,:),w(:,:),ei(:,:)

c..   local variables

      integer k,j
      real r3,n23,prtr,dre,dprre
      real drho,e0,v2,rv
      real gamav2,vmula2,vmuta2,vmud2,gamkap
      real uav2,vav2,wav2
      real s0av2,s1av2,s2av2,s3av2,s4av2,s5av2
      real s6av2,s7av2,s8av2,s9av2
      real du,dv,dw,dei

      allocate(s0(mdim,mdim),s1(mdim,mdim),s2(mdim,mdim))
      allocate(s3(mdim,mdim),s4(mdim,mdim),s5(mdim,mdim),s6(mdim,mdim))
      allocate(s7(mdim,mdim),s8(mdim,mdim),s9(mdim,mdim))
      allocate(u(mdim,mdim),v(mdim,mdim),w(mdim,mdim),ei(mdim,mdim))

      r3=1.d0/3.d0
      n23=-2.d0/3.d0
      prtr  = pr/prt
      dre   = 0.0625/rey
      dprre = 0.0625/(pr*rey)
      do k = 1,kmax
      do j = 1,jmax
         drho       = 1./q(j,k,l,1)
         u(j,k)     = q(j,k,l,2)*drho
         v(j,k)     = q(j,k,l,3)*drho
         w(j,k)     = q(j,k,l,4)*drho
         e0         = q(j,k,l,5)*drho
         v2         = 0.5*(u(j,k)**2 + v(j,k)**2 + w(j,k)**2)
         ei(j,k)    = e0 - v2
         rv         = 1./q(j,k,l,6)
         r3rv       = r3*rv
         s0(j,k)    = (xx(j,k,l)*yx(j,k,l) + xy(j,k,l)*yy(j,k,l)
     &               + xz(j,k,l)*yz(j,k,l))*rv
         s1(j,k)    = s0(j,k) + r3rv*xx(j,k,l)*yx(j,k,l)
         s2(j,k)    = s0(j,k) + r3rv*xy(j,k,l)*yy(j,k,l)
         s3(j,k)    = s0(j,k) + r3rv*xz(j,k,l)*yz(j,k,l)
         s4(j,k)    = (n23*xx(j,k,l)*yy(j,k,l) + xy(j,k,l)*yx(j,k,l))*rv
         s5(j,k)    = (n23*xx(j,k,l)*yz(j,k,l) + xz(j,k,l)*yx(j,k,l))*rv
         s6(j,k)    = (xx(j,k,l)*yy(j,k,l) + n23*xy(j,k,l)*yx(j,k,l))*rv
         s7(j,k)    = (n23*xy(j,k,l)*yz(j,k,l) + xz(j,k,l)*yy(j,k,l))*rv
         s8(j,k)    = (xx(j,k,l)*yz(j,k,l) + n23*xz(j,k,l)*yx(j,k,l))*rv
         s9(j,k)    = (xy(j,k,l)*yz(j,k,l) + n23*xz(j,k,l)*yy(j,k,l))*rv
      end do
      end do
c
c   compute j-direction fluxes.
c
      do k = 2,kmax-1
         kp       = k+1
         km       = k-1
         do j = 1,jmax-1
            jp       = j+1
            gamav2   = vgamma + vgamma
            vmula2   = vmul(jp,k,l) + vmul(j,k,l)
            vmuta2   = vmut(jp,k,l) + vmut(j,k,l)
            vmud2    = (vmula2 + vmuta2)*dre
            gamkap   = gamav2*(vmula2 + prtr*vmuta2)*dprre
            uav2     = u(jp,k) + u(j,k)
            vav2     = v(jp,k) + v(j,k)
            wav2     = w(jp,k) + w(j,k)
            s0av2    = s0(jp,k) + s0(j,k)
            s1av2    = s1(jp,k) + s1(j,k)
            s2av2    = s2(jp,k) + s2(j,k)
            s3av2    = s3(jp,k) + s3(j,k)
            s4av2    = s4(jp,k) + s4(j,k)
            s5av2    = s5(jp,k) + s5(j,k)
            s6av2    = s6(jp,k) + s6(j,k)
            s7av2    = s7(jp,k) + s7(j,k)
            s8av2    = s8(jp,k) + s8(j,k)
            s9av2    = s9(jp,k) + s9(j,k)
            du       = u(jp,kp) - u(jp,km) + u(j,kp) - u(j,km)
            dv       = v(jp,kp) - v(jp,km) + v(j,kp) - v(j,km)
            dw       = w(jp,kp) - w(jp,km) + w(j,kp) - w(j,km)
            dei      = ei(jp,kp) - ei(jp,km) + ei(j,kp) - ei(j,km)
c           fj(j,k,1)= 0.
            fj(j,k,2)= (s1av2*du + s4av2*dv + s5av2*dw)*vmud2
            fj(j,k,3)= (s6av2*du + s2av2*dv + s7av2*dw)*vmud2
            fj(j,k,4)= (s8av2*du + s9av2*dv + s3av2*dw)*vmud2
            fj(j,k,5)= 0.5*(fj(j,k,2)*uav2 + fj(j,k,3)*vav2
     &                    + fj(j,k,4)*wav2 + s0av2*dei*gamkap)
         enddo
      enddo
c
c   compute k-direction fluxes.
c
      do k = 1,kmax-1
         kp       = k+1
         do j = 2,jmax-1
            jp       = j+1
            jm       = j-1
            gamav2   = vgamma + vgamma
            vmula2   = vmul(j,kp,l) + vmul(j,k,l)
            vmuta2   = vmut(j,kp,l) + vmut(j,k,l)
            vmud2    = (vmula2 + vmuta2)*dre
            gamkap   = gamav2*(vmula2 + prtr*vmuta2)*dprre
            uav2     = u(j,kp) + u(j,k)
            vav2     = v(j,kp) + v(j,k)
            wav2     = w(j,kp) + w(j,k)
            s0av2    = s0(j,kp) + s0(j,k)
            s1av2    = s1(j,kp) + s1(j,k)
            s2av2    = s2(j,kp) + s2(j,k)
            s3av2    = s3(j,kp) + s3(j,k)
            s4av2    = s4(j,kp) + s4(j,k)
            s5av2    = s5(j,kp) + s5(j,k)
            s6av2    = s6(j,kp) + s6(j,k)
            s7av2    = s7(j,kp) + s7(j,k)
            s8av2    = s8(j,kp) + s8(j,k)
            s9av2    = s9(j,kp) + s9(j,k)
            du       = u(jp,kp) - u(jm,kp) + u(jp,k) - u(jm,k)
            dv       = v(jp,kp) - v(jm,kp) + v(jp,k) - v(jm,k)
            dw       = w(jp,kp) - w(jm,kp) + w(jp,k) - w(jm,k)
            dei      = ei(jp,kp) - ei(jm,kp) + ei(jp,k) - ei(jm,k)
c           fk(j,k,1)= 0.
            fk(j,k,2)= (s1av2*du + s6av2*dv + s8av2*dw)*vmud2
            fk(j,k,3)= (s4av2*du + s2av2*dv + s9av2*dw)*vmud2
            fk(j,k,4)= (s5av2*du + s7av2*dv + s3av2*dw)*vmud2
            fk(j,k,5)= 0.5*(fk(j,k,2)*uav2 + fk(j,k,3)*vav2
     &                    + fk(j,k,4)*wav2 + s0av2*dei*gamkap)
         enddo
      enddo
c
      deallocate(s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,u,v,w,ei)
      return
      end


      subroutine vfljl ( k,
     &                   rey,pr,prt,
     &                   q,vgamma,vmul,vmut,
     &                   xx,xy,xz,zx,zy,zz,
     &                   fj,fl,jmax,kmax,lmax,nd,mdim )
c
c   compute a jl-plane of jl (mixed derivative) viscous fluxes.
c   get fluxes at grid points j+1/2 in fj, l+1/2 in fl.
c
c   called by: rvjl2
c
c   calls: none
c
c   181 ops/pt.
c
c*******************************************************************
      
      integer k,jmax,kmax,lmax,mdim,nd

      real rey,pr,prt,vgamma
      real q(jmax,kmax,lmax,nd),vmul(jmax,kmax,lmax)
      real vmut(jmax,kmax,lmax)

      real xx(jmax,kmax,lmax),xy(jmax,kmax,lmax),xz(jmax,kmax,lmax)
      real zx(jmax,kmax,lmax),zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)

      real fj(mdim,mdim,5),fl(mdim,mdim,5)
!      real s0(mdim,mdim),s1(mdim,mdim),s2(mdim,mdim),
!     &     s3(mdim,mdim),s4(mdim,mdim),s5(mdim,mdim),
!     &     s6(mdim,mdim),s7(mdim,mdim),s8(mdim,mdim),
!     &     s9(mdim,mdim)
!      real u(mdim,mdim),v(mdim,mdim),w(mdim,mdim),ei(mdim,mdim)
      
      real, allocatable :: s0(:,:),s1(:,:),s2(:,:)
      real, allocatable :: s3(:,:),s4(:,:),s5(:,:),s6(:,:),s7(:,:)
      real, allocatable :: s8(:,:),s9(:,:)
      real, allocatable :: u(:,:),v(:,:),w(:,:),ei(:,:)


c..   local variables

      integer j,l
      real r3,n23,prtr,dre,dprre
      real drho,e0,v2,rv
      real gamav2,vmula2,vmuta2,vmud2,gamkap
      real uav2,vav2,wav2
      real s0av2,s1av2,s2av2,s3av2,s4av2,s5av2
      real s6av2,s7av2,s8av2,s9av2
      real du,dv,dw,dei

      allocate(s0(mdim,mdim),s1(mdim,mdim),s2(mdim,mdim))
      allocate(s3(mdim,mdim),s4(mdim,mdim),s5(mdim,mdim),s6(mdim,mdim))
      allocate(s7(mdim,mdim),s8(mdim,mdim),s9(mdim,mdim))
      allocate(u(mdim,mdim),v(mdim,mdim),w(mdim,mdim),ei(mdim,mdim))

       r3=1.d0/3.d0
       n23=-2.d0/3.d0
       prtr  = pr/prt
       dre   = 0.0625/rey
       dprre = 0.0625/(pr*rey)
c     
      do l = 1,lmax
      do j = 1,jmax
         drho       = 1./q(j,k,l,1)
         u(j,l)     = q(j,k,l,2)*drho
         v(j,l)     = q(j,k,l,3)*drho
         w(j,l)     = q(j,k,l,4)*drho
         e0         = q(j,k,l,5)*drho
         v2         = 0.5*(u(j,l)**2 + v(j,l)**2 + w(j,l)**2)
         ei(j,l)    = e0 - v2
         rv         = 1./q(j,k,l,6)
         r3rv       = r3*rv
         s0(j,l)    = (xx(j,k,l)*zx(j,k,l) + xy(j,k,l)*zy(j,k,l)
     &               + xz(j,k,l)*zz(j,k,l))*rv
         s1(j,l)    = s0(j,l) + r3rv*xx(j,k,l)*zx(j,k,l)
         s2(j,l)    = s0(j,l) + r3rv*xy(j,k,l)*zy(j,k,l)
         s3(j,l)    = s0(j,l) + r3rv*xz(j,k,l)*zz(j,k,l)
         s4(j,l)    = (n23*xx(j,k,l)*zy(j,k,l) + xy(j,k,l)*zx(j,k,l))*rv
         s5(j,l)    = (n23*xx(j,k,l)*zz(j,k,l) + xz(j,k,l)*zx(j,k,l))*rv
         s6(j,l)    = (xx(j,k,l)*zy(j,k,l) + n23*xy(j,k,l)*zx(j,k,l))*rv
         s7(j,l)    = (n23*xy(j,k,l)*zz(j,k,l) + xz(j,k,l)*zy(j,k,l))*rv
         s8(j,l)    = (xx(j,k,l)*zz(j,k,l) + n23*xz(j,k,l)*zx(j,k,l))*rv
         s9(j,l)    = (xy(j,k,l)*zz(j,k,l) + n23*xz(j,k,l)*zy(j,k,l))*rv
      enddo
      enddo
c
c   j-direction fluxes.
c
      do l = 2,lmax-1
         lp       = l+1
         lm       = l-1
         do j = 1,jmax-1
            jp       = j+1
            gamav2   = vgamma + vgamma
            vmula2   = vmul(jp,k,l) + vmul(j,k,l)
            vmuta2   = vmut(jp,k,l) + vmut(j,k,l)
            vmud2    = (vmula2 + vmuta2)*dre
            gamkap   = gamav2*(vmula2 + prtr*vmuta2)*dprre
            uav2     = u(jp,l) + u(j,l)
            vav2     = v(jp,l) + v(j,l)
            wav2     = w(jp,l) + w(j,l)
            s0av2    = s0(jp,l) + s0(j,l)
            s1av2    = s1(jp,l) + s1(j,l)
            s2av2    = s2(jp,l) + s2(j,l)
            s3av2    = s3(jp,l) + s3(j,l)
            s4av2    = s4(jp,l) + s4(j,l)
            s5av2    = s5(jp,l) + s5(j,l)
            s6av2    = s6(jp,l) + s6(j,l)
            s7av2    = s7(jp,l) + s7(j,l)
            s8av2    = s8(jp,l) + s8(j,l)
            s9av2    = s9(jp,l) + s9(j,l)
            du       = u(jp,lp) - u(jp,lm) + u(j,lp) - u(j,lm)
            dv       = v(jp,lp) - v(jp,lm) + v(j,lp) - v(j,lm)
            dw       = w(jp,lp) - w(jp,lm) + w(j,lp) - w(j,lm)
            dei      = ei(jp,lp) - ei(jp,lm) + ei(j,lp) - ei(j,lm)
c           fj(j,l,1)= 0.
            fj(j,l,2)= (s1av2*du + s4av2*dv + s5av2*dw)*vmud2
            fj(j,l,3)= (s6av2*du + s2av2*dv + s7av2*dw)*vmud2
            fj(j,l,4)= (s8av2*du + s9av2*dv + s3av2*dw)*vmud2
            fj(j,l,5)= 0.5*(fj(j,l,2)*uav2 + fj(j,l,3)*vav2
     &                    + fj(j,l,4)*wav2 + s0av2*dei*gamkap)
         enddo
      enddo
c
c   l-direction fluxes.
c
      do l = 1,lmax-1
         lp       = l+1
         do j = 2,jmax-1
            jp       = j+1
            jm       = j-1
            gamav2   = vgamma + vgamma
            vmula2   = vmul(j,k,lp) + vmul(j,k,l)
            vmuta2   = vmut(j,k,lp) + vmut(j,k,l)
            vmud2    = (vmula2 + vmuta2)*dre
            gamkap   = gamav2*(vmula2 + prtr*vmuta2)*dprre
            uav2     = u(j,lp) + u(j,l)
            vav2     = v(j,lp) + v(j,l)
            wav2     = w(j,lp) + w(j,l)
            s0av2    = s0(j,lp) + s0(j,l)
            s1av2    = s1(j,lp) + s1(j,l)
            s2av2    = s2(j,lp) + s2(j,l)
            s3av2    = s3(j,lp) + s3(j,l)
            s4av2    = s4(j,lp) + s4(j,l)
            s5av2    = s5(j,lp) + s5(j,l)
            s6av2    = s6(j,lp) + s6(j,l)
            s7av2    = s7(j,lp) + s7(j,l)
            s8av2    = s8(j,lp) + s8(j,l)
            s9av2    = s9(j,lp) + s9(j,l)
            du       = u(jp,lp) - u(jm,lp) + u(jp,l) - u(jm,l)
            dv       = v(jp,lp) - v(jm,lp) + v(jp,l) - v(jm,l)
            dw       = w(jp,lp) - w(jm,lp) + w(jp,l) - w(jm,l)
            dei      = ei(jp,lp) - ei(jm,lp) + ei(jp,l) - ei(jm,l)
c           fl(j,l,1)= 0.
            fl(j,l,2)= (s1av2*du + s6av2*dv + s8av2*dw)*vmud2
            fl(j,l,3)= (s4av2*du + s2av2*dv + s9av2*dw)*vmud2
            fl(j,l,4)= (s5av2*du + s7av2*dv + s3av2*dw)*vmud2
            fl(j,l,5)= 0.5*(fl(j,l,2)*uav2 + fl(j,l,3)*vav2
     &                    + fl(j,l,4)*wav2 + s0av2*dei*gamkap)
         enddo
      enddo
c
      deallocate(s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,u,v,w,ei)
      return
      end


      subroutine vflkl ( l,
     &                   rey,pr,prt,
     &                   q,vgamma,vmul,vmut,
     &                   yx,yy,yz,zx,zy,zz,
     &                   fk,fl,jmax,kmax,lmax,nd,mdim )
c
c   compute a kl-plane of kl (mixed derivative) viscous fluxes.
c   get k fluxes at k+1/2 in fk, l flux differences at l in fl.
c
c   note that we are doing a bunch of duplicate work because we want to
c   parallelize in l (and vectorize in j).  hopefully it is worth it.
c
c   called by: rvkl2
c
c   calls: none
c
c*******************************************************************

      integer l,jmax,kmax,lmax,mdim,nd

      real rey,pr,prt,vgamma
      real q(jmax,kmax,lmax,nd),vmul(jmax,kmax,lmax)
      real vmut(jmax,kmax,lmax)

      real yx(jmax,kmax,lmax),yy(jmax,kmax,lmax),yz(jmax,kmax,lmax)
      real zx(jmax,kmax,lmax),zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)

      real fk(mdim,mdim,5),fl(mdim,mdim,5)
!
!      real s0(mdim,mdim,3),s1(mdim,mdim,3),s2(mdim,mdim,3),
!     &     s3(mdim,mdim,3),s4(mdim,mdim,3),s5(mdim,mdim,3),
!     &     s6(mdim,mdim,3),s7(mdim,mdim,3),s8(mdim,mdim,3),
!     &     s9(mdim,mdim,3)
!      real u(mdim,mdim,3),v(mdim,mdim,3),w(mdim,mdim,3),ei(mdim,mdim,3)
 
      real, allocatable :: s0(:,:,:),s1(:,:,:),s2(:,:,:)
      real, allocatable :: s3(:,:,:),s4(:,:,:),s5(:,:,:)
      real, allocatable :: s6(:,:,:),s7(:,:,:),s8(:,:,:),s9(:,:,:)
      real, allocatable :: u(:,:,:),v(:,:,:),w(:,:,:),ei(:,:,:)

c..   local variables

      integer lx,k,j,lp,lm,lxp,lxm
      real r3,n23,prtr,dre,dprre
      real drho,e0,v2,rv
      real gamav2,vmula2,vmuta2,vmud2,gamkap
      real uav2,vav2,wav2
      real s0av2,s1av2,s2av2,s3av2,s4av2,s5av2
      real s6av2,s7av2,s8av2,s9av2
      real du,dv,dw,dei
      real fl1,fl2,fl3,fl4,fl5
      real fl1m,fl2m,fl3m,fl4m,fl5m

      allocate(s0(mdim,mdim,3),s1(mdim,mdim,3),s2(mdim,mdim,3))
      allocate(s3(mdim,mdim,3),s4(mdim,mdim,3),s5(mdim,mdim,3))
      allocate(s6(mdim,mdim,3),s7(mdim,mdim,3),s8(mdim,mdim,3))
      allocate(s9(mdim,mdim,3))
      allocate(u(mdim,mdim,3),v(mdim,mdim,3),w(mdim,mdim,3),
     &     ei(mdim,mdim,3))

      lp    = l+1
      lm    = l-1
c     
      r3=1.d0/3.d0
      n23=-2.d0/3.d0

      prtr  = pr/prt
      dre   = 0.0625/rey
      dprre = 0.0625/(pr*rey)
c
      do lx = 1,3
         if (lx.eq.1) then
            ll         = lm
         else if (lx.eq.2) then
            ll         = l
         else
            ll         = lp
         endif
        do k = 1,kmax
        do j = 1,jmax
            drho       = 1./q(j,k,ll,1)
            u(j,k,lx)  = q(j,k,ll,2)*drho
            v(j,k,lx)  = q(j,k,ll,3)*drho
            w(j,k,lx)  = q(j,k,ll,4)*drho
            e0         = q(j,k,ll,5)*drho
            v2         = 0.5*(u(j,k,lx)**2 + v(j,k,lx)**2
     &                                     + w(j,k,lx)**2)
            ei(j,k,lx) = e0 - v2
            rv         = 1./q(j,k,ll,6)
            r3rv       = r3*rv
            s0(j,k,lx) = ( yx(j,k,ll)*zx(j,k,ll) + yy(j,k,ll)*zy(j,k,ll)
     &                   + yz(j,k,ll)*zz(j,k,ll) )*rv
            s1(j,k,lx) = s0(j,k,lx) + r3rv*yx(j,k,ll)*zx(j,k,ll)
            s2(j,k,lx) = s0(j,k,lx) + r3rv*yy(j,k,ll)*zy(j,k,ll)
            s3(j,k,lx) = s0(j,k,lx) + r3rv*yz(j,k,ll)*zz(j,k,ll)
            s4(j,k,lx) = ( n23*yx(j,k,ll)*zy(j,k,ll)
     &                       + yy(j,k,ll)*zx(j,k,ll) )*rv 
            s5(j,k,lx) = ( n23*yx(j,k,ll)*zz(j,k,ll)
     &                       + yz(j,k,ll)*zx(j,k,ll) )*rv
            s6(j,k,lx) = (     yx(j,k,ll)*zy(j,k,ll)
     &                   + n23*yy(j,k,ll)*zx(j,k,ll) )*rv
            s7(j,k,lx) = ( n23*yy(j,k,ll)*zz(j,k,ll)
     &                       + yz(j,k,ll)*zy(j,k,ll) )*rv
            s8(j,k,lx) = (     yx(j,k,ll)*zz(j,k,ll)
     &                   + n23*yz(j,k,ll)*zx(j,k,ll) )*rv
            s9(j,k,lx) = (     yy(j,k,ll)*zz(j,k,ll)
     &                   + n23*yz(j,k,ll)*zy(j,k,ll) )*rv
         end do
         end do
      end do
c
      lx       = 2
      lxp      = 3
      lxm      = 1
c
c   k-direction fluxes.
c
      do k = 1,kmax-1
         kp       = k+1
         do j = 1,jmax
            gamav2   = vgamma + vgamma
            vmula2   = vmul(j,kp,l) + vmul(j,k,l)
            vmuta2   = vmut(j,kp,l) + vmut(j,k,l)
            vmud2    = (vmula2 + vmuta2)*dre
            gamkap   = gamav2*(vmula2 + prtr*vmuta2)*dprre
            uav2     = u(j,kp,lx) + u(j,k,lx)
            vav2     = v(j,kp,lx) + v(j,k,lx)
            wav2     = w(j,kp,lx) + w(j,k,lx)
            s0av2    = s0(j,kp,lx) + s0(j,k,lx)
            s1av2    = s1(j,kp,lx) + s1(j,k,lx)
            s2av2    = s2(j,kp,lx) + s2(j,k,lx)
            s3av2    = s3(j,kp,lx) + s3(j,k,lx)
            s4av2    = s4(j,kp,lx) + s4(j,k,lx)
            s5av2    = s5(j,kp,lx) + s5(j,k,lx)
            s6av2    = s6(j,kp,lx) + s6(j,k,lx)
            s7av2    = s7(j,kp,lx) + s7(j,k,lx)
            s8av2    = s8(j,kp,lx) + s8(j,k,lx)
            s9av2    = s9(j,kp,lx) + s9(j,k,lx)
            du       = u(j,kp,lxp) - u(j,kp,lxm)
     &               + u(j,k,lxp)  - u(j,k,lxm) 
            dv       = v(j,kp,lxp) - v(j,kp,lxm)
     &               + v(j,k,lxp)  - v(j,k,lxm)
            dw       = w(j,kp,lxp) - w(j,kp,lxm)
     &               + w(j,k,lxp)  - w(j,k,lxm)
            dei      = ei(j,kp,lxp) - ei(j,kp,lxm)
     &               + ei(j,k,lxp)  - ei(j,k,lxm)
c           fk(j,k,1)= 0.
            fk(j,k,2)= (s1av2*du + s4av2*dv + s5av2*dw)*vmud2
            fk(j,k,3)= (s6av2*du + s2av2*dv + s7av2*dw)*vmud2
            fk(j,k,4)= (s8av2*du + s9av2*dv + s3av2*dw)*vmud2
            fk(j,k,5)= 0.5*(fk(j,k,2)*uav2 + fk(j,k,3)*vav2
     &                    + fk(j,k,4)*wav2 + s0av2*dei*gamkap)
         enddo
      enddo
c
c   l-direction fluxes.
c
      do k = 2,kmax-1
         kp       = k+1
         km       = k-1
         do j = 1,jmax
c
c   l-direction fluxes at l+1/2.
c
            gamav2   = vgamma + vgamma
            vmula2   = vmul(j,k,lp) + vmul(j,k,l)
            vmuta2   = vmut(j,k,lp) + vmut(j,k,l)
            vmud2    = (vmula2 + vmuta2)*dre
            gamkap   = gamav2*(vmula2 + prtr*vmuta2)*dprre
            uav2     = u(j,k,lxp) + u(j,k,lx)
            vav2     = v(j,k,lxp) + v(j,k,lx)
            wav2     = w(j,k,lxp) + w(j,k,lx)
            s0av2    = s0(j,k,lxp) + s0(j,k,lx)
            s1av2    = s1(j,k,lxp) + s1(j,k,lx)
            s2av2    = s2(j,k,lxp) + s2(j,k,lx)
            s3av2    = s3(j,k,lxp) + s3(j,k,lx)
            s4av2    = s4(j,k,lxp) + s4(j,k,lx)
            s5av2    = s5(j,k,lxp) + s5(j,k,lx)
            s6av2    = s6(j,k,lxp) + s6(j,k,lx)
            s7av2    = s7(j,k,lxp) + s7(j,k,lx)
            s8av2    = s8(j,k,lxp) + s8(j,k,lx)
            s9av2    = s9(j,k,lxp) + s9(j,k,lx)
            du       = u(j,kp,lxp) - u(j,km,lxp)
     &               + u(j,kp,lx)  - u(j,km,lx)
            dv       = v(j,kp,lxp) - v(j,km,lxp)
     &               + v(j,kp,lx)  - v(j,km,lx)
            dw       = w(j,kp,lxp) - w(j,km,lxp)
     &               + w(j,kp,lx)  - w(j,km,lx)
            dei      = ei(j,kp,lxp) - ei(j,km,lxp)
     &               + ei(j,kp,lx)  - ei(j,km,lx)
c           fl1      = 0.
            fl2      = (s1av2*du + s6av2*dv + s8av2*dw)*vmud2
            fl3      = (s4av2*du + s2av2*dv + s9av2*dw)*vmud2
            fl4      = (s5av2*du + s7av2*dv + s3av2*dw)*vmud2
            fl5      = 0.5*(fl2*uav2 + fl3*vav2
     &                    + fl4*wav2 + s0av2*dei*gamkap)
c
c   l-direction fluxes at l-1/2.
c
            gamav2   = vgamma + vgamma
            vmula2   = vmul(j,k,l) + vmul(j,k,lm)
            vmuta2   = vmut(j,k,l) + vmut(j,k,lm)
            vmud2    = (vmula2 + vmuta2)*dre
            gamkap   = gamav2*(vmula2 + prtr*vmuta2)*dprre
            uav2     = u(j,k,lx) + u(j,k,lxm)
            vav2     = v(j,k,lx) + v(j,k,lxm)
            wav2     = w(j,k,lx) + w(j,k,lxm)
            s0av2    = s0(j,k,lx) + s0(j,k,lxm)
            s1av2    = s1(j,k,lx) + s1(j,k,lxm)
            s2av2    = s2(j,k,lx) + s2(j,k,lxm)
            s3av2    = s3(j,k,lx) + s3(j,k,lxm)
            s4av2    = s4(j,k,lx) + s4(j,k,lxm)
            s5av2    = s5(j,k,lx) + s5(j,k,lxm)
            s6av2    = s6(j,k,lx) + s6(j,k,lxm)
            s7av2    = s7(j,k,lx) + s7(j,k,lxm)
            s8av2    = s8(j,k,lx) + s8(j,k,lxm)
            s9av2    = s9(j,k,lx) + s9(j,k,lxm)
            du       = u(j,kp,lx)  - u(j,km,lx)
     &               + u(j,kp,lxm) - u(j,km,lxm)
            dv       = v(j,kp,lx)  - v(j,km,lx)
     &               + v(j,kp,lxm) - v(j,km,lxm)
            dw       = w(j,kp,lx)  - w(j,km,lx)
     &               + w(j,kp,lxm) - w(j,km,lxm)
            dei      = ei(j,kp,lx)  - ei(j,km,lx)
     &               + ei(j,kp,lxm) - ei(j,km,lxm)
c           fl1m     = 0.
            fl2m     = (s1av2*du + s6av2*dv + s8av2*dw)*vmud2
            fl3m     = (s4av2*du + s2av2*dv + s9av2*dw)*vmud2
            fl4m     = (s5av2*du + s7av2*dv + s3av2*dw)*vmud2
            fl5m     = 0.5*(fl2m*uav2 + fl3m*vav2
     &                    + fl4m*wav2 + s0av2*dei*gamkap)
c
c   flux differences.
c
c           fl(j,k,1)= fl1 - fl1m
            fl(j,k,2)= fl2 - fl2m
            fl(j,k,3)= fl3 - fl3m
            fl(j,k,4)= fl4 - fl4m
            fl(j,k,5)= fl5 - fl5m
         enddo
      enddo
c
      deallocate(s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,u,v,w,ei)
      return
      end
