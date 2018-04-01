c***********************************************************************
      subroutine visrhs(x,y,z,iblank,vnu,vnu0,vmul,turmu,q,s,xx,xy,xz,
     $			ug,yx,yy,yz,vg,zx,zy,zz,wg,tscale,kkr,kkp,grbc)
c
c  compute the viscous rhs 
c
c***********************************************************************
      
      use params_global
      use bcparam
      implicit none

      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)
      real q(jmax,kmax,lmax,nd), s(jmax,kmax,lmax,nv), 
     $     turmu(jmax,kmax,lmax),tscale(jmax,kmax,lmax)
      real xx(jmax,kmax,lmax), xy(jmax,kmax,lmax),xz(jmax,kmax,lmax),
     $     yx(jmax,kmax,lmax), yy(jmax,kmax,lmax),yz(jmax,kmax,lmax),
     $     zx(jmax,kmax,lmax), zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)
      real vnu(jmax,kmax,lmax),vmul(jmax,kmax,lmax)
      real vnu0(jmax,kmax,lmax)
      integer iblank(jmax,kmax,lmax)

      real ug(jmax,kmax,lmax),vg(jmax,kmax,lmax),wg(jmax,kmax,lmax)
      integer kkr(kmax),kkp(kmax)
      type(bc_t) :: grbc
c
c..   local variables

      real,allocatable :: f2(:),f3(:),f4(:),f5(:)
      real,allocatable :: s0(:),s1(:),s2(:),s3(:),s4(:),
     <    s5(:),s6(:),u(:),v(:),w(:),e(:),rr(:)

      integer ka,kb,j,k,l,l1

      real gkpr,prtr,dre,c2b,c2bp,ra,vmue,turm,vnu1,gkap,rj,tt
      real t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15
      real t16,du,dv,dw,dei,sqrt

      allocate(f2(mdim),f3(mdim),f4(mdim),f5(mdim))
      allocate(s0(mdim),s1(mdim),s2(mdim),s3(mdim),s4(mdim),
     <    s5(mdim),s6(mdim),u(mdim),v(mdim),w(mdim),e(mdim),rr(mdim))


c***********************************

      ka=  2-ksym
      kb = km + ksym
c
      if(.not. lamin) then
         if (iturb.eq.1) then
            call vmutur( x,y,z,q,s,turmu,xx,xy,xz,yx,yy,yz,
     <           zx,zy,zz,ug,vg,wg,kkr,kkp)
         else
            call vmu_sa(x,y,z,q,turmu,xx,xy,xz,yx,yy,yz,zx,zy,zz,
     <           ug,vg,wg,vnu,vnu0,vmul,tscale,iblank,grbc)
c
         endif
      endif

      gkpr = gamma/pr
      prtr = pr/0.9
      dre  = .5/rey
      c2b  =198.6/tinf
      c2bp = c2b +1.
      do 10 j = 2,jm
      do 10 k = ka,kb
        do 20 l = 1,lmax
          ra    = 1./q(j,k,l,1)
          u(l)  = q(j,k,l,2)*ra
          v(l)  = q(j,k,l,3)*ra
          w(l)  = q(j,k,l,4)*ra
          e(l)  = q(j,k,l,5)*ra-.5*(u(l)**2+v(l)**2+w(l)**2)
          tt    = ggm1*e(l)
          vmue  = c2bp*tt*sqrt(tt)/(c2b + tt)
          turm  = turmu(j,k,l)
          vnu1   = vmue+turm
          gkap  = vmue+prtr*turm
          rj    = 1./q(j,k,l,6)
          s0(l) = (zx(j,k,l)**2+zy(j,k,l)**2+zz(j,k,l)**2)*rj
          s1(l) = (s0(l)+zx(j,k,l)**2/3.*rj)*vnu1*dre
          s2(l) = (s0(l)+zy(j,k,l)**2/3.*rj)*vnu1*dre
          s3(l) = (s0(l)+zz(j,k,l)**2/3.*rj)*vnu1*dre
          s4(l) = (zx(j,k,l)*zy(j,k,l)/3.*rj)*vnu1*dre
          s5(l) = (zx(j,k,l)*zz(j,k,l)/3.*rj)*vnu1*dre
          s6(l) = (zy(j,k,l)*zz(j,k,l)/3.*rj)*vnu1*dre
          s0(l) = s0(l)*gkpr*gkap*dre
   20   continue
        do 30 l = 1,lm
          l1    = l+1
          t1    = s1(l1)+s1(l)
          t2    = s2(l1)+s2(l)
          t3    = s3(l1)+s3(l)
          t4    = s4(l1)+s4(l)
          t5    = s5(l1)+s5(l)
          t6    = s6(l1)+s6(l)
          t7    = u(l1)*s1(l1)+u(l)*s1(l)
          t8    = v(l1)*s2(l1)+v(l)*s2(l)
          t9    = w(l1)*s3(l1)+w(l)*s3(l)
          t10   = u(l1)*s4(l1)+u(l)*s4(l)
          t11   = u(l1)*s5(l1)+u(l)*s5(l)
          t12   = v(l1)*s4(l1)+v(l)*s4(l)
          t13   = v(l1)*s6(l1)+v(l)*s6(l)
          t14   = w(l1)*s5(l1)+w(l)*s5(l)
          t15   = w(l1)*s6(l1)+w(l)*s6(l)
          t16   = s0(l1)+s0(l)
          du    = u(l1)-u(l)
          dv    = v(l1)-v(l)
          dw    = w(l1)-w(l)
          dei   = e(l1)-e(l)
          f2(l) = t1*du+t4*dv+t5*dw
          f3(l) = t4*du+t2*dv+t6*dw
          f4(l) = t5*du+t6*dv+t3*dw
          f5(l) = (t7+t12+t14)*du+(t8+t10+t15)*dv+(t9+t11+t13)*dw
     *             +t16*dei
   30   continue
c
        do 40 l = 2,lm
c          s(j,k,l,1) = s(j,k,l,1) + 0.
Casitav...using max(iblank,0) instead of iblank
          s(j,k,l,2) = s(j,k,l,2) + (f2(l)-f2(l-1))*max(iblank(j,k,l),0)
          s(j,k,l,3) = s(j,k,l,3) + (f3(l)-f3(l-1))*max(iblank(j,k,l),0)
          s(j,k,l,4) = s(j,k,l,4) + (f4(l)-f4(l-1))*max(iblank(j,k,l),0)
          s(j,k,l,5) = s(j,k,l,5) + (f5(l)-f5(l-1))*max(iblank(j,k,l),0)
   40   continue
   10 continue
c
      return
      end


c***********************************************************************
      subroutine vmutur(x,y,z,q,s,turmu,xx,xy,xz,yx,yy,yz,
     $                  zx,zy,zz,ug,vg,wg,kkr,kkp)
c
c  turbulent eddy viscosity. model is algebraic model of baldwin
c  and lomax
c
c***********************************************************************

      use params_global
      
      implicit none
      
      real q(jmax,kmax,lmax,nd),s(jmax,kmax,lmax,nv),
     &     turmu(jmax,kmax,lmax)
      real x(jmax,kmax,lmax), y(jmax,kmax,lmax) ,z(jmax,kmax,lmax)
      real xx(jmax,kmax,lmax), xy(jmax,kmax,lmax) ,xz(jmax,kmax,lmax),
     $     yx(jmax,kmax,lmax), yy(jmax,kmax,lmax) ,yz(jmax,kmax,lmax),
     $     zx(jmax,kmax,lmax), zy(jmax,kmax,lmax) ,zz(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax), vg(jmax,kmax,lmax), wg(jmax,kmax,lmax)
      integer kkr(kmax),kkp(kmax)

c..   local variables

      real,allocatable :: fmax(:),smax(:),vmax(:),vmin(:),
     <          grdnrm(:),ra(:),snm(:),snp(:),vrtm(:),vrtp(:)
      real,allocatable :: sn(:,:),vort(:,:)
      real,allocatable :: qu(:,:),qv(:,:),qw(:,:),
     <          qup(:,:),qvp(:,:),qwp(:,:),
     <          qur(:,:),qvr(:,:),qwr(:,:)

      real c1,c2,c3,c4,c5
      DATA C1,C2,C3,C4,C5/0.4,26.0,0.01688,1.00,0.3/
      
      integer ka,kb,ledge,ledge1,j,k,l,kp,kr

      real re,dx2,dy2,dz2,wmu,ux,vx,wx,uy,vy,wy,uz,vz,wz
      real tx,ty,tz,rhowmuw,dx,dy,dz,fl,uu
      real flp,flm,dfm,dfp,dsm,dsp,am,bm,si,fli,dv2
      real t1,t2,fwake,t3,fkleb,slen,tmi,rhomuw

      allocate(fmax(mdim),smax(mdim),vmax(mdim),vmin(mdim),
     <          grdnrm(mdim),ra(mdim),snm(mdim),snp(mdim),vrtm(mdim),
     <          vrtp(mdim))
      allocate(sn(jmax,lmax),vort(jmax,lmax))
      allocate(qu(jmax,lmax),qv(jmax,lmax),qw(jmax,lmax),
     <          qup(jmax,lmax),qvp(jmax,lmax),qwp(jmax,lmax),
     <          qur(jmax,lmax),qvr(jmax,lmax),qwr(jmax,lmax))
c***********************************************************************
c      ka = 2 - ksym

      if (kroot.eq.1) then
         ka = 2 - ksym
      else
         ka=kroot
      endif

      kb = ktip
c
      re = rey
c
c..compute vorticity and normal distances from wall for entire flowfield
c
      ledge   = (3*lmax)/4
      ledge1  = ledge + 1
      dx2     = 0.5
      dy2     = 0.5
      dz2     = 0.5
c
c..zero out turmu
c
      do 12 l = 1,lmax
      do 12 k = 1,kmax
      do 12 j = 1,jmax
        turmu(j,k,l) = 0.0
   12 continue

      if (.not. is_wing) return

c
c..set wall molecular viscosity
c
      wmu = 1.0
      do 10 k = ka,kb
        kp = kkp(k)
        kr = kkr(k)
c
cgrs..for compatibility with turbulence model
c
        do 120 l = 1,ledge1+1
        do 120 j = jtail1-1,jtail2+1
          qu(j,l)  = q(j,k,l,2)  - q(j,k,l,1)*ug(j,k,l)
          qup(j,l) = q(j,kp,l,2) - q(j,kp,l,1)*ug(j,kp,l)
          qur(j,l) = q(j,kr,l,2) - q(j,kr,l,1)*ug(j,kr,l)
          qv(j,l)  = q(j,k,l,3)  - q(j,k,l,1)*vg(j,k,l)
          qvp(j,l) = q(j,kp,l,3) - q(j,kp,l,1)*vg(j,kp,l)
          qvr(j,l) = q(j,kr,l,3) - q(j,kr,l,1)*vg(j,kr,l)
          qw(j,l)  = q(j,k,l,4)  - q(j,k,l,1)*wg(j,k,l)
          qwp(j,l) = q(j,kp,l,4) - q(j,kp,l,1)*wg(j,kp,l)
          qwr(j,l) = q(j,kr,l,4) - q(j,kr,l,1)*wg(j,kr,l)
 120    continue
c
        l   = 1
        do 1 j = jtail1,jtail2
c
          ux  = (qu(j+1,l)/q(j+1,k,l,1) -
     *           qu(j-1,l)/q(j-1,k,l,1)) * dx2
          vx  = (qv(j+1,l)/q(j+1,k,l,1) -
     *           qv(j-1,l)/q(j-1,k,l,1)) * dx2
          wx  = (qw(j+1,l)/q(j+1,k,l,1) -
     *           qw(j-1,l)/q(j-1,k,l,1)) * dx2
          uy  = (qup(j,l)/q(j,kp,l,1) -
     *           qur(j,l)/q(j,kr,l,1)) * dy2
          vy  = (qvp(j,l)/q(j,kp,l,1) -
     *           qvr(j,l)/q(j,kr,l,1)) * dy2
          wy  = (qwp(j,l)/q(j,kp,l,1) -
     *           qwr(j,l)/q(j,kr,l,1)) * dy2
          uz  = -(3.0*qu(j,l)/q(j,k,l,1) -
     &            4.0*qu(j,l+1)/q(j,k,l+1,1) +
     &                qu(j,l+2)/q(j,k,l+2,1))*dz2
          vz  = -(3.0*qv(j,l)/q(j,k,l,1) -
     &            4.0*qv(j,l+1)/q(j,k,l+1,1) +
     &                qv(j,l+2)/q(j,k,l+2,1))*dz2
          wz  = -(3.0*qw(j,l)/q(j,k,l,1) -
     &            4.0*qw(j,l+1)/q(j,k,l+1,1) +
     &                qw(j,l+2)/q(j,k,l+2,1))*dz2
c     
          tx  =  xy(j,k,l)*ux -xx(j,k,l)*vx +yy(j,k,l)*uy
     *          -yx(j,k,l)*vy +zy(j,k,l)*uz -zx(j,k,l)*vz
          ty  =  xz(j,k,l)*vx -xy(j,k,l)*wx +yz(j,k,l)*vy
     *          -yy(j,k,l)*wy +zz(j,k,l)*vz -zy(j,k,l)*wz
          tz  =  xx(j,k,l)*wx -xz(j,k,l)*ux +yx(j,k,l)*wy
     *          -yz(j,k,l)*uy +zx(j,k,l)*wz -zz(j,k,l)*uz
c
          vort(j,l) = sqrt(tx**2 +ty**2 +tz**2)
          sn(j,l)  = 0.0
    1   continue
        l   = 1
        do 2 j = jtail1,jtail2
          rhomuw  = q(j,k,l,1)*q(j,k,l,6)/wmu
          ra(j)   = sqrt( re*rhomuw*vort(j,l) )/c2
c
cgrs..for moving blades
cgrs..   uu = sqrt(((uwal**2) + vwal**2) + wwal**2)
cgrs..  *             /q(j,k,l,1)
c
          uu = sqrt(qu(j,l)**2 +qv(j,l)**2 + qw(j,l)**2)
     &                                            /q(j,k,l,1)
c
          fmax(j)  = 1.e-3
          vmax(j)  = uu
          vmin(j)  = uu
          vrtp(j)  = 0.0
          vrtm(j)  = 0.0
          grdnrm(j) = sqrt(zx(j,k,l)**2 +zy(j,k,l)**2 +zz(j,k,l)**2)
    2   continue
c
        do 3 l = 2,ledge1
        do 3 j = jtail1,jtail2
          ux  = (qu(j+1,l)/q(j+1,k,l,1) -
     &           qu(j-1,l)/q(j-1,k,l,1)) * dx2
          vx  = (qv(j+1,l)/q(j+1,k,l,1) -
     &           qv(j-1,l)/q(j-1,k,l,1)) * dx2
          wx  = (qw(j+1,l)/q(j+1,k,l,1) -
     &           qw(j-1,l)/q(j-1,k,l,1)) * dx2
          uy  = (qup(j,l)/q(j,kp,l,1) -
     &           qur(j,l)/q(j,kr,l,1)) * dy2
          vy  = (qvp(j,l)/q(j,kp,l,1) -
     &           qvr(j,l)/q(j,kr,l,1)) * dy2
          wy  = (qwp(j,l)/q(j,kp,l,1) -
     &           qwr(j,l)/q(j,kr,l,1)) * dy2
          uz  = (qu(j,l+1)/q(j,k,l+1,1) -
     &           qu(j,l-1)/q(j,k,l-1,1)) * dz2
          vz  = (qv(j,l+1)/q(j,k,l+1,1) -
     &           qv(j,l-1)/q(j,k,l-1,1)) * dz2
          wz  = (qw(j,l+1)/q(j,k,l+1,1) -
     &           qw(j,l-1)/q(j,k,l-1,1)) * dz2
c
          tx  =  xy(j,k,l)*ux -xx(j,k,l)*vx +yy(j,k,l)*uy
     *          -yx(j,k,l)*vy +zy(j,k,l)*uz -zx(j,k,l)*vz
          ty  =  xz(j,k,l)*vx -xy(j,k,l)*wx +yz(j,k,l)*vy
     *          -yy(j,k,l)*wy +zz(j,k,l)*vz -zy(j,k,l)*wz
          tz  =  xx(j,k,l)*wx -xz(j,k,l)*ux +yx(j,k,l)*wy
     *          -yz(j,k,l)*uy +zx(j,k,l)*wz -zz(j,k,l)*uz
c
          vort(j,l) = sqrt(tx**2 +ty**2 +tz**2)
c
          dx  = x(j,k,l) -x(j,k,1)
          dy  = y(j,k,l) -y(j,k,1)
          dz  = z(j,k,l) -z(j,k,1)
c
          sn(j,l)  = (zx(j,k,1)*dx +zy(j,k,1)*dy +zz(j,k,1)*dz)/
     *                   grdnrm(j)
    3   continue
        l   = 2
        do 4 j = jtail1,jtail2
          snp(j)   = sn(j,l+1)
          smax(j)  = sn(j,l)
          snm(j)   = sn(j,l-1)
    4   continue
c
c..compute fmax, smax, vmax, vmin
c
        do 11 l = 2,ledge
        do 11 j = jtail1,jtail2
c
          fl       = sn(j,l)*vort(j,l)*(1.0 -exp(-ra(j)*sn(j,l)))
          fmax(j)  = amax1(fmax(j),fl)
ccray          smax(j)  = cvmgt(sn(j,l),smax(j),fmax(j).eq.fl)
ccray          snp(j)   = cvmgt(sn(j,l+1),snp(j),fmax(j).eq.fl)
ccray          snm(j)   = cvmgt(sn(j,l-1),snm(j),fmax(j).eq.fl)
ccray          vrtp(j)  = cvmgt(vort(j,l+1),vrtp(j),fmax(j).eq.fl)
ccray          vrtm(j)  = cvmgt(vort(j,l-1),vrtm(j),fmax(j).eq.fl)
          if(fmax(j).eq.fl) smax(j) = sn(j,l)
          if(fmax(j).eq.fl) snp(j) = sn(j,l+1)
          if(fmax(j).eq.fl) snm(j) = sn(j,l-1)
          if(fmax(j).eq.fl) vrtp(j) = vort(j,l+1)
          if(fmax(j).eq.fl) vrtm(j) = vort(j,l-1)
c
          uu = sqrt(qu(j,l)**2 +qv(j,l)**2 + qw(j,l)**2)
     <                                              /q(j,k,l,1)
c
cgrs..for moving blades
cgrs..    uu = sqrt(((uwal**2) + vwal**2) + wwal**2)
cgrs..   *             /q(j,k,l,1)
c
          vmax(j)  = amax1(uu,vmax(j))
          vmin(j)  = amin1(uu,vmin(j))
   11   continue
c
c..interpolation to improve estimate of fmax and smax
c
        do 21 j = jtail1,jtail2
          flp = snp(j)*vrtp(j)*(1.0 -exp(-ra(j)*snp(j)))
          flm = snm(j)*vrtm(j)*(1.0 -exp(-ra(j)*snm(j)))
          flp = amax1(flp,1.e-3)
          flm = amax1(flm,1.e-3)
          dfm = fmax(j) -flm
          dfp = flp -fmax(j)
          dsm = smax(j) -snm(j)
          dsp = snp(j) -smax(j)
          am  = (dsp**2*dfm +dsm**2*dfp)/(dsp*dsm*(dsp+dsm))
          bm  = (dsm*dfp -dsp*dfm)/(dsp*dsm*(dsp+dsm))
ccray          si  = cvmgt(smax(j) -0.5*am/(bm+1.e-21),smax(j),bm.lt.0.0 )
          si = smax(j)
          if(bm.lt.0.0) si  = smax(j) -0.5*am/(bm+1.e-21)
ccray          fli = cvmgt(fmax(j) -0.25*am**2/(bm+1.e-21),fmax(j),bm.lt.0.0)
          fli = fmax(j)
          if(bm.lt.0.0) fli  = fmax(j) -0.25*am**2/(bm+1.e-21)
ccray          fmax(j)  = cvmgt(fmax(j),fli, fli.lt.fmax(j)
ccray     *                       .or. si.lt.snm(j) .or. si.gt.snp(j))
          if(.not.(fli.lt.fmax(j) .or. si.lt.snm(j) .or. si.gt.snp(j))) 
     *                                                   fmax(j) = fli
ccray          smax(j)  = cvmgt(smax(j),si,fli.lt.fmax(j)
ccray     *                      .or. si.lt.snm(j) .or. si.gt.snp(j))
          if(.not.(fli.lt.fmax(j) .or. si.lt.snm(j) .or. si.gt.snp(j)))
     *                                                   smax(j) = si
   21   continue
c
c..compute outer eddy viscosity
c
        do 31 l = 2,ledge
        do 31 j = jtail1,jtail2
          dv2   = (vmax(j) -vmin(j))**2
          t1    = smax(j)*fmax(j)
          t2    = c4*smax(j)*dv2/fmax(j)
          fwake = amin1(t1,t2)
          t3    = (c5*sn(j,l)/smax(j))
          fkleb = 1.0 +5.5*amin1(1.e5,t3)**6
          turmu(j,k,l) = c3*q(j,k,l,1)*q(j,k,l,6)*fwake/fkleb
   31   continue
c
c..compute inner eddy viscosity and set final eddy viscosity
c
        do 41 l = 2,ledge
        do 41 j = jtail1,jtail2
          slen  = c1*sn(j,l)*(1.0 -exp(-ra(j)*sn(j,l)))
          tmi   = q(j,k,l,1)*q(j,k,l,6)*slen**2*vort(j,l)
          turmu(j,k,l) = amin1(turmu(j,k,l),tmi)*re
   41   continue
c
   10 continue
c
      return
      end




      



c*************************************************************************
      subroutine confine(q,s,xx,xy,xz,yx,yy,yz,zx,zy,zz)
c
c  add vorticity confinement term to rhs
c  see steinhoff's original paper for details
c
c*************************************************************************

      use params_global
      implicit none

      
      real q(jmax,kmax,lmax,nd), s(jmax,kmax,lmax,nv)
      
      real xx(jmax,kmax,lmax), xy(jmax,kmax,lmax) ,xz(jmax,kmax,lmax),
     $     yx(jmax,kmax,lmax), yy(jmax,kmax,lmax) ,yz(jmax,kmax,lmax),
     $     zx(jmax,kmax,lmax), zy(jmax,kmax,lmax) ,zz(jmax,kmax,lmax)
      

c..   local variables

      real,allocatable :: vort(:,:,:,:)
      integer l,lp1,lm1,k,kp1,km1,j,jp1,jm1
      
      real epsilon,rm1,fac,dux,dvx,dwx
      real duy,dvy,dwy,rp1,duz,dvz,dwz
      real dwdy,dvdz,dudz,dwdx,dvdx,dudy
      real dvorx,dvory,dvorz,snx,sny,snz,smag
      real tx,ty,tz,eterm


      allocate(vort(jmax,kmax,lmax,4))
      
      

c***********************************

c..   bridgeman et.al ahs tech spec meeting, jan 2002

      epsilon=0.005

c..   first compute vorticity

      do l=1,lmax

         lp1=min(l+1,lmax)
         lm1=max(l-1,1)

         do k=1,kmax

            kp1=min(k+1,kmax)
            km1=max(k-1,1)

            do j=1,jmax
               
               jp1=min(j+1,jmax)
               jm1=max(j-1,1)
               
c..   2nd order central in computational space

               rp1=1.0/q(jp1,k,l,1)
               rm1=1.0/q(jm1,k,l,1)

               fac=1.0/(jp1-jm1)

               dux=(q(jp1,k,l,2)*rp1-q(jm1,k,l,2)*rm1)*fac
               dvx=(q(jp1,k,l,3)*rp1-q(jm1,k,l,3)*rm1)*fac
               dwx=(q(jp1,k,l,4)*rp1-q(jm1,k,l,4)*rm1)*fac

               rp1=1.0/q(j,kp1,l,1)
               rm1=1.0/q(j,km1,l,1)
               fac=1.0/(kp1-km1)

               duy=(q(j,kp1,l,2)*rp1-q(j,km1,l,2)*rm1)*fac
               dvy=(q(j,kp1,l,3)*rp1-q(j,km1,l,3)*rm1)*fac
               dwy=(q(j,kp1,l,4)*rp1-q(j,km1,l,4)*rm1)*fac

               rp1=1.0/q(j,k,lp1,1)
               rm1=1.0/q(j,k,lm1,1)
               fac=1.0/(lp1-lm1)

               duz=(q(j,k,lp1,2)*rp1-q(j,k,lm1,2)*rm1)*fac
               dvz=(q(j,k,lp1,3)*rp1-q(j,k,lm1,3)*rm1)*fac
               dwz=(q(j,k,lp1,4)*rp1-q(j,k,lm1,4)*rm1)*fac
               
c..   transform to physical space

               dwdy=xy(j,k,l)*dwx+yy(j,k,l)*dwy+zy(j,k,l)*dwz
               dvdz=xz(j,k,l)*dvx+yz(j,k,l)*dvy+zz(j,k,l)*dvz

               dudz=xz(j,k,l)*dux+yz(j,k,l)*duy+zz(j,k,l)*duz
               dwdx=xx(j,k,l)*dwx+yx(j,k,l)*dwy+zx(j,k,l)*dwz

               dvdx=xx(j,k,l)*dvx+yx(j,k,l)*dvy+zx(j,k,l)*dvz
               dudy=xy(j,k,l)*dux+yy(j,k,l)*duy+zy(j,k,l)*duz
               
c..   vorticity and its magnitude

               vort(j,k,l,1)=dwdy-dvdz
               vort(j,k,l,2)=dudz-dwdx
               vort(j,k,l,3)=dvdx-dudy
               
               vort(j,k,l,4)=sqrt(vort(j,k,l,1)**2+vort(j,k,l,2)**2
     $              +vort(j,k,l,3)**2)
               
            enddo
         enddo
      enddo

          

               
c..   compute gradient vector of vorticity magnitude surface
c..   take cross product with vorticity vector and add it as source
c..   term

         do l=2,lmax-1
            
            lp1=l+1
            lm1=l-1
            
            do k=2,kmax-1
               
               kp1=k+1
               km1=k-1

               do j=2,jmax-1
                  
                  jp1=j+1
                  jm1=j-1

c..   gradient in computational space

                  dvorx=(vort(jp1,k,l,4)-vort(jm1,k,l,4))*0.5
                  dvory=(vort(j,kp1,l,4)-vort(j,km1,l,4))*0.5
                  dvorz=(vort(j,k,lp1,4)-vort(j,k,lm1,4))*0.5

c..   transform to physical space

                  snx=xx(j,k,l)*dvorx+yx(j,k,l)*dvory+zx(j,k,l)*dvorz
                  sny=xy(j,k,l)*dvorx+yy(j,k,l)*dvory+zy(j,k,l)*dvorz
                  snz=xz(j,k,l)*dvorx+yz(j,k,l)*dvory+zz(j,k,l)*dvorz

c..   normalize
                  smag=1.0/sqrt(snx**2+sny**2+snz**2)
                  
                  snx=snx*smag
                  sny=sny*smag
                  snz=snz*smag

c..   cross normal to center of vorticity and vorticity vector itself

                  tx=sny*vort(j,k,l,3)-snz*vort(j,k,l,2)
                  ty=snz*vort(j,k,l,1)-snx*vort(j,k,l,3)
                  tz=snx*vort(j,k,l,2)-sny*vort(j,k,l,1)

c..   add as source term in momentum equations (incompressible)
                
                  s(j,k,l,2)=s(j,k,l,2)-epsilon*tx*q(j,k,l,1)
                  s(j,k,l,3)=s(j,k,l,3)-epsilon*ty*q(j,k,l,1)
                  s(j,k,l,4)=s(j,k,l,4)-epsilon*tz*q(j,k,l,1)

c..   for compressible formulation (add dot product of body force and 
c..   velocity), body force is what we added in the momentum equations

                  eterm=-epsilon*(tx*q(j,k,l,2)+ty*q(j,k,l,3)
     &                 +tz*q(j,k,l,4))

                  s(j,k,l,5)=s(j,k,l,5)+eterm
                  
               enddo
            enddo
         enddo
                  
      return
      end


