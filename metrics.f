
c***********************************************************************
      subroutine mett( x,y,z,xx,xy,xz,yx,yy,yz,zx,zy,zz,
     &    ug,vg,wg)
c
c  time metrics here for unsteady calculations
c     
c  indicial.eq.1 -> step change in angle of attack
c  indicial.eq.2 -> step change in pitch rate
c
c  note: since the metrics are scaled by the jacobian, time metrics 
c        calculated below have in them this scaling.  do not scale
c        again by jacobian!
c
c     rf = reduced frequency
c
c***********************************************************************

      use params_global

      implicit none

      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax),vg(jmax,kmax,lmax),wg(jmax,kmax,lmax)
      
      real xx(jmax,kmax,lmax),xy(jmax,kmax,lmax),xz(jmax,kmax,lmax),
     <     yx(jmax,kmax,lmax),yy(jmax,kmax,lmax),yz(jmax,kmax,lmax),
     <     zx(jmax,kmax,lmax),zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)

c..   local variables
      
      real xpsi,fmm,decay,xr
      integer j,k,l
      
c..   first executable statement

      if (indicial.eq.1) then 

c..   step change in angle of attack

         xpsi=totime*rf
         do  l = 1,lmax
            do  k = 1,kmax
               fmm=fmtip*sqrt(x(jle,k,1)**2+y(jle,k,1)**2)/rartio+
     <              fsmach*sin(xpsi)
               do  j = 1,jmax
                  ug(j,k,l) = -rf * y(j,k,l)
                  vg(j,k,l) = rf * x(j,k,l)
                  wg(j,k,l) = -angmax*fmm
               enddo
            enddo
         enddo

      elseif (indicial.eq.2) then

c...  step change in pitch rate

         do l = 1,lmax
            do k = 1,kmax
               do j = 1,jmax
                  xr=sqrt((x(j,k,l)-0.25)**2+z(j,k,l)**2)
                  ug(j,k,l) = -rf * y(j,k,l)
     <                 + alpha_dot*decay(xr)*z(j,k,l)*fsmach
                  vg(j,k,l) = rf * x(j,k,l)
                  wg(j,k,l) = -alpha_dot*(x(j,k,l)-0.25)
     <                 *decay(xr)* fsmach
               enddo
            enddo
         enddo  
         
      else

         do  l = 1,lmax
            do  k = 1,kmax
               do j = 1,jmax
                  ug(j,k,l) = -rf * y(j,k,l)
                  vg(j,k,l) = rf * x(j,k,l)
                  wg(j,k,l) = 0.0 + xlam*fmtip
               enddo
            enddo
         enddo

      endif
         
      return
      end

c***********************************************************************
      subroutine rotate_grid(srot,x,y,z,xt2,yt2,zt2)

c***********************************************************************
      use params_global

      implicit none

      real srot
      real x(jmax,kmax,lmax), y(jmax,kmax,lmax), z(jmax,kmax,lmax)
      real xt2(jmax,kmax,lmax), yt2(jmax,kmax,lmax), zt2(jmax,kmax,lmax)

      !local
      integer j, k, l
      real xtmp, ytmp, ztmp
      real cs, ss

      ss = sin(srot)
      cs = cos(srot)

      do 81 l = 1,lmax
      do 81 k = 1,kmax
      do 81 j = 1,jmax

c..store grid config

        xt2(j,k,l) = x(j,k,l)
        yt2(j,k,l) = y(j,k,l)
        zt2(j,k,l) = z(j,k,l)

c..rotate to new grid config
        xtmp = x(j,k,l) * cs - y(j,k,l) * ss
        ytmp = x(j,k,l) * ss + y(j,k,l) * cs
        ztmp = z(j,k,l)

        x(j,k,l)  = xtmp
        y(j,k,l)  = ytmp
        z(j,k,l)  = ztmp

 81   continue

      return
      end

c***********************************************************************
      subroutine rotate(srot,x,y,z,xx,xy,xz,ug,yx,yy,yz,vg,zx,zy,zz,wg,
     &                  zx0,zy0,zz0,zt0)
c
c  this stores the old surface metrics as well as rotates the grid and 
c  updates the metrics.  
c  note: time derivates are an invariant of the transformation.
c
c***********************************************************************

      use params_global
      
      implicit none
      
      real xx(jmax,kmax,lmax),xy(jmax,kmax,lmax),xz(jmax,kmax,lmax),
     $     yx(jmax,kmax,lmax),yy(jmax,kmax,lmax),yz(jmax,kmax,lmax),
     $     zx(jmax,kmax,lmax),zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax), vg(jmax,kmax,lmax),wg(jmax,kmax,lmax)
      real x(jmax,kmax,lmax), y(jmax,kmax,lmax), z(jmax,kmax,lmax )
      real zx0(jmax,kmax),zy0(jmax,kmax),zz0(jmax,kmax),zt0(jmax,kmax)

c..   local variables

      integer j,k,l
      real cs,ss,srot
      real xtmp,ytmp,ztmp,xxt,xyt,yxt,yyt,zxt,zyt

c..   first executable statement

      ss = sin(srot)
      cs = cos(srot)

      do 450 j = 1,jmax
      do 450 k = 1,kmax
        zx0(j,k) = zx(j,k,1)
        zy0(j,k) = zy(j,k,1)
        zz0(j,k) = zz(j,k,1)
        zt0(j,k) = -ug(j,k,1)*zx(j,k,1)-vg(j,k,1)*zy(j,k,1)
     <             -wg(j,k,1)*zz(j,k,1)
  450 continue

      do 81 l = 1,lmax
      do 81 k = 1,kmax
      do 81 j = 1,jmax
c..grid
        xtmp= x(j,k,l) * cs - y(j,k,l) * ss
        ytmp= x(j,k,l) * ss + y(j,k,l) * cs
        ztmp= z(j,k,l)
        x(j,k,l)  = xtmp
        y(j,k,l)  = ytmp
        z(j,k,l)  = ztmp
c..xi metrics       
        xxt = xx(j,k,l) * cs - xy(j,k,l) * ss
        xyt = xx(j,k,l) * ss + xy(j,k,l) * cs
        xx(j,k,l) = xxt
        xy(j,k,l) = xyt
        xz(j,k,l) = xz(j,k,l)
c..eta metrics   
        yxt = yx(j,k,l) * cs - yy(j,k,l) * ss
        yyt = yx(j,k,l) * ss + yy(j,k,l) * cs
        yx(j,k,l) = yxt
        yy(j,k,l) = yyt
        yz(j,k,l) = yz(j,k,l)
c..zeta metrics  
        zxt = zx(j,k,l) * cs - zy(j,k,l) * ss
        zyt = zx(j,k,l) * ss + zy(j,k,l) * cs
        zx(j,k,l) = zxt
        zy(j,k,l) = zyt
        zz(j,k,l) = zz(j,k,l)
 81   continue

      return
      end



      function decay(r)

      real r
      real decay,angle,pi
      
      pi=acos(-1.)
      if (r.le.2.) then
         decay=1.
      elseif (r.gt.2.0.and.r.le.5.0) then
         angle=(r-2.0)*pi/3.0
         decay=0.5+0.5*cos(angle)
      elseif (r.gt.5.0) then
         decay=0.0
      endif
      
      return
      end


C***********************************************************************
      subroutine metfv(q,x,y,z,xx,xy,xz,yx,yy,yz,zx,zy,zz,a3,
     &            kkp,kkr,iprint)


      use params_global
      implicit none

C***  Prologue : ***********
C
C  Finite volume formulation         11/4/87  S.O.
C  Compute the metrics and the Jacobian for computational space
C  uniform computational space, DELTAS = 1 < averaged metrics >
C
C*** End Prologue: ****************************************************

      integer iprint
      
      real q(jmax,kmax,lmax,nd),x(jmax,kmax,lmax),y(jmax,kmax,lmax)
      real z(jmax,kmax,lmax),a3(jmax,kmax,lmax)

      real xx(jmax,kmax,lmax),xy(jmax,kmax,lmax),xz(jmax,kmax,lmax),
     <     yx(jmax,kmax,lmax),yy(jmax,kmax,lmax),yz(jmax,kmax,lmax),
     <     zx(jmax,kmax,lmax),zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)

      integer kkp(kmax),kkr(kmax)

      integer j,k,l,jp1,kp1,lp1
      integer ll,lm1,kk,km1,jj,jm1,nneg
      integer jjp(jmax),jjr(jmax),llp(lmax),llr(lmax)
      real dx11,dy11,dz11,dx22,dy22,dz22,sx1,sy1,sz1,qj
      real fac1,fac2

C********************************************************************

      jm      = jmax -1
      km      = kmax -1
      lm      = lmax -1
      fac1    = 0.5
      fac2    = 1./3.
c
      do 1  j = 1,jmax
       jjp(j) = j + 1
       jjr(j) = j - 1
    1 continue
      jjp(jmax)=jmax
      jjr(1   )=1
      do 2  l = 1,lmax
        llp(l) = l + 1
        llr(l) = l - 1
    2 continue
      llp(lmax)=lmax
      llr(1   )=1
c
c..xi derivatives
c
      do 10 l = 1,lm
      lp1 = llp(l)
      do 10 k = 1,km
      kp1 = kkp(k)
        do 11 j = 1,jmax
          dx11 = (x(j,kp1,lp1) -x(j,k,l))
          dy11 = (y(j,kp1,lp1) -y(j,k,l))
          dz11 = (z(j,kp1,lp1) -z(j,k,l))
          dx22 = (x(j,k,lp1) -x(j,kp1,l))
          dy22 = (y(j,k,lp1) -y(j,kp1,l))
          dz22 = (z(j,k,lp1) -z(j,kp1,l))
c
          xx(j,k,l) = fac1*( dy11*dz22 - dy22*dz11 )
          xy(j,k,l) = fac1*( dz11*dx22 - dz22*dx11 )
          xz(j,k,l) = fac1*( dx11*dy22 - dx22*dy11 )
   11   continue
   10 continue
c
c..eta derivatives
c
      do 20 l = 1,lm
      lp1 = llp(l)
      do 20 j = 1,jm
      jp1 = jjp(j)
        do 21 k = 1,kmax
          dx11 = (x(jp1,k,lp1) -x(j,k,l  ))
          dy11 = (y(jp1,k,lp1) -y(j,k,l  ))
          dz11 = (z(jp1,k,lp1) -z(j,k,l  ))
          dx22 = (x(jp1,k,l  ) -x(j,k,lp1))
          dy22 = (y(jp1,k,l  ) -y(j,k,lp1))
          dz22 = (z(jp1,k,l  ) -z(j,k,lp1))
c
          yx(j,k,l)  = fac1*( dy11*dz22 - dy22*dz11 )
          yy(j,k,l)  = fac1*( dz11*dx22 - dz22*dx11 )
          yz(j,k,l)  = fac1*( dx11*dy22 - dx22*dy11 )
   21   continue
   20 continue
c
c..zeta derivative
c
      do 30 k = 1,km
      kp1 = kkp(k)
      do 30 j = 1,jm
      jp1 = jjp(j)
        do 31 l = 1,lmax
          dx11 = (x(jp1,kp1,l) -x(j  ,k  ,l))
          dy11 = (y(jp1,kp1,l) -y(j  ,k  ,l))
          dz11 = (z(jp1,kp1,l) -z(j  ,k  ,l))
          dx22 = (x(j  ,kp1,l) -x(jp1,k  ,l))
          dy22 = (y(j  ,kp1,l) -y(jp1,k  ,l))
          dz22 = (z(j  ,kp1,l) -z(jp1,k  ,l))
c
          zx(j,k,l)  = fac1*( dy11*dz22 - dy22*dz11 )
          zy(j,k,l)  = fac1*( dz11*dx22 - dz22*dx11 )
          zz(j,k,l)  = fac1*( dx11*dy22 - dx22*dy11 )
   31   continue
   30 continue
      
c
c..compute cell volume
c
      do 40 l = 1,lm
      lp1 = llp(l)
      do 40 k = 1,km
      kp1 = kkp(k)
        do 41 j = 1,jm
          jp1 = jjp(j)
c
          sx1 = xx(j,k,l)+yx(j,k,l)+zx(j,k,l)
          sy1 = xy(j,k,l)+yy(j,k,l)+zy(j,k,l)
          sz1 = xz(j,k,l)+yz(j,k,l)+zz(j,k,l)
c
          dx11 = ( x(jp1,kp1,lp1)-x(j,k,l) )
          dy11 = ( y(jp1,kp1,lp1)-y(j,k,l) )
          dz11 = ( z(jp1,kp1,lp1)-z(j,k,l) )
c
          a3(j,k,l) = fac2*( sx1*dx11 + sy1*dy11 + sz1*dz11 )
   41 continue
   40 continue
c
c..finite-difference conventions
c
      do 42 ll = 1,lmax
      l  = min(ll,lm)
      lm1 = llr(ll)
      do 42 kk = 1,kmax
      k  = min(kk,km)
      km1 = kkr(kk)
        do 43 jj = 1,jmax
          j  = min(jj,jm)
          jm1 = jjr(jj)
          q(jj,kk,ll,6) = 8.0/( a3(j,k,l) +a3(jm1,km1,lm1)
     <                   +a3(jm1,k,l) +a3(j,km1,l) +a3(j,k,lm1)
     <                 +a3(jm1,km1,l) +a3(jm1,k,lm1) +a3(j,km1,lm1) )
   43   continue
   42 continue

c
      do 44 ll = lmax,1,-1
      l  = min(ll,lm)
      lm1 = llr(ll)
      do 44 kk = kmax,1,-1
      k  = min(kk,km)
      km1 = kkr(kk)
       do 45 j  = 1,jmax
       qj         = 0.25*q(j,kk,ll,6)
       xx(j,kk,ll)= qj*(xx(j,k,l)+xx(j,km1,l)+xx(j,k,lm1)+xx(j,km1,lm1))
       xy(j,kk,ll)= qj*(xy(j,k,l)+xy(j,km1,l)+xy(j,k,lm1)+xy(j,km1,lm1))
       xz(j,kk,ll)= qj*(xz(j,k,l)+xz(j,km1,l)+xz(j,k,lm1)+xz(j,km1,lm1))
   45  continue
   44 continue
c
      do 53 ll = lmax,1,-1
      l  = min(ll,lm)
      lm1 = llr(ll)
      do 53 jj = jmax,1,-1
      j  = min(jj,jm)
      jm1 = jjr(jj)
       do 54 k  = 1,kmax
       qj         = 0.25*q(jj,k,ll,6)
	if (kmax.eq.ktip .and.
     $      k.eq.kmax .and.
     $      ll.eq.1.and.jj.eq.((jmax+1)/2) .and.
     $      dom_type .ne. (DOM_WAKE.or.DOM_GROUNDWAKE)) then
c	dkarthik - physical corner point
	yx(jj,k,ll)=qj*4.*yx(j,k,l)
	yy(jj,k,ll)=qj*4.*yy(j,k,l)
	yz(jj,k,ll)=qj*4.*yz(j,k,l)
c      elseif (root_co.eq.1.and.k.eq.1.and.ll.eq.1.and.
c     &     jj.eq.((jmax+1)/2).and. 
c          dom_type .ne. (DOM_WAKE.or.DOM_GROUNDWAKE)) then
c        yx(jj,k,ll)=qj*4.*yx(j,k,l)
c        yy(jj,k,ll)=qj*4.*yy(j,k,l)
c        yz(jj,k,ll)=qj*4.*yz(j,k,l)
      else
       yx(jj,k,ll)= qj*(yx(j,k,l)+yx(jm1,k,l)+yx(j,k,lm1)+yx(jm1,k,lm1))
       yy(jj,k,ll)= qj*(yy(j,k,l)+yy(jm1,k,l)+yy(j,k,lm1)+yy(jm1,k,lm1))
       yz(jj,k,ll)= qj*(yz(j,k,l)+yz(jm1,k,l)+yz(j,k,lm1)+yz(jm1,k,lm1))
      endif

   54  continue
   53 continue
c
      do 63 kk = kmax,1,-1
      k  = min(kk,km)
      km1 = kkr(kk)
      do 63 jj = jmax,1,-1
      j  = min(jj,jm)
      jm1 = jjr(jj)
       do 64 l  = 1,lmax
       qj         = 0.25*q(jj,kk,l,6)
       zx(jj,kk,l)= qj*(zx(j,k,l)+zx(j,km1,l)+zx(jm1,k,l)+zx(jm1,km1,l))
       zy(jj,kk,l)= qj*(zy(j,k,l)+zy(j,km1,l)+zy(jm1,k,l)+zy(jm1,km1,l))
       zz(jj,kk,l)= qj*(zz(j,k,l)+zz(j,km1,l)+zz(jm1,k,l)+zz(jm1,km1,l))
   64  continue
   63  continue
c
c..check for negative jacobians
c
      nneg = 0
      do 900 l = 1, lmax
        do 910 k = 1, kmax
          do 920 j = 1, jmax
ccray            nneg = cvmgt(nneg+1,nneg,q(j,k,l,6).le.0.)
            if(q(j,k,l,6).le.0.) nneg = nneg + 1
 920      continue
 910    continue
 900  continue
c
      if(nneg .ne. 0) then
        write(6,*) nneg, ' negative jacobians in block'
        do 74 l = 1,lmax
        do 74 k = 1,kmax
        do 74 j = 1,jmax
          if( q(j,k,l,6).lt.0.0 ) then
            write(6,603) q(j,k,l,6), j, k, l
          end if
   74   continue
      endif
c
  603 format( ' ',10x,'negative jacobian = ',1p,e10.3,1x,'at j,k,l =',
     $                 3i5,5x)
c
      return
      end

