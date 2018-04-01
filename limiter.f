C************************************************************************
      subroutine 
     <  muscld_new(f,ql,qr,is,imax,im,th,qt,eps,fmin,fmax,ibmin,ibmax,
     <     mdim,iba)
C
C  Differentiable limiter for 3rd order accuracy
C  Modified by dkarthik(2003) for higher order at boundaries
C  Also added iblanking
C*************************************************************************

      implicit none
      
      integer is,imax,im,ibmin,ibmax,mdim
      
      real th,qt,eps
      real f(mdim,5),ql(mdim,5),qr(mdim,5)
      real fmin(5),fmax(5)
      real f2(mdim,5)
      integer iba(mdim)

c..   local variables

      integer i,n
      real thm,thp,f2i,f2i1,a1,a2,epsf,f3i,f3qt

c****************

      IF(QT.EQ.0.0)THEN
        DO 20 N=1,5
        DO 20 I=IS,IMAX
          QL(I,N)=F(I,N)
          QR(I,N)=F(I,N)
   20   CONTINUE
        RETURN
      ELSE
        THM    = 1.0 - TH
        THP    = 1.0 + TH
C	no point in blanking boundaries
C	iba(IS)=1
C	iba(IMAX)=1
        DO 1 N=1,5
          DO 11 I=IS,IM
            F2(I,N) = F(I+1,N) - F(I,N)
 11       CONTINUE
C
          DO 12 I=IS+1,IM
           F2I    = F2(I,N)
           F2I1   = F2(I-1,N)
           A1     = 3.0*(F2I*F2I1+EPS)
           A2     = 2.0*(F2I-F2I1)**2 + A1
c Outlawed!!    EPSF   = EPS*( 0.5+SIGN( 0.5,-ABS(A2) ) )
           F3I    =  A1/A2 
           F3QT   = QT*F3I
     >			*IBA(I)*IBA(I+1)*IBA(I-1)
           QL(I,N)= F(I,N)+F3QT*( THM*F2I1 + THP*F2I )
           QR(I,N)= F(I,N)-F3QT*( THP*F2I1 + THM*F2I )
   12     CONTINUE
C
C..BOUNDARIES AT PRESENT ONLY 3rd Order
	   QL(IS,N)=F(IS,N)
	   QR(IS,N)=F(IS,N)
	   QL(IMAX,N)=F(IMAX,N)
	   QR(IMAX,N)=F(IMAX,N)
	   if(ibmin.eq.-2) then
	   F2I =F(IS+1,N)-F(IS,N)
	   F2I1=F(IS+2,N)-F(IS+1,N)
	   QL(IS,N)=F(IS,N)+0.5*(F2I*4./3.-F2I1*1./3.)
     >			      *IBA(IS)*IBA(IS+1)*IBA(IS+2)
	   QR(IS,N)=F(IS,N)	  
	   endif

	   if(ibmax.eq.-2) then
	   F2I  =F(IMAX,N)-F(IMAX-1,N)
	   F2I1 =F(IMAX-1,N)-F(IMAX-2,N)
	   QR(IMAX,N)=F(IMAX,N)-0.5*(F2I*4./3.-F2I1*1./3.)
     >			      *IBA(IMAX)*IBA(IMAX-1)*IBA(IMAX-2)
	   QL(IMAX,N)=F(IMAX,N)	  
	   endif

    1   CONTINUE

c...Check for -ve pressure / density
	  if(ql(is,1).le.0..or.ql(is,5).le.0.) then
	  print*,'Boundary extrapolation: p or rho going negative
     <		  hence dropping to first order'
	  QL(IS,1)=F(IS,1)
	  QL(IS,2)=F(IS,2)
	  QL(IS,3)=F(IS,3)
	  QL(IS,4)=F(IS,4)
	  QL(IS,5)=F(IS,5)
	  endif	  

	  if(qr(imax,1).le.0..or.qr(imax,5).le.0.) then
	  print*,'Boundary extrapolation: p or rho going negative
     <		  hence dropping to first order'
	  QR(Imax,1)=F(Imax,1)
	  QR(Imax,2)=F(Imax,2)
	  QR(Imax,3)=F(Imax,3)
	  QR(Imax,4)=F(Imax,4)
	  QR(Imax,5)=F(Imax,5)
	  endif	  
      ENDIF
C
      RETURN
      END

      function ammd(fl,fr)
      real fl,fr,ammd

      ammd=0.5*(sign(1.,fl)+sign(1.,fr))*amin1(abs(fl),abs(fr))

      return
      end

c*************************************************************************
      subroutine sonica(f,ql,qr,is,ie,im,th,qt,eps,fmin,fmax,
     <                  ibmin,ibmax,mdim)
c
c sonic-a scheme of hyunh et al
c note: qr(is) and ql(ie) are never used.
c*************************************************************************

      implicit none
      
      integer is,ie,imax,im,ibmin,ibmax,mdim
      
      real th,qt,eps
      real f(mdim,5),ql(mdim,5),qr(mdim,5)
      real fmin(5),fmax(5)
      real f0(mdim),f1(mdim),f2m(mdim),f2(mdim)

c..   local variables

      integer i,n
      real thm,thp,f2i,f2i1,a1,a2,epsf,f3i,f3qt
      real f0bmin,f1bmin,f2bmin,f0bmax,f1bmax,f2bmax
      real at,at1,s1,t1,slope,ati1
      real ammd

c..this is just 1st order upwind

      if(qt.eq.0.0)then
        do n=1,5
          do i=is,ie
            ql(i,n)=f(i,n)
            qr(i,n)=f(i,n)
          enddo
        enddo
        return
      else
c
c..the m3-a interpolation scheme of hyunh follows
c
        do 1 n=1,5
c
c..let's load up a few difference arrays
c
          do i=is,ie
            f0(i) = f(i,n)
          enddo
c..1st difference at i+1/2
          do i=is,ie-1
            f1(i) = f0(i+1) - f0(i)
          enddo
c..2nd difference at i
          do i=is+1,ie-1
            f2(i)  = f1(i) - f1(i-1)
          enddo
c..extrapolate at the boundaries
          f2(is) = 2.*f2(is+1)-f2(is+2)
          f2(ie) = 2.*f2(ie-1)-f2(ie-2)
c..modify at boundaries, if needed
          if(ibmin.eq.2) then
            f0bmin = fmin(n)
            f1bmin = f0(is) - f0bmin
            f2(is) = f1(is) - f1bmin
            f2bmin = 2.*f2(is)-f2(is+1)
          endif
          if(ibmax.eq.2) then
            f0bmax = fmax(n)
            f1bmax = f0bmax - f0(ie)
            f2(ie) = f1bmax - f1(ie)
            f2bmax = 2.*f2(ie)-f2(ie-1)
          endif
c..limit 2nd difference to i+1/2
          do i = is,ie-1
            f2m(i) = ammd(f2(i),f2(i+1))
          enddo
c
c..now combine everything to get ql and qr
c
          do i=is+1,ie-1
c..include limited curvatures to calculate new slopes
            at    = f1(i)   - 0.5*f2m(i)
            ati1  = f1(i-1) + 0.5*f2m(i-1)
c..limit slopes at i
            s1    = ammd(f1(i),f1(i-1))
            t1    = ammd(at,ati1)
c..now find appropriate slope
            slope  =  sign(1.,t1)*
     <         amin1(0.5*abs(at+ati1),amax1(2.*abs(s1),abs(t1)))
c..use slope to calculate ql and qr
            ql(i,n) = f0(i) + 0.5*slope
            qr(i,n) = f0(i) - 0.5*slope
          enddo
c
c..first-order at boundary?
c
          slope = f1(is) - 0.5*f2m(is)
          slope = ammd(slope,2.*f1(is))
          ql(is,n) = f0(is) + 0.5*slope
c
          if(ibmin.eq.2) then
            i = is
            s1    = ammd(f1(i),f1bmin)
            at    = f1(i)  - 0.5*f2m(i)
            ati1  = f1bmin + 0.5*ammd(f2(i),f2bmin)
            t1    = ammd(at,ati1)
            slope  =  sign(1.,t1)*
     <         amin1(0.5*abs(at+ati1),amax1(2.*abs(s1),abs(t1)))
            ql(i,n) = f0(i) + 0.5*slope
          endif
c
          slope = f1(ie-1) + 0.5*f2m(ie-1)
          slope = ammd(slope,2.*f1(ie-1))
          qr(ie,n) = f0(ie) - 0.5*slope
c
          if(ibmax.eq.2) then
            i = ie
            s1    = ammd(f1bmax,f1(i-1))
            at    = f1bmax  - 0.5*ammd(f2(i-1),f2bmax)
            ati1  = f1(i-1) + 0.5*f2m(i-1)
            t1    = ammd(at,ati1)
            slope  =  sign(1.,t1)*amax1(2.*abs(s1),abs(t1))
            slope  =  ammd(0.5*(at+ati1),slope)
            qr(i,n) = f0(i) - 0.5*slope
          endif
    1   continue
      endif
c
      return
      end

      function amedian(b,c,d)

      real b,c,d,amedian

      amedian=b+0.5*(sign(1.,c-b)+sign(1.,d-b))*
     &     amin1(abs(c-b),abs(d-b))
      
      return
      end

c**************************************************************************
      subroutine quad(f,ql,qr,is,ie,im,th,qt,eps,fmin,fmax,
     <                  fmin1,fmax1,fmin2,fmax2,ibmin,ibmax,mdim)
c
c..piecewise quadratic reconstruction with 6th order compact evaluation of
c  nodal derivatives and new monotonicity-preserving constraint.
c
c**************************************************************************

      implicit none

      integer is,ie,im,ibmin,ibmax,mdim

      real th,qt,eps
      real f(mdim,5),ql(mdim,5),qr(mdim,5)
      real fmin(5),fmax(5),fmin1(5),fmax1(5),fmin2(5),fmax2(5)

c..   local variables

      integer i,is1,is2,ie2,ie1,istart,iend,n

      real f1i,f1i1,f1i2,f1i3,f1i4,f2i1,f2i2,f2i3,f2i,f1bmin,f2bmin
      real f2mi,f2mi1,f2mi2,at,sl,ati1,t1,slopee2,sm,cm,si,ci,sr,cr
      real f2mi3,at1,f1bmax,f2bmax,cl
      real slopes2

      real d1(mdim+2),d2(mdim),con1(mdim+2),con2(mdim+2),con3(mdim+2)
      real s1(mdim+2),sf(mdim+2),sc(mdim+2),sb(mdim+2)
      real amedian,ammd


c..   first executable statement

      if(qt.eq.0.)then
        do 10 n=1,5
        do 10 i=is,ie
          ql(i,n)=f(i,n)
          qr(i,n)=f(i,n)
   10   continue
        return
      else
        is1=is+1
        is2=is+2
        ie2=ie-2
        ie1=ie-1
c
        do 1 n=1,5
c
c..slope:
c
        if(ibmin.ne.3) then
          i=is
          f1i=f(i+1,n)-f(i,n)
          f1i1=f(i+2,n)-f(i+1,n)
          f1i2=f(i+3,n)-f(i+2,n)
          f1i3=f(i+4,n)-f(i+3,n)
          f2i1=f1i1-f1i
          f2i2=f1i2-f1i1
          f2i3=f1i3-f1i2
          f2i=2.*f2i1-f2i2
          if(ibmin.eq.2) then
            f1bmin=f(i,n)-fmin(n)
            f2i=f1i-f1bmin
            f2bmin=2.*f2i-f2i1
          endif
          f2mi=ammd(f2i,f2i1)
          f2mi1=ammd(f2i1,f2i2)
          f2mi2=ammd(f2i2,f2i3)
          at=f1i-0.5*f2mi
          d1(i)=ammd(at,2.*f1i)
          if(ibmin.eq.2) then
            sl=ammd(f1i,f1bmin)
            at=f1i-0.5*f2mi
            ati1=f1bmin+0.5*ammd(f2i,f2bmin)
            t1=ammd(at,ati1)
            d1(i)=sign(1.,t1)*amin1(0.5*abs(at+ati1),
     &            amax1(2.*abs(sl),abs(t1)))
          endif
c
c..sonica scheme at the boundary:
c
          ql(i,n)=f(i,n)+0.5*d1(i)
          qr(i,n)=f(i,n)-0.5*d1(i)
c
          i=is1
          at=f1i1-0.5*f2mi1
          ati1=f1i+0.5*f2mi
          sl=ammd(f1i1,f1i)
          t1=ammd(at,ati1)
          d1(i)=sign(1.,t1)*
     &          amin1(0.5*abs(at+ati1),amax1(2.*abs(sl),abs(t1)))
c
          ql(i,n)=f(i,n)+0.5*d1(i)
          qr(i,n)=f(i,n)-0.5*d1(i)
c
          i=is2
          at=f1i2-0.5*f2mi2
          ati1=f1i1+0.5*f2mi1
          sl=ammd(f1i2,f1i1)
          t1=ammd(at,ati1)
          slopes2=sign(1.,t1)*
     &            amin1(0.5*abs(at+ati1),amax1(2.*abs(sl),abs(t1)))
c
          istart=is2
        else
          i=is-1
          con1(i)=1./3.
          con2(i)=1.
          con3(i)=(f(i+2,n)+28.*f(i+1,n)-28.*fmin1(n)-fmin2(n))/36.
          i=is
          con1(i)=1./3.
          con2(i)=1.
          con3(i)=(f(i+2,n)+28.*f(i+1,n)-28.*fmin(n)-fmin1(n))/36.
          i=is1
          con1(i)=1./3.
          con2(i)=1.
          con3(i)=(f(i+2,n)+28.*f(i+1,n)-28.*f(i-1,n)-fmin(n))/36.
          istart=is-1
        endif
c
        do i=is2,ie2
          con1(i)=1./3.
          con2(i)=1.
          con3(i)=(f(i+2,n)+28.*f(i+1,n)-28.*f(i-1,n)-f(i-2,n))/36.
        enddo
c
        if(ibmax.eq.3) then
          i=ie1
          con1(i)=1./3.
          con2(i)=1.
          con3(i)=(fmax(n)+28.*f(i+1,n)-28.*f(i-1,n)-f(i-2,n))/36.
          i=ie
          con1(i)=1./3.
          con2(i)=1.
          con3(i)=(fmax1(n)+28.*fmax(n)-28.*f(i-1,n)-f(i-2,n))/36.
          i=ie+1
          con1(i)=1./3.
          con2(i)=1.
          con3(i)=(fmax2(n)+28.*fmax1(n)-28.*f(i-1,n)-f(i-2,n))/36.
          iend=ie+1
        else
          i=ie
          f1i1=f(i,n)-f(i-1,n)
          f1i2=f(i-1,n)-f(i-2,n)
          f1i3=f(i-2,n)-f(i-3,n)
          f1i4=f(i-3,n)-f(i-4,n)
          f2i1=f1i1-f1i2
          f2i2=f1i2-f1i3
          f2i3=f1i3-f1i4
          f2i=2.*f2i1-f2i2
          if(ibmax.eq.2) then
            f1bmax=fmax(n)-f(i,n)
            f2i=f1bmax-f1i1
            f2bmax=2.*f2i-f2i1
          endif
          f2mi1=ammd(f2i1,f2i)
          f2mi2=ammd(f2i2,f2i1)
          f2mi3=ammd(f2i3,f2i2)
          at1=f1i1+0.5*f2mi1
          d1(i)=ammd(at1,2.*f1i1)
          if(ibmax.eq.2) then
            sl=ammd(f1bmax,f1i1)
            at=f1bmax-0.5*ammd(f2i1,f2bmax)
            ati1=f1i1+0.5*f2mi1
            t1=ammd(at,ati1)
            d1(i)=sign(1.,t1)*amin1(0.5*abs(at+ati1),
     &            amax1(2.*abs(sl),abs(t1)))
          endif
c
          ql(i,n)=f(i,n)+0.5*d1(i)
          qr(i,n)=f(i,n)-0.5*d1(i)
c
          i=ie1
          at=f1i1-0.5*f2mi1
          ati1=f1i2+0.5*f2mi2
          sl=ammd(f1i1,f1i2)
          t1=ammd(at,ati1)
          d1(i)=sign(1.,t1)*
     &          amin1(0.5*abs(at+ati1),amax1(2.*abs(sl),abs(t1)))
c
          ql(i,n)=f(i,n)+0.5*d1(i)
          qr(i,n)=f(i,n)-0.5*d1(i)
c
          i=ie2
          at=f1i2-0.5*f2mi2
          ati1=f1i3+0.5*f2mi3
          sl=ammd(f1i2,f1i3)
          t1=ammd(at,ati1)
          slopee2=sign(1.,t1)*
     &            amin1(0.5*abs(at+ati1),amax1(2.*abs(sl),abs(t1)))
c
          iend=ie2
        endif
c
        call tridag(con1,con2,con1,con3,d1,istart,iend)
c
c..curvature:
c
        if(ibmin.ne.3) then
          i=is
          d2(i)=0.
          i=is1
          d2(i)=0.
        else
          i=is
          con1(i)=2./11.
          con2(i)=1.
          con3(i)=(3.*f(i+2,n)+48.*f(i+1,n)-102.*f(i,n)+48.*fmin(n)
     &             +3.*fmin1(n))/44.
c
          i=is1
          con1(i)=2./11.
          con2(i)=1.
          con3(i)=(3.*f(i+2,n)+48.*f(i+1,n)-102.*f(i,n)+48.*f(i-1,n)
     &             +3.*fmin(n))/44.
c
          istart=is
        endif
c
        do i=is2,ie2
          con1(i)=2./11.
          con2(i)=1.
          con3(i)=(3.*f(i+2,n)+48.*f(i+1,n)-102.*f(i,n)+48.*f(i-1,n)
     &             +3.*f(i-2,n))/44.
        enddo
c
        if(ibmax.eq.3) then
          i=ie1
          con1(i)=2./11.
          con2(i)=1.
          con3(i)=(3.*fmax(n)+48.*f(i+1,n)-102.*f(i,n)+48.*f(i-1,n)
     &             +3.*f(i-2,n))/44.
c
          i=ie
          con1(i)=2./11.
          con2(i)=1.
          con3(i)=(3.*fmax1(n)+48.*fmax(n)-102.*f(i,n)+48.*f(i-1,n)
     &             +3.*f(i-2,n))/44.
c
          iend=ie
        else
          i=ie1
          d2(i)=0.
          i=ie
          d2(i)=0.
        endif
        call tridag(con1,con2,con1,con3,d2,istart,iend)
c
c..correct the sign of the slope and curvature:
c
        do i=is1,ie
          s1(i)=f(i,n)-f(i-1,n)
        enddo
        do i=is,ie2
          sf(i)=(-f(i+2,n)+4.*f(i+1,n)-3.*f(i,n))/2.
        enddo
        do i=is1,ie1
          sc(i)=(f(i+1,n)-f(i-1,n))/2.
        enddo
        do i=is2,ie
          sb(i)=(f(i-2,n)-4.*f(i-1,n)+3.*f(i,n))/2.
        enddo
        if(ibmin.eq.3) then
          s1(is)=f(is,n)-fmin(n)
          i=is-1
          sf(i)=(-f(i+2,n)+4.*f(i+1,n)-3.*fmin(n))/2.
          i=is
          sc(i)=(f(i+1,n)-fmin(n))/2.
          i=is-1
          sc(i)=(f(i+1,n)-fmin1(n))/2.
          i=is1
          sb(i)=(fmin(n)-4.*f(i-1,n)+3.*f(i,n))/2.
          i=is
          sb(i)=(fmin1(n)-4.*fmin(n)+3.*f(i,n))/2.
          i=is-1
          sb(i)=(fmin2(n)-4.*fmin1(n)+3.*fmin(n))/2.
        endif
        if(ibmax.eq.3) then
          s1(ie+1)=fmax(n)-f(ie,n)
          i=ie1
          sf(i)=(-fmax(n)+4.*f(i+1,n)-3.*f(i,n))/2.
          i=ie
          sf(i)=(-fmax1(n)+4.*fmax(n)-3.*f(i,n))/2.
          i=ie+1
          sf(i)=(-fmax2(n)+4.*fmax1(n)-3.*fmax(n))/2.
          i=ie
          sc(i)=(fmax(n)-f(i-1,n))/2.
          i=ie+1
          sc(i)=(fmax1(n)-f(i-1,n))/2.
          i=ie+1
          sb(i)=(f(i-2,n)-4.*f(i-1,n)+3.*fmax(n))/2.
        endif
c
        if(ibmin.eq.3) then
          i=is-1
          sm=amedian(sf(i),sc(i),sb(i))
          if(sm.eq.sf(i)) cm=f(i+2,n)-2.*f(i+1,n)+fmin(n)
          if(sm.eq.sc(i)) cm=f(i+1,n)-2.*fmin(n)+fmin1(n)
          if(sm.eq.sb(i)) cm=fmin(n)-2.*fmin1(n)+fmin2(n)
          d1(i)=amedian(d1(i),sm,sc(i))
c
          i=is
          sm=amedian(sf(i),sc(i),sb(i))
          if(sm.eq.sf(i)) cm=f(i+2,n)-2.*f(i+1,n)+f(i,n)
          if(sm.eq.sc(i)) cm=f(i+1,n)-2.*f(i,n)+fmin(n)
          if(sm.eq.sb(i)) cm=f(i,n)-2.*fmin(n)+fmin1(n)
          d1(i)=amedian(d1(i),sm,sc(i))
          if(d1(i).eq.sm) d2(i)=cm
          if(d1(i).eq.sc(i)) d2(i)=f(i+1,n)-2.*f(i,n)+fmin(n)
c
          i=is1
          sm=amedian(sf(i),sc(i),sb(i))
          if(sm.eq.sf(i)) cm=f(i+2,n)-2.*f(i+1,n)+f(i,n)
          if(sm.eq.sc(i)) cm=f(i+1,n)-2.*f(i,n)+f(i-1,n)
          if(sm.eq.sb(i)) cm=f(i,n)-2.*f(i-1,n)+fmin(n)
          d1(i)=amedian(d1(i),sm,sc(i))
          if(d1(i).eq.sm) d2(i)=cm
          if(d1(i).eq.sc(i)) d2(i)=f(i+1,n)-2.*f(i,n)+f(i-1,n)
        endif
c
        do i=is2,ie2
          sm=amedian(sf(i),sc(i),sb(i))
          if(sm.eq.sf(i)) cm=f(i+2,n)-2.*f(i+1,n)+f(i,n)
          if(sm.eq.sc(i)) cm=f(i+1,n)-2.*f(i,n)+f(i-1,n)
          if(sm.eq.sb(i)) cm=f(i,n)-2.*f(i-1,n)+f(i-2,n)
          d1(i)=amedian(d1(i),sm,sc(i))
          if(d1(i).eq.sm) d2(i)=cm
          if(d1(i).eq.sc(i)) d2(i)=f(i+1,n)-2.*f(i,n)+f(i-1,n)
        enddo
c
        if(ibmax.eq.3) then
          i=ie-1
          sm=amedian(sf(i),sc(i),sb(i))
          if(sm.eq.sf(i)) cm=fmax(n)-2.*f(i+1,n)+f(i,n)
          if(sm.eq.sc(i)) cm=f(i+1,n)-2.*f(i,n)+f(i-1,n)
          if(sm.eq.sb(i)) cm=f(i,n)-2.*f(i-1,n)+f(i-2,n)
          d1(i)=amedian(d1(i),sm,sc(i))
          if(d1(i).eq.sm) d2(i)=cm
          if(d1(i).eq.sc(i)) d2(i)=f(i+1,n)-2.*f(i,n)+f(i-1,n)
c
          i=ie
          sm=amedian(sf(i),sc(i),sb(i))
          if(sm.eq.sf(i)) cm=fmax1(n)-2.*fmax(n)+f(i,n)
          if(sm.eq.sc(i)) cm=fmax(n)-2.*f(i,n)+f(i-1,n)
          if(sm.eq.sb(i)) cm=f(i,n)-2.*f(i-1,n)+f(i-2,n)
          d1(i)=amedian(d1(i),sm,sc(i))
          if(d1(i).eq.sm) d2(i)=cm
          if(d1(i).eq.sc(i)) d2(i)=fmax(n)-2.*f(i,n)+f(i-1,n)
c
          i=ie+1
          sm=amedian(sf(i),sc(i),sb(i))
          if(sm.eq.sf(i)) cm=fmax2(n)-2.*fmax1(n)+fmax(n)
          if(sm.eq.sc(i)) cm=fmax1(n)-2.*fmax(n)+f(i-1,n)
          if(sm.eq.sb(i)) cm=fmax(n)-2.*f(i-1,n)+f(i-2,n)
          d1(i)=amedian(d1(i),sm,sc(i))
          if(d1(i).eq.sm) d2(i)=cm
          if(d1(i).eq.sc(i)) d2(i)=fmax1(n)-2.*fmax(n)+f(i-1,n)
        endif
c
        do i=istart,iend
          si=d1(i)
          ci=d2(i)
c..impose monotonicity-preserving constraints:
          if(abs(d1(i)).gt.2.0*abs(s1(i+1)).and.abs(d1(i)).gt.
     &       2.0*abs(s1(i))) then
            if(abs(d1(i)).le.abs(d1(i+1))) then
              sr=d1(i)
              cr=2.*(s1(i+1)-d1(i))
            else
              sr=2.*s1(i+1)-d1(i+1)
              cr=2.*(-s1(i+1)+d1(i+1))
            endif
            if(abs(d1(i)).le.abs(d1(i-1))) then
              sl=d1(i)
              cl=2.*(-s1(i)+d1(i))
            else
              sl=2.*s1(i)-d1(i-1)
              cl=2.*(s1(i)-d1(i-1))
            endif
            si=ammd(sl,sr)
            ci=ammd(cl,cr)
          else
            if(abs(d1(i)).gt.2.0*abs(s1(i+1))) then
              if(abs(d1(i)).le.abs(d1(i+1))) then
                sr=d1(i)
                cr=2.*(s1(i+1)-d1(i))
              else
                sr=2.*s1(i+1)-d1(i+1)
                cr=2.*(-s1(i+1)+d1(i+1))
              endif
              si=ammd(d1(i),sr)
              ci=ammd(d2(i),cr)
            endif
            if(abs(d1(i)).gt.2.0*abs(s1(i))) then
              if(abs(d1(i)).le.abs(d1(i-1))) then
                sl=d1(i)
                cl=2.*(-s1(i)+d1(i))
              else
                sl=2.*s1(i)-d1(i-1)
                cl=2.*(s1(i)-d1(i-1))
              endif
              si=ammd(sl,d1(i))
              ci=ammd(cl,d2(i))
            endif
          endif
c..piecewise quadratic reconstruction:
          ql(i,n)=f(i,n)+0.5*si+ci/12.
          qr(i,n)=f(i,n)-0.5*si+ci/12.
        enddo
c
c..consistency with the sonica scheme:
c
        if(ibmin.ne.3) qr(is2,n)=f(is2,n)-0.5*slopes2
        if(ibmax.ne.3) ql(ie2,n)=f(ie2,n)+0.5*slopee2
 1      continue
      endif
c
      return
      end

      function ammd4(w,x,y,z)

      real w,x,y,z,ammd4

      ammd4= 0.125*(sign(1.,w)+sign(1.,x))*
     <            abs((sign(1.,w)+sign(1.,y))*(sign(1.,w)+sign(1.,z)))*
     <            amin1(abs(w),abs(x),abs(y),abs(z))

      return
      end

c*************************************************************************
      subroutine shmp5(f,ql,qr,is,ie,im,th,qt,eps,fmin,fmax,
     <                  ibmin,ibmax,mdim)
c
c mp5 scheme of suresh and hyunh
c note: qr(is) and ql(ie) are never used.
c*************************************************************************

      implicit none

      integer is,ie,im,ibmax,ibmin,mdim
      
      real th,qt,eps
      real f(mdim,5),fmin(5),fmax(5)
      real ql(mdim,5),qr(mdim,5)
      real f0(mdim),f1(mdim),f2(mdim),f2m(mdim)

c..   local variables

      integer i,n

      real f0bmin,f1bmin,f2min,f2bmax,at,at1
      real s1,ati1,t1,f2bmin,f0bmax,f1bmax
      real ammd,smp5,slope

c..this is just 1st order upwind
      if(qt.eq.0.0)then
        do n=1,5
          do i=is,ie
            ql(i,n)=f(i,n)
            qr(i,n)=f(i,n)
          enddo
        enddo
        return
      else
c
c..the mp5 interpolation scheme of suresh and hyunh follows
c
        do 1 n=1,5
c
c..let's load up a few difference arrays
c
          do i=is,ie
            f0(i) = f(i,n)
          enddo
c..1st difference at i+1/2
          do i=is,ie-1
            f1(i) = f0(i+1) - f0(i)
          enddo
c..2nd difference at i
          do i=is+1,ie-1
            f2(i)  = f1(i) - f1(i-1)
          enddo
c..extrapolate at the boundaries
          f2(is) = 2.*f2(is+1)-f2(is+2)
          f2(ie) = 2.*f2(ie-1)-f2(ie-2)
c..modify at boundaries, if needed
          if(ibmin.eq.2) then
            f0bmin = fmin(n)
            f1bmin = f0(is) - f0bmin
            f2(is) = f1(is) - f1bmin
            f2bmin = 2.*f2(is)-f2(is+1)
          endif
          if(ibmax.eq.2) then
            f0bmax = fmax(n)
            f1bmax = f0bmax - f0(ie)
            f2(ie) = f1bmax - f1(ie-1)
            f2bmax = 2.*f2(ie)-f2(ie-1)
          endif
c..limit 2nd difference to i+1/2
          do i = is,ie-1
            f2m(i) = ammd(f2(i),f2(i+1))
          enddo
c
c..now combine everything to get ql and qr
c
c
c..first-order at boundary?
c
          at    = f1(is) - 0.5*f2m(is)
          slope = ammd(at,2.*f1(is))
          ql(is,n) = f0(is) + 0.5*slope
c
          if(ibmin.eq.2) then
            i = is
            s1    = ammd(f1(i),f1bmin)
            at    = f1(i)  - 0.5*f2m(i)
            ati1  = f1bmin + 0.5*ammd(f2(i),f2bmin)
            t1    = ammd(at,ati1)
            slope  =  sign(1.,t1)*
     <         amin1(0.5*abs(at+ati1),amax1(2.*abs(s1),abs(t1)))
            ql(i,n) = f0(i) + 0.5*slope
          endif
c
          at1   = f1(ie-1) + 0.5*f2m(ie-1)
          slope = ammd(at1,2.*f1(ie-1))
          qr(ie,n) = f0(ie) - 0.5*slope
c
          if(ibmax.eq.2) then
            i = ie
            s1    = ammd(f1bmax,f1(i-1))
            at    = f1bmax  - 0.5*ammd(f2(i-1),f2bmax)
            ati1  = f1(i-1) + 0.5*f2m(i-1)
            t1    = ammd(at,ati1)
            slope =  sign(1.,t1)*
     <         amin1(0.5*abs(at+ati1),amax1(2.*abs(s1),abs(t1)))
            qr(i,n) = f0(i) - 0.5*slope
          endif
c
c..sonic-a near the boundary?
c
          do i=is+1,ie-1,ie-is-2
c..include limited curvatures to calculate new slopes
            at    = f1(i)   - 0.5*f2m(i)
            ati1  = f1(i-1) + 0.5*f2m(i-1)
c..limit slopes at i
            s1    = ammd(f1(i),f1(i-1))
            t1    = ammd(at,ati1)
c..now find appropriate slope
            slope  =  sign(1.,t1)*
     <         amin1(0.5*abs(at+ati1),amax1(2.*abs(s1),abs(t1)))
c..use slope to calculate ql and qr
            ql(i,n) = f0(i) + 0.5*slope
            qr(i,n) = f0(i) - 0.5*slope
          enddo
c
c..suresh at interior
c
          do i=is+2,ie-2
          ql(i,n) = smp5(f0(i-2),f0(i-1),f0(i),f0(i+1),f0(i+2))
          qr(i,n) = smp5(f0(i+2),f0(i+1),f0(i),f0(i-1),f0(i-2))
          enddo
c
c..sonic-a near the boundary?
c
          i=is+2
c..include limited curvatures to calculate new slopes
            at    = f1(i)   - 0.5*f2m(i)
            ati1  = f1(i-1) + 0.5*f2m(i-1)
c..limit slopes at i
            s1    = ammd(f1(i),f1(i-1))
            t1    = ammd(at,ati1)
c..now find appropriate slope
            slope  =  sign(1.,t1)*
     <         amin1(0.5*abs(at+ati1),amax1(2.*abs(s1),abs(t1)))
c..use slope to calculate ql and qr
            ql(i,n) = f0(i) + 0.5*slope
c
c..sonic-a near the boundary?
c
          i=ie-1
c..include limited curvatures to calculate new slopes
            at    = f1(i)   - 0.5*f2m(i)
            ati1  = f1(i-1) + 0.5*f2m(i-1)
c..limit slopes at i
            s1    = ammd(f1(i),f1(i-1))
            t1    = ammd(at,ati1)
c..now find appropriate slope
            slope  =  sign(1.,t1)*
     <         amin1(0.5*abs(at+ati1),amax1(2.*abs(s1),abs(t1)))
c..use slope to calculate ql and qr
            qr(i,n) = f0(i) - 0.5*slope
    1   continue
      endif
c
      return
      end

      function smp5(a,b,c,d,e)

      implicit none

      real a,b,c,d,e
      real ammd,ammd4
      
      real b1,b2,alpha,epsm,vor,vmp,djm1,dj,djp1,dm4jph,dm4jmh
      real vul,vav,vmd,vlc,vmin,vmax,smp5

      b1 = 1./60.
      b2 = 4./3.
      alpha = 4.
      epsm = 1.e-8
      vor = b1*(2*a-13*b+47*c+27*d-3*e)
      vmp = c + ammd(d-c,alpha*(c-b))
      if( (vor-c)*(vor-vmp) .le. epsm) then
        smp5=vor
      else
        djm1 = a-2.*b+c
        dj   = b-2.*c+d
        djp1 = c-2.*d+e
        dm4jph = ammd4(4.*dj-djp1,4.*djp1-dj,dj,djp1)
        dm4jmh = ammd4(4.*dj-djm1,4.*djm1-dj,dj,djm1)
        vul = c + alpha*(c-b)
        vav = 0.5*(c+d)
        vmd = vav - 0.5*dm4jph
        vlc = c + 0.5*(c-b) +b2*dm4jmh
        vmin = amax1(amin1(c,d,vmd),amin1(c,vul,vlc))
        vmax = amin1(amax1(c,d,vmd),amax1(c,vul,vlc))
        smp5=vor+ammd(vmin-vor,vmax-vor)
      endif
      return
      end


c*************************************************************************
      subroutine weno(f,ql,qr,is,ie,im,th,qt,eps,fmin,fmax,
     <                  ibmin,ibmax,mdim)
c
c 5th order weno scheme
c note: qr(is) and ql(ie) are never used.
c*************************************************************************
      
      implicit none
      
      integer is,ie,im,ibmin,ibmax,mdim
      
      real th,qt,eps
      real f(mdim,5),fmin(5),fmax(5)
      real ql(mdim,5),qr(mdim,5)

c..   local variables

      integer i,n

      real f0(mdim),f1(mdim),f2(mdim),f2m(mdim),slope(mdim)
      real ammd,ammd4,weno5

      real f0bmin,f1bmin,f2bmin,f0bmax,f1bmax,f2bmax
      real at,s1,ati1,t1,at1


c..this is just 1st order upwind
      if(qt.eq.0.0)then
        do n=1,5
          do i=is,ie
            ql(i,n)=f(i,n)
            qr(i,n)=f(i,n)
          enddo
        enddo
        return
      else
c
c 5th order weno scheme follows
c
        do 1 n=1,5
c
c..let's load up a few difference arrays
c
          do i=is,ie
            f0(i) = f(i,n)
          enddo
c..1st difference at i+1/2
          do i=is,ie-1
            f1(i) = f0(i+1) - f0(i)
          enddo
c..2nd difference at i
          do i=is+1,ie-1
            f2(i)  = f1(i) - f1(i-1)
          enddo
c..extrapolate at the boundaries
          f2(is) = 2.*f2(is+1)-f2(is+2)
          f2(ie) = 2.*f2(ie-1)-f2(ie-2)
c..modify at boundaries, if needed
          if(ibmin.eq.2) then
            f0bmin = fmin(n)
            f1bmin = f0(is) - f0bmin
            f2(is) = f1(is) - f1bmin
            f2bmin = 2.*f2(is)-f2(is+1)
          endif
          if(ibmax.eq.2) then
            f0bmax = fmax(n)
            f1bmax = f0bmax - f0(ie)
            f2(ie) = f1bmax - f1(ie-1)
            f2bmax = 2.*f2(ie)-f2(ie-1)
          endif
c..limit 2nd difference to i+1/2
          do i = is,ie-1
            f2m(i) = ammd(f2(i),f2(i+1))
          enddo
c
c..now combine everything to get ql and qr
c
c
c..first-order at boundary?
c
          at    = f1(is) - 0.5*f2m(is)
          slope(i) = ammd(at,2.*f1(is))
          ql(is,n) = f0(is) + 0.5*slope(i)
c
          if(ibmin.eq.2) then
            i = is
            s1    = ammd(f1(i),f1bmin)
            at    = f1(i)  - 0.5*f2m(i)
            ati1  = f1bmin + 0.5*ammd(f2(i),f2bmin)
            t1    = ammd(at,ati1)
            slope(i)  =  sign(1.,t1)*
     <         amin1(0.5*abs(at+ati1),amax1(2.*abs(s1),abs(t1)))
            ql(i,n) = f0(i) + 0.5*slope(i)
          endif
c
          at1   = f1(ie-1) + 0.5*f2m(ie-1)
          slope(i) = ammd(at1,2.*f1(ie-1))
          qr(ie,n) = f0(ie) - 0.5*slope(i)
c
          if(ibmax.eq.2) then
            i = ie
            s1    = ammd(f1bmax,f1(i-1))
            at    = f1bmax  - 0.5*ammd(f2(i-1),f2bmax)
            ati1  = f1(i-1) + 0.5*f2m(i-1)
            t1    = ammd(at,ati1)
            slope(i)  =  sign(1.,t1)*
     <         amin1(0.5*abs(at+ati1),amax1(2.*abs(s1),abs(t1)))
            qr(i,n) = f0(i) - 0.5*slope(i)
          endif
c
c..sonic-a near the boundary?
c
          do i=is+1,ie-1,ie-is-2
c..include limited curvatures to calculate new slopes
            at    = f1(i)   - 0.5*f2m(i)
            ati1  = f1(i-1) + 0.5*f2m(i-1)
c..limit slopes at i
            s1    = ammd(f1(i),f1(i-1))
            t1    = ammd(at,ati1)
c..now find appropriate slope
            slope(i)  =  sign(1.,t1)*
     <         amin1(0.5*abs(at+ati1),amax1(2.*abs(s1),abs(t1)))
c..use slope to calculate ql and qr
            ql(i,n) = f0(i) + 0.5*slope(i)
            qr(i,n) = f0(i) - 0.5*slope(i)
          enddo
c
c..suresh at interior
c
          do i=is+2,ie-2
          ql(i,n) = weno5(f0(i-2),f0(i-1),f0(i),f0(i+1),f0(i+2))
          qr(i,n) = weno5(f0(i+2),f0(i+1),f0(i),f0(i-1),f0(i-2))
          enddo
c
c..sonic-a near the boundary?
c
          i=is+2
c..include limited curvatures to calculate new slopes
            at    = f1(i)   - 0.5*f2m(i)
            ati1  = f1(i-1) + 0.5*f2m(i-1)
c..limit slopes at i
            s1    = ammd(f1(i),f1(i-1))
            t1    = ammd(at,ati1)
c..now find appropriate slope
            slope(i)  =  sign(1.,t1)*
     <         amin1(0.5*abs(at+ati1),amax1(2.*abs(s1),abs(t1)))
c..use slope to calculate ql and qr
            ql(i,n) = f0(i) + 0.5*slope(i)
c
c..sonic-a near the boundary?
c
          i=ie-1
c..include limited curvatures to calculate new slopes
            at    = f1(i)   - 0.5*f2m(i)
            ati1  = f1(i-1) + 0.5*f2m(i-1)
c..limit slopes at i
            s1    = ammd(f1(i),f1(i-1))
            t1    = ammd(at,ati1)
c..now find appropriate slope
            slope(i)  =  sign(1.,t1)*
     <         amin1(0.5*abs(at+ati1),amax1(2.*abs(s1),abs(t1)))
c..use slope to calculate ql and qr
            qr(i,n) = f0(i) - 0.5*slope(i)
    1   continue
      endif
c
      return
      end

c..   function weno5

      function weno5(a,b,c,d,e)

      implicit none

      real a,b,c,d,e
      real b1,b2,epsw,djm1,ejm1,dj,ej
      real djp1,ejp1,dis0,dis1,dis2,q30,q31,q32
      real d01,d02,a1ba0,a2ba0,w0,w1,w2,weno5
      
      b1 = 13./12. 
      b2 = 1./6.
      epsw = 1.e-6
      djm1 = a-2.*b+c
      ejm1 = a-4.*b+3.*c
      dj   = b-2.*c+d
      ej   = b-d
      djp1 = c-2.*d+e
      ejp1 = 3.*c-4.*d+e
      dis0 = b1*djm1*djm1+0.25*ejm1*ejm1+epsw
      dis1 = b1*dj*dj+0.25*ej*ej+epsw
      dis2 = b1*djp1*djp1+0.25*ejp1*ejp1+epsw
      q30 = 2.*a-7.*b+11.*c
      q31 = -b+5.*c+2.*d
      q32 = 2.*c+5.*d-e
      d01 = dis0/dis1
      d02 = dis0/dis2
      a1ba0 = 6.*d01
      a2ba0 = 3.*d02
      w0 = 1./(1.+a1ba0+a2ba0)
      w1 = a1ba0*w0
      w2 = 1.-w0-w1
      weno5 = b2*(w0*q30+w1*q31+w2*q32)
      
      return
      end



c_____________________________________________________________________________

      subroutine tridag(a,b,c,f,z,ni,nl)
c_____________________________________________________________________________
      
      implicit none

      integer ni,nl,j,nd,j1

      real a(nl),b(nl),c(nl),f(nl),z(nl)
      real w(nl),g(nl)
      real nipl,d,rd
c
      w(ni)=c(ni)/b(ni)
      g(ni)=f(ni)/b(ni)
      nipl=ni+1
      do 10 j=nipl,nl
      d=b(j)-a(j)*w(j-1)
      rd=1.0/d
      w(j)=c(j)*rd
      g(j)=(f(j)-a(j)*g(j-1))*rd
 10   continue
      z(nl)=g(nl)
      nd=nl-ni
      do 20 j1=1,nd
      j=nl-j1
      z(j)=g(j)-w(j)*z(j+1)
 20   continue
      return
      end
