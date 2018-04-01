c***********************************************************************
      subroutine step(x,y,z,iblank,q,qnewt,s,turmu,vnu,vnu0,
     <     xx,xy,xz,ug,yx,yy,yz,vg,zx,zy,zz,wg,zx0,zy0,zz0,zt0,bt,
     <     xaj,yaj,zaj,xak,yak,zak,xal,yal,zal,vnaj,vnak,vnal,
     <     kkp,kkr,wgrid,wvec,resmax,resrho,rsum,grbc)
c
c  note that the boundaries are updated explicitly
c
c***********************************************************************

      use params_global
      use arf_mod
      use dualtime
      use bcparam

      implicit none

      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)
      real q(jmax,kmax,lmax,nd),s(jmax,kmax,lmax,nv)
      real qnewt(jmax,kmax,lmax,nd)
      real resmax,resrho,rsum
      integer iblank(jmax,kmax,lmax)
      integer imesh

      real xx(jmax,kmax,lmax),xy(jmax,kmax,lmax),xz(jmax,kmax,lmax),
     <     yx(jmax,kmax,lmax),yy(jmax,kmax,lmax),yz(jmax,kmax,lmax),
     <     zx(jmax,kmax,lmax),zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)

      real xaj(jmax,kmax,lmax), yaj(jmax,kmax,lmax), zaj(jmax,kmax,lmax)
      real xak(jmax,kmax,lmax), yak(jmax,kmax,lmax), zak(jmax,kmax,lmax)
      real xal(jmax,kmax,lmax), yal(jmax,kmax,lmax), zal(jmax,kmax,lmax)
      real vnaj(jmax,kmax,lmax), vnak(jmax,kmax,lmax), vnal(jmax,kmax,lmax)

      real zx0(jmax,kmax),zy0(jmax,kmax),zz0(jmax,kmax),zt0(jmax,kmax)

      real turmu(jmax,kmax,lmax),vnu(jmax,kmax,lmax)
      real vnu0(jmax,kmax,lmax),bt(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax),vg(jmax,kmax,lmax),wg(jmax,kmax,lmax)
      real wgrid(jmax,kmax,lmax),wvec(mdim,25)
      integer kkp(kmax),kkr(kmax)
      type(bc_t) :: grbc

c...  local variables

      integer j,k,l,ka,kb,jd1,kd1,ld1,n,nd1,igrid
      real oat,smrs,smrs1,smrs2,smrs3,smrs4,smrs5,volum
      real, allocatable :: vmul(:,:,:), tscale(:,:,:)

      integer jc, kc, lc
      real gcl_l2norm, gcl_rmax

c***** first executable statement

      allocate(vmul(jmax,kmax,lmax))
      allocate(tscale(jmax,kmax,lmax))

      ka = 2 - ksym
      kb = km + ksym


c..   zero s array 

      do  10 l = 1,lmax
      do  10 k = 1,kmax
      do  10 j = 1,jmax
        s(j,k,l,1) = 0.
        s(j,k,l,2) = 0.
        s(j,k,l,3) = 0.
        s(j,k,l,4) = 0.
        s(j,k,l,5) = 0.
   10 continue

      oat = 1.0
      if (ntac.eq.2 .and. istep.gt.1) oat = 2./3.
c..   Now compute local time scaling according to dual time or Newton...
      !print*, "Computing tscale, idual = ", idual
      if(idual .eq. 0)then
        call newton_timescale(timeac, oat, h, q, iblank, tscale, bt, jmax, kmax, lmax)
      else
        call dualtime_timescale(oat, h, q, iblank, tscale, bt, iprecon, jmax, kmax, lmax)
      end if

c..   compute the right hand side and store in s array
      if (iunst.eq. FL_UNSTEADY) then

          if(usegcl)then
             print*, "CALLING RHSUP_GCL..."
             call rhsup_gcl(q,s,xaj,yaj,zaj,xak,yak,zak,xal,yal,zal,
     <            vnaj,vnak,vnal,kkr,kkp,iblank)
             print*, "...RHSUP_GCL RETURNED."
          else
             call rhsup1( q,s,xx,xy,xz,yx,yy,yz,zx,zy,zz,
     $            x,y,z,iblank,ug,vg,wg,kkr,kkp)
          end if

         if(arf_opt .eq. 1)then
           !Add source terms to rhs
           print*, "step: using ARF source terms"
           call add_arf_source(q,s)
         end if

      elseif (iunst.eq. FL_HOVER) then
         call rhsup1( q,s,xx,xy,xz,yx,yy,yz,zx,zy,zz,
     $        x,y,z,iblank,ug,vg,wg,kkr,kkp)
         call source( q,s)
      else
         call rhsup1( q,s,xx,xy,xz,yx,yy,yz,zx,zy,zz,
     $        x,y,z,iblank,ug,vg,wg,kkr,kkp)

c.. vorticity confinement
         if (iconfine.eq.1) then
            call confine(q,s,xx,xy,xz,yx,yy,yz,zx,zy,zz)
         endif

      endif

c..   compute viscous fluxes
      if( .not. invisc) then
         if (.not.ithin) then
            call fulvisrhs(vnu,vnu0,turmu,vmul,q,s,xx,xy,xz,yx,yy,yz,zx,zy,
     $                  zz,x,y,z,ug,vg,wg,tscale,kkr,kkp,iblank,grbc)
         else
            call visrhs(x,y,z,iblank,vnu,vnu0,vmul,turmu,q,s,xx,xy,xz,ug,
     $           yx,yy,yz,vg,zx,zy,zz,wg,tscale,kkr,kkp,grbc)
         endif
	     !write(STDOUT,*) 'in viscous ', iwm
      else
         do l=1,lmax
            do k=1,kmax
               do j=1,jmax
                  turmu(j,k,l)=0.
               enddo
            enddo
         enddo
      endif

	    !insurance
	    jm=jmax-1
	    km=kmax-1
	    lm=lmax-1

c..   make sure to account for time accuracy
     
c      do  15 l = 1,lmax
c      do  15 k = 1,kmax
c      do  15 j = 1,jmax
c        s(j,k,l,1) = s(j,k,l,1)*oat
c        s(j,k,l,2) = s(j,k,l,2)*oat
c        s(j,k,l,3) = s(j,k,l,3)*oat
c        s(j,k,l,4) = s(j,k,l,4)*oat
c        s(j,k,l,5) = s(j,k,l,5)*oat
c 15   continue
c
c..   start newton iteration here at each time step for convergence
c..   Still call newton if dual time stepping; same RHS (BS)
      
      if(itnmax.gt.1) call newton(q,qnewt,s,ka,kb)

      if(iprecon) call rhslom(q,s,2,jm,ka,kb,2,lm,bt)

      rsum = 0.
      resrho  = 0.0
      resmax  = 0.0
      !print*, "[in step] jm, ka, kb, lm", jm, ka, kb, lm
      do 23 l = 2,lm
      do 23 k = ka,kb

        do 22 n = 1,nv
        do 22 j = 2,jm
          smrs = s(j,k,l,n)*max(iblank(j,k,l),0)
          if (smrs .gt. resmax) then
          jd1 = j
          kd1 = k
          ld1 = l
          nd1 = n
          resmax = smrs
          endif
  22    continue

        do 24 j = 2,jm
          smrs1 = s(j,k,l,1)*max(iblank(j,k,l),0)
          smrs2 = s(j,k,l,2)*max(iblank(j,k,l),0)
          smrs3 = s(j,k,l,3)*max(iblank(j,k,l),0)
          smrs4 = s(j,k,l,4)*max(iblank(j,k,l),0)
          smrs5 = s(j,k,l,5)*max(iblank(j,k,l),0)
          resrho  = resrho + smrs1**2
          rsum = rsum + smrs1*smrs1 + smrs2*smrs2 + smrs3*smrs3
     <                + smrs4*smrs4 + smrs5*smrs5
          s(j,k,l,1) = smrs1*tscale(j,k,l)
          s(j,k,l,2) = smrs2*tscale(j,k,l)
          s(j,k,l,3) = smrs3*tscale(j,k,l)
          s(j,k,l,4) = smrs4*tscale(j,k,l)
          s(j,k,l,5) = smrs5*tscale(j,k,l)
 24    continue
 23   continue

      volum  = float((jm-1)*(km-1)*(lm-1))
      volum = float(jmax*kmax*lmax*5)

      resrho = sqrt(resrho/volum)
      rsum   = sqrt(rsum/volum)

      if( mod(istep,npnorm).eq.0 ) then
         write(STDOUT,102) jd1,kd1,ld1,nd1,resmax,resrho,rsum
	 write(reshst,'(I6,I4,4(x,E12.5))')istep0,itn,totime,rsum,resmax,resrho
	 call flush(6)
      endif
  102 format('  j,k,l,n,rmax,l2rho,l2 =',4i4,3(1x,e10.4))

c..   check on convergence

      if( rsum .gt.1000.0 ) then
          write(STDOUT,602)  rsum
  602     format(' ',10x,'norm is out of bounds,'
     $           ,1x,'l2norm = ',f16.12,1x,'solution suspended' )
          stop 'norm'
      end if

c..   now do the implicit part
!       call ilu3d_old(q,s,2,jm,ka,kb,2,lm
!     &      ,xx,xy,xz,yx,yy,yz,zx,zy,zz,ug,vg,wg,turmu,kkr,kkp,tscale/oat)
      if(ilhs.eq.1) then
        call ilu3d(q,s,2,jm,ka,kb,2,lm
     &       ,xx,xy,xz,yx,yy,yz,zx,zy,zz,ug,vg,wg,iblank,
     $      turmu,tscale,kkr,kkp)
      elseif (ilhs.eq.2 .or. ilhs.eq.3) then
         if (.not.iprecon) then
            call arc3d(q,s,2,jm,ka,kb,2,lm,xx,xy,xz,yx,yy,yz,zx,zy,zz,
     &                 ug,vg,wg,iblank,turmu,tscale,bt,kkr,kkp)
         else
            call arc3d_precon(q,s,2,jm,ka,kb,2,lm,xx,xy,xz,yx,yy,yz,
     &                 zx,zy,zz,ug,vg,wg,iblank,turmu,tscale,bt,kkr,kkp)
         endif
      endif
!        !insurance
	jm=jmax-1
	km=kmax-1
	lm=lmax-1
      

c..   update q with corrections
      
      do 31 l = 2,lm
      do 31 k = ka,kb
      do 31 j = 2,jm
        q(j,k,l,1) = q(j,k,l,1) + s(j,k,l,1)*max(iblank(j,k,l),0)
        q(j,k,l,2) = q(j,k,l,2) + s(j,k,l,2)*max(iblank(j,k,l),0)
        q(j,k,l,3) = q(j,k,l,3) + s(j,k,l,3)*max(iblank(j,k,l),0)
        q(j,k,l,4) = q(j,k,l,4) + s(j,k,l,4)*max(iblank(j,k,l),0)
        q(j,k,l,5) = q(j,k,l,5) + s(j,k,l,5)*max(iblank(j,k,l),0)
   31 continue

c..   update bc

c      call bc(x,y,z,q,xx,xy,xz,ug,yx,yy,yz,vg,zx,zy,zz,wg,
c     <     zx0,zy0,zz0,zt0,kkr,kkp)


!BS..  Let's properly clean up memory before exiting...

      deallocate(vmul)
      deallocate(tscale)

      return
      end

C***********************************************************************
      SUBROUTINE MONITOR(MSTOP,Q)
C
C  This subroutine checks for negative speed of sound
C
C***********************************************************************

      use params_global
      implicit none

      integer mstop
      real q(jmax,kmax,lmax,nd)

      real qtemp(6)
      integer l,k,j
      integer ja,jb,la,lb,ka,kb,nneg
      real asqmin,rho,qq,asq

C***********************************

      if (is_wing) then
       KA = 1 
       KB = KMAX
   
       JA=1
       JB=JMAX

       LA=1
       LB=LMAX
      else
	KA=2
	KB=KM
	JA=2
	JB=JM
	LA=2
	LB=LM
      endif
C
C..Negative speed of sound check
C
      MSTOP = 0
      ASQMIN = 0.00002
      NNEG = 0
      DO 900 L = LA,LB
       DO 910 K = KA,KB
        DO 920 J = JA,JB
          RHO = Q(J,K,L,1)
          QQ = Q(J,K,L,2)**2 +Q(J,K,L,3)**2 +Q(J,K,L,4)**2
          ASQ = GGM1*(Q(J,K,L,5) - .5*QQ/RHO)/RHO
CCRAY          NNEG = CVMGT(NNEG+1,NNEG,ASQ.LE.ASQMIN)
          IF(ASQ.LE.ASQMIN) NNEG = NNEG+1
 920    CONTINUE
 910   CONTINUE
 900  CONTINUE
C
      IF(NNEG .NE. 0) THEN

        WRITE(STDOUT,*) NNEG, ' NEGATIVE SPEED OF SOUND'
        write(STDERR,*) NNEG, ' negative speed of sound!!!'
        DO 74 L = LA,LB 
        DO 74 K = KA,KB
        DO 74 J = JA,JB
          RHO = Q(J,K,L,1)
          QQ = Q(J,K,L,2)**2 +Q(J,K,L,3)**2 +Q(J,K,L,4)**2
          ASQ = GGM1*(Q(J,K,L,5) - .5*QQ/RHO)/RHO
          IF( ASQ.LE.ASQMIN ) THEN
c$$$       	    write(*,*)Q(J,K,L,5),rho,.5*QQ/RHO,q(j,k,l,2),q(j,k,l,3),
c$$$     &                q(j,k,l,4)
            WRITE(STDOUT,601) J,K,L,ASQ
            write(STDOUT,*) 'q(j,k,l)=',q(j,k,l,:)
            write(STDOUT,*) 'fixing using previous value for now :'
            q(j,k,l,1:5)=qtemp(1:5)*qtemp(6)/q(j,k,l,6)
          ELSE
             QTEMP=Q(j,k,l,:)
          ENDIF
   74   CONTINUE

        if (nneg.gt.10) mstop=1
           
      END IF
C
 601  FORMAT(' J,K,L,ASQ FROM MONITOR ',3I5,F12.5)
C
      RETURN
      END



C***********************************************************************
      SUBROUTINE ILU3D_old(Q,S,JS,JE,KS,KE,LS,LE
     &      ,XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ,UG,VG,WG,TURMU,KKR,KKP,TSCALE)
C
C  Calculate the implicit inversion of the LHS
C  This involves two bidiagonal scalar inversions
C
C***********************************************************************

      use params_global
      implicit none


      integer js,je,ks,ke,ls,le

      real q(jmax,kmax,lmax,nd),s(jmax,kmax,lmax,nv),
     $     turmu(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax),vg(jmax,kmax,lmax),wg(jmax,kmax,lmax)
      real xx(jmax,kmax,lmax), xy(jmax,kmax,lmax) ,xz(jmax,kmax,lmax),
     $     yx(jmax,kmax,lmax), yy(jmax,kmax,lmax) ,yz(jmax,kmax,lmax),
     $     zx(jmax,kmax,lmax), zy(jmax,kmax,lmax) ,zz(jmax,kmax,lmax)
      real tscale(jmax,kmax,lmax)
      integer kkr(kmax),kkp(kmax)

c..   local variables

      real d(jmax,kmax,lmax)
      real a(mdim,5),c(mdim,5)
      real uv(jmax,lmax),vn(jmax,lmax),ge(jmax,lmax)
      real qx(jmax,lmax),qz(jmax,lmax),cx(jmax,lmax),cz(jmax,lmax)
      real b(jmax,lmax,5)
      real vn_j(jmax,lmax),vn_k(jmax,lmax),vn_l(jmax,lmax)
      
      real tcos,tsin,eps2,epsv,oat,dj,f1,scale,scale1,scale2
      real rhou,rhov,uu,vv,ww,er,uvw,cjkl,rj1,rj2,rj3
      real qq1,qqx,rr1,rl1,rl2,rl3,qq3,qqz,rr3,vnu,svt
      real ri1,ri2,ri3,ri4,qq,cc,qqy,sp1,sm2,spec
      real rk1,rk2,rk3,rr2,chkx,chky,chkz
      real s1,s2,s3,s4,s5,abcd,r2,r3,sv,di
      real c1, c2, c3, c4, c5, c6, c7, c8
      real visc_j,visc_k,visc_l
      real vnu_j,vnu_k,vnu_l,vnud
c     
      integer j,k,l,jj,llo,l1,m,lup,i,kr,j1
      integer ms(mdim*2),me(mdim*2)

C***********************************
C... viscous terms (this affects the vnu term)

      visc_j=1.0
      visc_k=1.0
      visc_l=1.0

      JM= JMAX - 1
      KM= KMAX - 1
      LM= LMAX - 1

      TCOS = COS(BLANG)
      TSIN = SIN(BLANG)
C
      EPS2 = EPSE*2.
      EPSV = 1. + EPS2
C
C..Set-up for hyper-plane loop
C
      DO 1 M=JS+LS,JE+LE
      I     = MAX(JS,M-LE)
      MS(M) = I
      ME(M) = MIN(JE,M-JS)
    1 CONTINUE
C
C..Set time-accuracy
C
      OAT = 1.0
      IF (NTAC.EQ.2 .AND. ISTEP.GT.1) OAT = 2./3.
C
C..Store DENSITY,U,V,W in nonconservative variables                  
C
      DO 2 L = LS-1, LE+1
      DO 2 K = KS-1, KE+1
      DO 2 J = JS-1, JE+1
      DJ         = 1.0 / Q(J,K,L,1)
      Q(J,K,L,2) = Q(J,K,L,2)*DJ
      Q(J,K,L,3) = Q(J,K,L,3)*DJ
      Q(J,K,L,4) = Q(J,K,L,4)*DJ
      Q(J,K,L,5) = Q(J,K,L,5)*DJ
      Q(J,K,L,1) = Q(J,K,L,1)*Q(J,K,L,6)
    2 CONTINUE
C..In the wake, average points above and below
C..for faster convergence?
C
C  REMOVED
C..for faster convergence at front and back?
C
C..Forward sweep
C..Setup D, the diagonal term and B contribution in K-direction
C
      DO 100 K = KS,KE
C
        DO 111 L = LS-1,LE+1
        DO 111 J = JS-1,JE+1
          UU  = Q(J,K,L,2)
          VV  = Q(J,K,L,3)
          WW  = Q(J,K,L,4)
          UVW = 0.5*( UU*UU + VV*VV + WW*WW )
          CJKL = SQRT( GGM1*( Q(J,K,L,5) - UVW ) ) !speed of sound ^ 2
C
          RJ1 = XX(J,K,L)
          RJ2 = XY(J,K,L)
          RJ3 = XZ(J,K,L)
          QQ1 = RJ1*UU + RJ2*VV + RJ3*WW
          QQX = ABS( -RJ1*UG(J,K,L)-RJ2*VG(J,K,L)-RJ3*WG(J,K,L) + QQ1 )
          RR1 = SQRT( RJ1**2 + RJ2**2 + RJ3**2 )
          RK1 = YX(J,K,L)
          RK2 = YY(J,K,L)
          RK3 = YZ(J,K,L)
          QQY = ABS( RK1*(UU-UG(J,K,L)) + RK2*(VV-VG(J,K,L))
     <             + RK3*(WW-WG(J,K,L)) )
          RR2 = SQRT( RK1**2 + RK2**2 + RK3**2 )
          RL1 = ZX(J,K,L)
          RL2 = ZY(J,K,L)
          RL3 = ZZ(J,K,L)
          QQ3 = RL1*UU + RL2*VV + RL3*WW
          QQZ = ABS( -RL1*UG(J,K,L)-RL2*VG(J,K,L)-RL3*WG(J,K,L) + QQ3 )
          RR3 = SQRT( RL1**2 + RL2**2 + RL3**2 )

          VNU_J = visc_j*2.0*RR1*RR1*(RMUE+TURMU(J,K,L))/(REY*Q(J,K,L,1)
     <		)

          VNU_K = visc_k*2.0*RR2*RR2*(RMUE+TURMU(J,K,L))/(REY*Q(J,K,L,1)
     <		)

          VNU_L = visc_l*2.0*RR3*RR3*(RMUE+TURMU(J,K,L))/(REY*Q(J,K,L,1)
     <		)

	  VNUD  = VNU_J+VNU_K+VNU_L

          SVT = OAT*TSCALE(J,K,L)
          D(J,K,L)=1.0/
     &     ( 1.0+SVT*((QQX+QQY+QQZ +CJKL*(RR1+RR2+RR3))*EPSV +VNUD) )
C
          UV(J,L) = UVW
          QX(J,L) = QQ1
          QZ(J,L) = QQ3
          CX(J,L) = CJKL*RR1
          CZ(J,L) = CJKL*RR3
          VN_J(J,L) = VNU_J
          VN_K(J,L) = VNU_K
          VN_L(J,L) = VNU_L
          GE(J,L) = GAMMA*Q(J,K,L,5) - GM1*UVW
  111   CONTINUE
C  
        KR=K-1
        DO 112 L = LS-1,LE+1
        DO 112 J = JS-1,JE+1
          UU  = Q(J,KR,L,2)
          VV  = Q(J,KR,L,3)
          WW  = Q(J,KR,L,4)
          ER  = Q(J,KR,L,5)
          UVW = 0.5*( UU*UU + VV*VV + WW*WW )
          RI1 = YX(J,KR,L)
          RI2 = YY(J,KR,L)
          RI3 = YZ(J,KR,L)
          RI4 = -RI1*UG(J,KR,L) - RI2*VG(J,KR,L) - RI3*WG(J,KR,L)
          RR2 = RI1**2 + RI2**2 + RI3**2
          QQ  = RI1*UU + RI2*VV + RI3*WW
          CC  = SQRT( GGM1*( ER-UVW )*RR2 )
          QQY = RI4 + QQ
          SP1 = ABS(QQY) + CC
          SM2 = EPS2*SP1
          CHKY= 0.5 + SIGN( 0.5, QQY+CC )
          SPEC= CHKY*( SP1 ) + VN_K(J,L) + SM2 
C
          S1 = S(J,KR,L,1)
          S2 = S(J,KR,L,2)
          S3 = S(J,KR,L,3)
          S4 = S(J,KR,L,4)
          S5 = S(J,KR,L,5)
C
          C1 = RI1*UU + RI2*VV + RI3*WW
          C2 = GM1*UVW
          C3 = RI4 + C1
          C4 = C1*S1
          C5 = UU*S2 + VV*S3 + WW*S4
          C6 = RI1*S2 + RI2*S3 + RI3*S4
          C7 = C6 - C4
          C8 = C2*S1 - GM1*(C5 - S5)
C 
          B(J,L,1) = SPEC*S1 + CHKY*(RI4*S1 + C6)
          B(J,L,2) = SPEC*S2 + CHKY*(C3*S2 + RI1*C8 + UU*C7)
          B(J,L,3) = SPEC*S3 + CHKY*(C3*S3 + RI2*C8 + VV*C7)
          B(J,L,4) = SPEC*S4 + CHKY*(C3*S4 + RI3*C8 + WW*C7)
          B(J,L,5) = SPEC*S5 + CHKY*( (RI4 + GAMMA*C1)*S5
     <                       + (GAMMA*ER - C2)*C7 + C2*C4 - GM1*C1*C5 )
  112   CONTINUE
C..Let's precondition at the wake cut and front overlap!
C
C..Loop on hyper-plane
C
      DO 120 M = JS+LS,JE+LE
C
C..Setup A contribution in J-direction
C
      DO 121 J = MS(M),ME(M)
        L  = M-J
        J1 = J-1
C
        UU  = Q(J1,K,L,2)
        VV  = Q(J1,K,L,3)
        WW  = Q(J1,K,L,4)
        UVW = UV(J1,L)
        RI1 = XX(J1,K,L)
        RI2 = XY(J1,K,L)
        RI3 = XZ(J1,K,L)
        RI4 = -RI1*UG(J1,K,L) - RI2*VG(J1,K,L) - RI3*WG(J1,K,L)
        QQ  = QX(J1,L)
        CC  = CX(J1,L)
        QQX = RI4 + QQ
        SP1 = ABS(QQX) + CC
        SM2 = EPS2*SP1
        CHKX= 0.5 + SIGN( 0.5, QQX+CC )
        SPEC= CHKX*( SP1 ) + VN_J(J1,L) + SM2
C
        S1 = S(J1,K,L,1)
        S2 = S(J1,K,L,2)
        S3 = S(J1,K,L,3)
        S4 = S(J1,K,L,4)
        S5 = S(J1,K,L,5)
C
        C1 = RI1*UU + RI2*VV + RI3*WW
        C2 = GM1*UVW
        C3 = RI4 + C1
        C4 = C1*S1
        C5 = UU*S2 + VV*S3 + WW*S4
        C6 = RI1*S2 + RI2*S3 + RI3*S4
        C7 = C6 - C4
        C8 = C2*S1 - GM1*(C5 - S5)
C 
        A(J,1) = SPEC*S1 + CHKX*(RI4*S1 + C6)
        A(J,2) = SPEC*S2 + CHKX*(C3*S2 + RI1*C8 + UU*C7)
        A(J,3) = SPEC*S3 + CHKX*(C3*S3 + RI2*C8 + VV*C7)
        A(J,4) = SPEC*S4 + CHKX*(C3*S4 + RI3*C8 + WW*C7)
        A(J,5) = SPEC*S5 + CHKX*( (RI4 + GAMMA*C1)*S5
     <                   + GE(J1,L)*C7 + C2*C4 - GM1*C1*C5 )
  121 CONTINUE
C
C..Setup C contribution in L-direction
C
      DO 122 J = MS(M),ME(M)
        L  = M-J
        L1 = L-1
C
        UU  = Q(J,K,L1,2)
        VV  = Q(J,K,L1,3)
        WW  = Q(J,K,L1,4)
        UVW = UV(J,L1)
        RI1 = ZX(J,K,L1)
        RI2 = ZY(J,K,L1)
        RI3 = ZZ(J,K,L1)
        RI4 = -RI1*UG(J,K,L1) - RI2*VG(J,K,L1) - RI3*WG(J,K,L1)
        QQ  = QZ(J,L1)
        CC  = CZ(J,L1)
        QQZ = RI4 + QQ
        SP1 = ABS(QQZ) + CC
        SM2 = EPS2*SP1
        CHKZ= 0.5 + SIGN( 0.5, QQZ+CC )
        SPEC= CHKZ*( SP1 ) + VN_L(J,L1) + SM2
C
        S1 = S(J,K,L1,1)
        S2 = S(J,K,L1,2)
        S3 = S(J,K,L1,3)
        S4 = S(J,K,L1,4)
        S5 = S(J,K,L1,5)
C
        C1 = RI1*UU + RI2*VV + RI3*WW
        C2 = GM1*UVW
        C3 = RI4 + C1
        C4 = C1*S1
        C5 = UU*S2 + VV*S3 + WW*S4
        C6 = RI1*S2 + RI2*S3 + RI3*S4
        C7 = C6 - C4
        C8 = C2*S1 - GM1*(C5 - S5)
C 
        C(J,1) = SPEC*S1 + CHKZ*(RI4*S1 + C6)
        C(J,2) = SPEC*S2 + CHKZ*(C3*S2 + RI1*C8 + UU*C7)
        C(J,3) = SPEC*S3 + CHKZ*(C3*S3 + RI2*C8 + VV*C7)
        C(J,4) = SPEC*S4 + CHKZ*(C3*S4 + RI3*C8 + WW*C7)
        C(J,5) = SPEC*S5 + CHKZ*( (RI4 + GAMMA*C1)*S5
     <                   + GE(J,L1)*C7 + C2*C4 - GM1*C1*C5 )
  122 CONTINUE
C
C..Bi-diagonal inversion, include source term for quasi-steady
C     
      IF(IUNST.EQ. FL_HOVER) THEN
Cdir$ ivdep
        DO 123 J = MS(M),ME(M)
          L  = M-J
          SVT = OAT*TSCALE(J,K,L)
          SV = SVT*0.5
          DI = D(J,K,L)
          S(J,K,L,1) = ( S(J,K,L,1) + SV*( A(J,1)+B(J,L,1)+C(J,1) ) )*DI
          R2 = ( S(J,K,L,2) + SV*(A(J,2)+B(J,L,2)+C(J,2) ) ) * DI
          R3 = ( S(J,K,L,3) + SV*(A(J,3)+B(J,L,3)+C(J,3) ) ) * DI
          ABCD = 1.0 + (2.*DI*SV*RF)**2
          S(J,K,L,2) = ( R2 + DI*2.*SV*RF*R3)/ABCD
          S(J,K,L,3) = ( R3 - DI*2.*SV*RF*R2)/ABCD
          S(J,K,L,4) = ( S(J,K,L,4) + SV*( A(J,4)+B(J,L,4)+C(J,4) ) )*DI
          S(J,K,L,5) = ( S(J,K,L,5) + SV*( A(J,5)+B(J,L,5)+C(J,5) ) )*DI
  123   CONTINUE
      ELSE
C
C..Bi-diagonal inversion
C     
Cdir$ ivdep
        DO 124 J = MS(M),ME(M)
          L  = M-J
          SVT = OAT*TSCALE(J,K,L)
          SV = SVT*0.5
          DI = D(J,K,L)
          S(J,K,L,1) = ( S(J,K,L,1) + SV*( A(J,1)+B(J,L,1)+C(J,1) ) )*DI
          S(J,K,L,2) = ( S(J,K,L,2) + SV*( A(J,2)+B(J,L,2)+C(J,2) ) )*DI
          S(J,K,L,3) = ( S(J,K,L,3) + SV*( A(J,3)+B(J,L,3)+C(J,3) ) )*DI
          S(J,K,L,4) = ( S(J,K,L,4) + SV*( A(J,4)+B(J,L,4)+C(J,4) ) )*DI
          S(J,K,L,5) = ( S(J,K,L,5) + SV*( A(J,5)+B(J,L,5)+C(J,5) ) )*DI
 124    CONTINUE
      ENDIF
C
  120 CONTINUE
  100 CONTINUE
C
C..Backward sweep
C..B contribution in K-direction
C
      DO 200 K = KE,KS,-1
C
        KR=K+1
        DO 211 L = LS-1,LE+1
        DO 211 J = JS-1,JE+1
          UU  = Q(J,K,L,2)
          VV  = Q(J,K,L,3)
          WW  = Q(J,K,L,4)
          UVW = 0.5*( UU*UU + VV*VV + WW*WW )
          CJKL = SQRT( GGM1*( Q(J,K,L,5) - UVW ) )
C
          RJ1 = XX(J,K,L)
          RJ2 = XY(J,K,L)
          RJ3 = XZ(J,K,L)
          QQ1 = RJ1*UU + RJ2*VV + RJ3*WW
C
          RR1 = SQRT( RJ1**2 + RJ2**2 + RJ3**2 )
          RL1 = ZX(J,K,L)
          RL2 = ZY(J,K,L)
          RL3 = ZZ(J,K,L)
          QQ3 = RL1*UU + RL2*VV + RL3*WW
C
          RR3 = SQRT( RL1**2 + RL2**2 + RL3**2 )
          
	  RK1 = YX(J,K,L)
          RK2 = YY(J,K,L)
          RK3 = YZ(J,K,L)
          RR2 = SQRT( RK1**2 + RK2**2 + RK3**2 )

          VNU_J = visc_j*2.0*RR1*RR1*(RMUE+TURMU(J,K,L))/(REY*Q(J,K,L,1)
     <		)
          VNU_K = visc_k*2.0*RR2*RR2*(RMUE+TURMU(J,K,L))/(REY*Q(J,K,L,1)
     <		)
          VNU_L = visc_l*2.0*RR3*RR3*(RMUE+TURMU(J,K,L))/(REY*Q(J,K,L,1)
     <		)
	  
          UV(J,L) = UVW
          QX(J,L) = QQ1
          QZ(J,L) = QQ3
          CX(J,L) = CJKL*RR1
          CZ(J,L) = CJKL*RR3
          VN_J(J,L) = VNU_J
          VN_K(J,L) = VNU_K
          VN_L(J,L) = VNU_L
          GE(J,L) = GAMMA*Q(J,K,L,5) - GM1*UVW
  211   CONTINUE
C
        DO 212 L = LS-1,LE+1
        DO 212 J = JS-1,JE+1
          UU  = Q(J,KR,L,2)
          VV  = Q(J,KR,L,3)
          WW  = Q(J,KR,L,4)
          ER  = Q(J,KR,L,5)
          UVW = 0.5*( UU*UU + VV*VV + WW*WW )
          RI1 = YX(J,KR,L)
          RI2 = YY(J,KR,L)
          RI3 = YZ(J,KR,L)
          RI4 = -RI1*UG(J,KR,L) - RI2*VG(J,KR,L) - RI3*WG(J,KR,L)
          RR2 = RI1**2 + RI2**2 + RI3**2
          QQ  = RI1*UU + RI2*VV + RI3*WW
          CC  = SQRT( GGM1*( ER-UVW )*RR2 )
          QQY = RI4 + QQ
          SP1 = ABS(QQY) + CC
          SM2 = EPS2*SP1
          CHKY= 0.5 - SIGN( 0.5, QQY-CC )
          SPEC= CHKY*( -SP1 ) - VN_K(J,L) - SM2
C
          S1 = S(J,KR,L,1)
          S2 = S(J,KR,L,2)
          S3 = S(J,KR,L,3)
          S4 = S(J,KR,L,4)
          S5 = S(J,KR,L,5)
C
          C1 = RI1*UU + RI2*VV + RI3*WW
          C2 = GM1*UVW
          C3 = RI4 + C1
          C4 = C1*S1
          C5 = UU*S2 + VV*S3 + WW*S4
          C6 = RI1*S2 + RI2*S3 + RI3*S4
          C7 = C6 - C4
          C8 = C2*S1 - GM1*(C5 - S5)
C 
          B(J,L,1) = SPEC*S1 + CHKY*(RI4*S1 + C6)
          B(J,L,2) = SPEC*S2 + CHKY*(C3*S2 + RI1*C8 + UU*C7)
          B(J,L,3) = SPEC*S3 + CHKY*(C3*S3 + RI2*C8 + VV*C7)
          B(J,L,4) = SPEC*S4 + CHKY*(C3*S4 + RI3*C8 + WW*C7)
          B(J,L,5) = SPEC*S5 + CHKY*( (RI4 + GAMMA*C1)*S5
     <                       + (GAMMA*ER - C2)*C7 + C2*C4 - GM1*C1*C5 )
  212   CONTINUE
C
C..Loop on hyper-plane
C
      DO 220 M = JE+LE,JS+LS,-1
C
C..Setup A contribution in J-direction
C
      DO 221 J = MS(M),ME(M)
        L  = M-J
        J1 = J+1
C
        UU  = Q(J1,K,L,2)
        VV  = Q(J1,K,L,3)
        WW  = Q(J1,K,L,4)
        UVW = UV(J1,L)
        RI1 = XX(J1,K,L)
        RI2 = XY(J1,K,L)
        RI3 = XZ(J1,K,L)
        RI4 = -RI1*UG(J1,K,L) - RI2*VG(J1,K,L) - RI3*WG(J1,K,L)
        QQ  = QX(J1,L)
        CC  = CX(J1,L)
        QQX = RI4 + QQ
        SP1 = ABS(QQX) + CC
        SM2 = EPS2*SP1
        CHKX= 0.5 - SIGN( 0.5, QQX-CC )
        SPEC= CHKX*( -SP1 ) - VN_J(J1,L) - SM2
C
        S1 = S(J1,K,L,1)
        S2 = S(J1,K,L,2)
        S3 = S(J1,K,L,3)
        S4 = S(J1,K,L,4)
        S5 = S(J1,K,L,5)
C
        C1 = RI1*UU + RI2*VV + RI3*WW
        C2 = GM1*UVW
        C3 = RI4 + C1
        C4 = C1*S1
        C5 = UU*S2 + VV*S3 + WW*S4
        C6 = RI1*S2 + RI2*S3 + RI3*S4
        C7 = C6 - C4
        C8 = C2*S1 - GM1*(C5 - S5)
C 
        A(J,1) = SPEC*S1 + CHKX*(RI4*S1 + C6)
        A(J,2) = SPEC*S2 + CHKX*(C3*S2 + RI1*C8 + UU*C7)
        A(J,3) = SPEC*S3 + CHKX*(C3*S3 + RI2*C8 + VV*C7)
        A(J,4) = SPEC*S4 + CHKX*(C3*S4 + RI3*C8 + WW*C7)
        A(J,5) = SPEC*S5 + CHKX*( (RI4 + GAMMA*C1)*S5
     <                   + GE(J1,L)*C7 + C2*C4 - GM1*C1*C5 )
  221 CONTINUE
C
C..Setup C contribution in L-direction
C
      DO 222 J = MS(M),ME(M)
        L  = M-J
        L1 = L+1
C
        UU  = Q(J,K,L1,2)
        VV  = Q(J,K,L1,3)
        WW  = Q(J,K,L1,4)
        UVW = UV(J,L1)
        RI1 = ZX(J,K,L1)
        RI2 = ZY(J,K,L1)
        RI3 = ZZ(J,K,L1)
        RI4 = -RI1*UG(J,K,L1) - RI2*VG(J,K,L1) - RI3*WG(J,K,L1)
        QQ  = QZ(J,L1)
        CC  = CZ(J,L1)
        QQZ = RI4 + QQ
        SP1 = ABS(QQZ) + CC
        SM2 = EPS2*SP1
        CHKZ= 0.5 - SIGN( 0.5, QQZ-CC )
        SPEC= CHKZ*( -SP1 ) - VN_L(J,L1) - SM2
C
        S1 = S(J,K,L1,1)
        S2 = S(J,K,L1,2)
        S3 = S(J,K,L1,3)
        S4 = S(J,K,L1,4)
        S5 = S(J,K,L1,5)
C
        C1 = RI1*UU + RI2*VV + RI3*WW
        C2 = GM1*UVW
        C3 = RI4 + C1
        C4 = C1*S1
        C5 = UU*S2 + VV*S3 + WW*S4
        C6 = RI1*S2 + RI2*S3 + RI3*S4
        C7 = C6 - C4
        C8 = C2*S1 - GM1*(C5 - S5)
C 
        C(J,1) = SPEC*S1 + CHKZ*(RI4*S1 + C6)
        C(J,2) = SPEC*S2 + CHKZ*(C3*S2 + RI1*C8 + UU*C7)
        C(J,3) = SPEC*S3 + CHKZ*(C3*S3 + RI2*C8 + VV*C7)
        C(J,4) = SPEC*S4 + CHKZ*(C3*S4 + RI3*C8 + WW*C7)
        C(J,5) = SPEC*S5 + CHKZ*( (RI4 + GAMMA*C1)*S5
     <                   + GE(J,L1)*C7 + C2*C4 - GM1*C1*C5 )
  222 CONTINUE
C
C..Bi-diagonal inversion
C     
Cdir$ ivdep
      DO 223 J = MS(M),ME(M)
        L  = M-J
        SVT = OAT*TSCALE(J,K,L)
        SV = SVT*0.5
        DI = D  (J,K,L)
        S(J,K,L,1) = S(J,K,L,1) - SV*( A(J,1)+B(J,L,1)+C(J,1) )*DI
        S(J,K,L,2) = S(J,K,L,2) - SV*( A(J,2)+B(J,L,2)+C(J,2) )*DI
        S(J,K,L,3) = S(J,K,L,3) - SV*( A(J,3)+B(J,L,3)+C(J,3) )*DI
        S(J,K,L,4) = S(J,K,L,4) - SV*( A(J,4)+B(J,L,4)+C(J,4) )*DI
        S(J,K,L,5) = S(J,K,L,5) - SV*( A(J,5)+B(J,L,5)+C(J,5) )*DI
  223 CONTINUE
C
  220 CONTINUE
  200 CONTINUE
C
C..Restore DENSITY , U , V , W in conservative variables             
C                                                                       
      DO 4 L = LS-1, LE+1
      DO 4 K = KS-1, KE+1
      DO 4 J = JS-1, JE+1
        Q(J,K,L,1) = Q(J,K,L,1) / Q(J,K,L,6)
        Q(J,K,L,2) = Q(J,K,L,2)*Q(J,K,L,1)
        Q(J,K,L,3) = Q(J,K,L,3)*Q(J,K,L,1)
        Q(J,K,L,4) = Q(J,K,L,4)*Q(J,K,L,1)
        Q(J,K,L,5) = Q(J,K,L,5)*Q(J,K,L,1)
    4 CONTINUE
C                                                                       
C    
      RETURN
      END

c************************************************************************
      subroutine source(q,s)
c
c  add source term to the rhs
c
c*************************************************************************

      use params_global
      implicit none

      real q(jmax,kmax,lmax,nd),s(jmax,kmax,lmax,nv)
      integer l,k,j

      do 10 l=2,lmax-1
      do 10 k=2,kmax-1
      do 10 j=2,jmax-1
        s(j,k,l,2) = s(j,k,l,2) + rf*q(j,k,l,3)
        s(j,k,l,3) = s(j,k,l,3) - rf*q(j,k,l,2)
  10  continue

      return
      end

c***********************************************************************
      subroutine rhsup1( q,s,xx,xy,xz,yx,yy,yz,zx,zy,zz,
     $                  x,y,z,iblank,ug,vg,wg,kkr,kkp)
c
c  muscl approach:
c  qt = 0                       1st-order  ; |irhsy| = 1
c  qt = 0.25
c    th = 0     upwind-biased   2nd-order  ; |irhsy| = 2
c    th = 1/3   upwind-biased   3rd-order  ; |irhsy| = 3
c
c  limter = 1   differentiable limiter for 3rd-order muscl
c  limter = 2   differentiable limiter for 2nd-order muscl
c  limter = 3   minmod limiter for 2nd and 3rd-order muscl
c
c***********************************************************************

      use params_global
      use weno_routines
      use domain_info, only: this_mesh
c      use work

      implicit none
      
      real q(jmax,kmax,lmax,nd),s(jmax,kmax,lmax,nv)
      integer iblank(jmax,kmax,lmax)
      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)
      real xx(jmax,kmax,lmax),xy(jmax,kmax,lmax),xz(jmax,kmax,lmax),
     <     yx(jmax,kmax,lmax),yy(jmax,kmax,lmax),yz(jmax,kmax,lmax),
     <     zx(jmax,kmax,lmax),zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax),vg(jmax,kmax,lmax),wg(jmax,kmax,lmax)
      integer kkr(kmax),kkp(kmax)
      
c..   local variables

      real, allocatable :: ppp(:,:,:)
      real tj(mdim),xa(mdim),ya(mdim),za(mdim)
      real f(mdim,5),ql(mdim,5),qr(mdim,5)
      real fmin(5),fmax(5),fmin1(5),fmax1(5),fmin2(5),fmax2(5)
      integer iba(mdim)
      integer irhs,limter,ka,kb,ibmin,ibmax,j,k,l
      integer jlo,jup,kout,kin,lout,kr,jj,llo,lup

      real th,qt,eps,tcos,tsin,rhou,rhov,rhoi,fac1,temp,foso
      real f1,q6j,q6j1,xdota,ydota,zdota,q6k,q6k1,q6l,q6l1
      real epsj,epsk,epsl

c************************************************************************
c..   first exectutable statement

      allocate(ppp(jmax,kmax,lmax))

      irhs   = iabs(irhsy)
      limter = 1
      if( irhs .eq. 2 ) limter = 2
      th     = real( irhs - 2 )/real( irhs )
      qt     = 0.25
      if( irhs .eq. 1 ) qt = 0.0
      eps    = 1.e-6
      
      epsj   = (10./real(this_mesh%jm-1))**3
      epsk   = (10./real(this_mesh%km-1))**3
      epsl   = (10./real(this_mesh%lm-1))**3
      
      tcos = cos(blang)
      tsin = sin(blang)
c
      ka = 2 - ksym
      kb = km + ksym

c  *** Compute pressure over interior of mesh (BS):

      do 1 l = 1, lmax
      do 1 k = 1, kmax
      do 1 j = 1, jmax
        ppp(j,k,l) = gm1*( q(j,k,l,5) -0.5*(((q(j,k,l,2)**2)
     $   +q(j,k,l,3)**2) +q(j,k,l,4)**2)/q(j,k,l,1) )
    1 continue

c
c..xi fluxes
c..nonconservative variables
c     
      do 13 l = 2,lm
      do 13 k = ka,kb
c
        do 15 j = 1,jmax
          rhoi   = 1.0/q(j,k,l,1)
          f(j,1) = q(j,k,l,1)*q(j,k,l,6)
          f(j,2) = q(j,k,l,2)*rhoi
          f(j,3) = q(j,k,l,3)*rhoi
          f(j,4) = q(j,k,l,4)*rhoi
          f(j,5) = ppp(j,k,l)*q(j,k,l,6)
          ql(j,1) = 0.25*( ug(j,k+1,l+1)+ug(j,k+1,l-1)
     <                    +ug(j,k-1,l+1)+ug(j,k-1,l-1) )
          ql(j,2) = 0.25*( vg(j,k+1,l+1)+vg(j,k+1,l-1)
     <                    +vg(j,k-1,l+1)+vg(j,k-1,l-1) )
          ql(j,3) = 0.25*( wg(j,k+1,l+1)+wg(j,k+1,l-1)
     <                    +wg(j,k-1,l+1)+wg(j,k-1,l-1) )
          iba(j)=abs(iblank(j,k,l))
   15   continue
c
        f1 = 0.125
        do 14 j = 1,jm
          q6j   = .5/q(j,k,l,6)
          q6j1  = .5/q(j+1,k,l,6)
          xa(j) = q6j*xx(j,k,l) + q6j1*xx(j+1,k,l)
          ya(j) = q6j*xy(j,k,l) + q6j1*xy(j+1,k,l)
          za(j) = q6j*xz(j,k,l) + q6j1*xz(j+1,k,l)
          xdota = 0.5*(ql(j,1)+ql(j+1,1))
          ydota = 0.5*(ql(j,2)+ql(j+1,2))
          zdota = 0.5*(ql(j,3)+ql(j+1,3))
          tj(j) = -xdota*xa(j) - ydota*ya(j) - zdota*za(j)
   14   continue

c..at boundaries

        ibmin = -2
        ibmax = -2
        

        if(ilim.eq.1)then
          call sonica(f,ql,qr,1,jmax,jm,th,qt,eps,fmin,fmax,ibmin,ibmax,
     &          mdim)

        elseif(ilim.eq.0)then
          !print*,'USING MUSCL'
          call muscld_new(f,ql,qr,1,jmax,jm,th,qt,epsj,fmin,fmax,
     &          ibmin,ibmax,mdim,iba)

        elseif(ilim.eq.3)then
          call shmp5(f,ql,qr,1,jmax,jm,th,qt,eps,fmin,fmax,ibmin,ibmax,
     &          mdim)

        elseif(ilim.eq.4)then
          !print*,'USING WENO'
          call weno(f,ql,qr,1,jmax,jm,th,qt,eps,fmin,fmax,ibmin,ibmax,mdim)
        elseif(ilim.eq.6)then
          !print*,'USING WENO (Debo) (Boundary treatment is different'
          call weno5_scalar(f,ql,qr,1,jmax,jm,th,qt,eps,fmin,fmax,ibmin,ibmax,mdim)
        elseif(ilim.eq.8)then
          !print*,'USING CRWENO'
          call crweno5_scalar(f,ql,qr,1,jmax,jm,th,qt,eps,fmin,fmax,ibmin,ibmax,mdim)

        elseif(ilim.eq.2)then
          call quad(f,ql,qr,1,jmax,jm,th,qt,eps,fmin,fmax,fmin1,fmax1,
     &              fmin2,fmax2,ibmin,ibmax,mdim)
        endif
c
c..compute the generalized numerical flux in roe's upwinding
c
        

        if (.not.iprecon) then
          call roeflx( f,ql,qr,xa,ya,za,tj,1,jm )
        else
          call roetrklflx( f,ql,qr,xa,ya,za,tj,1,jm)
        endif


c
        do 12 j = 2,jm
          s(j,k,l,1) = s(j,k,l,1) - ( f(j,1) - f(j-1,1) )
          s(j,k,l,2) = s(j,k,l,2) - ( f(j,2) - f(j-1,2) )
          s(j,k,l,3) = s(j,k,l,3) - ( f(j,3) - f(j-1,3) )
          s(j,k,l,4) = s(j,k,l,4) - ( f(j,4) - f(j-1,4) )
          s(j,k,l,5) = s(j,k,l,5) - ( f(j,5) - f(j-1,5) )
   12   continue
c
   13 continue
c
c..eta fluxes
c
      if(kmax.gt.3) then
c
      do 23 j = 2,jm
      do 23 l = 2,lm
c
        do 25 k = 1,kmax
          rhoi   = 1.0/q(j,k,l,1)
          f(k,1) = q(j,k,l,1)*q(j,k,l,6)
          f(k,2) = q(j,k,l,2)*rhoi
          f(k,3) = q(j,k,l,3)*rhoi
          f(k,4) = q(j,k,l,4)*rhoi
          f(k,5) = ppp(j,k,l)*q(j,k,l,6)
          ql(k,1) = 0.25*( ug(j+1,k,l+1)+ug(j+1,k,l-1)
     <                    +ug(j-1,k,l+1)+ug(j-1,k,l-1) )
          ql(k,2) = 0.25*( vg(j+1,k,l+1)+vg(j+1,k,l-1)
     <                    +vg(j-1,k,l+1)+vg(j-1,k,l-1) )
          ql(k,3) = 0.25*( wg(j+1,k,l+1)+wg(j+1,k,l-1)
     <                    +wg(j-1,k,l+1)+wg(j-1,k,l-1) )
          iba(k)=abs(iblank(j,k,l))
   25   continue
c
        do 24 k = 1,km
          q6k   = .5/q(j,k,l,6)
          q6k1  = .5/q(j,k+1,l,6)
          xa(k) = q6k*yx(j,k,l) + q6k1*yx(j,k+1,l)
          ya(k) = q6k*yy(j,k,l) + q6k1*yy(j,k+1,l)
          za(k) = q6k*yz(j,k,l) + q6k1*yz(j,k+1,l)
          xdota = 0.5*(ql(k,1)+ql(k+1,1))
          ydota = 0.5*(ql(k,2)+ql(k+1,2))
          zdota = 0.5*(ql(k,3)+ql(k+1,3))
          tj(k) = - xdota*xa(k) - ydota*ya(k) - zdota*za(k)
   24   continue

c
        ibmin = -2
        ibmax = -2
c
        if(ilim.eq.1)then
          call sonica(f,ql,qr,1,kmax,km,th,qt,eps,fmin,fmax,ibmin,ibmax,
     &          mdim)
        elseif(ilim.eq.0)then
          call muscld_new(f,ql,qr,1,kmax,km,th,qt,epsk,fmin,fmax,
     &          ibmin,ibmax,mdim,iba)
        elseif(ilim.eq.3)then
          call shmp5(f,ql,qr,1,kmax,km,th,qt,eps,fmin,fmax,ibmin,ibmax,
     &          mdim)
        elseif(ilim.eq.4)then
          call weno(f,ql,qr,1,kmax,km,th,qt,eps,fmin,fmax,ibmin,ibmax,
     &          mdim)
        elseif(ilim.eq.6)then
          call weno5_scalar(f,ql,qr,1,kmax,km,th,qt,eps,fmin,fmax,ibmin,ibmax,
     &          mdim)
        elseif(ilim.eq.8)then
          call crweno5_scalar(f,ql,qr,1,kmax,km,th,qt,eps,fmin,fmax,ibmin,ibmax,
     &          mdim)
        elseif(ilim.eq.2)then
          call quad(f,ql,qr,1,kmax,km,th,qt,eps,fmin,fmax,fmin1,fmax1,
     &              fmin2,fmax2,ibmin,ibmax,mdim)
        endif
c
c..compute the generalized numerical flux in roe's upwinding
c
        if (.not.iprecon) then
           call roeflx( f,ql,qr,xa,ya,za,tj,1,km )
        else
           call roetrklflx( f,ql,qr,xa,ya,za,tj,1,km)
        endif

c
        do 22 k = ka,kb
          kr = kkr(k)
          s(j,k,l,1) = s(j,k,l,1) - ( f(k,1) - f(kr,1) )
          s(j,k,l,2) = s(j,k,l,2) - ( f(k,2) - f(kr,2) )
          s(j,k,l,3) = s(j,k,l,3) - ( f(k,3) - f(kr,3) )
          s(j,k,l,4) = s(j,k,l,4) - ( f(k,4) - f(kr,4) )
          s(j,k,l,5) = s(j,k,l,5) - ( f(k,5) - f(kr,5) )
   22   continue
c
   23 continue
c     
      endif
c
c..zeta fluxes
c
      do 33 j = 2,jm
      do 33 k = ka,kb
c
        do 35 l = 1,lmax
          rhoi   = 1.0/q(j,k,l,1)
          f(l,1) = q(j,k,l,1)*q(j,k,l,6)
          f(l,2) = q(j,k,l,2)*rhoi
          f(l,3) = q(j,k,l,3)*rhoi
          f(l,4) = q(j,k,l,4)*rhoi
          f(l,5) = ppp(j,k,l)*q(j,k,l,6)
          ql(l,1) = 0.25*( ug(j+1,k+1,l)+ug(j+1,k-1,l)
     <                    +ug(j-1,k+1,l)+ug(j-1,k-1,l) )
          ql(l,2) = 0.25*( vg(j+1,k+1,l)+vg(j+1,k-1,l)
     <                    +vg(j-1,k+1,l)+vg(j-1,k-1,l) )
          ql(l,3) = 0.25*( wg(j+1,k+1,l)+wg(j+1,k-1,l)
     <                    +wg(j-1,k+1,l)+wg(j-1,k-1,l) )
          iba(l)=abs(iblank(j,k,l))
   35   continue
c
        do 34 l = 1,lm
          q6l   = .5/q(j,k,l,6)
          q6l1  = .5/q(j,k,l+1,6)
          xa(l) = q6l*zx(j,k,l) + q6l1*zx(j,k,l+1)
          ya(l) = q6l*zy(j,k,l) + q6l1*zy(j,k,l+1)
          za(l) = q6l*zz(j,k,l) + q6l1*zz(j,k,l+1)
          xdota = 0.5*(ql(l,1)+ql(l+1,1))
          ydota = 0.5*(ql(l,2)+ql(l+1,2))
          zdota = 0.5*(ql(l,3)+ql(l+1,3))
          tj(l) = - xdota*xa(l) - ydota*ya(l) - zdota*za(l)
   34   continue
c..at boundaries

        ibmin = -2
        ibmax = -2

        if(ilim.eq.1) then
          call sonica(f,ql,qr,1,lmax,lm,th,qt,eps,fmin,fmax,ibmin,ibmax,
     &          mdim)
        elseif(ilim.eq.0)then
          call muscld_new(f,ql,qr,1,lmax,lm,th,qt,epsl,fmin,fmax,
     &          ibmin,ibmax,mdim,iba)
        elseif(ilim.eq.3)then
          call shmp5(f,ql,qr,1,lmax,lm,th,qt,eps,fmin,fmax,ibmin,ibmax,
     &          mdim)
        elseif(ilim.eq.4)then
          call weno(f,ql,qr,1,lmax,lm,th,qt,eps,fmin,fmax,ibmin,ibmax,
     &          mdim)
        elseif(ilim.eq.6)then
          call weno5_scalar(f,ql,qr,1,lmax,lm,th,qt,eps,fmin,fmax,ibmin,ibmax,
     &          mdim)
        elseif(ilim.eq.8)then
          call crweno5_scalar(f,ql,qr,1,lmax,lm,th,qt,eps,fmin,fmax,ibmin,ibmax,
     &          mdim)
        elseif(ilim.eq.2)then
          call quad(f,ql,qr,1,lmax,lm,th,qt,eps,fmin,fmax,fmin1,fmax1,
     &              fmin2,fmax2,ibmin,ibmax,mdim)
        endif
c
c..compute the generalized numerical flux in roe's upwinding
c
        if (.not.iprecon) then
           call roeflx( f,ql,qr,xa,ya,za,tj,1,lm )
        else
           call roetrklflx( f,ql,qr,xa,ya,za,tj,1,lm)
        endif

c
        do 32 l = 2,lm
          s(j,k,l,1) = s(j,k,l,1) - ( f(l,1) - f(l-1,1) )
          s(j,k,l,2) = s(j,k,l,2) - ( f(l,2) - f(l-1,2) )
          s(j,k,l,3) = s(j,k,l,3) - ( f(l,3) - f(l-1,3) )
          s(j,k,l,4) = s(j,k,l,4) - ( f(l,4) - f(l-1,4) )
          s(j,k,l,5) = s(j,k,l,5) - ( f(l,5) - f(l-1,5) )
   32   continue
c
   33 continue
c
      return
      end

c***********************************************************************
      subroutine rhsup_gcl(q,s,xaj,yaj,zaj,xak,yak,zak,xal,yal,zal,
     <           vnaj,vnak,vnal,kkr,kkp,iblank)


c     based on gcl, calculates exact volumes and metrics
c     the refined mesh calculated in metfv
c
c  muscl approach:
c  qt = 0                       1st-order  ; |irhsy| = 1
c  qt = 0.25
c    th = 0     upwind-biased   2nd-order  ; |irhsy| = 2
c    th = 1/3   upwind-biased   3rd-order  ; |irhsy| = 3
c
c  limter = 1   differentiable limiter for 3rd-order muscl
c  limter = 2   differentiable limiter for 2nd-order muscl
c  limter = 3   minmod limiter for 2nd and 3rd-order muscl
c
c***********************************************************************
      
      use params_global
      use weno_routines

      implicit none
      
      real q(jmax,kmax,lmax,nd),s(jmax,kmax,lmax,nv)
      integer iblank(jmax,kmax,lmax)
      integer kkr(kmax),kkp(kmax)

      real xaj(jmax,kmax,lmax), yaj(jmax,kmax,lmax), zaj(jmax,kmax,lmax)
      real xak(jmax,kmax,lmax), yak(jmax,kmax,lmax), zak(jmax,kmax,lmax)
      real xal(jmax,kmax,lmax), yal(jmax,kmax,lmax), zal(jmax,kmax,lmax)
      real vnaj(jmax,kmax,lmax), vnak(jmax,kmax,lmax), vnal(jmax,kmax,lmax)

c..   local variables ! arrays

      real, allocatable :: ppp(:,:,:)
      real f(mdim,5),ql(mdim,5),qr(mdim,5)
      real tj(mdim),xa(mdim),ya(mdim),za(mdim)
      real fmin(5),fmax(5),fmin1(5),fmax1(5),fmin2(5),fmax2(5)

C..   local variables 

      integer iba(mdim)
      integer irhs,limter,ntac1,ka,kb,ibmin,ibmax
      integer jlo,jup,kout,kin,lout,kp2,km2,k21,lp2,lm2,l21,j21,jp2,jm2
      integer jx,kx,lx,j,k,l,i,kr,jj,llo,lup
      
      
      real th,qt,eps,tcos,tsin,rhou,rhov,rhoi,fac1,temp,foso
      real epsj, epsk, epsl

c************************************************************************

      allocate(ppp(jmax,kmax,lmax))

      irhs   = iabs(irhsy)
      limter = 1
      if( irhs .eq. 2 ) limter = 2
      th     = real( irhs - 2 )/real( irhs )
      qt     = 0.25
      if( irhs .eq. 1 ) qt = 0.0
      eps    = 1.e-6

      epsj   = (10./real(jm))**3
      epsk   = (10./real(km))**3
      epsl   = (10./real(lm))**3

      tcos = cos(blang)
      tsin = sin(blang)

      ntac1= ntac
c
      ka = 2 - ksym
      kb = km + ksym
c
      do 1 l = 1, lmax
      do 1 k = 1, kmax
      do 1 j = 1, jmax
        ppp(j,k,l) = gm1*( q(j,k,l,5) -0.5*(((q(j,k,l,2)**2)
     $   +q(j,k,l,3)**2) +q(j,k,l,4)**2)/q(j,k,l,1) )
    1 continue

c
c..xi fluxes
c..nonconservative variables
c     

      do 13 l = 2,lm
      do 13 k = ka,kb
c
        do 15 j = 1,jmax
          rhoi   = 1.0/q(j,k,l,1)
          f(j,1) = q(j,k,l,1)*q(j,k,l,6)
          f(j,2) = q(j,k,l,2)*rhoi
          f(j,3) = q(j,k,l,3)*rhoi
          f(j,4) = q(j,k,l,4)*rhoi
          f(j,5) = ppp(j,k,l)*q(j,k,l,6)
          iba(j) =abs(iblank(j,k,l))
 15    continue

c..at boundaries

        ibmin = -2
        ibmax = -2

        if(half.eq.1) then
          ibmax = 2
          fmax(1) = f(jmax-3,1)
          fmax(2) = f(jmax-3,2)
          fmax(3) = f(jmax-3,3)
          fmax(4) = -f(jmax-3,4)
          fmax(5) = f(jmax-3,5)
        endif

        if(iwake.eq.2) then
          ibmin = 3
          jlo = jle+1-l
c
          rhou = q(jlo,k,lmax-2,2)*tcos + q(jlo,k,lmax-2,3)*tsin
          rhov = -q(jlo,k,lmax-2,2)*tsin + q(jlo,k,lmax-2,3)*tcos
c
          rhoi   = 1.0/q(jlo,k,lmax-2,1)
          fmin(1) = q(jlo,k,lmax-2,1)*q(jlo,k,lmax-2,6)
          fmin(2) = rhou*rhoi
          fmin(3) = rhov*rhoi
          fmin(4) = q(jlo,k,lmax-2,4)*rhoi
          fmin(5) = ppp(jlo,k,lmax-2)*q(jlo,k,lmax-2,6)
ctl
          rhou = q(jlo,k,lmax-3,2)*tcos + q(jlo,k,lmax-3,3)*tsin
          rhov = -q(jlo,k,lmax-3,2)*tsin + q(jlo,k,lmax-3,3)*tcos
c
          rhoi   = 1.0/q(jlo,k,lmax-3,1)
          fmin1(1) = q(jlo,k,lmax-3,1)*q(jlo,k,lmax-3,6)
          fmin1(2) = rhou*rhoi
          fmin1(3) = rhov*rhoi
          fmin1(4) = q(jlo,k,lmax-3,4)*rhoi
          fmin1(5) = ppp(jlo,k,lmax-3)*q(jlo,k,lmax-3,6)
c
          rhou = q(jlo,k,lmax-4,2)*tcos + q(jlo,k,lmax-4,3)*tsin
          rhov = -q(jlo,k,lmax-4,2)*tsin + q(jlo,k,lmax-4,3)*tcos
c
          rhoi   = 1.0/q(jlo,k,lmax-4,1)
          fmin2(1) = q(jlo,k,lmax-4,1)*q(jlo,k,lmax-4,6)
          fmin2(2) = rhou*rhoi
          fmin2(3) = rhov*rhoi
          fmin2(4) = q(jlo,k,lmax-4,4)*rhoi
          fmin2(5) = ppp(jlo,k,lmax-4)*q(jlo,k,lmax-4,6)
        endif
        if(iwake.eq.2 .and. half.ne.1) then
          ibmax = 3
          jup = jle-1+l
c
          rhou = q(jup,k,lmax-2,2)*tcos + q(jup,k,lmax-2,3)*tsin
          rhov = -q(jup,k,lmax-2,2)*tsin + q(jup,k,lmax-2,3)*tcos
c
          rhoi   = 1.0/q(jup,k,lmax-2,1)
          fmax(1) = q(jup,k,lmax-2,1)*q(jup,k,lmax-2,6)
          fmax(2) = rhou*rhoi
          fmax(3) = rhov*rhoi
          fmax(4) = q(jup,k,lmax-2,4)*rhoi
          fmax(5) = ppp(jup,k,lmax-2)*q(jup,k,lmax-2,6)
ctl
          rhou = q(jup,k,lmax-3,2)*tcos + q(jup,k,lmax-3,3)*tsin
          rhov = -q(jup,k,lmax-3,2)*tsin + q(jup,k,lmax-3,3)*tcos
c
          rhoi   = 1.0/q(jup,k,lmax-3,1)
          fmax1(1) = q(jup,k,lmax-3,1)*q(jup,k,lmax-3,6)
          fmax1(2) = rhou*rhoi
          fmax1(3) = rhov*rhoi
          fmax1(4) = q(jup,k,lmax-3,4)*rhoi
          fmax1(5) = ppp(jup,k,lmax-3)*q(jup,k,lmax-3,6)
c
          rhou = q(jup,k,lmax-4,2)*tcos + q(jup,k,lmax-4,3)*tsin
          rhov = -q(jup,k,lmax-4,2)*tsin + q(jup,k,lmax-4,3)*tcos
c
          rhoi   = 1.0/q(jup,k,lmax-4,1)
          fmax2(1) = q(jup,k,lmax-4,1)*q(jup,k,lmax-4,6)
          fmax2(2) = rhou*rhoi
          fmax2(3) = rhov*rhoi
          fmax2(4) = q(jup,k,lmax-4,4)*rhoi
          fmax2(5) = ppp(jup,k,lmax-4)*q(jup,k,lmax-4,6)
        endif

c..set cell face area vectors and mean face velocity

        do j=1,jm
            xa(j) = xaj(j,k,l)
            ya(j) = yaj(j,k,l)
            za(j) = zaj(j,k,l)
            tj(j) = vnaj(j,k,l)
        end do

c..limit

        kout = ktip+3.*(kmax-ktip)/4
        kin  = ktip/2
        lout = lmax/2

        if(ilim.eq.1)then
          call sonica(f,ql,qr,1,jmax,jm,th,qt,eps,fmin,fmax,ibmin,ibmax,
     &          mdim)

        elseif(ilim.eq.0)then
          call muscld_new(f,ql,qr,1,jmax,jm,th,qt,epsj,fmin,fmax,
     &         ibmin,ibmax,mdim,iba)

        elseif(ilim.eq.3)then
          call shmp5(f,ql,qr,1,jmax,jm,th,qt,eps,fmin,fmax,ibmin,ibmax,
     &          mdim)

        elseif(ilim.eq.4)then
          call weno(f,ql,qr,1,jmax,jm,th,qt,eps,fmin,fmax,ibmin,ibmax,
     &          mdim)
        elseif(ilim.eq.6)then
          call weno5_scalar(f,ql,qr,1,jmax,jm,th,qt,eps,fmin,fmax,ibmin,ibmax,
     &          mdim)
        elseif(ilim.eq.8)then
          call crweno5_scalar(f,ql,qr,1,jmax,jm,th,qt,eps,fmin,fmax,ibmin,ibmax,
     &          mdim)

        elseif(ilim.eq.2)then
          call quad(f,ql,qr,1,jmax,jm,th,qt,eps,fmin,fmax,fmin1,fmax1,
     &              fmin2,fmax2,ibmin,ibmax,mdim)
        endif
        
c
c..compute the generalized numerical flux in roe's upwinding
c
        call roeflx( f,ql,qr,xa,ya,za,tj,1,jm )

c
        do 12 j = 2,jm
          s(j,k,l,1) = s(j,k,l,1) - ( f(j,1) - f(j-1,1) )
          s(j,k,l,2) = s(j,k,l,2) - ( f(j,2) - f(j-1,2) )
          s(j,k,l,3) = s(j,k,l,3) - ( f(j,3) - f(j-1,3) )
          s(j,k,l,4) = s(j,k,l,4) - ( f(j,4) - f(j-1,4) )
          s(j,k,l,5) = s(j,k,l,5) - ( f(j,5) - f(j-1,5) )
   12   continue
c
   13 continue

c
c..eta fluxes
c
      if(kmax.gt.3) then
c
      do 23 j = 2,jm
      do 23 l = 2,lm
c
        do 25 k = 1,kmax
          rhoi   = 1.0/q(j,k,l,1)
          f(k,1) = q(j,k,l,1)*q(j,k,l,6)
          f(k,2) = q(j,k,l,2)*rhoi
          f(k,3) = q(j,k,l,3)*rhoi
          f(k,4) = q(j,k,l,4)*rhoi
          f(k,5) = ppp(j,k,l)*q(j,k,l,6)
          iba(k) =abs(iblank(j,k,l))
   25   continue

        ibmin = -2
        ibmax = -2

        if(invisc) then
          ibmin = 2
          foso = 0.99
          fmin(1) = (1.+foso)*f(1,1)-foso*f(2,1)
          fmin(2) = (1.+foso)*f(1,2)-foso*f(2,2)
          fmin(3) = (1.+foso)*f(1,3)-foso*f(2,3)
          fmin(4) = (1.+foso)*f(1,4)-foso*f(2,4)
          fmin(5) = (1.+foso)*f(1,5)-foso*f(2,5)
        endif

c..set cell face area vectors and mean face velocity
        do k=1,km
            xa(k) = xak(j,k,l)
            ya(k) = yak(j,k,l)
            za(k) = zak(j,k,l)
            tj(k) = vnak(j,k,l)
        enddo

c.. limit      

        if(ilim.eq.1)then
          call sonica(f,ql,qr,1,kmax,km,th,qt,eps,fmin,fmax,ibmin,ibmax,
     &          mdim)
        elseif(ilim.eq.0)then
          call muscld_new(f,ql,qr,1,kmax,km,th,qt,epsk,fmin,fmax,
     &         ibmin,ibmax,mdim,iba)
        elseif(ilim.eq.3)then
          call shmp5(f,ql,qr,1,kmax,km,th,qt,eps,fmin,fmax,ibmin,ibmax,
     &          mdim)
        elseif(ilim.eq.4)then
          call weno(f,ql,qr,1,kmax,km,th,qt,eps,fmin,fmax,ibmin,ibmax,
     &          mdim)
        elseif(ilim.eq.6)then
          call weno5_scalar(f,ql,qr,1,kmax,km,th,qt,eps,fmin,fmax,ibmin,ibmax,
     &          mdim)
        elseif(ilim.eq.8)then
          call crweno5_scalar(f,ql,qr,1,kmax,km,th,qt,eps,fmin,fmax,ibmin,ibmax,
     &          mdim)
        elseif(ilim.eq.2)then
          call quad(f,ql,qr,1,kmax,km,th,qt,eps,fmin,fmax,fmin1,fmax1,
     &              fmin2,fmax2,ibmin,ibmax,mdim)
        endif
     
c     
c..   compute the generalized numerical flux in roe's upwinding
c     
        call roeflx( f,ql,qr,xa,ya,za,tj,1,km )
c     
        do 22 k = ka,kb
          kr = kkr(k)
          s(j,k,l,1) = s(j,k,l,1) - ( f(k,1) - f(kr,1) )
          s(j,k,l,2) = s(j,k,l,2) - ( f(k,2) - f(kr,2) )
          s(j,k,l,3) = s(j,k,l,3) - ( f(k,3) - f(kr,3) )
          s(j,k,l,4) = s(j,k,l,4) - ( f(k,4) - f(kr,4) )
          s(j,k,l,5) = s(j,k,l,5) - ( f(k,5) - f(kr,5) )
   22   continue
c
   23 continue
c
      endif
c
c..zeta fluxes
c

      do 33 j = 2,jm
      do 33 k = ka,kb
c
        do 35 l = 1,lmax
          rhoi   = 1.0/q(j,k,l,1)
          f(l,1) = q(j,k,l,1)*q(j,k,l,6)
          f(l,2) = q(j,k,l,2)*rhoi
          f(l,3) = q(j,k,l,3)*rhoi
          f(l,4) = q(j,k,l,4)*rhoi
          f(l,5) = ppp(j,k,l)*q(j,k,l,6)
          iba(l) =abs(iblank(j,k,l))
 35    continue

c..at boundaries

        ibmin = -2
        ibmax = -2

c..use symmetry if only half-plane, otherwise use other side
        if(half.eq.1) then
          if(k.gt.ktip) then
            ibmin = 2
            fmin(1) = f(2,1)
            fmin(2) = f(2,2)
            fmin(3) = f(2,3)
            fmin(4) = -f(2,4)
            fmin(5) = f(2,5)
          endif
          if(k.le.ktip .and. j.le.jtail1) then
            ibmin = 2
            fmin(1) = f(2,1)
            fmin(2) = f(2,2)
            fmin(3) = f(2,3)
            fmin(4) = -f(2,4)
            fmin(5) = f(2,5)
          endif
        else
          if(k.gt.ktip) then
            ibmin = 2
            jj = jmax-j+1
            rhoi   = 1.0/q(jj,k,2,1)
            fmin(1) = q(jj,k,2,1)*q(jj,k,2,6)
            fmin(2) = q(jj,k,2,2)*rhoi
            fmin(3) = q(jj,k,2,3)*rhoi
            fmin(4) = q(jj,k,2,4)*rhoi
            fmin(5) = ppp(jj,k,2)*q(jj,k,2,6) 
          endif
          if(k.le.ktip .and. (j.le.jtail1 .or. j.ge.jtail2) ) then 
            ibmin = 2
            jj = jmax-j+1
            rhoi   = 1.0/q(jj,k,2,1)
            fmin(1) = q(jj,k,2,1)*q(jj,k,2,6)
            fmin(2) = q(jj,k,2,2)*rhoi
            fmin(3) = q(jj,k,2,3)*rhoi
            fmin(4) = q(jj,k,2,4)*rhoi
            fmin(5) = ppp(jj,k,2)*q(jj,k,2,6) 
          endif
        endif
c..try to raise accuracy at wall if inviscid
c       if(invisc .and. k.le.ktip .and.
c    <          (j.gt.jtail1 .and. j.lt.jtail2) ) then
c         ibmin = 2
c         foso = 0.99
c         fmin(1) = (1.+foso)*f(1,1)-foso*f(2,1)
c         fmin(2) = (1.+foso)*f(1,2)-foso*f(2,2)
c         fmin(3) = (1.+foso)*f(1,3)-foso*f(2,3)
c         fmin(4) = (1.+foso)*f(1,4)-foso*f(2,4)
c         fmin(5) = (1.+foso)*f(1,5)-foso*f(2,5)
c       endif
c..limit
        if(iwake.eq.2 .and. (j.ge.(jle+1-lmax).and.j.le.jle)) then
          ibmax = 3
          llo = jle+1-j
c
          rhou = q(3,k,llo,2)*tcos - q(3,k,llo,3)*tsin
          rhov = q(3,k,llo,2)*tsin + q(3,k,llo,3)*tcos
c
          rhoi   = 1.0/q(3,k,llo,1)
          fmax(1) = q(3,k,llo,1)*q(3,k,llo,6)
          fmax(2) = rhou*rhoi
          fmax(3) = rhov*rhoi
          fmax(4) = q(3,k,llo,4)*rhoi
          fmax(5) = ppp(3,k,llo)*q(3,k,llo,6)
ctl
          rhou = q(4,k,llo,2)*tcos - q(4,k,llo,3)*tsin
          rhov = q(4,k,llo,2)*tsin + q(4,k,llo,3)*tcos
c
          rhoi   = 1.0/q(4,k,llo,1)
          fmax1(1) = q(4,k,llo,1)*q(4,k,llo,6)
          fmax1(2) = rhou*rhoi
          fmax1(3) = rhov*rhoi
          fmax1(4) = q(4,k,llo,4)*rhoi
          fmax1(5) = ppp(4,k,llo)*q(4,k,llo,6)
c
          rhou = q(5,k,llo,2)*tcos - q(5,k,llo,3)*tsin
          rhov = q(5,k,llo,2)*tsin + q(5,k,llo,3)*tcos
c
          rhoi   = 1.0/q(5,k,llo,1)
          fmax2(1) = q(5,k,llo,1)*q(5,k,llo,6)
          fmax2(2) = rhou*rhoi
          fmax2(3) = rhov*rhoi
          fmax2(4) = q(5,k,llo,4)*rhoi
          fmax2(5) = ppp(5,k,llo)*q(5,k,llo,6)
        endif
        if(iwake.eq.2 .and. (j.le.(jle-1+lmax).and.j.ge.jle) .and.
     >     half.ne.1) then
          ibmax = 3
          lup = -jle+1+j
c
          rhou = q(jmax-2,k,lup,2)*tcos - q(jmax-2,k,lup,3)*tsin
          rhov = q(jmax-2,k,lup,2)*tsin + q(jmax-2,k,lup,3)*tcos
c
          rhoi   = 1.0/q(jmax-2,k,lup,1)
          fmax(1) = q(jmax-2,k,lup,1)*q(jmax-2,k,lup,6)
          fmax(2) = rhou*rhoi
          fmax(3) = rhov*rhoi
          fmax(4) = q(jmax-2,k,lup,4)*rhoi
          fmax(5) = ppp(jmax-2,k,lup)*q(jmax-2,k,lup,6)
ctl
          rhou = q(jmax-3,k,lup,2)*tcos - q(jmax-3,k,lup,3)*tsin
          rhov = q(jmax-3,k,lup,2)*tsin + q(jmax-3,k,lup,3)*tcos
c
          rhoi   = 1.0/q(jmax-3,k,lup,1)
          fmax1(1) = q(jmax-3,k,lup,1)*q(jmax-3,k,lup,6)
          fmax1(2) = rhou*rhoi
          fmax1(3) = rhov*rhoi
          fmax1(4) = q(jmax-3,k,lup,4)*rhoi
          fmax1(5) = ppp(jmax-3,k,lup)*q(jmax-3,k,lup,6)
c
          rhou = q(jmax-4,k,lup,2)*tcos - q(jmax-4,k,lup,3)*tsin
          rhov = q(jmax-4,k,lup,2)*tsin + q(jmax-4,k,lup,3)*tcos
c
          rhoi   = 1.0/q(jmax-4,k,lup,1)
          fmax2(1) = q(jmax-4,k,lup,1)*q(jmax-4,k,lup,6)
          fmax2(2) = rhou*rhoi
          fmax2(3) = rhov*rhoi
          fmax2(4) = q(jmax-4,k,lup,4)*rhoi
          fmax2(5) = ppp(jmax-4,k,lup)*q(jmax-4,k,lup,6)
        endif

        do l=1,lm
            xa(l) = xal(j,k,l)
            ya(l) = yal(j,k,l)
            za(l) = zal(j,k,l)
            tj(l) = vnal(j,k,l)
        enddo
c
        if(ilim.eq.1) then
          call sonica(f,ql,qr,1,lmax,lm,th,qt,eps,fmin,fmax,ibmin,ibmax,
     &          mdim)
        elseif(ilim.eq.0)then
          call muscld_new(f,ql,qr,1,lmax,lm,th,qt,epsl,fmin,fmax,
     &         ibmin,ibmax,mdim,iba)
        elseif(ilim.eq.3)then
          call shmp5(f,ql,qr,1,lmax,lm,th,qt,eps,fmin,fmax,ibmin,ibmax,
     &          mdim)
        elseif(ilim.eq.4)then
          call weno(f,ql,qr,1,lmax,lm,th,qt,eps,fmin,fmax,ibmin,ibmax,
     &          mdim)
        elseif(ilim.eq.6)then
          call weno5_scalar(f,ql,qr,1,lmax,lm,th,qt,eps,fmin,fmax,ibmin,ibmax,
     &          mdim)
        elseif(ilim.eq.8)then
          call crweno5_scalar(f,ql,qr,1,lmax,lm,th,qt,eps,fmin,fmax,ibmin,ibmax,
     &          mdim)
        elseif(ilim.eq.2)then
          call quad(f,ql,qr,1,lmax,lm,th,qt,eps,fmin,fmax,fmin1,fmax1,
     &              fmin2,fmax2,ibmin,ibmax,mdim)
        endif
c
c..compute the generalized numerical flux in roe's upwinding
c
        call roeflx( f,ql,qr,xa,ya,za,tj,1,lm )
c
        do 32 l = 2,lm
          s(j,k,l,1) = s(j,k,l,1) - ( f(l,1) - f(l-1,1) )
          s(j,k,l,2) = s(j,k,l,2) - ( f(l,2) - f(l-1,2) )
          s(j,k,l,3) = s(j,k,l,3) - ( f(l,3) - f(l-1,3) )
          s(j,k,l,4) = s(j,k,l,4) - ( f(l,4) - f(l-1,4) )
          s(j,k,l,5) = s(j,k,l,5) - ( f(l,5) - f(l-1,5) )
   32   continue
c
   33 continue
c
      deallocate(ppp)

      return
      end


c***********************************************************************
      subroutine roeflx( f,ql,qr,xa,ya,za,tj,is,ie )
c
c  compute the generalized numerical flux in roe's upwinding
c  by s.o.                          
c  mod by jdb to incorporate smoother entropy check
c
c***********************************************************************

      use params_global
      
      implicit none

      integer is, ie

      real f(mdim,5)
      real tj(mdim),xa(mdim),ya(mdim),za(mdim)
      real ql(mdim,5),qr(mdim,5)

c..   local variables

      integer i,i1

      real eps
      real rlft,ulft,vlft,wlft,plft
      real rlfti,rulft,rvlft,rwlft,uvwl
      real elft,hlft,clft
      real rrht,urht,vrht,wrht,prht
      real rrhti,rurht,rvrht,rwrht,uvwr,erht,hrht,crht
      real rat,rati,rav,uav,vav,wav,hav,uvw,cav
      real aq1,aq2,aq3,aq4,aq5,ri1,ri2,ri3,ri4
      real rr2,rr,r0,r1,r2,r3,r4
      real uu,c2,c2i,auu,aupc,aumc,uulft,uurht
      real upclft,upcrht,umclft,umcrht
      real dauu,dauus,daupc,daupcs,daumc,daumcs
      real rcav,aquu,c2ih,ruuav,b1,b2,b3,b4,b5,b6,b7,b8
      real aj,plar,eplft,eprht,fssub

      eps = 1.e-6
      do 11 i = is,ie
c
       i1    = i + 1
       rlft = ql(i,1)
       ulft = ql(i,2)
       vlft = ql(i,3)
       wlft = ql(i,4)
       plft = ql(i,5)
       rlfti = 1.0/rlft
       rulft = rlft*ulft
       rvlft = rlft*vlft
       rwlft = rlft*wlft
       uvwl = 0.5*( ulft*ulft + vlft*vlft + wlft*wlft )
       elft = plft/gm1 + rlft*uvwl
       hlft = ( elft + plft )*rlfti
       clft = sqrt( gm1*( hlft - uvwl ) )
c
       rrht = qr(i1,1)
       urht = qr(i1,2)
       vrht = qr(i1,3)
       wrht = qr(i1,4)
       prht = qr(i1,5)
       rrhti = 1.0/rrht
       rurht = rrht*urht
       rvrht = rrht*vrht
       rwrht = rrht*wrht
       uvwr = 0.5*( urht*urht + vrht*vrht + wrht*wrht )
       erht = prht/gm1 + rrht*uvwr
       hrht = ( erht + prht )*rrhti
       crht = sqrt( gm1*( hrht - uvwr ) )
c
       rat  = sqrt( rrht*rlfti )
       rati = 1.0/( rat + 1. )
       rav  =   rat*rlft
       uav  = ( rat*urht + ulft )*rati
       vav  = ( rat*vrht + vlft )*rati
       wav  = ( rat*wrht + wlft )*rati
       hav  = ( rat*hrht + hlft )*rati
       uvw  = 0.5*( uav*uav + vav*vav + wav*wav )
       cav  = sqrt( gm1*( hav - uvw ) )
c
       aq1  = rrht - rlft
       aq2  = urht - ulft
       aq3  = vrht - vlft
       aq4  = wrht - wlft
       aq5  = prht - plft
c
       ri1 = xa(i)
       ri2 = ya(i)
       ri3 = za(i)
       ri4 = tj(i)
       rr2 = ri1*ri1 + ri2*ri2 + ri3*ri3
       rr  = sqrt( rr2 )
       r0  = 1.0 / rr
       r1  = ri1*r0
       r2  = ri2*r0
       r3  = ri3*r0
       r4  = ri4*r0
c
       uu  = r1*uav + r2*vav + r3*wav + r4
       c2  = cav*cav
       c2i = 1.0/c2
c
       auu   = abs( uu    )
       aupc  = abs( uu+cav )
       aumc  = abs( uu-cav )
c     
       uulft = r1*ulft + r2*vlft + r3*wlft + r4
       uurht = r1*urht + r2*vrht + r3*wrht + r4
       upclft= uulft + clft
       upcrht= uurht + crht
       umclft= uulft - clft
       umcrht= uurht - crht
c
       dauu = 4.*(uurht-uulft)+eps
       dauus = amax1(dauu,0.0)
ccray       auu = cvmgt(auu**2/dauu+0.25*dauu,auu,auu.le.0.5*dauus)
       if (auu.le.0.5*dauus) auu = auu**2/dauu+0.25*dauu
       daupc = 4.*(upcrht-upclft)+eps
       daupcs = amax1(daupc,0.0)
ccray       aupc = cvmgt(aupc**2/daupc+0.25*daupc,aupc,aupc.le.0.5*daupcs)
       if (aupc.le.0.5*daupcs) aupc = aupc**2/daupc+0.25*daupc
       daumc = 4.*(umcrht-umclft)+eps
       daumcs = amax1(daumc,0.0)
ccray       aumc = cvmgt(aumc**2/daumc+0.25*daumc,aumc,aumc.le.0.5*daumcs)
       if (aumc.le.0.5*daumcs) aumc = aumc**2/daumc+0.25*daumc
c     
       rcav = rav*cav
       aquu = uurht - uulft
       c2ih = 0.5*c2i
       ruuav= auu*rav
       b1   = auu*( aq1 - c2i*aq5 )
       b2   = c2ih*aupc*( aq5 + rcav*aquu )
       b3   = c2ih*aumc*( aq5 - rcav*aquu )
       b4   = b1 + b2 + b3
       b5   = cav*( b2 - b3 )
       b6   = ruuav*( aq2 - r1*aquu )
       b7   = ruuav*( aq3 - r2*aquu )
       b8   = ruuav*( aq4 - r3*aquu )
c
       aq1 = b4
       aq2 = uav*b4 + r1*b5 + b6
       aq3 = vav*b4 + r2*b5 + b7
       aq4 = wav*b4 + r3*b5 + b8
       aq5 = hav*b4 + ( uu-r4 )*b5 + uav*b6 + vav*b7 + wav*b8
     <       - c2*b1/gm1
c
       aj    = 0.5*rr
       plar  = plft + prht
       eplft = elft + plft
       eprht = erht + prht

       fssub = rr*r4
       if(usegcl) fssub = 0.0
       fssub=0.0

       f(i,1) = aj*(  rlft*uulft +  rrht*uurht           - aq1 )
     <                                          -fssub*rinf
       f(i,2) = aj*( rulft*uulft + rurht*uurht + r1*plar - aq2 )
     <                                          -fssub*rinf*uinf
       f(i,3) = aj*( rvlft*uulft + rvrht*uurht + r2*plar - aq3 )
     <                                          -fssub*rinf*vinf
       f(i,4) = aj*( rwlft*uulft + rwrht*uurht + r3*plar - aq4 )
     <                                          -fssub*rinf*winf
       f(i,5) = aj*( eplft*uulft + eprht*uurht - r4*plar - aq5 )
     <                                          -fssub*einf
   11 continue
c
      return
      end


c***********************************************************************
      subroutine roetrklflx( f,ql,qr,xa,ya,za,tj,is,ie)
c
c  compute the generalized numerical flux in roe's upwinding
c  by s.o.                          
c  mod by jdb to incorporate smoother entropy check
c
c***********************************************************************

      use params_global
      
      implicit none

      integer is, ie

      real f(mdim,5)
      real tj(mdim),xa(mdim),ya(mdim),za(mdim)
      real ql(mdim,5),qr(mdim,5)

c..   local variables

      integer i,i1

      real eps
      real rlft,ulft,vlft,wlft,plft
      real rlfti,rulft,rvlft,rwlft,uvwl
      real elft,hlft,clft
      real rrht,urht,vrht,wrht,prht
      real rrhti,rurht,rvrht,rwrht,uvwr,erht,hrht,crht
      real rat,rati,rav,uav,vav,wav,hav,uvw,cav
      real aq1,aq2,aq3,aq4,aq5,ri1,ri2,ri3,ri4
      real rr2,rr,r0,r1,r2,r3,r4
      real uu,c2,c2i,auu,aupc,aumc,uulft,uurht
      real upclft,upcrht,umclft,umcrht
      real dauu,dauus,daupc,daupcs,daumc,daumcs
      real aquu,ruuav,b1,b2,b3,b4,b5,b6,b7,b8
      real aj,plar,eplft,eprht,fssub
      real R,S,X,bSq

      eps = 1.e-6
      do 11 i = is,ie
c
       i1    = i + 1
       rlft = ql(i,1)
       ulft = ql(i,2)
       vlft = ql(i,3)
       wlft = ql(i,4)
       plft = ql(i,5)
       rlfti = 1.0/rlft
       rulft = rlft*ulft
       rvlft = rlft*vlft
       rwlft = rlft*wlft
       uvwl = 0.5*( ulft*ulft + vlft*vlft + wlft*wlft )
       elft = plft/gm1 + rlft*uvwl
       hlft = ( elft + plft )*rlfti
       clft = sqrt( gm1*( hlft - uvwl ) )
c
       rrht = qr(i1,1)
       urht = qr(i1,2)
       vrht = qr(i1,3)
       wrht = qr(i1,4)
       prht = qr(i1,5)
       rrhti = 1.0/rrht
       rurht = rrht*urht
       rvrht = rrht*vrht
       rwrht = rrht*wrht
       uvwr = 0.5*( urht*urht + vrht*vrht + wrht*wrht )
       erht = prht/gm1 + rrht*uvwr
       hrht = ( erht + prht )*rrhti
       crht = sqrt( gm1*( hrht - uvwr ) )
c
       rat  = sqrt( rrht*rlfti )
       rati = 1.0/( rat + 1. )
       rav  =   rat*rlft
       uav  = ( rat*urht + ulft )*rati
       vav  = ( rat*vrht + vlft )*rati
       wav  = ( rat*wrht + wlft )*rati
       hav  = ( rat*hrht + hlft )*rati
       uvw  = 0.5*( uav*uav + vav*vav + wav*wav )
       cav  = sqrt( gm1*( hav - uvw ) )
c
       aq1  = rrht - rlft
       aq2  = urht - ulft
       aq3  = vrht - vlft
       aq4  = wrht - wlft
       aq5  = prht - plft
c
       ri1 = xa(i)
       ri2 = ya(i)
       ri3 = za(i)
       ri4 = tj(i)
       rr2 = ri1*ri1 + ri2*ri2 + ri3*ri3
       rr  = sqrt( rr2 )
       r0  = 1.0 / rr
       r1  = ri1*r0
       r2  = ri2*r0
       r3  = ri3*r0
       r4  = ri4*r0
c
       uu  = r1*uav + r2*vav + r3*wav + r4
       c2  = cav*cav
       c2i = 1.0/c2
c
       bSq = Mp**2

       X = sqrt( (1.-bSq)*uu*(1.-bSq)*uu+4.*bSq*c2 )
       auu   = abs( uu    )
       aupc  = 0.5*abs( (1.+bSq)*uu + X )
       aumc  = 0.5*abs( (1.+bSq)*uu - X )
       R     = 0.5*((1.-bSq)*uu + X)
       S     = 0.5*((1.-bSq)*uu - X)
c
       uulft = r1*ulft + r2*vlft + r3*wlft + r4
       uurht = r1*urht + r2*vrht + r3*wrht + r4
       X = sqrt( (1.-bSq)*uulft*((1.-bSq)*uulft)+4.*bSq*clft*clft )
       upclft= 0.5*( (1.+bSq)*uulft + X )
       umclft= 0.5*( (1.+bSq)*uulft - X )
       X = sqrt( (1.-bSq)*uurht*((1.-bSq)*uurht)+4.*bSq*crht*crht )
       upcrht= 0.5*( (1.+bSq)*uurht + X )
       umcrht= 0.5*( (1.+bSq)*uurht - X )
c
       dauu = 4.*(uurht-uulft)+eps
       dauus = amax1(dauu,0.0)
ccray       auu = cvmgt(auu**2/dauu+0.25*dauu,auu,auu.le.0.5*dauus)
       if (auu.le.0.5*dauus) auu = auu**2/dauu+0.25*dauu
       daupc = 4.*(upcrht-upclft)+eps
       daupcs = amax1(daupc,0.0)
ccray       aupc = cvmgt(aupc**2/daupc+0.25*daupc,aupc,aupc.le.0.5*daupcs)
       if (aupc.le.0.5*daupcs) aupc = aupc**2/daupc+0.25*daupc
       daumc = 4.*(umcrht-umclft)+eps
       daumcs = amax1(daumc,0.0)
ccray       aumc = cvmgt(aumc**2/daumc+0.25*daumc,aumc,aumc.le.0.5*daumcs)
       if (aumc.le.0.5*daumcs) aumc = aumc**2/daumc+0.25*daumc
c
       aquu = uurht - uulft
       ruuav= auu*rav
       X    = sqrt( (1.-bSq)*uu*(1.-bSq)*uu+4.*bSq*c2 )
       b1   = auu*( aq1 - c2i*aq5 )
       b2   = aupc*( aq5/R + rav*aquu )/X
       b3   = aumc*(-aq5/S - rav*aquu )/X
       b4   = b1 + b2 + b3
       b5   = ( R*b2 + S*b3 )
       b6   = ruuav*( aq2 - r1*aquu )
       b7   = ruuav*( aq3 - r2*aquu )
       b8   = ruuav*( aq4 - r3*aquu )
c
       aq1 = b4
       aq2 = uav*b4 + r1*b5 + b6
       aq3 = vav*b4 + r2*b5 + b7
       aq4 = wav*b4 + r3*b5 + b8
       aq5 = hav*b4 + (uu-r4)*b5 + uav*b6 + vav*b7 + wav*b8 - c2*b1/gm1

c
       aj    = 0.5*rr
       plar  = plft + prht
       eplft = elft + plft
       eprht = erht + prht
c       fssub = rr*r4
       fssub = 0.0
       f(i,1) = aj*(  rlft*uulft +  rrht*uurht           - aq1 )
     <                                          -fssub*rinf
       f(i,2) = aj*( rulft*uulft + rurht*uurht + r1*plar - aq2 )
     <                                          -fssub*rinf*uinf
       f(i,3) = aj*( rvlft*uulft + rvrht*uurht + r2*plar - aq3 )
     <                                          -fssub*rinf*vinf
       f(i,4) = aj*( rwlft*uulft + rwrht*uurht + r3*plar - aq4 )
     <                                          -fssub*rinf*winf
       f(i,5) = aj*( eplft*uulft + eprht*uurht - r4*plar - aq5 )
     <                                          -fssub*einf
   11 continue
c
      return
      end

c***********************************************************************
      subroutine rhslom(q,s,js,je,ks,ke,ls,le,bt)
c
c***********************************************************************
      use params_global
      use weno_routines
c***********************************************************************
      implicit none
c*********************************************************************c

         real q(jmax,kmax,lmax,6),s(jmax,kmax,lmax,5)
         real bt(jmax,kmax,lmax)
         integer js,je,ks,ke,ls,le

      ! local variables
         real,allocatable :: tmp(:)

         integer j,k,l,n
         real u,v,w,ge,aSq,phiSq,bSq

         allocate(tmp(5))

c***  first executable statement

         do j = js,je
         do k = ks,ke
            do l = ls,le
               u = q(j,k,l,2)/q(j,k,l,1)
               v = q(j,k,l,3)/q(j,k,l,1)
               w = q(j,k,l,4)/q(j,k,l,1)
               phiSq = 0.5*( u*u + v*v + w*w )
               ge    = q(j,k,l,5)/q(j,k,l,1)*gamma-phiSq*gm1
               aSq   = (q(j,k,l,5)/q(j,k,l,1)-phiSq)*gm1*gamma
               do n = 1,5 
                 tmp(n)=s(j,k,l,n)
               enddo
               s(j,k,l,1)=phiSq*tmp(1)-u*tmp(2)-v*tmp(3)-w*tmp(4)+tmp(5)
               s(j,k,l,2)=u*s(j,k,l,1)
               s(j,k,l,3)=v*s(j,k,l,1)
               s(j,k,l,4)=w*s(j,k,l,1)
               s(j,k,l,5)=ge*s(j,k,l,1)
       
               bSq = Mp**2/(bt(j,k,l)-Mp**2*(bt(j,k,l)-1)) 

               do n = 1,5
                  s(j,k,l,n) = s(j,k,l,n)*(gm1*(bSq-1.)/aSq)
                  s(j,k,l,n) = s(j,k,l,n) + tmp(n)
               enddo

           enddo
        enddo
        enddo

      return
      end

c******************************************************************

c***********************************************************************
      subroutine ilu3d(q,s,js,je,ks,ke,ls,le
     &      ,xx,xy,xz,yx,yy,yz,zx,zy,zz,ug,vg,wg,iblank,turmu,tscale,
     &      kkr,kkp)
c
c  calculate the implicit inversion of the lhs
c  this involves two bidiagonal scalar inversions
c
c***********************************************************************

      use params_global
      implicit none


      integer js,je,ks,ke,ls,le

      real q(jmax,kmax,lmax,nd),s(jmax,kmax,lmax,nv),
     $     turmu(jmax,kmax,lmax),tscale(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax),vg(jmax,kmax,lmax),wg(jmax,kmax,lmax)
      real xx(jmax,kmax,lmax), xy(jmax,kmax,lmax) ,xz(jmax,kmax,lmax),
     $     yx(jmax,kmax,lmax), yy(jmax,kmax,lmax) ,yz(jmax,kmax,lmax),
     $     zx(jmax,kmax,lmax), zy(jmax,kmax,lmax) ,zz(jmax,kmax,lmax)
      integer kkr(kmax),kkp(kmax)
      integer iblank(jmax,kmax,lmax)

c..   local variables

      real,allocatable :: d(:,:,:)
      real,allocatable :: a(:,:),c(:,:)
      real,allocatable :: uv(:,:),vn(:,:),ge(:,:)
      real,allocatable :: qx(:,:),qz(:,:),cx(:,:),cz(:,:)
      real,allocatable :: b(:,:,:)
      real,allocatable :: vn_j(:,:),vn_k(:,:),vn_l(:,:)
      
      real tcos,tsin,eps2,epsv,oat,dj,f1,scale,scale1,scale2
      real rhou,rhov,uu,vv,ww,er,uvw,cjkl,rj1,rj2,rj3
      real qq1,qqx,rr1,rl1,rl2,rl3,qq3,qqz,rr3,vnu,svt
      real ri1,ri2,ri3,qq,cc,qqy,sp1,sm2,spec
      real rk1,rk2,rk3,rr2,chkx,chky,chkz
      real s1,s2,s3,s4,s5,a2,a5,abcd,r2,r3,sv,di
      real visc_j,visc_k,visc_l
      real vnu_j,vnu_k,vnu_l,vnud
c     
      integer,allocatable :: ms(:),me(:)
      integer j,k,l,jj,llo,l1,m,lup,i,kr,j1

      allocate(d(jmax,kmax,lmax))
      allocate(a(mdim,5),c(mdim,5))
      allocate(uv(jmax,lmax),vn(jmax,lmax),ge(jmax,lmax))
      allocate(qx(jmax,lmax),qz(jmax,lmax),cx(jmax,lmax),cz(jmax,lmax))
      allocate(b(jmax,lmax,5))
      allocate(vn_j(jmax,lmax),vn_k(jmax,lmax),vn_l(jmax,lmax))
      allocate(ms(mdim*2),me(mdim*2))
      
c***********************************
c... viscous terms (this affects the vnu term)

      visc_j=1.0
      visc_k=1.0
      visc_l=1.0

      jm= jmax - 1
      km= kmax - 1
      lm= lmax - 1

      tcos = cos(blang)
      tsin = sin(blang)
c
      eps2 = epse*2.
      epsv = 1. + eps2
c
c..set-up for hyper-plane loop
c
      do 1 m=js+ls,je+le
      i     = max(js,m-le)
      ms(m) = i
      me(m) = min(je,m-js)
    1 continue
c
c..set time-accuracy
c
      oat = 1.0
      if (ntac.eq.2 .and. istep.gt.1) oat = 2./3.
c
c..store density,u,v,w in nonconservative variables                  
c
      do 2 l = ls-1, le+1
      do 2 k = ks-1, ke+1
      do 2 j = js-1, je+1
      dj         = 1.0 / q(j,k,l,1)
      q(j,k,l,2) = q(j,k,l,2)*dj
      q(j,k,l,3) = q(j,k,l,3)*dj
      q(j,k,l,4) = q(j,k,l,4)*dj
      q(j,k,l,5) = q(j,k,l,5)*dj
      q(j,k,l,1) = q(j,k,l,1)*q(j,k,l,6)
    2 continue
c..in the wake, average points above and below
c..for faster convergence?
      if(half.ne.1.and.1.eq.0) then
        l  = 1
        l1 = l + 1
        f1 = 0.99!vinod
        do k = ks-1, ke+1
         if(k.ge.ktip) then
          do j = 1, (jmax+1)/2
            jj = jmax - j + 1
            scale = q(j,k,l1,6)/q(j,k,l,6)
            scale1= q(jj,k,l1,6)/q(j,k,l,6)
            scale2= q(j,k,l,6)/q(jj,k,l,6)
            s(j,k,l,1)= f1*0.5*(s(j,k,l1,1)*scale+s(jj,k,l1,1)*scale1)
            s(j,k,l,2)= f1*0.5*(s(j,k,l1,2)*scale+s(jj,k,l1,2)*scale1)
            s(j,k,l,3)= f1*0.5*(s(j,k,l1,3)*scale+s(jj,k,l1,3)*scale1)
            s(j,k,l,4)= f1*0.5*(s(j,k,l1,4)*scale+s(jj,k,l1,4)*scale1)
            s(j,k,l,5)= f1*0.5*(s(j,k,l1,5)*scale+s(jj,k,l1,5)*scale1)
            s(jj,k,l,1)= s(j,k,l,1)*scale2
            s(jj,k,l,2)= s(j,k,l,2)*scale2
            s(jj,k,l,3)= s(j,k,l,3)*scale2
            s(jj,k,l,4)= s(j,k,l,4)*scale2
            s(jj,k,l,5)= s(j,k,l,5)*scale2
          enddo
         else
          do j = 1, jtail1-1
            jj = jmax - j + 1
            scale = q(j,k,l1,6)/q(j,k,l,6)
            scale1= q(jj,k,l1,6)/q(j,k,l,6)
            scale2= q(j,k,l,6)/q(jj,k,l,6)
            s(j,k,l,1)= f1*0.5*(s(j,k,l1,1)*scale+s(jj,k,l1,1)*scale1)
            s(j,k,l,2)= f1*0.5*(s(j,k,l1,2)*scale+s(jj,k,l1,2)*scale1)
            s(j,k,l,3)= f1*0.5*(s(j,k,l1,3)*scale+s(jj,k,l1,3)*scale1)
            s(j,k,l,4)= f1*0.5*(s(j,k,l1,4)*scale+s(jj,k,l1,4)*scale1)
            s(j,k,l,5)= f1*0.5*(s(j,k,l1,5)*scale+s(jj,k,l1,5)*scale1)
            s(jj,k,l,1)= s(j,k,l,1)*scale2
            s(jj,k,l,2)= s(j,k,l,2)*scale2
            s(jj,k,l,3)= s(j,k,l,3)*scale2
            s(jj,k,l,4)= s(j,k,l,4)*scale2
            s(jj,k,l,5)= s(j,k,l,5)*scale2
          enddo
         endif
        enddo
      endif
c
c..for faster convergence at front and back?
c
c..forward sweep
c..setup d, the diagonal term and b contribution in k-direction
c
      do 100 k = ks,ke
c
        do 111 l = ls-1,le+1
        do 111 j = js-1,je+1
          uu  = q(j,k,l,2)
          vv  = q(j,k,l,3)
          ww  = q(j,k,l,4)
          uvw = 0.5*( uu*uu + vv*vv + ww*ww )
          cjkl = sqrt( ggm1*( q(j,k,l,5) - uvw ) )
c
          rj1 = xx(j,k,l)
          rj2 = xy(j,k,l)
          rj3 = xz(j,k,l)
          qq1 = rj1*uu + rj2*vv + rj3*ww
          qqx = abs( -rj1*ug(j,k,l)-rj2*vg(j,k,l)-rj3*wg(j,k,l) + qq1 )
          rr1 = sqrt( rj1**2 + rj2**2 + rj3**2 )
          rk1 = yx(j,k,l)
          rk2 = yy(j,k,l)
          rk3 = yz(j,k,l)
          qqy = abs( rk1*(uu-ug(j,k,l)) + rk2*(vv-vg(j,k,l))
     <             + rk3*(ww-wg(j,k,l)) )
          rr2 = sqrt( rk1**2 + rk2**2 + rk3**2 )
          rl1 = zx(j,k,l)
          rl2 = zy(j,k,l)
          rl3 = zz(j,k,l)
          qq3 = rl1*uu + rl2*vv + rl3*ww
          qqz = abs( -rl1*ug(j,k,l)-rl2*vg(j,k,l)-rl3*wg(j,k,l) + qq3 )
          rr3 = sqrt( rl1**2 + rl2**2 + rl3**2 )

          vnu_j = visc_j*2.0*rr1*rr1*(rmue+turmu(j,k,l))/(rey*q(j,k,l,1))

          vnu_k = visc_k*2.0*rr2*rr2*(rmue+turmu(j,k,l))/(rey*q(j,k,l,1))

          vnu_l = visc_l*2.0*rr3*rr3*(rmue+turmu(j,k,l))/(rey*q(j,k,l,1))

	  vnud  = vnu_j+vnu_k+vnu_l

          svt = tscale(j,k,l)
          d(j,k,l)=1.0/
     &     ( 1.0+svt*((qqx+qqy+qqz +cjkl*(rr1+rr2+rr3))*epsv +vnud) )
c
          uv(j,l) = uvw
          qx(j,l) = qq1
          qz(j,l) = qq3
          cx(j,l) = cjkl*rr1
          cz(j,l) = cjkl*rr3
          vn_j(j,l) = vnu_j
          vn_k(j,l) = vnu_k
          vn_l(j,l) = vnu_l
          ge(j,l) = gamma*q(j,k,l,5) - gm1*uvw
  111   continue
c  
        kr=k-1
        do 112 l = ls-1,le+1
        do 112 j = js-1,je+1
          uu  = q(j,kr,l,2)
          vv  = q(j,kr,l,3)
          ww  = q(j,kr,l,4)
          er  = q(j,kr,l,5)
          uvw = 0.5*( uu*uu + vv*vv + ww*ww )
          ri1 = yx(j,kr,l)
          ri2 = yy(j,kr,l)
          ri3 = yz(j,kr,l)
          rr2 = ri1**2 + ri2**2 + ri3**2
          qq  = ri1*uu + ri2*vv + ri3*ww
          cc  = sqrt( ggm1*( er-uvw )*rr2 )
          qqy = -ri1*ug(j,kr,l)-ri2*vg(j,kr,l)-ri3*wg(j,kr,l) + qq
          sp1 = abs(qqy) + cc
          sm2 = eps2*sp1
          chky= 0.5 + sign( 0.5, qqy+cc )
          spec= chky*( qqy + sp1 ) + vn_k(j,l) + sm2 
c
          s1 = s(j,kr,l,1)
          s2 = s(j,kr,l,2)
          s3 = s(j,kr,l,3)
          s4 = s(j,kr,l,4)
          s5 = s(j,kr,l,5)
c 
          a5 = chky*( ri1*s2 + ri2*s3 + ri3*s4 - qq*s1 )
          a2 = chky*gm1*( uvw*s1 - ( uu*s2 + vv*s3 + ww*s4 ) + s5 )
          b(j,l,1) = a5 + spec*s1
          b(j,l,2) = ri1*a2 + uu*a5 +spec*s2
          b(j,l,3) = ri2*a2 + vv*a5 +spec*s3
          b(j,l,4) = ri3*a2 + ww*a5 +spec*s4
          b(j,l,5) = qq*a2 + ( gamma*er-gm1*uvw  )*a5 +spec*s5
  112   continue
c..let's precondition at the wake cut and front overlap!
      if (1.eq.0) then
      l=1
      do j=js-1,je+1
        s(j,k,l,1) = d(j,k,l)*s(j,k,l,1)
        s(j,k,l,2) = d(j,k,l)*s(j,k,l,2)
        s(j,k,l,3) = d(j,k,l)*s(j,k,l,3)
        s(j,k,l,4) = d(j,k,l)*s(j,k,l,4)
        s(j,k,l,5) = d(j,k,l)*s(j,k,l,5)
      enddo
      l=lmax
      do j=js-1,je+1
        s(j,k,l,1) = d(j,k,l)*s(j,k,l,1)
        s(j,k,l,2) = d(j,k,l)*s(j,k,l,2)
        s(j,k,l,3) = d(j,k,l)*s(j,k,l,3)
        s(j,k,l,4) = d(j,k,l)*s(j,k,l,4)
        s(j,k,l,5) = d(j,k,l)*s(j,k,l,5)
      enddo
      endif
c
c..loop on hyper-plane
c
      do 120 m = js+ls,je+le
c
c..setup a contribution in j-direction
c
      do 121 j = ms(m),me(m)
        l  = m-j
        j1 = j-1
c
        uu  = q(j1,k,l,2)
        vv  = q(j1,k,l,3)
        ww  = q(j1,k,l,4)
        uvw = uv(j1,l)
        ri1 = xx(j1,k,l)
        ri2 = xy(j1,k,l)
        ri3 = xz(j1,k,l)
        qq  = qx(j1,l)
        cc  = cx(j1,l)
        qqx = -ri1*ug(j1,k,l)-ri2*vg(j1,k,l)-ri3*wg(j1,k,l) + qq
        sp1 = abs(qqx) + cc
        sm2 = eps2*sp1
        chkx= 0.5 + sign( 0.5, qqx+cc )
        spec= chkx*( qqx + sp1 ) + vn_j(j1,l) + sm2
c
        s1 = s(j1,k,l,1)
        s2 = s(j1,k,l,2)
        s3 = s(j1,k,l,3)
        s4 = s(j1,k,l,4)
        s5 = s(j1,k,l,5)
c
        a5 = chkx*( ri1*s2 + ri2*s3 + ri3*s4 - qq*s1 )
        a2 = chkx*gm1*( uvw*s1 - ( uu*s2 + vv*s3 + ww*s4 ) + s5 )
        a(j,1) = a5                  +spec*s1
        a(j,2) = ri1*a2 + uu*a5      +spec*s2
        a(j,3) = ri2*a2 + vv*a5      +spec*s3
        a(j,4) = ri3*a2 + ww*a5      +spec*s4
        a(j,5) = qq*a2 + ge(j1,l)*a5 +spec*s5
  121 continue
c
c..setup c contribution in l-direction
c
      do 122 j = ms(m),me(m)
        l  = m-j
        l1 = l-1
c
        uu  = q(j,k,l1,2)
        vv  = q(j,k,l1,3)
        ww  = q(j,k,l1,4)
        uvw = uv(j,l1)
        ri1 = zx(j,k,l1)
        ri2 = zy(j,k,l1)
        ri3 = zz(j,k,l1)
        qq  = qz(j,l1)
        cc  = cz(j,l1)
        qqz = -ri1*ug(j,k,l1)-ri2*vg(j,k,l1)-ri3*wg(j,k,l1) + qq
        sp1 = abs(qqz) + cc
        sm2 = eps2*sp1
        chkz= 0.5 + sign( 0.5, qqz+cc )
        spec= chkz*( qqz + sp1 ) + vn_l(j,l1) + sm2
c
        s1 = s(j,k,l1,1)
        s2 = s(j,k,l1,2)
        s3 = s(j,k,l1,3)
        s4 = s(j,k,l1,4)
        s5 = s(j,k,l1,5)
c
        a5 = chkz*( ri1*s2 + ri2*s3 + ri3*s4 - qq*s1 )
        a2 = chkz*gm1*( uvw*s1 - ( uu*s2 + vv*s3 + ww*s4 ) + s5 )
        c(j,1) = a5                  +spec*s1
        c(j,2) = ri1*a2 + uu*a5      +spec*s2
        c(j,3) = ri2*a2 + vv*a5      +spec*s3
        c(j,4) = ri3*a2 + ww*a5      +spec*s4
        c(j,5) = qq*a2 + ge(j,l1)*a5 +spec*s5
  122 continue
c
c..bi-diagonal inversion, include source term for quasi-steady
c     
      if(iunst.eq.2) then
cdir$ ivdep
        do 123 j = ms(m),me(m)
          l  = m-j
          svt = tscale(j,k,l)
          sv = svt*0.5
          di = d(j,k,l)
          s(j,k,l,1) = ( s(j,k,l,1) + sv*( a(j,1)+b(j,l,1)+c(j,1) ) )*di
          r2 = ( s(j,k,l,2) + sv*(a(j,2)+b(j,l,2)+c(j,2) ) ) * di
          r3 = ( s(j,k,l,3) + sv*(a(j,3)+b(j,l,3)+c(j,3) ) ) * di
          abcd = 1.0 + (2.*di*sv*rf)**2
          s(j,k,l,2) = ( r2 + di*2.*sv*rf*r3)/abcd
          s(j,k,l,3) = ( r3 - di*2.*sv*rf*r2)/abcd
          s(j,k,l,4) = ( s(j,k,l,4) + sv*( a(j,4)+b(j,l,4)+c(j,4) ) )*di
          s(j,k,l,5) = ( s(j,k,l,5) + sv*( a(j,5)+b(j,l,5)+c(j,5) ) )*di
  123   continue
      else
c
c..bi-diagonal inversion
c     
cdir$ ivdep
        do 124 j = ms(m),me(m)
          l  = m-j
          svt = tscale(j,k,l)
          sv = svt*0.5
          di = d(j,k,l)
          s(j,k,l,1) = ( s(j,k,l,1) + sv*( a(j,1)+b(j,l,1)+c(j,1) ) )*di
          s(j,k,l,2) = ( s(j,k,l,2) + sv*( a(j,2)+b(j,l,2)+c(j,2) ) )*di
          s(j,k,l,3) = ( s(j,k,l,3) + sv*( a(j,3)+b(j,l,3)+c(j,3) ) )*di
          s(j,k,l,4) = ( s(j,k,l,4) + sv*( a(j,4)+b(j,l,4)+c(j,4) ) )*di
          s(j,k,l,5) = ( s(j,k,l,5) + sv*( a(j,5)+b(j,l,5)+c(j,5) ) )*di
 124    continue
      endif
c
  120 continue
  100 continue
c
c..backward sweep
c..b contribution in k-direction
c
      do 200 k = ke,ks,-1
c
        kr=k+1
        do 211 l = ls-1,le+1
        do 211 j = js-1,je+1
          uu  = q(j,k,l,2)
          vv  = q(j,k,l,3)
          ww  = q(j,k,l,4)
          uvw = 0.5*( uu*uu + vv*vv + ww*ww )
          cjkl = sqrt( ggm1*( q(j,k,l,5) - uvw ) )
c
          rj1 = xx(j,k,l)
          rj2 = xy(j,k,l)
          rj3 = xz(j,k,l)
          qq1 = rj1*uu + rj2*vv + rj3*ww
c
          rr1 = sqrt( rj1**2 + rj2**2 + rj3**2 )
          rl1 = zx(j,k,l)
          rl2 = zy(j,k,l)
          rl3 = zz(j,k,l)
          qq3 = rl1*uu + rl2*vv + rl3*ww
c
          rr3 = sqrt( rl1**2 + rl2**2 + rl3**2 )
          
	  rk1 = yx(j,k,l)
          rk2 = yy(j,k,l)
          rk3 = yz(j,k,l)
          rr2 = sqrt( rk1**2 + rk2**2 + rk3**2 )

          vnu_j = visc_j*2.0*rr1*rr1*(rmue+turmu(j,k,l))/(rey*q(j,k,l,1))
          vnu_k = visc_k*2.0*rr2*rr2*(rmue+turmu(j,k,l))/(rey*q(j,k,l,1))
          vnu_l = visc_l*2.0*rr3*rr3*(rmue+turmu(j,k,l))/(rey*q(j,k,l,1))
	  
          uv(j,l) = uvw
          qx(j,l) = qq1
          qz(j,l) = qq3
          cx(j,l) = cjkl*rr1
          cz(j,l) = cjkl*rr3
          vn_j(j,l) = vnu_j
          vn_k(j,l) = vnu_k
          vn_l(j,l) = vnu_l
          ge(j,l) = gamma*q(j,k,l,5) - gm1*uvw
  211   continue
c
        do 212 l = ls-1,le+1
        do 212 j = js-1,je+1
          uu  = q(j,kr,l,2)
          vv  = q(j,kr,l,3)
          ww  = q(j,kr,l,4)
          er  = q(j,kr,l,5)
          uvw = 0.5*( uu*uu + vv*vv + ww*ww )
          ri1 = yx(j,kr,l)
          ri2 = yy(j,kr,l)
          ri3 = yz(j,kr,l)
          rr2 = ri1**2 + ri2**2 + ri3**2
          qq  = ri1*uu + ri2*vv + ri3*ww
          cc  = sqrt( ggm1*( er-uvw )*rr2 )
          qqy = -ri1*ug(j,kr,l)-ri2*vg(j,kr,l)-ri3*wg(j,kr,l) + qq
          sp1 = abs(qqy) + cc
          sm2 = eps2*sp1
          chky= 0.5 - sign( 0.5, qqy-cc )
          spec= chky*( qqy - sp1 ) - vn_k(j,l) - sm2
c
          s1 = s(j,kr,l,1)
          s2 = s(j,kr,l,2)
          s3 = s(j,kr,l,3)
          s4 = s(j,kr,l,4)
          s5 = s(j,kr,l,5)
c
          a5 = chky*( ri1*s2 + ri2*s3 + ri3*s4 - qq*s1 )
          a2 = chky*gm1*( uvw*s1 - ( uu*s2 + vv*s3 + ww*s4 ) + s5 )
          b(j,l,1) = a5 + spec*s1
          b(j,l,2) = ri1*a2 + uu*a5 +spec*s2
          b(j,l,3) = ri2*a2 + vv*a5 +spec*s3
          b(j,l,4) = ri3*a2 + ww*a5 +spec*s4
          b(j,l,5) = qq*a2 + ( gamma*er-gm1*uvw  )*a5 +spec*s5
  212   continue
c
c..loop on hyper-plane
c
      do 220 m = je+le,js+ls,-1
c
c..setup a contribution in j-direction
c
      do 221 j = ms(m),me(m)
        l  = m-j
        j1 = j+1
c
        uu  = q(j1,k,l,2)
        vv  = q(j1,k,l,3)
        ww  = q(j1,k,l,4)
        uvw = uv(j1,l)
        ri1 = xx(j1,k,l)
        ri2 = xy(j1,k,l)
        ri3 = xz(j1,k,l)
        qq  = qx(j1,l)
        cc  = cx(j1,l)
        qqx = -ri1*ug(j1,k,l)-ri2*vg(j1,k,l)-ri3*wg(j1,k,l) + qq
        sp1 = abs(qqx) + cc
        sm2 = eps2*sp1
        chkx= 0.5 - sign( 0.5, qqx-cc )
        spec= chkx*( qqx - sp1 ) - vn_j(j1,l) - sm2
c
        s1 = s(j1,k,l,1)
        s2 = s(j1,k,l,2)
        s3 = s(j1,k,l,3)
        s4 = s(j1,k,l,4)
        s5 = s(j1,k,l,5)
c
        a5 = chkx*( ri1*s2 + ri2*s3 + ri3*s4 - qq*s1 )
        a2 = chkx*gm1*( uvw*s1 - ( uu*s2 + vv*s3 + ww*s4 ) + s5 )
        a(j,1) = a5                  +spec*s1
        a(j,2) = ri1*a2 + uu*a5      +spec*s2
        a(j,3) = ri2*a2 + vv*a5      +spec*s3
        a(j,4) = ri3*a2 + ww*a5      +spec*s4
        a(j,5) = qq*a2 + ge(j1,l)*a5 +spec*s5
  221 continue
c
c..setup c contribution in l-direction
c
      do 222 j = ms(m),me(m)
        l  = m-j
        l1 = l+1
c
        uu  = q(j,k,l1,2)
        vv  = q(j,k,l1,3)
        ww  = q(j,k,l1,4)
        uvw = uv(j,l1)
        ri1 = zx(j,k,l1)
        ri2 = zy(j,k,l1)
        ri3 = zz(j,k,l1)
        qq  = qz(j,l1)
        cc  = cz(j,l1)
        qqz = -ri1*ug(j,k,l1)-ri2*vg(j,k,l1)-ri3*wg(j,k,l1) + qq
        sp1 = abs(qqz) + cc
        sm2 = eps2*sp1
        chkz= 0.5 - sign( 0.5, qqz-cc )
        spec= chkz*( qqz - sp1 ) - vn_l(j,l1) - sm2
c
        s1 = s(j,k,l1,1)
        s2 = s(j,k,l1,2)
        s3 = s(j,k,l1,3)
        s4 = s(j,k,l1,4)
        s5 = s(j,k,l1,5)
c
        a5 = chkz*( ri1*s2 + ri2*s3 + ri3*s4 - qq*s1 )
        a2 = chkz*gm1*( uvw*s1 - ( uu*s2 + vv*s3 + ww*s4 ) + s5 )
        c(j,1) = a5                  +spec*s1
        c(j,2) = ri1*a2 + uu*a5      +spec*s2
        c(j,3) = ri2*a2 + vv*a5      +spec*s3
        c(j,4) = ri3*a2 + ww*a5      +spec*s4
        c(j,5) = qq*a2 + ge(j,l1)*a5 +spec*s5
  222 continue
c
c..bi-diagonal inversion
c     
cdir$ ivdep
      do 223 j = ms(m),me(m)
        l  = m-j
        svt = tscale(j,k,l)
        sv = svt*0.5
        di = d  (j,k,l)
        s(j,k,l,1) = s(j,k,l,1) - sv*( a(j,1)+b(j,l,1)+c(j,1) )*di
        s(j,k,l,2) = s(j,k,l,2) - sv*( a(j,2)+b(j,l,2)+c(j,2) )*di
        s(j,k,l,3) = s(j,k,l,3) - sv*( a(j,3)+b(j,l,3)+c(j,3) )*di
        s(j,k,l,4) = s(j,k,l,4) - sv*( a(j,4)+b(j,l,4)+c(j,4) )*di
        s(j,k,l,5) = s(j,k,l,5) - sv*( a(j,5)+b(j,l,5)+c(j,5) )*di
  223 continue
c
  220 continue
  200 continue
c
c..restore density , u , v , w in conservative variables             
c                                                                       
      do 4 l = ls-1, le+1
      do 4 k = ks-1, ke+1
      do 4 j = js-1, je+1
        q(j,k,l,1) = q(j,k,l,1) / q(j,k,l,6)
        q(j,k,l,2) = q(j,k,l,2)*q(j,k,l,1)
        q(j,k,l,3) = q(j,k,l,3)*q(j,k,l,1)
        q(j,k,l,4) = q(j,k,l,4)*q(j,k,l,1)
        q(j,k,l,5) = q(j,k,l,5)*q(j,k,l,1)
    4 continue
                                                          
      return
      end

c***********************************************************************
      subroutine arc3d(q,s,js,je,ks,ke,ls,le,xx,xy,xz,yx,yy,yz,zx,zy,zz,
     &                 ug,vg,wg,iblank,turmu,tscale,bt,kkr,kkp)
c
c  calculate the implicit inversion of the lhs
c
c***********************************************************************
      
      use params_global
      implicit none

      integer js,je,ks,ke,ls,le

      real q(jmax,kmax,lmax,nd),s(jmax,kmax,lmax,nv),
     c     turmu(jmax,kmax,lmax),tscale(jmax,kmax,lmax),
     c     bt(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax),vg(jmax,kmax,lmax),wg(jmax,kmax,lmax)
      real xx(jmax,kmax,lmax), xy(jmax,kmax,lmax) ,xz(jmax,kmax,lmax),
     c     yx(jmax,kmax,lmax), yy(jmax,kmax,lmax) ,yz(jmax,kmax,lmax),
     c     zx(jmax,kmax,lmax), zy(jmax,kmax,lmax) ,zz(jmax,kmax,lmax)
      integer kkr(kmax),kkp(kmax)
      integer iblank(jmax,kmax,lmax)

c..   local variables

      real, allocatable :: p(:,:,:),vn(:,:,:)
      real, allocatable :: qx(:,:,:),qy(:,:,:),qz(:,:,:)
      real, allocatable :: cx(:,:,:),cy(:,:,:),cz(:,:,:)
      real, allocatable :: d2px(:,:,:),d2py(:,:,:),d2pz(:,:,:)
      real, allocatable :: tmp(:),uv(:,:,:)
      
      integer j,k,l,n
      real dj,uu,vv,ww,uvw,cjkl,ge
      real qq1,rj1,rj2,rj3,rr1
      real qq2,rk1,rk2,rk3,rr2
      real qq3,rl1,rl2,rl3,rr3
      real c2i
      real a1,a2,a3,a4,a5,a6
      real r2,r3,abcd

      real vnu
      real visc_j,visc_k,visc_l

c...   memory allocation block

       allocate(p(jmax,kmax,lmax),vn(jmax,kmax,lmax))
       allocate(qx(jmax,kmax,lmax),qy(jmax,kmax,lmax),
     c                             qz(jmax,kmax,lmax))
       allocate(cx(jmax,kmax,lmax),cy(jmax,kmax,lmax),
     c                             cz(jmax,kmax,lmax))
       allocate(d2px(jmax,kmax,lmax),d2py(jmax,kmax,lmax),
     c                               d2pz(jmax,kmax,lmax))
       allocate(tmp(nv),uv(jmax,kmax,lmax))

c***********************************************************************
c***  first executable statement
c
c... viscous terms (this affects the vnu term)
      visc_j = 1.0
      visc_k = 1.0
      visc_l = 1.0
c
c..store density,u,v,w in nonconservative variables                  
c

      do 2 l = ls-1, le+1
      do 2 k = ks-1, ke+1
      do 2 j = js-1, je+1
         dj  = 1.0 / q(j,k,l,1)
         q(j,k,l,2) = q(j,k,l,2)*dj
         q(j,k,l,3) = q(j,k,l,3)*dj
         q(j,k,l,4) = q(j,k,l,4)*dj
         q(j,k,l,5) = q(j,k,l,5)*dj
         q(j,k,l,1) = q(j,k,l,1)*q(j,k,l,6)
    2 continue
c
c   multiply by T_eta^inv
c
      do l = ls-1,le+1
      do k = ks-1,ke+1
         do j = js-1,je+1
            uu  = q(j,k,l,2)
            vv  = q(j,k,l,3)
            ww  = q(j,k,l,4)
         
            uvw  = 0.5*( uu*uu + vv*vv + ww*ww )
            cjkl = sqrt(ggm1*( q(j,k,l,5) - uvw ))
            c2i  = 1.0/(cjkl*cjkl)
            p(j,k,l) = ( q(j,k,l,5) - uvw )*gm1*q(j,k,l,1)
         
            rj1 = xx(j,k,l)
            rj2 = xy(j,k,l)
            rj3 = xz(j,k,l)
            qq1 = rj1*( uu - ug(j,k,l) ) + rj2*( vv - vg(j,k,l) )
     c                                   + rj3*( ww - wg(j,k,l) )

            rr1 = sqrt( rj1*rj1 + rj2*rj2 + rj3*rj3 )
            rj1 = rj1/rr1
            rj2 = rj2/rr1
            rj3 = rj3/rr1
         
            rk1 = yx(j,k,l)
            rk2 = yy(j,k,l)
            rk3 = yz(j,k,l)
            qq2 = rk1*( uu - ug(j,k,l) ) + rk2*( vv - vg(j,k,l) )
     c                                   + rk3*( ww - wg(j,k,l) )

            rr2 = sqrt( rk1*rk1 + rk2*rk2 + rk3*rk3 )
            rk1 = rk1/rr2
            rk2 = rk2/rr2
            rk3 = rk3/rr2
         
            rl1 = zx(j,k,l)
            rl2 = zy(j,k,l)
            rl3 = zz(j,k,l)
            qq3 = rl1*( uu - ug(j,k,l) ) + rl2*( vv - vg(j,k,l) )
     c                                   + rl3*( ww - wg(j,k,l) )

            rr3 = sqrt( rl1*rl1 + rl2*rl2 + rl3*rl3 )
            rl1 = rl1/rr3
            rl2 = rl2/rr3
            rl3 = rl3/rr3
            
            vnu = (rmue+turmu(j,k,l))/(rey*q(j,k,l,1))
 
            uv(j,k,l) = uvw
            qx(j,k,l) = qq1
            cx(j,k,l) = cjkl*rr1
            qy(j,k,l) = qq2
            cy(j,k,l) = cjkl*rr2
            qz(j,k,l) = qq3
            cz(j,k,l) = cjkl*rr3
            vn(j,k,l) = visc_k*vnu
          
            a1 = s(j,k,l,2)*uu + s(j,k,l,3)*vv + ww*s(j,k,l,4)
     c                         - s(j,k,l,5)
            a1 = a1*gm1*c2i
            a1 = a1 + s(j,k,l,1)*(1. - uvw*GM1*c2i)
         
            a2 = (rk2*ww-rk3*vv)*s(j,k,l,1) + rk3*s(j,k,l,3)
     c                                      - rk2*s(j,k,l,4)
            a3 = (rk3*uu-rk1*ww)*s(j,k,l,1) + rk1*s(j,k,l,4)
     c                                      - rk3*s(j,k,l,2)
            a4 = (rk1*vv-rk2*uu)*s(j,k,l,1) + rk2*s(j,k,l,2)
     c                                      - rk1*s(j,k,l,3)
            a5 = uvw*s(j,k,l,1) - s(j,k,l,2)*uu - s(j,k,l,3)*vv
     c                             - s(j,k,l,4)*ww + s(j,k,l,5)
            a5 = a5*gm1*c2i
         
            a6 = qq2*s(j,k,l,1)/rr2 - rk1*s(j,k,l,2)
     c          -rk2*s(j,k,l,3)     - rk3*s(j,k,l,4)
            a6 = a6
         
            s(j,k,l,1) = a1*rk1 + a2/q(j,k,l,1)
            s(j,k,l,2) = a1*rk2 + a3/q(j,k,l,1)
            s(j,k,l,3) = a1*rk3 + a4/q(j,k,l,1)
            s(j,k,l,4) =  0.5*(a5-a6/cjkl)
            s(j,k,l,5) =  0.5*(a5+a6/cjkl)
         enddo
      enddo
      enddo
         
      do l = ls-1,le+1
      do k = ks,ke
         do j = js-1,je+1
            d2py(j,k,l) = abs( p(j,k+1,l) - 2.*p(j,k,l) + p(j,k-1,l) ) 
            d2py(j,k,l) = d2py(j,k,l)/abs( p(j,k+1,l) + 2.*p(j,k,l) 
     c                                    + p(j,k-1,l) )
         enddo
      enddo
      enddo

      k = ks-1
      do l = ls-1,le+1
         do j = js-1,je+1
            d2py(j,k,l) = 0. 
         enddo
      enddo

      k = ke+1
      do l = ls-1,le+1
         do j = js-1,je+1
            d2py(j,k,l) = 0. 
         enddo
      enddo

c   K direction
      if (ilhs.eq.2) then
        call lhsinv3d(q,s,qy,cy,d2py,ks-1,ke+1,js-1,je+1,ls-1,le+1,2,
     c                              yx,yy,yz,vn,iblank,tscale,bt)
      elseif (ilhs.eq.3) then
        call lhsinv3d_up(q,s,qy,cy,d2py,ks-1,ke+1,js-1,je+1,ls-1,le+1,2,
     c                              yx,yy,yz,vn,iblank,tscale,bt)
      endif

c multiply by N^inv = T_si^(-1) T_eta 

      do l=ls-1,le+1
      do k=ks-1,ke+1
         do j=js-1,je+1

            cjkl = sqrt(ggm1*( q(j,k,l,5) - uv(j,k,l) ))
            rj1 = xx(j,k,l)
            rj2 = xy(j,k,l)
            rj3 = xz(j,k,l)
            rr1 = sqrt( rj1*rj1 + rj2*rj2 + rj3*rj3 )
            rj1 = rj1/rr1
            rj2 = rj2/rr1
            rj3 = rj3/rr1
        
            rk1 = yx(j,k,l)
            rk2 = yy(j,k,l)
            rk3 = yz(j,k,l)
            rr2 = sqrt( rk1*rk1 + rk2*rk2 + rk3*rk3 )
            rk1 = rk1/rr2
            rk2 = rk2/rr2
            rk3 = rk3/rr2

            vnu = (rmue+turmu(j,k,l))/(rey*q(j,k,l,1))
            vn(j,k,l) = visc_j*vnu

            a1 = ( rj1*rk1 + rj2*rk2 + rj3*rk3 )
            a2 = ( rj1*rk2 - rk1*rj2 )
            a3 = ( rj3*rk2 - rk3*rj2 )
            a4 = ( rj1*rk3 - rk1*rj3 )
            a5 = cjkl*( s(j,k,l,4) - s(j,k,l,5) )/(q(j,k,l,1))
            a6 = ( s(j,k,l,4) + s(j,k,l,5) )
        
            tmp(1) =   a1*s(j,k,l,1) + a2*s(j,k,l,2) + a4*s(j,k,l,3)
     c               + a5*a3
            tmp(2) = - a2*s(j,k,l,1) + a1*s(j,k,l,2) - a3*s(j,k,l,3)
     c               + a5*a4
            tmp(3) = - a4*s(j,k,l,1) + a3*s(j,k,l,2) + a1*s(j,k,l,3)
     c               - a5*a2
            tmp(4) = 0.5*(- q(j,k,l,1)*(a3*s(j,k,l,1) + a4*s(j,k,l,2)
     c                - a2*s(j,k,l,3) - a5*a1)/cjkl + a6)
            tmp(5) = 0.5*(  q(j,k,l,1)*(a3*s(j,k,l,1) + a4*s(j,k,l,2)
     c                - a2*s(j,k,l,3) - a5*a1)/cjkl + a6)

            do n = 1,5
               s(j,k,l,n) = tmp(n) 
            enddo
         enddo
      enddo
      enddo

      do l = ls-1,le+1
         do k = ks-1,ke+1
            do j = js,je
               d2px(j,k,l) =abs(p(j+1,k,l) - 2.*p(j,k,l) + p(j-1,k,l))
               d2px(j,k,l) = d2px(j,k,l)/abs( p(j+1,k,l)
     c                                     + 2.*p(j,k,l) + p(j-1,k,l))
            enddo
         enddo
      enddo

      j=js-1
      do l=ls-1,le+1
         do k=ks-1,ke+1
            d2px(j,k,l) = 0. 
         enddo
      enddo

      j=je+1
      do l=ls-1,le+1
         do k=ks-1,ke+1
            d2px(j,k,l) = 0. 
         enddo
      enddo

c  J direction

      if (ilhs.eq.2) then
        call lhsinv3d(q,s,qx,cx,d2px,js-1,je+1,ks-1,ke+1,ls-1,le+1,1,
     c                       xx,xy,xz,vn,iblank,tscale,bt)
      elseif (ilhs.eq.3) then
        call lhsinv3d_up(q,s,qx,cx,d2px,js-1,je+1,ks-1,ke+1,ls-1,le+1,1,
     c                       xx,xy,xz,vn,iblank,tscale,bt)
      endif
 
c  multiply by N^inv = T_zeta^(-1) T_si

      do l=ls-1,le+1
      do k=ks-1,ke+1
         do j=js-1,je+1

            cjkl = sqrt(ggm1*( q(j,k,l,5) - uv(j,k,l) ))
            rj1 = xx(j,k,l)
            rj2 = xy(j,k,l)
            rj3 = xz(j,k,l)
            rr1 = sqrt( rj1*rj1 + rj2*rj2 + rj3*rj3 )
            rj1 = rj1/rr1
            rj2 = rj2/rr1
            rj3 = rj3/rr1
        
            rl1 = zx(j,k,l)
            rl2 = zy(j,k,l)
            rl3 = zz(j,k,l)
            rr3 = sqrt( rl1*rl1 + rl2*rl2 + rl3*rl3 )
            rl1 = rl1/rr3
            rl2 = rl2/rr3
            rl3 = rl3/rr3
        
            vnu = (rmue+turmu(j,k,l))/(rey*q(j,k,l,1))
            vn(j,k,l) = visc_l*vnu

            a1 = (rj1*rl1+rj2*rl2+rj3*rl3)
            a2 = (rl1*rj2-rj1*rl2)
            a3 = (rl3*rj2-rj3*rl2)
            a4 = (rl1*rj3-rj1*rl3)
            a5 = cjkl*( s(j,k,l,4) - s(j,k,l,5) )/(q(j,k,l,1))
            a6 = ( s(j,k,l,4) + s(j,k,l,5) )
        
            tmp(1) =   a1*s(j,k,l,1) + a2*s(j,k,l,2) + a4*s(j,k,l,3)
     c               + a5*a3
            tmp(2) = - a2*s(j,k,l,1) + a1*s(j,k,l,2) - a3*s(j,k,l,3)
     c               + a5*a4
            tmp(3) = - a4*s(j,k,l,1) + a3*s(j,k,l,2) + a1*s(j,k,l,3)
     c               - a5*a2
            tmp(4) = 0.5*(- q(j,k,l,1)*(a3*s(j,k,l,1) + a4*s(j,k,l,2)
     c                - a2*s(j,k,l,3) - a5*a1)/cjkl + a6)
            tmp(5) = 0.5*(  q(j,k,l,1)*(a3*s(j,k,l,1) + a4*s(j,k,l,2)
     c                - a2*s(j,k,l,3) - a5*a1)/cjkl + a6)

            do n = 1,5
               s(j,k,l,n) = tmp(n) 
            enddo
         enddo
      enddo
      enddo

      do l = ls,le
      do k = ks-1,ke+1
         do j = js-1,je+1
            d2pz(j,k,l) = abs( p(j,k,l+1) - 2.*p(j,k,l) + p(j,k,l-1) ) 
            d2pz(j,k,l) = d2pz(j,k,l)/abs( p(j,k,l+1) + 2.*p(j,k,l) 
     c                                   + p(j,k,l-1) )
         enddo
      enddo
      enddo

      l = ls-1
      do k = ks-1,ke+1
         do j = js-1,je+1
            d2pz(j,k,l) = 0. 
         enddo
      enddo

      l = le+1
      do k = ks-1,ke+1
         do j = js-1,je+1
            d2pz(j,k,l) = 0. 
         enddo
      enddo
c   L direction
      if (ilhs.eq.2) then
        call lhsinv3d(q,s,qz,cz,d2pz,ls-1,le+1,js-1,je+1,ks-1,ke+1,3,
     c                             zx,zy,zz,vn,iblank,tscale,bt)
      elseif (ilhs.eq.3) then
        call lhsinv3d_up(q,s,qz,cz,d2pz,ls-1,le+1,js-1,je+1,ks-1,ke+1,3,
     c                             zx,zy,zz,vn,iblank,tscale,bt)
      endif
 
c  multiply by T_zeta 

      do l=ls-1,le+1
      do k=ks-1,ke+1
         do j=js-1,je+1
            uu  = q(j,k,l,2)
            vv  = q(j,k,l,3)
            ww  = q(j,k,l,4)
         
            uvw  = 0.5*( uu*uu + vv*vv + ww*ww )
            cjkl = sqrt(ggm1*( q(j,k,l,5) - uvw ))
            c2i  = 1.0/(cjkl*cjkl)
            ge    = gamma*q(j,k,l,5)- gm1*uvw
            qq3  = qz(j,k,l) 
        
            rl1 = zx(j,k,l)
            rl2 = zy(j,k,l)
            rl3 = zz(j,k,l)
            rr3 = sqrt( rl1*rl1 + rl2*rl2 + rl3*rl3 )
            rl1 = rl1/rr3
            rl2 = rl2/rr3
            rl3 = rl3/rr3
        
            a1 = rl1*s(j,k,l,1) + rl2*s(j,k,l,2) + rl3*s(j,k,l,3)
     c          +s(j,k,l,4) + s(j,k,l,5)
            a2 = cjkl*(s(j,k,l,4) - s(j,k,l,5))
            a3 = uvw*(rl1*s(j,k,l,1) + rl2*s(j,k,l,2) + rl3*s(j,k,l,3))
         
            tmp(1) = a1
            tmp(2) = a1*uu-q(j,k,l,1)*(rl3*s(j,k,l,2) - rl2*s(j,k,l,3))
            tmp(2) = tmp(2) + a2*rl1
            tmp(3) = a1*vv+q(j,k,l,1)*(rl3*s(j,k,l,1) - rl1*s(j,k,l,3))
            tmp(3) = tmp(3) + a2*rl2
            tmp(4) = a1*ww-q(j,k,l,1)*(rl2*s(j,k,l,1) - rl1*s(j,k,l,2))
            tmp(4) = tmp(4) + a2*rl3;
            tmp(5) = a3 + q(j,k,l,1)*((vv*rl3 - ww*rl2)*s(j,k,l,1)
     c                  + (ww*rl1 - uu*rl3)*s(j,k,l,2)
     c                  + (uu*rl2 - vv*rl1)*s(j,k,l,3))
            tmp(5) = tmp(5) + (ge+cjkl*qq3/(rr3))*s(j,k,l,4) 
     c                      + (ge-cjkl*qq3/(rr3))*s(j,k,l,5)

            do n = 1,5
               s(j,k,l,n) = tmp(n) 
            enddo
         enddo
      enddo
      enddo

      if (iunst.eq.2) then
         do l=ls-1,le+1
         do k=ks-1,ke+1
            do j=js-1,je+1
             r2 = s(j,k,l,2)
             r3 = s(j,k,l,3)
             abcd = 1.0 + (tscale(j,k,l)*rf)**2
             s(j,k,l,2) = ( r2 + tscale(j,k,l)*rf*r3)/abcd
             s(j,k,l,3) = ( r3 - tscale(j,k,l)*rf*r2)/abcd
            enddo
         enddo
         enddo
      endif
c
c..restore density , u , v , w in conservative variables             
c                                                                       
      do 4 l = ls-1, le+1
      do 4 k = ks-1, ke+1
      do 4 j = js-1, je+1
        q(j,k,l,1) = q(j,k,l,1) / q(j,k,l,6)
        q(j,k,l,2) = q(j,k,l,2)*q(j,k,l,1)
        q(j,k,l,3) = q(j,k,l,3)*q(j,k,l,1)
        q(j,k,l,4) = q(j,k,l,4)*q(j,k,l,1)
        q(j,k,l,5) = q(j,k,l,5)*q(j,k,l,1)
    4 continue

      return
      end


c***********************************************************************
      subroutine arc3d_precon(q,s,js,je,ks,ke,ls,le,xx,xy,xz,yx,yy,yz,
     &                 zx,zy,zz,ug,vg,wg,iblank,turmu,tscale,bt,kkr,kkp)
c
c  calculate the implicit inversion of the lhs
c
c***********************************************************************
      
      use params_global
      implicit none

      integer js,je,ks,ke,ls,le

      real q(jmax,kmax,lmax,nd),s(jmax,kmax,lmax,nv),
     c     turmu(jmax,kmax,lmax),tscale(jmax,kmax,lmax),
     c     bt(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax),vg(jmax,kmax,lmax),wg(jmax,kmax,lmax)
      real xx(jmax,kmax,lmax), xy(jmax,kmax,lmax) ,xz(jmax,kmax,lmax),
     c     yx(jmax,kmax,lmax), yy(jmax,kmax,lmax) ,yz(jmax,kmax,lmax),
     c     zx(jmax,kmax,lmax), zy(jmax,kmax,lmax) ,zz(jmax,kmax,lmax)
      integer kkr(kmax),kkp(kmax)
      integer iblank(jmax,kmax,lmax)

c..   local variables

      real, allocatable :: p(:,:,:),vn(:,:,:)
      real, allocatable :: qx(:,:,:),qy(:,:,:),qz(:,:,:)
      real, allocatable :: cx(:,:,:),cy(:,:,:),cz(:,:,:)
      real, allocatable :: d2px(:,:,:),d2py(:,:,:),d2pz(:,:,:)
      real, allocatable :: tmp(:),uv(:,:,:)
      
      integer j,k,l,n
      real dj,uu,vv,ww,uvw,cjkl,ge
      real qq1,rj1,rj2,rj3,rr1
      real qq2,rk1,rk2,rk3,rr2
      real qq3,rl1,rl2,rl3,rr3
      real c2i,X,Y,Z,X1,Y1,Z1,bSq
      real a1,a2,a3,a4,a5,a6
      real r2,r3,abcd

      real vnu
      real visc_j,visc_k,visc_l

c...   memory allocation block

       allocate(p(jmax,kmax,lmax),vn(jmax,kmax,lmax))
       allocate(qx(jmax,kmax,lmax),qy(jmax,kmax,lmax),
     c                             qz(jmax,kmax,lmax))
       allocate(cx(jmax,kmax,lmax),cy(jmax,kmax,lmax),
     c                             cz(jmax,kmax,lmax))
       allocate(d2px(jmax,kmax,lmax),d2py(jmax,kmax,lmax),
     c                               d2pz(jmax,kmax,lmax))
       allocate(tmp(nv),uv(jmax,kmax,lmax))

c***********************************************************************
c***  first executable statement
c
c... viscous terms (this affects the vnu term)
      visc_j = 1.0
      visc_k = 1.0
      visc_l = 1.0
c
c..store density,u,v,w in nonconservative variables                  
c
      do 2 l = ls-1, le+1
      do 2 k = ks-1, ke+1
      do 2 j = js-1, je+1
         dj  = 1.0 / q(j,k,l,1)
         q(j,k,l,2) = q(j,k,l,2)*dj
         q(j,k,l,3) = q(j,k,l,3)*dj
         q(j,k,l,4) = q(j,k,l,4)*dj
         q(j,k,l,5) = q(j,k,l,5)*dj
         q(j,k,l,1) = q(j,k,l,1)*q(j,k,l,6)
    2 continue
c
c   multiply by T_eta^inv
c
      do l = ls-1,le+1
      do k = ks-1,ke+1
         do j = js-1,je+1
            uu  = q(j,k,l,2)
            vv  = q(j,k,l,3)
            ww  = q(j,k,l,4)
         
            uvw  = 0.5*( uu*uu + vv*vv + ww*ww )
            cjkl = sqrt(ggm1*( q(j,k,l,5) - uvw ))
            c2i  = 1.0/(cjkl*cjkl)
            p(j,k,l) = ( q(j,k,l,5) - uvw )*gm1*q(j,k,l,1)
         
            rj1 = xx(j,k,l)
            rj2 = xy(j,k,l)
            rj3 = xz(j,k,l)
            qq1 = rj1*( uu - ug(j,k,l) ) + rj2*( vv - vg(j,k,l) )
     c                                   + rj3*( ww - wg(j,k,l) )
            rr1 = sqrt( rj1*rj1 + rj2*rj2 + rj3*rj3 )
            rj1 = rj1/rr1
            rj2 = rj2/rr1
            rj3 = rj3/rr1
         
            rk1 = yx(j,k,l)
            rk2 = yy(j,k,l)
            rk3 = yz(j,k,l)
            qq2 = rk1*( uu - ug(j,k,l) ) + rk2*( vv - vg(j,k,l) )
     c                                   + rk3*( ww - wg(j,k,l) )
            rr2 = sqrt( rk1*rk1 + rk2*rk2 + rk3*rk3 )
            rk1 = rk1/rr2
            rk2 = rk2/rr2
            rk3 = rk3/rr2
         
            rl1 = zx(j,k,l)
            rl2 = zy(j,k,l)
            rl3 = zz(j,k,l)
            qq3 = rl1*( uu - ug(j,k,l) ) + rl2*( vv - vg(j,k,l) )
     c                                   + rl3*( ww - wg(j,k,l) )
            rr3 = sqrt( rl1*rl1 + rl2*rl2 + rl3*rl3 )
            rl1 = rl1/rr3
            rl2 = rl2/rr3
            rl3 = rl3/rr3
            
            vnu = (rmue+turmu(j,k,l))/(rey*q(j,k,l,1))
 
            uv(j,k,l) = uvw
            qx(j,k,l) = qq1
            cx(j,k,l) = cjkl*rr1
            qy(j,k,l) = qq2
            cy(j,k,l) = cjkl*rr2
            qz(j,k,l) = qq3
            cz(j,k,l) = cjkl*rr3
            vn(j,k,l) = visc_k*vnu
          
            bSq = Mp**2/(bt(j,k,l)-Mp**2*(bt(j,k,l)-1))
        
            X = sqrt(((1.0 - bSq)*qq2)*((1.0 - bSq)*qq2)
     c         +4.*bSq*cy(j,k,l)*cy(j,k,l))/rr2
            Y = 0.5*((1.0 - bSq)*qq2/rr2 + X)
            Z = 0.5*((1.0 - bSq)*qq2/rr2 - X)
         
            a1 = s(j,k,l,2)*uu + s(j,k,l,3)*vv + ww*s(j,k,l,4)
     c                         - s(j,k,l,5)
            a1 = a1*gm1*c2i
            a1 = a1 + s(j,k,l,1)*(1. - uvw*GM1*c2i)
         
            a2 = (rk2*ww-rk3*vv)*s(j,k,l,1) + rk3*s(j,k,l,3)
     c                                      - rk2*s(j,k,l,4)
            a3 = (rk3*uu-rk1*ww)*s(j,k,l,1) + rk1*s(j,k,l,4)
     c                                      - rk3*s(j,k,l,2)
            a4 = (rk1*vv-rk2*uu)*s(j,k,l,1) + rk2*s(j,k,l,2)
     c                                      - rk1*s(j,k,l,3)
            a5 = uvw*s(j,k,l,1) - s(j,k,l,2)*uu - s(j,k,l,3)*vv
     c                             - s(j,k,l,4)*ww + s(j,k,l,5)
            a5 = a5*gm1*c2i
         
            a6 = qq2*s(j,k,l,1)/rr2 - rk1*s(j,k,l,2)
     c          -rk2*s(j,k,l,3)     - rk3*s(j,k,l,4)
            a6 = a6*bSq
         
            s(j,k,l,1) = a1*rk1 + a2/q(j,k,l,1)
            s(j,k,l,2) = a1*rk2 + a3/q(j,k,l,1)
            s(j,k,l,3) = a1*rk3 + a4/q(j,k,l,1)
            s(j,k,l,4) = -(Z*a5+a6)/X
            s(j,k,l,5) =  (Y*a5+a6)/X
         enddo
      enddo
      enddo
         
      do l = ls-1,le+1
      do k = ks,ke
         do j = js-1,je+1
            d2py(j,k,l) = abs( p(j,k+1,l) - 2.*p(j,k,l) + p(j,k-1,l) ) 
            d2py(j,k,l) = d2py(j,k,l)/abs( p(j,k+1,l) + 2.*p(j,k,l) 
     c                                    + p(j,k-1,l) )
         enddo
      enddo
      enddo

      k = ks-1
      do l = ls-1,le+1
         do j = js-1,je+1
            d2py(j,k,l) = 0. 
         enddo
      enddo

      k = ke+1
      do l = ls-1,le+1
         do j = js-1,je+1
            d2py(j,k,l) = 0. 
         enddo
      enddo

c   K direction
      if (ilhs.eq.2) then
        call lhsinv3d(q,s,qy,cy,d2py,ks-1,ke+1,js-1,je+1,ls-1,le+1,2,
     c                              yx,yy,yz,vn,iblank,tscale,bt)
      elseif (ilhs.eq.3) then
        call lhsinv3d_up(q,s,qy,cy,d2py,ks-1,ke+1,js-1,je+1,ls-1,le+1,2,
     c                              yx,yy,yz,vn,iblank,tscale,bt)
      endif

c multiply by N^inv = T_si^(-1) T_eta 

      do l=ls-1,le+1
      do k=ks-1,ke+1
         do j=js-1,je+1
            rj1 = xx(j,k,l)
            rj2 = xy(j,k,l)
            rj3 = xz(j,k,l)
            rr1 = sqrt( rj1*rj1 + rj2*rj2 + rj3*rj3 )
            rj1 = rj1/rr1
            rj2 = rj2/rr1
            rj3 = rj3/rr1
        
            rk1 = yx(j,k,l)
            rk2 = yy(j,k,l)
            rk3 = yz(j,k,l)
            rr2 = sqrt( rk1*rk1 + rk2*rk2 + rk3*rk3 )
            rk1 = rk1/rr2
            rk2 = rk2/rr2
            rk3 = rk3/rr2

            vnu = (rmue+turmu(j,k,l))/(rey*q(j,k,l,1))
            vn(j,k,l) = visc_j*vnu

            bSq = Mp**2/(bt(j,k,l)-Mp**2*(bt(j,k,l)-1))
        
            X = sqrt(((1.0 - bSq)*qy(j,k,l))*((1.0 - bSq)*qy(j,k,l))
     c          +4.*bSq*cy(j,k,l)*cy(j,k,l))/rr2
            Y = 0.5*((1.0 - bSq)*qy(j,k,l)/rr2 + X)
            Z = 0.5*((1.0 - bSq)*qy(j,k,l)/rr2 - X)

            X1 = sqrt(((1.0 - bSq)*qx(j,k,l))*((1.0 - bSq)*qx(j,k,l))
     c         +4.*bSq*cx(j,k,l)*cx(j,k,l))/rr1
            Y1 = 0.5*((1.0 - bSq)*qx(j,k,l)/rr1 + X1)
            Z1 = 0.5*((1.0 - bSq)*qx(j,k,l)/rr1 - X1)

            a1 = ( rj1*rk1 + rj2*rk2 + rj3*rk3 )
            a2 = ( rj1*rk2 - rk1*rj2 )
            a3 = ( rj3*rk2 - rk3*rj2 )
            a4 = ( rj1*rk3 - rk1*rj3 )
            a5 = ( Y*s(j,k,l,4) + Z*s(j,k,l,5) )/(bSq*q(j,k,l,1))
            a6 = ( s(j,k,l,4) + s(j,k,l,5) )
        
            tmp(1) =   a1*s(j,k,l,1) + a2*s(j,k,l,2) + a4*s(j,k,l,3)
     c               + a5*a3
            tmp(2) = - a2*s(j,k,l,1) + a1*s(j,k,l,2) - a3*s(j,k,l,3)
     c               + a5*a4
            tmp(3) = - a4*s(j,k,l,1) + a3*s(j,k,l,2) + a1*s(j,k,l,3)
     c               - a5*a2
            tmp(4) = (- bSq*q(j,k,l,1)*(a3*s(j,k,l,1) + a4*s(j,k,l,2)
     c                - a2*s(j,k,l,3) - a5*a1)-Z1*a6)/X1
            tmp(5) = (  bSq*q(j,k,l,1)*(a3*s(j,k,l,1) + a4*s(j,k,l,2)
     c                - a2*s(j,k,l,3) - a5*a1)+Y1*a6)/X1

            do n = 1,5
               s(j,k,l,n) = tmp(n) 
            enddo
         enddo
      enddo
      enddo

      do l = ls-1,le+1
         do k = ks-1,ke+1
            do j = js,je
               d2px(j,k,l) =abs(p(j+1,k,l) - 2.*p(j,k,l) + p(j-1,k,l))
               d2px(j,k,l) = d2px(j,k,l)/abs( p(j+1,k,l)
     c                                     + 2.*p(j,k,l) + p(j-1,k,l))
            enddo
         enddo
      enddo

      j=js-1
      do l=ls-1,le+1
         do k=ks-1,ke+1
            d2px(j,k,l) = 0. 
         enddo
      enddo

      j=je+1
      do l=ls-1,le+1
         do k=ks-1,ke+1
            d2px(j,k,l) = 0. 
         enddo
      enddo

c  J direction

      if (ilhs.eq.2) then
        call lhsinv3d(q,s,qx,cx,d2px,js-1,je+1,ks-1,ke+1,ls-1,le+1,1,
     c                       xx,xy,xz,vn,iblank,tscale,bt)
      elseif (ilhs.eq.3) then
        call lhsinv3d_up(q,s,qx,cx,d2px,js-1,je+1,ks-1,ke+1,ls-1,le+1,1,
     c                       xx,xy,xz,vn,iblank,tscale,bt)
      endif
 
c  multiply by N^inv = T_zeta^(-1) T_si

      do l=ls-1,le+1
      do k=ks-1,ke+1
         do j=js-1,je+1

            rj1 = xx(j,k,l)
            rj2 = xy(j,k,l)
            rj3 = xz(j,k,l)
            rr1 = sqrt( rj1*rj1 + rj2*rj2 + rj3*rj3 )
            rj1 = rj1/rr1
            rj2 = rj2/rr1
            rj3 = rj3/rr1
        
            rl1 = zx(j,k,l)
            rl2 = zy(j,k,l)
            rl3 = zz(j,k,l)
            rr3 = sqrt( rl1*rl1 + rl2*rl2 + rl3*rl3 )
            rl1 = rl1/rr3
            rl2 = rl2/rr3
            rl3 = rl3/rr3
        
            vnu = (rmue+turmu(j,k,l))/(rey*q(j,k,l,1))
            vn(j,k,l) = visc_l*vnu

            bSq = Mp**2/(bt(j,k,l)-Mp**2*(bt(j,k,l)-1))
        
            X = sqrt(((1.0 - bSq)*qx(j,k,l))*((1.0 - bSq)*qx(j,k,l))
     c          +4.*bSq*cx(j,k,l)*cx(j,k,l))/rr1
            Y = 0.5*((1.0 - bSq)*qx(j,k,l)/rr1 + X)
            Z = 0.5*((1.0 - bSq)*qx(j,k,l)/rr1 - X)
        
            X1 = sqrt(((1.0 - bSq)*qz(j,k,l))*((1.0 - bSq)*qz(j,k,l))
     c          +4.*bSq*cz(j,k,l)*cz(j,k,l))/rr3
            Y1 = 0.5*((1.0 - bSq)*qz(j,k,l)/rr3 + X1)
            Z1 = 0.5*((1.0 - bSq)*qz(j,k,l)/rr3 - X1)
        
            a1 = (rj1*rl1+rj2*rl2+rj3*rl3)
            a2 = (rl1*rj2-rj1*rl2)
            a3 = (rl3*rj2-rj3*rl2)
            a4 = (rl1*rj3-rj1*rl3)
            a5 = ( Y*s(j,k,l,4) + Z*s(j,k,l,5) )/(bSq*q(j,k,l,1))
            a6 = ( s(j,k,l,4) + s(j,k,l,5) )
        
            tmp(1) =   a1*s(j,k,l,1) + a2*s(j,k,l,2) + a4*s(j,k,l,3)
     c               + a5*a3
            tmp(2) = - a2*s(j,k,l,1) + a1*s(j,k,l,2) - a3*s(j,k,l,3)
     c               + a5*a4
            tmp(3) = - a4*s(j,k,l,1) + a3*s(j,k,l,2) + a1*s(j,k,l,3)
     c               - a5*a2
            tmp(4) = (- bSq*q(j,k,l,1)*(a3*s(j,k,l,1) + a4*s(j,k,l,2)
     c                - a2*s(j,k,l,3) - a5*a1)-Z1*a6)/X1
            tmp(5) = (  bSq*q(j,k,l,1)*(a3*s(j,k,l,1) + a4*s(j,k,l,2)
     c                - a2*s(j,k,l,3) - a5*a1)+Y1*a6)/X1

            do n = 1,5
               s(j,k,l,n) = tmp(n) 
            enddo
         enddo
      enddo
      enddo

      do l = ls,le
      do k = ks-1,ke+1
         do j = js-1,je+1
            d2pz(j,k,l) = abs( p(j,k,l+1) - 2.*p(j,k,l) + p(j,k,l-1) ) 
            d2pz(j,k,l) = d2pz(j,k,l)/abs( p(j,k,l+1) + 2.*p(j,k,l) 
     c                                   + p(j,k,l-1) )
         enddo
      enddo
      enddo

      l = ls-1
      do k = ks-1,ke+1
         do j = js-1,je+1
            d2pz(j,k,l) = 0. 
         enddo
      enddo

      l = le+1
      do k = ks-1,ke+1
         do j = js-1,je+1
            d2pz(j,k,l) = 0. 
         enddo
      enddo
c   L direction
      if (ilhs.eq.2) then
        call lhsinv3d(q,s,qz,cz,d2pz,ls-1,le+1,js-1,je+1,ks-1,ke+1,3,
     c                             zx,zy,zz,vn,iblank,tscale,bt)
      elseif (ilhs.eq.3) then
        call lhsinv3d_up(q,s,qz,cz,d2pz,ls-1,le+1,js-1,je+1,ks-1,ke+1,3,
     c                             zx,zy,zz,vn,iblank,tscale,bt)
      endif
 
c  multiply by T_zeta 

      do l=ls-1,le+1
      do k=ks-1,ke+1
         do j=js-1,je+1
            uu  = q(j,k,l,2)
            vv  = q(j,k,l,3)
            ww  = q(j,k,l,4)
         
            uvw  = 0.5*( uu*uu + vv*vv + ww*ww )
            cjkl = sqrt(ggm1*( q(j,k,l,5) - uvw ))
            c2i  = 1.0/(cjkl*cjkl)
            ge    = gamma*q(j,k,l,5)- gm1*uvw
            qq3  = qz(j,k,l) 
        
            rl1 = zx(j,k,l)
            rl2 = zy(j,k,l)
            rl3 = zz(j,k,l)
            rr3 = sqrt( rl1*rl1 + rl2*rl2 + rl3*rl3 )
            rl1 = rl1/rr3
            rl2 = rl2/rr3
            rl3 = rl3/rr3
        
            bSq = Mp**2/(bt(j,k,l)-Mp**2*(bt(j,k,l)-1))
        
            X = sqrt(((1.0 - bSq)*qq3)*((1.0 - bSq)*qq3)
     c          +4.*bSq*cz(j,k,l)*cz(j,k,l))/rr3
            Y = 0.5*((1.0 - bSq)*qz(j,k,l)/rr3 + X)
            Z = 0.5*((1.0 - bSq)*qz(j,k,l)/rr3 - X)
        
            a1 = rl1*s(j,k,l,1) + rl2*s(j,k,l,2) + rl3*s(j,k,l,3)
     c          +s(j,k,l,4) + s(j,k,l,5)
            a2 = (Y*s(j,k,l,4) + Z*s(j,k,l,5))/bSq
            a3 = uvw*(rl1*s(j,k,l,1) + rl2*s(j,k,l,2) + rl3*s(j,k,l,3))
         
            tmp(1) = a1
            tmp(2) = a1*uu-q(j,k,l,1)*(rl3*s(j,k,l,2) - rl2*s(j,k,l,3))
            tmp(2) = tmp(2) + a2*rl1
            tmp(3) = a1*vv+q(j,k,l,1)*(rl3*s(j,k,l,1) - rl1*s(j,k,l,3))
            tmp(3) = tmp(3) + a2*rl2
            tmp(4) = a1*ww-q(j,k,l,1)*(rl2*s(j,k,l,1) - rl1*s(j,k,l,2))
            tmp(4) = tmp(4) + a2*rl3;
            tmp(5) = a3 + q(j,k,l,1)*((vv*rl3 - ww*rl2)*s(j,k,l,1)
     c                  + (ww*rl1 - uu*rl3)*s(j,k,l,2)
     c                  + (uu*rl2 - vv*rl1)*s(j,k,l,3))
            tmp(5) = tmp(5) + (ge+Y*qq3/(rr3*bSq))*s(j,k,l,4) 
     c                      + (ge+Z*qq3/(rr3*bSq))*s(j,k,l,5)

            do n = 1,5
               s(j,k,l,n) = tmp(n) 
            enddo
         enddo
      enddo
      enddo

      if (iunst.eq.2) then
         do l=ls-1,le+1
         do k=ks-1,ke+1
            do j=js-1,je+1
             r2 = s(j,k,l,2)
             r3 = s(j,k,l,3)
             abcd = 1.0 + (tscale(j,k,l)*rf)**2
             s(j,k,l,2) = ( r2 + tscale(j,k,l)*rf*r3)/abcd
             s(j,k,l,3) = ( r3 - tscale(j,k,l)*rf*r2)/abcd
            enddo
         enddo
         enddo
      endif
c
c..restore density , u , v , w in conservative variables             
c                                                                       
      do 4 l = ls-1, le+1
      do 4 k = ks-1, ke+1
      do 4 j = js-1, je+1
        q(j,k,l,1) = q(j,k,l,1) / q(j,k,l,6)
        q(j,k,l,2) = q(j,k,l,2)*q(j,k,l,1)
        q(j,k,l,3) = q(j,k,l,3)*q(j,k,l,1)
        q(j,k,l,4) = q(j,k,l,4)*q(j,k,l,1)
        q(j,k,l,5) = q(j,k,l,5)*q(j,k,l,1)
    4 continue

      return
      end

c***********************************************************************
      subroutine lhsinv3d(q,s,qn,cn,d2p,ms,me,ns,ne,ls,le,idir,
     c                                    xn,yn,zn,vnu,iblank,tscale,bt)

      use params_global
      implicit none

      integer ms,me,ns,ne,ls,le,idir
      real q(jmax,kmax,lmax,nd),s(jmax,kmax,lmax,nv)
      real qn(jmax,kmax,lmax),cn(jmax,kmax,lmax),d2p(jmax,kmax,lmax)
      real xn(jmax,kmax,lmax),yn(jmax,kmax,lmax),zn(jmax,kmax,lmax)
      real vnu(jmax,kmax,lmax),tscale(jmax,kmax,lmax),bt(jmax,kmax,lmax)
      integer iblank(jmax,kmax,lmax)

      !local variables
      integer i,m,n,l
      real oat,dis2,dis4,X,bSq
      real c2,c2m,c4m2,c4m,c4,c4p
      real,allocatable :: jcbn(:),tscal(:)
      real,allocatable :: vistrm1(:),vistrm2(:),vistrm3(:),g(:)
      real,allocatable :: a(:),b(:),c(:),d(:),e(:),f(:)
      real,allocatable :: a1(:),b1(:),c1(:),d1(:),e1(:)
      real,allocatable :: diag(:,:)
      integer,allocatable :: iblnk(:)

      allocate(jcbn(me),tscal(me))
      allocate(vistrm1(me),vistrm2(me),vistrm3(me),g(me))
      allocate(a(me),b(me),c(me),d(me),e(me),f(me))
      allocate(a1(me),b1(me),c1(me),d1(me),e1(me))
      allocate(diag(me,nd))
      allocate(iblnk(me))
c
c..set time-accuracy
c
      oat = 1.0
      if (ntac.eq.2 .and. istep.gt.1) oat = 2.d0/3.d0
c
c..some initialization
c
      dis2 = 10.0
      dis4 = 0.1

      do l = ls,le
      do n = ns,ne
         if (idir.eq.1) then

            do m = ms,me
               tscal(m) = tscale(m,n,l)
               jcbn(m)  = q(m,n,l,6)
               g(m)     = d2p(m,n,l)
               iblnk(m) = max(iblank(m,n,l),0)

               bSq = Mp**2/(bt(m,n,l)-Mp**2*(bt(m,n,l)-1))

               X= sqrt(((1.-bSq)*qn(m,n,l))*((1.-bSq)*qn(m,n,l))
     c            +4.*bSq*cn(m,n,l)*cn(m,n,l))
               diag(m,1) = qn(m,n,l)
               diag(m,2) = qn(m,n,l)
               diag(m,3) = qn(m,n,l)
               diag(m,4) = 0.5*((bSq+1.)*qn(m,n,l)+X)
               diag(m,5) = 0.5*((bSq+1.)*qn(m,n,l)-X)
c              spectral radius divided by jacobian
               diag(m,6)= 0.5*((bSq+1.)*abs(qn(m,n,l))+X)/q(m,n,l,6)
            enddo

            do m = ms+1,me-1
               vistrm1(m) = (xn(m,n,l)*xn(m,n,l)+yn(m,n,l)*yn(m,n,l)
     c                      +zn(m,n,l)*zn(m,n,l))/q(m,n,l,6)
               vistrm1(m) = vistrm1(m) + (xn(m+1,n,l)*xn(m+1,n,l)
     c                     + yn(m+1,n,l)*yn(m+1,n,l) 
     c                     + zn(m+1,n,l)*zn(m+1,n,l))/q(m+1,n,l,6)
               vistrm1(m) = vistrm1(m)*q(m,n,l,6)*vnu(m,n,l)

               vistrm3(m) = (xn(m,n,l)*xn(m,n,l)+yn(m,n,l)*yn(m,n,l)
     c                      +zn(m,n,l)*zn(m,n,l))/q(m,n,l,6)
               vistrm3(m) = vistrm3(m) + (xn(m-1,n,l)*xn(m-1,n,l)
     c                     + yn(m-1,n,l)*yn(m-1,n,l)
     c                     + zn(m-1,n,l)*zn(m-1,n,l))/q(m-1,n,l,6)
               vistrm3(m) = vistrm3(m)*q(m,n,l,6)*vnu(m,n,l)

               vistrm2(m) = 0.5*(vistrm1(m) + vistrm3(m))
            enddo
            m = ms
            vistrm1(m) = 0
            vistrm2(m) = 0
            vistrm3(m) = 0
            m = me
            vistrm1(m) = 0
            vistrm2(m) = 0
            vistrm3(m) = 0

         elseif (idir.eq.2) then

            do m = ms,me
               tscal(m) = tscale(n,m,l)
               jcbn(m) = q(n,m,l,6)
               g(m)    = d2p(n,m,l)
               iblnk(m) = max(iblank(n,m,l),0)

               bSq = Mp**2/(bt(n,m,l)-Mp**2*(bt(n,m,l)-1))

               X= sqrt(((1.-bSq)*qn(n,m,l))*((1.-bSq)*qn(n,m,l))
     c            +4.*bSq*cn(n,m,l)*cn(n,m,l))
               diag(m,1) = qn(n,m,l)
               diag(m,2) = qn(n,m,l)
               diag(m,3) = qn(n,m,l)
               diag(m,4) = 0.5*((bSq+1.)*qn(n,m,l)+X)
               diag(m,5) = 0.5*((bSq+1.)*qn(n,m,l)-X)
c              spectral radius divided by jacobian
               diag(m,6)= 0.5*((bSq+1.)*abs(qn(n,m,l))+X)/q(n,m,l,6)
            enddo

            do m = ms+1,me-1
               vistrm1(m) = (xn(n,m,l)*xn(n,m,l)+yn(n,m,l)*yn(n,m,l)
     c                      +zn(n,m,l)*zn(n,m,l))/q(n,m,l,6)
               vistrm1(m) = vistrm1(m) + (xn(n,m+1,l)*xn(n,m+1,l)
     c                     + yn(n,m+1,l)*yn(n,m+1,l) 
     c                     + zn(n,m+1,l)*zn(n,m+1,l))/q(n,m+1,l,6)
               vistrm1(m) = vistrm1(m)*q(n,m,l,6)*vnu(n,m,l)

               vistrm3(m) = (xn(n,m,l)*xn(n,m,l)+yn(n,m,l)*yn(n,m,l)
     c                      +zn(n,m,l)*zn(n,m,l))/q(n,m,l,6)
               vistrm3(m) = vistrm3(m) + (xn(n,m-1,l)*xn(n,m-1,l)
     c                     + yn(n,m-1,l)*yn(n,m-1,l)
     c                     + zn(n,m-1,l)*zn(n,m-1,l))/q(n,m-1,l,6)
               vistrm3(m)= vistrm3(m)*q(n,m,l,6)*vnu(n,m,l)

               vistrm2(m) = 0.5*(vistrm1(m) + vistrm3(m))
            enddo
            m = ms
            vistrm1(m) = 0
            vistrm2(m) = 0
            vistrm3(m) = 0
            m = me
            vistrm1(m) = 0
            vistrm2(m) = 0
            vistrm3(m) = 0

         elseif (idir.eq.3) then

            do m = ms,me
               tscal(m) = tscale(n,l,m)
               jcbn(m) = q(n,l,m,6)
               g(m)    = d2p(n,l,m)
               iblnk(m) = max(iblank(n,l,m),0)

               bSq = Mp**2/(bt(n,l,m)-Mp**2*(bt(n,l,m)-1))

               X= sqrt(((1.-bSq)*qn(n,l,m))*((1.-bSq)*qn(n,l,m))
     c            +4.*bSq*cn(n,l,m)*cn(n,l,m))
               diag(m,1) = qn(n,l,m)
               diag(m,2) = qn(n,l,m)
               diag(m,3) = qn(n,l,m)
               diag(m,4) = 0.5*((bSq+1.)*qn(n,l,m)+X)
               diag(m,5) = 0.5*((bSq+1.)*qn(n,l,m)-X)
c              spectral radius divided by jacobian
               diag(m,6)= 0.5*((bSq+1.)*abs(qn(n,l,m))+X)/q(n,l,m,6)
            enddo

            do m = ms+1,me-1
               vistrm1(m) = (xn(n,l,m)*xn(n,l,m)+yn(n,l,m)*yn(n,l,m)
     c                      +zn(n,l,m)*zn(n,l,m))/q(n,l,m,6)
               vistrm1(m) = vistrm1(m) + (xn(n,l,m+1)*xn(n,l,m+1)
     c                     + yn(n,l,m+1)*yn(n,l,m+1) 
     c                     + zn(n,l,m+1)*zn(n,l,m+1))/q(n,l,m+1,6)
               vistrm1(m) = vistrm1(m)*q(n,l,m,6)*vnu(n,l,m)

               vistrm3(m) = (xn(n,l,m)*xn(n,l,m)+yn(n,l,m)*yn(n,l,m)
     c                      +zn(n,l,m)*zn(n,l,m))/q(n,l,m,6)
               vistrm3(m) = vistrm3(m) + (xn(n,l,m-1)*xn(n,l,m-1)
     c                     + yn(n,l,m-1)*yn(n,l,m-1)
     c                     + zn(n,l,m-1)*zn(n,l,m-1))/q(n,l,m-1,6)
               vistrm3(m)= vistrm3(m)*q(n,l,m,6)*vnu(n,l,m)

               vistrm2(m) = 0.5*(vistrm1(m) + vistrm3(m))
            enddo
            m = ms
            vistrm1(m) = 0
            vistrm2(m) = 0
            vistrm3(m) = 0
            m = me
            vistrm1(m) = 0
            vistrm2(m) = 0
            vistrm3(m) = 0

         endif

         do m = ms+3,me-3
c           second order dissipation
            c2   = 0.5*((g(m)+g(m+1)))*dis2
            c2m  = 0.5*((g(m)+g(m-1)))*dis2

            c4p  = (dis4-min(dis4,.5*dis2*(g(m+2)+g(m+1))))
            c4   = (dis4-min(dis4,1.*c2))
            c4m  = (dis4-min(dis4,1.*c2m))
            c4m2 = (dis4-min(dis4,.5*dis2*(g(m-2)+g(m-1))))

            a1(m) = (diag(m+1,nd)+diag(m+2,nd))*c4p
            b1(m) = -((diag(m+2,nd)+diag(m+1,nd))*c4p+
     c               (diag(m+1,nd)+diag(m,nd))*(3.*c4+c2))
            c1(m) = +((diag(m-1,nd)+diag(m,nd))*(3.*c4m+c2m)+
     c               (diag(m+1,nd)+diag(m,nd))*(3.*c4+c2))
            d1(m) = -((diag(m-2,nd)+diag(m-1,nd))*c4m2+
     c               (diag(m,nd)+diag(m-1,nd))*(3.*c4m+c2m))
            e1(m) = (diag(m-1,nd)+diag(m-2,nd))*c4m2
         enddo 

         m    = ms+2
         c2   = 0.5*((g(m)+g(m+1)))*dis2
         c2m  = 0.5*((g(m)+g(m-1)))*dis2
         c4p  = (dis4-min(dis4,.5*dis2*(g(m+2)+g(m+1))))
         c4   = (dis4-1.*min(dis4,c2))
         c4m  = (dis4-min(dis4,1.*c2m))
         c4m2 = (dis4-min(dis4,1.*dis2*(g(m-1))))
         a1(m) = (diag(m+1,nd)+diag(m+2,nd))*c4p
         b1(m) = -((diag(m+2,nd)+diag(m+1,nd))*c4p+
     c            (diag(m+1,nd)+diag(m,nd))*(3.*c4+c2))
         c1(m) = ((diag(m-1,nd)+diag(m,nd))*(3.*c4m+c2m)+
     c           (diag(m+1,nd)+diag(m,nd))*(3.*c4+c2))
         d1(m) = -((2.*diag(m-1,nd))*c4m2+
     c            (diag(m,nd)+diag(m-1,nd))*(3.*c4m+c2m))
         
         m    = ms+1
         c2   = 0.5*((g(m)+g(m+1)))*dis2
         c2m  = g(m)*dis2
         c4p  = (dis4-min(dis4,.5*dis2*(g(m+2)+g(m+1))))
         c4   = (dis4-min(dis4,1.*c2))
         c4m  = (dis4-min(dis4,1.*c2m))
         a1(m) = (diag(m+1,nd)+diag(m+2,nd))*c4p
         b1(m) = -((diag(m+2,nd)+diag(m+1,nd))*c4p+
     c            (diag(m+1,nd)+diag(m,nd))*(3.*c4+c2))
         c1(m) = (2.*(diag(m,nd))*(2.*c4m+c2m)+
     c            (diag(m+1,nd)+diag(m,nd))*(3.*c4+c2))
         m    = me-2
         c2   = 0.5*((g(m)+g(m+1)))*dis2
         c2m  = 0.5*((g(m)+g(m-1)))*dis2
         c4p  = (dis4-min(dis4,1.*dis2*(g(m+1))))
         c4   = (dis4-min(dis4,1.*c2))
         c4m  = (dis4-min(dis4,1.*c2m))
         c4m2 = (dis4-min(dis4,.5*dis2*(g(m-2)+g(m-1))))
         b1(m) = -(2.*(diag(m+1,nd))*c4p+
     c            (diag(m+1,nd)+diag(m,nd))*(3.*c4+c2))
         c1(m) = ((diag(m-1,nd)+diag(m,nd))*(3.*c4m+c2m)+
     c           (diag(m+1,nd)+diag(m,nd))*(2.*c4+c2))
         d1(m) = -((diag(m-2,nd)+diag(m-1,nd))*c4m2+
     c            (diag(m,nd)+diag(m-1,nd))*(3.*c4m+c2m))
         e1(m) = (diag(m-1,nd)+diag(m-2,nd))*c4m2
         
         m    = me-1  
         c2   = ((g(m)))*dis2
         c2m  = 0.5*((g(m)+g(m-1)))*dis2
         c4   = (dis4-min(dis4,1.*c2))
         c4m  = (dis4-min(dis4,1.*c2m))
         c4m2 = (dis4-min(dis4,.5*dis2*(g(m-2)+g(m-1))))
         c1(m) = +((diag(m-1,nd)+diag(m,nd))*(3.*c4m+c2m)+
     c            2.*(diag(m,nd))*(2.*c4+c2))
         d1(m) = -((diag(m-2,nd)+diag(m-1,nd))*c4m2+
     c            (diag(m,nd)+diag(m-1,nd))*(3.*c4m+c2m))
         e1(m) = (diag(m-1,nd)+diag(m-2,nd))*c4m2

         do m = ms,me
            a1(m) = a1(m)*jcbn(m)
            b1(m) = b1(m)*jcbn(m)
            c1(m) = c1(m)*jcbn(m)
            d1(m) = d1(m)*jcbn(m)
            e1(m) = e1(m)*jcbn(m)
c         viscous terms
            b1(m) = b1(m) - 0.5*vistrm1(m)
            c1(m) = c1(m) + vistrm2(m)
            d1(m) = d1(m) - 0.5*vistrm3(m)
         enddo

         do i = 1,nd-1
            do m = ms,me
               a(m) = a1(m)
               b(m) = b1(m)
               c(m) = c1(m)
               d(m) = d1(m)
               e(m) = e1(m)
               if (idir.eq.1) then
                  f(m)    = s(m,n,l,i)
               elseif (idir.eq.2) then
                  f(m)    = s(n,m,l,i)
               elseif (idir.eq.3) then
                  f(m)    = s(n,l,m,i)
               endif

c           euler terms

               b(m) = b(m) - diag(m,i)*0.5
               d(m) = d(m) + diag(m,i)*0.5

c         multiply by delta t
               if(m.le.me-2) then
                  a(m) = a(m)*tscal(m+2)*iblnk(m+2)
               else
                  a(m) = a(m)*tscal(m)*iblnk(m)
               endif
               if(m.le.me-1) then
                  b(m) = b(m)*tscal(m+1)*iblnk(m+1)
               else
                  b(m) = b(m)*tscal(m)*iblnk(m)
               endif
               if(m.ge.ms+1) then
                  d(m) = d(m)*tscal(m-1)*iblnk(m-1)
               else
                  d(m) = d(m)*tscal(m)*iblnk(m)
               endif
               if(m.ge.ms+2) then
                  e(m) = e(m)*tscal(m-2)*iblnk(m-2)
               else
                  e(m) = e(m)*tscal(m)*iblnk(m)
               endif
               c(m) = c(m)*tscal(m)*iblnk(m)

c          add identity to matrix
               c(m) = c(m) + 1.0
            enddo

            do m = ms+1,me-1
               a(m-ms) = a(m) 
               b(m-ms) = b(m) 
               c(m-ms) = c(m) 
               d(m-ms) = d(m) 
               e(m-ms) = e(m) 
               f(m-ms) = f(m) 
            enddo
            d(1) = 0
            e(1) = 0
            e(2) = 0
            b(me-ms-1) = 0
            a(me-ms-1) = 0
            a(me-ms-2) = 0

            call pentadag(a,b,c,d,e,f,me-ms-1)

            do m=ms+1,me-1
               if (idir.eq.1) then
                  s(m,n,l,i)=f(m-ms)
               elseif (idir.eq.2) then
                  s(n,m,l,i)=f(m-ms)
               elseif (idir.eq.3) then
                  s(n,l,m,i)=f(m-ms)
               endif
            enddo
         enddo
      enddo 
      enddo

      return
      end

c***********************************************************************
      subroutine lhsinv3d_up(q,s,qn,cn,d2p,ms,me,ns,ne,ls,le,idir,
     c                                    xn,yn,zn,vnu,iblank,tscale,bt)

      use params_global
      implicit none

      integer ms,me,ns,ne,ls,le,idir
      real q(jmax,kmax,lmax,nd),s(jmax,kmax,lmax,nv)
      real qn(jmax,kmax,lmax),cn(jmax,kmax,lmax),d2p(jmax,kmax,lmax)
      real xn(jmax,kmax,lmax),yn(jmax,kmax,lmax),zn(jmax,kmax,lmax)
      real vnu(jmax,kmax,lmax),tscale(jmax,kmax,lmax),bt(jmax,kmax,lmax)
      integer iblank(jmax,kmax,lmax)

      !local variables
      integer i,m,n,l
      real oat,dis2,dis4,X,bSq
      real c2,c2m,c4m2,c4m,c4,c4p
      real,allocatable :: jcbn(:),tscal(:)
      real,allocatable :: vistrm1(:),vistrm2(:),vistrm3(:),g(:)
      real,allocatable :: a(:),b(:),c(:),d(:),e(:),f(:)
      real,allocatable :: a1(:),b1(:),c1(:),d1(:),e1(:)
      real,allocatable :: diag_plus(:,:),diag_minus(:,:)
      integer,allocatable :: iblnk(:)
      real eig1,eig2,eig3,epsval,eps,fac

      allocate(jcbn(me),tscal(me))
      allocate(vistrm1(me),vistrm2(me),vistrm3(me),g(me))
      allocate(a(me),b(me),c(me),d(me),e(me),f(me))
      allocate(a1(me),b1(me),c1(me),d1(me),e1(me))
      allocate(diag_plus(me,nv),diag_minus(me,nv))
      allocate(iblnk(me))
c
c..set time-accuracy
c
      oat = 1.0
      if (ntac.eq.2 .and. istep.gt.1) oat = 2.d0/3.d0
c
c..some initialization
c
      epsval = 0.08
      fac    = 1.05

      do l = ls,le
      do n = ns,ne
         if (idir.eq.1) then

            do m = ms,me
               tscal(m) = tscale(m,n,l)
               jcbn(m) = q(m,n,l,6)
               g(m)    = d2p(m,n,l)
               iblnk(m) = max(iblank(m,n,l),0)

               bSq = Mp**2/(bt(m,n,l)-Mp**2*(bt(m,n,l)-1))

               X= sqrt(((1.-bSq)*qn(m,n,l))*((1.-bSq)*qn(m,n,l))
     c            +4.*bSq*cn(m,n,l)*cn(m,n,l))
               eig1 = qn(m,n,l)
               eig2 = 0.5*((bSq+1.)*qn(m,n,l)+X)
               eig3 = 0.5*((bSq+1.)*qn(m,n,l)-X)
               eps  = epsval*sqrt( xn(m,n,l)**2 + yn(m,n,l)**2 
     c                                          + zn(m,n,l)**2 )
               diag_plus(m,1)  = 0.5*(eig1 + fac*sqrt(eig1**2 + eps**2))
               diag_plus(m,2)  = 0.5*(eig1 + fac*sqrt(eig1**2 + eps**2))
               diag_plus(m,3)  = 0.5*(eig1 + fac*sqrt(eig1**2 + eps**2))
               diag_plus(m,4)  = 0.5*(eig2 + fac*sqrt(eig2**2 + eps**2))
               diag_plus(m,5)  = 0.5*(eig3 + fac*sqrt(eig3**2 + eps**2))
               diag_minus(m,1) = 0.5*(eig1 - fac*sqrt(eig1**2 + eps**2))
               diag_minus(m,2) = 0.5*(eig1 - fac*sqrt(eig1**2 + eps**2))
               diag_minus(m,3) = 0.5*(eig1 - fac*sqrt(eig1**2 + eps**2))
               diag_minus(m,4) = 0.5*(eig2 - fac*sqrt(eig2**2 + eps**2))
               diag_minus(m,5) = 0.5*(eig3 - fac*sqrt(eig3**2 + eps**2))
            enddo

            do m = ms+1,me-1
               vistrm1(m) = (xn(m,n,l)*xn(m,n,l)+yn(m,n,l)*yn(m,n,l)
     c                      +zn(m,n,l)*zn(m,n,l))/q(m,n,l,6)
               vistrm1(m) = vistrm1(m) + (xn(m+1,n,l)*xn(m+1,n,l)
     c                     + yn(m+1,n,l)*yn(m+1,n,l) 
     c                     + zn(m+1,n,l)*zn(m+1,n,l))/q(m+1,n,l,6)
               vistrm1(m) = vistrm1(m)*q(m,n,l,6)*vnu(m,n,l)

               vistrm3(m) = (xn(m,n,l)*xn(m,n,l)+yn(m,n,l)*yn(m,n,l)
     c                      +zn(m,n,l)*zn(m,n,l))/q(m,n,l,6)
               vistrm3(m) = vistrm3(m) + (xn(m-1,n,l)*xn(m-1,n,l)
     c                     + yn(m-1,n,l)*yn(m-1,n,l)
     c                     + zn(m-1,n,l)*zn(m-1,n,l))/q(m-1,n,l,6)
               vistrm3(m) = vistrm3(m)*q(m,n,l,6)*vnu(m,n,l)

               vistrm2(m) = 0.5*(vistrm1(m) + vistrm3(m))
            enddo
            m = ms
            vistrm1(m) = 0
            vistrm2(m) = 0
            vistrm3(m) = 0
            m = me
            vistrm1(m) = 0
            vistrm2(m) = 0
            vistrm3(m) = 0

         elseif (idir.eq.2) then

            do m = ms,me
               tscal(m) = tscale(n,m,l)
               jcbn(m) = q(n,m,l,6)
               g(m)    = d2p(n,m,l)
               iblnk(m) = max(iblank(n,m,l),0)

               bSq = Mp**2/(bt(n,m,l)-Mp**2*(bt(n,m,l)-1))

               X= sqrt(((1.-bSq)*qn(n,m,l))*((1.-bSq)*qn(n,m,l))
     c            +4.*bSq*cn(n,m,l)*cn(n,m,l))
               eig1 = qn(n,m,l)
               eig2 = 0.5*((bSq+1.)*qn(n,m,l)+X)
               eig3 = 0.5*((bSq+1.)*qn(n,m,l)-X)
               eps  = epsval*sqrt( xn(n,m,l)**2 + yn(n,m,l)**2 
     c                                          + zn(n,m,l)**2 )
               diag_plus(m,1)  = 0.5*(eig1 + fac*sqrt(eig1**2 + eps**2))
               diag_plus(m,2)  = 0.5*(eig1 + fac*sqrt(eig1**2 + eps**2))
               diag_plus(m,3)  = 0.5*(eig1 + fac*sqrt(eig1**2 + eps**2))
               diag_plus(m,4)  = 0.5*(eig2 + fac*sqrt(eig2**2 + eps**2))
               diag_plus(m,5)  = 0.5*(eig3 + fac*sqrt(eig3**2 + eps**2))
               diag_minus(m,1) = 0.5*(eig1 - fac*sqrt(eig1**2 + eps**2))
               diag_minus(m,2) = 0.5*(eig1 - fac*sqrt(eig1**2 + eps**2))
               diag_minus(m,3) = 0.5*(eig1 - fac*sqrt(eig1**2 + eps**2))
               diag_minus(m,4) = 0.5*(eig2 - fac*sqrt(eig2**2 + eps**2))
               diag_minus(m,5) = 0.5*(eig3 - fac*sqrt(eig3**2 + eps**2))
            enddo

            do m = ms+1,me-1
               vistrm1(m) = (xn(n,m,l)*xn(n,m,l)+yn(n,m,l)*yn(n,m,l)
     c                      +zn(n,m,l)*zn(n,m,l))/q(n,m,l,6)
               vistrm1(m) = vistrm1(m) + (xn(n,m+1,l)*xn(n,m+1,l)
     c                     + yn(n,m+1,l)*yn(n,m+1,l) 
     c                     + zn(n,m+1,l)*zn(n,m+1,l))/q(n,m+1,l,6)
               vistrm1(m) = vistrm1(m)*q(n,m,l,6)*vnu(n,m,l)

               vistrm3(m) = (xn(n,m,l)*xn(n,m,l)+yn(n,m,l)*yn(n,m,l)
     c                      +zn(n,m,l)*zn(n,m,l))/q(n,m,l,6)
               vistrm3(m) = vistrm3(m) + (xn(n,m-1,l)*xn(n,m-1,l)
     c                     + yn(n,m-1,l)*yn(n,m-1,l)
     c                     + zn(n,m-1,l)*zn(n,m-1,l))/q(n,m-1,l,6)
               vistrm3(m)= vistrm3(m)*q(n,m,l,6)*vnu(n,m,l)

               vistrm2(m) = 0.5*(vistrm1(m) + vistrm3(m))
            enddo
            m = ms
            vistrm1(m) = 0
            vistrm2(m) = 0
            vistrm3(m) = 0
            m = me
            vistrm1(m) = 0
            vistrm2(m) = 0
            vistrm3(m) = 0

         elseif (idir.eq.3) then

            do m = ms,me
               tscal(m) = tscale(n,l,m)
               jcbn(m) = q(n,l,m,6)
               g(m)    = d2p(n,l,m)
               iblnk(m) = max(iblank(n,l,m),0)

               bSq = Mp**2/(bt(n,l,m)-Mp**2*(bt(n,l,m)-1))

               X= sqrt(((1.-bSq)*qn(n,l,m))*((1.-bSq)*qn(n,l,m))
     c            +4.*bSq*cn(n,l,m)*cn(n,l,m))
               eig1 = qn(n,l,m)
               eig2 = 0.5*((bSq+1.)*qn(n,l,m)+X)
               eig3 = 0.5*((bSq+1.)*qn(n,l,m)-X)
               eps  = epsval*sqrt( xn(n,l,m)**2 + yn(n,l,m)**2 
     c                                          + zn(n,l,m)**2 )
               diag_plus(m,1)  = 0.5*(eig1 + fac*sqrt(eig1**2 + eps**2))
               diag_plus(m,2)  = 0.5*(eig1 + fac*sqrt(eig1**2 + eps**2))
               diag_plus(m,3)  = 0.5*(eig1 + fac*sqrt(eig1**2 + eps**2))
               diag_plus(m,4)  = 0.5*(eig2 + fac*sqrt(eig2**2 + eps**2))
               diag_plus(m,5)  = 0.5*(eig3 + fac*sqrt(eig3**2 + eps**2))
               diag_minus(m,1) = 0.5*(eig1 - fac*sqrt(eig1**2 + eps**2))
               diag_minus(m,2) = 0.5*(eig1 - fac*sqrt(eig1**2 + eps**2))
               diag_minus(m,3) = 0.5*(eig1 - fac*sqrt(eig1**2 + eps**2))
               diag_minus(m,4) = 0.5*(eig2 - fac*sqrt(eig2**2 + eps**2))
               diag_minus(m,5) = 0.5*(eig3 - fac*sqrt(eig3**2 + eps**2))
            enddo

            do m = ms+1,me-1
               vistrm1(m) = (xn(n,l,m)*xn(n,l,m)+yn(n,l,m)*yn(n,l,m)
     c                      +zn(n,l,m)*zn(n,l,m))/q(n,l,m,6)
               vistrm1(m) = vistrm1(m) + (xn(n,l,m+1)*xn(n,l,m+1)
     c                     + yn(n,l,m+1)*yn(n,l,m+1) 
     c                     + zn(n,l,m+1)*zn(n,l,m+1))/q(n,l,m+1,6)
               vistrm1(m) = vistrm1(m)*q(n,l,m,6)*vnu(n,l,m)

               vistrm3(m) = (xn(n,l,m)*xn(n,l,m)+yn(n,l,m)*yn(n,l,m)
     c                      +zn(n,l,m)*zn(n,l,m))/q(n,l,m,6)
               vistrm3(m) = vistrm3(m) + (xn(n,l,m-1)*xn(n,l,m-1)
     c                     + yn(n,l,m-1)*yn(n,l,m-1)
     c                     + zn(n,l,m-1)*zn(n,l,m-1))/q(n,l,m-1,6)
               vistrm3(m)= vistrm3(m)*q(n,l,m,6)*vnu(n,l,m)

               vistrm2(m) = 0.5*(vistrm1(m) + vistrm3(m))
            enddo
            m = ms
            vistrm1(m) = 0
            vistrm2(m) = 0
            vistrm3(m) = 0
            m = me
            vistrm1(m) = 0
            vistrm2(m) = 0
            vistrm3(m) = 0

         endif

         do m = ms,me
            a1(m) = 0
            b1(m) = 0
            c1(m) = 0
            d1(m) = 0
            e1(m) = 0

c         viscous terms
            b1(m) = b1(m) - 0.5*vistrm1(m)
            c1(m) = c1(m) + vistrm2(m)    
            d1(m) = d1(m) - 0.5*vistrm3(m)
         enddo

         do i = 1,nd-1
            do m = ms,me
               a(m) = a1(m)
               b(m) = b1(m)
               c(m) = c1(m)
               d(m) = d1(m)
               e(m) = e1(m)
               if (idir.eq.1) then
                  f(m)    = s(m,n,l,i)
               elseif (idir.eq.2) then
                  f(m)    = s(n,m,l,i)
               elseif (idir.eq.3) then
                  f(m)    = s(n,l,m,i)
               endif

c           euler terms

               b(m) = b(m) - diag_plus(m,i)
               d(m) = d(m) + diag_minus(m,i)
               c(m) = c(m) + diag_plus(m,i) - diag_minus(m,i)


c         multiply by delta t
               if(m.le.me-2) then
                  a(m) = a(m)*tscal(m+2)*iblnk(m+2)
               else
                  a(m) = a(m)*tscal(m)*iblnk(m)
               endif
               if(m.le.me-1) then
                  b(m) = b(m)*tscal(m+1)*iblnk(m+1)
               else
                  b(m) = b(m)*tscal(m)*iblnk(m)
               endif
               if(m.ge.ms+1) then
                  d(m) = d(m)*tscal(m-1)*iblnk(m-1)
               else
                  d(m) = d(m)*tscal(m)*iblnk(m)
               endif
               if(m.ge.ms+2) then
                  e(m) = e(m)*tscal(m-2)*iblnk(m-2)
               else
                  e(m) = e(m)*tscal(m)*iblnk(m)
               endif
                 c(m) = c(m)*tscal(m)*iblnk(m)

c          add identity to matrix
               c(m) = c(m) + 1.0

            enddo

            do m = ms+1,me-1
               a(m-ms) = a(m) 
               b(m-ms) = b(m) 
               c(m-ms) = c(m) 
               d(m-ms) = d(m) 
               e(m-ms) = e(m) 
               f(m-ms) = f(m) 
            enddo
            d(1) = 0
            e(1) = 0
            e(2) = 0
            b(me-ms-1) = 0
            a(me-ms-1) = 0
            a(me-ms-2) = 0

            call pentadag(a,b,c,d,e,f,me-ms-1)

            do m=ms+1,me-1
               if (idir.eq.1) then
                  s(m,n,l,i)=f(m-ms)
               elseif (idir.eq.2) then
                  s(n,m,l,i)=f(m-ms)
               elseif (idir.eq.3) then
                  s(n,l,m,i)=f(m-ms)
               endif
            enddo
         enddo
      enddo 
      enddo

      return
      end
c**********************************************************************
      subroutine pentadag(a,b,c,d,e,f,N)

c**********************************************************************
      use params_global
      implicit none

      integer N
      real a(N),b(N),c(N),d(N),e(N),f(N)
 
      !local variables
      integer i
      real d0,d1,d2

c     making diagonal 1.
      i = 1
      d0 = 1.0/c(i)
      d(i+1) = d(i+1)*d0
      e(i+2) = e(i+2)*d0
      f(i)   = f(i)*d0
      
c     making diagonal 1. and lower diagonal 0.
      i = 2
      d1 = b(i-1)
      d0 = 1.0/(c(i)-d1*d(i))
      f(i)   = (f(i)-d1*f(i-1))*d0
      d(i+1) = (d(i+1)-d1*e(i+1))*d0
      e(i+2) = e(i+2)*d0
      
c     making diaongal 1. and lower diagonals 0.
      do i = 3,N
         d1 = a(i-2)
         d2 = (b(i-1) - a(i-2)*d(i-1))
         d0 = 1./(c(i)-d2*d(i)-d1*e(i))
         f(i) = (f(i)-d2*f(i-1)-d1*f(i-2))*d0
         if (i.le.N-1) then
            d(i+1) = (d(i+1)-d2*e(i+1))*d0
            if (i.le.N-2) then
               e(i+2) = (e(i+2))*d0
            endif
         endif
      enddo 

c     backward sweep
      i  = N-1
      d1 = 0
      d2 = 0
      do i = N-1,1,-1
         d1 = d(i+1)
         if (i.lt.N-1) then
            d2 = e(i+2)
         endif
         f(i) = f(i) - d1*f(i+1) - d2*f(i+2)
      enddo

      return
      end

