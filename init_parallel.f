c***********************************************************************
      subroutine setup

c     setup the global physical variables
c     freestream direction, reynolds number
c     non-dimensional values for pressure, density and energy

      use params_global
      ! Math: add these to define zero and two (used if fmtip==0)
      use constants

      ! Math: add implicit none (check that zero and two are defined)
      implicit none

      real :: cs, ss, psi

c..   calculate a few things based on inputs

      dblad=2*pi/nblade
     
!     dblad=2*pi/numprocs
      
      theta0=theta0*pi/180.
      theta1c=theta1c*pi/180.
      theta1s=theta1s*pi/180.
            
      angmax = pi*angmax/180.
      if (invisc) rmue = 0.0
      
      blang = 2.*pi/nblade

      ! Math: modify rf definition
      if (fmtip > zero) then
         rf=fmtip/rartio
      else
         rf=two*fsmach
      end if

      if(dt.lt.0.0 .and. iunst .ne. FL_STEADY .and. iunst .ne. FL_PITCH 
     &             .and. iunst .ne. FL_FLAP) dt = abs(dt)*pi/abs(rf)/180.
      ! Math: add kinematics
      ! Updated by Camli 9/2014 

      if(dt.lt.0.0 .and. iunst .eq. FL_PITCH) dt = abs(dt)*pi/omega/180.
      if(dt.lt.0.0 .and. iunst .eq. FL_FLAP)  then
         if (fmtip > 0) then
         	rf = 2*kin_rf*fmtip
         else
            rf = 2*kin_rf*fsmach
         end if
         dt = abs(dt)*pi/abs(rf)/180. !abs(kin_dt)
      endif

c.. addad by camli July/31
      kin_xoff = xoff

      !if(dt.lt.0.0 .and. iunst .eq. FL_FLAP)  dt = abs(dt)*pi/abs(rf)/180.
      print*,'DT in init_parallel = ',dt
      ! Math: irot_dir hasn't been read yet so this is useless!!!
      ! we do rf=rf*irot_dir when reading mesh_azi.inp later so it should be fine
      !if(irot_dir.eq.-1) xoff = -xoff
      !if(irot_dir.eq.0) rf = 0.0

c..   verify inputs

      !call check_inputs

c..   calculate some more things

      h = dt
      hd = .5*dt

      gm1   = gamma -1.0
      ggm1  = gamma*gm1
      cs    = cos( pi*alfa/180.0 )
      ss    = sin( pi*alfa/180.0 )
      alf   = alfa
      psi   = rf*totime
      
      einf  = 1.0/ggm1 +0.5*fsmach**2
      pinf  = 1.0/gamma
      rinf  = 1.0
      ainf  = gamma*pinf/rinf
      htinf = gamma*einf/rinf - 0.5*gm1*fsmach**2

      ! Math: add kinematics
      if(iunst.eq. FL_STEADY .or. iunst.eq. FL_PITCH) then
        uinf = fsmach*cs
        vinf = uinf*tan(beta1*pi/180.)
      else
        uinf  = 0.0
        vinf  = fsmach*cs
      endif

      winf  = fsmach*ss
      rey   = rey/(fsmach + fmtip)

      return
      end

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! --Xuan 03/30/2018
! --Calculate wing twist
     do j = 1, n_Steps*nTimeSteps/2, int(nTimeSteps/2)
       do i = 1, nBladeElements
          wing_twist(d,i) = twisting*( ( (i/real(nBladeElements)- 1.0/(2*nBladeElements))
      &                                  *temp_twist_interp(j)) +temp_const_interp(j)) *pi/180)
 
       end do
       d = d+1
     end do
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


c*********************************************************************
      subroutine qzero(q,vnu)
c
c  set initial values of the flow variables to their free-stream
c  values for an impulsive start. otherwise the initial values
c  are read from a restart file.
c
c*********************************************************************
      use params_global

      implicit none

      real q(jmax,kmax,lmax,nd),vnu(jmax,kmax,lmax)
      real ruinf,rvinf,rwinf
      integer j,k,l
      
      ruinf = rinf*uinf
      rvinf = rinf*vinf
      rwinf = rinf*winf

      do 11 l = 1,lmax
         do 11 k = 1,kmax
            do 11 j = 1,jmax
               q(j,k,l,1) = rinf
               q(j,k,l,2) = ruinf
               q(j,k,l,3) = rvinf
               q(j,k,l,4) = rwinf
               q(j,k,l,5) = einf
               vnu(j,k,l) = 0.1 !standard practise for starting up
 11         continue

      return
      end

      subroutine findpretwist(x,y,z,twist)

      use ioutils
      use params_global
      use pyMOD
      implicit none

      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)
      real twist(kmax)

c..   local variables
      
      real dx,dy,dz,dnorm
      integer k
      integer un

      if (pythonMode.eq.1) then
         nsa=ktip-kroot+1
         nkmax=kmax
         allocate(spanDisPY(nsa),spandist(nkmax))
      endif

      do k=1,kmax

         dx=x(jle,k,1)-x(jtail1,k,1)
         dy=y(jle,k,1)-y(jtail1,k,1)
         dz=z(jle,k,1)-z(jtail1,k,1)   
         
         dnorm=sqrt(dx**2+dy**2+dz**2)
         twist(k)=asin(dz/dnorm)
         
         if ((k-kroot)*(k-ktip).le.0) then
            if (pythonMode.eq.1) then
               spandist(k)=y(jle,k,1)/rartio !y(jle,ktip,1)
               spandisPY(k-kroot+1)=y(jle,k,1)/rartio !y(jle,ktip,1)
            endif
         endif

      enddo

      return
      end
      
c********************************************************************
      subroutine restr2(q,vnu,rest_file)
c
c  read the restart file from unit 3 
c
c********************************************************************
      use ioutils
      use params_global

      implicit none

      real q(jmax,kmax,lmax,nd),vnu(jmax,kmax,lmax)
      character *40 rest_file

c..   local variables
      
      real tsin,tcos,fstipr,alfr,reyprr
      real runew,rvnew
      integer logq,jmaxr,kmaxr,lmaxr,jlo,k,j,l,n,jup

c..   first executable statement
            
      tsin = sin(blang)
      tcos = cos(blang)

      call set_basedir('restart/')
      logq=open_file(rest_file,form='unformatted')

      rewind logq
      if(kmax.eq.3) then

c..2-d grid, add extra plane before and after

        k = 2

        read(logq) jmaxr,lmaxr

        if(jmaxr.ne.jmax) stop ' jmax not the same as qfile dimension'
        if(lmaxr.ne.lmax) stop ' lmax not the same as qfile dimension'

        read(logq) fstipr,alfr,reyprr,totime
        read(logq) (((q(j,k,l,n),j=1,jmax),l=1,lmax),n=1,2),
     <             (((q(j,k,l,n),j=1,jmax),l=1,lmax),n=4,5)
        read(logq) ((vnu(j,k,l),j=1,jmax),l=1,lmax)

        read(logq) istep0

        do 10 j = 1,jmax
          do 10 l = 1,lmax
            q(j,1,l,3) = 0.0
            q(j,2,l,3) = 0.0
            q(j,3,l,3) = 0.0
            q(j,1,l,1) = q(j,2,l,1)
            q(j,1,l,2) = q(j,2,l,2)
            q(j,1,l,4) = q(j,2,l,4)
            q(j,1,l,5) = q(j,2,l,5)
            q(j,3,l,1) = q(j,2,l,1)
            q(j,3,l,2) = q(j,2,l,2)
            q(j,3,l,4) = q(j,2,l,4)
            q(j,3,l,5) = q(j,2,l,5)
 10     continue
      else

c..3-d grid

        read(logq) jmaxr,kmaxr,lmaxr
        if(kmaxr.ne.kmax) stop ' kmax not the same as qfile dimension'
        if(lmaxr.ne.lmax) stop ' lmax not the same as qfile dimension'
        if(iwake.ne.2) then
          if(jmaxr.ne.jmax) stop ' jmax not the same as qfile 
     <          dimension'

          read(logq) fstipr,alfr,reyprr,totime
          read(logq)((((q(j,k,l,n),j=1,jmax),k=1,kmax),l=1,lmax),
     <         n=1,nv)
          read(logq) (((vnu(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax)
          read(logq) istep0

        else
         if(half.eq.0) then

c..add extra j plane at beginning and end
    
         if(jmaxr+2.ne.jmax) stop 'jmax not the same as qfile 
     <           dimension'

          read(logq) fstipr,alfr,reyprr,totime
          read(logq)((((q(j,k,l,n),j=2,jm),k=1,kmax),l=1,lmax),n=1,nv)
          read(logq) istep0
          do 20 k=1,kmax
          do 20 l=1,lmax
            jlo = jle-l+1
            runew = q(jlo,k,lmax-1,2)*tcos + q(jlo,k,lmax-1,3)*tsin
            rvnew =-q(jlo,k,lmax-1,2)*tsin + q(jlo,k,lmax-1,3)*tcos
            q(1,k,l,1) = q(jlo,k,lmax-1,1)
            q(1,k,l,2) = runew
            q(1,k,l,3) = rvnew
            q(1,k,l,4) = q(jlo,k,lmax-1,4)
            q(1,k,l,5) = q(jlo,k,lmax-1,5)
            jup = jle+l-1
            runew = q(jup,k,lmax-1,2)*tcos + q(jup,k,lmax-1,3)*tsin
            rvnew =-q(jup,k,lmax-1,2)*tsin + q(jup,k,lmax-1,3)*tcos
            q(jmax,k,l,1) = q(jup,k,lmax-1,1)
            q(jmax,k,l,2) = runew
            q(jmax,k,l,3) = rvnew
            q(jmax,k,l,4) = q(jup,k,lmax-1,4)
            q(jmax,k,l,5) = q(jup,k,lmax-1,5)
  20      continue
         else

c..add extra j plane at beginning
    
         if(jmaxr+1.ne.jmax) stop 'jmax not the same as qfile 
     <           dimension'

          read(logq) fstipr,alfr,reyprr,totime
          read(logq)((((q(j,k,l,n),j=2,jmax),k=1,kmax),l=1,lmax),
     <         n=1,nv)
          read(logq) istep0
          do 30 k=1,kmax
          do 30 l=1,lmax
            jlo = jle-l+1
            runew = q(jlo,k,lmax-1,2)*tcos + q(jlo,k,lmax-1,3)*tsin
            rvnew =-q(jlo,k,lmax-1,2)*tsin + q(jlo,k,lmax-1,3)*tcos
            q(1,k,l,1) = q(jlo,k,lmax-1,1)
            q(1,k,l,2) = runew
            q(1,k,l,3) = rvnew
            q(1,k,l,4) = q(jlo,k,lmax-1,4)
            q(1,k,l,5) = q(jlo,k,lmax-1,5)
  30      continue
         endif
        endif
      endif

      rewind logq
      close(logq)
      call set_basedir(' ')

      return
      end

c***********************************************************************
      subroutine coarse(x,y,z,q,kkp,kkr)
c
c  coarsens the grid
c  necessary critical points are transformed to the new grid
c
c***********************************************************************

      use params_global

      implicit none

      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)
      real q(jmax,kmax,lmax,nd)
      integer kkp(kmax),kkr(kmax)

c..   local variables
      
      real tcos,tsin
      real xnew,ynew,znew
      real runew, rvnew

      integer jmaxo,kmaxo,lmaxo
      integer jtest,jt2c
      integer j,k,l
      integer jnew,knew,lnew
      integer jmaxc,jmaxf,kmaxc,kmaxf,lmaxc,lmaxf
      integer krootc,ktipc
      integer jlec,jt1c,jbotc
      integer jlo,jup

c..   first executable statement

      tcos = cos(blang)
      tsin = sin(blang)
      
      jmaxo = jmax
      kmaxo = kmax
      lmaxo = lmax
 
c..   upload  in l
      
      if(lint.gt.1)  then
        do 101 l=1,lmax,lint
        lnew = (l-1)/lint+1
        do 101 k=1,kmax
        knew = k
        do 101 j=1,jmax
          jnew = j
          x(jnew,knew,lnew) =   x(j,k,l)
          y(jnew,knew,lnew) =   y(j,k,l)
          z(jnew,knew,lnew) =   z(j,k,l)
          q(jnew,knew,lnew,1) = q(j,k,l,1)
          q(jnew,knew,lnew,2) = q(j,k,l,2)
          q(jnew,knew,lnew,3) = q(j,k,l,3)
          q(jnew,knew,lnew,4) = q(j,k,l,4)
          q(jnew,knew,lnew,5) = q(j,k,l,5)
101     continue
        lmaxc = (lmax-1)/lint+1
        lmaxf = lmax
        lmax = lmaxc
      end if

c..   upload  in k
      
      if(kint.gt.1)  then
        do 103 l=1,lmax
        lnew = l
        do 103 k=1,kmax,kint
        knew = (k-1)/kint+1
        do 103 j=1,jmax
          jnew = j
          x(jnew,knew,lnew) =   x(j,k,l)
          y(jnew,knew,lnew) =   y(j,k,l)
          z(jnew,knew,lnew) =   z(j,k,l)
          q(jnew,knew,lnew,1) = q(j,k,l,1)
          q(jnew,knew,lnew,2) = q(j,k,l,2)
          q(jnew,knew,lnew,3) = q(j,k,l,3)
          q(jnew,knew,lnew,4) = q(j,k,l,4)
          q(jnew,knew,lnew,5) = q(j,k,l,5)
103     continue

        kmaxc = (kmax-1)/kint+1
        kmaxf = kmax
        kmax = kmaxc

        krootc = (kroot-1)/kint+1
        ktipc  = (ktip-1)/kint+1
        kroot = krootc
        ktip  = ktipc
      end if

c..   upload  in j
c..   add symmetry points if half of c-plane used
c..   add to back of grid if iwake=2
      
      if(jint.gt.1)  then
        if(iwake.ne.2) then
          do 105 j=1,jmax,jint
          jnew = (j-1)/jint+1
          do 105 k=1,kmax
          knew = k
          do 105 l=1,lmax
            lnew = l
            x(jnew,knew,lnew) =   x(j,k,l)
            y(jnew,knew,lnew) =   y(j,k,l)
            z(jnew,knew,lnew) =   z(j,k,l)
            q(jnew,knew,lnew,1) = q(j,k,l,1)
            q(jnew,knew,lnew,2) = q(j,k,l,2)
            q(jnew,knew,lnew,3) = q(j,k,l,3)
            q(jnew,knew,lnew,4) = q(j,k,l,4)
            q(jnew,knew,lnew,5) = q(j,k,l,5)
105       continue

          jmaxc = (jmax-1)/jint+1
          if(half.eq.1) jmaxc = (jmax-2)/jint+2
          jlec  = (jle-1)/jint+1
          jt1c  = (jtail1-1)/jint+1
          jt2c  = (jtail2-1)/jint+1
          jbotc = (jbot-1)/jint+1
          jmaxf =  jmax
          jmax =   jmaxc
          jle    = jlec
          jtail1 = jt1c
          jtail2 = jt2c
          jbot   = jbotc
          if(half.eq.1) then
            do 115 k=1,kmax
            do 115 l=1,lmax
              x(jmax,k,l) = x(jmax-2,k,l)
              y(jmax,k,l) = y(jmax-2,k,l)
              z(jmax,k,l) = -z(jmax-2,k,l)
              q(jmax,k,l,1) = q(jmax-2,k,l,1)
              q(jmax,k,l,2) = q(jmax-2,k,l,2)
              q(jmax,k,l,3) = q(jmax-2,k,l,3)
              q(jmax,k,l,4) = -q(jmax-2,k,l,4)
              q(jmax,k,l,5) = q(jmax-2,k,l,5)
 115        continue
          endif
  
        else
          do 205 j=2,jmax-1,jint
          jnew = (j-2)/jint+2
          do 205 k=1,kmax
          knew = k
          do 205 l=1,lmax
            lnew = l
            x(jnew,knew,lnew) =   x(j,k,l)
            y(jnew,knew,lnew) =   y(j,k,l)
            z(jnew,knew,lnew) =   z(j,k,l)
            q(jnew,knew,lnew,1) = q(j,k,l,1)
            q(jnew,knew,lnew,2) = q(j,k,l,2)
            q(jnew,knew,lnew,3) = q(j,k,l,3)
            q(jnew,knew,lnew,4) = q(j,k,l,4)
            q(jnew,knew,lnew,5) = q(j,k,l,5)
205       continue

          jmaxc = (jmax-3)/jint+3
          jlec  = (jle-2)/jint+2
          jt1c  = (jtail1-2)/jint+2
          jt2c  = (jtail2-2)/jint+2
          jbotc = (jbot-2)/jint+2
          jmaxf =  jmax
          jmax =   jmaxc
          jle    = jlec
          jtail1 = jt1c
          jtail2 = jt2c
          jbot   = jbotc

          if(half.eq.0) then
            do 215 k=1,kmax
            do 215 l=1,lmax
              jlo = jle-l+1
              xnew = x(jlo,k,lmax-1)*tcos + y(jlo,k,lmax-1)*tsin
              ynew = -x(jlo,k,lmax-1)*tsin + y(jlo,k,lmax-1)*tcos
              runew = q(jlo,k,lmax-1,2)*tcos + q(jlo,k,lmax-1,3)*tsin
              rvnew = -q(jlo,k,lmax-1,2)*tsin + q(jlo,k,lmax-1,3)*tcos
              x(1,k,l) = xnew
              y(1,k,l) = ynew
              z(1,k,l) = z(jlo,k,lmax-1)
              q(1,k,l,1) = q(jlo,k,lmax-1,1)
              q(1,k,l,2) = runew
              q(1,k,l,3) = rvnew
              q(1,k,l,4) = q(jlo,k,lmax-1,4)
              q(1,k,l,5) = q(jlo,k,lmax-1,5)
              jup = jle+l-1
              xnew = x(jup,k,lmax-1)*tcos + y(jup,k,lmax-1)*tsin
              ynew = -x(jup,k,lmax-1)*tsin + y(jup,k,lmax-1)*tcos
              runew = q(jup,k,lmax-1,2)*tcos + q(jup,k,lmax-1,3)*tsin
              rvnew = -q(jup,k,lmax-1,2)*tsin + q(jup,k,lmax-1,3)*tcos
              x(jmax,k,l) = xnew
              y(jmax,k,l) = ynew
              z(jmax,k,l) = z(jup,k,lmax-1)
              q(jmax,k,l,1) = q(jup,k,lmax-1,1)
              q(jmax,k,l,2) = runew
              q(jmax,k,l,3) = rvnew
              q(jmax,k,l,4) = q(jup,k,lmax-1,4)
              q(jmax,k,l,5) = q(jup,k,lmax-1,5)
 215        continue
          else
            do 225 k=1,kmax
            do 225 l=1,lmax
              jlo = jle-l+1
              xnew = x(jlo,k,lmax-1)*tcos + y(jlo,k,lmax-1)*tsin
              ynew = -x(jlo,k,lmax-1)*tsin + y(jlo,k,lmax-1)*tcos
              runew = q(jlo,k,lmax-1,2)*tcos + q(jlo,k,lmax-1,3)*tsin
              rvnew = -q(jlo,k,lmax-1,2)*tsin + q(jlo,k,lmax-1,3)*tcos
              x(1,k,l) = xnew
              y(1,k,l) = ynew
              z(1,k,l) = z(jlo,k,lmax-1)
              q(1,k,l,1) = q(jlo,k,lmax-1,1)
              q(1,k,l,2) = runew
              q(1,k,l,3) = rvnew
              q(1,k,l,4) = q(jlo,k,lmax-1,4)
              q(1,k,l,5) = q(jlo,k,lmax-1,5)
              x(jmax,k,l) = x(jmax-2,k,l)
              y(jmax,k,l) = y(jmax-2,k,l)
              z(jmax,k,l) = -z(jmax-2,k,l)
              q(jmax,k,l,1) = q(jmax-2,k,l,1)
              q(jmax,k,l,2) = q(jmax-2,k,l,2)
              q(jmax,k,l,3) = q(jmax-2,k,l,3)
              q(jmax,k,l,4) = -q(jmax-2,k,l,4)
              q(jmax,k,l,5) = q(jmax-2,k,l,5)
 225        continue
          endif
        endif
      endif

      jm = jmax - 1
      km = kmax - 1
      lm = lmax - 1
      do 2 k=1,kmax
        kkp(k) = k + 1
        kkr(k) = k - 1
 2    continue
      kkp(kmax) = kmax +ksym*(1-kmax)
      kkr(1) = 1 + ksym*(kmax-1)

      write(6,600)
      write(6,601) jmax,kmax,lmax
      write(6,602) jmaxo,kmaxo,lmaxo
600   format(2x,70('*')/20x,'grid has been modified ',/
     $       2x,70('*'))
601   format(15x,'jmaxnew  =',i3,3x,'kmaxnew  =',i3,
     $     3x,'lmaxnew = ',i3)
602   format(15x,'jmaxold  = ',i3,3x,'kmaxold  = ',i3,3x,
     $ 'lmaxold  = ',i3)

      return
      end

c*************************************************************************
      subroutine number(z,z1,z2,loc1,loc2)
c
c  determine relationship between front and back for wake capturing
c
c*************************************************************************

      use params_global
      
      implicit none

      real z(jmax,kmax,lmax),z1(2*lmax,kmax),z2(2*lmax,kmax)
      integer loc1(2*lmax,kmax),loc2(2*lmax,kmax)

c..   local variables
      
      integer jtop,jb,je,kb,ke,lb,le,js
      integer jd,kd,ld
      integer j,k,l

c..   first executables statement
      
      jd=jmax
      kd=kmax
      ld=lmax

      js=2*lmax
      jtop = jbot + 2 * lmax - 2
      jb = 1
      je = 2 * lmax -1
      kb = 1
      ke = kmax
      lb = 1
      le = lmax

      call converz(z,z1,z2,jbot,jtop,kb,ke,lb,le,js,kd,ld,jmax,jd)
      call search(loc1,z1,z2,jb,je,kb,ke,js,kd)
      call search(loc2,z2,z1,jb,je,kb,ke,js,kd)

      return
      end

c************************************************************************
      subroutine search(loc1,y1,y2,ib,ie,kb,ke,jd,kd)
c     
c  search for locations     
c
c*************************************************************************
      
      implicit none

      integer ib,ie,kb,ke,jd,kd

      integer loc1(jd,kd)
      real y1(jd,kd)
      real y2(jd,kd)
c..   local variables

      integer lb,le
      integer i,j,k,l

c..   first executable statement

      LB = IB
      LE = IE

c..do not forget extra one since finite-volume like for conservative

      do 10 k = kb, ke
      loc1(lb,k) = lb
      loc1(le+1,k) = le+1
      do 10 l = lb+1, le
        do 20 i = ib, ie+1
          if(y2(l,k).le.y1(i,k).and.y2(l,k).gt.y1(i-1,k)) loc1(l,k)= i
   20   continue
   10 continue
c
      return
      end

c**************************************************************************
      subroutine converz(z,z1,z2,jbot,jtop,kb,ke,lb,le,js,kd,ld,jmax,jd)
c  input zc
c  output z1 (for upstream), z2 (for downstream from bottom to top)
c  *lb=1, le= lmax*
c**************************************************************************

      integer jbot,jtop,kb,ke,lb,le,js,kd,ld,jmax,jd
      real z(jd,kd,ld),z1(js,kd),z2(js,kd)

c..   local variables

      real ztemp(jmax)
      integer jj
      integer j,k,l

c..   first executable statement

      do 80 k = kb, ke
      do 80 j = jbot, jtop
        jj = j - jbot + 1
        z1(jj, k) = z(j,k,le)
   80 continue
      do 90 k = kb, ke
      do 90 l = lb, le
        z2(le+1-l,k) = z(1,k,l)
        z2(le-1+l,k) = z(jmax,k,l)
   90 continue

c..   now let us move grid so that q are at cell centers! 
c..   (like finite volume)

      je = 2*le - 1
      do 200 k = kb, ke
        z1(1,k) = z1(1,k)
        z2(1,k) = z2(1,k)
        z1(je+1,k) = z1(je,k)
        z2(je+1,k) = z2(je,k)
        do 150 j = 2,je
          ztemp(j) = 0.5*(z1(j,k)+z1(j-1,k))
  150   continue
        do 160 j = 2,je
          z1(j,k) = ztemp(j)
  160   continue
        do 170 j = 2,je
          ztemp(j) = 0.5*(z2(j,k)+z2(j-1,k))
  170   continue
        do 180 j = 2,je
          z2(j,k) = ztemp(j)
  180   continue
  200 continue

      return
      end

c**********************************************************************
      subroutine eigen( q,xx,xy,xz,yx,yy,yz,zx,zy,zz,ug,vg,wg)
c
c  if dt.lt.0 compute dt from input cnbr
c
c***********************************************************************

      use params_global

      implicit none

      real q(jmax,kmax,lmax,nd)
      real xx(jmax,kmax,lmax),xy(jmax,kmax,lmax),xz(jmax,kmax,lmax),
     <     yx(jmax,kmax,lmax),yy(jmax,kmax,lmax),yz(jmax,kmax,lmax),
     <     zx(jmax,kmax,lmax),zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax),vg(jmax,kmax,lmax),wg(jmax,kmax,lmax)

      real sig(jmax)

c..   local variables

      integer ka,kb,jsig,ksig,lsig,j,k,l
      real u,v,w,edr,a
      real uu,vv,ww,sigx,sigy,sigz,sigmax,dtmin,dtmax,dxmin
      
c..   first executable statement

      ka = 2 - ksym
      kb = km + ksym

      sigmax = -1.e35
      dxmin  = 1.0
      dtmin  = 1000.0
      dtmax  = 0.0

      do 13 l = 2,lm
      do 13 k = ka,kb
        do 11 j = 2,jm

c..the physical velocities u, v, & w and the energy are

          u    = q(j,k,l,2)/q(j,k,l,1)
          v    = q(j,k,l,3)/q(j,k,l,1)
          w    = q(j,k,l,4)/q(j,k,l,1)
          edr  = q(j,k,l,5)/q(j,k,l,1)

c..the speed of sound

          a    = sqrt( ggm1*(edr -0.5*(((u**2) +v**2) +w**2)) )

c..the contravariant velocites u, v, & w are

          uu   = +(u-ug(j,k,l))*xx(j,k,l)+(v-vg(j,k,l))*xy(j,k,l)
     <           +(w-wg(j,k,l))*xz(j,k,l)
          vv   = +(u-ug(j,k,l))*yx(j,k,l)+(v-vg(j,k,l))*yy(j,k,l)
     <           +(w-wg(j,k,l))*yz(j,k,l)
          ww   = +(u-ug(j,k,l))*zx(j,k,l)+(v-vg(j,k,l))*zy(j,k,l)
     <           +(w-wg(j,k,l))*zz(j,k,l)
c
          sigx = abs( uu ) +a*sqrt( ((xx(j,k,l)**2)
     $              +xy(j,k,l)**2) +xz(j,k,l)**2 )
          sigy = abs( vv ) +a*sqrt( ((yx(j,k,l)**2)
     $              +yy(j,k,l)**2) +yz(j,k,l)**2 )
          sigz = abs( ww ) +a*sqrt( ((zx(j,k,l)**2)
     $              +zy(j,k,l)**2) +zz(j,k,l)**2 )
          sig(j) = amax1( sigx,sigy,sigz )
     <             *( 1.0 + 0.002*(1.-timeac)*sqrt(q(j,k,l,6)))
     <             /( 1.0 + (1.-timeac)*sqrt(q(j,k,l,6)) )

   11   continue

        do 12 j = 2,jm
          if( sig(j).gt.sigmax ) then
            sigmax = sig(j)
            jsig   = j
            ksig   = k
            lsig   = l
          end if
          dtmin = amin1( dtmin,h
     <            *( 1.0 + 0.002*(1.-timeac)*sqrt(q(j,k,l,6)))
     <            /( 1.0 + (1.-timeac)*sqrt(q(j,k,l,6)) ) )
          dtmax = amax1( dtmax,h
     <            *( 1.0 + 0.002*(1.-timeac)*sqrt(q(j,k,l,6)))
     <            /( 1.0 + (1.-timeac)*sqrt(q(j,k,l,6)) ) )
   12   continue
   13 continue

      if( dt.lt.0.0 ) then
        dt = cnbr*dxmin/sigmax
      end if

      cnbr = dt*sigmax/dxmin
      h = dt
      hd = .5*dt

c      write(6,601)
c      write(6,602) cnbr,dt,sigmax,jsig,ksig,lsig
c      write(6,*) ' dtmin,dtmax =',dtmin,dtmax

  601 format( '0',5x,'cnbr, dt, sigmax, jsig, ksig, lsig')
  602 format( ' ',5x,3f10.5,2x,3i5,/)

      return
      end

c***********************************************************************
      subroutine qdivj(q)
c
c  divide q by jacobian
c
c***********************************************************************
      use params_global

      implicit none
      real q(jmax,kmax,lmax,nd)

c..   local variables

      integer j,k,l

c..   first executable statement

      do 71 l = 1,lmax
      do 71 k = 1,kmax
      do 71 j = 1,jmax
        q(j,k,l,1) = q(j,k,l,1)/q(j,k,l,6)
        q(j,k,l,2) = q(j,k,l,2)/q(j,k,l,6)
        q(j,k,l,3) = q(j,k,l,3)/q(j,k,l,6)
        q(j,k,l,4) = q(j,k,l,4)/q(j,k,l,6)
        q(j,k,l,5) = q(j,k,l,5)/q(j,k,l,6)
   71 continue

      return
      end

c***********************************************************************
      subroutine qmulj(q)
c
c  multiply q by jacobian
c
c***********************************************************************
      use params_global

      implicit none
      real q(jmax,kmax,lmax,nd)

c..   local variables
      integer j,k,l

c..   first executable statement
      do 71 l = 1,lmax
      do 71 k = 1,kmax
      do 71 j = 1,jmax
        q(j,k,l,1) = q(j,k,l,1)*q(j,k,l,6)
        q(j,k,l,2) = q(j,k,l,2)*q(j,k,l,6)
        q(j,k,l,3) = q(j,k,l,3)*q(j,k,l,6)
        q(j,k,l,4) = q(j,k,l,4)*q(j,k,l,6)
        q(j,k,l,5) = q(j,k,l,5)*q(j,k,l,6)
   71 continue

      return
      end



