

c*********************************************************************** 
      subroutine tim1(x, y, z )
c
c  write out data for acoustics
c  this writes grid in plane and dimensions for pressure and q-file
c
c*********************************************************************** 

      use params_global

      implicit none

      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)

c..   local variables
      
      integer j,k,l,j1
      integer ka,kb,ja,jb,la,lb
      integer jtime,ktime,itime2

      real rad(kmax),theta(jmax+lmax,kmax)
      real xxx,yyy

c..   first executable statement

      ka = (1+ktip)/2
      kb = kmax
      ja = 1
      jb = jle
      la = 1
      lb = lmax

      jtime = (jb-ja)+(lb-la)+1
      ktime = (kb-ka)+1

      itime2 = (nsteps-istep0+nacoup-1)/nacoup

      do 50 k = ka,kb
        rad(k) = sqrt(x(2,k,1)**2+y(2,k,1)**2)
        j = jb
        do 100 l = lb,la+1,-1
          j1 = lb-l+1
          xxx = x(j,k,l) + 0.25
          yyy = y(j,k,l)
          theta(j1,k) = atan(xxx/yyy)
          if (xxx.gt.0 .and. yyy.lt.0) theta(j1,k) = theta(j1,k) + pi
  100   continue

        l = la
        do 200 j = jb,ja,-1
          j1 = (lb-la)+jb-j+1
          xxx = x(j,k,l) + 0.25
          yyy = y(j,k,l)
          theta(j1,k) = atan(xxx/yyy)
c         if (xxx.gt.0 .and. yyy.lt.0) theta(j1,k) = theta(j1,k) + pi
          if (yyy.lt.0) theta(j1,k) = theta(j1,k) + pi
  200   continue

 50   continue

      return
      end
      

c**********************************************************************
      subroutine newton(q,qnewt,s,ka,kb)
c
c  take into account time terms on rhs for newton iterations
c
c**********************************************************************
      
      use params_global

      implicit none

      integer ka,kb
      real q(jmax,kmax,lmax,nd),qnewt(jmax,kmax,lmax,nd)
      real s(jmax,kmax,lmax,nv)
      
      integer j,k,l
      real oat,tac

c..set time-accuracy
      oat = 1.0
      if (ntac.eq.2 .and. istep.gt.1) oat = 2./3.
      
      do 10 l = 2, lm
      do 10 k = ka, kb
      do 10 j = 2, jm
c        tac = 1./( h*( 1.0 + 0.002*(1.-timeac)*sqrt(q(j,k,l,6)) ) 
c     <              /( 1.0 + (1.-timeac)*sqrt(q(j,k,l,6)) ) )
        tac=1.0D0/(oat*h)
        s(j,k,l,1) = s(j,k,l,1) - (q(j,k,l,1) - qnewt(j,k,l,1))*tac
        s(j,k,l,2) = s(j,k,l,2) - (q(j,k,l,2) - qnewt(j,k,l,2))*tac
        s(j,k,l,3) = s(j,k,l,3) - (q(j,k,l,3) - qnewt(j,k,l,3))*tac
        s(j,k,l,4) = s(j,k,l,4) - (q(j,k,l,4) - qnewt(j,k,l,4))*tac
        s(j,k,l,5) = s(j,k,l,5) - (q(j,k,l,5) - qnewt(j,k,l,5))*tac
   10 continue
c     
      return
      end


c***********************************************************************
      subroutine stqol(q,qnewt,qtn,vnu,vnu0)
c     
c  store data at previous time steps
c
c*********************************************************************** 

      use params_global

      implicit none

      integer l,k,j
      real q(jmax,kmax,lmax,nd),qnewt(jmax,kmax,lmax,nd),
     &     qtn(jmax,kmax,lmax,nd),vnu(jmax,kmax,lmax),
     &     vnu0(jmax,kmax,lmax) 

      do 40 l = 1, lmax 
      do 40 k = 1, kmax
        do 10 j = 1, jmax
          qnewt(j,k,l,1) = q(j,k,l,1)
          qnewt(j,k,l,2) = q(j,k,l,2)
          qnewt(j,k,l,3) = q(j,k,l,3)
          qnewt(j,k,l,4) = q(j,k,l,4)
          qnewt(j,k,l,5) = q(j,k,l,5)
          vnu0(j,k,l)    = vnu(j,k,l)
 10    continue
        if (ntac.eq.2 .and. istep.gt.1) then
          do 20 j = 1, jmax
            qnewt(j,k,l,1) = qnewt(j,k,l,1)+(q(j,k,l,1)-qtn(j,k,l,1))/3.
            qnewt(j,k,l,2) = qnewt(j,k,l,2)+(q(j,k,l,2)-qtn(j,k,l,2))/3.
            qnewt(j,k,l,3) = qnewt(j,k,l,3)+(q(j,k,l,3)-qtn(j,k,l,3))/3.
            qnewt(j,k,l,4) = qnewt(j,k,l,4)+(q(j,k,l,4)-qtn(j,k,l,4))/3.
            qnewt(j,k,l,5) = qnewt(j,k,l,5)+(q(j,k,l,5)-qtn(j,k,l,5))/3.
 20      continue
        endif
        do 30 j = 1, jmax
          qtn(j,k,l,1) = q(j,k,l,1)
          qtn(j,k,l,2) = q(j,k,l,2)
          qtn(j,k,l,3) = q(j,k,l,3)
          qtn(j,k,l,4) = q(j,k,l,4)
          qtn(j,k,l,5) = q(j,k,l,5)
 30    continue
 40   continue

      return
      end

      subroutine newton_timescale(timeac, oat, h, q, iblank, tscale, bt, jmax, kmax, lmax)
      !DESC: Computes local time scaling for Newton iterations. If timeac
      !      is 1.0 then a constant scaling is used everywhere.
  
         implicit none
         integer, intent(in) :: jmax, kmax, lmax
         real, intent(in) :: h, timeac, oat
         integer, dimension(jmax,kmax,lmax), intent(in) :: iblank
         real, dimension(jmax, kmax, lmax,6), intent(in) :: q
         real, dimension(jmax, kmax, lmax), intent(out) :: tscale, bt

         !local var
         integer :: j, k, l
         real :: s

         do l = 1, lmax
           do k = 1, kmax
             do j = 1, jmax

               s = ( 1.0 + 0.002*(1.0 - timeac)*sqrt(q(j,k,l,6)) )
               s = s/( 1.0 + (1.0 - timeac)*sqrt(q(j,k,l,6)) )

               tscale(j,k,l) = oat*h*s*max(iblank(j,k,l),0)
               bt(j,k,l) = 0.

             end do
           end do
         end do

      end subroutine newton_timescale
