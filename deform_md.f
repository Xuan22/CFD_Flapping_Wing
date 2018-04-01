! Math: add pitch/plunge kinematics
C*******************************************************************************
      subroutine pitch_plunge_wing(x,y,z,xg,yg,zg,ug,vg,wg,zx,zy,zz,
     &                             zx0,zy0,zz0,zt0,init)
C*******************************************************************************
C....Added by Ria and Vinod
C*****************************************************************************
      use params_global

      implicit none

      real zx(jmax,kmax,lmax),zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax), vg(jmax,kmax,lmax),wg(jmax,kmax,lmax)
      real xg(jmax,kmax,lmax), yg(jmax,kmax,lmax),zg(jmax,kmax,lmax)
      real x(jmax,kmax,lmax), y(jmax,kmax,lmax), z(jmax,kmax,lmax )
      real zx0(jmax,kmax),zy0(jmax,kmax),zz0(jmax,kmax),zt0(jmax,kmax)
      real xnew, ynew, znew
      real omega_t, alpha
      real flap_angle
      integer j,k,l
      integer init
C***********************************************************************

      if (init.eq.0) then
        do j = 1,jmax
          do k = 1,kmax
            zx0(j,k) = zx(j,k,1)
            zy0(j,k) = zy(j,k,1)
            zz0(j,k) = zz(j,k,1)
            zt0(j,k) = -ug(j,k,1)*zx(j,k,1)-vg(j,k,1)*zy(j,k,1)
     &                 -wg(j,k,1)*zz(j,k,1)
          enddo
        enddo
      else
        zx0 = 0. 
        zy0 = 0. 
        zz0 = 0. 
        zt0 = 0. 
      endif

      do l = 1,lmax
      do k = 1,kmax
      do j = 1,jmax
        ug(j,k,l) = 0
        vg(j,k,l) = 0
        wg(j,k,l) = 0
      enddo
      enddo
      enddo

      omega_t=omega*(totime)
      alpha = pitch_knot*pi/180+pitch_angle*pi/180*sin(omega_t+phase*pi/180)
      flap_angle = flap_amp*pi/180*(sin(omega_t))

c..end here

      do l = 1,lmax
      do k = 1,kmax
      do j = 1,jmax

       if (isroot.eq.1) then

c   if (j.eq.134 .and. k.eq.1 .and. l.eq.1 ) then
c 	 print *,'xg=',xg(j,k,l)
c 	 print *,'yg=',yg(j,k,l)
c	 print *,'zg=',zg(j,k,l)
c 	endif

          znew=(yg(j,k,l)+root)*sin(flap_angle)-(xg(j,k,l)-xac)*sin(alpha)
     &      +(zg(j,k,l)-zac)*cos(alpha)+zac

          xnew=(xg(j,k,l)-xac)*cos(alpha)
     &      +(zg(j,k,l)-zac)*sin(alpha)+xac

          ynew = (yg(j,k,l)+root)*cos(flap_angle)

       else

        xnew = (xg(j,k,l)-xac)*cos(alpha) + (zg(j,k,l)-zac)*sin(alpha) + xac
        ynew = y(j,k,l)
        znew = plunge_amp*sin(omega_t)-(xg(j,k,l)-xac)*sin(alpha)
     &                                +(zg(j,k,l)-zac)*cos(alpha) + zac

       endif

        if (init.eq.0) then
          ug(j,k,l) = (xnew-x(j,k,l))/dt
          vg(j,k,l) = (ynew-y(j,k,l))/dt
          wg(j,k,l) = (znew-z(j,k,l))/dt
        endif

        x(j,k,l)=xnew
        y(j,k,l)=ynew
        z(j,k,l)=znew

c   if (j.eq.134 .and. k.eq.1 .and. l.eq.1) then
c 		print *,'x=',x(j,k,l)
c 		print *,'y=',y(j,k,l)
c       print *,'z=',z(j,k,l)
c       print *,'Ug=',ug(j,k,l)
c       print *,'Vg=',vg(j,k,l)
c		print *,'Wg=',wg(j,k,l)
c	endif

      enddo
      enddo
      enddo

      end subroutine pitch_plunge_wing

C*******************************************************************************
      subroutine flap_pitch_wing_md(x,y,z,xg,yg,zg,ug,vg,wg,zx,zy,zz,
     &                             zx0,zy0,zz0,zt0,init)
C*******************************************************************************
C....Added by Xuan, 03/20/2018 for kinematics of 60% span of flapping wing
C....Added by Camli and Mathieu
C... 7/7/2014 Camli add stroke plane beta to kinematics
C... 3/2015 add smooth transition to new variables where i.e :
C...        beta=beta_old+d*(dbeta) and d=0.5*(1-cos(psi))
C*****************************************************************************
      use params_global

      implicit none

      real zx(jmax,kmax,lmax), zy(jmax,kmax,lmax), zz(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax), vg(jmax,kmax,lmax), wg(jmax,kmax,lmax)
      real xg(jmax,kmax,lmax), yg(jmax,kmax,lmax), zg(jmax,kmax,lmax)
      real x(jmax,kmax,lmax), y(jmax,kmax,lmax), z(jmax,kmax,lmax )
      real zx0(jmax,kmax), zy0(jmax,kmax), zz0(jmax,kmax), zt0(jmax,kmax)
      real xnew, ynew, znew
      integer j,k,l
      integer init
      integer :: ro=1
     
C***********************************************************************
      real :: flap1(3),flap2(3),pitch1(3),pitch2(3)
      real :: dphi,dalpha
      real :: temp
      real :: deg
      real :: p_new(3)
      real :: kin_span = 1.0
      real :: azm, azm1
      real :: d = 0.0
      real :: beta_old, phiOff_old, phiMax_old
      real :: delta_phiMax, delta_kin_phiOff, delta_kin_beta
      real :: alpha_g, delta_beta=0.0, delta_phiOff=0.0
      !real :: dt1, totime1

      deg = pi/180 ! deg to radians
      !dt1 = kin_dt
      !totime1 = totime/dt*dt1


      if (init.eq.0) then
        do j = 1,jmax
          do k = 1,kmax
            zx0(j,k) = zx(j,k,1)
            zy0(j,k) = zy(j,k,1)
            zz0(j,k) = zz(j,k,1)
            zt0(j,k) = -ug(j,k,1)*zx(j,k,1)-vg(j,k,1)*zy(j,k,1)
     &                 -wg(j,k,1)*zz(j,k,1)
          enddo
        enddo
      else
        zx0 = 0. 
        zy0 = 0. 
        zz0 = 0. 
        zt0 = 0. 
      endif

      do l = 1,lmax
      do k = 1,kmax
      do j = 1,jmax
        ug(j,k,l) = 0
        vg(j,k,l) = 0
        wg(j,k,l) = 0
      enddo
      enddo
      enddo
C  *********************************************************************
C  ********* added by Camli Badrya July/30 ***************
C  ********* define the relaxtion factor d ***********
C  ** need to define: d and azm1
C  ** phiMax_old
C  ** phiOff_old
C  ** beta_Old
C  *********************************************************************
      azm  = rf *totime*180/pi
      !kin_xoff1 = xoff

      !!! deefine ismth to switch the relaxation 
      IF (ismth .GT.0) THEN

         azm1 = abs(azm) - Ncyc*360

         IF ((azm1 .GT.0) .AND. (azm1 .LT.90) ) THEN
             d = 0.5 *(1- cos(2*azm1*pi/180))
         ELSE IF (azm1 .GE. 90) THEN
             d=1
         ELSE
             d=0
         END IF

         phiMax_old = kin_phiMax
         phiOff_old = kin_phiOff
         beta_old   = kin_beta

         delta_phiMax = kin_phiMax_new - phiMax_old
         delta_kin_phiOff = kin_phiOff_new - phiOff_old
         delta_kin_beta = kin_beta_new - beta_old
      
         IF (delta_phiMax .NE. 0) THEN
             kin_phiMax = phiMax_old + (delta_phiMax)*d
         END IF 

         IF (delta_kin_phiOff .NE. 0) THEN
             kin_phiOff = phiOff_old + delta_kin_phiOff*d
             delta_phiOff = kin_phiOff - phiOff_old
         END IF

         IF (delta_kin_beta .NE. 0) THEN
             kin_beta   = beta_old + (delta_kin_beta)*d
             delta_beta = kin_beta - beta_old
         END IF

      END IF

c.... define alpha_geomtric 
      alpha_g = kin_alphaInit - kin_beta

      !****----Flap----************************************************************
      !define rotation axis for flapping
      !point1: at origin (is at root of wing at point along chord specified
      !        by kin_rotAxis)

C... define the parameters for the the second wing rotation 
      if (irot_dir.eq.-1) then
          kin_xoff    = -xoff
         !kin_phiInit = -kin_phiInit
         !kin_phiOff  = -kin_phiOff
          ro = -1
      endif
C... done

      flap1(1)= 0.0 + kin_xoff
      flap1(2)= 0.0 + yoff
      flap1(3)= 0.0 + zoff
      !point2: origin +1 in z-direction (rotational axis is z unit vector)
      flap2(1)= 0.0 + kin_xoff
      flap2(2)=-sin(kin_beta*deg) + yoff
      flap2(3)= cos(kin_beta*deg) + zoff

      !if restarting, don't use offset angles, the mesh is already where it should be
      !if (iread.gt.0) then
         !kin_phiInit   =  0.0
         !kin_phiOff    =  0.0
         !kin_alphaInit =  0.0
         !alpha_g       =  0.0
      !endif

      !update flapping angle:
      if (init.eq.1) then
         dphi   = ro*kin_phiOff !kin_phiOff ! set grid to initial flap position
         kin_phi= ro*kin_phiInit+dphi
      else
         !1st order derivative of kinematics equation:
         dphi=dt *rf *kin_phiMax *sin(abs(rf) *totime 
     &                             + kin_phi_phase *deg) + ro*delta_phiOff
         kin_phi=kin_phi+dphi
      endif

      !****----PITCH----**************************************************************
      !define rotation axis for pitching:
      !point1: at origin (is at root of wing at point along chord specified
      !        by kin_rotAxis)
      pitch1(1)= 0 + kin_xoff
      pitch1(2)= 0 + yoff
      pitch1(3)= 0 + zoff
      !point2: at tip of wing, at same chord position as point1
      pitch2(1)= -kin_span *sin(kin_phi*deg) + kin_xoff
      pitch2(2)=  kin_span *cos(kin_phi*deg) *cos(kin_beta*deg) + yoff
      pitch2(3)=  kin_span *cos(kin_phi*deg) *sin(kin_beta*deg) + zoff
      
      !update pitching angle
      if (init.eq.1) then
         dalpha = ro*alpha_g		! set grid to initial pitch position
         kin_alpha = kin_alphaInit
      else
         !1st order derivative of kinematics equation:
         temp = abs(rf) *totime + kin_alpha_phase*deg
         dalpha = -dt*abs(rf) *kin_alphaMax *sin(temp) - ro *delta_beta
         kin_alpha = kin_alpha +dalpha
      endif
      
      write(STDOUT,101) kin_phi,kin_alpha,kin_beta,dt,totime,azm,d
      write(STDOUT,102) delta_phiOff, delta_beta, rf
	  !******************************************************************************
      do l = 1,lmax
		  do k = 1,kmax
		  do j = 1,jmax

		     !! simple rotation
		     !temp = x(j,k,l)
		     !p_new(1) = cos(dt*rf)*x(j,k,l) - sin(dt*rf)*y(j,k,l)
		     !p_new(2) = cos(dt*rf)*y(j,k,l) + sin(dt*rf)*temp
		     !p_new(3) = z(j,k,l)

		     p_new(1) = x(j,k,l)
		     p_new(2) = y(j,k,l)
		     p_new(3) = z(j,k,l)

		     !apply flap rotation to each grid point
		     call rotatePointArbitraryAxis(p_new, flap1, flap2, dphi*deg)

		     !apply pitch rotation to each grid point
		     call rotatePointArbitraryAxis(p_new, pitch1, pitch2, dalpha*deg)

		     if (init.eq.0) then
		        ug(j,k,l) = (p_new(1) -x(j,k,l))/dt
		        vg(j,k,l) = (p_new(2) -y(j,k,l))/dt
		        wg(j,k,l) = (p_new(3) -z(j,k,l))/dt
		     !else !set initial grid
		     !   xg(j,k,l) = p_new(1)
		     !   yg(j,k,l) = p_new(2)
		     !   zg(j,k,l) = p_new(3)
		     endif

		     x(j,k,l) = p_new(1)
		     y(j,k,l) = p_new(2)
		     z(j,k,l) = p_new(3)

		  enddo
		  enddo
      enddo
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! --Xuan 03/30/2018
! --add to consider wing twist effect
	  do l = 1,lmax
	    do k = 1,kmax
	    do j = 1,jmax
		if (t_index.eq.1) then   
			d_twist = wing_twist(1,x(j,k,l)) -wing_twist(timestep,x(j,k,l))
			if (k.eq.1 .and. j.eq.1) then
				d_twist = wing_twist(1,xLocs(2)) -wing_twist(timestep,xLocs(2))
			end if
		else
			d_twist = wing_twist(t_index,x(j,k,l)) -wing_twist(t_index-1,x(j,k,l))
			if (k.eq.1 .and. j.eq.1) then
				d_twist = wing_twist(t_index,xLocs(2)) -wing_twist(t_index-1,xLocs(2))
			end if
		end if
		end do
		end do
	  end do
	
      do l = 1,lmax
		  do k = 1,kmax
		  do j = 1,jmax
			p_new(1) = x(j,k,l)
			p_new(2) = y(j,k,l)
			p_new(3) = z(j,k,l)
			call rotatePointArbitraryAxis(p_new, pitch1, pitcht2, -d_twist)
		  end do
		  end do
	  end do
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

101   format(/,' Flapping / pitching wing: Phi = ',f10.5,', Alpha = ',f10.5,', Beta = ',f10.5,' (dt=',f10.5,
     &             ', totime=',f10.5,') ',', azm =' ,f10.5,', d = ',f10.5)
102   format(/,' Flapping / pitching wing: delta_phiOff = ',f10.5,', delta_beta = ',f10.5,', rf = ',f10.5)
      contains
         !----------------------------------------------------------------
         subroutine rotatePointArbitraryAxis(p_old, p1, p2, angle)
         !
         ! Rotate point about an arbitrary axis defined by 2 points. 
         ! The code is the full rotation/translation matrix of the 
         ! following procedure:
         !  1) translate rotation axis so that it passes through origin
         !  2) rotate space about xaxis to put rotation axis in xz plane
         !  3) rotate space about yaxis to put rotation axis along +z axis
         !  4) rotate point about z axis
         !  5) perform inverse of steps 3,2,1 to return rotation axis
         !
         !  p_old(3) -> x,y,z coordinates of point to be moved
         !  p_new(3) -> x,y,z coordinates of point's new location
         !  p1(3)    -> x,y,z coordinates of point on rotation axis
         !  p2(3)    -> x,y,z coordinates of point on rotation axis located
         !              in the positive axis direction from p1
         !  angle    -> rotation angle about the rotation axis in radians
         !----------------------------------------------------------------
         implicit none
         !----------------------------------------------------------------
         real, intent(inout) :: p_old(:)
         real, intent(in) :: p1(3), p2(3), angle
         ! Local variables
         real :: ux,uy,uz,length
         real :: cn,sn
         real :: p_new(3)

         cn = cos(angle)
         sn = sin(angle)

         ! calculate unit vector along rotation axis
         ux = p2(1) -p1(1)
         uy = p2(2) -p1(2)
         uz = p2(3) -p1(3)

         ! normalize unit vector along rotation axis
         length = 1.0d0/sqrt(ux*ux+uy*uy+uz*uz)
         ux = ux*length
         uy = uy*length
         uz = uz*length

         ! calculate vector from point p_old to p1
         p_old(1) = p_old(1) - p1(1)
         p_old(2) = p_old(2) - p1(2)
         p_old(3) = p_old(3) - p1(3)

         p_new(1) = p_old(1)*( cn+(1-cn)*ux*ux )
     &            + p_old(2)*( (1-cn)*ux*uy-uz*sn )
     &            + p_old(3)*( (1-cn)*ux*uz+uy*sn )
     &            + p1(1)

         p_new(2) = p_old(1)*( (1-cn)*ux*uy+uz*sn )
     &            + p_old(2)*( cn+(1-cn)*uy*uy )
     &            + p_old(3)*( (1-cn)*uz*uy-ux*sn )
     &            + p1(2)

         p_new(3) = p_old(1)*( (1-cn)*ux*uz-uy*sn )
     &            + p_old(2)*( (1-cn)*uy*uz+ux*sn )
     &            + p_old(3)*( cn+(1-cn)*uz*uz )
     &            + p1(3)

         p_old = p_new


         return
         end subroutine rotatePointArbitraryAxis

      end subroutine flap_pitch_wing

c***********************************************************************
      subroutine init_deform(srot,x,y,z,xx,xy,xz,ug,yx,yy,yz,
     &     vg,zx,zy,zz,wg,zx0,zy0,zz0,zt0,
     &     xg,yg,zg,xt2,yt2,zt2)

c***  Prlogue : ********
c
c     Subroutine for deforming the mesh according to the linear
c     and rotational deflections. This subroutine deforms and 
c     creates the initial mesh. i.e., no grid velocities are evaluated
c     here. The rest of the subroutine follows the same as the 
c     deform() subroutine below, see the documentation of that
c
c     last updated by jaina 07/10/2002
c***********************************************************************

      use mpi_wrapper, only: print_message
      use params_global
      use deflections
      use arf_mod

      implicit none

      real srot
      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)
      real xx(jmax,kmax,lmax),xy(jmax,kmax,lmax),xz(jmax,kmax,lmax)
      real yx(jmax,kmax,lmax),yy(jmax,kmax,lmax),yz(jmax,kmax,lmax)
      real zx(jmax,kmax,lmax),zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax),vg(jmax,kmax,lmax),wg(jmax,kmax,lmax)
      real xg(jmax,kmax,lmax),yg(jmax,kmax,lmax),zg(jmax,kmax,lmax)
      real zx0(jmax,kmax),zy0(jmax,kmax),zz0(jmax,kmax),zt0(jmax,kmax)
      real xt2(jmax,kmax,lmax),yt2(jmax,kmax,lmax),zt2(jmax,kmax,lmax)

c..   local variables

      real rt(3,3)
      real xpsi,dpsi,cs,ss,beta,cb,sb
      real epsdist,dmax,fac,fac1,dist
      real e1,e2,e3,e4,e5,e6
      real xe,ye,ze
      real a1,a2,a3
      real a11,a12,a13,a14
      real a21,a22,a23,a24
      real a31,a32,a33,a34
      real rt11,rt12,rt13,rt21,rt22,rt23,rt31,rt32,rt33
      real xtmp1,ytmp1,ztmp1,xtmp,ytmp,ztmp

      integer i,j,k,l   

c..   first exectuable statement
c     xpsi=totime*rf

      call print_message("Start of init_deform")

      xpsi=srot

      if(arf_opt .eq. 0)then
        dpsi=2*pi/iazimuth
        i=mod(nint(xpsi/dpsi),iazimuth)+1
      else
        i = 1
      end if

      cs=cos(xpsi)
      ss=sin(xpsi)

      beta=beta0+beta1c*cs+beta1s*ss
      cb=cos(beta*pi/180)
      sb=sin(beta*pi/180)

c..   decay parameters

      epsdist=1.5
      dmax=10.0
      fac=0.5

      do k=1,kmax

         xe=xyzref(1,k)
         ye=xyzref(2,k)
         ze=xyzref(3,k)
         
         do j=1,jmax
            do l=1,lmax

c..   turns coordinate system 90 degree clockwise of deformations

               a1=xg(j,k,l)+ye
               a2=yg(j,k,l)-xe
               a3=zg(j,k,l)-ze

               dist=sqrt(a1**2+a2**2+a3**2)
               dist=max(0.0,dist-epsdist)

               fac1=1.0
               if (idecay.eq.1) then
                   if (dist.gt.dmax) then
                       fac1=0.0
                   else
                       fac1=(1+cos(pi*dist/dmax))*0.5
                   endif
	       endif
               
               e4=defl_dat(4,k,i)
               e5=defl_dat(5,k,i)
               e6=defl_dat(6,k,i)

               if (idefl.eq.1) then
                  call rmat_extract(e4,e5,e6,rt)
               elseif (idefl.eq.2) then
                  call rmat_extract_rcas(e4,e5,e6,rt)
               else
                  e4=e4*fac1
                  e5=e5*fac1
                  e6=e6*fac1
                  call rmat_extract_tdu(e4,e5,e6,rt)
               endif

               rt11=rt(1,1);rt12=rt(1,2);rt13=rt(1,3);
               rt21=rt(2,1);rt22=rt(2,2);rt23=rt(2,3);
               rt31=rt(3,1);rt32=rt(3,2);rt33=rt(3,3);
               
               e1=defl_dat(1,k,i)*rartio*fac1
               e2=defl_dat(2,k,i)*rartio*fac1
               e3=defl_dat(3,k,i)*rartio*fac1

               a11=ss*rt12+cs*rt22
               a12=-ss*rt11-cs*rt21
               a13=-ss*rt13-cs*rt23
               a14=-(-ss*rt11-cs*rt21)*xe-(-ss*rt12-cs*rt22)*ye
     &              -(-ss*rt13-cs*rt23)*ze-ss*xe-cs*ye-ss*e1-cs*e2
               
               a21=-cs*rt12+ss*rt22
               a22=cs*rt11-ss*rt21
               a23=cs*rt13-ss*rt23
               a24=-(cs*rt11-ss*rt21)*xe-(cs*rt12-ss*rt22)*ye
     &              -(cs*rt13-ss*rt23)*ze+cs*xe-ss*ye+cs*e1-ss*e2
               
               a31=-rt32
               a32=rt31
               a33=rt33
               a34=-rt31*xe-rt32*ye-rt33*ze+ze+e3
         
c..   introduce rigid coning and flapping
               
               xtmp1=xg(j,k,l)
               ytmp1=(yg(j,k,l)-eflap)*cb-zg(j,k,l)*sb + eflap
               ztmp1=zg(j,k,l)*cb+(yg(j,k,l)-eflap)*sb
                 
               xtmp=xtmp1*a11+ytmp1*a12+ztmp1*a13+a14
               ytmp=xtmp1*a21+ytmp1*a22+ztmp1*a23+a24
               ztmp=xtmp1*a31+ytmp1*a32+ztmp1*a33+a34  
               
               x(j,k,l)=xtmp
               y(j,k,l)=ytmp
               z(j,k,l)=ztmp
               
            enddo
         enddo
         
      enddo

      return
      end



c***********************************************************************
      subroutine rmat_extract(e1,e2,e3,rmat)

c**   Prologue :
c
c     transformation matrices
c     for euler parameters
c***********************************************************************

      real rmat(3,3)

      e0= sqrt(1.-e1*e1-e2*e2-e3*e3)

      rmat(1,1)= 1.0 - 2.0*(e2*e2+e3*e3)
      rmat(1,2)= 2.0*(e1*e2-e3*e0)
      rmat(1,3)= 2.0*(e1*e3+e2*e0)
      rmat(2,1)= 2.0*(e1*e2+e3*e0)
      rmat(2,2)= 1.0 - 2.0*(e1*e1+e3*e3)
      rmat(2,3)= 2.0*(e2*e3-e1*e0)
      rmat(3,1)= 2.0*(e1*e3-e2*e0)
      rmat(3,2)= 2.0*(e2*e3+e1*e0)
      rmat(3,3)= 1.0 - 2.0*(e1*e1+e2*e2)

      return
      end

c***********************************************************************
      subroutine rmat_extract_rcas(e1,e2,e3,rmat)

c***  Prologue :
c
c     e1->dxang
c     e2->dyang
c     e3->dzang
c
c  -       -     -       -     -        -
c |1  0   0 |   |ct  0  st|   |cs -ss  0 |
c |0 cp -sp | X |0   1  0 | X |ss  cs  0 |
c |0 sp  cp |   |-st 0  ct|   |0   0   1 |
c  -       -     -       -     -        -
c***********************************************************************

      real rmat(3,3)

      cp=cos(e1)
      sp=sin(e1)

      ct=cos(e2)
      st=sin(e2)

      cs=cos(e3)
      ss=sin(e3)

      rmat(1,1)=ct*cs
      rmat(2,1)=ct*ss
      rmat(3,1)=-st

      rmat(1,2)=sp*st*cs-cp*ss
      rmat(2,2)=sp*st*ss+cp*cs
      rmat(3,2)=sp*ct

      rmat(1,3)=cp*st*cs+sp*ss
      rmat(2,3)=cp*st*ss-sp*cs
      rmat(3,3)=cp*ct

      return
      end

c***********************************************************************

      subroutine rmat_extract_tdu(e1,e2,e3,rmat)
      real rmat(3,3)

c**   Prologue :
c
c     tdu from umarc manual
c     UMARC format (u,v,w,v',w',phi)
c
c     e1 -> v'
c     e2 -> w'
c     e3 -> phi
c     
c***********************************************************************

      ct=cos(e3)
      st=sin(e3)

      rmat(1,1)=1-0.5*(e1**2+e2**2)
      rmat(2,1)=e1
      rmat(3,1)=e2
      
      rmat(1,2)=-(e1*ct+e2*st)
      rmat(2,2)=(1-e1**2*0.5)*ct-e1*e2*st
      rmat(3,2)=st*(1-e2**2*0.5)

      rmat(1,3)=e1*st-e2*ct
      rmat(2,3)=-(1-e1**2*0.5)*st-e1*e2*ct
      rmat(3,3)=ct*(1-e2**2*0.5)
      
      return
      end        


c***********************************************************************
      subroutine deform(srot,psi_rot,x,y,z,xx,xy,xz,ug,yx,yy,yz,
     &     vg,zx,zy,zz,wg,zx0,zy0,zz0,zt0,xt2,
     &     yt2,zt2,xg,yg,zg)

c***  Prologue:
c
c     Deform the blade according to flap-lag-torsion deflections
c     The linear and angular deformations (totally 6 in number) is
c     known at every section and every azimuthal station. These are
c     used to deform the mesh at every station. A cosine decay can
c     be applied such that the outer boundaries remain stationary if
c     necessary.
c
c     Three formats of deformation formats are accepted now
c
c
c     idefl = 0 -> (u,v,w,v',w',phi)
c                  where u,v,w are linear deflections in x,y and z
c                  primes represent derivatives and phi is the elastic 
c                  torsional deformation. Make sure the cyclic and 
c                  collective angles are included in the input file
c                  in this case.
c                  (UMARC follows this)
c
c     idefl = 1 -> euler parameters  (e1,e2,e3 -> linear, 
c                                     e4,e5,e6 -> rotational)
c                  (DYMORE follows this)
c
c     idefl = 2 -> (dx,dy,dz,dxang,dyang,dzang)
c                  Most general and straight forward format
c                  linear deflections in x,y and z and rotations
c                  about each axis
c                  (CAMRAD and RCAS follow this)
c
c
c     Last updated by jaina 07/13/04
c***  End prologue :************************************************


      use params_global
      use deflections
      use arf_mod

      implicit none

      real srot,psi_rot
      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)
      real xx(jmax,kmax,lmax),xy(jmax,kmax,lmax),xz(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax)
      real yx(jmax,kmax,lmax),yy(jmax,kmax,lmax),yz(jmax,kmax,lmax)
      real vg(jmax,kmax,lmax)
      real zx(jmax,kmax,lmax),zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)
      real wg(jmax,kmax,lmax)
      real zx0(jmax,kmax),zy0(jmax,kmax),zz0(jmax,kmax),zt0(jmax,kmax)
      real xt2(jmax,kmax,lmax),yt2(jmax,kmax,lmax),zt2(jmax,kmax,lmax)
      real xg(jmax,kmax,lmax),yg(jmax,kmax,lmax),zg(jmax,kmax,lmax)

c..   local variables

      real rt(3,3)
      real dpsi,cs,ss,beta,cb,sb,xpsi,fsm1,dtval
      real epsdist,dmax,fac,fac1,dist
      real e1,e2,e3,e4,e5,e6
      real xe,ye,ze
      real a1,a2,a3
      real a11,a12,a13,a14
      real a21,a22,a23,a24
      real a31,a32,a33,a34
      real rt11,rt12,rt13,rt21,rt22,rt23,rt31,rt32,rt33
      real xtmp1,ytmp1,ztmp1,xtmp,ytmp,ztmp
      
      integer i,j,k,l,n,nstp

C asitav (slow start of deformation)
      real    :: scal,tbar
      integer :: nstps

c..   first executable statement

      xpsi=psi_rot
      dpsi=2*pi/iazimuth

      if (arf_opt .eq. 0) then
         tbar=(istep-1)*dpsi*6.0/pi
         
         if (tbar .lt. 1.0D0) then
            scal=3*tbar**2-2*tbar**3
         else
            scal=1.0D0
         endif
         
         if (def_first) then
            scal=1.0D0
         endif
         
      endif
 8901 format(i5,i5,2(1x,f10.5))
         
      i=mod(nint(xpsi/dpsi),iazimuth)+1
      fsm1=fsmach

      if(arf_opt .eq. 0)then
        dpsi=2*pi/iazimuth
        i=mod(nint(xpsi/dpsi),iazimuth)+1
      else
        i = 2
      end if

      if (irotate.eq.0) then
         cs=cos(xpsi)
         ss=sin(xpsi)
      else
         cs=cos(totime0*rf)
         ss=sin(totime0*rf)
      endif

      dtval=RF*0.5/SROT

      beta=beta0+beta1c*cs+beta1s*ss
      cb=cos(beta*pi/180)
      sb=sin(beta*pi/180)

      epsdist=1.5
      fac=0.5
      dmax=10.0

      if (itn .eq. 1) then

        ! Store surface zeta metrics for wall pressure solver...
        ! (This should probably be factored out of this procedure)

        do j = 1,jmax
           do  k = 1,kmax

              zx0(j,k) = zx(j,k,1)
              zy0(j,k) = zy(j,k,1)
              zz0(j,k) = zz(j,k,1)
              zt0(j,k) = -ug(j,k,1)*zx(j,k,1)-vg(j,k,1)*zy(j,k,1)
     <                   -wg(j,k,1)*zz(j,k,1)

           end do
        end do

      end if

      ! Will now sweep from root to tip deforming k-planes...

      do k=1,kmax
         xe=xyzref(1,k)
         ye=xyzref(2,k)
         ze=xyzref(3,k)
         
         do j=1,jmax
            do l=1,lmax 
               
               a1=xg(j,k,l)+ye
               a2=yg(j,k,l)-xe
               a3=zg(j,k,l)-ze

               dist=sqrt(a1**2+a2**2+a3**2)
               dist=max(0.0,dist-epsdist)
               fac1=1.0

               if (idecay.eq.1) then
                   if (dist.gt.dmax) then
                       fac1=0.0
                   else
                       fac1=(1+cos(pi*dist/dmax))*0.5
                   endif
               endif

               if (arf_opt .gt. 0) then
                  e4=defl_dat(4,k,i)
                  e5=defl_dat(5,k,i)
                  e6=defl_dat(6,k,i)
               else
C asitav (slow start of deformation)
                  e4=scal*defl_dat(4,k,i)+(1.0D0-scal)*defl_dat_prev(4,k,i)
                  e5=scal*defl_dat(5,k,i)+(1.0D0-scal)*defl_dat_prev(5,k,i)
                  e6=scal*defl_dat(6,k,i)+(1.0D0-scal)*defl_dat_prev(6,k,i)
               endif

               if (idefl.eq.1) then
                  call rmat_extract(e4,e5,e6,rt)
               elseif (idefl.eq.2) then
                  call rmat_extract_rcas(e4,e5,e6,rt)
               else
                  e4=e4*fac1
                  e5=e5*fac1
                  e6=e6*fac1
                  call rmat_extract_tdu(e4,e5,e6,rt)
               endif

               xe=xyzref(1,k)
               ye=xyzref(2,k)
               ze=xyzref(3,k)
               
               rt11=rt(1,1);rt12=rt(1,2);rt13=rt(1,3);
               rt21=rt(2,1);rt22=rt(2,2);rt23=rt(2,3);
               rt31=rt(3,1);rt32=rt(3,2);rt33=rt(3,3);

               if (arf_opt .gt. 0) then
                  e1=defl_dat(1,k,i)
                  e2=defl_dat(2,k,i)
                  e3=defl_dat(3,k,i)
               else
C asitav (slow start of deformation)
                  e1=scal*defl_dat(1,k,i)+(1.0D0-scal)*defl_dat_prev(1,k,i)
                  e2=scal*defl_dat(2,k,i)+(1.0D0-scal)*defl_dat_prev(2,k,i)
                  e3=scal*defl_dat(3,k,i)+(1.0D0-scal)*defl_dat_prev(3,k,i)
               endif
               e1=e1*rartio*fac1
               e2=e2*rartio*fac1
               e3=e3*rartio*fac1
	
               a11=ss*rt12+cs*rt22
               a12=-ss*rt11-cs*rt21
               a13=-ss*rt13-cs*rt23
               a14=-(-ss*rt11-cs*rt21)*xe-(-ss*rt12-cs*rt22)*ye
     &              -(-ss*rt13-cs*rt23)*ze-ss*xe-cs*ye-ss*e1-cs*e2
               
               a21=-cs*rt12+ss*rt22
               a22=cs*rt11-ss*rt21
               a23=cs*rt13-ss*rt23
               a24=-(cs*rt11-ss*rt21)*xe-(cs*rt12-ss*rt22)*ye
     &              -(cs*rt13-ss*rt23)*ze+cs*xe-ss*ye+cs*e1-ss*e2
               
               a31=-rt32
               a32=rt31
               a33=rt33
               a34=-rt31*xe-rt32*ye-rt33*ze+ze+e3

c..   rigid flapping

               xtmp1=xg(j,k,l)
               ytmp1=(yg(j,k,l)-eflap)*cb-zg(j,k,l)*sb + eflap
               ztmp1=zg(j,k,l)*cb+(yg(j,k,l)-eflap)*sb
               
               xtmp=xtmp1*a11+ytmp1*a12+ztmp1*a13+a14
               ytmp=xtmp1*a21+ytmp1*a22+ztmp1*a23+a24
               ztmp=xtmp1*a31+ytmp1*a32+ztmp1*a33+a34  
               
               if (istep.gt.1) then
                  ug(j,k,l)=(3*xtmp-4*x(j,k,l)+xt2(j,k,l))*dtval
                  vg(j,k,l)=(3*ytmp-4*y(j,k,l)+yt2(j,k,l))*dtval
                  wg(j,k,l)=(3*ztmp-4*z(j,k,l)+zt2(j,k,l))*dtval
     &                 +xlam*fmtip
               else
                  ug(j,k,l)=(xtmp-x(j,k,l))*dtval*2
                  vg(j,k,l)=(ytmp-y(j,k,l))*dtval*2
                  wg(j,k,l)=(ztmp-z(j,k,l))*dtval*2+xlam*fmtip
               endif
               
               ug(j,k,l)=ug(j,k,l)+irotate*(-rf * y(j,k,l)
     $              -fsm1*sin(xpsi))
               vg(j,k,l)=vg(j,k,l)+irotate*( rf * x(j,k,l))


              ! Store prev grid config for next iteration:

               xt2(j,k,l)=x(j,k,l)
               yt2(j,k,l)=y(j,k,l)
               zt2(j,k,l)=z(j,k,l)

               
               ! Update deformed grid configuration:

               x(j,k,l)=xtmp
               y(j,k,l)=ytmp
               z(j,k,l)=ztmp
               
            enddo
         enddo
         
      enddo

      return
      end


c****************************************************************************
c..   Volume of a hexahedron defined by 8 nodes
c****************************************************************************

      function volume (index,xvert,yvert,zvert)

      implicit none
      
      integer index(3)
      real xvert(8),yvert(8),zvert(8)

c..   local variables
      real tet
      real volume
      
      volume=0.0
      volume=volume+tet(xvert,yvert,zvert,1,6,2,3)
      volume=volume+tet(xvert,yvert,zvert,1,5,6,3)
      volume=volume+tet(xvert,yvert,zvert,1,4,5,3)
      volume=volume+tet(xvert,yvert,zvert,4,8,5,7)
      volume=volume+tet(xvert,yvert,zvert,7,4,3,5)
      volume=volume+tet(xvert,yvert,zvert,7,3,6,5)
      

      return
      end

c****************************************************************************
c..  volume of a tetrahedron with nodes j1-j4
c****************************************************************************

      function tet (xvert,yvert,zvert,j1,j2,j3,j4)

      implicit none

      real xvert(8),yvert(8),zvert(8)
      integer j1,j2,j3,j4
      real dx1,dy1,dz1,dx2,dy2,dz2,dx3,dy3,dz3
      
c..   local variables

      real fac1
      real tet

c..   first executable statement

      fac1=1.0/6

      dx1=xvert(j2)-xvert(j1)
      dy1=yvert(j2)-yvert(j1)
      dz1=zvert(j2)-zvert(j1)

      dx2=xvert(j3)-xvert(j1)
      dy2=yvert(j3)-yvert(j1)
      dz2=zvert(j3)-zvert(j1)

      dx3=xvert(j4)-xvert(j1)
      dy3=yvert(j4)-yvert(j1)
      dz3=zvert(j4)-zvert(j1)

      tet= dx1*dy2*dz3-dx1*dz2*dy3-dx2*dy1*dz3+
     &     dx2*dz1*dy3+dx3*dy1*dz2-dx3*dz1*dy2

      tet=fac1*tet

      return
      end

C     ******************************************************************
C     Rotates the TEF of a 3d-wing
C     Changes the base grid (xg,yg,zg)
C     ********************************

      subroutine rigid_flap(x,y,z,psi_rot,psi_rot_old,init)
C     x->xg
C     y->yg
C     z->zg
C
C***********************************************************************
      
      use params_global

      implicit none
      
      real x(jmax,kmax,lmax), y(jmax,kmax,lmax), z(jmax,kmax,lmax )
      real psi_rot,psi_rot_old
      logical init

c ..   local variables

       integer j,k,l,ihar
       real cs,ss
       real xtmp,ytmp,ztmp,xxt,xyt,yxt,yyt,zxt,zyt
       real origin_,psi_new,theta_new,beta_new
       real psi_old,theta_old,beta_old
       real dth,e,cosp1,cosp2,sinp1,sinp2,cosb1,cosb2
       real sinb1,sinb2,sindth,cosdth
       real t11,t12,t13,t14,t21,t22,t23,t24,t31,t32,t33,t34
       real xo,yo,zo,xxo,xyo,xzo,yxo,yyo,yzo,zxo,zyo,zzo

C .....Asitav
       real cs1,ss1,cs2,ss2,theta_fnew,theta_fold
       real::xc,zc,f1,f2,f3,val,yplane,angle,angle0,theta,dist
       real ttef11,ttef12,ttef13,ttef14,ttef31,ttef32,ttef33,ttef34
       real::psitemp_new,thetatemp_new,betatemp_new,psitemp_old,
     $       thetatemp_old,betatemp_old,dthtemp,dtheta_f1,dtheta_f2,
     $       pzc,z1,z2,delx1,delx2,xt1,xt2,zt1,zt2,xsurf(jmax),zsurf(jmax)
C       real,parameter::pxc=0.80, pzc=0.0, ys=1.44, ye=2.31, rmax=2.5, rmin=0.2,
C     $                 dely = 0.2
       real,parameter::rmax=2.5, rmin=0.2

       real, save :: theta_prev

       pi = 4.0*atan(1.0)

       theta_fnew=theta_f0
       theta_fold=theta_f0
       
       do ihar=1,nharmflap
          theta_fnew=theta_fnew+ampFlap(ihar)
     &         *cos(ihar*psi_rot-phiFlap(ihar))
          theta_fold=theta_fold+ampFlap(ihar)
     &         *cos(ihar*psi_rot_old-phiFlap(ihar))
       enddo

c$$$       theta_flap= theta_fnew
c$$$       write(6,*) "theta_flap ",theta_flap*180/pi
c$$$       
c$$$       dtheta_f2 = theta_fnew - theta_f0         !default = 0.0
c$$$       dtheta_f1 = -(theta_fold - theta_f0)

       if (init) then
          dtheta_f2=theta_fnew
          dtheta_f1=0.0D0
       else
          dtheta_f2=theta_fnew
          dtheta_f1=-theta_prev
       endif
       write(6,'(A,F12.5,F12.5)') "TEF angles: ", dtheta_f2*180./pi,
     $      dtheta_f1*180.0/pi
       
       do k = 1,kmax
C.... test
        xt1 = x(jtail1,k,1); zt1 = z(jtail1,k,1)
        xt2 = x(jtail2,k,1); zt2 = z(jtail2,k,1)
        do j=1,jmax
         xsurf(j) = x(j,k,1); zsurf(j) = z(j,k,1);
        end do
     
        xc = (pxc*(x(jtail1,k,1) - x(jle,k,1)) + x(jle,k,1))
C       locate hinge
C       ------------
       hinge: do j=jtail1,jtail2
        delx1 = x(j,k,1)-xc
        delx2 = x(j+1,k,1)-xc
        if( (delx1*delx2 < 0.0) .and. (j < (jtail1+jtail2)/2 ) ) then
         z1 = 0.5*(z(j,k,1)+z(j+1,k,1))
        else if( (delx1*delx2 < 0.0) .and. (j > (jtail1+jtail2)/2 ) ) then
         z2 = 0.5*(z(j,k,1)+z(j+1,k,1))
         exit hinge
        else
         cycle hinge
        end if
       end do hinge
c
c check jaina
c
        z1=(z1+z2)*0.5 
        z2=z1

C     ------------
   
        yplane = y((jtail1+(jtail2-jtail1)/4),k,1 ) !semi-chord section
       
        if(yplane < flap_ys-dely .or. yplane > flap_ye+dely ) then
         f3 =0;
        elseif(yplane >flap_ys+dely .and. yplane <flap_ye-dely ) then
         f3=1.0;
        elseif(yplane >= flap_ys-dely .and. yplane <= flap_ys+dely) then
         val = yplane-flap_ys+dely;
         val = pi*val/(2.0*2.0*dely)
         f3 = sin(val);
        elseif(yplane >= flap_ye-dely .and. yplane <= flap_ye+dely) then
         val = yplane-flap_ye+dely;
         val = pi*val/(2.0*2.0*dely);
         f3 = cos(val);
        end if
  
C        xc = x(jtail1+ (jtail2-jtail1)/8,k,1);
C        zc = z(jtail1+ (jtail2-jtail1)/8,k,1);

  
        do l = 1,lmax
         do j = 1,jmax
C........ Point of rotation
          if(j<(jtail1+jtail2)/2) then
           zc = z1 
          else if(j>(jtail1+jtail2)/2) then
           zc = z2 
          end if
 
C          dist=sqrt((x(j,k,l)-x(j,k,1))**2+(z(j,k,l)-z(j,k,1))**2);
          dist=sqrt((x(j,k,l)-xsurf(j))**2+(z(j,k,l)-zsurf(j))**2);
          if (j > jtail2) then
C           dist=sqrt((x(j,k,l)-xt2)**2+(z(j,k,l)-zt2)**2);
           dist=sqrt((x(j,k,l)-xsurf(jtail2))**2+(z(j,k,l)-zsurf(jtail2))**2);
          else if (j < jtail1) then
C           dist=sqrt((x(j,k,l)-xt1)**2+(z(j,k,l)-zt1)**2);
           dist=sqrt((x(j,k,l)-xsurf(jtail1))**2+(z(j,k,l)-zsurf(jtail1))**2);
          end if
       
          if (dist > rmax) then
           f1=0.0;
          else if (dist < rmin) then
           f1=1.0;
          else
           val=dist-rmin;
           val=val/(rmax-rmin)*pi;
           f1=(1+cos(val))*0.5;
          end if
       
          angle =atan2( (z(j,k,l)-zc), (x(j,k,l)-xc) );
          angle0=atan2( (z(jtail1,k,1)-zc), (x(jtail1,k,1)-xc) );
          angle = angle - angle0;


          if (abs(angle) > pi*0.5+dela) then
           f2=0.0
          elseif (abs(angle)<pi*0.5-dela)  then
           f2=1.0
          else
           val=abs(angle)-pi*0.5+dela
           val=pi*val/(2.0*dela)
           f2=(1+cos(val))*0.5
          end if

          cs=cos( (dtheta_f1+dtheta_f2) *f3*f2*f1)
          ss=sin( (dtheta_f1+dtheta_f2) *f3*f2*f1)
          
C  .......The t-matrix for te flap
          ttef11 = cs
          ttef12 = 0.0
          ttef13 = ss
          ttef14 = xc*(1. - cs) - zc*ss
 
          ttef31 = -ss
          ttef32 = 0.0
          ttef33 = cs
          ttef34 = xc*ss + zc*(1-cs)
C  .......move grid
          xo=x(j,k,l) !- xc
          yo=y(j,k,l)
          zo=z(j,k,l) !- zc
          x(j,k,l) = xo*ttef11+zo*ttef13 + ttef14 
          y(j,k,l) = yo
          z(j,k,l) = xo*ttef31+zo*ttef33 + ttef34
  
C with no blade rotation...
C *****************************************************
C          cs=cos((theta_fnew - theta_fold)*f3*f2*f1)
C          ss=sin((theta_fnew - theta_fold)*f3*f2*f1)
C
C          xo=x(j,k,l) - xc
C          yo=y(j,k,l)
C          zo=z(j,k,l) - zc
C          x(j,k,l)=xo*cs-zo*ss+xc
C          z(j,k,l)=xo*ss+zo*cs+zc
C          y(j,k,l) = yo
C *****************************************************

C          if(x(j,k,l) .ge. 10.0 .or. z(j,k,l) .ge. 20.0 ) then
C           print*,"x ",x(j,k,l)," z ",z(j,k,l)
C          end if

         end do
        end do

       end do

       theta_prev=theta_fnew
      end subroutine rigid_flap 


C*******************************************************************************
C....Added by Asitav 
C*******************************************************************************
      subroutine rotate_rigid_flap(srot,x,y,z,xx,xy,xz,ug,yx,yy,yz,vg,
     &                 zx,zy,zz,wg,zx0,zy0,zz0,zt0,istart)
C
C***********************************************************************
      
      use params_global

      implicit none
      
      real xx(jmax,kmax,lmax),xy(jmax,kmax,lmax),xz(jmax,kmax,lmax),
     $     yx(jmax,kmax,lmax),yy(jmax,kmax,lmax),yz(jmax,kmax,lmax),
     $     zx(jmax,kmax,lmax),zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)
      real ug(jmax,kmax,lmax), vg(jmax,kmax,lmax),wg(jmax,kmax,lmax)
      real x(jmax,kmax,lmax), y(jmax,kmax,lmax), z(jmax,kmax,lmax )
      real zx0(jmax,kmax),zy0(jmax,kmax),zz0(jmax,kmax),zt0(jmax,kmax)
      integer istart

c..   local variables

      integer j,k,l
      real cs,ss,srot
      real xtmp,ytmp,ztmp,xxt,xyt,yxt,yyt,zxt,zyt
      real origin_,psi_new,theta_new,beta_new
      real psi_old,theta_old,beta_old
      real dth,e,cosp1,cosp2,sinp1,sinp2,cosb1,cosb2
      real sinb1,sinb2,sindth,cosdth
      real t11,t12,t13,t14,t21,t22,t23,t24,t31,t32,t33,t34
      real xo,yo,zo,xxo,xyo,xzo,yxo,yyo,yzo,zxo,zyo,zzo

C.....Asitav
      real,allocatable::xold(:,:,:),yold(:,:,:),zold(:,:,:)
      real theta_fnew,theta_fold
      real::xc,zc,f1,f2,f3,val,yplane,angle,angle0,theta,dist
      real ttef11,ttef12,ttef13,ttef14,ttef31,ttef32,ttef33,ttef34
      real::psitemp_new,thetatemp_new,betatemp_new,psitemp_old,thetatemp_old,
     $      betatemp_old,dthtemp,dtheta_f1,dtheta_f2,pzc,z1,z2,delx1,delx2
!      real,parameter::pxc=0.80, pzc=0.0, ys=1.44, ye=2.31, rmax=2.5, rmin=0.2,
!     $                dely = 0.2
      real,parameter::rmax=2.5, rmin=0.2

      pi = 4.0*atan(1.0)

c..   first executable statement
      
C.....Total displacement from (psi_1,theta_1,beta_1,theta_teflap_1) to
C     (psi_2,theta_2,beta_2,theta_teflap_2), 2=new, 1=old
C     ----------------------------------------------------------------
      origin_ = 0.0
      psi_new=totime*rf

C.....psi_1 co-ord
      allocate(xold(jmax,kmax,lmax),yold(jmax,kmax,lmax),zold(jmax,kmax,lmax))

      xold = x  
      yold = y
      zold = z
      
      theta_new=theta0+theta1c*cos(psi_new)+theta1s*sin(psi_new)
      beta_new=beta0+beta1c*cos(psi_new)+beta1s*sin(psi_new)

      if(istart.eq.1) then
         psi_old=0.
         theta_old=0.
         beta_old=0.
      else
         psi_old=psi_new-srot
         theta_old=theta0+theta1c*cos(psi_old)+theta1s*sin(psi_old)
         beta_old=beta0+beta1c*cos(psi_old)+beta1s*sin(psi_old)
      endif
C     ----------------------------------------------------------------

C      dth=theta_new-theta_old
C      e=eflap

C     psi_1 to psi_0 (=0.0 azimuth is zero)
C     *************************************

      origin_ = 0.0
      psitemp_new= 0.0 !azimuth = 0.0
      
      thetatemp_new=theta0+theta1c*cos(psitemp_new)+theta1s*sin(psitemp_new)
      betatemp_new=beta0+beta1c*cos(psitemp_new)+beta1s*sin(psitemp_new)

      psitemp_old=psi_old !psi_1
      thetatemp_old=theta0+theta1c*cos(psitemp_old)+theta1s*sin(psitemp_old)
      betatemp_old=beta0+beta1c*cos(psitemp_old)+beta1s*sin(psitemp_old)

      dthtemp=thetatemp_new-thetatemp_old
      e=eflap

      if (fmtip.eq.0) then
         cosp1=1.
         cosp2=1.
         sinp1=0.
         sinp2=0.
      else
         cosp1=cos(psitemp_old)
         cosp2=cos(psitemp_new)
         sinp1=sin(psitemp_old)
         sinp2=sin(psitemp_new)
      endif

      cosb1=cos(betatemp_old)
      cosb2=cos(betatemp_new)
      sinb1=sin(betatemp_old)
      sinb2=sin(betatemp_new)
      sindth=sin(dthtemp)
      cosdth=cos(dthtemp)
 
      t11=cosp2*(cosdth*cosp1+sindth*sinb1*sinp1)-sinp2
     &     *(-cosb2*cosb1*sinp1-sinb2*(-sindth*cosp1
     &     +cosdth*sinb1*sinp1))
      t12=cosp2*(cosdth*sinp1-sindth*sinb1*cosp1)-sinp2
     &     *(cosb2*cosb1*cosp1-sinb2*(-sindth*sinp1
     &     -cosdth*sinb1*cosp1))
      t13=cosp2*sindth*cosb1-sinp2*(cosb2*sinb1-sinb2
     &     *cosdth*cosb1)
      t14=cosp2*(-origin_*cosdth+sindth*sinb1*e+origin_)-sinp2
     &     *(-cosb2*cosb1*e-sinb2*(origin_*sindth
     &     +cosdth*sinb1*e)+e)

      t21=sinp2*(cosdth*cosp1+sindth*sinb1*sinp1)+cosp2
     &     *(-cosb2*cosb1*sinp1-sinb2*(-sindth*cosp1
     &     +cosdth*sinb1*sinp1))
      t22=sinp2*(cosdth*sinp1-sindth*sinb1*cosp1)+cosp2
     &     *(cosb2*cosb1*cosp1-sinb2*(-sindth*sinp1-cosdth
     &     *sinb1*cosp1))
      t23=sinp2*sindth*cosb1+cosp2*(cosb2*sinb1-sinb2
     &     *cosdth*cosb1)
      t24=sinp2*(-origin_*cosdth+sindth*sinb1*e+origin_)+cosp2
     &     *(-cosb2*cosb1*e-sinb2*(origin_*sindth+cosdth
     &     *sinb1*e)+e)

      t31=-sinb2*cosb1*sinp1+cosb2*(-sindth*cosp1
     &     +cosdth*sinb1*sinp1)
      t32=sinb2*cosb1*cosp1+cosb2*(-sindth*sinp1
     &     -cosdth*sinb1*cosp1)
      t33=sinb2*sinb1+cosb2*cosdth*cosb1
      t34=-sinb2*cosb1*e+cosb2*(origin_*sindth+cosdth*sinb1*e)
c
C      do j = 1,jmax
C         do k = 1,kmax
C            zx0(j,k) = zx(j,k,1)
C            zy0(j,k) = zy(j,k,1)
C            zz0(j,k) = zz(j,k,1)
C            zt0(j,k) = -ug(j,k,1)*zx(j,k,1)-vg(j,k,1)*zy(j,k,1)
C     <           -wg(j,k,1)*zz(j,k,1)
C         enddo
C      enddo
 
c..warning!!!!!!!!!the following is a temporary fix
c..till the exact grid velocities are calculated
c..by differentiating the t matrix
      
         do l = 1,lmax
            do k = 1,kmax
               do j = 1,jmax
c..  only grid
                  xo=x(j,k,l)
                  yo=y(j,k,l)
                  zo=z(j,k,l)
                  x(j,k,l)=t11*xo+t12*yo+t13*zo+t14
                  y(j,k,l)=t21*xo+t22*yo+t23*zo+t24
                  z(j,k,l)=t31*xo+t32*yo+t33*zo+t34
               enddo
            enddo
         enddo

C    TE-Flap displacement at psi = 0.0 
c    *********************************
C    The T-matrix of TE-flap involves un-flapping TE flap by current theta_teflap  
C    & flapping back to new theta_teflap from zero TE-flap position

      
      theta_fnew=angmax_tef*sin(rf_tef*istep*dt)    !theta_f2
      theta_fold=angmax_tef*sin(rf_tef*(istep-1)*dt)!theta_f1
      theta_flap=theta_fnew
      
      dtheta_f2 = theta_fnew - theta_f0         !default theta_f0= 0.0
      dtheta_f1 = -(theta_fold - theta_f0)

!      dela=30*pi/180
    
      do k = 1,kmax


       xc = (pxc*(x(jtail1,k,1) - x(jle,k,1)) + x(jle,k,1))
C     locate hinge
C     ------------
       hinge: do j=jtail1,jtail2
        delx1 = x(j,k,1)- xc !(pxc*(x(jtail1,k,1) - x(jle,k,1)) + x(jle,k,1))
        delx2 = x(j+1,k,1)- xc !(pxc*(x(jtail1,k,1) - x(jle,k,1)) + x(jle,k,1))
        if( (delx1*delx2 < 0.0) .and. (j < (jtail1+jtail2)/2 ) ) then
         z1 = 0.5*(z(j,k,1)+z(j+1,k,1))
        else if( (delx1*delx2 < 0.0) .and. (j > (jtail1+jtail2)/2 ) ) then
         z2 = 0.5*(z(j,k,1)+z(j+1,k,1))
         exit hinge
        else
         cycle hinge
        end if
       end do hinge
C     ------------

       yplane = y((jtail1+(jtail2-jtail1)/4),k,1 ) !semi-chord section
      
       if(yplane < flap_ys-dely .or. yplane > flap_ye+dely ) then
        f3 =0;
       elseif(yplane >flap_ys+dely .and. yplane <flap_ye-dely ) then
        f3=1.0;
       elseif(yplane >= flap_ys-dely .and. yplane <= flap_ys+dely) then
        val = yplane-flap_ys+dely;
        val = pi*val/(2.0*2.0*dely)
        f3 = sin(val);
       elseif(yplane >= flap_ye-dely .and. yplane <= flap_ye+dely) then
        val = yplane-flap_ye+dely;
        val = pi*val/(2.0*2.0*dely);
        f3 = cos(val);
       end if

       do l = 1,lmax
        do j = 1,jmax
        
C........Point of rotation
         if(j<(jtail1+jtail2)/2) then
          zc = z1 
         else if(j>(jtail1+jtail2)/2) then
          zc = z2 
         end if

         dist=sqrt((x(j,k,l)-x(j,k,1))**2+(z(j,k,l)-z(j,k,1))**2);
         if (j > jtail2) then
          dist=sqrt((x(j,k,l)-x(jtail2,k,1))**2+(z(j,k,l)-z(jtail2,k,1))**2);
         else if (j < jtail1) then
          dist=sqrt((x(j,k,l)-x(jtail1,k,1))**2+(z(j,k,l)-z(jtail1,k,1))**2);
         end if
      
         if (dist > rmax) then
          f1=0.0;
         else if (dist < rmin) then
          f1=1.0;
         else
          val=dist-rmin;
          val=val/(rmax-rmin)*pi;
          f1=(1+cos(val))*0.5;
         end if
      
         angle =atan2( (z(j,k,l)-zc), (x(j,k,l)-xc) );
         angle0=atan2( (z(jtail1,k,1)-zc), (x(jtail1,k,1)-xc) );
         angle = angle - angle0;
      
         if (abs(angle) > pi*0.5+dela) then
          f2=0.0
         elseif (abs(angle)<pi*0.5-dela)  then
          f2=1.0
         else
          val=abs(angle)-pi*0.5+dela
          val=pi*val/(2.0*dela)
          f2=(1+cos(val))*0.5
         end if
  
          cs=cos( (dtheta_f1+dtheta_f2) *f3*f2*f1)
          ss=sin( (dtheta_f1+dtheta_f2) *f3*f2*f1)
          
C  .......The t-matrix for te flap
          ttef11 = cs
          ttef12 = 0.0
          ttef13 = ss
          ttef14 = xc*(1. - cs) - zc*ss
 
          ttef31 = -ss
          ttef32 = 0.0
          ttef33 = cs
          ttef34 = xc*ss + zc*(1-cs)
C  .......move grid
          xo=x(j,k,l) !- xc
          yo=y(j,k,l)
          zo=z(j,k,l) !- zc
          x(j,k,l) = xo*ttef11+zo*ttef13 + ttef14 
          y(j,k,l) = yo
          z(j,k,l) = xo*ttef31+zo*ttef33 + ttef34
 
        end do
       end do

      end do

C     psi_0 to psi_2 
C     **************

      origin_ = 0.0
      psitemp_new = psi_new ! psi_2
      
      thetatemp_new=theta0+theta1c*cos(psitemp_new)+theta1s*sin(psitemp_new)
      betatemp_new=beta0+beta1c*cos(psitemp_new)+beta1s*sin(psitemp_new)

      psitemp_old = 0.0!psi = 0.0, zero azimuth
      thetatemp_old=theta0+theta1c*cos(psitemp_old)+theta1s*sin(psitemp_old)
      betatemp_old=beta0+beta1c*cos(psitemp_old)+beta1s*sin(psitemp_old)

      dthtemp=thetatemp_new-thetatemp_old
      e=eflap

      if (fmtip.eq.0) then
         cosp1=1.
         cosp2=1.
         sinp1=0.
         sinp2=0.
      else
         cosp1=cos(psitemp_old)
         cosp2=cos(psitemp_new)
         sinp1=sin(psitemp_old)
         sinp2=sin(psitemp_new)
      endif

      cosb1=cos(betatemp_old)
      cosb2=cos(betatemp_new)
      sinb1=sin(betatemp_old)
      sinb2=sin(betatemp_new)
      sindth=sin(dthtemp)
      cosdth=cos(dthtemp)
 
      t11=cosp2*(cosdth*cosp1+sindth*sinb1*sinp1)-sinp2
     &     *(-cosb2*cosb1*sinp1-sinb2*(-sindth*cosp1
     &     +cosdth*sinb1*sinp1))
      t12=cosp2*(cosdth*sinp1-sindth*sinb1*cosp1)-sinp2
     &     *(cosb2*cosb1*cosp1-sinb2*(-sindth*sinp1
     &     -cosdth*sinb1*cosp1))
      t13=cosp2*sindth*cosb1-sinp2*(cosb2*sinb1-sinb2
     &     *cosdth*cosb1)
      t14=cosp2*(-origin_*cosdth+sindth*sinb1*e+origin_)-sinp2
     &     *(-cosb2*cosb1*e-sinb2*(origin_*sindth
     &     +cosdth*sinb1*e)+e)

      t21=sinp2*(cosdth*cosp1+sindth*sinb1*sinp1)+cosp2
     &     *(-cosb2*cosb1*sinp1-sinb2*(-sindth*cosp1
     &     +cosdth*sinb1*sinp1))
      t22=sinp2*(cosdth*sinp1-sindth*sinb1*cosp1)+cosp2
     &     *(cosb2*cosb1*cosp1-sinb2*(-sindth*sinp1-cosdth
     &     *sinb1*cosp1))
      t23=sinp2*sindth*cosb1+cosp2*(cosb2*sinb1-sinb2
     &     *cosdth*cosb1)
      t24=sinp2*(-origin_*cosdth+sindth*sinb1*e+origin_)+cosp2
     &     *(-cosb2*cosb1*e-sinb2*(origin_*sindth+cosdth
     &     *sinb1*e)+e)

      t31=-sinb2*cosb1*sinp1+cosb2*(-sindth*cosp1
     &     +cosdth*sinb1*sinp1)
      t32=sinb2*cosb1*cosp1+cosb2*(-sindth*sinp1
     &     -cosdth*sinb1*cosp1)
      t33=sinb2*sinb1+cosb2*cosdth*cosb1
      t34=-sinb2*cosb1*e+cosb2*(origin_*sindth+cosdth*sinb1*e)
c .....?????.............
C      do j = 1,jmax
C         do k = 1,kmax
C            zx0(j,k) = zx(j,k,1)
C            zy0(j,k) = zy(j,k,1)
C            zz0(j,k) = zz(j,k,1)
C            zt0(j,k) = -ug(j,k,1)*zx(j,k,1)-vg(j,k,1)*zy(j,k,1)
C     <           -wg(j,k,1)*zz(j,k,1)
C         enddo
C      enddo
 
      if(istart.ne.1) then
         do l = 1,lmax
            do k = 1,kmax
               do j = 1,jmax
c..   grid
                  xo=x(j,k,l)
                  yo=y(j,k,l)
                  zo=z(j,k,l)
                  x(j,k,l)=t11*xo+t12*yo+t13*zo+t14
                  y(j,k,l)=t21*xo+t22*yo+t23*zo+t24
                  z(j,k,l)=t31*xo+t32*yo+t33*zo+t34
c..   grid velocities (x(psi_2) - x(psi_1)) / dt
                  ug(j,k,l)=(x(j,k,l)-xold(j,k,l))/dt
                  vg(j,k,l)=(y(j,k,l)-yold(j,k,l))/dt
                  wg(j,k,l)=(z(j,k,l)-zold(j,k,l))/dt
               enddo
            enddo
         enddo
      else
         do l = 1,lmax
            do k = 1,kmax
               do j = 1,jmax
c..   grid
                  xo=x(j,k,l)
                  yo=y(j,k,l)
                  zo=z(j,k,l)
                  x(j,k,l)=t11*xo+t12*yo+t13*zo+t14
                  y(j,k,l)=t21*xo+t22*yo+t23*zo+t24
                  z(j,k,l)=t31*xo+t32*yo+t33*zo+t34
               enddo
            enddo
         enddo
      endif

      deallocate(xold,yold,zold)
      return
      end subroutine rotate_rigid_flap

      subroutine rigid_flapoh(x,y,z,psi_rot,psi_rot_old,init)
C
C     x->xg
C     y->yg
C     z->zg
C
C     Modifies only the base grid w/o having to worry about storing the
C     old/new metrices.
C
C     *****Use it when deformation subroutine is called later*****
C
C     Written by Asitav Oct 1st, 2007
C***********************************************************************
      
      use params_global

      implicit none
      
      real x(jmax,kmax,lmax), y(jmax,kmax,lmax), z(jmax,kmax,lmax )
      real psi_rot,psi_rot_old
      logical init

c ..   local variables

      integer j,k,l,ihar
      real cs,ss
      real xtmp,ytmp,ztmp,xxt,xyt,yxt,yyt,zxt,zyt
      real origin_,psi_new,theta_new,beta_new
      real psi_old,theta_old,beta_old
      real dth,e,cosp1,cosp2,sinp1,sinp2,cosb1,cosb2
      real sinb1,sinb2,sindth,cosdth
      real t11,t12,t13,t14,t21,t22,t23,t24,t31,t32,t33,t34
      real xo,yo,zo,xxo,xyo,xzo,yxo,yyo,yzo,zxo,zyo,zzo

C ....Asitav
      real cs1,ss1,cs2,ss2,theta_fnew,theta_fold
      real f1,f2,f3,val,yplane,angle,angle0,theta,dtheta,dist,xc,zc
      real ttef11,ttef12,ttef13,ttef14,ttef31,ttef32,ttef33,ttef34
      real psitemp_new,thetatemp_new,betatemp_new,psitemp_old,
     $      thetatemp_old,betatemp_old,dthtemp,dtheta_f1,dtheta_f2,
     $      pzc,z1,z2,delx1,delx2,xt1,xt2,zt1,zt2
      real xsurf(jmax),zsurf(jmax)
      real,parameter ::  rmax=2.5, rmin=0.2

C.... Asitav (oh stuff)
      integer j1,j2,j3,j4
      real ang_temp,xdummy,zdummy
      real aa,aa_l,aa_u,angs,ange,xfil_s,xfil_e,
     $     xfil_sm,xfil_em,zfil_sm,zfil_em,
     $     xfil_el,xfil_eu,zfil_el,zfil_eu,xfil_m,zfil_m


      pi = 4.0*atan(1.0)

      theta_fnew=0.
      theta_fold=0.
      
      do ihar=1,nharmflap
         theta_fnew=theta_fnew+ampFlap(ihar)
     &        *cos(ihar*psi_rot-phiFlap(ihar))
         theta_fold=theta_fold+ampFlap(ihar)
     &        *cos(ihar*psi_rot_old-phiFlap(ihar))
      enddo

      dtheta = theta_fnew-theta_fold

      if (init) dtheta=theta_fnew

      j1 = joh1
      j2 = joh1+noh      !joh2
      j3 = jmax-(j2-1)   !joh3
      j4 = jmax-(joh1-1) !joh4

      do k = 1,kmax
       do j=1,jmax
        xsurf(j) = x(j,k,1); zsurf(j) = z(j,k,1)
       end do

       !define boundary of filler
       !xfil_e = xflph - (1.0-xflph)*flpovh
       !xfil_s = xfil_e - cgap
       !alternate definition
       xfil_s = xflph - (1.0-xflph)*flpovh -0.10*cgap
       xfil_e = xfil_s + 0.90*cgap

       !assumes chord=1.0
       hinge: do j=jtail1,jtail2
        delx1 = x(j,k,1)-xflph
        delx2 = x(j+1,k,1)-xflph
        if( (delx1*delx2 < 0.0) .and. (j < (jtail1+jtail2)/2 ) ) then
         z1 = 0.5*(z(j,k,1)+z(j+1,k,1))
        else if( (delx1*delx2 < 0.0) .and. (j > (jtail1+jtail2)/2 ) ) then
         z2 = 0.5*(z(j,k,1)+z(j+1,k,1))
         exit hinge
        else
         cycle hinge
        end if
       end do hinge
       
       zflph = (z1+z2)/2.0
        
       !*******************************************
       !following executions change every time step
       !*******************************************

       !find the angle decay fn 'aa'   (to be done every time step)
       !---------------------------- * ----------------------------

       xfil_sm = (xsurf(j2)+xsurf(j3))/2.
       zfil_sm = (zsurf(j2)+zsurf(j3))/2.
       xfil_em = (xsurf(j1)+xsurf(j4))/2.
       zfil_em = (zsurf(j1)+zsurf(j4))/2.

       xfil_m = (xfil_sm+xfil_em)/2.
       zfil_m = (zfil_sm+zfil_em)/2.

       xfil_eu = (xsurf(j2)*2+xsurf(j3)*7)/9.
       zfil_eu = (zsurf(j2)*2+zsurf(j3)*7)/9.
       xfil_el = (xsurf(j2)*7+xsurf(j3)*2)/9 
       zfil_el = (zsurf(j2)*7+zsurf(j3)*2)/9.

       !upper 'aa'
       angs = atan2( (zsurf(j3)-zfil_eu),(xsurf(j3)-xfil_eu) )
       ange = atan2( (zsurf(j4)-zfil_eu),(xsurf(j4)-xfil_eu) )
       aa_u = abs(ange - angs)
       !lower 'aa'
       angs = atan2( (zsurf(j2)-zfil_el),(xsurf(j2)-xfil_el))
       ange = atan2( (zsurf(j1)-zfil_el),(xsurf(j1)-xfil_el))
       aa_l = abs(ange - angs)
C     ------------
   
       yplane = y((jtail1+(jtail2-jtail1)/4),k,1 ) !semi-chord section
       
       if(yplane < flap_ys-dely .or. yplane > flap_ye+dely ) then
        f3 =0;
       elseif(yplane >flap_ys+dely .and. yplane <flap_ye-dely ) then
        f3=1.0;
       elseif(yplane >= flap_ys-dely .and. yplane <= flap_ys+dely) then
        val = yplane-flap_ys+dely;
        val = pi*val/(2.0*2.0*dely)
        f3 = sin(val);
       elseif(yplane >= flap_ye-dely .and. yplane <= flap_ye+dely) then
        val = yplane-flap_ye+dely;
        val = pi*val/(2.0*2.0*dely);
        f3 = cos(val);
       end if
 
       !====================================================================
       !treat the blade mesh
       !********************
       !====================================================================
       do l=1,lmax
        do j=1,jmax
        
         dist=sqrt((x(j,k,l)-xsurf(j))**2+(z(j,k,l)-zsurf(j))**2);
         if (j > jtail2) then
          dist=sqrt((x(j,k,l)-xsurf(jtail2))**2+(z(j,k,l)-zsurf(jtail2))**2);
         elseif (j < jtail1) then
          dist=sqrt((x(j,k,l)-xsurf(jtail1))**2+(z(j,k,l)-zsurf(jtail1))**2);
         end if

         if (dist > rmax) then
          f1=0.0;
         elseif (dist < rmin) then
          f1=1.0;
         else
          val=dist-rmin;
          val=val/(rmax-rmin)*pi;
          f1=(1+cos(val))*0.5;
         end if

         if(j < (jtail1+jtail2)/2) then
          aa = aa_l
          xdummy = xfil_el
          zdummy = zfil_el
         else if(j > (jtail1+jtail2)/2) then
          aa = aa_u
          xdummy = xfil_eu
          zdummy = zfil_eu
         end if

         xc = xdummy; zc = zdummy
         angle =atan2( (z(j,k,l)-zc), (x(j,k,l)-xc) );
         angle0=atan2( (zsurf(jtail1)-zc), (xsurf(jtail1)-xc));
         angle = angle - angle0;

         ang_temp = pi*0.5
         !if (abs(angle) > ang_temp+aa) then
         if (abs(angle) > ang_temp) then
          f2=0.
         elseif (abs(angle)<ang_temp-aa) then
          f2=1.
         else
          val=abs(angle)-ang_temp+aa
          val=pi*val/aa
          f2=(1+cos(val))*0.5
         end if
         
         cs = cos(f1*f2*f3*dtheta)
         ss = sin(f1*f2*f3*dtheta)

C  ..... The t-matrix for te flap
         ttef11 = cs
         ttef12 = 0.0
         ttef13 = ss
         ttef14 = xflph*(1. - cs) - zflph*ss
 
         ttef31 = -ss
         ttef32 = 0.0
         ttef33 = cs
         ttef34 = xflph*ss + zflph*(1-cs)
C  ..... move grid
         xo=x(j,k,l) 
         yo=y(j,k,l)
         zo=z(j,k,l) 
         x(j,k,l) = xo*ttef11+zo*ttef13 + ttef14 
         y(j,k,l) = yo
         z(j,k,l) = xo*ttef31+zo*ttef33 + ttef34
 
        end do
       end do !end blade mesh treatment

      end do !k=1,kmax
      end subroutine rigid_flapoh





