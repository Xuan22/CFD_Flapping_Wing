function line_eqn(x1,y1,num,den,x,y)
  implicit none

  real :: line_eqn
  real, intent(in) :: x1, y1, num,den, x, y
  
  line_eqn=(y-y1)*den-num*(x-x1)
end function line_eqn

integer function locate_yindex(ypos)
  use hqptef
  implicit none

  real, intent(in) :: ypos
  integer :: ki
  integer :: k
  
  if ((ypos .lt. yflap(1)) .or. (ypos .gt. yflap(kfmax))) then
     ki=-1
  else
     do k=1,kfmax
        if (yflap(k) .ge. ypos) then
           ki=k
	   exit
        end if
     end do
  end if
  locate_yindex=ki
end function locate_yindex

subroutine init_flap_params
  use ioutils
  use hqptef
  implicit none

  integer :: un

  un=open_file('tef.inp',form='formatted',status='old')
  read(un,*) kfmax

  allocate(yflap(kfmax),zupper(kfmax),zlower(kfmax),zte(kfmax))
  allocate(xupper(kfmax),xlower(kfmax),xte(kfmax))
  allocate(phiupmx(kfmax),phiupmi(kfmax))
  allocate(philowmx(kfmax),philowmi(kfmax))
  allocate(flphx(kfmax),flphz(kfmax))
  
  read(un,*) yflap, xupper, zupper, xlower, zlower, flphx, flphz,&
       & xte, zte, phiupmx, phiupmi, philowmx, philowmi

  close(un)

end subroutine init_flap_params

subroutine flap_decay_funcs(xco, yco, zco, decay, kp)
  use params_global, only: flap_ys, flap_ye, dely
  use hqptef
  implicit none

  real :: xco, yco, zco
  real :: decay
  integer :: kp
  
  real :: xfac,yfac,zfac
  real :: ymin, ymax, ypos, xmid, zmid, &
       &xup, zup, xlow, zlow, num, den, num1, den1
  real :: fx0, fx1, fx2,xfac1,xfac2,sl1,sl2
  real :: fz0, fz1,dist,angle
  real :: angmax, angmin
  real :: rmin, rmax, pi

  integer :: locate_yindex
  real :: line_eqn

  pi=4.0D0*atan(1.0D0)
  rmin=0.2D0
  rmax=2.0D0
  
  xfac=1.0D0
  yfac=1.0D0
  zfac=1.0D0
  
  kp=locate_yindex(yco)
  ypos=yco
  span: if (kp .lt. 0) then 
     !point not in the flap region spanwise location
     yfac=0.0D0
  else
     xlow=xlower(kp)
     zlow=zlower(kp)
     xup=xupper(kp)
     zup=zupper(kp)
     xmid=0.5D0*(xlow+xup)
     zmid=0.5D0*(zlow+zup)
     sl1=atan2((zte(kp)-zmid),xte(kp)-xmid)
     
     num=zup-zlow
     den=xup-xlow
     num1=zte(kp)-zmid
     den1=xte(kp)-xmid
     fx0=line_eqn(xup,zup,num,den,0.0D0,0.0D0)
     
     fx1=line_eqn(xup,zup,num,den,xco,zco)
     fx2=line_eqn(xte(kp),zte(kp),num,den,xco,zco)
     
     if (abs(fx1) .le. 1D-5) fx1=0.0D0
     chord: if ((fx0*fx1 .ge. 0.0D0) .and. (fx1*fx2 .ge. 0.0D0)) then
        !point in the blade section ahead of TEF LE
        xfac=0.0D0
     else
        if (fx1*fx2 .le. 0.0D0) then 
           !Within the flap section chordwise 
           dist=abs(zco-zte(kp))
        else
           dist=sqrt((xco-xte(kp))**2+(zco-zte(kp))**2)
        end if
        if (dist .gt. rmax) then
           zfac=0.0D0
        else if (dist .ge. rmin) then
           zfac=pi*(dist-rmin)/(rmax-rmin)
           zfac=0.5D0*(1.0D0+cos(zfac))
        end if
        
        if ((ypos .ge. yflap(1)) .and. (ypos .le. flap_ys+dely)) then
           yfac=ypos-flap_ys+dely
           yfac=0.5D0*pi*yfac/dely
           yfac=0.5D0*(1.0D0-cos(yfac))
        else if((ypos .ge. flap_ye-dely) .and. (ypos .le. yflap(kfmax))) then
           yfac=ypos-flap_ye+dely
           yfac=0.5D0*pi*yfac/dely
           yfac=0.5D0*(1.0D0+cos(yfac))
        end if
        
        fz0=line_eqn(xmid,zmid,num1,den1,xup,zup)
        fz1=line_eqn(xmid,zmid,num1,den1,xco,zco)
        sl2=atan2((zco-zmid),(xco-xmid))
        angle=abs(sl2-sl1)
        
        if (fz0*fz1 .gt. 0.0D0) then
           angmax=phiupmx(kp)
           angmin=phiupmi(kp)
        else
           angmax=philowmx(kp)
           angmin=philowmi(kp)
        end if
        
        if ((angle .ge. angmin) .and. (angle .le. angmax)) then
           xfac=angmax-angle
           xfac=pi*xfac/(angmax-angmin)
           xfac=0.5D0*(1.0D0-cos(xfac))
        end if
     end if chord
  end if span

  decay=xfac*yfac*zfac
end subroutine flap_decay_funcs

subroutine deflect_flap(x,y,z,psinew,psiold,init)
  use ioutils
  use params_global
  use hqptef, only: flphx, flphz
  implicit none

  real, dimension(jmax,kmax,lmax) :: x, y, z
  real :: psinew, psiold
  logical :: init

  integer :: j, k, l, ihar, kp
  real :: thnew, thold, theta
  real :: xco, yco, zco, decay,cphi,sphi
  
  thnew=0.0D0
  thold=0.0D0
  do ihar=1,nharmflap
     thnew=thnew+ampflap(ihar)*sin(ihar*psinew+phiflap(ihar))
     thold=thold+ampflap(ihar)*sin(ihar*psiold+phiflap(ihar))
  end do
  
  if (init) then 
     theta=thnew
  else
     theta=thnew-thold
  end if

  write(STDOUT,*) "Trailing Edge Flap information"
  write(STDOUT,1001) psinew*180.0/pi, thnew*180.0/pi, theta*180.0/pi
  do l=1,lmax
     do k=1,kmax
        do j=1,jmax
           xco=x(j,k,l)
           yco=y(j,k,l)
           zco=z(j,k,l)
           call flap_decay_funcs(xco,yco,zco,decay,kp)
           if (kp .gt. 0) then
              cphi=cos(theta*decay)
              sphi=sin(theta*decay)
              x(j,k,l)=flphx(kp)+cphi*(xco-flphx(kp))+sphi*(zco-flphz(kp))
              z(j,k,l)=flphz(kp)-sphi*(xco-flphx(kp))+cphi*(zco-flphz(kp))
           endif
        end do
     end do
  end do

1001 format(3E15.6)
end subroutine deflect_flap
