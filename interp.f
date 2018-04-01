

c***************************************************************************c
      subroutine twodinterp(cl,cl1,x,x1,nr,nz,nr1,nz1)
      
      real cl(nz,nr),x(nr)
      real cl1(nz1,nr1),x1(nr1)

      real clt(nz1,nr)
      real defl(nz),rdef(nr)
      real defl1(nz1),rdef1(nr1)

c..   interpolate azimuthally

      do k=1,nr
         do i=1,nz
            defl(i)=cl(i,k)
         enddo

         call fourier_extend(defl,defl1,nz,nz1)
         
         do i=1,nz1
            clt(i,k)=defl1(i)
         enddo

      enddo
      
c..   now interpolate radially

      do i=1,nz1

         do k=1,nr
            rdef(k)=clt(i,k)
         enddo

c         call interp1(x,rdef,x1,rdef1,nr,nr1)
         
         call spline(x,rdef,x1,rdef1,nr,nr1)
         
         do k=1,nr1
            cl1(i,k)=rdef1(k)
         enddo
         
      enddo

      return
      end



      subroutine fourier_extend1(a,b,n,m)
      real a(*),b(*)
      complex, allocatable :: fhat(:),ghat(:)

      allocate(fhat(n),ghat(m))

      status=sfft('R','C','F',a,fhat,n,1)

      do i=1,m
         ghat(i)=0
      enddo

      if (m.ge.n) then

c..      fourier extension

         if (mod(n,2).eq.0) then
            nby2=n/2
            ghat(1:nby2)=fhat(1:nby2)
            ghat(nby2+1)=fhat(nby2+1)
         else
            nmin2=(n-1)/2
            ghat(1:nmin2+1)=fhat(1:nmin2+1)
         endif
      else

c..   fourier truncation

         if (mod(m,2).eq.0) then
            m1=m/2
            ghat(1:m1)=fhat(1:m1)
            ghat(m1+1)=(fhat(m1+1)+fhat(n-m1+1))*0.5
         else
            m1=(m-1)/2
            ghat(1:m1+1)=fhat(1:m1+1)
         endif
           
      endif

      ghat=ghat*m/n

c..   low pass filter, use only first 15 frequencies

      do i=17,m
         ghat(i)=0.
      enddo

      status=sfft('C','R','B',b,ghat,m,1)

      return
      end


      subroutine fourier_extend(a,b,n,m)
      real a(*),b(*)
      complex, allocatable :: fhat(:),ghat(:)

      allocate(fhat(n),ghat(m))

      status=sfft('R','C','F',a,fhat,n,1)

      do i=1,m
         ghat(i)=0
      enddo

      if (m.ge.n) then

c..      fourier extension

         if (mod(n,2).eq.0) then
            nby2=n/2
            ghat(1:nby2)=fhat(1:nby2)
            ghat(nby2+1)=fhat(nby2+1)
         else
            nmin2=(n-1)/2
            ghat(1:nmin2+1)=fhat(1:nmin2+1)
         endif
      else

c..   fourier truncation

         if (mod(m,2).eq.0) then
            m1=m/2
            ghat(1:m1)=fhat(1:m1)
            ghat(m1+1)=(fhat(m1+1)+fhat(n-m1+1))*0.5
         else
            m1=(m-1)/2
            ghat(1:m1+1)=fhat(1:m1+1)
         endif
           
      endif

      ghat=ghat*m/n
      
c      do i=17,m
c         ghat(i)=0.
c      enddo

      status=sfft('C','R','B',b,ghat,m,1)

      return
      end




            
c..   Wrapper function for VFFTPACK
c..   remember temp and fhat stick same positions for forward and backward
c..   transforms

      function sfft(a,b,c,temp,fhat,n,m)
      
      character a,b,c
      real temp(n)
      complex fhat(n)
      
      real, allocatable::r(:,:),rt(:,:),wsave(:)

      mdimx=1

      allocate(r(mdimx,n),rt(mdimx,n))
      allocate(wsave(n+15))
      fac=sqrt(1.0*n)
      fac1=1.0/fac
      
      if (c.eq.'F') then

         do i=1,n
            r(1,i)=temp(i)
         enddo

         call vrffti(n,wsave)
         call vrfftf(1,n,r,rt,mdimx,wsave)
         
         aa=r(1,1)*fac
         bb=0.0
         
         fhat(1)=cmplx(aa,bb)

         do k=1,n/2-1
            aa=r(1,2*k)*fac
            bb=r(1,2*k+1)*fac
            
            fhat(k+1)=cmplx(aa,bb)
         enddo

         if (mod(n,2).eq.0) then
	   aa=r(1,n)*fac
	   bb=0.0
	   fhat(n/2+1)=cmplx(aa,bb)
	 endif

      elseif (c.eq.'B') then
         
         r(1,1)=real(fhat(1))*fac1
         
         do k=1,n/2-1
            r(1,2*k)=real(fhat(k+1))*fac1
            r(1,2*k+1)=imag(fhat(k+1))*fac1
         enddo
         
         if (mod(n,2).eq.0) then
	   r(1,n)=real(fhat(n/2+1))*fac1
         endif

         call vrffti(n,wsave)
         call vrfftb(1,n,r,rt,mdimx,wsave)
         
         do i=1,n
            temp(i)=r(1,i)
         enddo

      endif

      sfft=1
      return
      end


c------------------------------------------------------------------
c----------------------- linear interpolation ---------------------
      subroutine interp1(x,y,xk,yk,n,m)
      
      real x(*),y(*),xk(*),yk(*)
      integer found
            
      do i=1,m
c..   extrapolate if outside the limit

         if (xk(i).lt.x(1)) then
            yk(i)=(xk(i)-x(2))*(y(2)-y(1))/(x(2)-x(1))+y(2)
         elseif (xk(i).gt.x(n)) then
            yk(i)=(xk(i)-x(n))*(y(n)-y(n-1))/(x(n)-x(n-1))+y(n)
         else
c..   interpolate otherwise finding the closest points
            found=0
            j=2
            do while(found.eq.0.and.j.le.n)
               if ((x(j)-xk(i))*(x(j-1)-xk(i)).le.0) then
                  yk(i)=(y(j)*(xk(i)-x(j-1))+y(j-1)*(x(j)-xk(i)))/
     +                 (x(j)-x(j-1))
                  found=1
               endif
               j=j+1
            enddo
         endif
      enddo

      return
      end


      subroutine spline(x,y,xk,yk,n,m)

      implicit none

      integer n,m
      real x(n),y(n),xk(m),yk(m)
      real coef(2*n)

c..   local variables

      real sigma,xout,yout
      integer i

c..   extrapolate if outside the limit

      sigma=0.5
      call splico(x,y,coef,sigma,n)
      
      do i=1,m
         
         if (xk(i).lt.x(1)) then
            yk(i)=(xk(i)-x(2))*(y(2)-y(1))/(x(2)-x(1))+y(2)
         elseif (xk(i).gt.x(n)) then
            yk(i)=(xk(i)-x(n))*(y(n)-y(n-1))/(x(n)-x(n-1))+y(n)
         else
            xout=xk(i)
            call speval(x,y,coef,sigma,xout,yout,n)
            yk(i)=yout
         endif
         
      enddo
      
      return
      end
      
c...this subroutine computes the spline coefficients necessary to fit an
c   exponential spline to the input abscissa x and ordinate y 
c
c argument list:
c     x(n)       input   input data abscissa
c     y(n)       input   input data ordinate
c     coef(2*n)  output  spline coeficients
c     sigma      input   spline tension factor
c     n          input   input array size
c
      subroutine splico(x,y,coef,sigma,n)                                    
      dimension x(n),y(n),coef(n*2)                                            
c                                                                              
c--> end points                                                                
c                                                                              
      delx1= x(2)-x(1)                                                         
      dx1= (y(2)-y(1))/delx1                                                   
      delx2= x(3)-x(2)                                                         
      delx12= x(3)-x(1)                                                        
      c1= -(delx12+delx1)/delx12/delx1                                         
      c2= delx12/delx1/delx2                                                   
      c3= -delx1/delx12/delx2                                                  
      slpp1= c1*y(1)+c2*y(2)+c3*y(3)                                           
      deln= x(n)-x(n-1)                                                        
      delnm1= x(n-1)-x(n-2)                                                    
      delnn= x(n)-x(n-2)                                                       
      c1= (delnn+deln)/delnn/deln                                              
      c2= -delnn/deln/delnm1                                                   
      c3= deln/delnn/delnm1                                                    
      slppn= c3*y(n-2)+c2*y(n-1)+c1*y(n)                                       
      sigmap= sigma*(n-1)/(x(n)-x(1))                                          
      dels= sigmap*delx1                                                       
      exps= exp(dels)                                                          
      sinhs= (exps-1./exps)/2.                                                 
      sinhin= 1./(delx1*sinhs)                                                 
      diag1= sinhin*(dels*0.5*(exps+1./exps)-sinhs)                            
      diagin= 1./diag1                                                         
      coef(1)= diagin*(dx1-slpp1)                                              
      spdiag= sinhin*(sinhs-dels)                                              
      coef(n+1)= diagin*spdiag                                                 
c                                                                              
c--> interior points                                                           
c                                                                              
      do i=2,n-1                                                               
        delx2= x(i+1)-x(i)                                                     
        dx2= (y(i+1)-y(i))/delx2                                               
        dels= sigmap*delx2                                                     
        exps= exp(dels)                                                        
        sinhs= (exps-1./exps)/2.                                               
        sinhin= 1./(delx2*sinhs)                                               
        diag2= sinhin*(dels*0.5*(exps+1./exps)-sinhs)                          
        diagin= 1./(diag1+diag2-spdiag*coef(n+i-1))                            
        coef(i)= diagin*(dx2-dx1-spdiag*coef(i-1))                             
        spdiag= sinhin*(sinhs-dels)                                            
        coef(i+n)= diagin*spdiag                                               
        dx1= dx2                                                               
        diag1= diag2                                                           
      enddo                                                                    
      diagin= 1./(diag1-spdiag*coef(2*n-1))                                    
      coef(n)= diagin*(slppn-dx2-spdiag*coef(n-1))                             
      do i=n-1,1,-1                                                            
        coef(i)= coef(i)-coef(i+n)*coef(i+1)                                   
      enddo                                                                    
c                                                                              
      return                                                                   
      end                                                                      

                                                                              
c...this subroutine evaluates an exponential spline at the abscissa xout
c   and returns the value yout.
c
c argument list:
c     x(n)       input   input data abscissa
c     y(n)       input   input data ordinate
c     coef(2*n)  output  spline coeficients
c     sigma      input   spline tension factor
c     xout       input   abscissa at which spline evaluated
c     yout       input   spline value
c     n          input   input array size
c
      subroutine speval(x,y,coef,sigma,xout,yout,n)                          
      dimension x(n),y(n),coef(n*2)
chs                                                                            
      integer   index/0/,indm1/0/                                              
chs                                                                            
c                                                                              
c**spline evaluation                                                           
      s= x(n)-x(1)                                                             
      sigmap= sigma*(n-1)/s                                                    
      do i=1,n-1                                                               
        if(xout.ge.x(i) .and. xout.le.x(i+1)) then                             
          index= i+1                                                           
          indm1= i                                                             
          goto 100                                                             
        endif                                                                  
      enddo                                                                    
 100  del1= xout-x(indm1)                                                      
      del2= x(index)-xout                                                      
      dels= x(index)-x(indm1)                                                  
      exps1= exp(sigmap*del1)                                                  
      sinhd1= 0.5*(exps1-1./exps1)                                             
      exps= exp(sigmap*del2)                                                   
      sinhd2= 0.5*(exps-1./exps)                                               
      exps= exps1*exps                                                         
      sinhs= 0.5*(exps-1./exps)                                                
      yout= (coef(index)*sinhd1+coef(indm1)*sinhd2)/sinhs +                    
     &     ((y(index)-coef(index))*del1+(y(indm1)-coef(indm1))*del2)           
     &      /dels                                                              
c                                                                              
      return                                                                   
      end                                                                      

