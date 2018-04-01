 
c---- interpolate values from one mesh to the other based on connectivity

      subroutine  do_interpolations(q,vnu,jmax,kmax,lmax,
     <     q1,vnu1,jmax1,kmax1,lmax1,
     <     imesh,idonor,frac,ndonor,Mpts)

      real q(jmax,kmax,lmax,6)
      real q1(jmax1,kmax1,lmax1,6)
      real vnu(jmax,kmax,lmax)
      real vnu1(jmax1,kmax1,lmax1)
      integer imesh(Mpts,3),idonor(Mpts,3)
      real frac(Mpts,3),frac1(Mpts,3)

      do id=1,ndonor
         i=imesh(id,1)
         j=imesh(id,2)
         k=imesh(id,3)

         ii=idonor(id,1)
         jj=idonor(id,2)
         kk=idonor(id,3)

         iip=min(ii+1,jmax1)
         jjp=min(jj+1,kmax1)
         kkp=min(kk+1,lmax1)

         dj0=1.-frac(id,1)
         dk0=1.-frac(id,2)
         dl0=1.-frac(id,3)
         dj1=frac(id,1)
         dk1=frac(id,2)
         dl1=frac(id,3)
         
         w1   = (dj0*dk0)*dl0
         w2   = (dj1*dk0)*dl0
         w3   = (dj0*dk1)*dl0
         w4   = (dj1*dk1)*dl0
         w5   = (dj0*dk0)*dl1
         w6   = (dj1*dk0)*dl1
         w7   = (dj0*dk1)*dl1
         w8   = (dj1*dk1)*dl1
         
         do n=1,5
            q(i,j,k,n)= w1*q1(ii ,jj ,kk ,n)*q1(ii ,jj ,kk ,6) 
     &                + w2*q1(iip,jj ,kk ,n)*q1(iip,jj ,kk ,6)
     &                + w3*q1(ii ,jjp,kk ,n)*q1(ii ,jjp,kk ,6) 
     &                + w4*q1(iip,jjp,kk ,n)*q1(iip,jjp,kk ,6)
     &                + w5*q1(ii ,jj ,kkp,n)*q1(ii ,jj ,kkp,6) 
     &                + w6*q1(iip,jj ,kkp,n)*q1(iip,jj ,kkp,6)
     &                + w7*q1(ii ,jjp,kkp,n)*q1(ii ,jjp,kkp,6) 
     &                + w8*q1(iip,jjp,kkp,n)*q1(iip,jjp,kkp,6)
            q(i,j,k,n)=q(i,j,k,n)/q(i,j,k,6)
         enddo
         
         vnu(i,j,k)= w1*vnu1(ii ,jj ,kk )!/q1(ii ,jj ,kk ,6) 
     &        	   + w2*vnu1(iip,jj ,kk )!/q1(iip,jj ,kk ,6)
     &             + w3*vnu1(ii ,jjp,kk )!/q1(ii ,jjp,kk ,6) 
     &		   + w4*vnu1(iip,jjp,kk )!/q1(iip,jjp,kk ,6)
     &             + w5*vnu1(ii ,jj ,kkp)!/q1(ii ,jj ,kkp,6) 
     &		   + w6*vnu1(iip,jj ,kkp)!/q1(iip,jj ,kkp,6)
     &             + w7*vnu1(ii ,jjp,kkp)!/q1(ii ,jjp,kkp,6) 
     &		   + w8*vnu1(iip,jjp,kkp)!/q1(iip,jjp,kkp,6)


      enddo
      
      return
      end
      

c---- set iblanks to one --------------------------------------------c
      subroutine set_iblanks_to_one(iblank,jmax,kmax,lmax)

      integer iblank(jmax,kmax,lmax)

      do j=1,jmax
         do k=1,kmax
            do l=1,lmax
               iblank(j,k,l)=1.
            enddo
         enddo
      enddo
      
      return
      end


c---- build_octree for background mesh ----------------------------c      

      subroutine make_octree(x1,y1,z1,jmax1,kmax1,lmax1)

      use work
      use bg_octree

      real x1(jmax1,kmax1,lmax1),y1(jmax1,kmax1,lmax1),z1(jmax1,kmax1,lmax1)
      

      jmax2=jmax1-1
      kmax2=kmax1-1
      lmax2=lmax1-1

      max_sor=jmax2*kmax2*lmax2
      allocate(xs1(max_sor),ys1(max_sor),zs1(max_sor))
      allocate(pindex(max_sor))
c
c..   allocation block for work arrays in split_box
c
c      write(6,*) 'Allocating arrays for split box'
      allocate(xpp(max_sor),ypp(max_sor),zpp(max_sor))
      allocate(xppp(max_sor),yppp(max_sor),zppp(max_sor))
      allocate(kpb(max_sor),kpc(max_sor))
c      write(6,*) 'finished allocating arrays for split box'

      is=0
      do l=1,lmax2
          do k=1,kmax2
             do j=1,jmax2

                   is=is+1
                   xs1(is)=0
                   ys1(is)=0
                   zs1(is)=0
                   
                   do il=0,1
                      do ik=0,1
                         do ij=0,1
                            xs1(is)=xs1(is)+x1(j+ij,k+ik,l+il)
                            ys1(is)=ys1(is)+y1(j+ij,k+ik,l+il)
                            zs1(is)=zs1(is)+z1(j+ij,k+ik,l+il)
                         enddo
                      enddo
                   enddo
                   
                   xs1(is)=xs1(is)*0.125
                   ys1(is)=ys1(is)*0.125
                   zs1(is)=zs1(is)*0.125

                   pindex(is)=(l-1)*kmax2*jmax2+(k-1)*jmax2+j


             enddo
          enddo
       enddo


      nsor=is
c      print *,'nsor',nsor
      max_pts=10
      mbox=nsor/max_pts
      lvl=log(1.*mbox)/log(8.)
c      write(6,*) 'lvl=',lvl
      allocate(xcls(0:mbox),ycls(0:mbox),zcls(0:mbox))
      allocate(lks(max_sor,lvl),ncls(mbox),inxcbs(8,0:mbox),kns(mbox))
      allocate(nboxs(0:lvl),lboxs(lvl),nsups(mbox),sboxs(0:lvl))

      call build_octree(nsor,xs1,ys1,zs1,lvl,
     &     xcls,ycls,zcls,
     &     sboxs,
     &     lks,ncls,nboxs,lboxs,nsups,inxcbs,kns,
     &     nsor,lvl,mbox,max_pts)

c      write(6,*) 'octree build complete'

      deallocate(xpp,ypp,zpp)
      deallocate(xppp,yppp,zppp)
      deallocate(kpb,kpc)

      return
      end


      subroutine do_connect(x,y,z,jmax,kmax,lmax,
     &     x1,y1,z1,iblank1,jmax1,kmax1,lmax1,
     &     imesh,idonor,frac,
     &     imesh1,idonor1,frac1,ndonor,nfri,
     &     Mpts,psi_rot,root_co)
      
      use bg_octree

      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)
      real x1(jmax1,kmax1,lmax1),y1(jmax1,kmax1,lmax1),z1(jmax1,kmax1,lmax1)
      integer iblank1(jmax1,kmax1,lmax1)
      integer imesh(Mpts,3),idonor(Mpts,3)
      integer imesh1(Mpts,3),idonor1(Mpts,3)
      real frac(Mpts,3),frac1(Mpts,3)
      real psi_rot
      integer root_co

      integer, allocatable:: iblank(:,:,:)

      real, allocatable:: xs11(:),ys11(:),zs11(:)
      real, allocatable:: xs(:),ys(:),zs(:)
      integer, allocatable:: ioct(:),ioct1(:)

      real aa(3,3),eigenv(3)


c***  first executable statement
      
      allocate(iblank(jmax,kmax,lmax))

      cs1=cos(psi_rot)
      ss1=sin(psi_rot)
       
      !rotate mesh back to pzi=0

       do l=1,lmax
          do k=1,kmax
             do j=1,jmax
              
                xxx=x(j,k,l)*cs1+y(j,k,l)*ss1
                yyy=y(j,k,l)*cs1-x(j,k,l)*ss1
              
                x(j,k,l)=xxx
                y(j,k,l)=yyy
             enddo
          enddo
       enddo

       do l=1,lmax1
          do k=1,kmax1
             do j=1,jmax1
              
                xxx=x1(j,k,l)*cs1+y1(j,k,l)*ss1
                yyy=y1(j,k,l)*cs1-x1(j,k,l)*ss1
              
                x1(j,k,l)=xxx
                y1(j,k,l)=yyy
             enddo
          enddo
       enddo
      

      iblank=1

      ! find boundary pts

      jtemp=jmax
      ktemp=kmax
      ltemp=lmax

      if (root_co.eq.0) then
        Nbpts=jtemp*ktemp+jtemp*(ltemp-1)+2*ktemp*(ltemp-1)
      else
        Nbpts=jtemp*ktemp+2*ktemp*(ltemp-1)
      endif

!      allocate(xs(Nbpts),ys(Nbpts),zs(Nbpts),ioct(Nbpts))
      allocate(xs(Mpts),ys(Mpts),zs(Mpts),ioct(Mpts))

c      print *,'finding boundary points'

      call find_bndrypts(xs,ys,zs,ioct,x,y,z,Mpts,Nbpts,imesh,
     &     jtemp,ktemp,ltemp,jmax,kmax,lmax,root_co,0)


c      print *,'finding outside points'

      call find_outsidevol(xs,ys,zs,ioct,Mpts,Nbpts,
     &     x1,y1,z1,imesh,idonor,frac,iblank1,jmax1,kmax1,lmax1,1,
     &     max_sor,mbox,lvl,xs1,ys1,zs1,pindex,
     &     xcls,ycls,zcls,
     &     sboxs,
     &     lks,ncls,nboxs,lboxs,nsups,inxcbs,kns)

      ndonor=Nbpts

      ! define a smaller boundary for hole cut

      jtemp=jmax-30
      ktemp=kmax
      ltemp=lmax-20

      Nbpts1=jtemp*ktemp
      Nbpts2=Nbpts1

!      allocate(xs11(Nbpts1),ys11(Nbpts1),zs11(Nbpts1),ioct1(Nbpts1))
 
      allocate(xs11(Mpts),ys11(Mpts),zs11(Mpts),ioct1(Mpts))

c      print *, 'Boundary points 2'
      call find_bndrypts(xs11,ys11,zs11,ioct1,x,y,z,Mpts,Nbpts1,imesh1,
     &     jtemp,ktemp,ltemp,jmax,kmax,lmax,root_co,1)


c      print *, 'Hole cut'
      call hole_cut(xs11,ys11,zs11,ioct1,imesh1,Mpts,Nbpts1,
     &     x,y,z,jmax,kmax,lmax,
     &     x1,y1,z1,iblank1,jmax1,kmax1,lmax1)

c      print *, 'Finished hole cut'
      
      !jtemp=jtemp-2
      !ktemp=ktemp-2
      !ltemp=ltemp-2

      do j=(jmax-jtemp)/2+1,jmax-(jmax-jtemp)/2
         do k=kmax-ktemp+1,kmax
            do l=1,ltemp-1
               iblank(j,k,l)=0
            enddo
         enddo
      enddo

c      print *, 'Inside volume'
      call find_insidevol(xs11,ys11,zs11,ioct1,Mpts,Nbpts1,
     &     x,y,z,imesh1,idonor1,frac1,iblank,jmax,kmax,lmax,0)

c      print *, 'Processed ivol'
      nfri=Nbpts1

      deallocate(iblank)
      deallocate(xs11,ys11,zs11)
      deallocate(xs,ys,zs)
      deallocate(ioct,ioct1)

      ! rotate mesh back to original position


      do l=1,lmax
         do k=1,kmax
             do j=1,jmax
              
                xxx=x(j,k,l)*cs1-y(j,k,l)*ss1
                yyy=y(j,k,l)*cs1+x(j,k,l)*ss1
              
                x(j,k,l)=xxx
                y(j,k,l)=yyy
             enddo
          enddo
       enddo

       do l=1,lmax1
          do k=1,kmax1
             do j=1,jmax1
              
                xxx=x1(j,k,l)*cs1-y1(j,k,l)*ss1
                yyy=y1(j,k,l)*cs1+x1(j,k,l)*ss1
              
                x1(j,k,l)=xxx
                y1(j,k,l)=yyy
             enddo
          enddo
       enddo

       
       if (1.eq.0) then
      write(6,*) 'Writing files..'

      open(unit=21,file='mcon',status='unknown')
      do i=1,Nbpts
         write(21,1001) imesh(i,1),imesh(i,2),imesh(i,3),idonor(i,1),
     &        idonor(i,2),idonor(i,3),frac(i,1),frac(i,2),frac(i,3)
      enddo
      close(21)

      open(unit=21,file='mcon1',status='unknown')
      do i=1,Nbpts1
         write(21,1001) imesh1(i,1),imesh1(i,2),imesh1(i,3),idonor1(i,1),
     &        idonor1(i,2),idonor1(i,3),frac1(i,1),frac1(i,2),frac1(i,3)
      enddo
      close(21)

      write(2) Nbpts
      write(2) (imesh(i,1),i=1,Nbpts),
     <     (imesh(i,2),i=1,Nbpts),
     <     (imesh(i,3),i=1,Nbpts),
     <     (2,i=1,Nbpts),
     <     (idonor(i,1),i=1,Nbpts),
     <     (idonor(i,2),i=1,Nbpts),
     <     (idonor(i,3),i=1,Nbpts),
     <     (frac(i,1),i=1,Nbpts),
     <     (frac(i,2),i=1,Nbpts),
     <     (frac(i,3),i=1,Nbpts)

      write(2) Nbpts1
      write(2) (imesh1(i,1),i=1,Nbpts1),
     <     (imesh1(i,2),i=1,Nbpts1),
     <     (imesh1(i,3),i=1,Nbpts1),
     <     (1,i=1,Nbpts1),
     <     (idonor1(i,1),i=1,Nbpts1),
     <     (idonor1(i,2),i=1,Nbpts1),
     <     (idonor1(i,3),i=1,Nbpts1),
     <     (frac1(i,1),i=1,Nbpts1),
     <     (frac1(i,2),i=1,Nbpts1),
     <     (frac1(i,3),i=1,Nbpts1)

      endif

 1001 format(6(I5,X),X,3(F14.8,X))       


      end


c---  subroutine hole_cut

      subroutine hole_cut(xs,ys,zs,ioct,imesh,Mpts,Nbpts,
     &     x,y,z,jmax,kmax,lmax,
     &     x1,y1,z1,iblank,jmax1,kmax1,lmax1)

      real xs(Mpts),ys(Mpts),zs(Mpts)
      integer ioct(Mpts),imesh(Mpts,3)
      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)
      real x1(jmax1,kmax1,lmax1),y1(jmax1,kmax1,lmax1),z1(jmax1,kmax1,lmax1)
      integer iblank(jmax1,kmax1,lmax1)
      real aa(3,3),eigenv(3)

      ! find the principal axes first by  finding the least-square axes

      xbar=0.
      ybar=0.
      zbar=0.
      
      do inb=1,Nbpts
         xbar=xbar+xs(inb)
         ybar=ybar+ys(inb)
         zbar=zbar+zs(inb)
      enddo

      xbar=xbar/Nbpts
      ybar=ybar/Nbpts
      zbar=zbar/Nbpts
      
      aa=0.

      do inb=1,Nbpts
         aa(1,1)=aa(1,1)+(xs(inb)-xbar)*(xs(inb)-xbar)
         aa(1,2)=aa(1,2)+(xs(inb)-xbar)*(ys(inb)-ybar)
         aa(1,3)=aa(1,3)+(xs(inb)-xbar)*(zs(inb)-zbar)
         
         aa(2,1)=aa(1,2)+(ys(inb)-ybar)*(xs(inb)-xbar)
         aa(2,2)=aa(2,2)+(ys(inb)-ybar)*(ys(inb)-ybar)
         aa(2,3)=aa(2,3)+(ys(inb)-ybar)*(zs(inb)-zbar)
         
         aa(3,1)=aa(3,1)+(zs(inb)-zbar)*(xs(inb)-xbar)
         aa(3,2)=aa(3,2)+(zs(inb)-zbar)*(ys(inb)-ybar)
         aa(3,3)=aa(3,3)+(zs(inb)-zbar)*(zs(inb)-zbar)
      enddo

      ! kaisers method for finding eigen values of real,symmetric matrices

      call kaiser(aa,3,3,eigenv,trace,sume,ier)

      ! remember eigen vectors are principal axes (refer inertial bisection)
      ! now find the largest distance in the principal directions
      ! this defines the box that is used to hole cut

      xdistx=0
      ydistx=0
      zdistx=0

      xdisti=1e6
      ydisti=1e6
      zdisti=1e6

      do inb=1,Nbpts

         ! dot products with principal axes

         xd=(xs(inb)-xbar)*aa(1,1)+(ys(inb)-ybar)*aa(2,1)
     &        +(zs(inb)-zbar)*aa(3,1)
         yd=(xs(inb)-xbar)*aa(1,2)+(ys(inb)-ybar)*aa(2,2)
     &        +(zs(inb)-zbar)*aa(3,2)
         zd=(xs(inb)-xbar)*aa(1,3)+(ys(inb)-ybar)*aa(2,3)
     &        +(zs(inb)-zbar)*aa(3,3)

         if (xd.gt.xdistx) xdistx=xd
         if (yd.gt.ydistx) ydistx=yd
         if (zd.gt.zdistx) zdistx=zd

         if (xd.lt.xdisti) xdisti=xd
         if (yd.lt.ydisti) ydisti=yd
         if (zd.lt.zdisti) zdisti=zd

      enddo

      xdist=(xdistx-xdisti)*0.5
      ydist=(ydistx-ydisti)*0.5
      zdist=(zdistx-zdisti)*0.5

      xbar1=(xdistx+xdisti)*0.5
      ybar1=(ydistx+ydisti)*0.5
      zbar1=(zdistx+zdisti)*0.5

      xbar=xbar1*aa(1,1)+ybar1*aa(1,2)+zbar1*aa(1,3)+xbar
      ybar=xbar1*aa(2,1)+ybar1*aa(2,2)+zbar1*aa(2,3)+ybar
      zbar=xbar1*aa(3,1)+ybar1*aa(3,2)+zbar1*aa(3,3)+zbar

c      ydist=1.5*ydist




      ! now let us iblank
      jmi=jmax1
      kmi=kmax1
      lmi=lmax1

      jmx=1
      kmx=1
      lmx=1

      Ncpts=0

      do l=1,lmax1
         do k=1,kmax1
            do j=1,jmax1

!               iblank(j,k,l)=1
               xd=(x1(j,k,l)-xbar)*aa(1,1)+(y1(j,k,l)-ybar)*aa(2,1)
     &              +(z1(j,k,l)-zbar)*aa(3,1)
               yd=(x1(j,k,l)-xbar)*aa(1,2)+(y1(j,k,l)-ybar)*aa(2,2)
     &              +(z1(j,k,l)-zbar)*aa(3,2)
               zd=(x1(j,k,l)-xbar)*aa(1,3)+(y1(j,k,l)-ybar)*aa(2,3)
     &              +(z1(j,k,l)-zbar)*aa(3,3)

               if (abs(xd).lt.xdist.and.abs(yd).lt.ydist.
     &              and.abs(zd).lt.zdist) then
                  
                  iblank(j,k,l)=0
              
                  if (j.lt.jmi) jmi=j
                  if (j.gt.jmx) jmx=j

                  if (k.lt.kmi) kmi=k
                  if (k.gt.kmx) kmx=k

                  if (l.lt.lmi) lmi=l
                  if (l.gt.lmx) lmx=l
                  
                  Ncpts=Ncpts+1

               endif
           
      
            enddo
         enddo
      enddo


      ! find hole fringe points
      
      write(6,*) jmi,jmx,kmi,kmx,lmi,lmx

      nfringe=1

      Ncpts=0

      jmi=jmi-nfringe
      jmx=jmx+nfringe

      if (jmi.lt.1) jmi=1
      if (jmx.gt.jmax1) jmx=jmax1

      kmi=kmi-nfringe
      kmx=kmx+nfringe
      
      if (kmi.lt.1) kmi=1
      if (kmx.gt.kmax1) kmx=kmax1

      lmi=lmi-nfringe
      lmx=lmx+nfringe

      if (lmi.lt.1) lmi=1
      if (lmx.gt.lmax1) lmx=lmax1

      ! j-direction search

      do l=lmi,lmx
         do k=kmi,kmx

            ifound=0
            icount=0

            do j=jmi,jmx
               
               
               if (iblank(j,k,l).eq.0 ) then

                  ifound=1
                
                  if (icount.eq.0) then
                     do ii=1,nfringe
                        jm1=j-ii
                        if (jm1.lt.1) jm1=1
                        iblank(jm1,k,l)=-1
                        Ncpts=Ncpts+1
                     enddo
                     icount=1
                  endif

               else if (ifound.eq.1) then
                  
                     do ii=0,nfringe-1
                        jp1=j+ii
                        if (jp1.gt.jmax1) jp1=jmax1
                        iblank(jp1,k,l)=-1
                        Ncpts=Ncpts+1
                     enddo
                                    
                  icount=0
                  ifound=0
               endif

            enddo
         enddo
      enddo
      

      ! k-direction search

      do l=lmi,lmx
         do j=jmi,jmx

            ifound=0
            icount=0

            do k=kmi,kmx

               
               if (iblank(j,k,l).eq.0) then

                  ifound=1

                  if (icount.eq.0) then
                     do ii=1,nfringe
                        km1=k-ii
                        if (km1.lt.1) km1=1
                        iblank(j,km1,l)=-1
                        Ncpts=Ncpts+1
                     enddo
                     icount=1
                  endif

               else if (ifound.eq.1) then
                  
                     do ii=0,nfringe-1
                        kp1=k+ii
                        if (kp1.gt.kmax1) kp1=kmax1
                        iblank(j,kp1,l)=-1
                        Ncpts=Ncpts+1
                     enddo
                                    
                  icount=0
                  ifound=0
               endif

            enddo
         enddo
      enddo

      ! l-direction search

      do j=jmi,jmx
         do k=kmi,kmx

            ifound=0
            icount=0

            do l=lmi,lmx

               
               if ( iblank(j,k,l).eq.0) then

                  ifound=1
                  
                  if (icount.eq.0) then
                     do ii=1,nfringe
                        lm1=l-ii
                        if (lm1.lt.1) lm1=1
                        iblank(j,k,lm1)=-1
                        Ncpts=Ncpts+1
                     enddo
                     icount=1
                  endif

               else if (ifound.eq.1) then
                  
                     do ii=0,nfringe-1
                        lp1=l+ii
                        if (lp1.gt.lmax1) lp1=lmax1
                        iblank(j,k,lp1)=-1
                        Ncpts=Ncpts+1
                     enddo
                                    
                  icount=0
                  ifound=0
               endif

            enddo
         enddo
      enddo

      Ncpts=0
      ii=1
      ioct(ii)=1
      
      do l=lmi,lmx
         do j=jmi,jmx
            
            iflag=0
            do k=kmi,kmx
               if (iblank(j,k,l).eq.-1) then
                  Ncpts=Ncpts+1
                  xs(Ncpts)=x1(j,k,l)
                  ys(Ncpts)=y1(j,k,l)
                  zs(Ncpts)=z1(j,k,l)
                  
                  imesh(Ncpts,1)=j
                  imesh(Ncpts,2)=k
                  imesh(Ncpts,3)=l
!                  iblank(j,k,l)=1
                  iblank(j,k,l)=0
                  iflag=1
!               else
!                  iblank(j,k,l)=1
               endif
            enddo

            if (iflag.eq.1) then
               ii=ii+1
               ioct(ii)=Ncpts
            endif

         enddo
      enddo

      do i=ii+1,Ncpts
         ioct(ii)=0
      enddo

      Nbpts=Ncpts

      return
      end
      



c---- subroutine print_volume

      subroutine print_volume(x1,y1,z1,j,k,l,jmax,kmax,lmax)
      
      real x1(jmax,kmax,lmax),y1(jmax,kmax,lmax),z1(jmax,kmax,lmax)
      real xv(8),yv(8),zv(8)


      xv(1)=x1(j,k,l);yv(1)=y1(j,k,l);zv(1)=z1(j,k,l)
      xv(2)=x1(j+1,k,l);yv(2)=y1(j+1,k,l);zv(2)=z1(j+1,k,l)
      xv(3)=x1(j+1,k+1,l);yv(3)=y1(j+1,k+1,l);zv(3)=z1(j+1,k+1,l)
      xv(4)=x1(j,k+1,l);yv(4)=y1(j,k+1,l);zv(4)=z1(j,k+1,l)
      
      xv(5)=x1(j,k,l+1);yv(5)=y1(j,k,l+1);zv(5)=z1(j,k,l+1)
      xv(6)=x1(j+1,k,l+1);yv(6)=y1(j+1,k,l+1);zv(6)=z1(j+1,k,l+1)
      xv(7)=x1(j+1,k+1,l+1);yv(7)=y1(j+1,k+1,l+1);zv(7)=z1(j+1,k+1,l+1)
      xv(8)=x1(j,k+1,l+1);yv(8)=y1(j,k+1,l+1);zv(8)=z1(j,k+1,l+1)

      do i=1,8
         write(21,*) xv(i),yv(i),zv(i)
      enddo

      return
      end

c---- subroutine check_inside

      subroutine check_inside(xp,yp,zp,x,y,z,wface,chi,f)

      real x(8),y(8),z(8)
      real f,f1,f2,f3,f4,f5,f6
      integer wface(3),chi,wface1(3)
      real EPS

      EPS=0.0

      do i=1,3
c         wface1(i)=wface(i)
         wface1(i)=0
         wface(i)=0
      enddo

      if (wface1(3).ne.1) then
         f1= check_face(xp,yp,zp,1,2,3,4,x,y,z)
         if (f1.lt.EPS) then
c            write(6,*) 'f1=',f1
            wface(3)=-1
            chi=0
            f=f1
            return
         endif
      endif

      if (wface1(3).ne.-1) then
         f2= check_face(xp,yp,zp,5,8,7,6,x,y,z)
         if (f2.lt.EPS) then
c            write(6,*) 'f2=',f2
            wface(3)=1
            chi=0
            f=f2
            return
         endif
      endif

      if (wface1(2).ne.1) then
         f3= check_face(xp,yp,zp,1,5,6,2,x,y,z)
         if (f3.lt.EPS) then
            wface(2)=-1
            chi=0
            f=f3
            return
         endif
      endif

      if (wface1(2).ne.-1) then
         f4=check_face(xp,yp,zp,4,3,7,8,x,y,z)
         if (f4.lt.EPS) then
            wface(2)=1
            chi=0
            f=f4
            return
         endif
      endif

      
      if (wface1(1).ne.-1) then
         f5=check_face(xp,yp,zp,2,6,7,3,x,y,z)
         if (f5.lt.EPS) then
            wface(1)=1
            chi=0
            f=f5
            return
         endif
      endif

      if (wface1(1).ne.1) then
         f6=check_face(xp,yp,zp,1,4,8,5,x,y,z)
         if (f6.lt.EPS) then
            wface(1)=-1
            chi=0
            f=f6
            return
         endif
      endif

      chi=1

      return
      end


c--   subroutine check_face

      function check_face(xp,yp,zp,i,j,k,l,x,y,z)

      real x(8),y(8),z(8)
      real ax,ay,az,cx,cy,cz

      ax=x(k)-x(i)
      ay=y(k)-y(i)
      az=z(k)-z(i)

      cx=x(l)-x(j)
      cy=y(l)-y(j)
      cz=z(l)-z(j)

      px=ay*cz-az*cy
      py=az*cx-ax*cz
      pz=ax*cy-ay*cx

      pnorm=sqrt(px**2+py**2+pz**2)

      dist=(xp-x(i))*px+(yp-y(i))*py+(zp-z(i))*pz

      check_face=dist/pnorm

      return
      end


      function check_face1(xp,yp,zp,i,j,k,l,x,y,z)

      real x(8),y(8),z(8)
      real ax,ay,az,dx,dy,dz,cx,cy,cz
      real check_face

      ax=x(j)-x(i)
      ay=y(j)-y(i)
      az=z(j)-z(i)

      cx=x(k)-x(i)
      cy=y(k)-y(i)
      cz=z(k)-z(i)

      dx=x(l)-x(i)
      dy=y(l)-y(i)
      dz=z(l)-z(i)


      den=(ax*dy*cz-ax*cy*dz-ay*dx*cz+ay*cx*dz+az*dx*cy-az*cx*dy)
     
      if (abs(den).lt.1e-6) then
c         print *,ax,ay,az,dx,dy,dz,px,py,pz

         px=ay*dz-az*dy
         py=az*dx-ax*dz
         pz=ax*dy-ay*dx
         
         pnorm=sqrt(px*px+py*py+pz*pz)

         check_face1=((xp-x(i))*px+(yp-y(i))*py+(zp-z(i))*pz)
     &        /pnorm
      else

         den=1./den

         t11=(dy*cz-cy*dz+ay*dz-dy*az)*den
         t12=(dx*az-dx*cz+cx*dz-ax*dz)*den
         t13=(dx*cy-cx*dy+ax*dy-dx*ay)*den
             
         t21=(cy*az-ay*cz+ay*dz-dy*az)*den
         t22=(dx*az+ax*cz-cx*az-ax*dz)*den
         t23=(cx*ay+ax*dy-dx*ay-ax*cy)*den

         t31=(ay*dz-dy*az)
         t32=(dx*az-ax*dz)
         t33=(ax*dy-dx*ay)
         
         px=xp-x(i)
         py=yp-y(i)
         pz=zp-z(i)

         vx=t11*px + t12*py + t13*pz
         vy=t21*px + t22*py + t23*pz
         vz=t31*px + t32*py + t33*pz

         check_face1 = vz-vx*vy/den

      endif

      return
      end

c--   subroutine find_boundary points

      subroutine find_bndrypts(xs,ys,zs,ioct,x,y,z,Mpts,Nbpts,imesh,
     &     jtemp,ktemp,ltemp,jmax,kmax,lmax,root_co,itype)

      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)
      real xs(Mpts),ys(Mpts),zs(Mpts)
      integer root_co,itype
      integer imesh(Mpts,3)
      integer ioct(Mpts)
      
      jtail1=(jmax-jtemp)/2+1
      kbeg=kmax-ktemp+1
      ktemp=kmax

      i=0
      ii=0
      l=ltemp

      do j=jtail1,jmax-jtail1+1
         ii=ii+1
         ioct(ii)=i+1
         do k=kbeg,ktemp
            i=i+1
            xs(i)=x(j,k,l)
            ys(i)=y(j,k,l)
            zs(i)=z(j,k,l)
            imesh(i,1)=j
            imesh(i,2)=k
            imesh(i,3)=l
         enddo
      enddo

      if (itype.eq.1) then
         do i=ii+1,Nbpts
            ioct(ii)=0
         enddo
         return
      endif

      if (root_co.eq.0) then 
       k=kbeg
       do j=jtail1,jmax-jtail1+1
         ii=ii+1
         ioct(ii)=i+1
         do l=1,ltemp-1
            i=i+1
            xs(i)=x(j,k,l)
            ys(i)=y(j,k,l)
            zs(i)=z(j,k,l)
            imesh(i,1)=j
            imesh(i,2)=k
            imesh(i,3)=l
         enddo
      enddo
      endif

C       k=ktemp

C       do j=jtail1,jmax-jtail1+1
C          ii=ii+1
C          ioct(ii)=i+1
C          do l=1,ltemp-1
C             i=i+1
C             xs(i)=x(j,k,l)
C             ys(i)=y(j,k,l)
C             zs(i)=z(j,k,l)
C             imesh(i,1)=j
C             imesh(i,2)=k
C             imesh(i,3)=l
C          enddo
C       enddo

      j=jtail1
      do k=kbeg,ktemp
         ii=ii+1
         ioct(ii)=i+1
         do l=1,ltemp-1
            i=i+1
            xs(i)=x(j,k,l)
            ys(i)=y(j,k,l)
            zs(i)=z(j,k,l)
            imesh(i,1)=j
            imesh(i,2)=k
            imesh(i,3)=l
         enddo
      enddo

      j=jmax-jtail1+1
      do k=kbeg,ktemp
         ii=ii+1
         ioct(ii)=i+1
         do l=1,ltemp-1
            i=i+1
            xs(i)=x(j,k,l)
            ys(i)=y(j,k,l)
            zs(i)=z(j,k,l)
            imesh(i,1)=j
            imesh(i,2)=k
            imesh(i,3)=l
         enddo
      enddo

      do i=ii+1,Nbpts
         ioct(ii)=0
      enddo


      return
      end


c--   subroutine find outside volume

      subroutine find_insidevol(xs,ys,zs,ioct,Mpts,Nbpts,x1,y1,z1,
     &     imesh,idonor,frac,iblank,jmax1,kmax1,lmax1,itype)

      use work
      real xs(Mpts),ys(Mpts),zs(Mpts)
      integer ioct(Mpts)
      integer imesh(Mpts,3),idonor(Mpts,3)
      real x1(jmax1,kmax1,lmax1),y1(jmax1,kmax1,lmax1),z1(jmax1,kmax1,lmax1)
      integer iblank(jmax1,kmax1,lmax1)
      real xv(8),yv(8),zv(8)
      real frac(Mpts,3)
      integer wface(3),chi

      !local variables

      real, allocatable:: xs1(:),ys1(:),zs1(:)
      real, allocatable:: boxsize(:,:),boxcenter(:,:),sboxs(:)
      real, allocatable:: xcls(:),ycls(:),zcls(:)
      integer, allocatable::lks(:,:),ncls(:),inxcbs(:,:)
      integer, allocatable::kns(:),nboxs(:),lboxs(:),nsups(:)
      integer, allocatable::pindex(:)
      integer itimes
      real mindist
      
      jmax2=jmax1-1
      kmax2=kmax1-1
      lmax2=lmax1-1

      max_sor=jmax2*kmax2*lmax2
      allocate(xs1(max_sor),ys1(max_sor),zs1(max_sor),pindex(max_sor))
      !
      ! work arrays for split box
      !
      allocate(xpp(max_sor),ypp(max_sor),zpp(max_sor))
      allocate(xppp(max_sor),yppp(max_sor),zppp(max_sor))
      allocate(kpb(max_sor),kpc(max_sor))


      is=0

      do l=1,lmax2
          do k=1,kmax2
             do j=1,jmax2
              
                if (iblank(j,k,l).eq.1) then
                   is=is+1
                   xs1(is)=0
                   ys1(is)=0
                   zs1(is)=0
                   
                   do il=0,1
                      do ik=0,1
                         do ij=0,1
                            xs1(is)=xs1(is)+x1(j+ij,k+ik,l+il)
                            ys1(is)=ys1(is)+y1(j+ij,k+ik,l+il)
                            zs1(is)=zs1(is)+z1(j+ij,k+ik,l+il)
                         enddo
                      enddo
                   enddo
                   
                   xs1(is)=xs1(is)*0.125
                   ys1(is)=ys1(is)*0.125
                   zs1(is)=zs1(is)*0.125

                   pindex(is)=(l-1)*kmax2*jmax2+(k-1)*jmax2+j

                endif
             enddo
          enddo
       enddo

      nsor=is
c      print *,'forming octree for ', nsor,' points'
      max_pts=10
      mbox=nsor/max_pts
      lvl=log(1.*mbox)/log(8.)
      
      allocate(xcls(0:mbox),ycls(0:mbox),zcls(0:mbox))
      allocate(lks(nsor,lvl),ncls(mbox),inxcbs(8,0:mbox),kns(mbox))
      allocate(nboxs(0:lvl),lboxs(lvl),nsups(mbox),sboxs(0:lvl))

      call build_octree(nsor,xs1,ys1,zs1,lvl,
     &     xcls,ycls,zcls,
     &     sboxs,
     &     lks,ncls,nboxs,lboxs,nsups,inxcbs,kns,
     &     nsor,lvl,mbox,max_pts)
      
      !
      ! deallocate the work arrays for split box
      ! 
      deallocate(xpp,ypp,zpp)
      deallocate(xppp,yppp,zppp)
      deallocate(kpb,kpc)

c      write(6,*) nboxs
c      write(6,*) lboxs

c      print *,Nbpts
      ii=1
      i1=0
      iflag=0
      itimes=0

      do inb=1,Nbpts

         xx=xs(inb)
         yy=ys(inb)
         zz=zs(inb)

         jii=imesh(inb,1)
         kii=imesh(inb,2)
         lii=imesh(inb,3)
         
         iflag=0

 100     continue

         if (inb.eq.ioct(ii).or.iflag.eq.1) then

            if (iflag.eq.0) ii=ii+1

            l=0
            nbxs=0

            xi=(xx-xcls(nbxs)+sboxs(l))/sboxs(l)
            yi=(yy-ycls(nbxs)+sboxs(l))/sboxs(l)
            zi=(zz-zcls(nbxs)+sboxs(l))/sboxs(l)
            fi=xi*yi*zi
            
            if ( xi*(xi-2).ge.0.or.yi*(yi-2).ge.0.or.zi*(zi-2).ge.0) then
               print *,'out of domain'
               goto 200
            endif
            
            ix=abs(xi)
            iy=abs(yi)
            iz=abs(zi)
            
            ich=4*iz+2*iy+ix+1
            ichbx=inxcbs(ich,nbxs)
            
            do l=1,lvl-1
               
               nbxs=ichbx
               
               ix=(xx-xcls(nbxs)+sboxs(l))/sboxs(l)
               iy=(yy-ycls(nbxs)+sboxs(l))/sboxs(l)
               iz=(zz-zcls(nbxs)+sboxs(l))/sboxs(l)
               
               ich=4*iz+2*iy+ix+1
               if (inxcbs(ich,nbxs).eq.0) goto 30
               ichbx=nboxs(l)+inxcbs(ich,nbxs)
               
            enddo

 30         continue


            mindist=1e6
            mpt=0
            
            do ipts=1,ncls(ichbx)
               
               ip=lks(nsups(ichbx)+ipts,l)
               ip1=ip
               
               dist=(xx-xs1(ip1))**2+(yy-ys1(ip1))**2+(zz-zs1(ip1))**2
               if (dist.lt.mindist) then
                  mindist=dist
                  mpt=ip1
               endif
            enddo
            
            ip=pindex(mpt)
            ip=ip-1
            j=mod(ip,jmax2)+1
            ip=ip/jmax2
            k=mod(ip,kmax2)+1
            ip=ip/kmax2
            l=mod(ip,lmax2)+1

            jx=j
            kx=k
            lx=l

         endif  ! if (inb.eq.ioct(ii)) then
         
         wface(1)=0
         wface(2)=0
         wface(3)=0
         
         chi=0
         iter=0

         do while (chi.eq.0.and.iter.lt.100)
            
            j=j+wface(1)
            k=k+wface(2)
            l=l+wface(3)         

            call meshbc(j,k,l,jmax2,kmax2,lmax2,itype)

            xv(1)=x1(j,k,l);yv(1)=y1(j,k,l);zv(1)=z1(j,k,l)
            xv(2)=x1(j+1,k,l);yv(2)=y1(j+1,k,l);zv(2)=z1(j+1,k,l)
            xv(3)=x1(j+1,k+1,l);yv(3)=y1(j+1,k+1,l);zv(3)=z1(j+1,k+1,l)
            xv(4)=x1(j,k+1,l);yv(4)=y1(j,k+1,l);zv(4)=z1(j,k+1,l)
            
            xv(5)=x1(j,k,l+1);yv(5)=y1(j,k,l+1);zv(5)=z1(j,k,l+1)
            xv(6)=x1(j+1,k,l+1);yv(6)=y1(j+1,k,l+1);zv(6)=z1(j+1,k,l+1)
            xv(7)=x1(j+1,k+1,l+1);yv(7)=y1(j+1,k+1,l+1);zv(7)=z1(j+1,k+1,l+1)
            xv(8)=x1(j,k+1,l+1);yv(8)=y1(j,k+1,l+1);zv(8)=z1(j,k+1,l+1)
            
            call check_inside(xx,yy,zz,xv,yv,zv,wface,chi,f)
            
            iter=iter+1
         enddo

         if (chi.eq.0) then
            if (iflag.eq.0) then
               iflag=1
               itimes=itimes+1
               goto 100
            else
               i1=i1+1
               ! degenerate fractions to the closest point
               ! ideally it shouldn't come here, but this is fail safe

               write(23,'(I10,x,3(E10.4,x))') inb,xx,yy,zz

               imesh(i1,1)=jii
               imesh(i1,2)=kii
               imesh(i1,3)=lii
               
               idonor(i1,1)=jx
               idonor(i1,2)=kx
               idonor(i1,3)=lx

               frac(i1,1)=1.0
               frac(i1,2)=1.0
               frac(i1,3)=1.0
            endif
         else
            i1=i1+1
            imesh(i1,1)=jii
            imesh(i1,2)=kii
            imesh(i1,3)=lii

            idonor(i1,1)=j
            idonor(i1,2)=k
            idonor(i1,3)=l

            call find_fractions(f1,f2,f3,xx,yy,zz,xv,yv,zv)

            frac(i1,1)=f1
            frac(i1,2)=f2
            frac(i1,3)=f3
         endif

 200     continue
c        if (mod(inb,1000).eq.1) print *,inb
      enddo

c      write(6,*) 'resorted to octree search ',itimes,' times'

      Nbpts=i1

      return
      end


      subroutine find_fractions(f1,f2,f3,xp,yp,zp,x,y,z)

      real x(8),y(8),z(8)

      z1=check_face(xp,yp,zp,1,2,3,4,x,y,z)
      z2=check_face(xp,yp,zp,5,8,7,6,x,y,z)

      y1=check_face(xp,yp,zp,1,5,6,2,x,y,z)
      y2=check_face(xp,yp,zp,4,3,7,8,x,y,z)

      x1=check_face(xp,yp,zp,1,4,8,5,x,y,z)
      x2=check_face(xp,yp,zp,2,6,7,3,x,y,z)
 
      f1=x1/(x1+x2)
      f2=y1/(y1+y2)
      f3=z1/(z1+z2)

      return
      end



      subroutine meshbc(j,k,l,jmax2,kmax2,lmax2,itype)

      if (itype.eq.1) then
         
         if (j.gt.jmax2) j=1
         if (j.lt.1) j=jmax2
         
         if (k.gt.kmax2) k=kmax2
         if (k.lt.1) k=1
         
         if (l.gt.lmax2) l=lmax2
         if (l.lt.1) l=1

      else

         ! c-o mesh

         ! c-cut

         if (l.lt.1.and.k.gt.kmax2) then
            l=1
            k=kmax2
            j=jmax2-j+1
         else if (l.lt.1) then
            l=1
            j=jmax2-j+1
         else if (k.gt.kmax2) then
            k=kmax2
            j=jmax2-j+1
         endif

         if (k.lt.1) k=1
         if (j.lt.1) j=1
         if (j.gt.jmax2) j=jmax2
         if (l.gt.lmax2) l=lmax2

      endif

      return
      end



c---- subroutine check_inside

      subroutine check_inside1(xp,yp,zp,x,y,z,wface,chi,fx,imin)

      real x(8),y(8),z(8)
      real f(6)
      integer wface(3),chi,wface1(3)
      real EPS

      EPS=-1e-6

      do i=1,3
         wface(i)=0
      enddo

      f(1)=check_face1(xp,yp,zp,1,2,3,4,x,y,z)
      f(2)=check_face1(xp,yp,zp,5,8,7,6,x,y,z)
      f(3)=check_face1(xp,yp,zp,1,5,6,2,x,y,z)
      f(4)=check_face1(xp,yp,zp,4,3,7,8,x,y,z)
      f(5)=check_face1(xp,yp,zp,2,6,7,3,x,y,z)
      f(6)=check_face1(xp,yp,zp,1,4,8,5,x,y,z)

      fmin=f(1)
      imin=1

      do i=2,6
         if (f(i).lt.fmin) then
            fmin=f(i)
            imin=i
         endif
      enddo

      if (fmin.lt.EPS) then
         chi=0
         fx=fmin
         if (imin.eq.1) wface(3)=-1
         if (imin.eq.2) wface(3)=1
         if (imin.eq.3) wface(2)=-1
         if (imin.eq.4) wface(2)=1
         if (imin.eq.5) wface(1)=1
         if (imin.eq.6) wface(1)=-1
      else
         chi=1
      endif

      return
      end



c--   subroutine find outside volume


      subroutine find_outsidevol(xs,ys,zs,ioct,Mpts,Nbpts,x1,y1,z1,
     &     imesh,idonor,frac,iblank,jmax1,kmax1,lmax1,itype,
     &     max_sor,mbox,lvl,xs1,ys1,zs1,pindex,
     &     xcls,ycls,zcls,
     &     sboxs,
     &     lks,ncls,nboxs,lboxs,nsups,inxcbs,kns)

      real xs(Mpts),ys(Mpts),zs(Mpts)
      integer ioct(Mpts)
      integer imesh(Mpts,3),idonor(Mpts,3)
      real x1(jmax1,kmax1,lmax1),y1(jmax1,kmax1,lmax1),z1(jmax1,kmax1,lmax1)
      integer iblank(jmax1,kmax1,lmax1)
      real xv(8),yv(8),zv(8)
      real frac(Mpts,3)
      integer wface(3),chi

      !local variables

      real xs1(max_sor),ys1(max_sor),zs1(max_sor)
      real xcls(0:mbox),ycls(0:mbox),zcls(0:mbox)
      integer inxcbs(8,0:mbox),ncls(mbox),lks(max_sor,lvl)
      integer kns(mbox),pindex(max_sor)
      integer nboxs(0:lvl),lboxs(lvl),nsups(mbox)
      real sboxs(0:lvl)

      real mindist

c***  first executable statement

      jmax2=jmax1-1
      kmax2=kmax1-1
      lmax2=lmax1-1

c      print *,Nbpts

      ii=1
      i1=0
      iflag=0
      itimes=0

      do inb=1,Nbpts

         xx=xs(inb)
         yy=ys(inb)
         zz=zs(inb)

         jii=imesh(inb,1)
         kii=imesh(inb,2)
         lii=imesh(inb,3)
         
         iflag=0

 100     continue

         if (inb.eq.ioct(ii).or.iflag.eq.1) then

            if (iflag.eq.0) ii=ii+1

            l=0
            nbxs=0

            xi=(xx-xcls(nbxs)+sboxs(l))/sboxs(l)
            yi=(yy-ycls(nbxs)+sboxs(l))/sboxs(l)
            zi=(zz-zcls(nbxs)+sboxs(l))/sboxs(l)
            fi=xi*yi*zi
            
            if ( xi*(xi-2).ge.0.or.yi*(yi-2).ge.0.or.zi*(zi-2).ge.0) then
               print *,'out of domain'
               goto 200
            endif
            
            ix=abs(xi)
            iy=abs(yi)
            iz=abs(zi)
            
            ich=4*iz+2*iy+ix+1
            ichbx=inxcbs(ich,nbxs)
            
            do l=1,lvl-1
               
               nbxs=ichbx
               
               ix=(xx-xcls(nbxs)+sboxs(l))/sboxs(l)
               iy=(yy-ycls(nbxs)+sboxs(l))/sboxs(l)
               iz=(zz-zcls(nbxs)+sboxs(l))/sboxs(l)
               
               ich=4*iz+2*iy+ix+1
               if (inxcbs(ich,nbxs).eq.0) goto 30
               ichbx=nboxs(l)+inxcbs(ich,nbxs)
               
            enddo

 30         continue

            mindist=1e6
            mpt=0
            
            do ipts=1,ncls(ichbx)
               
               ip=lks(nsups(ichbx)+ipts,l)
               ip1=ip
               
               dist=(xx-xs1(ip1))**2+(yy-ys1(ip1))**2+(zz-zs1(ip1))**2
               if (dist.lt.mindist) then
                  mindist=dist
                  mpt=ip1
               endif
            enddo
            
            ip=pindex(mpt)
            ip=ip-1
            j=mod(ip,jmax2)+1
            ip=ip/jmax2
            k=mod(ip,kmax2)+1
            ip=ip/kmax2
            l=mod(ip,lmax2)+1

            jx=j
            kx=k
            lx=l

         endif  ! if (inb.eq.ioct(ii)) then
         
         wface(1)=0
         wface(2)=0
         wface(3)=0
         
         chi=0
         iter=0

         do while (chi.eq.0.and.iter.lt.100)
            
            j=j+wface(1)
            k=k+wface(2)
            l=l+wface(3)         

            call meshbc(j,k,l,jmax2,kmax2,lmax2,itype)

            xv(1)=x1(j,k,l);yv(1)=y1(j,k,l);zv(1)=z1(j,k,l)
            xv(2)=x1(j+1,k,l);yv(2)=y1(j+1,k,l);zv(2)=z1(j+1,k,l)
            xv(3)=x1(j+1,k+1,l);yv(3)=y1(j+1,k+1,l);zv(3)=z1(j+1,k+1,l)
            xv(4)=x1(j,k+1,l);yv(4)=y1(j,k+1,l);zv(4)=z1(j,k+1,l)
            
            xv(5)=x1(j,k,l+1);yv(5)=y1(j,k,l+1);zv(5)=z1(j,k,l+1)
            xv(6)=x1(j+1,k,l+1);yv(6)=y1(j+1,k,l+1);zv(6)=z1(j+1,k,l+1)
            xv(7)=x1(j+1,k+1,l+1);yv(7)=y1(j+1,k+1,l+1);zv(7)=z1(j+1,k+1,l+1)
            xv(8)=x1(j,k+1,l+1);yv(8)=y1(j,k+1,l+1);zv(8)=z1(j,k+1,l+1)
            
            call check_inside(xx,yy,zz,xv,yv,zv,wface,chi,f)
            
            iter=iter+1
         enddo

         if (chi.eq.0) then
            if (iflag.eq.0) then
               iflag=1
               itimes=itimes+1
               goto 100
            else
               i1=i1+1
               ! degenerate fractions to the closest point
               ! ideally it shouldn't come here, but this is fail safe

               write(23,'(I10,x,3(E10.4,x))') inb,xx,yy,zz

               imesh(i1,1)=jii
               imesh(i1,2)=kii
               imesh(i1,3)=lii
               
               idonor(i1,1)=jx
               idonor(i1,2)=kx
               idonor(i1,3)=lx

               frac(i1,1)=1.0
               frac(i1,2)=1.0
               frac(i1,3)=1.0
            endif
         else
            i1=i1+1
            imesh(i1,1)=jii
            imesh(i1,2)=kii
            imesh(i1,3)=lii

            idonor(i1,1)=j
            idonor(i1,2)=k
            idonor(i1,3)=l

            call find_fractions(f1,f2,f3,xx,yy,zz,xv,yv,zv)

            frac(i1,1)=f1
            frac(i1,2)=f2
            frac(i1,3)=f3
         endif

 200     continue
c        if (mod(inb,1000).eq.1) print *,inb
      enddo

      write(6,*) 'resorted  to octree search ',itimes,' times'
                     
      Nbpts=i1

      return
      end
