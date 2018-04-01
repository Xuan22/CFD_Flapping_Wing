cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     octree.f - generate the octree data structure
c
c     This subroutine generate cubic tree structure of the boxes and
c     sorts the points into the boxes at different levels of the tree
c     structure.
c
c     inputs:
c     [xyz]sor(MAX_TARGETS)	location of points (either field or source pts)
c     nsor		        number of points
c     lvl		         suggested number of levels of the tree 
c     
c     outputs:
c     [xyz]cls(MAX_BOXES)	        center of box, for all boxes
c     sboxs(MAX_LEVEL)	                diagonal side of a box of given level
c     [abc]boxs(MAX_BOXES)	        dimensions of box
c     lks(MAX_TARGETS,MAX_LEVEL)	global index of point (given level point index, level)
c     ncls(MAX_BOXES)	                number of points in given box
c     inxcbs(8,MAX_BOXES)	        box index of given child of given box
c     kns(MAX_BOXES)		        number of non-empty children of given box
c     nboxs(0:MAX_LEVEL)	        number of boxes in given level AND all rootward levels
c			Note that this starts at level 0, with 0 boxes,
c			then at level 1, it reads 1-8 boxes (depending);
c                       this is because the rootmost parent box is never empty
c     lboxs(1:MAX_LEVEL)	        number of boxes in given level
c     nsups(MAX_BOXES)	                total number of points in boxes 1 to given box - 1
c     
c     local only:
c     [xyz]c(MAX_TARGETS)	reordered point locations
c     kpa(MAX_TARGETS)	integer pointers to all points
c     [xyz]m[ix]l(MAX_BOXES)	min, max bounds of all boxes
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine build_octree(nsor,xsor,ysor,zsor,
     &   lvl,xcls,ycls,zcls,sboxs,
c     &   aboxs,bboxs,cboxs,
     &   lks,ncls,nboxs,lboxs,nsups,
     &   inxcbs,kns,MAX_TARGETS,MAX_LEVEL,MAX_BOXES,MAX_PTS)
      
      implicit none
      
      integer MAX_TARGETS,MAX_LEVEL,MAX_BOXES,MAX_PTS
      real xsor(MAX_TARGETS), ysor(MAX_TARGETS), zsor(MAX_TARGETS)
      integer nsor,lvl
      real xcls(0:MAX_BOXES), ycls(0:MAX_BOXES), zcls(0:MAX_BOXES)
      real sboxs(0:MAX_LEVEL)
c      real aboxs(MAX_BOXES),bboxs(MAX_BOXES),cboxs(MAX_BOXES)
      integer lks(MAX_TARGETS,MAX_LEVEL), ncls(MAX_BOXES),
     &        inxcbs(8,0:MAX_BOXES), kns(MAX_BOXES),
     &        nboxs(0:MAX_LEVEL),  lboxs(MAX_LEVEL), nsups(MAX_BOXES)


c..   local variables

      real, allocatable :: xmil(:),xmxl(:),ymil(:),ymxl(:)
      real, allocatable :: zmil(:),zmxl(:)
      real, allocatable :: xp(:),yp(:),zp(:)
      real, allocatable :: xs1(:),ys1(:),zs1(:)
      real, allocatable :: xs2(:),ys2(:),zs2(:)
      real, allocatable :: xc(:),yc(:),zc(:)
      integer, allocatable :: kc(:),kpa(:)

      integer nc(8),chid(8)
      integer ip,ipts,ichbx

c---- center, max and min of the child boxes

      real xcc(8),ycc(8),zcc(8)
      real xmaxc(8),ymaxc(8),zmaxc(8)
      real xminc(8),yminc(8),zminc(8)

      real xmi,xmx,ymi,ymx,zmi,zmx,xmm,ymm,zmm,rmm,drr
      real xsmi,xsmx,ysmi,ysmx,zsmi,zsmx

      integer i,j,k,l,kne,kp2,mii,j8,kp1,kl,np,nbxs
      real xminb,xmaxb,yminb,ymaxb,zminb,zmaxb

c=======================================================================
c---- find the bdry of the parent box, use 1st point as initial max/min

      allocate(xmil(MAX_BOXES), xmxl(MAX_BOXES), ymil(MAX_BOXES),
     &          ymxl(MAX_BOXES), zmil(MAX_BOXES), zmxl(MAX_BOXES))
      allocate(xp(MAX_TARGETS),yp(MAX_TARGETS),zp(MAX_TARGETS))
      allocate(xs1(MAX_TARGETS),ys1(MAX_TARGETS),zs1(MAX_TARGETS))
      allocate(xs2(MAX_TARGETS),ys2(MAX_TARGETS),zs2(MAX_TARGETS))
      allocate(xc(MAX_TARGETS),yc(MAX_TARGETS),zc(MAX_TARGETS))
      allocate( kc(MAX_TARGETS),kpa(MAX_TARGETS))

      xmi=xsor(1)
      xmx=xsor(1)
      ymi=ysor(1)
      ymx=ysor(1)
      zmi=zsor(1)
      zmx=zsor(1)
      do i=1,nsor
         if(xsor(i).lt.xmi) xmi=xsor(i)
         if(xsor(i).gt.xmx) xmx=xsor(i)
         if(ysor(i).lt.ymi) ymi=ysor(i)
         if(ysor(i).gt.ymx) ymx=ysor(i)
         if(zsor(i).lt.zmi) zmi=zsor(i)
         if(zsor(i).gt.zmx) zmx=zsor(i)
      enddo

c---- make sure it is a square cube.

      xmm=xmx-xmi
      ymm=ymx-ymi
      zmm=zmx-zmi
      rmm=xmm
      if(ymm.gt.rmm) rmm=ymm
      if(zmm.gt.rmm) rmm=zmm
      xmx=0.5*(xmi+xmx+rmm)
      xmi=xmx-rmm
      ymx=0.5*(ymi+ymx+rmm)
      ymi=ymx-rmm
      zmx=0.5*(zmi+zmx+rmm)
      zmi=zmx-rmm

c---- add epsilon to the outer boundaries, but do it so unevenly
c     so that fld pts on a pure axis are all in the same quadrant

      drr=amax1(rmm/1000., 1.0e-5)

      xsmi=xmi-1.5*drr
      xsmx=xmx+0.5*drr
      ysmi=ymi-1.5*drr
      ysmx=ymx+0.5*drr
      zsmi=zmi-1.5*drr
      zsmx=zmx+0.5*drr

      goto 20
 10   lvl=lvl-1
 20   continue

c---- set kpa pointers to initial point locations (in array)
      do k=1,nsor
         kpa(k)=k
      enddo

c---- set the diagonal side of a box of each level

      xcls(0)=xsmi+(xsmx-xsmi)*0.5
      ycls(0)=ysmi+(ysmx-ysmi)*0.5
      zcls(0)=zsmi+(zsmx-zsmi)*0.5

      do l=0,lvl
c         sboxs(l)=sqrt(3.)*(xsmx-xsmi)/2**l
         sboxs(l)=(xsmx-xsmi)/2**(l+1)
         ! write(*,*) 'sboxs(',l,') = ',sboxs(l)
      enddo

c---- set the number of full boxes at level 0 to 0

      nboxs(0)=0

c----------------------------------------------------------------------------
c---- run the routine through root level once
      l=1
c---- split the rootmost parent box into 8 child boxes
      call split_box(xsmi,xsmx,ysmi,ysmx,zsmi,zsmx,
     +     xsor,ysor,zsor,kpa,nsor,
     +     xc,yc,zc,kc,nc,kne,
     +     xminc,xmaxc,yminc,ymaxc,zminc,zmaxc,xcc,ycc,zcc,
     +     chid,MAX_TARGETS)

      kp2=0
      mii=0
c---- loop over the number of non-empty child boxes (0 up to 8)
      do j8=1,kne
c------- for each point in the box, store locations and pointer in a temp array
         do i=1,nc(j8)
            kp2=kp2+1
            mii=mii+1
            xs1(kp2)=xc(mii)
            ys1(kp2)=yc(mii)
            zs1(kp2)=zc(mii)
            lks(kp2,1)=kc(mii)
         enddo
         inxcbs(chid(j8),0)=j8
      enddo

c---- essentially, we just copied xc to xs1, kc to lks
c---- now, add the bounds, centers, sizes, and point counts of each
c        of the non-empty boxes to the global box arrays
      do j8=1,kne
         xcls(j8)=xcc(j8)
         ycls(j8)=ycc(j8)
         zcls(j8)=zcc(j8)
         ncls(j8)=nc(j8)

         xmil(j8)=xminc(j8)
         ymil(j8)=yminc(j8)
         zmil(j8)=zminc(j8)
         xmxl(j8)=xmaxc(j8)
         ymxl(j8)=ymaxc(j8)
         zmxl(j8)=zmaxc(j8)
c         aboxs(j8)=xmxl(j8)-xmil(j8)
c         bboxs(j8)=ymxl(j8)-ymil(j8)
c         cboxs(j8)=zmxl(j8)-zmil(j8)
      enddo

c---- and, since this is the first time we split, the number of boxes
c        in level 1 is 8 (or slightly less), and the total number of
c        boxes in the simulation is the same
      nboxs(1)=kne
      lboxs(1)=kne
         
c----------------------------------------------------------------------------
c     and for each level up the tree, do a similar thing
     
      do l=2,lvl
         kp1=0
         kp2=0
         kl=0

c------- for each of the non-empty boxes in the previous level,
c           treat it as a parent, and split it

         do j=1,lboxs(l-1)
c---------- so, for l=2, nboxs(0)=0, and np=ncls(j)
            np=ncls(nboxs(l-2)+j)
            if (np.gt.max_pts) then
               nbxs=nboxs(l-2)+j
c----------grab the dimensions of this box
               xminb=xmil(nbxs)
               yminb=ymil(nbxs)
               zminb=zmil(nbxs)
               xmaxb=xmxl(nbxs)
               ymaxb=ymxl(nbxs)
               zmaxb=zmxl(nbxs)
c----------and create the temporary array of point locations and pointer
c     for all of the points in this particular box (xs1 are the
c     reordered points, lks are reordered pointers)
               do i=1,np
                  kp1=kp1+1 
                  xp(i)=xs1(kp1)
                  yp(i)=ys1(kp1)
                  zp(i)=zs1(kp1)
                  kpa(i)=lks(kp1,l-1)
               enddo
               
               call split_box(xminb,xmaxb,yminb,ymaxb,zminb,zmaxb,
     +              xp,yp,zp,kpa,np,
     +              xc,yc,zc,kc,nc,kne,
     +              xminc,xmaxc,yminc,ymaxc,zminc,zmaxc,xcc,ycc,zcc,
     +              chid,MAX_TARGETS)
               
c----------just like before, copy all of the reordered locations and
c     pointers from the subroutine's arrays into the 
               mii=0
c----------loop over the total number of non-empty boxes returned by split_box
               do j8=1,kne
                  do i=1,nc(j8)
                     kp2=kp2+1
                     mii=mii+1
                     xs2(kp2)=xc(mii)
                     ys2(kp2)=yc(mii)
                     zs2(kp2)=zc(mii)
                     lks(kp2,l)=kc(mii)
                  enddo
               enddo
               
c----------do the same for the non-empty child box's dimensions and point count
               do j8=1,kne
c-------------kl is the counter for the number of boxes in level l
                  kl=kl+1
c-------------important: nbxs is the pointer to this new box
                  nbxs=nboxs(l-1)+kl
                  
                  if (nbxs.gt.MAX_BOXES) then
                   write(*,*) 'WARNING (octree.f): trying to allocate ',
     +                    'box',nbxs,' when MAX_BOXES is',MAX_BOXES
                     goto 10
                                ! stop
                  endif

                  xcls(nbxs)=xcc(j8)
                  ycls(nbxs)=ycc(j8)
                  zcls(nbxs)=zcc(j8)
                  ncls(nbxs)=nc(j8)
                  
                  xmil(nbxs)=xminc(j8)
                  ymil(nbxs)=yminc(j8)
                  zmil(nbxs)=zminc(j8)
                  xmxl(nbxs)=xmaxc(j8)
                  ymxl(nbxs)=ymaxc(j8)
                  zmxl(nbxs)=zmaxc(j8)
c                  aboxs(nbxs)=xmxl(nbxs)-xmil(nbxs)
c                  bboxs(nbxs)=ymxl(nbxs)-ymil(nbxs)
c                  cboxs(nbxs)=zmxl(nbxs)-zmil(nbxs)

c------------- important: set the child box pointers of the split box
                  inxcbs(chid(j8),nboxs(l-2)+j)=kl
               enddo
c----------set the array of non-empty children of the split box
               kns(nboxs(l-2)+j)=kne

            endif               ! end if loop (np.gt.max_pts) 
    
         enddo                  ! end loop over number of boxes in level l-1
         
c-------number of valid boxes in this level
         lboxs(l)=kl
c-------number of valid boxes in this and all rootward levels
         nboxs(l)=nboxs(l-1)+kl
         
c-------copy reordered array of locations
         do k=1,nsor
            xs1(k)=xs2(k)
            ys1(k)=ys2(k)
            zs1(k)=zs2(k)
         enddo
         
      enddo                     ! end loop l=2 to lvl
      
c-----------------------------------------------------------------------------

c---- nsups(i) is the pointer to the first point in box i, regardless of level
c     more correctly, nsups(i) is the number of points in boxes 1 to i-1

      do l=1,lvl
         nsups(nboxs(l-1)+1)=0

         do i=2,lboxs(l)
            nsups(nboxs(l-1)+i)=nsups(nboxs(l-1)+i-1)+
     +                          ncls(nboxs(l-1)+i-1)

         enddo
      enddo

      return
      end


c ============================================================================
c
c     This subroutine divides a box (parent) into 8 smaller boxes (children)
c     and sorts a set of points into the children boxes. 
c
c     inputs:
c     [xyz]m[in|ax]	outer bounds of the parent box
c     [xyz]p(MAX_SOURCES)	locations of ONLY the points in parent box
c     kpa(MAX_SOURCES)		integer pointers to the points in parent box
c     np		number of points in the parent box
c     
c     outputs:
c     [xyz]c(MAX_SOURCES)	*reordered* locations of ONLY the points in parent box
c     kc(MAX_SOURCES)		*reordered* integer pointers of the points in parent box
c     nc(8)		number of points in the child boxes
c     kne		number of non-empty children of this parent (up to 8)
c     [xyz]m[in|ax]c(8)	outer bounds of the child boxes
c     [xyz]cc(8)	center of the child boxes
c     

      subroutine split_box(xmin,xmax,ymin,ymax,zmin,zmax,
     +      xp,yp,zp,kpa,np,
     +      xc,yc,zc,kc,nc,kne,
     +      xminc,xmaxc,yminc,ymaxc,zminc,zmaxc,xcc,ycc,zcc,chid,
     +     MAX_SOURCES)

      use work
      implicit none

      integer MAX_SOURCES
      real xmin,xmax,ymin,ymax,zmin,zmax
      real xp(MAX_SOURCES),yp(MAX_SOURCES),zp(MAX_SOURCES)
      integer kpa(MAX_SOURCES),np
      real xc(MAX_SOURCES),yc(MAX_SOURCES),zc(MAX_SOURCES)
      integer kc(MAX_SOURCES),nc(8),kne,chid(8)
      real xminc(8),yminc(8),zminc(8)
      real xmaxc(8),ymaxc(8),zmaxc(8)
      real xcc(8),ycc(8),zcc(8)

c---- working arrays


      integer nco(8)
      integer i,j,k,ii,npp,mii,nin,nout
c---- temporary arrays to store all 8 child boxes dimensions
      real xmxc(8),ymxc(8),zmxc(8)
      real xmic(8),ymic(8),zmic(8)
      real xcco(8),ycco(8),zcco(8)

c----------------------------------------------------------------------------
c---- boundaries of children boxes


      xmic(1)=xmin
      xmxc(1)=xmin+0.5*(xmax-xmin)
      ymic(1)=ymin
      ymxc(1)=ymin+0.5*(ymax-ymin)
      zmic(1)=zmin
      zmxc(1)=zmin+0.5*(zmax-zmin)

      xmic(2)=xmxc(1)
      xmxc(2)=xmax
      ymic(2)=ymic(1)
      ymxc(2)=ymxc(1)
      zmic(2)=zmic(1)
      zmxc(2)=zmxc(1) 

      xmic(3)=xmic(1)
      xmxc(3)=xmxc(1)
      ymic(3)=ymxc(1)
      ymxc(3)=ymax
      zmic(3)=zmic(1)
      zmxc(3)=zmxc(1)

      xmic(4)=xmxc(1)
      xmxc(4)=xmax
      ymic(4)=ymxc(1)
      ymxc(4)=ymax
      zmic(4)=zmic(1)
      zmxc(4)=zmxc(1)

c      xmic(3)=xmic(2)
c      xmxc(3)=xmxc(2)
c      ymic(3)=ymxc(1)
c      ymxc(3)=ymax
c      zmic(3)=zmic(1)
c      zmxc(3)=zmxc(1) 
 
c      xmic(4)=xmic(1)
c      xmxc(4)=xmxc(1)
c      ymic(4)=ymic(3)
c      ymxc(4)=ymxc(3)
c      zmic(4)=zmic(1)
c      zmxc(4)=zmxc(1)

      xmic(5)=xmic(1)
      xmxc(5)=xmxc(1)
      ymic(5)=ymic(1)
      ymxc(5)=ymxc(1)
      zmic(5)=zmxc(1)
      zmxc(5)=zmax

      xmic(6)=xmic(2)
      xmxc(6)=xmxc(2)
      ymic(6)=ymic(2)
      ymxc(6)=ymxc(2)
      zmic(6)=zmic(5)
      zmxc(6)=zmxc(5)

      xmic(7)=xmic(3)
      xmxc(7)=xmxc(3)
      ymic(7)=ymic(3)
      ymxc(7)=ymxc(3)
      zmic(7)=zmic(5)
      zmxc(7)=zmxc(5)

      xmic(8)=xmic(4)
      xmxc(8)=xmxc(4)
      ymic(8)=ymic(4)
      ymxc(8)=ymxc(4)
      zmic(8)=zmic(5)
      zmxc(8)=zmxc(5)
      
c---- center of children boxes
      do i=1,8
         xcco(i)=0.5*(xmic(i)+xmxc(i))
         ycco(i)=0.5*(ymic(i)+ymxc(i))
         zcco(i)=0.5*(zmic(i)+zmxc(i))
      enddo

c----------------------------------------------------------------------------
c---- sorting of points

c---- first, put all points from the parent box into temporary arrays
      npp=np
      do j=1,np
         xpp(j)=xp(j)
         ypp(j)=yp(j)
         zpp(j)=zp(j)
         kpb(j)=kpa(j)
      enddo

      mii=0	! number of parent's points that have found child boxes
c---- then, for each child box...
      do 10 i=1,8
         nin=0
         nout=0
c------- loop through all parent's points
         do j=1,npp
            if(xpp(j).gt.xmic(i) .and. xpp(j).le.xmxc(i) .and.
     +         ypp(j).gt.ymic(i) .and. ypp(j).le.ymxc(i) .and.
     +         zpp(j).gt.zmic(i) .and. zpp(j).le.zmxc(i)) then
c------------- if the point is in the child box, put it in a reordered array
               nin=nin+1
               mii=mii+1
               xc(mii)=xpp(j)
               yc(mii)=ypp(j)
               zc(mii)=zpp(j)
               kc(mii)=kpb(j)
            else
c------------- if the point is not in the box, put it in a temporary array
               nout=nout+1
               xppp(nout)=xpp(j)
               yppp(nout)=ypp(j)
               zppp(nout)=zpp(j)
               kpc(nout)=kpb(j)
            endif
         enddo

c------- move all unassigned points to the "remaining" arrays
         npp=nout
         do k=1,npp
            xpp(k)=xppp(k)
            ypp(k)=yppp(k)
            zpp(k)=zppp(k)
            kpb(k)=kpc(k)
         enddo

c------- number of points in child box i
         nco(i)=nin

10    continue		! end loop over 8 children

c---- reorder the child box arrays to exclude any boxes with 0 points
      ii=0
      do 20 i=1,8
         if(nco(i).eq.0) goto 20
         ii=ii+1
         nc(ii)=nco(i)
         xcc(ii)=xcco(i)
         ycc(ii)=ycco(i)
         zcc(ii)=zcco(i)
         xminc(ii)=xmic(i)
         yminc(ii)=ymic(i)
         zminc(ii)=zmic(i)
         xmaxc(ii)=xmxc(i)
         ymaxc(ii)=ymxc(i)
         zmaxc(ii)=zmxc(i)
         chid(ii)=i
20    continue		! end loop over 8 children

      kne=ii	! sets kne equal to the number of populated child boxes

      return
      end
