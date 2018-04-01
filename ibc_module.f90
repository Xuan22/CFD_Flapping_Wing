MODULE IBC_MODULE

   USE immersedBoundVars
   USE domain_info

   IMPLICIT NONE

   PRIVATE :: getImmersedBody,calcBounds,hole_cut,getBoundaryNodes,calcNearestNorm
   PRIVATE :: inTri_bary,calcBoundInterp,trilinearInterp_ibc,deallocate_ibcvars
   PRIVATE :: assemble_ibc_mesh,send_iblk_ibc,send_boundnodes,update_ibc_block_iblanks

   PUBLIC :: init_ibc

   CONTAINS

      !----------------------------------------------------------------
      subroutine init_ibc(x,y,z,iblank)
      !
      ! setup immersed boundarys.
      !
      !----------------------------------------------------------------
      implicit none
      !----------------------------------------------------------------
		! incoming variables

		REAL(KIND=rdp),DIMENSION(jmax,kmax,lmax),INTENT(IN) :: x,y,z
		INTEGER(KIND=idp),DIMENSION(jmax,kmax,lmax),INTENT(INOUT) :: iblank

		! local variables

		INTEGER(KIND=idp) :: imin_mesh,jmin_mesh,kmin_mesh,imax_mesh,jmax_mesh,kmax_mesh,num_gpt

		INTEGER :: i,j,k,l

		INTEGER(KIND=idp) :: nbound_zone,n,nn,iint,jint,kint

		REAL(KIND=rdp),ALLOCATABLE,DIMENSION(:) :: xmesh,ymesh,zmesh

		REAL(KIND=rdp) :: dx,dy,dz

                real(kind=rdp) :: half=0.5_rdp

                integer :: idx

                character*32 :: filename

		REAL(KIND=rdp),ALLOCATABLE,DIMENSION(:,:,:) :: xg,yg,zg ! entire IBC grid for master proc
                INTEGER(KIND=idp),ALLOCATABLE,DIMENSION(:,:,:) :: iblk

		! First executable statement

                if (this_mesh%dom_type.ne.DOM_IBC) goto 100

		! Assemble IBC mesh on IBC master proc
                allocate(xg(this_mesh%jm,this_mesh%km,this_mesh%lm))
                allocate(yg(this_mesh%jm,this_mesh%km,this_mesh%lm))
                allocate(zg(this_mesh%jm,this_mesh%km,this_mesh%lm))
                call assemble_ibc_mesh(x,y,z,xg,yg,zg)

                if (.not.this_isMaster) goto 100 ! do IBC only on IBC master proc

                allocate(iblk(this_mesh%jm,this_mesh%km,this_mesh%lm))

		num_gpt = 3!num_ghost_points

      if(.not. allocated(bodyFacet)) call getImmersedBody

      if(.not. ibc_present) return

			imin_mesh = 1
			jmin_mesh = 1
			kmin_mesh = 1
			imax_mesh = this_mesh%jm
			jmax_mesh = this_mesh%km
			kmax_mesh = this_mesh%lm

			! If using cell centerd formulation: Set number of cells on entire mesh (num_cells = num_nodes-1)
			!imax_mesh = imax_mesh - 1
			!jmax_mesh = jmax_mesh - 1
			!kmax_mesh = kmax_mesh - 1

			! Calculate x,y,z of mesh (nodes or cell centers coordinates)
			allocate(xmesh(imax_mesh),ymesh(jmax_mesh),zmesh(kmax_mesh))
			do i=1,imax_mesh
				xmesh(i)=xg(i,1,1)
				!xmesh(i)=half(xg(i,1,1) + xg(i+1,1,1))
			end do
			do j=1,jmax_mesh
				ymesh(j)=yg(1,j,1)
				!ymesh(j)=half(yg(1,j,1) + yg(1,j+1,1))
			end do
			do k=1,kmax_mesh
				zmesh(k)=zg(1,1,k)
				!zmesh(k)=half(zg(1,1,k) + zg(1,1,k+1))
			end do

                        if (xmesh(1)<xb_min.and.xmesh(imax_mesh)>xb_max) print*,'STL body bounded by IBC mesh in X direction on PID:',PID 
                        if (ymesh(1)<yb_min.and.ymesh(jmax_mesh)>yb_max) print*,'STL body bounded by IBC mesh in Y direction on PID:',PID 
                        if (zmesh(1)<zb_min.and.zmesh(kmax_mesh)>zb_max) print*,'STL body bounded by IBC mesh in Z direction on PID:',PID 

                        ! Compute new iblank arrays
			call hole_cut(xmesh,ymesh,zmesh,iblk,imin_mesh,jmin_mesh,kmin_mesh,imax_mesh,jmax_mesh,kmax_mesh)
                        print*,'1 IBC: hole cutting done on PID:',PID

			! determine which iblanked points are immersed boundary points
			call getBoundaryNodes(iblk,nbound,imin_mesh,jmin_mesh,kmin_mesh,imax_mesh,jmax_mesh,kmax_mesh,num_gpt)
                        print*,'2 IBC: boundary nodes found on PID:',PID
			print*,'       number of bound nodes = ',nbound

			! initialize velocities associated with boundary nodes
			! This is not really necessary except at the first iteration.
			boundNode(:)%vel  = 0

			! for ib points, find the nearest normal facet
			call calcNearestNorm(nbound,boundNode,xmesh,ymesh,zmesh,imin_mesh,jmin_mesh,kmin_mesh,imax_mesh,jmax_mesh,kmax_mesh,num_gpt)
                        print*,'3 IBC: nearest norms found on PID:',PID

			call calcBoundInterp(iblk,nbound,boundNode,xmesh,ymesh,zmesh,imin_mesh,jmin_mesh,kmin_mesh,imax_mesh,jmax_mesh,kmax_mesh,num_gpt)
                        print*,'4 IBC: interpolation info found on PID:',PID

			deallocate(xmesh,ymesh,zmesh)

100   call barrier_mpi

      ! Update IBC block iblank (only needed if no overset, otherwise iblk_ibc)
      if (this_mesh%dom_type.eq.DOM_IBC) call update_ibc_block_iblanks(iblk,iblank)

      ! Send IBC iblanks to overset master for IHC routine
      call send_iblk_ibc(iblk)
      call barrier_mpi

      ! Send boundNodes to their respective IBC blocks/procs
      call send_boundnodes()
      call barrier_mpi

      if (this_mesh%dom_type.eq.DOM_IBC) deallocate(xg,yg,zg)
      if (this_mesh%dom_type.eq.DOM_IBC.and.this_isMaster) deallocate(iblk)

      return

      end subroutine init_ibc

      !----------------------------------------------------------------
      subroutine getImmersedBody
      !
      ! read in body stl file, allocate and calculate some bodyFacet 
      ! values
      !
      !----------------------------------------------------------------
      implicit none
      !----------------------------------------------------------------
      ! Local variables
      logical               :: file_exists
      character(LEN = 80)   :: text
      integer(kind=idp) :: idx1,idx2
      real(kind=rdp)   :: sa,sb,sc,ss,x1,x2,x3,y1,y2,y3,z1,z2,z3
      real(kind=rdp)   :: x21,x31,x32,y21,y31,y32,z21,z31,z32
      real(kind=rdp)   :: third
      integer :: bodytype
      ! Math: add stuff
      real(kind=rdp) :: nx,ny,nz,norm_normal

      iforce = zero
      itestMove = zero
      imove = zero


!      if(.not. ibound) return


      ! check to see if bodyfile exists
      bodytype=0 ! 0: generic file, 1: sphere_special
      if (bodytype==0) inquire(file='surface.stl',exist=file_exists)
      if (bodytype==1) inquire(file='sphere_special.stl',exist=file_exists)

      if(.not. file_exists) then
        if (bodytype==0) write(*,*) '*** Input file:  "surface.stl"  not found for immersed boundaries'
        if (bodytype==1) write(*,*) '*** Input file:  "sphere_special.stl"  not found for immersed boundaries'
        ibc_present = .false.
        !return
        stop
      else
         ibc_present = .true.
      endif

      write(*,*) ' *** Reading STL file ...'

      ! open bodyfile
      if (bodytype==0) open(unit=11,file='surface.stl',form='formatted')
      if (bodytype==1) open(unit=11,file='sphere_special.stl',form='formatted')

      ! count facets in STL file     
      nfacet = 0
      do while (text.ne.'endsolid')
        read(11,*) text
        if(text.eq.'facet') nfacet = nfacet+1
      enddo
      rewind(11)

      ! allocate bodyFacet and force arrays
      allocate(bodyFacet(nfacet))

      ! read-in facet data
      read(11,*) text
      do idx1=1,nfacet

        read(11,*) text,text,bodyFacet(idx1)%norm(1),  &
                             bodyFacet(idx1)%norm(2),  &
                             bodyFacet(idx1)%norm(3)
        read(11,*) text
        do idx2=1,3
          read(11,*) text,bodyFacet(idx1)%vertex(idx2,1), &
                          bodyFacet(idx1)%vertex(idx2,2), &
                          bodyFacet(idx1)%vertex(idx2,3)
        enddo
        read(11,*) text
        read(11,*) text

      enddo

      ! close bodyfile
      close(11)


      write(*,*) '     STL file number of facets:',nfacet
      !write(*,*) ''

      third     = 1d0/3d0
      totalArea = 0
      do idx1=1,nfacet
        ! facet coordinates

        if (bodytype==1) then ! resize and move sphere_special
           print*,'Warning: resizing and moving sphere_special'
           bodyFacet(idx1)%vertex(1,1) = bodyFacet(idx1)%vertex(1,1)*8
           bodyFacet(idx1)%vertex(1,2) = bodyFacet(idx1)%vertex(1,2)*8-4
           bodyFacet(idx1)%vertex(1,3) = bodyFacet(idx1)%vertex(1,3)*8-8
           bodyFacet(idx1)%vertex(2,1) = bodyFacet(idx1)%vertex(2,1)*8
           bodyFacet(idx1)%vertex(2,2) = bodyFacet(idx1)%vertex(2,2)*8-4
           bodyFacet(idx1)%vertex(2,3) = bodyFacet(idx1)%vertex(2,3)*8-8
           bodyFacet(idx1)%vertex(3,1) = bodyFacet(idx1)%vertex(3,1)*8
           bodyFacet(idx1)%vertex(3,2) = bodyFacet(idx1)%vertex(3,2)*8-4
           bodyFacet(idx1)%vertex(3,3) = bodyFacet(idx1)%vertex(3,3)*8-8
        endif

        x1 = bodyFacet(idx1)%vertex(1,1)
        y1 = bodyFacet(idx1)%vertex(1,2)
        z1 = bodyFacet(idx1)%vertex(1,3)
        x2 = bodyFacet(idx1)%vertex(2,1)
        y2 = bodyFacet(idx1)%vertex(2,2)
        z2 = bodyFacet(idx1)%vertex(2,3)
        x3 = bodyFacet(idx1)%vertex(3,1)
        y3 = bodyFacet(idx1)%vertex(3,2)
        z3 = bodyFacet(idx1)%vertex(3,3)
        ! Math: check normals
        nx = bodyFacet(idx1)%norm(1)
        ny = bodyFacet(idx1)%norm(2)
        nz = bodyFacet(idx1)%norm(3)    

        ! Make sure normals are unit vectors
        norm_normal = sqrt(nx*nx + ny*ny + nz*nz)
        if (norm_normal>1.05.or.norm_normal<0.95) then
           print*,'Rescaling normal vector:',idx1
           bodyFacet(idx1)%norm(1) = nx/norm_normal
           bodyFacet(idx1)%norm(2) = ny/norm_normal
           bodyFacet(idx1)%norm(3) = nz/norm_normal
        endif

        x21 = x2-x1
        x31 = x3-x1
        x32 = x3-x2
        y21 = y2-y1
        y31 = y3-y1
        y32 = y3-y2
        z21 = z2-z1
        z31 = z3-z1
        z32 = z3-z2

        sa = sqrt(x21**2+y21**2+z21**2)
        sb = sqrt(x31**2+y31**2+z31**2)
        sc = sqrt(x32**2+y32**2+z32**2)
        ss = 0.5d0*(sa+sb+sc)

        ! area of facet
        bodyFacet(idx1)%area = sqrt(ss*(ss-sa)*(ss-sb)*(ss-sc))
        !if (bodyFacet(idx1)%area<10e-5) print*,'Panel area too small:',idx1,bodyFacet(idx1)%area
        totalArea = totalArea + bodyFacet(idx1)%area

        ! centroid of facet
        bodyFacet(idx1)%center(1) = third*(x1+x2+x3)
        bodyFacet(idx1)%center(2) = third*(y1+y2+y3)
        bodyFacet(idx1)%center(3) = third*(z1+z2+z3)

      enddo
      write(*,*) '     STL body surface area =',totalArea


      call calcBounds

      return
      end subroutine getImmersedBody

      !----------------------------------------------------------------
      subroutine calcBounds
      !
      ! Calculate the geometric bounds of the immersed body and store
      ! in xb_min,xb_max,yb_min,etc which are in the immersedBoundVars
      ! module
      !
      !----------------------------------------------------------------
      implicit none
      !----------------------------------------------------------------
      ! Local variables
      integer(kind=idp) :: idx1,idx2

      ! find bounds of the immersed body
      xb_max = bodyFacet(1)%vertex(1,1)
      xb_min = xb_max
      yb_max = bodyFacet(1)%vertex(1,2)
      yb_min = yb_max
      zb_max = bodyFacet(1)%vertex(1,3)
      zb_min = zb_max

      do idx1=1,nfacet
        do idx2=1,3
          xb_max = max(xb_max,bodyFacet(idx1)%vertex(idx2,1))
          xb_min = min(xb_min,bodyFacet(idx1)%vertex(idx2,1))
          yb_max = max(yb_max,bodyFacet(idx1)%vertex(idx2,2))
          yb_min = min(yb_min,bodyFacet(idx1)%vertex(idx2,2))
          zb_max = max(zb_max,bodyFacet(idx1)%vertex(idx2,3))
          zb_min = min(zb_min,bodyFacet(idx1)%vertex(idx2,3))
        enddo
      enddo

      return
      end subroutine calcBounds

      !----------------------------------------------------------------
      subroutine hole_cut(x,y,z,ib,jmin,kmin,lmin,jmax,kmax,lmax)
      !
      !     core routine to perform the hole-cutting
      !
      !     Author  :: Jay Sitaraman
      !     Started :: 09/14/06
      !
      !----------------------------------------------------------------
      !  prologue:
      ! 
      !  {x,y,z}    :: 3-D arrays of cutteee mesh, dim=(jmax,kmax,lmax)
      !  {xt,yt,zt} :: coordinates of vertices of triangles, dim=(3,nfacet)
      !        ib   :: iblank values, dim=(jmax,kmax,lmax)       
      !
      !----------------------------------------------------------------
      implicit none
      !----------------------------------------------------------------
      !
      ! subroutine arguments
      !
      integer(kind=idp)  :: jmin,kmin,lmin,jmax,kmax,lmax
      !real(kind=rdp)  :: ib(jmin:,kmin:,lmin:)
      integer(kind=idp)  :: ib(jmin:,kmin:,lmin:)
      real(kind=rdp)    :: x(jmin:),y(kmin:),z(lmin:)
      !
      ! local variables
      !
      integer(kind=idp)  :: i,itr,j,k,l,ii,jj,inside,maxtri
      integer(kind=idp)  :: js,je,ks,ke,ls,le
      real(kind=rdp)    :: l1,l2,l3
      real(kind=rdp)    :: tmp,yy,zz,eps,den
      real(kind=rdp)    :: fac,xlen,xcm,ylen,ycm,zlen,zcm
      real(kind=rdp)    :: xmax,xmin,ymax,ymin,zmax,zmin
      real(kind=rdp), allocatable  :: cof(:,:),xd(:),dist(:)
      real(kind=rdp)	:: dist_tmp,xp,yp,zp,xf,yf,zf,nx,ny,nz
      !
      ! first executable statement
      !
      !
      ! increase the boundaries by 5 percent for preventing
      ! round off errors
      !
      fac=1.05
      xlen=(xb_max-xb_min)*0.5*fac ! Half length of body (+xx%)
      xcm=(xb_max+xb_min)*0.5      ! Mid point of body
      xmin=xcm-xlen                ! Min point of body
      xmax=xcm+xlen                ! Max point of body

      ylen=(yb_max-yb_min)*0.5*fac
      ycm=(yb_max+yb_min)*0.5
      ymin=ycm-ylen
      ymax=ycm+ylen

      zlen=(zb_max-zb_min)*0.5*fac
      zcm=(zb_max+zb_min)*0.5
      zmin=zcm-zlen
      zmax=zcm+zlen

      !
      ! now find bounds of the coordinates which are contained
      ! inside the bounding box
      !
      js=jmax
      je=1
      ks=kmax
      ke=1
      ls=lmax
      le=1

      do l=1,lmax
         do k=1,kmax
            do j=1,jmax
               
               if ((x(j)-xmin)*(x(j)-xmax).le.0 .and.    &
                   (y(k)-ymin)*(y(k)-ymax).le.0 .and.    &
                   (z(l)-zmin)*(z(l)-zmax).le.0) then

                  je=max(j,je)
                  js=min(j,js)
                  ke=max(k,ke)
                  ks=min(k,ks)
                  le=max(l,le)
                  ls=min(l,ls)

               endif
               
            enddo
         enddo
      enddo

      !
      ! preprocess the projected areas in Y-Z plane
      ! note :: have to change to 
      !         generalize and pick the plane with maximum
      !         exposure
      ! 
      eps=1e-12
      maxtri=500
      allocate(cof(5,nfacet),xd(maxtri),dist(maxtri))

      do i=1,nfacet

         cof(1,i)=bodyFacet(i)%vertex(1,2)-bodyFacet(i)%vertex(3,2)
         cof(2,i)=bodyFacet(i)%vertex(2,2)-bodyFacet(i)%vertex(3,2)
         cof(3,i)=bodyFacet(i)%vertex(1,3)-bodyFacet(i)%vertex(3,3)
         cof(4,i)=bodyFacet(i)%vertex(2,3)-bodyFacet(i)%vertex(3,3)
         den=cof(1,i)*cof(4,i)-cof(2,i)*cof(3,i)
         !
         ! if zero projected area then don't check this triangle
         !
         if (abs(den).ge.eps) then
            cof(5,i)=1./den
         else
            cof(5,i)=0
         endif
      enddo
      !
      ! set iblanks to 1 everywhere first
      !      
      ib=1
      !
      ! now perform the actual iblanking
      !
      do l=ls,le
         do k=ks,ke
            
            itr=0
            do i=1,nfacet
               if (cof(5,i).ne.0) then
                  yy=(bodyFacet(i)%vertex(3,2)-y(k))
                  zz=(bodyFacet(i)%vertex(3,3)-z(l))
                  !     
                  !     area coordinates
                  !     
                  l1=(cof(2,i)*zz-cof(4,i)*yy)*cof(5,i)
                  l2=(cof(3,i)*yy-cof(1,i)*zz)*cof(5,i)
                  l3=1-l1-l2
                  !     
                  !     if convex set the point is inside
                  !
                  if (l1.gt.-eps .and. l2.gt.-eps .and. l3.gt.-eps) then

                     itr=itr+1
                     xd(itr)=l1*bodyFacet(i)%vertex(1,1)  &
                            +l2*bodyFacet(i)%vertex(2,1)  &
                            +l3*bodyFacet(i)%vertex(3,1)  
                     dist(itr)=abs(x(js)-xd(itr))

                  endif
               endif
            enddo

            if (itr.gt.1) then
            !
            !     sort the distance vector
            !
               do ii=1,itr
                  do jj=ii+1,itr
                     if (dist(jj).lt.dist(ii)) then
                        tmp=dist(ii)
                        dist(ii)=dist(jj)
                        dist(jj)=tmp
                        tmp=xd(ii)
                        xd(ii)=xd(jj)
                        xd(jj)=tmp
                     endif
                  enddo
               enddo
            !
            !     cull out same distances
            !
               ii=2
               do while(ii.le.itr)
                  if (abs(xd(ii)-xd(ii-1)).lt.eps) then
                     do jj=ii+1,itr
                        xd(jj-1)=xd(jj)
                        dist(jj-1)=dist(jj)
                     enddo
                     itr=itr-1
                     ii=ii-1
                  endif
                  ii=ii+1
               enddo
                  
               ii=1
               j=js
               inside=0
            !
            !     iblank if inside
            !
               do while(j.le.je)
                  if ((x(j)-xd(ii))*(x(j)-xd(ii+1)).le.eps) then
                     inside=1
                     ib(j,k,l)=0
                  else
                     if (inside.eq.1) then
                        if (ii+2.lt.itr) ii=ii+2
                        inside=0
                     endif
                  endif
                  j=j+1
               enddo
            endif

         enddo
      enddo

      deallocate(cof,xd,dist)

      return
      end subroutine hole_cut

      !----------------------------------------------------------------
      subroutine getBoundaryNodes(ib,nb,imin,jmin,kmin,imax,jmax,kmax,num_gpt)
      !
      ! tag boundary nodes based on the iblank array, ib. number of
      ! boundary nodes is returned in nb. 
      !
      !----------------------------------------------------------------
      implicit none
      !----------------------------------------------------------------
      integer(kind=idp), intent(out) :: nb
      integer(kind=idp), intent(in)  :: imin,jmin,kmin
      integer(kind=idp), intent(in)  :: imax,jmax,kmax,num_gpt
      integer(kind=idp), intent(in)  :: ib(imin:,jmin:,kmin:)
      ! Local variables    
      integer(kind=idp) :: i,j,k,idx,jdx,kdx
      logical               :: test
      real(kind=rdp), allocatable :: bnodes_tmp(:,:)
      ! Math: add aa
      integer(kind=idp) :: aa

      test = .false.
      nb   = 0

		if(num_gpt/=3) then
			print*,'message from getBoundaryNodes in ibc_module.f90 - number of ghost points must be set to 3. terminating simulation'
			stop
		endif
      
     ! note that nodes within the last two rows of the mesh are not
     ! considered boundary nodes. These rows are used as buffer rows
     ! and, aside from block mesh solutions, the body should not be 
     ! allowed to occupy these rows.

      ! Math: add offset for boundary node search
      aa = 0
      do k=kmin+num_gpt+aa,kmax-num_gpt-aa
        do j=jmin+num_gpt+aa,jmax-num_gpt-aa
          do i=imin+num_gpt+aa,imax-num_gpt-aa
        
            if(ib(i,j,k).eq.0) then 

              loopk: do kdx=k-num_gpt,k+num_gpt
                loopj: do jdx=j-num_gpt,j+num_gpt
                  loopi: do idx=i-num_gpt,i+num_gpt

                    if(ib(idx,jdx,kdx).eq.1) then
                      test = .true.
                      exit loopk
                    endif

                  enddo loopi
                enddo loopj
              enddo loopk 

              if(test) then
                test        = .false.
                nb          = nb+1                
              endif          
            endif
            
          enddo
        enddo
      enddo

      allocate(bnodes_tmp(nb,num_gpt))

      nb = 0

      do k=kmin+num_gpt+aa,kmax-num_gpt-aa
        do j=jmin+num_gpt+aa,jmax-num_gpt-aa
          do i=imin+num_gpt+aa,imax-num_gpt-aa
        
            if(ib(i,j,k).eq.0) then  
                
              loopkk: do kdx=k-num_gpt,k+num_gpt
                loopjj: do jdx=j-num_gpt,j+num_gpt
                  loopii: do idx=i-num_gpt,i+num_gpt

                    if(ib(idx,jdx,kdx).eq.1) then
                      test = .true.
                      exit loopkk
                    endif

                  enddo loopii
                enddo loopjj
              enddo loopkk 

              if(test) then
                nb = nb + 1
                test        = .false.
                bnodes_tmp(nb,1) = i
                bnodes_tmp(nb,2) = j
                bnodes_tmp(nb,3) = k
              endif          
            endif
            
          enddo
        enddo
      enddo

      ! allocate boundNode array to correct number of boundary nodes
      ! and store boundary node indices

       if(allocated(boundNode)) deallocate(boundNode)
       allocate(boundNode(nb))
       do idx=1,nb
         boundNode(idx)%i = bnodes_tmp(idx,1)      
         boundNode(idx)%j = bnodes_tmp(idx,2)            
         boundNode(idx)%k = bnodes_tmp(idx,3)
       enddo

      deallocate(bnodes_tmp)

      return
      end subroutine getBoundaryNodes

      !----------------------------------------------------------------
      subroutine calcNearestNorm(nb,bnodes,x,y,z,imin,jmin,kmin,imax,jmax,kmax,num_gpt)
      !
      ! For each boundary node, identify the nearest normal distance 
      ! to one of the facets in the immersed body. 
      !
      !----------------------------------------------------------------
      implicit none
      !----------------------------------------------------------------
      integer(kind=idp) :: nb,imin,jmin,kmin,imax,jmax,kmax,num_gpt
      real(kind=rdp)   :: x(imin:),y(jmin:),z(kmin:)
      type(BoundaryNode)    :: bnodes(1:)
      ! Local variables    
      integer(kind=idp) :: i,j,k,idx,ii,ii_min
      real(kind=rdp)    :: dx,dy,dz,ds,zptmp,yptmp,xptmp,nx,ny,nz
      real(kind=rdp)    :: xp,yp,zp,x1,y1,z1,x2,z2,y2,x3,z3,y3
      real(kind=rdp)    :: eps,dist,dist_tmp
      logical           :: test

      ! Math: add control point
      real(kind=rdp)    :: xcp,ycp,zcp,dist_norm
      integer(kind=idp) :: nacc


      eps            = 10d0**(-13)    
      nacc = 0

      forEachBoundNode: do idx=1,nb   
      
        ! boundary node mesh indices
        i = bnodes(idx)%i
        j = bnodes(idx)%j
        k = bnodes(idx)%k       
              
        ! boundary node coords
        xp = x(i)
        yp = y(j)
        zp = z(k)  

        ! calculate local grid spacing at boundary node
        dx = (x(i+1)-x(i-1))
        dy = (y(j+1)-y(j-1))  
        dz = (z(k+1)-z(k-1))  
        ds = 0.5*max(dx,dy,dz)            

        ! some initial large distance
        dist      = max(x(imax)-x(1),y(jmax)-y(1),z(kmax)-z(1)) 
        ii_min    = -1
        

        forEachFacet: do ii=1,nfacet
          
          ! Math: check panel area, disregard if too small (normals are likely to be wrong)
          if (bodyFacet(ii)%area>10e-6) then
          ! facet node coords
          x1 = bodyFacet(ii)%vertex(1,1)
          y1 = bodyFacet(ii)%vertex(1,2)
          z1 = bodyFacet(ii)%vertex(1,3)  
          x2 = bodyFacet(ii)%vertex(2,1)
          y2 = bodyFacet(ii)%vertex(2,2)
          z2 = bodyFacet(ii)%vertex(2,3) 
          x3 = bodyFacet(ii)%vertex(3,1)
          y3 = bodyFacet(ii)%vertex(3,2)
          z3 = bodyFacet(ii)%vertex(3,3) 
          nx = bodyFacet(ii)%norm(1)
          ny = bodyFacet(ii)%norm(2)
          nz = bodyFacet(ii)%norm(3)    

          ! Math: find control point (center of panel)
          xcp = bodyFacet(ii)%center(1) 
          ycp = bodyFacet(ii)%center(2) 
          zcp = bodyFacet(ii)%center(3) 

          ! Math: use control point to compute distance as some panels may be skewed
          dist_tmp = sqrt((xp-xcp)**2+(yp-ycp)**2+(zp-zcp)**2) 
          !dist_tmp = nx*(xp-xcp)+ny*(yp-ycp)+nz*(zp-zcp) 


          ! since all tested points are inside body and normals point
          ! outwards, dist should always be negative. error if dist>0.
!	This is true only for very simple, convex shapes... not true for more
!	complex bodies... so over-ruled.
          !if(dist_tmp.gt.eps) then
             !write(*,*) '*** Error. Bound node found OUTSIDE of body',dist_tmp
             !stop
     !!!       cycle
          !endif

          
          ! bound node must have a neighbor that's outside of body 
          ! thus the triangle should be within a grid space in the 
          ! direction of the facet normal  !! Not completely tested !!
          !if(abs(dist_tmp) .gt. &
          !   0.5d0*(abs(dx*nx)+abs(dy*ny)+abs(dz*nz)) ) cycle            


          ! if dist to facet ii is less than current min dist, determine
          ! if normal from boundary node to plane intersects the facet
          ! within facet bounds
          if(abs(dist_tmp).lt.abs(dist)) then
                 
            dist_norm = nx*(xp-xcp)+ny*(yp-ycp)+nz*(zp-zcp) 
            xptmp = xp + abs(dist_norm)*nx  ! define intersection of normal line
            yptmp = yp + abs(dist_norm)*ny  ! from point to facet
            zptmp = zp + abs(dist_norm)*nz  !

            ! is [xptmp,yptmp,zptmp] inside facet?
            if(abs(nx).gt.eps) then
              call inTri_bary(zptmp,yptmp,z1,y1,z2,y2,z3,y3,test)
            elseif(abs(nz).gt.eps) then                  
              call inTri_bary(xptmp,yptmp,x1,y1,x2,y2,x3,y3,test)
            else
              call inTri_bary(zptmp,xptmp,z1,x1,z2,x2,z3,x3,test) 
            endif

            if(test) then ! If normal projection is inside facet
              dist   = dist_tmp
              ii_min = ii
            endif     
            
          endif             
          
          if(abs(dist_norm).le.eps) then
             print*,'distance to surface smaller than 10-13, exiting loop',dist
             exit forEachFacet
          endif
          endif

        enddo forEachFacet
        
        ! check for dist that are more than 10X max local grid spacing
        ! This absolutely SHOULD NOT happen.
        if(abs(dist).gt.10.0*ds .and. ii_min.ne.-1) then
           !print*,'Distance to great, resetting:',idx
           ii_min = -1
        endif
        if (ii_min.ne.-1) nacc = nacc + 1
        ! check for bound nodes where no normal is found
        ! just use closest facet, even if normal projection isn't inside facet

	if(ii_min.eq.-1) then
              !print*,'bound node:',idx,', no normal found, fixing without checking bary'
		! some initial large distance
		dist = max(x(imax)-x(1),y(jmax)-y(1),z(kmax)-z(1)) 

		do ii=1,nfacet
                   ! Math: check panel area, disregard if too small
                   if (bodyFacet(ii)%area>10e-6) then
			! facet node coords
			x1 = bodyFacet(ii)%vertex(1,1)
			y1 = bodyFacet(ii)%vertex(1,2)
			z1 = bodyFacet(ii)%vertex(1,3)  
			x2 = bodyFacet(ii)%vertex(2,1)
			y2 = bodyFacet(ii)%vertex(2,2)
			z2 = bodyFacet(ii)%vertex(2,3) 
			x3 = bodyFacet(ii)%vertex(3,1)
			y3 = bodyFacet(ii)%vertex(3,2)
			z3 = bodyFacet(ii)%vertex(3,3) 
			nx = bodyFacet(ii)%norm(1)
			ny = bodyFacet(ii)%norm(2)
			nz = bodyFacet(ii)%norm(3)    

                        ! Math: find control point (center of panel)
                        xcp = bodyFacet(ii)%center(1) 
                        ycp = bodyFacet(ii)%center(2) 
                        zcp = bodyFacet(ii)%center(3) 

			! find normal dist from boundary node to facet ii
                        ! Math: use control point for distance
			dist_tmp = sqrt((xp-xcp)**2+(yp-ycp)**2+(zp-zcp)**2) 
                        !dist_tmp = nx*(xp-xcp)+ny*(yp-ycp)+nz*(zp-zcp) 

			if(abs(dist_tmp).lt.abs(dist)) then
                                dist_norm = nx*(xp-xcp)+ny*(yp-ycp)+nz*(zp-zcp) 
                                xptmp = xp + abs(dist_norm)*nx  ! define intersection of normal line
                                yptmp = yp + abs(dist_norm)*ny  ! from point to facet
                                zptmp = zp + abs(dist_norm)*nz  !
			        dist   = dist_tmp
			        ii_min = ii
			endif
                   endif  
          	enddo             
          
	endif     
 
	if(ii_min.eq.-1) print*,'bound node:',idx,', stil no normal found' ! should not happen, remove this check later
                    
        ! for each boundary point, store nearest normal and dist
        bnodes(idx)%norm(1) = bodyfacet(ii_min)%norm(1)
        bnodes(idx)%norm(2) = bodyfacet(ii_min)%norm(2)
        bnodes(idx)%norm(3) = bodyfacet(ii_min)%norm(3)

        ! point of normal intersection with panel
        bnodes(idx)%surf(1) = xptmp
        bnodes(idx)%surf(2) = yptmp
        bnodes(idx)%surf(3) = zptmp
        
        bnodes(idx)%normFacet = ii_min 

        ! if point is on surface, set dist = -0.25*min(dx,dy,dz)  (arbitrary negative val)
        if(abs(dist_norm).gt.eps) then
          bnodes(idx)%d2surf  = abs(dist_norm)
        else  
          print*,'Bound node is on surface, fixing normal distance:',idx
          bnodes(idx)%d2surf  = -0.25*min(dx,dy,dz)
        endif

      enddo forEachBoundNode

      print*,'IBC solution',real(nacc)/real(nb)*100.0,'% accurate (based on facet intersection)',nacc

      return      
      end subroutine calcNearestNorm

      !----------------------------------------------------------------
      subroutine inTri_bary(zp,yp,z1,y1,z2,y2,z3,y3,test)
      ! given a point [pz,py], determine if it lies inside the 2D 
      ! triangle defined by the points [v1z,v1y], [v2z,v2y], [v3z,v3y].
      !----------------------------------------------------------------
      !use constants
      implicit none      
      !----------------------------------------------------------------
      real(kind=rdp)   :: zp,yp,z1,y1,z2,y2,z3,y3,eps_dp
      logical               :: test
      ! local variables      
      real(kind=rdp)   :: c1,c2,c3,c4,c5,b1,b2
      real(kind=rdp)   :: L1,L2,L3
      
      
      ! vector components and cross product
      b1 = z3-zp
      b2 = y3-yp

      c1 = y1-y3
      c2 = y2-y3
      c3 = z1-z3
      c4 = z2-z3
      c5 = 1./(c4*c1-c2*c3)
      
      ! barycentric coordinates
      L1 = (c2*b1 - c4*b2)*c5
      L2 = (c3*b2 - c1*b1)*c5
      L3 = 1-L1-L2
       
      eps_dp = 10**(-14.d0)
      ! if convex set, point is inside
      if(L1.ge.-eps_dp .and. L2.ge.-eps_dp .and. L3.ge.-eps_dp) then
        test = .true.
      else
        test = .false.
      endif

      return
      end subroutine inTri_bary

      !----------------------------------------------------------------
      subroutine calcBoundInterp(ib,nb,bnodes,x,y,z,imin,jmin,kmin,imax,jmax,kmax,num_gpt)
      !
      ! For each boundary node, calculate location for field 
      ! interpolation node and calculate the interpolation coeffs of 
      ! surrounding 8 points using trilinear interpolation
      !
      !----------------------------------------------------------------
      implicit none
      !----------------------------------------------------------------
      integer(kind=idp) :: nb,imin,jmin,kmin,imax,jmax,kmax,num_gpt
      integer(kind=idp) :: ib(imin:,jmin:,kmin:)
      real(kind=rdp)   :: x(imin:),y(jmin:),z(kmin:)
      type(BoundaryNode)    :: bnodes(1:)
      ! Local variables   
      integer(kind=idp) :: i,j,k,bi,bj,bk,idx,nflow
	integer(kind=idp) :: itip,jtip,ktip
	integer(kind=idp) :: idiff, jdiff,kdiff
      real(kind=rdp)   :: nx,ny,nz,xb,yb,zb,ds,dz,dy,dx,dist
      real(kind=rdp)   :: xtip,ytip,ztip,dtip,fact      
      logical          :: skip,verboseWarning
	logical	:: is_it_in
	integer :: temp1

	print*, 'PID:',PID,'Number of boundary nodes = ',nb

	forEachBoundNode: do idx=1,nb
		i = bnodes(idx)%i
		j = bnodes(idx)%j
		k = bnodes(idx)%k
		xb = x(i)
		yb = y(j)
		zb = z(k)
		nx = bnodes(idx)%norm(1)
		ny = bnodes(idx)%norm(2)
		nz = bnodes(idx)%norm(3)
		dist = bnodes(idx)%d2surf
		dx = 0.5 * (x(i+1)-x(i-1))
		dy = 0.5 * (y(j+1)-y(j-1))
		dz = 0.5 * (z(k+1)-z(k-1))

		ds = 0.01d0*min(dx,dy,dz)

		xtip = xb + abs(dist)*nx
		ytip = yb + abs(dist)*ny
		ztip = zb + abs(dist)*nz
		dtip = abs(dist)


		if (dtip .lt. 0.0_8) then
			print *, 'ERROR: Bug in code - dtip should not be negative at this point!'
			stop
		endif

		is_it_in = .false.
		isitin:	do while (.not. is_it_in)
			itip = i
			jtip = j
			ktip = k

			if (xtip > x(i)) then
				do while (x(itip) < xtip)
					itip = itip+1
				enddo
			else
				do while (x(itip-1) > xtip)
					itip = itip-1
				enddo
			endif

			if (ytip > y(j)) then
				do while (y(jtip) < ytip)
					jtip=jtip+1
				enddo
			else
				do while (y(jtip-1) > ytip)
					jtip=jtip-1
				enddo
			endif

			if (ztip > z(k)) then
				do while (z(ktip) < ztip)
					ktip = ktip+1
				enddo
			else
				do while (z(ktip-1) > ztip)
					ktip = ktip-1
				enddo
			endif

			nflow = 0
			nflow = nflow + ib(itip,jtip,ktip)
			nflow = nflow + ib(itip-1,jtip,ktip)
			nflow = nflow + ib(itip,jtip-1,ktip)
			nflow = nflow + ib(itip,jtip,ktip-1)
			nflow = nflow + ib(itip-1,jtip-1,ktip)
			nflow = nflow + ib(itip,jtip-1,ktip-1)
			nflow = nflow + ib(itip-1,jtip,ktip-1)
			nflow = nflow + ib(itip-1,jtip-1,ktip-1)

			if (nflow.eq.8) then
				is_it_in = .true.
			elseif (nflow.lt.8) then
				is_it_in = .false.
				xtip = xtip + nx*abs(ds)
				ytip = ytip + ny*abs(ds)
				ztip = ztip + nz*abs(ds)
				dtip = dtip + abs(ds)
			else
				print*, 'ERROR: Bug in code - counting interior points surrounding probe tip. at (j,k,l) = ',i,j,k,'PID: ',PID
                                stop
			endif
		enddo isitin

		!Now it's in
		bnodes(idx)%iint = itip-1
		bnodes(idx)%jint = jtip-1
		bnodes(idx)%kint = ktip-1


		if (boundNode(idx)%iint < imin-num_gpt) then
			write(21,*) 'WARNING: Found IBC Interpolation cells that are out of bounds. Forcefully readjusting.', idx
			boundNode(idx)%iint = imin-num_gpt
		elseif (boundNode(idx)%iint+1 > imax+num_gpt) then
			write(21,*) 'WARNING: Found IBC Interpolation cells that are out of bounds. Forcefully readjusting.', idx
			boundNode(idx)%iint = imax+num_gpt-1
		endif

		if (boundNode(idx)%jint < jmin-num_gpt) then
			write(21,*) 'WARNING: Found IBC Interpolation cells that are out of bounds. Forcefully readjusting.', idx
			boundNode(idx)%jint = jmin-num_gpt
		elseif (boundNode(idx)%jint+1 > jmax+num_gpt) then
			write(21,*) 'WARNING: Found IBC Interpolation cells that are out of bounds. Forcefully readjusting.', idx
			boundNode(idx)%jint = jmax+num_gpt-1
		endif

		if (boundNode(idx)%kint < kmin-num_gpt) then
			write(21,*) 'WARNING: Found IBC Interpolation cells that are out of bounds. Forcefully readjusting.', idx
			boundNode(idx)%kint = kmin-num_gpt
		elseif (boundNode(idx)%kint+1 > kmax+num_gpt) then
			write(21,*) 'WARNING: Found IBC Interpolation cells that are out of bounds. Forcefully readjusting.', idx
			boundNode(idx)%kint = kmax+num_gpt-1
		endif

		if(dist.gt.0) then
			bnodes(idx)%d2int_ratio  = dtip/(dtip+dist)
			bnodes(idx)%d2surf_ratio = dist/(dtip+dist)
		else
			bnodes(idx)%d2int_ratio  = 1
			bnodes(idx)%d2surf_ratio = 0
		endif


		call trilinearInterp_ibc(x(itip-1),y(jtip-1),z(ktip-1),x(itip),y(jtip),z(ktip), &
                             xtip,ytip,ztip,bnodes(idx)%interp(:))
	enddo forEachBoundNode

      return 
      end subroutine calcBoundInterp

!----------------------------------------------------------------------
      subroutine trilinearInterp_ibc(x1,y1,z1,x2,y2,z2,x,y,z,coeff)
!
! This subroutine takes the coordinates of 2 points in which define a
! cell in a Cartesian mesh and returns the 8 interpolation coefficients
! needed to approximate the value of the point [x,y,z]. Components of
! point 2 should be greater than those of point 1.
!
! If the point x1,y1,z1 is indexed i,j,k then the coeffs returned are:
!   coeff(1) =>  ( i,j,k )      
!   coeff(2) =>  ( i+1,j,k )    
!   coeff(3) =>  ( i,j+1,k )    
!   coeff(4) =>  ( i+1,j+1,k )  
!   coeff(5) =>  ( i,j,k+1 )    
!   coeff(6) =>  ( i+1,j,k+1 )  
!   coeff(7) =>  ( i,j+1,k+1 )  
!   coeff(8) =>  ( i+1,j+1,k+1 )
!           
!----------------------------------------------------------------------
      !use Precision_Def
      implicit none
!----------------------------------------------------------------------
      real(kind=rdp), intent(in)  :: x1,y1,z1,x2,y2,z2,x,y,z
      real(kind=rdp), intent(out) :: coeff(8)
      ! Local variables
      real(kind=rdp) :: voli,dx1,dx2,dy1,dy2,dz1,dz2

      voli = 1.d0/( (x2-x1)*(y2-y1)*(z2-z1) )
      dx1  = x-x1
      dx2  = x2-x
      dy1  = y-y1
      dy2  = y2-y
      dz1  = z-z1
      dz2  = z2-z
                                   ! interpolation for node:
      coeff(1) = dz2*dy2*dx2*voli  ! ( i,j,k )
      coeff(2) = dz2*dy2*dx1*voli  ! ( i+1,j,k )     
      coeff(3) = dz2*dy1*dx2*voli  ! ( i,j+1,k )    
      coeff(4) = dz2*dy1*dx1*voli  ! ( i+1,j+1,k )   
      coeff(5) = dz1*dy2*dx2*voli  ! ( i,j,k+1 )   
      coeff(6) = dz1*dy2*dx1*voli  ! ( i+1,j,k+1 ) 
      coeff(7) = dz1*dy1*dx2*voli  ! ( i,j+1,k+1 )
      coeff(8) = dz1*dy1*dx1*voli  ! ( i+1,j+1,k+1 )


      return
      end subroutine trilinearInterp_ibc


      subroutine deallocate_ibcvars()

         deallocate(bodyFacet)

      end subroutine deallocate_ibcvars

   ! Math: assemble IBC block on IBC master
   subroutine assemble_ibc_mesh(x,y,z,xg,yg,zg)

      implicit none

      real(kind=rdp), dimension(jmax,kmax,lmax), intent(in) :: x, y, z
      real(kind=rdp), dimension(this_mesh%jm,this_mesh%km,this_mesh%lm), intent(inout) :: xg, yg, zg

      integer :: j, k, l, bufsize, i, tag, n
      integer :: jj, kk, ll
      real(kind=rdp), dimension(:), allocatable :: xbuffer
      integer, dimension(:), pointer :: bids
      integer, dimension(:,:), pointer :: bst, bend

      nullify(bids,bst,bend)
      tag=101
      bufsize=b_jm*b_km*b_lm*3
      allocate(xbuffer(bufsize))

      call barrier_mpi
      master_check: if (this_isMaster) then 
         bids=>this_mesh%block_id
         bst=>this_mesh%start
         bend=>this_mesh%end
         mesh_blocks: do n=1,this_mesh%nblocks
            if (bids(n) == (PID+1)) then 
               !Local so copy data
               do l=bst(LDIR,n),bend(LDIR,n)
                  ll=send_lims(LDIR)+(l-bst(LDIR,n))
                  do k=bst(KDIR,n),bend(KDIR,n)
                     kk=send_lims(KDIR)+(k-bst(KDIR,n))
                     do j=bst(JDIR,n),bend(JDIR,n)
                        jj=send_lims(JDIR)+(j-bst(JDIR,n))
                        xg(j,k,l)=x(jj,kk,ll)
                        yg(j,k,l)=y(jj,kk,ll)
                        zg(j,k,l)=z(jj,kk,ll)
                     end do
                  end do
               end do
            else ! Receive data from remote processor
               call mpi_recv(xbuffer,bufsize,REAL_TYPE,bids(n)-1,tag,&
                     DEFAULT_COMM,stats_mpi,ierr)
               i=0
               do l=bst(LDIR,n),bend(LDIR,n)
                  do k=bst(KDIR,n),bend(KDIR,n)
                     do j=bst(JDIR,n),bend(JDIR,n)
                        xg(j,k,l)=xbuffer(i+1)
                        yg(j,k,l)=xbuffer(i+2)
                        zg(j,k,l)=xbuffer(i+3)
                        i=i+3
                     end do
                  end do
              end do
            end if
         end do mesh_blocks
      else  master_check ! We will send information
         i=0
         do l=send_lims(LDIR),send_lims(LDIR+3)
            do k=send_lims(KDIR),send_lims(KDIR+3)
               do j=send_lims(JDIR),send_lims(JDIR+3)
                  xbuffer(i+1)=x(j,k,l)
                  xbuffer(i+2)=y(j,k,l)
                  xbuffer(i+3)=z(j,k,l)
                  i=i+3
               end do
            end do
         end do
         call mpi_bsend(xbuffer,bufsize,REAL_TYPE,(this_mesh%master_id-1),&
               tag,DEFAULT_COMM,ierr)
      end if master_check

      call barrier_mpi
      deallocate(xbuffer)
   end subroutine assemble_ibc_mesh

   ! Math: send iblank array from IBC master to overset master
   subroutine send_iblk_ibc(iblk)

      use domain_connectivity

      implicit none

      INTEGER(KIND=idp),DIMENSION(:,:,:), intent(in) :: iblk
      integer :: ii, jj, kk, mm, mid, o_master

      integer :: bufsize, i, j, k, l, tag
      real(kind=rdp), dimension(:), allocatable :: xbuffer

      call barrier_mpi
      if (.not. this_isMaster) go to 100

      tag=201
      bufsize=m_jm*m_km*m_lm
      allocate(xbuffer(bufsize))

      ovset_groups: do ii=1,o_ngroups
         jj=o_grpptr(ii)
         kk=this_ovset(ii)
         o_master=ovset(kk)%master_id-1 ! overset master

         master_proc: if (this_mesh%dom_type==DOM_IBC) then
            ! I am the IBC mesh master, I should send iblk to overset master (if I'm not it)
            ovset_master: if (PID.ne.o_master) then
               !! Copy iblk into the MPI send buffer
               i=0
               do l=1,this_mesh%lm
                  do k=1,this_mesh%km
                     do j=1,this_mesh%jm
                        xbuffer(i+1)=iblk(j,k,l)
                        i=i+1
                     end do
                  end do
               end do
               bufsize=this_mesh%jm*this_mesh%km*this_mesh%lm
               call mpi_bsend(xbuffer,bufsize,INT_TYPE,o_master,tag,&
                     DEFAULT_COMM,ierr)
            else ovset_master
               print*,'setting iblk on overset master (==IBC master)'
               if(allocated(iblk_ibc)) deallocate(iblk_ibc)
               allocate(iblk_ibc(this_mesh%jm,this_mesh%km,this_mesh%lm))
               iblk_ibc=iblk
               print*,'set iblk on overset master (==IBC master)'
            end if ovset_master
         else master_proc
            ovset_master2: if (PID==o_master) then ! I'm the overset master but not IBC mesh master, I have to receive iblk of IBC mesh
               group_meshes: do mm=1,ovset(kk)%nmesh
                  mid=ovset(kk)%mesh_ids(mm)
                  o_master=meshes(mid)%master_id-1
                  ibc_mesh: if (meshes(mid)%dom_type==DOM_IBC) then ! this is the IBC mesh master, receive iblk from it
                     bufsize=meshes(mid)%jm*meshes(mid)%km*meshes(mid)%lm
                     call mpi_recv(xbuffer,bufsize,INT_TYPE,o_master,tag, &
                           DEFAULT_COMM,stats_mpi,ierr)
                     !! Copy iblk from the MPI send buffer
                     if(allocated(iblk_ibc)) deallocate(iblk_ibc)
                     allocate(iblk_ibc(meshes(mid)%jm,meshes(mid)%km,meshes(mid)%lm))
                     i=0
                     do l=1,meshes(mid)%lm
                        do k=1,meshes(mid)%km
                           do j=1,meshes(mid)%jm
                              iblk_ibc(j,k,l)=xbuffer(i+1)
                              i=i+1
                           end do
                        end do
                     end do
                  end if ibc_mesh
               end do group_meshes
            end if ovset_master2
         end if master_proc
      end do ovset_groups

      deallocate(xbuffer)
100   call barrier_mpi

   end subroutine send_iblk_ibc

   ! Math: send boundNode from IBC master to IBC blocks if doing parallel IBC
   subroutine send_boundnodes()

      !use domain_connectivity

      implicit none

      integer :: idx,bi,bj,bk,i,j,k,n
      integer, dimension(:,:), pointer :: bst, bend
      integer, dimension(:), pointer :: bids
      integer, allocatable, dimension(:,:) :: blk_lim

      integer(kind=idp), dimension(:), allocatable :: ibuffer
      real(kind=rdp), dimension(:), allocatable :: xbuffer

      nullify(bids,bst,bend)

      ! Determine if interpolation nodes are within bounds of their respective blocks
      ! If yes then can do parallel IBC
      ! Else do IBC interpolation on IBC master at every time step
      parallel_ibc = -1
      if (this_mesh%dom_type.ne.DOM_IBC) goto 100
      if (.not.this_isMaster) go to 101
      if (nbound==0) go to 101 ! If no bound nodes were found, parallel_ibc is -1

      ! Get block indices (without gost points)
      bids=>this_mesh%block_id
      bst=>this_mesh%start
      bend=>this_mesh%end
      allocate(blk_lim(6,this_mesh%nblocks))
      call determine_send_limits

      do idx=1,nbound ! Go through boundary nodes
 
         ! boundary node indices
         bi = boundNode(idx)%i
         bj = boundNode(idx)%j
         bk = boundNode(idx)%k
         ! leading interpolation node indices
         i = boundNode(idx)%iint
         j = boundNode(idx)%jint
         k = boundNode(idx)%kint
 
         do n=1,this_mesh%nblocks ! go through each block and check if boundnode belongs to that block
            if ((bi>=bst(JDIR,n).and.bi<=bend(JDIR,n)) .and. &
                (bj>=bst(KDIR,n).and.bj<=bend(KDIR,n)) .and. &
                (bk>=bst(LDIR,n).and.bk<=bend(LDIR,n))) then
               if ((i<blk_lim(JDIR,n).or.i+1>blk_lim(JDIR+3,n)) .or. &
                   (j<blk_lim(KDIR,n).or.j+1>blk_lim(KDIR+3,n)) .or. &
                   (k<blk_lim(LDIR,n).or.k+1>blk_lim(LDIR+3,n))) then
                  parallel_ibc = 0 ! Interpolation node out of bound of block so can't do parallel IBC
                  cycle
               end if
            end if
         end do
         if (parallel_ibc==0) cycle
      end do
      if (parallel_ibc==0) then
         print*,'!!!!! Parallel IBC not possible (interpolation node out of bound of block)'
      else
         parallel_ibc = 1
         print*,'!!!!! Parallel IBC possible!!!'
      end if
               
101   call barrier_mpi

      ! Send parallel_ibc (and boundNode if parallel) to blocks
      master_check: if (this_isMaster) then ! IBC master, send info 
         if (parallel_ibc==1) then ! Need to send boundNode, so create buffer
            allocate(ibuffer(nbound*6),xbuffer(nbound*10))
            do idx=1,nbound ! Go through boundary nodes
               ibuffer(6*(idx-1)+1) = boundNode(idx)%i
               ibuffer(6*(idx-1)+2) = boundNode(idx)%j
               ibuffer(6*(idx-1)+3) = boundNode(idx)%k
               ibuffer(6*(idx-1)+4) = boundNode(idx)%iint
               ibuffer(6*(idx-1)+5) = boundNode(idx)%jint
               ibuffer(6*(idx-1)+6) = boundNode(idx)%kint
               do n=1,8
                  xbuffer(10*(idx-1)+n) = boundNode(idx)%interp(n)
               end do
               xbuffer(10*(idx-1)+9)  = boundNode(idx)%d2int_ratio
               xbuffer(10*(idx-1)+10) = boundNode(idx)%d2surf_ratio
            end do
         end if
            
         mesh_blocks: do n=1,this_mesh%nblocks
            if (bids(n).ne.(PID+1)) then 
               ! Send parallel_ibc
               call mpi_bsend(parallel_ibc,1,INT_TYPE,bids(n)-1,&
                     101,DEFAULT_COMM,ierr)

               if (parallel_ibc==1) then
                  call mpi_bsend(nbound,1,INT_TYPE,bids(n)-1,&
                        102,DEFAULT_COMM,ierr)
                  call mpi_bsend(ibuffer,nbound*6,INT_TYPE,bids(n)-1,&
                        103,DEFAULT_COMM,ierr)
                  call mpi_bsend(xbuffer,nbound*10,REAL_TYPE,bids(n)-1,&
                        104,DEFAULT_COMM,ierr)
               end if

            end if
         end do mesh_blocks
      else  master_check ! IBC blocks, receive info
         call mpi_recv(parallel_ibc,1,INT_TYPE,this_mesh%master_id-1,&
               101,DEFAULT_COMM,stats_mpi,ierr)

         if (parallel_ibc==1) then
            call mpi_recv(nbound,1,INT_TYPE,this_mesh%master_id-1,&
                  102,DEFAULT_COMM,stats_mpi,ierr)
            allocate(ibuffer(nbound*6),xbuffer(nbound*10))
            call mpi_recv(ibuffer,nbound*6,INT_TYPE,this_mesh%master_id-1,&
                  103,DEFAULT_COMM,stats_mpi,ierr)
            call mpi_recv(xbuffer,nbound*10,REAL_TYPE,this_mesh%master_id-1,&
                  104,DEFAULT_COMM,stats_mpi,ierr)
            if(allocated(boundNode)) deallocate(boundNode)
            allocate(boundNode(nbound))
            do idx=1,nbound ! Set boundNode info from buffer
               boundNode(idx)%i = ibuffer(6*(idx-1)+1)
               boundNode(idx)%j = ibuffer(6*(idx-1)+2)
               boundNode(idx)%k = ibuffer(6*(idx-1)+3)
               boundNode(idx)%iint = ibuffer(6*(idx-1)+4)
               boundNode(idx)%jint = ibuffer(6*(idx-1)+5)
               boundNode(idx)%kint = ibuffer(6*(idx-1)+6)
               do n=1,8
                  boundNode(idx)%interp(n) = xbuffer(10*(idx-1)+n)
               end do
               boundNode(idx)%d2int_ratio  = xbuffer(10*(idx-1)+9)
               boundNode(idx)%d2surf_ratio = xbuffer(10*(idx-1)+10)
            end do
         end if

      end if master_check

      if (parallel_ibc==1) deallocate(ibuffer,xbuffer)

100   call barrier_mpi

   contains
      subroutine determine_send_limits
         !! Determine the lower and upper j,k,l limits of each block
         !! taking ghost points into account. This information isn't
         !! stored in the master processor, so we need to calculate
         !! them.
         integer, dimension(6) :: extents
         integer, dimension(3) :: maxes

         maxes=(/ this_mesh%jm, this_mesh%km, this_mesh%lm /)

         do n=1,this_mesh%nblocks
            extents=(/ bst(JDIR,n), bst(KDIR,n), bst(LDIR,n), &
                  bend(JDIR,n), bend(KDIR,n), bend(LDIR,n) /)

            do i=1,3
               if (extents(i) == 1) then 
                  blk_lim(i,n) = extents(i)
               else
                  blk_lim(i,n) = extents(i)-NGHOST
               end if
            end do
            
            ! The max faces
            do i=4,6
               if (extents(i) == maxes(i-3)) then 
                  blk_lim(i,n) = extents(i)
               else
                  blk_lim(i,n) = extents(i)+NGHOST
               end if
            end do
         end do
      end subroutine determine_send_limits
   end subroutine send_boundnodes

   subroutine update_ibc_block_iblanks(iblk,iblank)
      !! Relay per-block iblank information to the individual blocks
      !! from the mesh masters. Here we also send the iblank
      !! information for the ghost points as well.
      integer, dimension(this_mesh%jm,this_mesh%km,this_mesh%lm), intent(in) :: iblk
      integer, dimension(jmax,kmax,lmax), intent(inout) :: iblank

      integer :: j, k, l, bufsize, i, tag, n
      integer :: j1, k1, l1
      integer, dimension(:), allocatable :: iblbuf
      integer, dimension(:), pointer :: bids
      integer, dimension(:,:), pointer :: bst, bend
      integer, dimension(6) :: blk_lim

      tag=101
      bufsize=b_jm*b_km*b_lm
      allocate(iblbuf(bufsize))

      call barrier_mpi
      master_check: if (this_isMaster) then
         bids=>this_mesh%block_id
         bst=>this_mesh%start
         bend=>this_mesh%end
         mesh_blocks: do n=1,this_mesh%nblocks
            if (bids(n) == (PID+1)) then
               do l=jkl_lims(LDIR),jkl_lims(LDIR+3)
                  l1=l-jkl_lims(LDIR)+1
                  do k=jkl_lims(KDIR),jkl_lims(KDIR+3)
                     k1=k-jkl_lims(KDIR)+1
                     do j=jkl_lims(JDIR),jkl_lims(JDIR+3)
                        j1=j-jkl_lims(JDIR)+1
                        iblank(j1,k1,l1)=iblk(j,k,l)
                     end do
                  end do
               end do
            else
               call determine_send_limits
               i=0
               do l=blk_lim(LDIR),blk_lim(LDIR+3)
                  do k=blk_lim(KDIR),blk_lim(KDIR+3)
                     do j=blk_lim(JDIR),blk_lim(JDIR+3)
                        i=i+1
                        iblbuf(i)=iblk(j,k,l)
                     end do
                  end do
               end do
               call mpi_bsend(iblbuf,bufsize,INT_TYPE,bids(n)-1,tag, &
                     DEFAULT_COMM,ierr)
            end if
         end do mesh_blocks
      else master_check
         call mpi_recv(iblbuf,bufsize,INT_TYPE,(this_mesh%master_id-1),tag,&
               DEFAULT_COMM,stats_mpi,ierr)
         i=0
         do l=1,lmax
            do k=1,kmax
               do j=1,jmax
                  i=i+1
                  iblank(j,k,l)=iblbuf(i)
               end do
            end do
         end do
      end if master_check
      call barrier_mpi
      deallocate(iblbuf)
   contains
      subroutine determine_send_limits
         !! Determine the lower and upper j,k,l limits of each block
         !! taking ghost points into account. This information isn't
         !! stored in the master processor, so we need to calculate
         !! them.
         integer, dimension(6) :: extents
         integer, dimension(3) :: maxes

         extents=(/ bst(JDIR,n), bst(KDIR,n), bst(LDIR,n), &
               bend(JDIR,n), bend(KDIR,n), bend(LDIR,n) /)
         maxes=(/ this_mesh%jm, this_mesh%km, this_mesh%lm /)

         do i=1,3
            if (extents(i) == 1) then 
               blk_lim(i) = extents(i)
            else
               blk_lim(i) = extents(i)-NGHOST
            end if
         end do
         
         ! The max faces
         do i=4,6
            if (extents(i) == maxes(i-3)) then 
               blk_lim(i) = extents(i)
            else
               blk_lim(i) = extents(i)+NGHOST
            end if
         end do
      end subroutine determine_send_limits
   end subroutine update_ibc_block_iblanks

END MODULE IBC_MODULE



!                        write (filename, 10) PID
!10                      format ('fus_box_ibc_', I0, '.p3d')
!		        open(unit=30,file=filename,form='unformatted',status='unknown')
!		        write(30) imax_mesh,jmax_mesh,kmax_mesh
!			write(30) (((xg(i,j,k),i=1,imax_mesh),j=1,jmax_mesh),k=1,kmax_mesh),&
!				  (((yg(i,j,k),i=1,imax_mesh),j=1,jmax_mesh),k=1,kmax_mesh),&
!				  (((zg(i,j,k),i=1,imax_mesh),j=1,jmax_mesh),k=1,kmax_mesh),&
!				  (((iblk(i,j,k),i=1,imax_mesh),j=1,jmax_mesh),k=1,kmax_mesh)
!		        close(30)


!                        write (filename, 11) PID
!11                      format ('bound_nodes_ibc_', I0, '.p3d')
!		        open(unit=30,file=filename,form='unformatted',status='unknown')
!		        write(30) nbound,1,1
!			write(30) (((xg(boundNode(i)%i,boundNode(i)%j,boundNode(i)%k),i=1,nbound),j=1,1),k=1,1),&
!			          (((yg(boundNode(i)%i,boundNode(i)%j,boundNode(i)%k),i=1,nbound),j=1,1),k=1,1),&
!			          (((zg(boundNode(i)%i,boundNode(i)%j,boundNode(i)%k),i=1,nbound),j=1,1),k=1,1),&
!			          (((1,i=1,nbound),j=1,1),k=1,1)
!		        close(30)


!                        ! write normal intersection points
!                        write (filename, 13) PID
!13                      format ('surf_nodes_ibc_', I0, '.p3d')
!		        open(unit=30,file=filename,form='unformatted',status='unknown')
!		        write(30) nbound,1,1
!			write(30) (boundNode(i)%surf(1),i=1,nbound),&
!			          (boundNode(i)%surf(2),i=1,nbound),&
!			          (boundNode(i)%surf(3),i=1,nbound),&
!			          (1,i=1,nbound)
!		        close(30)

