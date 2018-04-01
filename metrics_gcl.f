

c***********************************************************************
      subroutine metgcl(ntac_loc,kkp,kkr,x,y,z,xt2,yt2,zt2,svj,svk,svl,
     <                  xx,xy,xz,yx,yy,yz,zx,zy,zz,
     <                  xaj,yaj,zaj,xak,yak,zak,xal,yal,zal,
     <                  q,a3,vnaj,vnak,vnal,iprint)

c  Finite volume formulation  consistent with gcl
c
c  Compute the metrics and the jacobian for computational space
c  uniform computational space, deltas = 1 < averaged metrics >
c
c  Also compute finite-volume quantities for right-hand-side:
c  cell face normal velocities (FV analog of time metrics) and 
c  cell face area vectors (FV analog of grad(xi) * inv(J)). Note,
c  that these quantities correspond to cell faces (e.g j => j+1/2),
c  whereas the "metrics" correspond to cell centers. Hence, to avoid
c  confusion, the data corresponding to cell faces are not named 
c  according to their FD metric analogs.
c
c  Historically, FV cell data has been computed on the fly using stored
c  doubly refined grids at n-1, n, and n+1 time levels. However, this 
c  results in a significant increase in the amount of memory consumed
c  and redundant computation of cell geometry (e.g. cell face normals,
c  and swept volumes) every sub-iteration of the solver. Thus, the
c  doubly refined meshes are now only used as temporary/intermediate 
c  data for generating the FD data required by the RHS procedures.
c
c  INPUT:
c  ntac_loc...... local time accuracy (either 1 or 2)
c  kkp,kkr ...... k-plane pointer arrays
c  x,y,z......... Cartesian coordinates of deformed grid at n+1
c  xt2,yt2,zt2... Cartesian coordinates of deformed grid at n
c  iprint........... Option flag for writting doubly refined mesh to file;
c                    if iprint = 1, doubly refined mesh is written to 
c                    plot3d formated grid file. Used mostly for debugging.
c
c  IN/OUT:
c  svj,svk,svl... Cell face swept volumes corresponding to [n-1,n] interval
c                 svj: j+1/2 faces; svk: k+1/2 faces; svl: l+1/2 faces
c
c  OUTPUT:
c  xx,xy,xz...... xi gradient components corresponding to grid at n+1
c  yx,yy,yz...... eta gradient components corresponding to grid at n+1
c  zx,zy,zz...... zeta gradient components corresponding to grid at n+1
c  q............. Solution vector at n+1 with updated Jacobian (see note 1)
c  a3............ cell volumes corresponding to grid at n+1 (see note 2)
c  xaj,yaj,zaj... cell face area vector Cartesian components, j+1/2 faces
c  xak,yak,zak... cell face area vector Cartesian components, k+1/2 faces
c  xal,yal,zal... cell face area vector Cartesian components, l+1/2 faces
c  vnaj,vnak,vnal... mean cell face velocity dotted with cell face area
c
c  NOTES:
c  1.   Only the Jacobian is updated in the solution vector, q; i.e.
c       q(:,:,:,6) vector. Flow variables are not modified. This is
c       a bit of a kludge, but it is consistent with the metfv 
c       subroutine behavior; for the time being, will keep things 
c       consistent.
c  2.   q(:,:,:,6) and a3 contain redundant data: q(:,:,:,6) = 1.0/a3.
c       Can we get rid of one of these?
c
c***********************************************************************

      use params_global

      implicit none
      
      integer iprint
      integer ntac_loc
      
      real q(jmax,kmax,lmax,nd)
      real x(jmax,kmax,lmax),y(jmax,kmax,lmax),z(jmax,kmax,lmax)
      real xt2(jmax,kmax,lmax),yt2(jmax,kmax,lmax),zt2(jmax,kmax,lmax)

      real a3(jmax,kmax,lmax)

      real xx(jmax,kmax,lmax),xy(jmax,kmax,lmax),xz(jmax,kmax,lmax),
     <     yx(jmax,kmax,lmax),yy(jmax,kmax,lmax),yz(jmax,kmax,lmax),
     <     zx(jmax,kmax,lmax),zy(jmax,kmax,lmax),zz(jmax,kmax,lmax)

      real svj(jmax,kmax,lmax), svk(jmax,kmax,lmax), svl(jmax,kmax,lmax)
      real vnaj(jmax,kmax,lmax), vnak(jmax,kmax,lmax), vnal(jmax,kmax,lmax)
      real xaj(jmax,kmax,lmax), yaj(jmax,kmax,lmax), zaj(jmax,kmax,lmax)
      real xak(jmax,kmax,lmax), yak(jmax,kmax,lmax), zak(jmax,kmax,lmax)
      real xal(jmax,kmax,lmax), yal(jmax,kmax,lmax), zal(jmax,kmax,lmax)

      integer kkp(kmax),kkr(kmax)

c..   local variables

      logical :: use_refmesh_metrics
      integer :: jbig, kbig, lbig
      integer :: j, k, l
      real, allocatable, dimension(:,:,:) :: xbig, ybig, zbig
      real, allocatable, dimension(:,:,:) :: xold, yold, zold

c..   Will only use doubly refined mesh for FV cell data; space metrics
c..   will be constructed from baseline mesh. Use of doubly refined mesh
c..   metrics was found to cause slow convergence of the Newton 
c..   iterations with LU-SGS factorization. In principle, solver accuracy 
c..   is dictated by the RHS. So, provided that the RHS is treated in 
c..   a manner consistent with GCL, the metrics used in the LHS can be
c..   formulated such that fast and stable convergence is realized.

      use_refmesh_metrics = .false.

c..   Set up temporary data structure

      jbig = 2*jmax-1
      kbig = 2*kmax-1
      lbig = 2*lmax-1

      allocate(xbig(jbig,kbig,lbig), 
     <         ybig(jbig,kbig,lbig), 
     <         zbig(jbig,kbig,lbig))
      allocate(xold(jbig,kbig,lbig), 
     <         yold(jbig,kbig,lbig), 
     <         zold(jbig,kbig,lbig))


c..   Construct doubly refined meshes

      call build_refmesh(xt2,yt2,zt2,xold,yold,zold,jmax,kmax,lmax)
      call build_refmesh(x,y,z,xbig,ybig,zbig,jmax,kmax,lmax)

      if (iprint.eq.1) then

         write(2) 2*jmax-1,2*kmax-1,2*lmax-1
         write(2)(((xbig(j,k,l),j=1,2*jmax-1),k=1,2*kmax-1),
     $        l=1,2*lmax-1),
     1        (((ybig(j,k,l),j=1,2*jmax-1),k=1,2*kmax-1),l=1,2*lmax-1),
     2        (((zbig(j,k,l),j=1,2*jmax-1),k=1,2*kmax-1),l=1,2*lmax-1)
      endif

      if(.not.use_refmesh_metrics)then
c..     Compute metrics using baseline mesh
        call metfv(q,x,y,z,xx,xy,xz,yx,yy,yz,zx,zy,zz,a3,kkp,kkr,iprint)
      end if

c.. *** COMPUTE FINITE VOLUME DATA FOR RHS ***

c.. Compute the cell volumes following a dual refined cell approach
c.. ( see zhang et al. computers fluids vol.22, no.1 1993 ) 

      call calc_cell_volumes(xbig,ybig,zbig,a3,jmax,kmax,lmax)

c.. Compute cell face areas and associated swept volumes
c.. This will ensure that the RHS satisfies the GCL

      call calc_cell_face(dt,xbig,ybig,zbig,xold,yold,zold,
     <     svj,svk,svl,xaj,yaj,zaj,xak,yak,zak,xal,yal,zal,
     <     vnaj,vnak,vnal,jmax,kmax,lmax,ntac_loc)

c.. *** COMPUTE FINITE DIFFERENCE METRICS FOR LHS ***

c.. Set Jacobians in solution vector
c.. Here we set Jacobians to inverse of cell volume to satsify GCL

      call calc_jacobians(a3,q,jmax,kmax,lmax,nd)

c.. let's actually use the refined mesh in finite volume like way
c   rather than use slopes form metslope

      if(use_refmesh_metrics)then
        call compute_metrics(q,xbig,ybig,zbig,xx,xy,xz,yx,yy,yz,
     <       zx,zy,zz,jmax,kmax,lmax,nd)
      end if

c..   clean up memeory
      deallocate(xbig,ybig,zbig)
      deallocate(xold,yold,zold)

      end

c***********************************************************************
      subroutine build_refmesh(x,y,z,xbig,ybig,zbig,jmax,kmax,lmax)
c***********************************************************************

      implicit none

      integer, intent(in) :: jmax, kmax, lmax
      real, dimension(jmax,kmax,lmax), intent(in) ::  x, y, z
      real, dimension(2*jmax-1, 2*kmax-1, 2*lmax-1), intent(out) ::
     <  xbig, ybig, zbig

      !local var
      integer :: j, k, l
      integer :: k21, l21
      
c..let's define the refined mesh for metric calculations

      do l=1,lmax
         do k=1,kmax
            do j=1,jmax
               xbig(2*j-1,2*k-1,2*l-1)=x(j,k,l)
               ybig(2*j-1,2*k-1,2*l-1)=y(j,k,l)
               zbig(2*j-1,2*k-1,2*l-1)=z(j,k,l)
            enddo
         enddo
      enddo

c..   interpolating using metslope to find midpoints (slope not used)
c..   first find all points on each k-plane

      do k=1,kmax
         k21=2*k-1
         
         do l=1,lmax
            l21=2*l-1
            do j=1,jmax-1
               xbig(2*j,k21,l21)=(x(j,k,l)+x(j+1,k,l))*0.5
               ybig(2*j,k21,l21)=(y(j,k,l)+y(j+1,k,l))*0.5
               zbig(2*j,k21,l21)=(z(j,k,l)+z(j+1,k,l))*0.5
            enddo
         enddo

         do j=1,2*jmax-1
            do l=1,lmax-1
               xbig(j,k21,2*l)=(xbig(j,k21,2*l-1)+xbig(j,k21,2*l+1))*0.5
               ybig(j,k21,2*l)=(ybig(j,k21,2*l-1)+ybig(j,k21,2*l+1))*0.5
               zbig(j,k21,2*l)=(zbig(j,k21,2*l-1)+zbig(j,k21,2*l+1))*0.5
            enddo
         enddo
      enddo

      do j=1,2*jmax-1
         do l=1,2*lmax-1
            do k=1,kmax-1
               xbig(j,2*k,l) = (xbig(j,2*k-1,l)+xbig(j,2*k+1,l))*0.5
               ybig(j,2*k,l) = (ybig(j,2*k-1,l)+ybig(j,2*k+1,l))*0.5
               zbig(j,2*k,l) = (zbig(j,2*k-1,l)+zbig(j,2*k+1,l))*0.5
            enddo
         enddo
      enddo 
      end


c***********************************************************************
      subroutine calc_cell_volumes(xbig,ybig,zbig,a3,jmax,kmax,lmax)

c.. Compute the cell volumes following a dual refined cell approach
c.. ( see zhang et al. computers fluids vol.22, no.1 1993 ) 

c.. each grid point uses all the hexahedra in the refined mesh which
c.. contains it, so essentially for a inside node there are volume of
c.. 8 hexahedra to be computed

c***********************************************************************

      implicit none

      integer, intent(in) :: jmax, kmax, lmax
      real, dimension(2*jmax-1,2*kmax-1,2*lmax-1), intent(in) :: 
     <  xbig,ybig,zbig
      real, dimension(jmax,kmax,lmax), intent(out) :: a3

      !local var
      integer :: i
      integer :: nneg
      integer :: j, k, l
      integer :: jx,kx,lx
      integer :: k2, j2, l2
      integer :: k21, j21, l21

      integer jvert(8),kvert(8),lvert(8),index(3)

      real xv1,yv1,zv1
      real xvert(8),yvert(8),zvert(8)
      real volume
      
c..   initializing to zero for peace of mind
      
      do l=1,lmax
         do k=1,kmax
            do j=1,jmax
               a3(j,k,l) = 0.0
            enddo
         enddo
      enddo
  
c..   ubr piece

      do lx=0,1
         do kx=0,1
            do jx=0,1
               do l=2-lx,lmax-lx
                  l21=2*l-2
                  l2=2*l-1
                  
                  do k=2-kx,kmax-kx
                     k21=2*k-1
                     k2=2*k-2
                     
                     do j=2-jx,jmax-jx
                        j21=2*j-1
                        j2=2*j-2
                        
                        jvert(1)=j21; kvert(1)=k2;    lvert(1)=l21;
                        jvert(2)=j21; kvert(2)=k21;   lvert(2)=l21;  
                        jvert(3)=j2;  kvert(3)=k21;   lvert(3)=l21;  
                        jvert(4)=j2;  kvert(4)=k2;    lvert(4)=l21;
                        jvert(5)=j21; kvert(5)=k2;     lvert(5)=l2;
                        jvert(6)=j21; kvert(6)=k21;    lvert(6)=l2;  
                        jvert(7)=j2;  kvert(7)=k21;    lvert(7)=l2;  
                        jvert(8)=j2;  kvert(8)=k2;     lvert(8)=l2;
                        
                        do i=1,8
                           jvert(i)=jvert(i)+jx
                           kvert(i)=kvert(i)+kx
                           lvert(i)=lvert(i)+lx
                        enddo

                        index(1)=jvert(4)
                        index(2)=kvert(4)
                        index(3)=lvert(4)
                        
                        i=1
                        xv1=xbig(jvert(i),kvert(i),lvert(i))
                        yv1=ybig(jvert(i),kvert(i),lvert(i))
                        zv1=zbig(jvert(i),kvert(i),lvert(i))

                        xvert(i)=0.
                        yvert(i)=0.
                        zvert(i)=0.
                        
                        do i=2,8
                           xvert(i)=xbig(jvert(i),kvert(i),lvert(i))
     $                          -xv1
                           yvert(i)=ybig(jvert(i),kvert(i),lvert(i))
     $                          -yv1
                           zvert(i)=zbig(jvert(i),kvert(i),lvert(i))
     $                          -zv1
                        enddo
                        
                        a3(j,k,l)=a3(j,k,l)+volume(index,xvert,yvert,
     $                       zvert)
                        
                     enddo
                  enddo
               enddo
               
            enddo
         enddo
      enddo

c..check for negative jacobians

      nneg=0
      do 75 l = 1,lmax
         do 75 k = 1,kmax
            do 75 j = 1,jmax
           if( a3(j,k,l).le.0.0 ) then
            write(6,603) a3(j,k,l), j, k, l
            nneg=nneg+1
          end if
 75    continue

        if (nneg.gt.0) then
        print *,'nneg=',nneg
        stop
        endif

  603 format( ' ',10x,'negative jacobian = ',1p,e10.3,1x,'at j,k,l =',
     $                 3i5,5x)

      end

c***********************************************************************
      subroutine calc_jacobians(a3,q,jmax,kmax,lmax,nd)
c***********************************************************************

      implicit none

      integer, intent(in) :: jmax, kmax, lmax, nd
      real, dimension(jmax,kmax,lmax), intent(in) :: a3
      real, dimension(jmax,kmax,lmax,nd), intent(inout) :: q

      !local
      integer :: j, k, l

c.. setting the jacobians

      do j=1,jmax
         do k=1,kmax
            do l=1,lmax
               q(j,k,l,6)=1.0/a3(j,k,l)
            enddo
         enddo
      enddo

c.. setting jacobians on the boundary planes right

      do l=1,lmax,lmax-1
         do k=1,kmax
            do j=1,jmax
               q(j,k,l,6)=0.5*q(j,k,l,6)
            enddo
         enddo
      enddo
                  
      do j=1,jmax,jmax-1
         do l=1,lmax
            do k=1,kmax
               q(j,k,l,6)=0.5*q(j,k,l,6)
            enddo
         enddo
      enddo

      do k=1,kmax,kmax-1
         do j=1,jmax
            do l=1,lmax
               q(j,k,l,6)=0.5*q(j,k,l,6)
            enddo
         enddo
      enddo

      end 


c***********************************************************************
      subroutine compute_metrics(q,xbig,ybig,zbig,xx,xy,xz,yx,yy,yz,
     <           zx,zy,zz,jmax,kmax,lmax,nd)
c***********************************************************************

      implicit none

      integer, intent(in) :: jmax, kmax, lmax, nd
      real, dimension(2*jmax-1,2*kmax-1,2*lmax-1), intent(in) ::
     <  xbig, ybig, zbig
      real, dimension(jmax,kmax,lmax,nd), intent(in) :: q
      real, dimension(jmax,kmax,lmax), intent(out) :: xx, xy, xz,
     <                                                yx, yy, yz,
     <                                                zx, zy, zz

      ! local var
      integer :: j, k, l
      integer :: j21, k21, l21
      integer :: jp1, kp1, lp1
      integer :: jm1, km1, lm1
      integer :: jm, km, lm
      real :: qj
      real :: fac1,fac2
      real :: dx2,dy2,dz2,dx0,dy0,dz0

      integer :: kkp1(kmax),kkr1(kmax)
      integer :: jjp(jmax),jjr(jmax),llp(lmax),llr(lmax)

      ! set up temporary data structure

      jm      = jmax -1
      km      = kmax -1
      lm      = lmax -1

      fac1    = 0.5
      fac2    = 1./3.

      do 1  j = 1,jmax
       jjp(j) = 2*j
       jjr(j) = 2*j - 2
    1 continue
      jjp(jmax)=2*jmax-1
      jjr(1   )=1

      do 2  l = 1,lmax
         llp(l) = 2*l
         llr(l) = 2*l - 2
    2 continue
      llp(lmax)=2*lmax-1
      llr(1   )=1
      
      do 3  k = 1,kmax
         kkp1(k) = 2*k
         kkr1(k) = 2*k - 2
 3    continue
      kkp1(kmax)=2*kmax - 1
      kkr1(1   )=1

c.. Compute space metrics

c..xi derivatives

      do l=1,lmax

         lp1=llp(l)
         lm1=llr(l)

         do k=1,kmax

            kp1=kkp1(k)
            km1=kkr1(k)

            do j=1,jmax
               j21=2*j-1
               dx0=xbig(j21,kp1,lp1)-xbig(j21,km1,lm1)
               dy0=ybig(j21,kp1,lp1)-ybig(j21,km1,lm1)
               dz0=zbig(j21,kp1,lp1)-zbig(j21,km1,lm1)

               dx2=xbig(j21,km1,lp1)-xbig(j21,kp1,lm1)
               dy2=ybig(j21,km1,lp1)-ybig(j21,kp1,lm1)
               dz2=zbig(j21,km1,lp1)-zbig(j21,kp1,lm1)

               xx(j,k,l)=fac1*( dy0*dz2 - dy2*dz0 )
               xy(j,k,l)=fac1*( dz0*dx2 - dz2*dx0 )
               xz(j,k,l)=fac1*( dx0*dy2 - dx2*dy0 )
            enddo
         enddo
      enddo
            
      do j=1,jmax
         do k=1,kmax,kmax-1
            do l=1,lmax
               xx(j,k,l)=xx(j,k,l)*2
               xy(j,k,l)=xy(j,k,l)*2
               xz(j,k,l)=xz(j,k,l)*2
            enddo
         enddo
         
         do l=1,lmax,lmax-1
            do k=1,kmax
               xx(j,k,l)=xx(j,k,l)*2
               xy(j,k,l)=xy(j,k,l)*2
               xz(j,k,l)=xz(j,k,l)*2
            enddo
         enddo
      enddo

c..eta derivatives

      do l=1,lmax

         lp1=llp(l)
         lm1=llr(l)

         do j=1,jmax

            jp1=jjp(j)
            jm1=jjr(j)

            do k=1,kmax
               k21=2*k-1
               dx2=xbig(jp1,k21,lp1)-xbig(jm1,k21,lm1)
               dy2=ybig(jp1,k21,lp1)-ybig(jm1,k21,lm1)
               dz2=zbig(jp1,k21,lp1)-zbig(jm1,k21,lm1)

               dx0=xbig(jm1,k21,lp1)-xbig(jp1,k21,lm1)
               dy0=ybig(jm1,k21,lp1)-ybig(jp1,k21,lm1)
               dz0=zbig(jm1,k21,lp1)-zbig(jp1,k21,lm1)

               yx(j,k,l)=fac1*( dy0*dz2 - dy2*dz0 )
               yy(j,k,l)=fac1*( dz0*dx2 - dz2*dx0 )
               yz(j,k,l)=fac1*( dx0*dy2 - dx2*dy0 )
            enddo
         enddo
      enddo
            
      do k=1,kmax
         do l=1,lmax,lmax-1
            do j=1,jmax
               yx(j,k,l)=yx(j,k,l)*2
               yy(j,k,l)=yy(j,k,l)*2
               yz(j,k,l)=yz(j,k,l)*2
            enddo
         enddo
         
         do j=1,jmax,jmax-1
            do l=1,lmax
               yx(j,k,l)=yx(j,k,l)*2
               yy(j,k,l)=yy(j,k,l)*2
               yz(j,k,l)=yz(j,k,l)*2
            enddo
         enddo
      enddo

c..   zeta derivative

      do k=1,kmax
         
         kp1=kkp1(k)
         km1=kkr1(k)
         
         do j=1,jmax

            jp1=jjp(j)
            jm1=jjr(j)

            do l=1,lmax
               l21=2*l-1
               dx0=xbig(jp1,kp1,l21)-xbig(jm1,km1,l21)
               dy0=ybig(jp1,kp1,l21)-ybig(jm1,km1,l21)
               dz0=zbig(jp1,kp1,l21)-zbig(jm1,km1,l21)

               dx2=xbig(jm1,kp1,l21)-xbig(jp1,km1,l21)
               dy2=ybig(jm1,kp1,l21)-ybig(jp1,km1,l21)
               dz2=zbig(jm1,kp1,l21)-zbig(jp1,km1,l21)

               zx(j,k,l)=fac1*( dy0*dz2 - dy2*dz0 )
               zy(j,k,l)=fac1*( dz0*dx2 - dz2*dx0 )
               zz(j,k,l)=fac1*( dx0*dy2 - dx2*dy0 )
            enddo
         enddo
      enddo
            
      do l=1,lmax
         do k=1,kmax,kmax-1
            do j=1,jmax
               zx(j,k,l)=zx(j,k,l)*2
               zy(j,k,l)=zy(j,k,l)*2
               zz(j,k,l)=zz(j,k,l)*2
            enddo
         enddo
         
         do j=1,jmax,jmax-1
            do k=1,kmax
               zx(j,k,l)=zx(j,k,l)*2
               zy(j,k,l)=zy(j,k,l)*2
               zz(j,k,l)=zz(j,k,l)*2
            enddo
         enddo
      enddo

c..   scaling the metrics by the jacobian

      do j=1,jmax
         do k=1,kmax
            do l=1,lmax
               qj=q(j,k,l,6)
               
               xx(j,k,l)=xx(j,k,l)*qj
               xy(j,k,l)=xy(j,k,l)*qj
               xz(j,k,l)=xz(j,k,l)*qj

               yx(j,k,l)=yx(j,k,l)*qj
               yy(j,k,l)=yy(j,k,l)*qj
               yz(j,k,l)=yz(j,k,l)*qj

               zx(j,k,l)=zx(j,k,l)*qj
               zy(j,k,l)=zy(j,k,l)*qj
               zz(j,k,l)=zz(j,k,l)*qj
               
            enddo
         enddo
      enddo

      end 
  
c***********************************************************************
      subroutine calc_cell_face(dt,xbig,ybig,zbig,xold,yold,zold,
     <           svj,svk,svl,xaj,yaj,zaj,xak,yak,zak,xal,yal,zal,
     <           vnaj,vnak,vnal,jmax,kmax,lmax,ntac_loc)

c  (May want to decompose this into two subroutines: one for face areas
c   and one for swept volumes?)

c***********************************************************************

      implicit none

      integer :: jmax, kmax, lmax, ntac_loc
      real :: dt
      real, dimension(2*jmax-1, 2*kmax-1, 2*lmax-1), intent(in) ::
     < xbig, ybig, zbig, xold, yold, zold
      real, dimension(jmax, kmax, lmax), intent(inout) :: svj, svk, svl
      real, dimension(jmax, kmax, lmax), intent(out) :: xaj, yaj, zaj
      real, dimension(jmax, kmax, lmax), intent(out) :: xak, yak, zak
      real, dimension(jmax, kmax, lmax), intent(out) :: xal, yal, zal
      real, dimension(jmax, kmax, lmax), intent(out) :: vnaj, vnak, vnal

c..   local variables ! arrays


C..   local variables 

      integer :: i
      integer :: j, k, l
      integer :: jm, km, lm
      integer :: jm2, km2, lm2
      integer :: jp2, kp2, lp2
      integer :: j21, k21, l21
      integer :: jx, kx, lx

      integer jvert(8),kvert(8),lvert(8),index(3)
      
      real :: fac1, temp
      real :: dx0,dy0,dz0,dx2,dy2,dz2
      real :: volume

      real xvert(8),yvert(8),zvert(8)

c************************************************************************

      jm = jmax - 1
      km = kmax - 1
      lm = lmax - 1

      fac1=0.5

c..xi faces 

      do 13 l = 2, lm
      do 13 k = 2, km

        kp2=2*k
        km2=2*k-2
        lp2=2*l
        lm2=2*l-2
        k21=2*k-1
        l21=2*l-1

        do j=1,jm

c..        Initialize jth cell face area and swept volume to zero

           xaj(j,k,l) = 0.0
           yaj(j,k,l) = 0.0
           zaj(j,k,l) = 0.0

           temp=0.0

c..        Compute each of the four facets on the xi-face and the volume
c..        swept by each facet over the interval [t_n, t_{n+1}]

           do kx=0,1
              do lx=0,1

c..              Set vertex indicies of a quadralaterial contained 
c..              in the jth face

                 jvert(1)=2*j ; kvert(1)=km2+kx; lvert(1)=lm2+lx
                 jvert(2)=2*j ; kvert(2)=k21+kx; lvert(2)=lm2+lx
                 jvert(3)=2*j ; kvert(3)=k21+kx; lvert(3)=l21+lx
                 jvert(4)=2*j ; kvert(4)=km2+kx; lvert(4)=l21+lx

c..              *** Compute contribution to cell face surface integral ***

                 dx0=xbig(jvert(3),kvert(3),lvert(3))
     <              -xbig(jvert(1),kvert(1),lvert(1))
                 dy0=ybig(jvert(3),kvert(3),lvert(3))
     <              -ybig(jvert(1),kvert(1),lvert(1))
                 dz0=zbig(jvert(3),kvert(3),lvert(3))
     <              -zbig(jvert(1),kvert(1),lvert(1))
                 
                 dx2=xbig(jvert(4),kvert(4),lvert(4))
     <              -xbig(jvert(2),kvert(2),lvert(2))
                 dy2=ybig(jvert(4),kvert(4),lvert(4))
     <              -ybig(jvert(2),kvert(2),lvert(2))
                 dz2=zbig(jvert(4),kvert(4),lvert(4))
     <              -zbig(jvert(2),kvert(2),lvert(2))
                 
                 xaj(j,k,l)= xaj(j,k,l) + fac1*( dy0*dz2 - dy2*dz0 )
                 yaj(j,k,l)= yaj(j,k,l) + fac1*( dz0*dx2 - dz2*dx0 )
                 zaj(j,k,l)= zaj(j,k,l) + fac1*( dx0*dy2 - dx2*dy0 )

c..              *** Compute contribution to swept volume ***

c..              Set the corresponding cartesian coordinates
c..              of the verticies at n and n+1.

                 do i=1,4

                    xvert(i)=xold(jvert(i),kvert(i),lvert(i))
                    yvert(i)=yold(jvert(i),kvert(i),lvert(i))
                    zvert(i)=zold(jvert(i),kvert(i),lvert(i))

                    xvert(i+4)=xbig(jvert(i),kvert(i),lvert(i))
                    yvert(i+4)=ybig(jvert(i),kvert(i),lvert(i))
                    zvert(i+4)=zbig(jvert(i),kvert(i),lvert(i))

                 enddo
           
c..              Note that 'index' array is not used in volume -- must be
c..              an artifact from an old implementation (09/21/2008).

                 index(1)=jvert(1)
                 index(2)=kvert(1)
                 index(3)=lvert(1)
                 
c..              Add swept volume contribution -- note that storing the
c..              verticies of xold in xvert(1:4), yvert(1:4), and zvert(1:4)
c..              and verticies of xbig in xvert(5:8), yvert(5:8), and zvert(5:8)
c..              will result in a negative volume for positive swept volumes.
c..              Note that this sign convention corresponds to the 
c..              sign associated with the 'time-metric' analog at j+1/2.

                 temp=temp-volume(index,xvert,yvert,zvert)

              enddo
           enddo

           ! Compute negative of mean face velocity dotted with ndA
           if(ntac_loc.eq.2)then 
                vnaj(j,k,l) = 1.5*temp/dt - 0.5*svj(j,k,l)/dt
           else
                vnaj(j,k,l)=temp/dt
           end if

           ! store swept volume for next iter
           svj(j,k,l) = temp

        enddo ! loop over j
   13 continue ! loop over k and l

c
c..eta fluxes
c
      if(kmax.gt.3) then
c
      do 23 j = 2,jm
      do 23 l = 2,lm

        jp2=2*j
        jm2=2*j-2
        lp2=2*l
        lm2=2*l-2
        j21=2*j-1
        l21=2*l-1

        do k=1,km

c..        Initialize kth cell face area and swept volume to zero

           xak(j,k,l) = 0.0
           yak(j,k,l) = 0.0
           zak(j,k,l) = 0.0
           
c..        Compute each of the four facets on the xi-face and the volume
c..        swept by each facet over the interval [t_n, t_{n+1}]

           temp=0.0

           do jx=0,1
              do lx=0,1

c..              Set vertex indicies of a quadralaterial contained 
c..              in the jth face

                 jvert(1)=jp2-jx ; kvert(1)=k*2; lvert(1)=l21+lx
                 jvert(2)=jp2-jx ; kvert(2)=k*2; lvert(2)=lm2+lx
                 jvert(3)=j21-jx ; kvert(3)=k*2; lvert(3)=lm2+lx
                 jvert(4)=j21-jx ; kvert(4)=k*2; lvert(4)=l21+lx

c..              *** Compute contribution to cell face surface integral ***

                 dx0=xbig(jvert(4),kvert(4),lvert(4))
     <              -xbig(jvert(2),kvert(2),lvert(2))
                 dy0=ybig(jvert(4),kvert(4),lvert(4))
     <              -ybig(jvert(2),kvert(2),lvert(2))
                 dz0=zbig(jvert(4),kvert(4),lvert(4))
     <              -zbig(jvert(2),kvert(2),lvert(2))

                 dx2=xbig(jvert(1),kvert(1),lvert(1))
     <              -xbig(jvert(3),kvert(3),lvert(3))
                 dy2=ybig(jvert(1),kvert(1),lvert(1))
     <              -ybig(jvert(3),kvert(3),lvert(3))
                 dz2=zbig(jvert(1),kvert(1),lvert(1))
     <              -zbig(jvert(3),kvert(3),lvert(3))
                                 
                 xak(j,k,l)= xak(j,k,l) + fac1*( dy0*dz2 - dy2*dz0 )
                 yak(j,k,l)= yak(j,k,l) + fac1*( dz0*dx2 - dz2*dx0 )
                 zak(j,k,l)= zak(j,k,l) + fac1*( dx0*dy2 - dx2*dy0 )
                 
c..              *** Compute contribution to swept volume ***

                 do i=1,4

                    xvert(i)=xold(jvert(i),kvert(i),lvert(i))
                    yvert(i)=yold(jvert(i),kvert(i),lvert(i))
                    zvert(i)=zold(jvert(i),kvert(i),lvert(i))

                    xvert(i+4)=xbig(jvert(i),kvert(i),lvert(i))
                    yvert(i+4)=ybig(jvert(i),kvert(i),lvert(i))
                    zvert(i+4)=zbig(jvert(i),kvert(i),lvert(i))
                 enddo

                 index(1)=jvert(2)
                 index(2)=kvert(2)
                 index(3)=lvert(2)
                 
                 temp=temp-volume(index,xvert,yvert,zvert)
              enddo
           enddo

           ! Compute negative of mean face velocity dotted with ndA
           if(ntac_loc.eq.2)then 
                vnak(j,k,l) = 1.5*temp/dt - 0.5*svk(j,k,l)/dt
           else
                vnak(j,k,l)=temp/dt
           end if

           ! store swept volume for next iter
           svk(j,k,l) = temp

        end do ! loop over k

   23 continue ! loop over j and l
c     
      endif

c
c..zeta fluxes
c
      do 33 j = 2,jm
      do 33 k = 2,km
c
        jp2=2*j
        jm2=2*j-2
        kp2=2*k
        km2=2*k-2
        j21=2*j-1
        k21=2*k-1

        do l=1,lm

c..        Initialize jth cell face area and swept volume to zero

           xal(j,k,l) = 0.0
           yal(j,k,l) = 0.0
           zal(j,k,l) = 0.0

           temp=0.0

c..        Compute each of the four facets on the xi-face and the volume
c..        swept by each facet over the interval [t_n, t_{n+1}]

           do jx=0,1
              do kx=0,1

c..              Set vertex indicies of a quadralaterial contained 
c..              in the lth face

                 jvert(1)=jp2-jx ; kvert(1)=k21-kx; lvert(1)=l*2
                 jvert(2)=jp2-jx ; kvert(2)=kp2-kx; lvert(2)=l*2
                 jvert(3)=j21-jx ; kvert(3)=kp2-kx; lvert(3)=l*2
                 jvert(4)=j21-jx ; kvert(4)=k21-kx; lvert(4)=l*2

c..              *** Compute contribution to cell face surface integral ***

                 dx0=xbig(jvert(3),kvert(3),lvert(3))
     <              -xbig(jvert(1),kvert(1),lvert(1))
                 dy0=ybig(jvert(3),kvert(3),lvert(3))
     <              -ybig(jvert(1),kvert(1),lvert(1))
                 dz0=zbig(jvert(3),kvert(3),lvert(3))
     <              -zbig(jvert(1),kvert(1),lvert(1))
                 
                 dx2=xbig(jvert(4),kvert(4),lvert(4))
     <              -xbig(jvert(2),kvert(2),lvert(2))
                 dy2=ybig(jvert(4),kvert(4),lvert(4))
     <              -ybig(jvert(2),kvert(2),lvert(2))
                 dz2=zbig(jvert(4),kvert(4),lvert(4))
     <              -zbig(jvert(2),kvert(2),lvert(2))
                 
                 xal(j,k,l)= xal(j,k,l) + fac1*( dy0*dz2 - dy2*dz0 )
                 yal(j,k,l)= yal(j,k,l) + fac1*( dz0*dx2 - dz2*dx0 )
                 zal(j,k,l)= zal(j,k,l) + fac1*( dx0*dy2 - dx2*dy0 )

c..              *** Compute contribution to swept volume ***
                 
                 do i=1,4

                    xvert(i)=xold(jvert(i),kvert(i),lvert(i))
                    yvert(i)=yold(jvert(i),kvert(i),lvert(i))
                    zvert(i)=zold(jvert(i),kvert(i),lvert(i))

                    xvert(i+4)=xbig(jvert(i),kvert(i),lvert(i))
                    yvert(i+4)=ybig(jvert(i),kvert(i),lvert(i))
                    zvert(i+4)=zbig(jvert(i),kvert(i),lvert(i))
                 enddo
           
                 index(1)=jvert(4)
                 index(2)=kvert(4)
                 index(3)=lvert(4)
                 
                 temp=temp-volume(index,xvert,yvert,zvert)
              enddo
           enddo

           ! Compute negative of mean face velocity dotted with ndA
           if(ntac_loc.eq.2)then 
                vnal(j,k,l) = 1.5*temp/dt - 0.5*svl(j,k,l)/dt
           else
                vnal(j,k,l)=temp/dt
           end if

           ! store swept volume for next iter
           svl(j,k,l) = temp

        end do ! loop over l
   33 continue ! loop over j and k

      end

