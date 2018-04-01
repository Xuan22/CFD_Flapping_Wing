!********************************************************************
module octree
  implicit none
!********************************************************************

  integer :: lvlmin,lvlmax
  integer,private :: nbox
  real,pointer :: xs(:,:),xe(:,:),scale(:,:)
  logical,pointer :: empty_flag(:,:)
  integer,pointer :: lvlpoint(:,:)
  integer,pointer :: istart(:,:,:),iend(:,:,:)

  integer,allocatable,private :: check_alloc(:)

  private :: init_octree,localboxindex,index_ltog,interleave,parent
contains

!********************************************************************
  subroutine init_octree(ngrids,llmin,llmax)
   
    integer :: ngrids,llmin,llmax,lvl
  
    lvlmin = llmin
    lvlmax = llmax

    nbox = 0
    do lvl = 0,lvlmax
      nbox = nbox + ishft(1,3*lvl)
    enddo
   
    nullify(xs,xe,scale,empty_flag,lvlpoint,istart,iend)
    if(allocated(check_alloc)) deallocate(check_alloc)

    allocate(xs(3,ngrids),xe(3,ngrids),scale(3,ngrids))
    allocate(empty_flag(nbox,ngrids),lvlpoint(nbox,ngrids))
    allocate(istart(nbox,3,ngrids),iend(nbox,3,ngrids))

    allocate(check_alloc(1))
  
  end subroutine init_octree

!********************************************************************
  subroutine make_octree(x,y,z,jmax,kmax,lmax,ngrids,llmin,llmax,init)

!********************************************************************

    integer,optional :: init
    integer :: ngrids,llmin,llmax
    integer :: jmax(:),kmax(:),lmax(:)
    real    :: x(:,:,:,:),y(:,:,:,:),z(:,:,:,:)

!.. local variables
    integer :: j,k,l,n
    integer :: lvl,nout,npout,nmout,nbox,npbox,nsbox,nebox
    integer :: maxbox,maxlev
    real    :: eps_x,eps_y,eps_z
    real    :: xbar(3)

!...allocate variables under user request

    if(present(init)) then
      if (init.eq.1) call init_octree(ngrids,llmin,llmax)
    endif

!...allocate variables if not initiated yet

    if (.not.allocated(check_alloc)) call init_octree(ngrids,llmin,llmax)

    empty_flag = .true.
    istart = 1
    iend = 1

!... scaling

    do n = 1,ngrids
      xs(1,n) = minval(x(1:jmax(n),1:kmax(n),1:lmax(n),n))
      xe(1,n) = maxval(x(1:jmax(n),1:kmax(n),1:lmax(n),n))

      xs(2,n) = minval(y(1:jmax(n),1:kmax(n),1:lmax(n),n))
      xe(2,n) = maxval(y(1:jmax(n),1:kmax(n),1:lmax(n),n))

      xs(3,n) = minval(z(1:jmax(n),1:kmax(n),1:lmax(n),n))
      xe(3,n) = maxval(z(1:jmax(n),1:kmax(n),1:lmax(n),n))

      eps_x=(xe(1,n)-xs(1,n))*0.001
      eps_y=(xe(2,n)-xs(2,n))*0.001
      eps_z=(xe(3,n)-xs(3,n))*0.001
      
      xs(1,n) = xs(1,n) - eps_x
      xs(2,n) = xs(2,n) - eps_y
      xs(3,n) = xs(3,n) - eps_z
      xe(1,n) = xe(1,n) + eps_x
      xe(2,n) = xe(2,n) + eps_y
      xe(3,n) = xe(3,n) + eps_z

      scale(1,n) = xe(1,n) - xs(1,n)
      scale(2,n) = xe(2,n) - xs(2,n)
      scale(3,n) = xe(3,n) - xs(3,n)

      maxlev = 10

      do j=1,jmax(n)
      do k=1,kmax(n)
      do l=1,lmax(n)

        xbar(1) = (x(j,k,l,n)-xs(1,n))/scale(1,n)
        xbar(2) = (y(j,k,l,n)-xs(2,n))/scale(2,n)
        xbar(3) = (z(j,k,l,n)-xs(3,n))/scale(3,n)

        lvl = lvlmax
        call boxindex(xbar,lvl,nout,3)

        if (empty_flag(nout,n)) then
          istart(nout,1,n) = j 
          istart(nout,2,n) = k 
          istart(nout,3,n) = l 
          call boxindex(xbar,maxlev,nsbox,3)

          iend(nout,1,n) = j 
          iend(nout,2,n) = k 
          iend(nout,3,n) = l 
          nebox = nsbox
          empty_flag(nout,n) = .false.
        else
          call boxindex(xbar,maxlev,nmout,3)
          if (nmout.lt.nsbox) then
            istart(nout,1,n) = j 
            istart(nout,2,n) = k 
            istart(nout,3,n) = l 
            nsbox = nmout
          elseif (nmout.gt.nebox) then
            iend(nout,1,n) = j 
            iend(nout,2,n) = k 
            iend(nout,3,n) = l 
            nebox = nmout
          endif
        endif

      enddo
      enddo
      enddo

      do lvl = lvlmax,1,-1
        maxbox = 2**(3*lvl)
        do nbox = 0,maxbox-1
          call parent(nbox,npbox,3)
          call index_ltog(nbox,nout,lvl,3)
          call index_ltog(npbox,npout,lvl-1,3)

          if (.not.empty_flag(nout,n)) then
            if (empty_flag(npout,n)) then
              istart(npout,1,n) = istart(nout,1,n) 
              istart(npout,2,n) = istart(nout,2,n) 
              istart(npout,3,n) = istart(nout,3,n) 
            endif
            empty_flag(npout,n) = .false.
            iend(npout,1,n) = iend(nout,1,n) 
            iend(npout,2,n) = iend(nout,2,n) 
            iend(npout,3,n) = iend(nout,3,n) 
          endif
        enddo
      enddo

      maxbox = 2**(3*lvlmax)

      do nbox = 0,maxbox-1
        call index_ltog(nbox,nout,lvlmax,3)
        lvlpoint(nout,n) = lvlmax

        if (empty_flag(nout,n)) then
          npbox = nbox

          do lvl = lvlmax-1,lvlmin,-1
            call parent(npbox,npbox,3)
            call index_ltog(npbox,npout,lvl,3)
            if (empty_flag(npout,n)) then
              cycle
            else
              lvlpoint(nout,n) = lvl
              exit
            endif
          enddo 
        endif

      enddo

    enddo

  end subroutine make_octree

!*********************************************************************
  subroutine boxindex(x,l,n,d)
! finding the global box index given a point
!*********************************************************************
    integer :: n,l,d
    real :: x(d)

    call localboxindex(x,l,n,d)
    call index_ltog(n,n,l,d)

    end subroutine boxindex

!*********************************************************************
  subroutine localboxindex(x,l,n,d)
! finding the box index given a point
! This gives per level box number
!*********************************************************************
    integer :: n,l,d
    real :: x(d)

!.. local variables
    integer :: dim
    integer :: ncoord(d)

    do dim = 1,d
      ncoord(dim) = (ishft(1,l)*x(dim))
    enddo
    call interleave(ncoord,n,d)

    end subroutine localboxindex

!*********************************************************************
  subroutine index_ltog(n,m,l,d)
! convert local box index to global 
!*********************************************************************
    integer :: n,m,l,d
    
    m = n + 1
    do l = 0,l-1
      m = m + ishft(1,d*l)
    enddo
    
  end subroutine index_ltog

!*********************************************************************
  subroutine interleave(n,nout,d)
! subroutine for interleaving
!*********************************************************************
    integer :: nout,d
    integer :: n(d)

!.. local variables
    integer :: ntemp(d)
    integer :: j,dim,count,bit

    ntemp = n
    nout = 0
    count = 0

    do while (maxval(ntemp).ne.0)
      do dim = d,1,-1
        bit = ntemp(dim) - 2*ishft(ntemp(dim),-1)
        nout = nout + ishft(bit,count)
        count = count+1
        ntemp(dim) = ishft(ntemp(dim),-1)
      enddo
    enddo

  end subroutine interleave

!*********************************************************************
  subroutine parent(n,m,d)
!
! finding the parent
!*********************************************************************
   integer :: n,m,d

   m = ishft(n,-d)

   end subroutine parent

end module octree

!********************************************************************
module ihc
  implicit none
!********************************************************************

  integer,private :: ngrids

  type ptr_to_integer3D
     integer, pointer :: arr(:,:,:)
  end type ptr_to_integer3D

  type ptr_to_real3D
     real, pointer :: arr(:,:,:)
  end type ptr_to_real3D

  type(ptr_to_integer3D),allocatable,private :: cell_bc(:),ichk(:)
  type(ptr_to_real3D),allocatable,private    :: vold(:),volr(:)

  integer,pointer,private :: jmax(:),kmax(:),lmax(:)
  real,pointer,private    :: x(:,:,:,:),y(:,:,:,:),z(:,:,:,:)
  integer,pointer,private :: iblank(:,:,:,:)

  character*40,pointer,private :: bc_file(:)

  integer,allocatable,private :: check_alloc(:)

  real,private :: overlap_factor

  private :: init_ihc,read_bc,boundconnect,recvconnect,holeconnect,metfv, &
             finddonor,find_fractions,trilinearinterp,matrixinv3x3

contains

!********************************************************************
  subroutine init_ihc(jmx,kmx,lmx,ngrds,bcfile)
!********************************************************************
    integer :: ngrds
    integer :: jmx(:),kmx(:),lmx(:)
    character*40 :: bcfile(:)

!...local variables
    integer :: n

    ngrids = ngrds

    if(allocated(cell_bc)) deallocate(cell_bc)
    if(allocated(ichk)) deallocate(ichk)
    if(allocated(vold)) deallocate(vold)
    if(allocated(volr)) deallocate(volr)
    if(allocated(check_alloc)) deallocate(check_alloc)

    nullify(jmax,kmax,lmax)
    nullify(bc_file)
    nullify(x,y,z)

    allocate(jmax(ngrids),kmax(ngrids),lmax(ngrids))
    allocate(bc_file(ngrids))

    jmax = jmx
    kmax = kmx
    lmax = lmx
    bc_file = bcfile

    allocate(vold(ngrids),volr(ngrids))
    allocate(cell_bc(ngrids),ichk(ngrids))

    do n = 1,ngrids
      allocate(vold(n)%arr(jmax(n),kmax(n),lmax(n)))
      allocate(volr(n)%arr(jmax(n),kmax(n),lmax(n)))
      allocate(cell_bc(n)%arr(jmax(n),kmax(n),lmax(n)))
      allocate(ichk(n)%arr(jmax(n),kmax(n),lmax(n)))
    enddo

    overlap_factor = 0.8

    !allocate(x(jm,km,lm,ngrids),y(jm,km,lm,ngrids),z(jm,km,lm,ngrids))
    !allocate(iblank(jm,km,lm,ngrids))

    allocate(check_alloc(1))

    call read_bc
 
  end subroutine init_ihc

!********************************************************************
  subroutine read_bc
!..bc types
!...1 --> not a receiver 
!...2 --> iblank = 0  
!...3 --> is a receiver 
!...4 --> iblank.ne.0
!...5 --> wall point 
!********************************************************************
    integer :: j,k,l,n,ib
    integer :: js,je,ks,ke,ls,le
    integer :: nbc,ibtyp(25),ibdir(25)
    integer :: jbcs(25),kbcs(25),lbcs(25)
    integer :: jbce(25),kbce(25),lbce(25)
    integer :: ibproc(25)

    namelist/ihcbcinp/ nbc,ibtyp,ibdir,jbcs,jbce,kbcs,kbce,lbcs,lbce,ibproc

    do n = 1,ngrids
      cell_bc(n)%arr = 0
    enddo

    do n = 1,ngrids
      open(unit=21,file=bc_file(n),status='unknown')
      read(21,ihcbcinp)
      close(21)

      do ib = 1,nbc
        js = jbcs(ib)
        je = jbce(ib)
        ks = kbcs(ib)
        ke = kbce(ib)
        ls = lbcs(ib)
        le = lbce(ib)
        if(js.lt.0) js = jmax(n) + js + 1
        if(ks.lt.0) ks = kmax(n) + ks + 1
        if(ls.lt.0) ls = lmax(n) + ls + 1
        if(je.lt.0) je = jmax(n) + je + 1
        if(ke.lt.0) ke = kmax(n) + ke + 1
        if(le.lt.0) le = lmax(n) + le + 1

        do j = js,je
        do k = ks,ke
        do l = ls,le
          cell_bc(n)%arr(j,k,l) = ibtyp(ib)
        enddo
        enddo
        enddo
      enddo 
    enddo

  end subroutine read_bc

!*********************************************************************
  subroutine do_connectihc(xyz,iblnk,jmx,kmx,lmx,imesh,idonor,frac,&
              imap,ndonor,nfringe,idsize3,Njmx,Nkmx,Nlmx, &
              ngrds,mibc,init) ! Math: add IBC

    integer,optional :: init
    integer :: ngrds,Njmx,Nkmx,Nlmx,idsize3,mibc ! Math: add IBC
    integer :: jmx(:),kmx(:),lmx(:)
    integer :: ndonor(:),nfringe(:)
    real,target    :: xyz(:,:,:,:,:)
    integer,target :: iblnk(:,:,:,:)
    integer :: imesh(:,:,:),idonor(:,:,:)
    real    :: frac(:,:,:) 
    integer :: imap(:,:,:)
    
    integer :: n,ng,ngd,nf,nd,is,ie
    character*40 :: integer_string

    real,pointer :: x1(:,:,:,:),y1(:,:,:,:),z1(:,:,:,:)
    integer,allocatable :: imesh1(:,:,:),idonor1(:,:,:)
    integer,allocatable :: ibc1(:,:)
    real,allocatable :: frac1(:,:,:)
    integer,allocatable :: iisptr(:),iieptr(:)
    character*40,allocatable :: bcfile(:)

    allocate(imesh1(idsize3,3,ngrds),idonor1(idsize3,3,ngrds))
    allocate(ibc1(idsize3,ngrds),frac1(idsize3,3,ngrds))
    allocate(iisptr(ngrds),iieptr(ngrds))
    allocate(bcfile(ngrds))

    x1 => xyz(1,:,:,:,:)
    y1 => xyz(2,:,:,:,:)
    z1 => xyz(3,:,:,:,:)

    do n=1,ngrds
      write(integer_string,*) n-1
      integer_string = adjustl(integer_string)
      bcfile(n) = 'ihcbc.'//trim(integer_string)
    enddo

    call connect3d(x1,y1,z1,iblnk,idonor1,frac1,imesh1,ibc1,&
              iisptr,iieptr,ndonor,nfringe,jmx,kmx,lmx,ngrds,&
              bcfile,mibc,init) ! Math: add IBC

    do ng = 1,ngrds
      do nf = 1,nfringe(ng)
        do n = 1,3
          imesh(n,nf,ng) = imesh1(nf,n,ng)
        enddo 
        do ngd = 1,ngrds
          is = iisptr(ngd)
          ie = iieptr(ngd)
          if (ibc1(nf,ng).ge.is.and.ibc1(nf,ng).le.ie) then
            imap(1,nf,ng) = ibc1(nf,ng) - is + 1
            imap(2,nf,ng) = ngd
            exit
          endif
        enddo
      enddo

      do nd = 1,ndonor(ng)
        do n = 1,3
          idonor(n,nd,ng) = idonor1(nd,n,ng)
          frac(n,nd,ng) = frac1(nd,n,ng)
        enddo 
      enddo
    enddo 

  end subroutine do_connectihc

!*********************************************************************
  subroutine connect3d(x1,y1,z1,iblnk,idonor,frac,imesh,ibc,&
              iisptr,iieptr,ndonor,nfringe,jmx,kmx,lmx,ngrds,bcfile,mibc,init) ! Math: add IBC
    use octree
    use immersedBoundVars ! Math: iblk_ibc
!********************************************************************

    integer,optional :: init
    integer :: ngrds,mibc ! Math: add IBC
    integer :: jmx(:),kmx(:),lmx(:)
    integer :: ndonor(:),nfringe(:)
    integer :: iisptr(:),iieptr(:)
    real,target :: x1(:,:,:,:),y1(:,:,:,:),z1(:,:,:,:)
    integer,target :: iblnk(:,:,:,:)
    integer :: imesh(:,:,:),idonor(:,:,:)
    real    :: frac(:,:,:) 
    integer :: ibc(:,:)
    character*40 :: bcfile(:)

    integer jkl
    integer j,k,l,j1,k1,l1,m,n,nn,n1
    integer idsize,count,llmin,llmax
    integer,allocatable :: nholept(:),iholept(:,:,:)
    integer,allocatable :: ibctmp(:,:),idontmp(:,:)
    real,allocatable :: diholept(:,:,:)
    logical,allocatable :: ioverlap(:,:)

!...allocate variables under user request

    if(present(init)) then
      if (init.eq.1) call init_ihc(jmx,kmx,lmx,ngrds,bcfile)
    endif

!...allocate variables if not initiated yet

    if(.not.allocated(check_alloc)) call init_ihc(jmx,kmx,lmx,ngrds,bcfile)
     
!...set pointers

    x => x1 ! x,y,z set as pointers to avoid additional memory usage
    y => y1 ! have to be careful not to change x,y,z in any of the
    z => z1 ! subroutines associated with this module 

    iblank => iblnk ! iblank values do change in this module

!...some allocatation

    idsize = 0
    do n = 1,ngrids
      jkl = jmax(n)*kmax(n)*lmax(n)
      if(jkl.gt.idsize) idsize = jkl 
    enddo

    allocate(nholept(ngrids))
    allocate(iholept(idsize,7,ngrids))
    allocate(diholept(idsize,3,ngrids))
    allocate(ibctmp(idsize,ngrids),idontmp(idsize,ngrids))

    allocate(ioverlap(ngrids,ngrids))

!...calculating volume of the cells

    call metfv

!...reassign cell volume based on the bc

    do n = 1,ngrids
      do j = 1,jmax(n)
      do k = 1,kmax(n)
      do l = 1,lmax(n)
        if (cell_bc(n)%arr(j,k,l).eq.1) volr(n)%arr(j,k,l) = 0.
        if (cell_bc(n)%arr(j,k,l).eq.1) vold(n)%arr(j,k,l) = 0.
        if (cell_bc(n)%arr(j,k,l).eq.2) volr(n)%arr(j,k,l) = -1.e30
        if (cell_bc(n)%arr(j,k,l).eq.2) vold(n)%arr(j,k,l) = 1.e30
        if (cell_bc(n)%arr(j,k,l).eq.3) volr(n)%arr(j,k,l) = 1.e30
        if (cell_bc(n)%arr(j,k,l).eq.3) vold(n)%arr(j,k,l) = 1.e30
        if (cell_bc(n)%arr(j,k,l).eq.5) volr(n)%arr(j,k,l) = 0.
        if (cell_bc(n)%arr(j,k,l).eq.5) vold(n)%arr(j,k,l) = 0.
      enddo
      enddo
      enddo
    enddo

!...generate octree data
    llmin = 4
    llmax = 5
    call make_octree(x,y,z,jmax,kmax,lmax,ngrids,llmin,llmax)

!...initialize iblank array

    do n = 1,ngrids
      do j = 1,jmax(n)
      do k = 1,kmax(n)
      do l = 1,lmax(n)
         ! Math: add IBC
         if (n==mibc) then
            iblank(j,k,l,n) = iblk_ibc(j,k,l)
         else
            iblank(j,k,l,n) = 1
         endif
         if (cell_bc(n)%arr(j,k,l).eq.2) iblank(j,k,l,n) = 0
      enddo
      enddo
      enddo
    enddo

    do n = 1,ngrids
      ichk(n)%arr = 0
    enddo

!...check for overlapping meshes 

    ioverlap = .false.

    do n = 1,ngrids
      call boundconnect(n,ioverlap)
    enddo

!...do connectivity for forced receivers first

    nholept = 0

    do n = 1,ngrids
      call recvconnect(n,ioverlap,nholept,iholept,diholept,idsize)
    enddo

!...do connectivity for all other points

    do n = 1,ngrids
      call holeconnect(n,ioverlap,nholept,iholept,diholept,idsize)
    enddo

!... get connectivity info....
    ndonor = 0
    nfringe = 0

    do n=1,ngrids
     do m=1,nholept(n)
       n1 = iholept(m,4,n)

       j = iholept(m,1,n)
       k = iholept(m,2,n)
       l = iholept(m,3,n)
       ! Math: add IBC
       j1 = iholept(m,5,n)
       k1 = iholept(m,6,n)
       l1 = iholept(m,7,n)
       if (n1==mibc.and.iblank(j1,k1,l1,n1)==0) iblank(j,k,l,n)=0

       if (iblank(j,k,l,n).ne.0) then
         j1 = iholept(m,5,n)
         k1 = iholept(m,6,n)
         l1 = iholept(m,7,n)
         nfringe(n) = nfringe(n) + 1
         imesh(nfringe(n),1,n) = j
         imesh(nfringe(n),2,n) = k
         imesh(nfringe(n),3,n) = l

         ndonor(n1) = ndonor(n1) + 1
         ibctmp(nfringe(n),n) = ndonor(n1)
         idontmp(nfringe(n),n) = n1

         idonor(ndonor(n1),1,n1) = j1
         idonor(ndonor(n1),2,n1) = k1
         idonor(ndonor(n1),3,n1) = l1
         frac(ndonor(n1),1,n1) = diholept(m,1,n)
         frac(ndonor(n1),2,n1) = diholept(m,2,n)
         frac(ndonor(n1),3,n1) = diholept(m,3,n)
         iblank(j,k,l,n) = -1
       endif
     enddo
    enddo
    
    !global donor cell ptrs
    iisptr = 0
    iieptr = 0
    iisptr(1) = 1
    iieptr(1) = ndonor(1)
    do n=2,ngrids
      iisptr(n) = iisptr(n-1) + ndonor(n-1) 
      iieptr(n) = iisptr(n) + ndonor(n) -1 
    end do 

    do n=1,ngrids
     do j=1,nfringe(n)
       n1 = idontmp(j,n)
       count = 0
       if(n1.gt.1) then
        do nn=2,n1
         count = count+ndonor(nn-1)
        end do
       end if
       ibc(j,n) = count+ibctmp(j,n)
     end do
    end do

  end subroutine connect3d

!********************************************************************
  subroutine boundconnect(ng,ioverlap)
    use octree
!********************************************************************

    integer :: ng
    logical :: ioverlap(ngrids,ngrids)

    integer :: j,k,l,n,nface
    integer :: lvl,lvlcheck,nlvl,nout
    integer :: jl,kl,ll,js,ks,ls,jstep,kstep,lstep

    real    :: xbar(3),xp(3),frac1(3)
    real,allocatable :: xg(:,:,:,:)

    integer :: nstop,ntime
    integer,allocatable :: jst(:),kst(:),lst(:)
    integer,allocatable :: jet(:),ket(:),let(:)
    integer,allocatable :: box(:)

    allocate(jst(2),kst(2),lst(2))
    allocate(jet(2),ket(2),let(2))
    allocate(box(2))

    do n = 1,ngrids

    if (ioverlap(n,ng)) ioverlap(ng,n) = .true.

    if (n.ne.ng.and..not.ioverlap(ng,n)) then
      allocate(xg(jmax(n),kmax(n),lmax(n),3))
      xg(:,:,:,1) = x(:,:,:,n)
      xg(:,:,:,2) = y(:,:,:,n)
      xg(:,:,:,3) = z(:,:,:,n)

      jl = 2; kl = 2; ll =  2
      js = -1; ks = -1; ls = -1

      j = 1; k = 1; l = 1
      jstep = 1; kstep = 1; lstep = 1
      
      nface = 1
      do while (1) !loop for all boundary points 

        if (nface.eq.1) j = 1
        if (nface.eq.2) j = jmax(ng)
        if (nface.eq.3) k = 1
        if (nface.eq.4) k = kmax(ng)
        if (nface.eq.5) l = 1
        if (nface.eq.6) l = lmax(ng)

        xbar(1) = (x(j,k,l,ng)-xs(1,n))/scale(1,n)
        xbar(2) = (y(j,k,l,ng)-xs(2,n))/scale(2,n)
        xbar(3) = (z(j,k,l,ng)-xs(3,n))/scale(3,n)

        if (xbar(1).ge.0.and.xbar(2).ge.0.and.xbar(3).ge.0.and. &
            xbar(1).lt.1.and.xbar(2).lt.1.and.xbar(3).lt.1) then

          call boxindex(xbar,lvlmin,nout,3)
          if (empty_flag(nout,n)) go to 10

          xp(1) = x(j,k,l,ng)
          xp(2) = y(j,k,l,ng)
          xp(3) = z(j,k,l,ng)

          js = jl
          ks = kl
          ls = ll
          call finddonor(xp,xg,js,ks,ls,jmax(n),kmax(n),lmax(n),100)

          if (js.le.0) then

            call boxindex(xbar,lvlmax,nout,3)
            lvlcheck = lvlpoint(nout,n)
            nlvl = 2
            if (lvlcheck .eq. lvlmin) nlvl = 1

            do lvl = 1,nlvl 
              call boxindex(xbar,lvlcheck-lvl+1,nout,3)
              box(lvl) = nout
              jst(lvl) = istart(nout,1,n)
              kst(lvl) = istart(nout,2,n)
              lst(lvl) = istart(nout,3,n)
              jet(lvl) = iend(nout,1,n)
              ket(lvl) = iend(nout,2,n)
              let(lvl) = iend(nout,3,n)
            enddo

            nstop = 0
            ntime = 0

            do while(nstop.lt.2*nlvl.and.ntime.lt.1)

              do lvl = 1,nlvl 

                nout = box(lvl)

                if (js.le.0.and.kst(lvl).gt.0) then
                  js = abs(jst(lvl))
                  ks = abs(kst(lvl))
                  ls = abs(lst(lvl))
                  call finddonor(xp,xg,js,ks,ls,jmax(n),kmax(n),lmax(n),100)
                  jst(lvl) = js
                  kst(lvl) = ks
                  lst(lvl) = ls

                  if (js.le.0.and.ks.le.0) then
                    nstop = nstop + 1
                  endif
                else
                  nstop = nstop + 1
                endif

                if (js.le.0.and.ket(lvl).gt.0) then
                  js = abs(jet(lvl))
                  ks = abs(ket(lvl))
                  ls = abs(let(lvl))
                  call finddonor(xp,xg,js,ks,ls,jmax(n),kmax(n),lmax(n),100)
                  jet(lvl) = js
                  ket(lvl) = ks
                  let(lvl) = ls

                  if (js.le.0.and.ks.le.0) then
                    nstop = nstop + 1
                  endif
                else
                  nstop = nstop + 1
                endif

              enddo
              ntime = ntime + 1

            enddo
          endif

          if (js.gt.0) then
            jl = js
            kl = ks
            ll = ls
            ioverlap(ng,n) = .true.
            ioverlap(n,ng) = .true.
            exit
          endif
        endif

 10     if (nface.eq.1.or.nface.eq.2) then
          l = l + lstep 
          if (l.gt.lmax(ng).or.l.lt.1) then
            lstep = -1*lstep
            l = l + lstep
            k = k + kstep
            if (k.gt.kmax(ng).or.k.lt.1) then
              nface = nface + 1
              j = 1; k = 1; l = 1
              jstep = 1; kstep = 1; lstep = 1
              cycle
            endif
          endif
        endif

        if (nface.eq.3.or.nface.eq.4) then
          j = j + jstep 
          if (j.gt.jmax(ng).or.j.lt.1) then
            jstep = -1*jstep
            j = j + jstep
            l = l + lstep
            if (l.gt.lmax(ng).or.l.lt.1) then
              nface = nface + 1
              j = 1; k = 1; l = 1
              jstep = 1; kstep = 1; lstep = 1
              cycle
            endif
          endif
        endif

        if (nface.eq.5.or.nface.eq.6) then
          k = k + kstep 
          if (k.gt.kmax(ng).or.k.lt.1) then
            kstep = -1*kstep
            k = k + kstep
            j = j + jstep
            if (j.gt.jmax(ng).or.j.lt.1) then
              nface = nface + 1
              j = 1; k = 1; l = 1
              jstep = 1; kstep = 1; lstep = 1
              cycle
            endif
          endif
        endif

        if (nface.eq.7) exit
      enddo

      deallocate(xg)
    endif

    enddo

  end subroutine boundconnect

!********************************************************************
  subroutine recvconnect(ng,ioverlap,nholept,iholept,diholept,idsize)
    use octree
!********************************************************************

    integer :: ng,idsize
    integer :: nholept(ngrids)
    integer :: iholept(idsize,7,ngrids)
    real    :: diholept(idsize,3,ngrids)
    logical :: ioverlap(ngrids,ngrids)

    integer :: j,k,l,jm1,km1,lm1,jp1,kp1,lp1,n,nc
    integer :: lvl,nlvl,lvlcheck,nout,ncount,nbodycross,nblank
    integer :: jl,kl,ll,js,ks,ls,jstep,kstep,lstep
    real    :: voldonor
    real    :: scalei1,scalei2,scalei3

    real    :: xbar(3),xp(3),frac1(3)
    integer,allocatable :: idonorptr(:,:,:)
    real,allocatable :: xg(:,:,:,:)
    real,allocatable :: volrecv(:,:,:)
 
    integer :: nstop,ntime
    integer,allocatable :: jst(:),kst(:),lst(:)
    integer,allocatable :: jet(:),ket(:),let(:)
    integer,allocatable :: box(:)

    allocate(jst(2),kst(2),lst(2))
    allocate(jet(2),ket(2),let(2))
    allocate(box(2))

    allocate(volrecv(jmax(ng),kmax(ng),lmax(ng)))
    allocate(idonorptr(jmax(ng),kmax(ng),lmax(ng)))

    volrecv(:,:,:) = overlap_factor*volr(ng)%arr(:,:,:)

    idonorptr = 0

    ncount = nholept(ng)

    do n = 1,ngrids
    if (n.ne.ng.and.ioverlap(ng,n)) then
      allocate(xg(jmax(n),kmax(n),lmax(n),3))
      xg(:,:,:,1) = x(:,:,:,n)
      xg(:,:,:,2) = y(:,:,:,n)
      xg(:,:,:,3) = z(:,:,:,n)

      jl = 2; kl = 2; ll =  2
      js = -1; ks = -1; ls = -1

      j = 1; k = 1; l = 1
      jstep = 1; kstep = 1; lstep = 1
      
      scalei1 = 1./scale(1,n)
      scalei2 = 1./scale(2,n)
      scalei3 = 1./scale(3,n)

      do while (1) !loop for all points 

        nbodycross = 0

        if (cell_bc(ng)%arr(j,k,l).ne.3) go to 10

        xbar(1) = (x(j,k,l,ng)-xs(1,n))*scalei1
        xbar(2) = (y(j,k,l,ng)-xs(2,n))*scalei2
        xbar(3) = (z(j,k,l,ng)-xs(3,n))*scalei3

        if (volrecv(j,k,l).ge.0..and. &
            xbar(1).ge.0.and.xbar(2).ge.0.and.xbar(3).ge.0.and. &
            xbar(1).lt.1.and.xbar(2).lt.1.and.xbar(3).lt.1) then

          call boxindex(xbar,lvlmin,nout,3)
          if (empty_flag(nout,n)) go to 10

          xp(1) = x(j,k,l,ng)
          xp(2) = y(j,k,l,ng)
          xp(3) = z(j,k,l,ng)

          js = jl
          ks = kl
          ls = ll
          call finddonor(xp,xg,js,ks,ls,jmax(n),kmax(n),lmax(n),100)

          if (js.le.0) then
            if (ks.le.0.and.cell_bc(n)%arr(abs(js),abs(ks),abs(ls)).eq.5) then
              nbodycross = nbodycross + 1
            endif

            call boxindex(xbar,lvlmax,nout,3)
            lvlcheck = lvlpoint(nout,n)
            nlvl = 2
            if (lvlcheck .eq. lvlmin) nlvl = 1

            do lvl = 1,nlvl
              call boxindex(xbar,lvlcheck-lvl+1,nout,3)
              box(lvl) = nout
              jst(lvl) = istart(nout,1,n)
              kst(lvl) = istart(nout,2,n)
              lst(lvl) = istart(nout,3,n)
              jet(lvl) = iend(nout,1,n)
              ket(lvl) = iend(nout,2,n)
              let(lvl) = iend(nout,3,n)
            enddo

            nstop = 0
            ntime = 0

            do while(nstop.lt.2*nlvl.and.ntime.lt.5)

              do lvl = 1,nlvl

                nout = box(lvl)

                if (js.le.0.and.kst(lvl).gt.0) then
                  js = abs(jst(lvl))
                  ks = abs(kst(lvl))
                  ls = abs(lst(lvl))
                  call finddonor(xp,xg,js,ks,ls,jmax(n),kmax(n),lmax(n),100)
                  jst(lvl) = js
                  kst(lvl) = ks
                  lst(lvl) = ls

                  if (js.le.0.and.ks.le.0) then
                    if (cell_bc(n)%arr(abs(js),abs(ks),abs(ls)).eq.5) then
                      nbodycross = nbodycross + 1
                    endif
                    nstop = nstop + 1
                  endif
                else
                  nstop = nstop + 1
                endif

                if (js.le.0.and.ket(lvl).gt.0) then
                  js = abs(jet(lvl))
                  ks = abs(ket(lvl))
                  ls = abs(let(lvl))
                  call finddonor(xp,xg,js,ks,ls,jmax(n),kmax(n),lmax(n),100)
                  jet(lvl) = js
                  ket(lvl) = ks
                  let(lvl) = ls

                  if (js.le.0.and.ks.le.0) then
                    if (cell_bc(n)%arr(abs(js),abs(ks),abs(ls)).eq.5) then
                      nbodycross = nbodycross + 1
                    endif
                    nstop = nstop + 1
                  endif
                else
                  nstop = nstop + 1
                endif

              enddo
              ntime = ntime + 1

            enddo 
          endif

          if (js.gt.0) then

            voldonor = vold(n)%arr(js,ks,ls)
            if (cell_bc(n)%arr(js  , ks  , ls  ).eq.3.or. &
                cell_bc(n)%arr(js  , ks  , ls+1).eq.3.or. &
                cell_bc(n)%arr(js  , ks+1, ls  ).eq.3.or. &
                cell_bc(n)%arr(js  , ks+1, ls+1).eq.3.or. &
                cell_bc(n)%arr(js+1, ks  , ls  ).eq.3.or. &
                cell_bc(n)%arr(js+1, ks  , ls+1).eq.3.or. &
                cell_bc(n)%arr(js+1, ks+1, ls  ).eq.3.or. &
                cell_bc(n)%arr(js+1, ks+1, ls+1).eq.3) voldonor = 1.e30

            if (voldonor.lt.volrecv(j,k,l)) then
              volrecv(j,k,l) = voldonor
              if (idonorptr(j,k,l).eq.0) then
                ncount = ncount + 1
                nc = ncount
                nholept(ng) = ncount
                idonorptr(j,k,l) = ncount
              else
                nc = idonorptr(j,k,l)
              endif
              iholept(nc,1,ng) = j
              iholept(nc,2,ng) = k
              iholept(nc,3,ng) = l
              iholept(nc,4,ng) = n
              iholept(nc,5,ng) = js
              iholept(nc,6,ng) = ks
              iholept(nc,7,ng) = ls
              jm1 = max(j-1,1)
              km1 = max(k-1,1)
              lm1 = max(l-1,1)
              vold(ng)%arr(jm1, km1, lm1) = 1.e30
              vold(ng)%arr(jm1, km1, l  ) = 1.e30
              vold(ng)%arr(jm1, k  , lm1) = 1.e30
              vold(ng)%arr(jm1, k  , l  ) = 1.e30
              vold(ng)%arr(j  , km1, lm1) = 1.e30
              vold(ng)%arr(j  , km1, l  ) = 1.e30
              vold(ng)%arr(j  , k  , lm1) = 1.e30
              vold(ng)%arr(j  , k  , l  ) = 1.e30
              call find_fractions(xp,xg,frac1,js,ks,ls, &
                                          jmax(n),kmax(n),lmax(n))
              diholept(nc,1,ng) = frac1(1)
              diholept(nc,2,ng) = frac1(2)
              diholept(nc,3,ng) = frac1(3)
            endif

            jl = js
            kl = ks
            ll = ls
          endif

          if (js.le.0.and.nbodycross.gt.2*nlvl) then
            if (cell_bc(ng)%arr(j,k,l).ne.4) then
              iblank(j,k,l,ng) = 0
              volrecv(j,k,l) = -1.e30
            endif
          endif

          ichk(ng)%arr(j,k,l) = 1
        endif

 10     l = l + lstep 
        if (l.gt.lmax(ng).or.l.lt.1) then
          lstep = -1*lstep
          l = l + lstep
          k = k + kstep
          if (k.gt.kmax(ng).or.k.lt.1) then
            kstep = -1*kstep
            k = k + kstep
            j = j + jstep
          endif
        endif
        if (j.gt.jmax(ng)) exit
      enddo

      deallocate(xg)
    endif
    enddo

    do nc = 1,ncount
      n  = iholept(nc,4,ng)
      js = iholept(nc,5,ng)
      ks = iholept(nc,6,ng)
      ls = iholept(nc,7,ng)
      jp1 = min(js+1,jmax(n))
      kp1 = min(ks+1,kmax(n))
      lp1 = min(ls+1,lmax(n))
      volr(n)%arr(js , ks , ls ) = 0.
      volr(n)%arr(js , ks , lp1) = 0.
      volr(n)%arr(js , kp1, ls ) = 0.
      volr(n)%arr(js , kp1, lp1) = 0.
      volr(n)%arr(jp1, ks , ls ) = 0.
      volr(n)%arr(jp1, ks , lp1) = 0.
      volr(n)%arr(jp1, kp1, ls ) = 0.
      volr(n)%arr(jp1, kp1, lp1) = 0.
    enddo

  end subroutine recvconnect

!********************************************************************
  subroutine holeconnect(ng,ioverlap,nholept,iholept,diholept,idsize)
    use octree
!********************************************************************

    integer :: ng,idsize
    integer :: nholept(ngrids)
    integer :: iholept(idsize,7,ngrids)
    real    :: diholept(idsize,3,ngrids)
    logical :: ioverlap(ngrids,ngrids)

    integer :: j,k,l,jm1,km1,lm1,jp1,kp1,lp1,n,nc
    integer :: lvl,nlvl,lvlcheck,nout,ncount,nbodycross,nblank
    integer :: jl,kl,ll,js,ks,ls,jstep,kstep,lstep
    real    :: voldonor
    real    :: scalei1,scalei2,scalei3

    real    :: xbar(3),xp(3),frac1(3)
    integer,allocatable :: idonorptr(:,:,:)
    real,allocatable :: xg(:,:,:,:)
    real,allocatable :: volrecv(:,:,:)

    integer :: nstop,ntime
    integer,allocatable :: jst(:),kst(:),lst(:)
    integer,allocatable :: jet(:),ket(:),let(:)
    integer,allocatable :: box(:)

    allocate(jst(2),kst(2),lst(2))
    allocate(jet(2),ket(2),let(2))
    allocate(box(2))

    allocate(volrecv(jmax(ng),kmax(ng),lmax(ng)))
    allocate(idonorptr(jmax(ng),kmax(ng),lmax(ng)))

    volrecv(:,:,:) = overlap_factor*volr(ng)%arr(:,:,:)

    idonorptr = 0

    ncount = nholept(ng)

    do n = 1,ngrids
    if (n.ne.ng.and.ioverlap(ng,n)) then
      allocate(xg(jmax(n),kmax(n),lmax(n),3))
      xg(:,:,:,1) = x(:,:,:,n)
      xg(:,:,:,2) = y(:,:,:,n)
      xg(:,:,:,3) = z(:,:,:,n)

      jl = 2; kl = 2; ll =  2
      js = -1; ks = -1; ls = -1

      j = 1; k = 1; l = 1
      jstep = 1; kstep = 1; lstep = 1
      
      scalei1 = 1./scale(1,n)
      scalei2 = 1./scale(2,n)
      scalei3 = 1./scale(3,n)

      do while (1) !loop for all points 

        nbodycross = 0

        if (volrecv(j,k,l).lt.0..or.ichk(ng)%arr(j,k,l).eq.1) go to 10

        xbar(1) = (x(j,k,l,ng)-xs(1,n))*scalei1
        xbar(2) = (y(j,k,l,ng)-xs(2,n))*scalei2
        xbar(3) = (z(j,k,l,ng)-xs(3,n))*scalei3

        if (xbar(1).ge.0.and.xbar(2).ge.0.and.xbar(3).ge.0.and. &
            xbar(1).lt.1.and.xbar(2).lt.1.and.xbar(3).lt.1) then

          call boxindex(xbar,lvlmin,nout,3)
          if (empty_flag(nout,n)) go to 10

          xp(1) = x(j,k,l,ng)
          xp(2) = y(j,k,l,ng)
          xp(3) = z(j,k,l,ng)

          js = jl
          ks = kl
          ls = ll
          call finddonor(xp,xg,js,ks,ls,jmax(n),kmax(n),lmax(n),1000)

          if (js.le.0) then

            if (ks.le.0.and.cell_bc(n)%arr(abs(js),abs(ks),abs(ls)).eq.5) then
              nbodycross = nbodycross + 1
            endif

            call boxindex(xbar,lvlmax,nout,3)
            lvlcheck = lvlpoint(nout,n)
            nlvl = 2
            if (lvlcheck .eq. lvlmin) nlvl = 1

            do lvl = 1,nlvl 
              call boxindex(xbar,lvlcheck-lvl+1,nout,3)
              box(lvl) = nout
              jst(lvl) = istart(nout,1,n)
              kst(lvl) = istart(nout,2,n)
              lst(lvl) = istart(nout,3,n)
              jet(lvl) = iend(nout,1,n)
              ket(lvl) = iend(nout,2,n)
              let(lvl) = iend(nout,3,n)
            enddo

            nstop = 0
            ntime = 0

            do while(nstop.lt.2*nlvl.and.ntime.lt.5)

              do lvl = 1,nlvl 

                nout = box(lvl)

                if (js.le.0.and.kst(lvl).gt.0) then
                  js = abs(jst(lvl))
                  ks = abs(kst(lvl))
                  ls = abs(lst(lvl))
                  call finddonor(xp,xg,js,ks,ls,jmax(n),kmax(n),lmax(n),100)
                  jst(lvl) = js
                  kst(lvl) = ks
                  lst(lvl) = ls

                  if (js.le.0.and.ks.le.0) then
                    if (cell_bc(n)%arr(abs(js),abs(ks),abs(ls)).eq.5) then
                      nbodycross = nbodycross + 1
                    endif
                    nstop = nstop + 1
                  endif
                else
                  nstop = nstop + 1
                endif

                if (js.le.0.and.ket(lvl).gt.0) then
                  js = abs(jet(lvl))
                  ks = abs(ket(lvl))
                  ls = abs(let(lvl))
                  call finddonor(xp,xg,js,ks,ls,jmax(n),kmax(n),lmax(n),100)
                  jet(lvl) = js
                  ket(lvl) = ks
                  let(lvl) = ls

                  if (js.le.0.and.ks.le.0) then
                    if (cell_bc(n)%arr(abs(js),abs(ks),abs(ls)).eq.5) then
                      nbodycross = nbodycross + 1
                    endif
                    nstop = nstop + 1
                  endif
                else
                  nstop = nstop + 1
                endif

              enddo
              ntime = ntime + 1

            enddo
          endif

          if (js.gt.0) then

            voldonor = vold(n)%arr(js,ks,ls)
            if (cell_bc(n)%arr(js  , ks  , ls).eq.3.or. &
               cell_bc(n)%arr(js  , ks  , ls+1).eq.3.or. &
               cell_bc(n)%arr(js  , ks+1, ls  ).eq.3.or. &
               cell_bc(n)%arr(js  , ks+1, ls+1).eq.3.or. &
               cell_bc(n)%arr(js+1, ks  , ls  ).eq.3.or. &
               cell_bc(n)%arr(js+1, ks  , ls+1).eq.3.or. &
               cell_bc(n)%arr(js+1, ks+1, ls  ).eq.3.or. &
               cell_bc(n)%arr(js+1, ks+1, ls+1).eq.3) voldonor = 1.e30

            if (voldonor.lt.volrecv(j,k,l)) then
              volrecv(j,k,l) = voldonor
              if (idonorptr(j,k,l).eq.0) then
                ncount = ncount + 1
                nc = ncount
                nholept(ng) = ncount
                idonorptr(j,k,l) = ncount
              else
                nc = idonorptr(j,k,l)
              endif
              iholept(nc,1,ng) = j
              iholept(nc,2,ng) = k
              iholept(nc,3,ng) = l
              iholept(nc,4,ng) = n
              iholept(nc,5,ng) = js
              iholept(nc,6,ng) = ks
              iholept(nc,7,ng) = ls
              jm1 = max(j-1,1)
              km1 = max(k-1,1)
              lm1 = max(l-1,1)
              vold(ng)%arr(jm1, km1, lm1) = 1.e30
              vold(ng)%arr(jm1, km1, l  ) = 1.e30
              vold(ng)%arr(jm1, k  , lm1) = 1.e30
              vold(ng)%arr(jm1, k  , l  ) = 1.e30
              vold(ng)%arr(j  , km1, lm1) = 1.e30
              vold(ng)%arr(j  , km1, l  ) = 1.e30
              vold(ng)%arr(j  , k  , lm1) = 1.e30
              vold(ng)%arr(j  , k  , l  ) = 1.e30
              call find_fractions(xp,xg,frac1,js,ks,ls, &
                                      jmax(n),kmax(n),lmax(n))
              diholept(nc,1,ng) = frac1(1)
              diholept(nc,2,ng) = frac1(2)
              diholept(nc,3,ng) = frac1(3)
            endif

            jl = js
            kl = ks
            ll = ls
          endif

          if (js.le.0.and.nbodycross.gt.2*nlvl) then
            if (cell_bc(ng)%arr(j,k,l).ne.4) then
              iblank(j,k,l,ng) = 0
              volrecv(j,k,l) = -1.e30
            endif
          endif
        endif

  10    l = l + lstep 
        if (l.gt.lmax(ng).or.l.lt.1) then
          lstep = -1*lstep
          l = l + lstep
          k = k + kstep
          if (k.gt.kmax(ng).or.k.lt.1) then
            kstep = -1*kstep
            k = k + kstep
            j = j + jstep
          endif
        endif
        if (j.gt.jmax(ng)) exit
      enddo

      deallocate(xg)
    endif
    enddo

    print*,'Number of fringe points in mesh',ng,'is', ncount

    do nc = 1,ncount
      n  = iholept(nc,4,ng)
      js = iholept(nc,5,ng)
      ks = iholept(nc,6,ng)
      ls = iholept(nc,7,ng)
      jp1 = min(js+1,jmax(n))
      kp1 = min(ks+1,kmax(n))
      lp1 = min(ls+1,lmax(n))
      volr(n)%arr(js , ks , ls ) = 0.
      volr(n)%arr(js , ks , lp1) = 0.
      volr(n)%arr(js , kp1, ls ) = 0.
      volr(n)%arr(js , kp1, lp1) = 0.
      volr(n)%arr(jp1, ks , ls ) = 0.
      volr(n)%arr(jp1, ks , lp1) = 0.
      volr(n)%arr(jp1, kp1, ls ) = 0.
      volr(n)%arr(jp1, kp1, lp1) = 0.
    enddo

  end subroutine holeconnect

!********************************************************************
   subroutine metfv
!********************************************************************
      
    integer :: j,k,l,n,jm,km,lm,jp1,kp1,lp1,jmx,kmx,lmx
    integer :: ll,lm1,kk,km1,jj,jm1,nneg
    real :: dx11,dy11,dz11,dx22,dy22,dz22,sx1,sy1,sz1
    real :: fac1,fac2

    integer,allocatable :: jjp(:),jjr(:),kkp(:),kkr(:),llp(:),llr(:)
    real,allocatable :: xx(:,:,:),xy(:,:,:),xz(:,:,:)
    real,allocatable :: yx(:,:,:),yy(:,:,:),yz(:,:,:)
    real,allocatable :: zx(:,:,:),zy(:,:,:),zz(:,:,:)
    real,allocatable :: a3(:,:,:)

    jm = 0
    km = 0
    lm = 0
    do n = 1,ngrids
      if(jmax(n).gt.jm) jm = jmax(n) 
      if(kmax(n).gt.km) km = kmax(n) 
      if(lmax(n).gt.lm) lm = lmax(n) 
    enddo

    allocate(jjp(jm),jjr(jm),kkp(km),kkr(km),llp(lm),llr(lm))
    allocate(xx(jm,km,lm),xy(jm,km,lm),xz(jm,km,lm))
    allocate(yx(jm,km,lm),yy(jm,km,lm),yz(jm,km,lm))
    allocate(zx(jm,km,lm),zy(jm,km,lm),zz(jm,km,lm))
    allocate(a3(jm,km,lm))

    do n = 1,ngrids
      jmx = jmax(n)
      kmx = kmax(n)
      lmx = lmax(n)

      fac1    = 0.5
      fac2    = 1./3.

      do  j = 1,jmx
        jjp(j) = j + 1
        jjr(j) = j - 1
      enddo
      jjp(jmx) = jmx
      jjr(1) = 1

      do  k = 1,kmx
        kkp(k) = k + 1
        kkr(k) = k - 1
      enddo
      kkp(kmx) = kmx
      kkr(1) = 1

      do l = 1,lmx
        llp(l) = l + 1
        llr(l) = l - 1
      enddo
      llp(lmx) = lmx
      llr(1) = 1
!
!..xi derivatives
!
      do l = 1,lmx-1
      lp1 = llp(l)
        do k = 1,kmx-1
          kp1 = kkp(k)
          do j = 1,jmx
            dx11 = (x(j,kp1,lp1,n) -x(j,k,l,n))
            dy11 = (y(j,kp1,lp1,n) -y(j,k,l,n))
            dz11 = (z(j,kp1,lp1,n) -z(j,k,l,n))
            dx22 = (x(j,k,lp1,n) -x(j,kp1,l,n))
            dy22 = (y(j,k,lp1,n) -y(j,kp1,l,n))
            dz22 = (z(j,k,lp1,n) -z(j,kp1,l,n))

            xx(j,k,l) = fac1*( dy11*dz22 - dy22*dz11 )
            xy(j,k,l) = fac1*( dz11*dx22 - dz22*dx11 )
            xz(j,k,l) = fac1*( dx11*dy22 - dx22*dy11 )
          enddo
        enddo
      enddo
!
!..eta derivatives
!
      do l = 1,lmx-1
        lp1 = llp(l)
        do j = 1,jmx-1
          jp1 = jjp(j)
          do k = 1,kmx
            dx11 = (x(jp1,k,lp1,n) - x(j,k,l,n))
            dy11 = (y(jp1,k,lp1,n) - y(j,k,l,n))
            dz11 = (z(jp1,k,lp1,n) - z(j,k,l,n))
            dx22 = (x(jp1,k,l,n) - x(j,k,lp1,n))
            dy22 = (y(jp1,k,l,n) - y(j,k,lp1,n))
            dz22 = (z(jp1,k,l,n) - z(j,k,lp1,n))

            yx(j,k,l) = fac1*( dy11*dz22 - dy22*dz11 )
            yy(j,k,l) = fac1*( dz11*dx22 - dz22*dx11 )
            yz(j,k,l) = fac1*( dx11*dy22 - dx22*dy11 )
          enddo
        enddo
      enddo
!
!..zeta derivative
!
      do k = 1,kmx-1
        kp1 = kkp(k)
        do j = 1,jmx-1
          jp1 = jjp(j)
          do l = 1,lmx
            dx11 = (x(jp1,kp1,l,n) - x(j,k,l,n))
            dy11 = (y(jp1,kp1,l,n) - y(j,k,l,n))
            dz11 = (z(jp1,kp1,l,n) - z(j,k,l,n))
            dx22 = (x(j,kp1,l,n) - x(jp1,k,l,n))
            dy22 = (y(j,kp1,l,n) - y(jp1,k,l,n))
            dz22 = (z(j,kp1,l,n) - z(jp1,k,l,n))

            zx(j,k,l) = fac1*( dy11*dz22 - dy22*dz11 )
            zy(j,k,l) = fac1*( dz11*dx22 - dz22*dx11 )
            zz(j,k,l) = fac1*( dx11*dy22 - dx22*dy11 )
          enddo
        enddo
      enddo
    
!
!..compute cell volume
!
      do l = 1,lmx-1
        lp1 = llp(l)
        do k = 1,kmx-1
          kp1 = kkp(k)
          do j = 1,jmx-1
            jp1 = jjp(j)

            sx1 = xx(j,k,l)+yx(j,k,l)+zx(j,k,l)
            sy1 = xy(j,k,l)+yy(j,k,l)+zy(j,k,l)
            sz1 = xz(j,k,l)+yz(j,k,l)+zz(j,k,l)

            dx11 = (x(jp1,kp1,lp1,n) - x(j,k,l,n))
            dy11 = (y(jp1,kp1,lp1,n) - y(j,k,l,n))
            dz11 = (z(jp1,kp1,lp1,n) - z(j,k,l,n))

            a3(j,k,l) = fac2*( sx1*dx11 + sy1*dy11 + sz1*dz11 )
          enddo
        enddo
      enddo
!
!..finite-difference conventions
!
      do ll = 1,lmx
        l  = min(ll,lmx-1)
        lm1 = llr(ll)
        do kk = 1,kmx
          k  = min(kk,kmx-1)
          km1 = kkr(kk)
          do jj = 1,jmx
            j  = min(jj,jmx-1)
            jm1 = jjr(jj)
            vold(n)%arr(jj,kk,ll) = 0.125*( a3(j,k,l) +a3(jm1,km1,lm1) &
                         +a3(jm1,k,l) +a3(j,km1,l) +a3(j,k,lm1) &
                         +a3(jm1,km1,l) +a3(jm1,k,lm1) +a3(j,km1,lm1) )
            volr(n)%arr(jj,kk,ll) = vold(n)%arr(jj,kk,ll)
          enddo
        enddo
      enddo

!
!..check for negative jacobians
!
      nneg = 0
      do l = 1, lmx
      do k = 1, kmx
      do j = 1, jmx
        if(vold(n)%arr(j,k,l).le.0.) nneg = nneg + 1
      enddo
      enddo
      enddo

      if(nneg .ne. 0) then
        write(6,*) nneg, ' negative jacobians in block'
        do l = 1,lmx
        do k = 1,kmx
        do j = 1,jmx
          if( vold(n)%arr(j,k,l).lt.0.0 ) then
            write(6,603) vold(n)%arr(j,k,l), j, k, l
          endif
        enddo
        enddo
        enddo
      endif
    enddo


603 format( ' ',10x,'negative jacobian = ',1p,e10.3,1x,'at j,k,l =', &
                    3i5,5x)

    end subroutine metfv

!********************************************************************
      subroutine finddonor(xp,x,js,ks,ls,jmax,kmax,lmax,maxsearch)
!*********************************************************************
      integer :: js,ks,ls,jmax,kmax,lmax,maxsearch
      real :: x(jmax,kmax,lmax,3)
      real :: xp(3)

      integer :: jj,kk,ll,j,k,l,m,n,nsearch
      real    :: dum
      logical :: Inside,Outside,alter,notoutside
      integer :: movej,movek,movel
      integer :: movejp,movekp,movelp
      integer,allocatable :: i_(:)
      real,allocatable :: xc(:,:)
      real,allocatable :: p(:),q(:),r(:),pqr(:)
      integer :: mo(6),ma(6),mb(6)

      allocate(i_(3))
      allocate(xc(8,3))
      allocate(p(8),q(8),r(8),pqr(2*8)) ! 2**Ndim

      DATA mo(1),mo(2),mo(3),mo(4),mo(5),mo(6) /1,5,5,2,3,6/
      DATA ma(1),ma(2),ma(3),ma(4),ma(5),ma(6) /2,1,6,6,4,5/
      DATA mb(1),mb(2),mb(3),mb(4),mb(5),mb(6) /3,7,1,4,7,8/

!..Get starting cell index (previous cell) j,k,l for the search
        j = min(js,jmax-1)
        k = min(ks,kmax-1)
        l = min(ls,lmax-1)

        Inside = .FALSE.                    !assumed false to enter the loop
        Outside = .FALSE.                   !assumed false initially
        notoutside = .FALSE.                !assumed false for maxsearch
        nsearch = 0
        if (maxsearch.eq.0) maxsearch = 60

        movejp = 0; movekp = 0; movelp = 0
        do while (.not.Inside.and..not.Outside)
          Inside = .TRUE.                   !assumed true intially
          jj = j
          kk = k
          ll = l

          DO m=0,1
            do n = 1,3
              xc(4*m+1,n) = x(jj,  kk,  ll+m,n)
              xc(4*m+2,n) = x(jj+1,kk,  ll+m,n)
              xc(4*m+3,n) = x(jj,  kk+1,ll+m,n)
              xc(4*m+4,n) = x(jj+1,kk+1,ll+m,n)
            enddo
          ENDDO

          movej  = 0; movek  = 0; movel  = 0
!..Cross+dot product pqr=(OPxOQ).OR for in/outside cell test of point xp
          DO m=1,6                       !2*Ndim = number of cell faces
            DO n=1,3
              dum = xc(mo(m),n)
              p(n) = xc(ma(m),n) - dum
              q(n) = xc(mb(m),n) - dum
              r(n) = xp(n) - dum
            ENDDO
            pqr(m) = (p(2)*q(3)-p(3)*q(2))*r(1) &
                  + (p(3)*q(1)-p(1)*q(3))*r(2) &
                  + (p(1)*q(2)-p(2)*q(1))*r(3)
            IF(pqr(m) .lt. 0.) THEN  !If outside, get neighboring cell index
              Inside = .FALSE.
              IF(m .eq. 1) then
                l = l-1
                movel = movel-1
              endif
              IF(m .eq. 2) then
                j = j-1
                movej = movej-1
              endif
              IF(m .eq. 3) then
                k = k-1
                movek = movek-1
              endif
              IF(m .eq. 4) then
                j = j+1
                movej = movej+1
              endif
              IF(m .eq. 5) then
                k = k+1
                movek = movek+1
              endif
              IF(m .eq. 6) then
                l = l+1
                movel = movel+1
              endif
            ENDIF
          ENDDO

          alter = .false.
          if ((j.lt.1.or.j.ge.jmax).and.(movek.ne.0.or.movel.ne.0)) then
            j = j - movej
            movej = 0
            alter = .true.
          endif
          if ((k.lt.1.or.k.ge.kmax).and.(movej.ne.0.or.movel.ne.0)) then
            k = k - movek
            movek = 0
            alter = .true.
          endif
          if ((l.lt.1.or.l.ge.lmax).and.(movej.ne.0.or.movek.ne.0)) then
            l = l - movel
            movel = 0
            alter = .true.
          endif

          if ((movejp+movej).eq.0.and.(movekp+movek).eq.0.and. &
              (movelp+movel).eq.0.and.alter) then
            movej = 0
            movek = 0
            movel = 0
          endif

          if (movej.eq.0.and.movek.eq.0.and.movel.eq.0.and..not.inside) &
            outside = .true.

          movejp = movej
          movekp = movek
          movelp = movel

          if (j.lt.1.or.k.lt.1.or.l.lt.1.or.j.ge.jmax.or.k.ge.kmax.or. &
              l.ge.lmax) Outside = .TRUE.

          if(nsearch .ge. 3) then
            if(3*int((nsearch-1)/3) .eq. nsearch-1) then
              i_(1) = j
              i_(2) = k
              i_(3) = l
            else
              if(i_(1).eq.j .and. i_(2).eq.k .and. i_(3).eq.l) &
                inside = .true.
            endif
          endif
          nsearch = nsearch + 1
          if (nsearch.gt.maxsearch) then
            outside = .true.
            notoutside = .true.
          endif

        enddo

        if (Inside.and..not.outside) then
          js = j
          ks = k
          ls = l
        else
          js = -abs(j-min(movej,0))
          ks = -abs(k-min(movek,0))
          ls = -abs(l-min(movel,0))
          if (notoutside) ks = -ks !make ks positive if stopped by maxsearch
        endif

      end subroutine finddonor

!********************************************************************
      subroutine find_fractions(xp,x,frac,jd,kd,ld,jmax,kmax,lmax)
!..Interpolate coordinates frac in cell jd,kd,ld by using Newton iteration
!********************************************************************
        integer :: jd,kd,ld,jmax,kmax,lmax
        real :: x(jmax,kmax,lmax,3),frac(3)
        real :: xp(3)

        integer :: j,k,l,jj,kk,ll,m,n,niter
        real    :: a2,a3,a4,a5,a6,a7,a8,resid
        integer,allocatable :: i_(:)
        real,allocatable :: xx(:,:,:,:),xc(:,:),xnew(:),B(:,:)

        allocate(i_(3))
        allocate(xx(2,2,2,3),xc(8,3),xnew(3),B(3,4))

!..get coordinates of stencil points for interpolation

        do ll = 1,2
        do kk = 1,2
        do jj = 1,2
          i_(1) = jd - 1 + jj
          i_(2) = kd - 1 + kk
          i_(3) = ld - 1 + ll
          do n=1,3
            xx(jj,kk,ll,n) = x(i_(1),i_(2),i_(3),n)
          enddo
        enddo
        enddo
        enddo

        do m=0,1
          do n = 1,3
            xc(4*m+1,n) = x(jd,  kd,  ld+m,n)
            xc(4*m+2,n) = x(jd+1,kd,  ld+m,n)
            xc(4*m+3,n) = x(jd,  kd+1,ld+m,n)
            xc(4*m+4,n) = x(jd+1,kd+1,ld+m,n)
          enddo
        enddo

!..initial frac guess (linear with no cross terms)

        do n=1,3
          b(n,1) = xc(2,n) - xc(1,n)
          b(n,2) = xc(3,n) - xc(1,n)
          b(n,3) = xc(5,n) - xc(1,n)
          b(n,4) = xp(n) - xc(1,n)
        enddo

        call matrixinv3x3(b)

        do n=1,3
          frac(n) = b(n,4)
          if((frac(n)-2)*(frac(n)+1) .gt. 0.) frac(n) = .5
        enddo

!..second frac guess (linear with cross terms)

        do n=1,3
          a2 = xc(2,n) - xc(1,n)
          a3 = xc(3,n) - xc(1,n)
          a4 = xc(5,n) - xc(1,n)
          a5 = xc(1,n) - xc(2,n) - xc(3,n) + xc(4,n)
          a6 = xc(1,n) - xc(2,n) - xc(5,n) + xc(6,n)
          a7 = xc(1,n) - xc(3,n) - xc(5,n) + xc(7,n)
          a8 = xc(2,n) - xc(1,n) + xc(3,n) - xc(4,n) &
          + xc(5,n) - xc(6,n) - xc(7,n) + xc(8,n)
          b(n,1) = a2 + a5*frac(2) + a6*frac(3) + a8*frac(2)*frac(3)
          b(n,2) = a3 + a5*frac(1) + a7*frac(3) + a8*frac(1)*frac(3)
          b(n,3) = a4 + a6*frac(1) + a7*frac(2) + a8*frac(1)*frac(2)
          b(n,4) = xc(1,n)-xp(n) + a2*frac(1) + a3*frac(2) + a4*frac(3) &
               + a5*frac(1)*frac(2) + a6*frac(1)*frac(3) &
               + a7*frac(2)*frac(3) + a8*frac(1)*frac(2)*frac(3)
        enddo
        call matrixinv3x3(b)
        do n=1,3
          frac(n) = frac(n) - b(n,4)
          if((frac(n)-1.5)*(frac(n)+.5) .gt. 0.) frac(n) = .5
        enddo

        niter = 0
        resid = 1.e10

        do while (niter.lt.10.and.resid.gt.1e-8)

!..solve f(frac) =0 =func(frac)-xp for frac. [dfdxi(frac)]*deltafrac = -f(frac)
          do n=1,3
            do m=1,3
              call trilinearinterp(m,frac,n, xx, b(n,m))
            enddo
            call trilinearinterp(0,frac,n, xx, xnew(n))
            b(n,4) = xnew(n) - xp(n)
          enddo

          call matrixinv3x3(b)

          do n=1,3
            frac(n) = frac(n) - b(n,4)
          enddo

          resid = 0.
          do n=1,3
            resid = resid + (xnew(n)-xp(n))**2
          enddo

          niter = niter + 1
        enddo

!...fix fractions if they don't fall in [0-1] range.

        if((frac(1)+.0001)*(frac(1)-1.0001).gt.0. .or. &
          (frac(2)+.0001)*(frac(2)-1.0001).gt.0. .or. &
          (frac(3)+.0001)*(frac(3)-1.0001).gt.0. ) then
          if(frac(1) .lt. -.0001) frac(1) = 0.01
          if(frac(2) .lt. -.0001) frac(2) = 0.01
          if(frac(3) .lt. -.0001) frac(3) = 0.01
          if(frac(1) .gt. 1.0001) frac(1) = 0.99
          if(frac(2) .gt. 1.0001) frac(2) = 0.99
          if(frac(3) .gt. 1.0001) frac(3) = 0.99
        end if

      end subroutine find_fractions

!********************************************************************
        subroutine trilinearinterp(ideriv,dpsi,n,f,fint)

!!  trilinear interpolation
!!  (returns fint)

!********************************************************************
        real dpsi(3),f(2,2,2,3),fint
        integer j,k,l,m, n,ideriv
        real g(2,3)

        do m=1,3
        if(m .ne. ideriv) then
          g(1,m) = 1 - dpsi(m)
          g(2,m) = dpsi(m)
        else
          g(1,m) = -1
          g(2,m) = 1
        endif
        enddo

        fint = 0.
        do l=1,2
        do k=1,2
        do j=1,2
          fint = fint + f(j,k,l,n) * g(j,1)*g(k,2)*g(l,3)
        enddo
        enddo
        enddo

        end subroutine trilinearinterp

!********************************************************************
        subroutine matrixinv3x3(a)
!
!  cramer's method to solve ax=b
!  (both utilizes & returns a)
!!
!! a(3,4) --- matrix to be inverted
!! 
!! a(:,4) contains the right hand side vector. the final solution is also stored
!! in a(:,4).
!!
!********************************************************************
        real a(3,4),x,y,z,det
!
        det = a(1,1) * (a(2,2)*a(3,3) - a(2,3)*a(3,2)) &
            - a(1,2) * (a(2,1)*a(3,3) - a(2,3)*a(3,1)) &
            + a(1,3) * (a(2,1)*a(3,2) - a(2,2)*a(3,1))
        x = a(1,4)/det
        y = a(2,4)/det
        z = a(3,4)/det
        a(1,4) = x * (a(2,2)*a(3,3) - a(2,3)*a(3,2)) &
               + y * (a(1,3)*a(3,2) - a(1,2)*a(3,3)) &
               + z * (a(1,2)*a(2,3) - a(1,3)*a(2,2))
        a(2,4) = x * (a(2,3)*a(3,1) - a(2,1)*a(3,3)) &
               + y * (a(1,1)*a(3,3) - a(1,3)*a(3,1)) &
               + z * (a(1,3)*a(2,1) - a(1,1)*a(2,3))
        a(3,4) = x * (a(2,1)*a(3,2) - a(2,2)*a(3,1)) &
               + y * (a(1,2)*a(3,1) - a(1,1)*a(3,2)) &
               + z * (a(1,1)*a(2,2) - a(1,2)*a(2,1))
!
        end subroutine matrixinv3x3

!********************************************************************

end module ihc
