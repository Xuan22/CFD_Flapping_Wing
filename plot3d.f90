!!! plot3d - Plot 3D processing utilities
!!! 
!!! Provides utilities for reading, writing p3d files, as well as
!!! extracting surface details from a given grid.
!!! 
!!! Defines 4 new types to handle grid and q for single and double
!!! precision file formats:
!!!     p3d  - plot3d grid type (double precision)
!!!     p3dq - plot3d solution type (double precision)
!!!     p3dsingle - plot3d grid type (single precision)
!!!     p3dqsp - plot3d solution type (single precision)
!!!
!!! There are separate subroutines defined for each type. However,
!!! they are all wrapped up in interface definitions, so that the same
!!! subroutine names work for all precision types.
!!!
!!! Limitations: Currently handles only a single block
!!! definition. Needs to be extended to handle multi-block grid types.
!!!

!!! BP_EDIT When running on deepthought2 (SLURM system) OVERTURNS reading and writing grid process
!!! was very slow due to implied do loops So they were rewritten to manually collapse them
!!! This may or may not also be an issue for speed in other parts of the code, but this is where it was most egregious

module plot3d
  use kind_param
  use ioutils

  implicit none

  integer, parameter :: QDIM=5  !Variables in solution

  type p3d
     integer :: jm,km,lm
     real(kind=double), dimension(:,:,:), allocatable :: x,y,z
  end type p3d

  type p3dq
     integer :: jm,km,lm
     real(kind=double) :: fstip, alf, reypr, totime
     real(kind=double), dimension(:,:,:,:), allocatable :: qsol
  end type p3dq

  type p3dsingle
     integer :: jm,km,lm
     real(kind=float), dimension(:,:,:), allocatable :: x,y,z
  end type p3dsingle

  type p3dqsp
     integer :: jm,km,lm
     real(kind=float) :: fstip, alf, reypr, totime
     real(kind=float), dimension(:,:,:,:), allocatable :: qsol
  end type p3dqsp

  !Begin wrapper function definitions

  !Allocate/free memory for P3D grid type
  interface allocate_grid_arrays
     module procedure allo_grd_arr, allo_grd_arr_single
  end interface
  interface free_grid_arrays
     module procedure free_grd_arr, free_grd_arr_single
  end interface

  ! read_p3d(grid,filename) - Read a p3d file from disk
  ! write_p3d(grid,filename) - Write out a p3d file
  ! read_p3d_sol(soln,filename) - Read a q file from disk
  ! write_p3d_sol(soln,filename) - Write out a q file
  !
  !    grid - p3d or p3dsingle
  !    soln - p3dq or p3dqsp
  !    Filename - string denoting the filename on disk
  !
  ! Aborts if the file is not found
  interface read_p3d
     module procedure read_p3d_double, read_p3d_single
  end interface
  interface write_p3d
     module procedure write_p3d_double, write_p3d_single
  end interface
  interface read_p3d_sol
     module procedure read_p3dsol_double, read_p3dsol_single
  end interface
  interface write_p3d_sol
     module procedure write_p3dsol_double, write_p3dsol_single
  end interface

  ! extract_surface(grid1, grid2, js, je, ks, ke, ls, le)
  ! extract_solution(qsol1,qsol2, js, je, ks, ke, ls, le)
  ! 
  ! Given a range of j,k,l values, extract the surface or the chunk
  !specified by those ranges
  !
  !  grid1, grid2 - objects of p3d or p3dsingle type
  !  qsol1, qsol2 - objects of p3dq or p3dqsp type
  !  js, je - limits in J direction (wraparound)
  !  ks, ke - limits in K direction (spanwise)
  !  ls, le - limits in L direction (normal)
  interface extract_surface
     module procedure extract_surface_double, extract_surface_single
  end interface
  interface extract_solution
     module procedure extract_solution_double, extract_solution_single
  end interface

contains
  subroutine read_p3d_double(grid,flnam)
    type(p3d), intent(out) :: grid
    character(len=*), optional, intent(in) :: flnam

    integer :: j,k,l,jm,km,lm, un
    character(len=30) :: surffile

    if (present(flnam)) then
       surffile=flnam
    else
       surffile='surf.p3d'
    end if
    
    un=open_file(surffile,form='unformatted',status='old')
    read(un) jm,km,lm
    
    grid%jm=jm
    grid%km=km
    grid%lm=lm
    
    allocate(grid%x(jm,km,lm),grid%y(jm,km,lm),grid%z(jm,km,lm))
! BP_EDIT deepthought2 hada problem with speed manually collapsed implied do loop
    read(un) grid%x,grid%y,grid%z
!   read(un) (((grid%x(j,k,l),j=1,jm),k=1,km),l=1,lm),&
 !        &(((grid%y(j,k,l),j=1,jm),k=1,km),l=1,lm),&
  !       &(((grid%z(j,k,l),j=1,jm),k=1,km),l=1,lm)

    close(un)
  end subroutine read_p3d_double

  subroutine read_p3d_single(grid,flnam)
    type(p3dsingle), intent(out) :: grid
    character(len=*), optional, intent(in) :: flnam

    integer :: j,k,l,jm,km,lm, un
    character(len=30) :: surffile

    if (present(flnam)) then
       surffile=flnam
    else
       surffile='surf.p3d'
    end if
    
    un=open_file(surffile,form='unformatted',status='old')
    read(un) jm,km,lm
    
    grid%jm=jm
    grid%km=km
    grid%lm=lm
    
    allocate(grid%x(jm,km,lm),grid%y(jm,km,lm),grid%z(jm,km,lm))
! BP_EDIT deepthought2 hada problem with speed manually collapsed implied do loop
     read(un) grid%x,grid%y,grid%z

   ! read(un) (((grid%x(j,k,l),j=1,jm),k=1,km),l=1,lm),&
    !     &(((grid%y(j,k,l),j=1,jm),k=1,km),l=1,lm),&
     !    &(((grid%z(j,k,l),j=1,jm),k=1,km),l=1,lm)

    close(un)
  end subroutine read_p3d_single

  subroutine extract_surface_double(grd1,surf,js,je,ks,ke,ls,le)
    type(p3d), intent(in) :: grd1
    type(p3d), intent(out) :: surf
    integer, intent(in) :: js,je,ks,ke,ls,le

    integer :: j, k, l, j1, k1, l1

    surf%jm=je-js+1
    surf%km=ke-ks+1
    surf%lm=le-ls+1
    
    allocate(surf%x(surf%jm,surf%km,surf%lm))
    allocate(surf%y(surf%jm,surf%km,surf%lm))
    allocate(surf%z(surf%jm,surf%km,surf%lm))
    
    do l=ls,le
       l1=l-ls+1
       do k=ks,ke
          k1=k-ks+1
          do j=js,je
             j1=j-js+1
             surf%x(j1,k1,l1)=grd1%x(j,k,l)
             surf%y(j1,k1,l1)=grd1%y(j,k,l)
             surf%z(j1,k1,l1)=grd1%z(j,k,l)
          end do
       end do
    end do
  end subroutine extract_surface_double
  
  subroutine extract_surface_single(grd1,surf,js,je,ks,ke,ls,le)
    type(p3dsingle), intent(in) :: grd1
    type(p3dsingle), intent(out) :: surf
    integer, intent(in) :: js,je,ks,ke,ls,le

    integer :: j, k, l, j1, k1, l1

    surf%jm=je-js+1
    surf%km=ke-ks+1
    surf%lm=le-ls+1
    
    allocate(surf%x(surf%jm,surf%km,surf%lm))
    allocate(surf%y(surf%jm,surf%km,surf%lm))
    allocate(surf%z(surf%jm,surf%km,surf%lm))
    
    do l=ls,le
       l1=l-ls+1
       do k=ks,ke
          k1=k-ks+1
          do j=js,je
             j1=j-js+1
             surf%x(j1,k1,l1)=grd1%x(j,k,l)
             surf%y(j1,k1,l1)=grd1%y(j,k,l)
             surf%z(j1,k1,l1)=grd1%z(j,k,l)
          end do
       end do
    end do
  end subroutine extract_surface_single
  
  subroutine write_p3d_double(grd,flnam)
    type(p3d), intent(in) :: grd
    character(len=*), intent(in), optional :: flnam
    
    integer :: j, k, l, un
    character(len=30) :: grdfile
    
    if (present(flnam)) then
       grdfile=flnam
    else
       grdfile='out.p3d'
    end if

    un=open_file(grdfile,form='unformatted')
    
    write(un) grd%jm,grd%km,grd%lm
! BP_EDIT deepthought2 hada problem with speed manually collapsed implied do loop 	
   write(un) grd%x,grd%y,grd%z
!  write(un) (((grd%x(j,k,l),j=1,grd%jm),k=1,grd%km),l=1,grd%lm),&
 !        &(((grd%y(j,k,l),j=1,grd%jm),k=1,grd%km),l=1,grd%lm),&
  !       &(((grd%z(j,k,l),j=1,grd%jm),k=1,grd%km),l=1,grd%lm)

    close(un)
  end subroutine write_p3d_double

  subroutine write_p3d_single(grd,flnam)
    type(p3dsingle), intent(in) :: grd
    character(len=*), intent(in), optional :: flnam
    
    integer :: j, k, l, un
    character(len=30) :: grdfile
    
    if (present(flnam)) then
       grdfile=flnam
    else
       grdfile='out.p3d'
    end if

    un=open_file(grdfile,form='unformatted')
    
    write(un) grd%jm,grd%km,grd%lm
! BP_EDIT deepthought2 hada problem with speed manually collapsed implied do loop
     write(un) grd%x,grd%y,grd%z
!    write(un) (((grd%x(j,k,l),j=1,grd%jm),k=1,grd%km),l=1,grd%lm),&
 !        &(((grd%y(j,k,l),j=1,grd%jm),k=1,grd%km),l=1,grd%lm),&
  !       &(((grd%z(j,k,l),j=1,grd%jm),k=1,grd%km),l=1,grd%lm)

    close(un)
  end subroutine write_p3d_single

  subroutine allo_grd_arr(grd)
    type(p3d), intent(inout) :: grd

    allocate(grd%x(grd%jm,grd%km,grd%lm))
    allocate(grd%y(grd%jm,grd%km,grd%lm))
    allocate(grd%z(grd%jm,grd%km,grd%lm))
  end subroutine allo_grd_arr
  subroutine allo_grd_arr_single(grd)
    type(p3dsingle), intent(inout) :: grd

    allocate(grd%x(grd%jm,grd%km,grd%lm))
    allocate(grd%y(grd%jm,grd%km,grd%lm))
    allocate(grd%z(grd%jm,grd%km,grd%lm))
  end subroutine allo_grd_arr_single
  
  subroutine free_grd_arr(grd)
    type(p3d), intent(inout) :: grd
    
    deallocate(grd%x,grd%y,grd%z)
  end subroutine free_grd_arr
  subroutine free_grd_arr_single(grd)
    type(p3dsingle), intent(inout) :: grd
    
    deallocate(grd%x,grd%y,grd%z)
  end subroutine free_grd_arr_single

  subroutine grid_double_to_single(grd1,grd2)
    type(p3d), intent(in) :: grd1
    type(p3dsingle), intent(out) :: grd2

    integer :: j, k, l

    grd2%jm=grd1%jm
    grd2%km=grd1%km
    grd2%lm=grd1%lm

    call allocate_grid_arrays(grd2)
    do l=1,grd1%lm
       do k=1,grd1%km
          do j=1,grd1%jm
             grd2%x(j,k,l)=grd1%x(j,k,l)
             grd2%y(j,k,l)=grd1%y(j,k,l)
             grd2%z(j,k,l)=grd1%z(j,k,l)
          end do
       end do
    end do
  end subroutine grid_double_to_single

  subroutine read_p3dsol_double(qvals,flnam)
    type(p3dq), intent(out) :: qvals
    character(len=*), optional, intent(in) :: flnam

    integer :: j,k,l,jm,km,lm,n, un
    character(len=50) :: qfile
    
    if (present(flnam)) then
       qfile=trim(flnam)
    else
       qfile='q.p3d'
    end if
    
    un=get_free_unit()
    open(un,file=qfile,form='unformatted',status='old',err=100)
    print *, 'Reading solution file: '//qfile
    read(un) jm,km,lm
    
    qvals%jm=jm
    qvals%km=km
    qvals%lm=lm
    
    allocate(qvals%qsol(jm,km,lm,QDIM))

    read(un) qvals%fstip, qvals%alf, qvals%reypr, qvals%totime
	! BP_EDIT deepthought2 hada problem with speed manually collapsed implied do loop
     read(un) qvals%qsol	
    ! read(un) ((((qvals%qsol(j,k,l,n),j=1,jm),k=1,km),l=1,lm),n=1,QDIM)
	close(un)
    
    return
100 print *, 'Cannot open file for reading: '// qfile
    stop
  end subroutine read_p3dsol_double

  subroutine read_p3dsol_single(qvals,flnam)
    type(p3dqsp), intent(out) :: qvals
    character(len=*), optional, intent(in) :: flnam

    integer :: j,k,l,jm,km,lm,n, un
    character(len=50) :: qfile
    
    if (present(flnam)) then
       qfile=trim(flnam)
    else
       qfile='q.p3d'
    end if
    
    un=get_free_unit()
    open(un,file=qfile,form='unformatted',status='old',err=100)
    print *, 'Reading solution file: '//qfile
    read(un) jm,km,lm
    
    qvals%jm=jm
    qvals%km=km
    qvals%lm=lm
    
    allocate(qvals%qsol(jm,km,lm,QDIM))

    read(un) qvals%fstip, qvals%alf, qvals%reypr, qvals%totime
    ! BP_EDIT deepthought2 hada problem with speed manually collapsed implied do loop
    read(un) qvals%qsol
	! read(un) ((((qvals%qsol(j,k,l,n),j=1,jm),k=1,km),l=1,lm),n=1,QDIM)
    close(un)
    
    return
100 print *, 'Cannot open file for reading: '//qfile
    stop
  end subroutine read_p3dsol_single

  subroutine write_p3dsol_double(qvals,flnam)
    type(p3dq), intent(in) :: qvals
    character(len=*), optional, intent(in) :: flnam

    integer :: j,k,l,jm,km,lm,n, un
    character(len=50) :: qfile
    
    if (present(flnam)) then
       qfile=trim(flnam)
    else
       qfile='q.p3d'
    end if
    
    un=get_free_unit()
    open(un,file=qfile,form='unformatted',err=100)
    write(un) qvals%jm,qvals%km,qvals%lm
    
    write(un) qvals%fstip, qvals%alf, qvals%reypr, qvals%totime
! BP_EDIT deepthought2 hada problem with speed manually collapsed implied do loop	
     write(un) qvals%qsol
    ! write(un) ((((qvals%qsol(j,k,l,n),j=1,qvals%jm),&
    !     &k=1,qvals%km),l=1,qvals%lm),n=1,QDIM)  
	close(un)
    
    print *, 'Double precision soln output to: '//qfile
    return
100 print *, 'Cannot open file for writing: '//qfile
    stop
  end subroutine write_p3dsol_double

  subroutine write_p3dsol_single(qvals,flnam)
    type(p3dqsp), intent(in) :: qvals
    character(len=*), optional, intent(in) :: flnam

    integer :: j,k,l,jm,km,lm,n, un
    character(len=50) :: qfile
    
    if (present(flnam)) then
       qfile=trim(flnam)
    else
       qfile='q.p3d'
    end if
    
    un=get_free_unit()
    open(un,file=qfile,form='unformatted',err=100)
    write(un) qvals%jm,qvals%km,qvals%lm
    

    write(un) qvals%fstip, qvals%alf, qvals%reypr, qvals%totime
! BP_EDIT deepthought2 hada problem with speed manually collapsed implied do loop
    write(un) qvals%qsol
  !  write(un) ((((qvals%qsol(j,k,l,n),j=1,qvals%jm),&
   !      &k=1,qvals%km),l=1,qvals%lm),n=1,QDIM)
	close(un)
    
    print *, 'Single precision soln output to: '//qfile
    return
100 print *, 'Cannot open file for writing: '//qfile
    stop
  end subroutine write_p3dsol_single

  subroutine extract_solution_double(qval1,qval2,js,je,ks,ke,ls,le)
    type(p3dq), intent(in) :: qval1
    type(p3dq), intent(out) :: qval2
    integer, intent(in) :: js,je,ks,ke,ls,le

    integer :: j, k, l, j1, k1, l1, i

    qval2%jm=je-js+1
    qval2%km=ke-ks+1
    qval2%lm=le-ls+1
    
    allocate(qval2%qsol(qval2%jm,qval2%km,qval2%lm,QDIM))

    qval2%fstip=qval1%fstip
    qval2%alf=qval1%alf
    qval2%reypr=qval1%reypr
    qval2%totime=qval1%totime

    do l=ls,le
       l1=l-ls+1
       do k=ks,ke
          k1=k-ks+1
          do j=js,je
             j1=j-js+1
             do i=1,QDIM
                qval2%qsol(j1,k1,l1,i)=qval1%qsol(j,k,l,i)
             end do
          end do
       end do
    end do
  end subroutine extract_solution_double

  subroutine extract_solution_single(qval1,qval2,js,je,ks,ke,ls,le)
    type(p3dqsp), intent(in) :: qval1
    type(p3dqsp), intent(out) :: qval2
    integer, intent(in) :: js,je,ks,ke,ls,le

    integer :: j, k, l, j1, k1, l1, i

    qval2%jm=je-js+1
    qval2%km=ke-ks+1
    qval2%lm=le-ls+1
    
    allocate(qval2%qsol(qval2%jm,qval2%km,qval2%lm,QDIM))

    qval2%fstip=qval1%fstip
    qval2%alf=qval1%alf
    qval2%reypr=qval1%reypr
    qval2%totime=qval1%totime

    do l=ls,le
       l1=l-ls+1
       do k=ks,ke
          k1=k-ks+1
          do j=js,je
             j1=j-js+1
             do i=1,QDIM
                qval2%qsol(j1,k1,l1,i)=qval1%qsol(j,k,l,i)
             end do
          end do
       end do
    end do
  end subroutine extract_solution_single

  subroutine soln_double_to_single(grd1,grd2)
    type(p3dq), intent(in) :: grd1
    type(p3dqsp), intent(out) :: grd2

    integer :: j, k, l,i

    grd2%jm=grd1%jm
    grd2%km=grd1%km
    grd2%lm=grd1%lm

    allocate(grd2%qsol(grd2%jm,grd2%km,grd2%lm,QDIM))
    do l=1,grd1%lm
       do k=1,grd1%km
          do j=1,grd1%jm
             do i=1,QDIM
                grd2%qsol(j,k,l,i)=grd1%qsol(j,k,l,i)
             end do
          end do
       end do
    end do
  end subroutine soln_double_to_single

end module plot3d

