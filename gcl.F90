
module gcl

#ifdef USEIO
    use io_base
    use cplot3d
#endif

    implicit none

    integer :: nj, nk, nl
    real(kind=8), allocatable, dimension(:,:,:) :: cvol, cvol0
    real(kind=8), allocatable, dimension(:,:,:) :: svol
    real(kind=8), allocatable, dimension(:,:,:) :: res
    real(kind=8), allocatable, dimension(:) :: tvol

contains

subroutine init(jmax, kmax, lmax)
    integer, intent(in) :: jmax, kmax, lmax

    integer :: mdim

    nj = jmax
    nk = kmax
    nl = lmax

    mdim = jmax
    if(kmax > jmax) mdim = kmax
    if(lmax > kmax) mdim = lmax

    allocate(cvol(nj, nk, nl))
    allocate(cvol0(nj, nk, nl))
    allocate(svol(nj, nk, nl))
    allocate(tvol(mdim))

    cvol = 0.0
    cvol0 = 0.0
    svol = 0.0

end subroutine init

subroutine gcl_delete()

    deallocate(cvol)
    deallocate(cvol0)
    deallocate(svol)
    if(allocated(res)) deallocate(res)

end subroutine gcl_delete

subroutine set_jacobians(jac)
    real(kind=8), dimension(nj,nk,nl), intent(in) :: jac

    !local
    integer :: j, k, l

    forall(j=1:nj, k=1:nk, l=1:nl)

        svol(j,k,l) = 0.0

        cvol0(j,k,l) = cvol(j,k,l)
        cvol(j,k,l) = 1.0/jac(j,k,l)

    end forall

end subroutine set_jacobians

subroutine set_cell_volumes(cellvol)
    real(kind=8), dimension(nj,nk,nl), intent(in) :: cellvol

    !local
    integer :: j, k, l

    forall(j=1:nj, k=1:nk, l=1:nl)

        svol(j,k,l) = 0.0

        cvol0(j,k,l) = cvol(j,k,l)
        cvol(j,k,l) = cellvol(j,k,l)

    end forall

end subroutine set_cell_volumes

subroutine calc_residual(l2norm, rmax, jc, kc, lc)
    real(kind=8), intent(out) :: l2norm, rmax
    integer, intent(out) :: jc, kc, lc

    ! local
    integer :: j, k, l
    real(kind=8) :: temp

    if(.not.allocated(res)) allocate(res(nj, nk, nl))

    res = 0.0
    forall(j=2:nj-1, k=2:nk-1, l=2:nl-1)
        res(j,k,l) = cvol(j,k,l) - cvol0(j,k,l) - svol(j,k,l)
        res(j,k,l) = res(j,k,l)/cvol(j,k,l)
    end forall

    l2norm = 0.0

    rmax = abs(res(2,2,2))
    jc = 2
    kc = 2
    lc = 2

    do l = 2, nl-1
        do k = 2, nk-1
            do j = 2, nj-1

                temp = abs(res(j,k,l))

                l2norm = l2norm + temp*temp

                if(temp > rmax)then
                    rmax = temp
                    jc = j
                    kc = k
                    lc = l
                end if

            end do
        end do
    end do

    l2norm = sqrt(l2norm/(nj*nk*nl))

end subroutine calc_residual

#ifdef USEIO
subroutine write_resfield(filename)
    character(len=128), intent(in) :: filename

    !local
    integer :: un
    
    un = open_binary_write(filename)
    call write_func(un, res, nj, nk, nl)
    close(un)

end subroutine write_resfield
#endif

end module gcl
