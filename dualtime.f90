module dualtime

    use ioutils
    use circlist_class
    use params_global, only: dtpseudo

    implicit none

    type(circlist), private :: dtpseudo_range

contains

subroutine dualtime_init()
    !DESC:  Loads list of pseudo times from ASCII input file into
    !       dtpseudo_range. Currently assuming input file name
    !       is dtpseudo.inp.

    !local var
    logical :: inuse
    integer :: n, un
    character(len=2) :: key
    character(len=*), parameter :: ifile = "dtpseudo.inp"
    real :: val

    ! initialize empty circlist...
    call new(dtpseudo_range)

    un=open_file(ifile,status='old')
    do
        read(un, '(a2)', advance='no', eor=99, end=100) key
        !print*, "key = "//key
        if(key == ">>")then
            read(un, *, end=100, err=100) val
            !print*, "dtpseudo val = ", val
            call append(dtpseudo_range, val)
        end if
99      continue
    end do
100 close(un)

    call reset(dtpseudo_range)

end subroutine dualtime_init

subroutine dualtime_del()
    call del(dtpseudo_range)
end subroutine dualtime_del

subroutine dualtime_reset_sequence()
    dtpseudo = first(dtpseudo_range)
end subroutine dualtime_reset_sequence

subroutine dualtime_next_dtpseudo()
    dtpseudo = next(dtpseudo_range)
end subroutine dualtime_next_dtpseudo

subroutine dualtime_timescale(oat, h, q, iblank, tscale, bt, iprecon, jmax, kmax, lmax)
      ! DESC: Computes local time scaling for dual time stepping
      ! INPUTS:
      !   oat... Factor dependent on physical time difference operator
      !           e.g. oat = 2./3. for 2nd order backward
      !                oat = 1.    for Euler implicit
      !   h.......... Physical (non-dim) time step size
      !   dtpseudo... Pseudo time step size
      !   q.......... Conserved flow var
      !   iblank..... I-blank array 
      !   jmax, kmax, lmax... dimensions of current grid
      !
      ! OUTPUT:
      !   tscale... array of local time scale factors
        integer, intent(in) :: jmax, kmax, lmax
        real, intent(in) :: oat, h
        integer, dimension(jmax,kmax,lmax), intent(in) :: iblank
        real, dimension(jmax, kmax, lmax,6), intent(in) :: q
        real, dimension(jmax, kmax, lmax), intent(out) :: tscale, bt
          
        !local var
        integer :: j, k, l
        real :: s

        logical :: iprecon

        !For ease of diffing output..
        !write(STDOUT,*) "[in dualtime_timescale] dtpseudo = ", dtpseudo

        ! Reset indices from 1 to MAX for consistency with Vinod's code.
        do l = 1, lmax !- 1
            do k = 1, kmax !- 1
                do j = 1, jmax !- 1
                    s = ( 1.0 + 0.002*sqrt(q(j,k,l,6)) )
                    s = s/( 1.0 + sqrt(q(j,k,l,6)) )

                    if (iprecon) then !vinod...
                     s = 80.*( 1.0 + 0.01*sqrt(q(j,k,l,6)))
                     s = s/( 1. + sqrt(q(j,k,l,6)))
                     if (s.gt.10.) s = 10.
                    endif
                    s = s*dtpseudo

                    ! Changed term in denominator slightly  
                    ! Set timescales to zero if fringe/hole point
                    tscale(j,k,l) = max(iblank(j,k,l),0)*s/(1. + s/h/oat)
                    if (iprecon) then
                      bt(j,k,l) = max(iblank(j,k,l),0)*1.0/(1. + s/h/oat)
                    endif
                end do
            end do
        end do

end subroutine dualtime_timescale

end module dualtime
