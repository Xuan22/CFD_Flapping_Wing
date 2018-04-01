!!! kind_param - Module for controlling precision of primitive types
!!!
!!! For real kinds, defines 3 different kinds: float, double and
!!! rdp. float by default is single precision, double is double
!!! precision and the precision of rdp depends on compile time
!!! directive SINGLE_PRECISION. If SINGLE_PRECISION is defined at
!!! compile time, rdp is equivalent to float, else it is equivalent to
!!! double.
!!!
!!! By default, use real(kind=rdp) for real declarations. However, if
!!! you want more control over the precision, then use
!!! real(kind=double) to have exact precision. See plot3d module for
!!! such an example.

module kind_param
  implicit none

  integer, parameter :: intknd=kind(1)  
  integer, parameter :: float=kind(1.0) 
  integer, parameter :: double=kind(1.0D0)

  ! Math: add IBC
  INTEGER, PARAMETER :: isp = SELECTED_INT_KIND(4), &
                        idp = SELECTED_INT_KIND(8), &
                        rsp = SELECTED_REAL_KIND(6)
#ifdef USE_SINGLE_PRECISION
  integer, parameter :: rdp=float
#else
  integer, parameter :: rdp=double
#endif

end module kind_param

