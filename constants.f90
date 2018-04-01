!!! constants - Define useful constants 
!!! 
!!! Defines often used constants for the code. Precision of the
!!! constants is set to kind rdp, so that it can be controlled during
!!! compile time.

module constants
  use kind_param
  implicit none
  save
  
  real(kind=rdp), parameter :: m_pi=3.14159265358979323846_rdp
  
  real(kind=rdp), parameter :: zero=0.0_rdp
  real(kind=rdp), parameter :: one=1.0_rdp
  real(kind=rdp), parameter :: two=2.0_rdp
  real(kind=rdp), parameter :: three=3.0_rdp
  real(kind=rdp), parameter :: four=4.0_rdp
  real(kind=rdp), parameter :: five=5.0_rdp
  real(kind=rdp), parameter :: six=6.0_rdp
  real(kind=rdp), parameter :: seven=7.0_rdp
  real(kind=rdp), parameter :: eight=8.0_rdp
  real(kind=rdp), parameter :: nine=9.0_rdp

  real(kind=rdp), parameter :: half=0.5_rdp
  real(kind=rdp), parameter :: third=one/three
  real(kind=rdp), parameter :: fourth=0.25_rdp
  real(kind=rdp), parameter :: eighth=0.125_rdp

  real(kind=rdp), parameter :: twopi=two*m_pi
  real(kind=rdp), parameter :: piInv=one/m_pi
  real(kind=rdp), parameter :: d2r=m_pi/180.0_rdp
  real(kind=rdp), parameter :: r2d=one/d2r

  complex(kind=rdp), parameter :: icmplx=cmplx(zero,one)
  complex(kind=rdp), parameter :: czero=cmplx(zero,zero)

  ! Special characters
  character(len=1), parameter :: cTab=achar(9)
  character(len=1), parameter :: cSpace=achar(32)
  character(len=1), parameter :: cNewline=achar(10)
  character(len=1), parameter :: cReturn=achar(13)

end module constants

