!----------------------------------------------------------------------
!module constants
!----------------------------------------------------------------------
!use Precision_Def
!implicit none
!----------------------------------------------------------------------

!  real(kind=rdp), parameter :: eps_sp = 10**(-7.d0)
!  real(kind=rdp), parameter :: eps_dp = 10**(-14.d0)
!  real(kind=rdp), parameter :: third  = 1.d0/3.d0


!end module constants

!----------------------------------------------------------------------

!----------------------------------------------------------------------
module immersedBoundVars
!----------------------------------------------------------------------
!use Precision_Def
use kind_param
implicit none
!----------------------------------------------------------------------

  ! input options
  logical               :: ibound, iforce, itestMove,ibc_present
  character(LEN = 80)   :: bodyfile
  integer(kind=idp) :: imove

  ! initial body position
  character(LEN = 80)   :: initial_x,initial_y,initial_z
  
  ! facets of immersed body
  type Facet                             
    real(kind=rdp)   :: vertex(3,3)
    real(kind=rdp)   :: norm(3),center(3)
    real(kind=rdp)   :: area
  end type

  type(Facet), allocatable :: bodyFacet(:)
  integer(kind=idp)    :: nfacet
  real(kind=rdp)      :: totalArea
  real(kind=rdp)      :: xb_min,xb_max,yb_min,yb_max,zb_min,zb_max
  

  ! force calculation parameters
  type FacetVelIdx                             
    integer(kind=idp) :: index
    real(kind=rdp)   :: centVel(3)
  end type

  type(Facet), allocatable         :: forceFacet(:)
  type(FacetVelIdx), allocatable   :: velIndex(:)
  integer(kind=idp)            :: nfacetForce  


  real(kind=rdp), allocatable :: bodyVarsAll(:,:)
  logical                          :: bodySurfaceVar

  ! boundary nodes
  type BoundaryNode                             
    integer(kind=idp) :: i,j,k
    integer(kind=idp) :: iint,jint,kint
    integer(kind=idp) :: normFacet
    real(kind=rdp)   :: d2surf,d2surf_ratio,d2int_ratio
    real(kind=rdp)   :: norm(3),surf(3)
    real(kind=rdp)   :: vel
    real(kind=rdp)   :: interp(8)
  end type
  
  type(BoundaryNode), allocatable :: boundNode(:)
  integer(kind=idp)           :: nbound 

  ! Math: iblank array for ibc mesh
  integer(kind=idp),allocatable :: iblk_ibc(:,:,:)
  ! Math: flag for parallel IBC
  integer(kind=idp) :: parallel_ibc

end module immersedBoundVars
