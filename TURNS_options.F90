module TURNS_options
   implicit none

   !! Directions
   integer, parameter :: JDIR=1
   integer, parameter :: KDIR=2
   integer, parameter :: LDIR=3
   
   !! Boundary condition 
   integer, parameter :: BC_INVISCWALL=4
   integer, parameter :: BC_WALL=5
   integer, parameter :: BC_GROUND=6
   integer, parameter :: BC_EXTRAPOLATE=10
   integer, parameter :: BC_SYMMETRIC=11
	
   !Math: periodic BC
   integer, parameter :: BC_PERIODIC=12
   integer, parameter :: BC_AVERAGE=14
   integer, parameter :: BC_AXISYM=22
   integer, parameter :: BC_FARFIELD=47
   integer, parameter :: BC_HOVER=48
   integer, parameter :: BC_WAKEAVG=51
   integer, parameter :: BC_INTERFACE=87
   integer, parameter :: BC_INTERNAL=101

   !! Turbulence options
   integer, parameter :: BALDWIN_LOMAX=1
   integer, parameter :: SPALART_ALLMARAS=2

   !! Spatial order of accuracy
   integer, parameter :: RHS_FIRST=-1
   integer, parameter :: RHS_SECOND=-2
   integer, parameter :: RHS_THRID=-3

   !! Type of limiting
   integer, parameter :: LIM_ALL=1
   integer, parameter :: LIM_JDIR=-1
   integer, parameter :: LIM_JTIP=0
   integer, parameter :: LIM_NONE=2

   !! Field velocity wake evaluations
   integer, parameter :: FV_BEDDOES=1
   integer, parameter :: FV_BIOTSAVART=2
   integer, parameter :: FV_FAST=3
   integer, parameter :: FV_USER=4

   !! Deformation
   integer, parameter :: NO_DEFORM=0
   integer, parameter :: DEFORMING=1

   !! Type of unsteady flow
   integer, parameter :: FL_STEADY=0
   integer, parameter :: FL_UNSTEADY=1
   integer, parameter :: FL_HOVER=2
   integer, parameter :: FL_QUASI=3
   ! Math: add pitch/plunge kinematics
   integer, parameter :: FL_PITCH=5
   integer, parameter :: FL_FLAP=6

   !! Time accuracy
   integer, parameter :: NTAC_FIRST=1
   integer, parameter :: NTAC_SECOND=2

   !! Computational domain types
   integer, parameter :: DOM_WAKE=0
   integer, parameter :: DOM_BLADE=1
   integer, parameter :: DOM_SLAT=2
   integer, parameter :: DOM_GROUNDWAKE=3

   ! Math: add IBC
   integer, parameter :: DOM_IBC=4 ! box mesh around fuselage for IBC

   ! Math: adding domain for fuselage mesh (not IBC)
   integer, parameter :: DOM_FUS=5

   !! Blade deformation types
   integer, parameter :: DEFL_TDU=0
   integer, parameter :: DEFL_EULER=1
   integer, parameter :: DEFL_RCAS=2

   !! Connectivity types
   integer, parameter :: CONN_JAINA=0
   integer, parameter :: CONN_IHC=1

   !! Connectivity modes
   integer, parameter :: CONN_STATIC=0
   integer, parameter :: CONN_DYNAMIC=1

   !! IHC mesh type options
   integer, parameter :: IHC_NO_WALL=1
   integer, parameter :: IHC_OH=2
   integer, parameter :: IHC_CH=3
   integer, parameter :: IHC_1WALL=4
   integer, parameter :: IHC_CO=5
   integer, parameter :: IHC_FARFIELD=-1

   !! Mesh splitting strategies
   integer, parameter :: SPLIT_AUTO=0
   integer, parameter :: SPLIT_MANUAL_3D=1
end module TURNS_options


