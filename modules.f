C*******************************************************************
c..   global parameters module
C*******************************************************************
      module params_global

      use TURNS_options
      use ioutils
      
      implicit none

C..   number of meshes and the current mesh under process
      integer :: nmesh,iwm

c...  offset between blades
      real    :: dblad

C..Parameters for outputting acoustic information
C..NACOU and NACOUP set by input
      integer :: NACOU,NACOUP,NAKP1,NAKP2,NAKP3,
     &     NAKP4,NAKP5,NWW,NWWQ,NACOUQ ! acoustic output

C..Parameters that define critical points for C-H mesh
C..JTAIL1,KTIP,HALF set by input
C..KROOT always set to 1
C..JLE,JTAIL2,JBOT calculated based on input in INITIA
      integer :: JLE,JTAIL1, JTAIL2, JBOT, KROOT, KTIP, HALF ! c-mesh

      integer :: ksym

C..Parameters for space and time
C..DX1,DY1,DZ1 always set to 1
C..H,HD calculated based on input in INITIA
      real :: DX1, DY1, DZ1, H, HD !dels

C..Parameters that define the flow
C..ALF,FSMACH,REY set by input
C..GAMMA,GM1,GGM1,PR,TINF set to fixed values in INITIA
      real :: ALF, GM1, GAMMA, GGM1, FSMACH, PR, REY, TINF !floprm
      real :: ALFA

C..Free stream Parameters
      real :: EINF,HTINF,PINF,RINF, UINF, VINF, WINF,AINF !frstrm

C..Parameters that define perturbation flow
      integer :: IPERT,INITPE
      real    :: XVOR, ZVOR, RCORE, VORGAM !perprm

C..Parameters that define grid spatial dimensions and time steps
C..JMAX,KMAX,LMAX,NSTEPS set by input
C..ISTEP0 is initial time step read in with Q-file
C..JM,KM,LM calculated based on input in INITIA
      integer :: JM, JMAX, KM, KMAX, LM, LMAX, ISTEP0, 
     &     NSTEPS,MDIM,NV,ND          !grdprm

C.. Dimension variables
C
C  Parameters for type of limiting and pressure bc
C..ILIM set by input
C..IBCWP set to 1 in INITIA
      integer :: ilim,ibcwp

c..Parameters for preconditioning
C..IPRECON,MP set by input
      logical :: iprecon
      real    :: MP

C..Parameters for symmetry in k-direction (not used)
C..KSYM set to 0 in INITIA
C..KKR and KKP set in INITIA
      integer :: ksymm

C..Parameters used for grid coarsening
C..JINT,KINT,LINT set by input
      integer :: jint,kint,lint

C..Parameters used for spatial order of accuracy
C..IRHSY set by input
      integer :: irhsy
C
C  Parameters for describing rotor conditions
C..FMTIP,RARTIO and TOTIME set by input
C..RF calculated based on input in INITIA
      real :: fmtip,rartio,totime,rf,angmax,alpha_dot,totime0

C..Parameters for restart, storing out restart and writing residual
C..IREAD,NREST and NPNORM set by input
      integer :: iread,nrest,npnorm,nmovie,imovie
      integer :: restart_counter

C..Parameters for damping left-hand-side
C..EPSE set by input
      real :: epse

C..Parameter for left-hand-side inversion
C..ILHS set by input
      integer :: ilhs

C..Parameters for time 
C..DT,TIMEAC,IUNST,NTAC,ITNMAX set by input
C..CNBR set to 1 in INITIA (not used)
C..ISTEP is current time step in the run
      real :: cnbr, dt, timeac, dtpseudo
      integer :: iunst,itnmax,ntac,istep,itn,indicial,idual

C..Parameters for viscous flow
C..INVISC,LAMIN set by input
C..RMUE set in INITIA
      real :: rmue
      logical :: invisc,lamin,ithin
      integer iturb

C..Parameters for ground plane simulation
C..ZGROUND set by input
      real :: zground

C..GCL option flag 
C..If true will use GCL satisfying algorithms -- more expensive
C..If false will use free-stream subtraction -- less expensive
      logical :: usegcl

C..Test GCL option flag
C..If true all mesh bc's set to free-stream.
      logical :: testgcl

C..Parameters for wake capturing
C..NBLADE,CTINP,IWAKE set by input
C..BLANG calulated based on input in INITIA
      real :: blang,ctinp
      integer :: iwake,nblade, nblade_tot
      
C.. parameters determining wake coupling
      integer :: imth,iconfine,irefine

C.. parameters defining deformation and deformation scheme
      integer :: ideform,idefl,irotate, idecay
      
c.. controls
      real :: theta0,theta1c,theta1s,xlam

C.. coning
      real :: beta0,beta1c,beta1s,eflap

      real :: beta1

C..   PI
      real :: pi

CC..   MPI stuff
C DEPRECATED!!!!
C
C      integer :: myid, numprocs

c..   c-o mesh at the root also ?
      integer :: root_co

c..   parameters for trailing edge flap 
      real :: pxc,flap_ys,flap_ye,dely,dela,rf_tef,angmax_tef
      real :: theta_flap,theta_f0
      integer :: nharmflap
      real :: ampFlap(10),phiFlap(10)
      integer :: iteflap, flpdef

c     Parameters for the overhang flap
      real :: xflph, zflph, flpovh, cgap
      integer :: noh, joh1

c     Misc parameters
      integer :: ivec, ntot, nrev

c     Grid size parameters
      integer gsize,ssize,qsize,dsize,dsize3,size2d,size1d,size2dp

c     dom_type - Domain type
c     irot_dir - Direction of rotation (if zero no rotation)
c     is_wing - Is the mesh a wing, blade/slat?
c     is_ground - Is the mesh a background with ground plane?
      integer :: dom_type, irot_dir
      logical :: is_wing,is_ground, is_body
c     Connectivity type
c     Connectivity mode: Static or dynamic connectivity?
      integer :: conn_type, conn_mode
c     Splitting strategy
      integer :: split_type

c     File units for various outputs
c     reshst - Time history of residuals
      integer :: reshst
      
! Math: Added by Ria for pitch_plunge_wing
! Camli : 3/2015 added deltab (kinematics smooth transthion)
      real :: isroot, root, flap_amp, pitch_knot, pitch_angle, phase, plunge_amp, omega, xac, zac
      real :: kin_rf, kin_freq, kin_phiMax, kin_phiOff, kin_phi_phase
      real :: kin_alphaMax, kin_alphaOff, kin_alpha_phase
      real :: kin_phi, kin_alpha, xoff, yoff, zoff
      real :: kin_beta, kin_phiInit, kin_alphaInit 
      real :: kin_beta_new, kin_phiOff_new, kin_phimax_new
      real :: kin_xoff

!*****************************
c     Added by James Lankford to save restart files for a given interval (2/18/2014)
      integer :: nstart, nstop, ninterval

c     Added by Camli.B declare the cg coordinates
      real :: xcg, ycg, zcg
      integer :: isolno
      integer :: ismth, Ncyc
	
c     Added for subroutine (read_kine_def) by James Lankford (03.02.18)
      integer :: read_flap_kine, read_pitch_kine
      real, allocatable :: flap_kine_data(:), pitch_kine_data(:)

C*******************************************************************
C      Default values of variables
C*******************************************************************
      DATA NMESH /1/
      DATA DX1,DY1,DZ1 / 3*1.0 /
      DATA GAMMA,PR,RMUE,TINF /1.4,0.72,1.0,540./
      DATA KSYM,KROOT /0,1/
      DATA IBCWP /1/
      DATA CNBR /1./
      DATA IREAD /0/
      DATA JMAX,KMAX,LMAX,JTAIL1,KTIP,HALF /109,36,31,19,27,0/
      DATA NSTEPS,NREST,NPNORM,NMOVIE,IMOVIE /500,1000,25,4000,0/
      DATA RESTART_COUNTER /0/
      DATA FSMACH,ALFA,REY,INVISC,LAMIN,ITHIN,ITURB
     <           /0.16,0.0,3900000.,.FALSE.,.FALSE.,.FALSE.,1/
      DATA IUNST,NTAC,ITNMAX,DT,TIMEAC,INDICIAL /3,1,1,.1,1.,0/
      DATA IDUAL, DTPSEUDO /0, 1.0/
      DATA ZGROUND /1.e10/
      DATA USEGCL /.FALSE./
      DATA TESTGCL /.FALSE./
      DATA EPSE,IRHSY,ILIM,TOTIME,ANGMAX /0.01,-3,0,0.0,0.0/
      DATA ILHS,IPRECON,MP /1,.FALSE.,1./
      DATA FMTIP,RARTIO,ALPHA_DOT,BETA1 /0.8,7.0,0.0,0.0/
      DATA IWAKE,CTINP,NBLADE /0,0.00666,1/
      DATA NACOU,NACOUP,NAKP1,NAKP2,NAKP3,NAKP4,NAKP5 /100,20,5*1/
      DATA NWW,NWWQ /1,1000/
      DATA JINT,KINT,LINT /1,1,1/
      DATA IPERT,INITPE,XVOR,ZVOR,RCORE,VORGAM /0,1,-5.,-0.26,.05,.20/
      data imth,ivec/1,0/
      data ntot,nrev/72,3/
      data theta0,theta1c,theta1s,xlam/0.0,0.0,0.0,0.0/
      data irotate,totime0/0,0.0/
      data ideform/0/
      data iconfine,irefine/0,0/
      data idefl,idecay/0,0/
      data beta0,beta1c,beta1s,eflap/0.0,0.0,0.0,0.0/
      data root_co /0/
      data pxc,flap_ys,flap_ye,dely /0.8,1.44,2.31,0.2/
      data dela,rf_tef,angmax_tef /30,1.0,5.0/
      data theta_flap,theta_f0 /0.,0./
      data iteflap,flpdef/0,0/
      data ampFlap/10*0./
      data phiFlap/10*0./
      data xflph,zflph,flpovh,cgap/0.,0.,0.,0./
      data noh,joh1/0,0/
      data conn_type,conn_mode/CONN_JAINA,CONN_DYNAMIC/
      data split_type/SPLIT_AUTO/
      data is_wing,is_ground,is_body/.false.,.false.,.false./
      data irot_dir/1/
      
      ! Math: add kinematics
      data isroot, root, flap_amp, pitch_knot, pitch_angle, phase, plunge_amp, omega, xac, zac/0,0.0,0.0,0,0.0,0.0,0.0,0.0,0.0,0.0/
      data kin_rf, kin_freq, kin_phiMax, kin_phiOff, kin_phi_phase/0.00001,0.0,0.0,0.0,0.0/
      data kin_alphaMax, kin_alphaOff, kin_alpha_phase/0.0,0.0,0.0/
      data xoff, yoff, zoff/0.0,0.0,0.0/
      data kin_beta,kin_phiInit, kin_alphaInit/0.0,0.0,0.0/
      data kin_beta_new, kin_phiOff_new, kin_phimax_new/0.0,0.0,0.0/

      ! Camli : add  nstart, nstop, ninterval and xcg, ycg, zcg
      data nstart, nstop, ninterval/0,-1,1/
      data xcg, ycg, zcg/0.0,0.0,0.0/
      
C*******************************************************************
C      End definitions
C*******************************************************************      

      end module params_global

C*******************************************************************
c..   deflections module
C*******************************************************************

      module deflections

      real, allocatable :: defl_dat(:,:,:),xyzref(:,:)
      real, allocatable :: defl_dat_prev(:,:,:)
      integer iazimuth,icom,def_alloc_flag
      logical def_first

      data icom,def_alloc_flag /1,0/
      data def_first/.true./

      end module deflections

c******************************************************************
c..   twice as refined mesh
c..   used for GCL purposes
c*****************************************************************

      module refmesh

      real, allocatable:: xbig(:,:,:),ybig(:,:,:),zbig(:,:,:),
     <                    xold(:,:,:),yold(:,:,:),zold(:,:,:),
     <                    xole(:,:,:),yole(:,:,:),zole(:,:,:),
     <                     ugp(:,:,:),vgp(:,:,:),wgp(:,:,:),
     <                     ug1(:,:,:),vg1(:,:,:),wg1(:,:,:)

      end module refmesh


c******************************************************************
c..   free wake geometry
c*****************************************************************

      module freewake

      real, allocatable:: pcx(:,:,:),pcy(:,:,:),
     &     pcz(:,:,:),circ(:,:),rrc(:,:)
      real dzeta,ft

      integer iadim,izdim,nr,nw,np,nz,fw_alloc

      data dzeta,fw_alloc /5.,0/

      end module freewake

c********************************************************************
c..	work arrays
c..     more to be added here
c********************************************************************

	module work

	real, allocatable :: uge(:),vge(:),wge(:)
	real, allocatable :: ppp(:),d(:)
        real, allocatable :: xpp(:),ypp(:),zpp(:)
        real, allocatable :: xppp(:),yppp(:),zppp(:)
        integer, allocatable :: kpb(:),kpc(:)

	end module work


c..    background grid octree
c..    made only once, so made sense to keep it as common block

        module bg_octree

	integer max_sor,mbox,lvl

        real, allocatable:: xs1(:),ys1(:),zs1(:)
        real, allocatable:: sboxs(:)
        real, allocatable:: xcls(:),ycls(:),zcls(:)
        integer, allocatable::lks(:,:),ncls(:),inxcbs(:,:)
        integer, allocatable::kns(:),nboxs(:),lboxs(:),nsups(:)
        integer, allocatable::pindex(:)
        
        end module bg_octree

C     Movie module moved to src/overturns/movie_utils.f90

c********************************************************************
c..     parameters for root execution in python
c
c********************************************************************

        module rootVariables

        use TURNS_options
        use bcparam
        use mpi_wrapper
        
        real, allocatable :: s(:),q(:)
        real, allocatable:: qnewt(:),qtn(:)
        real, allocatable:: x(:),y(:),z(:)
        real, allocatable:: bt(:)
        real, allocatable:: xg(:),yg(:),zg(:)
        real, allocatable:: xt2(:),yt2(:),zt2(:)
        real, allocatable:: z1(:),z2(:)
        integer, allocatable :: iblank(:)
        
        integer, allocatable :: loc1(:),loc2(:)
        
        real, allocatable :: xx(:),xy(:),xz(:)
        real, allocatable :: yx(:),yy(:),yz(:)
        real, allocatable :: zx(:),zy(:),zz(:)

        real, allocatable :: svj(:), svk(:), svl(:)
        real, allocatable :: vnaj(:), vnak(:), vnal(:)
        real, allocatable :: xaj(:), yaj(:), zaj(:)
        real, allocatable :: xak(:), yak(:), zak(:)
        real, allocatable :: xal(:), yal(:), zal(:)
        
        real, allocatable :: ug(:),vg(:),wg(:)
        real, allocatable :: turmu(:),wgrid(:),vnu(:),vnu0(:)
        real, allocatable :: zx0(:),zy0(:),zz0(:),zt0(:)
        real, allocatable :: wvec(:,:)
        real, allocatable :: pwork(:),awork(:),bwork(:),cwork(:)
        
        integer, allocatable :: kkr(:),kkp(:)
        integer, allocatable :: imesh(:),idonor(:)
        integer, allocatable :: imesh1(:),idonor1(:)
        
        !integer, allocatable :: ipointer(:)
        !character *40 , allocatable :: filenames(:)
        !integer, allocatable :: jgmx(:),kgmx(:),lgmx(:)
        !integer, allocatable :: jlemx(:),jtail1mx(:),jtail2mx(:),jbotmx(:)
        !integer, allocatable :: ktipmx(:),krootmx(:)
        integer :: ndonor,nfringe
        
        real, allocatable    :: buffer(:), twist(:)
        real, allocatable :: frac(:), frac1(:)
        
        real psi,srot,resrho,rsum,reschk,resmax,x0,y0,z0
        real ct,cq,figofmerit,psi_rot
        integer nrc2,mstop,mstopNode,mstopTmp,im
        !integer ig,ig1,ig2,ig21,igq,igs,igd,igdl,igql,igl
        !character *40 grid_file,soln_file,rest_file, movie_file
        character *128 arg,outfile,integer_string
        real t1,t2,dummy
        integer turnsInit
        
        data turnsInit /0/
        !integer ireq(4),mpistatus(MPI_STATUS_SIZE,4)
        !integer rc,iproc,tag,status(MPI_STATUS_SIZE),ifile
        
        type(bc_t) :: gridbc
        
        end module rootVariables
        
        module pyMOD
        
        integer :: nsa,nts,nkmax
        real,allocatable :: loads(:,:,:)
        real,allocatable :: spandisPY(:),spandist(:)
        integer :: pythonMode
        data pythonMode /0/
        end module pyMOD

        module hqptef
        
        integer :: kfmax
        real,dimension(:),allocatable :: yflap
        real,dimension(:),allocatable :: xupper, xlower
        real,dimension(:),allocatable :: zupper, zlower
        real,dimension(:),allocatable :: xte, flphx
        real,dimension(:),allocatable :: zte, flphz
        real,dimension(:),allocatable :: phiupmx, philowmx
        real,dimension(:),allocatable :: phiupmi, philowmi
        
        end module hqptef

        module tightcouple
        real, allocatable, dimension(:) :: xt2ref, yt2ref, zt2ref
        real, allocatable, dimension(:) :: svjref, svkref, svlref
        end module tightcouple

        module io_filenames
        !!! Variables used for storing input filenames
        ! Modified by James Lankford 03/02/18
        ! Added inputs for reading in flap and pitch kinematics

        integer, parameter, private :: fname_len=128

        ! mesh_inp - The mesh input file
        ! overset_inp - Overset connectivity input file
        ! mesh_positions - Mesh orientation, rotation direction
        ! movie_inp - Movie input file
        ! def_file - Deflections file (for standalone runs)
        ! soln_dir - Directory where solutions must be output to
        ! rest_dir - Directory where restart files must be read from
        ! tef_inp - Trailing edge flap parameters input file
        ! flap_kine_inp - Prescribed flap kinematics
        ! pitch_kine_inp - Prescribed pitch kinematics
        character(len=fname_len) :: mesh_inp
        character(len=fname_len) :: overset_inp
        character(len=fname_len) :: mesh_positions
        character(len=fname_len) :: movie_inp, tef_inp
        character(len=fname_len) :: soln_dir, rest_dir
        character(len=fname_len) :: def_file
        character(len=fname_len) :: palloc_file
        character(len=fname_len) :: ihc_inp
        character(len=fname_len) :: flap_kine_file
        character(len=fname_len) :: pitch_kine_file

        data mesh_inp/'grids/mesh.inp'/
        data overset_inp/'grids/overset.inp'/
        data mesh_positions/'grids/mesh_azi.inp'/
        data soln_dir/'soln/'/
        data rest_dir/'restart/'/
        data movie_inp/'movie.inp'/
        data tef_inp/'tef.inp'/
        data def_file/'defs.mod'/
        data palloc_file/'grids/proc_alloc.inp'/
        data ihc_inp/'grids/ihc.inp'/
        data flap_kine_file/'flap_kine_data.txt'/
        data pitch_kine_file/'pitch_kine_data.txt'/
        
        end module io_filenames


