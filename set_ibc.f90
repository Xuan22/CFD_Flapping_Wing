MODULE SET_IMMERSED_BOUNDARY

  USE immersedBoundVars
  use domain_info
  use mpi_wrapper

  IMPLICIT NONE

  PRIVATE :: assemble_q!, relay_q_to_blocks

  PUBLIC :: SET_IBC

  CONTAINS

   SUBROUTINE SET_IBC(q)

      !----------------------------------------------------------------
      !
      ! applies immersed boundary conditions for velocity
      !
      !----------------------------------------------------------------
      implicit none
      !----------------------------------------------------------------
      ! Incoming Variables
      real(kind=rdp),intent(inout),dimension(jmax,kmax,lmax,nd) :: q

      ! Local variables
      integer(kind=idp) :: idx,bi,bj,bk,i,j,k,l,nb,comp
      integer(kind=idp) :: jmx,kmx,lmx,in
      real(kind=rdp)   :: c1,c2,c3,c4,c5,c6,c7,c8,cnode,ctip
      real(kind=rdp)   :: phi1,phi2,phi3,phi4
      real(kind=rdp)   :: phi5,phi6,phi7,phi8,phitip,velb,V_Total_Sqr,inv_rho,gamma_minus_one
      REAL(KIND=rdp),ALLOCATABLE,DIMENSION(:,:,:,:) :: Qtmp
      REAL(KIND=rdp),ALLOCATABLE,DIMENSION(:,:,:,:) :: qg

      !integer :: logg
      !character*32 :: flname

      ! First executable statement

      if (this_mesh%dom_type.ne.DOM_IBC) goto 100 ! Already checked in turns_api
      if (parallel_ibc==-1) goto 100 ! No bound nodes were found

      ! If not doing parallel IBC: assemble q variables on the IBC master proc
      if (parallel_ibc==0) then
         jmx = this_mesh%jm
         kmx = this_mesh%km
         lmx = this_mesh%lm
         allocate(qg(jmx,kmx,lmx,nv))
         call assemble_q(q,qg)
      elseif (parallel_ibc==1) then
         jmx = jmax
         kmx = kmax
         lmx = lmax
         allocate(qg(jmx,kmx,lmx,nv))
         do l=1,lmx
            do k=1,kmx
               do j=1,jmx
                  do in=1,nv
                     qg(j,k,l,in)=q(j,k,l,in)
                  end do
               end do
            end do
         end do
      else
         print*,'ERROR: parrallel_ibc should have been set to 0 or 1, check!!!'
         stop
      end if

      ! do the interpolation on the IBC master proc
      master_check: if ((parallel_ibc==0.and.this_isMaster).or.parallel_ibc==1) then
      print*,'Doing IBC interpolations on PID:',PID,nbound
      gamma_minus_one = 1.4-1.0

      ALLOCATE(Qtmp(jmx,kmx,lmx,nv))

      ! Compute and save primitive variables (rho,u,v,w,p) to Qtmp
      do k = 1,lmx
      do j = 1,kmx
      do i = 1,jmx
         inv_rho = one / qg(i,j,k,1)
         Qtmp(i,j,k,1) = qg(i,j,k,1)
         Qtmp(i,j,k,2) = qg(i,j,k,2)*inv_rho
         Qtmp(i,j,k,3) = qg(i,j,k,3)*inv_rho
         Qtmp(i,j,k,4) = qg(i,j,k,4)*inv_rho
         V_Total_Sqr = Qtmp(i,j,k,2)*Qtmp(i,j,k,2) + Qtmp(i,j,k,3)*Qtmp(i,j,k,3) + Qtmp(i,j,k,4)*Qtmp(i,j,k,4)
         Qtmp(i,j,k,5) = gamma_minus_one *( qg(i,j,k,5) - 0.5_rdp*Qtmp(i,j,k,1)*V_Total_Sqr )
      enddo
      enddo
      enddo


         do idx=1,nbound
 
            ! Boundary node global indices
            bi = boundNode(idx)%i
            bj = boundNode(idx)%j
            bk = boundNode(idx)%k
            if (parallel_ibc==1) then ! Indices have to be converted to local indices
               bi = bi - jkl_lims(JDIR) + 1 
               bj = bj - jkl_lims(KDIR) + 1 
               bk = bk - jkl_lims(LDIR) + 1 
            end if

            if (parallel_ibc==0.or.(parallel_ibc==1 .and. &
               (bi>=send_lims(JDIR).and.bi<=send_lims(JDIR+3)) .and. &
               (bj>=send_lims(KDIR).and.bj<=send_lims(KDIR+3)) .and. &
               (bk>=send_lims(LDIR).and.bk<=send_lims(LDIR+3)))) then

               ! Leading interpolation node global indices
               i = boundNode(idx)%iint
               j = boundNode(idx)%jint
               k = boundNode(idx)%kint
               if (parallel_ibc==1) then ! Indices have to be converted to local indices
                  i = i - jkl_lims(JDIR) + 1 
                  j = j - jkl_lims(KDIR) + 1 
                  k = k - jkl_lims(LDIR) + 1 
               end if
 
               ! trilinear interpolation coefficients
               c1 = boundNode(idx)%interp(1)
               c2 = boundNode(idx)%interp(2)
               c3 = boundNode(idx)%interp(3)
               c4 = boundNode(idx)%interp(4)
               c5 = boundNode(idx)%interp(5)
               c6 = boundNode(idx)%interp(6)        
               c7 = boundNode(idx)%interp(7)
               c8 = boundNode(idx)%interp(8)
       
               ! linear interp coeffs between interpolated and boundary nodes
               ctip  = boundNode(idx)%d2int_ratio
               cnode = boundNode(idx)%d2surf_ratio
  
               do comp = 1,nv
                  ! 8 velocity nodes used in trilinear interpolation
                  phi1 = Qtmp(i,j,k,comp)
                  phi2 = Qtmp(i+1,j,k,comp)
                  phi3 = Qtmp(i,j+1,k,comp)
                  phi4 = Qtmp(i+1,j+1,k,comp)
                  phi5 = Qtmp(i,j,k+1,comp)
                  phi6 = Qtmp(i+1,j,k+1,comp)
                  phi7 = Qtmp(i,j+1,k+1,comp)
                  phi8 = Qtmp(i+1,j+1,k+1,comp)
    
                  ! interpolated value
                  phitip = phi1*c1+phi2*c2+phi3*c3+phi4*c4     &
                        +  phi5*c5+phi6*c6+phi7*c7+phi8*c8
               
                  if (comp==1.or.comp==5) then
                     Qtmp(bi,bj,bk,comp) = phitip
                  else
                     Qtmp(bi,bj,bk,comp) = -cnode*phitip/ctip
                  endif
               enddo
 
               ! update q
               qg(bi,bj,bk,1) = Qtmp(bi,bj,bk,1)
               qg(bi,bj,bk,2) = Qtmp(bi,bj,bk,1)*Qtmp(bi,bj,bk,2)
               qg(bi,bj,bk,3) = Qtmp(bi,bj,bk,1)*Qtmp(bi,bj,bk,3)
               qg(bi,bj,bk,4) = Qtmp(bi,bj,bk,1)*Qtmp(bi,bj,bk,4)
               V_Total_Sqr = Qtmp(bi,bj,bk,2)*Qtmp(bi,bj,bk,2) + Qtmp(bi,bj,bk,3)*Qtmp(bi,bj,bk,3) &
                             + Qtmp(bi,bj,bk,4)*Qtmp(bi,bj,bk,4)
               qg(bi,bj,bk,5) = Qtmp(bi,bj,bk,5)/gamma_minus_one + 0.5_rdp*Qtmp(bi,bj,bk,1)*V_Total_Sqr
            endif
         enddo

         deallocate(Qtmp)
      end if master_check
     
      ! Set q variables from qg
      if (parallel_ibc==0) then
         call relay_q_to_blocks(qg,q)
      elseif (parallel_ibc==1) then
         do l=1,lmx
            do k=1,kmx
               do j=1,jmx
                  do in=1,nv
                     q(j,k,l,in)=qg(j,k,l,in)
                  end do
               end do
            end do
         end do
      end if
      deallocate(qg)

100   call barrier_mpi

      !print*,'STOPING ON ALL PROCS'
      !stop

      return

      end subroutine SET_IBC

   subroutine assemble_q(q,qg)
      real(kind=rdp), dimension(jmax,kmax,lmax,nd),intent(in) :: q
      real(kind=rdp), dimension(this_mesh%jm,this_mesh%km,this_mesh%lm,nv),intent(inout) :: qg

      integer :: j, k, l, bufsize, i, ii, in, n
      integer :: jj, kk, ll
      real(kind=rdp), dimension(:), allocatable :: qbuf
      integer, dimension(:), pointer :: bids
      integer, dimension(:,:), pointer :: bst, bend

      nullify(bids,bst,bend)
      bufsize=b_jm*b_km*b_lm*nv
      allocate(qbuf(bufsize))

      !call barrier_mpi

      master_check: if (this_isMaster) then
         bids=>this_mesh%block_id
         bst=>this_mesh%start
         bend=>this_mesh%end
         mesh_blocks: do n=1,this_mesh%nblocks
            if (bids(n) == (PID+1)) then 
               !Local so copy data
               do l=bst(LDIR,n),bend(LDIR,n)
                  ll=send_lims(LDIR)+(l-bst(LDIR,n))
                  do k=bst(KDIR,n),bend(KDIR,n)
                     kk=send_lims(KDIR)+(k-bst(KDIR,n))
                     do j=bst(JDIR,n),bend(JDIR,n)
                        jj=send_lims(JDIR)+(j-bst(JDIR,n))
                        do in=1,nv
                           qg(j,k,l,in)=q(jj,kk,ll,in)*q(jj,kk,ll,6)
                        end do
                     end do
                  end do
               end do
            else ! Receive data from remote processor
               call mpi_recv(qbuf,bufsize,REAL_TYPE,bids(n)-1,100,&
                     DEFAULT_COMM,stats_mpi,ierr)
               !print*,'received: from PID:',bids(n)-1,'to',PID
               ii=0
               do l=bst(LDIR,n),bend(LDIR,n)
                  do k=bst(KDIR,n),bend(KDIR,n)
                     do j=bst(JDIR,n),bend(JDIR,n)
                        do in=1,nv
                           qg(j,k,l,in)=qbuf(ii+in)
                        end do
                        ii=ii+nv
                     end do
                  end do
              end do
            end if
         end do mesh_blocks
      else  master_check ! We will send information
         ii=0
         do l=send_lims(LDIR),send_lims(LDIR+3)
            do k=send_lims(KDIR),send_lims(KDIR+3)
               do j=send_lims(JDIR),send_lims(JDIR+3)
                  do in=1,nv
                     qbuf(ii+in)=q(j,k,l,in)*q(j,k,l,6)
                  end do
                  ii=ii+nv
               end do
            end do
         end do
         call mpi_bsend(qbuf,bufsize,REAL_TYPE,(this_mesh%master_id-1),&
               100,DEFAULT_COMM,ierr)
         !print*,'sent: from PID:',PID,'to',this_mesh%master_id-1
      end if master_check

      !call barrier_mpi
      deallocate(qbuf)
   end subroutine assemble_q
   
   subroutine relay_q_to_blocks(qg,q)
      real(kind=rdp), dimension(this_mesh%jm,this_mesh%km,this_mesh%lm,nv),intent(in) :: qg
      real(kind=rdp), dimension(jmax,kmax,lmax,nd),intent(inout) :: q

      integer :: j, k, l, bufsize, i, tag, n, in
      integer :: j1, k1, l1
      real, dimension(:), allocatable :: qbuf
      integer, dimension(:), pointer :: bids
      integer, dimension(:,:), pointer :: bst, bend
      integer, dimension(6) :: blk_lim

      tag=101
      bufsize=b_jm*b_km*b_lm*nv
      allocate(qbuf(bufsize))

      !call print_message("Updating block iblanks")
      !call barrier_mpi
      master_check: if (this_isMaster) then
         bids=>this_mesh%block_id
         bst=>this_mesh%start
         bend=>this_mesh%end
         mesh_blocks: do n=1,this_mesh%nblocks
            if (bids(n) == (PID+1)) then
               do l=jkl_lims(LDIR),jkl_lims(LDIR+3)
                  l1=l-jkl_lims(LDIR)+1
                  do k=jkl_lims(KDIR),jkl_lims(KDIR+3)
                     k1=k-jkl_lims(KDIR)+1
                     do j=jkl_lims(JDIR),jkl_lims(JDIR+3)
                        j1=j-jkl_lims(JDIR)+1
                        do in=1,nv
                           q(j1,k1,l1,in)=qg(j,k,l,in)/q(j1,k1,l1,6)
                        end do
                     end do
                  end do
               end do
            else
               call determine_send_limits
               i=0
               do l=blk_lim(LDIR),blk_lim(LDIR+3)
                  do k=blk_lim(KDIR),blk_lim(KDIR+3)
                     do j=blk_lim(JDIR),blk_lim(JDIR+3)
                        do in=1,nv
                           i=i+1
                           qbuf(i)=qg(j,k,l,in)
                        end do
                        if (qbuf(i-4)==0.0) print*,'problem on PID:',PID,j,k,l,qbuf(i-3),qbuf(i-2)
                     end do
                  end do
               end do
               call mpi_bsend(qbuf,bufsize,REAL_TYPE,bids(n)-1,tag, &
                     DEFAULT_COMM,ierr)
               !print*,'sent: from PID:',PID,'to',bids(n)-1,bufsize,i
            end if
         end do mesh_blocks
      else master_check
         call mpi_recv(qbuf,bufsize,real_type,(this_mesh%master_id-1),tag,&
               default_comm,stats_mpi,ierr)
         i=0
         do l=1,lmax
            do k=1,kmax
               do j=1,jmax
                  do in=1,nv
                     i=i+1
                     q(j,k,l,in)=qbuf(i)/q(j,k,l,6)
                  end do
               end do
            end do
         end do
         !print*,'received: from PID:',this_mesh%master_id-1,'to',PID,bufsize,i
      end if master_check

      !call barrier_mpi
      deallocate(qbuf)
   contains
      subroutine determine_send_limits
         !! Determine the lower and upper j,k,l limits of each block
         !! taking ghost points into account. This information isn't
         !! stored in the master processor, so we need to calculate
         !! them.
         integer, dimension(6) :: extents
         integer, dimension(3) :: maxes

         extents=(/ bst(JDIR,n), bst(KDIR,n), bst(LDIR,n), &
               bend(JDIR,n), bend(KDIR,n), bend(LDIR,n) /)
         maxes=(/ this_mesh%jm, this_mesh%km, this_mesh%lm /)

         do i=1,3
            if (extents(i) == 1) then 
               blk_lim(i) = extents(i)
            else
               blk_lim(i) = extents(i)-NGHOST
            end if
         end do
         
         ! The max faces
         do i=4,6
            if (extents(i) == maxes(i-3)) then 
               blk_lim(i) = extents(i)
            else
               blk_lim(i) = extents(i)+NGHOST
            end if
         end do
      end subroutine determine_send_limits
   end subroutine relay_q_to_blocks

END MODULE SET_IMMERSED_BOUNDARY
