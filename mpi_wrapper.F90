module mpi_wrapper
!!! A wrapper around MPI definitions so that this can be included in
!!! files and still be able to compile in sequential mode. If MPI
!!! support is not detected, then it uses the dummy_mpi module to
!!! provide empty definitions for most MPI functions used in the
!!! code. Do not include mpif.h or use mpi/dummy_mpi modules directly
!!! in your source code. Use this wrapper module.

#ifdef HAVE_MPIF_H
   use mpi
#else
   use dummy_mpi
#endif 
   use ioutils
   
   implicit none

   !! Wrapper around MPI_REAL or MPI_DOUBLE_PRECISION depending on
   !! compile time delcaration.
#ifdef USE_SINGLE_PRECISION
   integer, parameter :: REAL_TYPE=MPI_REAL
#else
   integer, parameter :: REAL_TYPE=MPI_DOUBLE_PRECISION
#endif

   ! Wrappers around integer and character types
   integer, parameter :: INT_TYPE=MPI_INTEGER
   integer, parameter :: CHAR_TYPE=MPI_CHARACTER

   ! MPI status array
   integer :: stats_mpi(MPI_STATUS_SIZE)

   ! Generalizing if (procid .eq. 0) kind of tests
   integer, parameter :: MASTER=0
   integer, parameter :: FRM_MASTER=0
   integer, parameter :: FRM_SLAVE=1

   !! DEFAULT_COMM :: The current communicator object
   !! PID :: rank of the current process in communicator
   !! NPROCS :: size of the communicator
   !! ierr :: Error value from most MPI calls
   !! isMASTER :: Is rank == zero?
   !! isNthProc :: Is rank == (size-1)?
   !! isParallel :: Are we running in parallel mode? 

   integer :: PID, NPROCS
   integer :: DEFAULT_COMM
   integer :: ierr
   logical :: isMASTER, isNthProc, isParallel

contains
   subroutine initialize_mpi_standalone
#ifdef HAVE_MPIF_H
      call mpi_init(ierr)
      call initialize_mpi(MPI_COMM_WORLD)
#else
      call initialize_mpi(0)
#endif
   end subroutine initialize_mpi_standalone

   subroutine initialize_mpi(comm)
      integer, intent(in) :: comm

      DEFAULT_COMM=comm
#ifdef HAVE_MPIF_H
      call mpi_comm_rank(DEFAULT_COMM,PID,ierr)
      call mpi_comm_size(DEFAULT_COMM,NPROCS,ierr)

      if (PID .eq. MASTER) then
         isMASTER = .true.
      else
         isMASTER = .false.
      end if
      
      if (PID .eq. (NPROCS-1)) then
         isNthProc = .true.
      else
         isNthProc = .false.
      end if

      isParallel=.true.
#else   
      PID=0
      NPROCS=1
      isMASTER=.true.
      isNthProc=.true.
      isParallel=.false.
#endif
   end subroutine initialize_mpi

   subroutine barrier_mpi
      call mpi_barrier(DEFAULT_COMM,ierr)
   end subroutine barrier_mpi

   function open_file_mpi(flnam,form,status) result(un)
      !! Wrapper function to handle opening of files. Opens file
      !! ${BASEDIR}/flnam. See set_basedir for more details. This
      !! function is essentially identially to the open_file function,
      !! with the minor difference that it opens the file with the
      !! processor id suffix to the filename.
      !!
      !! INPUTS:
      !!    flnam - Name of the file to be opened
      !!    form  - FORMATTED or UNFORMATTED (default: formatted)
      !!    status- (default: unknown)
      !!
      character(len=*) :: flnam
      character(len=*), optional :: status, form
      integer :: un

      character(len=50) :: flstat, flform
      character(len=FLNLEN) :: filename
      character(len=10) :: intstr

      if (present(status)) then
         flstat=status
      else
         flstat='unknown'
      end if

      if (present(form)) then
         flform=form
      else
         flform='formatted'
      end if

#ifdef HAVE_MPIF_H
      write(intstr,'(I10)') PID
      filename=trim(adjustl(flnam))//'.'//trim(adjustl(intstr))
#else
      filename=flnam
#endif
      un=open_file(filename,form=flform,status=flstat)
   end function open_file_mpi

   subroutine print_message(msg)
      !! Print the given message with the processor ID prefixed to it.
      character(len=*), intent(in) :: msg

#ifdef HAVE_MPIF_H
      write(STDERR,101) pid,msg
101   format('Proc ',I4,': ',A)
#else
      write(STDERR,'(A)') msg
#endif
   end subroutine print_message

end module mpi_wrapper

