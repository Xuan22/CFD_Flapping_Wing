program overturns

   use mpi_wrapper
   use params_global, only: nsteps
   use turns_api

   implicit none

   character*128 fname

   ! get the input file from command line
   call getarg(1,fname)
   fname=trim(adjustl(fname))

   ! Call the init routine
   call init(fname)

   ! March till the prescribed number of timesteps
   call run(nsteps)

   ! Clean up after ourselves
   call mpi_finalize(ierr)
end program overturns

