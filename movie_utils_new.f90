module movie_utils
!!! Implement movie output utilities capable of handling multiple
!!! meshes across parallelized computational blocks. 

   use io_filenames
   use mesh_types
   use domain_info

   implicit none

   ! Defined type to aggregate movie generation metadata for one mesh.
   type movie_type
      ! mesh_id - Mesh ID for the movie info
      ! mg - File unit for grid output (*.3001)
      ! mq - File unit for q output (*.3003)
      integer :: mesh_id
      integer :: mg, mq

      ! nplanes - Number of surfaces/chunk where movie is desired
      integer :: nplanes

      ! jloc, kloc, lloc - {start,end,skip} for each plane.
      integer, dimension(:,:), pointer :: jloc, kloc, lloc
   end type movie_type

   ! movie_list - Global list of movie generation metadata for all
   !              meshes during a computational run.
   ! this_movie - Pointer to the movie metadata for the mesh in the
   !              current processor.
   ! has_movie - Does this mesh have movie metadata associated?
   type(movie_type), dimension(:), pointer, private :: movie_list
   type(movie_type), pointer, private :: this_movie
   logical, private :: has_movie

contains
   subroutine load_movie_info
      !! Read the file, currently called movie.inp, and load up movie
      !! metadata information. Parse through this data and determine
      !! which mesh_masters need to be concerned with outputting movie
      !! files. Finally, open <mesh_name>.3001, <mesh_name>.3003 files
      !! for movie output.
      integer :: un, nlist, n, np, i
      character(len=128) :: mname

      has_movie=.false.
      nullify(this_movie)
      
      un=open_file(movie_inp)
      read(un,*) nlist

      if (nlist == 0) then
         close(un)
         return
      end if

      allocate(movie_list(nlist))

      do n=1,nlist
         read(un,*) mname, movie_list(n)%nplanes

         movie_list(n)%mesh_id=get_global_mesh_id(meshes,mname)

         if (movie_list(n)%nplanes > 0) then
            allocate(movie_list(n)%jloc(3,movie_list(n)%nplanes), &
                  movie_list(n)%kloc(3,movie_list(n)%nplanes), &
                  movie_list(n)%lloc(3,movie_list(n)%nplanes))
            
            do i=1,movie_list(n)%nplanes
               read(un,*) movie_list(n)%jloc(1,i),movie_list(n)%jloc(2,i), &
                     movie_list(n)%jloc(3,i), movie_list(n)%kloc(1,i), &
                     movie_list(n)%kloc(2,i), movie_list(n)%kloc(3,i), &
                     movie_list(n)%lloc(1,i), movie_list(n)%lloc(2,i), &
                     movie_list(n)%lloc(3,i)
            end do

            if (movie_list(n)%mesh_id == this_mesh%mesh_id) then
               has_movie=.true.
               this_movie=>movie_list(n)
            end if
         end if
      end do
      close(un)
      if (has_movie .and. this_isMaster) then
         !!this_movie%mg=open_file(trim(adjustl(this_mesh%mesh_name))//'.3001',&
         !!      form='unformatted')
         !!this_movie%mq=open_file(trim(adjustl(this_mesh%mesh_name))//'.3003', &
         !!      form='unformatted')
         !!close(this_movie%mg)
         !!close(this_movie%mq)
         open(unit=19,file=trim(adjustl(this_mesh%mesh_name))//'.3001',&
               form='unformatted')
         open(unit=23,file=trim(adjustl(this_mesh%mesh_name))//'.3003',&
               form='unformatted')
         this_movie%mg=19
         this_movie%mq=23
      else
         has_movie=.false.
      end if
      
   end subroutine load_movie_info

   subroutine store_movie
      !! Once the global mesh and the corresponding solution has been
      !! assembled, output the relevant sections of the mesh data to
      !! the movie files for post-processing.
      !!
      !! NOTE: This outputs the q-variables by removing Jacobian
      !! scaling. No need to duplicate it in the post-processing step.
      integer :: i, j, k, l, n
      integer :: mg, mq

      real(kind=rdp) :: reypr
      integer, dimension(:,:), pointer :: jloc, kloc, lloc
      
      ! Either not a mesh_master or no movie requested for this mesh.
      if (.not. has_movie) return

      nullify(jloc,kloc,lloc)
      reypr=rey*(fsmach+fmtip)
      jloc=>this_movie%jloc
      kloc=>this_movie%kloc
      lloc=>this_movie%lloc

      !!this_movie%mg=open_file(trim(adjustl(this_mesh%mesh_name))//'.3001',&
      !!      form='unformatted')
      !!this_movie%mq=open_file(trim(adjustl(this_mesh%mesh_name))//'.3003', &
      !!      form='unformatted')
      
      mq=this_movie%mq
      mg=this_movie%mg

      do i=1,this_movie%nplanes

         write(mg)
         write(mg) (((m_grid%x(j,k,l),j=jloc(1,i),jloc(2,i),jloc(3,i)),&
                      k=kloc(1,i),kloc(2,i),kloc(3,i)), &
                      l=lloc(1,i),lloc(2,i),lloc(3,i)), &
                   (((m_grid%y(j,k,l),j=jloc(1,i),jloc(2,i),jloc(3,i)),&
                      k=kloc(1,i),kloc(2,i),kloc(3,i)), &
                      l=lloc(1,i),lloc(2,i),lloc(3,i)), &
                   (((m_grid%z(j,k,l),j=jloc(1,i),jloc(2,i),jloc(3,i)),&
                      k=kloc(1,i),kloc(2,i),kloc(3,i)), &
                      l=lloc(1,i),lloc(2,i),lloc(3,i))

         write(mq) fsmach, alf, reypr, totime
         write(mq) ((((m_q(j,k,l,n),j=jloc(1,i),jloc(2,i),jloc(3,i)),&
               k=kloc(1,i),kloc(2,i),kloc(3,i)), &
               l=lloc(1,i),lloc(2,i),lloc(3,i)), &
               n=1,nv)
      end do
      !!close(mg)
      !!close(mq)
   end subroutine store_movie
end module movie_utils

               
         
