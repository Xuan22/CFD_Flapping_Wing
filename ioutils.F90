!!! ioutils - Utility functions for File I/O in Fortran 90
!!!
!!! Functions defined:
!!! 
!!!     get_free_unit - Get the first unused unit number for file
!!!         operations.
!!!     open_file - wrapper around FORTRAN open
!!!     gobble_comment - Skip upto comment close statement
!!!     is_comment - Checks whether the line is start of comment
!!!

module ioutils

#ifdef HAS_ISO_FORTRAN_ENV   
   use, intrinsic :: iso_fortran_env
#endif
   
   implicit none

   !Limits of file units
   integer, parameter, private :: start=9
   integer, parameter, private :: end=999

#ifdef HAS_ISO_FORTRAN_ENV   
   integer, parameter :: STDOUT=OUTPUT_UNIT
   integer, parameter :: STDERR=ERROR_UNIT 
#else
   integer, parameter :: STDOUT=6
   integer, parameter :: STDERR=0
#endif
   
   ! File name string lengths
   integer, parameter :: FLNLEN=128
   ! File comment handling
   integer, parameter :: COMLEN=2, RECLEN=72
   character(len=2), parameter :: COMSTART='#|', COMEND='|#'

   character(len=FLNLEN), private :: BASEDIR=' '

contains
   function get_free_unit(hint) result(un)
      !! Return a free unit that can be used to safely open files.
      !!
      !! Inputs:
      !!    hint - (optional) if present, then use this as the starting
      !!    unit to search for

      integer, optional, intent(in) :: hint
      integer :: un
      logical :: found
      logical :: fexist, fopen

      found=.false.
      if (present(hint)) then
         un=hint
      else
         un=start
      end if

      do while ((.not. found) .and. (un .le. end))
         un=un+1
         inquire(unit=un,exist=fexist,opened=fopen)

         if (fexist .and. (.not. fopen)) then
            found=.true.
         end if
      end do

      !Something is seriously wrong if we get here
      if (.not. found) then
         write(STDERR,*) "ERROR: Cannot find a free file unit"
         stop 'get_free_unit'
      end if
   end function get_free_unit

   subroutine set_basedir(dirnam)
      !! Set the base directory where files will be opened/closed. Must
      !! be absolute and must end with the file/directory separation
      !! character.
      !!
      character(len=*), intent(in) :: dirnam

      BASEDIR=trim(adjustl(dirnam))
   end subroutine set_basedir

   function open_file(flnam,form,status) result(un)
      !! Wrapper function to handle opening of files. Opens file
      !! ${BASEDIR}/flnam. See set_basedir for more details.
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

      filename=trim(adjustl(BASEDIR))//trim(adjustl(flnam))

      un=get_free_unit()
      open(un,file=filename,form=flform,status=flstat)

      write(STDERR,100) trim(adjustl(filename)), un

100   format('# Opened file: ',a,', unit ',I3)
   end function open_file

   subroutine gobble_comment(un)
      !! Given a file unit position on or after the comment start
      !! string consumes until it matches the END_OF_COMMENT string,
      !! see is_comment for testing whether a region should be parsed
      !! as comment.
      integer, intent(in) :: un
      character(len=COMLEN) :: cs

      do 
         read(un,*) cs
         if (cs == COMEND) exit
      end do
   end subroutine gobble_comment

   function is_comment(cs) result(tf)
      !! Given a string tests whether it starts with a comment start
      !! character
      character(len=*) :: cs
      logical :: tf

      if (cs(1:COMLEN) == COMSTART) then
         tf=.true.
      else
         tf=.false.
      end if
   end function is_comment

   subroutine stop_execution(flname,msg)
      !! A simple wrapper around stop to print out a message before
      !! existing on fatal errors.
      character(len=*), intent(in) :: flname, msg

      write(STDERR,1000) trim(adjustl(flname)),msg

      stop
1000  format(//'Fatal error in ',a,/,5x,a,/,'Run aborted'//)
   end subroutine stop_execution
end module ioutils

