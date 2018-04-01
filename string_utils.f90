module string_utils
!!! This module provides some simple string manipulation utilities. 
!!!
!!! Functions defined:
!!!     upcase - convert string to all uppercase
!!!     downcase - convert string to all lower case
!!!

   implicit none

   ! ASCII character set 'a'-''A'=32
   integer, parameter, private :: Aoffset=ichar('a')-ichar('A')

contains
   pure function upcase(str) result(upper)
      character(len=*), intent(in) :: str
      character(len=len(str)) :: upper

      integer :: i

      do i=1,len(str)
         if (str(i:i) >= "a" .and. str(i:i) <= "z") then
            upper(i:i)=achar(iachar(str(i:i))-Aoffset)
         else
            upper(i:i)=str(i:i)
         end if
      end do
   end function upcase

   pure function downcase(str) result(lower)
      character(len=*), intent(in) :: str
      character(len=len(str)) :: lower
      
      integer :: i

      do i=1,len(str)
         if (str(i:i) >= "A" .and. str(i:i) <= "Z") then 
            lower(i:i)=achar(iachar(str(i:i))+Aoffset)
         else
            lower(i:i)=str(i:i)
         end if
      end do
   end function downcase
end module string_utils

