module parser

   private :: read_integer,read_real,read_char
   interface read_data
       module procedure read_integer,read_real,read_char
   end interface read_data

contains

   subroutine readline (string,ifile)
   ! ===========================
   ! read real from string
   ! ===========================
       use prec
       implicit none

       integer(kind=I4P),intent(in) :: ifile
       character(len=*),intent(out) :: string

       read(ifile,'(a)') string
       string = to_lower(string)

   end subroutine readline
   
   pure function to_upper (str) result (string)
   !   ==============================
   !   Changes a string to upper case
   !   ==============================
       use prec
       implicit none

       character(*), intent(in) :: str
       character(len=len(str))  :: string
       integer(kind=I4P)        :: ic, i
   
       character(len=26), parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
       character(len=26), parameter :: low = 'abcdefghijklmnopqrstuvwxyz'
   
       string = str
       do i = 1, len_trim(str)
           ic = index(low, str(i:i))
           if (ic > 0) string(i:i) = cap(ic:ic)
       end do
   
   end function to_upper
   
   pure function to_lower (str) result (string)
   !   ==============================
   !   Changes a string to lower case
   !   ==============================
       use prec
       implicit none

       character(*), intent(in) :: str
       character(len=len(str))  :: string
       integer(kind=I4P)        :: ic, i
   
       character(len=26), parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
       character(len=26), parameter :: low = 'abcdefghijklmnopqrstuvwxyz'
   
       string = str
       do i = 1, len_trim(str)
           ic = index(cap, str(i:i))
           if (ic > 0) string(i:i) = low(ic:ic)
       end do
   
   end function to_lower

   subroutine read_integer(input_line,var_name,m)
   ! ===========================
   ! read real from string
   ! ===========================
       use prec
       implicit none
       character(*), intent(in)         :: input_line,var_name
       integer(kind=I4P),intent(in out) :: m

       integer(kind=I4P)                :: ieq,ivar

       ieq = index(input_line, "=")
       if (ieq == 0) return
       ivar = index(input_line,trim(var_name))
       if (ivar == 0 .or. ivar > ieq ) return

       read(input_line(ieq+1:len_trim(input_line)),*) m
   end subroutine read_integer

   subroutine read_real(input_line,var_name,r)
   ! ===========================
   ! read real from string
   ! ===========================
       use prec
       implicit none
       character(*), intent(in)         :: input_line,var_name
       real(kind=R8P),intent(in out)    :: r

       integer(kind=I4P)                :: ieq,ivar

       ieq = index(input_line, "=")
       if (ieq == 0) return
       ivar = index(input_line,trim(var_name))
       if (ivar == 0 .or. ivar > ieq ) return

       read(input_line(ieq+1:len_trim(input_line)),*) r
   end subroutine read_real

   subroutine read_char(input_line,var_name,str)
   ! ===========================
   ! read real from string
   ! ===========================
       use prec
       implicit none
       character(*),intent(in)         :: input_line,var_name
       character(*),intent(in out)     :: str

       integer(kind=I4P)                :: ieq,ivar

       ieq = index(input_line, "=")
       if (ieq == 0) return
       ivar = index(input_line,trim(var_name))
       if (ivar == 0 .or. ivar > ieq ) return

       read(input_line(ieq+1:len_trim(input_line)),'(a)') str
       str = trim(adjustl(str))
   end subroutine read_char

end module parser
