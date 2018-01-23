module file_set
    use prec
    implicit none

    logical            :: frst1,frst2,frst3,frst4
    integer(kind=I4P)  :: fres
    character(len=80)  :: filevo  ! time evolution
    character(len=80)  :: filres  ! residual
    character(len=80)  :: filrst  ! restart
    character(len=80)  :: filout  ! solution
    character(len=80)  :: filgrd  ! grid
    character(len=80)  :: filecc  ! topology and boundary condition
    character(len=80)  :: filini  ! initial condition

contains

   !-----------------------------------------------------------------------
   !  returns a free logic unit to be used in open
   !-----------------------------------------------------------------------
   subroutine getunit(Free_Unit)
      use prec
   
      implicit none
      integer(kind=I4P),intent(out) :: Free_Unit ! Free logic unit.
      integer(kind=I4P) :: n1        ! Counter.
      integer(kind=I4P) :: ios       ! Inquiring flag.
      logical::            lopen     ! Inquiring flag.

      Free_Unit = -1      ! initializing free logic unit
      n1=1                ! initializing counter
   
      do
   
         if ((n1/=5).and.(n1/=6).and.(n1/=9)) then    ! reserved units
   
            inquire (unit=n1,opened=lopen,iostat=ios) ! verify logic units
            if (ios==0) then
               if (.not.lopen) then
                  Free_Unit = n1                      ! assignment of free unit
                  return
               endif
            endif
   
         endif
         n1=n1+1                                      ! updating counter
   
      enddo
   
      return
      !-----------------------------------------------------
   end subroutine getunit
end module file_set
