module motion_par
   use prec
   implicit none
    
   logical           :: glomob    ! global motion
   logical           :: glodina
   logical           :: sysrot

   integer(kind=I4P) :: ngruppi,ncorpi

   real(kind=R8P)                 :: OMEGA_2
   real(kind=R8P),dimension(3)    :: OMEGA

   ! matrici per corpi mobili
   real(kind=R8P),allocatable,dimension(:,:)  :: velgr,omeabs,xcr

contains

   !.......................................................................
   subroutine gridve(ig,xc,ugr)
   !.......................................................................
      use prec
   
      integer(KIND=I4P),intent(in)  :: ig
      real(KIND=R8P),intent(in)     :: xc(3)
      real(KIND=R8P),intent(out)    :: ugr(3)

      real(KIND=R8P)  ::  xo(3)
   
      ! vettore posizione relativa
      xo = xc-xcr(1:3,ig)  
   
      ! velocita' griglia nel riferimento assoluto
      ugr = velgr(1:3,ig) + omeabs(1:3,ig).cross.xo
   
      return
   end subroutine gridve
end module motion_par
