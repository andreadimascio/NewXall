module concon
   use prec
   use var_type
   implicit none
   save; public

   ! lista condizioni al contorno

   integer(kind=I4P),parameter  :: Parete          = -1  ! condizioni naturali
   integer(kind=I4P),parameter  :: Simmetria       = -2
   integer(kind=I4P),parameter  :: Ingresso        = -3
   integer(kind=I4P),parameter  :: Uscita          = -4
   integer(kind=I4P),parameter  :: IngrAss         = -5
   integer(kind=I4P),parameter  :: PressAss        = -6
   integer(kind=I4P),parameter  :: UnormAss        = -7
   integer(kind=I4P),parameter  :: InvRiemAss      = -8
   integer(kind=I4P),parameter  :: Estrapol        = -9
   integer(kind=I4P),parameter  :: PareteTrasl     = -10
   integer(kind=I4P),parameter  :: PareteNonAt     = -11
   integer(kind=I4P),parameter  :: PareteIsoterma  = -12
   integer(kind=I4P),parameter  :: EstrAlter       = -19

   integer(kind=I4P),parameter  :: FacXf        = 20 ! facce chimera - centro faccia
   integer(kind=I4P),parameter  :: FacXf_i0     = 21
   integer(kind=I4P),parameter  :: FacXf_in     = 22
   integer(kind=I4P),parameter  :: FacXf_j0     = 23
   integer(kind=I4P),parameter  :: FacXf_jn     = 24
   integer(kind=I4P),parameter  :: FacXf_k0     = 25
   integer(kind=I4P),parameter  :: FacXf_kn     = 26

   integer(kind=I4P),parameter  :: CellChim     = 27 ! cella interna chimera
   integer(kind=I4P),parameter  :: CellParInt   = 28 ! parete interna

   integer(kind=I4P),parameter  :: FacXc        = 40 ! facce chimera - centro cella
   integer(kind=I4P),parameter  :: FacXc_i0     = 41 
   integer(kind=I4P),parameter  :: FacXc_in     = 42
   integer(kind=I4P),parameter  :: FacXc_j0     = 43
   integer(kind=I4P),parameter  :: FacXc_jn     = 44
   integer(kind=I4P),parameter  :: FacXc_k0     = 45
   integer(kind=I4P),parameter  :: FacXc_kn     = 46

   integer(kind=I4P),parameter  :: FacAd        = 60 ! adiacenza
   integer(kind=I4P),parameter  :: FacAd_i0     = 61
   integer(kind=I4P),parameter  :: FacAd_in     = 62
   integer(kind=I4P),parameter  :: FacAd_j0     = 63
   integer(kind=I4P),parameter  :: FacAd_jn     = 64
   integer(kind=I4P),parameter  :: FacAd_k0     = 65
   integer(kind=I4P),parameter  :: FacAd_kn     = 66

   integer(kind=I4P),parameter  :: Spigolo      = 80

   integer(kind=I4P),parameter  :: Search       = 1000

   type(gen_var)                :: qwal,qinf
   real(kind=R8P)               :: uinf(3)

   type cc_patch
      integer(kind=I4P) ::  bloc
      integer(kind=I4P) ::  imin
      integer(kind=I4P) ::  imax
      integer(kind=I4P) ::  jmin
      integer(kind=I4P) ::  jmax
      integer(kind=I4P) ::  kmin
      integer(kind=I4P) ::  kmax
      integer(kind=I4P) ::  tycc
   end type cc_patch

   integer(kind=I4P)                        :: npatch
   type (cc_patch),allocatable,dimension(:) :: patch

end module concon
