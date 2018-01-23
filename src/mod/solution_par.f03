module solution_par
   use prec
   implicit none

   logical                                 :: finest,solvin
   logical                                 :: level_set

   real(kind=R8P)                          :: resmx,epsv
!
! integrazione temporale
!
   logical                                 :: unsteady
   logical                                 :: implicito
   integer(kind=I4P)                       :: itime
   real(kind=R8P)                          :: tempo,dtempo
   real(kind=R8P)                          :: temposg,tempo_max
   real(kind=R8P),allocatable,dimension(:) :: stab
!
! integrazione pseudo-temporale
!
   real(kind=R8P)                          :: sqbeta,beta
!
!  coefficienti Runge kutta
!
   integer(kind=I4P)                       :: nrk
   real(kind=R8P),allocatable,dimension(:) :: ark,brk,crk

!
! iterazione Multigrid
!
   integer(kind=I4P)                       :: sgr
   real(kind=R8P),allocatable,dimension(:) :: npre,npost
!
! residui
!
   real(kind=R8P),allocatable,dimension(:) :: resq
!
!  parametri memorizzazione
!
   integer(kind=I4P)                       :: ifsto

!
! matrici per c.c. variabili
!
   real(kind=R8P),allocatable,dimension(:,:,:) :: matcc

end module solution_par
