module mpi_var
    use prec
    implicit none

    ! Parallelizzazione
    integer(kind=I4P) ::  nproc,myrank

    ! mappe di comunicazione
    type map_comm
       integer(kind=I4P)                            :: npt_from,npt_to
       integer(kind=I4P),allocatable,dimension(:,:) :: i_from,i_to
    contains
       procedure, pass(self) :: create_from        ! Crea mappa di ricezione
       procedure, pass(self) :: destroy_from       ! Distrugge mappa di ricezione
       procedure, pass(self) :: create_to          ! Crea mappa di spedizione
       procedure, pass(self) :: destroy_to         ! Distrugge mappa di spedizione
       procedure, pass(self) :: reduce_size_from   ! Restringe la memoria allocata
       procedure, pass(self) :: increase_size_from ! Aumenta la memoria allocata
    end type map_comm

    type(map_comm), allocatable,dimension(:,:)       :: exchange

    integer(kind=I4P),allocatable,dimension(:) :: iblloc
    integer(kind=I4P),allocatable,dimension(:) :: iblglo

    integer(kind=I4P),parameter :: ANNULLATO=-1

contains

   pure subroutine create_from(self, n)
      ! Crea mappa di ricezione
      class(map_comm)  ,intent(in out) :: self         ! Mappa
      integer(kind=I4P),intent(in out) :: n            ! Dimensioni

      call self%destroy_from

      allocate(self%i_from(4,n)) !  indici punti da ricevere
   
   end subroutine create_from
   
   elemental subroutine destroy_from(self)
      ! Distrugge mappa di ricezione
      class(map_comm)  ,intent(in out) :: self         ! Mappa
    
      if (allocated(self%i_from)) deallocate(self%i_from)
      
   end subroutine destroy_from

   pure subroutine create_to(self, n)
      ! Crea mappa di spedizione
      class(map_comm)  ,intent(in out) :: self         ! Mappa
      integer(kind=I4P),intent(in out) :: n            ! Dimensioni

      call self%destroy_to

      allocate(self%i_to(4,n))  !  indici punti da spedire
   
   end subroutine create_to
   
   elemental subroutine destroy_to(self)
      ! Distrugge mappa di spedizione
      class(map_comm)  ,intent(in out) :: self         ! Mappa
    
      if (allocated(self%i_to)) deallocate(self%i_to)
      
   end subroutine destroy_to

   pure subroutine reduce_size_from(self, n)
      ! Crea mappa di spedizione
      class(map_comm)  ,intent(in out) :: self     ! Mappa
      integer(kind=I4P),intent(in out) :: n        ! Dimensioni

      integer(kind=I4P),allocatable,dimension(:,:)  :: tmp         ! Mappa locale

      if (allocated(tmp)) deallocate(tmp)
      allocate(tmp(4,n))

      tmp(1:4,1:n) = self%i_from(1:4,1:n)
   
      call self%create_from(n)
      self%i_from(1:4,1:n) = tmp(1:4,1:n)

      deallocate(tmp)

   end subroutine reduce_size_from

   pure subroutine increase_size_from(self,n,chunk)
      ! Crea mappa di spedizione
      class(map_comm)  ,intent(in out) :: self   ! Mappa
      integer(kind=I4P),intent(in) :: n          ! Dimensioni occorrenti
      integer(kind=I4P),intent(in) :: chunk      ! incremento

      integer(kind=I4P),allocatable,dimension(:,:) :: tmp  ! Mappa locale
      integer(kind=I4P) :: n_new

      if (size(self%i_from,2) > n) return ! c'e' memoria sufficiente

      if (allocated(tmp)) deallocate(tmp)
      allocate(tmp(4,n))
      tmp(1:4,1:n) = self%i_from(1:4,1:n)
   
      n_new = n+chunk
      call self%create_from(n_new)
      self%i_from(1:4,1:n) = tmp(1:4,1:n)

      deallocate(tmp)

   end subroutine increase_size_from

end module mpi_var
