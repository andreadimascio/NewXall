module chimera_type

   use prec
   implicit none

   type block_map ! Mappa del blocco
      integer(kind=I4P)           :: priority
      integer(kind=I4P)           :: ni,nj,nk
      real(kind=R8P),dimension(3) :: low,hig
      integer(kind=I4P),allocatable,dimension(:,:,:) :: imin,imax,jmin,jmax,kmin,kmax
   contains
       procedure, pass(self) :: create_M     ! Crea mappa
       procedure, pass(self) :: destroy_M    ! Distrugge mappa
   end type block_map

   type X_point ! punto tipo chimera
      integer(kind=I4P)           :: b ,i ,j ,k  ! identificativo locale
      integer(kind=I4P)           :: bd,id,jd,kd ! indici donatore
      integer(kind=I4P)           :: tycc        ! tipo cella
      real(kind=R8P)              :: dsize  ! dimensione donatore
      real(kind=R8P),dimension(3) :: xc     ! coordinate richiedente
   end type X_point

   type X_array ! array di punti chimera
      integer(kind=I4P)  :: n
      type(X_point),allocatable,dimension(:) :: P
   contains
       procedure, pass(self) :: create_XA     ! Crea array chimera
       procedure, pass(self) :: destroy_XA    ! Distrugge array chimera
       procedure, pass(self) :: increase_size_XA ! Aumenta dimensioni
       procedure, pass(self) :: reduce_size_XA   ! Riduce dimensioni
   end type X_array
   interface assignment (=) 
       procedure X_eq_X
   end interface assignment (=) 

   type (block_map),allocatable,dimension(:) :: mappa
   type (X_array),allocatable,dimension(:)   :: Xpt_req,Xpt_send
   real(kind=R8P)                            :: chimera_bl

contains

   pure subroutine X_eq_X(p1,p2)
       type (X_point),intent(out) :: p1
       type (X_point),intent(in)  :: p2

       p1%b     = p2%b
       p1%i     = p2%i
       p1%j     = p2%j
       p1%k     = p2%k

       p1%bd    = p2%bd
       p1%id    = p2%id
       p1%jd    = p2%jd
       p1%kd    = p2%kd

       p1%tycc  = p2%tycc

       p1%dsize = p2%dsize

       p1%xc    = p2%xc
   end subroutine X_eq_X

   pure subroutine create_M(self, ni, nj, nk)
      ! Crea blocco
      class(block_map) ,intent(in out) :: self         ! Mappa del blocco
      integer(kind=I4P),intent(in    ) :: ni, nj, nk   ! Dimensioni

      call self%destroy_M
      self%ni = ni
      self%nj = nj
      self%nk = nk

      self%low = 0d0   ! coordinate minimo
      self%hig = 0d0   ! coordinate massimo

      allocate(self%imin(ni,nj,nk))
      allocate(self%imax(ni,nj,nk))
      allocate(self%jmin(ni,nj,nk))
      allocate(self%jmax(ni,nj,nk))
      allocate(self%kmin(ni,nj,nk))
      allocate(self%kmax(ni,nj,nk))
   
   end subroutine create_M
   
   elemental subroutine destroy_M(self)
      ! Distrugge blocco
      class(block_map) ,intent(in out) :: self         ! Mappa del blocco
    
      if (allocated(self%imin)) deallocate(self%imin)
      if (allocated(self%imax)) deallocate(self%imax)
      if (allocated(self%jmin)) deallocate(self%jmin)
      if (allocated(self%jmax)) deallocate(self%jmax)
      if (allocated(self%kmin)) deallocate(self%kmin)
      if (allocated(self%kmax)) deallocate(self%kmax)
      
   end subroutine destroy_M

   pure subroutine create_XA(self, n)
      ! Crea blocco
      class(X_array)   ,intent(in out) :: self         ! Array chimera
      integer(kind=I4P),intent(in    ) :: n            ! Dimensioni

      call self%destroy_XA
      self%n = n
      allocate(self%P(n))
   
   end subroutine create_XA
   
   elemental subroutine destroy_XA(self)
      ! Distrugge blocco
      class(X_array) ,intent(in out) :: self         ! Mappa del blocco
    
      if (allocated(self%P)) deallocate(self%P)
   end subroutine destroy_XA

   pure subroutine reduce_size_XA(self, n)
      ! Crea mappa di spedizione
      class(X_array)   ,intent(in out) :: self     ! Punto chimera
      integer(kind=I4P),intent(in    ) :: n        ! Dimensioni

      type(X_array) :: tmp  ! array locale 

      if (self%n == n) return ! dimensioni corrette

      call tmp%create_XA(n)
      tmp%P(1:n) = self%P(1:n)
   
      call self%create_XA(n)
      self%P(1:n) = tmp%P(1:n)

   end subroutine reduce_size_XA

   pure subroutine increase_size_XA(self, n, chunk)
      ! Crea mappa di spedizione
      class(X_array)   ,intent(in out) :: self     ! Punto chimera
      integer(kind=I4P),intent(in)     :: n        ! Dimensioni
      integer(kind=I4P),intent(in)     :: chunk      ! incremento

      type(X_array) :: tmp  ! array locale 

      if (self%n >= n) return ! memoria sufficiente

      call tmp%create_XA(n)
      tmp%P(1:n) = self%P(1:n)
   
      call self%create_XA(n+chunk)
      self%P(1:n) = tmp%P(1:n)

   end subroutine increase_size_XA

end module chimera_type
