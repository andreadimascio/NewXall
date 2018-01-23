module block_type

   use prec
   use var_type
   use cell_type
 
   implicit none
 
   public
 
   integer(kind=I4P) :: nbltot,nbltlo,nbl,ngr
   integer(kind=I4P) :: nvar      ! numero di variabili di stato
   integer(kind=I4P) :: naux      ! numero variabili ausiliarie
 
   type gen_block
       integer(kind=I4P) :: ni,nj,nk  ! dimensioni
 
       type(gen_cell),allocatable,dimension(:,:,:)   :: cell
       type(gen_var),allocatable,dimension(:,:,:)   :: q,qo,dq,sq
       type(gen_var),allocatable,dimension(:,:,:)   :: vaux
   contains
        procedure, pass(self) :: create_BMP   ! Crea blocco per multiprocessore
        procedure, pass(self) :: destroy_BMP  ! Distrugge blocco per multiprocessore

        ! block metric procedures
        procedure, pass(self) :: centroid           ! calcolo centri cella
        procedure, pass(self) :: faces              ! calcolo centri cella
        procedure, pass(self) :: volume             ! calcolo centri cella
        procedure, pass(self) :: tens               ! calcolo centri cella

        ! block variable procedures
        procedure, pass(self) :: null_res           ! inizializzazione res = 0
        procedure, pass(self) :: ini_sorg           ! inizializzazione res = -sorg
        procedure, pass(self) :: mem_sorg           ! memorizzazione sorgenti
        procedure, pass(self) :: euler_flux_balance ! bilancio flussi euleriani
        procedure, pass(self) :: visc_flux_balance  ! bilancio flussi viscosi
        procedure, pass(self) :: grid_vel_face      ! velocità griglia sulle facce
        procedure, pass(self) :: level_set_ext      ! estrapolazione level set fuori dominio
        procedure, pass(self) :: add_dert           ! aggiunta derivate temporali
        procedure, pass(self) :: res_dt             ! residui*dtempo/(1+dt/vol*fchi)
        procedure, pass(self) :: update_rk          ! bilancio flussi euleriani
   end type gen_block

   interface 

      ! block metric procedures
      module subroutine centroid(self)
         class(gen_block),intent(in out) :: self         ! Blocco
      end subroutine centroid
      module subroutine faces(self)
         class(gen_block),intent(in out) :: self    
      end subroutine faces
      module subroutine volume(self)
         class(gen_block),intent(in out) :: self    
      end subroutine volume
      module subroutine tens(self)
         class(gen_block),intent(in out) :: self    
      end subroutine tens

      ! block variable procedures
      module subroutine null_res(self)
         class(gen_block),intent(in out) :: self         ! Blocco
      end subroutine null_res

      module subroutine ini_sorg(self)
         class(gen_block),intent(in out) :: self   
      end subroutine ini_sorg

      module subroutine mem_sorg(self)
         class(gen_block),intent(in out) :: self    
      end subroutine mem_sorg

      module subroutine euler_flux_balance(self)
         class(gen_block), intent(in out) :: self 
      end subroutine euler_flux_balance

      module subroutine visc_flux_balance(self)
         class(gen_block), intent(in out) :: self  
      end subroutine visc_flux_balance

      module subroutine grid_vel_face(self,igruppo)
         class(gen_block), intent(in out) :: self   
         integer(kind=I4P),intent(in)     :: igruppo
      end subroutine grid_vel_face

      module subroutine level_set_ext(self)
         class(gen_block), intent(in out) :: self   
      end subroutine level_set_ext

      module subroutine add_dert(self,dtempo)
         class(gen_block),intent(in out) :: self         ! Blocco
         real(kind=R8P)  ,intent(in)     :: dtempo       ! passo temporale
      end subroutine add_dert

      module subroutine res_dt(self,dtempo)
         class(gen_block),intent(in out) :: self      
         real(kind=R8P)  ,intent(in)     :: dtempo    
      end subroutine res_dt

      module subroutine update_rk(self,ark,brk,crk)
         class(gen_block),intent(in out) :: self         ! Blocco
         real(kind=R8P)  ,intent(in)     :: ark,brk,crk  ! coeff. RK
      end subroutine update_rk

   end interface 

   type mg_block
        integer(kind=I4P) :: gruppo    
        integer(kind=I4P) :: corpo
        type(gen_block),allocatable,dimension(:) :: gr
   contains
        procedure, pass(self) :: create_B     ! Crea blocco.
        procedure, pass(self) :: destroy_B    ! Distrugge blocco.

        ! multigrid procedures
        procedure, pass(self) :: prolong  ! prolongation
        procedure, pass(self) :: restrict ! restriction
        procedure, pass(self) :: nullsorg ! set source term to zero
        procedure, pass(self) :: collres  ! collect residual into source
        procedure, pass(self) :: calcvar  ! calculate variation
        procedure, pass(self) :: prolvar  ! prolong variation
   end type mg_block

   interface
       module subroutine prolong(self,igr)
         class(mg_block),intent(in out) :: self         ! MG block
         integer(kind=I4P),intent(in)   :: igr          ! working grid level
       end subroutine prolong
       module subroutine restrict(self,igr)
         class(mg_block),intent(in out) :: self  
         integer(kind=I4P),intent(in)   :: igr  
       end subroutine restrict

       module subroutine nullsorg(self,igr)
         class(mg_block),intent(in out) :: self  
         integer(kind=I4P),intent(in)   :: igr  
       end subroutine nullsorg

       module subroutine collres(self,igr)
         class(mg_block),intent(in out) :: self  
         integer(kind=I4P),intent(in)   :: igr  
       end subroutine collres

       module subroutine calcvar(self,igr)
         class(mg_block),intent(in out) :: self 
         integer(kind=I4P),intent(in)   :: igr  
       end subroutine calcvar

       module subroutine prolvar(self,igr,ls)
         class(mg_block),intent(in out) :: self         ! MG block
         integer(kind=I4P),intent(in)   :: igr          ! working grid level
         logical                        :: ls           ! level set control
       end subroutine prolvar
   end interface


contains

   pure subroutine create_B(self, ni, nj, nk, gruppo, corpo)
      ! Crea blocco
      class(mg_block)   ,intent(in out) :: self         ! Blocco
      integer(kind=I4P),intent(in out) :: ni, nj, nk   ! Dimensioni
      integer(kind=I4P),optional,intent(in out) :: gruppo, corpo

      integer(kind=I4P)               :: i, j, k      ! Indici
      integer(kind=I4P)               :: igr 
      
      call self%destroy_B
      self%gruppo = 0
      self%corpo = 0
      if (present(gruppo)) self%gruppo = gruppo
      if (present(corpo )) self%corpo  = corpo

      allocate(self%gr(ngr)) ! allocazione livelli

      do igr=1,ngr


         self%gr(igr)%ni = ni
         self%gr(igr)%nj = nj
         self%gr(igr)%nk = nk

         allocate(self%gr(igr)%cell(-2:ni+2,-2:nj+2,-2:nk+2))
         allocate(self%gr(igr)%q(   -1:ni+2,-1:nj+2,-1:nk+2))
         allocate(self%gr(igr)%qo(  -1:ni+2,-1:nj+2,-1:nk+2))
         allocate(self%gr(igr)%dq(  -1:ni+2,-1:nj+2,-1:nk+2))
         allocate(self%gr(igr)%sq(  -1:ni+2,-1:nj+2,-1:nk+2))
         allocate(self%gr(igr)%vaux(-1:ni+2,-1:nj+2,-1:nk+2))
         do k=-1,nk+2
           do j=-1,nj+2
             do i=-1,ni+2
                allocate(self%gr(igr)%q(   i,j,k)%v(1:nvar))
                allocate(self%gr(igr)%qo(  i,j,k)%v(1:nvar))
                allocate(self%gr(igr)%dq(  i,j,k)%v(1:nvar))
                allocate(self%gr(igr)%sq(  i,j,k)%v(1:nvar))
                allocate(self%gr(igr)%vaux(i,j,k)%v(1:naux))
             end do
           end do
         end do

         ni = ni/2 ! dimensioni griglia successiva
         nj = nj/2
         nk = nk/2
      end do
   
   end subroutine create_B
   
   elemental subroutine destroy_B(self)
      ! Distrugge blocco
      class(mg_block), intent(in out) :: self 
      integer(kind=I4P)               :: igr 
    
      do igr=1,ngr
         if (allocated(self%gr(igr)%cell)) deallocate(self%gr(igr)%cell)
         if (allocated(self%gr(igr)%q   )) deallocate(self%gr(igr)%q   )
         if (allocated(self%gr(igr)%qo  )) deallocate(self%gr(igr)%q   )
         if (allocated(self%gr(igr)%dq  )) deallocate(self%gr(igr)%dq  )
         if (allocated(self%gr(igr)%sq  )) deallocate(self%gr(igr)%sq  )
         if (allocated(self%gr(igr)%vaux)) deallocate(self%gr(igr)%sq  )
      end do

      if (allocated(self%gr)) deallocate(self%gr)
      
   end subroutine destroy_B

   pure subroutine create_BMP(self,n)
      ! Crea blocco
      class(gen_block) ,intent(in out) :: self ! Blocco
      integer(kind=I4P),intent(in    ) :: n    ! Dimensioni

      integer(kind=I4P)               :: i     ! Indici
      
      call self%destroy_BMP
      self%ni = n
      self%nj = 1
      self%nk = 1

      allocate(self%q(1:n,1,1))
      do i=1,n
         allocate(self%q(i,1,1)%v(1:nvar))
      end do
   
   end subroutine create_BMP
   
   elemental subroutine destroy_BMP(self)
      ! Distrugge blocco
      class(gen_block), intent(in out) :: self 
    
      if (allocated(self%q   )) deallocate(self%q   )
   end subroutine destroy_BMP

end module block_type
