module cell_type
   use prec
   implicit none

   type gen_chi
       ! tipo celle e c.c.
       integer(kind=I4P) :: tycc  ! tipo c.c.
       integer(kind=I4P) :: npcc  ! numero donatori
       integer(kind=I4P),allocatable,dimension(:) :: bd,id,jd,kd
       real(kind=R8P),   allocatable,dimension(:) :: we
   contains
        procedure, pass(self) :: add_X    ! Crea cella chimera
        procedure, pass(self) :: remove_X ! Distrugge cella chimera
   end type gen_chi

   ! cella generica
   type gen_cell

       ! geometria
       real(kind=R8P) :: Center(3)
       real(kind=R8P) :: Vertex(3)
       ! metrica
       real(kind=R8P) :: Vol
       real(kind=R8P) :: Area(3)
       real(kind=R8P) :: Sn(3,3)
       real(kind=R8P) :: Csi(3,3)

       ! cc e chimera
       real(kind=R8P) :: fchi
       real(kind=R8P) :: dwall
       type(gen_chi)  :: Chi

       ! dt locale
       real(kind=R8P) :: dtmat
   end type gen_cell

contains

   pure subroutine add_X(self, npcc)
     ! Crea cella chimera
     class(gen_chi), intent(inout) :: self     ! Cella
     integer(kind=I4P),intent(in)  :: npcc     ! N. donatori
     integer(kind=I4P)             :: ierr
     
     call self%remove_X
     self%npcc = npcc
     if (npcc > 0) then
        allocate(self%we(npcc),STAT=ierr)
        allocate(self%bd(npcc),STAT=ierr)
        allocate(self%id(npcc),STAT=ierr)
        allocate(self%jd(npcc),STAT=ierr)
        allocate(self%kd(npcc),STAT=ierr)
     end if
   end subroutine add_X
 
   elemental subroutine remove_X(self)
      !< Distrugge chimera
      class(gen_chi), intent(inout) :: self !< Cella
    
      self%npcc = 0
      if (allocated(self%we)) deallocate(self%we)
      if (allocated(self%bd)) deallocate(self%bd)
      if (allocated(self%id)) deallocate(self%id)
      if (allocated(self%jd)) deallocate(self%jd)
      if (allocated(self%kd)) deallocate(self%kd)
   end subroutine remove_X

end module cell_type
