submodule (block_type) block_procedures
   use proc_pointers
   use motion_par

contains
   module subroutine update_rk(self,ark,brk,crk)
   
      class(gen_block),intent(in out) :: self         ! Blocco
      real(kind=R8P)  ,intent(in)     :: ark,brk,crk  ! coeff. RK
 
      integer(kind=I4P) ::  i, j, k
     !==========================================================
     !$omp parallel do collapse(3) schedule(static) &
     !$omp default(none) &
     !$omp shared(ark,brk,crk,self) &
     !$omp private(i,j,k)
     !==========================================================
      do k=1,self%nk
         do j=1,self%nj
            do i=1,self%ni
               self%q(i,j,k) = ark*self%qo(i,j,k)  &
                             + brk*self%q( i,j,k)  &
                             - crk*self%dq(i,j,k)
            end do
         end do
      end do
      !$omp end parallel do 
   end subroutine  update_rk

   module subroutine add_dert(self,dtempo)
      class(gen_block),intent(in out) :: self         ! Blocco
      real(kind=R8P)  ,intent(in)     :: dtempo       ! passo temporale

      integer(kind=I4P)             :: i,j,k
      real(kind=R8P)                ::  voldt
      !==============================================================================
      !$omp parallel do collapse(3) schedule(static) &
      !$omp default(none) &
      !$omp shared(self,dtempo) &
      !$omp private(i,j,k,voldt)
      !==============================================================================
      do k=1,self%nk
         do j=1,self%nj
            do i=1,self%ni
               voldt = self%cell(i,j,k)%Vol/dtempo
               self%dq(i,j,k) = self%dq(i,j,k) + voldt*(self%q(i,j,k)-self%qo(i,j,k))
            end do
         end do
      end do
      !$omp end parallel do ! nowait
   end subroutine add_dert

   module subroutine res_dt(self,dtempo)
      class(gen_block),intent(in out) :: self         ! Blocco
      real(kind=R8P)  ,intent(in)     :: dtempo       ! passo temporale

      integer(kind=I4P)             :: i,j,k
      real(kind=R8P)                :: voldt
      !==============================================================================
      !$omp parallel do collapse(3) schedule(static) &
      !$omp default(none) &
      !$omp shared(self,dtempo) &
      !$omp private(i,j,k,voldt)
      !==============================================================================
      do k=1,self%nk
         do j=1,self%nj
            do i=1,self%ni
               voldt = self%cell(i,j,k)%Vol/dtempo
               if (self%cell(i,j,k)%Chi%tycc /= 0) voldt = voldt/(1d0+voldt*self%cell(i,j,k)%fchi)
               self%dq(i,j,k) = self%dq(i,j,k)*voldt 
            end do
         end do
      end do
      !$omp end parallel do ! nowait
   end subroutine res_dt

   module subroutine euler_flux_balance(self)
      class(gen_block), intent(in out) :: self         ! Blocco

      integer(kind=I4P)    ::  i, j, k     ! indici
      integer(kind=I4P)    :: ni,nj,nk     ! dimensioni
      integer(kind=I4P)    :: ijkmax       ! massima dimensione lineare
      integer(kind=I4P)    :: ierr

      type(gen_var),allocatable,dimension(:) :: fq,vl,vr

      ni = self%ni
      nj = self%nj
      nk = self%nk
      ijkmax = max(ni,nj,nk)
      allocate(fq(0:ijkmax  ),STAT=ierr); call errore('All eul fq',ierr)
      allocate(vl(0:ijkmax+1),STAT=ierr); call errore('All eul vl',ierr)
      allocate(vr(0:ijkmax+1),STAT=ierr); call errore('All eul vr',ierr)

      !=================================================================
      !$omp parallel &
      !$omp default(none) &
      !$omp shared(lrstat,actual_flux,ni,nj,nk,self) &
      !$omp private(i,j,k,vr,vl,fq)
      !=================================================================

      !=================================================================
      !  flussi k
      !=================================================================
      !$omp do collapse (2) schedule(static)
      do j=1,nj
         do i=1,ni
      
            call lrstat(nk,self%q(i,j,-1:nk+2),vl,vr)
            do k=0,nk
               fq(k) = actual_flux(3,self%cell(i,j,k),self%vaux(i,j,k),vl(k),vr(k+1))
            end do
            do k=1,nk
               self%dq(i,j,k) = self%dq(i,j,k) + fq(k)-fq(k-1)
            end do
      
         end do
      end do
      !$omp end do ! nowait
      !====================================================================
      !  flussi j
      !====================================================================
      !$omp do collapse (2) schedule(static)
      do k=1,nk
         do i=1,ni
      
            call lrstat(nj,self%q(i,-1:nj+2,k),vl,vr)
            do j=0,nj
               fq(j) = actual_flux(2,self%cell(i,j,k),self%vaux(i,j,k),vl(j),vr(j+1))
            end do
            do j=1,nj
               self%dq(i,j,k) = self%dq(i,j,k) + fq(j)-fq(j-1)
            end do
      
         end do
      end do
      !$omp end do ! nowait
      !==============================================================================
      !  flussi i
      !==============================================================================
      !$omp do collapse (2) schedule(static)
      do k=1,nk
         do j=1,nj
      
            call lrstat(ni,self%q(-1:ni+2,j,k),vl,vr)
            do i=0,ni
               fq(k) = actual_flux(1,self%cell(i,j,k),self%vaux(i,j,k),vl(i),vr(i+1))
            end do
            do i=1,ni
               self%dq(i,j,k) = self%dq(i,j,k) + fq(i)-fq(i-1)
            end do
      
         end do
      end do
      !$omp end do ! nowait
      !$omp end parallel

      deallocate(  fq,STAT=ierr)  ; call errore('Deall eul fq',ierr)
      deallocate(  vr,STAT=ierr)  ; call errore('Deall eul vr',ierr)
      deallocate(  vl,STAT=ierr)  ; call errore('Deall eul vl',ierr)

   end subroutine euler_flux_balance

   module subroutine ini_sorg(self)
      class(gen_block),intent(in out) :: self         ! Blocco

      integer(kind=I4P)             :: i,j,k
      !==============================================================================
      !$omp parallel do collapse(3) schedule(static) &
      !$omp default(none) &
      !$omp shared(self) &
      !$omp private(i,j,k)
      !==============================================================================
      do k=1,self%nk
         do j=1,self%nj
            do i=1,self%ni
               self%dq(i,j,k) = -self%sq(i,j,k)
            end do
         end do
      end do
      !$omp end parallel do ! nowait
   end subroutine ini_sorg

   module subroutine null_res(self)
      class(gen_block),intent(in out) :: self         ! Blocco

      integer(kind=I4P)             :: i,j,k
      !==============================================================================
      !$omp parallel do collapse(3) schedule(static) &
      !$omp default(none) &
      !$omp shared(self) &
      !$omp private(i,j,k)
      !==============================================================================
      do k=1,self%nk
         do j=1,self%nj
            do i=1,self%ni
               self%dq(i,j,k) = 0d0
            end do
         end do
      end do
      !$omp end parallel do ! nowait
   end subroutine null_res


   module subroutine mem_sorg(self)
      class(gen_block),intent(in out) :: self         ! Blocco

      integer(kind=I4P)             :: i,j,k
      !==============================================================================
      !$omp parallel do collapse(3) schedule(static) &
      !$omp default(none) &
      !$omp shared(self) &
      !$omp private(i,j,k)
      !==============================================================================
      do k=1,self%nk
         do j=1,self%nj
            do i=1,self%ni
               self%sq(i,j,k) = self%dq(i,j,k)
            end do
         end do
      end do
      !$omp end parallel do ! nowait
   end subroutine mem_sorg

   module subroutine grid_vel_face(self,igruppo)
      class(gen_block), intent(in out) :: self         ! Blocco
      integer(kind=I4P),intent(in)     :: igruppo

      integer(kind=I4P)    ::  i, j, k     ! indici
      integer(kind=I4P)    :: ni,nj,nk     ! dimensioni
      integer(kind=I4P)    :: idir
      real(kind=R8P),dimension(3) :: xf,uf

      ni = self%ni
      nj = self%nj
      nk = self%nk

      !=================================================================
      !$omp parallel &
      !$omp default(none) &
      !$omp shared(ni,nj,nk,igruppo,idir,self) &
      !$omp private(i,j,k,xf,uf)
      !=================================================================
      idir = 3
      !=================================================================
      !  velocita' k
      !=================================================================
      !$omp do collapse (3) schedule(static)
      do k=0,nk
         do j=1,nj
            do i=1,ni
               xf = 0.25d0*(self%Cell(i  ,j  ,k)%vertex &
                           +self%Cell(i  ,j-1,k)%vertex &
                           +self%Cell(i-1,j  ,k)%vertex &
                           +self%Cell(i-1,j-1,k)%vertex )
               call gridve(igruppo,xf,uf)
               self%vaux(i,j,k)%v(vaux_un(idir)) = uf.dot.self%cell(i,j,k)%Sn(1:3,idir)
            end do
         end do
      end do
      !$omp end do ! nowait
      !====================================================================
      !  velocita' j
      !====================================================================
      idir = 2
      !$omp do collapse (3) schedule(static)
      do k=1,nk
         do j=0,nj
            do i=1,ni
               xf = 0.25d0*(self%Cell(i  ,j,k  )%vertex &
                           +self%Cell(i  ,j,k-1)%vertex &
                           +self%Cell(i-1,j,k  )%vertex &
                           +self%Cell(i-1,j,k-1)%vertex )
               call gridve(igruppo,xf,uf)
               self%vaux(i,j,k)%v(vaux_un(idir)) = uf.dot.self%cell(i,j,k)%Sn(1:3,idir)
            end do
         end do
      end do
      !$omp end do ! nowait
      !==============================================================================
      !  velocita' i
      !==============================================================================
      idir = 1
      !$omp do collapse (3) schedule(static)
      do k=1,nk
         do j=1,nj
            do i=1,ni
               xf = 0.25d0*(self%Cell(i,j  ,k  )%vertex &
                           +self%Cell(i,j  ,k-1)%vertex &
                           +self%Cell(i,j-1,k  )%vertex &
                           +self%Cell(i,j-1,k-1)%vertex )
               call gridve(igruppo,xf,uf)
               self%vaux(i,j,k)%v(vaux_un(idir)) = uf.dot.self%cell(i,j,k)%Sn(1:3,idir)
            end do
         end do
      end do
      !$omp end do ! nowait
      !$omp end parallel
   end subroutine grid_vel_face

end submodule block_procedures
