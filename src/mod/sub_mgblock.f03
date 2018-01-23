submodule (block_type) MGblock_procedures

contains

   module subroutine prolong(self,igr)
      implicit none
      class(mg_block),intent(in out) :: self         ! MG block
      integer(kind=I4P),intent(in)   :: igr          ! working grid level
   
      integer(KIND=I4P)  :: ico,jco,kco
      integer(KIND=I4P)  :: ifi,jfi,kfi
   
      real(KIND=R8P)  ::  csi,eta,zeta
      real(KIND=R8P)  ::  umcsi,umeta,umzeta

      !==============================================================================
      !$omp parallel do &
      !$omp  default(none) &
      !$omp  shared(igr,self) &
      !$omp  private(ifi,jfi,kfi,ico,jco,kco, &
      !$omp          csi,umcsi,eta,umeta,zeta,umzeta)
      !==============================================================================
      do kfi=1,self%gr(igr)%nk
         kco = kfi/2+1
         zeta = dble(1+2*mod(kfi,2))/4.d0
         umzeta = 1.0d0 - zeta

         do jfi=1,self%gr(igr)%nj
            jco = jfi/2+1
            eta = dble(1+2*mod(jfi,2))/4.d0
            umeta = 1.0d0 - eta

            do ifi=1,self%gr(igr)%ni
               ico = ifi/2+1
               csi = dble(1+2*mod(ifi,2))/4.d0
               umcsi = 1.0d0 - csi

               self%gr(igr)%q(ifi,jfi,kfi) = &
                 umcsi*( ( self%gr(igr+1)%q(ico-1,jco-1,kco-1)*umeta &
                         + self%gr(igr+1)%q(ico-1,jco  ,kco-1)*  eta)*umzeta &
                       + ( self%gr(igr+1)%q(ico-1,jco-1,kco  )*umeta &
                         + self%gr(igr+1)%q(ico-1,jco  ,kco  )*  eta)*  zeta &
                       ) + &
                   csi*( ( self%gr(igr+1)%q(ico  ,jco-1,kco-1)*umeta &
                         + self%gr(igr+1)%q(ico  ,jco  ,kco-1)*  eta)*umzeta &
                       + ( self%gr(igr+1)%q(ico  ,jco-1,kco  )*umeta &
                         + self%gr(igr+1)%q(ico  ,jco  ,kco  )*  eta)*  zeta &
                       )
            end do
         end do
      end do
      !$omp end parallel do
      return
   end subroutine prolong

   module subroutine restrict(self,igr)
      implicit none
      class(mg_block),intent(in out) :: self         ! MG block
      integer(kind=I4P),intent(in)   :: igr          ! working grid level
   
      integer(KIND=I4P)  :: ico,jco,kco
      integer(KIND=I4P)  :: i,j,k
      integer(KIND=I4P)  :: grf
      grf = igr - 1
      !==============================================================================
      !$omp parallel do collapse(3) schedule(static) &
      !$omp  default(none) &
      !$omp  shared(self,igr,grf) &
      !$omp  private(ico,jco,kco,i,j,k)
      !==============================================================================
      do kco=1,self%gr(igr)%nk
         do jco=1,self%gr(igr)%nj
            do ico=1,self%gr(igr)%ni

               ! per il collapse(3) OMP
               k = 2*(kco-1)+1
               j = 2*(jco-1)+1
               i = 2*(ico-1)+1
   
               self%gr(igr)%q(ico,jco,kco) = &
                 self%gr(grf)%q(i  ,j  ,k  )*self%gr(grf)%cell(i  ,j  ,k  )%Vol &
               + self%gr(grf)%q(i  ,j+1,k  )*self%gr(grf)%cell(i  ,j+1,k  )%Vol &
               + self%gr(grf)%q(i  ,j  ,k+1)*self%gr(grf)%cell(i  ,j  ,k+1)%Vol &
               + self%gr(grf)%q(i  ,j+1,k+1)*self%gr(grf)%cell(i  ,j+1,k+1)%Vol &
               + self%gr(grf)%q(i+1,j  ,k  )*self%gr(grf)%cell(i+1,j  ,k  )%Vol &
               + self%gr(grf)%q(i+1,j+1,k  )*self%gr(grf)%cell(i+1,j+1,k  )%Vol &
               + self%gr(grf)%q(i+1,j  ,k+1)*self%gr(grf)%cell(i+1,j  ,k+1)%Vol &
               + self%gr(grf)%q(i+1,j+1,k+1)*self%gr(grf)%cell(i+1,j+1,k+1)%Vol 
   
               self%gr(igr)%q(ico,jco,kco) =   self%gr(igr)%q(ico,jco,kco) &
                                           /self%gr(igr)%cell(ico,jco,kco)%Vol
               ! store in qo
               self%gr(igr)%qo(ico,jco,kco) = self%gr(igr)%q(ico,jco,kco)
            end do
         end do
      end do
      !$omp end parallel do
      return
   end subroutine restrict

   module subroutine nullsorg(self,igr)
      implicit none
      class(mg_block),intent(in out) :: self         ! MG block
      integer(kind=I4P),intent(in)   :: igr          ! working grid level
   
      integer(KIND=I4P)  :: i,j,k
      !==============================================================================
      !$omp parallel do collapse(3) schedule(static) &
      !$omp  default(none) &
      !$omp  shared(self,igr) &
      !$omp  private(i,j,k)
      !==============================================================================
      do k=1,self%gr(igr)%nk
         do j=1,self%gr(igr)%nj
            do i=1,self%gr(igr)%ni
               self%gr(igr)%sq(i,j,k) = 0d0
            end do
         end do
      end do
      !$omp end parallel do
      return
   end subroutine nullsorg

   module subroutine collres(self,igr)
      implicit none
      class(mg_block),intent(in out) :: self         ! MG block
      integer(kind=I4P),intent(in)   :: igr          ! working grid level
   
      integer(KIND=I4P)  :: ico,jco,kco
      integer(KIND=I4P)  :: i,j,k
      integer(KIND=I4P)  :: grf
      grf = igr - 1
      !==============================================================================
      !$omp parallel do collapse(3) schedule(static) &
      !$omp  default(none) &
      !$omp  shared(self,igr,grf) &
      !$omp  private(ico,jco,kco,i,j,k)
      !==============================================================================
      do kco=1,self%gr(igr)%nk
         do jco=1,self%gr(igr)%nj
            do ico=1,self%gr(igr)%ni

               ! per il collapse(3) OMP
               k = 2*(kco-1)+1
               j = 2*(jco-1)+1
               i = 2*(ico-1)+1
   
               self%gr(igr)%sq(ico,jco,kco) = self%gr(grf)%dq(i  ,j  ,k  ) &
                                            + self%gr(grf)%dq(i  ,j+1,k  ) &
                                            + self%gr(grf)%dq(i  ,j  ,k+1) &
                                            + self%gr(grf)%dq(i  ,j+1,k+1) &
                                            + self%gr(grf)%dq(i+1,j  ,k  ) &
                                            + self%gr(grf)%dq(i+1,j+1,k  ) &
                                            + self%gr(grf)%dq(i+1,j  ,k+1) &
                                            + self%gr(grf)%dq(i+1,j+1,k+1) 
            end do
         end do
      end do
      !$omp end parallel do
      return
   end subroutine collres

   !=======================================================================
   !  calculation of variation on coarse grid
   !=======================================================================
   module subroutine calcvar(self,igr)
      implicit none
      class(mg_block),intent(in out) :: self         ! MG block
      integer(kind=I4P),intent(in)   :: igr          ! working grid level
   
      integer(KIND=I4P)  :: i,j,k
   !==============================================================================
   !$omp parallel  &
   !$omp  default(none) &
   !$omp  shared(self,igr) &
   !$omp  private(i,j,k)
   !==============================================================================
   !$omp do collapse(3) schedule(static)
      do k=1,self%gr(igr)%nk
         do j=1,self%gr(igr)%nj
            do i=1,self%gr(igr)%ni
               self%gr(igr)%dq(i,j,k) = self%gr(igr)%q(i,j,k) &
                                      -self%gr(igr)%qo(i,j,k)
            end do
         end do
      end do
   !$omp end do

   !----------------------------------------------------------
   ! simplified c.c.
   !----------------------------------------------------------

   !$omp do collapse(2) schedule(static)
      do k=1,self%gr(igr)%nk
         do j=1,self%gr(igr)%nj
            self%gr(igr)%dq(0,j,k) = self%gr(igr)%dq(1,j,k)
            self%gr(igr)%dq(self%gr(igr)%ni+1,j,k) = self%gr(igr)%dq(self%gr(igr)%ni,j,k)
         end do
      end do
   !$omp end do
   !$omp do collapse(2) schedule(static)
      do k=1,self%gr(igr)%nk
         do i=0,self%gr(igr)%ni+1
            self%gr(igr)%dq(i,0,k) = self%gr(igr)%dq(i,1,k)
            self%gr(igr)%dq(i,self%gr(igr)%nj+1,k) = self%gr(igr)%dq(i,self%gr(igr)%nj,k)
         end do
      end do
   !$omp end do

   !$omp do collapse(2) schedule(static)
      do j=0,self%gr(igr)%nj+1
         do i=0,self%gr(igr)%ni+1
            self%gr(igr)%dq(i,j,0) = self%gr(igr)%dq(i,j,1)
            self%gr(igr)%dq(i,j,self%gr(igr)%nk+1) = self%gr(igr)%dq(i,j,self%gr(igr)%nk)
         end do
      end do
   !$omp end do
   !$omp end parallel 
   
      return
   end subroutine calcvar

   module subroutine prolvar(self,igr,ls)
      implicit none
      class(mg_block),intent(in out) :: self         ! MG block
      integer(kind=I4P),intent(in)   :: igr          ! working grid level
      logical                        :: ls           ! level set control
   
      integer(KIND=I4P)  :: ico,jco,kco
      integer(KIND=I4P)  :: ifi,jfi,kfi
   
      real(KIND=R8P)  ::  csi,eta,zeta
      real(KIND=R8P)  ::  umcsi,umeta,umzeta

      if (ls) then
         !==============================================================================
         !$omp parallel do &
         !$omp  default(none) &
         !$omp  shared(igr,self) &
         !$omp  private(ifi,jfi,kfi,ico,jco,kco, &
         !$omp          csi,umcsi,eta,umeta,zeta,umzeta)
         !==============================================================================
         do kfi=1,self%gr(igr)%nk
            kco = kfi/2+1
            zeta = dble(1+2*mod(kfi,2))/4.d0
            umzeta = 1.0d0 - zeta
       
            do jfi=1,self%gr(igr)%nj
               jco = jfi/2+1
               eta = dble(1+2*mod(jfi,2))/4.d0
               umeta = 1.0d0 - eta
       
               do ifi=1,self%gr(igr)%ni
                  if (self%gr(igr)%vaux(ifi,jfi,kfi)%v(vaux_ls) > 0d0) cycle ! point in air

                  ico = ifi/2+1
                  csi = dble(1+2*mod(ifi,2))/4.d0
                  umcsi = 1.0d0 - csi
       
                  self%gr(igr)%q(ifi,jfi,kfi) =  self%gr(igr)%q(ifi,jfi,kfi) + &
                    umcsi*( ( self%gr(igr+1)%dq(ico-1,jco-1,kco-1)*umeta &
                            + self%gr(igr+1)%dq(ico-1,jco  ,kco-1)*  eta)*umzeta &
                          + ( self%gr(igr+1)%dq(ico-1,jco-1,kco  )*umeta &
                            + self%gr(igr+1)%dq(ico-1,jco  ,kco  )*  eta)*  zeta &
                          ) + &
                      csi*( ( self%gr(igr+1)%dq(ico  ,jco-1,kco-1)*umeta &
                            + self%gr(igr+1)%dq(ico  ,jco  ,kco-1)*  eta)*umzeta &
                          + ( self%gr(igr+1)%dq(ico  ,jco-1,kco  )*umeta &
                            + self%gr(igr+1)%dq(ico  ,jco  ,kco  )*  eta)*  zeta &
                          )
               end do
            end do
         end do
         !$omp end parallel do
      else
         !==============================================================================
         !$omp parallel do &
         !$omp  default(none) &
         !$omp  shared(igr,self) &
         !$omp  private(ifi,jfi,kfi,ico,jco,kco, &
         !$omp          csi,umcsi,eta,umeta,zeta,umzeta)
         !==============================================================================
         do kfi=1,self%gr(igr)%nk
            kco = kfi/2+1
            zeta = dble(1+2*mod(kfi,2))/4.d0
            umzeta = 1.0d0 - zeta
       
            do jfi=1,self%gr(igr)%nj
               jco = jfi/2+1
               eta = dble(1+2*mod(jfi,2))/4.d0
               umeta = 1.0d0 - eta
       
               do ifi=1,self%gr(igr)%ni
                  ico = ifi/2+1
                  csi = dble(1+2*mod(ifi,2))/4.d0
                  umcsi = 1.0d0 - csi
       
                  self%gr(igr)%q(ifi,jfi,kfi) =  self%gr(igr)%q(ifi,jfi,kfi) + &
                    umcsi*( ( self%gr(igr+1)%dq(ico-1,jco-1,kco-1)*umeta &
                            + self%gr(igr+1)%dq(ico-1,jco  ,kco-1)*  eta)*umzeta &
                          + ( self%gr(igr+1)%dq(ico-1,jco-1,kco  )*umeta &
                            + self%gr(igr+1)%dq(ico-1,jco  ,kco  )*  eta)*  zeta &
                          ) + &
                      csi*( ( self%gr(igr+1)%dq(ico  ,jco-1,kco-1)*umeta &
                            + self%gr(igr+1)%dq(ico  ,jco  ,kco-1)*  eta)*umzeta &
                          + ( self%gr(igr+1)%dq(ico  ,jco-1,kco  )*umeta &
                            + self%gr(igr+1)%dq(ico  ,jco  ,kco  )*  eta)*  zeta &
                          )
               end do
            end do
         end do
         !$omp end parallel do
      end if
      return
   end subroutine prolvar

end submodule MGblock_procedures
