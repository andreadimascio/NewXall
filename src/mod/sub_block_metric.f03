submodule (block_type) block_metric
   use cell_type
   use block_type
contains
   module subroutine centroid(self)
      class(gen_block),intent(in out) :: self         ! Blocco
      integer(kind=I4P)             :: i,j,k
      !==============================================================================
      !$omp parallel do collapse(3) schedule(static) &
      !$omp default(none) &
      !$omp shared(self) &
      !$omp private(i,j,k)
      !==============================================================================
      do k=-1,self%nk+2
         do j=-1,self%nj+2
            do i=-1,self%ni+2
               self%cell(i,j,k)%Center = self%cell(i  ,j  ,k  )%Vertex &
                                       + self%cell(i-1,j  ,k  )%Vertex &
                                       + self%cell(i  ,j-1,k  )%Vertex &
                                       + self%cell(i-1,j-1,k  )%Vertex &
                                       + self%cell(i  ,j  ,k-1)%Vertex &
                                       + self%cell(i-1,j  ,k-1)%Vertex &
                                       + self%cell(i  ,j-1,k-1)%Vertex &
                                       + self%cell(i-1,j-1,k-1)%Vertex 
               self%cell(i,j,k)%Center = self%cell(i  ,j  ,k  )%Center*0.125d0
            end do
         end do
      end do
      !$omp end parallel do ! nowait
   end subroutine centroid

   module subroutine faces(self)
      class(gen_block),intent(in out) :: self         ! Blocco
      integer(kind=I4P)             :: i,j,k
      real(kind=R8P),dimension(3)   :: P1,P2,P3,sign_S

      ! TEST ORIENTAMENTO
      i = max(self%ni/2,1)
      j = max(self%nj/2,1)
      k = max(self%nk/2,1)

      ! direzione i
      P1 = self%cell(i,j,k)%Vertex - self%cell(i,j-1,k-1)%Vertex
      P2 = self%cell(i,j-1,k)%Vertex - self%cell(i,j,k-1)%Vertex
      P3 = P1.cross.P2

      P1 = self%cell(i,j,k)%Vertex - self%cell(i-1,j,k)%Vertex
      sign_S(1) = sign(1d0,P1.dot.P3)  ! orientamento normale

      ! direzione j
      P1 = self%cell(i,j,k)%Vertex - self%cell(i-1,j,k-1)%Vertex
      P2 = self%cell(i-1,j,k)%Vertex - self%cell(i,j,k-1)%Vertex
      P3 = P1.cross.P2

      P1 = self%cell(i,j,k)%Vertex - self%cell(i,j-1,k)%Vertex
      sign_S(2) = sign(1d0,P1.dot.P3)  ! orientamento normale

      ! direzione k
      P1 = self%cell(i,j,k)%Vertex - self%cell(i-1,j-1,k)%Vertex
      P2 = self%cell(i-1,j,k)%Vertex - self%cell(i,j-1,k)%Vertex
      P3 = P1.cross.P2

      P1 = self%cell(i,j,k)%Vertex - self%cell(i,j,k-1)%Vertex
      sign_S(3) = sign(1d0,P1.dot.P3)  ! orientamento normale

      ! N.B.
      sign_S = 0.5d0*sign_S

      !==============================================================================
      !$omp parallel  &
      !$omp default(none) &
      !$omp shared(self,sign_S) &
      !$omp private(i,j,k,P1,P2,P3)
      !==============================================================================

      ! direzione i

      !$omp do collapse(3) schedule(static)
      do k=-1,self%nk+2
         do j=-1,self%nj+2
            do i=-2,self%ni+2
               P1 = self%cell(i,j,k)%Vertex - self%cell(i,j-1,k-1)%Vertex
               P2 = self%cell(i,j-1,k)%Vertex - self%cell(i,j,k-1)%Vertex
               P3 = P1.cross.P2

               self%cell(i,j,k)%Sn(1:3,1) = sign_S(1)*P3
               self%cell(i,j,k)%Area(1)   = norm2(P3)
            end do
         end do
      end do
      !$omp end do 

      ! direzione j

      !$omp do collapse(3) schedule(static)
      do k=-1,self%nk+2
         do j=-2,self%nj+2
            do i=-1,self%ni+2
               P1 = self%cell(i,j,k)%Vertex - self%cell(i-1,j,k-1)%Vertex
               P2 = self%cell(i-1,j,k)%Vertex - self%cell(i,j,k-1)%Vertex
               P3 = P1.cross.P2

               self%cell(i,j,k)%Sn(1:3,2) = sign_S(2)*P3
               self%cell(i,j,k)%Area(2)   = norm2(P3)
            end do
         end do
      end do
      !$omp end do 

      ! direzione k

      !$omp do collapse(3) schedule(static)
      do k=-2,self%nk+2
         do j=-1,self%nj+2
            do i=-1,self%ni+2
               P1 = self%cell(i,j,k)%Vertex - self%cell(i-1,j-1,k)%Vertex
               P2 = self%cell(i-1,j,k)%Vertex - self%cell(i,j-1,k)%Vertex
               P3 = P1.cross.P2

               self%cell(i,j,k)%Sn(1:3,3) = sign_S(3)*P3
               self%cell(i,j,k)%Area(3)   = norm2(P3)
            end do
         end do
      end do
      !$omp end do 
      !$omp end parallel 
   end subroutine faces

   module subroutine volume(self)
      class(gen_block),intent(in out) :: self         ! Blocco
      integer(kind=I4P)             :: i,j,k
      real(kind=R8P),dimension(3)   :: P1
      !==============================================================================
      !$omp parallel do collapse(3) schedule(static) &
      !$omp default(none) &
      !$omp shared(self) &
      !$omp private(i,j,k,P1)
      !==============================================================================
      do k=-1,self%nk+2
         do j=-1,self%nj+2
            do i=-1,self%ni+2
               ! direzione i
               P1 = self%cell(i  ,j,k  )%Vertex+self%cell(i  ,j-1,k  )%Vertex &
                  + self%cell(i  ,j,k-1)%Vertex+self%cell(i  ,j-1,k-1)%Vertex
               self%cell(i,j,k)%Vol =                          P1.dot.self%cell(i,j,k)%Sn(1:3,1)
                                                                                               
               P1 = self%cell(i-1,j,k  )%Vertex+self%cell(i-1,j-1,k  )%Vertex &                
                  + self%cell(i-1,j,k-1)%Vertex+self%cell(i-1,j-1,k-1)%Vertex                  
               self%cell(i,j,k)%Vol = self%cell(i,j,k)%Vol - P1.dot.self%cell(i-1,j,k)%Sn(1:3,1)

               ! direzione j
               P1 = self%cell(i  ,j  ,k  )%Vertex+self%cell(i-1,j  ,k  )%Vertex &
                  + self%cell(i  ,j  ,k-1)%Vertex+self%cell(i-1,j  ,k-1)%Vertex
               self%cell(i,j,k)%Vol = self%cell(i,j,k)%Vol +   P1.dot.self%cell(i,j,k)%Sn(1:3,2)
                                                                                               
               P1 = self%cell(i  ,j-1,k  )%Vertex+self%cell(i-1,j-1,k  )%Vertex &              
                  + self%cell(i  ,j-1,k-1)%Vertex+self%cell(i-1,j-1,k-1)%Vertex                
               self%cell(i,j,k)%Vol = self%cell(i,j,k)%Vol - P1.dot.self%cell(i,j-1,k)%Sn(1:3,2)

               ! direzione k
               P1 = self%cell(i  ,j  ,k  )%Vertex+self%cell(i-1,j  ,k  )%Vertex &
                  + self%cell(i  ,j-1,k  )%Vertex+self%cell(i-1,j-1,k  )%Vertex
               self%cell(i,j,k)%Vol = self%cell(i,j,k)%Vol +   P1.dot.self%cell(i,j,k)%Sn(1:3,3)
                                                                                               
               P1 = self%cell(i  ,j  ,k-1)%Vertex+self%cell(i-1,j  ,k-1)%Vertex &              
                  + self%cell(i  ,j-1,k-1)%Vertex+self%cell(i-1,j-1,k-1)%Vertex                
               self%cell(i,j,k)%Vol = self%cell(i,j,k)%Vol - P1.dot.self%cell(i,j,k-1)%Sn(1:3,3)

               self%cell(i,j,k)%Vol = self%cell(i,j,k)%Vol / 12d0 ! sarebbe /4 /3
            end do
         end do
      end do
      !$omp end parallel do ! nowait
   end subroutine volume

   module subroutine tens(self)
      class(gen_block),intent(in out) :: self         ! Blocco
      integer(kind=I4P)             :: i,j,k
      !==============================================================================
      !$omp parallel do collapse(3) schedule(static) &
      !$omp default(none) &
      !$omp shared(self) &
      !$omp private(i,j,k)
      !==============================================================================
      do k=-1,self%nk+2
         do j=-1,self%nj+2
            do i=-1,self%ni+2
               self%cell(i,j,k)%Csi(1:3,1) = self%cell(i,j,k)%Sn(1:3,1) + self%cell(i-1,j,k)%Sn(1:3,1)
               self%cell(i,j,k)%Csi(1:3,2) = self%cell(i,j,k)%Sn(1:3,2) + self%cell(i,j-1,k)%Sn(1:3,2)
               self%cell(i,j,k)%Csi(1:3,3) = self%cell(i,j,k)%Sn(1:3,3) + self%cell(i,j,k-1)%Sn(1:3,3)
               self%cell(i,j,k)%Csi = 0.5d0*self%cell(i,j,k)%Csi/self%cell(i,j,k)%Vol
            end do
         end do
      end do
      !$omp end parallel do ! nowait
   end subroutine tens
end submodule block_metric
