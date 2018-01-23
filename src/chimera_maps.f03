subroutine build_block_map
   use prec
   use chimera_type
   use work_var
   use mpi_var
   use mpi
   implicit none

   integer(kind=I4P) :: ibl,igr
   integer(kind=I4P) :: jbl
   integer(kind=I4P) :: jrank
   integer(kind=I4P) :: ni,nj,nk
   integer(kind=I4P) :: ib,jb,kb
   integer(kind=I4P) :: i,j,k,l
   integer(kind=I4P) :: nfrac_map
   integer(kind=I4P) :: ierr,itag

   integer(kind=I4P),dimension(3) :: dimblo 
   integer(kind=I4P),allocatable,dimension(:) :: stts
   integer(kind=I4P),allocatable,dimension(:,:,:,:) :: indici

   real(kind=R8P),dimension(3)    :: DX,Xhig,Xlow
   real(kind=R8P),dimension(6)    :: block_size

   ! allocazione mappe per tutti i nbltot blocchi
   if (allocated(mappa)) deallocate(mappa)
   allocate(mappa(nbltot),stat=ierr)

   nfrac_map = 8  ! leggere da input

   igr = 1  ! mappa i blocchi solo sulla più fitta

   !****************************************************
   ! costruzione mappe per i blocchi locali
   !****************************************************
   do ibl=1,nbl
      jbl = iblglo(ibl) ! blocco globale

      ni = max(1,block(ibl)%gr(igr)%ni/nfrac_map)
      nj = max(1,block(ibl)%gr(igr)%nj/nfrac_map)
      nk = max(1,block(ibl)%gr(igr)%nk/nfrac_map)

      ! allocazione mappa per blocchi locali
      call mappa(jbl)%create_M(ni,nj,nk)

      ! estremi del blocco
      mappa(jbl)%low =  Infinity
      mappa(jbl)%hig = -Infinity
      do k=0,block(ibl)%gr(igr)%nk
      do j=0,block(ibl)%gr(igr)%nj
      do i=0,block(ibl)%gr(igr)%ni,block(ibl)%gr(igr)%ni
         do l=1,3
            mappa(jbl)%low(l) = min(mappa(jbl)%low(l),block(ibl)%gr(igr)%cell(i,j,k)%Vertex(l))
            mappa(jbl)%hig(l) = max(mappa(jbl)%hig(l),block(ibl)%gr(igr)%cell(i,j,k)%Vertex(l))
         end do
      end do
      end do
      end do
      do k=0,block(ibl)%gr(igr)%nk
      do j=0,block(ibl)%gr(igr)%nj,block(ibl)%gr(igr)%nj
      do i=0,block(ibl)%gr(igr)%ni
         do l=1,3
            mappa(jbl)%low(l) = min(mappa(jbl)%low(l),block(ibl)%gr(igr)%cell(i,j,k)%Vertex(l))
            mappa(jbl)%hig(l) = max(mappa(jbl)%hig(l),block(ibl)%gr(igr)%cell(i,j,k)%Vertex(l))
         end do
      end do
      end do
      end do
      do k=0,block(ibl)%gr(igr)%nk,block(ibl)%gr(igr)%nk
      do j=0,block(ibl)%gr(igr)%nj
      do i=0,block(ibl)%gr(igr)%ni
         do l=1,3
            mappa(jbl)%low(l) = min(mappa(jbl)%low(l),block(ibl)%gr(igr)%cell(i,j,k)%Vertex(l))
            mappa(jbl)%hig(l) = max(mappa(jbl)%hig(l),block(ibl)%gr(igr)%cell(i,j,k)%Vertex(l))
         end do
      end do
      end do
      end do

      ! parcellizzazione
      DX = mappa(jbl)%hig - mappa(jbl)%low
      DX(1) = DX(1)/float(ni)
      DX(2) = DX(2)/float(nj)
      DX(3) = DX(3)/float(nk)
      do kb=1,nk
      do jb=1,nj
      do ib=1,ni

         ! estremi parallelepipedo
         Xhig(1) = float(ib)*DX(1)
         Xhig(2) = float(jb)*DX(2)
         Xhig(3) = float(kb)*DX(3)

         Xhig    = Xhig + mappa(jbl)%low  ! + offset
         Xlow    = Xhig - DX

         mappa(jbl)%imax(ib,jb,kb) = 0
         mappa(jbl)%jmax(ib,jb,kb) = 0
         mappa(jbl)%kmax(ib,jb,kb) = 0

         mappa(jbl)%imin(ib,jb,kb) = ni
         mappa(jbl)%jmin(ib,jb,kb) = nj
         mappa(jbl)%kmin(ib,jb,kb) = nk


         do k=0,block(ibl)%gr(igr)%nk
         do j=0,block(ibl)%gr(igr)%nj
         do i=0,block(ibl)%gr(igr)%ni

            coord: do l=1,3
               if (block(ibl)%gr(igr)%cell(i,j,k)%Vertex(l) > Xhig(l) ) exit coord
               if (block(ibl)%gr(igr)%cell(i,j,k)%Vertex(l) < Xlow(l) ) exit coord

               ! il punto cade nel parallelepipedo
               mappa(jbl)%imax(ib,jb,kb) = max(mappa(jbl)%imax(ib,jb,kb),i)
               mappa(jbl)%jmax(ib,jb,kb) = max(mappa(jbl)%jmax(ib,jb,kb),j)
               mappa(jbl)%kmax(ib,jb,kb) = max(mappa(jbl)%kmax(ib,jb,kb),k)

               mappa(jbl)%imin(ib,jb,kb) = min(mappa(jbl)%imin(ib,jb,kb),i)
               mappa(jbl)%jmin(ib,jb,kb) = min(mappa(jbl)%jmin(ib,jb,kb),j)
               mappa(jbl)%kmin(ib,jb,kb) = min(mappa(jbl)%kmin(ib,jb,kb),k)

            end do coord
            
         end do
         end do
         end do

      end do
      end do
      end do
   end do

   if (nproc == 1) return
   ! matrici per mpi
   if (allocated(stts)) deallocate(stts)
   allocate(stts(MPI_STATUS_SIZE),STAT=ierr)

   !****************************************************
   ! trasmissione mappe - ricezione in attesa
   !****************************************************
   do jbl=1,nbltot  ! blocco globale
      jrank = iblloc(jbl)  ! blocco locale

      if (jrank < 0 ) then  
         jrank = - jrank ! processo del blocco remoto

         ! dimensioni blocco dal processo jrank
         itag = jbl
         call MPI_RECV(dimblo,3,MPI_INTEGER, &
                       jrank,itag,MPI_COMM_WORLD,stts,ierr)

         ! estremi del blocco
         itag = itag + nbltot
         call MPI_RECV(block_size,6,MPI_REAL8, &
                       jrank,itag,MPI_COMM_WORLD,stts,ierr)
         mappa(jbl)%low = block_size(1:3)
         mappa(jbl)%hig = block_size(4:6)

         ! allocazione mappa per blocchi remoti
         call mappa(jbl)%create_M(dimblo(1),dimblo(2),dimblo(3))
      end if

   end do

   !****************************************************
   ! trasmissione mappe - spedizione
   !****************************************************
   do jbl=1,nbltot  ! blocco globale
      jrank = iblloc(jbl)  ! processo remoto
      if (jrank > 0) then
         jrank = myrank  ! blocco locale
      else
         jrank = -jrank
      end if

      if (jrank /= myrank) then 

         ! dimensioni blocco al processo jrank
         itag = jbl
         dimblo(1) = mappa(jbl)%ni
         dimblo(2) = mappa(jbl)%nj
         dimblo(3) = mappa(jbl)%nk
         call MPI_SEND(dimblo,3,MPI_INTEGER, &
                       jrank,itag,MPI_COMM_WORLD,ierr)

         ! estremi del blocco
         itag = itag + nbltot
         block_size(1:3) = mappa(jbl)%low
         block_size(4:6) = mappa(jbl)%hig
         call MPI_SEND(block_size,6,MPI_REAL8, &
                       jrank,itag,MPI_COMM_WORLD,ierr)
      end if
   end do

   !***********************************************************************
   ! le mappe sono dimensionate e allocate --> scambio mappe
   !***********************************************************************
   ! trasmissione mappe - ricezione in attesa
   do jbl=1,nbltot  ! blocco globale
      jrank = iblloc(jbl)  ! blocco locale

      if (jrank < 0 ) then  
         jrank = - jrank ! processo del blocco remoto

         ni = mappa(jbl)%ni
         nj = mappa(jbl)%nj
         nk = mappa(jbl)%nk
         if (allocated(indici)) deallocate(indici)
         allocate(indici(ni,nj,nk,6))

         ! mappa dal processo jrank
         itag = jbl
         call MPI_RECV(indici,(ni*nj*nk*6),MPI_INTEGER, &
                       jrank,itag,MPI_COMM_WORLD,stts,ierr)
         mappa(jbl)%imin = indici(1:ni,1:nj,1:nk,1)
         mappa(jbl)%imax = indici(1:ni,1:nj,1:nk,2)
         mappa(jbl)%jmin = indici(1:ni,1:nj,1:nk,3)
         mappa(jbl)%jmax = indici(1:ni,1:nj,1:nk,4)
         mappa(jbl)%kmin = indici(1:ni,1:nj,1:nk,5)
         mappa(jbl)%kmax = indici(1:ni,1:nj,1:nk,6)
      end if

   end do

   ! spedizione
   do jbl=1,nbltot  ! blocco globale
      jrank = iblloc(jbl)  ! processo remoto
      if (jrank > 0) then
         jrank = myrank  ! blocco locale
      else
         jrank = -jrank
      end if

      if (jrank /= myrank) then 
         ni = mappa(jbl)%ni
         nj = mappa(jbl)%nj
         nk = mappa(jbl)%nk
         if (allocated(indici)) deallocate(indici)
         allocate(indici(ni,nj,nk,6))

         ! mappa al processo jrank
         itag = jbl
         indici(1:ni,1:nj,1:nk,1) =  mappa(jbl)%imin
         indici(1:ni,1:nj,1:nk,2) =  mappa(jbl)%imax
         indici(1:ni,1:nj,1:nk,3) =  mappa(jbl)%jmin
         indici(1:ni,1:nj,1:nk,4) =  mappa(jbl)%jmax
         indici(1:ni,1:nj,1:nk,5) =  mappa(jbl)%kmin
         indici(1:ni,1:nj,1:nk,6) =  mappa(jbl)%kmax

         call MPI_SEND(indici,(ni*nj*nk*6),MPI_INTEGER, &
                       jrank,itag,MPI_COMM_WORLD,ierr)
      end if
   end do
   !*************************************************
   ! mappe complete
   !*************************************************
   return
end subroutine build_block_map
!=======================================================================
