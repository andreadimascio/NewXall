!=======================================================================
subroutine build_comm_map(igr)
!=======================================================================
! Costruzione mappe di comunicazione 
!=======================================================================
   use prec
   use work_var
   use mpi_var
   use mpi
   implicit none
    
   integer(kind=I4P),intent(in) :: igr

   logical           :: not_present

   integer(kind=I4P) :: i,j,k,n
   integer(kind=I4P) :: ipcc,npcc,tycc
   integer(kind=I4P) :: ibl,irank
   integer(kind=I4P) :: jbl,jrank
   integer(kind=I4P) :: ierr,iloc
   integer(kind=I4P) :: itag,nptin
   integer(kind=I4P) :: i_from(4)
   integer(kind=I4P),allocatable,dimension(:) :: ipt,npt,stts

   integer(kind=I4P),parameter :: chunk = 1000
   integer(kind=I4P),parameter :: npt0  = 1000

   if (allocated(stts)) deallocate(stts)
   allocate(stts(MPI_STATUS_SIZE),STAT=ierr)

   if (allocated(ipt)) deallocate(ipt)
   if (allocated(npt)) deallocate(npt)

   allocate(ipt(0:nproc-1),stat=ierr)
   allocate(npt(0:nproc-1),stat=ierr)

   ! allocazione temporanea mappe ricezione (allocazione struttura in input)
   do jrank=0,nproc-1
      ipt(jrank) = 0
      npt(jrank) = npt0
      call exchange(jrank,igr)%create_from(npt(jrank))
   end do

   ! loop con memorizzazione ordinata indici
   ! e costruzione mappe locali
   do ibl=1,nbl     ! Loop sui blocchi locali
      do k=-1,block(ibl)%gr(igr)%nk+2
      do j=-1,block(ibl)%gr(igr)%nj+2
      do i=-1,block(ibl)%gr(igr)%ni+2

         tycc = block(ibl)%gr(igr)%cell(i,j,k)%Chi%tycc

         if (tycc < 0) then  ! c.c. naturale --> costruzione mappa locale

            block(ibl)%gr(igr)%cell(i,j,k)%Chi%npcc = 1
            call block(ibl)%gr(igr)%cell(i,j,k)%Chi%add_X(1)
            ! costruzione donatori per le naturali 
            block(ibl)%gr(igr)%cell(i,j,k)%Chi%id(1) = i
            block(ibl)%gr(igr)%cell(i,j,k)%Chi%jd(1) = j
            block(ibl)%gr(igr)%cell(i,j,k)%Chi%kd(1) = k
            if (i<1) then
               block(ibl)%gr(igr)%cell(i,j,k)%Chi%bd(1) = 1  ! faccia
               block(ibl)%gr(igr)%cell(i,j,k)%Chi%id(1) = 1-i
            else if (i>block(ibl)%gr(igr)%ni) then
               block(ibl)%gr(igr)%cell(i,j,k)%Chi%bd(1) = 2
               block(ibl)%gr(igr)%cell(i,j,k)%Chi%id(1) = 2*block(ibl)%gr(igr)%ni-i+1
            else if (j<1) then
               block(ibl)%gr(igr)%cell(i,j,k)%Chi%bd(1) = 3
               block(ibl)%gr(igr)%cell(i,j,k)%Chi%jd(1) = 1-j
            else if (j>block(ibl)%gr(igr)%nj) then
               block(ibl)%gr(igr)%cell(i,j,k)%Chi%bd(1) = 4
               block(ibl)%gr(igr)%cell(i,j,k)%Chi%jd(1) = 2*block(ibl)%gr(igr)%nj-j+1
            else if (k<1) then
               block(ibl)%gr(igr)%cell(i,j,k)%Chi%bd(1) = 5
               block(ibl)%gr(igr)%cell(i,j,k)%Chi%kd(1) = 1-k
            else if (k>block(ibl)%gr(igr)%nk) then
               block(ibl)%gr(igr)%cell(i,j,k)%Chi%bd(1) = 6
               block(ibl)%gr(igr)%cell(i,j,k)%Chi%kd(1) = 2*block(ibl)%gr(igr)%nk-k+1
            else ! DEBUG
               write(*,*) ' ?????????? '
               stop
            end if

         else if (tycc > 0) then  ! punto chimera generico

            npcc = block(ibl)%gr(igr)%cell(i,j,k)%Chi%npcc
            if (npcc < 1) cycle  ! DEBUG inutile

            ! Loop sui donatori
            do ipcc=1,npcc
               jbl = block(ibl)%gr(igr)%cell(i,j,k)%Chi%bd(ipcc)! Blocco (globale) del donatore
               jrank = iblloc(jbl)                              ! Processo del blocco in questione
               
               if (jrank > 0) then ! Il donatore e' locale

                  block(ibl)%gr(igr)%cell(i,j,k)%Chi%bd(ipcc) = jrank  ! numerazione locale

               else ! Il blocco donatore e' non locale

                  jrank = -jrank    ! Processo del blocco in questione
                  i_from(1) = block(ibl)%gr(igr)%cell(i,j,k)%Chi%bd(ipcc)
                  i_from(2) = block(ibl)%gr(igr)%cell(i,j,k)%Chi%id(ipcc)
                  i_from(3) = block(ibl)%gr(igr)%cell(i,j,k)%Chi%jd(ipcc)
                  i_from(4) = block(ibl)%gr(igr)%cell(i,j,k)%Chi%kd(ipcc)

                  ! eliminazione duplicati
                  not_present = .true.
                  do n=1,ipt(jrank) 
                     if (i_from(1) /= exchange(jrank,igr)%i_from(n,1)) cycle
                     if (i_from(2) /= exchange(jrank,igr)%i_from(n,2)) cycle
                     if (i_from(3) /= exchange(jrank,igr)%i_from(n,3)) cycle
                     if (i_from(4) /= exchange(jrank,igr)%i_from(n,4)) cycle

                     ! già considerato
                     not_present = .false.
                     iloc = n  ! indice numerazione locale
                     exit
                  end do

                  if (not_present) then
                     ipt(jrank) = ipt(jrank) + 1 ! dimensioni correnti

                     call exchange(jrank,igr)%increase_size_from(ipt(jrank),chunk)
                     iloc = ipt(jrank)  ! indice numerazione locale

                     exchange(jrank,igr)%i_from(iloc,1) = i_from(1)
                     exchange(jrank,igr)%i_from(iloc,2) = i_from(2)
                     exchange(jrank,igr)%i_from(iloc,3) = i_from(3)
                     exchange(jrank,igr)%i_from(iloc,4) = i_from(4)
                  end if


                  block(ibl)%gr(igr)%cell(i,j,k)%Chi%bd(ipcc) = nbl + jrank + 1  ! blocco locale
                  block(ibl)%gr(igr)%cell(i,j,k)%Chi%id(ipcc) = iloc  ! indice numerazione locale
                  block(ibl)%gr(igr)%cell(i,j,k)%Chi%jd(ipcc) = 1
                  block(ibl)%gr(igr)%cell(i,j,k)%Chi%kd(ipcc) = 1

               end if
            end do

         end if

      end do
      end do
      end do
   end do

   do jrank=0,nproc-1 ! punti da ricevere
      npt(jrank) = ipt(jrank)
      exchange(jrank,igr)%npt_from = npt(jrank)
      ! riallocazione mappe in ricezione
      call exchange(jrank,igr)%reduce_size_from(npt(jrank))
   end do

   ! Trasmissione del numero di punti da spedire
   do jrank=0,nproc-1
      call MPI_SCATTER(npt(0)                    ,1,MPI_INTEGER, &
                       exchange(jrank,igr)%npt_to,1,MPI_INTEGER, &
                       jrank,MPI_COMM_WORLD,ierr)
   end do

   ! allocazione mappe di spedizione (allocazione struttura in input)
   do jrank=0,nproc-1
      call exchange(jrank,igr)%create_to(npt(jrank))
   end do

   !======================================================================
   ! Trasmissione degli indici dei punti
   !=========================================================================
   do jrank=0,nproc-1

      if (myrank == jrank) then    ! Ricezione da tutti ...

         do irank=0,nproc-1
           if(irank /= myrank) then    ! ...tranne che da myrank
             itag = nproc*(irank+1)
             nptin = 4*exchange(irank,igr)%npt_to
             call MPI_RECV(exchange(irank,igr)%i_to(1,1),nptin,MPI_INTEGER, &
                           irank,itag,MPI_COMM_WORLD,stts,ierr)
           end if
         end do

      else                            ! Spedizione al processo jrank

         itag = nproc*(myrank+1)
         nptin = 4*exchange(jrank,igr)%npt_from
         call MPI_SEND(exchange(jrank,igr)%i_from(1,1),nptin,MPI_INTEGER, &
                       jrank,itag,MPI_COMM_WORLD,ierr)

      end if

   end do

   do jrank=0,nproc-1
      ibl = nbl+jrank+1
      ! Creazione blocchi per connessione fra processori
      call block(ibl)%gr(igr)%create_BMP(exchange(jrank,igr)%npt_from)
      ! Localizzazione indici blocchi da spedire
      do iloc=1,exchange(jrank,igr)%npt_to

         jbl = exchange(jrank,igr)%i_to(iloc,1)
         jbl = iblloc(jbl)

         if (jbl <= 0) then  ! DEBUG
            write(*,*) ' ??????????????????'
            stop
         end if
         exchange(jrank,igr)%i_to(iloc,1) = jbl

      end do
   end do

   return
end subroutine build_comm_map
!=======================================================================

!=======================================================================
subroutine comm_var(igr)
!=======================================================================
! Trasferimento dati tra i processori 
!=======================================================================
   use prec
   use work_var
   use mpi_var
   use mpi 
   implicit none

   integer(kind=I4P),intent(in) :: igr

   integer(kind=I4P) :: ipt,npt
   integer(kind=I4P) :: irank,jrank
   integer(kind=I4P) :: bd,id,jd,kd
   integer(kind=I4P) :: ibl,ivar
   integer(kind=I4P) :: ierr,itag,icnt

   integer(kind=I4P),allocatable,dimension(:,:) ::stts_ar
   integer(kind=I4P),allocatable,dimension(:)   ::req
   real(kind=R8P),allocatable,dimension(:)   :: qq
      
   if (allocated(stts_ar)) deallocate(stts_ar)
   if (allocated(req    )) deallocate(req)
   if (allocated(qq     )) deallocate(qq)

   allocate(stts_ar(MPI_STATUS_SIZE,2*nproc),STAT=ierr)
   allocate(    req(                2*nproc),STAT=ierr)

   npt = 0
   do irank=0,nproc-1
      npt = max(exchange(irank,igr)%npt_from,npt)
      npt = max(exchange(irank,igr)%npt_to  ,npt)
   end do
   allocate(qq(npt),stat=ierr)

   do ivar=1,nvar
      !=======================================================================
      ! Ricezione da tutti i processori tranne che da se stesso
      !=======================================================================
      icnt = 0
      do irank=0,nproc-1            ! Ricezione da tutti i processori
   
        if (irank /= myrank) then    ! tranne che da se stesso
          itag = nproc*(irank+1)
          icnt = icnt+1
     
          ibl = nbl+irank+1
          npt = exchange(irank,igr)%npt_from 
   
          call MPI_IRECV(qq(1),npt,MPI_REAL8, &
                         irank,itag,MPI_COMM_WORLD,req(icnt),ierr)
          do ipt=1,npt
             block(ibl)%gr(igr)%q(ipt,1,1)%v(ivar) = qq(ipt)
          end do
        end if
   
      end do
   
      !=======================================================================
      ! Spedizione a tutti i processori tranne che a se stesso
      !=======================================================================
      do jrank=0,nproc-1    ! Spedizione a tutti i processori
         if(jrank /= myrank) then        ! tranne che a se stesso
   
             itag = nproc*(myrank+1)
             icnt = icnt+1
   
             npt = exchange(jrank,igr)%npt_to
             do ipt=1,npt
                bd = exchange(jrank,igr)%i_to(ipt,1)
                id = exchange(jrank,igr)%i_to(ipt,2)
                jd = exchange(jrank,igr)%i_to(ipt,3)
                kd = exchange(jrank,igr)%i_to(ipt,4)
                qq(ipt) = block(bd)%gr(igr)%q(id,jd,kd)%v(ivar)
             end do

             call MPI_ISEND(qq(1),npt,MPI_REAL8, &
                            jrank,itag,MPI_COMM_WORLD,req(icnt),ierr)
   
         end if
      end do
      !=======================================================================
      call MPI_WAITALL(icnt,req,stts_ar,ierr)
      !=======================================================================
   end do

   deallocate(stts_ar,stat=ierr)
   deallocate(    req,stat=ierr)

   return
end subroutine comm_var
!=======================================================================
