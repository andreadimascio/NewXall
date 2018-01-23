!=======================================================================
subroutine search_principal_donor(igr)
   use prec
   use chimera_type
   use concon, only: Search
   use work_var
   use mpi_var
   use mpi

   implicit none
   integer(kind=I4P),intent(in) :: igr

   integer(kind=I4P) :: ibl,jbl,kbl
   integer(kind=I4P) :: i,j,k
   integer(kind=I4P) :: ib,jb,kb
   integer(kind=I4P) :: id,jd,kd
   integer(kind=I4P) :: ni,nj,nk
   integer(kind=I4P) :: jrank,tycc
   integer(kind=I4P) :: iloc
   integer(kind=I4P) :: nptmax
   integer(kind=I4P) :: itag,ierr
   integer(kind=I4P),allocatable,dimension(:) :: npt
   integer(kind=I4P),allocatable,dimension(:) :: stts
   integer(kind=I4P),allocatable,dimension(:,:) :: iblo
   integer(kind=I4P),parameter :: chunk = 1000
   integer(kind=I4P),parameter :: npt0  = 1000

   real(kind=R8P)              :: dsize
   real(kind=R8P),dimension(3) :: xc,DX
   real(kind=R8P),allocatable,dimension(:,:) :: rblo

   ! allocazione strutture
   if (allocated(npt)) deallocate(npt)
   allocate(npt(0:nproc-1))

   if (allocated(Xpt_req)) deallocate(Xpt_req)
   allocate(Xpt_req(0:nproc-1))
   if (allocated(Xpt_send)) deallocate(Xpt_send)
   allocate(Xpt_send(0:nproc-1))
   ! allocazione strutture  - fine

   do jrank=0,nproc-1  ! dimensionamento iniziale
      call Xpt_req(jrank)%create_XA(npt0)
   end do
   
   npt = 0

   loc_blocks: do ibl=1,nbl ! loop blocchi locali
      kbl = iblglo(ibl) ! numerazione globale

      ! ricerca possibili donatori per tutti i punti
      glo_blocks: do jbl=1,nbltot

         if (kbl == jbl) cycle glo_blocks! il blocco è lo stesso
         if (mappa(kbl)%priority > mappa(jbl)%priority) cycle glo_blocks ! priorità jbl più bassa

         if (mappa(kbl)%hig .smaller. mappa(jbl)%low) cycle glo_blocks
         if (mappa(kbl)%low .larger.  mappa(jbl)%hig) cycle glo_blocks

         ! c'e' possibile sovrapposizione
         jrank = iblloc(jbl) ! rango del blocco
         if (jrank > 0) then
            jrank = myrank  ! blocco locale
         else
            jrank = -jrank
         end if

         ni = mappa(jbl)%ni
         nj = mappa(jbl)%nj
         nk = mappa(jbl)%nk
         
         ! parcellizzazione
         DX = mappa(jbl)%hig - mappa(jbl)%low
         DX(1) = DX(1)/float(ni)
         DX(2) = DX(2)/float(nj)
         DX(3) = DX(3)/float(nk)

         do k=-1,block(ibl)%gr(igr)%nk+2
         do j=-1,block(ibl)%gr(igr)%nj+2
         do i=-1,block(ibl)%gr(igr)%ni+2


            xc = block(ibl)%gr(igr)%cell(i,j,k)%Center ! punto in analisi

            call check_patch(kbl,ibl,igr,i,j,k, &
                             id,jd,kd,tycc)
            if (tycc < 0) cycle  ! c.c. naturale

            ! sottoblocco donatore
            ib = ceiling( (xc(1)-mappa(jbl)%low(1))/DX(1))
            if (ib < 1) cycle ; if (ib > ni) cycle

            jb = ceiling( (xc(2)-mappa(jbl)%low(2))/DX(2))
            if (jb < 1) cycle ; if (jb > nj) cycle

            kb = ceiling( (xc(3)-mappa(jbl)%low(3))/DX(3))
            if (kb < 1) cycle ; if (kb > nk) cycle

            ! è un possibile buco: donatore su jrank
            npt(jrank) = npt(jrank) + 1

            iloc = npt(jrank)
            call Xpt_req(jrank)%increase_size_XA(iloc,chunk)

            Xpt_req(jrank)%P(iloc)%xc = xc
            ! in ingresso : distanza dalla (eventuale) parete
            Xpt_req(jrank)%P(iloc)%dsize = block(ibl)%gr(igr)%cell(i,j,k)%dwall

            Xpt_req(jrank)%P(iloc)%bd = jbl  ! in andata: circoscrizione donatore
            Xpt_req(jrank)%P(iloc)%id = ib
            Xpt_req(jrank)%P(iloc)%jd = jb
            Xpt_req(jrank)%P(iloc)%kd = kb

            Xpt_req(jrank)%P(iloc)%b  = ibl  ! identificativo locale
            Xpt_req(jrank)%P(iloc)%i  = i 
            Xpt_req(jrank)%P(iloc)%j  = j 
            Xpt_req(jrank)%P(iloc)%k  = k 

            ! crea chimera temporanea con un solo donatore 
            call block(ibl)%gr(igr)%cell(i,j,k)%Chi%add_X(1)
            block(ibl)%gr(igr)%cell(i,j,k)%Chi%tycc  = Search
            block(ibl)%gr(igr)%cell(i,j,k)%Chi%we(1) = zero ! test dimensioni

         end do
         end do
         end do

      end do glo_blocks
   end do loc_blocks

   do jrank=0,nproc-1  ! allocazione esatta
      call Xpt_req(jrank)%reduce_size_XA(npt(jrank))
   end do

   !************************************************************
   ! scambio informazioni
   !************************************************************
   if (allocated(stts)) deallocate(stts)
   allocate(stts(MPI_STATUS_SIZE),STAT=ierr)

   ! trasmissione numero di punti
   do jrank=0,nproc-1
      call MPI_SCATTER(npt(0)           ,1,MPI_INTEGER, &
                       Xpt_send(jrank)%n,1,MPI_INTEGER, &
                       jrank,MPI_COMM_WORLD,ierr)
   end do

   nptmax = 0
   do jrank=0,nproc-1  ! dimensionamento buffet di spedizione
      call Xpt_send(jrank)%create_XA(npt(jrank))
      nptmax = max(nptmax,npt(jrank))
   end do

   ! matrici ausiliarie
   if (allocated(iblo)) deallocate(iblo)
   allocate(iblo(nptmax,4))
   if (allocated(rblo)) deallocate(rblo)
   allocate(rblo(nptmax,4))

   ! ricezione punti per i quali bisogna cercare un donatore
   ! sul processo myrank e rispedirli al mittente
   do jrank=0,nproc-1
      if (jrank == myrank) then
         ! non c'è spedizione se la ricerca è locale
         Xpt_send(jrank)%P(1:npt(jrank)) = Xpt_req(myrank)%P(1:npt(jrank))
      else
         itag = jbl
         call MPI_RECV(iblo,4*npt(jrank),MPI_INTEGER, &
                       jrank,itag,MPI_COMM_WORLD,stts,ierr)
         Xpt_send(jrank)%P(1:npt(jrank))%bd = iblo(1:npt(jrank),1)
         Xpt_send(jrank)%P(1:npt(jrank))%id = iblo(1:npt(jrank),2)
         Xpt_send(jrank)%P(1:npt(jrank))%jd = iblo(1:npt(jrank),3)
         Xpt_send(jrank)%P(1:npt(jrank))%kd = iblo(1:npt(jrank),4)

         itag = jbl+nbltot
         call MPI_RECV(rblo,4*npt(jrank),MPI_INTEGER, &
                       jrank,itag,MPI_COMM_WORLD,stts,ierr)
         Xpt_send(jrank)%P(1:npt(jrank))%xc(1) = rblo(1:npt(jrank),1)
         Xpt_send(jrank)%P(1:npt(jrank))%xc(2) = rblo(1:npt(jrank),2)
         Xpt_send(jrank)%P(1:npt(jrank))%xc(3) = rblo(1:npt(jrank),3)
         Xpt_send(jrank)%P(1:npt(jrank))%dsize = rblo(1:npt(jrank),4)

      end if
   end do

   ! spedizione punti per i quali bisogna ricevere un donatore
   ! dal processo jrank
   do jrank=0,nproc-1
      if (jrank /= myrank) then
         itag = jbl
         iblo(1:npt(jrank),1)  = Xpt_req(jrank)%P(1:npt(jrank))%bd 
         iblo(1:npt(jrank),2)  = Xpt_req(jrank)%P(1:npt(jrank))%id 
         iblo(1:npt(jrank),3)  = Xpt_req(jrank)%P(1:npt(jrank))%jd 
         iblo(1:npt(jrank),4)  = Xpt_req(jrank)%P(1:npt(jrank))%kd 
         call MPI_SEND(iblo,4*npt(jrank),MPI_INTEGER, &
                       jrank,itag,MPI_COMM_WORLD,ierr)

         itag = jbl+nbltot
         rblo(1:npt(jrank),1) = Xpt_req(jrank)%P(1:npt(jrank))%xc(1)
         rblo(1:npt(jrank),2) = Xpt_req(jrank)%P(1:npt(jrank))%xc(2)
         rblo(1:npt(jrank),3) = Xpt_req(jrank)%P(1:npt(jrank))%xc(3)
         rblo(1:npt(jrank),4) = Xpt_req(jrank)%P(1:npt(jrank))%dsize
         call MPI_SEND(rblo,4*npt(jrank),MPI_INTEGER, &
                       jrank,itag,MPI_COMM_WORLD,ierr)

      end if
   end do

! ***************************************************
! ricerca locale donatori
! ***************************************************

   do jrank=0,nproc-1  ! punti da rispedire a jrank
      do iloc = 1,Xpt_send(jrank)%n

         xc  = Xpt_send(jrank)%P(iloc)%xc  ! candidato

         jbl = Xpt_send(jrank)%P(iloc)%bd  ! blocco donatore
         ibl = iblloc(jbl) 
         if (ibl < 0 ) then ! DEBUG
            write(*,*) ' ???????????? '
            stop
         end if

         id = Xpt_send(jrank)%P(iloc)%id  ! in ingresso : sottoblocco
         jd = Xpt_send(jrank)%P(iloc)%jd
         kd = Xpt_send(jrank)%P(iloc)%kd

         ! in ingresso: distanza dalla parete propria
         dsize = Xpt_send(jrank)%P(iloc)%dsize
         call find_donor(xc,mappa(jbl),block(ibl)%gr(igr), &
                         id,jd,kd,tycc,dsize)

         Xpt_send(jrank)%P(iloc)%id = id ! in uscita : indici donatore
         Xpt_send(jrank)%P(iloc)%jd = jd
         Xpt_send(jrank)%P(iloc)%kd = kd
         Xpt_send(jrank)%P(iloc)%kd = tycc  ! tipo cella (se trovato)
         Xpt_send(jrank)%P(iloc)%dsize = dsize
                                               
      end do
   end do

   !************************************************************
   ! scambio inverso di informazioni
   !************************************************************
   ! matrici ausiliarie
   !  if (allocated(iblo)) deallocate(iblo)
   !  allocate(iblo(nptmax,4))  ha le stesse dimensioni di prima
   if (allocated(rblo)) deallocate(rblo)
   allocate(rblo(nptmax,1))

   ! ricezione punti per i quali è stato cercato 
   ! un possibile donatore su jrank
   do jrank=0,nproc-1
      if (jrank == myrank) then
         Xpt_req(jrank)%P(1:npt(jrank)) = Xpt_send(myrank)%P(1:npt(jrank))
      else
         itag = jbl
         call MPI_RECV(iblo,4*npt(jrank),MPI_INTEGER, &
                       jrank,itag,MPI_COMM_WORLD,stts,ierr)
         Xpt_req(jrank)%P(1:npt(jrank))%id = iblo(1:npt(jrank),1)
         Xpt_req(jrank)%P(1:npt(jrank))%jd = iblo(1:npt(jrank),2)
         Xpt_req(jrank)%P(1:npt(jrank))%kd = iblo(1:npt(jrank),3)
         Xpt_req(jrank)%P(1:npt(jrank))%tycc = iblo(1:npt(jrank),4)

         itag = jbl+nbltot
         call MPI_RECV(rblo,npt(jrank),MPI_INTEGER, &
                       jrank,itag,MPI_COMM_WORLD,stts,ierr)
         Xpt_req(jrank)%P(1:npt(jrank))%dsize = rblo(1:npt(jrank),1)
      end if
   end do

   ! spedizione punti per i quali è stato cercato un donatore
   ! per il processo jrank
   do jrank=0,nproc-1
      if (jrank /= myrank) then
         itag = jbl
         iblo(1:npt(jrank),1)  = Xpt_req(jrank)%P(1:npt(jrank))%id 
         iblo(1:npt(jrank),2)  = Xpt_req(jrank)%P(1:npt(jrank))%jd 
         iblo(1:npt(jrank),3)  = Xpt_req(jrank)%P(1:npt(jrank))%kd 
         iblo(1:npt(jrank),4)  = Xpt_req(jrank)%P(1:npt(jrank))%tycc
         call MPI_SEND(iblo,4*npt(jrank),MPI_INTEGER, &
                       jrank,itag,MPI_COMM_WORLD,ierr)

         itag = jbl+nbltot
         rblo(1:npt(jrank),1) = Xpt_req(jrank)%P(1:npt(jrank))%dsize
         call MPI_SEND(rblo,npt(jrank),MPI_INTEGER, &
                       jrank,itag,MPI_COMM_WORLD,ierr)

      end if
   end do

   deallocate(rblo)
   deallocate(iblo)
   deallocate(stts)

   ! costruzione celle chimera
   ! analisi dei possibili donatori ricevuti
   do jrank=0,nproc-1
      npt(jrank) = Xpt_req(jrank)%n
      do iloc=1,npt(jrank)

         jbl = Xpt_req(jrank)%P(iloc)%b  ! identificativo
         ibl = iblloc(jbl)  ! blocco locale

         i   = Xpt_req(jrank)%P(iloc)%i
         j   = Xpt_req(jrank)%P(iloc)%j
         k   = Xpt_req(jrank)%P(iloc)%k

         ! se la cella non è stata esaminata, ciclo
         if (.not.allocated(block(ibl)%gr(igr)%cell(i,j,k)%Chi%we)) cycle
         if (block(ibl)%gr(igr)%cell(i,j,k)%Chi%tycc == Search) cycle

         ! se il possibile donatore è più grande, non lo considera
         if (block(ibl)%gr(igr)%cell(i,j,k)%Vol < Xpt_req(jrank)%P(iloc)%dsize) cycle

         ! altrimenti si tiene il più grande di tutti quelli trovati
         if (block(ibl)%gr(igr)%cell(i,j,k)%Chi%we(1) < Xpt_req(jrank)%P(iloc)%dsize) then

            block(ibl)%gr(igr)%cell(i,j,k)%Chi%bd(1) = Xpt_req(jrank)%P(iloc)%bd
            block(ibl)%gr(igr)%cell(i,j,k)%Chi%id(1) = Xpt_req(jrank)%P(iloc)%id
            block(ibl)%gr(igr)%cell(i,j,k)%Chi%jd(1) = Xpt_req(jrank)%P(iloc)%jd
            block(ibl)%gr(igr)%cell(i,j,k)%Chi%kd(1) = Xpt_req(jrank)%P(iloc)%kd
            block(ibl)%gr(igr)%cell(i,j,k)%Chi%we(1) = Xpt_req(jrank)%P(iloc)%dsize
            block(ibl)%gr(igr)%cell(i,j,k)%Chi%tycc  = Xpt_req(jrank)%P(iloc)%tycc

         end if
         
      end do
   end do

   ! rimozione dei tentativi residui
   do ibl=1,nbl
      do k=-1,block(ibl)%gr(igr)%nk+2
      do j=-1,block(ibl)%gr(igr)%nj+2
      do i=-1,block(ibl)%gr(igr)%ni+2
         if (.not.allocated(block(ibl)%gr(igr)%cell(i,j,k)%Chi%we)) cycle
         if (block(ibl)%gr(igr)%cell(i,j,k)%Chi%tycc <= 0 ) cycle ! c.c.naturale

         if (block(ibl)%gr(igr)%cell(i,j,k)%Chi%tycc == Search) then
            call block(ibl)%gr(igr)%cell(i,j,k)%Chi%remove_X
         else
            block(ibl)%gr(igr)%cell(i,j,k)%Chi%we(1) = 1d0 ! un donatore, peso 1d0
         end if
       
      end do
      end do
      end do
   end do

   ! a questo punto tutte le celle sono marcate come regolari, contorni naturali,
   ! celle chimera o parete interna e c'è solo il donatore principale

   deallocate(npt)

contains

   subroutine check_patch(kbl,ibl,igr,i,j,k,id,jd,kd,tycc)
      use prec
      use concon
      use work_var
      implicit none 

      integer(kind=I4P),intent(in)  :: kbl,ibl,igr,i,j,k
      integer(kind=I4P),intent(out) :: id,jd,kd,tycc

      logical           :: internal
      integer(kind=I4P) :: ipatch
      integer(kind=I4P) :: jfcc
      integer(kind=I4P) :: ni,nj,nk

      ni = block(ibl)%gr(igr)%ni
      nj = block(ibl)%gr(igr)%nj
      nk = block(ibl)%gr(igr)%nk

      internal = i*j*k*(ni-i+1)*(nj-j+1)*(nk-k+1) > 0
      tycc = 0

      if (internal) return
      
      do ipatch=1,npatch
         if (kbl /= patch(ipatch)%bloc) cycle

         if (  i  < patch(ipatch)%imin) cycle
         if (  i  > patch(ipatch)%imax) cycle
         if (  j  < patch(ipatch)%jmin) cycle
         if (  j  > patch(ipatch)%jmax) cycle
         if (  k  < patch(ipatch)%kmin) cycle
         if (  k  > patch(ipatch)%kmax) cycle

         tycc = patch(ipatch)%tycc
         call block(ibl)%gr(igr)%cell(i,j,k)%Chi%add_X(1)

         id = i 
         jd = j
         kd = k
         if (i<1) then
            jfcc = 1
            id = -i+1
         else if (i>ni) then
            jfcc = 2
            id = 2*ni-i+1
         else if (j<1) then
            jfcc = 3
            jd = -j+1
         else if (j>nj) then
            jfcc = 4
            jd = 2*nj-j+1
         else if (k<1) then
            jfcc = 5
            kd = -k+1
         else if (k>nk) then
            jfcc = 6
            kd = 2*nk-k+1
         end if

         block(ibl)%gr(igr)%cell(i,j,k)%Chi%tycc  = tycc
         block(ibl)%gr(igr)%cell(i,j,k)%Chi%bd(1) = jfcc 
         block(ibl)%gr(igr)%cell(i,j,k)%Chi%id(1) = id   
         block(ibl)%gr(igr)%cell(i,j,k)%Chi%jd(1) = jd   
         block(ibl)%gr(igr)%cell(i,j,k)%Chi%kd(1) = kd   
         return

      end do

   end subroutine check_patch
   
   subroutine find_donor(xc,mp,bb,id,jd,kd,tycc,dsize)
      use prec
      use block_type
      use chimera_type
      use concon
      implicit none

      integer(kind=I4P), intent(in out) :: id,jd,kd
      integer(kind=I4P), intent(out)    :: tycc
      real(kind=R8P)   , intent(in)     :: xc(1:3)
      real(kind=R8P)   , intent(out)    :: dsize

      type(gen_block),intent(in) :: bb
      type(block_map),intent(in) :: mp

      logical                     :: trovato
      integer(kind=I4P)           :: imax,imin,jmax,jmin,kmax,kmin
      integer(kind=I4P)           :: ni,nj,nk
      integer(kind=I4P)           :: istep,jstep,kstep
      integer(kind=I4P)           :: dirwall
      integer(kind=I4P)           :: internal_check
      real(kind=R8P)              :: d,d_old
      real(kind=R8P)              :: d1,d2,d3
      real(kind=R8P),dimension(3) :: xp,dp
      real(kind=R8P),dimension(3) :: econtr

      imin = mp%imin(id,jd,kd)
      imax = mp%imax(id,jd,kd)
      jmin = mp%jmin(id,jd,kd)
      jmax = mp%jmax(id,jd,kd)
      kmin = mp%kmin(id,jd,kd)
      kmax = mp%kmax(id,jd,kd)

      ! dimensioni blocco
      ni = bb%ni
      nj = bb%nj
      nk = bb%nk

      ! stima iniziale
      id = max(1,(imax+imin)/2) ;  istep=max(1,(imax-imin)/2)
      jd = max(1,(jmax+jmin)/2) ;  jstep=max(1,(jmax-jmin)/2)
      kd = max(1,(kmax+kmin)/2) ;  kstep=max(1,(kmax-kmin)/2)

      i = id  ! stima iniziale
      j = jd
      k = kd

      xp = bb%cell(i,j,k)%Center
      dp = xc - xp
      d_old = dp .dot. dp ! distanza ^2

      !=========================================
      ! ricerca cartesiana
      !=========================================
      do 
         ! *************************************
         ! ricerca in i
         ! *************************************
         i = id + istep ; i = min(i,ni+1)
         xp = bb%cell(i,j,k)%Center
         dp = xc - xp
         d  = dp .dot. dp ! distanza ^2
         if (d < d_old) then
            d = d_old
            id = i
         else
            i = id - istep ; i = max(i,0)
            xp = bb%cell(i,j,k)%Center
            dp = xc - xp
            d  = dp .dot. dp ! distanza ^2
            if (d < d_old) then
               d = d_old
               id = i
            end if
         end if

         ! *************************************
         ! ricerca in j
         ! *************************************
         j = jd + jstep ; j = min(j,nj+1)
         xp = bb%cell(i,j,k)%Center
         dp = xc - xp
         d  = dp .dot. dp ! distanza ^2
         if (d < d_old) then
            d = d_old
            jd = j
         else
            j = jd - jstep ; j = max(j,0)
            xp = bb%cell(i,j,k)%Center
            dp = xc - xp
            d  = dp .dot. dp ! distanza ^2
            if (d < d_old) then
               d = d_old
               jd = j
            end if
         end if
         ! *************************************
         ! ricerca in k
         ! *************************************
         k = kd + kstep ; k = min(k,nk+1)
         xp = bb%cell(i,j,k)%Center
         dp = xc - xp
         d  = dp .dot. dp ! distanza ^2
         if (d < d_old) then
            d = d_old
            kd = k
         else
            k = kd - kstep ; k = max(k,0)
            xp = bb%cell(i,j,k)%Center
            dp = xc - xp
            d  = dp .dot. dp ! distanza ^2
            if (d < d_old) then
               d = d_old
               kd = k
            end if
         end if
  
         ! dimezzamento del passo
         istep = istep / 2
         jstep = jstep / 2
         kstep = kstep / 2
         if ( (istep+jstep+kstep) == 0) exit

         istep = max(1,istep)
         jstep = max(1,jstep)
         kstep = max(1,kstep)
      end do
      !=========================================
      ! ricerca cartesiana - fine
      !=========================================

      ! direzione parete locale
      dirwall = 0
      if (dsize < chimera_bl) then
         d1 = abs(bb%cell(id+1,jd,kd)%dwall - bb%cell(id-1,jd,kd)%dwall)
         d2 = abs(bb%cell(id,jd+1,kd)%dwall - bb%cell(id,jd-1,kd)%dwall)
         d3 = abs(bb%cell(id,jd,kd+1)%dwall - bb%cell(id,jd,kd-1)%dwall)
         if ((d1.gt.d2).and.(d1.gt.d3)) then
            if (minval(bb%cell(imin:imax,jd,kd)%dwall) < chimera_bl) &
            dirwall = 1
         else if ((d2.gt.d1).and.(d2.gt.d3)) then
            if (minval(bb%cell(id,jmin:jmax,kd)%dwall) < chimera_bl) &
            dirwall = 2
         else 
            if (minval(bb%cell(id,kd,kmin:kmax)%dwall) < chimera_bl) &
            dirwall = 3
         end if
      end if
         
      !=========================================
      ! ricerca covariante (non normalizzata)
      !=========================================
      
      ! direzione i
      ! vettore controvariante
      if (dirwall /= 1) then
         econtr = bb%cell(id,jd,kd)%Sn(1:3,1)
         
         xp = bb%cell(id,jd,kd)%Center
         dp = xc - xp
         ! distanza covariante 
         d_old  = abs(dp.dot.econtr)
         do     ! ricerca in i
            trovato = .true.
            i = id + 1 ; if (i > ni) exit  ! fuori blocco
            xp = bb%cell(i,j,k)%Center
            dp = xc - xp
            d  = abs(dp.dot.econtr)
            if (d < d_old) then
               trovato = (i==id)
               d = d_old
               id = i
            else
               i = id - 1; if (i < 1) exit  ! fuori blocco
               xp = bb%cell(i,j,k)%Center
               dp = xc - xp
               d  = abs(dp.dot.econtr)
               if (d < d_old) then
                  trovato = (i==id)
                  d = d_old
                  id = i
               end if
            end if
            if (trovato) exit
         end do
      else
         ! distanza dalla parete
         d_old  = abs(bb%cell(id,jd,kd)%dwall - dsize)
         do     ! ricerca in i
            trovato = .true.
            i  = min(id + 1,imax)
            d  = abs(bb%cell(i,j,k)%dwall - dsize)
            if (d < d_old) then
               trovato = (i==id)
               d = d_old
               id = i
            else
               i = max(id - 1,imin)
               d  = abs(bb%cell(i,j,k)%dwall - dsize)
               if (d < d_old) then
                  trovato = (i==id)
                  d = d_old
                  id = i
               end if
            end if
            if (trovato) exit
         end do
      end if
      
      ! direzione j
      ! vettore controvariante
      if (dirwall /= 2) then
         econtr = bb%cell(id,jd,kd)%Sn(1:3,2)
        
         xp = bb%cell(id,jd,kd)%Center
         dp = xc - xp
         ! distanza covariante 
         d_old  = abs(dp.dot.econtr)
         do     ! ricerca in j
            trovato = .true.
            j = jd + 1 ; if (j > nj) exit  ! fuori blocco
            xp = bb%cell(i,j,k)%Center
            dp = xc - xp
            d  = abs(dp.dot.econtr)
            if (d < d_old) then
               trovato = (j==jd)
               d = d_old
               jd = j
            else
               j = jd - 1; if (j < 1) exit  ! fuori blocco
               xp = bb%cell(i,j,k)%Center
               dp = xc - xp
               d  = abs(dp.dot.econtr)
               if (d < d_old) then
                  trovato = (j==jd)
                  d = d_old
                  jd = j
               end if
            end if
            if (trovato) exit
         end do
      else
         ! distanza dalla parete
         d_old  = abs(bb%cell(id,jd,kd)%dwall - dsize)
         do     ! ricerca in j
            trovato = .true.
            j  = min(jd + 1,jmax)
            d  = abs(bb%cell(i,j,k)%dwall - dsize)
            if (d < d_old) then
               trovato = (j==jd)
               d = d_old
               jd = j
            else
               j = max(jd - 1,jmin)
               d  = abs(bb%cell(i,j,k)%dwall - dsize)
               if (d < d_old) then
                  trovato = (j==jd)
                  d = d_old
                  jd = j
               end if
            end if
            if (trovato) exit
         end do
      end if

      ! direzione k
      ! vettore controvariante
      if (dirwall /= 3) then
         econtr = bb%cell(id,jd,kd)%Sn(1:3,3)
        
         xp = bb%cell(id,jd,kd)%Center
         dp = xc - xp
         ! distanza covariante 
         d_old  = abs(dp.dot.econtr)
         do     ! ricerca in k
            trovato = .true.
            k = kd + 1 ; if (k > nk) exit  ! fuori blocco
            xp = bb%cell(i,j,k)%Center
            dp = xc - xp
            d  = abs(dp.dot.econtr)
            if (d < d_old) then
               trovato = (k==kd)
               d = d_old
               kd = k
            else
               k = kd - 1; if (k < 1) exit  ! fuori blocco
               xp = bb%cell(i,j,k)%Center
               dp = xc - xp
               d  = abs(dp.dot.econtr)
               if (d < d_old) then
                  trovato = (k==kd)
                  d = d_old
                  kd = k
               end if
            end if
            if (trovato) exit
         end do
      else
         ! distanza dalla parete
         d_old  = abs(bb%cell(id,jd,kd)%dwall - dsize)
         do     ! ricerca in k
            trovato = .true.
            k  = min(kd + 1,kmax)
            d  = abs(bb%cell(i,j,k)%dwall - dsize)
            if (d < d_old) then
               trovato = (k==kd)
               d = d_old
               kd = k
            else
               k = max(kd - 1,kmin)
               d  = abs(bb%cell(i,j,k)%dwall - dsize)
               if (d < d_old) then
                  trovato = (k==kd)
                  d = d_old
                  kd = k
               end if
            end if
            if (trovato) exit
         end do
      end if

      !=========================================
      ! ricerca controvariante - fine
      !=========================================

      ! verifica punto interno
      internal_check = id * (ni-id+1) &
                     * jd * (nj-jd+1) &
                     * kd * (nk-kd+1) 
      if (internal_check > 0 ) then
         ! donatore chimera
         dsize = bb%cell(id,jd,kd)%Vol
         tycc  = CellChim
      else
         ! possibile parete interna
         call check_wall()  ! verifica la possibilità 
                            ! di parete interna
      end if

   end subroutine find_donor

end subroutine search_principal_donor
subroutine check_wall()
    write(*,*) ' check_wall da implementare '
    stop
end subroutine check_wall

