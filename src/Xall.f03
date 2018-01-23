!=======================================================================
!     main program
!=======================================================================
program Xall
   use prec
   use solution_par
   use motion_par
   use compute_par
   use work_var
   use mpi_var
   use file_set
   use mpi

   implicit none

   logical           :: marcia
   integer(kind=I4P) :: igr,ibl,ifgr,n,ierr,iter,nmaxsg
!$ integer(kind=I4P) :: omp_get_num_threads, omp_get_thread_num
   real(kind=R8P)    :: tempoi

!=======================================================================
!  Inizializzazione MPI
!=======================================================================
   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc,ierr)
   write(*,'(a,i3,a,i3,a)') 'Myrank: ', myrank, ' of ', nproc, ' is alive'

!$omp parallel   &
!$omp default(none) &
!$omp shared(myrank,n)
!$    n = omp_get_num_threads()
!$    write(*,'(a,i2,a,i2,a,i2,a)')  &
!$         'Myrank: ', myrank, ' Thread: ', omp_get_thread_num(), &
!$         ' of ', omp_get_num_threads(), ' is alive'
!$omp end parallel 

   if (myrank==0) then
     write(*,'(a,i3,a,i3)') 'Myrank:', myrank, &
                           ' numero processi (MPI):', nproc
!$   write(*,'(a,i3,a,i3)') 'Myrank:', myrank, &
!$                         ' numero  threads (OMP):', n
     write(*,*)
   end if

   cpu_io    = 0.0d0
   cpu_rcc   = 0.0d0
   cpu_com   = 0.0d0
   cpu_residui = 0.0d0
   cpu_start = MPI_Wtime()    ! Tempo di riferimento

!
!  input parametri
!
   if (myrank==0) write(*,*) ' Lettura parametri del problema '
   call input_par()
   call input_proc()
!
!  lettura e creazione blocchi 
!
   if (myrank==0) write(*,*) ' Lettura blocchi '
   call input_block()
!
!  lettura topologia e condizioni al contorno
!
   if (myrank==0) write(*,*) ' Lettura topologia e condizioni al contorno'
   call input_cc()
!
!  lettura condizioni iniziali
!
   if (myrank==0) write(*,*) ' Lettura condizioni iniziali'
   do igr=1,sgr
      call input_init(igr)
   end do
!
!  lettura gruppi e corpi per il moto griglia
!
   if (myrank==0) write(*,*) ' Lettura gruppi e corpi'
   call input_motion()

   if (myrank==0) write(*,*) ' Fine lettura input '
!=======================================================================
      cpu_io = cpu_io + (MPI_Wtime() - cpu_start)
!==============================================================================
! Matrici di comunicazione e allocazione della memoria per
! i punti "non strutturati"
   if (nproc/=1) then   ! Non occorre per seriale
      do igr=1,ngr        ! L'ordine del ciclo non va cambiato!!!!
         call build_comm_map(igr)
      end do
   end if
!=======================================================================

!-----------------------------------------------------------------------
!  soluzione  Full MultiGrid  o ciclo V semplice
!-----------------------------------------------------------------------
   n = 0
   iter  = 0
   work  = 0d0

   workgr(1) = 1d0
   do igr=2,ngr
      workgr(igr) = 0.125d0*workgr(igr-1)
   end do

   do ifgr=sgr,1,-1

!.......................................................................
!  calcolo termini metrici
!.......................................................................
      do igr=ifgr,ngr
         if (igr==ifgr .or. ifgr==sgr) then
            call metr_mb(igr)
            if (myrank==0) &
            write(*,'(a,i3,a,i2.2)') ' calcolo metrica livello...',igr
         end if
      end do
      if (myrank==0) &
      write(*,'(a)') ' calcolo metrica - fatto'

      igr = ifgr
!.......................................................................
!
!  stima iniziale della soluzione a partire dalla griglia rada
!
      if (igr<sgr) then
         do ibl=1,nbl
            call block(ibl)%restrict(igr)
         end do
      end if

      ! Evoluzione della soluzione nel tempo
      tempoi = temposg*(sgr-igr+1)
      if (igr==1) tempoi = tempo_max
      marcia = .true.
      itime = 0

      do while (marcia)

         cpu_tmp = MPI_Wtime()
         itime = itime + 1 ! indice passo (pseudo)temporeale

         if (unsteady) then
            finest = .true.  ! calcolo passo temporale DTEMPO
            call step_mb(ifgr,.true.)

            tempo = tempo + dtempo
            n = n + 1

            call updinf()   ! aggiornamento velocita' all'infinito
            call updmcc()   ! aggiornamento (eventuale) matcc

            ! memorizzazione soluzione al passo n
            call mem_old_sol(ifgr)

            ! PRE SOLUZIONE
            call pre_step(ifgr)

            if (implicito) then
               ! memorizzazione soluzione al passo n
               ! sulle griglie rade
               do igr=ifgr+1,ngr
                  call mem_old_sol(igr)
               end do

               ! integrazione implicita con MG
               call MGsolve(iter,ifgr,nmaxsg)  
            else
               finest = .true.
               ! integrazione Runge-Kutta
               call RungeKutta(ifgr)
            end if

            ! POST SOLUZIONE
            call post_step(ifgr)

            ! .........................................................
            ! calcolo e stampa forze 
            ! .........................................................
            call for_mb(ifgr)
            if (myrank==0) then  ! scrittura forze
               call wrievol('u',ifgr,ngruppi)
               call wrievol('h',ifgr,ncorpi)
               frst2 = .false.
            end if

         else  ! stazionario

             
            if (implicito) then
               call MGsolve(iter,ifgr,nmaxsg)  
            else
               do iter=1,nmaxsg
                  ! calcolo residui
                     call RungeKutta(ifgr)
               
                  call l2_norm_res(iter,ifgr)
                  if (resmx < epsv) exit
                  if (    mod(iter,ifsto)==0 ) &  ! mem. soluzione
                  call memsol_freq(ifgr)
               
               end do
            end if

         end if
         cpu_residui = cpu_residui + (MPI_Wtime()-cpu_tmp)

         marcia = unsteady.and.(tempo<=tempoi)

!        if (runavr) call Running(ifgr)  ! Running Averages

! scrittura ultima soluzione calcolata

         if (    (unsteady.and.mod(itime,ifsto)==0) &   ! mem. soluzione
             .or.(.not.marcia)) &  ! Stazionario - ultima soluzione
         call memsol_freq(ifgr)

         ! estrazione soluzione
         call extsol_mb(ifgr)
         ! Stampa tempi di calcolo
         if(unsteady) then
            cpu_tot = MPI_Wtime() - cpu_start
            call wricputime(frst4,filevo, &
                            itime,tempo,cpu_tot, &
                            cpu_io,cpu_rcc,cpu_com,cpu_residui)
         end if

      end do

!-----------------------------------------------------------------------
!  fine soluzione  Full MultiGrid
!-----------------------------------------------------------------------
   end do

   if (myrank==0) write(*,'(a,i3,a,i3,a)') ' Fine esecuzione'

   cpu_tot = MPI_Wtime() - cpu_start
   call wricputime(frst4,filevo, &
                   itime,tempo,cpu_tot, &
                   cpu_io,cpu_rcc,cpu_com,cpu_residui)
!=======================================================================

   call MPI_FINALIZE(ierr)     ! FINE
   stop

contains

   subroutine memsol_freq(ifgr)
      use prec
      implicit none
      integer(kind=I4P)    :: ifgr,igr

      if (nproc/=1) then   ! Non occorre per seriale
         do igr = ifgr,ngr
           call comm_invmodrcc(igr) ! Globalizzazione di rcc
         end do
         call memsol_mb(ifgr)
         do igr = ifgr,ngr
            call comm_locrcc(igr)    ! Localizzazione di rcc
         end do
      else
         call memsol_mb(ifgr)
      end if
   end subroutine memsol_freq

end program Xall
!=======================================================================
subroutine wricputime(frst4,filevo,i,t,tot,io,rcc,com,res)
!.......................................................................
! Stampa cpu time
!.......................................................................
   use prec
   use mpi_var, only:myrank
   use mpi
   implicit none

   logical                     :: frst4
   character(len=*),intent(in) :: filevo

   integer(kind=I4P)      :: l
   integer(kind=I4P)      :: i
   integer(kind=I4P)      :: fcpu

   real(kind=R8P)         :: t,tot,io,rcc,com,res

   l = len_trim(filevo)

! Il processo 0 (quello che scrive) riceve da tutti gli altri le
! informazioni sui tempi di calcolo

   if (myrank==0) then    ! Ricezione da tutti

      call getunit(fcpu)
      if (frst4) then
         open(fcpu,file=filevo(1:l)//'cpu.dat',form='formatted', &
                  status='unknown')
         rewind(fcpu)
         frst4 = .false.
      else
         open(fcpu,file=filevo(1:l)//'cpu.dat',form='formatted',  &
                  status='old',position='append')
      end if
      write(fcpu,'(1x,i8,6(e20.12))') i,t,tot,io,rcc,com,res
      close(fcpu)
   end if

   return
end subroutine wricputime
!.......................................................................
