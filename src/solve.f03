!=======================================================================
!     calcolo variabili di stato 
!=======================================================================
subroutine solve(nmaxl,igr,ifin,ncf)
!.......................................................................
      use prec
      use motion_par
      use solution_par
      implicit none

      integer(kind=I4P),intent(in)     :: igr,ifin,nmaxl
      integer(kind=I4P),intent(in out) :: ncf

      logical,parameter :: init_res = .true.,add_der_temp = .true.
      integer(kind=I4P) :: n

      if (nmaxl.le.0) return

      n     = 0
!
!  inizio iterazione
!
      do while (n.lt.nmaxl)
         n = n+1

         if (glodina.and.finest) then
            call din_mb(.false.)    ! calcolo dinamica del corpo
            call din_vfi()   ! aggiornamento velocita' figli
         end if

         ! calcolo residui
         call res_mb(igr,init_res,add_der_temp)

         !  calcolo passo pseudo-temporale
         call step_mb(igr,.false.)

         ! norma l2 dei residui  
         if (ncf.eq.0.and.finest) call l2res_mb(igr,ncf)

         ! aggiornamento variabili
         call agg_mb(igr,ifin)

!        CHECK  !!!!!!
!        ! limite var conservate
!        call limvar_mb(igr)
!
!        ! calcolo variabili ausiliarie
!        call varaux_mb(igr)

      end do

      return
end subroutine solve

!=======================================================================
subroutine l2res_mb(igr,ncf)
!.......................................................................
      use prec
      use var_type
      use work_var
      use solution_par, only: resmx,resq
      use mpi_var, only: nproc
      use mpi
      implicit none

      integer(kind=I4P),intent(in    ) :: igr
      integer(kind=I4P),intent(in out) :: ncf
      integer(kind=I4P) :: ibl,ierr
      real(kind=R8P)    :: sumloc(nvar),sumglo(nvar) ! automatiche per MPI
      real(kind=R8P)    :: volloc,volglo


      sumloc = 0d0
      volloc = 0.0d0
      do ibl=1,nbl
         call residui(ibl,igr,sumloc,volloc)
      end do

      if (nproc.ne.1) then
         call MPI_ALLREDUCE(volloc,volglo,1,MPI_REAL8, &
                            MPI_SUM,MPI_COMM_WORLD,ierr)
         volglo =  max(volglo,1d-20)

         call MPI_ALLREDUCE(sumloc,sumglo,nvar,MPI_REAL8, &
                            MPI_SUM,MPI_COMM_WORLD,ierr)
         resq   = sqrt(sumglo/volglo)
      else
         resq   = sqrt(sumloc/volloc)
      end if

      ncf = 1
      resmx = maxval(resq)

      return
contains
   !=======================================================================
   subroutine residui(ibl,igr,locresq,voltot)
   !.......................................................................
         use prec 
         use var_type
         use work_var
         implicit none

         integer(kind=I4P),intent(in)   :: ibl,igr
         real(kind=R8P),intent(in out)  :: locresq(1:nvar)
         real(kind=R8P),intent(in out)  :: voltot
   
         integer(kind=I4P) :: i,j,k
         real(kind=R8P)  :: den,vol
   
   !==============================================================================
   !     cpu_tmp = MPI_Wtime()       ! DEBUG
   !==============================================================================
   !$omp parallel do collapse(3) schedule(static) &
   !$omp  default(none) &
   !$omp  shared(ibl,igr,nvar,block) &
   !$omp  private(i,j,k,den,vol) &
   !$omp  reduction(+: locresq,voltot)
   !==============================================================================
         do k=1,block(ibl)%gr(igr)%nk
            do j=1,block(ibl)%gr(igr)%nj
               do i=1,block(ibl)%gr(igr)%ni
                  vol = block(ibl)%gr(igr)%cell(i,j,k)%Vol
                  if (block(ibl)%gr(igr)%cell(i,j,k)%Chi%tycc == 0) then
                     den = vol
                  else
                     den = max(1d0,block(ibl)%gr(igr)%cell(i,j,k)%fchi**2) &
                         / vol
                  end if
                  den = 1d0/den
                  locresq(1:nvar) = locresq(1:nvar)  &
                                  + den*block(ibl)%gr(igr)%dq(i,j,k)%v(1:nvar)**2
                  voltot  = voltot  + vol
               end do
            end do
         end do
   !$omp end parallel do
   !==============================================================================
   !     cpu_residui = cpu_residui + (MPI_Wtime()-cpu_tmp)      ! DEBUG
   !==============================================================================
         return
   end subroutine residui
end subroutine l2res_mb
!=======================================================================
subroutine agg_mb(igr,ifin)
!.......................................................................
      use prec
      use work_var
      use solution_par, only: stab
      implicit none

      integer(kind=I4P),intent(in) :: igr,ifin
      integer(kind=I4P)            :: ibl

      do ibl=1,nbl     ! AGGIORNAMENTO VARIABILI
         call aggvar(ibl,igr,stab(ifin))
      end do

      return
contains
   !.......................................................................
   subroutine aggvar(ibl,igr,stab)
   !.......................................................................
         use prec
         use work_var
         use motion_par
         implicit none
   
         integer(kind=I4P),intent(in) :: ibl,igr
         real(kind=R8P)   ,intent(in) :: stab

         integer(kind=I4P) :: i,j,k
         real(kind=R8P)    :: dt
         real(kind=R8P)    :: omi(3,3),dd(3),rr(3)
   !==============================================================================
   !     cpu_tmp = MPI_Wtime()       ! DEBUG
   !==============================================================================
   !$omp parallel   & 
   !$omp  default(none) & 
   !$omp  shared(sysrot,OMEGA,OMEGA_2,ibl,igr,block) & 
   !$omp  private(i,j,k,dt,omi,dd,rr)
   !==============================================================================
   !$omp do collapse(3) schedule(static)
         do k=1,block(ibl)%gr(igr)%nk
            do j=1,block(ibl)%gr(igr)%nj
               do i=1,block(ibl)%gr(igr)%ni
                  dt = block(ibl)%gr(igr)%cell(i,j,k)%dtmat
                  block(ibl)%gr(igr)%dq(i,j,k) = dt*block(ibl)%gr(igr)%dq(i,j,k)
               end do
            end do
         end do
   !$omp enddo !nowait
         if (sysrot) then
            !$omp do collapse(3) schedule(static)
            do k=1,block(ibl)%gr(igr)%nk
               do j=1,block(ibl)%gr(igr)%nj
                  do i=1,block(ibl)%gr(igr)%ni
                     dt = block(ibl)%gr(igr)%cell(i,j,k)%dtmat
                     call matrotin(dt,block(ibl)%gr(igr)%cell(i,j,k)%Vol,OMEGA_2,OMEGA,omi)
                     dd(1) = block(ibl)%gr(igr)%dq(i,j,k)%v(q_u)
                     dd(2) = block(ibl)%gr(igr)%dq(i,j,k)%v(q_v)
                     dd(3) = block(ibl)%gr(igr)%dq(i,j,k)%v(q_w)
                     rr = matmul(omi,dd)
                     block(ibl)%gr(igr)%dq(i,j,k)%v(q_u) = rr(1)
                     block(ibl)%gr(igr)%dq(i,j,k)%v(q_v) = rr(2)
                     block(ibl)%gr(igr)%dq(i,j,k)%v(q_w) = rr(3)
                  end do
               end do
            end do
            !$omp enddo !nowait
         end if
   !$omp end parallel
   !==============================================================================
   !     cpu_residui = cpu_residui + (MPI_Wtime()-cpu_tmp)      ! DEBUG
   !==============================================================================
         if (stab.ge.0.499d0) call implsm(ibl,igr)
   !==============================================================================
   !     cpu_tmp = MPI_Wtime()       ! DEBUG
   !==============================================================================
   !$omp parallel &
   !$omp  default(none)  &
   !$omp  shared(ibl,igr,block) &
   !$omp  private(i,j,k)
   !==============================================================================
   
   !$omp do collapse(3) schedule(static)
         do k=1,block(ibl)%gr(igr)%nk
            do j=1,block(ibl)%gr(igr)%nj
               do i=1,block(ibl)%gr(igr)%ni
                  block(ibl)%gr(igr)%q(i,j,k) = block(ibl)%gr(igr)%q(i,j,k) &
                                              - block(ibl)%gr(igr)%dq(i,j,k)
               end do
            end do
         end do
   !$omp enddo !nowait
   !$omp end parallel
   !==============================================================================
   !     cpu_residui = cpu_residui + (MPI_Wtime()-cpu_tmp)      ! DEBUG
   !==============================================================================
         return
   end subroutine aggvar
   subroutine matrotin(dt,vol,OMEGA_2,OMEGA,omi)
      use prec
      implicit none
   
      real(kind=R8P),intent(in)                 :: dt,vol,OMEGA_2
      real(kind=R8P),dimension(3),intent(in)    :: OMEGA
      real(kind=R8P),dimension(3,3),intent(out) :: omi

      real(kind=R8P)                :: den,locdt
      real(kind=R8P),dimension(3)   :: omdt
   
      locdt = dt*vol
      omdt = OMEGA*locdt
   
      locdt = OMEGA_2*locdt**2
      den   = 1d0/(1d0+locdt)
   
      omi(1,1) = 1d0+omdt(1)**2
      omi(1,2) =  omdt(3)+omdt(1)*omdt(2)
      omi(1,3) = -omdt(2)+omdt(1)*omdt(3)
   
      omi(2,1) = -omdt(3)+omdt(2)*omdt(1)
      omi(2,2) = 1d0+omdt(2)**2
      omi(2,3) =  omdt(1)+omdt(2)*omdt(3)
   
      omi(3,1) =  omdt(2)+omdt(3)*omdt(1)
      omi(3,2) = -omdt(1)+omdt(3)*omdt(2)
      omi(3,3) = 1d0+omdt(3)**2
   
      omi = omi*den
   end subroutine matrotin
end subroutine agg_mb
