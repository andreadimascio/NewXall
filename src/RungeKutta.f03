subroutine RungeKutta(igr)
   use prec
   use work_var
   use solution_par, only: finest,ark,brk,crk,nrk,dtempo

   implicit none

   integer(kind=I4P), intent(in) :: igr

   logical,parameter :: add_der_temp = .false.
   logical,parameter :: init_res     = .true.
   integer(kind=I4P) :: irk,ibl

   ! Runge Kutta TVD (assegnato in input)
   !
   !  RK TDV 1^o ordine (Eulero) : 
   !              nrk = 1
   !              ark(1) = 1d0     ; brk(1) = 0d0     ; crk(1)*dt  = dt
   !
   !  RK TDV 2^o ordine : 
   !              nrk = 2
   !              ark(1) = 1d0     ; brk(1) = 0d0     ; crk(1)*dt  = dt
   !              ark(2) = 1d0/2d0 ; brk(2) = 1d0/2d0 ; crk(2)*dt  = dt/2d0
   !
   !  RK TDV 3^o ordine : 
   !              nrk = 3
   !              ark(1) = 1d0     ; brk(1) = 0d0      ; crk(1)*dt  = dt
   !              ark(2) = 3d0/4d0 ; brk(2) = 1d0/4d0  ; crk(2)*dt  = dt/4d0
   !              ark(2) = 1d0/3d0 ; brk(2) = 2d0/3d0  ; crk(3)*dt  = 2d0*dt/3d0
   !
   !  RK NON TDV 2^o ordine (per i centrati) : 
   !      nrk = 3
   !      ark(1) = 1d0 ; brk(1) = 0d0  ; crk(1) = 0.5d0
   !      ark(2) = 1d0 ; brk(2) = 0d0  ; crk(2) = 0.5d0
   !      ark(3) = 1d0 ; brk(3) = 0d0  ; crk(3) = 1d0

   finest = .true.
   do irk=1,nrk

      call res_mb(igr,init_res,add_der_temp)

      do ibl=1,nbl
         call block(ibl)%gr(igr)%res_dt(dtempo)
         ! N. B. in uscita = dru/dt*vol : servono per le forze
         ! residui * dt /(1d0+dt*kfor/delta)
      end do

      do ibl=1,nbl ! update variables
         call block(ibl)%gr(igr)%update_rk(ark(irk),brk(irk),crk(irk))
      end do
   end do

end subroutine RungeKutta
