subroutine res_mb(igr,init_res,add_der_temp)
   use prec
   use work_var
   use motion_par  , only: glomob
   use physical_par, only: rey
   use solution_par, only: level_set,finest,dtempo

   implicit none

   logical          ,intent(in) :: init_res,add_der_temp
   integer(kind=I4P),intent(in) :: igr

   integer(kind=I4P)            :: ibl

   if (init_res) then
      if (finest) then
         do ibl=1,nbl ! set residual to zero
            call block(ibl)%gr(igr)%null_res()
         end do
      else
         do ibl=1,nbl ! set residual equal to source term
            call block(ibl)%gr(igr)%ini_sorg()
         end do
      end if
   end if

   call cont_mb(igr)   ! condizioni al contorno + forzanti chimera

   ! Calcolo velocità di griglia
   if (glomob) then
      do ibl=1,nbl
         call block(ibl)%gr(igr)%grid_vel_face(block(ibl)%gruppo)
      end do
   end if

   ! calcolo termini euleriani
   do ibl=1,nbl
      call block(ibl)%gr(igr)%euler_flux_balance()
   end do

   ! calcolo termini euleriani
   if (rey < infinity) then
      do ibl=1,nbl
         call block(ibl)%gr(igr)%visc_flux_balance()
      end do
   end if

   ! derivate temporali
   if (finest.and.add_der_temp) then
      do ibl=1,nbl
         call block(ibl)%gr(igr)%add_dert(dtempo)
      end do
   end if
 
   ! estrapolazione level-set
   if (level_set) then
      do ibl=1,nbl
         call block(ibl)%gr(igr)%level_set_ext()
      end do
   end if

   ! memorizzazioni sorgenti
   if (.not.init_res) then
      do ibl=1,nbl
         call block(ibl)%gr(igr)%mem_sorg()
      end do
   end if

end subroutine res_mb
