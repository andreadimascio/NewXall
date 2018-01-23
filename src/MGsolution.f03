subroutine MGsolution(iter,ifgr,nmax)
   use prec
   use compute_par
   use solution_par,     only: resmx,epsv

   implicit none

   integer(kind=I4P)               ,intent(in out) :: iter
   integer(kind=I4P)               ,intent(in)     :: ifgr,nmax

   integer(kind=I4P) :: icycle

   do icycle=1,nmax
      work_old = work
      call Vcycle(ifgr)

      iter = iter + 1
      call monitor(ifgr,iter,work_old)

      if (resmx<epsv) exit
   end do

   return

contains
   
   subroutine Vcycle(ifgr)
      use prec
      use work_var
      use solution_par
   
      implicit none
   
      integer(kind=I4P),intent(in) :: ifgr
   
      logical,parameter :: add_der_temp = .true.
      logical           :: init_res
      integer(kind=I4P) :: igr,ibl,ncf
      integer(kind=I4P) :: ifin
   
      ncf = 0
      !.......................................................
      !  DESCENT
      !.......................................................
      do igr=ifgr,ngr
         finest = (igr==ifgr)
         ifin   = igr-ifgr+1

         solvin = .false.
         if (.not.finest) then
            ! restriction on coarse grid
            do ibl=1,nbl
               call block(ibl)%restrict(igr)
            end do
         end if

         if (igr/=ifgr) then ! source term calculation
   
            ! set souce terms to zero
            do ibl=1,nbl
              call block(ibl)%nullsorg(igr)
            end do
   
            ! source term on fine grid
            finest   = (igr-1)==ifgr
            init_res = .true.
            call res_mb(igr-1,init_res,add_der_temp)
            work = work + workgr(igr-1)
   
            ! collect residuals
            do ibl=1,nbl
               call block(ibl)%collres(igr)
            end do
   
            ! residui on coarse grid and sum into source
            finest   = .false.
            init_res = .true.
            call res_mb(igr,init_res,add_der_temp)
            work = work + workgr(igr)
   
         end if
   
         solvin = .true.
         call solve(npre(ifin),igr,ifin,ncf)
   
         ! residual norm computed only at cycle start
         if (finest.and.(resmx<epsv)) return
   
         work = work + workgr(igr)*npre(ifin)
   
      end do
   
   !..........................................................
   !  ASCENT
   !..........................................................
   
      do igr=ngr-1,ifgr,-1
   
         finest = igr==ifgr
         ifin   = igr-ifgr+1

         !  coarse --> fine correction
         ! 1) variations on coarse grid
         do ibl=1,nbl
            call block(ibl)%calcvar(igr)
         end do
         ! 2) correction on fine grid
         do ibl=1,nbl
            call block(ibl)%prolvar(igr,level_set)
         end do
   
         solvin = .true.
         call solve(npost(ifin),igr,ifin,ncf)
   
         work = work + workgr(igr)*npost(ifin)*nrk
      end do

      !..........................................................
      !  END OF V-CYCLE
      !..........................................................
   
   end subroutine Vcycle
   
   
   subroutine monitor(igr,iter,work)
      use prec
      use solution_par  , only: resq
      use work_var, only: nvar
      use file_set

      implicit none
   
      integer(kind=I4P),intent(in) :: igr,iter
      integer(kind=I4P)            :: l,ivar
      real(kind=R8P)               :: work
   
      l = len_trim(filevo)
   
      if (frst1) then
         call getunit(fres)
         open(fres,file=filevo(1:l)//'res.dat',form='formatted', &
              status='unknown',position='rewind')
         frst1 = .false.
      else
         call getunit(fres)
         open(fres,file=filevo(1:l)//'res.dat',form='formatted', &
              status='old',position='append')
      end if
   
      write(fres,'(i2,i9,10e16.8)') igr,iter,work,(resq(ivar),ivar=1,nvar)
      close(fres)
   
   end subroutine monitor
   
   
end subroutine MGsolution
