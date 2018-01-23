!***********************************************************************  
! calcolo stati destro e sinistro per il problema di Riemann
!***********************************************************************  
!
!
!***********************************************************************  
!***********************************************************************  
! N.B.     UL(i) : stato a SINISTRA della faccia i+1/2
!          UR(i) : stato a DESTRA   della faccia i-1/2
!***********************************************************************  
!***********************************************************************  
!
!
!.......................................................................
! GODUNOV 1 ORDINE
!.......................................................................
subroutine lrgo1(n,q,ql,qr)
!.......................................................................
   use prec
   use var_type
   implicit none

   integer(kind=I4P),intent(in)  :: n
   type(gen_var),intent(in)     :: q(-1:*)
   type(gen_var),intent(out)    :: ql(0:*),qr(0:*)

   integer(kind=I4P)  :: i

   do i=0,n+1
      ql(i) = q(i)
      qr(i) = q(i)
   end do

   return
end subroutine lrgo1
!.......................................................................
! TVD SECONDO ORDINE
!.......................................................................
subroutine lrtvd(n,q,ql,qr)
   use prec
   use var_type
   implicit none

   integer(kind=I4P),intent(in) :: n
   type(gen_var),intent(in)     :: q(-1:*)
   type(gen_var),intent(out)    :: ql(0:*),qr(0:*)

   integer(kind=I4P)                        :: i,ierr
   type(gen_var)                            :: mmdiff
   type(gen_var),allocatable,dimension(:)   :: diff
   !===============================================
   interface
      elemental function minmod(x,y) result(mm)
         use prec
         use var_type
         use block_type, only: nvar
         implicit none

         type(gen_var),intent(in)  :: x,y
         type(gen_var)             :: mm
      end function minmod
   end interface
   !===============================================

   allocate(diff(0:n+2),STAT=ierr) ; call errore('ALLocate diff',ierr)

   do i=0,n+2
      diff(i) = q(i)-q(i-1)
   end do

   do i=0,n+1
      mmdiff = 0.5d0*minmod(diff(i),diff(i+1))
      ql(i) = q(i) + mmdiff
      qr(i) = q(i) - mmdiff
   end do
   deallocate(diff,STAT=ierr)  ; call errore('DEALLocate diff',ierr)

   return
end subroutine lrtvd
!.......................................................................
! WENO  TERZO ORDINE
!.......................................................................
subroutine lr3we(n,q,ql,qr)
!.......................................................................
   use prec
   use var_type
   implicit none

   integer(kind=I4P),intent(in) :: n
   type(gen_var),intent(in)     :: q(-1:*)
   type(gen_var),intent(out)    :: ql(0:*),qr(0:*)

   integer(kind=I4P)                       :: i,ierr
   type(gen_var)                           :: umk4,upk4
   type(gen_var)                           :: den,sumIS
   type(gen_var),allocatable,dimension(:)  :: diff,IS

   allocate(diff(0:n+2),STAT=ierr) ; call errore('ALLocate diff',ierr)
   allocate(  IS(0:n+2),STAT=ierr) ; call errore('ALLocate IS  ',ierr)

   do i=0,n+2
      diff(i) = q(i)-q(i-1)
      IS(i)   = diff(i)**4
   end do
   
   do i=0,n+1
      sumIS = IS(i+1)+IS(i)+zero 
   
      den   = 1d0/(sumIS+IS(i))
      umk4  = IS(i+1)*den
      upk4  = 2d0*IS(i)*den
      ql(i) = q(i) + 0.5d0*(upk4*diff(i+1)+umk4*diff(i  ))
      
      den   = 1d0/(sumIS+IS(i+1))
      umk4  = IS(i)*den
      upk4  = 2d0*IS(i+1)*den
      qr(i) = q(i) - 0.5d0*(upk4*diff(i  )+umk4*diff(i+1))
   end do

   deallocate(diff,STAT=ierr)  ; call errore('DEALLocate diff',ierr)
   deallocate(  IS,STAT=ierr)  ; call errore('DEALLocate   IS',ierr)

   return
end subroutine lr3we
!.......................................................................
!  minmod classico
!.......................................................................
!
elemental function minmod(x,y) result(mm)
   use prec
   use var_type
   use block_type, only: nvar
   implicit none

   type(gen_var),intent(in)  :: x,y
   type(gen_var)             :: mm

   integer(kind=I4P)         :: ivar
   do ivar=1,nvar
      mm%v(ivar) = sign(min(abs(x%v(ivar)),abs(y%v(ivar))),x%v(ivar))   ! classico
      if ((x%v(ivar)*y%v(ivar))<=0.d0) then
         mm%v(ivar) = 0.d0
      end if
   end do

!  minmod = 0.0d0
!  if ((x*y) <= 0.0d0)  return
!  minmod=2.0d0*x*y/(x+y)

   return
end function minmod
