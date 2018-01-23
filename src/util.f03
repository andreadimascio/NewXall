subroutine errore(Stringa,ierr)
   use prec
   use mpi
   implicit none

   character(len=*),intent(in)     :: Stringa
   integer(kind=I4P),intent(in)    :: ierr
   integer(kind=I4P)               :: ierrl
   
   if (ierr==0) return

   write(*,*) trim(Stringa)
   call MPI_FINALIZE(ierrl)     ! FINE
   stop

end subroutine errore
