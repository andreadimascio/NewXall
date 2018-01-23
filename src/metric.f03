subroutine metr_mb(igr)
   use prec
   use work_var
   implicit none

   integer(kind=I4P),intent(in) :: igr
   integer(kind=I4P)            :: ibl

   do ibl=1,nbl
      call block(ibl)%gr(igr)%Centroid
      call block(ibl)%gr(igr)%faces
      call block(ibl)%gr(igr)%volume
      call block(ibl)%gr(igr)%tens
   end do
end subroutine metr_mb
