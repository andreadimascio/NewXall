!=======================================================================
! routine di transito per scelta schema
!=======================================================================
subroutine reseul_mb(igr)
   use prec
   use var_type
   use work_var
   use proc_pointers
   use solution_par, only:finest

   implicit none
   
   integer(kind=I4P),intent(in)   ::  igr
   integer(kind=I4P)              ::  ibl

   if (finest) then
      lrstat => lrstat_finest
   else
      lrstat => lrgo1
   end if

   do ibl=1,nbl
      call block(ibl)%gr(igr)%euler_flux_balance()
   end do

   if (finest) call chifor_mb(igr) ! forzanti chimera
   return
end subroutine reseul_mb
!-----------------------------------------------------------------------
subroutine chifor_mb(igr)
!.......................................................................
   use prec
   use cell_type
   use block_type
   use work_var
   use concon
   implicit none


   integer(kind=I4P),intent(in)  :: igr

   integer(kind=I4P)  :: ibl,jbl,i1,j1,k1
   integer(kind=I4P)  :: i,j,k,ipcc,npcc,tycc
   integer(kind=I4P)  :: jgruppo
   type(gen_var),allocatable,dimension(:)     ::  qdon

   do ibl=1,nbl
!==============================================================================
!$omp parallel &
!$omp  default(none) &
!$omp  shared(nvar,nbl,igr,ibl,block) &
!$omp  private(i,j,k,ipcc,tycc,npcc,jbl,jgruppo,i1,j1,k1,qdon)
!==============================================================================
!$omp do
      do k=1,block(ibl)%gr(igr)%nk
      do j=1,block(ibl)%gr(igr)%nj
      do i=1,block(ibl)%gr(igr)%ni

         tycc = block(ibl)%gr(igr)%cell(i,j,k)%Chi%tycc
         if (tycc == 0) cycle ! cella regolare

         ! chimera o parete interna

         jgruppo = 0
         if (tycc==CellChim) then
            ! numero donatori nel punto
            npcc = block(ibl)%gr(igr)%cell(i,j,k)%Chi%npcc

            if (allocated(qdon)) deallocate(qdon)
            allocate(qdon(npcc))

            do ipcc=1,npcc
               jbl = block(ibl)%gr(igr)%cell(i,j,k)%Chi%bd(ipcc)
               i1  = block(ibl)%gr(igr)%cell(i,j,k)%Chi%id(ipcc)
               j1  = block(ibl)%gr(igr)%cell(i,j,k)%Chi%jd(ipcc)
               k1  = block(ibl)%gr(igr)%cell(i,j,k)%Chi%kd(ipcc)

               call ass_var( &
               block(jbl)%gr(igr)%ni,block(jbl)%gr(igr)%nj,block(jbl)%gr(igr)%nk, &
               i1,j1,k1, block(jbl)%gr(igr)%q,qdon(ipcc))
            end do
         else  ! parete interna
            npcc = 1

            if (allocated(qdon)) deallocate(qdon)
            allocate(qdon(npcc))

            ipcc = 1
            i1  = i
            j1  = j
            k1  = k
            jbl = ibl
            call ass_var( &
               block(jbl)%gr(igr)%ni,block(jbl)%gr(igr)%nj,block(jbl)%gr(igr)%nk, &
               i1,j1,k1, block(jbl)%gr(igr)%q,qdon(ipcc))

            ! gruppo della parete interna  
            jgruppo = nint(block(ibl)%gr(igr)%cell(i,j,k)%Chi%we(ipcc))
         end if

         call reschi(tycc,jgruppo,i,j,k, &
              block(ibl)%gr(igr)%ni,block(ibl)%gr(igr)%nj,block(ibl)%gr(igr)%nk, &
              npcc,qdon,block(ibl)%gr(igr)%cell(i,j,k), &
              block(ibl)%gr(igr)%q,block(ibl)%gr(igr)%dq)

      end do
      end do
      end do
!$omp end do !nowait
!$omp end parallel
   end do
   return

contains 
!***********************************************************************
  subroutine reschi(tycc,jgruppo,i,j,k,ni,nj,nk,npcc, &
                    qdon,cell,q,dq)
!.......................................................................
      use prec
      use var_type
      use motion_par
      implicit none

      integer(kind=I4P),intent(in)     :: tycc,jgruppo,i,j,k,ni,nj,nk
      integer(kind=I4P),intent(in)     :: npcc
      type(gen_cell)   ,intent(in)     :: cell
      type(gen_var)    ,intent(in)     :: qdon(*)
      type(gen_var)    ,intent(in)     ::  q(-1:ni+2,-1:nj+2,-1:nk+2)
      type(gen_var)    ,intent(in out) :: dq(-1:ni+2,-1:nj+2,-1:nk+2)


      integer(kind=I4P)  :: ii,jj,kk,ipcc
      real(kind=R8P)     :: ff,rul(3)
      type(gen_var)      :: fq


!
! termini forzanti per i chimera interni
!
      if (tycc==CellChim) then

         fq = q(i,j,k)
         do ipcc=1,npcc
            fq = fq - qdon(ipcc)*cell%Chi%we(ipcc)
         end do
         dq(i,j,k) = dq(i,j,k) + cell%fchi*fq

      else 

         fq = 0d0
         do kk=k-1,k+1
         do jj=j-1,j+1
         do ii=i-1,i+1
            ff = ff + 1d0
            fq = fq + q(ii,jj,kk)
         end do
         end do
         end do
         fq = q(i,j,k)-(fq-q(i,j,k))/(ff-1d0)

         ! correzione velocità
         if (glomob.and.(jgruppo/=0)) then
             call gridve(jgruppo, &
                         block(ibl)%gr(igr)%cell(i,j,k)%Center, &
                         rul)! velocita' parete interna
            rul = rul*q(i,j,k)%v(q_r)
            fq%v(q_u) = q(i,j,k)%v(q_u) - rul(1)
            fq%v(q_v) = q(i,j,k)%v(q_v) - rul(2)
            fq%v(q_w) = q(i,j,k)%v(q_w) - rul(3)
         else
            fq%v(q_u) = q(i,j,k)%v(q_u)
            fq%v(q_v) = q(i,j,k)%v(q_v)
            fq%v(q_w) = q(i,j,k)%v(q_w)
         end if

         dq(i,j,k) = dq(i,j,k) + cell%fchi*fq
      end if

      return
   end subroutine reschi

end subroutine chifor_mb

