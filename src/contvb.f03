!=======================================================================
!
!    condizioni al contorno
!
!         Se condizione fisica :
!
!         rcc(icc(i,j,k)) = 1      parete
!         rcc(icc(i,j,k)) = 2      scorrimento o simmetria
!         rcc(icc(i,j,k)) = 3      ingresso
!         rcc(icc(i,j,k)) = 4      ingresso o uscita
!         rcc(icc(i,j,k)) = 5      ingresso assegnato
!         rcc(icc(i,j,k)) = 6      pressione assegnata
!         rcc(icc(i,j,k)) = 7      velocita'' normale entrante assegnata
!         rcc(icc(i,j,k)) = 8      Inv. Riemann assegnato
!         rcc(icc(i,j,k)) = 9      tutto estrapolato
!         rcc(icc(i,j,k)) = 10     parete traslante con vel. Uinf
!         rcc(icc(i,j,k)) = 11     parete non attiva
!         rcc(icc(i,j,k)) = 19     tutto estrapolato (alternativa)
!
!=======================================================================
subroutine cont_mb(igr)
!.......................................................................
    use prec
    use work_var
    use mpi_var

    implicit none

    integer(kind=I4P),intent(in)   :: igr

    integer(kind=I4P)   :: ibl,i,j,k,i1,j1,k1

!.......................................................................
! Comunicazione informazioni tra processori
!.......................................................................
    if (nproc/=1) call comm_var(igr)
!.......................................................................
    do ibl=1,nbl

       do k=block(ibl)%gr(igr)%nk/2,-1,-1
          do j=block(ibl)%gr(igr)%nj/2,-1,-1
             do i=block(ibl)%gr(igr)%ni/2,-1,-1
                call generic_boundary(i,j,k,ibl,igr) 
             
                k1 = block(ibl)%gr(igr)%nk - k + 1
                j1 = block(ibl)%gr(igr)%nj - j + 1
                i1 = block(ibl)%gr(igr)%ni - i + 1

                call generic_boundary(i1,j1,k1,ibl,igr)
             end do
          end do
       end do
    end do

    return
end subroutine cont_mb

!-----------------------------------------------------------------------
subroutine generic_boundary(i,j,k,ibl,igr)
!-----------------------------------------------------------------------
    use prec
    use work_var
    use proc_pointers
    use solution_par
    use motion_par
    use concon

    implicit none

    integer(kind=I4P),intent(in)   :: i,j,k,igr,ibl

    logical             :: ccnat
    integer(kind=I4P)   :: ipcc,tycc,npcc
    integer(kind=I4P)   :: jbl,id,jd,kd,jfcc
    integer(kind=I4P)   :: jgruppo
    real(kind=R8P)      :: we,u(3)
    type(gen_var)       :: vd

    ! condizioni  contorno nel punto
    tycc = block(ibl)%gr(igr)%cell(i,j,k)%Chi%tycc

    if (tycc == 0) return ! punto regolare

    ccnat = tycc < 0  ! condizione al contorno naturale

    if (ccnat) then ! condizioni al contorno naturali
                    ! MODIFICARE IN INPUT

       jfcc = block(ibl)%gr(igr)%cell(i,j,k)%Chi%bd(1)
       id   = block(ibl)%gr(igr)%cell(i,j,k)%Chi%id(1)
       jd   = block(ibl)%gr(igr)%cell(i,j,k)%Chi%jd(1)
       kd   = block(ibl)%gr(igr)%cell(i,j,k)%Chi%kd(1)

       call natural_boundary(tycc,ibl,igr,jfcc,i,j,k,id,jd,kd)

    else ! chimera generico

       npcc = block(ibl)%gr(igr)%cell(i,j,k)%Chi%npcc

       vd = 0d0
       do ipcc=1,npcc
          jbl = block(ibl)%gr(igr)%cell(i,j,k)%Chi%bd(1)
          id  = block(ibl)%gr(igr)%cell(i,j,k)%Chi%id(1)
          jd  = block(ibl)%gr(igr)%cell(i,j,k)%Chi%jd(1)
          kd  = block(ibl)%gr(igr)%cell(i,j,k)%Chi%kd(1)
          we  = block(ibl)%gr(igr)%cell(i,j,k)%Chi%we(1)

          vd = vd + block(jbl)%gr(igr)%q(id,jd,kd)*we
       end do

       if ((tycc == CellChim  ) .or. &
           (tycc == CellParInt) ) then  ! chimera interna
          if (.not.finest) return

          if (tycc == CellParInt) then ! parete interna
             ! gruppo della parete interna  
             jgruppo = nint(block(ibl)%gr(igr)%cell(i,j,k)%Chi%we(ipcc))
             call gridve(jgruppo, &
                         block(ibl)%gr(igr)%cell(i,j,k)%Center, &
                         u)! velocita' parete interna
             vd%v(q_u) = u(1)
             vd%v(q_v) = u(2)
             vd%v(q_w) = u(3)  ! ricorda densità se compressibile o bifase
          end if

          block(ibl)%gr(igr)%dq(i,j,k) = block(ibl)%gr(igr)%cell(i,j,k)%fchi &
                                       * (block(ibl)%gr(igr)%q(i,j,k)-vd)

       else ! chimera di contorno

          block(ibl)%gr(igr)%q(i,j,k) = vd
          ! (SI DEVONO CAMBIARE DONATORI PER ELIMINARE L'IF)
          ! if ((tycc>=FacXf_i0).and.(tycc<=FacXf_kn))  then

       end if

    end if

    return
end subroutine generic_boundary

