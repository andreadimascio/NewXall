function incompressible_flux(idir,cell,vaux,ql,qr) result(F_inc)
   use prec
   use work_var
   use solution_par,only:sqbeta

   implicit none

    integer(kind=I4P),intent(in) :: idir
    type(gen_cell)   ,intent(in) :: cell
    type(gen_var)    ,intent(in) :: vaux
    type(gen_var)    ,intent(in) :: ql,qr
    type(gen_var)                :: F_inc

    real(kind=R8P)               :: unl,unr,pl,pr
    real(kind=R8P)               :: ugc
    real(kind=R8P)               :: unm,pm
    
    ! caso incompressibile 
    ! q(1) = pressione
    ! q(2) = u
    ! q(3) = v
    ! q(4) = w
    ! q(5...n) = variabili turbolente


    unl  = ql%v(2)*cell%Sn(1,idir) &
          +ql%v(3)*cell%Sn(2,idir) &
          +ql%v(4)*cell%Sn(3,idir)
    unr  = qr%v(2)*cell%Sn(1,idir) &
          +qr%v(3)*cell%Sn(2,idir) &
          +qr%v(4)*cell%Sn(3,idir)

    pr   = qr%v(1)
    pl   = ql%v(1)

    ! pressione e velocità interfaccia
    pm   = 0.5d0*( pl +  pr - sqbeta*(unr-unl)/cell%Area(idir))
    unm  = 0.5d0*(unl + unr - (pr-pl)/sqbeta*cell%Area(idir))

    ugc  = vaux%v( vaux_un(idir) )
    ugc  = unm - ugc ! ugc in ingresso contiene la velocità di griglia

    F_inc = 0d0

    ! divergenza 
    F_inc%v(1) = unm  ! o ugc
    F_inc%v(2) = pm*cell%Sn(1,idir)
    F_inc%v(3) = pm*cell%Sn(2,idir)
    F_inc%v(4) = pm*cell%Sn(3,idir)

    if (ugc > 0d0) then
       F_inc = F_inc + ugc*ql
    else
       F_inc = F_inc + ugc*qr
    end if
end function incompressible_flux

subroutine natural_BC_incompressible(tycc,ibl,igr,jfcc,i0,j0,k0,i1,j1,k1)
    use prec
    use work_var
    use mpi_var
    use solution_par
    use concon
    implicit none

    integer(kind=I4P), intent(in out) :: tycc
    integer(kind=I4P), intent(in)     :: ibl,igr,jfcc,i0,j0,k0,i1,j1,k1

    integer(kind=I4P) :: id,jd,kd,idir
    real(kind=R8P)    :: nor(3),nnorm,unorm,penf,uenf

    if (tycc==Uscita.or.tycc==Ingresso) then ! uscita o ingresso
       id = (i0+i1)/2
       jd = (j0+j1)/2
       kd = (k0+k1)/2
       idir = (jfcc+1)/2

       nor(1:3) = block(ibl)%gr(igr)%cell(id,jd,kd)%Sn(1:3,idir)
       unorm = uinf.dot.nor
              
       if (unorm<zero) then
          tycc = Ingresso
       else
          tycc = Uscita
       end if
    end if 
!.......................................................................
! Condizioni al Contorno
!.......................................................................
    select case (tycc)
       !................................................................
       case (Parete,PareteTrasl,PareteNonAt)
       !................................................................
          ! qwal da calcolare all'inizio
          ! if (mobile(ibl) ) ricalcola uwal

          block(ibl)%gr(igr)%q(i0,j0,k0) = qwal - block(ibl)%gr(igr)%q(i1,j1,k1)
          ! pressione
          block(ibl)%gr(igr)%q(i0,j0,k0)%v(q_p) = block(ibl)%gr(igr)%q(i1,j1,k1)%v(q_p)

          if (q_eps > 0) &  ! two equation model - epsilon
          block(ibl)%gr(igr)%q(i0,j0,k0)%v(q_eps) = block(ibl)%gr(igr)%q(i1,j1,k1)%v(q_eps)

       !................................................................
       case (Simmetria)
       !................................................................
          block(ibl)%gr(igr)%q(i0,j0,k0) = block(ibl)%gr(igr)%q(i1,j1,k1)

          ! velocity correction
          id = (i0+ i1)/2
          jd = (j0+ j1)/2
          kd = (k0+ k1)/2
          idir = (jfcc+1)/2

          nor(1:3) = block(ibl)%gr(igr)%cell(id,jd,kd)%Sn(1:3,idir)
          nnorm = norm2(nor) + zero
          nor = nor/nnorm

          unorm = 2d0*(block(ibl)%gr(igr)%q(i1,j1,k1)%v(q_u)*nor(1) &
                      +block(ibl)%gr(igr)%q(i1,j1,k1)%v(q_v)*nor(2) &
                      +block(ibl)%gr(igr)%q(i1,j1,k1)%v(q_w)*nor(3) &

                                                -qwal%v(q_u)*nor(1) &
                                                -qwal%v(q_v)*nor(2) &
                                                -qwal%v(q_w)*nor(3)) 

          block(ibl)%gr(igr)%q(i0,j0,k0)%v(q_u) = block(ibl)%gr(igr)%q(i1,j1,k1)%v(q_u)-nor(1)*unorm
          block(ibl)%gr(igr)%q(i0,j0,k0)%v(q_v) = block(ibl)%gr(igr)%q(i1,j1,k1)%v(q_v)-nor(2)*unorm
          block(ibl)%gr(igr)%q(i0,j0,k0)%v(q_w) = block(ibl)%gr(igr)%q(i1,j1,k1)%v(q_w)-nor(3)*unorm
       !................................................................
       case (Ingresso)
       !................................................................
          block(ibl)%gr(igr)%q(i0,j0,k0) = qinf
          ! pressione
          block(ibl)%gr(igr)%q(i0,j0,k0)%v(q_p) = block(ibl)%gr(igr)%q(i1,j1,k1)%v(q_p)
       !................................................................
       case (Uscita)
       !................................................................
          block(ibl)%gr(igr)%q(i0,j0,k0) = block(ibl)%gr(igr)%q(i1,j1,k1)
          ! pressione
          block(ibl)%gr(igr)%q(i0,j0,k0)%v(q_p) = -block(ibl)%gr(igr)%q(i1,j1,k1)%v(q_p)
       !................................................................
       case (IngrAss,InvRiemAss)
       !................................................................
            return ! non va cambiato niente
       !................................................................
       case (PressAss)
          block(ibl)%gr(igr)%q(i0,j0,k0) = block(ibl)%gr(igr)%q(i1,j1,k1)

          penf = matcc(1,jfcc,iblglo(ibl))
          block(ibl)%gr(igr)%q(i0,j0,k0)%v(q_p) = penf
       !................................................................
       case (UnormAss)
       !................................................................
          block(ibl)%gr(igr)%q(i0,j0,k0) = block(ibl)%gr(igr)%q(i1,j1,k1)

          ! velocity correction
          id = (i0+ i1)/2
          jd = (j0+ j1)/2
          kd = (k0+ k1)/2
          idir = (jfcc+1)/2

          nor(1:3) = block(ibl)%gr(igr)%cell(id,jd,kd)%Sn(1:3,idir)

          nnorm = norm2(nor) + zero
          nor = -nor/nnorm ! perché Sn è la  normale uscente

          ! uenf > 0 ==> flusso entrante nel dominio
          uenf = matcc(1,jfcc,iblglo(ibl))
          unorm = 2d0*(uenf                     +qwal%v(q_u)*nor(1) &
                                                +qwal%v(q_v)*nor(2) &
                                                +qwal%v(q_w)*nor(3) &
                      -block(ibl)%gr(igr)%q(i1,j1,k1)%v(q_u)*nor(1) &
                      -block(ibl)%gr(igr)%q(i1,j1,k1)%v(q_v)*nor(2) &
                      -block(ibl)%gr(igr)%q(i1,j1,k1)%v(q_w)*nor(3)) 
          block(ibl)%gr(igr)%q(i0,j0,k0)%v(q_u) = block(ibl)%gr(igr)%q(i1,j1,k1)%v(q_u)+nor(1)*unorm
          block(ibl)%gr(igr)%q(i0,j0,k0)%v(q_v) = block(ibl)%gr(igr)%q(i1,j1,k1)%v(q_v)+nor(2)*unorm
          block(ibl)%gr(igr)%q(i0,j0,k0)%v(q_w) = block(ibl)%gr(igr)%q(i1,j1,k1)%v(q_w)+nor(3)*unorm
       !................................................................
       case (Estrapol,EstrAlter)
       !................................................................
          block(ibl)%gr(igr)%q(i0,j0,k0) = block(ibl)%gr(igr)%q(i1,j1,k1)
       !................................................................
       case default
       !................................................................
           write(*,*) ' ATTENZIONE !!!!!'
           write(*,*) ' Condizione sconosciuta : ',tycc
    end select

end subroutine natural_BC_incompressible

