function compressible_flux(idir,cell,vaux,ql,qr) result(F_comp)
   use prec
   use work_var
   use physical_par
   use var_type
   use eqstat_comp

   implicit none

    integer(kind=I4P),intent(in) :: idir
    type(gen_cell)   ,intent(in) :: cell
    type(gen_var)    ,intent(in) :: vaux
    type(gen_var)    ,intent(in) :: ql,qr
    type(gen_var)                :: F_comp

    real(kind=R8P)               :: unl,unr
    real(kind=R8P)               :: unm,ugc
    real(kind=R8P)               :: rr,rl
    real(kind=R8P)               :: a2l,a2r
    real(kind=R8P)               :: sqrl,sqrr,aal,aar
    real(kind=R8P)               :: erl,err
    real(kind=R8P)               :: hl,hr,hm
    real(kind=R8P)               :: u2m
    real(kind=R8P)               :: am,am2,ar,al,amA
    real(kind=R8P)               :: Sr,Sl,Srl,DS
    type(gen_var)                :: fl,fr
    
    ! caso compressibile 
    ! q(1) = densità
    ! q(2) = ru
    ! q(3) = rv
    ! q(4) = rw
    ! q(e) = energia totale per unità di volume
    ! q(6...n) = variabili turbolente


    unl  = ql%v(q_u)*cell%Sn(1,idir) &
          +ql%v(q_v)*cell%Sn(2,idir) &
          +ql%v(q_w)*cell%Sn(3,idir)
    unr  = qr%v(q_u)*cell%Sn(1,idir) &
          +qr%v(q_v)*cell%Sn(2,idir) &
          +qr%v(q_w)*cell%Sn(3,idir)

    rr   = qr%v(q_r)
    rl   = ql%v(q_r)

    ! flussi di volume normali u*area
    unl = unl/rl
    unr = unr/rr

    ! aggiunta velocità di griglia
    ugc  = vaux%v( vaux_un(idir) ) ! ugc in ingresso: la velocità di griglia
    unl  = unl - ugc               ! normale alla faccia * A
    unr  = unr - ugc

    ! velocità del suono ^ 2
    a2l = velsuo_2(gam,gamm1,ql)
    a2r = velsuo_2(gam,gamm1,qr)

    sqrl = sqrt(rl)
    sqrr = sqrt(rr)

    aal = sqrl/(sqrl+sqrr)
    aar = 1d0 - aal
   
    a2l = a2l*invgam   ! a^2/gamma
    a2r = a2r*invgam

    erl = qr%v(q_e)/rl  ! e/rho
    err = ql%v(q_e)/rr
      
    hl  = a2l + erl    ! H
    hr  = a2r + err    ! H

    ! entalpia totale media
    hm  = aal*hl + aar*hr

    ! velocita' ^2 /2 media
    ! u^2 = e/rho - p/rho/(gamma-1) 
    !     = e/eho - (gamma*p/rho)/gamma/(gamma-1)
    !     = e/eho -      (a^2/gamma)   /(gamma-1)
    u2m  = aal*(erl-a2l*invgamm1) + aar*(err-a2r*invgamm1)
    u2m  = max(0d0,u2m)

    ! velocità del suono ^2 media
    am2 = gamm1*max(1d-8,hm-0.5d0*u2m)
    am  = sqrt(am2)

    ar = sqrt(a2r)
    al = sqrt(a2l)

    unm = aal*unl + aar*unr
    amA = am*cell%Area(idir)
    Sr = max(0d0,unm+amA,unr+ar*cell%Area(idir))
    Sl = min(0d0,unm-amA,unl-al*cell%Area(idir))

    F_comp = 0d0

    if (Sl.gt.0d0) then
       call Fcomp_face(cell%Sn(1:3,idir),unl,ql,F_comp)
    else if (Sr.lt.0d0) then
       call Fcomp_face(cell%Sn(1:3,idir),unr,qr,F_comp)
    else
       DS = 1d0/max(zero,(Sr-Sl))
       aal =  Sr*DS
       aar = -Sl*DS
       Srl = Sr*Sl*DS
       call Fcomp_face(cell%Sn(1:3,idir),unl,ql,fl)
       call Fcomp_face(cell%Sn(1:3,idir),unr,qr,fr)

       F_comp = aal*fl + aar*fl + Srl*(qr-ql)
    end if
end function compressible_flux

subroutine natural_BC_compressible(tycc,ibl,igr,jfcc,i0,j0,k0,i1,j1,k1)
    use prec
    use work_var
    use mpi_var
    use solution_par
    use concon
    use physical_par
    use eqstat_comp
    implicit none

    integer(kind=I4P), intent(in out) :: tycc
    integer(kind=I4P), intent(in)     :: ibl,igr,jfcc,i0,j0,k0,i1,j1,k1

    integer(kind=I4P)             :: id,jd,kd,idir
    real(kind=R8P)                :: nor(3),nnorm,unorm,uenf,eenf,mach

    if (tycc==Uscita.or.tycc==Ingresso) then ! uscita o ingresso
       id = (i0+i1)/2
       jd = (j0+j1)/2
       kd = (k0+k1)/2
       idir = (jfcc+1)/2

       nor(1:3) = block(ibl)%gr(igr)%cell(id,jd,kd)%Sn(1:3,idir)
       nnorm = norm2(nor) + zero
       nor = nor/nnorm
       unorm = uinf.dot.nor
       unorm =  unorm/block(ibl)%gr(igr)%q(i1,j1,k1)%v(q_r)
              
       if (unorm<(-1d-20)) then
          tycc = Ingresso
       else
          tycc = Uscita
       end if

       mach = unorm/sqrt(velsuo_2(gam,gamm1,block(ibl)%gr(igr)%q(i1,j1,k1)))
    else 
       mach = 0d0
    end if 
!.......................................................................
! Condizioni al Contorno
!.......................................................................
    select case (tycc)
       !................................................................
       case (Parete,PareteTrasl,PareteNonAt)  ! Adiabatiche
       !................................................................
          ! qwal da calcolare all'inizio
          ! if (mobile(ibl) ) ricalcola uwal

          block(ibl)%gr(igr)%q(i0,j0,k0) = qwal - block(ibl)%gr(igr)%q(i1,j1,k1)

          ! densità
          block(ibl)%gr(igr)%q(i0,j0,k0)%v(q_r) = block(ibl)%gr(igr)%q(i1,j1,k1)%v(q_r)
          ! energia totale
          block(ibl)%gr(igr)%q(i0,j0,k0)%v(q_e) = block(ibl)%gr(igr)%q(i1,j1,k1)%v(q_e)

          if (q_eps > 0) &  ! two equation model - epsilon
          block(ibl)%gr(igr)%q(i0,j0,k0)%v(q_eps) = block(ibl)%gr(igr)%q(i1,j1,k1)%v(q_eps)
       
       !................................................................
       case (PareteIsoterma)
       !................................................................

          ! qwal da calcolare all'inizio
          ! if (mobile(ibl) ) ricalcola uwal

          block(ibl)%gr(igr)%q(i0,j0,k0) = qwal - block(ibl)%gr(igr)%q(i1,j1,k1)
          ! densità
          block(ibl)%gr(igr)%q(i0,j0,k0)%v(q_r) = block(ibl)%gr(igr)%q(i1,j1,k1)%v(q_r)

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
       case (Ingresso)  ! il check INGRESSO/USCITA si può fare una volta per tutte
       !................................................................
          block(ibl)%gr(igr)%q(i0,j0,k0) = qinf
          ! energia totale
          if (mach > -1d0) then
             block(ibl)%gr(igr)%q(i0,j0,k0)%v(q_e) = block(ibl)%gr(igr)%q(i1,j1,k1)%v(q_e)
          end if
       !................................................................
       case (Uscita)
       !................................................................
          block(ibl)%gr(igr)%q(i0,j0,k0) = block(ibl)%gr(igr)%q(i1,j1,k1)
          ! energia totale
          if (mach < 1d0) then
             block(ibl)%gr(igr)%q(i0,j0,k0)%v(q_e) = qinf%v(q_e)
          end if
       !................................................................
       case (IngrAss,InvRiemAss)
       !................................................................
            return ! non va cambiato niente
       !................................................................
       case (PressAss)
          block(ibl)%gr(igr)%q(i0,j0,k0) = block(ibl)%gr(igr)%q(i1,j1,k1)

          eenf = matcc(1,jfcc,iblglo(ibl))
          block(ibl)%gr(igr)%q(i0,j0,k0)%v(q_e) = eenf
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
          uenf = uenf*block(ibl)%gr(igr)%q(i0,j0,k0)%v(q_r)

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

end subroutine natural_BC_compressible

