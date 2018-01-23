module eqstat_comp
contains

    subroutine Fcomp_face(Sn,unA,q,f)
       use prec
       use var_type
       use physical_par
       implicit none

       real(kind=R8P), intent(in) :: unA
       real(kind=R8P), intent(in) :: Sn(3)
       type(gen_var) , intent(in) :: q
       type(gen_var) , intent(out):: f

       real(kind=R8P)             :: p,press

       p = press(gamm1,q)
       f = unA*q  ! convezione

       f%v(q_u) = f%v(q_u) + p*Sn(1)
       f%v(q_v) = f%v(q_v) + p*Sn(2)
       f%v(q_w) = f%v(q_w) + p*Sn(3)
       f%v(q_e) = f%v(q_e) + p*unA
    end subroutine Fcomp_face

    function press(gamm1,q) result(p)
       use prec
       use var_type
       implicit none
    
       real(kind=R8P),intent(in) :: gamm1
       type(gen_var) ,intent(in) :: q
       real(kind=R8P)            :: p
    
       p = gamm1*(q%v(q_e) - 0.5d0*(q%v(q_u)*q%v(q_u) &
                                  + q%v(q_v)*q%v(q_v) &
                                  + q%v(q_w)*q%v(q_w)) &
                                   /q%v(q_r))
       p = abs(p) 
    end function press

    function velsuo_2(gam,gamm1,q) result(a2)
       use prec
       use var_type
       implicit none
    
       real(kind=R8P),intent(in) :: gam,gamm1
       type(gen_var) ,intent(in) :: q
       real(kind=R8P) :: a2
    
       real(kind=R8P) :: r,p
       real(kind=R8P) :: press
    
       r = abs(q%v(q_r))+zero
       p = press(gamm1,q)
       a2 = gam*p/r
    end function velsuo_2
end module eqstat_comp
