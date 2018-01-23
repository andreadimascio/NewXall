module proc_pointers
    use prec
    use cell_type
    use var_type

    implicit none

    procedure(flux), pointer :: actual_flux => null()
    procedure(flux)          :: incompressible_flux
    procedure(flux)          ::   compressible_flux
    interface 
       function flux(idir,cell,vaux,q1,q2) result(F)
          import I4P,gen_cell,gen_var
          integer(kind=I4P),intent(in) :: idir
          type(gen_cell)   ,intent(in) :: cell
          type(gen_var)    ,intent(in) :: vaux
          type(gen_var)    ,intent(in) :: q1,q2
          type(gen_var)                :: F
       end function flux
    end interface

    procedure(left_right_states), pointer :: lrstat => null()
    procedure(left_right_states), pointer :: lrstat_finest => null()
    procedure(left_right_states)          :: lrgo1,lrtvd,lr3up,lr3we,lr4ce

    interface
       subroutine left_right_states(n,q,ql,qr)
           import I4P,gen_var
           integer(kind=I4P),intent(in ) :: n
           type(gen_var)    ,intent(in ) :: q(-1:)
           type(gen_var)    ,intent(out) :: ql(0:),qr(0:)
       end subroutine left_right_states
    end interface

    procedure(natural_bc), pointer :: natural_boundary => null()
    procedure(natural_bc)          :: natural_BC_incompressible 
    procedure(natural_bc)          :: natural_BC_compressible 
    interface 
       subroutine natural_bc(tycc,ibl,igr,jfcc,i,j,k,id,jd,kd)
           import I4P
          integer(kind=I4P),intent(in) :: tycc,ibl,igr,jfcc,i,j,k,id,jd,kd
       end subroutine natural_bc
    end interface
 
end module proc_pointers
