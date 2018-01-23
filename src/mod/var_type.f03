module var_type
   use prec
   implicit none

   integer(kind=I4P) :: ivar1,ivar2

   ! indici associazione
   integer(kind=I4P),parameter    :: q_r = 1
   integer(kind=I4P),parameter    :: q_p = 1
   integer(kind=I4P),parameter    :: q_u = 2
   integer(kind=I4P),parameter    :: q_v = 3
   integer(kind=I4P),parameter    :: q_w = 4
   integer(kind=I4P),parameter    :: q_e = 5
   integer(kind=I4P)              :: q_k = 0
   integer(kind=I4P)              :: q_eps = 0
   integer(kind=I4P)              :: q_vitl = 0

   integer(kind=I4P),parameter    :: vaux_ls = 1
   integer(kind=I4P),parameter,dimension(3)    :: vaux_un = [2,3,4]

   type gen_var
       real(kind=R8P),allocatable,dimension(:) :: v
   end type gen_var

   interface operator (+) 
       procedure add_gen_var,add_gen_real,add_real_gen
   end interface operator (+) 

   interface operator (-) 
       procedure sub_gen_var,unary_gen_var
   end interface operator (-) 

   interface operator (/) 
       procedure div_gvScal,div_Scalgv
   end interface operator (/) 

   interface operator (*) 
       procedure mult_gvXgv, mult_realXgv, mult_gvXreal
   end interface operator (*) 

   interface operator (**) 
       procedure power
   endinterface operator (**) 

   interface assignment (=) 
       procedure gv_eq_scal,gv_eq_gv
   end interface assignment (=) 

contains

   pure function add_gen_var(q1,q2) result (sum_q)
       type (gen_var), intent(in) :: q1,q2
       type (gen_var)             :: sum_q
       sum_q%v(ivar1:ivar2) = q1%v(ivar1:ivar2) + q2%v(ivar1:ivar2)
   end function add_gen_var

   pure function add_gen_real(q1,q2) result (sum_q)
       type (gen_var), intent(in) :: q1
       real(kind=R8P), intent(in) :: q2
       type (gen_var)             :: sum_q
       sum_q%v(ivar1:ivar2) = q1%v(ivar1:ivar2) + q2
   end function add_gen_real
   pure function add_real_gen(q1,q2) result (sum_q)
       real(kind=R8P), intent(in) :: q1
       type (gen_var), intent(in) :: q2
       type (gen_var)             :: sum_q
       sum_q%v(ivar1:ivar2) = q2%v(ivar1:ivar2) + q1
   end function add_real_gen

   pure function sub_gen_var(q1,q2) result (diff_q)
       type (gen_var), intent(in) :: q1,q2
       type (gen_var)             :: diff_q
       diff_q%v(ivar1:ivar2) = q1%v(ivar1:ivar2) - q2%v(ivar1:ivar2)
   end function sub_gen_var

   pure function unary_gen_var(q1) result (mq)
       type (gen_var), intent(in) :: q1
       type (gen_var)             :: mq
       mq%v(ivar1:ivar2) = -q1%v(ivar1:ivar2)
   end function unary_gen_var

   pure function mult_gvXgv(q1,q2) result (qXq)
       type (gen_var), intent(in) :: q1,q2
       type (gen_var)             :: qXq
       qXq%v(ivar1:ivar2) = q1%v(ivar1:ivar2) * q2%v(ivar1:ivar2)
   end function mult_gvXgv

   pure function mult_realXgv(q1,q2) result (rXq)
       real(kind=R8P), intent(in) :: q1
       type (gen_var), intent(in) :: q2
       type (gen_var)             :: rXq
       rXq%v(ivar1:ivar2) = q1 * q2%v(ivar1:ivar2)
   end function mult_realXgv

   pure function mult_gvXreal(q1,q2) result (qXr)
       type (gen_var), intent(in) :: q1
       real(kind=R8P), intent(in) :: q2
       type (gen_var)             :: qXr
       qXr%v(ivar1:ivar2) = q2 * q1%v(ivar1:ivar2)
   end function mult_gvXreal

   pure function power(q1,n) result (q_at)
       type (gen_var)   ,intent(in) :: q1
       integer(kind=I4P),intent(in) :: n
       type (gen_var)               :: q_at
       q_at%v(ivar1:ivar2) = q1%v(ivar1:ivar2)**n
   end function power

   pure function div_gvScal(q1,q2) result (qXr)
       type (gen_var), intent(in) :: q1
       real(kind=R8P) , intent(in) :: q2

       type (gen_var)             :: qXr
       qXr%v(ivar1:ivar2) = q1%v(ivar1:ivar2)/q2
   end function div_gvScal

   pure function div_Scalgv(q1,q2) result (qXr)
       real(kind=R8P), intent(in) :: q1
       type (gen_var), intent(in) :: q2

       type (gen_var)             :: qXr
       qXr%v(ivar1:ivar2) = q1/q2%v(ivar1:ivar2)
   end function div_Scalgv

   pure subroutine gv_eq_scal(ql,qr)
       real(kind=R8P),intent(in)  :: qr
       type (gen_var),intent(out) :: ql

       ql%v(ivar1:ivar2) = qr
   end subroutine gv_eq_scal

   pure subroutine gv_eq_gv(ql,qr)
       type (gen_var),intent(in)  :: qr
       type (gen_var),intent(out) :: ql

       ql%v(ivar1:ivar2) = qr%v(ivar1:ivar2)
   end subroutine gv_eq_gv

end module var_type
