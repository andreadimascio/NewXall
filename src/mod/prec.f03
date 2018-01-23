module prec
   implicit none
   save; public
   
   integer, parameter:: R8P  = selected_real_kind(15,307)
   integer, parameter:: R4P  = selected_real_kind( 6,37)
   integer, parameter:: I4P  = selected_int_kind(9)
   integer, parameter:: I1P  = selected_int_kind(1)

   real(kind=R8P),parameter :: Infinity=1d50, Zero=1d-50

   interface operator(.smaller.)
       module procedure smaller_vec
   end interface
   interface operator(.larger.)
       module procedure larger_vec
   end interface

   interface operator(.dot.)
       module procedure scalar_prod
   end interface
   interface operator(.cross.)
       module procedure cross_prod
   end interface

contains
   function scalar_prod(A,B) result(C)
      implicit none
      real(kind=R8P),dimension(3),intent(in)  :: A,B
      real(kind=R8P)                          :: C

      C = A(1)*B(1) + A(2)*B(2) + A(3)*B(3)
   end function scalar_prod

   function cross_prod(A,B) result(C)
      implicit none
      real(kind=R8P),dimension(3),intent(in)  :: A,B
      real(kind=R8P),dimension(3)             :: C

      C(1) = A(2)*B(3) - A(3)*B(2) 
      C(2) = A(3)*B(1) - A(1)*B(3) 
      C(3) = A(1)*B(2) - A(2)*B(1) 
   end function cross_prod

   function smaller_vec(A,B) result(C)
      implicit none
      real(kind=R8P),dimension(3),intent(in)  :: A,B
      logical                                 :: C

      C = (A(1) < B(1)) .or. &
          (A(2) < B(2)) .or. &
          (A(3) < B(3)) 
   end function smaller_vec

   function larger_vec(A,B) result(C)
      implicit none
      real(kind=R8P),dimension(3),intent(in)  :: A,B
      logical                                 :: C

      C = (A(1) > B(1)) .or. &
          (A(2) > B(2)) .or. &
          (A(3) > B(3)) 
   end function larger_vec

end module prec
