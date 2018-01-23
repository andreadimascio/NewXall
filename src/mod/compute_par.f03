module compute_par
   use prec

   ! Lavoro in termini di CPU
   real(kind=R8P)    :: work,work_old
   real(kind=R8P)    ::  cpu_tot,cpu_io,cpu_rcc,cpu_com  &
                        ,cpu_residui,cpu_start,cpu_tmp
   real(kind=R8P),allocatable,dimension(:)    ::  workgr

end module compute_par
