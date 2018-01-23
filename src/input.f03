subroutine input_par()
   use prec
   use parser
   use physical_par
   use file_set
   use solution_par
   use proc_pointers

   implicit none
   integer(kind=I4P)  :: ifile_par,ifile_echo
   integer(kind=I4P)  :: ierr
   character(len=130) :: input_line
   character(len=80 ) :: flow_type
   character(len=3  ) :: dis


   ! read input
   call getunit(ifile_par)
   open(unit=ifile_par,file='Xall.input',form='formatted', &
        status='old',action='read',position='rewind')

   ! defaults
   filevo = 'EVOL/'       ! evolution directory
   filgrd = 'CC/cc'       ! grid
   filecc = 'CC/cc'       ! boundary conditions 
   filini = 'INI/ini'     ! initial conditions
   filrst = 'RST/rst'     ! restart file
   filout = 'OUT/EXsol'   ! extraction file

   dis = "tvd"        ! discretization
   rey = infinity
   flow_type = "incompressible"
   gam       = 1.4d0  ! specific heat ratio
   sqbeta    = 1d0

   do 
      read(ifile_par,'(a)',iostat=ierr) input_line
      if (ierr /= 0) exit  ! end file reached

      call readline(input_line,ifile_par)

      !****************************************************
      ! Physical parameters
      !****************************************************
      ! equation type
      call read_data(input_line,"flow_type",flow_type)
      ! Reynolds number
      call read_data(input_line,"reynolds",rey)
      ! Specific heat ratio
      call read_data(input_line,"gamma",gam)

      !****************************************************
      ! Computation parameter
      !****************************************************
      ! pseudo sound speed
      call read_data(input_line,"pseudo sound",sqbeta)
      ! discretizzazione
      call read_data(input_line,"discretization",dis)

      !****************************************************
      ! File Set
      !****************************************************
      ! root for grid files
      call read_data(input_line,"grid",filgrd)
      ! root for boundary condition files
      call read_data(input_line,"boundary",filecc)
      ! root for initial condition
      call read_data(input_line,"initial",filini)
      ! root for restart file
      call read_data(input_line,"restart",filrst)
      ! root for restart file
      call read_data(input_line,"extraction",filout)
      ! directory for history files
      call read_data(input_line,"evolution",filevo)
      
   end do

   close(ifile_par)

   ! *******************************************************
   ! valori acquisiti
   ! *******************************************************
   rrey = 1_R8P/rey

   gamm1    = gam - 1d0
   invgamm1 = 1d0/gamm1
   invgam   = 1d0/gam

   flow_type = to_lower(flow_type)
   if (flow_type == "incompressible" ) then
      actual_flux => incompressible_flux
      natural_boundary => natural_BC_incompressible
   else if (flow_type == "compressible" ) then
      actual_flux => compressible_flux
      natural_boundary => natural_BC_compressible
   end if

   ! assegnazione puntatori a procedure
   if (dis=="god") then
      lrstat => lrgo1
   else if (dis=="tvd") then
      lrstat => lrtvd
   else if (dis=="up3") then
      lrstat => lr3we   ! lr3up  DA IMPLEMENTARE
   else if (dis=="we3") then
      lrstat => lr3we
   else if (dis=="ce4") then
      lrstat => lr3we  ! lr4ce DA IMPLEMENTARE
   end if
   lrstat_finest => lrstat

   ! echo file
   call getunit(ifile_echo)
   open(unit=ifile_echo,file='Xall.input.echo',form='formatted', &
        status='unknown',action='write',position='rewind')

   write(ifile_echo,'(a      )') ' ======================'
   write(ifile_echo,'(a      )') ' Physical parameter   :'
   write(ifile_echo,'(a      )') ' ===================== '
   write(ifile_echo,'(a      )') '                       '
   write(ifile_echo,'(a,a    )') ' Flow type           = ',flow_type
   write(ifile_echo,'(a,g20.8)') ' reynolds number     = ',rey
   write(ifile_echo,'(a,g20.8)') ' Specific heat ratio = ',gam

   write(ifile_echo,'(a      )') '                       '
   write(ifile_echo,'(a      )') ' ======================'
   write(ifile_echo,'(a      )') ' Discretization terms :'
   write(ifile_echo,'(a      )') ' ======================'
   write(ifile_echo,'(a      )') '                       '
   write(ifile_echo,'(a,a    )') ' discretization      = ',dis
   write(ifile_echo,'(a,a    )') ' pseudo sound speed  = ',sqbeta

   write(ifile_echo,'(a      )') '                       '
   write(ifile_echo,'(a      )') ' ======================'
   write(ifile_echo,'(a      )') ' File set             :'
   write(ifile_echo,'(a      )') ' ======================'
   write(ifile_echo,'(a      )') '                       '
   write(ifile_echo,'(a,a)'    ) ' Grid            = ',filgrd
   write(ifile_echo,'(a,a)'    ) ' Boundary Cond.  = ',filecc
   write(ifile_echo,'(a,a)'    ) ' Initial  Cond.  = ',filini
   write(ifile_echo,'(a,a)'    ) ' Restart         = ',filrst
   write(ifile_echo,'(a,a)'    ) ' Evolution       = ',filevo
   write(ifile_echo,'(a,a)'    ) ' Extraction      = ',filout

   close(ifile_echo)
end subroutine input_par

!***********************************************************************
! read grid dimension and block allocation
!***********************************************************************
subroutine input_blocks()
   use prec
   use file_set
   use compute_par
   use mpi_var
   use work_var
   use mpi

   implicit none

   integer(kind=I4P)  :: fgrd
   integer(kind=I4P)  :: ierr,idum
   integer(kind=I4P)  :: ibl,igr,ni,nj,nk

   character(len=130) :: finame
   character(len=3)   :: nprc
!
!  lettura numero e dimensioni blocchi dal file reticolo fitto
!
   write(finame,'(a)') trim(filecc)//'.01.grd'

   if (nproc/=1) then
      write(nprc,'(i3.3)') myrank
      finame = trim(finame)//'.p'//nprc
   end if

   call getunit(fgrd)
   open(fgrd,file=finame,status='old',form='unformatted', &
        action='read',position='rewind')

   read(fgrd) nbl

   idum = nbl
   if (nproc/=1) call MPI_ALLREDUCE(nbl,idum,1,MPI_INTEGER, &
                            MPI_SUM,MPI_COMM_WORLD,ierr)
   if (idum/=nbltot) then
      if (myrank==0) then
         write(*,*) 'Attenzione numero blocchi totali'
         write(*,*) 'nei file di reticolo non coincide'
         write(*,*) 'con quanto letto in proc.input'
      end if
      call MPI_FINALIZE(ierr)
      stop
   end if

   ! lettura dimensioni e allocazione blocchi
   allocate(block(nbltot))
   do ibl = 1,nbl
      read(fgrd) ni,nj,nk
      call block(ibl)%create_B(ni,nj,nk)
   end do

   ! lettura griglia fitta
   igr = 1
   do ibl=1,nbl
      call input_x(myrank,fgrd, &
           block(ibl)%gr(igr)%ni,block(ibl)%gr(igr)%ni,block(ibl)%gr(igr)%nk, &
           block(ibl)%gr(igr)%cell)
   end do

contains

   subroutine input_x(myrank,fgrd,ni,nj,nk,cell)
      use prec
      use var_type
      use mpi
      implicit none

      integer(kind=I4P),intent(in)     :: fgrd
      integer(kind=I4P),intent(in)     :: ni,nj,nk,myrank
      integer(kind=I4P)                :: ierr1,ierr2,ierr3
      integer(kind=I4P)                :: i,j,k
      type(gen_cell),intent(in out)    :: cell(-2:ni+2,-2:nj+2,-2:nk+2)

      ierr1 = 0
      ierr2 = 0
      ierr3 = 0
      read(fgrd,iostat=ierr1)(((cell(i,j,k)%Vertex(1),i=-2,ni+2),j=-2,nj+2),k=-2,nk+2)
      read(fgrd,iostat=ierr2)(((cell(i,j,k)%Vertex(2),i=-2,ni+2),j=-2,nj+2),k=-2,nk+2)
      read(fgrd,iostat=ierr3)(((cell(i,j,k)%Vertex(3),i=-2,ni+2),j=-2,nj+2),k=-2,nk+2)

      write(*,*) 'Myid = ', myrank, ' Errore lettura in input_x : ' &
                 ,ierr1,ierr2,ierr3
      call MPI_ABORT(MPI_COMM_WORLD,ierr1,ierr2)

   end subroutine input_x

end subroutine input_blocks
!***********************************************************************
subroutine input_cc()
!***********************************************************************
   use prec
   use file_set
   use compute_par
   use work_var
   use mpi_var
   use mpi

   implicit none

   integer(kind=I4P)  :: ficc
   integer(kind=I4P)  :: ibl,igr
   integer(kind=I4P)  :: i,j,k
   integer(kind=I4P)  :: ndim_rcc
   integer(kind=I4P)  :: ierr
   character(len=130) :: str
   real(kind=R4P),allocatable,dimension(:) :: rcc

   do igr=1,ngr

      write(str,'(a,i2.2,a,i3.3)') trim(filecc)//'.',igr
      if (nproc/=1) then
         write(str,'(a,i2.2,a,i3.3)') trim(filecc)//'.',igr,'.p',myrank
      end if

      call getunit(ficc)
      open(ficc,file=str,status='old',form='unformatted')
      rewind(ficc)

      ! lettura (dummy) del numero blocchi locali
      read(ficc) i

      ! lettura (dummy) delle dimensioni dei blocchi
      do ibl=1,nbl
         read(ficc) i,j,k
      end do

      do ibl=1,nbl
         call input_chi(ficc, &
              block(ibl)%gr(igr)%ni,block(ibl)%gr(igr)%ni,block(ibl)%gr(igr)%nk, &
              block(ibl)%gr(igr)%cell)
      end do

      read(ficc) ndim_rcc
      if (allocated(rcc)) deallocate(rcc,stat=ierr)
      allocate(rcc(ndim_rcc),stat=ierr)
      read(ficc,iostat=ierr) (rcc(i),i=1,ndim_rcc)

      call input_rcc(rcc(1), &
           block(ibl)%gr(igr)%ni,block(ibl)%gr(igr)%ni,block(ibl)%gr(igr)%nk, &
           block(ibl)%gr(igr)%cell)

      deallocate(rcc,stat=ierr)


      close(ficc)
   end do

   return

contains

   subroutine input_chi(ficc,ni,nj,nk,cell)
      use prec
      use cell_type

      implicit none

      integer(kind=I4P),intent(in)  :: ficc
      integer(kind=I4P),intent(in)  :: ni,nj,nk
      type(gen_cell),intent(in out) :: cell(-1:ni+2,-1:nj+2,-1:nk+2)

      integer(kind=I4P)             :: i,j,k

      read(ficc) (((cell(i,j,k)%Chi%tycc,i=-1,ni+2),j=-1,nj+2),k=-1,nk+2)
      ! su tycc temporanemente il vecchio icc

      return
   end subroutine input_chi

   subroutine input_rcc(rcc,ni,nj,nk,cell)
      use prec
      use cell_type

      implicit none

      real(kind=R4P)   ,intent(in)  :: rcc(*)
      integer(kind=I4P),intent(in)  :: ni,nj,nk
      type(gen_cell),intent(in out) :: cell(-1:ni+2,-1:nj+2,-1:nk+2)

      integer(kind=I4P)             :: i,j,k
      integer(kind=I4P)             :: pvec,tycc,npcc
      integer(kind=I4P)             :: n,offset

      do k=-1,nk+2
      do j=-1,nj+2
      do i=-1,ni+2

         pvec = cell(i,j,k)%Chi%tycc  ! puntatore a rcc
         if (pvec < 1) cycle

         tycc = nint(rcc(pvec))
         pvec = pvec + 1
         npcc = nint(rcc(pvec))
         pvec = pvec + 1

         cell(i,j,k)%Chi%tycc = tycc  ! tipo cella
         cell(i,j,k)%Chi%npcc = npcc  ! numero donatori
         if (npcc < 1) cycle

         call cell(i,j,k)%Chi%add_X(npcc)
         do n = 1,npcc
            offset = pvec + 5*(n-1)
            cell(i,j,k)%Chi%bd(n) = nint(rcc(offset+1))
            cell(i,j,k)%Chi%id(n) = nint(rcc(offset+2))
            cell(i,j,k)%Chi%jd(n) = nint(rcc(offset+3))
            cell(i,j,k)%Chi%kd(n) = nint(rcc(offset+4))
            cell(i,j,k)%Chi%we(n) =      rcc(offset+5)
         end do
    
      end do
      end do
      end do


      return
   end subroutine input_rcc

end subroutine input_cc
