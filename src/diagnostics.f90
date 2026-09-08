module diagnostics

  use mpi
  use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
  use, intrinsic :: iso_c_binding, only : c_int
  implicit none
  private

  integer, parameter, public :: DIAG_INFO = 0
  integer, parameter, public :: DIAG_WARNING = 1
  integer, parameter, public :: DIAG_ERROR = 2
  integer, parameter, public :: DIAG_FATAL = 3

  type, public :: diagnostic_t
     integer :: severity = DIAG_INFO
     character(len=32) :: code = ''
     character(len=64) :: section = ''
     character(len=64) :: field = ''
     character(len=256) :: message = ''
     character(len=256) :: suggestion = ''
     integer :: block_id = 0
     integer :: source_rank = -1
  end type diagnostic_t

  type, public :: diagnostic_list_t
     type(diagnostic_t), allocatable :: items(:)
   contains
     procedure :: add => add_diagnostic
     procedure :: clear => clear_diagnostics
     procedure :: count => count_diagnostics
     procedure :: has_errors => diagnostics_have_errors
     procedure :: print_summary => print_diagnostic_summary
  end type diagnostic_list_t

  public :: warn_once, fatal_local, terminate_collectively

contains

  subroutine add_diagnostic(this, severity, code, message, section, field, &
                            suggestion, block_id, source_rank)
    class(diagnostic_list_t), intent(inout) :: this
    integer, intent(in) :: severity
    character(*), intent(in) :: code, message
    character(*), intent(in), optional :: section, field, suggestion
    integer, intent(in), optional :: block_id, source_rank

    type(diagnostic_t), allocatable :: expanded(:)
    integer :: n

    if (allocated(this%items)) then
       n = size(this%items)
       allocate(expanded(n+1))
       if (n > 0) expanded(1:n) = this%items
       call move_alloc(expanded, this%items)
    else
       n = 0
       allocate(this%items(1))
    end if

    this%items(n+1)%severity = severity
    this%items(n+1)%code = trim(code)
    this%items(n+1)%message = trim(message)
    if (present(section)) this%items(n+1)%section = trim(section)
    if (present(field)) this%items(n+1)%field = trim(field)
    if (present(suggestion)) this%items(n+1)%suggestion = trim(suggestion)
    if (present(block_id)) this%items(n+1)%block_id = block_id
    if (present(source_rank)) this%items(n+1)%source_rank = source_rank
  end subroutine add_diagnostic


  subroutine clear_diagnostics(this)
    class(diagnostic_list_t), intent(inout) :: this
    if (allocated(this%items)) deallocate(this%items)
  end subroutine clear_diagnostics


  integer function count_diagnostics(this, severity) result(n)
    class(diagnostic_list_t), intent(in) :: this
    integer, intent(in), optional :: severity
    integer :: i

    n = 0
    if (.not.allocated(this%items)) return
    if (.not.present(severity)) then
       n = size(this%items)
       return
    end if
    do i = 1, size(this%items)
       if (this%items(i)%severity == severity) n = n + 1
    end do
  end function count_diagnostics


  logical function diagnostics_have_errors(this) result(found)
    class(diagnostic_list_t), intent(in) :: this
    integer :: i

    found = .false.
    if (.not.allocated(this%items)) return
    do i = 1, size(this%items)
       if (this%items(i)%severity >= DIAG_ERROR) then
          found = .true.
          return
       end if
    end do
  end function diagnostics_have_errors


  subroutine print_diagnostic_summary(this, title)
    class(diagnostic_list_t), intent(in) :: this
    character(*), intent(in), optional :: title
    integer :: i

    if (present(title)) write(output_unit,'(A)') trim(title)
    if (.not.allocated(this%items)) then
       write(output_unit,'(A)') 'Summary: 0 errors, 0 warnings.'
       return
    end if

    do i = 1, size(this%items)
       call print_diagnostic(this%items(i), output_unit)
    end do
    write(output_unit,'(A,I0,A,I0,A)') 'Summary: ', &
         this%count(DIAG_ERROR)+this%count(DIAG_FATAL), ' errors, ', &
         this%count(DIAG_WARNING), ' warnings.'
    flush(output_unit)
  end subroutine print_diagnostic_summary


  subroutine print_diagnostic(item, unit)
    type(diagnostic_t), intent(in) :: item
    integer, intent(in) :: unit
    character(len=8) :: label

    select case (item%severity)
    case (DIAG_INFO)
       label = 'INFO'
    case (DIAG_WARNING)
       label = 'WARNING'
    case (DIAG_ERROR)
       label = 'ERROR'
    case default
       label = 'FATAL'
    end select

    write(unit,'(/,A,2X,A)') trim(label), trim(item%code)
    if (len_trim(item%section) > 0) write(unit,'(A,A)') '  Section: ', trim(item%section)
    if (len_trim(item%field) > 0) write(unit,'(A,A)') '  Field: ', trim(item%field)
    if (item%block_id > 0) write(unit,'(A,I0)') '  Block: ', item%block_id
    if (item%source_rank >= 0) write(unit,'(A,I0)') '  Source rank: ', item%source_rank
    write(unit,'(A,A)') '  ', trim(item%message)
    if (len_trim(item%suggestion) > 0) &
         write(unit,'(A,A)') '  Suggested fix: ', trim(item%suggestion)
  end subroutine print_diagnostic


  subroutine warn_once(code, message, routine)
    character(*), intent(in) :: code, message
    character(*), intent(in), optional :: routine
    integer :: ierr, world_rank
    logical :: initialized

    call MPI_Initialized(initialized, ierr)
    world_rank = 0
    if (initialized) call MPI_Comm_rank(MPI_COMM_WORLD, world_rank, ierr)
    if (world_rank /= 0) return

    write(output_unit,'(/,A,2X,A)') 'WARNING', trim(code)
    if (present(routine)) write(output_unit,'(A,A)') '  Routine: ', trim(routine)
    write(output_unit,'(A,A)') '  ', trim(message)
    flush(output_unit)
  end subroutine warn_once


  subroutine fatal_local(code, message, routine, block_id, root_only)
    character(*), intent(in) :: code, message
    character(*), intent(in), optional :: routine
    integer, intent(in), optional :: block_id
    logical, intent(in), optional :: root_only
    integer :: ierr, world_rank
    logical :: initialized, print_here

    call MPI_Initialized(initialized, ierr)
    world_rank = 0
    if (initialized) call MPI_Comm_rank(MPI_COMM_WORLD, world_rank, ierr)

    print_here = .true.
    if (present(root_only)) print_here = .not.root_only .or. world_rank == 0
    if (print_here) then
       write(error_unit,'(/,A,2X,A)') 'FATAL', trim(code)
       write(error_unit,'(A,I0)') '  Rank: ', world_rank
       if (present(block_id)) write(error_unit,'(A,I0)') '  Block: ', block_id
       if (present(routine)) write(error_unit,'(A,A)') '  Routine: ', trim(routine)
       write(error_unit,'(A,A)') '  ', trim(message)
       flush(error_unit)
    end if

    if (initialized) then
       call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    else
       error stop 1
    end if
  end subroutine fatal_local


  subroutine terminate_collectively(exit_code)
    integer, intent(in), optional :: exit_code
    interface
      subroutine c_exit(status) bind(C, name='exit')
        import :: c_int
        integer(c_int), value :: status
      end subroutine c_exit
    end interface

    integer :: code, ierr, world_rank
    logical :: initialized

    code = 1
    if (present(exit_code)) code = max(1, exit_code)
    call MPI_Initialized(initialized, ierr)
    if (initialized) then
       call MPI_Comm_rank(MPI_COMM_WORLD, world_rank, ierr)
       if (world_rank == 0) flush(output_unit)
       call MPI_Barrier(MPI_COMM_WORLD, ierr)
       call MPI_Finalize(ierr)
       call c_exit(int(code, c_int))
    else
       call c_exit(int(code, c_int))
    end if
  end subroutine terminate_collectively

end module diagnostics
