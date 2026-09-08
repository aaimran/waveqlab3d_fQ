module elastic

  !> elastic module has routine to set rates using elastic wave equation

  !> this routine could be in fields.f90, but it is instead separated out in order to have separate modules
  !> for each physics: elastic.f90, acoustic.f90, etc., instead of a general rates.f90

  implicit none

contains

  !> set rates in block using finite difference discretization of elastic wave equation
  subroutine set_rates_elastic(F, G, M, type_of_mesh)

    use datatypes, only: block_type, block_grid_t, block_material
    use RHS_Interior, only : RHS_Center, RHS_near_boundaries
    implicit none

    type(block_type), intent(inout):: F
    type(block_grid_t), intent(in) :: G
    type(block_material), intent(inout) :: M

    character(len=64), intent(in) :: type_of_mesh

    call RHS_Center(F, G, M, type_of_mesh)

    call RHS_Near_Boundaries(F, G, M)

  end subroutine set_rates_elastic


end module elastic
