! Copyright (C) 2022 National Center for Atmospheric Research,
! National Technology & Engineering Solutions of Sandia, LLC (NTESS),
! and the U.S. Environmental Protection Agency (USEPA)
!
! SPDX-License-Identifier: Apache-2.0
!
module aero_model

  use aero_constants,              only : real_kind

  implicit none
  private

  public :: model_t, model_ptr

  type, abstract :: model_t
  contains
    procedure(model_name),     deferred :: name
    procedure(create_state),   deferred :: create_state
    procedure(optics_grid),    deferred :: optics_grid
    procedure(compute_optics), deferred :: compute_optics
  end type model_t

  type :: model_ptr
    class(model_t), pointer :: ptr_ => null( )
  contains
    final :: model_ptr_finalize
  end type model_ptr

interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the name of the aerosol model/package
  function model_name( this )

    import :: model_t

    !> Unique name for the aerosol model
    character(len=:), allocatable :: model_name
    !> Aerosol model
    class(model_t), intent(in) :: this

  end function model_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns an newly created aerosol state
  function create_state( this ) result( state )

    use aero_state,                    only : state_t
    import :: model_t

    !> New aerosol state
    class(state_t), pointer    :: state
    !> Aerosol model
    class(model_t), intent(in) :: this

  end function create_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the aerosol optics grid, discretized in wavenumber space
  function optics_grid( this )

    use aero_grid,                     only : grid_t
    import :: model_t

    !> Optical property wave number grid
    type(grid_t) :: optics_grid
    !> Aerosol model
    class(model_t), intent(in) :: this

  end function optics_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculates optical property data, given an aerosol state and destination
  !! arrays
  subroutine compute_optics( this, state, od, od_ssa, od_asym )

    use aero_array,                    only : array_t
    use aero_state,                    only : state_t
    import :: model_t

    !> Aerosol model
    class(model_t),    intent(inout) :: this
    !> Aerosol state
    class(state_t),    intent(inout) :: state
    !> Aerosol optical depth [m]
    class(array_t),    intent(inout) :: od
    !> Aerosol scattering optical depth [m]
    class(array_t),    intent(inout) :: od_ssa
    !> Aerosol asymmetric scattering optical depth [m]
    class(array_t),    intent(inout) :: od_asym

  end subroutine compute_optics

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Frees resources assocaited with a model pointer
  subroutine model_ptr_finalize( this )

    type(model_ptr), intent(inout) :: this

    if( associated( this%ptr_ ) ) deallocate( this%ptr_ )

  end subroutine model_ptr_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module aero_model
