! Copyright (C) 2022 National Center for Atmospheric Research,
! National Technology & Engineering Solutions of Sandia, LLC (NTESS),
! and the U.S. Environmental Protection Agency (USEPA)
!
! SPDX-License-Identifier: Apache-2.0
!


module aero_oslo_optics

  use aero_constants,                  only : rk => real_kind
  use aero_grid,                       only : grid_t
  use aero_model,                      only : model_t
  use aero_state,                      only : state_t
  use pmxsub

  implicit none
  private

  public :: aero_oslo_t

  !> Aerosol model parameters
  type, extends(model_t) :: aero_oslo_t
    private
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Put here configuration data for your aerosol package, parameters !!
    !! that do not vary during the course of a simulation.              !!
    !!                                                                  !!
    !! In this simplified example we include the wave number grid for   !!
    !! aerosol optical properties, as well as optical properties for    !!
    !! "mixed-type" aerosols.                                           !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
    !> Optics grid in wave number [m-1]
    type(grid_t) :: grid_
    !> Aerosol optical depth [m]
    real(kind=rk), allocatable :: tau_(:,:,:)
    !> Single scattering albedo [-]
    real(kind=rk), allocatable :: omega_(:,:,:)
    !> Asymmetry parameter [-]
    real(kind=rk), allocatable :: g_(:,:,:)

    

    integer :: &
              nmodes_,   & ! Number of modes
              pver_,     & ! number of vertical levels
              nbands_,   & ! number of bands
              ncol_,     & ! number of columns
              pcols_,    &
              lchnk_,     &


  contains
    procedure :: name => model_name
    procedure :: create_state
    procedure :: optics_grid
    procedure :: compute_optics
  end type aero_oslo_t
  !> Aerosol state specific to this model
  type, extends(state_t) :: aero_oslo_state_t
    private
    real(kind=rk), allocatable :: mixed_type_
    real(kind=rk), allocatable :: od_work_(:)
    real(kind=rk), allocatable :: Nnatk_(:,:,:) ! number concentration
  end type aero_oslo_state_t




  interface aero_oslo_t
    module procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Creates and configure the aerosol model
  function constructor( description_file ) result( model )

    use aero_array,                    only : array_t
    use aero_util,                     only : assert_msg
    use commondefinitions,             only : nmodes, nbmodes                         
#ifdef AERO_USE_NETCDF
    use netcdf,                        only : nf90_open, nf90_close,          &
                                              NF90_NOWRITE, NF90_NOERR
#endif
    type(aero_oslo_t), pointer    :: model
    character(len=*), intent(in) :: description_file

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Set parameters/configuration data for the aerosol package here !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, parameter :: nCol = 1
    ! Initialize the aerosol grid with wavelength data pulled from
    ! https://acp.copernicus.org/articles/18/7815/2018/acp-18-7815-2018-f03.pdf
    real(kind=rk) :: wavelengths(4) = & ! [nm]
      (/ 440.0_rk, 675.0_rk, 870.0_rk, 1020.0_rk /)
    real(kind=rk) :: wave_numbers(4) ! [m-1]
    class(array_t), pointer :: interfaces

    integer :: i, netcdf_file, k
#ifdef AERO_USE_NETCDF
    ! access NetCDF data
    if( len_trim( description_file ) > 0 ) then
      call assert_msg( 724306399,                                             &
          nf90_open( trim( description_file ), NF90_NOWRITE, netcdf_file )    &
            == NF90_NOERR,                                                    &
          "Error opening NetCDF file '"//trim( description_file )//"'" )
      call assert_msg( 299736143, nf90_close( netcdf_file ) == NF90_NOERR,    &
          "Error closing NetCDF file '"//trim( description_file )//"'" )
    end if
#endif
    ! Convert to wave numbers for the grid's interfaces [m-1]
    do i = 1, 4
      wave_numbers(i) = 1.0e-9_rk / wavelengths(5-i)
    end do
    interfaces => array_t( wave_numbers )

    allocate( model )
    model%grid_ = grid_t( interfaces )
    model%nmodes_ = nmodes
    model%pver_ = 1
    model%pcols_ = nCol 
    model%lchnk_ = 1
    model%nbands_ = 14
    ! Load the averaged optical properties from
    ! https://acp.copernicus.org/articles/18/7815/2018/acp-18-7815-2018-f03.pdf
    allocate( model%tau_(  model%pcols_, model%pver_ ,interfaces%size()   ) )
    allocate( model%omega_( model%pcols_, model%pver_ ,interfaces%size()  ) )
    allocate( model%g_(    model%pcols_, model%pver_ ,interfaces%size()  ) )
    ! TODO: The wavelenghts need to match the ones used by pmxsub_light...
    do k=1, model%pver_
      do i=1, model%pcols_
        model%tau_(i,k,:)   = (/ 0.27_rk,  0.35_rk,   0.5_rk, 0.75_rk /)
        model%omega_(i,k,:) = (/ 0.88_rk, 0.895_rk, 0.905_rk, 0.88_rk /)
        model%g_(i,k,:)     = (/  0.3_rk, 0.035_rk, 0.045_rk, 0.09_rk /)
      end do
    end do  
  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the name of the aerosol model/package
  function model_name( this )

    character(len=:), allocatable :: model_name
    class(aero_oslo_t), intent(in) :: this

    model_name = "aero oslo"

  end function model_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns a newly created aerosol state
  function create_state( this ) result( state )

    class(state_t),    pointer    :: state
    class(aero_oslo_t), intent(in) :: this

    allocate( aero_oslo_state_t  :: state )
    select type( state )
    class is( aero_oslo_state_t )

      !! create a working array for use in calculating optical properties
      allocate( state%od_work_( size( this%tau_ ) ) )

      !! Set some intial state (in a real simulation this would evolve over
      !! time)
      state%mixed_type_ = 0.92_rk

    end select

  end function create_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the aerosol optics grid, discretized in wavenumber space
  function optics_grid( this )

    !> Copy of optical property wave number grid
    type(grid_t) :: optics_grid
    !> My aerosol model
    class(aero_oslo_t), intent(in) :: this

    optics_grid = this%grid_%clone( )

  end function optics_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes optical property data, given an aerosol state and destination
  !! arrays
  subroutine compute_optics( this, state, od, od_ssa, od_asym )

    use aero_array,                    only : array_t

    !> My aerosol model
    class(aero_oslo_t), intent(inout) :: this
    !> Aerosol state
    class(state_t),    intent(inout) :: state
    !> Aerosol optical depth [m]
    class(array_t),    intent(inout) :: od
    !> Aerosol scattering optical depth [m]
    class(array_t),    intent(inout) :: od_ssa
    !> Aerosol asymmetric scattering optical depth [m]
    class(array_t),    intent(inout) :: od_asym

    select type( state )
    class is( aero_oslo_state_t )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Calculate optical properties for the current aerosol state here !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! aerosol optical depth
      call pmxsub_light(this%lchnk_,this%pcols_,this%pver_,this%pcols_,state%Nnatk_,this%tau_, this%omega_, this%g_)

    end select

  end subroutine compute_optics

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module aero_oslo_optics
