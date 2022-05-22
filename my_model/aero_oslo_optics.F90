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
    real(kind=rk), allocatable :: tau_(:)
    !> Single scattering albedo [-]
    real(kind=rk), allocatable :: omega_(:)
    !> Asymmetry parameter [-]
    real(kind=rk), allocatable :: g_(:)

    

    integer :: &
              nmodes_,   & ! Number of modes
              pver_,     & ! number of vertical levels
              nbands_,   & ! number of bands
              ncol_,     & ! number of columns
              pcols_,    &
              lchnk_    


  contains
    procedure :: name => model_name
    procedure :: create_state
    procedure :: optics_grid
    procedure :: compute_optics
  end type aero_oslo_t
  !> Aerosol state specific to this model
  type, extends(state_t) :: aero_oslo_state_t
    private
    real(kind=rk), allocatable :: number_conc_
    !real(kind=rk), allocatable :: Nnatk_(:,:,:) ! number concentration
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
    real(kind=rk) :: wavelengths(14) = & ! [nm]
    (/ 200.0_rk, 263.0_rk, 345.0_rk ,442.0_rk, 625.0_rk,778.0_rk,& 
    1242.0_rk, 1299.0_rk, 1626.0_rk, 1942.0_rk, 2150.0_rk, &
    2500.0_rk, 3077.0_rk,3846.0_rk /)
    real(kind=rk) :: wave_numbers(14) ! [m-1]
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
    do i = 1, 14
      wave_numbers(i) = 1.0e-9_rk / wavelengths(15-i)
    end do
    interfaces => array_t( wave_numbers )
    !print *, interfaces
    allocate( model )
    model%grid_ = grid_t( interfaces )
    model%nmodes_ = nmodes
    model%pver_ = 1
    model%pcols_ = nCol 
    model%lchnk_ = 1
    model%nbands_ = 14
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

      !! Set some intial state (in a real simulation this would evolve over
      !! time)
      state%number_conc_ = 1.0e7_rk
      

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
    use aero_constants,                only : rk => real_kind

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

    ! local variables:
    real(rk) :: per_tau(this%pcols_, 0:this%pver_, this%nbands_) ! aerosol extinction optical depth
    real(rk) :: per_tau_w(this%pcols_, 0:this%pver_, this%nbands_) ! aerosol single scattering albedo * tau
    real(rk) :: per_tau_w_g(this%pcols_, 0:this%pver_, this%nbands_) ! aerosol assymetry parameter * w * tau
    real(rk) :: od_temp( this%nbands_ )
    real(rk) :: Nnatk(this%pcols_, this%pver_, 0:this%nmodes_)
    integer :: i,j,k
    select type( state )
    class is( aero_oslo_state_t )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Calculate optical properties for the current aerosol state here !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do i=1, this%pcols_
        do j=1, this%pver_
          do k=0, this%nmodes_
            if (k.eq.6.or.k.eq.7) then
              Nnatk(i,j,k) = state%number_conc_ 
            else
              Nnatk(i,j,k) = 0.0_rk
            endif

          end do 
        end do
      end do
      ! aerosol optical depth
      print *, 'calling pmxsub'
      call pmxsub_light(this%lchnk_,this%pcols_,this%pver_,this%pcols_,Nnatk,per_tau, per_tau_w, per_tau_w_g)
      !state%od_work_(:)=
      call od%copy_in( per_tau(1,1,:) )
      !state%od_work_(:)=per_tau_w(1,1,:)
      call od_ssa%copy_in( per_tau_w(1,1,:) )
      !state%od_work_(:)=per_tau_w_g(1,1,:)
      call od_asym%copy_in( per_tau_w_g(1,1,:) )

    end select

  end subroutine compute_optics

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module aero_oslo_optics
