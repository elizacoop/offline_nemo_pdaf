!> Initialize PDAF
!!
!! This module contains the initialization routine for PDAF
!! `init_pdaf`. Here the ensemble is initialized and distributed
!! and the statevector and state variable information are computed.
!!
!! Contributors:
!! The coupling and PDAF user codes for NEMO contains code parts
!! from different contributors. Namely
!! - Wibke Duesterhoeft-Wriggers, BSH, Germany
!! - Nicholas Byrne, NCEO and University of Reading, UK
!! - Yumeng Chen, NCEO and University of Reading, UK
!! - Yuchen Sun, AWI, Germany
!! - Lars Nerger, AWI, Germany
!!
module initialize_pdaf
   implicit none

contains

  !> This routine collects the initialization of variables for PDAF.
  !!
  !! The initialization routine `PDAF_init` is called
  !! such that the internal initialization of PDAF is performed.
  !! The initialization is used to set-up local domain and filter options
  !! such as the filter type, inflation, and localization radius.
  !!
  !! This variant is for the offline mode of PDAF.
  !!
  !! Much of the initialisation is read from a PDAF-specific namelist.
  !! This is performed in `read_config_pdaf`.
  !!
   subroutine init_pdaf()

      use mod_kind_pdaf
      use PDAF, only: PDAF_init, PDAF_set_iparam, PDAF_set_rparam
      use parallel_pdaf, only: n_modeltasks, task_id, COMM_model, COMM_filter, &
                               COMM_couple, mype_ens, filterpe, abort_parallel
      use config_pdaf, only: screen, step_null, filtertype, subtype, dim_ens, &
                             type_forget, forget, type_trans, type_sqrt, locweight
      use nemo_pdaf, only: set_nemo_grid
      use statevector_pdaf, only: dim_state, dim_state_p, setup_statevector
      use utils_pdaf, only: init_info_pdaf
      use timer, only: timeit, time_temp

      implicit none

      ! *** Local variables ***
      integer :: filter_param_i(2) ! Integer parameter array for filter
      real    :: filter_param_r(1) ! Real parameter array for filter
      integer :: status_pdaf       ! PDAF status flag

      ! *** External subroutines ***
      external :: init_ens_pdaf    ! Ensemble initialization

      ! ***************************
      ! ***   Initialize PDAF   ***
      ! ***************************

      call timeit(2,'old')
      call timeit(3,'new')

      if (mype_ens == 0) then
         write (*, '(/a,1x,a)') 'NEMO-PDAF', 'INITIALIZE PDAF'
      end if

      ! ******************************************
      ! *** Namelist reading and screen output ***
      ! ******************************************
      ! Screen output for PDAF parameters
      if (mype_ens == 0) call init_info_pdaf()
      ! ************************************************
      ! *** Specify state vector and state dimension ***
      ! ************************************************
      ! Initialize dimension information for NEMO grid
      call set_nemo_grid(screen)
      ! Setup state vector
      call setup_statevector(screen)
      ! *****************************************************
      ! *** Call PDAF initialization routine on all PEs.  ***
      ! *****************************************************
      ! *** All filters except LKNETF/EnKF/LEnKF ***
      filter_param_i(1) = dim_state_p ! State dimension
      filter_param_i(2) = dim_ens     ! Size of ensemble
      filter_param_r(1) = forget      ! Forgetting factor
      call PDAF_init(filtertype, subtype, step_null, &
            filter_param_i, 2, &
            filter_param_r, 1, &
            COMM_model, COMM_filter, COMM_couple, &
            task_id, 1, filterpe, init_ens_pdaf, &
            screen, status_pdaf)
      call PDAF_set_iparam(5, type_forget, status_pdaf) ! Type of forgetting factor
      call PDAF_set_iparam(6, type_trans, status_pdaf)  ! Type of ensemble transformation
      call PDAF_set_iparam(7, type_sqrt, status_pdaf)   ! Type of transform square-root (SEIK-sub2/ESTKF)

    ! *** Check whether initialization of PDAF was successful ***
    if (status_pdaf /= 0) then
       write (*, '(/1x,a6,i3,a43,i4,a1/)') &
            'ERROR ', status_pdaf, &
            ' in initialization of PDAF - stopping! (PE ', mype_ens, ')'
       call abort_parallel()
    end if

    call timeit(3,'old')
    call timeit(4,'new')

  end subroutine init_pdaf

end module initialize_pdaf
