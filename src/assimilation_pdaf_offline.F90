!> Assimilation Parameters
!!
!! This module provides variables needed for the
!! assimilation.
!!
!! See `initialize_pdaf` for where many of these
!! variables are initialised.
!!
module assimilation_pdaf
   use mod_kind_pdaf

   implicit none
   save

   !< Type of coupling between NEMO and PDAF
   !< offline: 'rest', 'incr', 'ieoi'
   character(len=4)   :: coupling_nemo = 'ieoi'
   ! *** Model- and data specific variables ***
   integer :: dim_state     !< Global model state dimension
   integer :: dim_state_p   !< Model state dimension for PE-local domain
   ! Settings for time stepping - available as namelist read-in
   !< initial time step of assimilation
   integer :: step_null = 0
   ! *** Below are the generic variables used for configuring PDAF ***
   ! *** Their values are set in init_PDAF                         ***
   ! Settings for observations - available as command line options
   integer :: delt_obs         !< time step interval between assimilation steps
   !< (1) use global obs.; (0) use domain-reduced full obs.
   integer :: use_global_obs = 1

   ! General control of PDAF - available as command line options
   integer :: screen       !< Control verbosity of PDAF
                           !< * (0) no outputs
                           !< * (1) progress info
                           !< * (2) add timings
                           !< * (3) debugging output
   integer :: dim_ens      !< Size of ensemble
   integer :: filtertype   !< Select filter algorithm:
                           !<   * SEEK (0), SEIK (1), EnKF (2), LSEIK (3), ETKF (4)
                           !<   LETKF (5), ESTKF (6), LESTKF (7), NETF (9), LNETF (10)
                           !<   PF (12), GENOBS (100), 3DVAR (200)
   integer :: subtype      !< Subtype of filter algorithm
                           !<   * SEEK:
                           !<     (0) evolve normalized modes
                           !<     (1) evolve scaled modes with unit U
                           !<     (2) fixed basis (V); variable U matrix
                           !<     (3) fixed covar matrix (V,U kept static)
                           !<   * SEIK:
                           !<     (0) ensemble forecast; new formulation
                           !<     (1) ensemble forecast; old formulation
                           !<     (2) fixed error space basis
                           !<     (3) fixed state covariance matrix
                           !<     (4) SEIK with ensemble transformation
                           !<   * EnKF:
                           !<     (0) analysis for large observation dimension
                           !<     (1) analysis for small observation dimension
                           !<   * LSEIK:
                           !<     (0) ensemble forecast;
                           !<     (2) fixed error space basis
                           !<     (3) fixed state covariance matrix
                           !<     (4) LSEIK with ensemble transformation
                           !<   * ETKF:
                           !<     (0) ETKF using T-matrix like SEIK
                           !<     (1) ETKF following Hunt et al. (2007)
                           !<       There are no fixed basis/covariance cases, as
                           !<       these are equivalent to SEIK subtypes 2/3
                           !<   * LETKF:
                           !<     (0) LETKF using T-matrix like SEIK
                           !<     (1) LETKF following Hunt et al. (2007)
                           !<       There are no fixed basis/covariance cases, as
                           !<       these are equivalent to LSEIK subtypes 2/3
                           !<   * ESTKF:
                           !<     (0) standard ESTKF
                           !<       There are no fixed basis/covariance cases, as
                           !<       these are equivalent to SEIK subtypes 2/3
                           !<   * LESTKF:
                           !<     (0) standard LESTKF
                           !<       There are no fixed basis/covariance cases, as
                           !<       these are equivalent to LSEIK subtypes 2/3
                           !<   * NETF:
                           !<     (0) standard NETF
                           !<   * LNETF:
                           !<     (0) standard LNETF
                           !<   * PF:
                           !<     (0) standard PF
                           !<   * 3D-Var:
                           !<     (0) parameterized 3D-Var
                           !<     (1) 3D Ensemble Var using LESTKF for ensemble update
                           !<     (4) 3D Ensemble Var using ESTKF for ensemble update
                           !<     (6) hybrid 3D-Var using LESTKF for ensemble update
                           !<     (7) hybrid 3D-Var using ESTKF for ensemble update
   integer :: incremental  !< Perform incremental updating in LSEIK
   integer :: dim_lag      !< Number of time instances for smoother

   ! Filter settings - available as command line options
   !    ! General
   integer   :: type_forget  !< Type of forgetting factor
   real(pwp) :: forget       !< Forgetting factor for filter analysis
   integer   :: dim_bias     !< dimension of bias vector

   !    ! ENKF
   integer   :: rank_analysis_enkf  !< Rank to be considered for inversion of HPH

   !    ! SEIK/ETKF/ESTKF/LSEIK/LETKF/LESTKF
   integer :: type_trans    !< Type of ensemble transformation
                              !< * SEIK/LSEIK:
                              !< (0) use deterministic omega
                              !< (1) use random orthonormal omega orthogonal to (1,...,1)^T
                              !< (2) use product of (0) with random orthonormal matrix with
                              !<     eigenvector (1,...,1)^T
                              !< * ETKF/LETKF with subtype=4:
                              !< (0) use deterministic symmetric transformation
                              !< (2) use product of (0) with random orthonormal matrix with
                              !<     eigenvector (1,...,1)^T
                              !< * ESTKF/LESTKF:
                              !< (0) use deterministic omega
                              !< (1) use random orthonormal omega orthogonal to (1,...,1)^T
                              !< (2) use product of (0) with random orthonormal matrix with
                              !<     eigenvector (1,...,1)^T
                              !< * NETF/LNETF:
                              !< (0) use random orthonormal transformation orthogonal to (1,...,1)^T
                              !< (1) use identity transformation

   !    ! LSEIK/LETKF/LESTKF/LNETF
   integer :: locweight     !< Type of localizing weighting of observations
                     !<   * (0) constant weight of 1
                     !<   * (1) exponentially decreasing with SRANGE
                     !<   * (2) use 5th-order polynomial
                     !<   * (3) regulated localization of R with mean error variance
                     !<   * (4) regulated localization of R with single-point error variance
   !    ! SEIK-subtype4/LSEIK-subtype4/ESTKF/LESTKF
   integer :: type_sqrt     !< Type of the transform matrix square-root
                     !<   * (0) symmetric square root
                     !<   * (1) Cholesky decomposition


   !< Indices of local state vector in global vector
   integer, allocatable :: id_lstate_in_pstate(:)
   !< Coordinates of local analysis domain
   real(pwp) :: domain_coords(2)


!$OMP THREADPRIVATE(domain_coords, id_lstate_in_pstate)

contains
   !> Performing the Assimilation Step
   !!
   !! This routine is called to perform the analysis step in
   !! offline mode.
   !!
   subroutine assimilate_pdaf()
      use PDAF, only: PDAF3_assim_offline
      use parallel_pdaf, only: mype_ens, abort_parallel
      implicit none
      ! *** Local variables ***
      integer :: status_pdaf         ! PDAF status flag
      ! *** External subroutines ***
      ! Interface between model and PDAF, and prepoststep
      external :: collect_state_pdaf, &  ! Collect a state vector from model fields
                  prepoststep_pdaf       ! User supplied pre/poststep routine

      ! Localization of state vector
      external :: init_n_domains_pdaf, & ! Provide number of local analysis domains
            init_dim_l_pdaf, &            ! Initialize state dimension for local analysis domain
            g2l_state_pdaf, &             ! Get state on local analysis domain from global state
            l2g_state_pdaf                ! Update global state from state on local analysis domain

      ! Interface to PDAF-OMI for local and global filters
      external :: &
            init_dim_obs_pdafomi, &       ! Get dimension of full obs. vector for PE-local domain
            obs_op_pdafomi, &             ! Obs. operator for full obs. vector for PE-local domain
            init_dim_obs_l_pdafomi        ! Get dimension of obs. vector for local analysis domain

      ! *********************************
      ! *** Call assimilation routine ***
      ! *********************************
      call PDAF3_assim_offline(init_dim_obs_pdafomi, obs_op_pdafomi, &
            init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi, &
            prepoststep_pdaf, status_pdaf)

      ! Check for errors during execution of PDAF
      if (status_pdaf /= 0) then
         write (*, '(/1x,a6,i3,a43,i4,a1/)') &
               'ERROR ', status_pdaf, &
               ' in PDAF3_assimilate - stopping! (PE ', mype_ens, ')'
         call abort_parallel()
      end if

   end subroutine assimilate_pdaf

end module assimilation_pdaf
