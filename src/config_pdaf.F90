module config_pdaf
   use mod_kind_pdaf
   implicit none
   ! *** Below are the generic variables used for configuring PDAF ***
   ! Settings for time stepping - available as namelist read-in
   !< initial time step of assimilation
   integer :: step_null = 0
   ! General control of PDAF - available as command line options
   integer :: screen = 2       !< Control verbosity of PDAF
                           !< * (0) no outputs
                           !< * (1) progress info
                           !< * (2) add timings
                           !< * (3) debugging output
   integer :: dim_ens      !< Size of ensemble
   integer :: filtertype = 7   !< Select filter algorithm:
                           !<   * SEEK (0), SEIK (1), EnKF (2), LSEIK (3), ETKF (4)
                           !<   LETKF (5), ESTKF (6), LESTKF (7), NETF (9), LNETF (10)
                           !<   PF (12), GENOBS (100), 3DVAR (200)
   integer :: subtype = 0      !< Subtype of filter algorithm
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

   ! Filter settings - available as command line options
   !    ! General
   integer   :: type_forget = 0  !< Type of forgetting factor
   real(pwp) :: forget = 1.0       !< Forgetting factor for filter analysis
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
   integer :: locweight = 0     !< Type of localizing weighting of observations
                     !<   * (0) constant weight of 1
                     !<   * (1) exponentially decreasing with SRANGE
                     !<   * (2) use 5th-order polynomial
                     !<   * (3) regulated localization of R with mean error variance
                     !<   * (4) regulated localization of R with single-point error variance
   !    ! SEIK-subtype4/LSEIK-subtype4/ESTKF/LESTKF
   integer :: type_sqrt = 0     !< Type of the transform matrix square-root
                     !<   * (0) symmetric square root
                     !<   * (1) Cholesky decomposition

   namelist /pdaf_nml/ screen, dim_ens, filtertype, subtype, &
                       type_trans, type_sqrt, type_forget, forget, locweight

contains
   !> Print Assimilation Configuration
   !! This routine prints the assimilation configuration
   !! to the standard output.
   SUBROUTINE print_pdaf_configuration()
      implicit none
      ! *** Print configuration ***
      write (*, '(a,3x,a)') 'NEMO-PDAF','[pdaf_nml]:'
      write (*, '(a,5x,a,i10)') 'NEMO-PDAF','screen       ', screen
      write (*, '(a,5x,a,i10)') 'NEMO-PDAF','filtertype   ', filtertype
      write (*, '(a,5x,a,i10)') 'NEMO-PDAF','subtype      ', subtype
      write (*, '(a,5x,a,i10)') 'NEMO-PDAF','type_trans   ', type_trans
      write (*, '(a,5x,a,i10)') 'NEMO-PDAF','type_sqrt    ', type_sqrt
      write (*, '(a,5x,a,i10)') 'NEMO-PDAF','type_forget  ', type_forget
      write (*, '(a,5x,a,f10.3)') 'NEMO-PDAF','forget       ', forget
      write (*, '(a,5x,a,i10)') 'NEMO-PDAF','locweight    ', locweight
   end SUBROUTINE print_pdaf_configuration
end module config_pdaf