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
      ! Interface for prepoststep
      external :: prepoststep_pdaf       ! User supplied pre/poststep routine
      ! Localization of state vector
      external :: init_n_domains_pdaf, & ! Provide number of local analysis domains
                  init_dim_l_pdaf        ! Initialize state dimension for local analysis domain
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
