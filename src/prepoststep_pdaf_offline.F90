!> Pre- and Post-Processing of the PDAF output
!!
!!  - For global filters (e.g. ESTKF), the routine is called
!! before the analysis and after the ensemble transformation.
!!  - For local filters (e.g. LESTKF), the routine is called
!! before and after the loop over all local analysis
!! domains.
!!
!! The routine provides full access to the state
!! ensemble to the user.
!! Thus, user-controlled pre- and poststep
!! operations can be performed here. For example
!! the forecast and the analysis states and ensemble
!! covariance matrix can be analyzed, e.g. by
!! computing the estimated variances.
!! For the offline mode, this routine is the place
!! in which the writing of the analysis ensemble
!! can be performed.
!!
!! If a user considers to perform adjustments to the
!! estimates (e.g. for balances), this routine is
!! the right place for it.
!!
!! **Calling Sequence**
!!
!! Called by: PDAF_init_forecst and PDAF3_assimilate
!!
subroutine prepoststep_pdaf(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
     state_p, Uinv, ens_p, flag)
   use config_pdaf, only: screen
   use mod_memcount_pdaf, only: memcount
   use mod_kind_pdaf
   use io_pdaf, only: write_asmdin_mv, write_asminc_mv
   use parallel_pdaf, only: mype=>mype_filter
   use transforms_pdaf, only: transform_field_mv
   implicit none

   ! *** Arguments ***
   integer, intent(in) :: step           !< Current time step (negative for call after forecast)
   integer, intent(in) :: dim_p          !< PE-local state dimension
   integer, intent(in) :: dim_ens        !< Size of state ensemble
   integer, intent(in) :: dim_ens_p      !< PE-local size of ensemble
   integer, intent(in) :: dim_obs_p      !< Dimension of observation vector
   real, intent(inout) :: state_p(dim_p) !< PE-local forecast/analysis state
   !< (The array 'state_p' is not generally not initialized in the case of SEIK.
   !< It can be used freely here.)
   real, intent(inout) :: Uinv(dim_ens-1, dim_ens-1) !< Inverse of matrix U
   real, intent(inout) :: ens_p(dim_p, dim_ens)      !< PE-local state ensemble
   integer, intent(in) :: flag           !< PDAF status flag

   ! *** local variables ***
   integer :: member              ! counters
   integer :: verbose             ! control verbosity of transform_field_mv
   real, save, allocatable :: ens_f_p(:,:) ! Store forecast ensemble for increment file writing (Ensemble mode)
   logical, save :: first = .true. ! Flag for first call to this routine


   ! **********************
   ! *** output increment files ***
   ! **********************
   if (first) then
      if (mype==0) write (*,'(a, 5x,a)') 'NEMO-PDAF', 'Analyze forecast state ensemble'
      ! For using NEMO's asminc module in offline mode: store forecast
      ! Store full forecast ensemble
      allocate(ens_f_p(dim_p, dim_ens))
      ens_f_p = ens_p
      call memcount(2, 'r', dim_p*dim_ens)
      ! *** Transform fields
      do member = 1 , dim_ens
         if (mype==0 .and. member==1) then
            verbose = screen
         else
            verbose = 0
         end if
         call transform_field_mv(1, ens_p(:,member), 11, verbose)
      end do
      first = .false.
   else
      if (mype==0) write (*,'(a, 5x,a)') 'NEMO-PDAF', 'Analyze assimilated state ensemble'
      ! Ensemble KF - store an ensemble of increment files
      if (mype == 0) write (*,'(a,5x,a)') 'NEMO-PDAF', '--- Write ensemble of increments'
      do member = 1, dim_ens
         ! Store member of analysis ensemble
         state_p = ens_p(:,member)
         call write_asmdin_mv(member, ens_f_p(:,member))
         call write_asminc_mv(member, state_p, ens_f_p(:,member))
      end do
      deallocate(ens_f_p)
   end if
end subroutine prepoststep_pdaf
