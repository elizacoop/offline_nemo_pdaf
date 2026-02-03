!> Ensemble Initialization
!!
!! This routine calls the routines for initializing the ensemble.
!!
!! The routine is called when the filter is initialized in
!! `PDAF_init`.
!!
!! The routine is called by all filter processes and
!! initializes the ensemble for the process-local domain.
!!
!!  **Calling Sequence**
!!
!!  - Called from: `init_pdaf/PDAF_init` (PDAF module)
!!
subroutine init_ens_pdaf(filtertype, dim_p, dim_ens, state_p, Uinv, &
     ens_p, flag)
   use mod_kind_pdaf
   use parallel_pdaf, only: mype_filter
   use io_pdaf, only: read_restart
   use transforms_pdaf, only: transform_field_mv
   implicit none
   ! *** Arguments ***
   integer, intent(in) :: filtertype                     !< Type of filter to initialize
   integer, intent(in) :: dim_p                          !< PE-local state dimension
   integer, intent(in) :: dim_ens                        !< Size of ensemble
   real(pwp), intent(inout) :: state_p(dim_p)            !< PE-local model state
   !< It is not necessary to initialize the array 'state_p' for LETKF/LESTKF.
   !< It is available here only for convenience and can be used freely.
   real(pwp), intent(inout) :: Uinv(dim_ens-1,dim_ens-1) !< Array not referenced for SEIK
   real(pwp), intent(out)   :: ens_p(dim_p, dim_ens)     !< PE-local state ensemble
   integer, intent(inout) :: flag                        !< PDAF status flag

   ! *** Local variables ***
   integer :: i, member              ! Counters
   real(pwp) :: inv_dim_ens          ! Inverse ensemble size
   integer :: verbose                ! Control verbosity

   ! ************************************
   ! *** Generate ensemble from files ***
   ! ************************************
   ! Read ensemble states as model snapshots from separate files

   if (mype_filter==0) write (*,'(a,1x,a)') 'NEMO-PDAF', 'Initialize ensemble from list of output files'

   if (mype_filter==0 .and. member==1) then
      verbose = 1
   else
      verbose = 0
   end if

   do member = 1 , dim_ens
      call read_restart(member, ens_p(:, member))
      ! *** Transform fields
      call transform_field_mv(1, ens_p(:,member), 11, verbose)
   end do

   inv_dim_ens = 1.0_pwp/real(dim_ens, kind=pwp)
   ! Scale ensemble perturbations - either all using 'ensscale' or field-wise
   do i = 1, dim_p
!$OMP PARALLEL DO private(member)
      do member = 1, dim_ens
         state_p(i) = state_p(i) + inv_dim_ens*ens_p(i, member)
      end do
!$OMP END PARALLEL DO
   end do

end subroutine init_ens_pdaf
