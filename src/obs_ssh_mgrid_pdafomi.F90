!> PDAF-OMI observation module for ssh observations (on model grid)
!!
!! Observation type: SSH on model grid
!!
!! The subroutines in this module are for the particular handling of
!! ssh observations available on the model grid.
!!
!! The routines are called by the different call-back routines of PDAF.
!!
!! The module uses two derived data type (obs_f and obs_l), which contain
!! all information about the full and local observations. Only variables
!! of the type obs_f need to be initialized in this module. The variables
!! in the type obs_l are initialized by the generic routines from `PDAFomi`.
!!
!! Author: Nicholas Byrne, NCEO & University of Reading, UK
!!
module obs_ssh_mgrid_pdafomi
   use mod_kind_pdaf
   use parallel_pdaf, only: mype_filter, abort_parallel
   use PDAF, only: obs_f, obs_l
   use nemo_pdaf, only: i0, j0
   use netcdf
   implicit none
   save

   !> Whether to assimilate this data type
   logical :: assim_ssh_mgrid = .false.
   !> (1) use global obs.; (0) use domain-reduced full obs.
   integer :: use_global_obs = 1
   !> Index in sfields array
   integer :: id_sfields = 1
   !> Observation error standard deviation (for constant errors)
   real(pwp) :: rms = 0.1
   !> Localization cut-off radius
   real(pwp) :: lradius = 1.0
   !> Support radius for weight function
   real(pwp) :: sradius = 1.0
   !> Path to observation file
   character(lc) :: file_ssh_mgrid = 'ssh_obs.nc'
   !> Name of SSH variable in the observation file
   character(lc) :: varname_ssh_mgrid = 'zos'
   !> Time step index to read from the observation file
   integer :: nc_step = 1

   namelist /mgrid_ssh_nml/ assim_ssh_mgrid, use_global_obs, id_sfields, rms, &
                            lradius, sradius, file_ssh_mgrid, varname_ssh_mgrid, &
                            nc_step

   !> Instance of full observation data type - see `PDAFomi` for details.
   type(obs_f), target, public :: thisobs
   !> Instance of local observation data type - see `PDAFomi` for details.
   type(obs_l), target, public :: thisobs_l

!$OMP THREADPRIVATE(thisobs_l)

contains

   !> Print configuration of this observation type to screen
   subroutine print_ssh_mgrid_configuration()
      write (*, '(a,3x,a)')    'NEMO-PDAF', '[obs_ssh_mgrid_nml]:'
      write (*, '(a,5x,a,l)')  'NEMO-PDAF', 'assim_ssh_mgrid  ', assim_ssh_mgrid
      if (assim_ssh_mgrid) then
         write (*, '(a,5x,a,a)')    'NEMO-PDAF', 'file_ssh_mgrid   ', trim(file_ssh_mgrid)
         write (*, '(a,5x,a,a)')    'NEMO-PDAF', 'varname_ssh_mgrid', trim(varname_ssh_mgrid)
         write (*, '(a,5x,a,i0)')   'NEMO-PDAF', 'nc_step          ', nc_step
         write (*, '(a,5x,a,f12.4)') 'NEMO-PDAF', 'rms              ', rms
         write (*, '(a,5x,a,i0)')   'NEMO-PDAF', 'use_global_obs   ', use_global_obs
         write (*, '(a,5x,a,es12.4)') 'NEMO-PDAF', 'lradius          ', lradius
         write (*, '(a,5x,a,es12.4)') 'NEMO-PDAF', 'sradius          ', sradius
      end if
   end subroutine print_ssh_mgrid_configuration

   !> Initialize information on the observation
   !!
   !! The routine is called by each filter process at the beginning of
   !! the analysis step before the loop through all local analysis domains.
   !!
   !! It counts the number of process-local and full observations, initializes
   !! the vector of observations and their inverse variances, the coordinate
   !! array, and the index array for observed elements of the state vector.
   !!
   !! The following variables must be initialized:
   !! - **thisobs%doassim**    - Whether to assimilate
   !! - **thisobs%disttype**   - Type of distance computation for localization
   !! - **thisobs%ncoord**     - Number of coordinates for distance computation
   !! - **thisobs%id_obs_p**   - Index of observation in PE-local state vector
   !!
   subroutine init_dim_obs_ssh_mgrid(step, dim_obs)
      use PDAF, only: PDAFomi_gather_obs
      use statevector_pdaf, only: sfields
      use nemo_pdaf, only: nwet, wet_pts, glamt, gphit, jpiglo, jpjglo,i0,j0
      use io_pdaf, only: check

      integer, intent(in)    :: step    !< Current time step
      integer, intent(inout) :: dim_obs !< Dimension of full observation vector

      integer :: i                              !< Counter
      integer :: dim_obs_p                      !< Number of process-local observations
      integer :: ncid, id_var                   !< NetCDF file and variable IDs
      integer :: pos(3), cnt(3)                 !< NetCDF hyperslab start/count
      real(pwp), allocatable :: obs_global(:,:) !< Global SSH field read from file
      real(pwp), allocatable :: obs_p(:)        !< PE-local observation vector
      real(pwp), allocatable :: ivar_obs_p(:)   !< PE-local inverse obs. error variance
      real(pwp), allocatable :: ocoord_p(:,:)   !< PE-local observation coordinates
      real(pwp), parameter   :: rad_conv = 3.141592653589793_pwp / 180.0_pwp

      ! *****************************
      ! *** Global setting config ***
      ! *****************************
      if (mype_filter == 0) &
         write (*, '(a,4x,a)') 'NEMO-PDAF', 'Assimilate observations - obs_ssh_mgrid'

      if (assim_ssh_mgrid) thisobs%doassim = 1

      thisobs%disttype    = 3            ! 3=Haversine
      thisobs%ncoord      = 2
      thisobs%use_global_obs = use_global_obs
      thisobs%inno_omit   = 1.0e6_pwp   ! Discard obs with innovations above threshold

      ! **********************************
      ! *** Read observations from file ***
      ! **********************************
      if (mype_filter == 0) then
         write (*, '(a,4x,a,a)')  'NEMO-PDAF', '--- SSH obs file: ', trim(file_ssh_mgrid)
         write (*, '(a,4x,a,a)')  'NEMO-PDAF', '--- SSH variable: ', trim(varname_ssh_mgrid)
         write (*, '(a,4x,a,i0)') 'NEMO-PDAF', '--- reading nc_step: ', nc_step
      end if

      allocate(obs_global(jpjglo, jpiglo))

      call check( nf90_open(trim(file_ssh_mgrid), NF90_NOWRITE, ncid) )
      call check( nf90_inq_varid(ncid, trim(varname_ssh_mgrid), id_var) )

      
      pos = (/nc_step, 1, 1/)
      cnt = (/1, jpjglo, jpiglo/)
      call check( nf90_get_var(ncid, id_var, obs_global, start=pos, count=cnt) )
      !write(*,*) 'obs_global min/max:', minval(obs_global), maxval(obs_global)
      !write(*,*) 'obs_global(10,10), (50,50), (100,100):', &
      !       obs_global(10,10), obs_global(50,50), obs_global(100,100)
      
      
      
      !pos = (/nc_step, 1, 1/)
      !cnt = (/jpiglo, jpjglo, 1/)
      !call check( nf90_get_var(ncid, id_var, obs_global, start=pos, count=cnt) )
      call check( nf90_close(ncid) )

      ! ***********************************************************
      ! *** Initialize index and coordinate arrays on wet points ***
      ! ***********************************************************
      dim_obs_p = nwet

      if (dim_obs_p > 0) then
         allocate(obs_p(dim_obs_p))
         allocate(ocoord_p(2, dim_obs_p))
         allocate(thisobs%id_obs_p(1, dim_obs_p))
         allocate(ivar_obs_p(dim_obs_p))

         do i = 1, nwet

            if (i <= 5) then
               write(*,*) 'wet_pts global:', wet_pts(1,i), wet_pts(2,i), &
                 'local:', wet_pts(6,i), wet_pts(7,i), &
                 'obs:', obs_global(wet_pts(2,i), wet_pts(1,i)), &
                 'glamt:', glamt(wet_pts(6,i), wet_pts(7,i))
             end if
            ! Read SSH at this wet point from the global field.
            ! wet_pts(6,:) and wet_pts(7,:) hold the global i/j indices.
            !obs_p(i) = obs_global(wet_pts(6, i), wet_pts(7, i))
            !obs_p(i) = obs_global(wet_pts(7, i), wet_pts(6, i))
            !obs_p(i) = obs_global(wet_pts(7, i) + j0 - 1, wet_pts(6, i) + istart - i0)
            obs_p(i) = obs_global(wet_pts(2, i), wet_pts(1, i)) !!**IS THIS RIGHT??

            ! Observation coordinates in radians (required by PDAFOMI Haversine)
            ocoord_p(1, i) = glamt(wet_pts(6, i), wet_pts(7, i)) * rad_conv
            ocoord_p(2, i) = gphit(wet_pts(6, i), wet_pts(7, i)) * rad_conv

            ! Map to state vector index
            thisobs%id_obs_p(1, i) = i + sfields(id_sfields)%off
         end do
      else
         ! No wet points on this PE - set dummy arrays
         allocate(obs_p(1))
         allocate(ocoord_p(2, 1))
         allocate(thisobs%id_obs_p(1, 1))
         allocate(ivar_obs_p(1))
         obs_p(1)      = thisobs%inno_omit
         ivar_obs_p(1) = epsilon(ivar_obs_p)
         ocoord_p      = 0.0_pwp
         thisobs%id_obs_p(1, 1) = 1
      end if


     !if (mype_filter == 48) then
   write (*, '(a,4x,a)') 'NEMO-PDAF', '--- Sample SSH observations:'
   do i = 1, min(3,dim_obs_p)
      write (*, '(a,4x,a,i5,a,f12.6)') 'NEMO-PDAF', '    obs(', i, ') = ', obs_p(i)
   end do
      !end if


      deallocate(obs_global)

      ! ****************************************************************
      ! *** Define observation errors for process-local observations ***
      ! ****************************************************************
      ivar_obs_p(:) = 1.0_pwp / (rms * rms)

      if (mype_filter == 0) &
         write (*, '(a,4x,a,i7)') 'NEMO-PDAF', '--- number of SSH observations: ', dim_obs_p

      ! ****************************************
      ! *** Gather global observation arrays ***
      ! ****************************************
      call PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
                              thisobs%ncoord, lradius, dim_obs)
      write(*,*) 'dim_obs_p:', dim_obs_p
      ! ********************
      ! *** Finishing up ***
      ! ********************
      deallocate(obs_p, ocoord_p, ivar_obs_p)

   end subroutine init_dim_obs_ssh_mgrid

   !> Apply the observation operator for SSH observations
   !!
   !! This routine applies the full observation operator.
   !! It is called by all filter processes.
   !!
   subroutine obs_op_ssh_mgrid(dim_p, dim_obs, state_p, ostate)
      use PDAF, only: PDAFomi_obs_op_gridpoint

      integer,   intent(in)    :: dim_p           !< PE-local state dimension
      integer,   intent(in)    :: dim_obs          !< Dimension of full observed state
      real(pwp), intent(in)    :: state_p(dim_p)  !< PE-local model state
      real(pwp), intent(inout) :: ostate(dim_obs) !< Full observed state

      if (thisobs%doassim == 1) then
         call PDAFomi_obs_op_gridpoint(thisobs, state_p, ostate)
      end if

   end subroutine obs_op_ssh_mgrid

   !> Initialize local observation information for one local analysis domain
   !!
   !! Called during the loop over all local analysis domains.
   !!
   subroutine init_dim_obs_l_ssh_mgrid(domain_p, step, dim_obs, dim_obs_l)
      use PDAF, only: PDAFomi_init_dim_obs_l
      use assimilation_pdaf, only: domain_coords
      use config_pdaf, only: locweight

      integer, intent(in)  :: domain_p  !< Index of current local analysis domain
      integer, intent(in)  :: step      !< Current time step
      integer, intent(in)  :: dim_obs   !< Full dimension of observation vector
      integer, intent(out) :: dim_obs_l !< Local dimension of observation vector

      call PDAFomi_init_dim_obs_l(thisobs_l, thisobs, domain_coords, &
           locweight, lradius, sradius, dim_obs_l)

   end subroutine init_dim_obs_l_ssh_mgrid

end module obs_ssh_mgrid_pdafomi
