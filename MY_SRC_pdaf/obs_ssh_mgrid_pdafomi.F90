!> PDAF-OMI observation module for ssh observations (on model grid)
!!
!! Observation type: SSH on model grid
!!
!! The subroutines in this module are for the particular handling of
!! ssh observations available on the model grid. The observation module
!! also allows to perform a twin experiment in which model output
!! is read and used as observations after adding noise to the values.
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
   use parallel_pdaf, &
      only: mype_filter, abort_parallel
   use PDAF, &
      only: obs_f, obs_l
   use nemo_pdaf, &
        only: i0, j0, istart, jstart
   use netcdf

   implicit none
   save

   !> Whether to assimilate this data type
   logical :: assim_ssh_mgrid = .false.
   !> Observation error standard deviation (for constant errors)
   real(pwp) :: rms_ssh_mgrid = 0.1 !1
   !> Localization cut-off radius
   real(pwp) :: lradius_ssh_mgrid = 1.0
   !> Support radius for weight function
   real(pwp) :: sradius_ssh_mgrid = 1.0
   !> Whether to perform an identical twin experiment
   logical :: twin_exp_ssh_mgrid = .false.
   !> Standard deviation for Gaussian noise in twin experiment
   real(pwp) :: noise_amp_ssh_mgrid = 1
   ! NetCDF file holding observations
   character(lc) :: file_ssh_mgrid = 'my_nemo_ssh_file.nc'
   ! Name of SSH variable in the observation file
   character(lc) :: varname_ssh_mgrid = 'zos'

   !> Instance of full observation data type - see `PDAFomi` for details.
   type(obs_f), target, public :: thisobs
   !> Instance of local observation data type - see `PDAFomi` for details.
   type(obs_l), target, public :: thisobs_l

!$OMP THREADPRIVATE(thisobs_l)

contains

   !>Initialize information on the observation
   !!
   !! The routine is called by each filter process.
   !! at the beginning of the analysis step before
   !! the loop through all local analysis domains.
   !!
   !! It has to count the number of process-local and full
   !! observations, initialize the vector of observations
   !! and their inverse variances, initialize the coordinate
   !! array and index array for indices of observed elements
   !! of the state vector.
   !!
   !! The following four variables have to be initialized in this routine:
   !!
   !! - **thisobs%doassim** - Whether to assimilate ssh
   !! - **thisobs%disttype** - type of distance computation for localization
   !! with ssh
   !! - **thisobs%ncoord** - number of coordinates used for distance
   !! computation
   !! - **thisobs%id_obs_p** - index of module-type observation in PE-local state
   !! vector
   !!
   !!
   !! Optional is the use of:
   !! - **thisobs%icoeff_p**       - Interpolation coefficients for obs. operator (only if interpolation is used)
   !! - **thisobs%domainsize**     - Size of domain for periodicity for *disttype=1* (<0 for no periodicity)
   !! - **thisobs%obs_err_type**   - Type of observation errors for particle filter and NETF (default: 0=Gaussian)
   !! - **thisobs%use_global_obs** - Whether to use global observations or restrict the observations
   !!                               to the relevant ones (default: *.true.* i.e use global full observations)
   !!
   !! Further variables are set when the routine PDAFomi_gather_obs is called.
   !!
   subroutine init_dim_obs_ssh_mgrid(step, dim_obs)

      use PDAF, &
         only: PDAFomi_gather_obs
      use assimilation_pdaf, &
         only: filtertype, delt_obs, use_global_obs
      use statevector_pdaf, &
           only: id, sfields
      use io_pdaf, &
           only: check
      use nemo_pdaf, &
         only: ni_p, nj_p, jpiglo, jpjglo, glamt, gphit, ndastp

      integer, intent(in)    :: step    !< Current time step
      integer, intent(inout) :: dim_obs !< Dimension of full observation vector

      integer :: i, j                          !> Counters
      integer :: cnt_p, cnt0_p                 !> Counters
      integer :: ncid_in                       !> ID for NetCDF file
      integer :: nc_step = 0                   !> Step for observations in NetCDF file
      integer :: id_var                        !> IDs for fields
      integer :: pos(3), cnt(3)                !> NetCDF position arrays for 3D field
      integer :: dim_obs_p                     !> Number of process-local observations
      integer :: i_obs, j_obs                  !> Global gridbox coordinates of observations
      real(pwp), allocatable :: obs(:, :, :)   !> Global observation field
      real(pwp), allocatable :: obs_p(:)       !> PE-local observation vector
      real(pwp), allocatable :: ivar_obs_p(:)  !> PE-local inverse observation error variance
      real(pwp), allocatable :: ocoord_p(:, :) !> PE-local observation coordinates
      real(pwp) :: rad_conv = 3.141592653589793/180.0 !> Degree to radian conversion


      ! *****************************
      ! *** Global setting config ***
      ! *****************************

      if (mype_filter == 0) &
         write (*, '(a,4x,a)') 'NEMO-PDAF', 'Assimilate observations - obs_ssh_mgrid'

      ! Store whether to assimilate this observation type (used in routines
      ! below)
      if (assim_ssh_mgrid) thisobs%doassim = 1

      ! Specify type of distance computation
      thisobs%disttype = 3   ! 3=Haversine

      ! Number of coordinates used for distance computation.
      ! The distance compution starts from the first row
      thisobs%ncoord = 2

      ! SEt to use limited full observations
      thisobs%use_global_obs = use_global_obs

      ! **********************************
      ! *** Read PE-local observations ***
      ! **********************************

      ! Format of ndastp is YYYYMMDD
      if (mype_filter == 0) then
         write (*, '(a, 4x, a, 1x, i8)') &
              'NEMO-PDAF', '--- obs_ssh_mgrid current date:', ndastp
         write (*, '(a,4x,a,a)') 'NEMO-PDAF', '--- name of SSH file variable: ', trim(varname_ssh_mgrid)
      end if

      call check( nf90_open(file_ssh_mgrid, NF90_NOWRITE, ncid_in) )

      call check( nf90_inq_varid(ncid_in, trim(varname_ssh_mgrid), id_var) )

      allocate (obs(jpiglo, jpjglo, 1))
      ! Increment time in NetCDF file so correct obs read
      nc_step = nc_step + delt_obs

      nc_step = 40
if (mype_filter == 0) write (*,*) 'NEMO-PDAF:    Warning: reading step ', nc_step, 'is hard-coded'

      pos = (/1, 1, nc_step/)
      cnt = (/jpiglo, jpjglo, 1/)

      call check( nf90_get_var(ncid_in, id_var, obs, start=pos, count=cnt) )

      call check( nf90_close(ncid_in) )


      ! ***********************************************************
      ! *** Count available observations for the process domain ***
      ! *** and initialize index and coordinate arrays.         ***
      ! ***********************************************************

      cnt_p = 0

      do j = 1, nj_p
         do i = 1, ni_p
            ! Convert to global coordinates
            cnt_p = cnt_p + 1
         end do
      end do

      ! Set number of local observations
      dim_obs_p = cnt_p

      if (cnt_p == 0) write (*, '(/9x, a, i3, 3x, a, i4)') &
         'WARNING: No ssh_mgrid observations on PE:', mype_filter, &
         'NetCDF file step=', nc_step

      obs_nonzero: if (dim_obs_p > 0) then
         ! Vector of observations on the process sub-domain
         allocate (obs_p(dim_obs_p))
         ! Coordinate array of observations on the process sub-domain
         allocate (ocoord_p(2, dim_obs_p))
         ! Coordinate array for observation operator
         allocate (thisobs%id_obs_p(1, dim_obs_p))
         allocate (ivar_obs_p(dim_obs_p))

         cnt_p = 0
         cnt0_p = 0

         do j = 1, nj_p
            do i = 1, ni_p
               ! State vector index counter for observation operator.
               cnt0_p = cnt0_p + 1

               ! Convert to global coordinates.
               i_obs = istart + i - 1
               j_obs = jstart + j - 1

               cnt_p = cnt_p + 1
               obs_p(cnt_p) = obs(i_obs, j_obs, 1)

               ! Observation coordinates - must be in radians for PDAFOMI
               ocoord_p(1, cnt_p) = glamt(i + i0, j + j0)*rad_conv
               ocoord_p(2, cnt_p) = gphit(i + i0, j + j0)*rad_conv

               ! Coordinates for observation operator (gridpoint)
               thisobs%id_obs_p(1, cnt_p) = cnt0_p + sfields(id%ssh)%off
            end do
         end do
      else
         ! No observations on PE, create dummy arrays to pass to PDAFOMI
         allocate (obs_p(1))
         allocate (ivar_obs_p(1))
         allocate (ocoord_p(2, 1))
         allocate (thisobs%id_obs_p(1, 1))
         obs_p = -999999.0
         ivar_obs_p = epsilon(ivar_obs_p)
         ocoord_p = 0
         thisobs%id_obs_p = 1
      end if obs_nonzero

      ! ****************************************************************
      ! *** Define observation errors for process-local observations ***
      ! ****************************************************************

      ! Set inverse observation error variances
      ivar_obs_p(:) = 1.0/(rms_ssh_mgrid*rms_ssh_mgrid)

      ! *********************************************************
      ! *** For twin experiment: Read synthetic observations  ***
      ! *********************************************************

      if (twin_exp_ssh_mgrid) then
         if (dim_obs_p > 0) call add_noise(dim_obs_p, obs_p)
      end if

      ! ****************************************
      ! *** Gather global observation arrays ***
      ! ****************************************

      call PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
                              thisobs%ncoord, lradius_ssh_mgrid, dim_obs)

      ! ********************
      ! *** Finishing up ***
      ! ********************

      ! Deallocate all local arrays
      deallocate (obs)
      deallocate (obs_p, ocoord_p, ivar_obs_p)

      ! Arrays in THISOBS have to be deallocated after the analysis step
      ! by a call to deallocate_obs() in prepoststep_pdaf.

   end subroutine init_dim_obs_ssh_mgrid

   !>###Implementation of observation operator
   !>
   !>This routine applies the full observation operator
   !>for the ssh observations.
   !>
   !>The routine is called by all filter processes.
   !>
   subroutine obs_op_ssh_mgrid(dim_p, dim_obs, state_p, ostate)

      use PDAF, &
         only: PDAFomi_obs_op_gridpoint

      !> PE-local state dimension
      integer, intent(in) :: dim_p
      !> Dimension of full observed state (all observed fields)
      integer, intent(in) :: dim_obs
      !> PE-local model state
      real(pwp), intent(in) :: state_p(dim_p)
      !> Full observed state
      real(pwp), intent(inout) :: ostate(dim_obs)

      ! ******************************************************
      ! *** Apply observation operator H on a state vector ***
      ! ******************************************************

      if (thisobs%doassim == 1) then
         call PDAFomi_obs_op_gridpoint(thisobs, state_p, ostate)
      end if

   end subroutine obs_op_ssh_mgrid

   !>###Initialize local information on the module-type observation
   !>
   !>The routine is called during the loop over all local
   !>analysis domains. It has to initialize the information
   !>about local ssh observations.
   !>
   !>This routine calls the routine `PDAFomi_init_dim_obs_l`
   !>for each observation type. The call allows to specify a
   !>different localization radius and localization functions
   !>for each observation type and local analysis domain.
   !>
   subroutine init_dim_obs_l_ssh_mgrid(domain_p, step, dim_obs, dim_obs_l)

      use PDAF, &
         only: PDAFomi_init_dim_obs_l
      use assimilation_pdaf, &
         only: domain_coords, locweight

      !> Index of current local analysis domain
      integer, intent(in)  :: domain_p
      !> Current time step
      integer, intent(in)  :: step
      !> Full dimension of observation vector
      integer, intent(in)  :: dim_obs
      !> Local dimension of observation vector
      integer, intent(out) :: dim_obs_l

      ! **********************************************
      ! *** Initialize local observation dimension ***
      ! **********************************************

      call PDAFomi_init_dim_obs_l(thisobs_l, thisobs, domain_coords, &
           locweight, lradius_ssh_mgrid, sradius_ssh_mgrid, dim_obs_l)

   end subroutine init_dim_obs_l_ssh_mgrid

   !>###Routine to add model error.
   !>
   subroutine add_noise(dim_obs_p, obs)

      !> Number of process-local observations
      integer, intent(in) :: dim_obs_p
      !> Process-local observations
      real, intent(inout) :: obs(dim_obs_p)

      !> Random noise
      real, allocatable :: noise(:)
      !> Seed for random number generator
      integer, save :: iseed(4)
      !> Flag for first call
      logical, save :: firststep = .true.

      ! Seeds taken from PDAF Lorenz96 routine
      if (firststep) then
         write (*, '(9x, a)') '--- Initialize seed for ssh_mgrid noise'
         iseed(1) = 2*220 + 1
         iseed(2) = 2*100 + 5
         iseed(3) = 2*10 + 7
         iseed(4) = 2*30 + 9
         firststep = .false.
      end if

      ! Generate random Gaussian noise
      allocate (noise(dim_obs_p))
      call dlarnv(3, iseed, dim_obs_p, noise)

      obs = obs + (noise_amp_ssh_mgrid*noise)

      deallocate (noise)

   end subroutine add_noise

 end module obs_ssh_mgrid_pdafomi
