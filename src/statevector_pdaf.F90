!> Building the state vector
!!
!! This module provides variables & routines for
!! defining the state vector.
!!
!! The module contains three routines
!! - **init_id** - initialize the array `id`
!! - **init_sfields** - initialize the array `sfields`
!! - **setup_statevector** - generic routine controlling the initialization
!!
!! The declarations of **id** and **sfields** as well as the
!! routines **init_id** and **init_sfields** might need to be
!! adapted to a particular modeling case. However, for most
!! parts also the configruation using the namelist is possible.
!!
module statevector_pdaf

   use mod_kind_pdaf
   use nemo_pdaf, only: use_wet_state, numcat
   implicit none
   save

   ! Declare Fortran type holding the definitions for model fields
   type state_field
      integer :: ndims = 0                    !< Number of field dimensions (2 or 3)
      integer :: dim = 0                      !< Dimension of the field
      integer :: off = 0                      !< Offset of field in state vector
      character(len=20)  :: variable = ''     !< Name of field
      character(len=30)  :: name_incr = ''    !< Name of field in increment file
      character(len=20)  :: name_rest_n = ''  !< Name of field in restart file (n-field)
      character(len=30)  :: name_bkg_din = '' !< Name of field in background file for direct initialization
                                              !< NEMO reads only ocean physics variables
                                              !< other variables are still be output for sanity checks
      character(len=256) :: rst_file = ''     !< Name of restart file
      character(len=3)   :: k_name = 'lev'    !< Name of field in vertical level direction
                                              !< 'lev' for vertical levels
                                              !< or 'cat' for sea ice categories
                                              !< luckily, all variables with numcat
                                              !< dimension seem to be 2D.
                                              !< This could be different for other sea ice models.
                                              !< but let's keep it the case for now.
      character(len=20)  :: operation = ''    !< operations to form state vector
                                              !< 'sum_over_cat_ice': summation over sea ice categories
                                              !< 'select_cat_ice': selection of a specific sea ice category
      integer            :: ice_cat = 0       !< sea ice category to select if operation is 'select_cat_ice'
      character(len=20)  :: unit = ''         !< Unit of variable
      integer :: transform = 0                !< Type of variable transformation
      real(pwp) :: trafo_shift = 0.0_pwp      !< Constant to shift value in transformation
      integer :: limit = 0                    !< Whether to limit the value of the variable
                                              !< 0: no limits, 1: lower limit, 2: upper limit, 3: both limits
      real(pwp) :: max_limit = 0.0_pwp        !< Upper limit of variable
      real(pwp) :: min_limit = 0.0_pwp        !< Lower limit of variable
   end type state_field

   ! Declare Fortran type holding the definitions for local model fields
   ! This is separate from state_field to support OpenMP
   type state_field_l
      integer :: dim = 0                    !< Dimension of the field
      integer :: off = 0                    !< Offset of field in state vector
   end type state_field_l

   ! Type variable holding the definitions of model fields
   type(state_field), allocatable :: sfields(:)
   ! Type variable holding the definitions of local model fields
   ! This is separate from sfields to support OpenMP
   type(state_field_l), allocatable :: sfields_l(:)

   integer :: dim_state     !< Global model state dimension
   integer :: dim_state_p   !< Model state dimension for PE-local domain

   ! Variables to handle multiple fields in the state vector
   integer :: n_fields          !< number of fields in state vector

   namelist /sv_nml/ n_fields, use_wet_state, numcat
   namelist /sfields_nml/ sfields

!$OMP THREADPRIVATE(sfields_l)

contains
   !> This initializes the array sfields
   !!
   !! This routine initializes the sfields array with specifications
   !! of the fields in the state vector.
   !!
   subroutine init_sfields()
      use mod_kind_pdaf
      use nemo_pdaf, only: sdim2d, sdim3d
      use parallel_pdaf, only: abort_parallel
      implicit none
      ! *** Local variables ***
      integer :: id_var            ! Index of a variable in state vector
      ! *** Specifications for each model field in state vector ***
      ! Some examples:
      ! SSH
      ! sfields(id_var)%ndims = 2
      ! sfields(id_var)%variable = 'zos'
      ! sfields(id_var)%name_incr = 'bckineta'
      ! sfields(id_var)%name_rest_n = 'sshn'
      ! sfields(id_var)%rst_file = 'restart_in.nc'
      ! sfields(id_var)%unit = 'm'

      ! Temperature
      ! sfields(id_var)%ndims = 3
      ! sfields(id_var)%variable = 'thetao'
      ! sfields(id_var)%name_incr = 'bckint'
      ! sfields(id_var)%name_rest_n = 'tn'
      ! sfields(id_var)%rst_file = 'restart_in.nc'
      ! sfields(id_var)%unit = 'degC'

      ! Salinity
      ! sfields(id_var)%ndims = 3
      ! sfields(id_var)%variable = 'so'
      ! sfields(id_var)%name_incr = 'bckins'
      ! sfields(id_var)%name_rest_n = 'sn'
      ! sfields(id_var)%rst_file = 'restart_in.nc'
      ! sfields(id_var)%unit = 'psu'
      ! sfields(id_var)%transform = 0
      ! sfields(id_var)%trafo_shift = 0.0
      ! sfields(id_var)%limit = 0
      ! sfields(id_var)%min_limit = 0.000001

      ! U-velocity
      ! sfields(id_var)%ndims = 3
      ! sfields(id_var)%variable = 'uo'
      ! sfields(id_var)%name_incr = 'bckinu'
      ! sfields(id_var)%name_rest_n = 'un'
      ! sfields(id_var)%rst_file = 'restart_in.nc'
      ! sfields(id_var)%unit = 'm/s'

      ! V-velocity
      ! sfields(id_var)%ndims = 3
      ! sfields(id_var)%variable = 'vo'
      ! sfields(id_var)%name_incr = 'bckinv'
      ! sfields(id_var)%name_rest_n = 'vn'
      ! sfields(id_var)%rst_file = 'restart_in.nc'
      ! sfields(id_var)%unit = 'm/s'

      open (20,file='namelist_cfg.pdaf')
      read (20,NML=sfields_nml)
      close (20)

      do id_var = 1, n_fields
         if (sfields(id_var)%variable == '') then
            write (*, '(a,i2,a)') 'NEMO-PDAF: variable name for field', id_var, ' is not specified.'
            call abort_parallel()
         end if
         if (sfields(id_var)%name_incr == '') then
            write (*, '(a,i2,a)') 'NEMO-PDAF: name of field in increment file for field', id_var, ' is not specified.'
            call abort_parallel()
         end if
         if (sfields(id_var)%rst_file == '') then
            write (*, '(a,i2,a)') 'NEMO-PDAF: name of restart file for field', id_var, ' is not specified.'
            call abort_parallel()
         end if
         if (sfields(id_var)%ndims == 2) then
            sfields(id_var)%dim = sdim2d
         else if (sfields(id_var)%ndims == 3) then
            sfields(id_var)%dim = sdim3d
         else
            write (*, '(a,i2,a)') 'NEMO-PDAF: cannot handle', &
                                  sfields(id_var)%ndims, ' number of dimensions.'
            call abort_parallel()
         end if
      end do

   end subroutine init_sfields
   ! ===================================================================================

   !> Calculate the dimension of the process-local statevector.
   !!
   !! This routine is generic. case-specific adaptions should only
   !! by done in the routines init_id and init_sfields.
   !!
   subroutine setup_statevector()
      use config_pdaf, only: screen
      use mod_kind_pdaf
      use parallel_pdaf, &
            only: mype=>mype_ens, comm_ensemble, &
            comm_model, MPI_SUM, MPI_INTEGER, MPIerr
      implicit none
      ! *** Local variables ***
      integer :: i                 ! Counters
      ! ***********************************
      ! *** Initialize the state vector ***
      ! ***********************************
      ! *** Initialize array `sfields` ***
      allocate(sfields(n_fields))
      call init_sfields()
      ! *** Compute offsets ***
      ! Define offsets in state vector
      sfields(1)%off = 0
      do i = 2, n_fields
         sfields(i)%off = sfields(i-1)%off + sfields(i-1)%dim
      end do
      ! *** Set state vector dimension ***
      dim_state_p = sum(sfields(:)%dim)
      ! *** Write information about the state vector ***
      if (mype==0) then
         write (*,'(/a,2x,a)') 'NEMO-PDAF', '*** Setup of state vector ***'
         write (*,'(a,5x,a,i5)') 'NEMO-PDAF', '--- Number of fields in state vector:', n_fields
         write (*,'(a,a4,3x,a2,3x,a8,6x,a5,7x,a3,7x,a6,4x,a6)') &
               'NEMO-PDAF','pe','ID', 'variable', 'ndims', 'dim', 'offset'
      end if

      if (mype==0 .or. screen>2) then
         do i = 1, n_fields
            write (*,'(a, i4, i5,3x,a10,2x,i5,3x,i10,3x,i10,4x,l)') 'NEMO-PDAF', &
                  mype, i, sfields(i)%variable, sfields(i)%ndims, sfields(i)%dim, sfields(i)%off
         end do
      end if


      if (screen>1 .or. mype==0) &
            write (*,'(a,2x,a,1x,i4,2x,a,1x,i10)') &
            'NEMO-PDAF', 'PE', mype, 'PE-local full state dimension: ',dim_state_p

      call MPI_Reduce(dim_state_p, dim_state, 1, MPI_INTEGER, MPI_SUM, 0, COMM_model, MPIerr)
      if (mype==0) &
         write (*,'(a,2x,a,1x,i10)') 'NEMO-PDAF', 'Global state dimension: ',dim_state
      call MPI_Barrier(comm_ensemble, MPIerr)

   end subroutine setup_statevector

end module statevector_pdaf
