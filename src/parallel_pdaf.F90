!> Setup parallelisation
!!
!! This module provides variables for the MPI parallelization
!! to be shared between model-related routines. There are variables
!! that are used in the model even without PDAF, and additional variables
!! that are only used if data assimilaion with PDAF is performed.
!! The initialization of communicators for execution with PDAF is
!! performed in `init_parallel_pdaf`.
!!
module parallel_pdaf
   use mpi
   use mod_kind_pdaf

   implicit none
   save

   ! Basic variables for model state integrations
   integer :: COMM_model         ! MPI communicator for model tasks
   integer :: mype_model         ! Rank in COMM_model
   integer :: npes_model         ! Size of COMM_model

   integer :: COMM_ensemble      ! Communicator for entire ensemble
   integer :: mype_ens           ! Rank in COMM_ensemble
   integer :: npes_ens           ! Size of COMM_ensemble

   ! Additional variables for use with PDAF
   integer :: n_modeltasks = 1   ! Number of parallel model tasks

   integer :: COMM_filter        ! MPI communicator for filter PEs
   integer :: mype_filter        ! Rank in COMM_filter
   integer :: npes_filter        ! Size of COMM_filter

   integer :: COMM_couple        ! MPI communicator for coupling filter and model
   integer :: mype_couple        ! Rank in COMM_couple
   integer :: npes_couple        ! Size in COMM_couple

   integer :: mype_world         ! Rank in MPI_COMM_WORLD
   integer :: npes_world         ! Size in MPI_COMM_WORLD

   logical :: modelpe            ! Whether we are on a PE in a COMM_model
   logical :: filterpe = .true.  ! Whether we are on a PE in a COMM_filter
   integer :: MPIerr             ! Error flag for MPI
   integer :: MPIstatus(MPI_STATUS_SIZE)       ! Status array for MPI

contains

   !-------------------------------------------------------------------------------
   !> Initialize MPI communicators for PDAF
   !!
   subroutine init_parallel_pdaf(screen)
      use PDAF, only: PDAF3_set_parallel
      use timer, only: timeit
      implicit none
      integer, intent(in) :: screen ! Control verbosity of PDAF (see config_pdaf)
      integer :: my_color, color_couple !< Variables for communicator-splitting
      integer :: flag                   !< Flag for PDAF communicator setup

      call timeit(5,'ini')
      call timeit(5,'new')
      call timeit(1,'new')

      ! Online in case of online mode: dim_ens = number of parallel model tasks
      n_modeltasks = 1
      ! ***              COMM_ENSEMBLE                ***
      ! *** Generate communicator for ensemble runs   ***
      ! *** only used to generate model communicators ***
      COMM_ensemble = MPI_COMM_WORLD
      call MPI_Comm_Size(COMM_ensemble, npes_ens, MPIerr)
      call MPI_Comm_Rank(COMM_ensemble, mype_ens, MPIerr)
      ! Initialize communicators for ensemble evaluations
      if (mype_ens == 0) &
         write (*, '(/a, 2x, a)') 'PDAF', 'Initialize communicators for assimilation with PDAF'
      ! ***              COMM_MODEL               ***
      ! *** Generate communicators for model runs ***
      my_color = 1
      call MPI_Comm_split(COMM_ensemble, my_color, mype_ens, COMM_model, MPIerr)
      ! Re-initialize PE information according to model communicator
      call MPI_Comm_Size(COMM_model, npes_model, MPIerr)
      call MPI_Comm_Rank(COMM_model, mype_model, MPIerr)
      if (screen > 1) &
         write (*, *) 'PDAF-MODEL: mype(w)= ', mype_ens, '; model task: ', 1, &
            '; mype(m)= ', mype_model, '; npes(m)= ', npes_model
      ! ***         COMM_FILTER                 ***
      ! *** Generate communicator for filter    ***
      call MPI_Comm_split(COMM_ensemble, my_color, mype_ens, &
                          COMM_filter, MPIerr)
      ! Initialize PE information according to filter communicator
      call MPI_Comm_Size(COMM_filter, npes_filter, MPIerr)
      call MPI_Comm_Rank(COMM_filter, mype_filter, MPIerr)
      ! ***              COMM_COUPLE                 ***
      ! *** Generate communicators for communication ***
      ! *** between model and filter PEs             ***
      color_couple = mype_model + 1
      call MPI_Comm_split(COMM_ensemble, color_couple, mype_ens, &
                          COMM_couple, MPIerr)
      ! Initialize PE information according to coupling communicator
      call MPI_Comm_Size(COMM_couple, npes_couple, MPIerr)
      call MPI_Comm_Rank(COMM_couple, mype_couple, MPIerr)
      if (screen > 0) then
         if (mype_ens == 0) then
            write (*, '(/18x, a)') 'PE configuration:'
            write (*, '(a, 2x, a6, a9, a10, a14, a13, /a, 2x, a5, a9, a7, a7, a7, a7, a7, /a, 2x, a)') &
               'Pconf', 'world', 'filter', 'model', 'couple', 'filterPE', &
               'Pconf', 'rank', 'rank', 'task', 'rank', 'task', 'rank', 'T/F', &
               'Pconf', '----------------------------------------------------------'
         end if
         call MPI_Barrier(COMM_ensemble, MPIerr)
         write (*, '(a, 2x, i4, 4x, i4, 4x, i3, 4x, i3, 4x, i3, 4x, i3, 5x, l3)') &
            'Pconf', mype_ens, mype_filter, 1, mype_model, color_couple, &
            mype_couple, filterpe
         call MPI_Barrier(COMM_ensemble, MPIerr)
         if (mype_ens == 0) write (*, '(/a)') ''
      end if
      ! *****************************************************
      ! *** Set communicator within which PDAF operates.  ***
      ! *****************************************************
      CALL PDAF3_set_parallel(COMM_ensemble, COMM_model, COMM_filter, COMM_couple, &
                              1, n_modeltasks, filterpe, flag)
      call timeit(1,'old')
      call timeit(2,'new')

   end subroutine init_parallel_pdaf

   !-------------------------------------------------------------------------------
   !> Initialize the MPI execution environment.
   !!
   subroutine init_parallel()
      implicit none

      integer :: i

      call MPI_INIT(i);
      call MPI_Comm_Size(MPI_COMM_WORLD,npes_world,i)
      call MPI_Comm_Rank(MPI_COMM_WORLD,mype_world,i)

      ! Initialize model communicator, its size and the process rank
      ! Here the same as for MPI_COMM_WORLD
      comm_model = MPI_COMM_WORLD
      npes_model = npes_world
      mype_model = mype_world

   end subroutine init_parallel

   !-------------------------------------------------------------------------------
   !> Finalize the MPI execution environment.
   !!
   subroutine finalize_parallel()

      implicit none

      call  MPI_Barrier(MPI_COMM_WORLD,MPIerr)
      call  MPI_Finalize(MPIerr)

   end subroutine finalize_parallel

   !-------------------------------------------------------------------------------
   !> Terminate the MPI execution environment.
   !!
   subroutine abort_parallel()

      call MPI_Abort(MPI_COMM_WORLD, 1, MPIerr)

   end subroutine abort_parallel

end module parallel_pdaf
