!> Module holding IO operations for NEMO-PDAF
!!
!! This code bases in wide parts on the implementation
!! by Wibke Duesterhoeft-Wriggers, BSH, Germany for the
!! CMEMS Baltic Monitoring and forecasting center
!!
module io_pdaf
   use mod_kind_pdaf
   use mpi

   implicit none
   save

   integer :: verbose_io=0   ! Set verbosity of IO routines (0,1,2,3)

   ! Control of IO
   logical :: save_ens_sngl=.false.           ! write set of files holding ensemble of selected field
   logical :: do_deflate=.false.              ! Deflate variables in NC files (this seems to fail for parallel nc)
   character(len=3) :: sgldbl_io='sgl'        ! Write PDAF output in single (sgl) or double (dbl) precision

   character(len=256) :: fname_dom               ! Name of domain file
   character(len=256) :: path_dom                ! Path for NEMO file holding dimensions
   character(len=256) :: f_basename_rst          ! Name of domain file
   character(len=256) :: path_rst                ! Path for NEMO file holding dimensions

   namelist /io_nml/ verbose_io, sgldbl_io, &
                      path_dom, fname_dom, path_rst, f_basename_rst

contains
   !> Print configuration of IO module
   !!
   SUBROUTINE print_io_configuration()
      implicit none
      ! Print PDAF IO configuration to screen
      write (*, '(a,3x,a)') 'NEMO-PDAF','[io_nml]:'
      write (*, '(a,5x,a,i10)') 'NEMO-PDAF','verbose_io ', verbose_io
      write (*, '(a,5x,a,6x,a)')'NEMO-PDAF','sgldbl_io   ', trim(sgldbl_io)
      write (*, '(a,5x,a,6x,a)')'NEMO-PDAF','path to model domain file ', trim(path_dom)
      write (*, '(a,5x,a,6x,a)')'NEMO-PDAF','model domain filename   ', trim(fname_dom)
      write (*, '(a,5x,a,6x,a)')'NEMO-PDAF','path to restart files    ', trim(path_rst)
      write (*, '(a,5x,a,6x,a)')'NEMO-PDAF','basename of restart files', trim(f_basename_rst)
      write (*, '(a,1x,a)') 'NEMO-PDAF','-- End of PDAF IO configuration overview --'
   END SUBROUTINE print_io_configuration

   !> Read domain local information from restart files
   !!
   SUBROUTINE read_local_domain()
      use mpi
      USE netcdf
      use mod_memcount_pdaf, only: memcount
      use parallel_pdaf, only: mype_model, npes_model, comm_model,MPIerr
      use nemo_pdaf, only: i0, j0, ni_p, nj_p, nav_lat, nav_lon
      IMPLICIT NONE
      ! Local variables
      integer :: w                ! domain index
      character(len=256) :: fname ! file name
      INTEGER :: ncid             ! netCDF file identifier
      INTEGER :: varid            ! variable identifier
      INTEGER :: ierr             ! error status

      integer :: dom_size_local(2)
      integer :: dom_pos_first(2)

      call add_slash(path_rst)
      write(fname, '(a,i4.4)') TRIM(f_basename_rst)//'_', mype_model
      !!----------------------------------------------------------------------
      if (mype_model == 0 .and. verbose_io>0) then
         WRITE(*,*)
         WRITE(*,'(/a,1x,a)') 'read_local_domain : Reading local domain information from file'
         WRITE(*,'(a,1x,a/)') '   Input file: ', TRIM(fname)//'.nc'
      end if
      ! Open the NetCDF file
      call check(nf90_open( trim(path_rst)//trim(fname)//'.nc', NF90_NOWRITE, ncid ))
      ! Read 2D grid variables for local domain
      ! dom_size_local
      call check(nf90_get_att( ncid, NF90_GLOBAL, 'DOMAIN_size_local', dom_size_local ))
      ! dom_pos_first
      call check(nf90_get_att( ncid, NF90_GLOBAL, 'DOMAIN_position_first', dom_pos_first ))
      ! nav_lon
      allocate( nav_lon(ni_p, nj_p) )
      allocate( nav_lat(ni_p, nj_p) )
      call memcount(1, 'r', 2*ni_p*nj_p)
      call check(nf90_inq_varid( ncid, 'nav_lon', varid ))
      call check(nf90_get_var( ncid, varid, nav_lon, [1, 1], [ni_p, nj_p] ))
      ! nav_lat
      call check(nf90_inq_varid( ncid, 'nav_lat', varid ))
      call check(nf90_get_var( ncid, varid, nav_lat, [1, 1], [ni_p, nj_p] ))
      ! Close the NetCDF files
      call check (nf90_close( ncid ))
      i0 = dom_pos_first(1)
      j0 = dom_pos_first(2)
      ni_p = dom_size_local(1)
      nj_p = dom_size_local(2)
      ! Screen output
      if (npes_model>1 .and. verbose_io>0) then
         if (mype_model == 0) then
            write (*,'(/a,3x,a)') 'NEMO-PDAF','Grid decomposition:'
            write (*,'(a, 8x,a,2x,a,a,2x,a,a,1x,a,6(1x,a))') &
               'NEMO-PDAF','rank ', 'istart', '  iend', 'jstart', '  jend', '  idim', '  jdim'
         end if
         call MPI_Barrier(comm_model, MPIerr)
         write (*,'(a,2x, a,i6,1x,2i7,2i7,2i7/)') 'NEMO-PDAF', 'RANK', mype_model, i0, i0+ni_p-1, j0, j0+nj_p-1, ni_p, nj_p
      end if
   END SUBROUTINE read_local_domain

   !> Read global domain information from restart files
   !!
   SUBROUTINE read_global_domain()
      use mpi
      USE netcdf
      use mod_memcount_pdaf, only: memcount
      use parallel_pdaf, only: mype_model, npes_model, comm_model, MPIerr
      use nemo_pdaf, only: jpiglo, jpjglo, jpk, i0, j0, ni_p, nj_p, nk_p, &
                           glamt, glamu, glamv, gphit, gphiu, gphiv, &
                           gdept_1d, tmask
      IMPLICIT NONE
      ! Local variables
      INTEGER :: ncid       ! netCDF file identifier
      INTEGER :: varid      ! variable identifier
      INTEGER :: dimid_x    ! dimension id for x
      INTEGER :: dimid_y    ! dimension id for y
      INTEGER :: dimid_t    ! dimension id for t
      INTEGER :: jpi_loc    ! local i-dimension
      INTEGER :: jpj_loc    ! local j-dimension
      INTEGER :: jpt_loc    ! local t-dimension
      INTEGER :: ierr       ! error status
      INTEGER :: i, j, iktop, ikbot ! counter
      integer, allocatable :: k_top(:, :), k_bot(:, :) ! top and bottom wet levels
      !!----------------------------------------------------------------------
      if (mype_model == 0 .and. verbose_io>0) then
         WRITE(*,*)
         WRITE(*,*) 'NEMO-PDAF', 'read_grid_variables : Reading grid variables from file'
         WRITE(*,*) 'NEMO-PDAF', '   Input file: ', TRIM(fname_dom)
      end if
      call add_slash(path_dom)
      ! Open the NetCDF file
      call check(nf90_open( trim(path_dom)//trim(fname_dom), NF90_NOWRITE, ncid ))
      ! Read scalar variables (stored as scalars or 0D variables)
      ! jpiglo
      call check(nf90_inq_varid( ncid, 'jpiglo', varid ))
      call check(nf90_get_var( ncid, varid, jpiglo ))
      ! jpjglo
      call check(nf90_inq_varid( ncid, 'jpjglo', varid ))
      call check(nf90_get_var( ncid, varid, jpjglo ))
      ! jpk
      call check(nf90_inq_varid( ncid, 'jpkglo', varid ))
      call check(nf90_get_var( ncid, varid, jpk ))
      nk_p = jpk
      ! Allocate arrays with dimensions (time, y, x)
      ALLOCATE( glamt(ni_p, nj_p) )
      ALLOCATE( glamu(ni_p, nj_p) )
      ALLOCATE( glamv(ni_p, nj_p) )
      ALLOCATE( gphit(ni_p, nj_p) )
      ALLOCATE( gphiu(ni_p, nj_p) )
      ALLOCATE( gphiv(ni_p, nj_p) )
      ALLOCATE( gdept_1d(nk_p) )
      call memcount(1, 'r', 6*ni_p*nj_p + nk_p)
      ! Read 3D grid variables
      ! glamt
      call check(nf90_inq_varid( ncid, 'glamt', varid ))
      call check(nf90_get_var( ncid, varid, glamt, [i0, j0, 1], [ni_p, nj_p, 1] ))
      ! glamu
      call check(nf90_inq_varid( ncid, 'glamu', varid ))
      call check(nf90_get_var( ncid, varid, glamu, [i0, j0, 1], [ni_p, nj_p, 1] ))
      ! glamv
      call check(nf90_inq_varid( ncid, 'glamv', varid ))
      call check(nf90_get_var( ncid, varid, glamv, [i0, j0, 1], [ni_p, nj_p, 1] ))
      ! gphit
      call check(nf90_inq_varid( ncid, 'gphit', varid ))
      call check(nf90_get_var( ncid, varid, gphit, [i0, j0, 1], [ni_p, nj_p, 1] ))
      ! gphiu
      call check(nf90_inq_varid( ncid, 'gphiu', varid ))
      call check(nf90_get_var( ncid, varid, gphiu, [i0, j0, 1], [ni_p, nj_p, 1] ))
      call check(nf90_get_var( ncid, varid, gphiu, [i0, j0, 1], [ni_p, nj_p, 1] ))
      ! gphiv
      call check(nf90_inq_varid( ncid, 'gphiv', varid ))
      call check(nf90_get_var( ncid, varid, gphiv, [i0, j0, 1], [ni_p, nj_p, 1] ))
      ! gdept_1d
      call check(nf90_inq_varid( ncid, 'gdept_1d', varid ))
      call check(nf90_get_var( ncid, varid, gdept_1d) )
      ! calculate t_mask
      allocate( k_top(ni_p, nj_p) )
      allocate( k_bot(ni_p, nj_p) )
      allocate( tmask(ni_p, nj_p, nk_p) )
      call memcount(1, 'r', ni_p*nj_p*nk_p )
      ! k_top
      call check(nf90_inq_varid( ncid, 'top_level', varid ))
      call check(nf90_get_var( ncid, varid, k_top, [i0, j0], [ni_p, nj_p] ))
      ! k_bot
      call check(nf90_inq_varid( ncid, 'bottom_level', varid ) )
      call check(nf90_get_var( ncid, varid, k_bot, [i0, j0], [ni_p, nj_p] ))
      ! k_top and k_bot
      tmask(:,:,:) = 0._pwp
      DO j = 1, nj_p
         DO i = 1, ni_p
            iktop = k_top(i,j)
            ikbot = k_bot(i,j)
            IF( iktop /= 0 ) THEN       ! water in the column
               tmask(i, j, iktop:ikbot  ) = 1._pwp
            ENDIF
         END DO
      END DO
      ! deallocate k_top and k_bot
      deallocate( k_top, k_bot )
      ! Close the NetCDF file
      call check (nf90_close( ncid ))
      ! *** Screen output ***
      if (mype_model==0 .and. verbose_io>0) then
         write (*,'(/a,5x,a)') 'NEMO-PDAF', '*** NEMO: grid dimensions ***'
         write(*,'(a,3x,2(6x,a),9x,a)') 'NEMO-PDAF', 'jpiglo','jpjglo','jpk'
         write(*,'(a,3x,3i12)') 'NEMO-PDAF', jpiglo, jpjglo, jpk
         write(*,'(a,5x,a,i12)') 'NEMO-PDAF', 'Dimension of global 3D grid box', jpiglo*jpjglo*jpk
         write(*,'(a,5x,a,i12)') 'NEMO-PDAF', 'Number of global surface points', jpiglo*jpjglo
      end if
      !
      call MPI_Barrier(comm_model, MPIerr)
      if (npes_model>1 .and. verbose_io>1) then
         write(*,'(a,2x,a,1x,i4,1x,a,i12)') &
               'NEMO-PDAF', 'PE', mype_model, 'Dimension of local 3D grid box', ni_p*nj_p*nk_p
         write(*,'(a,2x,a,1x,i4,1x,a,i12)') &
               'NEMO-PDAF', 'PE', mype_model, 'Number of local surface points', ni_p * nj_p
      end if
   END SUBROUTINE read_global_domain
   !============================================================================
   !> Check status of NC operation
   !!
   subroutine check(status)
      use netcdf
      use parallel_pdaf, only: abort_parallel
      ! *** Arguments ***
      ! Reading status
      integer, intent ( in) :: status
      ! end program with error message if status is not nf90_noerr
      if(status /= nf90_noerr) then
         print *, trim(nf90_strerror(status))
         call abort_parallel()
      end if
   end subroutine check
   ! ===========================================================================
   !> Add a trailing slash to a path string
   !!
   !! This routine ensures that a string defining a path
   !! has a trailing slash.
   !!
   subroutine add_slash(path)
      implicit none
      ! *** Arguments ***
      !< String holding the path
      character(len=100) :: path
      ! *** Local variables ***
      integer :: strlength
      ! *** Add trailing slash ***
      strlength = len_trim(path)
      if (path(strlength:strlength) /= '/') then
         path = trim(path) // '/'
      end if
   end subroutine add_slash
   !============================================================================
   !> Convert an integer to a strong of length 4
   !!
   character(len=4) function str(k)
      implicit none
      !< number
      integer, intent(in) :: k
      ! string representation
      write (str, '(i4.4)') k
   end function str
   !============================================================================
   !> Check whether a file exists
   !!
   function file_exists(filename) result(res)
      implicit none
      !< File name
      character(len=*),intent(in) :: filename
      !< Status of file
      logical                     :: res
      ! Check if the file exists
      inquire( file=trim(filename), exist=res )
   end function file_exists

end module io_pdaf
