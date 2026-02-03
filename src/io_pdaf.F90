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
   character(len=256) :: f_basename_rst          ! Name of restart file
   ! the parth to restart files is constructed by the following variables
   ! path_rst_root//ens_prefix[1-N]//path_rst_suffix
   character(len=256) :: path_rst_root           ! Path for NEMO file holding dimensions
   character(len=256) :: ens_prefix = 'ens_'
   character(len=256) :: path_rst_suffix

   namelist /io_nml/ verbose_io, sgldbl_io, path_dom, fname_dom, &
                     path_rst_root, ens_prefix, path_rst_suffix,f_basename_rst

contains
   !> Print configuration of IO module
   !!
   SUBROUTINE print_io_configuration()
      implicit none
      character(len=256) :: path_rst ! path to restart files
      ! Construct restart file path
      call add_slash(path_rst_root)
      write(path_rst, '(a,i0,"/",a)') TRIM(path_rst_root)//TRIM(ens_prefix),1,TRIM(path_rst_suffix)
      call add_slash(path_rst)
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
      character(len=256) :: path_rst ! path to restart files
      INTEGER :: ncid             ! netCDF file identifier
      INTEGER :: varid            ! variable identifier
      INTEGER :: ierr             ! error status

      integer :: dom_size_local(2)
      integer :: dom_pos_first(2)

      ! Construct restart file path
      call add_slash(path_rst_root)
      write(path_rst, '(a,i0,"/",a)') TRIM(path_rst_root)//TRIM(ens_prefix),1,TRIM(path_rst_suffix)
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

   !> Read restart file to form state vector
   !!
   SUBROUTINE read_restart(ens_member, state_p)
      use mpi
      USE netcdf
      use mod_memcount_pdaf, only: memcount
      use parallel_pdaf, only: mype_model, npes_model, comm_model, MPIerr
      use transforms_pdaf, only: field2state
      use nemo_pdaf, only: i0, j0, ni_p, nj_p, nk_p, tmp_4d
      use statevector_pdaf, only: sfields, n_fields
      IMPLICIT NONE
      !*** Arguments ***
      integer, intent(in) :: ens_member !< Ensemble member index
      real(pwp), intent(inout) :: state_p(:) !< State vector
      ! Local variables
      character(len=256) :: fname ! file name
      character(len=256) :: path_rst ! path to restart files
      INTEGER :: ncid             ! netCDF file identifier
      INTEGER :: varid            ! variable identifier
      integer :: i                ! counter
      ! Construct restart file path
      call add_slash(path_rst_root)
      write(path_rst, '(a,i0,"/",a)') TRIM(path_rst_root)//TRIM(ens_prefix), &
                                      ens_member,TRIM(path_rst_suffix)
      call add_slash(path_rst)

      if (verbose_io>0 .and. mype_model==0) &
            write(*,'(a,4x,a)') 'NEMO-PDAF','*** Ensemble: Reading model restart file'

      if (.not. allocated(tmp_4d)) allocate(tmp_4d(ni_p, nj_p, nk_p, 1))

      ! Initialize state
      state_p = 0.0_pwp

      do i = 1, n_fields
         write(fname, '(a,i4.4)') TRIM(sfields(i)%rst_file)//'_', mype_model
         if (verbose_io>1 .and. mype_model==0) then
            write(*,'(a,2x,a)') 'NEMO-PDAF', 'Reading: '//trim(path_rst)//trim(fname)//'.nc'
            write (*,'(a,i5,1x,a,a,a,i10)') &
                  'NEMO-PDAF', i, 'Variable: ',trim(sfields(i)%variable), ',  offset', sfields(i)%off
         end if
         ! Open the file
         call check( nf90_open(trim(path_rst)//trim(fname)//'.nc', nf90_nowrite, ncid) )
         !  Read field
         call check( nf90_inq_varid(ncid, trim(sfields(i)%name_rest_n), varid) )
         if (sfields(i)%ndims == 3) then
            call check( nf90_get_var(ncid, varid, tmp_4d, &
                  start=(/1, 1, 1, 1/), count=(/ni_p, nj_p, nk_p, 1/)) )
         else
            call check( nf90_get_var(ncid, varid, tmp_4d(:,:,1,1), &
                  start=(/1, 1, 1/), count=(/ni_p, nj_p, 1/)) )
         end if
         ! Close the file
         call check( nf90_close(ncid) )
         ! Convert field to state vector
         call field2state(tmp_4d, state_p, sfields(i)%off, sfields(i)%ndims)
      end do

      if (verbose_io>2) then
         do i = 1, n_fields
            write(*,*) 'Min and max for ',trim(sfields(i)%variable),' :     ', &
                  minval(state_p(sfields(i)%off+1:sfields(i)%off+sfields(i)%dim)), &
                  maxval(state_p(sfields(i)%off+1:sfields(i)%off+sfields(i)%dim))
         enddo
      end if
   END SUBROUTINE read_restart

   !===========================================================================
   !> Write a state vector as model fields into a file
   !!
   subroutine write_increment_mv(state, state_f, filename, &
         attime, nsteps, step, transform)
      use netcdf
      use nemo_pdaf, only: tmp_4d, ni_p, nj_p, nk_p
      implicit none
      ! *** Arguments ***
      real(pwp),        intent(inout) :: state(:)   ! Analysis state vector
      real(pwp),        intent(inout) :: state_f(:) ! Forecast state vector
      character(len=*), intent(in) :: filename      ! File name
      real(pwp),        intent(in) :: attime        ! Time attribute
      integer(4),       intent(in) :: transform     ! Whether to transform fields
      ! *** Local variables ***
      integer(4) :: ncid
      integer(4) :: dimids_field(4)
      integer(4) :: i
      integer(4) :: dimid_time, dimid_lvls, dimid_lat, dimid_lon
      integer(4) :: id_dateb, id_datef
      integer(4) :: id_lat, id_lon, id_lev, id_time, id_incr
      integer(4) :: startC(2), countC(2)
      integer(4) :: startt(4), countt(4)
      integer(4) :: startz(1), countz(1)
      integer(4) :: nf_prec      ! Precision for netcdf output of model fields
      real(pwp)  :: fillval
      integer(4) :: verbose      ! Control verbosity
      ! **********************
      ! *** Initialization ***
      ! **********************
      if (verbose_io>0 .and. mype==0) write (*,'(8x,a)') '--- Write increment file'
      ! *** Set increment times ***
      ! Time for direct initialisation in Nemo (time of restart file which is used for adding to increment file)
      timeInIncr(1)=incrTime
      ! Start time of interval on which increment is valid (later for time ramp initialisation of increment)
      if (incrTime>0.0 .and. startIncrTime==0.0) then
         bgnTimeInterv(1)=incrTime
      else
         bgnTimeInterv(1)=startIncrTime
      end if
      ! End time of interval on which increment is valid (later for time ramp initialisation of increment)
      if (incrTime>0.0 .and. endIncrTime==0.0) then
         finTimeInterv(1)=incrTime
      else
         finTimeInterv(1)=endIncrTime
      end if
      ! Prepare file writing
      if (.not. allocated(tmp_4d)) allocate(tmp_4d(ni_p, nj_p, nk_p, 1))
      nf_prec = NF90_DOUBLE
      fillval = 1.0e20_pwp
      ! *****************************
      ! *** Create and write file ***
      ! *****************************
      ! *** Create file ***
      if (verbose_io>0 .and. mype==0) &
            write (*,'(a,1x,a,a)') 'NEMO-PDAF', 'Create file: ', trim(filename)

      if (npes==1) then
         call check( NF90_CREATE(trim(filename),NF90_NETCDF4,ncid))
      else
         call check( NF90_CREATE_PAR(trim(filename), NF90_NETCDF4, comm_filter, MPI_INFO_NULL, ncid))
      end if
      call check( NF90_PUT_ATT(ncid,  NF90_GLOBAL, 'title', &
                               'Increment for NEMO-PDAF data assimilation'))

      ! define dimensions for NEMO-input file
      if (npes==1) then
         call check( NF90_DEF_DIM(ncid,'t', NF90_UNLIMITED, dimid_time))
      else
         call check( NF90_DEF_DIM(ncid,'t', nsteps, dimid_time))
      end if
      call check( NF90_DEF_DIM(ncid, 'z', nk_p, dimid_lvls))
      call check( NF90_DEF_DIM(ncid, 'y', nj_p, dimid_lat) )
      call check( NF90_DEF_DIM(ncid, 'x', ni_p, dimid_lon) )

      dimids_field(4)=dimid_time
      dimids_field(3)=dimid_lvls
      dimids_field(2)=dimid_lat
      dimids_field(1)=dimid_lon

      ! define variables
      call check( NF90_DEF_VAR(ncid, 'time', NF90_DOUBLE, id_time))
      call check( NF90_DEF_VAR(ncid, 'z_inc_dateb', NF90_DOUBLE, id_dateb))
      call check( NF90_DEF_VAR(ncid, 'z_inc_datef', NF90_DOUBLE, id_datef))
      call check( NF90_DEF_VAR(ncid, 'nav_lat', NF90_FLOAT, dimids_field(1:2), id_lat))
      call check( NF90_DEF_VAR(ncid, 'nav_lon', NF90_FLOAT, dimids_field(1:2), id_lon))
      call check( NF90_DEF_VAR(ncid, 'nav_lev', NF90_FLOAT, dimids_field(3), id_lev))
      if (do_deflate) then
         call check( NF90_def_var_deflate(ncid, id_lat, 0, 1, 1) )
         call check( NF90_def_var_deflate(ncid, id_lon, 0, 1, 1) )
         call check( NF90_def_var_deflate(ncid, id_lev, 0, 1, 1) )
      end if

      do i = 1, n_fields
         if (sfields(i)%update) then
            if (sfields(i)%ndims==3) then
               dimids_field(3)=dimid_lvls
               call check( NF90_DEF_VAR(ncid, trim(sfields(i)%name_incr), nf_prec, dimids_field(1:4), id_incr) )
            else
               dimids_field(3)=dimid_time
               call check( NF90_DEF_VAR(ncid, trim(sfields(i)%name_incr), nf_prec, dimids_field(1:3), id_incr) )
            end if
            if (do_deflate) &
                  call check( NF90_def_var_deflate(ncid, id_incr, 0, 1, 1) )

            call check( nf90_put_att(ncid, id_incr, "long_name", trim(sfields(i)%variable)//'_'//trim('Increment')) )
            call check( nf90_put_att(ncid, id_incr, "units", trim(sfields(i)%unit)) )
            call check( nf90_put_att(ncid, id_incr, "coordinates", "nav_lat nav_lon") )
            call check( nf90_put_att(ncid, id_incr, "_FillValue", fillval) )
            call check( nf90_put_att(ncid, id_incr, "missing_value", fillval) )
         end if
      end do

      ! End define mode
      call check( NF90_ENDDEF(ncid) )

      ! write coordinates
      startz(1)=1
      countz(1)=nlvls

      startC(1) = istart
      countC(1) = ni_p
      startC(2) = jstart
      countC(2) = nj_p

      call check( nf90_put_var(ncid, id_lon, lons, startC, countC))
      call check( nf90_put_var(ncid, id_lat, lats, startC, countC))

      if (mype==0) then
         call check( nf90_put_var(ncid,id_lev,depths,startz,countz))

         call check( nf90_put_var(ncid, id_time, timeInIncr, start=(/1/), count=(/1/)))
         call check( nf90_put_var(ncid, id_dateb, bgnTimeInterv, start=(/1/), count=(/1/)))
         call check( nf90_put_var(ncid, id_datef, finTimeInterv, start=(/1/), count=(/1/)))
      end if


      ! *** Write fields

      call check( nf90_inq_varid(ncid, 'time', id_time) )
   !    call check( nf90_VAR_PAR_ACCESS(NCID, id_time, NF90_COLLECTIVE) )

      ! Backwards transformation of state fields
      if (mype==0) then
         verbose = 1
      else
         verbose = 0
      end if
      if (transform==1) then
         call transform_field_mv(2, state, 0, verbose)
         call transform_field_mv(2, state_f, 0, verbose)
      end if

      ! Compute increment
      state = state - state_f

      ! Write each updated field
      do i = 1, n_fields

         if (sfields(i)%update) then

            tmp_4d = fillval
            call state2field(state, tmp_4d, sfields(i)%off, sfields(i)%ndims, tmask)

            if (verbose_io>1 .and. mype==0) &
                  write (*,'(a,1x,a,a)') 'NEMO-PDAF', '--- write variable: ', trim(sfields(i)%variable)

            call check( nf90_inq_varid(ncid, trim(sfields(i)%name_incr), id_incr) )
   !       call check( nf90_VAR_PAR_ACCESS(NCID, id_field, NF90_COLLECTIVE) )

            startt(1) = istart
            countt(1) = ni_p
            startt(2) = jstart
            countt(2) = nj_p
            startt(3) = 1
            countt(3) = nlvls
            startt(4) = step
            countt(4) = 1

            if (sfields(i)%ndims==3) then
               startt(3) = 1
               countt(3) = nlvls

               call check( nf90_put_var(ncid, id_incr, tmp_4d, startt, countt))
            else
               startt(3) = step
               countt(3) = 1

               call check( nf90_put_var(ncid, id_incr, tmp_4d, startt(1:3), countt(1:3)))
            end if
         end if
      end do

      ! *** close file with state sequence ***
      call check( NF90_CLOSE(ncid) )

   end subroutine write_increment_mv
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
