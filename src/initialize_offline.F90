module initialize_offline
   implicit none

   character(len=20)  :: vname_mask              ! Name of variable to read to determine mask
   character(len=256) :: fname_dom               ! Name of domain file
   character(len=256) :: path_dom                ! Path for NEMO file holding dimensions
   character(len=256) :: f_basename_rst          ! Name of domain file
   character(len=256) :: path_rst                ! Path for NEMO file holding dimensions
   integer :: screen=1                           ! Verbosity flag

contains

   !> Initialize parameters for PDAF offline implementation
   !!
   !! This routine reads in the pdaf offline namelist to
   !! initialize parameters for the offline implementation.
   !! The routine afterwards calls the routine that initializes
   !! the model grid information.
   !!
   subroutine initialize
      use mod_kind_pdaf
      use parallel_pdaf, only: mype_ens
      use assimilation_pdaf, only: step_null
      use nemo_pdaf, only: path_dims, file_dims, jptra, ndastp, use_wet_state, type_limcoords
#if defined key_top
      use nemo_pdaf, only: sn_tracer
#endif
      use io_pdaf, only: verbose_io, check, add_slash, coupling_nemo, incrTime, startIncrTime, endIncrTime
      implicit none
      ! *** Local variables ***
      integer   :: day, year, month
      real(pwp) :: rdate
      namelist /pdaf_offline/ fname_dom

      ! *** Read namelist file for PDAF-offline ***
      open (500,file='pdaf_offline.nml')
      read (500,NML=pdaf_offline)
      close (500)

      call add_slash(path_dims)


#if defined key_top
      if (jptra>0) allocate(sn_tracer(jptra))
#endif

      ! Print PDAF parameters to screen
      showconf: if (mype_ens == 0) then
         write (*, '(/a,1x,a)') 'NEMO-PDAF','-- Overview of PDAF-offline configuration --'
         write (*, '(a,3x,a)') 'NEMO-PDAF','[pdaf_offline]:'
         write (*, '(a,5x,a,i10)') 'NEMO-PDAF','screen       ', screen
         write (*, '(a,5x,a,i10)') 'NEMO-PDAF','use_wet_state', use_wet_state
#if defined key_top
         write (*, '(a,5x,a,i10)') 'NEMO-PDAF','jptra        ', jptra
#endif
         write (*, '(a,5x,a,6x,a)')'NEMO-PDAF','path_dims    ', trim(path_dims)
         write (*, '(a,5x,a,6x,a)')'NEMO-PDAF','file_dims    ', trim(file_dims)
         write (*, '(a,5x,a,6x,a)')'NEMO-PDAF','varname      ', trim(varname)
         write (*, '(a,5x,a,1x,i10)') 'NEMO-PDAF','ndastp       ', ndastp
         write (*, '(a,5x,a,f12.2)') 'NEMO-PDAF','incrTime       ', incrTime
         write (*, '(a,5x,a,f12.2)') 'NEMO-PDAF','startIncrTime  ', startincrTime
         write (*, '(a,5x,a,f12.2)') 'NEMO-PDAF','endIncrTime    ', endincrTime
         write (*, '(a,1x,a/)') 'NEMO-PDAF','-- End of PDAF-offline configuration overview --'
      end if showconf


      ! *** Initialize model grid information ***
      call read_local_domain()
      call read_global_domain()

      ! *** Set initialization of lim_coords to using min/max of glamt/gphit
      type_limcoords = 2

   end subroutine initialize

   !> Read domain local information from restart files
   !!
   SUBROUTINE read_local_domain()
      USE netcdf
      use parallel_pdaf, only: mype_model, npes_model
      use nemo_pdaf, only: i0, j0, ni_p, nj_p
      use io_pdaf, only: check
      IMPLICIT NONE
      ! Local variables
      integer :: w          ! domain index
      character(len=256) :: fname ! file name
      INTEGER :: ncid       ! netCDF file identifier
      INTEGER :: varid      ! variable identifier
      INTEGER :: ierr       ! error status

      integer :: dom_size_local(2)
      integer :: dom_pos_first(2)

      write(fname, '(a,i4.4)') TRIM(f_basename_rst), mype_model
      !!----------------------------------------------------------------------
      if (mype_model == 0 .and. screen>0) then
         WRITE(*,*)
         WRITE(*,*) 'read_local_domain : Reading local domain information from file'
         WRITE(*,*) '~~~~~~~~~~~~~~~~~~~'
         WRITE(*,*) '   Input file: ', TRIM(fname)//'.nc'
      end if
      ! Open the NetCDF file
      ierr = check(nf90_open( trim(path_rst)//trim(fname)//'.nc', NF90_NOWRITE, ncid ))
      ! Read 2D grid variables for local domain
      ! dom_size_local
      ierr = check(nf90_get_att( ncid, NF90_GLOBAL, 'DOMAIN_size_local', dom_size_local ))
      ! dom_pos_first
      ierr = check(nf90_get_att( ncid, NF90_GLOBAL, 'DOMAIN_position_first', dom_pos_first ))
      ! Close the NetCDF files
      ierr = check (nf90_close( ncid ))
      i0 = dom_pos_first(1)
      j0 = dom_pos_first(2)
      ni = dom_size_local(1)
      nj = dom_size_local(2)
      ! Screen output
      if (npes_model>1 .and. screen>0) then
         write (*,'(/a,3x,a)') 'NEMO-PDAF','Grid decomposition:'
         write (*,'(a, 8x,a,2x,a,a,2x,a,a,1x,a,6(1x,a))') &
               'NEMO-PDAF','rank ', 'istart', '  iend', 'jstart', '  jend', '  idim', '  jdim'
         write (*,'(a,2x, a,i6,1x,2i7,2i7,2i7)') 'NEMO-PDAF', 'RANK', i, i0, i0+ni-1, j0, j0+nj-1, ni, nj
      end if
   END SUBROUTINE read_local_domain

   !> Read global domain information from restart files
   !!
   SUBROUTINE read_global_domain()
      USE netcdf
      use io_pdaf, only: check
      use parallel_pdaf, only: mype_model
      use nemo_pdaf, only: jpiglo, jpjglo, jpk, &
                           i0, j0, &
                           ni_p, nj_p, nk_p, &
                           glamt, glamu, glamv, &
                           gphit, gphiu, gphiv, &
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

      integer, allocatable :: k_top(:, :), k_bot(:, :)
      !!----------------------------------------------------------------------
      if (mype_model == 0 .and. screen>0) then
         WRITE(*,*)
         WRITE(*,*) 'NEMO-PDAF', 'read_grid_variables : Reading grid variables from file'
         WRITE(*,*) 'NEMO-PDAF', '   Input file: ', TRIM(fname_dom)
      end if
      ! Open the NetCDF file
      ierr = check(nf90_open( trim(path_dom)//trim(fname_dom), NF90_NOWRITE, ncid ))
      ! Read scalar variables (stored as scalars or 0D variables)
      ! jpiglo
      ierr = check(nf90_inq_varid( ncid, 'jpiglo', varid ))
      ierr = check(nf90_get_var( ncid, varid, jpiglo ))
      ! jpjglo
      ierr = check(nf90_inq_varid( ncid, 'jpjglo', varid ))
      ierr = check(nf90_get_var( ncid, varid, jpjglo ))
      ! jpk
      ierr = check(nf90_inq_varid( ncid, 'jpkglo', varid ))
      ierr = check(nf90_get_var( ncid, varid, jpk ))
      nk_p = jpk
      ! Screen output
      if (mype_model == 0 .and. screen>0) then
         WRITE(*,*) 'NEMO-PDAF', '   Grid dimensions:'
         WRITE(*,*) 'NEMO-PDAF', '      jpiglo = ', jpiglo
         WRITE(*,*) 'NEMO-PDAF', '      jpjglo = ', jpjglo
         WRITE(*,*) 'NEMO-PDAF', '      jpk = ', jpk
      end if
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
      ierr = check(nf90_inq_varid( ncid, 'glamt', varid ))
      ierr = check(nf90_get_var( ncid, varid, glamt, [i0, j0, 1], [ni_p, nj_p, 1] ))
      ! glamu
      ierr = check(nf90_inq_varid( ncid, 'glamu', varid ))
      ierr = check(nf90_get_var( ncid, varid, glamu, [i0, j0, 1], [ni_p, nj_p, 1] ))
      ! glamv
      ierr = check(nf90_inq_varid( ncid, 'glamv', varid ))
      ierr = check(nf90_get_var( ncid, varid, glamv, [i0, j0, 1], [ni_p, nj_p, 1] ))
      ! gphit
      ierr = check(nf90_inq_varid( ncid, 'gphit', varid ))
      ierr = check(nf90_get_var( ncid, varid, gphit, [i0, j0, 1], [ni_p, nj_p, 1] ))
      ! gphiu
      ierr = check(nf90_inq_varid( ncid, 'gphiu', varid ))
      ierr = check(nf90_get_var( ncid, varid, gphiu, [i0, j0, 1], [ni_p, nj_p, 1] ))
      ierr = nf90_get_var( ncid, varid, gphiu, [i0, j0, 1], [ni_p, nj_p, 1] )
      ! gphiv
      ierr = check(nf90_inq_varid( ncid, 'gphiv', varid ))
      ierr = check(nf90_get_var( ncid, varid, gphiv, [i0, j0, 1], [ni_p, nj_p, 1] ))
      ! gdept_1d
      ierr = check(nf90_inq_varid( ncid, 'deptht_1d', varid ))
      ierr = check(nf90_get_var( ncid, varid, gdept_1d) )
      ! calculate t_mask
      allocate( k_top(ni_p, nj_p) )
      allocate( k_bot(ni_p, nj_p) )
      allocate( tmask(ni_p, nj_p, nk_p) )
      call memcount(1, 'a', ni_p*nj_p*nk_p )
      ! k_top
      ierr = check(nf90_inq_varid( ncid, 'top_level', varid ))
      ierr = check(nf90_get_var( ncid, varid, k_top, [i0, j0], [ni_p, nj_p] ))
      ! k_bot
      ierr = check(nf90_inq_varid( ncid, 'bottom_level', varid ) )
      ierr = check(nf90_get_var( ncid, varid, k_bot, [i0, j0], [ni_p, nj_p] ))
      ! k_top and k_bot
      tmask(:,:,:) = 0._wp
      DO j = 1, nj_p
         DO i = 1, ni_p
            iktop = k_top(i,j)
            ikbot = k_bot(i,j)
            IF( iktop /= 0 ) THEN       ! water in the column
               tmask(i, j, iktop:ikbot  ) = 1._wp
            ENDIF
         END DO
      END DO
      deallocate( k_top, k_bot )
      ! Close the NetCDF file
      ierr = check (nf90_close( ncid ))
      ! *** Screen output ***
      if (mype==0 .and. screen>0) then
         write (*,'(/a,5x,a)') 'NEMO-PDAF', '*** NEMO: grid dimensions ***'
         write(*,'(a,3x,2(6x,a),9x,a)') 'NEMO-PDAF', 'jpiglo','jpjglo','jpk'
         write(*,'(a,3x,3i12)') 'NEMO-PDAF', jpiglo, jpjglo, jpk
         write(*,'(a,5x,a,i12)') 'NEMO-PDAF', 'Dimension of global 3D grid box', jpiglo*jpjglo*jpk
         write(*,'(a,5x,a,i12)') 'NEMO-PDAF', 'Number of global surface points', dim_2d
      end if
      !
      if (npes>1 .and. screen>1) then
         write(*,'(a,2x,a,1x,i4,1x,a,i12)') &
               'NEMO-PDAF', 'PE', mype, 'Dimension of local 3D grid box', ni_p*nj_p*nz
         write(*,'(a,2x,a,1x,i4,1x,a,i12)') &
               'NEMO-PDAF', 'PE', mype, 'Number of local surface points', dim_2d_p
      end if
   END SUBROUTINE read_global_domain

end module initialize_offline
