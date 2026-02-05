module initialize_offline
   implicit none
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
      use io_pdaf, only: read_local_domain, read_global_domain
      implicit none

      ! *** Initialize model grid information ***
      call read_local_domain()
      call read_global_domain()

   end subroutine initialize

end module initialize_offline
