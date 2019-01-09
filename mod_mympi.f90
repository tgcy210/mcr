!for mpi variables

!!-----------------------------------!!
! preprocessor
#define _USE_MPI  0
!!-----------------------------------!!

module mod_mympi
  
#if _USE_MPI
   use mpi
#endif

   integer ::  myid=0,totp=1,merr=0

end module  mod_mympi

