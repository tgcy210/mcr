!for mpi variables

module mod_mympi
  
#ifdef _USE_MPI
   use mpi
#endif

   integer ::  myid=0,totp=1,merr=0

end module  mod_mympi

