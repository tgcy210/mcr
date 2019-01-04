!!!-----------------------------------------------------------------!!!
!  Monte Carlo simulation of light scattering off an spherical particle
!  
!  Author: Ce Yi
!  Last Updated: 10.2018
!!!-----------------------------------------------------------------!!!

program mcrad
   
   use mod_optic
   use mod_utest
   use mpi
   
   implicit none
   !include 'mpif.h'
   
   integer i,idx
   real :: rval
   integer myid,totp,merr
      
   real(R_KD) :: tlen_lim, tlen, r_n, rv1
   logical :: inSphere

   real(R_KD) :: dir_n(3),dir_r(3), dot_ni

!   call ut_flec
!   stop
   
   call MPI_INIT(merr)

   call MPI_COMM_SIZE(MPI_COMM_WORLD, totp, merr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, myid, merr)

   c_scatter=0
   c_absorb=0
   c_sca_tot=0
   c_abs_tot=0
   do i=1, num_p
     call GetInitPnt
     call convert_c2s(pos_xyz, pos_rtp)
     
     
     r_n=n_frac(1)/n_frac(2)
     dir_n=pos_xyz

     if (IsFlec(dir_n, r_n)) then
        call add_scatter(dir_cos)
        cycle
     else
        inSphere=.true.
!        write(*,"('i=',I0,' dir_cos=',f10.4)") i,dir_cos(1)**2+dir_cos(2)**2+dir_cos(3)**2
        
        call RANDOM_NUMBER(rval)
        tlen_lim=-log(rval)/siga(2)
        tlen=0d0
     endif

     do while ( inSphere)

        r_n=n_frac(2)/n_frac(1)        
        !calculate track length
        dot_ni=0d0
        do idx=1, 3
           dot_ni=dot_ni+dir_cos(idx)*pos_xyz(idx)
        enddo
        if (dot_ni .gt. 0) then
            write(*,"('ERROR: lost particle i=', I0,' dot=',f10.4)") i,dot_ni
            write(*,"('dir_cos=', 3(f10.4,x) )") dir_cos
            write(*,"('pos_xyz=', 3(f10.4,x) )") pos_xyz
            stop
        endif

        rv1=-2d0*dot_ni    
        tlen=tlen+rv1
        if (tlen .gt. tlen_lim) then
           c_absorb=c_absorb+1
           inSphere=.false.
           cycle
        endif

!DEBUG section begins
!        write(*,"('dir_cos=', 3(f10.4,x) )") dir_cos
!        write(*,"('pos_xyz=', 3(f10.4,x) )") pos_xyz
!        write(*,"('dot_ni=',f10.4)") dot_ni

!        write(*,"('before: dir_cos=',f10.4)") dir_cos(1)**2+dir_cos(2)**2+dir_cos(3)**2
!        write(*,"('before: pos_r=',f10.4)") pos_xyz(1)**2+pos_xyz(2)**2+pos_xyz(3)**2
!        write(*,"('rv1=', f10.4)") rv1
!DEBUG section ends

        pos_xyz=pos_xyz+rv1*dir_cos
        rv1= pos_xyz(1)**2+pos_xyz(2)**2+pos_xyz(3)**2
        pos_xyz=1d0/dsqrt(rv1)*pos_xyz

!ERROR trapping section begins
        if (rv1 .lt. 0.98) then
           write(*,"('i=',I0, ' rv1=', f10.4)") i, rv1
           write(*,"('pos_xyz=', 3(f10.4,x) )") pos_xyz
           stop
        endif
!ERROR trapping  section ends

        dir_n=-pos_xyz
        inSphere=IsFlec(dir_n, r_n)
        if (.not. inSphere) call add_scatter(dir_cos)
     enddo
     
   enddo
  
   
   call MPI_REDUCE(c_scatter,c_sca_tot,n_tr,MPI_INTEGER,MPI_SUM,0, MPI_COMM_WORLD,merr)
   call MPI_REDUCE(c_absorb,c_abs_tot,1,MPI_INTEGER, MPI_SUM,0, MPI_COMM_WORLD,merr)
    
   if (myid .eq. 0) then 
      write(*,"('number of absorbed  =',I0  )" ) c_abs_tot
      write(*,"('number of scattered =',I0 )" ) sum(c_sca_tot)
      rv1=-1d0
      do i=1, n_tr
         write(*,"('   scattered to [',f8.3,2x,f8.3,']: ', I0)") rv1, rv1+bin_size,c_sca_tot(i)
         rv1=rv1+bin_size
      enddo
   endif
   call MPI_FINALIZE(merr)
end program mcrad
