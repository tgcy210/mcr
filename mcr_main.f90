!!!-----------------------------------------------------------------!!!
!  Monte Carlo simulation of light scattering off an spherical particle
!  
!  Author: Ce Yi
!  Last Updated: 10.2018
!!!-----------------------------------------------------------------!!!

program mcrad
   
   use mod_optic
   use mod_ioset
   use mod_utest
   use mod_mympi
   
   implicit none
   !include 'mpif.h'
   
   integer i,idx
   real :: rval,sval(0:10)
      
   real(R_KD) :: tlen_lim, tlen, r_n, rv1
   logical :: inSphere

   real(R_KD) :: dir_n(3),dir_r(3), dot_ni
   integer rk

#ifdef _USE_MPI   
      call MPI_INIT(merr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, totp, merr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, merr)
      ! write(*,"('myid=', I3)") myid
#endif


   call ReadInp
   c_scatter=0
   c_absorb=0
   c_sca_tot=0
   c_abs_tot=0
   
   !perform unit test
   if (nutdbg .ne. 0) then
      write(disp_msg,"('perform testing routines ... ')") 
      call ShowMsg(1)
      call mcr_test      
   endif

   do i=1, num_p
     call GetInitPnt
     call convert_c2s(pos_xyz, pos_rtp)
    
     rk=0 
     
     r_n=n_frac(1)/n_frac(2)
     dir_n=pos_xyz/r_sph

     if (IsFlec(dir_n, r_n)) then
        call add_scatter(dir_cos,rk)
        cycle
     else
        inSphere=.true.
!        write(*,"('i=',I0,' dir_cos=',f10.4)") i,dir_cos(1)**2+dir_cos(2)**2+dir_cos(3)**2
        
        call RANDOM_NUMBER(rval)
        tlen_lim=-log(rval)/siga(2)
        tlen=0d0
     endif

     do while ( inSphere)
        rk=rk+1
        r_n=n_frac(2)/n_frac(1)        
        !calculate track length
        dot_ni=0d0
        do idx=1, 3
           dot_ni=dot_ni+dir_cos(idx)*pos_xyz(idx)/r_sph
        enddo
        if (dot_ni .gt. 0) then
            write(*,"('ERROR: lost particle i=', I0,' dot=',f10.4)") i,dot_ni
            write(*,"('dir_cos=', 3(f10.4,x) )") dir_cos
            write(*,"('pos_xyz=', 3(f10.4,x) )") pos_xyz
            stop
        endif

        rv1=-2d0*dot_ni*r_sph    
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
        pos_xyz=1d0*r_sph/dsqrt(rv1)*pos_xyz

!ERROR trapping section begins
        if (rv1 .lt. 0.98*r_sph) then
           write(*,"('i=',I0, ' rv1=', f10.4)") i, rv1
           write(*,"('pos_xyz=', 3(f10.4,x) )") pos_xyz
           stop
        endif
!ERROR trapping  section ends

        dir_n=-pos_xyz/r_sph
        inSphere=IsFlec(dir_n, r_n)
        if (.not. inSphere) call add_scatter(dir_cos,rk)
     enddo
     
   enddo
  
   
#ifdef _USE_MPI    
      call MPI_REDUCE(c_scatter,c_sca_tot,n_tr,MPI_INTEGER,MPI_SUM,0, MPI_COMM_WORLD,merr)
      call MPI_REDUCE(c_tracker,c_trk_tot,n_tr*n_rk,MPI_INTEGER,MPI_SUM,0, MPI_COMM_WORLD,merr)
      call MPI_REDUCE(c_absorb,c_abs_tot,1,MPI_INTEGER, MPI_SUM,0, MPI_COMM_WORLD,merr)
#else
      c_sca_tot=c_scatter
      c_trk_tot=c_tracker
      c_abs_tot=c_absorb
#endif
   sval=0.0d0 
   if (myid .eq. 0) then 
      write(*,"('number of absorbed  =',I0  )" ) c_abs_tot
      write(*,"('number of scattered =',I0 )" ) sum(c_sca_tot)
      rv1=-1d0
      do i=1, n_tr
         if (i .ge. n_tr*0.76)  then 
            sval(0)=sval(0)+c_sca_tot(i)
         else 
            sval(10)=sval(10)+c_sca_tot(i)
         endif
         if (i .gt. n_tr/2) then
            sval(2)=sval(2)+c_sca_tot(i)
         else
            sval(1)=sval(1)+c_sca_tot(i)
         endif 
         write(*,"('   scattered to [',f8.3,2x,f8.3,']: ', I0)") rv1, rv1+bin_size,c_sca_tot(i)
         rv1=rv1+bin_size
      enddo
      rval=sval(2)/num_p
      sval(3)=(sval(10)+c_absorb)/sval(0)
      sval(4)=sval(10)/sval(0)
      sval(5)=1.0d0 - sval(1)/sval(2)

      sval(6)=sum(c_trk_tot(:,1))
      sval(7)=sum(c_trk_tot)-sval(6)
!      sval(3)=(sval(7)+c_absorb)/sval(6)
!      sval(4)=sval(7)/sval(6)      
      
      write(*,"('   extiction  eff.: ', ES12.5 )") 2.0
      write(*,"('   scattering eff.: ', ES12.5 )") sval(4)
      !write(*,"('   asym. factor   : ', f8.3 )") sval(5)
      
   
   endif

#ifdef _USE_MPI  
     call MPI_FINALIZE(merr)
#endif

end program mcrad
