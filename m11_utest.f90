!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  unit testing module for mod_optic
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mod_utest
   use mod_optic
   use mod_ioset   


   !for 3-mirror test and 1-slab test
   real(R_KD) n_mr(3,3)
   complex(R_KD) E_vec_in(3,2)
   
contains

    subroutine mcr_test

       select case (nutdbg)

          case(1)
             call ut_flec
          case(2)
             call ut_3_mirror
          case(3)
             call ut_1_slab
      end select
      call StopMe
     
    end subroutine mcr_test
   
    subroutine ut_flec
       
       implicit none
       real(R_KD) dir_i_ut(3), dir_r_ut(3),dir_n_ut(3)
       real(R_KD) r_n_ut
       
       dir_n_ut=[0d0,0d0,1d0]
       dir_cos=[0d0,0d0,-1d0]

       
       call GetFlecDir(dir_n_ut)
       
       write(*,"('dir_r_ut=', 3(f8.4,x))") dir_cos

       dir_cos=[dsin(pi/3d0), 0d0, -dcos(pi/3d0)]
       call GetFlecDir(dir_n_ut)
       write(*,"('dir_r_ut=', 3(f8.4,x))") dir_cos
   
       r_n_ut=0.8d0
       dir_cos=[dsin(pi/3d0), 0d0, -dcos(pi/3d0)]
       call GetFracDir(dir_n_ut, r_n_ut)
       write(*,"('dir_r_ut=', 3(f8.4,x))") dir_cos

   
    end subroutine

    subroutine ut_3_mirror

        integer i, j, k
        real(R_KD) rval,r_n
        logical :: is_reflected

        rval=dsqrt(2.0d0)/2d0

        !set-up the mirror normal direction
        n_mr(:,1)=[rval,0d0, -rval]
        n_mr(:,2)=[-rval, rval,0d0]
        n_mr(:,3)=[0d0,-rval, rval]
 
        !incident light
        dir_cos=[0d0,0d0,1d0]
        E_vec_in(:,1)=[1d0,0d0,0d0]  !x_polarized
        E_vec_in(:,2)=[0d0,1d0,0d0]  !y_polarized

        
        r_n=0d0
        
        do i=1, 2  ! two tests: x and y polarized incident lights
           E_vec=E_vec_in(:,i)
           do j=1,3  !3-mirror
              is_reflected=IsFlec(n_mr(:,j),r_n)
              if (.not. is_reflected) then
                 write(disp_msg(1), "('Utest failure: 3-mirror test light not reflected')")
                 call ShowMsg(1)                  
              endif
           enddo
           write(disp_msg(1),"('3-mirror test Part ', I0)") i           
           write(disp_msg(2),"('Incident polarization E=', 3('(',ES12.5,x,ES12.5,')') )") &
                    (realpart(E_vec_in(k,i)), imagpart(E_vec_in(k,i)), k=1,3)     
           write(disp_msg(3),"('Outgoing polarization E=', 3('(',ES12.5,x,ES12.5,')') )") &
                    (realpart(E_vec(k)), imagpart(E_vec(k)), k=1,3)     
           write(disp_msg(4),"('     ')")
           call ShowMsg(4)
        enddo
        
    end subroutine ut_3_mirror

    subroutine ut_1_slab

        integer i, j, k,ei
        real(R_KD) rval,r_n
        integer :: counter(3)=0
        logical :: is_reflected, inSlab

        rval=dsqrt(2.0d0)/2d0

        !set-up the mirror normal direction
        n_mr(:,1)=[0d0,0d0, -1d0]
        n_mr(:,2)=[0d0,0d0,1d0]
        n_mr(:,3)=[0d0,-rval, rval]
 
        !incident light
        E_vec_in(:,1)=[1d0,0d0,0d0]  !x_polarized
        E_vec_in(:,2)=[0d0,1d0,0d0]  !y_polarized

        ei=1
        counter=0
        do i=1, num_p   
           E_vec=E_vec_in(:,ei)
           r_n=n_frac(1)/n_frac(2)
           dir_cos=[0d0,rval ,rval]
           inSlab=.false.
           is_reflected=IsFlec(n_mr(:,1),r_n)
           if (is_reflected) then
               counter(3)=counter(3)+1
           else
               inSlab=.true.
           endif
           k=0
           do while (inSlab)
             j=mod(k,2)+1  !j=1: upper surface; j=2: lower surface
             r_n=n_frac(2)/n_frac(1)
             inSlab=IsFlec(n_mr(:,j), r_n)
             if (.not. inSlab) then
                counter(j)=counter(j)+1
             endif
             k=k+1
           enddo
        enddo
        
        write(disp_msg(1),"('1-slab test: tot_particle=',I0)") num_p           
        write(disp_msg(2),"('Incident polarization E=', 3('(',ES12.5,x,ES12.5,')') )") &
                    (realpart(E_vec_in(k,ei)), imagpart(E_vec_in(k,ei)), k=1,3)     
        write(disp_msg(3),"('  Reflected   : ', I0 )") counter(3)
        write(disp_msg(4),"('  backscatted : ', I0 )") counter(2)
        
        write(disp_msg(5),"('  Passed      : ', I0 )") counter(1)
        call ShowMsg(5)
    end subroutine ut_1_slab


end module mod_utest
