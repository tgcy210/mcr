!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  unit testing module for mod_optic
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mod_utest
   use mod_optic
   use mod_ioset   


   !for 3-mirror test
   real(R_KD) n_mr(3,3)
   complex(R_KD) E_vec_in(3,2)
contains

    subroutine mcr_test

       select case (nutdbg)

          case(1)
             call ut_flec
          case(2)
             call ut_3_mirror
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


end module mod_utest
