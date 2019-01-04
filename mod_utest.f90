!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  unit testing module for mod_optic
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mod_utest
   use mod_optic


contains
   
    subroutine ut_flec
       
       implicit none
       real(R_KD) dir_i_ut(3), dir_r_ut(3),dir_n_ut(3)
       real(R_KD) r_n_ut
       
       dir_n_ut=[0d0,0d0,1d0]
       dir_i_ut=[0d0,0d0,-1d0]
       
       call GetFlecDir(dir_n_ut,dir_i_ut, dir_r_ut)
       
       write(*,"('dir_r_ut=', 3(f8.4,x))") dir_r_ut

       dir_i_ut=[dsin(pi/3d0), 0d0, -dcos(pi/3d0)]
       call GetFlecDir(dir_n_ut,dir_i_ut,dir_r_ut)
       write(*,"('dir_r_ut=', 3(f8.4,x))") dir_r_ut
   
       r_n_ut=0.8d0
       call GetFracDir(dir_n_ut, dir_i_ut,r_n_ut,dir_r_ut)
       write(*,"('dir_r_ut=', 3(f8.4,x))") dir_r_ut

   
    end subroutine


end module mod_utest
