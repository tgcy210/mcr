!!!-----------------------------------------------------------------!!!
!  optical ray tracing in a unit sphere
!  for simulation of light scattering off an spherical particle
!  
!  Ce Yi 10.2018 
!
!!!-----------------------------------------------------------------!!!

module mod_optic

    !constants 
    integer, parameter :: R_KD = selected_real_kind(15,307)  !real*8
    real(R_KD), parameter ::  pi=dacos(-1.0d0)
    
    !pos_xyz: x,y,z coordinates
    !pos_rtp: r,cos(theta),phi, spherical cooridate
    real(R_KD) ::  pos_xyz(3)=0, pos_rtp(3)=0
    
    !light direction: mu_x, mu_y, mu_z (cosines with the x,y,z axies)
    real(R_KD) ::  dir_cos(3)=[0,0,1]   
    
    !material params, assuming incident light in air (m1), sphere made of water (m2)
    real(R_KD) :: n_frac(2)=[1.0d0, 1.50d0]
    !absorption cross section measured by 1/r  
    real(R_KD) :: siga(2)=[0d0, 0.0d0]
    real(R_KD) :: r_size=10.0d0  !radius measured by lambda/2pi
    
    !simulation params
    !number of particles    
    integer :: num_p=100000
    !number of cosine bins for counter
    integer :: n_tr=8
    real(R_KD) ::  bin_size 
    !counters arrays
    integer :: c_absorb=0, c_abs_tot=0
    integer, allocatable :: c_scatter(:), c_sca_tot(:)

    !for debugging
    integer :: n_out=0, npsdbg=12
    integer :: nutdbg=0  !unit test index

    !polarization
    complex(kind=8) :: E_ini(3)=[(1d0,1d0),(1d0,2d0), (0d0,0d0) ]
    complex(kind=8) :: E_vec(3)=[1.0d0,1.0d0,0.0d0]
    real(R_KD) :: k_wave, w_wave
    real(R_KD) :: k_zero, c_zero
contains 
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine GetInitPnt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  uniformly sample a point inside a unit circle (on XY plane)          
!  as the incident light beam initial intersection with the unit sphere
!  assuming incident light travelling along z+
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      real :: rval(3), rad
      
      pos_rtp(1)=1.0d0
      pos_rtp(2)=1.0  !cos(theta)
      call RANDOM_NUMBER(rval)
      pos_rtp(3)=2*pi*rval(1)-pi
      
      rad=rval(2)+rval(3)
      if (rad .gt. 1) rad=2-rad
 
      pos_xyz(1)=rad*dcos(pos_rtp(3))*r_size
      pos_xyz(2)=rad*dsin(pos_rtp(3))*r_size
      pos_xyz(3)=-dsqrt(r_size**2-pos_xyz(1)**2-pos_xyz(2)**2)
      
      dir_cos=[0d0,0d0,1d0]     

      !initialize polarization state
      k_wave =15000*(2*pi)  !cm^-1
      w_wave = 0.0d0 
      E_vec=E_ini

   end subroutine GetInitPnt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine GetFlecDir(d_n)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  calculate a unit vector pointing the reflection direction
!  Input: d_n: surface normal direction, pointing to incident light
!
!  local variables: 
!     d_i: incident light, pointing to surface
!     d_flec: reflected light direction
!  global variable(s):
!     dir_cos: incident dir on input, reflected dir on output
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit none
      real(R_KD), intent(in) :: d_n(3)

      real(R_KD) :: d_i(3), d_flec(3)
      real(R_KD) :: dot_ni
      integer i

      d_i=dir_cos
      dot_ni=0d0
      do i=1,3
         dot_ni=dot_ni+d_n(i)*d_i(i)
      enddo
     
      if (dot_ni .ge. 0) then
         write(*,"('ERROR: dot_ni is greater or equal zero in GetFlecDir:', f10.4 )") dot_ni
         write(*,"('dir_n=', 3(f8.4,x))") d_n
         write(*,"('dir_i=', 3(f8.4,x))") d_i
         stop
      endif

      d_flec=d_i(:)-2*dot_ni*d_n(:)
      dir_cos=d_flec

   end subroutine GetFlecDir

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine GetFracDir(d_n,r_n,IsTir)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  calculate a unit vector pointing the refraction direction
!  Input: d_n(3): surface normal direction, pointing to incident light
!         r_n : ratio of refraction index n1/n2, where n1 is incident medium
!  Output:
!         IsTir: logical, true if total internal reflection,
!                 if not present, assuming false always
!  local variables
!         d_i(3): incident light, pointing to surface
!         d_frac(3): refractioned light direction if IsTir=false
!                    reflectioned direction if IsTir=true
!  global variable(s)
!         dir_cos: incoming direction on input, refracted dir on output
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit none
      real(R_KD), intent(in) :: d_n(3), r_n
      logical,optional,intent(out) :: IsTir

      real(R_KD) :: d_i(3), d_frac(3)
      real(R_KD) :: dot_ni, rv1,rv2
      integer i
      logical :: IsTir_loc=.false.
     
      IsTir_loc=.false. 
      d_i=dir_cos      
      dot_ni=0d0
      do i=1,3
         dot_ni=dot_ni+d_n(i)*d_i(i)
      enddo
     
      if (dot_ni .ge. 0) then
         write(*,"('ERROR: dot_ni is less or equal zero in GetFracDir')")
         stop
      endif

      rv1=dsqrt(1.0d0-dot_ni**2)
      if (rv1*r_n .gt. 1) IsTir_loc=.true.    
     
      if (IsTir_loc) then
         call GetFlecDir(d_n)
         if (present(IsTir)) IsTir=IsTir_loc
         return
      endif
     
      d_frac=r_n*d_i(:)-(r_n*dot_ni+dsqrt(1d0-(r_n*rv1)**2))*d_n(:)
      dir_cos=d_frac

   end subroutine GetFracDir

   logical function IsFlec(d_n,r_n)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  input:
!     d_n: surface normal dir, pointing to the side of incident light 
!     r_n: refaction index ratio n1/n2  
!  local variables:
!     d_i: incident dir, k axis for incident light
!     d_s: s axis in local coordiate (s p k) where k is light moving dir
!          perpendicular to incident plane, 
!          remains unchanged after reflection/refraction
!     d_p: p axis in local coordiate (s p k) for incident light 
!          o the incident plane
!     d_pp: aka p', p axis for ougoing light
!
!  global variable:
!     dir_cos: incident dir on input, aka k
!              outgoing dir on output, aka k'
!  
!  References:                 
!     Fresnel equation:  https://en.wikipedia.org/wiki/Fresnel_equations 
!     Polarization tracking:
!          Garam Yun et al, Three-dimensional polarization ray-tracing 
!          calculus I: definition and diattenuation
!     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit none
      real(R_KD), intent(in) :: d_n(3) ,r_n

      real(R_KD) :: d_i(3), d_s(3), d_p(3), d_pp(3)
      real, parameter :: rflec(2)=[0.1,0.8]
      real rv1, prob_r
      complex (kind=8) :: csp(2),cv1
      integer i
      real(R_KD) :: r_cos, r_sin2, rsp(2), tsp(2),E2_sp(3)      
      complex (kind=8) :: E_sp(2)
      
      d_i=dir_cos
      ! calculate reflection rate for s, p components
      r_cos=0d0
      do i=1, 3
         r_cos=r_cos+d_n(i)*d_i(i)
      enddo
      r_sin2=(1.0d0-r_cos*r_cos)
      
      rv1=1.0-r_sin2*r_n**2
      cv1=zsqrt(complex(rv1,0.0d0))
      
      rv1=-r_n*r_cos
      csp(1)=(rv1-cv1)/(rv1+cv1)
      rsp(1)=zabs(csp(1))**2

      !rv1=-r_cos/r_n
      rv1=-r_cos
      cv1=cv1*r_n
      csp(2)=(cv1-rv1)/(cv1+rv1)
      rsp(2)=zabs(csp(2))**2

      !refraction 
      do i=1,2
         tsp(i)=1.0d0-rsp(i)
      enddo

      !incident light global to local coordiate transformation
      !local cooridate (s p k)
      call vec_x_prod(d_i, d_n,d_s)
      d_s=d_s/vec_mag(d_s)

      call vec_x_prod(d_i, d_s, d_p)
      d_p=d_p/vec_mag(d_p)
     
      E_sp=0d0
      do i=1,3
         E_sp(1)=E_sp(1)+d_s(i)*E_vec(i)
         E_sp(2)=E_sp(2)+d_p(i)*E_vec(i)
      enddo  
       
      do i=1,2
         E2_sp(i)=E_sp(i)*conjg(E_sp(i))
      enddo
      E2_sp(3)=E2_sp(1)+E2_sp(2)

      prob_r=(E2_sp(1)*rsp(1)+E2_sp(2)*rsp(2))/E2_sp(3)

      call RANDOM_NUMBER(rv1)
      IsFlec=(rv1 .le. prob_r)
      if (IsFlec) then
         call GetFlecDir(d_n)
      else 
         call GetFracDir(d_n,r_n,IsFlec)
      endif

      !update E for polarization state
      call vec_x_prod(d_i, dir_cos,d_s)
      d_s=d_s/vec_mag(d_s)

      call vec_x_prod(d_i, d_s, d_p)
      d_p=d_p/vec_mag(d_p)
     
      call vec_x_prod(dir_cos, d_s, d_pp)
      d_pp=d_pp/vec_mag(d_p)
      
      E_sp=0d0
      do i=1,3
         E_sp(1)=E_sp(1)+d_s(i)*E_vec(i)
         E_sp(2)=E_sp(2)+d_p(i)*E_vec(i)
      enddo  
      if (IsFlec) then
         E_sp(1:2)=E_sp(1:2)*rsp(1:2)
      else
         E_sp(1:2)=E_sp(1:2)*(1.0d0-rsp(1:2))
      endif

      do i=1,3
         E_vec(i)=d_s(i)*E_sp(1)+d_pp(i)*E_sp(2)
      enddo
         
      !!debug section begins: output 
      if (n_out .lt. npsdbg) then
         n_out=n_out+1
         !write(*,"('cv1=', 2(ES12.5,x),' rv1=', ES12.5  )") cv1, rv1
         !write(*,"('csp(1)=',2(ES12.5,2x),'  csp(2)=', 2(ES12.5,2x) )") (realpart(csp(i)), imagpart(csp(i)), i=1,2)
         
         !write(*,"('r_n=',ES12.5,'  reflectance_sp=', 2(ES12.5,2x) )") r_n,rsp(1:2)
         !write(*,"('p_r=',ES12.5,'  E_sp=', 4(ES12.5,2x) )") prob_r, E_sp(1:2)
         write(*,"('polarization E=', 3('(',ES12.5,x,ES12.5,')') )") (realpart(E_vec(i)), imagpart(E_vec(i)), i=1,3)
         cv1=0d0
         do i=1,3
            cv1=cv1+E_vec(i)*dir_cos(i)
         enddo
         write(*,"('verifying  E \dot k =', 3('(',ES12.5,x,ES12.5,')') )") realpart(cv1), imagpart(cv1)

      endif      
      !!debug section ends
      
      
      !if (r_n .lt. 1) then
      !   IsFlec=(rv1 .le. rflec(1) )
      !else
      !   IsFlec=(rv1 .le. rflec(2) )
      !endif

      return

   end function IsFlec

!--------------------------utilities--------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine convert_c2s(xyz,rtp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  convert Cartesian to spherical coordinates
!  Input: xyz(3)=(x, y, z)
!  Output: rtp=(r, cos(t), phi), where -pi<phi<pi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real*8, intent(in) :: xyz(3)
      real*8, intent(out):: rtp(3)

      rtp=0d0
  
      rtp(1)=dsqrt(xyz(1)**2+xyz(2)**2+xyz(3)**2)
      if(rtp(1) .eq. 0) return    
      
      rtp(2)=xyz(3)/rtp(1)
      rtp(3)=datan2(xyz(2),xyz(1))

   end subroutine convert_c2s

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine add_scatter(d_r)     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  counting outgoing light into cosine bins
!  Input: d_r: outgoing direction
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(R_KD) :: d_r(3)
      integer ibin
       
      ibin=int((1d0+d_r(3))/bin_size)+1
      if (ibin .gt. n_tr .or. ibin .lt. 1) then
         write(*,"('ERROR: in add_scatter, d_r(3)=', f10.4)") d_r(3)
      else
         c_scatter(ibin)=c_scatter(ibin)+1
      endif
   end subroutine add_scatter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine vec_x_prod(v1,v2,v_out)     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  calculate vector cross product: v_out=v1 x v2
!  Input: v1, v2  (vectors)
!  output: v_out
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(R_KD),intent(in) :: v1(3), v2(3)
      real(R_KD),intent(out) :: v_out(3)
      
      v_out(1)= v1(2)*v2(3)-v1(3)*v2(2)
      v_out(2)= -(v1(1)*v2(3)-v1(3)*v2(1))
      v_out(3)= v1(1)*v2(2)-v1(2)*v2(1)
 
   end subroutine vec_x_prod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function  vec_mag(v1) result(mag)     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  calculate magnitude of a vector
!  Input: v1, v2  (vectors)
!  output: v_out
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit none
      real(R_KD),intent(in) :: v1(3)
      real(R_KD) :: mag
      
      integer i

      mag=0d0
      do i=1,3
         mag=mag+v1(i)**2
      enddo
      mag=dsqrt(mag)
 
   end function vec_mag


end module  mod_optic
