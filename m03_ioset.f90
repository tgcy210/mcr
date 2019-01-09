! I/O modules

module mod_ioset
   
   integer, parameter :: MAXLEN=4001
   integer, parameter :: MAXUN=20
   character (LEN=MAXLEN) :: in4line=''
   character (LEN=MAXLEN) :: lineseg=''

   integer :: lineflag=0

   character (LEN=3) ::  cmtchar="/!#"
   character (LEN=1) ::  secmark="$"
   character (LEN=3) ::  omitchar=",="
   character (LEN=54) :: letterchar='/!abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
   character (LEN=13)   :: digchar='0123456789.+-'
   character (LEN=2)    :: scichar='eE'
   character (LEN=4)  ::  fidochar='RrQqZz'
   real, allocatable :: fidobank(:)
   integer :: num_need=0, num_read=0

   logical :: ex_inp=.false.
   !file info
   character(len=50) :: file_inp="params_mcr.inp"
   integer :: U_inp=21, U_log=31
   
   !display/error/warning
   character(len=400) :: disp_msg(5)=""
   character(len=80) :: cur_var=""
   
contains

subroutine ReadInp

   implicit none
   integer ierr,i
   character(len=3) :: head_char="" 

   inquire(file=file_inp, exist=ex_inp)
   if (.not. ex_inp) then  
      write(disp_msg(1),"('input file not found:', A)") trim(file_inp)
      write(disp_msg(2),"('using default parameters')") 
      call ShowMsg(2)
      return
   endif      

   open(file=trim(file_inp), unit=U_inp)
   ierr=0
   i=0
   
   do while(ierr .eq. 0)
      i=i+1
      read(U_inp, "(A)", iostat=ierr ) head_char
      if (head_char(1:1) .eq. secmark .and. ierr .eq. 0) then
         if( index('0123456789',head_char(3:3)) .eq. 0) head_char(3:3)=' '
         select case (head_char(2:3))
            case ('0 ')
               call ProcInp_sec0
            case default
               write(disp_msg(1),"('section number is not supported on line:', I3)") i
               call ShowMsg(1)
         end select
      endif
   enddo

end subroutine ReadInp

!*********************************************************
!  processing input section 0
!*********************************************************

subroutine ProcInp_sec0
   use mod_optic 
   
   implicit none

   integer, parameter :: num_keyword =3
   integer len_eff
   integer i,j,k
   integer ierr
   integer :: cur_key=0, next_key=0, pos_key=0
   integer :: pos_start=1, pos_end=1

   character (LEN=7) list_keyword(num_keyword)

   data list_keyword /'numips', 'Evecxy', 'npsdbg'/
   integer :: cout_keyword(num_keyword)=0

   real(R_KD)  rval(10) 

   ierr=0

   loop_line : do while (ierr .eq. 0)
      in4line=''
      read(U_inp,"(A)",iostat=ierr) in4line
      if(ierr .ne. 0) then      !reading failed
        if (cur_key .eq. 0) then  
           exit loop_line        !ignore reading error if currently no keywords on deck
        else
           write(disp_msg(1),"('error when reading keyword: ', A )") &
              trim(list_keyword(cur_key))
           call ShowMsg(1)
           stop
        endif
      endif

      in4line=adjustl(in4line)
      len_eff=len_trim(in4line)

      !check blank line
      if(len_eff .eq. 0) cycle loop_line

      !check comment line
      if( index(cmtchar, in4line(1:1)) .ne. 0) cycle loop_line

      !check section mark
      if ( index(secmark, in4line(1:1)) .ne. 0) then
         backspace(U_inp)
         exit loop_line
      endif

      !pre-processing the line
      do i=1, len_eff
         if (index(omitchar,in4line(i:i)) .ne. 0) in4line(i:i)=' '
         if (index(cmtchar, in4line(i:i)) .ne. 0 .and. i .gt. 1 .and. in4line(max(i-1,1):max(i-1,1)) .eq. ' ') then
            len_eff=i-1
            exit
         endif 
      enddo

      pos_start=1
      pos_end=len_eff

      !if a fresh line, try to locate the first keyword
      if (cur_key .eq. 0) then 
         do i=1, num_keyword
            pos_key=index(in4line(pos_start:len_eff), trim(list_keyword(i)) )
            if (pos_key  .ne. 0) then
               cur_key=i
               pos_start=pos_key+len_trim(list_keyword(i))
               exit
            endif
         enddo
      endif

      !if failed to find a keyword, unwanted line, continue on next line
      if (cur_key .eq. 0) then
         write(disp_msg(1),"(' Warning: unrecognised line in Section 0: ')") 
         call ShowMsg(1)
         write(*,"(2x, A)") trim(in4line)
         cycle loop_line
    endif

    !begin line scan
    line_scan : do while( pos_start .le. len_eff)

      !check for another keyword
      !if found, pos_end=pos_key-1
      !if not, pos_end=len_eff
      next_key=0
      pos_end=len_eff
      do i=1, num_keyword
         pos_key=index(in4line(pos_start:len_eff), trim(list_keyword(i)) )
         if (pos_key  .ne. 0) then
            pos_end=pos_start+pos_key-2
            next_key=i
            exit
         endif
      enddo

      if(pos_end .eq. 0) then
         write(disp_msg(1),"('error while reading ', A )") &
           trim(list_keyword(cur_key)) 
         call ShowMsg(1)
       endif
      lineseg=in4line(pos_start:pos_end)

      select case (cur_key)

        case (1)
          cur_var='numips (number of groups for the the finer group structure)'
          rval=0d0
          read(lineseg,*,err=1002,end=1002) rval(1)
          !error check
          if(rval(1) .ge. 1) then
            num_p=int(rval(1))
            cur_key=0
            cout_keyword(1)= cout_keyword(1)+1
          else
            !to be refine error message
            write(disp_msg(1),"('invalid numips value ')") 
            call ShowMsg(1)
          endif

          write(disp_msg(1),"(' numisp (number of trial incident light per processor) ',' : ', I10 )" )  num_p
          call ShowMsg(1)
        case (2)
          !tgtgrp : number of groups of target group structure
          cur_var='Evecxy (x and y components of incident E )'
          read(lineseg,*,err=1002,end=1002) rval(1:4)
          E_ini(1)=dcmplx(rval(1), rval(2))
          E_ini(2)=dcmplx(rval(3), rval(4))
          E_ini(3)=0d0
          cur_key=0
          cout_keyword(2)= cout_keyword(2)+1
          write(disp_msg(1),"('E_xy=', 4(ES12.5,2x) )") rval(1:4)
          call ShowMsg(1)
        case (3)
          cur_var='npsdbg (debug output for first npsdbg reflection/refraction)'
          rval=0d0
          read(lineseg,*,err=1002,end=1002) rval(1)
          !error check
          if(rval(1) .ge. -1.0) then
            npsdbg=int(rval(1))
            cur_key=0
            cout_keyword(3)= cout_keyword(3)+1
          else
            !to be refine error message
            write(disp_msg(1),"('invalid npsdbg value ')") 
            call ShowMsg(1)
          endif

        case default
           write(disp_msg(1),"('invalid keyword ')") 
           call ShowMsg(1)
      end select

      if(next_key .ne. 0) then
        pos_start=pos_end+len_trim(list_keyword(next_key))+1
        cur_key=next_key
      else
        pos_start=len_eff+1
      endif

    enddo line_scan

  enddo loop_line


  !To do: check section read status
  do i=1, 0
    if(cout_keyword(i) .eq. 0) then
      write(disp_msg,"('required keyword not found in Section 0', A)") &
        trim(list_keyword(i))
      call ShowMsg(1)
    endif
  enddo

  
  return
  1002 write(disp_msg(1),"('error while processing input file: ', A)") trim(file_inp) 
  call ShowMsg(1)
  stop


end subroutine ProcInp_sec0

!*********************************************************
!utility routines
!*********************************************************   
subroutine ShowMsg(num_line)
   use mod_mympi , only: myid

   implicit none
 
   integer num_line
   
   integer i
   if (myid .eq. 0) then 
     do i=1, num_line
        write(*,"(A)") trim(disp_msg(i))
     enddo
   endif

end subroutine ShowMsg

end module mod_ioset

