      PROGRAM REFORMAT_B
! Rewrite *.b file created by archv2data3z.f
! to match format required in relaxv.f 
! to process z-level fields in *.a and *.b files
!
!
      IMPLICIT NONE

      character(240) :: flnm,frmt*80,cline, fout,&
                        cline_new
      character(79)  :: cline1, cline2, cline3, &
                        cline4, cline5, cc1, cc2

      real :: hmin, hmax, zlev

      integer :: ios, k, i, j

      read(*,'(A)') flnm
      print*,' input file: ',flnm
      read(*,'(A)') cline1
      read(*,'(A)') cline2
      read(*,'(A)') cline3
      read(*,'(A)') cline4
      read(*,'(A)') cline5

      fout=trim(flnm)//'F'

      open(11,file=trim(flnm),&
           action='read', form='formatted',&
           iostat=ios)
      if (ios>0) then
        print*,'    *** ERR:  ERROR opening ',trim(flnm)
        STOP
      endif

      open(41,file=trim(fout),&
           action='write', form='formatted',&
           iostat=ios)
      if (ios>0) then
        print*,'    *** ERR:  ERROR opening ',trim(fout)
        STOP
      endif

! Add header lines
      write(41,'(A)') trim(cline1)
      write(41,'(A)') trim(cline2)
      write(41,'(A)') trim(cline3)
      write(41,'(A)') trim(cline4)
      write(41,'(A)') trim(cline5)

      DO k=1,9999
        read(11,'(A)', end=100) cline
        I = index(cline,'=')
        read(cline(I+1:I+8),*) zlev
        cc1 = cline(1:I-1)//',range ='
        J = index(cline,':')
        read(cline(J+1:),*) hmin,hmax
        J = index(cline,'[')
        cc2 = trim(cline(11:J-2))
        write(cline_new,110) trim(cc2),trim(cc1),zlev,hmin,hmax
        print*,trim(cline_new)
        write(41,'(A)') trim(cline_new)
      ENDDO

  100 CONTINUE
      close(11)
      close(41)

 110  format (5x,A,3x,A,f8.2,1x,2g15.6)


      STOP
      END PROGRAM REFORMAT_B
