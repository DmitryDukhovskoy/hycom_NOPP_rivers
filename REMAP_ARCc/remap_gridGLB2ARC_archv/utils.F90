      MODULE UTILS
! 
! 
!
      USE all_variables

      IMPLICIT NONE

! -----------------------
      contains
! -----------------------

      SUBROUTINE GET_REMAP_INDX

      integer :: ntoto, nrec, nrecl, npad
      integer :: i, j
      real(SP), allocatable :: toto(:), dmm(:), dmm1d(:)
 
      print*,'GET_REMAP_INDX: ', trim(finGMa)

      npad = 4096-mod(IJDMA,4096)
      nrec = (IJDMA+npad)*4
      allocate(toto(npad), dmm(IJDMA+npad), dmm1d(1:IJDMA))

      print*,'npad = ',npad,' IJDMA=',IJDMA
      print*,'size dmm=',size(dmm)

      inquire (iolength=nrecl) dmm
      print*,'Determined Record length = ',nrecl

      open(11,file=trim(finGMa),&
              action='read',form='unformatted',&
              access='direct',recl=nrecl, iostat=ios)
      if (ios>0) then
        print*,'    *** ERR:  ERROR opening ',trim(finGMa)
        STOP
      endif

      read(11,rec=1, iostat=ios) dmm
      if (ios<0) STOP('READING HIT EOF UNIT=11 *.a');
      if (ios>0) STOP('READING ERROR UNIT=11 *.a')
      dmm1d = dmm(1:IJDMA)
      xmap = reshape(dmm1d,(/IDMA,JDMA/))
!      print*,'xmap =',xmap(1:5,1)

      read(11,rec=2, iostat=ios) dmm
      if (ios<0) STOP('READING HIT EOF UNIT=11 *.a');
      if (ios>0) STOP('READING ERROR UNIT=11 *.a')
      dmm1d = dmm(1:IJDMA)
      ymap = reshape(dmm1d,(/IDMA,JDMA/))

      print*,' Remap indices: GLBb -> ARCc'
      DO i=1,2
      DO j=1,3
        print*,'i=',i,' j=',j,'xmap =',xmap(i,j),&
               &' ymap =',ymap(i,j)
      ENDDO
      ENDDO
      print*,'    '

      close(11)      

      END SUBROUTINE GET_REMAP_INDX


      END MODULE UTILS
