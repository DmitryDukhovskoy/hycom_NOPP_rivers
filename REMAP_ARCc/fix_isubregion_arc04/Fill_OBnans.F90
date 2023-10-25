      PROGRAM FILL_OB_NANS
! isubregion generates -999 along the OBs
! for U and V, which is probably leftover from
! unsuccessful bi-linear interpolation near the OBs
!
! Fill nan/-999 values at the OB
! in the nest file
! by copying the next row/column over
!
! Name of PARAM??.dat file is provided
! at the call of the program
! in case several processes are running at the
! same time to avoid overriding PARAM.dat
! ./fillnans.x PARAM107.dat
! If not specified, PARAM.dat is used
!
      USE all_variables
 !     USE utils

      IMPLICIT NONE

      CHARACTER(60) :: fmat1
      CHARACTER(8)  :: aa
      CHARACTER(80) :: str, smb*1
      CHARACTER(80) :: str1, str2

      INTEGER :: i,j, kk, cnt
!      INTEGER :: nrecIN, nrecOUT
      INTEGER :: ich, jch, imm, k, &
                 irec, irec2, &
                 llist, lstr, chck, dii, djj, &
                 i0, j0, np, iM

      REAL*4 :: mday, dens, amin, amax, bmin, bmax
      REAL*4 :: amino, amaxo, dmm1, dmm2, derr
 
      
      CALL READ_PARAM

      CALL READ_TOPO
      print*,'Topo read: HH(1000,2500)=',HH(1000,2500)
! Read the nest file

! HYCOM *.b format:
!      fmat1='(A28,1x,A1,1x,I3,1x,I2,4x,se14.7,2x,se14.7)'
      fmat1='(A8,1x,A1,2x,I9,1x,f10.3,1x,I2,1x,f6.3,2x,se14.7,2x,se14.7)'


      open(11, file=trim(fina), action='read', &
           form='unformatted', access='direct', &
           recl=nrecL1, iostat=ios)
      if (ios>0) then
        print*,'ERROR openning HYCOM input *a'
        STOP
      endif
      open(12, file=trim(finb), action='read', &
               form='formatted', iostat=ios)
      if (ios>0) then
        print*,'ERROR openning HYCOM input *b'
        STOP
      endif

      open(13, file=trim(fouta), action='write', &
           form='unformatted', access='direct', &
           recl=nrecL1, iostat=ios)
      if (ios>0) then
        print*,'ERROR openning output for writing *a'
        STOP
      endif
      open(14, file=trim(foutb), action='write', &
               form='formatted', iostat=ios)
      if (ios>0) then
        print*,'ERROR openning output for writing *b'
        STOP
      endif
!
! Read header from *.b and write it to the new file
      print*,'Reading header: '
      DO i=1,7
        read(12,'(A)') str
        write(14,'(A)') trim(str)
        print*,'i=',i,'  ',trim(str)
      ENDDO
      read(12,'(I5,A)') ich,str1
      print*,'I dimensions, input file ich=',ich
      read(12,'(I5,A)') jch, str2
      print*,'J dimensions, input file jch=',jch
      write(14,'(I5,A)') IDM,trim(str1)
      write(14,'(I5,A)') JDM,trim(str2)
      read(12,'(A)') str
      print*,str
      write(14,'(A)') trim(str)

! Read output fields:
      print*,trim(fmat1)
      irec=0
      irec2=0
      DO
        read(12,trim(fmat1),iostat=ios) &
                  aa,smb,imm,mday,k,dens,bmin,bmax
         print*,'*b: aa=',aa,' mday=',mday,' dens=',dens,&
               ' bmin=',bmin,' bmax=',bmax
!        pause
        if (ios<0) exit ! EOF
        if (ios>0) STOP('*** ERR: READING ERROR UNIT=12')
        if (trim(aa)=='') exit
        irec=irec+1
        read(11, rec=irec, iostat=ios) FinT
!        read(11, iostat=ios) dmm
        if (ios<0) STOP('READING HIT EOF UNIT=11 *.a')
        if (ios>0) STOP('READING ERROR UNIT=11 *.a')
        fin1d = FinT(1:IJDM)
        fin2d = reshape(fin1d,(/IDM,JDM/))

        amin=minval(fin1d, mask = fin1d .lt. 1.e20)
        amax=maxval(fin1d, mask = fin1d .lt. 1.e20)
        print*,'Input *.a:',aa,' k=',k,'min=',amin,'max=',amax
!
        if (aa(2:5) .eq. '-vel' .or. &
            aa(2:5) .eq. '_btr') then
          CALL FILL_NANS
! Now put land mask values for U,V variables
          if ( aa(2:5) .ne. '_btr' ) then
            fout2d(:,1)   = hg
            fout2d(:,JDM) = hg
            fout2d(1,:)   = hg
            fout2d(IDM,:) = hg 
          endif
        else
          fout2d = fin2d ! no changes in the input field
        endif

! Write out
        irec2  = irec2+1
        fout1d = reshape(fout2d,(/IJDM/))
        FoutT(1:IJDM)=fout1d
        FoutT(IJDM+1:IJDM+npad)=2.**100
        amino=minval(fout1d, mask = fout1d .lt. 1.e20)
        amaxo=maxval(fout1d, mask = fout1d .lt. 1.e20)
        print*,'Write ',aa,' Min=',amino,' Max=',amaxo
!        pause
        write(14,fmat1) aa,smb,imm,mday,k,dens,amino,amaxo
        write(13, rec=irec2) FoutT
      ENDDO



      close(11)
      close(12)
      close(13)
      close(14)

      print*,'------------ ALL DONE ---------'
 

      END PROGRAM FILL_OB_NANS
