      PROGRAM REMAP_GLB_TO_ARC
! Code remaps fields from HYCOM mean archive (archm) 
!  or restart file
! saved on a GLBb tripolar grid onto Arctic grid ARCc
! In archm u-, v-vel are total velocities 
! in *.b files - need to delete several lines that
! are in archm files (kebtrop, kemix, ...)
! also the first field 'montg1' k = 0 and dens=0
! in archv k = 6 dens=34.00 (ref. dens)
! U and V fields: In the ARCc - the half that is 
! rotated from the GLBb grid (up from ARCc j= 1249)
! all vectors (U and V components) need 
! to be rotated by 180degrees
! to match + direction in the ARCc grid
!
! access can be 'stream' then no record lenght is needed
! HYCOM files are 'direct', then record length is needed
! all records have to be same lengths
! note that record length is platform-dependent
! measure of the record
!
! Dmitry Dukhovskoy, COAPS FSU, April 2016
! Jan 2017 - fixed bugs with rotating U,V in the upper half of ARCc domain
!          
!
      USE all_variables
      USE utils

      IMPLICIT NONE

      CHARACTER(60) :: fmat1
      CHARACTER(8)  :: aa
      CHARACTER(80) :: str, smb*1
      CHARACTER(80) :: str1, str2

      INTEGER :: i,j,iG,jG
      INTEGER :: npadg, npada, nrecIN, nrecOUT
      INTEGER :: ich, jch, imm, k, &
                 nrecl1, irec, nrecl2, irec2, &
                 llist, lstr

      REAL*4 :: mday, dens, amin, amax, bmin, bmax
      REAL*4 :: amino, amaxo
      REAL*4, allocatable :: FinT(:), fin1d(:), fin2d(:,:), &
                             FoutT(:), fout1d(:), fout2d(:,:), &
                             ubtrop(:,:), vbtrop(:,:)


      CALL READ_PARAM

! Read remapping indices
      CALL GET_REMAP_INDX

!
! HYCOM *.b format:
      fmat1='(A8,1x,A1,2x,I9,1x,f10.3,1x,I2,1x,f6.3,2x,se14.7,2x,se14.7)'
!
! Define padding HYCOM file GLBb:
      npadg = 4096-mod(IJDMG,4096)
! Define size/padding for output file, ARCc:
      npada = 4096-mod(IJDMA,4096)
! Arrays for reading/writing full record with npad
! for direct-access file
      allocate(FinT(IJDMG+npadg),FoutT(IJDMA+npada))
! Arrays with no padding 1D and 2D
      allocate(fin1d(IJDMG), fin2d(IDMG,JDMG))
      allocate(fout1d(IJDMA), fout2d(IDMA,JDMA), &
               ubtrop(IDMA,JDMA), vbtrop(IDMA,JDMA))
      inquire (iolength=nrecl1) FinT
      inquire (iolength=nrecl2) FoutT

      ubtrop = hg
      vbtrop = hg

      fin1d=0.
      print*,'  '
      print*,' GLBb: size 1 rec=',size(fin1d),' npad GLB=',npadg
      print*,' ARCc: size 1 rec=',size(fout1d),' npad ARC=',npada
      print*,' Record lengths, GLBb=',nrecl1,'  ARC=', nrecl2

      open(11, file=trim(fina), action='read', &
           form='unformatted', access='direct', &
           recl=nrecl1, iostat=ios)
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
           recl=nrecl2, iostat=ios)
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
        if (i==4) &
         str='depth_GLBb0.08_11; GLBb0.08 GFOS-3.1 remapped to ARCc'
        write(14,'(A)') trim(str)
        print*,'i=',i,'  ',trim(str)
      ENDDO
      read(12,'(I5,A)') ich,str1
      print*,'I dimensions, input file ich=',ich
      read(12,'(I5,A)') jch, str2
      print*,'J dimensions, input file jch=',jch
      write(14,'(I5,A)') IDMA,trim(str1)
      write(14,'(I5,A)') JDMA,trim(str2)
      read(12,'(A)') str
      print*,str
      str = 'field       time step  model day  k  dens        min              max'
      write(14,'(A)') trim(str)

! Read output fields:
! in *b: Skip extra fields from archm not needed in archv
! and make nonzero k and dens for montg1
      print*,trim(fmat1)
      irec=0
      irec2=0
      DO
        read(12,trim(fmat1),iostat=ios) &
                  aa,smb,imm,mday,k,dens,bmin,bmax
!        print*,'*b: aa=',aa,' mday=',mday,' k=',k,' dens=',dens,&
!               ' bmin=',bmin,' bmax=',bmax
!        pause
        if (ios<0) exit ! EOF
        if (ios>0) STOP('*** ERR: READING ERROR UNIT=12')
        if (trim(aa)=='') exit
        print*,'*b: aa=',aa,' mday=',mday,' dens=',dens,&
               ' bmin=',bmin,' bmax=',bmax

        irec=irec+1
        read(11, rec=irec, iostat=ios) FinT
!        read(11, iostat=ios) dmm
        if (ios<0) STOP('READING HIT EOF UNIT=11 *.a')
        if (ios>0) STOP('READING ERROR UNIT=11 *.a')

        if ( aa(1:5) == 'montg' ) then
          k = 6
          dens = 34.0
        endif

        if ( aa(1:5) == 'tmix ' .or. &
             aa(1:5) == 'smix ' .or. &
             aa(1:5) == 'thmix' .or. &
             aa(1:5) == 'umix ' .or. &
             aa(1:5) == 'vmix ' .or. &
             aa(1:5) == 'kemix' .or. &
             aa(1:5) == 'kebtr' .or. &
             aa(1:5) == 'k.e. ' .or. &
             aa(1:5) == 'densi') then
          print*,'skipping  ',aa(1:9)
          cycle 
        endif

        fin1d = FinT(1:IJDMG)
        fin2d = reshape(fin1d,(/IDMG,JDMG/))

        amin=minval(fin1d, mask = fin1d .lt. 1.e20)
        amax=maxval(fin1d, mask = fin1d .lt. 1.e20)
        print*,'Input *.a:',aa,' k=',k,'min=',amin,'max=',amax
!        print*,'Check: Fin=',FinT(1800:1801)
!        pause
!
!  Remap GLBb -> ARCc
!  rotate U and V
! for u-vel, v-vel - subtract depth-average U,V (btrop)
! transition occurs at j=1249
! all i indices are the same (in columns) for j=1:1248
! and then they switch  to another value for 1249:end
! do not change land values 
        DO i=1,IDMA
        DO j=1,JDMA
          iG = int(xmap(i,j))
          jG = int(ymap(i,j))
          fout2d(i,j) = fin2d(iG,jG)
          if ( aa(2:5) .eq. '_btr' .or. &
               aa(2:5) .eq. '-vel' ) then
            if ( abs(iG-int(xmap(i,1))) > 1 .and. &
                fin2d(iG,jG) < 1.e20) fout2d(i,j) = -fin2d(iG,jG)
          endif

! Checking
!          if (i.eq.1006 .and. j.eq.1009) then
!            print*,'Checking i/j= ',i,j
!            print*,'iG =',iG,' jG=',jG
!            print*,'fout2d=',fout2d(i,j)
!          endif

        ENDDO
        ENDDO

        if ( aa(1:5) .eq. 'u_btr' ) ubtrop = fout2d
        if ( aa(1:5) .eq. 'v_btr' ) vbtrop = fout2d
        if ( aa(1:5) .eq. 'u-vel' ) fout2d = fout2d-ubtrop
        if ( aa(1:5) .eq. 'v-vel' ) fout2d = fout2d-vbtrop

! Write out
        irec2  = irec2+1
        fout1d = reshape(fout2d,(/IJDMA/))
        FoutT(1:IJDMA)=fout1d
        FoutT(IJDMA+1:IJDMA+npada)=2.**100
        amino=minval(fout1d, mask = fout1d .lt. 1.e20)
        amaxo=maxval(fout1d, mask = fout1d .lt. 1.e20)
        print*,'Rec: ',irec,' Write ',aa,' Min=',amino,' Max=',amaxo
!        pause
        write(14,fmat1) aa,smb,imm,mday,k,dens,amino,amaxo
        write(13, rec=irec2) FoutT
      ENDDO



      close(11)
      close(12)
      close(13)
      close(14)

      print*,'Completed: archv from GLBb archm '
      print*,trim(fouta)
      print*,trim(foutb)

      STOP
      END PROGRAM REMAP_GLB_TO_ARC
