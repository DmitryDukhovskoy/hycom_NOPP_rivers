      PROGRAM TRACERS2REST
! Uses restart from ARCc0.08 with all fields including tracers
! and interpolated fields except tracers into ARCc0.04
! Reads in all fields before tracers from existing ARCc0.04 - no tracers ->
! ARCc0.04 with tracers
! Reads tracers from ARCc0.08, interpolates into ARCc0.04 and
! writes them into restart ARCc0.04
!
      USE utils
      USE interp

      IMPLICIT NONE

      CHARACTER(60) :: fmat1
      CHARACTER(8)  :: aa, aar*28
      CHARACTER(80) :: str, smb*1
      CHARACTER(80) :: str1, str2

      INTEGER :: kk, lc, ktstp
      INTEGER :: irec08, cnt, &
                 irec04, cnt04
!      INTEGER :: nrecIN, nrecOUT
      INTEGER :: ich, jch, imm, &
                 irec, irec2, &
                 llist, lstr, chck, dii, djj, &
                 i0, j0, np, iM, itr, ilr
      INTEGER :: irec_out

      REAL*4 :: mday, dens, amin, amax, bmin, bmax, tstp
      REAL*4 :: amino, amaxo, dmm1, dmm2, derr
 
      fmat1='(A28,1x,A1,1x,I3,1x,I2,4x,se14.7,2x,se14.7)'
      smb = '='
      
      CALL READ_PARAM

      CALL READ_TOPO_004
      print*,'1. Topo read ARCc0.04: HH(1000,2500)=',Hout(1000,2500),&
             'LMSK=',LMSKout(1000,2500)
      print*,'2. Topo read ARCc0.04: HH(536,2678)=',Hout(536,2678),&
             'LMSK=',LMSKout(536,2678)

      CALL READ_TOPO_008
      print*,'1. Topo read ARCc0.08: HH(500,1250)=',Hin(500,1250),&
             'LMSK=',LMSKin(500,1250)
      print*,'2. Topo read ARCc0.08: HH(268,1339)=',Hin(268,1339),&
             'LMSK=',LMSKin(268,1339)

! Open ARCc0.08 restart and wind it until 
! tracer # 1
      open(11, file=trim(fina08), action='read', &
           form='unformatted', access='direct', &
           recl=nrecLin, iostat=ios)
      if (ios.ne.0) CALL ERR_STOP(fina08,1,ios)
      open(12, file=trim(finb08), action='read', &
               form='formatted', iostat=ios)
      if (ios.ne.0) CALL ERR_STOP(finb08,1,ios)

! Read header from *.b and write it to the new file
      print*,'Reading header  ',trim(frst08),'.b:'
      read(12,'(A)') str
      print*,trim(str)
      read(12,'(A)') str
      print*,trim(str)

! Read output fields from old restart ARCc0.08:
! until hit tracers
! also check if *a and *b correspond 
      irec08  = 0 
      cnt = 0
      DO 
        cnt = cnt+1
!        read(12, trim(fmat1), iostat=ios) &
!                  aar,smb,k,tstp,bmin,bmax
        read(12,'(A)',iostat=ios) cline
        if (ios.ne.0) CALL ERR_STOP(finb08,2,ios)
        lc = index(cline,'=') 
        read(cline(1:10),'(A)') aar
        read(cline(lc+2:),*) k,tstp,bmin,bmax
        print*,'restart 0.08 Rec# ',cnt,' ',trim(aar),& 
               ' layer=',k,' t.step=',tstp,&
               ' bmin=',bmin,' bmax=',bmax

        if (aar(1:6) .eq. 'tracer') exit

        irec08 = irec08+1
        read(11, rec=irec08, iostat=ios) FinT
        if (ios.ne.0) CALL ERR_STOP(fina08,2,ios)
        fin1d = FinT(1:IJDMin)
        amino=minval(fin1d, mask = fin1d .lt. 1.e20)
        amaxo=maxval(fin1d, mask = fin1d .lt. 1.e20)

        if (abs(bmin-amino)>1.e-3) then
          print*, trim(frst08),' ERR:  *b and *a do not match:'
          print*,'min *a/*b:',amino,' ',bmin
          STOP
        endif

        if (abs(bmax-amaxo)>1.e-3) then
          print*, trim(frst08),'ERR:  *b and *a do not match:'
          print*,'max *a/*b:',amaxo,' ',bmax
          STOP
        endif

      ENDDO

      print*,'# ARCc0.08 Records before Tracer =',irec08
      print*,' -------------------------------'
! Get to the right record in *b 0.08:
      rewind(12)
      read(12,'(A)') str
      read(12,'(A)') str
      DO irec=1,irec08
        read(12,'(A)',iostat=ios) cline
        if (ios.ne.0) CALL ERR_STOP(finb08,2,ios)
      ENDDO

!      pause
! Open existing ARCc0.04 files for reading
! and copying over existing fields, units 13-14
!      if (.not.file_exists(trim(fina04)) then
!        print*,'File ',trim(fina04),' should exist, not found'
!        STOP
!      endif
!      if (.not.file_exists(trim(finb)) then
!        print*,'File ',trim(finb),' should exist, not found'
!        STOP
!      endif
      open(13, file=trim(fina04), action='read', &
           form='unformatted', access='direct', &
           recl=nrecLout, iostat=ios)
      if (ios .ne. 0) CALL ERR_STOP(fina04,1,ios)

      open(14, file=trim(finb04),  action='read', &
               form='formatted', iostat=ios)
      if (ios .ne. 0) CALL ERR_STOP(finb04,1,ios)
!
!  Open output file for writing:
      open(15, file=trim(fouta), action='write', &
           form='unformatted', access='direct', &
           recl=nrecLout, iostat=ios)
      if (ios .ne. 0) CALL ERR_STOP(fouta,1,ios)

      open(16, file=trim(foutb), action='write', &
                form='formatted', iostat=ios)
      if (ios .ne. 0) CALL ERR_STOP(foutb,1,ios)

!
! Copy all interpoalted to 0.04
! restart fields before tracers:
! There should be irec08 records
      read(14,'(A)') str
      write(16,'(A)') trim(str)
      read(14,'(A)') str
      write(16,'(A)') trim(str)
      irec_out = 0
!      print*,'!!!!!!!!!!!!!!!!  '
!      print*,' For debugging skipping DO irec loop'
!      print*,' where existing 0.04 restart fields are copied over'
!      print*,' need to uncomment DO irec_out=1,irec08 in tracers2restart.F90'
!      DO irec_out=irec08+1,irec08
      DO irec_out=1,irec08
        read(14,'(A)',iostat=ios) str
        if (ios .ne. 0) CALL ERR_STOP(finb04,2,ios)
        write(16,'(A)',iostat=ios) trim(str)
        if (ios .ne. 0) CALL ERR_STOP(foutb,3,ios)
        print*,'Restart ARCc0.04: ',trim(str)
        print*,'==== irec_out=',irec_out

        read(13, rec=irec_out, iostat=ios) FoutT 
        if (ios .ne. 0) CALL ERR_STOP(fina04,2,ios)
        write(15, rec=irec_out, iostat=ios) FoutT
        if (ios .ne. 0) CALL ERR_STOP(fouta,3,ios)
        
! Check *a and *b:
        lc = index(str,'=') 
        read(str(1:28),'(A)') aar
        read(str(lc+2:),*) k,tstp,bmin,bmax
        amin = minval(FoutT, mask = FoutT .lt. 1.e20)
        amax = maxval(FoutT, mask = FoutT .lt. 1.e20)

        if (abs(bmin-amin)>1.e-3 .or. &
            abs(bmax-amax)>1.e-3) then
          print*,'ERR: Reading existing restart 0.04:'
          print*,'ERR: Min/max values in *b/*a do not match:'
          print*,'min *b, *a',bmin,amin
          print*,'max *b, *a',bmax,amax
          STOP
        endif

      ENDDO

      irec_out = irec08

! Now use restart 0.08
! Loop through all tracers and interpolate to 0.04
! Note Land mask mismatch
      DO itr=1,nTr
!        read_tracers from old restart 0.08
       DO ktstp=1,2  ! time levels 1,2
        DO ilr=1,nlyrs
          print*,'Tracer=',itr,' Layer=',ilr
!          pause
          read(12,'(A)',iostat=ios) cline
          if (ios.ne.0) CALL ERR_STOP(finb08,2,ios)
          lc = index(cline,'=') 
          read(cline(1:28),'(A)') aar
          read(cline(lc+2:),*) k,tstp,bmin,bmax
          print*,'0.08 *b: v.layer=',k,' t.step=',tstp,&
                 ' bmin=',bmin,' bmax=',bmax

          irec08 = irec08+1
          read(11, rec=irec08, iostat=ios) FinT
          if (ios.ne.0) CALL ERR_STOP(fina08,2,ios)
          fin1d = FinT(1:IJDMin)
          amin = minval(fin1d, mask = fin1d .lt. 1.e20)
          amax = maxval(fin1d, mask = fin1d .lt. 1.e20)
          fin2d = reshape(fin1d,(/IDMin,JDMin/))
          print*,'0.08 *.a: Rec #',irec08,&
            ' v.lev=',ilr,'min=',amin,'max=',amax
! Check:
          if (ilr .ne. k) then
            print*,'ERR: V. layers in *b/*a do not match: ',k,ilr
            STOP
          endif
          if (ktstp .ne. tstp) then
            print*,'ERR: Time levels in *b/*a do not match', ktstp, tstp
            STOP
          endif
          if (abs(bmin-amin)>1.e-3 .or. &
              abs(bmax-amax)>1.e-3) then
            print*,'ERR: Min/max values in *b/*a do not match:'
            print*,'min *b, *a',bmin,amin
            print*,'max *b, *a',bmax,amax
            STOP
          endif

          if (ldebug>1) &
          print*,'Tr(60,766)=',fin2d(60,766),&
                 'Tr(543:544,607)=',fin2d(543:544,607)

!        interpolate: fin2d -> fout2d
          CALL INTERP_COARSE2FINE

!        write interpolated ARCc0.04:
          irec_out = irec_out+1
          fout1d = reshape(fout2d,(/IJDMout/))
          FoutT(1:IJDMout) = fout1d
          FoutT(IJDMout+1:IJDMout+NPADout)=2.**100
!          amino=minval(fout1d, mask = fout1d .lt. 1.e20)
!          amaxo=maxval(fout1d, mask = fout1d .lt. 1.e20)
          amino=minval(FoutT, mask = FoutT .lt. 1.e20)
          amaxo=maxval(FoutT, mask = FoutT .lt. 1.e20)

!        aar = "tracer  : layer,tlevel,range "
          print*,'Rec# ',irec_out,' Write to ARCc0.04 restart *b: ',&
                  trim(aar), ilr, ktstp, amino, amaxo
          write(16, fmat1, iostat=ios) aar,smb,ilr,ktstp,amino,amaxo
          if (ios .ne. 0) CALL ERR_STOP(foutb,3,ios)

          write(15, rec=irec_out, iostat=ios) FoutT
          if (ios .ne. 0) CALL ERR_STOP(fouta,3,ios)

        ENDDO ! v. layers
       ENDDO  ! time levels =1,2
      ENDDO   ! tracers

      close(11)
      close(12)
      close(13)
      close(14)
      close(15)
      close(16)

      print*,'------------ ALL DONE ---------'
 

      END PROGRAM TRACERS2REST
