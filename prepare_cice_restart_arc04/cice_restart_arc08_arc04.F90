      PROGRAM CICE_RESTART_ARC04
! Rewrite existing cice restart 
! from ARCc0.08 -> ARCc0.04
! All fields are linearly interpoalted
! into the new grid
      USE all_variables
      USE interp_cice

      IMPLICIT NONE

      integer :: nu_dump, &
                 istep1,  &
                 n, nflds

      real*8 :: rday, & !restart day(sec)
                fday    !forc.day (sec)

      CALL READ_PARAM

      CALL READ_TOPO_004
      print*,'Topo read ARCc0.04: HH(1000,2500)=',Hout(1000,2500),&
             'LMSK=',LMSKout(1000,2500)
      print*,'Topo HH(1000,2500) should be 206.67'

      CALL READ_TOPO_008
      print*,'Topo read ARCc0.08: HH(500,1250)=',Hin(500,1250),&
             'LMSK=',LMSKin(500,1250)

!
! Open old cice restart file
      open(11, file=trim(fina), action='read', &
           form='unformatted', iostat=ios)
      if (ios>0) call err_stop(fina,1,ios)
!
! Open new restart file
      open(12, file=trim(fouta), action='write', &
           form='unformatted', iostat=ios)
      if (ios>0) call err_stop(fouta,1,ios)

      read(11) istep1,rday,fday
      print*,'<-Readin CICE Restart ARCc0.08, header:'
      print*,'<-istep1 =',istep1,' time=',rday/86400,' ftime=',fday/86400
      if (rdayN>0) rday=rdayN*86400
      write(12) istep1,rday,fday
      print*,'->Write out CICE Restart ARCc0.04, header:'
      print*,'->istep1 =',istep1,' time=',rday/86400,' ftime=',fday/86400

      nflds = 0 ! count fields in the restart
! ---------------
! State variables
! ---------------
      print*,'State variables ...'
      do n=1,ncat
        print*,'Cat =',n

        read(11, iostat=ios) fin2d   ! aicen
        if (ios.ne.0) CALL ERR_STOP(fina,2,ios)
        CALL INTERP_FLD
        aicen(n,:,:) = fout2d ! keep area
!        i=171
!        j=1733
!        print*,'=== Check area cat:  i=',i,' j=',j,' aicen=',fout2d(i,j)
        CALL AGGREGATE_AREA(n)
!        i=171
!        j=1733
!        print*,'=== Check area cat:  i=',i,' j=',j,' aicen=',fout2d(i,j)
!        print*,'  Calling AGGREGATE AREA 2nd time'
!        CALL AGGREGATE_AREA(n)
        write(12) fout2d 
        nflds = nflds+1

        read(11, iostat=ios) fin2d   ! vicen
        if (ios.ne.0) CALL ERR_STOP(fina,2,ios)
        CALL INTERP_FLD
        write(12) fout2d 
        nflds = nflds+1

        read(11, iostat=ios) fin2d   ! vsnon
        if (ios.ne.0) CALL ERR_STOP(fina,2,ios)
        CALL INTERP_FLD
        write(12) fout2d 
        nflds = nflds+1

        read(11, iostat=ios) fin2d   ! Tracer # 1 - surf T
        if (ios.ne.0) CALL ERR_STOP(fina,2,ios)
        CALL INTERP_FLD
        write(12) fout2d 
        nflds = nflds+1
        print*,'nfields =',nflds
      enddo

      print*,'Ice enthalpy ...'
      do k=1,ntilyr
        read(11, iostat=ios) fin2d   ! eicen
        if (ios.ne.0) CALL ERR_STOP(fina,2,ios)
        CALL INTERP_FLD
        write(12) fout2d 
        nflds = nflds+1
      enddo

      print*,'Snow enthalpy ...'
      do k=1,ntslyr
        read(11, iostat=ios) fin2d   ! esnon
        if (ios.ne.0) CALL ERR_STOP(fina,2,ios)
        CALL INTERP_FLD
        write(12) fout2d 
        nflds = nflds+1
      enddo
      print*,'nfields =',nflds

! ----------------
! Velocity
! ----------------
      print*,'Velocity ...'
      read(11, iostat=ios) fin2d   ! uvel
      if (ios.ne.0) CALL ERR_STOP(fina,2,ios)
      CALL INTERP_FLD
      write(12) fout2d 
      nflds = nflds+1

      read(11, iostat=ios) fin2d   ! vvel
      if (ios.ne.0) CALL ERR_STOP(fina,2,ios)
      CALL INTERP_FLD
      write(12) fout2d 
      nflds = nflds+1
      print*,'nfields =',nflds

! ------------------
! Radiation fields
! -----------------
      print*,'Radiation Fields ...'
      read(11, iostat=ios) fin2d   ! scale_factor
      if (ios.ne.0) CALL ERR_STOP(fina,2,ios)
      CALL INTERP_FLD
      write(12) fout2d 
      nflds = nflds+1
      print*,'nfields =',nflds

!
      read(11, iostat=ios) fin2d   ! swvdr
      if (ios.ne.0) CALL ERR_STOP(fina,2,ios)
      CALL INTERP_FLD
      write(12) fout2d 
      nflds = nflds+1
      print*,'nfields =',nflds

      read(11, iostat=ios) fin2d   ! swvdf
      if (ios.ne.0) CALL ERR_STOP(fina,2,ios)
      CALL INTERP_FLD
      write(12) fout2d
      nflds = nflds+1
 
      read(11, iostat=ios) fin2d   ! swidr
      if (ios.ne.0) CALL ERR_STOP(fina,2,ios)
      CALL INTERP_FLD
      write(12) fout2d 
      nflds = nflds+1

      read(11, iostat=ios) fin2d   ! swidf
      if (ios.ne.0) CALL ERR_STOP(fina,2,ios)
      CALL INTERP_FLD
      write(12) fout2d 
      nflds = nflds+1
      print*,'nfields =',nflds

! ------------------
! Ocean Stress
! -----------------
      print*,'Ocean stress ...'
      read(11, iostat=ios) fin2d   ! strocnxT
      if (ios.ne.0) CALL ERR_STOP(fina,2,ios)
      CALL INTERP_FLD
      write(12) fout2d 
      nflds = nflds+1

      read(11, iostat=ios) fin2d   ! strocnyT
      if (ios.ne.0) CALL ERR_STOP(fina,2,ios)
      CALL INTERP_FLD
      write(12) fout2d 
      nflds = nflds+1

! ------------------
! Internal stress
! -----------------
      print*,'Internal Stress ...'
      read(11, iostat=ios) fin2d   ! stressp_1
      if (ios.ne.0) CALL ERR_STOP(fina,2,ios)
      CALL INTERP_FLD
      write(12) fout2d 
      nflds = nflds+1

      read(11, iostat=ios) fin2d   ! stressp_3
      if (ios.ne.0) CALL ERR_STOP(fina,2,ios)
      CALL INTERP_FLD
      write(12) fout2d 
      nflds = nflds+1

      read(11, iostat=ios) fin2d   ! stressp_2
      if (ios.ne.0) CALL ERR_STOP(fina,2,ios)
      CALL INTERP_FLD
      write(12) fout2d 
      nflds = nflds+1

      read(11, iostat=ios) fin2d   ! stressp_4
      if (ios.ne.0) CALL ERR_STOP(fina,2,ios)
      CALL INTERP_FLD
      write(12) fout2d 
      nflds = nflds+1

      read(11, iostat=ios) fin2d   ! stressm_1
      if (ios.ne.0) CALL ERR_STOP(fina,2,ios)
      CALL INTERP_FLD
      write(12) fout2d 
      nflds = nflds+1

      read(11, iostat=ios) fin2d   ! stressm_3
      if (ios.ne.0) CALL ERR_STOP(fina,2,ios)
      CALL INTERP_FLD
      write(12) fout2d 
      nflds = nflds+1

      read(11, iostat=ios) fin2d   ! stressm_2
      if (ios.ne.0) CALL ERR_STOP(fina,2,ios)
      CALL INTERP_FLD
      write(12) fout2d 
      nflds = nflds+1

      read(11, iostat=ios) fin2d   ! stressm_4
      if (ios.ne.0) CALL ERR_STOP(fina,2,ios)
      CALL INTERP_FLD
      write(12) fout2d 
      nflds = nflds+1

      read(11, iostat=ios) fin2d   ! stresss12_1
      if (ios.ne.0) CALL ERR_STOP(fina,2,ios)
      CALL INTERP_FLD
      write(12) fout2d 
      nflds = nflds+1

      read(11, iostat=ios) fin2d   ! stresss12_3
      if (ios.ne.0) CALL ERR_STOP(fina,2,ios)
      CALL INTERP_FLD
      write(12) fout2d 
      nflds = nflds+1

      read(11, iostat=ios) fin2d   ! stresss12_2
      if (ios.ne.0) CALL ERR_STOP(fina,2,ios)
      CALL INTERP_FLD
      write(12) fout2d 
      nflds = nflds+1

      read(11, iostat=ios) fin2d   ! stresss12_4
      if (ios.ne.0) CALL ERR_STOP(fina,2,ios)
      CALL INTERP_FLD
      write(12) fout2d 
      nflds = nflds+1

!----------------------
! Ice Mask for dynamics
! ---------------------
      read(11, iostat=ios) fin2d   ! work1
      if (ios.ne.0) CALL ERR_STOP(fina,2,ios)
      CALL INTERP_FLD
      write(12) fout2d 
      nflds = nflds+1
      print*,'nfields =',nflds


      close(11)
      close(12)


      END PROGRAM CICE_RESTART_ARC04
