! ---------------
! State variables
! ---------------
      do n=1,ncat
        read(11) fin2d   ! aicen
        CALL INTERP_FLD
        write(12) fout2d 

        read(11) fin2d   ! vicen
        CALL INTERP_FLD
        write(12) fout2d 

        read(11) fin2d   ! vsnon
        CALL INTERP_FLD
        write(12) fout2d 

        read(11) fin2d   ! Tracer # 1 - surf T
        CALL INTERP_FLD
        write(12) fout2d 
      enddo

      do k=1,ntilyr
        read(11) fin2d   ! eicen
        CALL INTERP_FLD
        write(12) fout2d 
      enddo

      do k=1,ntslyr
        read(11) fin2d   ! esnon
        CALL INTERP_FLD
        write(12) fout2d 
      enddo

! ----------------
! Velocity
! ----------------
      read(11) fin2d   ! uvel
      CALL INTERP_FLD
      write(12) fout2d 

      read(11) fin2d   ! vvel
      CALL INTERP_FLD
      write(12) fout2d 

! ------------------
! Radiation fields
! -----------------
      read(11) fin2d   ! scale_factor
      CALL INTERP_FLD
      write(12) fout2d 

!
      read(11) fin2d   ! swvdr
      CALL INTERP_FLD
      write(12) fout2d 

      read(11) fin2d   ! swvdf
      CALL INTERP_FLD
      write(12) fout2d
 
      read(11) fin2d   ! swidr
      CALL INTERP_FLD
      write(12) fout2d 

      read(11) fin2d   ! swidf
      CALL INTERP_FLD
      write(12) fout2d 

! ------------------
! Ocean Stress
! -----------------
      read(11) fin2d   ! strocnxT
      CALL INTERP_FLD
      write(12) fout2d 

      read(11) fin2d   ! strocnyT
      CALL INTERP_FLD
      write(12) fout2d 

! ------------------
! Internal stress
! -----------------
      read(11) fin2d   ! stressp_1
      CALL INTERP_FLD
      write(12) fout2d 

      read(11) fin2d   ! stressp_3
      CALL INTERP_FLD
      write(12) fout2d 

      read(11) fin2d   ! stressp_2
      CALL INTERP_FLD
      write(12) fout2d 

      read(11) fin2d   ! stressp_4
      CALL INTERP_FLD
      write(12) fout2d 

      read(11) fin2d   ! stressm_1
      CALL INTERP_FLD
      write(12) fout2d 

      read(11) fin2d   ! stressm_3
      CALL INTERP_FLD
      write(12) fout2d 

      read(11) fin2d   ! stressm_2
      CALL INTERP_FLD
      write(12) fout2d 

      read(11) fin2d   ! stressm_4
      CALL INTERP_FLD
      write(12) fout2d 

      read(11) fin2d   ! stresss12_1
      CALL INTERP_FLD
      write(12) fout2d 

      read(11) fin2d   ! stresss12_3
      CALL INTERP_FLD
      write(12) fout2d 

      read(11) fin2d   ! stresss12_2
      CALL INTERP_FLD
      write(12) fout2d 

      read(11) fin2d   ! stresss12_4
      CALL INTERP_FLD
      write(12) fout2d 

!----------------------
! Ice Mask for dynamics
! ---------------------
      read(11) fin2d   ! work1
      CALL INTERP_FLD
      write(12) fout2d 
