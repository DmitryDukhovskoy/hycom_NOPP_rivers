! -------------------------
! Bi-linear interpolation
! of CICE fields
! -------------------------

      MODULE INTERP_CICE

      USE all_variables
      IMPLICIT NONE

      CONTAINS

      SUBROUTINE INTERP_FLD
!
! Interpolate CICE fields into higher-res grid 
! from a coarser-res. grid 
! bilinear interpolation
!     Approach: for (i,j) in ARCc0.04 find corresponding
!     i8,j8 and other vertices around this node
!     interpolate 
!
      INTEGER :: i4,j4,i8,j8,i8p1,j8p1
      INTEGER :: L11, L21, L12, L22, nlnd

      REAL :: dmm, dmm1, dmm2, dx, dy
      REAL :: F11,F12,F22,F21, &
                Fxy1, Fxy2, Fi
      REAL :: dltX, dltY

      print*,'Interpolating CICE field ...'
      
      Fout2d = 0.0 
      dx = 1./float(nfctr)
      dy = 1./float(nfctr)
      DO i=1,IDMout
! Find vertices in ARCc0.08
        i8 = floor(float(i+nfctr-1)/float(nfctr))
        i8p1 = i8+1
        if (i8p1>IDMin) i8p1=IDMin
! Find corresponding indx in ARCc0.04
        i4 = (i8-1)*nfctr+1  
        dltX = abs(float(i-i4)*dx) ! Should be <=1
        if (dltX>1.000001) then
          print*,'dltX > 1', dltX
          STOP
        endif

        DO j=1,JDMout
          j8 = floor(float(j+nfctr-1)/float(nfctr))
          j8p1 = j8+1
          if (j8p1>JDMin) j8p1=JDMin
          j4 = (j8-1)*nfctr+1
          dltY = abs(float(j-j4)*dy)
          if (dltY>1.000001) then
             print*,'dltX > 1', dltX
             STOP
          endif

          F11 = fin2d(i8,j8)
          F12 = fin2d(i8,j8p1)
          F21 = fin2d(i8p1,j8)
          F22 = fin2d(i8p1,j8p1)

! Check land
          L11 = LMSKin(i8,j8)
          L12 = LMSKin(i8,j8+1)
          L21 = LMSKin(i8+1,j8)
          L22 = LMSKin(i8+1,j8+1)
          if (L11+L12+L21+L22 == 0 .or. &
              LMSKout(i,j) == 0) then
            Fout2d(i,j) = 0.0
            cycle
          else ! partial land
            nlnd = L11+L12+L21+L22
            dmm = F11*float(L11) + &
                  F12*float(L12) + &
                  F21*float(L21) + &
                  F22*float(L22)
            dmm = dmm/float(nlnd)

            if (L11==0) F11=dmm
            if (L12==0) F12=dmm
            if (L21==0) F21=dmm
            if (L22==0) F22=dmm
          endif
!
! Bilinear interpolation:
          Fxy1 = (1.-dltX)*F11+dltX*F21
          Fxy2 = (1.-dltX)*F21+dltX*F22
          Fi   = (1.-dltY)*Fxy1+dltY*Fxy2
!
! Checking:
          if (ldebug>1 .and. F11>1.e-8) then
            print*,'       '
            print*,'Interp subroutine:'
            print*,'ARCc0.04  i=',i,' j=',j
            print*,'Found coresponding ARCc0.08 indices:'
            print*,'  i8=',i8,' j8=',j8,' i8+1=',i8p1,' j8+1=',j8p1
            print*,' Depth in ARCc0.04 = ',Hout(i,j),' LMSK=',LMSKout(i,j)
            print*,' Depth in ARCc0.08 = ',Hin(i8,j8),' LMSK=',LMSKin(i8,j8)
            print*,' LMSK 0.08: i:i+1,j:j+1 ', L11,L21,L22,L12
            print*,' Input Values: i:i+1,j:j+1 ', F11,F21,F22,F12
            print*,' Interpoalted: ',Fi
!            pause
          endif

          if (abs(Fi)>abs(F11) .and. &
              abs(Fi)>abs(F12) .and. &
              abs(Fi)>abs(F21) .and. &
              abs(Fi)>abs(F22)) then
              print*,'ERR interpolation: abs(F itnerp) > interpolants'
              print*,'F interp = ',Fi
              print*,'F surrounding: ',F11,F12,F21,F22
              print*,'ARCc0.04i=',i,' j=',j
              print*,'Found coresponding ARCc0.08 indices:'
              print*,'  i8=',i8,' j8=',j8,' i8+1=',i8p1,' j8+1=',j8p1
              print*,' Depth in ARCc0.04 = ',Hout(i,j),' LMSK=',LMSKout(i,j)
              print*,' Depth in ARCc0.08 = ',Hin(i8,j8),' LMSK=',LMSKin(i8,j8)
              print*,' LMSK 0.08: i:i+1,j:j+1 ', L11,L21,L22,L12
              print*,' Input Values: i:i+1,j:j+1 ', F11,F21,F22,F12
              STOP
          endif
          Fout2d(i,j) = Fi
           
        ENDDO
      ENDDO

      END SUBROUTINE INTERP_FLD

      SUBROUTINE AGGREGATE_AREA(n)
! CICE checks aggreatge area
! in the restart field, 
! It has to be <=1
! see: ice_itd.F90
! Check that the area is <=1 and adjust 
! the area where needed in the last category
! that has been read in (n)
!
      INTEGER, intent(in) :: n
      INTEGER :: ict
      REAL*8 :: aice(IDMout,JDMout)
      REAL*8 :: armax(1), aa, an, da

      aice=0.
      DO ict=1,n
        aice=aice+aicen(ict,:,:)
      ENDDO  ! ice categories
      armax = maxval(aice)
      print*,'Cat ',n,'  Max aggregated area: ',armax(1)
 
      if (armax(1) > 1.0_DP) then
        DO i=1,IDMout
        DO j=1,JDMout
          aa = aice(i,j)
          if (aa <= 1.0_DP) cycle
          an = aicen(n,i,j)
          da = aa-1.0
          
          if (da<0.0) then
            print*,'Negative adjustment of ice area, da= ',da
            STOP
          endif

          if (da > an) then
            print*,'Category =',n
            print*,'Adjustment > area_category, da=',da,' an=',an
            print*,'Cannot adjust ice area in category ',n
            print*,'Need to adjust through previous categories'
            STOP
          endif
       
          aicen(n,i,j) = an-da
        ENDDO
        ENDDO

! Check corrected ice area:
        aice=0.
        DO ict=1,n
          aice=aice+aicen(ict,:,:)
        ENDDO  ! ice categories
        armax = maxval(aice)
        print*,'After adjustment, max aggregated ice area =',armax(1)
        print*,'  '
        
        if (armax(1) > 1.0_DP) then
          print*,'Adjustment of aggregated area failed ...'
          STOP
        endif

        fout2d = aicen(n,:,:)
      else
        print*,'Aggregated area <= 1, nothing to adjust'
      endif
          
      i=171
      j=1733
      print*,'Check area:  i=',i,' j=',j,' aice=',aice(i,j)

      END SUBROUTINE AGGREGATE_AREA

      END MODULE INTERP_CICE
