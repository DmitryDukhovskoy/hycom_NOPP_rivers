! ---------------------------------
! Interpolated field from
! coarser grid (ARCc0.08)
! to finer grid (ARCc0.04)
! Grid resolution is integer time
! the old grid res. 
! ---------------------------------

      MODULE INTERP

      USE utils

      IMPLICIT NONE

! ------------------------------
      CONTAINS
! ------------------------------

      SUBROUTINE INTERP_COARSE2FINE
!
! Interpolation is bilinear
!
!     Approach: for (i,j) in ARCc0.04 find corresponding
!     i8,j8 and other vertices around this node
!     interpolate 
!
      INTEGER :: i4,j4,i8,j8,i8p1,j8p1
      INTEGER :: L11, L21, L12, L22, nlnd

      REAL*4 :: dmm, dmm1, dmm2, dx, dy
      REAL*4 :: F11,F12,F22,F21, &
                Fxy1, Fxy2, Fi
      REAL*4 :: dltX, dltY

      print*,'Interpolating ...'

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
! Note that land masks in 0.04 and 0.08 do not match
! there can be water in ARCc0.04 with corresponding land in ARCc0.08
! Here, make 0 concentration everywhere where land in 0.08 or 0.04
! better approach - if LMSKout!=0, interpolate from adjacent 
! points with non-zero concentration
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

          Fxy1 = (1.-dltX)*F11+dltX*F21
          Fxy2 = (1.-dltX)*F21+dltX*F22
          Fi   = (1.-dltY)*Fxy1+dltY*Fxy2
!
! Checking:
          if (ldebug>1 .and. F11>1.e-8) then
!          if (j==2678 .and. i==536) then
            print*,'Interp subroutine:'
            print*,'ARCc0.04i=',i,' j=',j
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



      END SUBROUTINE INTERP_COARSE2FINE

      END MODULE INTERP

