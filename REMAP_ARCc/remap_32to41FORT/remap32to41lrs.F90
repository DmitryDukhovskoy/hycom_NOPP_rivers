      PROGRAM REMAP2NEW_LAYERS
! Code remaps fields from HYCOM archive or restart file
! saved on a GLBb tripolar grid onto Arctic grid ARCc
!
! access can be 'stream' then no record lenght is needed
! HYCOM files are 'direct', then record length is needed
! all records have to be same lengths
! note that record length is platform-dependent
! measure of the record
!
! Dmitry Dukhovskoy, COAPS FSU, April 2016
!
      USE init_variables
      USE utils
      USE interp

      IMPLICIT NONE

!      CHARACTER(60) :: 
      CHARACTER(8)  :: cdm3, vname
      CHARACTER(90) :: aa1, aa

      INTEGER :: ii, jj, kk, &
                 nlevU, &
                 nlevV, &
                 nlevD, &
                 nlevT, &
                 nlevS

      REAL(SP) :: h1, h0

!      REAL*4, allocatable ::

      CALL READ_PARAM

      CALL READ_TOPO
      print*,'Topo, HH(800,1500)=',HH(800,1500)
      print*,'Topo, HH(1,1)=',HH(1,1)
!
! Get min z-level spacing for the new grid:
      CALL READ_VLAYERS

! Read new layer thicknesses from relax/template file
      fld = 'thknss'
      DPnew = 0.*DPnew
      CALL READ_HYCOM(finRa, finRb, fld, DPnew)
! Check
      i=800
      j=1500
      h1=sum(DPnew(:,i,j))
      h0=HH(i,j)
      print*,'RELAX file: Depths i=',i,' j=',j,'from DP/Topo: ',h1,h0
      if (abs(1.-h0/h1)>1.e-3) then
        print*,'Depth = ',h0
        print*,'Depth from relax DPnew=',h1
        print*,'Do not match'
        STOP
      endif
!
      CALL LAYER_INTRF_DPTH(DPnew,ZZin,IDM,JDM,KDM) 

! Get old layer thicknesses
      fld = 'thknss'
      DPold = 0.*DPold
      CALL READ_HYCOM(fina,finb,fld,DPold)
!      print*,'DPold(:,800,1500):'
!      print*,DPold(:,800,1500)
!
! Check
      i=ichck
      j=jchck
      h1=sum(DPold(:,i,j))
      h0=HH(i,j)
      print*,'Depths i=',i,' j=',j,'from DP/Topo: ',h1,h0
      if (abs(1.-h0/h1)>1.e-3) then
        print*,'Depth = ',h0
        print*,'Depth from DP=',h1
        print*,'Do not match'
        STOP
      endif
!
      CALL LAYER_INTRF_DPTH(DPold,ZZo,IDM,JDM,KDMold) 
! Assumption: all new layers added to 32-layer grid
! are above Lx_old layer (which corresponds to Lx_new in 41-layer)
! and these new layers are fixed-depth layers (z)
! Below Lx_new, layers match old 32-layer grid
! except for the bottom most (one target density is eliminated)
      CALL ADJUST_VLAYERS 

!      deallocate(ZZin)

      DO kk=1,nvars
        Fnew = -999.9
        vname = VARS(kk)
        print*,' '
        if (trim(vname) == 'thknss') then
          print*,kk,' Updating thickness for new grid ',trim(vname)
          CALL LAYER_DEPTH2dP(ZZn,Fnew)
          Fint_ALL(kk,:,:,:) = Fnew
          cycle
        endif
        print*,kk,' Interpolating into new grid: ',trim(vname)

        CALL READ_HYCOM(fina,finb,vname,Fold)
        CALL INTERP_FLD32to41(vname,Fold,Fnew)
        Fint_ALL(kk,:,:,:) = Fnew
      ENDDO

! Write interpolated field:
      CALL WRITE_HYCOM(fina,finb,fouta,foutb)

      print*,' ==== Done ============='

!      STOP
      END PROGRAM REMAP2NEW_LAYERS
