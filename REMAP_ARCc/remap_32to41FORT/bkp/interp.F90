      MODULE INTERP
!
!
!
      USE init_variables
      USE utils

      IMPLICIT NONE


! ----------------------
      contains
! ----------------------

      SUBROUTINE ADJUST_VLAYERS
! First: update new vertical grid (ZZn)
! above Layer Lx_new - fixed depths, as given in 
!       dp0k -1D depth array of minimum z-thickness in deep ocean
! below Lx_new - layers identically match the 32-layer grid
!         for layer >=Lx_old
! Target densities of the layers Lx_new and Lx_old 
! are the same
! In some locations, Lx_old (=14) happens to be
! shallower than Lx_new(=24) (although should match)
! In this case, stick with the new grid depth
! the old grid will be assumed at the same depth
! during field interpolation
!
! the bottom-most layer combines 2 bottom layers of 32 grid
! Bathymetry for both vertical grids should match!
! 
! ZZo (3D array of interf. depths) -> ZZn (3D with 41+1 interf. depths)
! both have negative depths
!% 41-layer fields:
!% DP and thicknesses need to be changed
!% in the deep ocean to make z-level fixed
!% for layers < Lx_new
!% layers >=Lx_new should match the old grid
!% (layers >=Lx_old)
      real :: h0, sgmH, zXnew, zXold
      real :: Zold(KDMold+1), Znew(KDM+1), &
              zn(KDM+1), dHb(KDM+1), imm(1)
      real :: sgmDP
      integer :: Ib, k1, n1, ii, jj
      integer :: i0, j0

      logical :: lchck
!
      lchck = .false. ! check info flag

      print*,'  '
      print*,'===  adjust_vlevels :'
      ZZn = ZZin ! all sigma (shallow) regions not changed
                 ! adjust only deep ocean points
      i0 = ichck
      j0 = jchck
! z-layer spacing is read from blkdat.input
      DO ii=1,IDM
      DO jj=1,JDM
        h0   = -HH(ii,jj)
        Zold = ZZo(:,ii,jj)
        Znew = ZZin(:,ii,jj)
        zn   = Znew

        if (lchck) then
        if (ii.eq.i0 .and. jj.eq.j0) then
          print*,'Check: i, j=', ii, jj,' h0=',h0
!          print*,'Znew=',Znew
!          print*,'Zold=',Zold
        endif
        endif
! Bottom intrf of Lx_old should be collocated with the bottom intrf 
! of Lx_new (ZZn(Lx_new+1))
! in some   cases, this is not true
!% Check that bottom interface Lx_old is not shallower than layer Lx_new
!% if it is, then ignore Lx_old and match layers (>Lx_old) in both 
!% grid only for those Old Layers (32 lr) that are deeper or same depth
!% as the Lx_new
!% then Layers in the new grid are added as z-level grid until
!% the corresponding layer in the old grid becomes
!% deeper, after that match the rest of the layers
        sgmDP= Znew(nsigma+1)
        sgmH = abs(1.-abs(sgmDP)/abs(h0)) ! sigma-coord points
        if (sgmH>1.e-2) then ! deep ocean, not sigma coord.
          zn = 0.
          dHb = abs(abs(Znew)-abs(h0))
!          imm = maxloc(dHb, (dHb>1.e-2)) ! bottom index
          DO k=1,KDM
            if (dHb(k) <= 1.e-2) exit
          ENDDO
          Ib = k-1
          if (Ib<Lx_new) then ! bottom shallower than Lx_new
            DO i=2,Ib+1
              zn(i) = zn(i-1)-dp0k(i-1)
            ENDDO
            DO i=Ib+1,KDM+1
              zn(i)=h0
            ENDDO
            k1 = KDM+1
            n1 = KDMold+1
            zXold = h0-100.
            zXnew = h0
          else    ! bottom below Lx_new, 
            DO i=2,Lx_new
              zn(i) = zn(i-1)-dp0k(i-1) ! intrf dpths above Lx_new
            ENDDO
            k1 = Lx_new
            n1 = Lx_old
            zXold = Zold(Lx_old+1) ! bottom intrfc Lx_old lr
            zXnew = zn(Lx_new)-dp0k(Lx_new)
            zXnew = max(zXnew,h0)
          endif ! bottom vs Lx_new

          if (lchck) then
          if (ii.eq.i0 .and. jj.eq.j0) then
            print*,'Check 1: i, j=', ii, jj
            print*,'Ib=',Ib,' Lx_new=',Lx_new
            print*,'k1=',k1,' n1=',n1
            print*,'zXnew=',zXnew,' zXold=',zXold
            print*,'zn=',zn
          endif
          endif

!% zXold should = zXnew, if zXold shallower zXnew
!% This means that the targ. dens sigma41(Lx_new)< density at this depth
!% Start adding z-levels using z-level spacing dp0k (dPdeep)
!% specified in the blkdat until the new (41) hycom layer
!% becomes shallower than the old (32) hycom layer
          DO WHILE (nint(zXold)>nint(zXnew))
            k1 = k1+1
            n1 = n1+1
            if (k1>KDM) exit
            zn(k1) = zn(k1-1)-dp0k(k1-1)
            if (abs(zn(k1)-h0)<0.01 .or. &
                zn(k1)<h0) then    ! hit bottom
              zn(k1:KDM+1) = h0
              k1 = KDM+1
              n1 = KDMold+1
              exit
            endif
            zXold = Zold(n1)
            zXnew = zn(k1)
          ENDDO  ! WHILE loop
!
          if (lchck) then
          if (ii.eq.i0 .and. jj.eq.j0) then
            print*,'Check 2: i, j=', ii, jj
            print*,'k1=',k1,' n1=',n1
            print*,'zXnew=',zXnew,' zXold=',zXold
            print*,'With no deep lrs: zn=',zn
          endif
          endif
! Match deeped layers if needed
! Note that last deep layer in old grid (32)
! disappeared in 41-layer version
          if (k1<=KDM) zn(k1+1:KDM+1)=Zold(n1+1:KDMold)
          zn(KDM+1) = Zold(KDMold+1) ! match the bottom

          ZZn(:,ii,jj) = zn
        endif  ! deep ocean

      ENDDO
      ENDDO ! i/j loops

      if (lchck) then
        print*,'Adjusted V.Layers: '
        print*,'i=',i0,' j=',j0
        print*,'ZZn(:,i,j) =',ZZn(:,i0,j0)
      endif
          
      END SUBROUTINE ADJUST_VLAYERS
!
!
!
      SUBROUTINE INTERP_FLD32to41(vname, Fold,FF)
! Interpolate 3D fields from 32 ---> 41 layers
!Interpolation is "mass" weighted, where mass is the layer 
!% thickness
!% ZM*, ZZ* - mean depths and interface depths of vertical layers
!% Lx_new & Lx_old  - layer # in 41-layer grid where target densities
!%        start to match the 32 layer grid = tdens(Lx_old)
!%  -------------------------------------------------- 
!%  For each interface in new t.dens. find
!%  nearest upper and bottom interfaces in old t.dens.
!%
!%   Old v. layers             New v. layers
!%   -------------
!%       Z(kold) *           ----------   
!%                            Zn(knew)   *
!%   -------------
!%
!%     Z(kold+1) *           ----------   
!%                            Zn(knew+1) *
!%   -------------
      character(*), intent(in) :: vname
      real, intent(in) :: Fold(:,:,:)
      real, intent(inout) :: FF(:,:,:)

      logical :: lchck
      integer :: ii, jj, k1, k2, kk, &
                 nZ, cc, kf, cc1, cpp, &
                 rpp, cpp_old, k3
      real :: sgmDP, sgmH, dHn, &
              z1_new, z2_new, pU1, pU2, &
              dZS1, dZS2, dZfull, Ufull, &
              dZ1, Utot_new, Utot_old, &
              derr, dff, h0, dZtot, dZZ, &
              dZold, dZnew, zz0
      real :: ffold(KDMold), Zold(KDMold+1), &
              Znew(KDM+1), ffnew(KDM)

      lchck = .false.
      FF = hg
      cc1  = 0
      cpp_old = 0
      print*,'Interpolating 32 ----> 41 lr : ',trim(vname)

      DO ii=1,IDM
      DO jj=1,JDM
        cc1 = cc1+1
        cpp = nint(real(cc1)/real(IJDM)*100.)
        rpp = mod(cpp,5)
        if (rpp == 0 .and. cpp_old.ne.cpp) then
          cpp_old = cpp
          print*,'==  ',trim(vname),' done ',cpp,'%'
        endif
!        print*,'cpp=',cpp,' rpp=',rpp,' cc1=',cc1

        h0   = -HH(ii,jj)
! Only ocean points, skip land
        if (abs(h0)>0.01*hg) cycle

        Zold = ZZo(:,ii,jj)
        ffold = Fold(:,ii,jj)
! Adjust truncation error of bottom depth
        DO k3 = 1,KDMold+1
          zz0 = abs(Zold(k3)-h0)
          if (zz0<0.01) Zold(k3)=h0-1.e-7
        ENDDO
          
!        print*,'h0',h0
! 41-lr field:
        Znew = ZZin(:,ii,jj)
!        sgmDP= Znew(nsigma+1)
!        sgmH = abs(1.-abs(sgmDP)/abs(h0)) ! sigma-coord points
!% In theory, interpolation is needed
!% only in the new layers added above Lx_new
!% and all layers below should be identical
!% to the old 32-layer grid (below Lx_old)
!% however, in several locations there was
!% some mismatch, so have to interpolate
!% over all layers
        ffnew = ffold(KDMold)
        ffnew(Lx_new:KDM)=ffold(Lx_old:KDMold-1)
!        if (ii==ichck .and. jj==jchck) then
!          print*,'cc1=',cc1,' cpp=',cpp,' rpp=',rpp
!          print*,'Check 1: ffnew=',ffnew
!        print*,'Check 1: ffold=',ffold
!        endif

        DO k=1,KDM
          dHn = abs(Znew(k+1)-Znew(k)) !Lr. thickness
          if (dHn<1.e-2) cycle
          z1_new = Znew(k)
          z2_new = Znew(k+1)
! Fix problem finding k2 near the bottom
! due to truncation error
!          z2_new = floor(abs(z2_new*100.))*(-0.01) !trunc error
          zz0 = abs(z2_new-h0)
          if (zz0<0.01) z2_new=h0+1.e-6

! Find corresponding layer surfaces in old grid (32lr)
! that contain new surfaces
          k1 = 0
          k2 = 0
          DO kk=1,KDMold+1
            if (Zold(kk)<z1_new .and. k1==0) k1=kk-1
            if (Zold(kk)<=z2_new .and. k2==0) then
              k2=kk-1
              exit
            endif

          ENDDO
!          if (ii==184 .and. jj==464 .and. k==9) then
!            print*,'kk=',kk,' z1_new=',z1_new,' z2_new=',z2_new
!            print*,'dHn=',dHn,' Znew(k+1)=',Znew(k+1),' k=',k
!          endif

!          if (ii==177 .and. jj==616) then
!            print*,vname,': i=',ii,' j=',jj, 'k=',k,' h0=',h0
!            print*,'k1=',k1,' k2=',k2
!            print*,'Pause 1'
!            pause
!          endif

          if (k1==0 .or. k2==0) & 
          CALL ERR_MSG01(k1,k2,vname,ii,jj,k,h0,Znew,&
                         dHn,z1_new,z2_new,Zold)
 ! how many complete old z-levels inside the new layer:
          nZ = k2-k1-1

          if (k1>KDMold+1 .or. k2>KDMold+1) &
          CALL ERR_MSG02(ii,jj,k1,k2)

          if (nZ<0) then ! new layer is completely within old lr
            ffnew(k)=ffold(k1)
          else
!% top and bottom of new layer are in different layers of the old grid
!% nZ = can be 0 or >0 complete "old" layers within the new layer
!% Delta thickness of partial layers      
            dZS1 = abs(z1_new-Zold(k1+1))
            dZS2 = abs(z2_new-Zold(k2))
            pU1  = ffold(k1)*dZS1
            if (k2<=KDMold) then
              pU2 = ffold(k2)*dZS2
            else
              pU2 = 0.
            endif
            
            dZfull = 0.
            Ufull  = 0.
            cc     = 0
            DO kf=k1+1,k2-1
              cc  = cc+1
              dZ1 = abs(Zold(kf)-Zold(kf+1))
              dZfull = dZfull+dZ1
              Ufull = Ufull+ffold(kf)*dZ1
            ENDDO  ! kf
! CHeck total thickness of the layer
            dZtot = dZS1+dZfull+dZS2
            derr  = abs(abs(dZtot/(z2_new-z1_new))-1)
            if (derr>1.e-2) &
              CALL ERR_DERR (nZ,dZtot,z1_new,z2_new,&
                          k1,k2,dZS1,dZS2,ii,jj,k)

! Interpolated value in the new layer:
            ffnew(k) = (pU1+Ufull+pU2)/dZtot
          endif ! if nZ<0

!          if (ii==184 .and. jj==464 .and. k<=9) then
!            print*,':: k1=',k1,' k2=',k2,' k=',k
!            print*,':: Zold(k1)=',Zold(k1),' Zold(k2)=',Zold(k2)
!            print*,'dZS1=',dZS1,' dZS2=',dZS2,' dZfull=',dZfull
!            PAUSE
!          endif
        ENDDO   ! k - vert.layers

! ---------------------------
! Check interpolation:
! ---------------------------
        Utot_old = 0.
        Utot_new = 0.
        CALL CHECK_TOTAL(Zold, Znew, ffold, &
                         ffnew, Utot_old, &
                         Utot_new, dZold, dZnew)
! If interpolated and old profiles differ > eps
! try to adjust it
! if it is too big, something is wrong - exit
        derr = abs(1.-Utot_new/Utot_old)
        dff  = abs(Utot_new-Utot_old)
        if (derr>0.5 .and. dff>0.5) then
          print*,'ERR: Utot big derr=',derr,' dff=',dff
          CALL ERR_UTOT(vname,ii,jj,h0,Utot_old, &
                          Utot_new, dZnew, dZold, ffold, &
                          ffnew, Znew, Zold)
        endif

        if (derr>0.1 .and. dff>1.e-2) then
          print*,'Interp Adjustment: derr=',derr,' dff=',dff
          CALL ADJUST_INTERP(Zold,Znew,ffold,ffnew, &
                             Utot_old,Utot_new,h0)
!
! Check if the error has been reduced
          CALL CHECK_TOTAL(Zold,Znew,ffold, &
                           ffnew, Utot_old, &
                           Utot_new, dZold, dZnew)
          derr = abs(1.-Utot_new/Utot_old)
          dff  = abs(Utot_new-Utot_old)
          print*,'====  After Adjusted: derr=',derr,' dff=',dff
!          pause
        endif 

        if (derr>0.1 .and. dff>1.e-2) &
          CALL ERR_UTOT(vname,ii,jj,h0,Utot_old, &
                        Utot_new,dZnew, dZold, ffold, &
                        ffnew, Znew, Zold)


        FF(:,ii,jj) = ffnew

        if (lchck .and. ii==ichck .and. jj==jchck) then
          print*,'Check interpolation i, j =',ii,jj,' ',vname
          print*,'interpolated FF=',FF(:,ii,jj)
        endif

      ENDDO ! jj
      ENDDO ! i

              


      END SUBROUTINE INTERP_FLD32to41

      END MODULE INTERP
