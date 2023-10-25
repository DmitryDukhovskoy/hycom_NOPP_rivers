      MODULE UTILS
!
!
!
      USE init_variables

      IMPLICIT NONE

! -----------------------
      contains
! -----------------------


      SUBROUTINE READ_TOPO

      integer :: ntoto, nrecl, npad

      print*,'READ_TOPO: ', trim(ftopa)

      inquire (iolength=nrecl) dmm
      print*,'Determined Record length = ',nrecl

      open(11,file=trim(ftopa),&
              action='read',form='unformatted',&
              access='direct',recl=nrecl, iostat=ios)
      if (ios.ne.0) call ERR_STOP(ftopa,1,ios)
      
      read(11, rec=1, iostat=ios) dmm
      if (ios.ne.0) call ERR_STOP(ftopa,2,ios)
      dmm1d = dmm(1:IJDM)
      HH = reshape(dmm1d,(/IDM,JDM/))
      close(11)

      END SUBROUTINE READ_TOPO
!
!
      SUBROUTINE READ_VLAYERS
! Read dp0k - min z-level layer thicknesses
! same as dPdeep in sub_adjust_vlayers.m
! and Targ. densities for new grid
! from blkdat.input
      character(100) :: cheader
      integer :: nn, mm, kk, ix
      real*4 :: rmm

      print*,' ' 
      print*, '         READING z-spacing   ',fblk41
      open(11,file=trim(fblk41), &
              action='read',form='formatted', &
              iostat=ios)
      if (ios.ne.0) call ERR_STOP(fblk41,1,ios)
! Read the header lines:
      DO i=1,6 
        read(11,'(A)',iostat=ios) cheader
        if (ios.ne.0) call ERR_STOP(fblk41,2,ios)
        print*,trim(cheader)
      ENDDO

      read(11,*,iostat=ios) nn, cheader
      if (ios.ne.0) call ERR_STOP(fblk41,2,ios)
      read(11,*,iostat=ios) mm, cheader
      if (ios.ne.0) call ERR_STOP(fblk41,2,ios)

      print*,'IDM=',nn,' JDM=',mm
      if (nn.ne.IDM .or. mm.ne.JDM) then
        print*,'ERR: I,J=',nn,mm,' in blkdat.input do not match ', IDM,JDM
        STOP
      endif

      read(11,'(A)') cheader
      print*,trim(cheader)
      read(11,'(A)') cheader
      print*,trim(cheader)
      read(11,*) kk,cheader ! number of layers
      read(11,'(A)') cheader      ! N. hybrid lrs
      print*,trim(cheader)
      read(11,*) nsigma, cheader  
      print*,nsigma,' ',trim(cheader)

      if (kk.ne.KDM) then
        print*,'ERR: v.layers =',kk,' expected =',KDM
        STOP
      endif

! 1D array allocate:
      allocate(dp0k(KDM), ds0k(nsigma), TDENS(KDM))
! 
! Read deep z-level spacing min thickness, m
      DO i=1,KDM
        read(11,*) rmm,cheader
        dp0k(i)=rmm
      ENDDO
!
! Read shallow z-level spacing (sigma-coord. layers):
      DO i=1,nsigma
        read(11,*) rmm,cheader
        ds0k(i)=rmm
      ENDDO
      print*,'min/max deep    z-level',minval(dp0k),maxval(dp0k)
      print*,'min/max shallow z-level',minval(ds0k),maxval(ds0k)

! Read New target densities
      DO
        read(11,'(A)',iostat=ios) cheader
        if (ios.ne.0) call ERR_STOP(fblk41,2,ios)
        ix = index(cheader,"vsigma")
        if (ix>0) exit
      ENDDO

      DO k=1,KDM
        read(11,*) rmm, cheader
        TDENS(k)=rmm
      ENDDO
      print*,'min/max TDENS = ',minval(TDENS),maxval(TDENS)

      close(11)

      END SUBROUTINE READ_VLAYERS
!
      SUBROUTINE READ_HYCOM(fina,finb,fld,FF)
!  reads hycom binary archive files (model output), 
!  returns specified field 'fld'
!  FF is a 3D array(z,x,y) has to be 
!  correct size
      character(*), intent(in) :: fina, finb
      character(*), intent(in) :: fld
      real*4, intent(inout) :: FF(:,:,:)

      logical :: lthck

      character(len=85) :: ch1
      integer :: mm, nn, ll, ix, jx, irec, &
                 nrecl, &
                 a1, a2, a3
      
      lthck = .false.
      if (fld == 'thknss') lthck = .true.
      a1 = size(FF,1)
      a2 = size(FF,2)
      a3 = size(FF,3)

!      print*,'a1=',a1,' a2=',a2,' a3=',a3
      print*,' ::: Reading ',fld,' from ',trim(fina)

      inquire (iolength=nrecl) dmm
      open(unit=20, file=trim(fina), &
           action='read',form='unformatted',&
           access='direct',recl=nrecl, iostat=ios)
      if (ios.ne.0) call ERR_STOP(fina,1,ios)

      open(unit=21, file=finb, &
           action='read', status='old', &
           form='formatted',iostat=ios)
      if (ios.ne.0) call ERR_STOP(finb,1,ios)

! Read the header lines:
      print*,'  '
!      print*,trim(finb)
      DO 
        read(21,'(A)',iostat=ios) cheader
        if (ios.ne.0) call ERR_STOP(finb,2,ios)
        print*,trim(cheader)
        ix = index(cheader,"time step")
        if (ix>0) exit
        ix = index(cheader,"idm ")
        jx = index(cheader,"jdm ")
        if (ix>0) read(cheader,*) nn,ch1
        if (jx>0) read(cheader,*) mm,ch1
      ENDDO

      if (nn.ne.a2 .or. mm.ne.a3) then
        print*,'Size of F idm=',a2, ' jdm=',a3
        print*,'Check idm ',nn,' & jdm ',mm,' from ',finb
        STOP
      endif
      
      ll = 0
      irec = 0
      DO 
        read(21,*,iostat=ios) fld0,ch1
        if (ios<0) exit
        if (ios>0) call ERR_STOP(finb,2,ios)
        irec=irec+1
        ix = index(fld0,trim(fld)) ! find field in *b file
        if (ix>0) then
          ll=ll+1
!          print*,'ll=',ll,' ',fld0
          read(20, rec=irec, iostat=ios) dmm
          if (ios.ne.0) call ERR_STOP(fina,2,ios)
          dmm1d = dmm(1:IJDM)
          FF(ll,:,:) = reshape(dmm1d,(/IDM,JDM/))
        endif
      ENDDO

      if (ll.ne.a1) then
        print*,'# of v.levels in F=',a1
        print*,'does not match v.lev. in file =',ll
        print*,'Check ',finb
        STOP
      endif
      
      if (lthck) FF=FF/rg ! convert Pa -> m

      close(20)
      close(21)

      END SUBROUTINE READ_HYCOM


      SUBROUTINE LAYER_INTRF_DPTH(dPr0,ZZ,n,m,l)
! Calculate interface depths of the layers
! Given layer thickness (m) dPr0
! Interface depths - negative, 0 - on land
      real, intent(in) :: dPr0(:,:,:)
      real, intent(inout) :: ZZ(:,:,:)
      real :: dP(l,n,m)
      real :: zr(l), dp0(l), dp1(l), zmm(n,m)
      real :: h0, hz0, drr, hmin, hmax, zmin, zmax
      integer :: n,m,l, i0, j0
      dP=dPr0
      print*,'Layer Interface Depths'
!
! Convert all land dP to 0
      zr=zr*0.
      DO i=1,n
      DO j=1,m
        dp0=dP(:,i,j)
        dp1=pack(dp0, (dp0<1.e20), zr)
        dP(:,i,j)=dp1
!        print*,'dP0=',dp0
!        print*,'dP=',dp1
!         pause
      ENDDO
      ENDDO
!  
! Interf depths
      ZZ=ZZ*0.
      DO k=1,l
        zmm=abs(dP(k,:,:))
        ZZ(k+1,:,:)=ZZ(k,:,:)-zmm
      ENDDO

! Check
      print*,'Layers Intrf. Depths: checking ...'
      DO i=1,n
      DO j=1,m
        h0=abs(HH(i,j))
        if (h0>1.e20) h0=0. ! land
        hz0=abs(ZZ(l+1,i,j))
!        print*,'hz0=',hz0,' ZZ=',ZZ(l+1,i,j),' l=',l
!        print*,'ZZ=',ZZ(:,i,j)
!        pause
        if (abs(hz0)<1.e-8) then
          if (abs(h0)>1.e-8) then
            print*,'ZZ=',ZZ(:,i,j)
            print*,'DP32=', dP(:,i,j)
            print*,'ERR: Intrf depths: i=',i,' j=',j
            print*,'LAND?   Z btm=',hz0,' H topo=',h0
            pause
          endif
        else
          drr=abs(1.-h0/hz0)
          if (drr>1.e-6) then
            print*,'ZZ=',ZZ(:,i,j)
            print*,'DP32=', dP(:,i,j)
            print*,'ERR: Intrf depths: i=',i,' j=',j
            print*,'   Z btm=',hz0,' H topo=',h0
            pause
          endif
        endif
      ENDDO
      ENDDO

      zmin=minval(abs(ZZ(l,:,:)), (abs(ZZ(l,:,:))>0))
      zmax=maxval(abs(ZZ(l,:,:)), (abs(ZZ(l,:,:))>0))
      hmin=minval(abs(HH), (abs(HH)<1.e20))
      hmax=maxval(abs(HH), (abs(HH)<1.e20))

      print*,'Interf Dpeths:  ZZ check OK:'
      print*,'min/max ZZ bttm  =  ',zmin,' ',zmax
      print*,'min/max HH, Topo =  ',hmin,' ',hmax
      print*,'    '

!      i0=800
!      j0=20
!      print*,'Checking ZZin, i=',i0,' j=',j0
!      print*,'ZZin=',ZZin(:,i0,j0)
!      pause

      END SUBROUTINE LAYER_INTRF_DPTH
!
!  
      SUBROUTINE LAYER_DEPTH2dP(ZZ,DP)
! Calculate layer thickness (in pressure units, dP)
! From layer interface depths (m)
      real, intent(in) :: ZZ(:,:,:)
      real, intent(inout) :: DP(:,:,:)

      integer :: klev, kk, ii, jj
      real  :: h0
      real :: dmm(IDM,JDM)

      klev = size(DP,1)
      DO kk=1,klev
        dmm = abs(ZZ(kk+1,:,:)-ZZ(kk,:,:))
        DO ii=1,IDM
        DO jj=1,JDM
          h0 = HH(ii,jj)
          if (dmm(ii,jj)<1.e-2 .or. h0>0.01*hg) then
            dmm(ii,jj)=hg
          else
            dmm(ii,jj)=dmm(ii,jj)*rg
          endif
        ENDDO
        ENDDO

        DP(kk,:,:) = dmm
      ENDDO  ! kk - v. layers

      END SUBROUTINE LAYER_DEPTH2dP


      SUBROUTINE WRITE_HYCOM(finaOLD, finbOLD,&
                             finaNEW, finbNEW)
! Write interpolated fields into
! HYCOM binary
! Use old *a, *b files as templates
! NOTE: horizontal dimensions in OLD and NEW fields 
! are the same IDM x JDM
! only vertical grid is different
      character(*), intent(in) :: finaOLD,finbOLD,&
                                  finaNEW, finbNEW

      character(60) :: fmat1 ! *b IO format
      character(90) :: cmm
      character(8)  :: aa, smb*1, vname
      character(80) :: str1, str2

      integer :: nl, imm, jmm, cnt, kk, &
                 irec, irec2, nrecl, iv, &
                 nlr
      integer :: NLEV(nvars)

      real(SP) :: dens, bmin, bmax, td0, mday
      real(SP) :: FinT(IJDM+npad), & !read/write 1D field+padding
                  fout1d(IJDM),     & ! 1D field w/o padding
                  fout2d(IDM,JDM)     ! 2D field
!
      inquire (iolength=nrecl) FinT
      print*,'HYCOM 2D+padding record length =',nrecl
!      if (nrecl.ne.nrecl0) then
!        print*,'ERR: write_hycom: Record length ??'
!        print*,'I/O record length=',nrecl
!        print*,'Expected  record length=',nrecl0
!        STOP
!      endif

      NLEV = 0 ! layer counter for interpolated 3D fields

! Input files:
      open(unit=20, file=trim(finaOLD), &
           action='read',form='unformatted',&
           access='direct',recl=nrecl, iostat=ios)
      if (ios.ne.0) call ERR_STOP(finaOLD,1,ios)

      open(unit=21, file=trim(finbOLD), &
           action='read',form='formatted',&
           iostat=ios)
      if (ios.ne.0) call ERR_STOP(finbOLD,1,ios)
! Output files:
      open(unit=40, file=trim(finaNEW), &
           action='write',form='unformatted',&
           access='direct',recl=nrecl, iostat=ios)
      if (ios.ne.0) call ERR_STOP(finaNEW,1,ios)

      open(unit=41, file=trim(finbNEW), &
           action='write',form='formatted',&
           iostat=ios)
      if (ios.ne.0) call ERR_STOP(finbNEW,1,ios)


! Read/write Header
      DO nl=1,10
        read(21,'(A)',iostat=ios) ccline
!        read(21,'(A)',iostat=ios) str1
        if (ios.ne.0) call ERR_STOP(finbOLD,2,ios)
        write(*,'(A)'), trim(ccline)
!        print*,trim(str1)

        if (nl==4) ccline=cheadr
        if (nl==8) read(ccline,*) imm, cmm     
        if (nl==9) read(ccline,*) jmm, cmm     

        write(41,'(A)') trim(ccline)
      ENDDO    

      if (imm.ne.IDM .or. jmm.ne.JDM) then
        print*,'IDM and JDM in input *.b not consistent'
        print*,'IDM =',IDM,' JDM=',JDM
        print*,'File ',finbOLD
        print*,'I dim=',imm,' J dim=',jmm
      endif;

! Read/write fields by layers
      fmat1='(A8,1x,A1,2x,I9,1x,f10.3,1x,I2,1x,f5.2,2x,se14.7,2x,se14.7)'

! Write 2D/barotropic fields that do not need
! interpolation     
      irec  = 0
      irec2 = 0
      cnt   = 0
      DO
        cnt=cnt+1
        read(21,trim(fmat1),iostat=ios) &
             vname,smb,imm,mday,kk,dens,bmin,bmax
        if (ios<0) exit ! EOF
        if (ios>0) call ERR_STOP(finbOLD,2,ios)
        if (trim(vname)=='') exit
        print*,cnt,'*b: var=',vname,' mday=',mday,' dens=',dens,&
               ' bmin=',bmin,' bmax=',bmax

        if (kk>0 .and. trim(vname).ne.'montg1') exit
        if (trim(vname) == 'u-vel.') then
          print*,'ERR: write_hycom: '
          print*,'unexpected 3D field found in *.b:',vname
          STOP
        endif

        irec=irec+1
        read(20, rec=irec, iostat=ios) FinT
        if (ios/=0) call ERR_STOP(finaOLD,2,ios)
!        fin1d = FinT(1:IJDM)
        write(40, rec=irec, iostat=ios) FinT
        if (ios/=0) call ERR_STOP(finaNEW,3,ios)
        print*,'Write ',vname,bmin,bmax
        write(41,fmat1) vname,smb,imm,mday,kk,dens,bmin,bmax

      ENDDO  ! barotropic part
 
! Baroclinic part, interpolated fields
! NOTE: order of the fields in VARS should match *.b
      DO nlr=1,KDM     ! v.layers
        td0=TDENS(nlr)
        DO iv=1,nvars  ! variables
          vname = VARS(iv)
          NLEV(iv)=NLEV(iv)+1
          irec = irec+1
          fout2d = Fint_ALL(iv,nlr,:,:)
          fout1d = reshape(fout2d,(/IJDM/))
!          FinT = -999999.9
          FinT(1:IJDM) = fout1d
          FinT(IJDM+1:IJDM+npad) = hg
          bmin=minval(fout1d, mask = fout1d .lt. 1.e20)
          bmax=maxval(fout1d, mask = fout1d .lt. 1.e20)

          write(40, rec=irec, iostat=ios) FinT
          if (ios/=0) call ERR_STOP(finaNEW,3,ios)
          write(41,fmat1) vname,smb,imm,mday,nlr,td0,bmin,bmax
        ENDDO
      ENDDO

      close(20)
      close(21)
      close(40)
      close(41)
      print*,' Saved VARS, layers=',NLEV
      write(*,'(2A)'),'::Fields saved: ',trim(fouta)

      END SUBROUTINE WRITE_HYCOM
!
!
      SUBROUTINE CHECK_TOTAL(Zold,Znew,ffold,ffnew,Utot_old,Utot_new, &
                             dZold, dZnew)
! Check depth-integrated old and new interpolated fields
      real, intent(in) :: Zold(:), Znew(:), ffold(:), ffnew(:)
      real, intent(out) :: Utot_old, Utot_new
      real, intent(out) :: dZold, dZnew
      real :: dZZ
      integer :: k

      Utot_old = 0.
      Utot_new = 0.
      dZold = 0.
      DO k=1,KDMold
        dZZ = abs(Zold(k+1)-Zold(k))
        dZold = dZold+dZZ
        Utot_old = Utot_old+ffold(k)*dZZ
      ENDDO

      dZnew = 0.
      DO k=1,KDM
        dZZ = abs(Znew(k+1)-Znew(k))
        dZnew = dZnew+dZZ
        Utot_new = Utot_new+ffnew(k)*dZZ
      ENDDO

      END SUBROUTINE CHECK_TOTAL
!
!
      SUBROUTINE ADJUST_INTERP(Zold,Znew,ffold,ffnew, &
                             Utot_old,Utot_new,h0)
! For small errors in depth-integrated interpolated profiles
! and old profiles, adjust interpolation by redistributing
! error over the profile
      real, intent(in) :: Zold(:), Znew(:), ffold(:), &
                          Utot_old, Utot_new, h0
      real, intent(inout) :: ffnew(:)
      real :: derr_dh, dZZ, adjF
      integer :: k

      derr_dh = (Utot_old-Utot_new)/abs(h0)
      DO k=1,KDM
        dZZ=abs(Znew(k+1)-Znew(k))
        if (dZZ<1.e-4) exit
        adjF = derr_dh*dZZ
        ffnew(k) = ffnew(k)+adjF
      ENDDO

      END SUBROUTINE ADJUST_INTERP
!
!
!
      SUBROUTINE ERR_MSG01(k1,k2,vname,ii,jj,k,h0,Znew,&
                           dHn,z1_new,z2_new,Zold)
! k1 and k2 - layer interface numbers in the old grid
! k1 - interface above the top interface of the new layer
! k2 - interface above the bottom interface of the new layer
! k1 or k2 cannot be 0
      integer, intent(in) :: k1,k2,ii,jj,k
      character, intent(in) :: vname
      real, intent(in) :: Znew(:), dHn, z1_new, z2_new,&
                          Zold(:), h0
      print*,'      '
      print*,'=======  ERR_MSG01: '
      print*,'interp  ERR: k1 or k2 is zero, k1, k2=',k1,k2
      print*,'Couldnot locate old grid layer interface' 
      print*,vname,': i=',ii,' j=',jj, 'k=',k,' h0=',h0
      print*,'Znew(k+1)=',Znew(k+1),' dHn=',dHn
      print*,'z1_new=',z1_new,' z2_new=',z2_new
      print*,'Zold=',Zold
      STOP

      END SUBROUTINE ERR_MSG01
!
!
      SUBROUTINE ERR_MSG02(ii,jj,k1,k2)
! Layer Interface number (k1 or k2) exceeds
! The total number of interfaces KDMold+1
      integer, intent(in) :: ii, jj, k1, k2

      print*,'   '
      print*,'==== ERR_MSG02: '
      print*,'interp: ERR: k1 or k2>KDMold+1, error layer index'
      print*,'i=',ii,'j=',jj
      print*,'k1=',k1,' k2=',k2,' KDMold+1=',KDMold+1
      STOP

      END SUBROUTINE ERR_MSG02
!
!
      SUBROUTINE ERR_DERR(nZ,dZtot,z1_new,z2_new,&
                          k1,k2,dZS1,dZS2,ii,jj,k)
! CHecking layer thickness in the new (41lr) dZtot
! it shoud be identical to the summation of the partial
! dZS1 (top) and bottom (dZS2) and full layers of the old (32)
! grid that are within the new layer
      integer, intent(in) :: nZ,k1,k2,ii,jj,k
      real, intent(in) :: dZtot, z1_new, z2_new, &
                          dZS1, dZS2

      print*,'   '
      print*,'==== ERR_DERR: '
      print*,'interp: ERR: check total layer'
      print*,' thickness over ',nZ,' nZ'
      print*,' dZtot = ',dZtot,' m'
      print*,'z1_new=',z1_new,' z2_new=',z2_new
      print*,'k1=',k1,' k2=',k2,' dZS1=',dZS1,' dZS2=',dZS2
      print*,'i=',ii,' j=',jj,' k=',k
      STOP

      END SUBROUTINE ERR_DERR
!
!
      SUBROUTINE ERR_UTOT(vname,ii,jj,h0,Utot_old, &
                          Utot_new,dZnew, dZold, ffold, &
                          ffnew, Znew, Zold)
! CHecking interplolated and original fields at point ii,jj
! by integrating over the total depth, should be same
      character, intent(in) :: vname
      integer, intent(in) :: ii,jj
      real, intent(in) :: h0, Utot_old, Utot_new, &
                          dZnew, dZold, ffold(:), &
                          ffnew(:), Znew(:), Zold(:)

      print*,'ERR: interp: Total value '
      print*,' is not conserved after interpolation'
      print*,vname, ' i=',ii,' j=',jj,' h0=',h0
      print*,'Utot_old=',Utot_old,' Utot_new=',Utot_new
      print*,'Sum dZnew=',dZnew,' dZold=',dZold
      print*,'ffold =',ffold
      print*,'ffnew=',ffnew
      print*,'Znew=',Znew
      print*,'Zold=',Zold
      STOP

      END SUBROUTINE ERR_UTOT

      END MODULE UTILS
