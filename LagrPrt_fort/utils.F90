! --------------------------------------------
! Read Param, prepare filenames, etc
! --------------------------------------------

      MODULE UTILS

      IMPLICIT NONE

      integer, parameter :: DP=SELECTED_REAL_KIND(14,50) ! Double precision
      integer, parameter :: SP=SELECTED_REAL_KIND(6,30)  !  single precis.
      integer, parameter :: ldebug = 0    ! <2 - low output info
!      real(SP), parameter  :: rg = 9806., hg=2.**100
!      real(SP), parameter :: R=6371.0e3, pi=3.14159265358979
      real*4, parameter  :: rg = 9806., hg=2.**100
      real*4, parameter :: R=6371.0e3, pi=3.14159265358979

      character(len=90) :: pthin, pthout, &
                           pthTOPO
      character(len=30) :: flnma, flnmb, flnmout, &
                           fnmTOPO
      character(len=130) :: fina, finb, &
                            fout, &
                            ftopa, ftopb, &
                            fgrida, fgridb

      character(len=8)  :: RGNM       ! region name HYCOM simulation
      character(len=7)  :: RegSeed    ! region name particles are seeded 
      character(len=72) :: cline
      character(len=15) :: param_file
! All fields for ARCc0.08 have "in" - those that need to be interpolated
! ARCc0.04 fields use  "out"
      integer :: i, j, k, expt, &
                 yrS, nyrs, yrE, imm
      integer :: IDM, JDM, IJDM, ios, NPAD, &
                 nrecL ! ARCc0.08
      integer :: lr1, lr2      ! min/max HYCOM layers for seeding particles
      integer :: NPmax, NPlr   ! # of particles max, 1 layer
      integer :: nlyrs, nTr, & ! # v. layers , # of passive tracers
                 nfctr         ! factor new grid resolution/old grid - integer
      integer :: Nprt

      integer, allocatable :: LMSK(:,:)


      real*4 :: dTday, NuDiff
      real*4 :: rday  ! read day

      real*4, allocatable :: fin1d(:), fout1d(:), &
                               fin2d(:,:), fout2d(:,:), &
                               FinT(:), & 
                               HH(:,:), LON(:,:), LAT(:,:),  &
                               DX(:,:), DY(:,:)
      logical :: yrpt  ! repeat same years N times


! ----------------------
!  Particles structures
! ----------------------
!
      TYPE T_PRTCL
        logical, pointer :: alive(:)

        integer, pointer :: kzlev(:)
        integer, pointer :: Nprt

        real*4, pointer :: time(:)
        real*4, pointer :: iindx(:)
        real*4, pointer :: jindx(:)
        real*4, pointer :: xcrd(:)
        real*4, pointer :: ycrd(:)
      END TYPE T_PRTCL
      TYPE(T_PRTCL), allocatable :: PRTCL(:)
     


! -----------------------
      contains
! -----------------------
      SUBROUTINE READ_PARAM

      character :: cdm1*9,cdm2*2
      integer :: ll

      call get_command_argument(1,param_file)
      if (len_trim(param_file) == 0) then
        print*,'   Param file name not specified, will use PARAM.dat'
        param_file = 'PARAM.dat'
      endif

      print*, 'Reading ',trim(param_file)
      open(unit=11, file=trim(param_file), &
               action='read', status='old', &
               form='formatted',iostat=ios)

      if (ios>0) then
        print*,'    *** ERR:  ERROR opening ', trim(param_file)
        STOP
      endif

      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),*) expt 

      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),'(A)') RGNM
      print*,'RGNM:  ',trim(RGNM)

      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),'(A)') pthin
      print*,'pthin:  ',trim(pthin)

      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),'(A)') pthout
      print*,'pthout:  ',trim(pthin)

      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),'(A)') pthTOPO

      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),'(A)') fnmTOPO

      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),*) IDM

      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),*) JDM

      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),*) nlyrs  ! # of vert. levels 
      print*,'nlyrs=',nlyrs

      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),*) nTr  ! # of tracers

      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),*) yrS  ! Start year

      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),*) nyrs  ! # of years for Lagr. tracking
      print*,'N of years=',nyrs

      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),*) imm  ! repeat N times
      yrpt=.true.
      IF (imm==0) yrpt=.false.
 
      read(11,'(A)') cline
      ll = index(cline,'=') 
       read(cline(ll+2:),*) dTday  ! time step for time integration dX/dt, days
      print*,'dTday=',dTday

      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),*) NPmax  ! max # of part

      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),*) NPlr ! # of part in 1 layer

      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),*) NuDiff  ! diffusivity, random walk m2/s

      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),*) lr1 ! 1st layer to seed prtc

      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),*) lr2 ! last layer to seed prtcl

      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),*) RegSeed ! name of seeding region

      print*,'All input read'
      close(11)

! ---------------------------------
! Define dimensions and allocate arrays
! ---------------------------------
!      fina  = trim(pthin)//trim(frst08)//'.a'
!      finb  = trim(pthf08)//trim(frst08)//'.b'
!      fout  = trim(pthfout)//trim(fnstout)//'.a'
      ftopa  = trim(pthTOPO)//trim(fnmTOPO)//'.a'
      ftopb  = trim(pthTOPO)//trim(fnmTOPO)//'.b'
      fgrida = trim(pthTOPO)//'regional.grid.a'
      fgridb = trim(pthTOPO)//'regional.grid.b'

      print*,' Domain Dim 008: i=',IDM,' j=',JDM
      print*,' Input  Dir: ',trim(pthin)
      print*,' Output Dir: ',trim(pthout)
      print*,' Topo file : ',trim(ftopa)
      print*,' Topo file : ',trim(ftopb)

      IJDM  = IDM*JDM
! Define padding HYCOM file GLBb:
      NPAD  = 4096-mod(IJDM,4096)

! Arrays with no padding 1D and 2D
      allocate(fin2d(IDM,JDM))
      allocate(fin1d(IJDM),HH(IDM,JDM), LON(IDM,JDM), LAT(IDM,JDM))
      allocate(LMSK(IDM,JDM), DX(IDM,JDM), DY(IDM,JDM))
! Arrays for reading/writing full record with npad
! for direct-access file
      allocate(FinT(IJDM+NPAD))
      inquire (iolength=nrecL) FinT

      print*,' ARCc: input rec = ',size(fin1d),' npad=',NPAD
      print*,' Input Record lengths = ',nrecL
      print*,' # of vertical levels in restart = ',nlyrs
      print*,' # of Passive Tracers = ',nTr
      print*,' Lagrangian tracking options:  '
      print*,' Start year:',yrS
      if (yrpt) then
        print*,'  repeated ',nyrs,' times'
      else
        print*,' End year: ',yrS+nyrs-1
      endif

!
! # of Lagr particles
      Nprt=(lr2-lr1+1)*NPlr
      print*,'N Particles requested = ',Nprt
      IF (Nprt.gt.NPmax) THEN
        print*,'N Partilces exceeds max:',Nprt,NPmax
        STOP
      ENDIF

      allocate (PRTCL(1))
      allocate (PRTCL(1) % alive(Nprt))
      allocate (PRTCL(1) % kzlev(Nprt))
      allocate (PRTCL(1) % time(Nprt))
      allocate (PRTCL(1) % iindx(Nprt))
      allocate (PRTCL(1) % jindx(Nprt))
      allocate (PRTCL(1) % xcrd(Nprt))
      allocate (PRTCL(1) % ycrd(Nprt))

      print*,'ARRAYS ALLOCATED  '
      print*,'  '


      END SUBROUTINE READ_PARAM
!
!
!      SUBROUTINE READ_TOPO_004
!
!      real(SP) :: dmm
!
!      print*,'READ_TOPO: ', trim(ftopa)
!
!      open(11,file=trim(ftopa),&
!              action='read',form='unformatted',&
!              access='direct',recl=nrecLout, iostat=ios)
!      if (ios.ne.0) call ERR_STOP(ftopa,1,ios)
!      
!      read(11, rec=1, iostat=ios) FoutT
!      if (ios.ne.0) call ERR_STOP(ftopa,2,ios)
!      fout1d = FoutT(1:IJDMout)
!      Hout = reshape(fout1d,(/IDMout,JDMout/))
!      close(11)
! Create Land Mask
!      DO i=1,IDMout
!      DO j=1,JDMout
!        if (Hout(i,j)<1.e20) then
!          LMSKout(i,j)=1
!        else
!          LMSKout(i,j)=0
!        endif
!      ENDDO
!      ENDDO
!
!      END SUBROUTINE READ_TOPO_004
!
      SUBROUTINE READ_TOPO_008

      real*4 :: dmm, HNg(IDM,JDM)

      print*,'READ_TOPO: ', trim(ftopa)

      open(11,file=trim(ftopa),&
              action='read',form='unformatted',&
              access='direct',recl=nrecL, iostat=ios)
      if (ios.ne.0) call ERR_STOP(ftopa,1,ios)
      
      read(11, rec=1, iostat=ios) FinT
      if (ios.ne.0) call ERR_STOP(ftopa,2,ios)
      fin1d = FinT(1:IJDM)
      HH = reshape(fin1d,(/IDM,JDM/))
      close(11)
!
! Create Land Mask
      HNg=HH*0.0
      DO i=1,IDM
      DO j=1,JDM
        if (HH(i,j)<1.e20) then
          LMSK(i,j)=1
          HNg(i,j)=-HH(i,j)
        else
          LMSK(i,j)=0
          HNg(i,j)=100.0 
        endif
      ENDDO
      ENDDO
      HH=HNg  ! depths are negative
!

      END SUBROUTINE READ_TOPO_008

      SUBROUTINE READ_LON_LAT
! Read plon, plat 
      open(11,file=trim(fgrida),                  &
              action='read',form='unformatted',   &
              access='direct',recl=nrecL, iostat=ios)
      if (ios.ne.0) call ERR_STOP(fgrida,1,ios)

      print*,'Reading Lon/Lat'
      read(11, rec=1, iostat=ios) FinT
      if (ios.ne.0) call ERR_STOP(ftopa,2,ios)
      fin1d = FinT(1:IJDM)
      LON = reshape(fin1d,(/IDM,JDM/))

      read(11, rec=2, iostat=ios) FinT
      if (ios.ne.0) call ERR_STOP(ftopa,2,ios)
      fin1d = FinT(1:IJDM)
      LAT = reshape(fin1d,(/IDM,JDM/))

!      print*,'file=',trim(fgrida)
!      print*,'IDM=',IDM, 'JDM=',JDM,'recL=',nrecL
!      print*,'LAT=',LAT(10,20),'LON=',LON(10,20)
      close(11)

!      print*,'LON/LAT read'
!      CALL wait()


      END SUBROUTINE READ_LON_LAT

!
      SUBROUTINE DISTANCE_SPHERIC_Y(y0d,x0d,Yd,Xd,ny,D)
!
! Use Vincenty formula 
! Code for 1D array, distances along Y, x is fixed
!
!      REAL(DP), parameter :: R=6371.0e3, pi=3.14159265358979
      integer, intent(in) :: ny
      REAL*4, intent(in) :: Xd, Yd(ny), x0d, y0d
      REAL*4, intent(out) :: D(ny)

      integer :: i, j
      REAL :: cf, dX, X, x0, y0, dmm1
      REAL, dimension(ny) :: dmm2, dmm3, dsgm, Y, dY

! Convert to radians
      cf=pi/180.
      y0=cf*y0d
      x0=cf*x0d

      Y=cf*Yd
      X=cf*Xd
      dY=Y-y0
      dX=X-x0
! Central angle between 2 pnts:
      dmm1=(cos(y0)*sin(dX))**2
      dmm2=(cos(Y)*sin(y0)-sin(Y)*cos(y0)*cos(dX))**2
      dmm3=abs(sin(Y)*sin(y0)+cos(Y)*cos(y0)*cos(dX))

      dsgm = atan(sqrt((dmm1+dmm2)/dmm3))
! The great-circle distance:
      D = R*dsgm

      END SUBROUTINE DISTANCE_SPHERIC_Y

       SUBROUTINE DISTANCE_SPHERIC(y0d,x0d,Yd,Xd,D)
!
! Use Vincenty formula 
! Code for a distance given y(lat) and x(lon) is fixed
!
!      REAL(SP), parameter :: R=6371.0e3, pi=3.14159265358979
       REAL*4, intent(in) :: Xd, Yd, x0d, y0d
       REAL*4, intent(out) :: D

       integer :: i, j
       REAL :: cf, dY, Y, x0, y0
       REAL :: dmm1,dmm2, dmm3, dsgm, X, dX

! Convert to radians
       cf=pi/180.0
       y0=cf*y0d
       x0=cf*x0d

       Y=cf*Yd
       X=cf*Xd
       dY=Y-y0
       dX=X-x0
! Central angle between 2 pnts:
       dmm1=(cos(y0)*sin(dX))**2
       dmm2=(cos(Y)*sin(y0)-sin(Y)*cos(y0)*cos(dX))**2
       dmm3=abs(sin(Y)*sin(y0)+cos(Y)*cos(y0)*cos(dX))

      dsgm = atan(sqrt((dmm1+dmm2)/dmm3))
! The great-circle distance:
      D = R*dsgm
!       print*,'R=',R,' dsgm=',dsgm,' D=',D

      END SUBROUTINE DISTANCE_SPHERIC
!     
!
      SUBROUTINE SEED_INI(alive, kzlev, time, &
                          iindx, jindx, xcrd, &
                          ycrd, Nprt)
! Initial seeding 
! Different for different regions
! Currently only Gr Shelf exists
      USE ifport

      logical, intent(inout) :: alive

      integer, intent(inout) :: kzlev, Nprt

      real*4, intent(inout) :: time, iindx, jindx
      real*4, intent(inout) :: xcrd, ycrd

      integer :: i1, i2, j1, j2, icc
      integer :: ix, irr, iloc, npp, IL, ipp, jpp
      integer :: II1(IJDM), Ipr(NPlr)

!
! Find points on the SW Gr Shelf
      II1=0
      npp=0
      i1=510
      i2=605
      j1=401
      j2=625
      icc=0
      DO i=i1,i2
        DO j=j1,j2
          IF (HH(i,j).ge.-600.0 .and. HH(i,j).lt.-50.0) THEN
            icc=icc+1
            CALL sub2ind(IDM,JDM,i,j,ix)
            II1(icc)=ix
          ENDIF
        ENDDO
      ENDDO
!
! Subset particles to have no more than NPlr
! If more than needed - randomly pick poins
! avoid repetition
      IF (icc>NPlr) THEN
        npp=0
        DO i=1,NPlr
          irr=1+rand(0)*(NPlr-1)
          IF (i.gt.1) THEN   ! avoid repetition of random indices
            iloc=findloc(Ipr,II1(irr))
            DO WHILE (iloc.gt.0)
              irr=1+random_number(1)*(NPlr-1)
              iloc=findloc(Ipr,II1(irr))
            ENDDO 
            Ipr(i)=II1(irr)
          ELSE
            Ipr(i)=II1(irr)
          ENDIF
          npp=npp+1
        ENDDO
      ELSE
        Ipr(1:icc)=II1(1:icc)
        npp=icc
      ENDIF

      Nprt=0
      DO k=lr1,lr2
      DO i=1,npp
        Nprt=Nprt+1
        IL=Ipr(i)
        CALL IND2SUB(IDM,JDM,ipp,jpp,IL)
        iindx(i)=ipp
        jindx(i)=jpp
        alive(i)=.true.
        xcrd(i)=LON(ipp,jpp)
        ycrd(i)=LAT(ipp,jpp)
        kzlev(i)=k
      ENDDO
      ENDDO        
     
      END SUBROUTINE SEED_INI


      SUBROUTINE SUB2IND(n,m,i,j,LI)
!  
! Convert Matrix indices to linear index
!
      integer, intent(in) :: n,m,i,j
      integer, intent(inout) :: LI

      LI=(i-1)*m+j

      END SUBROUTINE SUB2IND

      SUBROUTINE IND2SUB(n,m,i,j,LI)
! 
! Convert Linear index to matrix i,j indices
      integer, intent(in) :: n,m,LI
      integer, intent(inout) :: i,j     
      i=floor(LI/m)+1
      j=LI-(i-1)*m

      END SUBROUTINE IND2SUB

!
      SUBROUTINE ERR_STOP(fnm,fmode,ios)

      character(*), intent(in) :: fnm
      integer, intent(in) :: fmode, ios
! 1 - open
! 2 - reading
! 3 - writing
      if (fmode==1) then
        write(*,'(2A)'),'    *** ERROR opening File: ',trim(fnm)
      elseif (fmode==2) then
        if (ios>0) then 
          write(*,'(2A)'),'    *** ERROR reading: check input ',trim(fnm)
        elseif (ios<0) then
          write(*,'(2A)'),'    *** ERROR reading: E-o-F ',trim(fnm)
        endif
      else
        write(*,'(2A)'),'    *** ERROR writing to ',trim(fnm)
      endif
      write(*,'(A, I)'),' IOSTAT = ',ios

      STOP

      END SUBROUTINE ERR_STOP


       SUBROUTINE WAIT()
!  Use instead of PAUSE command - obsolete in Fortran 90/95
       implicit none
       character(1) Key
       write (*,*) 'press any key to continue, q to exit'

       read *,Key

       if (key.eq.'q') STOP

       END SUBROUTINE WAIT

       END MODULE UTILS

