! ------------------------------
! Declare variables, paths, etc.
! ------------------------------

      MODULE ALL_VARIABLES
 
      IMPLICIT NONE

      integer, parameter :: SP=SELECTED_REAL_KIND(6,30)  !  single precis.
      real*4, parameter  :: rg = 9806., hg=2.**100

      character(len=90) :: pthfin, pthfout, &
                           pthTOPO
      character(len=30) :: fnstin, fnstout, &
                           fnmTOPO
      character(len=130) :: fina, finb, &
                            fouta, foutb, &
                            ftopa, ftopb
      character(len=70) :: cline
      character(len=15) :: param_file
      integer :: IDM, JDM, IJDM, ios, npad, &
                 nrecL1

!      real(SP) :: rnan  ! nan value that needs to be replaced 
      real(SP), allocatable :: fin1d(:), fout1d(:), &
                               fin2d(:,:), fout2d(:,:), &
                               FinT(:), FoutT(:), HH(:,:)


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
      read(cline(ll+2:),'(A)') pthfin
      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),'(A)') pthfout
      print*,'cline:  ',cline
!      print*,'ll=',ll
!      print*,'cline:  ',cline(ll+2:)
!      print*,'pthfout=',pthfout
      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),'(A)') fnstin
      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),'(A)') fnstout
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

      fina  = trim(pthfin)//trim(fnstin)//'.a'
      finb  = trim(pthfin)//trim(fnstin)//'.b'
      fouta = trim(pthfout)//trim(fnstout)//'.a'
      foutb = trim(pthfout)//trim(fnstout)//'.b'
      ftopa = trim(pthTOPO)//trim(fnmTOPO)//'.a'
      ftopb = trim(pthTOPO)//trim(fnmTOPO)//'.b'

      print*,' Domain Dim: i=',IDM,' j=',JDM
      print*,'Input files : ',trim(fina)
      print*,'Input files : ',trim(finb)
      print*,'Output files: ',trim(fouta)
      print*,'Output files: ',trim(foutb)
      print*,'Topo file:    ',trim(ftopa)
      print*,'Topo file:    ',trim(ftopb)

!      print*,'Replaced values: ',rnan

      IJDM = IDM*JDM
! Define padding HYCOM file GLBb:
      npad = 4096-mod(IJDM,4096)

! Arrays with no padding 1D and 2D
      allocate(fin2d(IDM,JDM),fout2d(IDM,JDM),HH(IDM,JDM))
      allocate(fin1d(IJDM),fout1d(IJDM))
! Arrays for reading/writing full record with npad
! for direct-access file
      allocate(FinT(IJDM+npad),FoutT(IJDM+npad))
      inquire (iolength=nrecL1) FinT

      print*,' ARCc: size 1 rec=',size(fout1d),' npad ARC=',npad
      print*,' Record lengths, ARCc=',nrecL1

      close(11)


      END SUBROUTINE READ_PARAM
!
      SUBROUTINE READ_TOPO

      real(SP) :: dmm

      print*,'READ_TOPO: ', trim(ftopa)

      open(11,file=trim(ftopa),&
              action='read',form='unformatted',&
              access='direct',recl=nrecL1, iostat=ios)
      if (ios.ne.0) call ERR_STOP(ftopa,1,ios)
      
      read(11, rec=1, iostat=ios) FinT
      if (ios.ne.0) call ERR_STOP(ftopa,2,ios)
      fin1d = FinT(1:IJDM)
      HH = reshape(fin1d,(/IDM,JDM/))
      close(11)

      END SUBROUTINE READ_TOPO
!
!
      SUBROUTINE FILL_NANS
! Fill nans/missing values near OBs
      REAL*4  :: dmm1, dmm2, derr
      INTEGER :: i,j, kk, cnt, cnan
      INTEGER :: ich, jch, imm, k, &
                 irec, irec2, &
                 llist, lstr, chck, dii, djj, &
                 i0, j0, np, iM, mm
      INTEGER :: min1(1)
      INTEGER, allocatable :: Iindx(:), Jindx(:), DST(:)
      INTEGER, allocatable :: Inan(:), Jnan(:)

!      print*,'==== Subroutine: Fill_NaNs: 1'
      cnt = 0
      cnan = 0
      DO i=1,IDM
      DO j=1,JDM
        if (fin2d(i,j) .lt. 1.e20 .and. &
            fin2d(i,j) .gt. -100. .and. &
            HH(i,j) .lt. 1.e20) cnt = cnt+1
        if (fin2d(i,j) < -100.) cnan=cnan+1
      ENDDO
      ENDDO

      if (allocated(Iindx)) deallocate(Iindx)
      if (allocated(Jindx)) deallocate(Jindx)
      if (allocated(DST)) deallocate(DST)
      if (allocated(Inan)) deallocate(Inan)
      if (allocated(Jnan)) deallocate(Jnan)
      allocate(Iindx(cnt), Jindx(cnt), DST(cnt))
      allocate(Inan(cnan), Jnan(cnan))

!      print*,'===== Subroutine: Fill_NaNs: 2'
      kk = 0
      mm = 0
      DO i=1,IDM
        DO j=1,JDM
          if (fin2d(i,j) .lt. 1.e20 .and. &
              fin2d(i,j) .gt. -100. .and. &
              HH(i,j) .lt. 1.e20) then
            kk = kk+1
            Iindx(kk) = i
            Jindx(kk) = j
          endif
          if (fin2d(i,j) < -100.) then
            mm = mm+1
            Inan(mm) = i
            Jnan(mm) = j
          endif 
        ENDDO
      ENDDO

      if (kk .ne. cnt) then
        print*,'ERR: mismatch Number of non-nans ...',cnt, kk
        STOP
      endif
      if (mm .ne. cnan) then
        print*,'ERR: mismatch Number of NaNs ...',cnan, mm
        STOP
      endif

!      print*,'==== Subroutine: Fill_NaNs: 3'
      print*,'Found ',cnan,' bad points to be replaced'
!  Fill -999 (nans) at the northern OB
! with closest not-nan
! some values are half of -999 - bilenear interpolated
      fout2d = fin2d

      if (cnan > 0) then
        DO kk=1,cnan
!          if (mod(kk,1000) == 0) &
!          print*,'----  kk =  ',kk

          i = Inan(kk)
          j = Jnan(kk)
          dmm1 = fout2d(i,j)
          if (dmm1>-100.) cycle
          if (HH(i,j) > 1.e20) then
            fout2d(i,j) = hg
            cycle
          endif
!
! Find closest point: 
          DST=sqrt(float((Iindx-i)**2+(Jindx-j)**2))
          min1 = minloc(DST)
          iM = min1(1)
          i0 = Iindx(iM)
          j0 = Jindx(iM)
          derr = DST(iM)
          if (derr > 10) then
            print*,'kk=',kk,' Identified closest point is too far away:'
            print*,'Pnt i=',i,' j=',j
            print*,'Closest point: i0=',i0,' j0=',j0
!            pause
            STOP 
          endif
          fout2d(i,j) = fin2d(i0,j0)
        ENDDO
      endif
 !     print*,'==== Subroutine: Fill_NaNs: 4 DONE'

      END SUBROUTINE FILL_NANS


      SUBROUTINE ERR_STOP(fnm,fmode,ios)

      character(*), intent(in) :: fnm
      integer, intent(in) :: fmode, ios
! 1 - open
! 2 - reading
! 3 - writing
      if (fmode==1) then
        write(*,'(2A)'),'    *** ERROR opening ',trim(fnm)
      elseif (fmode==2) then
        if (ios>0) then 
          write(*,'(2A)'),'    *** ERROR reading: check input ',trim(fnm)
        elseif (ios<0) then
          write(*,'(2A)'),'    *** ERROR reading: E-o-F ',trim(fnm)
        endif
      else
        write(*,'(2A)'),'    *** ERROR writing ',trim(fnm)
      endif
      write(*,'(A, I)'),' IOSTAT = ',ios
        
      STOP

      END SUBROUTINE ERR_STOP


      END MODULE ALL_VARIABLES
