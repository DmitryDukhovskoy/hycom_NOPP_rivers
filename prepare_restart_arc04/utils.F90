! --------------------------------------------
! Read Param, prepare filenames, etc
! --------------------------------------------

      MODULE UTILS

      IMPLICIT NONE

      integer, parameter :: SP=SELECTED_REAL_KIND(6,30)  !  single precis.
      integer, parameter :: ldebug = 0    ! <2 - low output info
      real*4, parameter  :: rg = 9806., hg=2.**100

      character(len=90) :: pthfin, pthfout, &
                           pthTOPO, pthf08, &
                           pthTOPO8
      character(len=30) :: frst08, fnstin, fnstout, &
                           fnmTOPO, fnmTOPO8
      character(len=130) :: fina08, finb08, &
                            fina04, finb04, &
                            fouta, foutb, &
                            ftopa, ftopb, &
                            ftop8a, ftop8b

      character(len=72) :: cline
      character(len=15) :: param_file
! All fields for ARCc0.08 have "in" - those that need to be interpolated
! ARCc0.04 fields use  "out"
      integer :: i, j, k
      integer :: IDMin, JDMin, IJDMin, ios, NPADin, &
                 nrecLin ! ARCc0.08
      integer :: IDMout, JDMout, IJDMout,  NPADout, &
                 nrecLout ! ARCc0.04
      integer :: nlyrs, nTr, & ! # v. layers , # of passive tracers
                 nfctr         ! factor new grid resolution/old grid - integer
      integer, allocatable :: LMSKin(:,:), LMSKout(:,:)


      real(SP) :: rstday  ! restart day # to be placed in new *.b
                          ! if =0 - not changed
      real(SP), allocatable :: fin1d(:), fout1d(:), &
                               fin2d(:,:), fout2d(:,:), &
                               FinT(:), FoutT(:), & 
                               Hin(:,:), Hout(:,:)

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
      read(cline(ll+2:),'(A)') pthf08

      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),'(A)') pthfin

      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),'(A)') pthfout
      print*,'cline:  ',cline

      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),'(A)') frst08

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
      read(cline(ll+2:),'(A)') pthTOPO8

      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),'(A)') fnmTOPO8

      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),*) rstday

      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),*) IDMin

      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),*) JDMin

      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),*) IDMout

      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),*) JDMout

      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),*) nlyrs  ! # of vert. levels 

      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),*) nTr  ! # of tracers

      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),*) nfctr  ! resolutin factor

      fina08= trim(pthf08)//trim(frst08)//'.a'
      finb08= trim(pthf08)//trim(frst08)//'.b'
      fina04= trim(pthfin)//trim(fnstin)//'.a'
      finb04= trim(pthfin)//trim(fnstin)//'.b'
      fouta = trim(pthfout)//trim(fnstout)//'.a'
      foutb = trim(pthfout)//trim(fnstout)//'.b'
      ftopa = trim(pthTOPO)//trim(fnmTOPO)//'.a'
      ftopb = trim(pthTOPO)//trim(fnmTOPO)//'.b'
      ftop8a= trim(pthTOPO8)//trim(fnmTOPO8)//'.a'
      ftop8b= trim(pthTOPO8)//trim(fnmTOPO8)//'.b'

      print*,' Domain Dim 008: i=',IDMin,' j=',JDMin
      print*,' Domain Dim 004: i=',IDMout,' j=',JDMout
      print*,'Input files 008 : ',trim(fina08)
      print*,'Input files 008 : ',trim(finb08)
      print*,'Input files 004 : ',trim(fina04)
      print*,'Input files 004 : ',trim(finb04)
      print*,'Output files: ',trim(fouta)
      print*,'Output files: ',trim(foutb)
      print*,'Topo file 0.04:    ',trim(ftopa)
      print*,'Topo file 0.04:    ',trim(ftopb)
      print*,'Topo file 0.08:    ',trim(ftop8a)
      print*,'Topo file 0.08:    ',trim(ftop8b)

      IJDMin  = IDMin*JDMin
      IJDMout = IDMout*JDMout
! Define padding HYCOM file GLBb:
      NPADin  = 4096-mod(IJDMin,4096)
      NPADout = 4096-mod(IJDMout,4096)

! Arrays with no padding 1D and 2D
      allocate(fin2d(IDMin,JDMin),fout2d(IDMout,JDMout),Hout(IDMout,JDMout))
      allocate(fin1d(IJDMin),fout1d(IJDMout), Hin(IDMin,JDMin))
      allocate(LMSKin(IDMin,JDMin), LMSKout(IDMout,JDMout))
! Arrays for reading/writing full record with npad
! for direct-access file
      allocate(FinT(IJDMin+NPADin),FoutT(IJDMout+NPADout))
      inquire (iolength=nrecLin) FinT
      inquire (iolength=nrecLout) FoutT

      print*,' ARCc: output rec = ',size(fout1d),' npad=',NPADout
      print*,' Output Record lengths = ',nrecLout
      print*,' # of vertical levels in restart = ',nlyrs
      print*,' # of Passive Tracers = ',nTr
      print*,'          '

      close(11)

      END SUBROUTINE READ_PARAM
!
!
      SUBROUTINE READ_TOPO_004

      real(SP) :: dmm

      print*,'READ_TOPO: ', trim(ftopa)

      open(11,file=trim(ftopa),&
              action='read',form='unformatted',&
              access='direct',recl=nrecLout, iostat=ios)
      if (ios.ne.0) call ERR_STOP(ftopa,1,ios)
      
      read(11, rec=1, iostat=ios) FoutT
      if (ios.ne.0) call ERR_STOP(ftopa,2,ios)
      fout1d = FoutT(1:IJDMout)
      Hout = reshape(fout1d,(/IDMout,JDMout/))
      close(11)
! Create Land Mask
      DO i=1,IDMout
      DO j=1,JDMout
        if (Hout(i,j)<1.e20) then
          LMSKout(i,j)=1
        else
          LMSKout(i,j)=0
        endif
      ENDDO
      ENDDO

      END SUBROUTINE READ_TOPO_004
!
      SUBROUTINE READ_TOPO_008

      real(SP) :: dmm

      print*,'READ_TOPO: ', trim(ftop8a)

      open(11,file=trim(ftop8a),&
              action='read',form='unformatted',&
              access='direct',recl=nrecLin, iostat=ios)
      if (ios.ne.0) call ERR_STOP(ftop8a,1,ios)
      
      read(11, rec=1, iostat=ios) FinT
      if (ios.ne.0) call ERR_STOP(ftop8a,2,ios)
      fin1d = FinT(1:IJDMin)
      Hin = reshape(fin1d,(/IDMin,JDMin/))
      close(11)
!
! Create Land Mask
      DO i=1,IDMin
      DO j=1,JDMin
        if (Hin(i,j)<1.e20) then
          LMSKin(i,j)=1
        else
          LMSKin(i,j)=0
        endif
      ENDDO
      ENDDO

      END SUBROUTINE READ_TOPO_008
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



      END MODULE UTILS

