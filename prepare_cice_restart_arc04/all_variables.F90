! --------------------------------
! Declare variables, paths, etc.
! --------------------------------

      MODULE ALL_VARIABLES

      IMPLICIT NONE

      integer, parameter :: SP=SELECTED_REAL_KIND(6,30) !single precis.
      integer, parameter :: DP=SELECTED_REAL_KIND(18,30)! double precis
      integer, parameter :: ldebug = 0 ! <2 - less info

      real*4, parameter :: rg = 9806., hg = 2.**100

      character(len=90) :: pthfin, pthfout, &
                           pthTOPO,   &
                           pthTOPO8
      character(len=30) :: fnstin, fnstout, &
                           fnmTOPO, fnmTOPO8
      character(len=130) :: fina, finb, &
                            fouta, foutb, &
                            ftopa, ftopb, &
                            ftop8a, ftop8b

      character(len=70) :: cline
      character(len=15) :: param_file

! All fields for ARCc0.08 have "in" - those that need to be interpolated
! ARCc0.04 fields use  "out"
      integer :: i, j, k, ios
      integer :: IDMin, JDMin, IJDMin, NPADin, &
                 nrecLin ! ARCc0.08
      integer :: IDMout, JDMout, IJDMout,  NPADout, &
                 nrecLout ! ARCc0.04
      integer :: ncat,  & ! # of ice categories
                 nilyr, & ! # of ice layers in each cat
                 nslyr, & ! # of snow layers
                 ntilyr,& ! # of total ice layers for all cat
                 ntslyr,& ! # of total snow layers for all cat
                 ntrcr, & ! # of tracers
                 nfctr    ! factor new grid resolution/old grid - integer

      integer, allocatable :: LMSKin(:,:), LMSKout(:,:)

      real*8 :: rdayN ! new run day (>0 if need to replace) HYCOM conven. 1901
      real(SP), allocatable :: Hin(:,:), & ! HYCOM topo ARCc0.08
                               Hout(:,:), &
                               fin1d(:), &
                               fout1d(:), &
                               FinT(:), FoutT(:)

      real*8, allocatable :: fin2d(:,:), fout2d(:,:), aicen(:,:,:)
                
! -------------------------
      CONTAINS
! -------------------------
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

      if (ios>0) call err_stop(param_file,1,ios)

      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),'(A)') pthfin

      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),'(A)') pthfout
!      print*,'cline:  ',cline

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
      read(cline(ll+2:),*) rdayN

      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),*) ncat

      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),*) nilyr

      read(11,'(A)') cline
      ll = index(cline,'=') 
      read(cline(ll+2:),*) nslyr

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
      read(cline(ll+2:),*) nfctr  ! resolutin factor


      fina  = trim(pthfin)//trim(fnstin)
      fouta = trim(pthfout)//trim(fnstout)
      ftopa = trim(pthTOPO)//trim(fnmTOPO)//'.a'
      ftopb = trim(pthTOPO)//trim(fnmTOPO)//'.b'
      ftop8a= trim(pthTOPO8)//trim(fnmTOPO8)//'.a'
      ftop8b= trim(pthTOPO8)//trim(fnmTOPO8)//'.b'

      ntilyr = ncat*nilyr ! # of ice layers in all categories
      ntslyr = ncat*nslyr ! # of snow layers in all cat.
      ntrcr  = 3          ! # of tracers, 1 - T surf

      print*,' Domain Dim: i=',IDMout,' j=',JDMout
      print*,'Input files : ',trim(fina)
      print*,'Output files: ',trim(fouta)
      print*,'Topo file 0.04:    ',trim(ftopa)
      print*,'Topo file 0.04:    ',trim(ftopb)
      print*,'Topo file 0.08:    ',trim(ftop8a)
      print*,'Topo file 0.08:    ',trim(ftop8b)
      print*,'# of ice cat: ',ncat
      print*,'# of ice layers: ',nilyr
      print*,'# of snow layers: ',nslyr
      if (rdayN > 0) then
        print*,'Change restart day: ',rdayN
      else
        print*,'Restart day wont be changed'
      endif

      IJDMin  = IDMin*JDMin
      IJDMout = IDMout*JDMout
! Define padding HYCOM file GLBb:
      NPADin  = 4096-mod(IJDMin,4096)
      NPADout = 4096-mod(IJDMout,4096)

! Arrays with no padding 1D and 2D
      allocate(fin2d(IDMin,JDMin),fout2d(IDMout,JDMout),Hout(IDMout,JDMout))
      allocate(fin1d(IJDMin),fout1d(IJDMout))
      allocate(Hin(IDMin,JDMin))
      allocate(LMSKin(IDMin,JDMin), LMSKout(IDMout,JDMout))
      allocate(aicen(ncat,IDMout,JDMout))
! Arrays for reading/writing full record with npad
! for direct-access file
      allocate(FinT(IJDMin+NPADin),FoutT(IJDMout+NPADout))
      inquire (iolength=nrecLin) FinT
      inquire (iolength=nrecLout) FoutT

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
!      
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

