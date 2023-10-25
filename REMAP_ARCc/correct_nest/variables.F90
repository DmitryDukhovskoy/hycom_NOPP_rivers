! ------------------------------
! Define variables, pathnames
! ------------------------------

      MODULE VARIABLES

      IMPLICIT NONE

      integer, parameter :: SP=SELECTED_REAL_KIND(6,30)
      real*4, parameter  :: rg = 9806., hg=2.**100

      logical :: lexst
      integer :: i, j, k, ios
      integer :: IDM, JDM, IJDM, & ! horizontal dim 
                 npad,           & ! padding
                 KDM               ! vertical lrs

!      character(len=*), parameter :: 
      character(len=8)   :: fld, fld0, fld1, fld2, fld3
      character(len=90)  :: ccline, &   ! dummy line for *.b files
                            chdr_write  ! info line for *.b
      character(len=40)  :: fnmIN, fnmOUT, fnmGMAP
      character(len=100) :: pthIN, pthOUT, pthGMAP
      character(len=100) :: cheader
      character(len=140) :: fina, finb,   &
                            fouta, foutb, &
                            finGMa, finGMb
!
! Allocatables
      real (SP), allocatable :: dmm(:),     & ! HYCOM IO: 1D array record + padding 
                                FinT(:),    & ! read/write 1D fld+padding
                                fin1d(:),   & ! input 1D field w/o padding
                                fin2d(:,:), & ! input 2D
                                fout1d(:),  & ! output 1D field w/o padding
                                fout2d(:,:),& ! output 2D field
                                toto(:),    & ! padding
                         xmap(:,:), ymap(:,:) ! GLBb->ARCc remap indices

! -----------------

      contains

! -----------------

      subroutine read_param
      character :: cdm1*11, cdm2*2, cdm3*8, cdm4*30

      inquire(file='PARAM.dat', exist=lexst)
      if (.not.lexst) then
        print*,'PARAM.dat does not exist'
        STOP
      endif
      print*,'Reading PARAM.dat'
      open(unit=11, file='PARAM.dat', &
               action='read', status='old', &
               form='formatted',iostat=ios)
      if (ios.ne.0) call ERR_STOP('PARAM.dat',1,ios)
!
! Read the header lines:
      print*,'Reading header: '
      DO i=1,10000
        read(11,'(A)',iostat=ios) cheader
        if (ios.ne.0) call ERR_STOP('PARAM.dat',2,ios)
        print*,trim(cheader)
        if (cheader(1:10).eq.'==========') exit
      ENDDO

      read(11,'(3(A))') cdm1, cdm2, pthIN
      read(11,'(3(A))') cdm1, cdm2, pthOUT
      read(11,'(3(A))') cdm1, cdm2, fnmIN
      read(11,'(3(A))') cdm1, cdm2, fnmOUT
      read(11,'(3(A))') cdm1, cdm2, pthGMAP
      read(11,'(3(A))') cdm1, cdm2, fnmGMAP
      read(11,'(2(A),I)') cdm1, cdm2, IDM
      read(11,'(2(A),I)') cdm1, cdm2, JDM
      read(11,'(2(A),I)') cdm1, cdm2, KDM
      read(11,'(3(A))') cdm1, cdm2, chdr_write
 
      close(11)
 
      fina  = trim(pthIN)//trim(fnmIN)//'.a'     ! error file
      finb  = trim(pthIN)//trim(fnmIN)//'.b'
      fouta = trim(pthOUT)//trim(fnmOUT)//'.a'   ! corrected
      foutb = trim(pthOUT)//trim(fnmOUT)//'.b'
      finGMa= trim(pthGMAP)//trim(fnmGMAP)//'.a'
      finGMb= trim(pthGMAP)//trim(fnmGMAP)//'.b'

      print*,'Input file: ',trim(fina)
      print*,'Input file: ',trim(finb)
      print*,'Output file: ',trim(fouta)
      print*,'Output file: ',trim(foutb)
      print*,'Remap index file:',trim(finGMa)
      print*,'Remap index file:',trim(finGMb)

! Dimensions of HYCOM binary:
      IJDM = IDM*JDM
      npad = 4096-mod(IJDM,4096)

! Allocate arrays:
      allocate (toto(npad), dmm(IJDM+npad), &
                FinT(IJDM+npad), fin1d(IJDM), &
                fout1d(IJDM), fin2d(IDM,JDM),&
                fout2d(IDM,JDM), xmap(IDM,JDM), &
                ymap(IDM,JDM))

      end subroutine read_param

      SUBROUTINE ERR_STOP(fnm,fmode,ios)
!
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


      END MODULE VARIABLES
