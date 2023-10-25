! ------------------------------
! Declare variables, paths, etc.
! ------------------------------

      MODULE ALL_VARIABLES
 
      IMPLICIT NONE

      integer, parameter :: SP=SELECTED_REAL_KIND(6,30)  !  single precis.
      real*4, parameter :: hg=2**100
!      character(len=*), parameter :: 
      character(len=80) ctitle   ! string to put in *.b
      character(len=90) pthG, pthA, pthGMAP
      character(len=30) fnmGLB, fnmARC, fnmGMAP
      character(len=120) fina, finb, fouta, foutb, &
                         finGMa, finGMb
      integer ::  IDMG,JDMG,IJDMG, &
                  IDMA,JDMA,IJDMA, &
                  ios 


      real(SP), allocatable :: xmap(:,:), ymap(:,:)
!      real(SP), allocatable :: fglb(:,:), farc(:,:)
!      real(SP), allocatable :: xmap(:), ymap(:)
      real(SP), allocatable :: fglb(:), farc(:)

! ----------------------
      contains
! ----------------------
      SUBROUTINE READ_PARAM

      character :: cdm1*9,cdm2*2

      print*,'Start reading PARAM.dat'
!      open(11,file='dd.dat',action='write')
!      write(11,*) pthG
!      close(11)
      open(unit=11, file='PARAM.dat', &
               action='read', status='old', &
               form='formatted',iostat=ios)

      if (ios>0) then
        print*,'    *** ERR:  ERROR opening PARAM.dat'
        STOP
      endif

      read(11,'(3(A))') cdm1, cdm2, pthG
      read(11,'(3(A))') cdm1, cdm2, pthA
      read(11,'(3(A))') cdm1, cdm2, fnmGLB
      read(11,'(3(A))') cdm1, cdm2, fnmARC
      read(11,'(3(A))') cdm1, cdm2, pthGMAP
      read(11,'(3(A))') cdm1, cdm2, fnmGMAP
      read(11,'(3(A))') cdm1, cdm2, ctitle
      read(11,'(2(A),I)') cdm1, cdm2, IDMG
      read(11,'(2(A),I)') cdm1, cdm2, JDMG
      read(11,'(2(A),I)') cdm1, cdm2, IDMA
      read(11,'(2(A),I)') cdm1, cdm2, JDMA

      fina  = trim(pthG)//trim(fnmGLB)//'.a'
      finb  = trim(pthG)//trim(fnmGLB)//'.b'
      fouta = trim(pthA)//trim(fnmARC)//'.a'
      foutb = trim(pthA)//trim(fnmARC)//'.b'
      finGMa= trim(pthGMAP)//trim(fnmGMAP)//'.a'
      finGMb= trim(pthGMAP)//trim(fnmGMAP)//'.b'

      print*,' Global Dim: i=',IDMG,' j=',JDMG
      print*,' ARC Dim:    i=',IDMA,' j=',JDMA
      print*,'Global input files:',trim(fina)
      print*,'Global input files:',trim(finb)
      print*,'ARCc output files:',trim(fouta)
      print*,'ARCc output files:',trim(foutb)
      print*,'Remap index files:',trim(finGMa)
      print*,'Remap index files:',trim(finGMb)

      IJDMG = IDMG*JDMG
      IJDMA = IDMA*JDMA

      allocate(xmap(IDMA,JDMA),ymap(IDMA,JDMA))
!      allocate(fglb(IDMG,JDMG),farc(IDMA,JDMA))
!      allocate(xmap(IJDMA),ymap(IJDMA))
      allocate(fglb(IJDMG),farc(IJDMA))

      close(11)

      END SUBROUTINE READ_PARAM

      END MODULE ALL_VARIABLES
