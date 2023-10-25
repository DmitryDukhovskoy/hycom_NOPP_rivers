! ------------------------------
! Declare variables, paths, etc.
! ------------------------------

      MODULE ALL_VARIABLES
 
      IMPLICIT NONE

      integer, parameter :: SP=SELECTED_REAL_KIND(6,30)  !  single precis.

!      character(len=*), parameter :: 
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

      logical :: irst ! restart-file flag
! ----------------------
      contains
! ----------------------
      SUBROUTINE READ_PARAM

      integer :: clen, clenL, i1, i2
      character :: cdm1*9,cdm2*2, cline*7

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

! check if restart file
      irst = .false.
      cline= 'restart'
      clenL= len_trim(cline)-1
      clen = len_trim(fnmGLB)
      i1=1
      i2=i1+clenL
      DO WHILE (i2<=clen)
        irst = (fnmGLB(i1:i2).eq.cline)
        if (irst) exit
        i1=i1+1
        i2=i1+clenL
      ENDDO

      if (irst) then
        print*,'    '
        print*,' ====  RESTART FILES ========'
        print*,'    '
      endif



      END SUBROUTINE READ_PARAM

      END MODULE ALL_VARIABLES
