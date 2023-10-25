! ------------------------------
! Declare variables, paths, etc.
! ------------------------------

      MODULE INIT_VARIABLES

      IMPLICIT NONE

      integer, parameter :: SP=SELECTED_REAL_KIND(6,30), & !single precis.
                            ichck = 800, jchck = 10   ! check points
      real*4, parameter  :: rg = 9806., hg=2.**100

      logical :: lexst
      integer :: i, j, k, ios
      integer :: IDM, JDM, IJDM, & 
                 npad,           &
                 KDM, KDMold,    &
                 nsigma,         &! # of sigma-coord.layers
                 Lx_new,         &!Layer 41-lr grid t.dens matches Lx_old
                 Lx_old,         &!Layer 32-lr grid T.dens matches Lx_new
                 nvars            !# of 3D fields interpolated

!      character(len=*), parameter :: 
      character(len=8)   :: fld, fld0, fld1, fld2, fld3
      character(len=90)  :: ccline, &   ! dummy line for *.b files
                            cheadr      ! info line for *.b
      character(len=40)  :: fnmTOPO, fnmIN, &
                            fnmOUT, fnmRLX
      character(len=100) :: pthTOPO, pthIN, &
                            pthOUT, pthBLKDT, &
                            pthrlx
      character(len=100) :: cheader
      character(len=140) :: fina, finb,   &
                            fouta, foutb, &
                            ftopa,ftopb,  &
                            finRa, finRb, & ! 41lr relax/template files
                            fblk32, fblk41 
 
      real(SP), allocatable :: HH(:,:), & ! ARCc topography T11
                               toto(:), & ! padding for binary file
                               dmm(:),  & ! array record+padding read in
                               dmm1d(:), & ! 1D for field
                               dp0k(:), & ! deep z-lev spac. min thck(m)
                               ds0k(:), & ! shall. z-lev spac. min thck(m)
                               TDENS(:),& ! target dens. for 41lr
                                 zmo(:),& !1D depth, 32lr
                                 zmn(:),& !1D depth, 41lr
                             ZZo(:,:,:),& !Old(32+1) vert intrf depths (<0)
                            ZZin(:,:,:),& !New(41+1) example, template
                             ZZn(:,:,:),& !New(41+1) vert intrf depths
                           DPold(:,:,:),& !Old layer thickness (m)
                           DPnew(:,:,:),& !New lr thickn. from template
                            Fold(:,:,:),& !dummy array to read 32lr
                            Fnew(:,:,:),& !interpolated field 41lr
                      Fint_ALL(:,:,:,:)   ! interp.41lr, for writing
                              
      character(len=8), allocatable :: VARS(:)

! ----------------------
      contains
! ----------------------
      SUBROUTINE READ_PARAM
! Read in input from PARAM.dat
      character :: cdm1*11,cdm2*2, cdm3*8, cdm4*30
      integer :: ix, lx, lxc

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
      DO i=1,10 
        read(11,'(A)') cheader
        print*,trim(cheader)
      ENDDO

      read(11,'(3(A))') cdm1, cdm2, pthTOPO
      read(11,'(3(A))') cdm1, cdm2, pthIN
      read(11,'(3(A))') cdm1, cdm2, pthOUT
      read(11,'(3(A))') cdm1, cdm2, pthBLKDT
      read(11,'(3(A))') cdm1, cdm2, pthRLX
      read(11,'(3(A))') cdm1, cdm2, fnmTOPO
      read(11,'(3(A))') cdm1, cdm2, fnmIN
      read(11,'(3(A))') cdm1, cdm2, fnmOUT
      read(11,'(3(A))') cdm1, cdm2, fnmRLX
      read(11,'(2(A),I)') cdm1, cdm2, IDM
      read(11,'(2(A),I)') cdm1, cdm2, JDM
      read(11,'(2(A),I)') cdm1, cdm2, KDM
      read(11,'(2(A),I)') cdm1, cdm2, KDMold
      read(11,'(2(A),I)') cdm1, cdm2, Lx_old
      read(11,'(2(A),I)') cdm1, cdm2, Lx_new
      read(11,'(2(A),I)') cdm1, cdm2, nvars
      print*,'Reading vars =',nvars
      allocate(VARS(nvars))
      DO k=1,nvars
        read(11,*,iostat=ios) cdm1,cdm2,cdm3
        if (ios>0) call ERR_STOP('PARAM.dat',2,ios)
        if (ios<0) exit
        print*,'cdm1= ',cdm1,' cdm2= ',cdm2,' cdm3= ',cdm3
!        ix = index(cdm4,"vars")
!        if (ix==0) exit
!        ix = index(cdm4,"=")
!        lx = len_trim(cdm4)
!        lxc = lx-(ix+1)
        lxc = len_trim(cdm3)
        print*,'cdm3=',cdm3
        VARS(k)=trim(cdm3)
      ENDDO
      read(11,'(3(A),I)') cdm1,cdm2,cheadr

 
      ftopa = trim(pthTOPO)//trim(fnmTOPO)//'.a' ! topo *a
      ftopb = trim(pthTOPO)//trim(fnmTOPO)//'.b' ! topo *b
      fina  = trim(pthIN)//trim(fnmIN)//'.a'     ! 32-lr *a file
      finb  = trim(pthIN)//trim(fnmIN)//'.b'
      fouta = trim(pthOUT)//trim(fnmOUT)//'.a'   ! 41-lr *a file
      foutb = trim(pthOUT)//trim(fnmOUT)//'.b'
      finRa = trim(pthRLX)//trim(fnmRLX)//'.a'   ! 41-lr relax *a file
      finRb = trim(pthRLX)//trim(fnmRLX)//'.b'
      fblk32= trim(pthBLKDT)//'blkdat.input_ARCc0.08_32lev'
      fblk41= trim(pthBLKDT)//'blkdat.input_ARCc0.08_41lev'


! Dimensions of the HYCOM binary file 
      IJDM = IDM*JDM
      npad = 4096-mod(IJDM,4096)
!      nrecl0 = (IJDM+npad)*4
      

! Allocate arrays
      allocate(toto(npad), dmm(IJDM+npad), dmm1d(IJDM))
      allocate(HH(IDM,JDM))
      allocate(ZZo(KDMold+1,IDM,JDM), &
               ZZin(KDM+1,IDM,JDM), &
               ZZn(KDM+1,IDM,JDM), &
               DPnew(KDM,IDM,JDM), &
               Fnew(KDM,IDM,JDM), &
               Fint_ALL(nvars,KDM,IDM,JDM),&
               DPold(KDMold,IDM,JDM), &
               Fold(KDMold,IDM,JDM))
      allocate(zmo(KDMold),zmn(KDM))
!
      toto = 9999.99
!
      print*,' ARCc Dim:    i=',IDM,' j=',JDM
      print*,' Interpolating from k=',KDMold,' -> k=',KDM
      write(*,'(2A)'),'TOPO files:',trim(ftopa)
      write(*,'(2A)'),'TOPO files:',trim(ftopb)
      write(*,'(2A)'),'Old input files:',trim(fina)
      write(*,'(2A)'),'Old input files:',trim(finb)
      write(*,'(2A)'),'New output files:',trim(fouta)
      write(*,'(2A)'),'New output files:',trim(foutb)
      write(*,'(2A)'),'Relax template files:',trim(finRa)
      write(*,'(2A)'),'Relax template files:',trim(finRb)
      write(*,'(2A)'),'BLKDAT.INPUT 32lrs:',trim(fblk32)
      write(*,'(2A)'),'BLKDAT.INPUT 41lrs:',trim(fblk41)
      write(*,'(2(A,I))'),'npad = ',npad,' IJDM=',IJDM
      write(*,'(A,I)'),'size dmm=',size(dmm)
      write(*,'(A,I)'),'Matching t.dens, Laer 32grid=',Lx_old
      write(*,'(A,I)'),'Matching t.dens, Laer 41grid=',Lx_new
      write(*,'(A)'),'Variables to Interpolate :'
      write(*,*),VARS
      write(*,'(A)'),trim(cheadr)
      write(*,*),' ================================'
      write(*,*),'  '
   
      close(11)
!      pause
      END SUBROUTINE READ_PARAM


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



      


      END MODULE INIT_VARIABLES
