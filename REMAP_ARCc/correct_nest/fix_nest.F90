      PROGRAM FIX_NEST
! Quick fix of the U,V, barotrop U errors 
! in the nest fields
! Original remapping GLBb -> ARCc0.08 didn't
! rotate U,V components in the upper part of the ARCc
! In order not to regenirate all nest fields
! This quick bug fix corrects the existing nest fields

      USE variables
      USE utils

      IMPLICIT NONE

      integer :: nrecl, nl, ix, jx, ixc, &
                 nn, mm, irec, irec2, &
                 cnt, imm, jmm, kk, nlr

      real :: bmin, bmax
      real(SP) :: dens, td0, mday

      character(60) :: fmat1 ! *b IO format
      character(8)  :: smb*1, vname
      character(80) :: str1, str2, ch1

      call read_param

! Read remapping indices
      print*,'Reading Remap indices '
      call GET_REMAP_INDX

! Read/Write fields
      inquire (iolength=nrecl) FinT
! Input files:
      open (unit=20, file=trim(fina), &
            action='read', form='unformatted', &
            access='direct', recl=nrecl, iostat=ios)
      if (ios.ne.0) call ERR_STOP(fina,1,ios)

      open(unit=21, file=trim(finb), &
           action='read', status='old', &
           form='formatted',iostat=ios)
      if (ios.ne.0) call ERR_STOP(finb,1,ios)

! Output files:
      open(unit=40, file=trim(fouta), &
           action='write',form='unformatted',&
           access='direct',recl=nrecl, iostat=ios)
      if (ios.ne.0) call ERR_STOP(fouta,1,ios)

      open(unit=41, file=trim(foutb), &
           action='write',form='formatted',&
           iostat=ios)
      if (ios.ne.0) call ERR_STOP(foutb,1,ios)


! Read the header lines:
      print*,'Reading *b:  '
!      print*,trim(finb)
      DO nl=1,1000
        read(21,'(A)',iostat=ios) ccline
        if (ios.ne.0) call ERR_STOP(finb,2,ios)
        write(*,'(A)'), trim(ccline)
        ixc = index(ccline,"time step")
        ix = index(ccline,"idm ")
        jx = index(ccline,"jdm ")
        if (ix>0) read(ccline,*) nn,ch1
        if (jx>0) read(ccline,*) mm,ch1
        if (nl == 4) ccline = trim(chdr_write)

        write(41,'(A)') trim(ccline)
        if (ixc>0) exit

      ENDDO

      print*,'header done   '

      if (nn.ne.IDM .or. mm.ne.JDM) then
        print*,'Size of F idm=',IDM, ' jdm=',JDM
        print*,'Check idm ',nn,' & jdm ',mm,' from ',finb
        STOP
      endif
      
! Read/write fields by layers
      fmat1='(A8,1x,A1,2x,I9,1x,f10.3,1x,I2,1x,f5.2,2x,se14.7,2x,se14.7)'

! Write all fields following *.b 
! all U/V fields - correct
      irec  = 0
      irec2 = 0
      cnt   = 0
      DO
        cnt = cnt+1
        read(21, trim(fmat1), iostat=ios) &
             vname,smb,imm,mday,kk,dens,bmin,bmax

        if (ios < 0) exit ! EOF
        if (ios > 0) call ERR_STOP(finb,2,ios)
        if (cnt > 1000) call ERR_STOP(finb,2,cnt)
        if (trim(vname)=='') exit

        print*,cnt,'*b: var=',vname,' mday=',mday,' dens=',dens,&
               ' bmin=',bmin,' bmax=',bmax
! Break after barotropic/surface fields:
!        if (kk>0 .and. trim(vname).ne.'montg1') exit
!        if (trim(vname) == 'u-vel.') then
!          print*,'ERR: fix_nest reading *b: '
!          print*,'unexpected 3D field found in *.b:',vname
!          STOP
!        endif

        irec=irec+1
        read(20, rec=irec, iostat=ios) FinT
        if (ios.ne.0) call ERR_STOP(fina,2,ios)

        if ( vname(2:5) .eq. '_btr' .or. &
             vname(2:5) .eq. '-vel') then
           print*,'Fixing vel fields ...'
!           print*,'Max val=',maxval(FinT(1:IJDM))
!           print*,'FinT(1)=',FinT(1)
            call FIX_ARC_NEST(FinT, bmin, bmax)        
        endif

        write(40, rec=irec, iostat=ios) FinT
        if (ios.ne.0) call ERR_STOP(fouta,3,ios)
        print*,'Write ',vname,bmin,bmax
        write(41,fmat1) vname,smb,imm,mday,kk,dens,bmin,bmax

      ENDDO  !
 
      close(20)
      close(21)
      close(40)
      close(41)


      END PROGRAM FIX_NEST

