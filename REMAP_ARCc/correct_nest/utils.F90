      MODULE UTILS


      USE variables

      IMPLICIT NONE

! -----------------
      contains
! -----------------
      subroutine GET_REMAP_INDX

      integer :: ntoto, nrec, nrecl, npad
      integer :: i, j
      real(SP) :: dmm1d(IJDM)
 
      print*,'GET_REMAP_INDX: ', trim(finGMa)

      inquire (iolength=nrecl) dmm
      print*,'Determined Record length = ',nrecl

      open(11,file=trim(finGMa),&
              action='read',form='unformatted',&
              access='direct',recl=nrecl, iostat=ios)
      if (ios>0) then
        print*,'    *** ERR:  ERROR opening ',trim(finGMa)
        STOP
      endif

      read(11, rec=1, iostat=ios) dmm
      if (ios<0) STOP('READING HIT EOF UNIT=11 *.a');
      if (ios>0) STOP('READING ERROR UNIT=11 *.a')
      dmm1d = dmm(1:IJDM)
      xmap = reshape(dmm1d,(/IDM,JDM/))
!      print*,'xmap =',xmap(1:5,1)

      read(11,rec=2, iostat=ios) dmm
      if (ios<0) STOP('READING HIT EOF UNIT=11 *.a');
      if (ios>0) STOP('READING ERROR UNIT=11 *.a')
      dmm1d = dmm(1:IJDM)
      ymap = reshape(dmm1d,(/IDM,JDM/))

      print*,' Remap indices: GLBb -> ARCc'
      DO i=1,2
      DO j=1,3
        print*,'i=',i,' j=',j,'xmap =',xmap(i,j),&
               &' ymap =',ymap(i,j)
      ENDDO
      ENDDO
      print*,' ----------------   '

      close(11)      

      end subroutine GET_REMAP_INDX
!
!
      subroutine FIX_ARC_NEST(FinT, bmin, bmax)
! Correct velocity fields in an ARCc nest file 
! in the upper part of the grid
! Rotate by 180 degrees
! Double check is rotation is needed to avoid
! accidental correction of corrected files
! at the cut-line of GLBb U and V are opposite if not rotated
!  rotate U and V
! transition occurs at j=1249
! Use GLBb->ARCc remap indices (xmap, ymap):
! all i indices are the same (in columns) for j=1:1248
! and then they switch  to another value for 1249:end
! do not change land values 
      real(SP), intent(inout) :: FinT(:)
      real(SP), intent(out) :: bmin, bmax
      real :: rsum1, rsum2, d1 , d2

      integer :: i, j, jL, cnt, iG, jG

! get data w/o padding:
      fin1d = FinT(1:IJDM)
      fin2d = reshape(fin1d,(/IDM,JDM/))
      FinT  = hg
      i = 1
      jL = 0
      do j=2,JDM
        iG = int(xmap(i,j))
        if ( abs(iG-int(xmap(i,1)))>1 ) then
          jL = j
          exit
        endif
      enddo

!      print*,'GLBb cut line j=',jL
!
! Check if this is an error file      
! Land can be either huge number or 0 
      cnt = 0
      rsum1 = 0
      do i=1,IDM
        d1 = fin2d(i,jL-1)
        d2 = fin2d(i,jL)
        if ( d1<1.e20 .and. d2<1.e20 &
             .and. abs(d1)>1.e-5 .and. abs(d2)>1.e-5) then
          cnt = cnt+1
          d1 = d1/abs(d1)
          d2 = d2/abs(d2)
          rsum1 = rsum1+abs(d1-d2)
!          print*,'i=',i,'  d1=',d1,'  d2=',d2
!          print*,'    rsum1=',rsum1
        endif
      enddo
      rsum1 = rsum1/float(cnt)
!      print*,'rsum1=',rsum1, 'cnt=',cnt

      if (rsum1<0.5) then
        print*,'Existing U/V fields seem to be consisten '
        print*,'across the cut-line j=',jL
        print*,'Double check if the corrected fields need any fix of U,V'
        print*,'Stat of mean(abs(sign(u1)-sugn(u2))) across the line j=',jL
        print*,'Should be close to 2 for screwed up fields and 0 for correct'
        print*,'Here: rsum1=',rsum1, '# of points : cnt=',cnt
        print*,'Stopping ...'
        STOP
      else
        print*,'Checked: correction is needed'
        print*,'==>  Stat of mean(abs(sign(u1)-sugn(u2))) across the line j=',jL
        print*,'==>  Should be close to 2 for error fields and 0 for correct'
        print*,'==>  Here: statistics=',rsum1, '# of points : cnt=',cnt
      endif

!  rotate U and V
! transition occurs at j=jL
! Use GLBb remap indices to find where U/V need to be rotated
! all i indices are the same (in columns) for j=1:1248
! and then they switch  to another value for 1249:end
! do not change land values 
      DO i=1,IDM
      DO j=1,JDM
!        iG=int(xmap(i,j)) ! remap indices from GLBb
!        jG=int(ymap(i,j))
!        if ( abs(iG-int(xmap(i,1)))>1 .and. &
 !             fin2d(i,j)<0.1*hg ) then 
        if ( j>=jL .and. abs(fin2d(i,j))<0.01*hg ) then
          fout2d(i,j) = -fin2d(i,j)
        else
          fout2d(i,j) = fin2d(i,j)
        endif
      ENDDO
      ENDDO

      fout1d = reshape(fout2d,(/IJDM/))
      FinT(1:IJDM) = fout1d

      bmin=minval(fout1d, mask = fout1d .lt. 1.e25)
      bmax=maxval(fout1d, mask = fout1d .lt. 1.e25)

!      pause

      end subroutine FIX_ARC_NEST

      END MODULE UTILS      
