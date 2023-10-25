      PROGRAM LagrPart
! Track Lagrangian particles
! seeded within some region
! Use Runge-Kutta method
! to solve advection eq.
      USE utils

      IMPLICIT NONE

!      real, allocatable(:) :: lat1, lon1
      real*4 :: lat1, lon1, lat2, lon2, dmm

      CALL READ_PARAM

      CALL READ_TOPO_008
      print*,'Topo read ARCc0.08: HH(500,1250)=',HH(500,1250),&
             'LMSK=',LMSK(500,1250)

!     pause
      CALL READ_LON_LAT

      print*,'LAT(10,11)=',LAT(10,11)
!      CALL wait()
!      print*,'LON=',LON(10,11)
      
!
! Calculate cartesian distances
! convert geogr -> cart. coord
      print*,'Calculating DX, DY ...' 
!      allocate(lat1(JDM), lon1(IDM))
     
      DO i=1,IDM-1
        DO j=1,JDM-1
          lat1=LAT(i,j)
          lat2=LAT(i+1,j)
          lon1=LON(i,j)
          lon2=LON(i+1,j)
!          print*,'lat1=',lat1,' lon1=',lon1
!          print*,'lat2=',lat2,' lon2=',lon2
!          CALL wait()
          CALL DISTANCE_SPHERIC(lat1,lon1,lat2,lon2,dmm)
          DX(i,j)=dmm

          lat2=LAT(i,j+1)
          lon2=LON(i,j+1)
          CALL DISTANCE_SPHERIC(lat1,lon1,lat2,lon2,dmm)
          DY(i,j)=dmm

!          print*,'i=',i,'j=',j,'dx=',DX(i,j),'DY=',DY(i,j)
!          CALL wait()

        ENDDO
      ENDDO

! Initialize
      CALL SEED_INI(PRTCL(1)%alive,   &
                    PRTCL(1)%kzlev,   &
                    PRTCL(1)%time,    &
                    PRTCL(1)%iindx,   &
                    PRTCL(1)%jindx,   &
                    PRTCL(1)%xcrd,    &
                    PRTCL(1)%ycrd,    &
                    PRTCL(1)%Nprt)

      print*,'------------ ALL DONE ---------'

      END PROGRAM LagrPart

