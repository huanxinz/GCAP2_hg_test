module geom_giss_mod
  
  USE CMN_SIZE_MOD

  implicit none
  PRIVATE

  PUBLIC :: GEOM_B
  PUBLIC :: POLWT, KMAXJ, IDIJ, IDJJ, RAPJ, COSIV, SINIV
  PUBLIC :: bySN, DRAT, NMIN
  PUBLIC :: DXP, DYP, BYDYP, DXV, DXYP, IMAXJ !, ZATMO
  PUBLIC :: REL_AREA, GEOFAC, GEOFAC_PC, XCOLMASS_FIX

  ! Dimensions
#if defined( GRIDM23 )
  integer, parameter :: IM = 72, JM = 46, LM = 23
#elif defined( GRIDF40 )
  integer, parameter :: IM = 144, JM = 90, LM = 40
#else
  integer, parameter :: IM = 0, JM = 0, LM = 0
#endif

  ! Constants
  real*8,  parameter :: xAVRX            = 1.d0  ! 2nd Order Method
  integer, parameter :: IMH=IM/2 ! Half num. of lat boxes, from MODEL_COM.f
  real*8,  parameter :: FIM=IM, byIM = 1./FIM ! from GEOM_B.f
  real*8,  parameter :: PI = 3.1415926535897932d0 ! from CONST.f
  real*8,  parameter :: TWOPI = 2d0 * PI ! from CONST.f
  real*8,  parameter :: RADIAN = PI / 180d0 ! from CONST.f
  real*8,  parameter :: g0 = 9.80665d0 ! [m s-2]
  real*8,  parameter :: D_LON = TWOPI*byIM ! grid spacing in longitude (deg), GEOM_B.f
                        ! Renamed from DLON, since conflict in CMN_SIZE_MOD.F
  real*8,  parameter :: teeny = 1d-30

  ! Values for global M23 simulation?
#if defined( GRIDM23 )
  integer, parameter :: I_0H = 1
  integer, parameter :: I_1H = IM
  integer, parameter :: J_0  = 1       ! J_STRT
  integer, parameter :: J_1  = JM      ! J_STOP
  integer, parameter :: J_0S = 2       ! J_STRT_SKP
  integer, parameter :: J_1S = JM-1    ! J_STOP_SKP
  integer, parameter :: J_0H = 1       ! J_STRT_HALO
  integer, parameter :: J_1H = JM      ! J_STOP_HALO
  integer, parameter :: J_0STG = 2     ! J_STRT_STGR
  integer, parameter :: J_1STG = JM    ! J_STOP_STGR
  integer, parameter :: do_polefix = 1 ! from MODEL_COM.f
  logical, parameter :: HAVE_SOUTH_POLE = .true.
  logical, parameter :: HAVE_NORTH_POLE = .true.
  real*8,  parameter :: D_LAT_DG = 180./(JM-1), D_LAT = D_LAT_DG * radian
#elif defined( GRIDF40 )
  integer, parameter :: I_0H = 1
  integer, parameter :: I_1H = IM
  integer, parameter :: J_0  = 1       ! J_STRT
  integer, parameter :: J_1  = JM      ! J_STOP
  integer, parameter :: J_0S = 2       ! J_STRT_SKP
  integer, parameter :: J_1S = JM-1    ! J_STOP_SKP
  integer, parameter :: J_0H = 1       ! J_STRT_HALO
  integer, parameter :: J_1H = JM      ! J_STOP_HALO
  integer, parameter :: J_0STG = 2     ! J_STRT_STGR
  integer, parameter :: J_1STG = JM    ! J_STOP_STGR
  integer, parameter :: do_polefix = 1 ! from MODEL_COM.f
  logical, parameter :: HAVE_SOUTH_POLE = .true.
  logical, parameter :: HAVE_NORTH_POLE = .true.
  real*8,  parameter :: D_LAT_DG = 180./(JM-1), D_LAT = D_LAT_DG * radian
#else ! Define for your simulation
  integer, parameter :: I_0H = -9999
  integer, parameter :: I_1H = -9999
  integer, parameter :: J_0  = -9999      
  integer, parameter :: J_1  = -9999      
  integer, parameter :: J_0S = -9999      
  integer, parameter :: J_1S = -9999    
  integer, parameter :: J_0H = -9999    
  integer, parameter :: J_1H = -9999    
  integer, parameter :: J_0STG = -9999  
  integer, parameter :: J_1STG = -9999 
  integer, parameter :: do_polefix = -9999
  logical, parameter :: HAVE_SOUTH_POLE = .true.
  logical, parameter :: HAVE_NORTH_POLE = .true.
  real*8,  parameter :: D_LAT_DG = -9999, D_LAT = D_LAT_DG * radian
#endif
  integer, parameter :: xstride = 1, ystride = IM, zstride = IM*(J_1H-J_0H+1)

  ! Variables set in subroutine GEOM_B used by AFLUX/AVRX
  !real*8, dimension(IM,JM) :: ZATMO ! geopotential height [m2/s2]
  real*8 :: POLWT
  real*8, dimension(JM) :: IMAXJ
  real*8, dimension(IM,JM) :: rel_area
  real*8, dimension(JM)    :: geofac
  real*8 :: geofac_pc
  real*8, dimension(IM,JM) :: xcolmass_fix

  real*8, dimension(JM)  :: DRAT
  integer, dimension(JM) :: NMIN
  real*8, dimension(IMH) :: BYSN

  !@var  DXP,DYP,BYDXP,BYDYP distance between points on primary grid
  !@+     (+inverse)
  REAL*8, DIMENSION(JM) :: BYDXP,DXP,DYP,BYDYP ! declared in module header
  !@var  DXV,DYV distance between velocity points (secondary grid)
  REAL*8, DIMENSION(JM) :: DYV,DXV ! declared in module header
  !@var  DXYP,BYDXYP area of grid box (+inverse) (m^2)
  !**** Note that this is not the exact area, but is what is required for
  !**** some B-grid conservation quantities
  REAL*8, DIMENSION(JM) :: BYDXYP, DXYP ! declared in module header
  !@var  DXYN,DXYS half box areas to the North,South of primary grid point
  REAL*8, DIMENSION(JM) :: DXYS,DXYN
  !@var AREAG global integral of area (m^2)
  REAL*8 :: AREAG
  !@var  RAPVS,RAPVN,RAVPS,RAVPN area scalings for primary and sec. grid
  REAL*8, DIMENSION(JM) :: RAPVS,RAPVN,RAVPS,RAVPN
  !@var DXYV,BYDXYV area of grid box around velocity point (recip.)(m^2)
  REAL*8, DIMENSION(JM) :: DXYV,BYDXYV
  !@var WTJ area weighting used in JLMAP, JKMAP (for hemispheric means)
  REAL*8, DIMENSION(JM,2,2) :: WTJ
  !@var  FCOR latitudinally varying coriolis parameter
  REAL*8, DIMENSION(JM) :: FCOR
  !@var SINIV,COSIV,SINIP,COSIP longitud. sin,cos for wind,pressure grid
  REAL*8, DIMENSION(IM) :: SINIV,COSIV,SINIP,COSIP
  !@var  RAVJ scaling for A grid U/V to B grid points (func. of lat. j)
  !@var  RAPJ scaling for B grid -> A grid conversion (1/4,1/im at poles)
  REAL*8, DIMENSION(IM,JM) :: RAPJ,RAVJ
  !@var  KMAXJ varying number of adjacent velocity points
  INTEGER, DIMENSION(JM) :: KMAXJ
  !@var  IDJJ J index of adjacent U/V points for A grid (func. of lat. j)
  INTEGER, DIMENSION(IM,JM) :: IDJJ
  !@var  IDIJ I index of adjacent U/V points for A grid (func. of lat/lon)
  INTEGER, DIMENSION(IM,IM,JM) :: IDIJ

  REAL*8, PARAMETER :: FJEQ=.5*(1+JM)    ! Equatorial value of J

  REAL*8, PARAMETER :: RADIUS = 6375000. ! Radius of spherical earth, same volume [m]
  REAL*8, PARAMETER :: EDPERD = 1.
  REAL*8, PARAMETER :: EDPERY = 365.
  REAL*8, PARAMETER :: SDAY = 86400. ! Seconds per day
  ! Earth's rotation rate (7.29 s^-1)
  REAL*8, PARAMETER :: OMEGA = TWOPI*(EDPERD+EDPERY)/(EDPERD*EDPERY*SDAY)

  real*8 :: LAT(JM), LAT_DG(JM,2), LON(IM), LON_DG(IM,2)

contains

  SUBROUTINE GEOM_B

    !USE DAO_MOD, ONLY : PHIS
    IMPLICIT NONE

    INTEGER :: I, J, K, N, IM1, JMHALF, JVPO
    REAL*8  :: RAVPO, LAT1, COSP1, DXP1

    REAL*8, DIMENSION(JM) :: SINP, COSP, COSV

    REAL*8 :: ACOR, ACOR2

    INTEGER :: ij
    REAL*8  :: dp           ! spacing in latitude (rad)
    REAL*8  :: ri2_gl
    REAL*8  :: rj2m1

    ! Assign ZATMO, the surface geopotential height [m2/s2]
    !ZATMO = State_Met%PHIS

    ! From GEOM_B.f
    ! Varying number of used longitudes
    IMAXJ(:)  = IM
    IMAXJ(1)  = 1
    IMAXJ(JM) = 1

    ! Calculate spherical geometry for B grid (from SUBROUTINE GEOM_B)
    LAT(1)  = -.25*TWOPI
    LAT(JM) = -LAT(1)
    SINP(1)  = -1.
    SINP(JM) = 1.
    COSP(1)  = 0.
    COSP(JM) = 0.
    DXP(1)  = 0.
    DXP(JM) = 0.
    DO J=2,JM-1
       LAT(J)  = D_LAT*(J-FJEQ)
       SINP(J) = SIN(LAT(J))
       COSP(J) = COS(LAT(J))
       DXP(J)  = RADIUS*D_LON*COSP(J)
    END DO
    BYDXP(2:JM-1) = 1.D0/DXP(2:JM-1)
    LAT1    = D_LAT*(1.-FJEQ)
    COSP1   = COS(LAT1)
    DXP1    = RADIUS*D_LON*COSP1
    DO J=2,JM
       COSV(J) = .5*(COSP(J-1)+COSP(J))
       DXV(J)  = .5*(DXP(J-1)+DXP(J))
       DYV(J)  = RADIUS*(LAT(J)-LAT(J-1))
       !**** The following corrections have no effect for half polar boxes
       !**** but are important for full and quarter polar box cases.
       IF (J.eq.2) THEN
          polwt = cosv(j)
          COSV(J) = .5*(COSP1+COSP(J))
          DXV(J)  = .5*(DXP1+DXP(J))
       END IF
       IF (J.eq.JM) THEN
          COSV(J) = .5*(COSP(J-1)+COSP1)
          DXV(J)  = .5*(DXP(J-1)+DXP1)
       END IF
       !****
    END DO
    DYP(1)  = RADIUS*(LAT(2)-LAT(1)-0.5*D_LAT)
    DYP(JM) = RADIUS*(LAT(JM)-LAT(JM-1)-0.5*D_LAT)
    DXYP(1) = .5*DXV(2)*DYP(1)
    BYDXYP(1) = 1./DXYP(1)
    DXYP(JM)= .5*DXV(JM)*DYP(JM)
    BYDXYP(JM) = 1./DXYP(JM)
    DXYS(1)  = 0.
    DXYS(JM) = DXYP(JM)
    DXYN(1)  = DXYP(1)
    DXYN(JM) = 0.
    polwt = (cosv(3)-cosv(2))/(cosv(3)-polwt)
    AREAG = DXYP(1)+DXYP(JM)

    DO J=2,JM-1
       DYP(J)  =  radius*D_LAT !.5*(DYV(J)+DYV(J+1))
       DXYP(J) = .5*(DXV(J)+DXV(J+1))*DYP(J)
       BYDXYP(J) = 1./DXYP(J)
       DXYS(J) = .5*DXYP(J)
       DXYN(J) = .5*DXYP(J)
       AREAG = AREAG+DXYP(J)
    END DO
    BYDYP(:) = 1.D0/DYP(:)
    AREAG = AREAG*FIM
    RAVPS(1)  = 0.
    RAPVS(1)  = 0.
    RAVPN(JM) = 0.
    RAPVN(JM) = 0.
    DO J=2,JM
       DXYV(J) = DXYN(J-1)+DXYS(J)
       BYDXYV(J) = 1./DXYV(J)
       RAPVS(J)   = .5*DXYS(J)/DXYV(J)
       RAPVN(J-1) = .5*DXYN(J-1)/DXYV(J)
       RAVPS(J)   = .5*DXYS(J)/DXYP(J)
       RAVPN(J-1) = .5*DXYN(J-1)/DXYP(J-1)
    END DO
    acor = dxyv(2)/(.5*dxp(2)*dyv(2)) ! gridbox area correction factor
    acor2 = dxyv(2)/(dxv(2)*dyv(2))
    !**** LONGITUDES (degrees); used in ILMAP
    LON_DG(1,1) = -180.+360./(2.*FLOAT(IM))
    LON_DG(1,2) = -180.+360./    FLOAT(IM)
    DO I=2,IM
       LON_DG(I,1) = LON_DG(I-1,1)+360./FLOAT(IM)
       LON_DG(I,2) = LON_DG(I-1,2)+360./FLOAT(IM)
    END DO
    !**** LATITUDES (degrees); used extensively in the diagn. print routines
    LAT_DG(1,1:2)=-90.
    LAT_DG(JM,1)=90.
    DO J=2,JM-1
       LAT_DG(J,1)=D_LAT_DG*(J-FJEQ)    ! primary (tracer) latitudes
    END DO
    DO J=2,JM
       LAT_DG(J,2)=D_LAT_DG*(J-JM/2-1)  ! secondary (velocity) latitudes
    END DO
    !**** WTJ: area weighting for JKMAP, JLMAP hemispheres
    JMHALF= JM/2
    DO J=1,JM
       WTJ(J,1,1)=1.
       WTJ(J,2,1)=2.*FIM*DXYP(J)/AREAG
    END DO
    DO J=2,JM
       WTJ(J,1,2)=1.
       WTJ(J,2,2)=2.*FIM*DXYV(J)/AREAG
    END DO
    !gsfc      WTJ(JMHALF+1,1,2)=.5
    !gsfc      WTJ(JMHALF+1,2,2)=WTJ(JMHALF+1,2,2)/2.
    WTJ(1,1,2)=0.
    WTJ(1,2,2)=0.
    !**** CALCULATE CORIOLIS PARAMETER
    !      OMEGA = TWOPI*(EDPERD+EDPERY)/(EDPERD*EDPERY*SDAY)
    FCOR(1)  = -OMEGA*DXV(2)*DXV(2)/D_LON
    FCOR(JM) = OMEGA*DXV(JM)*DXV(JM)/D_LON
    DO J=2,JM-1
       FCOR(J) = OMEGA*(DXV(J)*DXV(J)-DXV(J+1)*DXV(J+1))/D_LON
    END DO

    !**** Set indexes and scalings for the influence of A grid points on
    !**** adjacent velocity points

    !**** Calculate relative directions of polar box to nearby U,V points
    DO I=1,IM
       SINIV(I)=SIN((I-1)*D_LON)
       COSIV(I)=COS((I-1)*TWOPI*BYIM) ! D_LON)
       LON(I)=D_LON*(I-.5)
       SINIP(I)=SIN(LON(I))
       COSIP(I)=COS(LON(I))
    END DO

    !**** Conditions at the poles
    DO J=1,JM,JM-1
       IF(J.EQ.1) THEN
          JVPO=2
          RAVPO=2.*RAPVN(1)
       ELSE
          JVPO=JM
          RAVPO=2.*RAPVS(JM)
       END IF
       KMAXJ(J)=IM
       IMAXJ(J)=1
       RAVJ(1:KMAXJ(J),J)=RAVPO
       RAPJ(1:KMAXJ(J),J)=BYIM
       IDJJ(1:KMAXJ(J),J)=JVPO
       DO K=1,KMAXJ(J)
          IDIJ(K,1:IM,J)=K
       END DO
    END DO
    !**** Conditions at non-polar points
    DO J=2,JM-1
       KMAXJ(J)=4
       IMAXJ(J)=IM
       DO K=1,2
          RAVJ(K,J)=RAPVS(J)
          RAPJ(K,J)=RAVPS(J)    ! = .25
          IDJJ(K,J)=J
          RAVJ(K+2,J)=RAPVN(J)
          RAPJ(K+2,J)=RAVPN(J)  ! = .25
          IDJJ(K+2,J)=J+1
       END DO
       IM1=IM
       DO I=1,IM
          IDIJ(1,I,J)=IM1
          IDIJ(2,I,J)=I
          IDIJ(3,I,J)=IM1
          IDIJ(4,I,J)=I
          IM1=I
       END DO
    END DO

    !print*,'Primary (Tracer) grid: '
    !print*,LAT_DG(:,1)
    !print*,LON_DG(:,1)

    !print*,'Secondary (Velocity) grid: '
    !print*,LAT_DG(2:,2)
    !print*,LON_DG(2:,2)

    DO J=1,JM
       rel_area(:,j) = dxyp(j) / ( SUM(dxyp) * dble(IM) )
    ENDDO

    dp = PI / ( JM - 1 )
    do ij = 1, JM
       geofac(ij) = dp / (2.0d0 * SUM(rel_area(:,ij)))
    end do
    geofac_pc = &
         dp / (2.0d0 * SUM(rel_area(:,1)))



    ! Now that the GEOM_B.f code is finished, calculate the parameters from
    ! AVRX in ATMDYN.f
    BYSN = 0d0
    DRAT = 0d0
    NMIN = 0
    DO N=1,IMH
       BYSN(N)=xAVRX/SIN(.5*D_LON*N)
    ENDDO
    DO J=1,JM
       DRAT(J) = DXP(J)*BYDYP(3)
       DO N=IMH,1,-1
          IF(BYSN(N)*DRAT(J) .GT.1.) THEN
             NMIN(J) = N+1
             EXIT
          ENDIF
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE GEOM_B

end module geom_giss_mod
