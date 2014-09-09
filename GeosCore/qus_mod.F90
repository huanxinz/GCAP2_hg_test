module qus_mod

  USE CMN_SIZE_MOD
  USE GEOM_GISS_MOD
  USE ERROR_MOD, ONLY : GEOS_CHEM_STOP

  implicit none

  PRIVATE

  !PUBLIC :: advect_x, advect_y, advect_z                     ! Advection subroutines
  !PUBLIC :: S0, SM, SX, SY, SZ, SXX, SYY, SZZ, SXY, SXZ, SYZ ! Moments
  PUBLIC :: NMOM, TrMom, TrMPrev
  PUBLIC :: AVRX, AFLUX, CALC_PIJL, CALC_AMP, AADVT, PFIX_QUS
  PUBLIC :: PU, PV, SD, PIT, MB, MA

  ! Dimensions
#if defined( GRIDM23 )
  integer, parameter :: IM = 72, JM = 46, LM = 23
#elif defined( GRIDF40 )
!  integer, parameter :: IM = 144, JM = 91, LM = 40
#else
   integer, parameter :: IM=0, JM=0, LM=0
#endif

  logical, SAVE :: FIRST=.TRUE.

  ! Global moments
  !real*8 ::  S0(IM,JM,LM),  SM(IM,JM,LM)
  !real*8 ::  SX(IM,JM,LM),  SY(IM,JM,LM),  SZ(IM,JM,LM)
  !real*8 :: SXX(IM,JM,LM), SYY(IM,JM,LM), SZZ(IM,JM,LM)
  !real*8 :: SXY(IM,JM,LM), SXZ(IM,JM,LM), SYZ(IM,JM,LM)

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Define moments as in GISS  !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  real*8,  parameter :: xAVRX            = 1.d0  ! 2nd Order Method
  integer, parameter :: prather_limits   =  0    ! GISS default
  integer, parameter :: flux_negative    = -1
  integer, parameter :: flux_nonnegative = +1

  integer, parameter :: nmom=9, mx=1,  my=2,  mz=3
  integer, parameter :: mxx=4, myy=5, mzz=6, mxy=7, mzx=8, myz=9
  
  ! moments with vertical component
  integer, dimension(4), parameter :: zmoms=(/MZ,MZZ,MYZ,MZX/)
  ! moments with no vertical component
  integer, dimension(5), parameter :: xymoms=(/MX,MY,MXX,MXY,MYY/)
  ! moments with a horizontal component
  integer, dimension(7), parameter :: ihmoms=(/MX,MY,MXX,MYY,MXY,MYZ,MZX/)
  ! moments with no horizontal component
  integer, dimension(2), parameter :: zomoms=(/MZ,MZZ/)
  ! x-x, x-y, x-z switches
  integer, dimension(nmom), parameter :: xdir=(/mx,my,mz,mxx,myy,mzz,mxy,myz,mzx/)
  integer, dimension(nmom), parameter :: ydir=(/my,mx,mz,myy,mxx,mzz,mxy,mzx,myz/)
  integer, dimension(nmom), parameter :: zdir=(/mz,my,mx,mzz,myy,mxx,myz,mxy,mzx/)

  real*8, allocatable :: TrMom(:,:,:,:,:)
  real*8, allocatable :: TrMPrev(:,:,:,:)

  ! Time step limits
  INTEGER, PARAMETER :: NCMAX = 10
  INTEGER :: NSTEPX1(JM,LM,NCMAX)
  INTEGER :: NSTEPX2(JM,LM,NCMAX)
  INTEGER :: NSTEPY(LM,NCMAX) 
  INTEGER :: NSTEPZ(IM*JM,NCMAX)
  INTEGER :: NCYC

  ! What are these values for global M23 simulation?
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
#else
  integer, parameter :: I_0H = -99
  integer, parameter :: I_1H = -99
  integer, parameter :: J_0  = -99       ! J_STRT
  integer, parameter :: J_1  = -99      ! J_STOP
  integer, parameter :: J_0S = -99       ! J_STRT_SKP
  integer, parameter :: J_1S = -99    ! J_STOP_SKP
  integer, parameter :: J_0H = 1       ! J_STRT_HALO
  integer, parameter :: J_1H = JM     ! J_STOP_HALO
  integer, parameter :: J_0STG = -99     ! J_STRT_STGR
  integer, parameter :: J_1STG = -99    ! J_STOP_STGR
  integer, parameter :: do_polefix = -99 ! from MODEL_COM.f
  logical, parameter :: HAVE_SOUTH_POLE = .true.
  logical, parameter :: HAVE_NORTH_POLE = .true.
  integer, parameter :: LS1 = 1
  real*8, parameter  :: PSFMPT = 984
#endif
  integer, parameter :: xstride = 1, ystride = IM, zstride = IM*(J_1H-J_0H+1)

  ! Velocity field across eastern, northern, and upper edges [kg/s]
  real*8 :: UFLUX(IM,JM,LM), VFLUX(IM,JM,LM), WFLUX(IM,JM,LM)
  
  ! Mass Flux
  real*8, dimension(IM,JM,LM)         :: PU, PV    ! Horizontal mass flux [mb m^2 s-1]
  real*8, dimension(IM,JM,LM)         :: SPA       ! ?
  real*8, dimension(IM,JM,LM), target :: CONV
  real*8, dimension(IM,JM,LM)         :: MB, MA    ! Fluid mass before/after
  real*8, pointer                     :: SD(:,:,:) ! Sigma dot [mb m^2 s-1]
  real*8, pointer                     :: PIT(:,:)  ! Pressure tendency [mb m^2 s-1]

  real*8, dimension(JM,LM)            :: SCF, SCM, SFCM
  real*8, dimension(JM,LM)            :: SBF, SBM, SFBM
  
  integer, parameter :: IMH=IM/2 ! Half num. of lat boxes, from MODEL_COM.f
  real*8,  parameter :: FIM=IM, byIM = 1./FIM ! from GEOM_B.f
  real*8,  parameter :: PI = 3.1415926535897932d0 ! from CONST.f
  real*8,  parameter :: TWOPI = 2d0 * PI ! from CONST.f
  real*8,  parameter :: RADIAN = PI / 180d0 ! from CONST.f
  real*8,  parameter :: g0 = 9.80665d0 ! [m s-2]
  real*8,  parameter :: D_LON = TWOPI*byIM ! grid spacing in longitude (deg), GEOM_B.f
  real*8,  parameter :: teeny = 1d-30
#if defined( GRIDM23 )
  real*8,  parameter :: D_LAT_DG = 180./(JM-1), D_LAT = D_LAT_DG * radian
#else
!  real*8,  parameter :: D_LAT_DG = -99., D_LAT=-99.
#endif

contains

  SUBROUTINE AFLUX( ZATMO, U, V, PIJL )!, State_Met )
    ! Modified from the subroutine in Model E's ATMDYN.f for GEOS-Chem. 
    ! Added pressure fixer.

    !@sum  AFLUX Calculates horizontal/vertical air mass fluxes
    !@+    Input: U,V velocities, PIJL pressure (surface press <= LS1, PSMPT otherwise)
    !@+    Output: PIT  pressure tendency (mb m^2/s)
    !@+            SD   sigma dot (mb m^2/s)
    !@+            PU,PV horizontal mass fluxes (mb m^2/s)
    !@+            CONV  horizontal mass convergence (mb m^2/s)
    !@+            SPA
    !@auth Original development team
    !@ver  1.0
    
    USE PRESSURE_MOD, ONLY : SIGE, DSIG, byDSIG
    USE DIAG_MOD,     ONLY : MASSFLEW, MASSFLNS, MASSFLUP
    !USE GIGC_State_Met_Mod, ONLY : MetState

    IMPLICIT NONE
    
    !**** CONSTANT PRESSURE AT L=LS1 AND ABOVE, PU,PV CONTAIN DSIG
    !@var U,V input velocities (m/s)
    !@var PIJL input 3-D pressure field (mb) (no DSIG)

    ! Arguments
    REAL*8, INTENT(INOUT), DIMENSION(IM,JM,LM) :: U,V,PIJL
    REAL*8, INTENT(IN), DIMENSION(IM,JM) :: ZATMO
    !TYPE(MetState), INTENT(IN) :: State_Met  ! Meteorology State object
    
    ! Local variables
    REAL*8, DIMENSION(IM) :: DUMMYS,DUMMYN

    INTEGER I,J,L,IP1,IM1,IPOLE
    REAL*8 PUS,PUN,PVS,PVN,PBS,PBN
    REAL*8 WT,DXDSIG,DYDSIG,PVSA(LM),PVNA(LM),xx,twoby3
    real*8, dimension(IM,2,LM) :: usv0,vsv0
    integer :: jvs,jvn,jv
    real*8 :: wts

    IF ( FIRST ) THEN
       ! Assign pointers
       SD  => CONV(:,:,2:)
       PIT => CONV(:,:,1)
       FIRST = .FALSE.
    ENDIF

!****
!**** BEGINNING OF LAYER LOOP
!****
!$OMP PARALLEL DO PRIVATE (I,J,L,IM1,IP1,DXDSIG,DYDSIG,DUMMYS,DUMMYN,PUS,PUN,PVS,PVN,PBS,PBN)
      DO 2000 L=1,LM
!****
!**** COMPUTATION OF MASS FLUXES     P,T  PU     PRIMARY GRID ROW
!**** ARAKAWA'S SCHEME B             PV   U,V    SECONDARY GRID ROW
!****
!**** COMPUTE PU, THE WEST-EAST MASS FLUX, AT NON-POLAR POINTS
      DO 2154 J=2,JM-1
      DO 2154 I=1,IM
 2154 SPA(I,J,L)=U(I,J,L)+U(I,J+1,L)
      !CALL AVRX (SPA(1,1,L))
      I=IM
      DO 2166 J=2,JM-1
      DYDSIG = 0.25D0*DYP(J)*DSIG(L)
      DO 2165 IP1=1,IM
      PU(I,J,L)=DYDSIG*SPA(I,J,L)*(PIJL(I,J,L)+PIJL(IP1,J,L))
 2165 I=IP1
 2166 CONTINUE
!**** COMPUTE PV, THE SOUTH-NORTH MASS FLUX
      IM1=IM
      DO 2172 J=2,JM
      DXDSIG = 0.25D0*DXV(J)*DSIG(L)
      DO 2170 I=1,IM
      PV(I,J,L)=DXDSIG*(V(I,J,L)+V(IM1,J,L))*(PIJL(I,J,L)+PIJL(I,J-1,L))

 2170 IM1=I
 2172 CONTINUE
!**** COMPUTE PU*3 AT THE POLES
      PUS=0.
      PUN=0.
      PVS=0.
      PVN=0.
      DO 1110 I=1,IM
      PUS=PUS+U(I,2,L)
      PUN=PUN+U(I,JM,L)
      PVS=PVS+PV(I,2,L)
 1110 PVN=PVN+PV(I,JM,L)
      PUS=.25*DYP(2)*PUS*PIJL(1,1,L)*BYIM
      PUN=.25*DYP(JM-1)*PUN*PIJL(1,JM,L)*BYIM
      PVS=PVS*BYIM
      PVN=PVN*BYIM
      PVSA(L)=PVS
      PVNA(L)=PVN
      DUMMYS(1)=0.
      DUMMYN(1)=0.
      DO 1120 I=2,IM
      DUMMYS(I)=DUMMYS(I-1)+(PV(I,2,L) -PVS)*BYDSIG(L)
 1120 DUMMYN(I)=DUMMYN(I-1)+(PV(I,JM,L)-PVN)*BYDSIG(L)
      PBS=0.
      PBN=0.
      DO 1130 I=1,IM
      PBS=PBS+DUMMYS(I)
 1130 PBN=PBN+DUMMYN(I)
      PBS=PBS*BYIM
      PBN=PBN*BYIM
      DO 1140 I=1,IM
      PV(I,1,L)=0.
      SPA(I,1,L)=4.*(PBS-DUMMYS(I)+PUS)/(DYP(2)*PIJL(1,1,L))
      SPA(I,JM,L)=4.*(DUMMYN(I)-PBN+PUN)/(DYP(JM-1)*PIJL(1,JM,L))
      PU(I,1,L)=3.*(PBS-DUMMYS(I)+PUS)*DSIG(L)
 1140 PU(I,JM,L)=3.*(DUMMYN(I)-PBN+PUN)*DSIG(L)
!****
!**** CONTINUITY EQUATION
!****
!**** COMPUTE CONV, THE HORIZONTAL MASS CONVERGENCE
!     DO 1510 J=2,JM-1
!     IM1=IM
!     DO 1510 I=1,IM
!     CONV(I,J,L)=(PU(IM1,J,L)-PU(I,J,L)+PV(I,J,L)-PV(I,J+1,L))
!1510 IM1=I
!     CONV(1,1,L)=-PVS
!     CONV(1,JM,L)=PVN
 2000 CONTINUE
!$OMP  END PARALLEL DO
!****
!**** END OF HORIZONTAL ADVECTION LAYER LOOP
!****
!
! modify uphill air mass fluxes around steep topography
      do 2015 j=2,jm-1
      i = im
      do 2010 ip1=1,im
         xx = zatmo(ip1,j)-zatmo(i,j)
         if(xx.eq.0.0)  go to 2007
         DO 2005 L=1,LS1-1
!cc         if((zatmo(ip1,j)-zatmo(i,j))*pu(i,j,l).gt.0.) then
            if(xx*pu(i,j,l).gt.0.) then
               if(pu(i,j,l).gt.0.) then
                  wt = (pijl(ip1,j,l)/pijl(i,j,l)-sige(l+1))/dsig(l)
               else
                  wt = (pijl(i,j,l)/pijl(ip1,j,l)-sige(l+1))/dsig(l)
               endif
               if(wt.le.0.) then
                  pu(i,j,l+1) = pu(i,j,l+1) + pu(i,j,l)
                  pu(i,j,l) = 0.
               else
                  go to 2007
               endif
            endif
 2005    CONTINUE
 2007    CONTINUE
         i = ip1
 2010 CONTINUE
 2015 CONTINUE
!cc   do j=2,jm
      do 2035 j=3,jm-1
      do 2035 i=1,im
         xx = zatmo(i,j)-zatmo(i,j-1)
         if(xx.eq.0.0)  go to 2035
         DO 2020 L=1,LS1-1
!cc         if((zatmo(i,j)-zatmo(i,j-1))*pv(i,j,l).gt.0.) then
            if(xx*pv(i,j,l).gt.0.) then
               if(pv(i,j,l).gt.0.) then
                  wt = (pijl(i,j,l)/pijl(i,j-1,l)-sige(l+1))/dsig(l)
               else
                  wt = (pijl(i,j-1,l)/pijl(i,j,l)-sige(l+1))/dsig(l)
               endif
               if(wt.le.0.) then
                  pv(i,j,l+1) = pv(i,j,l+1) + pv(i,j,l)
                  pv(i,j,l) = 0.
               else
                  go to 2035
               endif
            endif
 2020    CONTINUE
 2035 CONTINUE
!
!     Now Really Do  CONTINUITY EQUATION
!
!     COMPUTE CONV, THE HORIZONTAL MASS CONVERGENCE
!
!$OMP  PARALLEL DO PRIVATE (I,J,L,IM1)
      DO 2400 L=1,LM
      DO 1510 J=2,JM-1
      IM1=IM
      DO 1510 I=1,IM
      CONV(I,J,L)=(PU(IM1,J,L)-PU(I,J,L)+PV(I,J,L)-PV(I,J+1,L))
 1510 IM1=I
      CONV(1,1,L)=-PVSA(L)
      CONV(1,JM,L)=PVNA(L)
 2400 CONTINUE
!$OMP  END PARALLEL DO
!
!**** COMPUTE PIT, THE PRESSURE TENDENCY
!C    PIT(I,J)=CONV(I,J,1)
!     DO 2420 L=LM,2,-1
!     PIT(1,1)=PIT(1,1)+CONV(1,1,L)
!     PIT(1,JM)=PIT(1,JM)+CONV(1,JM,L)
!     DO 2420 J=2,JM-1
!     DO 2420 I=1,IM
!2420 PIT(I,J)=PIT(I,J)+CONV(I,J,L)
!$OMP  PARALLEL DO PRIVATE(I,J,L)
      DO 2420 J=1,JM
         DO 2410 I=1,IMAXJ(J)
         DO 2410 L=LM-1,1,-1
            PIT(I,J)=PIT(I,J)+SD(I,J,L)
 2410    CONTINUE
 2420 CONTINUE
!$OMP  END PARALLEL DO
!**** COMPUTE SD, SIGMA DOT                        -------
!     SD(1, 1,LM-1)=CONV(1, 1,LM)                     |
!     SD(1,JM,LM-1)=CONV(1,JM,LM)             completely wasteful
!     DO 2430 J=2,JM-1                                |
!     DO 2430 I=1,IM                                  |
!2430 SD(I,J,LM-1)=CONV(I,J,LM)                    -------
!     DO 2435 L=LM-2,LS1-1,-1
!     SD(1, 1,L)=SD(1, 1,L+1)+CONV(1, 1,L+1)
!     SD(1,JM,L)=SD(1,JM,L+1)+CONV(1,JM,L+1)
!     DO 2435 J=2,JM-1
!     DO 2435 I=1,IM
!     SD(I, J,L)=SD(I, J,L+1)+CONV(I, J,L+1)
!2435 CONTINUE
!$OMP  PARALLEL DO PRIVATE(I,J,L)
      DO 2435 J=1,JM
         DO 2430 I=1,IMAXJ(J)
         DO 2430 L=LM-2,LS1-1,-1
            SD(I,J,L)=SD(I,J,L+1)+SD(I,J,L)
 2430    CONTINUE
 2435 CONTINUE
!$OMP  END PARALLEL DO
!     DO 2440 L=LS1-2,1,-1
!     SD(1, 1,L)=SD(1, 1,L+1)+CONV(1, 1,L+1)-DSIG(L+1)*PIT(1, 1)
!     SD(1,JM,L)=SD(1,JM,L+1)+CONV(1,JM,L+1)-DSIG(L+1)*PIT(1,JM)
!     DO 2440 J=2,JM-1
!     DO 2440 I=1,IM
!     SD(I, J,L)=SD(I, J,L+1)+CONV(I, J,L+1)-DSIG(L+1)*PIT(I, J)
!2440 CONTINUE
!$OMP  PARALLEL DO PRIVATE(I,J,L)
      DO 2440 J=1,JM
         DO 2438 I=1,IMAXJ(J)
         DO 2438 L=LS1-2,1,-1
            SD(I,J,L)=SD(I,J,L+1)+SD(I,J,L)-DSIG(L+1)*PIT(I,J)
 2438    CONTINUE
 2440 CONTINUE
!$OMP  END PARALLEL DO
      DO 2450 L=1,LM-1
      DO 2450 I=2,IM
      SD(I,1,L)=SD(1,1,L)
 2450 SD(I,JM,L)=SD(1,JM,L)
!**** temporary fix for CLOUDS module
!      SD_CLOUDS(:,:,1)    = PIT
!!$OMP PARALLEL DO PRIVATE (L)
!      DO L=2,LM
!        SD_CLOUDS(:,:,L) = SD(:,:,L-1)
!      END DO
!!$OMP END PARALLEL DO
!****

      !WRITE(6,*) IM, JM, LM     
      !WRITE(6,*) I, J, L, IP1, IM1, IPOLE
      !WRITE(6,*) PUS, PUN, PVS, PVN, PBS, PBN
      !WRITE(6,*) WT, DXDSIG, DYDSIG, XX, TWOBY3
      !WRITE(6,*) JVS, JVS, JV
      !WRITE(6,*) WTS

      !I=54
      !WRITE(6,'(46(F5.1","))') U(I,:,1)
      !WRITE(6,'(46(F5.1","))') V(I,:,1)
      !WRITE(6,'(46(F5.0","))') PIJL(I,:,LS1-1)
      !WRITE(6,'(46(E9.2","))') PU(I,:,1)
      !WRITE(6,'(46(E9.2","))') PV(I,:,1)
      !WRITE(6,'(46(E9.2","))') SD(I,:,1)
      !CALL GEOS_CHEM_STOP

      RETURN
  END SUBROUTINE AFLUX
  
  SUBROUTINE AVRX(X,jrange)
    !@sum  AVRX Smoothes zonal mass flux and geopotential near the poles
    !@auth Original development team
    !@ver  1.0

    USE ERROR_MOD, ONLY : GEOS_CHEM_STOP

    IMPLICIT NONE

    REAL*8, INTENT(INOUT), optional :: &
         X(IM,J_0H:J_1H)
    Integer, Intent(In), optional :: jrange(2)
    REAL*8, ALLOCATABLE, SAVE  :: DRAT(:)
    REAL*8, SAVE ::  BYSN(IMH)
    REAL*8, DIMENSION(0:IMH) :: AN,BN
    INTEGER, ALLOCATABLE, SAVE :: NMIN(:)
    INTEGER J,N
    LOGICAL, SAVE :: init = .false.

    INTEGER :: J0, J1

    if ( present(X) ) goto 1000
    IF (.NOT. init) THEN
       init = .true.
       CALL FFT0(IM)
       NMIN=0
       j0 = MAX(1,J_0H)
       j1 = MIN(JM,J_1H)
       ALLOCATE(DRAT(j0:j1), NMIN(j0:j1))
       DO N=1,IMH
          BYSN(N)=xAVRX/SIN(.5*D_LON*N)
       END DO
       DO J=j0,j1
          DRAT(J) = DXP(J)*BYDYP(3)
          DO N=IMH,1,-1
             IF(BYSN(N)*DRAT(J) .GT.1.) THEN
                NMIN(J) = N+1
                EXIT
             ENDIF
          END DO
       END DO
    END IF
    RETURN
    !****
!!!      ENTRY AVRX (X)
1000 continue
    !****

    CALL FFT0(IM)
    If (Present(jrange)) Then
       j0 = jrange(1)
       j1 = jrange(2)
    Else
       j0=J_0S
       j1=J_1S
    End If

    DO J=j0,j1
       IF (DRAT(J).GT.1) CYCLE
       CALL FFT (X(1,J),AN,BN)
       DO N=NMIN(J),IMH-1
          AN(N)=BYSN(N)*DRAT(J) * AN(N)
          BN(N)=BYSN(N)*DRAT(J) * BN(N)
       END DO
       AN(IMH) = BYSN(IMH)*DRAT(J) * AN(IMH)
       CALL FFTI(AN,BN,X(1,J))
    END DO

    RETURN
  END SUBROUTINE AVRX

  SUBROUTINE CALC_PIJL(lmax,p,pijl)
!@sum  CALC_PIJL Fills in P as 3-D
!@auth Jean Lerner
!@ver  1.0
    
    USE CMN_SIZE_MOD
    
    implicit none
    
    REAL*8, dimension(IM,JM),    intent(IN)  :: p    ! Surface pressure
    REAL*8, dimension(IM,JM,LM), intent(OUT) :: pijl
    INTEGER, intent(IN)                      :: lmax
    INTEGER :: L

    do l=1,ls1-1
       pijl(:,:,l) = p(:,:) - PTOP
    enddo
    do l=ls1,lmax
       pijl(:,:,l) = PSFMPT
    enddo
    return
  end subroutine calc_pijl


  SUBROUTINE CALC_AMP(p,amp)
    !@sum  CALC_AMP Calc. AMP: kg air*grav/100, incl. const. pressure strat
    !@auth Jean Lerner/Max Kelley
    !@ver  1.0
    USE CMN_SIZE_MOD
    USE PRESSURE_MOD, ONLY : dSIG
    
    implicit none
    REAL*8, dimension(im,jm)    :: p ! Surface pressure
    REAL*8, dimension(im,jm,lm) :: amp
    integer :: j,l
    logical, save :: first = .true.

    if ( first ) then
       CALL GEOM_B
       first = .false.
    endif

!$OMP PARALLEL DO        &
!$OMP DEFAULT( SHARED   )&
!$OMP PRIVATE( J, L )
    DO L=1,LM
       IF(L.LT.LS1) THEN
          do j=1,JM
             amp(:,j,l) = (p(:,j)-ptop)*dxyp(j)*dsig(l)
          enddo
       ELSE
          do j=1,JM
             amp(:,j,l) = (psfmpt)*dxyp(j)*dsig(l)
          enddo
       END IF
    enddo
!$OMP END PARALLEL DO

    return
  end subroutine calc_amp
    
  subroutine adv1d(s,smom, f,fmom, mass,dm, nx,qlimit,stride,dir,ierr,nerr)
      !@sum  adv1d implements the quadratic upstream scheme in one dimension
      !@auth G. Russell, modified by Maxwell Kelley
      !--------------------------------------------------------------
      ! adv1d advects tracers in x-direction using the qus
      ! the order of the moments in dir is: x,y,z,xx,yy,zz,xy,yz,zx
      !--------------------------------------------------------------
      implicit none
      !
      !@var s      mean tracer amount (kg or J)
      !@var smom   qus tracer moments (kg or J)
      !@var f      tracer flux (diagnostic output) (kg or J)
      !@var fmom   tracer moment flux (diagnostic output) (kg or J)
      !@var mass   mass field (kg)
      !@var dm     mass flux (kg)
      !@var nx     length of 1D vector
      !@var qlimit true if negative tracer is to be avoided
      !@var stride spacing in s array between elements of relevant 1D array
      !@var dir    direction switch (equals one of xdir ydir or zdir)
      !@var ierr, nerr error codes
      integer, intent(in) :: nx,stride
      logical, intent(in) :: qlimit
      REAL*8, dimension(nx) :: dm, f
      REAL*8, dimension(nx*stride) :: s,mass
      REAL*8, dimension(nmom,nx*stride) :: smom
      REAL*8, dimension(nmom,nx) :: fmom
      integer, dimension(nmom) :: dir
      integer :: mx,my,mz,mxx,myy,mzz,mxy,myz,mzx
      integer :: n,np1,nm1,nn,ns
      integer,intent(out) :: ierr,nerr
      REAL*8 :: fracm,frac1,bymnew,mnew,dm2,tmp
      ! qlimit variables
      REAL*8 :: an, anm1, fn, fnm1, sn, sxn, sxxn
      !
      ierr=0 ; nerr=0
      mx  = dir(1)
      my  = dir(2)
      mz  = dir(3)
      mxx = dir(4)
      myy = dir(5)
      mzz = dir(6)
      mxy = dir(7)
      myz = dir(8)
      mzx = dir(9)
      !-----------------------------------------------------------
      ! calculate tracer mass flux f
      !-----------------------------------------------------------
      n = nx
      do np1=1,nx

         if(dm(n).lt.0.) then ! air mass flux is negative
            nn=np1
            frac1=+1.
         else                 ! air mass flux is positive
            nn=n
            frac1=-1.
         endif
         ns=1+(nn-1)*stride
         fracm=dm(n)/mass(ns)
         if(mass(ns).le.0.d0) fracm=0.d0
         frac1=fracm+frac1
         f(n)=fracm*(s(ns)-frac1*(smom(mx,ns)-&
              (frac1+fracm)*smom(mxx,ns)))
         ! temporary storage of fracm in fx, to be used below
         fmom(mx,n)=fracm
         ! temporary storage of frac1 in fxx, to be used below
         fmom(mxx,n)=frac1
         !
         n = np1
      enddo
      if(qlimit) then
         nm1 = nx
         do n=1,nx
            ns=1+(n-1)*stride
            an = fmom(mx,n)      ! reading fracm which was stored in fx
            anm1 = fmom(mx,nm1)
            fn = f(n)
            fnm1 = f(nm1)
            sn = s(ns)
            sxn = smom(mx,ns)
            sxxn = smom(mxx,ns)
            call limitq(anm1,an,fnm1,fn,sn,sxn,sxxn,ierr)
            if (ierr.gt.0) then
               nerr=n
               if (ierr.eq.2) return
            end if
            f(n) = fn
            f(nm1) = fnm1
            smom(mx,ns) = sxn
            smom(mxx,ns) = sxxn
            nm1 = n
         enddo
      endif
      !--------------------------------------------------------------------
      ! calculate tracer fluxes of slopes and curvatures
      !--------------------------------------------------------------------
      n = nx
      do np1=1,nx
         if(dm(n).lt.0.) then ! air mass flux is negative
            nn=np1
         else                 ! air mass flux is positive
            nn=n
         endif
         ns=1+(nn-1)*stride
         ! retrieving fracm, which was stored in fx
         fracm=fmom(mx,n)
         ! retrieving frac1, which was stored in fxx
         frac1=fmom(mxx,n)
         !
         fmom(mx,n)=dm(n)*(fracm*fracm*(smom(mx,ns) &
              -3.*frac1*smom(mxx,ns))-3.*f(n))
         fmom(mxx,n)=dm(n)*(dm(n)*fracm**3 *smom(mxx,ns) &
              -5.*(dm(n)*f(n)+fmom(mx,n)))
         ! cross moments
         fmom(my,n)  = fracm*(smom(my,ns)-frac1*smom(mxy,ns))
         fmom(mxy,n) = dm(n)*(fracm*fracm*smom(mxy,ns)-3.*fmom(my,n))
         fmom(mz,n)  = fracm*(smom(mz,ns)-frac1*smom(mzx,ns))
         fmom(mzx,n) = dm(n)*(fracm*fracm*smom(mzx,ns)-3.*fmom(mz,n))
         fmom(myy,n) = fracm*smom(myy,ns)
         fmom(mzz,n) = fracm*smom(mzz,ns)
         fmom(myz,n) = fracm*smom(myz,ns)
         n = np1
      enddo
      !-------------------------------------------------------------------
      ! update tracer mass, moments of tracer mass, air mass distribution
      !-------------------------------------------------------------------
      nm1 = nx
      do n=1,nx
         ns=1+(n-1)*stride
         tmp=mass(ns)+dm(nm1)
         mnew=tmp-dm(n)
         !     mnew=mass(ns)+dm(nm1)-dm(n)
         bymnew = 1./mnew
         dm2=dm(nm1)+dm(n)
         tmp=s(ns)+f(nm1)
         s(ns)=tmp-f(n)
         !     s(ns)=s(ns)+f(nm1)-f(n)
         !
         smom(mx,ns)=(smom(mx,ns)*mass(ns)-3.*(-dm2*s(ns) &
              +mass(ns)*(f(nm1)+f(n)))+(fmom(mx,nm1)-fmom(mx,n)))*bymnew
         smom(mxx,ns) = (smom(mxx,ns)*mass(ns)*mass(ns) &
              +2.5*s(ns)*(mass(ns)*mass(ns)-mnew*mnew-3.*dm2*dm2) &
              +5.*(mass(ns)*(mass(ns)*(f(nm1)-f(n))-fmom(mx,nm1) &
              -fmom(mx,n))+dm2*smom(mx,ns)*mnew) &
              +(fmom(mxx,nm1)-fmom(mxx,n))) * (bymnew*bymnew)
         ! cross moments
         smom(my,ns)=smom(my,ns)+fmom(my,nm1)-fmom(my,n)
         smom(mxy,ns)=(smom(mxy,ns)*mass(ns)-3.*(-dm2*smom(my,ns) + &
              mass(ns)*(fmom(my,nm1)+fmom(my,n))) + &
              (fmom(mxy,nm1)-fmom(mxy,n)))*bymnew
         smom(mz,ns)=smom(mz,ns)+fmom(mz,nm1)-fmom(mz,n)
         smom(mzx,ns)=(smom(mzx,ns)*mass(ns)-3.*(-dm2*smom(mz,ns) + &
              mass(ns)*(fmom(mz,nm1)+fmom(mz,n))) + &
              (fmom(mzx,nm1)-fmom(mzx,n)))*bymnew
         !
         smom(myy,ns)=smom(myy,ns)+fmom(myy,nm1)-fmom(myy,n)
         smom(mzz,ns)=smom(mzz,ns)+fmom(mzz,nm1)-fmom(mzz,n)
         smom(myz,ns)=smom(myz,ns)+fmom(myz,nm1)-fmom(myz,n)
         !
         !------------------------------------------------------------------
         mass(ns) = mnew
         if(mass(ns).le.0.) then
            s(ns)=0.
            smom(:,ns)=0.
         endif
         if (qlimit .and. prather_limits.eq.1) then ! force Prather limits
            smom(mx,ns)=min(1.5*s(ns),max(-1.5*s(ns),smom(mx,ns)))
            smom(mxx,ns)=&
              min(2.*s(ns)-abs(smom(mx,ns))/3.,max(abs(smom(mx,ns))-s(ns),smom(mxx,ns)))
            smom(mxy,ns)=min(s(ns),max(-s(ns),smom(mxy,ns)))
            smom(mzx,ns)=min(s(ns),max(-s(ns),smom(mzx,ns)))
         end if
         !-----------------------------------------------------------------
         nm1 = n
      enddo
      return
    end subroutine adv1d
    
    subroutine limitq(anm1,an,fnm1,fn,sn,sx,sxx,ierr)
      !@sum  limitq adjusts moments to maintain non-neg. tracer means/fluxes
      !@auth G. Russell, modified by Maxwell Kelley
      implicit none
      REAL*8 :: anm1,an,fnm1,fn,sn,sx,sxx

      ! local variables
      REAL*8 :: sl,sc,sr, frl,frl1, frr,frr1, gamma,g13ab, fr,fr1, fsign,su,sd
      integer, intent(out) ::ierr
      ierr=0

!****
!**** modify the tracer moments so that the tracer mass in each
!**** division is non-negative
!****
!**** no air leaving the box
        if(anm1.ge.0. .and. an.le.0.) return
!**** air is leaving through both the left and right edges
        if(anm1.lt.0. .and. an.gt.0.) then
           sl = -fnm1
           sr = +fn
           sc = sn - sl
           sc = sc - sr
!**** all three divisions are non-negative
           if(sl.ge.0. .and. sr.ge.0. .and. sc.ge.0.) return
!**** first check for the cases when only one of the three is negative
           frl = anm1
           frl1 = frl+1.
           frr = an
           frr1 = frr-1.
           if(sl.ge.0. .and. sr.ge.0.) then ! center division
              gamma = 1.+(frl-frr)
              g13ab = gamma*gamma - 1. + 3.*(frl+frr)**2
              sxx = sxx - sc*10.*g13ab / &
                  (gamma*(12.*(frl+frr)**2 + 5.*g13ab*g13ab))
              sx = sx + sc*12.*(frl+frr) / &
                  (gamma*(12.*(frl+frr)**2 + 5.*g13ab*g13ab))
              sl = -frl*(sn-frl1*(sx-(frl+frl1)*sxx))
              sr = sn-sl
           else if(sr.ge.0.) then           ! leftmost division
              sxx = sxx + sl*(frl+frl1)/(frl*frl1*(.6d0+(frl+frl1)**2))
              sx = sn/frl1 + (frl+frl1)*sxx
              sr = frr*(sn-frr1*(sx-(frr+frr1)*sxx))
              sl = 0.
           else if(sl.ge.0.) then           ! rightmost division
              sxx = sxx - sr*(frr+frr1)/(frr*frr1*(.6d0+(frr+frr1)**2))
              sx = sn/frr1 + (frr+frr1)*sxx
              sl = -frl*(sn-frl1*(sx-(frl+frl1)*sxx))
              sr = 0.
           endif
           sc = sn - sl
           sc = sc - sr
! check for the cases where two of the three divisions are nonpositive
! these cases arise due to adjustments made when only one of the three
! divisions was negative
           if(sl.le.0. .and. sr.le.0.) then ! two outer divisions
              gamma = 1.+(frl-frr)
              sxx = sn*(1.+gamma)/(2.*gamma*frl1*frr1)
              sx  = sn*(frl+frr)*(1.+2.*gamma)/(2.*gamma*frl1*frr1)
              sl = 0.
              sr = 0.
           else if(sl.le.0. .and. sc.le.0.) then ! center/left divs
              sxx = sn/(2.*frr*frl1)
              sx  = sn*(frl+frr+.5)/(frr*frl1)
              sl = 0.
              sr = sn
           else if(sr.le.0. .and. sc.le.0.) then ! center/right divs
              sxx = sn/(2.*frl*frr1)
              sx  = sn*(frl+frr-.5)/(frl*frr1)
              sl = sn
              sr = 0.
           endif
           fnm1 = -sl
           fn   = +sr
        else
!**** air is leaving only through one edge
           if(an.gt.0.)  then ! right edge
              fr=an
              sd=fn
              fsign=-1.
           else                  ! left edge
              fr=anm1
              sd=-fnm1
              fsign=1.
           endif
           if(abs(fr).gt.1.)  then
!**** give warnings if fractional mass loss > 1
             ierr=1
             write(6,*) "limitq warning: abs(a)>1",fr,sd
             write(6,*) "limitq input: anm1,an,fnm1,fn,sn,sx,sxx",anm1,an,fnm1,fn,sn,sx,sxx
!**** only stop if net tracer mass is negative
             if (sn+fnm1-fn.lt.0) then
               ierr=2
               write(6,*) "limitq error: new sn < 0",sn,sn+fnm1-fn
               return
             end if
           end if
           su = sn-sd
           if(sd.ge.0. .and. su.ge.0.) return
           fr1=fr+fsign
           if(sd.lt.0.)  then
!**** downstream division is negative
              sxx = sxx +fsign*sd*(fr+fr1)/(fr*fr1*(.6d0+(fr+fr1)**2))
              sx = sn/fr1 + (fr+fr1)*sxx
              su = sn
           else
!**** upstream division is negative
              sxx = sxx -fsign*su*(fr+fr1)/(fr*fr1*(.6d0+(fr+fr1)**2))
              sx = sn/fr + (fr+fr1)*sxx
              su = 0.
           endif
           sd = sn - su
           if(an.gt.0.) then
              fn=sd
           else
              fnm1=-sd
           endif
        endif
        return
      end subroutine limitq

    SUBROUTINE PFIX_QUS( Ps_init, Ps_mid, Ps_end, dT )
      ! Implements the Cameron Smith LLNL pressure fixer into QUS

      USE DIAG_MOD, ONLY     : MASSFLEW, MASSFLNS, MASSFLUP
      USE PRESSURE_MOD, ONLY : dSIG
      IMPLICIT NONE

      REAL*8, INTENT(IN) :: Ps_init(IM,JM) ! Pressure at start of dynamic timestep [hPa]
      REAL*8, INTENT(IN) :: Ps_mid(IM,JM)  ! Pressure at middle of dynamic timestep [hPa]
      REAL*8, INTENT(IN) :: Ps_end(IM,JM)  ! Pressure at end of dynamic timestep [hPa]
      REAL*8, INTENT(IN) :: dT             ! Dynamic timestep [s]

      ! Arrays to hold adjusted fluxes
      REAL*8, DIMENSION(IM,JM,LM)         :: PU_FIXED, PV_FIXED
      REAL*8, DIMENSION(IM,JM,LM), TARGET :: CONV_FIXED
      REAL*8, DIMENSION(:,:,:), POINTER   :: SD_FIXED
      REAL*8, DIMENSION(:,:), POINTER     :: PIT_FIXED

      ! Pressure fixer variables
      !REAL*8, DIMENSION(IM,JM) :: PMET2, PCTM1, DPS_CTM
      !REAL*8                   :: DGPRESS

      REAL*8, DIMENSION(IM,JM) :: DPS_CTM ! Sum over vertical from original mass fluxes [hPa]
      REAL*8, DIMENSION(IM,JM) :: DPS     ! Change of surface pressure from met field pressure [hPa]
      REAL*8, DIMENSION(IM,JM) :: dDPS    ! Difference

      REAL*8, DIMENSION(LM) :: dbk

      ! Scalars
      INTEGER :: i, j, l, im1

      REAL*8  :: dgpress
      REAL*8  :: fxmean
      REAL*8  :: ri2

      ! Arrays
      REAL*8  :: fxintegral(im+1)
      REAL*8  :: mmfd(jm)
      REAL*8  :: mmf(jm)

      ! Only distribute mass in the sigma portion of the vertical grid to avoid upsetting
      ! vertical levels.
      dbk = 0d0
      dbk(1:LS1-1) = dSig(1:LS1-1)

      ! Initialize
      SD_FIXED  => CONV_FIXED(:,:,2:)
      PIT_FIXED => CONV_FIXED(:,:,1)

      ! Pressure tendency predicted from AFLUX (PIT)
      ! Convert from [mb m^2 s^-1] -> [hPa]
      DO J=1,JM
         DPS_CTM(:,J) = PIT(:,J) * dT / DXYP(J)
      ENDDO
      ! Actual pressure tendency from interpolated met fields [hPa]
      DPS     = ( Ps_end - Ps_init )

      !massflew(:,:,1,2) = massflew(:,:,1,2) + DPS_CTM
      !massflew(:,:,2,2) = massflew(:,:,2,2) + DPS

      ! Calculate difference between GCM and CTM predicted pressure
      dDPS    = DPS - DPS_CTM

      ! Calculate global pressure discrepancy
      dgpress = SUM( dDPS * rel_area )

      ! Calculate mean meridional flux divergence (df/dy)
      mmfd(1)    = -( dDPS(1,1)               - dgpress ) ! South Pole
      do J = 2, JM-1
         mmfd(J) = -( sum(dDPS(:,J))/DBLE(IM) - dgpress ) 
      enddo
      mmfd(JM)   = -( dDPS(1,JM)              - dgpress ) ! North Pole

      ! Calculate mean meridional fluxes (cos(e)*fy)
      mmf(2) = mmfd(1) / geofac_pc
      do J = 2, JM-1
         mmf(J+1) = mmf(J) + mmfd(J) / geofac(J)
      enddo

      ! Fix latitude bands
      do J = 2, JM-1
         fxintegral(:) = 0.0d0
         do I = 1, IM
            fxintegral(I+1) = &
                 fxintegral(I) - (dDPS(I,J) - dgpress) - mmfd(J)
         end do
         fxmean = SUM(fxintegral(2:IM+1)) / DBLE(IM)
         do I = 1, IM
            xcolmass_fix(I,J) = fxintegral(I) - fxmean
         end do
      end do
      
      ! Convert flux adjustments from hPa to mb m2 s-1
      do J=1,JM
         xcolmass_fix(:,j) = xcolmass_fix(:,j) * dXYP(j) / dT
         mmf(j) = mmf(j) * (dXYP(j)) / dT
      enddo

      ! Distribute fixed mass in the vertical
      do L = 1, LS1-1      
         do J = 2, JM-1
            do I = 1, IM
               pu_fixed(I,J,L) = pu(I,J,L) + xcolmass_fix(I,J) * dbk(L)
            end do
         end do
      end do

      do L = 1, LS1-1      
         do J = 2, JM
            do I = 1, IM
               pv_fixed(I,J,L) = pv(I,J,L) + mmf(J) * dbk(L)
            end do
         end do
      end do
      
      ! Recalculate pressure tendencies and vertical fluxes for the adjusted PU and PV values,
      ! and see if now consistent with model met
      DO L=1,LM
         DO J=J_0S,J_1S
            IM1=IM
            DO I=1,IM
               CONV_FIXED(I,J,L)=(PU_FIXED(IM1,J,L)-PU_FIXED(I,  J,L)+&
                                  PV_FIXED(  I,J,L)-PV_FIXED(I,J+1,L))
               IM1=I
            ENDDO
         ENDDO
         !CONV(1,1,L)=-PVSA(L) ! Incorporate this later!
         !CONV(1,JM,L)=PVNA(L) ! Incorporate this later!
      ENDDO
      DO J=J_0,J_1
         DO I=1,IMAXJ(J)
            DO L=LM-1,1,-1
               PIT_FIXED(I,J)=PIT_FIXED(I,J)+SD_FIXED(I,J,L)
            ENDDO
         ENDDO
      ENDDO
      DO J=J_0,J_1
         DO I=1,IMAXJ(J)
            DO L=LM-2,LS1-1,-1
               SD_FIXED(I,J,L)=SD_FIXED(I,J,L+1)+SD_FIXED(I,J,L)
            ENDDO
         ENDDO
      ENDDO
      DO J=J_0,J_1
         DO I=1,IMAXJ(J)
            DO L=LS1-2,1,-1
               SD_FIXED(I,J,L)=SD_FIXED(I,J,L+1)+SD_FIXED(I,J,L)-dSIG(L+1)*PIT_FIXED(I,J)
            ENDDO
         ENDDO
      ENDDO
      DO L=1,LM-1
         DO I=2,IM
            SD_FIXED(I, 1,L)=SD_FIXED(1, 1,L)
            SD_FIXED(I,JM,L)=SD_FIXED(1,JM,L)
         ENDDO
      ENDDO

      DO J=1,JM
         DPS_CTM(:,J) = PIT_FIXED(:,J) * dT / DXYP(J)
      ENDDO

      !massflew(:,:,:,1) = pu
      !massflns(:,:,:,1) = pv

      ! To compare with original
      !massflew(:,:,3,2) = massflew(:,:,3,2) + DPS_CTM

      !massflew(:,:,4,2) = massflew(:,:,4,2) + PIT
      !massflew(:,:,5,2) = massflew(:,:,5,2) + PIT_FIXED

      ! Update with corrections
      PU  = 0d0 !PU_FIXED
      PV  = 0d0 !PV_FIXED
      SD  = 0d0 !SD_FIXED
      PIT = 0d0 !PIT_FIXED
      
    END SUBROUTINE PFIX_QUS

    SUBROUTINE AADVT (MA,RM,RMOM,SD,PU,PV,DT,QLIMIT,FQU,FQV,tracer)
!@sum  AADVT advection driver
!@auth G. Russell, modified by Maxwell Kelley
!****
!**** AADVT advects tracers using the Quadradic Upstream Scheme.
!****
!**** input:
!****  pu,pv,sd (kg/s) = east-west,north-south,vertical mass fluxes
!****      qlimit = whether moment limitations should be used
!****         DT (s) = time step
!****
!**** input/output:
!****     rm = tracer concentration
!****   rmom = moments of tracer concentration
!****     ma (kg) = fluid mass
!****
      USE DIAG_MASS_FLUX_QUS, ONLY : DIAG_PU_TR, DIAG_PV_TR, DIAG_SD_TR, DIAG_STE_TR, DIAG_MFLUX_TR_CNT, DIAG_MFLUX
      !USE TROPOPAUSE_MOD, ONLY : GET_TPAUSE_LEVEL
      !USE DIAG_MOD, ONLY     : MASSFLEW, MASSFLNS, MASSFLUP
      !USE QUSDEF
      !USE QUSCOM, ONLY : IM,JM,LM, MFLX
      IMPLICIT NONE

      REAL*8, dimension(im,jm,lm) :: rm,ma,mflx
      REAL*8, dimension(NMOM,IM,JM,LM) :: rmom

      INTEGER, INTENT(IN) :: tracer
      REAL*8, INTENT(IN) :: DT
      REAL*8, dimension(im,jm,lm), intent(in) :: pu,pv
      REAL*8, dimension(im,jm,lm-1), intent(in) :: sd
      LOGICAL, INTENT(IN) :: QLIMIT

      REAL*8, dimension(im,jm), intent(inout) :: fqu,fqv
      REAL*8, DIMENSION(IM,JM,LM) :: hFQU, hFQV
      REAL*8, DIMENSION(IM,JM,LM-1) :: hFSD

      INTEGER :: I,J,L,N
      REAL*8 :: BYMA

!**** Initialise diagnostics
      FQU=0.  ; FQV=0.

! MOVE THIS TO TRANSPORT_MOD
!!**** Fill in values at the poles
!!$OMP  PARALLEL DO PRIVATE(I,L,N)
!      DO L=1,LM
!         DO I=2,IM
!           RM(I,1 ,L) =   RM(1,1 ,L)
!           RM(I,JM,L) =   RM(1,JM,L)
!           DO N=1,NMOM
!             RMOM(N,I,1 ,L) =  RMOM(N,1,1 ,L)
!             RMOM(N,I,JM,L) =  RMOM(N,1,JM,L)
!         enddo
!         enddo
!      enddo
!!$OMP  END PARALLEL DO
!!****
!!**** convert from concentration to mass units
!!****
!!$OMP  PARALLEL DO PRIVATE(I,J,L)
!      DO L=1,LM
!      DO J=1,JM
!      DO I=1,IM
!         RM(I,J,L)=RM(I,J,L)*MA(I,J,L)
!         RMOM(:,I,J,L)=RMOM(:,I,J,L)*MA(I,J,L)
!      enddo
!      enddo
!      enddo
!!$OMP  END PARALLEL DO

!****
!**** Advect the tracer using the quadratic upstream scheme
!****
!C    mflx(:,:,:)=pu(:,:,:)*(.5*dt)
!$OMP  PARALLEL DO PRIVATE(L)
       DO L=1,LM
          mflx(:,:,l)=pu(:,:,l)*(.5*dt)
          mflx(:,1,l)=0d0
          mflx(:,jm,l)=0d0
       ENDDO
!$OMP  END PARALLEL DO
      CALL AADVTX (RM,RMOM,MA,MFLX,QLIMIT,FQU,hFQU)
      IF ( DIAG_MFLUX ) DIAG_PU_TR(:,:,:,tracer) = DIAG_PU_TR(:,:,:,tracer) + hFQU
!C    mflx(:,1:jm-1,:)=pv(:,2:jm,:)*dt
!C    mflx(:,jm,:)=0.
!$OMP  PARALLEL DO PRIVATE(L)
       DO L=1,LM
          mflx(:,1:jm-1,l)=pv(:,2:jm,l)*dt
          mflx(:,jm,l)=0.
       ENDDO
!$OMP  END PARALLEL DO
      CALL AADVTY (RM,RMOM,MA,MFLX,QLIMIT,FQV,hFQV)
      IF ( DIAG_MFLUX ) DIAG_PV_TR(:,:,:,tracer) = DIAG_PV_TR(:,:,:,tracer) + hFQV
!C    mflx(:,:,1:lm-1)=sd(:,:,1:lm-1)*(-dt)
!C    mflx(:,:,lm)=0.
!$OMP  PARALLEL DO PRIVATE(L)
      DO L=1,LM
         IF(L.NE.LM)  THEN
            MFLX(:,:,L)=SD(:,:,L)*(-DT)
         ELSE
            MFLX(:,:,L)=0.
         END IF
      ENDDO
!$OMP  END PARALLEL DO
      CALL AADVTZ (RM,RMOM,MA,MFLX,QLIMIT,hFSD)
      IF ( DIAG_MFLUX ) THEN
         DIAG_SD_TR(:,:,:,tracer) = DIAG_SD_TR(:,:,:,tracer) + hFSD
!         DO I=1,IIPAR
!         DO J=1,JJPAR
!            DIAG_STE_TR(I,J,TRACER) =  DIAG_STE_TR(I,J,TRACER) + hFSD(I,J,GET_TPAUSE_LEVEL(I,J))
!         ENDDO
!         ENDDO
      ENDIF
!C    mflx(:,:,:)=pu(:,:,:)*(.5*dt)
!$OMP  PARALLEL DO PRIVATE(L)
       DO L=1,LM
          mflx(:,:,l)=pu(:,:,l)*(.5*dt)
          mflx(:,1,l)=0d0
          mflx(:,jm,l)=0d0
       ENDDO
!$OMP  END PARALLEL DO
      CALL AADVTX (RM,RMOM,MA,MFLX,QLIMIT,FQU,hFQU)
      IF ( DIAG_MFLUX ) DIAG_PU_TR(:,:,:,tracer) = DIAG_PU_TR(:,:,:,tracer) + hFQU

      IF ( DIAG_MFLUX .and. TRACER .eq. 1 ) DIAG_MFLUX_TR_CNT = DIAG_MFLUX_TR_CNT + 1d0

! MOVED TO TRANSPORT_MOD
!!****
!!**** convert from mass to concentration units
!!****
!!$OMP  PARALLEL DO PRIVATE(I,J,L,BYMA)
!      DO L=1,LM
!      DO J=1,JM
!      DO I=1,IM
!         BYMA = 1.D0/MA(I,J,L)
!         RM(I,J,L)=RM(I,J,L)*BYMA
!         RMOM(:,I,J,L)=RMOM(:,I,J,L)*BYMA
!      enddo
!      enddo
!      enddo
!!$OMP  END PARALLEL DO
      RETURN
    END SUBROUTINE AADVT

      subroutine aadvtx(rm,rmom,mass,mu,qlimit,fqu,hfqu)
!@sum  AADVTX advection driver for x-direction
!@auth Maxwell Kelley
!****
!**** aadvtx advects tracers in the west to east direction using the
!**** quadratic upstream scheme.  if qlimit is true, the moments are
!**** limited to prevent the mean tracer from becoming negative.
!****
!**** input:
!****     mu (kg) = west-east mass flux, positive eastward
!****      qlimit = whether moment limitations should be used
!****
!**** input/output:
!****     rm (kg) = tracer mass
!****   rmom (kg) = moments of tracer mass
!****   mass (kg) = fluid mass
!****
      !use QUSDEF
      !use QUSCOM, only : im,jm,lm, xstride
      implicit none
      REAL*8, dimension(im,jm,lm) :: rm,mass,mu!,hfqu
      REAL*8, dimension(NMOM,IM,JM,LM) :: rmom
      logical ::  qlimit
      REAL*8, INTENT(OUT), DIMENSION(IM,JM) :: FQU
      REAL*8, INTENT(OUT), DIMENSION(IM,JM,LM) :: hFQU
      REAL*8  AM(IM), F_I(IM), FMOM_I(NMOM,IM)
      integer :: i,j,l,ierr,nerr,ICKERR
!**** loop over layers and latitudes
      ICKERR=0
!$OMP PARALLEL DO PRIVATE(J,L,AM,F_I,FMOM_I,IERR,NERR)&
!$OMP DEFAULT( SHARED   )&
!!$OMP SHARED(IM,QLIMIT,XSTRIDE)&
!$OMP REDUCTION(+:ICKERR)
      do l=1,lm
      do j=2,jm-1
      am(:) = mu(:,j,l)
!****
!**** call 1-d advection routine
!****
      call adv1d(rm(1,j,l),rmom(1,1,j,l), f_i,fmom_i, mass(1,j,l),&
           am, im, qlimit,xstride,xdir,ierr,nerr)
      if (ierr.gt.0) then
        write(6,*) "Error in aadvtx: i,j,l=",nerr,j,l
        if (ierr.eq.2) then
          write(0,*) "Error in qlimit: abs(a) > 1"
!CC       call stop_model('Error in qlimit: abs(a) > 1',11)
          ICKERR=ICKERR+1
        end if
      end if
!****
!**** store tracer flux in fqu array
!****
!CC   fqu(:,j)  = fqu(:,j) + f_i(:)
      hfqu(:,j,l)  = f_i(:)
      enddo ! j
      enddo ! l
!$OMP  END PARALLEL DO
!
!     now sum into fqu
!
      do l=1,lm
      do j=2,jm-1
         fqu(:,j)  = fqu(:,j) + hfqu(:,j,l)
      enddo ! j
      enddo ! l
!
      IF(ICKERR.GT.0)  CALL geos_chem_stop !stop_model('Stopped in aadvtx',11)
!
      return
!****
      end subroutine aadvtx

      subroutine aadvty(rm,rmom,mass,mv,qlimit,fqv,hfqv)
!@sum  AADVTY advection driver for y-direction
!@auth Maxwell Kelley
!****
!**** aadvty advects tracers in the south to north direction using the
!**** quadratic upstream scheme.  if qlimit is true, the moments are
!**** limited to prevent the mean tracer from becoming negative.
!****
!**** input:
!****     mv (kg) = north-south mass flux, positive northward
!****      qlimit = whether moment limitations should be used
!****
!**** input/output:
!****     rm (kg) = tracer mass
!****   rmom (kg) = moments of tracer mass
!****   mass (kg) = fluid mass
!****
      !use QUSDEF
      !use QUSCOM, only : im,jm,lm, ystride,               byim
      implicit none
      REAL*8, dimension(im,jm,lm) :: rm,mass,mv
      REAL*8, dimension(NMOM,IM,JM,LM) :: rmom
      logical ::  qlimit
      REAL*8, intent(out), dimension(im,jm) :: fqv
      REAL*8, INTENT(OUT), DIMENSION(IM,JM,LM) :: hFQV
      REAL*8  BM(JM),F_J(JM),FMOM_J(NMOM,JM)
      integer :: i,j,l,ierr,nerr,ICKERR
      REAL*8 :: m_sp,m_np,rm_sp,rm_np,rzm_sp,rzm_np,rzzm_sp,rzzm_np
!**** loop over layers
      ICKERR=0
!$OMP  PARALLEL DO PRIVATE(I,L,M_SP,M_NP,RM_SP,RM_NP,RZM_SP,RZZM_SP,F_J,FMOM_J,RZM_NP,RZZM_NP,BM,IERR,NERR)&
!$OMP DEFAULT( SHARED   )&
!!$OMP SHARED(JM,QLIMIT,YSTRIDE)&
!$OMP REDUCTION(+:ICKERR)
      do l=1,lm
!**** scale polar boxes to their full extent
      mass(:,1:jm:jm-1,l)=mass(:,1:jm:jm-1,l)*im
      m_sp = mass(1,1 ,l)
      m_np = mass(1,jm,l)
      rm(:,1:jm:jm-1,l)=rm(:,1:jm:jm-1,l)*im
      rm_sp = rm(1,1 ,l)
      rm_np = rm(1,jm,l)
      do i=1,im
         rmom(:,i,1 ,l)=rmom(:,i,1 ,l)*im
         rmom(:,i,jm,l)=rmom(:,i,jm,l)*im
      enddo
      rzm_sp  = rmom(mz ,1,1 ,l)
      rzzm_sp = rmom(mzz,1,1 ,l)
      rzm_np  = rmom(mz ,1,jm,l)
      rzzm_np = rmom(mzz,1,jm,l)
!**** loop over longitudes
      do i=1,im
!****
!**** load 1-dimensional arrays
!****
      bm   (:) = mv(i,:,l) !/nstep
      bm(jm)= 0.
      rmom(ihmoms,i,1 ,l) = 0.! horizontal moments are zero at pole
      rmom(ihmoms,i,jm,l) = 0.
!****
!**** call 1-d advection routine
!****
      call adv1d(rm(i,1,l),rmom(1,i,1,l), f_j,fmom_j, mass(i,1,l),&
           bm, jm,qlimit,ystride,ydir,ierr,nerr)
      if (ierr.gt.0) then
        write(6,*) "Error in aadvty: i,j,l=",i,nerr,l
        if (ierr.eq.2) then
          write(0,*) "Error in qlimit: abs(b) > 1"
!cc       call stop_model('Error in qlimit: abs(b) > 1',11)
          ICKERR=ICKERR+1
        endif
      end if
!**** store tracer flux in fqv array
!cc   fqv(i,:) = fqv(i,:) + f_j(:)
!cc   fqv(i,jm) = 0.   ! play it safe
      hfqv(i,:,l) = f_j(:)
      rmom(ihmoms,i,1 ,l) = 0.! horizontal moments are zero at pole
      rmom(ihmoms,i,jm,l) = 0.
      enddo ! end loop over longitudes
!**** average and unscale polar boxes
      mass(:,1 ,l) = (m_sp + sum(mass(:,1 ,l)-m_sp))*byim
      mass(:,jm,l) = (m_np + sum(mass(:,jm,l)-m_np))*byim
      rm(:,1 ,l) = (rm_sp + sum(rm(:,1 ,l)-rm_sp))*byim
      rm(:,jm,l) = (rm_np + sum(rm(:,jm,l)-rm_np))*byim
      rmom(mz ,:,1 ,l) = (rzm_sp  + sum(rmom(mz ,:,1 ,l)-rzm_sp ))*byim
      rmom(mzz,:,1 ,l) = (rzzm_sp + sum(rmom(mzz,:,1 ,l)-rzzm_sp))*byim
      rmom(mz ,:,jm,l) = (rzm_np  + sum(rmom(mz ,:,jm,l)-rzm_np ))*byim
      rmom(mzz,:,jm,l) = (rzzm_np + sum(rmom(mzz,:,jm,l)-rzzm_np))*byim
      enddo ! end loop over levels
!$OMP  END PARALLEL DO
!
!     sum into fqv
!
      do l=1,lm
      do i=1,im
         fqv(i,:) = fqv(i,:) + hfqv(i,:,l)
         fqv(i,jm) = 0.
      enddo
      enddo
!
      IF(ICKERR.GT.0)  CALL geos_chem_stop !stop_model('Stopped in aadvty',11)
!
      return
!****
      end subroutine aadvty

      subroutine aadvtz(rm,rmom,mass,mw,qlimit,hfsd)
!@sum  AADVTZ advection driver for z-direction
!@auth Maxwell Kelley
!****
!**** aadvtz advects tracers in the upward vertical direction using the
!**** quadratic upstream scheme.  if qlimit is true, the moments are
!**** limited to prevent the mean tracer from becoming negative.
!****
!**** input:
!****     mw (kg) = vertical mass flux, positive upward
!****      qlimit = whether moment limitations should be used
!****
!**** input/output:
!****     rm (kg) = tracer mass
!****   rmom (kg) = moments of tracer mass
!****   mass (kg) = fluid mass
!****
      !use QUSDEF
      !use QUSCOM, only : im,jm,lm, zstride
      implicit none
      REAL*8, dimension(im,jm,lm) :: rm,mass,mw
      REAL*8, dimension(NMOM,IM,JM,LM) :: rmom
      logical ::  qlimit
      REAL*8  CM(LM),F_L(LM),FMOM_L(NMOM,LM)
      REAL*8, DIMENSION(IM,JM,LM-1), INTENT(OUT) :: hfsd
      integer :: i,j,l,ierr,nerr,ICKERR
!**** loop over latitudes and longitudes
      ICKERR=0
!$OMP PARALLEL DO PRIVATE(I,J,CM,F_L,FMOM_L,IERR,NERR)&
!$OMP DEFAULT( SHARED   )&
!!$OMP SHARED(LM,QLIMIT,ZSTRIDE)&
!$OMP REDUCTION(+:ICKERR)
      do j=1,jm
      do i=1,im
      cm(:) = mw(i,j,:)
      cm(lm)= 0.
!****
!**** call 1-d advection routine
!****
      call adv1d(rm(i,j,1),rmom(1,i,j,1),f_l,fmom_l,mass(i,j,1),&
           cm,lm,qlimit,zstride,zdir,ierr,nerr)
      if (ierr.gt.0) then
        write(6,*) "Error in aadvtz: i,j,l=",i,j,nerr
        if (ierr.eq.2) then
          write(0,*) "Error in qlimit: abs(c) > 1"
!cc       call stop_model('Error in qlimit: abs(c) > 1',11)
          ICKERR=ICKERR+1
        endif
      end if
      hfsd(i,j,:) = hfsd(i,j,:) + F_L(1:(LM-1))
      enddo ! i
      enddo ! j
!$OMP  END PARALLEL DO
!
      IF(ICKERR.GT.0) call geos_chem_stop !stop_model('Stopped in aadvtz',11)
      return
!****
      end subroutine aadvtz

  end module qus_mod
