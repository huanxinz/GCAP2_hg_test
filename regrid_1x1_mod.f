! $Id: regrid_1x1_mod.f,v 1.10 2007/02/22 18:26:29 bmy Exp $
      MODULE REGRID_1x1_MOD
!
!******************************************************************************
!  Module REGRID_1x1_MOD does online regridding of data on the GEOS-Chem 1x1 
!  grid to 1x1, 2x25, or 4x5 GEOS/GCAP grids. (bdf, bmy, 10/24/05, 2/20/07)
!
!  Module Variables:
!  ============================================================================
!  (1 ) A1x1      (REAL*8)      : Surface areas [m2] of 1x1 GEOS-Chem grid 
!  (2 ) A_GEN_1x1 (REAL*8)      : Surface areas [m2] of 1x1 GENERIC   grid 
!
!  Module Routines:
!  ============================================================================
!  (1 ) DO_REGRID_G2G_1x1       : Regrids GENERIC 1x1 to GEOS-Chem 1x1 GRID
!  (2 ) DO_REGRID_1x1_R4        : Passes 3D, REAL*4 array to DO_THE_REGRIDDING
!  (3 ) DO_REGRID_1x1_R4_2D     : Passes 2D, REAL*4 array to DO_THE_REGRIDDING
!  (4 ) DO_REGRID_1x1_R8        : Passes 3D, REAL*8 array to DO_THE_REGRIDDING
!  (5 ) DO_REGRID_1x1_R8_2D     : Passes 2D, REAL*8 input to DO_THE_REGRIDDING
!  (6 ) DO_THE_REGRIDDING       : Driver routine for regridding from 1x1
!  (7 ) ITS_CONCENTRATION_DATA  : Returns TRUE if it's concentration data
!  (8 ) REGRID_CONC_TO_4x5_GCAP : Regrids conc data from GEOS 1x1 -> GCAP 4x5
!  (9 ) REGRID_MASS_TO_4x5_GCAP : Regrids mass data from GEOS 1x1 -> GCAP 4x5
!  (10) REGRID_CONC_TO_4x5      : Regrids conc data from GEOS 1x1 -> GEOS 4x5
!  (11) REGRID_MASS_TO_4x5      : Regrids mass data from GEOS 1x1 -> GEOS 4x5
!  (12) REGRID_CONC_TO_2x25     : Regrids conc data from GEOS 1x1 -> GEOS 2x25
!  (13) REGRID_MASS_TO_2x25     : Regrids mass data from GEOS 1x1 -> GEOS 2x25
!  (14) REGRID_CONC_TO_1x125    : Regrids conc data from GEOS 1x1 -> GEOS 1x125
!  (15) REGRID_MASS_TO_1x125    : Regrids mass data from GEOS 1x1 -> GEOS 1x125
!  (16) INIT_REGRID_1x1         : Initializes all module variables
!  (17) CLEANUP_REGRID_1x1      : Deallocates all module variables
! 
!  GEOS-Chem modules referenced by "regrid_1x1_mod.f"
!  ============================================================================
!  (1 ) charpak_mod.f           : Module w/ string handling routines
!  (2 ) grid_mod.f              : Module w/ horizontal grid information
!
!  NOTES:  
!  (1 ) Added DO_REGRID_G2G_1x1 to regrid from GENERIC 1x1 to GEOS 1x1 grid.
!        (psk, bmy, 4/18/06)
!  (2 ) Added routines REGRID_CONC_TO_1x125 and REGRID_MASS_TO_1x125 to regrid
!        1x1 data to the GEOS-Chem 1x1.25 grid. (bdf, bmy, 8/2/06)
!  (3 ) DO_REGRID_G2G_1x1 now takes UNIT via the arg list (bmy, 8/9/06)
!  (4 ) Bug fix in REGRID_MASS_TO_4x5 (tw, bmy, 2/20/07)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "regrid_1x1_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE
 
      ! ... except these routines
      PUBLIC :: CLEANUP_REGRID_1x1
      PUBLIC :: DO_REGRID_1x1
      PUBLIC :: DO_REGRID_G2G_1x1

      !=================================================================
      ! MODULE VARIABLES 
      !=================================================================

      ! Arrays
      REAL*8, ALLOCATABLE :: A1x1(:)
      REAL*8, ALLOCATABLE :: A_GEN_1x1(:)

      !=================================================================
      ! MODULE INTERFACES -- "bind" two or more routines with different
      ! argument types or # of arguments under one unique name
      !================================================================= 
      INTERFACE DO_REGRID_1x1
         MODULE PROCEDURE DO_REGRID_1x1_R4
         MODULE PROCEDURE DO_REGRID_1x1_R4_2D
         MODULE PROCEDURE DO_REGRID_1x1_R8
         MODULE PROCEDURE DO_REGRID_1x1_R8_2D
      END INTERFACE

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE DO_REGRID_G2G_1x1( UNIT, GEN1x1, GEOS1x1 )
!
!******************************************************************************
!  Subroutine DO_REGRID_G2G_1x1 regrids 2-D data on the GENERIC 1x1
!  grid (1st box edged at -180, -90) to the GEOS-Chem 1x1 grid. 
!  (psk, bmy, 4/5/06, 8/9/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) UNIT    (CHARACTER) : Unit of the data to be regridded 
!  (2 ) GEN1x1  (REAL*4   ) : Data array on the GENERIC 1x1 grid
!
!  Arguments as Output:
!  ============================================================================
!  (3 ) GEOS1x1 (REAL*4   ) : Data array on the GEOS 1x1 grid
! 
!  NOTES:
!  (1 ) Now pass UNIT via the arg list and pass that to ITS_CONCENTRATION_DATA
!        to determine if the data to be regridded is concentration data
!        or mass data.  This is now consistent with routine DO_REGRID_1x1.
!        (bmy, 8/9/06)
!******************************************************************************
!
#     include "CMN_SIZE"             ! Size parameters
#     include "CMN_GCTM"             ! Physical constants

      ! Arguments
      CHARACTER(LEN=*), INTENT(IN)  :: UNIT
      REAL*8,           INTENT(IN)  :: GEN1x1(I1x1,J1x1-1)
      REAL*8,           INTENT(OUT) :: GEOS1x1(I1x1,J1x1)

      ! Local variables
      LOGICAL, SAVE                 :: FIRST = .TRUE.
      LOGICAL                       :: ITS_PER_UNIT_AREA
      INTEGER                       :: I, J, IX, JX, IE(2), JE(2)
      REAL*8                        :: RE_cm, LAT, DLAT2
      REAL*8                        :: A_GEN(J1x1-1)
      REAL*8                        :: A_GEOS(J1x1)

      !=================================================================
      ! DO_REGRID_G2G_1x1 begins here!      
      !=================================================================

      ! Initialize on first call (if necessary)
      IF ( FIRST ) THEN
         CALL INIT_REGRID_1x1
         FIRST = .FALSE.
      ENDIF

      ! Is this concentration data?  
      ITS_PER_UNIT_AREA = ITS_CONCENTRATION_DATA( UNIT )

      ! Surface area on generic grid [m2]
      A_GEN(:)  = A_GEN_1x1(:)

      ! Surface area on GEOS-Chem grid [m2]
      A_GEOS(:) = A1x1(:)

      !-----------------------------------
      ! Regrid quantity in mass units
      ! from GENERIC to GEOS-Chem grid
      !-----------------------------------
     
      ! Loop over GEOS-Chem latitudes
      DO J = 1, J1x1

         ! Set limits
         JE(1) = J - 1
         JE(2) = J

         ! Special case for South Pole
         IF ( J == 1 ) THEN
            JE(1) = 1
            JE(2) = 1
         ENDIF

         ! Special case for North Pole
         IF ( J == J1x1 ) THEN
            JE(1) = J1x1-1
            JE(2) = J1x1-1
         ENDIF

         ! Loop over GEOS-Chem longitudes
         DO I = 1, I1x1

            ! Zero quantity on GEOS-Chem 1x1 GRID
            GEOS1x1(I,J) = 0d0

            ! Set limits
            IE(1) = I - 1
            IE(2) = I

            ! Date line
            IF ( I == 1 ) THEN
               IE(1) = I1x1
               IE(2) = 1
            ENDIF
                
            ! Save into GEOS 1x1 grid
            IF ( ITS_PER_UNIT_AREA ) THEN

               ! Data on GENERIC 1x1 grid is per unit area
               ! We have to multiply by the generic grid area (A_GEN)
               DO JX = 1, 2
               DO IX = 1, 2
                  GEOS1x1(I,J) = GEOS1x1(I,J) +
     &               0.25d0 * GEN1x1( IE(IX), JE(JX) ) * A_GEN( JE(JX) )
               ENDDO 
               ENDDO

            ELSE

               ! Data on GENERIC 1x1 grid is a mass quantity
               ! We do not have to multiply by the generic grid area
               DO JX = 1, 2
               DO IX = 1, 2
                  GEOS1x1(I,J) = GEOS1x1(I,J) +
     &               0.25d0 * GEN1x1( IE(IX), JE(JX) ) 
               ENDDO 
               ENDDO 

            ENDIF

         ENDDO
      ENDDO

      ! If the data on the GENERIC 1x1 grid is per unit area...we also
      ! want to return data on the GEOS 1x1 grid as per unit area.
      ! Thus, we have to divide by the GEOS 1x1 area (A_GEOS).
      IF ( ITS_PER_UNIT_AREA ) THEN
         DO J = 1, J1x1
         DO I = 1, I1x1
            GEOS1x1(I,J) = GEOS1x1(I,J) / A_GEOS(J)
         ENDDO
         ENDDO
      ENDIF

      ! Return to calling program
      END SUBROUTINE DO_REGRID_G2G_1x1

!------------------------------------------------------------------------------

      SUBROUTINE DO_REGRID_1x1_R4( L1x1, UNIT, INDATA, OUTDATA )
!
!******************************************************************************
!  Subroutine DO_REGRID_1x1_R4 is a wrapper routine for DO_THE_REGRIDDING.
!  It takes a REAL*4 array as input and returns a 3-D REAL*8 array as output. 
!  (bdf, bmy, 10/24/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) L1x1    (INTEGER ) : Level dimension for INDATA and OUTDATA
!  (2 ) UNIT    (CHAR*(*)) : String containing the units of INDATA & OUTDATA
!  (3 ) INDATA  (REAL*4  ) : Input data array on 1x1 grid
!  
!  Arguments as Output:
!  ============================================================================
!  (4 ) OUTDATA (REAL*8  ) : Output data array 
!
!  NOTES:
!******************************************************************************
!
#     include "CMN_SIZE"             ! Size parameters
      
      ! Arguments
      INTEGER,          INTENT(IN)  :: L1x1
      REAL*4,           INTENT(IN)  :: INDATA(I1x1,J1x1,L1x1)
      REAL*8,           INTENT(OUT) :: OUTDATA(IIPAR,JJPAR,L1x1)
      CHARACTER(LEN=*), INTENT(IN)  :: UNIT

      !=================================================================
      ! DO_REGRID_1x1_R4 begins here
      !=================================================================
      
      ! Regrid data
      CALL DO_THE_REGRIDDING( L1x1, UNIT, DBLE( INDATA ), OUTDATA )

      ! Return to calling program
      END SUBROUTINE DO_REGRID_1x1_R4

!------------------------------------------------------------------------------

      SUBROUTINE DO_REGRID_1x1_R8( L1x1, UNIT, INDATA, OUTDATA )
!
!******************************************************************************
!  Subroutine DO_REGRID_1x1_R8 is a wrapper routine for DO_THE_REGRIDDING.
!  It takes a REAL*8 array as input and returns a 3-D REAL*8 array as output.
!  (bdf, bmy, 10/24/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) L1x1    (INTEGER ) : Level dimension for INDATA and OUTDATA
!  (2 ) UNIT    (CHAR*(*)) : String containing the units of INDATA & OUTDATA
!  (3 ) INDATA  (REAL*8  ) : Input data array on 1x1 grid
!  
!  Arguments as Output:
!  ============================================================================
!  (4 ) OUTDATA (REAL*8  ) : Output data array 
!
!  NOTES:
!******************************************************************************
!
#     include "CMN_SIZE"             ! Size parameters
      
      ! Arguments
      INTEGER,          INTENT(IN)  :: L1x1
      REAL*8,           INTENT(IN)  :: INDATA(I1x1,J1x1,L1x1)
      REAL*8,           INTENT(OUT) :: OUTDATA(IIPAR,JJPAR,L1x1)
      CHARACTER(LEN=*), INTENT(IN)  :: UNIT

      !=================================================================
      ! DO_REGRID_1x1_R8 begins here
      !=================================================================

      ! Regrid data
      CALL DO_THE_REGRIDDING( L1x1, UNIT, INDATA, OUTDATA )

      ! Return to calling program
      END SUBROUTINE DO_REGRID_1x1_R8

!------------------------------------------------------------------------------

      SUBROUTINE DO_REGRID_1x1_R4_2D( UNIT, INDATA, OUTDATA )
!
!******************************************************************************
!  Subroutine DO_REGRID_1x1_R4 is a wrapper routine for DO_THE_REGRIDDING.
!  It takes a REAL*4 array as input and saves a 2-D REAL*8 array as output.
!  (bdf, bmy, 10/24/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) L1x1    (INTEGER ) : Level dimension for INDATA and OUTDATA
!  (2 ) UNIT    (CHAR*(*)) : String containing the units of INDATA & OUTDATA
!  (3 ) INDATA  (REAL*4  ) : Input data array on 1x1 grid
!  
!  Arguments as Output:
!  ============================================================================
!  (4 ) OUTDATA (REAL*8  ) : Output data array 
!
!  NOTES:
!******************************************************************************
!
#     include "CMN_SIZE"             ! Size parameters
      
      ! Arguments
      REAL*4,           INTENT(IN)  :: INDATA(I1x1,J1x1,1)
      REAL*8,           INTENT(OUT) :: OUTDATA(IIPAR,JJPAR)
      CHARACTER(LEN=*), INTENT(IN)  :: UNIT

      ! Local variables
      REAL*8                        :: TMP_OUT(IIPAR,JJPAR,1)

      !=================================================================
      ! DO_REGRID_1x1_R4 begins here
      !=================================================================

      ! Regrid data
      CALL DO_THE_REGRIDDING( 1, UNIT, DBLE( INDATA ), TMP_OUT )

      ! Save output data to a 2D array 
      OUTDATA(:,:) = TMP_OUT(:,:,1)

      ! Return to calling program
      END SUBROUTINE DO_REGRID_1x1_R4_2D

!------------------------------------------------------------------------------

      SUBROUTINE DO_REGRID_1x1_R8_2D( UNIT, INDATA, OUTDATA )
!
!******************************************************************************
!  Subroutine DO_REGRID_1x1_R8_2D is a wrapper routine for DO_THE_REGRIDDING.
!  It takes a REAL*8 array as input and saves to a 2-D REAL*8 array as output.
!  (bdf, bmy, 10/24/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) L1x1    (INTEGER ) : Level dimension for INDATA and OUTDATA
!  (2 ) UNIT    (CHAR*(*)) : String containing the units of INDATA & OUTDATA
!  (3 ) INDATA  (REAL*8  ) : Input data array on 1x1 grid
!  
!  Arguments as Output:
!  ============================================================================
!  (4 ) OUTDATA (REAL*8  ) : Output data array 
!
!  NOTES:
!******************************************************************************
!
#     include "CMN_SIZE"             ! Size parameters
      
      ! Arguments
      REAL*8,           INTENT(IN)  :: INDATA(I1x1,J1x1,1)
      REAL*8,           INTENT(OUT) :: OUTDATA(IIPAR,JJPAR)
      CHARACTER(LEN=*), INTENT(IN)  :: UNIT

      ! Local variables
      REAL*8                        :: TMP_OUT(IIPAR,JJPAR,1)

      !=================================================================
      ! DO_REGRID_1x1_R8 begins here
      !=================================================================

      ! Regrid data
      CALL DO_THE_REGRIDDING( 1, UNIT, INDATA, TMP_OUT )

      ! Copy output data to a 2D array 
      OUTDATA(:,:) = TMP_OUT(:,:,1) 

      ! Return to calling program
      END SUBROUTINE DO_REGRID_1x1_R8_2D

!------------------------------------------------------------------------------

      SUBROUTINE DO_THE_REGRIDDING( L1x1, UNIT, INDATA, OUTDATA )
!
!******************************************************************************
!  Subroutine DO_THE_REGRIDDING is the driver routine for the regridding from 
!  the GEOS-Chem 1x1 grid to other CTM grids. (bmy, 10/24/05, 8/2/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) L1x1    (INTEGER ) : Level dimension for INDATA and OUTDATA
!  (2 ) UNIT    (CHAR*(*)) : String containing the units of INDATA & OUTDATA
!  (3 ) INDATA  (REAL*8  ) : Input data array on 1x1 grid
!  
!  Arguments as Output:
!  ============================================================================
!  (4 ) OUTDATA (REAL*8  ) : Output data array 
!
!  NOTES:
!  (1 ) Added #if block for 1 x 1.25 grid (bdf, bmy, 8/2/06)
!******************************************************************************
!
#     include "CMN_SIZE"             ! Size parameters
      
      ! Arguments
      INTEGER,          INTENT(IN)  :: L1x1
      REAL*8,           INTENT(IN)  :: INDATA(I1x1,J1x1,L1x1)
      REAL*8,           INTENT(OUT) :: OUTDATA(IIPAR,JJPAR,L1x1)
      CHARACTER(LEN=*), INTENT(IN)  :: UNIT

      ! Local variables
      LOGICAL, SAVE                 :: FIRST = .TRUE.
      LOGICAL                       :: IS_CONC

      !=================================================================
      ! DO_REGRID_1x1 begins here!
      !=================================================================

      ! Initialize on first call (if necessary)
      IF ( FIRST ) THEN
         CALL INIT_REGRID_1x1
         FIRST = .FALSE.
      ENDIF

      ! Is this concentration data?
      IS_CONC = ITS_CONCENTRATION_DATA( UNIT )

#if   defined( GCAP ) 

      !--------------------------------------------
      ! Regrid GEOS 1x1 grid to GCAP 4x5 grid
      !--------------------------------------------

      IF ( IS_CONC ) THEN

         ! Regrid concentration field to GCAP 4x5
         CALL REGRID_CONC_TO_4x5_GCAP( I1x1,  J1x1, L1x1, INDATA,
     &                                 IIPAR, JJPAR,      OUTDATA )

      ELSE

         ! Regrid mass field to GCAP 4x5
         CALL REGRID_MASS_TO_4x5_GCAP( I1x1,  J1x1, L1x1, INDATA,
     &                                 IIPAR, JJPAR,      OUTDATA )

      ENDIF

#elif defined( GRID4x5 )

      !--------------------------------------------
      ! Regrid GEOS 1x1 grid to GEOS 4x5 grid
      !--------------------------------------------

      IF ( IS_CONC ) THEN

         ! Regrid concentration field to 4x5
         CALL REGRID_CONC_TO_4x5( I1x1,  J1x1, L1x1, INDATA,
     &                            IIPAR, JJPAR,      OUTDATA )

      ELSE

         ! Regrid mass field to 4x5
         CALL REGRID_MASS_TO_4x5( I1x1,  J1x1, L1x1, INDATA,
     &                            IIPAR, JJPAR,      OUTDATA )

      ENDIF

#elif defined( GRID2x25 )

      !-------------------------------------------
      ! Regrid GEOS 1x1 grid to GEOS 2x2.5 grid
      !-------------------------------------------

      IF ( IS_CONC ) THEN

         ! Regrid concentration field to 4x5
         CALL REGRID_CONC_TO_2x25( I1x1,  J1x1, L1x1, INDATA,
     &                             IIPAR, JJPAR,      OUTDATA )

      ELSE

         ! Regrid mass field to 4x5
         CALL REGRID_MASS_TO_2x25( I1x1,  J1x1, L1x1, INDATA,
     &                             IIPAR, JJPAR,      OUTDATA )

      ENDIF

#elif defined( GRID1x125 )

      !--------------------------------------------
      ! Regrid GEOS 1x1 grid to GEOS 1x1.25 GRID
      !--------------------------------------------

      IF ( IS_CONC ) THEN

         ! Regrid concentration field to 1x125
         CALL REGRID_CONC_TO_1X125( I1x1,  J1x1,  L1x1, INDATA,
     &                              IIPAR, JJPAR,       OUTDATA )

      ELSE

         ! Regrid mass field to 1x125
         CALL REGRID_MASS_TO_1X125( I1x1,  J1x1,  L1x1, INDATA,
     &                              IIPAR, JJPAR,       OUTDATA )

      ENDIF

#elif defined( GRID1x1 ) && defined( NESTED_CH )

      !--------------------------------------------
      ! Regrid GEOS 1x1 grid to nested China grid
      !--------------------------------------------

      ! China nested grid has corners (70E,11S) and (150E,55N)
      ! which corresponds to 1x1 indices (251,80) and (331,146)
      OUTDATA = INDATA( 251:331, 80:146, : )

#elif defined( GRID1x1 ) && defined( NESTED_NA )

      !--------------------------------------------
      ! Regrid GEOS 1x1 grid to nested N. Am. grid
      !--------------------------------------------

      ! N. Am. nested grid has corners (10N,140W) and (60N,40W)
      ! which corresponds to 1x1 indices (41,101) and (141,151)
      OUTDATA = INDATA( 41:141, 101:151, : )

#elif defined( GRID1x1 )

      !--------------------------------------------
      ! GEOS 1x1 grid (no regridding necessary)
      !--------------------------------------------

      ! Copy data array
      OUTDATA = INDATA

#endif

      ! Return to calling program
      END SUBROUTINE DO_THE_REGRIDDING

!------------------------------------------------------------------------------

      FUNCTION ITS_CONCENTRATION_DATA( UNIT ) RESULT( IS_CONC )
!
!******************************************************************************
!  Subroutine ITS_CONCENTRATION_DATA returns TRUE if UNIT is a concentration
!  (i.e. is per unit area such as molec/cm2/s or is a ratio such as kg/kg).
!  (bmy, 10/24/05, 8/9/06)
!  
!  Arguments as Input:
!  ============================================================================
!  (1 ) UNIT (CHAR*(*)) : String with unit of data
!
!  NOTES:
!  (1 ) Added kg/s, kg/month, kg/season to CASE statement (bmy, 8/9/06)
!******************************************************************************
!
      ! References to F90 modules
      USE CHARPAK_MOD, ONLY : STRSQUEEZE
      USE ERROR_MOD,   ONLY : ERROR_STOP

      ! Arguments
      CHARACTER(LEN=* )    :: UNIT

      ! Local variables
      LOGICAL              :: IS_CONC
      CHARACTER(LEN=40)    :: THISUNIT
      CHARACTER(LEN=255)   :: MSG, LOC
      
      !=================================================================
      ! ITS_CONCENTRATION_DATA begins here!
      !=================================================================

      ! Copy UNIT to local variable
      THISUNIT = UNIT

      ! Remove all leading/trailing blanks
      CALL STRSQUEEZE( THISUNIT )

      ! Test if UNIT is a concentration unit (i.e. per unit area or a ratio)
      SELECT CASE ( TRIM( THISUNIT ) )     

         ! Concentration units
         CASE ( 'gC/m2/s'      )
            IS_CONC = .TRUE.
         CASE ( 'unitless'     ) 
            IS_CONC = .TRUE.
         CASE ( 'molec/cm2'    )
            IS_CONC = .TRUE.
         CASE ( 'molec/cm2/s'  )
            IS_CONC = .TRUE.
         CASE ( 'molec C/cm2/s')
            IS_CONC = .TRUE.
         CASE ( 'atom C/cm2/s' )
            IS_CONC = .TRUE.
         CASE ( 'atoms C/cm2/s' )
            IS_CONC = .TRUE.
         CASE ( 's-1'          )
            IS_CONC = .TRUE.
         CASE ( 'K'            ) 
            IS_CONC = .TRUE.
         CASE ( 'kg/kg'        )
            IS_CONC = .TRUE.
         CASE ( 'factor'       )
            IS_CONC = .TRUE.
         CASE ( 'm2/m2'        )
            IS_CONC = .TRUE.
         CASE ( 'cm2/cm2'      )
            IS_CONC = .TRUE.
         CASE ( 'DU'           )
            IS_CONC = .TRUE.
         CASE ( 'DU/day'       )
            IS_CONC = .TRUE.
         CASE ( 'mg C/m2/hr'   )
            IS_CONC = .TRUE.
         CASE ( 'ug C/m2/hr'   )
            IS_CONC = .TRUE.

         ! Mass units
         CASE ( 'kg'           )
            IS_CONC = .FALSE.
         CASE ( 'kg/s'         )
            IS_CONC = .FALSE.
         CASE ( 'kg/month'     )
            IS_CONC = .FALSE.
         CASE ( 'kg/season'    )
            IS_CONC = .FALSE.
         CASE ( 'kg/yr'        )
            IS_CONC = .FALSE.
         CASE ( 'kgC/yr'       )
            IS_CONC = .FALSE.
         CASE ( 'kg C/yr'      )
            IS_CONC = .FALSE.
         CASE ( 'kgN'          )
            IS_CONC = .FALSE.
         CASE ( 'kgEC'         )
            IS_CONC = .FALSE.
         CASE ( 'kgOC'         )
            IS_CONC = .FALSE.
         CASE ( 'kgEC/yr'      )
            IS_CONC = .FALSE.
         CASE ( 'kgOC/yr'      )
            IS_CONC = .FALSE.

         ! Unit not recognized
         CASE DEFAULT

            ! Set IS_CONC to false
            IS_CONC = .FALSE.

            ! Error msg
            MSG = TRIM( UNIT ) // ' is an unrecognized unit, ' //
     &            ' it must be added to the CASE statement!'
           
            ! Location of error
            LOC = 'IS_CONCENTRATION_DATA ("regrid_mod.f")'

            ! Stop run w/ error
            CALL ERROR_STOP( MSG, LOC )

      END SELECT

      ! Return to calling program
      END FUNCTION ITS_CONCENTRATION_DATA

!------------------------------------------------------------------------------

      SUBROUTINE REGRID_CONC_TO_4x5_GCAP( I1, J1, L1, IN, I4, J4, OUT )
!
!******************************************************************************
!  Subroutine REGRID_CONC_TO_4x5_GCAP regrids concentration data from the 
!  GEOS-Chem 1x1 grid to the GEOS-Chem 4x5 GCAP grid. (bdf, bmy, 10/24/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I1  (INTEGER) : 1x1 longitude dimension of IN array
!  (2 ) J1  (INTEGER) : 1x1 latitude  dimension of IN array
!  (3 ) L1  (INTEGER) : 1x1 altitude  dimension of IN array
!  (4 ) IN  (REAL*8 ) : Array containing input data on GEOS-Chem 1x1 grid
!  (5 ) I4  (INTEGER) : 4x5 longitude dimension of OUT array
!  (6 ) J4  (INTEGER) : 4x5 latitude  dimension of OUT array
!
!  Arguments as Output:
!  ============================================================================
!  (7 ) OUT (REAL*8 ) : Array containing output data on GCAP 4x5 grid
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_AREA_M2

      ! Arguments
      INTEGER, INTENT(IN)  :: I1, J1, L1, I4, J4
      REAL*8,  INTENT(IN)  :: IN(I1,J1,L1)
      REAL*8,  INTENT(OUT) :: OUT(I4,J4,L1)

      ! Local variables
      INTEGER              :: I, J, L, W, E, S, N
      REAL*8               :: M_TOT
      
      !==================================================================
      ! REGRID_CONC_TO_4x5_GCAP begins here!
      !==================================================================

      ! Loop over levels
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, W, E, M_TOT, S, N )
      DO L = 1, L1

         !-----------------------
         ! S and N Poles
         !-----------------------
         DO I = 1, I4

            ! 1x1 lon index at W edge of 4x5 box
            W           = MOD( 5 * ( I - 1 ) - 1 + I1, I1 )
 
            ! 1x1 lon index at E edge of 4x5 box
            E           = 5 * ( I - 1 ) + 3

            ! Total mass of 1x1 boxes w/in to 4x5 S pole box
            M_TOT       = SUM( IN( W  :W+1, 1,    L ) ) * A1x1(1) + 
     &                    SUM( IN( W  :W+1, 2,    L ) ) * A1x1(2) + 
     &                    SUM( IN( W  :W+1, 3,    L ) ) * A1x1(3) + 
     &                    SUM( IN( W  :W+1, 4,    L ) ) * A1x1(4) + 
     &              0.5d0*SUM( IN( W  :W+1, 5,    L ) ) * A1x1(5) + 
     &                    SUM( IN( E-2:E,   1,    L ) ) * A1x1(1) +
     &                    SUM( IN( E-2:E,   2,    L ) ) * A1x1(2) + 
     &                    SUM( IN( E-2:E,   3,    L ) ) * A1x1(3) + 
     &                    SUM( IN( E-2:E,   4,    L ) ) * A1x1(4) + 
     &              0.5d0*SUM( IN( E-2:E,   5,    L ) ) * A1x1(5)

            ! Output field at 4x5 S pole box
            OUT(I,1,L)  = M_TOT / ( 5d0 * ( A1x1(1) + A1x1(2) +
     &                                      A1x1(3) + A1x1(4) + 
     &                                        0.5d0 * A1x1(5) ) )

            ! Total mass of 1x1 boxes w/in to 4x5 N pole box
            M_TOT       = SUM( IN( W  :W+1, J1,    L ) ) * A1x1(J1  ) + 
     &                    SUM( IN( W  :W+1, J1-1,  L ) ) * A1x1(J1-1) + 
     &                    SUM( IN( W  :W+1, J1-2,  L ) ) * A1x1(J1-2) + 
     &                    SUM( IN( W  :W+1, J1-3,  L ) ) * A1x1(J1-3) + 
     &              0.5d0*SUM( IN( W  :W+1, J1-4,  L ) ) * A1x1(J1-4) + 
     &                    SUM( IN( E-2:E,   J1,    L ) ) * A1x1(J1  ) +
     &                    SUM( IN( E-2:E,   J1-1,  L ) ) * A1x1(J1-1) + 
     &                    SUM( IN( E-2:E,   J1-2,  L ) ) * A1x1(J1-2) + 
     &                    SUM( IN( E-2:E,   J1-3,  L ) ) * A1x1(J1-3) + 
     &              0.5d0*SUM( IN( E-2:E,   J1-4,  L ) ) * A1x1(J1-4)

            ! Output field at 4x5 N pole box
            OUT(I,J4,L) = M_TOT / ( 5d0 * ( A1x1(J1)   + A1x1(J1-1) +
     &                                      A1x1(J1-2) + A1x1(J1-3) +
     &                                           0.5d0 * A1x1(J1-4) ) )
         ENDDO

         !-----------------------
         ! Non-polar latitudes
         !-----------------------
         DO J = 2, J4-1 

            ! 1x1 lat index at S edge of 4x5 box
            S     = ( 4 * ( J - 1 ) ) + 1

            ! 1x1 lat index at N edge of 4x5 box
            N     = ( J * 4 ) + 1
      
         DO I = 1, I4
          
            ! 1x1 lon index at W edge of 4x5 box
            W     = MOD( 5*( I - 1 ) - 1 + I1, I1 )

            ! 1x1 lon index at E edge of 4x5 box
            E     = 5*( I -1 ) + 3

            ! Total mass w/in the 4x5 box at (I,J,L)
            M_TOT = 0.5d0 * SUM( IN( W  :W+1, S,   L ) ) * A1x1(S  ) + 
     &              0.5d0 * SUM( IN( E-2:E,   S,   L ) ) * A1x1(S  ) + 
     &                      SUM( IN( W  :W+1, S+1, L ) ) * A1x1(S+1) + 
     &                      SUM( IN( E-2:E,   S+1, L ) ) * A1x1(S+1) + 
     &                      SUM( IN( W  :W+1, S+2, L ) ) * A1x1(S+2) + 
     &                      SUM( IN( E-2:E,   S+2, L ) ) * A1x1(S+2) + 
     &                      SUM( IN( W  :W+1, S+3, L ) ) * A1x1(S+3) + 
     &                      SUM( IN( E-2:E,   S+3, L ) ) * A1x1(S+3) + 
     &              0.5d0 * SUM( IN( W  :W+1, N,   L ) ) * A1x1(N  ) +
     &              0.5d0 * SUM( IN( E-2:E,   N,   L ) ) * A1x1(N  )
            
            ! 4x5 output field at (I,J,L)
            OUT(I,J,L) = M_TOT / GET_AREA_M2( J )
         ENDDO 
         ENDDO  

      ENDDO 
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE REGRID_CONC_TO_4x5_GCAP

!------------------------------------------------------------------------------

      SUBROUTINE REGRID_MASS_TO_4x5_GCAP( I1, J1, L1, IN, I4, J4, OUT )
!
!******************************************************************************
!  Subroutine REGRID_MASS_TO_4x5_GCAP regrids mass data from the 
!  GEOS-Chem 1x1 grid to the GEOS-Chem 4x5 GCAP grid. (bdf, bmy, 10/24/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I1  (INTEGER) : 1x1 longitude dimension of IN array
!  (2 ) J1  (INTEGER) : 1x1 latitude  dimension of IN array
!  (3 ) L1  (INTEGER) : 1x1 altitude  dimension of IN array
!  (4 ) IN  (REAL*8 ) : Array containing input data on GEOS-Chem 1x1 grid
!  (5 ) I4  (INTEGER) : 4x5 longitude dimension of OUT array
!  (6 ) J4  (INTEGER) : 4x5 latitude  dimension of OUT array
!
!  Arguments as Output:
!  ============================================================================
!  (7 ) OUT (REAL*8 ) : Array containing output data on GEOS-Chem 4x5 grid
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN)  :: I1, J1, L1, I4, J4
      REAL*8,  INTENT(IN)  :: IN(I1,J1,L1)
      REAL*8,  INTENT(OUT) :: OUT(I4,J4,L1)

      ! Local variables
      INTEGER              :: I, J, L, W, E, S, N
      
      !==================================================================
      ! REGRID_MASS_TO_4x5_GCAP begins here!
      !==================================================================

      ! Loop over levels
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, W, E, S, N )
      DO L = 1, L1

         !-----------------------
         ! S and N Poles
         !-----------------------
         DO I = 1, I4

            ! 1x1 lon index at W edge of 4x5 box
            W           = MOD( 5 * ( I - 1 ) - 1 + I1, I1 )
 
            ! 1x1 lon index at E edge of 4x5 box
            E           = 5 * ( I - 1 ) + 3

            ! Total mass of 1x1 boxes w/in to 4x5 S pole box
            OUT(I,1,L)  = SUM( IN( W  :W+1, 1,     L ) ) + 
     &                    SUM( IN( W  :W+1, 2,     L ) ) + 
     &                    SUM( IN( W  :W+1, 3,     L ) ) + 
     &                    SUM( IN( W  :W+1, 4,     L ) ) + 
     &            0.5d0 * SUM( IN( W  :W+1, 5,     L ) ) + 
     &                    SUM( IN( E-2:E,   1,     L ) ) +
     &                    SUM( IN( E-2:E,   2,     L ) ) + 
     &                    SUM( IN( E-2:E,   3,     L ) ) + 
     &                    SUM( IN( E-2:E,   4,     L ) ) + 
     &            0.5d0 * SUM( IN( E-2:E,   5,     L ) )

            ! Total mass of 1x1 boxes w/in to 4x5 N pole box
            OUT(I,J4,L) = SUM( IN( W  :W+1, J1,    L ) ) +
     &                    SUM( IN( W  :W+1, J1-1,  L ) ) +
     &                    SUM( IN( W  :W+1, J1-2,  L ) ) +
     &                    SUM( IN( W  :W+1, J1-3,  L ) ) +
     &            0.5d0 * SUM( IN( W  :W+1, J1-4,  L ) ) +
     &                    SUM( IN( E-2:E,   J1,    L ) ) +
     &                    SUM( IN( E-2:E,   J1-1,  L ) ) +
     &                    SUM( IN( E-2:E,   J1-2,  L ) ) +
     &                    SUM( IN( E-2:E,   J1-3,  L ) ) +
     &            0.5d0 * SUM( IN( E-2:E,   J1-4,  L ) )

         ENDDO

         !-----------------------
         ! Non-polar latitudes
         !-----------------------
         DO J = 2, J4-1 

            ! 1x1 lat index at S edge of 4x5 box
            S          = ( 4 * ( J - 1 ) ) + 1

            ! 1x1 lat index at N edge of 4x5 box
            N          = ( J * 4 ) + 1
      
         DO I = 1, I4
          
            ! 1x1 lon index at W edge of 4x5 box
            W          = MOD( 5*( I - 1 ) - 1 + I1, I1 )

            ! 1x1 lon index at E edge of 4x5 box
            E          = 5*( I -1 ) + 3

            ! Total mass w/in the 4x5 box at (I,J,L)
            OUT(I,J,L) = 0.5d0 * SUM( IN( W  :W+1, S,   L ) ) +
     &                   0.5d0 * SUM( IN( E-2:E,   S,   L ) ) +
     &                           SUM( IN( W  :W+1, S+1, L ) ) +
     &                           SUM( IN( E-2:E,   S+1, L ) ) +
     &                           SUM( IN( W  :W+1, S+2, L ) ) +
     &                           SUM( IN( E-2:E,   S+2, L ) ) +
     &                           SUM( IN( W  :W+1, S+3, L ) ) +
     &                           SUM( IN( E-2:E,   S+3, L ) ) +
     &                   0.5d0 * SUM( IN( W  :W+1, N,   L ) ) +
     &                   0.5d0 * SUM( IN( E-2:E,   N,   L ) )

         ENDDO 
         ENDDO  

      ENDDO 
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE REGRID_MASS_TO_4x5_GCAP

!------------------------------------------------------------------------------

      SUBROUTINE REGRID_CONC_TO_4x5( I1, J1, L1, IN, I4, J4, OUT )
!
!******************************************************************************
!  Subroutine REGRID_CONC_TO_4x5 regrids concentration data from the 
!  GEOS-Chem 1x1 grid to the GEOS_CHEM 4x5 grid. (bdf, bmy, 10/24/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I1  (INTEGER) : 1x1 longitude dimension of IN array
!  (2 ) J1  (INTEGER) : 1x1 latitude  dimension of IN array
!  (3 ) L1  (INTEGER) : 1x1 altitude  dimension of IN array
!  (4 ) IN  (REAL*8 ) : Array containing input data on GEOS-Chem 1x1 grid
!  (5 ) I4  (INTEGER) : 4x5 longitude dimension of OUT array
!  (6 ) J4  (INTEGER) : 4x5 latitude  dimension of OUT array
!
!  Arguments as Output:
!  ============================================================================
!  (7 ) OUT (REAL*8 ) : Array containing output data on GEOS-Chem 4x5 grid
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_AREA_M2

      ! Arguments
      INTEGER, INTENT(IN)  :: I1, J1, L1, I4, J4
      REAL*8,  INTENT(IN)  :: IN(I1,J1,L1)
      REAL*8,  INTENT(OUT) :: OUT(I4,J4,L1)

      ! Local variables
      INTEGER              :: I, J, L, W, E, S, N
      REAL*8               :: M_TOT
      
      !==================================================================
      ! REGRID_CONC_TO_4x5 begins here!
      !==================================================================

      ! Loop over levels
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, W, E, M_TOT, S, N )
      DO L = 1, L1

         !-----------------------
         ! S and N Poles
         !-----------------------
         DO I = 1, I4

            ! 1x1 lon index at W edge of 4x5 box
            W           = MOD( 5 * ( I - 1 ) - 1 + I1, I1 )
 
            ! 1x1 lon index at E edge of 4x5 box
            E           = 5 * ( I - 1 ) + 3

            ! Total mass of 1x1 boxes w/in to 4x5 S pole box
            M_TOT       = SUM( IN( W  :W+1, 1,    L ) ) * A1x1(1) + 
     &                    SUM( IN( W  :W+1, 2,    L ) ) * A1x1(2) + 
     &                    SUM( IN( E-2:E,   1,    L ) ) * A1x1(1) +
     &                    SUM( IN( E-2:E,   2,    L ) ) * A1x1(2) + 
     &              0.5d0*SUM( IN( W  :W+1, 3,    L ) ) * A1x1(3) +
     &              0.5d0*SUM( IN( E-2:E,   3,    L ) ) * A1x1(3)

            ! Output field at 4x5 S pole box
            OUT(I,1,L)  = M_TOT /
     &                    ( 5d0* ( A1x1(1) + A1x1(2) + 0.5d0*A1x1(3) ) )

            ! Total mass of 1x1 boxes w/in to 4x5 N pole box
            M_TOT       = SUM( IN( W  :W+1, J1,   L ) ) * A1x1(J1  ) + 
     &                    SUM( IN( W  :W+1, J1-1, L ) ) * A1x1(J1-1) + 
     &                    SUM( IN( E-2:E,   J1,   L ) ) * A1x1(J1  ) + 
     &                    SUM( IN( E-2:E,   J1-1, L ) ) * A1x1(J1-1) + 
     &              0.5d0*SUM( IN( W  :W+1, J1-2, L ) ) * A1x1(J1-2) +
     &              0.5d0*SUM( IN( E-2:E,   J1-2, L ) ) * A1x1(J1-2)

            ! Output field at 4x5 N pole box
            OUT(I,J4,L) = M_TOT/
     &                  ( 5d0* ( A1x1(J1) + A1x1(J1-1)+ 0.5*A1x1(J1-2)))

         ENDDO

         !-----------------------
         ! Non-polar latitudes
         !-----------------------
         DO J = 2, J4-1 

            ! 1x1 lat index at S edge of 4x5 box
            S     = ( 4 * ( J - 1 ) ) - 1

            ! 1x1 lat index at N edge of 4x5 box
            N     = ( J * 4 ) - 1
      
         DO I = 1, I4
          
            ! 1x1 lon index at W edge of 4x5 box
            W     = MOD( 5*( I - 1 ) - 1 + I1, I1 )

            ! 1x1 lon index at E edge of 4x5 box
            E     = 5*( I -1 ) + 3

            ! Total mass w/in the 4x5 box at (I,J,L)
            M_TOT = 0.5d0*SUM( IN( W  :W+1, S,   L ) ) * A1x1(S  ) + 
     &              0.5d0*SUM( IN( E-2:E,   S,   L ) ) * A1x1(S  ) + 
     &                    SUM( IN( W  :W+1, S+1, L ) ) * A1x1(S+1) + 
     &                    SUM( IN( E-2:E,   S+1, L ) ) * A1x1(S+1) + 
     &                    SUM( IN( W  :W+1, S+2, L ) ) * A1x1(S+2) + 
     &                    SUM( IN( E-2:E,   S+2, L ) ) * A1x1(S+2) + 
     &                    SUM( IN( W  :W+1, S+3, L ) ) * A1x1(S+3) + 
     &                    SUM( IN( E-2:E,   S+3, L ) ) * A1x1(S+3) + 
     &              0.5d0*SUM( IN( W  :W+1, N,   L ) ) * A1x1(N  ) +
     &              0.5d0*SUM( IN( E-2:E,   N,   L ) ) * A1x1(N  )
            
            ! 4x5 output field at (I,J,L)
            OUT(I,J,L) = M_TOT / GET_AREA_M2( J )

         ENDDO 
         ENDDO  
      ENDDO 
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE REGRID_CONC_TO_4x5

!------------------------------------------------------------------------------

      SUBROUTINE REGRID_MASS_TO_4x5( I1, J1, L1, IN, I4, J4, OUT )
!
!******************************************************************************
!  Subroutine REGRID_MASS_TO_4x5 regrids mass data from the GEOS-Chem 1x1 
!  grid to the GEOS_CHEM 4x5 grid. (bdf, bmy, 10/24/05, 2/20/07)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I1  (INTEGER) : 1x1 longitude dimension of IN array
!  (2 ) J1  (INTEGER) : 1x1 latitude  dimension of IN array
!  (3 ) L1  (INTEGER) : 1x1 altitude  dimension of IN array
!  (4 ) IN  (REAL*8 ) : Array containing input data on GEOS-Chem 1x1 grid
!  (5 ) I4  (INTEGER) : 4x5 longitude dimension of OUT array
!  (6 ) J4  (INTEGER) : 4x5 latitude  dimension of OUT array
!
!  Arguments as Output:
!  ============================================================================
!  (7 ) OUT (REAL*8 ) : Array containing output data on GEOS-Chem 4x5 grid
!
!  NOTES:
!  (1 ) Bug fix: the lat index should be N, not S in the last 2 lines of the
!        non-polar latitude regridding.  (tw, bmy, 2/20/07)
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN)  :: I1, J1, L1, I4, J4
      REAL*8,  INTENT(IN)  :: IN(I1,J1,L1)
      REAL*8,  INTENT(OUT) :: OUT(I4,J4,L1)

      ! Local variables
      INTEGER              :: I, J, L, W, E, S, N
      REAL*8               :: M_TOT

      !=================================================================
      ! REGRID_MASS_TO_4x5 begins here!
      !=================================================================

      ! Loop over levels
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, W, E, S, N )
      DO L = 1, L1

         !-----------------------
         ! S and N Poles
         !-----------------------
         DO I = 1, I4

            ! 1x1 lon index at W edge of 4x5 box
            W           = MOD( 5 * ( I - 1 ) - 1 + I1, I1 )
 
            ! 1x1 lon index at E edge of 4x5 box
            E           = 5 * ( I - 1 ) + 3

            ! Output field at 4x5 S Pole box
            OUT(I,1,L)  =      SUM( IN( W  :W+1, 1:2,     L ) ) + 
     &                         SUM( IN( E-2:E,   1:2,     L ) ) + 
     &                   0.5d0*SUM( IN( W  :W+1, 3,       L ) ) +
     &                   0.5d0*SUM( IN( E-2:E,   3,       L ) )

            ! Output field at 4x5 N pole box
            OUT(I,J4,L) =      SUM( IN( W  :W+1, J1-1:J1, L ) ) + 
     &                         SUM( IN( E-2:E,   J1-1:J1, L ) ) + 
     &                   0.5d0*SUM( IN( W  :W+1, J1-2,    L ) ) +
     &                   0.5d0*SUM( IN( E-2:E,   J1-2,    L ) )
         ENDDO

         !-----------------------
         ! Non-polar latitudes
         !-----------------------
         DO J = 2, J4-1 

            ! 1x1 lat index at S edge of 4x5 box
            S          = ( 4 * ( J - 1 ) ) - 1

            ! 1x1 lat index at Northern edge of 4x5 box
            N          = ( J * 4 ) - 1
          
         DO I = 1, I4

            ! 1x1 lon index at W edge of the 4x5 box
            W          = MOD( 5 * ( I - 1 ) - 1 + I1, I1 )
         
            ! 1x1 lon index at E edge of 4x5 box
            E          = 5 * ( I - 1 ) + 3

            ! Output value for 4x5 grid box (I,J,L)
            OUT(I,J,L) = 0.5d0*SUM( IN( W  :W+1, S,       L ) ) + 
     &                   0.5d0*SUM( IN( E-2:E,   S,       L ) ) + 
     &                         SUM( IN( W  :W+1, S+1:N-1, L ) ) + 
     &                         SUM( IN( E-2:E,   S+1:N-1, L ) ) + 
!-----------------------------------------------------------------------------
! Prior to 2/20/07:
! Bug fix: the lat-index should be N, not S in the last 2 lines.
! (mb, bmy, 2/20/07)
!     &                   0.5d0*SUM( IN( W  :W+1, S,       L ) ) +
!     &                   0.5d0*SUM( IN( E-2:E,   S,       L ) )
!-----------------------------------------------------------------------------
     &                   0.5d0*SUM( IN( W  :W+1, N,       L ) ) +
     &                   0.5d0*SUM( IN( E-2:E,   N,       L ) )

         ENDDO
         ENDDO

      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE REGRID_MASS_TO_4x5

!------------------------------------------------------------------------------

      SUBROUTINE REGRID_CONC_TO_2x25( I1, J1, L1, IN, I2, J2, OUT )
!
!******************************************************************************
!  Subroutine REGRID_CONC_TO_2x25 regrids concentration data from the 
!  GEOS-Chem 1x1 grid to the GEOS_CHEM 2x25 grid. (bdf, bmy, 10/24/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I1  (INTEGER) : 1x1  longitude dimension of IN array
!  (2 ) J1  (INTEGER) : 1x1  latitude  dimension of IN array
!  (3 ) L1  (INTEGER) : 1x1  altitude  dimension of IN array
!  (4 ) IN  (REAL*8 ) : Array containing input data on GEOS-Chem 1x1 grid
!  (5 ) I2  (INTEGER) : 2x25 longitude dimension of OUT array
!  (6 ) J2  (INTEGER) : 2x25 latitude  dimension of OUT array
!
!  Arguments as Output:
!  ============================================================================
!  (7 ) OUT (REAL*8 ) : Array containing output data on GEOS-Chem 2x25 grid
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_AREA_M2

      ! Arguments
      INTEGER, INTENT(IN)  :: I1, J1, L1, I2, J2
      REAL*8,  INTENT(IN)  :: IN(I1,J1,L1)
      REAL*8,  INTENT(OUT) :: OUT(I2,J2,L1)

      ! Local variables
      INTEGER              :: I, J, L, W, E, S, N
      REAL*8               :: M_TOT

      !=================================================================
      ! REGRID_CONC_TO_2x25 begins here!
      !=================================================================
      
      ! Loop over levels
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, W, E, M_TOT, S, N )
      DO L = 1, L1

         !-----------------------
         ! S and N Poles
         !-----------------------
         DO I = 1, I2

            ! 1x1 lon index at W edge of 2x25 box
            W = FLOOR( 2.5d0 * ( I - 1 ) )
            IF ( W == 0 ) W = 360

            ! 1x1 lon index at E edge of 2x25 box
            E = FLOOR( 2.5d0 * I )

            ! Test for 3 or 4 contributing 1x1 longitude boxes
            IF ( MOD( I, 2 ) == 1 ) THEN

               !---------------------------------------
               ! 3 contributing 1x1 lon boxes at poles
               !---------------------------------------

               ! Total mass (w/ 3 contributing 1x1 boxes) at S Pole
               M_TOT       = 0.75d0  * IN( W,   1, L ) * A1x1(1) + 
     &                       0.375d0 * IN( W,   2, L ) * A1x1(2) +
     &                                 IN( E-1, 1, L ) * A1x1(1) +
     &                       0.5d0   * IN( E-1, 2, L ) * A1x1(2) +
     &                       0.75d0  * IN( E,   1, L ) * A1x1(1) +
     &                       0.375d0 * IN( E,   2, L ) * A1x1(2)

               ! Output field at 2 x 2.5 S pole box
               OUT(I,1,L)  = M_TOT / 
     &                     ( 2.5d0 * ( A1x1(1) + 0.5d0*A1x1(2) ) )

               ! Total mass (w/ 3 contributing 1x1 lon boxes) at N pole
               M_TOT       = 0.75d0  * IN( W,   J1,   L ) * A1x1(J1  ) + 
     &                       0.375d0 * IN( W,   J1-1, L ) * A1x1(J1-1) +
     &                                 IN( E-1, J1,   L ) * A1x1(J1  ) +
     &                       0.5d0   * IN( E-1, J1-1, L ) * A1x1(J1-1) +
     &                       0.75d0  * IN( E,   J1,   L ) * A1x1(J1  ) +
     &                       0.375d0 * IN( E,   J1-1, L ) * A1x1(J1-1)

               ! Output field at 2 x 2.5 N pole box
               OUT(I,J2,L) = M_TOT/
     &                     ( 2.5d0 * ( A1x1(J1) + 0.5d0*A1x1(J1-1) ) )

            ELSE

               !---------------------------------------
               ! 4 contributing 1x1 lon boxes at poles
               !---------------------------------------

               ! Total mass (w/ 4 contributing 1x1 lon boxes) at S pole
               M_TOT       = 
     &           0.25d0  *      IN( W,       1, L )   * A1x1(1) + 
     &           0.125d0 *      IN( W,       2, L )   * A1X1(2) +
     &                     SUM( IN( W+1:E-1, 1, L ) ) * A1x1(1) +
     &           0.5d0   * SUM( IN( W+1:E-1, 2, L ) ) * A1x1(2) +
     &           0.25d0  *      IN( E,       1, L )   * A1x1(1) +
     &           0.125d0 *      IN( E,       2, L )   * A1x1(2)

               ! Output field at 2 x 2.5 S pole box
               OUT(I,1,L)  = M_TOT/
     &                     ( 2.5d0* ( A1x1(1) + 0.5d0*A1x1(2) ) )

               ! Total mass (w/ 4 contributing 1x1 lon boxes) at N pole
               M_TOT =       
     &           0.25d0  *      IN( W,       J1,   L )   * A1x1(J1  ) + 
     &           0.125d0 *      IN( W,       J1-1, L )   * A1x1(J1-1) +
     &                     SUM( IN( W+1:E-1, J1,   L ) ) * A1x1(J1  ) +
     &           0.5d0   * SUM( IN( W+1:E-1, J1-1, L ) ) * A1x1(J1-1) +
     &           0.25d0  *      IN( E,       J1,   L )   * A1x1(J1  ) +
     &           0.125d0 *      IN( E,       J1-1, L )   * A1x1(J1-1)

               ! Output field at 2 x 2.5 N pole box
               OUT(I,J2,L) = M_TOT/
     &                     ( 2.5d0* ( A1x1(J1) + 0.5d0*A1x1(J1-1) ) )

            ENDIF
         ENDDO

         !-----------------------
         ! Non-polar latitudes
         !-----------------------
         DO J = 2, J2-1 
            
            ! 1x1 lat index at S edge of 2 x 2.5 box
            S = 2 * ( J - 1 )

            ! 1x1 lat index at N edge of 2 x 2.5 box
            N = 2 * J 

         DO I = 1, I2

            ! 1x1 lon index at W edge of 2 x 2.5 box
            W = FLOOR( 2.5d0 * ( I - 1 ) )
            IF ( W == 0 ) W = 360

            ! 1x1 lon index at E edge of 2 x 2.5 box
            E = FLOOR( 2.5d0 * I )

            ! Test for 3 or 4 contributing 1x1 lon boxes
            IF ( MOD( I, 2 ) == 1 )  THEN

               !------------------------------
               ! 3 contributing 1x1 lon boxes
               !------------------------------

               ! Total mass (w/ 3 contributing 1x1 lon boxes) in 2 x 2.5 box
               M_TOT       = 0.375d0 * IN(W,  S,  L) * A1x1(S  ) + 
     &                       0.75d0  * IN(W,  S+1,L) * A1x1(S+1) +
     &                       0.375d0 * IN(W,  N,  L) * A1x1(N  ) +
     &                       0.5d0   * IN(E-1,S,  L) * A1x1(S  ) +
     &                                 IN(E-1,S+1,L) * A1x1(S+1) +
     &                       0.5d0   * IN(E-1,N,  L) * A1x1(N  ) +
     &                       0.375d0 * IN(E,  S,  L) * A1x1(S  ) +
     &                       0.75d0  * IN(E,  S+1,L) * A1x1(S+1) +
     &                       0.375d0 * IN(E,  N,  L) * A1x1(N  )

               ! 2 x 2.5 output field at (I,J,L)
               OUT(I,J,L)  = M_TOT / GET_AREA_M2( J )

            ELSE

               !------------------------------
               ! 4 contributing 1x1 lon boxes
               !------------------------------

               ! Total mass (w/ 4 contributing 1x1 lon boxes) in 2 x 2.5 box
               M_TOT       = 
     &           0.125d0 *      IN( W,       S,  L )   * A1x1(S  ) + 
     &           0.25d0  *      IN( W,       S+1,L )   * A1x1(S+1) +
     &           0.125d0 *      IN( W,       N,  L )   * A1x1(N  ) +
     &           0.5d0   * SUM( IN( W+1:E-1, S,  L ) ) * A1x1(S  ) +
     &                     SUM( IN( W+1:E-1, S+1,L ) ) * A1x1(S+1) +
     &           0.5d0   * SUM( IN( W+1:E-1, N,  L ) ) * A1x1(N  ) +
     &           0.125d0 *      IN( E,       S,  L )   * A1x1(S  ) +
     &           0.25d0  *      IN( E,       S+1,L )   * A1x1(S+1) +
     &           0.125d0 *      IN( E,       N,  L )   * A1X1(N  )

               ! 2 x 2.5 output field at (I,J,L)
               OUT(I,J,L) = M_TOT / GET_AREA_M2( J )
            ENDIF
         ENDDO
         ENDDO

      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE REGRID_CONC_TO_2x25

!------------------------------------------------------------------------------

      SUBROUTINE REGRID_MASS_TO_2x25( I1, J1, L1, IN, I2, J2, OUT )
!
!******************************************************************************
!  Subroutine REGRID_CONC_TO_2x25 regrids mass data from the GEOS-Chem 1x1 
!  grid to the GEOS_CHEM 2x25 grid. (bdf, bmy, 10/24/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I1  (INTEGER) : 1x1  longitude dimension of IN array
!  (2 ) J1  (INTEGER) : 1x1  latitude  dimension of IN array
!  (3 ) L1  (INTEGER) : 1x1  altitude  dimension of IN array
!  (4 ) IN  (REAL*8 ) : Array containing input data on GEOS-Chem 1x1 grid
!  (5 ) I2  (INTEGER) : 2x25 longitude dimension of OUT array
!  (6 ) J2  (INTEGER) : 2x25 latitude  dimension of OUT array
!
!  Arguments as Output:
!  ============================================================================
!  (7 ) OUT (REAL*8 ) : Array containing output data on GEOS-Chem 2x25 grid
!
!  NOTES:
!  (1 ) Fixed typo: J should be J1 in "4 contrib boxes at poles" section.
!        (bmy, 4/18/06)
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN)  :: I1, J1, L1, I2, J2
      REAL*8,  INTENT(IN)  :: IN(I1,J1,L1)
      REAL*8,  INTENT(OUT) :: OUT(I2,J2,L1)

      ! Local variables
      INTEGER              :: I, J, L, W, E, S, N
      REAL*8               :: M_TOT

      !=================================================================
      ! REGRID_MASS_TO_2x25 begins here!
      !=================================================================

      ! Loop over levels
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, W, E, S, N )
      DO L = 1, L1

         !-----------------------
         ! S and N Poles
         !-----------------------

         DO I = 1, I2

            ! 1x1 lon index at W edge of 2x25 box
            W = FLOOR( 2.5d0 * ( I - 1 ) )
            IF ( W == 0 ) W = 360

            ! 1x1 lon index at E edge of 2x25 box
            E = FLOOR( 2.5d0 * I )

            ! Test for 3 or 4 contributing 1x1 longitude boxes
            IF ( MOD( I, 2 ) == 1 ) THEN

               !---------------------------------------
               ! 3 contributing 1x1 lon boxes at poles
               !---------------------------------------

               ! Output field at 2 x 2.5 S Pole box
               OUT(I,1,L)  = 0.75d0  * IN( W,   1,    L ) + 
     &                       0.375d0 * IN( W,   2,    L ) +
     &                                 IN( E-1, 1,    L ) +
     &                       0.5d0   * IN( E-1, 2,    L ) +
     &                       0.75d0  * IN( E,   1,    L ) +
     &                       0.375d0 * IN( E,   2,    L )

               ! Output field at 2 x 2.5 N pole box
               OUT(I,J2,L) = 0.75d0  * IN( W,   J1,   L ) + 
     &                       0.375d0 * IN( W,   J1-1, L ) +
     &                                 IN( E-1, J1,   L ) +
     &                       0.5d0   * IN( E-1, J1-1, L ) +
     &                       0.75d0  * IN( E,   J1,   L ) +
     &                       0.375d0 * IN( E,   J1-1, L )
            ELSE

               !---------------------------------------
               ! 4 contributing 1x1 lon boxes at poles
               !---------------------------------------

               ! Output field at 2 x 2.5 S Pole box
               OUT(I,1,L)  = 0.25d0  *      IN( W,       1,    L )   + 
     &                       0.125d0 *      IN( W,       2,    L )   +
     &                                 SUM( IN( W+1:E-1, 1,    L ) ) +
     &                       0.5d0   * SUM( IN( W+1:E-1, 2,    L ) ) +
     &                       0.25d0  *      IN( E,       1,    L )   +
     &                       0.125d0 *      IN( E,       2,    L )

               ! Output field at 2 x 2.5 N pole box
               OUT(I,J2,L) = 0.25d0  *      IN( W,       J1,   L )   + 
     &                       0.125d0 *      IN( W,       J1-1, L )   +
     &                                 SUM( IN( W+1:E-1, J1,   L ) ) +
     &                       0.5d0   * SUM( IN( W+1:E-1, J1-1, L ) ) +
     &                       0.25d0  *      IN( E,       J1,   L )   +
     &                       0.125d0 *      IN( E,       J1-1, L )
            ENDIF

         ENDDO
         
         !-----------------------
         ! Non-polar latitudes
         !-----------------------
         DO J = 2, J2-1    

            ! 1x1 lat index at S edge of 2x25 box
            S = 2 * ( J - 1 )

            ! 1x1 lat index at N edge of 2x25 box
            N = 2 * J

         DO I = 1, J1

            ! 1x1 lon index at W edge of 2x25 box
            W = FLOOR( 2.5d0 * ( I - 1 ) )
            IF ( W == 0 ) W = 360

            ! 1x1 lon index at E edge of 2x25 box
            E = FLOOR( 2.5d0 * I )

            ! Test for 3 or 4 contributing 1x1 lon boxes
            IF ( MOD( I, 2 ) == 1 ) THEN

               !------------------------------
               ! 3 contributing 1x1 lon boxes 
               !------------------------------

               ! Output value at 2x25 box (I,J,L)
               OUT(I,J,L) = 0.375d0 * IN( W,   S,   L ) + 
     &                      0.75d0  * IN( W,   S+1, L ) +
     &                      0.375d0 * IN( W,   N,   L ) +
     &                      0.5d0   * IN( E-1, S,   L ) +
     &                                IN( E-1, S+1, L ) +
     &                      0.5d0   * IN( E-1, N,   L ) +
     &                      0.375d0 * IN( E,   S,   L ) +
     &                      0.75d0  * IN( E,   S+1, L ) +
     &                      0.375d0 * IN( E,   N,   L )
            ELSE

               !------------------------------
               ! 4 contributing 1x1 lon boxes 
               !------------------------------

               ! Output value at 2 x 2.5 box (I,J,L)
               OUT(I,J,L) = 0.125d0  *      IN( W,       S,   L )   +
     &                      0.25d0   *      IN( W,       S+1, L )   +
     &                      0.125d0  *      IN( W,       N,   L )   +
     &                      0.5d0    * SUM( IN( W+1:E-1, S,   L ) ) +
     &                                 SUM( IN( W+1:E-1, S+1, L ) ) +
     &                      0.5d0    * SUM( IN( W+1:E-1, N,   L ) ) +
     &                      0.125d0  *      IN( E,       S,   L )   +
     &                      0.25d0   *      IN( E,       S+1, L )   +
     &                      0.125d0  *      IN( E,       N,   L )
            ENDIF
         ENDDO
         ENDDO

      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE REGRID_MASS_TO_2x25

!------------------------------------------------------------------------------

      SUBROUTINE REGRID_CONC_TO_1x125( I1, J1, L1, IN, I2, J2, OUT )
!
!******************************************************************************
!  Subroutine REGRID_CONC_TO_1x125 regrids conc data from the GEOS-Chem 
!  1x1 grid to the GEOS_CHEM 1x125 grid. (bdf, bmy, 8/2/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I1  (INTEGER) : 1x1  longitude dimension of IN array
!  (2 ) J1  (INTEGER) : 1x1  latitude  dimension of IN array
!  (3 ) L1  (INTEGER) : 1x1  altitude  dimension of IN array
!  (4 ) IN  (REAL*8 ) : Array containing input data on GEOS-Chem 1x1 grid
!  (5 ) I2  (INTEGER) : 1x125 longitude dimension of OUT array
!  (6 ) J2  (INTEGER) : 1x125 latitude  dimension of OUT array
!
!  Arguments as Output:
!  ============================================================================
!  (7 ) OUT (REAL*8 ) : Array containing output data on GEOS-Chem 1x125 grid
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD,    ONLY : GET_AREA_M2

      ! Arguments
      INTEGER, INTENT(IN)  :: I1, J1, L1, I2, J2
      REAL*8,  INTENT(IN)  :: IN(I1,J1,L1)
      REAL*8,  INTENT(OUT) :: OUT(I2,J2,L1)

      ! Local variables
      INTEGER              :: I, J, L, W, E, C, OFFSET, PLACE
      REAL*8               :: M_TOT

      !=================================================================
      ! REGRID_CONC_TO_1x125 begins here!
      !=================================================================

      ! Loop over levels
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, OFFSET, W, C, E, PLACE, M_TOT )
      DO L = 1, L1

         ! Poles can be done at same time, no latitude differences
         DO J = 1, J2
         DO I = 1, I2

            ! 1x1 offset.  there is 1 extra 1x1 box for every 4 1x125 boxes
            OFFSET = FLOOR( ( I - 1 ) / 4d0 )

            ! West, center, east longitude indices
            W      = I - 1 + OFFSET
            C      = I     + OFFSET
            E      = I + 1 + OFFSET

            ! Special handling 
            IF ( W == 0 ) W = 360

            ! There are 4 possible cases for overlap
            PLACE  = MOD( ( I - 1 ), 4 )

            !-----------------------------------------------------------
            ! Pick the right case for the overlap
            !
            ! Because there is no difference in the latitude coordinates
            ! between 1x1 and 1x1.25 grids, the concentration ratio is: 
            !
            !    [ area(1x1) / (1.24 or 1.26)*area(1x1) ]
            ! 
            ! The 1.24 or 1.26 depends on how much overlap the 1x125 
            ! grid has with the 1x1 grid.  
            !-----------------------------------------------------------

            SELECT CASE ( PLACE )

               !----------------------------------------------
               ! CASE 0: 1x1 and 1x125 are centered the same
               !----------------------------------------------
               CASE( 0 )
                  M_TOT = 0.12d0 * IN(W,J,L) + 
     &                             IN(C,J,L) + 
     &                    0.12d0 * IN(E,J,L)
               
                  ! Overlap factor for CASE 0 is 1.24
                  OUT(I,J,L) = M_TOT / 1.24d0

               !----------------------------------------------
               ! CASE 1: one to the right of a centered box
               !----------------------------------------------
               CASE ( 1 )
                  M_TOT = 0.88d0 * IN(C,J,L) +
     &                    0.38d0 * IN(E,J,L)  

                  ! Overlap factor for CASE 1 is 1.26
                  OUT(I,J,L) = M_TOT / 1.26d0

               !----------------------------------------------
               ! CASE 2: 1x1 and 1x125 are edged the same
               !----------------------------------------------
               CASE ( 2 )
                  M_TOT = 0.62d0 * IN(C,J,L) +
     &                    0.62d0 * IN(E,J,L)

                  ! Overlap factor for CASE 2 is 1.24
                  OUT(I,J,L) = M_TOT / 1.24d0

               !----------------------------------------------
               ! CASE 3: one to the left of a centered box
               !----------------------------------------------
               CASE ( 3 )
                  M_TOT = 0.38d0 * IN(C,J,L) +
     &                    0.88d0 * IN(E,J,L)

                  ! Overlap factor for CASE 3 is 1.26
                  OUT(I,J,L) = M_TOT / 1.26d0

            END SELECT
         ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE REGRID_CONC_TO_1x125

!------------------------------------------------------------------------------

      SUBROUTINE REGRID_MASS_TO_1x125( I1, J1, L1, IN, I2, J2, OUT )
!
!******************************************************************************
!  Subroutine REGRID_MASS_TO_1x125 regrids mass data from the 
!  GEOS-Chem 1x1 grid to the GEOS-Chem 1x125 grid. (bdf, bmy, 8/2/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I1  (INTEGER) : 1x1  longitude dimension of IN array
!  (2 ) J1  (INTEGER) : 1x1  latitude  dimension of IN array
!  (3 ) L1  (INTEGER) : 1x1  altitude  dimension of IN array
!  (4 ) IN  (REAL*8 ) : Array containing input data on GEOS-Chem 1x1 grid
!  (5 ) I2  (INTEGER) : 1x125 longitude dimension of OUT array
!  (6 ) J2  (INTEGER) : 1x125 latitude  dimension of OUT array
!
!  Arguments as Output:
!  ============================================================================
!  (7 ) OUT (REAL*8 ) : Array containing output data on GEOS-Chem 1x125 grid
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_AREA_M2

      ! Arguments
      INTEGER, INTENT(IN)  :: I1, J1, L1, I2, J2
      REAL*8,  INTENT(IN)  :: IN(I1,J1,L1)
      REAL*8,  INTENT(OUT) :: OUT(I2,J2,L1)

      ! Local variables
      INTEGER              :: I, J, L, W, E, C, OFFSET, PLACE
      REAL*8               :: M_TOT

      !=================================================================
      ! REGRID_MASS_TO_1x125 begins here!
      !=================================================================

      ! Loop over levels
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, OFFSET, W, C, E, PLACE )
      DO L = 1, L1

         ! Poles can be done at same time, no latitude differences
         DO J = 1, J2
         DO I = 1, I2

            ! 1x1 offset.  there is 1 extra 1x1 box for every 4 1x125 boxes
            OFFSET = FLOOR( ( I - 1 ) / 4d0 )

            ! West, center, east longitude indices
            W      = I - 1 + OFFSET
            C      = I     + OFFSET
            E      = I + 1 + OFFSET

            ! Special handling
            IF ( W == 0 ) W = 360

            ! There are 4 possible casses for overlap
            PLACE = MOD( ( I -1 ), 4 )

            SELECT CASE ( PLACE )

               !----------------------------------------------
               ! CASE 0: 1x1 and 1x125 are centered the same
               !----------------------------------------------
               CASE( 0 )
                  OUT(I,J,L) = 0.12d0 * IN(W,J,L) + 
     &                                  IN(C,J,L) + 
     &                         0.12d0 * IN(E,J,L)
               
               !----------------------------------------------
               ! CASE 1: one to the right of a centered box
               !----------------------------------------------
               CASE ( 1 )
                  OUT(I,J,L) = 0.88d0 * IN(C,J,L) +
     &                         0.38d0 * IN(E,J,L)

               !----------------------------------------------
               ! CASE 2: 1x1 and 1x125 are edged the same
               !----------------------------------------------
               CASE ( 2 )
                  OUT(I,J,L) = 0.62d0 * IN(C,J,L) +
     &                         0.62d0 * IN(E,J,L)

               !----------------------------------------------
               ! CASE 3: one to the left of a centered box
               !----------------------------------------------
               CASE ( 3 )
                  OUT(I,J,L) = 0.38d0 * IN(C,J,L) +
     &                         0.88d0 * IN(E,J,L)

            END SELECT
         ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE REGRID_MASS_TO_1x125

!------------------------------------------------------------------------------

      SUBROUTINE INIT_REGRID_1x1
!
!******************************************************************************
!  Subroutine INIT_REGRID_1x1 initializes module arrays 
!  (bdf, bmy, 10/24/05, 4/18/06)
!
!  NOTES:
!  (1 ) Now exit if we have already initialized (bmy, 4/18/06)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_GCTM"  ! Physical constants

      ! Local variables
      LOGICAL, SAVE      :: IS_INIT = .FALSE.
      INTEGER            :: AS, J
      REAL*8             :: S,  N, RLAT, YEDGE(J1x1+1)

      !=================================================================
      ! INIT_REGRID_1x1 begins here!
      !=================================================================

      ! Return if we have already initialized
      IF ( IS_INIT ) RETURN

      !---------------------------------------
      ! Surface area on GEOS-Chem 1x1 grid
      ! Uses same algorithm from "grid_mod.f"
      !---------------------------------------

      ! Allocate array
      ALLOCATE( A1x1( J1x1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'A1x1' )

      ! Initialize
      YEDGE(:) = 0d0

      ! 1x1 latitude edges
      DO J = 2, J1x1
         YEDGE(J) = -90.5d0 + ( J - 1 )
      ENDDO

      ! Special cases at poles
      YEDGE(1)      = -90.0d0
      YEDGE(2)      = -89.5d0
      YEDGE(J1x1+1) =  90.0d0

      ! Compute 1x1 surface area
      DO J = 1, J1x1

         ! Lat at S and N edges of 1x1 box [radians]
         S       = PI_180 * YEDGE(J  )
         N       = PI_180 * YEDGE(J+1) 

         ! S to N extent of grid box [unitless]
         RLAT    = SIN( N ) - SIN( S )

         ! 1x1 surface area [m2] (see "grid_mod.f" for algorithm)
         A1x1(J) = 2d0 * PI * Re * Re / DBLE( I1x1 ) * RLAT
      ENDDO

      !---------------------------------------
      ! Surface area on GENERIC 1x1 grid
      ! Uses same algorithm from "grid_mod.f"
      !---------------------------------------

      ! Initialize
      YEDGE(:) = 0d0

      ! Allocate array
      ALLOCATE( A_GEN_1x1( J1x1-1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'A_GEN_1x1' )

      ! 1x1 latitude edges
      DO J = 1, J1x1
         YEDGE(J) = -90d0 + ( J - 1 )
      ENDDO

      ! Compute 1x1 surface area
      DO J = 1, J1x1-1

         ! Lat at S and N edges of 1x1 box [radians]
         S            = PI_180 * YEDGE(J  )
         N            = PI_180 * YEDGE(J+1) 

         ! S to N extent of grid box [unitless]
         RLAT         = SIN( N ) - SIN( S )

         ! 1x1 surface area [m2] (see "grid_mod.f" for algorithm)
         A_GEN_1x1(J) = 2d0 * PI * Re * Re / DBLE( I1x1 ) * RLAT
      ENDDO

      ! We have now initialized
      IS_INIT = .TRUE.

      ! Return to calling program
      END SUBROUTINE INIT_REGRID_1x1

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_REGRID_1x1
!
!******************************************************************************
!  Subroutine CLEANUP_REGRID_1x1 deallocates all module arrays.
!  (bdf, bmy, 10/24/05)
!
!  NOTES:
!******************************************************************************
!      
      !=================================================================
      ! CLEANUP_REGRID_1x1 begins here!
      !=================================================================
      IF ( ALLOCATED( A1x1      ) ) DEALLOCATE( A1x1 )
      IF ( ALLOCATED( A_GEN_1x1 ) ) DEALLOCATE( A_GEN_1x1 )

      ! Return to calling program 
      END SUBROUTINE CLEANUP_REGRID_1x1

!------------------------------------------------------------------------------

      ! End of module
      END MODULE REGRID_1x1_MOD
