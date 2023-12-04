#include "GCHEM_OPTIONS.h"

CBOP
C !ROUTINE: GCHEM_FIELDS_LOAD

C !INTERFACE: ==========================================================
      SUBROUTINE GCHEM_FIELDS_LOAD (
     I           myTime, myIter, myThid )

C !DESCRIPTION:
C  calls routines which read in fields needed for any tracer experiment
C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_FIELDS.h"
#include "GCHEM.h"
#include "GCHEM_SIZE.h"
#include "GCHEM_FIELDS.h"

C !INPUT PARAMETERS: ===================================================
C  myTime               :: current time
C  myIter               :: current iteration
C  myThid               :: thread number
      _RL myTime
      INTEGER myIter
      INTEGER myThid

C !LOCAL VARIABLES: ====================================================
C =================
C Here starts the modification by Alan
C  gchemSource of all tracers is read from the same file
C  filename1
C 24 OCT 2023
C =================
C     i,j,k,bi,bj,jTr  :: loop indices
      INTEGER i,j,k,bi,bj,jTr
C  tmpSourceFile        :: temporary source file name to read in
      CHARACTER*(MAX_LEN_FNAM) tmpSourceFile
CEOP

C     gchemSource is initialized here
      DO jTr = 1, GCHEM_tendTr_num
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO k = 1, Nr
          DO j = 1-OLy, sNy+OLy
           DO i = 1-OLx, sNx+OLx
            gchemSource(i,j,k,bi,bj,jTr) = 0. _d 0
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
    
C     Read the source file into gchemSource.
C     But only one filename is used here.
C     Maybe find some way to loop over more filenames in the future.
      IF (nIter0.EQ.PTRACERS_Iter0) THEN
       DO jTr = 1, GCHEM_tendTr_num
C     read up to three files
C     last argument of gchemSource is the tracer number
        IF (jTr .EQ. 1) THEN
            tmpSourceFile = Filename1
        ELSEIF ((jTr.EQ.2).AND.(GCHEM_tendTr_num.GE.2)) THEN
            tmpSourceFile = Filename2
        ELSEIF ((jTr.EQ.3).AND.(GCHEM_tendTr_num.GE.3)) THEN
            tmpSourceFile = Filename3
        ENDIF
        IF ( tmpSourceFile .NE. ' ' ) THEN
         CALL READ_FLD_XYZ_RL(tmpSourceFile,' ',
     &        gchemSource(1-OLx,1-OLy,1,1,1,jTr),0,myThid)
        ENDIF

       ENDDO
      ENDIF

C     Apply mask
      DO jTr = 1, GCHEM_tendTr_num
       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
         DO k=1,Nr
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            IF (maskC(i,j,k,bi,bj).EQ.0.)
     &           gchemSource(i,j,k,bi,bj,jTr)=0. _d 0
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

C =========
C Here ends the modification.
C Hope this works :P
C =========

#ifdef ALLOW_GCHEM

cccccccccccccccccccccccccc
c load external data     c
cccccccccccccccccccccccccc
#ifdef ALLOW_DIC
      IF ( useDIC ) THEN
       CALL DIC_FIELDS_LOAD( myTime, myIter, myThid )
      ENDIF
#endif

#ifdef ALLOW_BLING
      IF ( useBLING ) THEN
       CALL BLING_FIELDS_LOAD( myTime, myIter, myThid )
      ENDIF
#endif

#ifdef ALLOW_CFC
      IF ( useCFC ) THEN
       CALL CFC_FIELDS_LOAD( myTime, myIter, myThid )
      ENDIF
#endif

#ifdef ALLOW_DARWIN
      IF ( useDARWIN ) THEN
       CALl DARWIN_FIELDS_LOAD(myIter,myTime,myThid)
      ENDIF
#endif

#endif /* ALLOW_GCHEM */

      RETURN
      END

