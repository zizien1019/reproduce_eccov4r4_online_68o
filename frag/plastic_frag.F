#include "GCHEM_OPTIONS.h"

CBOP
C     !ROUTINE: PLASTIC_FRAG
C     !INTERFACE:
      SUBROUTINE PLASTIC_FRAG(
     I                          pTr_plt1,
     U                          gplt1, gplt2,
     I                          bi, bj, iMin, iMax, jMin, jMax,
     I                          myTime, myIter, myThid )

C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE PLASTIC_FRAG
C     | o Calculate the changes
C     |   from larger to smaller plastic particles
C     |   due to fragmentation.
C     |   Smaller particles have lower sinking velocity.
C     *==========================================================*

C     !USES:
C     == GLobal variables ==
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


C     !INPUT/OUTPUT PARAMETERS:
C     pTr_plt1  :: ocean larger plastic concentration
C     gplt1     :: larger plastic tendency
C     gplt2     :: smaller plastic tendency
C     i,j       :: loop index
C     k         :: vertical level
C     bi, bj     :: current tile indices
C     iMin,iMax  :: computation domain, 1rst index bounds
C     jMin,jMax  :: computation domain, 2nd  index bounds
C     myTime     :: current time in simulation
C     myIter     :: current iteration number
C     myThid     :: my Thread Id number
      _RL  pTr_plt1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  gplt1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  gplt2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      INTEGER i,j,k
      INTEGER bi, bj
      INTEGER iMin, iMax, jMin, jMax
      _RL  myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef ALLOW_PTRACERS

      DO k = 1, Nr
       DO j = 1-OLy, sNy+OLy
        DO i = 1-OLx, sNx+OLx
            gplt1(i,j,k) = 
     &          - pTr_plt1(i,j,k)
     &            * PTRACERS_rFrag / 3600
            gplt2(i,j,k) =
     &            pTr_plt1(i,j,k)
     &            * PTRACERS_rFrag / 3600
C     this subroutine calculates the tendency so
C     do not   * dt   !!!!
        ENDDO
       ENDDO
      ENDDO

#endif /* ALLOW_PTRACERS */

       RETURN
       END
