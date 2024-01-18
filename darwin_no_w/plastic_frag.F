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
      IMPLICIT NONE
C     == GLobal variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "PTRACERS_FIELDS.h"
# include "PTRACERS_SIZE.h"
# include "PTRACERS_PARAMS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     pTr_plt1  :: ocean larger plastic concentration
C     gplt1     :: larger plastic tendency
C     gplt2     :: smaller plastic tendency
C     bi, bj     :: current tile indices
C     iMin,iMax  :: computation domain, 1rst index bounds
C     jMin,jMax  :: computation domain, 2nd  index bounds
C     myTime     :: current time in simulation
C     myIter     :: current iteration number
C     myThid     :: my Thread Id number
      _RL  pTr_plt1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  gplt1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  gplt2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      INTEGER bi, bj
      INTEGER iMin, iMax, jMin, jMax
      _RL  myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef ALLOW_PTRACERS

      DO jTr = 1, GCHEM_tendTr_num
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO k = 1, Nr
          DO j = 1-OLy, sNy+OLy
           DO i = 1-OLx, sNx+OLx
            glpt1(i,j,k,bi,bj,jTr) = 
     &            pTr_plt1(i,j,k,bi,bj,jTr)
     &            * PTRACERS_rFrag * delTrac / 3600
            glpt2(i,j,k,bi,bj,jTr) =
     &            pTr_plt1(i,j,k,bi,bj,jTr)
     &            * PTRACERS_rFrag * delTrac / 3600
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

#endif /* ALLOW_PTRACERS */

       RETURN
       END
