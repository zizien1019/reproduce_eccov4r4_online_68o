C    Modified 15 Oct 2023 by Alan
C    Mixing of tracer is now done from bottom to top
C    For negatively buoyant tracers
C    And from top to bottom for positively buoyant tracers.
C-          Modified so that this routine only updates ptracers convectively
C-          To include the effect of buoyancy and settling velocity
C-          Modify convective_weights.F
C-          So that the coefficients A, B correspond to the correct physics.

#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_SIZE.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif


CBOP
C     !ROUTINE: CONVECTIVE_ADJUSTMENT
C     !INTERFACE:
      SUBROUTINE CONVECTIVE_ADJUSTMENT(
     I                      bi, bj, myTime, myIter, myThid )
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE CONVECTIVE_ADJUSTMENT
C     | o Driver for vertical mixing or similar parameterization
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
C     == Global data ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DYNVARS.h"
#include "GRID.h"
#ifdef ALLOW_TIMEAVE
#include "TIMEAVE_STATV.h"
#endif
#ifdef ALLOW_AUTODIFF_TAMC
#include "tamc.h"
#endif /* ALLOW_AUTODIFF_TAMC */

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     bi,bj  :: tile indices
C     myTime :: Current time in simulation
C     myIter :: Current iteration in simulation
C     myThid :: My Thread Id number
      INTEGER bi,bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid

#ifdef INCLUDE_CONVECT_CALL

C     !FUNCTIONS:
      EXTERNAL DIFFERENT_MULTIPLE
      LOGICAL  DIFFERENT_MULTIPLE

C     !LOCAL VARIABLES:
C     == Local variables ==
C     iMin,iMax,jMin,jMax :: computation domain
C     i,j,k        :: Loop counters
C     rhoKm1, rhoK :: Density at adjacent levels (common ref. level)
C     ConvectCount :: Convection mixing freq. counter.
      INTEGER iMin,iMax,jMin,jMax
      INTEGER i, j, k, kTop, kBottom, kDir, deltaK
#ifdef ALLOW_AUTODIFF_TAMC
      INTEGER tkey, kkey
#endif
      _RL rhoKm1(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL rhoK  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL ConvectCount(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL weightA(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL weightB(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
CEOP

C--   Check to see if should convect now
      IF ( DIFFERENT_MULTIPLE(cAdjFreq,myTime,deltaTClock)
     &   ) THEN

C--   Define computation domain
        iMin=1-OLx
        iMax=sNx+OLx
        jMin=1-OLy
        jMax=sNy+OLy

C--   Initialise counters
        kTop    = 0
        kBottom = 0
        kDir    = 0
        deltaK  = 0

C-      Initialisation of Convection Counter
        DO k=1,Nr
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           ConvectCount(i,j,k) = 0.
          ENDDO
         ENDDO
        ENDDO

#ifdef ALLOW_AUTODIFF_TAMC
        tkey = bi + (bj-1)*nSx + (ikey_dynamics-1)*nSx*nSy
#endif /* ALLOW_AUTODIFF_TAMC */

        IF ( rkSign*gravitySign .GT. 0. ) THEN
C-          <=> usingZCoords:
         kTop    =  2
         kBottom = Nr
         kDir    =  1
         deltaK  = -1
        ENDIF

C-          Density of k-1 layer (above W(k)) reference to k-1 T-level
            CALL FIND_RHO_2D(
     &         iMin, iMax, jMin, jMax, k+deltaK,
     &         theta(1-OLx,1-OLy,k-1,bi,bj),
     &         salt (1-OLx,1-OLy,k-1,bi,bj),
     &         rhoKm1,
     &         k-1, bi, bj, myThid )
   
C-          Density of k layer (below W(k)) reference to k-1 T-level.
            CALL FIND_RHO_2D(
     &         iMin, iMax, jMin, jMax, k+deltaK,
     &         theta(1-OLx,1-OLy,k,bi,bj),
     &         salt (1-OLx,1-OLy,k,bi,bj),
     &         rhoK,
     &         k, bi, bj, myThid )

#ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE rhoKm1 = comlev1_bibj_k, key = kkey, kind = isbyte
CADJ STORE rhoK   = comlev1_bibj_k, key = kkey, kind = isbyte
#endif /* ALLOW_AUTODIFF_TAMC */

        IF (PTRACERS_wRise .GT. 0. _d 0) THEN
C--       If positively buoyant loop from top to bottom
C--       Loop over all *interior* layers
                DO k=kTop,kBottom,kDir

#ifdef ALLOW_AUTODIFF_TAMC
         kkey = k + (tkey-1)*Nr
C     It is important that the two k-levels of these fields are stored
C     in one statement because otherwise taf will only store one, which
C     is wrong (i.e. was wrong in previous versions).
CADJ STORE theta(:,:,k-1,bi,bj), theta(:,:,k,bi,bj) =
CADJ &     comlev1_bibj_k, key = kkey, kind = isbyte
CADJ STORE salt(:,:,k-1,bi,bj), salt(:,:,k,bi,bj) =
CADJ &     comlev1_bibj_k, key = kkey, kind = isbyte
#endif /* ALLOW_AUTODIFF_TAMC */


C-          Pre-calculate mixing weights for interface k
                CALL CONVECTIVE_WEIGHTS(
     I        bi,bj,k,rhoKm1,rhoK,
     O        weightA,weightB,ConvectCount,
     I        myThid)

C-          Convectively mix heat across interface k
C         CALL CONVECTIVELY_MIXTRACER(
C     I                               bi,bj,k,weightA,weightB,
C     U                               theta,
C     I                               myThid)

C-          Convectively mix salt across interface k
C         CALL CONVECTIVELY_MIXTRACER(
C     I                               bi,bj,k,weightA,weightB,
C     U                               salt,
C     I                               myThid)

#ifdef ALLOW_PTRACERS
C-          Convectively mix passive tracers across interface k
                IF ( usePTRACERS ) THEN
                CALL PTRACERS_CONVECT(
     I                          bi,bj,k,weightA,weightB,myThid)
                ENDIF
#endif /* ALLOW_PTRACERS */

C--       End DO k=2,Nr
                ENDDO
        ELSE
C--       If negatively buoyany loop from bottom to top
C--       Loop over all *interior* layers
                DO k=kBottom,kTop,-kDir

#ifdef ALLOW_AUTODIFF_TAMC
         kkey = k + (tkey-1)*Nr
C     It is important that the two k-levels of these fields are stored
C     in one statement because otherwise taf will only store one, which
C     is wrong (i.e. was wrong in previous versions).
CADJ STORE theta(:,:,k-1,bi,bj), theta(:,:,k,bi,bj) =
CADJ &     comlev1_bibj_k, key = kkey, kind = isbyte
CADJ STORE salt(:,:,k-1,bi,bj), salt(:,:,k,bi,bj) =
CADJ &     comlev1_bibj_k, key = kkey, kind = isbyte
#endif /* ALLOW_AUTODIFF_TAMC */


C-          Pre-calculate mixing weights for interface k
                CALL CONVECTIVE_WEIGHTS(
     I        bi,bj,k,rhoKm1,rhoK,
     O        weightA,weightB,ConvectCount,
     I        myThid)

C-          Convectively mix heat across interface k
C         CALL CONVECTIVELY_MIXTRACER(
C     I                               bi,bj,k,weightA,weightB,
C     U                               theta,
C     I                               myThid)

C-          Convectively mix salt across interface k
C         CALL CONVECTIVELY_MIXTRACER(
C     I                               bi,bj,k,weightA,weightB,
C     U                               salt,
C     I                               myThid)

#ifdef ALLOW_PTRACERS
C-          Convectively mix passive tracers across interface k
                IF ( usePTRACERS ) THEN
                CALL PTRACERS_CONVECT(
     I                          bi,bj,k,weightA,weightB,myThid)
                ENDIF
#endif /* ALLOW_PTRACERS */

C--       End DO k=Nr,2
                ENDDO
C--       Endif (wRise<=0)
        ENDIF

#ifdef ALLOW_TIMEAVE
        IF (myIter.NE.nIter0 .AND. taveFreq.GT.0.) THEN
         CALL TIMEAVE_CUMUL_1T(ConvectCountTave, ConvectCount,
     I                         Nr, deltaTClock, bi, bj, myThid)
        ENDIF
#endif /* ALLOW_TIMEAVE */

#ifdef ALLOW_DIAGNOSTICS
        IF ( myIter.NE.nIter0 .AND. useDiagnostics ) THEN
         CALL DIAGNOSTICS_FILL( ConvectCount, 'CONVADJ ',
     I                          0, Nr, 2, bi, bj, myThid )
        ENDIF
#endif /*  ALLOW_DIAGNOSTICS */

C--   End IF (DIFFERENT_MULTIPLE)
      ENDIF

#endif /* INCLUDE_CONVECT_CALL */

      RETURN
      END
