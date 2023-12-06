#include "GCHEM_OPTIONS.h"

CBOP
C !ROUTINE: GCHEM_CALC_TENDENCY

C !INTERFACE: ==========================================================
      SUBROUTINE GCHEM_CALC_TENDENCY(
     I     myTime, myIter, myThid )

C !DESCRIPTION:
C In the case of GCHEM_SEPARATE_FORCING not defined,
C this subroutine computes the tracer tendencies due to a
C bio-geogchemistry or ecosystem model and stores them on an array
C gchemTendency, that will be incorporated into regular timestepping in
C in ptracers_intergrate.F
C The current example uses the CFC package, but his is meant to
C be replaced by anything that the user provides.

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
#ifdef ALLOW_CFC
# include "CFC_SIZE.h"
#endif
#ifdef ALLOW_SPOIL
# include "SPOIL_SIZE.h"
#endif

C !INPUT PARAMETERS: ===================================================
C  myTime               :: current time
C  myIter               :: current timestep
C  myThid               :: thread number
C     d1,d2       :: Layer thicknesses
C     w_r         :: rise velocity
C     weightA     :: weight for tracer @ level k-1
C     weightB     :: weight for tracer @ level  k
C     delTrac     :: Difference between tracer in each layer
C                    or more generally tracer in one specific cell
C                    when using an upwind scheme like in this routine
      _RL myTime
      _RL w_r, d1, d2
      _RL weightA(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL weightB(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL delTrac
      INTEGER myIter, myThid
CEOP

#ifdef GCHEM_ADD2TR_TENDENCY
C !LOCAL VARIABLES: ====================================================
C  i,j                  :: loop indices
C  iTr                  :: ptracer number
C  jTr                  :: corresponding tracer tendency number
C  bi,bj                :: tile indices
C  k                    :: vertical level
      INTEGER i,j
      INTEGER iTr, jTr
      INTEGER bi,bj,iMin,iMax,jMin,jMax,k

C     ========================
C     Here starts the modification:
C     ======================== 
C     19 July 2023, by Alan:
C     In my case I'm currently not wanting the tendency to be recalculated
C     so I'll just comment out the reinitialize part and calculate part.
C     So basically it does nothing lol.
C     =
C     24 OCT 2023, by Alan:
C     Implement the scheme as in "darwin_sinking.F"
C     An upwind adjustment with respect to the rising/sinking velocity
C     ** Not sure if it is OK to put an upwind scheme here **

C     PTRACERS_wRise    :: the rising/sinking velocity
C      w_r = PTRACERS_wRise
C
      IF (useGCHEM) THEN

       DO jTr = 1, GCHEM_tendTr_num
        DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)
C     gchemTendency is re-initialized here
          DO k = 1, Nr
            DO j = 1-OLy, sNy+OLy
             DO i = 1-OLx, sNx+OLx
              gchemTendency(i,j,k,bi,bj,jTr) = 0. _d 0
             ENDDO
            ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO

C     gchemTendency is computed here
       IF (GCHEM_tendTr_num .GE. 2) THEN
        DO jTr = 1, GCHEM_tendTr_num-1
         DO bj=myByLo(myThid),myByHi(myThid)
          DO bi=myBxLo(myThid),myBxHi(myThid)
           CALL PLASTIC_FRAG(
     I                       pTracer     (1-OLx,1-OLy,1,bi,bj,jTr),
     U                       gchemTendency(1-OLx,1-OLy,1,bi,bj,jTr),
     U                       gchemTendency(1-OLx,1-OLy,1,bi,bj,jTr+1),
     I                       bi, bj, iMin, iMax, jMin, jMax,
     I                       myTime, myIter, myThid )

          ENDDO
         ENDDO
        ENDDO
       ENDIF

       DO jTr = 1, GCHEM_tendTr_num
        DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)
          DO k = 1, Nr
            DO j = 1-OLy, sNy+OLy
             DO i = 1-OLx, sNx+OLx
              gchemTendency(i,j,k,bi,bj,jTr) = 
     &           gchemTendency(i,j,k,bi,bj,jTr)
     &           + gchemSource(i,j,k,bi,bj,jTr)
             ENDDO
            ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO

C     Apply mask
       DO jTr = 1, GCHEM_tendTr_num
        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
          DO k=1,Nr
           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
             IF (maskC(i,j,k,bi,bj).EQ.0.)
     &           gchemTendency(i,j,k,bi,bj,jTr)=0. _d 0
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO

      ENDIF
C     ========================
C     Here ends the modification.
C     ======================== 

#ifdef ALLOW_AUTODIFF
      IF ( useGCHEM ) THEN
#endif
C     Here is the place for code to compute bio-geochemical
C     tendency terms (sometimes referred to as source-minus-sink
C     terms). The tendencies are stored on gchemTendency, as show
C     in the CFC-example.
C
C     loop over tiles
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
C     define horizontal loop ranges
        iMin = 1
        iMax = sNx
        jMin = 1
        jMax = sNy
c       DO jTr = 1, GCHEM_tendTr_num
c        iTr = jTr + gchem_sepFTr_num
c        DO k = 1, Nr
c         DO j = 1-OLy, sNy+OLy
c          DO i = 1-OLx, sNx+OLx
c           gchemTendency(i,j,k,bi,bj,jTr) = your specific model
c          ENDDO
c         ENDDO
c        ENDDO
c       ENDDO

C------------------------
C chemical forcing      |
C------------------------
#ifdef ALLOW_SPOIL
        IF ( useSPOIL ) THEN
          CALL SPOIL_CALC_TENDENCY( bi, bj,
     I                              myTime, myIter, myThid )
        ENDIF
#endif /* ALLOW_SPOIL */

#ifdef ALLOW_CFC
        IF ( useCFC ) THEN
         iTr = CFC_pTr_i1
         jTr = iTr - gchem_sepFTr_num
         CALL CFC11_FORCING(
     I                       pTracer     (1-OLx,1-OLy,1,bi,bj,iTr),
     U                       gchemTendency(1-OLx,1-OLy,1,bi,bj,jTr),
     I                       bi, bj, iMin, iMax, jMin, jMax,
     I                       myTime, myIter, myThid )
         iTr = CFC_pTr_i1 + 1
         jTr = iTr - gchem_sepFTr_num
         CALL CFC12_FORCING(
     I                       pTracer     (1-OLx,1-OLy,1,bi,bj,iTr),
     U                       gchemTendency(1-OLx,1-OLy,1,bi,bj,jTr),
     I                       bi, bj, iMin, iMax, jMin, jMax,
     I                       myTime, myIter, myThid )
        ENDIF
#endif /* ALLOW_CFC */

C     end of tile-loops
       ENDDO
      ENDDO
#ifdef ALLOW_AUTODIFF
      ENDIF
#endif
#endif /* GCHEM_ADD2TR_TENDENCY */

      RETURN
      END
