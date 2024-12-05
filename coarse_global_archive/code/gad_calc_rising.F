#include "GAD_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

CBOP
C !ROUTINE: GAD_CALC_RHS

C !INTERFACE: ==========================================================
      SUBROUTINE GAD_CALC_RISING(
     I           bi,bj,iMin,iMax,jMin,jMax,k,kM1,kUp,kDown,
     I           xA, yA, maskUp, uFld, vFld, wFld,
     I           uTrans, vTrans, rTrans, rTransKp1,
     I           diffKh, diffK4, KappaR, diffKr4, TracerN, TracAB,
     I           deltaTLev, trIdentity,
     I           advectionSchArg, vertAdvecSchArg,
     I           calcAdvection, implicitAdvection, applyAB_onTracer,
     I           trUseDiffKr4, trUseGMRedi, trUseKPP, trUseSmolHack,
     O           fZon, fMer,
     U           fVerT, gTracer,
     I           myTime, myIter, myThid )

C !DESCRIPTION:
C Calculates the tendency of a tracer due to advection and diffusion.
C It calculates the fluxes in each direction indepentently and then
C sets the tendency to the divergence of these fluxes. The advective
C fluxes are only calculated here when using the linear advection schemes
C otherwise only the diffusive and parameterized fluxes are calculated.
C
C Contributions to the flux are calculated and added:
C \begin{equation*}
C {\bf F} = {\bf F}_{adv} + {\bf F}_{diff} +{\bf F}_{GM} + {\bf F}_{KPP}
C \end{equation*}
C
C The tendency is the divergence of the fluxes:
C \begin{equation*}
C G_\theta = G_\theta + \nabla \cdot {\bf F}
C \end{equation*}
C
C The tendency is assumed to contain data on entry.

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SURFACE.h"
#include "GAD.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_PARAMS.h"
#endif /* ALLOW_AUTODIFF */

C !INPUT PARAMETERS: ===================================================
C bi, bj           :: tile indices
C iMin, iMax       :: for called routines, to get valid output "gTracer"
C jMin, jMax       ::                      over this range of indices
C k                :: vertical index
C kM1              :: =k-1 for k>1, =1 for k=1
C kUp              :: index into 2 1/2D array, toggles between 1|2
C kDown            :: index into 2 1/2D array, toggles between 2|1
C xA, yA           :: areas of X and Y face of tracer cells
C maskUp           :: 2-D array for mask at W points
C uFld, vFld, wFld :: Local copy of velocity field (3 components)
C uTrans, vTrans   :: 2-D arrays of volume transports at U,V points
C rTrans           :: 2-D arrays of volume transports at W points
C rTransKp1        :: 2-D array of volume trans at W pts, interf k+1
C diffKh           :: horizontal diffusion coefficient
C diffK4           :: horizontal bi-harmonic diffusion coefficient
C KappaR           :: 2-D array for vertical diffusion coefficient, interf k
C diffKr4          :: 1-D array for vertical bi-harmonic diffusion coefficient
C TracerN          :: tracer field @ time-step n (Note: only used
C                     if applying AB on tracer field rather than on tendency gTr)
C TracAB           :: current tracer field (@ time-step n if applying AB on gTr
C                     or extrapolated fwd in time to n+1/2 if applying AB on Tr)
C trIdentity       :: tracer identifier (required for KPP,GM)
C advectionSchArg  :: advection scheme to use (Horizontal plane)
C vertAdvecSchArg  :: advection scheme to use (Vertical direction)
C calcAdvection    :: =False if Advec computed with multiDim scheme
C implicitAdvection:: =True if vertical Advec computed implicitly
C applyAB_onTracer :: apply Adams-Bashforth on Tracer (rather than on gTr)
C trUseDiffKr4     :: true if this tracer uses vertical bi-harmonic diffusion
C trUseGMRedi      :: true if this tracer uses GM-Redi
C trUseKPP         :: true if this tracer uses KPP
C trUseSmolHack    :: true if this tracer uses Smolarkiewicz-Hack to remain > 0
C myTime           :: current time
C myIter           :: iteration number
C myThid           :: thread number
      INTEGER bi,bj,iMin,iMax,jMin,jMax
      INTEGER k,kUp,kDown,kM1
      _RS xA    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RS yA    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RS maskUp(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL uFld  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL vFld  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL wFld  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL uTrans(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL vTrans(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL rTrans(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL rTransKp1(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL diffKh, diffK4
      _RL KappaR(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL diffKr4(Nr)
      _RL TracerN(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL TracAB (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL deltaTLev(Nr)
      INTEGER trIdentity
      INTEGER advectionSchArg, vertAdvecSchArg
      LOGICAL calcAdvection
      LOGICAL implicitAdvection, applyAB_onTracer
      LOGICAL trUseDiffKr4, trUseGMRedi, trUseKPP, trUseSmolHack
      _RL     myTime
      INTEGER myIter, myThid

C !OUTPUT PARAMETERS: ==================================================
C gTracer          :: tendency array
C fZon             :: zonal flux
C fMer             :: meridional flux
C fVerT            :: 2 1/2D arrays for vertical advective flux
      _RL gTracer(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL fZon  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL fMer  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL fVerT (1-OLx:sNx+OLx,1-OLy:sNy+OLy,2)

C !FUNCTIONS:       ====================================================
#ifdef ALLOW_DIAGNOSTICS
      CHARACTER*4 GAD_DIAG_SUFX
      EXTERNAL    GAD_DIAG_SUFX
#endif /* ALLOW_DIAGNOSTICS */

C !LOCAL VARIABLES: ====================================================
C i,j              :: loop indices
C advectionScheme  :: local copy of routine argument advectionSchArg
C vertAdvecScheme  :: local copy of routine argument vertAdvecSchArg
C df4              :: used for storing del^2 T for bi-harmonic term
C af               :: advective flux
C df               :: diffusive flux
C localT           :: local copy of tracer field
C locABT           :: local copy of (AB-extrapolated) tracer field
      INTEGER i,j
      INTEGER advectionScheme, vertAdvecScheme
      _RS maskLocW(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RS maskLocS(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL df4   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL af    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL df    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL localT(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL locABT(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL advFac, rAdvFac
#ifdef GAD_SMOLARKIEWICZ_HACK
      _RL outFlux, trac, gTrFac
#endif
#ifdef ALLOW_DIAGNOSTICS
      CHARACTER*8 diagName
      CHARACTER*4 diagSufx
#endif
CEOP

C     make local copies to be tampered with if necessary
      advectionScheme = advectionSchArg
      vertAdvecScheme = vertAdvecSchArg
#ifdef ALLOW_AUTODIFF
#ifdef ALLOW_AUTODIFF_TAMC
CADJ INIT gad_local_tape = COMMON, 1
C     This store directive just suppresses a recomputation warning.
C     TAF generates an extra field with or without this directive.
CADJ STORE fvert = gad_local_tape
#endif
C--   only the kUp part of fverT is set in this subroutine
C--   the kDown is still required
      fVerT(1,1,kDown) = fVerT(1,1,kDown)
C
      IF ( inAdMode .AND. useApproxAdvectionInAdMode ) THEN
C     In AD-mode, we change non-linear, potentially unstable AD advection
C     schemes to linear schemes with more stability. So far only DST3 with
C     flux limiting is replaced by DST3 without flux limiting, but any
C     combination is possible.
       IF ( advectionSchArg.EQ.ENUM_DST3_FLUX_LIMIT )
     &      advectionScheme = ENUM_DST3
       IF ( vertAdvecSchArg.EQ.ENUM_DST3_FLUX_LIMIT )
     &      vertAdvecScheme = ENUM_DST3
C     here is room for more advection schemes as this becomes necessary
      ENDIF
#endif /* ALLOW_AUTODIFF */

#ifdef ALLOW_DIAGNOSTICS
C--   Set diagnostic suffix for the current tracer
      IF ( useDiagnostics ) THEN
        diagSufx = GAD_DIAG_SUFX( trIdentity, myThid )
      ENDIF
#endif

      advFac  = 0. _d 0
      IF (calcAdvection) advFac = 1. _d 0
      rAdvFac = rkSign*advFac
      IF (implicitAdvection) rAdvFac = rkSign

      DO j=1-OLy,sNy+OLy
       DO i=1-OLx,sNx+OLx
        fZon(i,j)      = 0. _d 0
        fMer(i,j)      = 0. _d 0
        fVerT(i,j,kUp) = 0. _d 0
        df(i,j)        = 0. _d 0
        df4(i,j)       = 0. _d 0
       ENDDO
      ENDDO

C--   Make local copy of tracer array
      IF ( applyAB_onTracer ) THEN
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          localT(i,j)=TracerN(i,j,k)
          locABT(i,j)= TracAB(i,j,k)
         ENDDO
        ENDDO
      ELSE
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          localT(i,j)=TracerN(i,j,k)
          locABT(i,j)=TracerN(i,j,k)
         ENDDO
        ENDDO
      ENDIF

C--   Compute vertical flux fVerT(kUp) at interface k (between k-1 & k):
C-    Advective flux in R
#ifdef ALLOW_AIM
C- a hack to prevent Water-Vapor vert.transport into the stratospheric level Nr
      IF (calcAdvection .AND. .NOT.implicitAdvection .AND. k.GE.2 .AND.
     &     (.NOT.useAIM .OR. trIdentity.NE.GAD_SALINITY .OR. k.LT.Nr)
     &   ) THEN
#else
      IF (calcAdvection .AND. .NOT.implicitAdvection .AND. k.GE.2) THEN
#endif
       IF ( applyAB_onTracer ) THEN
C-    Compute vertical advective flux in the interior using TracAB:
        IF ( vertAdvecScheme.EQ.ENUM_CENTERED_2ND ) THEN
           CALL GAD_C2_ADV_R( bi,bj,k, rTrans, TracAB, af, myThid )
        ELSEIF ( vertAdvecScheme.EQ.ENUM_UPWIND_1RST
     &          .OR. vertAdvecScheme.EQ.ENUM_DST2 ) THEN
           CALL GAD_DST2U1_ADV_R( bi,bj,k,vertAdvecScheme,deltaTLev(k),
     I              rTrans, wFld, TracAB,
     O              af, myThid )
        ELSEIF ( vertAdvecScheme.EQ.ENUM_FLUX_LIMIT ) THEN
           CALL GAD_FLUXLIMIT_ADV_R( bi,bj,k, deltaTLev(k),
     I              rTrans, wFld, TracAB,
     O              af, myThid )
        ELSEIF ( vertAdvecScheme.EQ.ENUM_UPWIND_3RD ) THEN
           CALL GAD_U3_ADV_R( bi,bj,k, rTrans, TracAB, af, myThid )
        ELSEIF (vertAdvecScheme.EQ.ENUM_CENTERED_4TH) THEN
           CALL GAD_C4_ADV_R( bi,bj,k, rTrans, TracAB, af, myThid )
        ELSEIF ( vertAdvecScheme.EQ.ENUM_DST3 ) THEN
           CALL GAD_DST3_ADV_R( bi,bj,k, deltaTLev(k),
     I              rTrans, wFld, TracAB,
     O              af, myThid )
        ELSEIF ( vertAdvecScheme.EQ.ENUM_DST3_FLUX_LIMIT ) THEN
           CALL GAD_DST3FL_ADV_R( bi,bj,k, deltaTLev(k),
     I              rTrans, wFld, TracAB,
     O              af, myThid )
#ifndef ALLOW_AUTODIFF
        ELSEIF ( vertAdvecScheme.EQ.ENUM_OS7MP ) THEN
           CALL GAD_OS7MP_ADV_R( bi,bj,k, deltaTLev(k),
     I              rTrans, wFld, TracAB,
     O              af, myThid )
#endif
        ELSE
          STOP 'GAD_CALC_RHS: Bad vertAdvecScheme (R)'
        ENDIF
       ELSE
C-    Compute vertical advective flux in the interior using TracerN:
        IF ( vertAdvecScheme.EQ.ENUM_CENTERED_2ND ) THEN
           CALL GAD_C2_ADV_R( bi,bj,k, rTrans, TracerN, af, myThid )
        ELSEIF ( vertAdvecScheme.EQ.ENUM_UPWIND_1RST
     &          .OR. vertAdvecScheme.EQ.ENUM_DST2 ) THEN
           CALL GAD_DST2U1_ADV_R( bi,bj,k,vertAdvecScheme,deltaTLev(k),
     I              rTrans, wFld, TracerN,
     O              af, myThid )
        ELSEIF ( vertAdvecScheme.EQ.ENUM_FLUX_LIMIT ) THEN
           CALL GAD_FLUXLIMIT_ADV_R( bi,bj,k, deltaTLev(k),
     I              rTrans, wFld, TracerN,
     O              af, myThid )
        ELSEIF ( vertAdvecScheme.EQ.ENUM_UPWIND_3RD ) THEN
           CALL GAD_U3_ADV_R( bi,bj,k, rTrans, TracerN, af, myThid )
        ELSEIF ( vertAdvecScheme.EQ.ENUM_CENTERED_4TH ) THEN
           CALL GAD_C4_ADV_R( bi,bj,k, rTrans, TracerN, af, myThid )
        ELSEIF( vertAdvecScheme.EQ.ENUM_DST3 ) THEN
           CALL GAD_DST3_ADV_R( bi,bj,k, deltaTLev(k),
     I              rTrans, wFld, TracerN,
     O              af, myThid )
        ELSEIF ( vertAdvecScheme.EQ.ENUM_DST3_FLUX_LIMIT ) THEN
           CALL GAD_DST3FL_ADV_R( bi,bj,k, deltaTLev(k),
     I              rTrans, wFld, TracerN,
     O              af, myThid )
#ifndef ALLOW_AUTODIFF
        ELSEIF ( vertAdvecScheme.EQ.ENUM_OS7MP ) THEN
           CALL GAD_OS7MP_ADV_R( bi,bj,k, deltaTLev(k),
     I              rTrans, wFld, TracerN,
     O              af, myThid )
#endif
        ELSE
          STOP 'GAD_CALC_RHS: Bad vertAdvecScheme (R)'
        ENDIF
       ENDIF
C-     add the advective flux to fVerT
       DO j=1-OLy,sNy+OLy
        DO i=1-OLx,sNx+OLx
          fVerT(i,j,kUp) = fVerT(i,j,kUp) + af(i,j)*maskInC(i,j,bi,bj)
        ENDDO
       ENDDO
#ifdef ALLOW_DIAGNOSTICS
       IF ( useDiagnostics ) THEN
          diagName = 'ADVr'//diagSufx
          CALL DIAGNOSTICS_FILL( af, diagName, k,1, 2,bi,bj, myThid )
C- note: needs to explicitly increment the counter since DIAGNOSTICS_FILL
C        does it only if k=1 (never the case here)
          IF ( k.EQ.2 ) CALL DIAGNOSTICS_COUNT(diagName,bi,bj,myThid)
#ifdef ALLOW_LAYERS
          IF ( useLayers ) THEN
            CALL LAYERS_FILL(af,trIdentity,'AFR',k,1,2,bi,bj,myThid)
          ENDIF
#endif /* ALLOW_LAYERS */
       ENDIF
#endif
      ENDIF


C--   Divergence of fluxes
C--   This routine by Zih-En (Alan) Tseng
C--   27/Apr/2024
C--   for vertical flux due to rising/sinking

C     Anelastic: scale vertical fluxes by rhoFac and leave Horizontal fluxes unchanged
C     for Stevens OBC: keep only vertical diffusive contribution on boundaries
      DO j=1-OLy,sNy+OLy-1
       DO i=1-OLx,sNx+OLx-1
        gTracer(i,j,k) = gTracer(i,j,k)
     &   -_recip_hFacC(i,j,k,bi,bj)*recip_drF(k)
     &   *recip_rA(i,j,bi,bj)*recip_deepFac2C(k)*recip_rhoFacC(k)
     &   *( 
     &     +(fVerT(i,j,kDown)-fVerT(i,j,kUp))*rkSign
     &    )
       ENDDO
      ENDDO

#ifdef ALLOW_DEBUG
      IF ( debugLevel .GE. debLevC
     &   .AND. trIdentity.EQ.GAD_TEMPERATURE
     &   .AND. k.EQ.2 .AND. myIter.EQ.1+nIter0
     &   .AND. nPx.EQ.1 .AND. nPy.EQ.1
     &   .AND. useCubedSphereExchange ) THEN
        CALL DEBUG_CS_CORNER_UV( ' fZon,fMer from GAD_CALC_RHS',
     &             fZon,fMer, k, standardMessageUnit,bi,bj,myThid )
      ENDIF
#endif /* ALLOW_DEBUG */

      RETURN
      END
