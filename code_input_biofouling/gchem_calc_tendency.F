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

C     tr_projA    :: tracer values projected on A-boundaries
C     k_bolz      :: bolzmann constant == 1.3806e-23
C     D_pl(Nr)    :: Brownian motion related diffusivity
C     D_A(Nr)     :: Brownian motion related diffusivity
C     b_brown(Nr) :: Brownian motion encounter kernel rate
C     w_r(Nr)     :: falling velocity of a particle
C     r_A         :: algae cell radius
C     Vbf         :: biofilm volume on one particle
C     Chla(Nr)    :: ambient Chl-a concentration
C     Chla_by_C(Nr)  :: ambient Carbon concentration
C     A_A(Nr)     :: ambient algae cell concentration
C     I_mag(Nr)   :: magnitude of radiation
C     Inc(Nr)     :: amount of incident radiation
C     mu_opt(Nr)  :: optimal growth rate in (1/day)
C     mu(Nr)      :: growth rate in (1/day)
C     phi(Nr)     :: temperature influence

      _RL myTime
      _RL d1, d2
      _RL weightA(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL weightB(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL delTrac

      _RL tr_projA, k_bolz, r_A, Vbf
      _RL D_pl(Nr), D_A(Nr), b_brown(Nr), w_r(Nr)
      _RL Chla_by_C(Nr), Chla(Nr), A_A(Nr)
      _RL I_mag(Nr), Inc(Nr)
      _RL mu_opt(Nr), mu(Nr), phi(Nr)

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

#ifdef ALLOW_AUTODIFF
      IF ( useGCHEM ) THEN
#endif

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
C     =
C     May 2024, by Alan:
C     sinking is now fully incorporated in ptracers package
C     for more detail see gad_calc_rhs and gad_calc_rising
C
      DO jTr = 1, GCHEM_tendTr_num
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
C     gchemTendency is re-initialized here
         DO k = 1, Nr
          DO j = 1-OLy, sNy+OLy
           DO i = 1-OLx, sNx+OLx
             gchemTendency(i,j,k,bi,bj,jTr) = 0. _d 0
             bio_dAdt(i,j,k,bi,bj,jTr) = 0. _d 0
             bio_flux(i,j,k,bi,bj,jTr) = 0. _d 0
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      DO jTr = 1, 2
C     Only release clean particles
C     jTr = 1     :: all clean without biofouling
C     jTr = 2:end :: biofouling & reacting tracers
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO k = 1, Nr
          DO j = 1-OLy, sNy+OLy
           DO i = 1-OLx, sNx+OLx
C     source term from specified source file
            IF (gchemSource(i,j,k,bi,bj,jTr).NE.0) THEN
            gchemTendency(i,j,k,bi,bj,jTr) = 
     &        gchemTendency(i,j,k,bi,bj,jTr)
     &        + gchemSource(i,j,k,bi,bj,jTr)
            ENDIF
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

C     Here is the place for code to compute bio-geochemical
C     tendency terms (sometimes referred to as source-minus-sink
C     terms). The tendencies are stored on gchemTendency, as show
C     in the CFC-example.
C     =
C     Feb 2025, by Alan:
C     Biofouling based on (KOOI 2017)

! C   % Boltzman constant           :: k_bolz
! C   % radius of an algae cell     :: r_A
      k_bolz = 1.3806E-23
      r_A = 3.6E-6

C     loop over tiles
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
C     define horizontal loop ranges
        iMin = 1
        iMax = sNx
        jMin = 1
        jMax = sNy

        DO jTr = 2, GCHEM_tendTr_num
C     Only do for tracers used with biofouling
         iTr = jTr + gchem_sepFTr_num

! C         % species specific          :: Vbf
         Vbf = 2.E-16 * (jTr * 1. _d 0 - 2.) * PTRACERS_delA 
     &      * 4. * PI * (PTRACERS_d_0/2.)**2

         DO j = 1-OLy, sNy+OLy
          DO i = 1-OLx, sNx+OLx

           DO k = 1, Nr
! C         % profiles

! C         % vertical coordinate       :: z_c == rC(Nr)
! C         % potential temperature     :: theta(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
! C           !!!! actual temperature = potential temperature + tref !!!  
! C         % salinity                  :: salt(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
! C         % water density             :: rhoInSitu(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
! C         % water dynamic viscosity   :: use 1000
! C         % Chl-a                     :: Chla(Nr)
            Chla(k) = 0.
            IF (rC(k) .GT. -99.9) THEN
            Chla(k) = 0.151 * (
     &            0.533 + 1.72E-3 * rC(k)
     &            + 1.194 * exp(-((rC(k)+92.01)/43.46)**2)
     &       )
            ENDIF

! C         % magnitude of radiation    :: I_mag(Nr)
! C         % radiation                 :: Inc(Nr)
            I_mag(k) = 1.2E8 * exp( (0.2+0.02*Chla(k)) * rC(k))
!
! C         % reaction step             :: sA == (jTr * 1. _d 0)
! C         % encounter kernel rate
! C         % diffusivity               :: D_pl(Nr)  
! C         % diffusivity               :: D_A(Nr)
! C         % velocity                  :: w_r(Nr)
! C         % kernel rate               :: b_brown(Nr)
! C         % ambient calbon            :: C(Nr)   (mg/m^3)
! C         % ambient algae             :: A_A(Nr) (#cell/m^3) 
            D_pl(k) = k_bolz * (theta(i,j,k,bi,bj)+tRef(k)+273.16) 
     &            / (6. * PI * 0.001 * PTRACERS_d_p(jTr) / 2.)
            D_A(k) = k_bolz * (theta(i,j,k,bi,bj)+tRef(k)+273.16) 
     &            / (6. * PI * 0.001 * r_A)
            w_r(k) = abs(PTRACERS_rho_p(jTr) - rhoInSitu(i,j,k,bi,bj))
     &            * 9.81 * PTRACERS_d_p(jTr)**2 / 18. / 0.001
            b_brown(k) = 4. * PI * (D_pl(k) + D_A(k)) * 
     &            (PTRACERS_d_p(jTr) / 2. + r_A)
     &            + 0.5 * PI * (PTRACERS_d_p(jTr)/2)**2 * w_r(k)
     &            + 1.3 * (1.7E5 / 86400.) * 
     &            (PTRACERS_d_p(jTr)/2. + r_A)**3
! C         % ambient algae concentration
            Inc(k) = I_mag(k) 
     &            * DCOS( 2. * PI * DBLE(myTime) / 86400. )
     &            * DCOS( yC(i,j,bi,bj) / 180. * PI )
            IF (Inc(k) .LT. 0) THEN
                  Inc(k) = 0.
            ENDIF
            Chla_by_C(k) = 0.003 
     &            + 0.0154 * exp(0.05 * (theta(i,j,k,bi,bj)+tRef(k))
     &            - 0.059 * Inc(k) / 1.E6)
            A_A(k) = Chla(k) / Chla_by_C(k) / (2726.E-9)
            
! C         % mu_opt(Nr) :: optimal growth rate in (1/day)
! C         % mu(Nr)     :: growth rate in (1/day)
! C         % phi(Nr)    :: temperature influence
            mu_opt(k) = 1.85 * Inc(k) / 86400. / 
     &            (Inc(k) / 86400. + 1.85 / 0.12)

            phi(k) = (theta(i,j,k,bi,bj)+tRef(k) - 33.3) 
     &            * (theta(i,j,k,bi,bj)+tRef(k) - 0.2)**2 / 26.5
     &            / ( 26.5 * (theta(i,j,k,bi,bj)+tRef(k) - 26.7) 
     &            + 6.6 * (33.5 - 2. * (theta(i,j,k,bi,bj)+tRef(k))))
            IF ( ((theta(i,j,k,bi,bj)+tRef(k)) .GT. 33.3) .OR.
     &           ((theta(i,j,k,bi,bj)+tRef(k)) .LT. 0.2 ) ) THEN
                  phi(k) = 0.
            ENDIF
            mu(k) = mu_opt(k) * phi(k)
            

            bio_dAdt(i,j,k,bi,bj,jTr) = 
     &       b_brown(k) * A_A(k) / (4. * PI * (PTRACERS_d_0/2.)**2)
     &       + mu(k) / 86400. * (jTr * 1. _d 0 - 2.) * PTRACERS_delA
     &       - 0.39 / 86400. * (jTr * 1. _d 0 - 2.) * PTRACERS_delA
     &       - 2. _d 0 **((theta(i,j,k,bi,bj)+tRef(k))/10. - 2.)
     &       * 0.1 / 86400. * (jTr * 1. _d 0 - 2.) * PTRACERS_delA

           ENDDO
          ENDDO
         ENDDO
        ENDDO


        DO jTr = 2, GCHEM_tendTr_num-1
! C   % reactive flux due to biofouling
! C   % located at A-boundaries
! C   % with zero-flux boundary condition
! C   % there's only need to compute through 2 to end-1
         iTr = jTr + gchem_sepFTr_num
         DO j = 1-OLy, sNy+OLy
          DO i = 1-OLx, sNx+OLx
           DO k = 1, Nr

            IF ((bio_dAdt(i,j,k,bi,bj,jTr)
     &         + bio_dAdt(i,j,k,bi,bj,jTr+1)) .GT. 0.) THEN
             tr_projA = pTracer(i,j,k,bi,bj,jTr)
            ELSE
             tr_projA = pTracer(i,j,k,bi,bj,jTr+1)
            ENDIF

            bio_flux(i,j,k,bi,bj,jTr) = 
     &   (bio_dAdt(i,j,k,bi,bj,jTr)+bio_dAdt(i,j,k,bi,bj,jTr+1))
     &   * tr_projA / 2.

           ENDDO
          ENDDO
         ENDDO
        ENDDO

        DO jTr = 2, GCHEM_tendTr_num
         iTr = jTr + gchem_sepFTr_num
         DO j = 1-OLy, sNy+OLy
          DO i = 1-OLx, sNx+OLx
           DO k = 1, Nr
      
               gchemTendency(i,j,k,bi,bj,jTr) =
     &      gchemTendency(i,j,k,bi,bj,jTr) -
     &      (bio_flux(i,j,k,bi,bj,jTr) - bio_flux(i,j,k,bi,bj,jTr-1))
     &      / PTRACERS_delA

           ENDDO
          ENDDO
         ENDDO
        ENDDO

C     ========================
C     Here ends the modification.
C     ========================


C------------------------
C chemical forcing      |
C------------------------
C they are now outside jTr loop

! #ifdef ALLOW_SPOIL
!         IF ( useSPOIL ) THEN
!           CALL SPOIL_CALC_TENDENCY( bi, bj,
!      I                              myTime, myIter, myThid )
!         ENDIF
! #endif /* ALLOW_SPOIL */

! #ifdef ALLOW_CFC
!         IF ( useCFC ) THEN
!          iTr = CFC_pTr_i1
!          jTr = iTr - gchem_sepFTr_num
!          CALL CFC11_FORCING(
!      I                       pTracer     (1-OLx,1-OLy,1,bi,bj,iTr),
!      U                       gchemTendency(1-OLx,1-OLy,1,bi,bj,jTr),
!      I                       bi, bj, iMin, iMax, jMin, jMax,
!      I                       myTime, myIter, myThid )
!          iTr = CFC_pTr_i1 + 1
!          jTr = iTr - gchem_sepFTr_num
!          CALL CFC12_FORCING(
!      I                       pTracer     (1-OLx,1-OLy,1,bi,bj,iTr),
!      U                       gchemTendency(1-OLx,1-OLy,1,bi,bj,jTr),
!      I                       bi, bj, iMin, iMax, jMin, jMax,
!      I                       myTime, myIter, myThid )
!         ENDIF
! #endif /* ALLOW_CFC */


C     end of tile-loops
       ENDDO
      ENDDO
#ifdef ALLOW_AUTODIFF
      ENDIF
#endif
#endif /* GCHEM_ADD2TR_TENDENCY */

      RETURN
      END
