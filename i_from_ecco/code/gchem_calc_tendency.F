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
CCCCCCCC include exf to gain access to swflux (negative downward)  
#include "EXF_FIELDS.h"
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
C     r_A         :: plankton cell radius
C     Vbf         :: biofilm volume on one particle

C     Chla(Nr)    :: ambient Chl-a concentration
C     Chla_by_C(Nr)  :: ambient Carbon concentration

C     plankton_ambient :: ambient plankton cell concentration

C     I_surf      :: radiation received at surface
C                    in uE/m2-day
C     Inc(Nr)     :: amount of incident radiation
C           note: I_surf from ECCOv4r5, exponential decay with depth

C     W_to_uE_day :: conversion factor from W/m2 to uE/m2-day
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
      _RL A_A(Nr)
      _RL I_surf, Inc(Nr), W_to_uE_day
      _RL Topt, Tmin, Tmax
      _RL mu_opt(Nr), mu(Nr), phi(Nr)

      INTEGER myIter, myThid
CEOP

#ifdef GCHEM_ADD2TR_TENDENCY
C !LOCAL VARIABLES: ====================================================
C  i,j                  :: loop indices
C  i_d0, i_A            :: loop indices related to biofouling
C  iTr                  :: ptracer number
C  jTr                  :: corresponding tracer tendency number
C  bi,bj                :: tile indices
C  k                    :: vertical level
      INTEGER i, j, k, bi, bj, iTr, jTr
      INTEGER i_d0, i_A, iMin, iMax, jMin, jMax

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
C     for more detail see ptracers_integrate, gad_calc_rhs and gad_calc_rising
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

      DO i_d0=1,PTRACERS_n_d_0
C     Only release clean particles
C     i_A=1                   :: all clean without biofouling
C     i_A=2:PTRACERS_n_A      :: biofouling & reacting tracers
       DO i_A=1,2
        jTr = (i_d0-1)*PTRACERS_n_A+i_A
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
      ENDDO
C     Here is the place for code to compute bio-geochemical
C     tendency terms (sometimes referred to as source-minus-sink
C     terms). The tendencies are stored on gchemTendency, as show
C     in the CFC-example.
C     =
C     Feb 2025, by Alan:
C     Biofouling based on (KOOI 2017)
C     https://doi.org/10.1016/j.marpol.2017.01.001

C     Define the conversion factor from W/m2 to uE/m2-day
C     assuming average photon wavelength of 550 nm
C     each photon has  3.61e-19 J  energy  
      W_to_uE_day = 397000.

! C   % Boltzman constant           :: k_bolz
! C   % radius of an plankton cell     :: r_A
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

C     Compute dA/dt between bins allocated at A-cell boundaries
C     so only do i_A = 2 : nA-1
C     i_A = 1 & nA serve as boundaries in the dimension of A
C     so fluxes there are zero, no need to update 
        DO i_d0 = 1, PTRACERS_n_d_0
        DO i_A = 2, PTRACERS_n_A-1
C     Only do for tracers used with biofouling
         jTr = (i_d0-1)*PTRACERS_n_A+i_A
         iTr = jTr + gchem_sepFTr_num

! C         % species specific          :: Vbf
         Vbf = 2.E-16 * (i_A * 1. _d 0 - 1.5) * PTRACERS_delA 
     &      * 4. * PI * (PTRACERS_d_0(i_d0)/2.)**2

         DO j = 1-OLy, sNy+OLy
          DO i = 1-OLx, sNx+OLx

            I_surf = - swflux(i,j,bi,bj) * W_to_uE_day


           DO k = 1, Nr
! C         % profiles

! C         % vertical coordinate       :: z_c == rC(Nr)
! C         % potential temperature     :: theta(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
! C           !!!! actual temperature = potential temperature + tref !!!  
! C         % salinity                  :: salt(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
! C         % water density             :: rhoInSitu(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
! C           !!!! actual rho = rhoInSitu + rhoConst
! C         % water dynamic viscosity   :: use 0.001
! C
! C         % reaction step             :: sA == (jTr * 1. _d 0)
! C         % encounter kernel rate
! C         % diffusivity               :: D_pl(Nr)  
! C         % diffusivity               :: D_A(Nr)
! C         % velocity                  :: w_r(Nr)
! C         % kernel rate               :: b_brown(Nr)
! C         % ambient plankton             :: A_A(Nr) (#cell/m^3) 
            D_pl(k) = k_bolz* (theta(i,j,k,bi,bj)+273.16) 
     &            / (6. * PI * 0.001 * PTRACERS_d_p(iTr) / 2.)
            D_A(k)  = k_bolz* (theta(i,j,k,bi,bj)+273.16) 
     &            / (6. * PI * 0.001 * r_A)
            w_r(k)=abs((PTRACERS_rho_p(iTr)+PTRACERS_rho_p(iTr+1))/2.
     &            - rhoConst - rhoInSitu(i,j,k,bi,bj))
     &            * 9.81 * PTRACERS_d_p(iTr)**2 / 18. / 0.001
            b_brown(k) = 4. * PI * (D_pl(k) + D_A(k)) *
     &            (PTRACERS_d_p(iTr) / 2. + r_A)
     &            + 0.5 * PI * (PTRACERS_d_p(iTr)/2)**2 * w_r(k)
     &            + 1.3 * (1.7E5 / 86400.) *
     &              (PTRACERS_d_p(iTr)/2. + r_A)**3
            A_A(k) = plankton_ambient(i,j,k,bi,bj)


! C         % radiation at surface   :: I_surf in uE/m2-day
! C         % radiation at depth r   :: Inc(Nr)
! C         % light exdinction       :: exponential decay with depth
            Inc(k) = I_surf * exp( 0.2 * rC(k))
            IF (Inc(k) .LT. 0) THEN
                  Inc(k) = 0.
            ENDIF

            
! C         % mu_opt(Nr) :: optimal growth rate in (1/day)
! C         % mu(Nr)     :: growth rate in (1/day)
! C         % phi(Nr)    :: temperature influence
! C         % Irradiance needs conversion into uE/m2-sec 
            mu_opt(k) = 1.85 * Inc(k)/86400 / 
     &            (Inc(k)/86400 + 1.85/0.21)
            Tmin = -3.7
            Topt = 27.5
            Tmax = 34.2
            phi(k) = (theta(i,j,k,bi,bj) - Tmax) 
     &         * (theta(i,j,k,bi,bj) - Tmin)**2 
     &         / (Topt-Tmin) / ( 
     &    (Topt-Tmin)*(theta(i,j,k,bi,bj) - Topt) 
     &  - (Topt-Tmax)*(Topt+Tmin-2.*theta(i,j,k,bi,bj)) 
     &        )
            IF ( (theta(i,j,k,bi,bj) .GT. Tmax) .OR.
     &           (theta(i,j,k,bi,bj) .LT. Tmin) ) THEN
                  phi(k) = 0.
            ENDIF
            mu(k) = mu_opt(k) * phi(k)
            

            bio_dAdt(i,j,k,bi,bj,jTr) = 
     &      b_brown(k)*A_A(k) / (4.*PI*(PTRACERS_d_0(i_d0)/2.)**2)
     &     + mu(k)/ 86400. * (i_A * 1. _d 0 - 1.5) * PTRACERS_delA
     &     - 0.05 / 86400. * (i_A * 1. _d 0 - 1.5) * PTRACERS_delA
     &     - 1.88**(theta(i,j,k,bi,bj)/10. - 2.)*0.06/86400.
     &       * (i_A * 1. _d 0 - 1.5) * PTRACERS_delA

           ENDDO
          ENDDO
         ENDDO
        ENDDO
        ENDDO

        DO i_d0 = 1, PTRACERS_n_d_0
        DO i_A = 2, PTRACERS_n_A-1
! C   % reactive flux due to biofouling
! C   % located at A-boundaries
! C   % with zero-flux boundary condition
! C   % there's only need to compute through 2 to end-1
         jTr = (i_d0-1)*PTRACERS_n_A+i_A
         iTr = jTr + gchem_sepFTr_num
         DO j = 1-OLy, sNy+OLy
          DO i = 1-OLx, sNx+OLx
           DO k = 1, Nr

            IF ((bio_dAdt(i,j,k,bi,bj,jTr)
     &         + bio_dAdt(i,j,k,bi,bj,jTr+1)) .GT. 0.) THEN
             tr_projA = pTracer(i,j,k,bi,bj,iTr)
            ELSE
             tr_projA = pTracer(i,j,k,bi,bj,iTr+1)
            ENDIF

            bio_flux(i,j,k,bi,bj,jTr) = 
     &   bio_dAdt(i,j,k,bi,bj,jTr) * tr_projA

           ENDDO
          ENDDO
         ENDDO
        ENDDO
        ENDDO

        DO i_d0 = 1, PTRACERS_n_d_0
C     i_A = 1 is reserved for particles not subject to biofouling
C     only compute chemical tendencies for i_A > 1
        DO i_A = 2, PTRACERS_n_A
         jTr = (i_d0-1)*PTRACERS_n_A+i_A
         iTr = jTr + gchem_sepFTr_num
         DO j = 1-OLy, sNy+OLy
          DO i = 1-OLx, sNx+OLx
           DO k = 1, Nr
      
               gchemTendency(i,j,k,bi,bj,jTr) =
     &      gchemTendency(i,j,k,bi,bj,jTr) -
     & (bio_flux(i,j,k,bi,bj,jTr) - bio_flux(i,j,k,bi,bj,jTr-1))
     &      / PTRACERS_delA

           ENDDO
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
