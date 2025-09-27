#ifdef ALLOW_GCHEM
CBOP
C    !ROUTINE: GCHEM_FIELDS.h
C    !INTERFACE:

C    !DESCRIPTION:
C Contains tracer fields specifically for chemical tracers.
C
C  gchemTendency :: 3DxPTRACER_num field that store the tendencies due
C                   to the bio-geochemical model
C  gchemSource   :: 3DxPTRACER_num field that store the tendencies due
C                   to a fixed input
C  average1      :: 3D field that store the average of the shortwave
C                   radiation at the ocean surface,
C                   * assuming the magnitude is UNITY.


#ifdef GCHEM_ADD2TR_TENDENCY
      _RL gchemTendency(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,
     &                  GCHEM_tendTr_num)
      _RL bio_dAdt(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,
     &                  GCHEM_tendTr_num)
      _RL bio_flux(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,
     &                  GCHEM_tendTr_num)
      _RL gchemSource(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,
     &                  GCHEM_tendTr_num)
      _RL algae_ambient(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)

      COMMON /GCHEM_FIELDS/
     &     gchemTendency,
     &     bio_dAdt,
     &     bio_flux,
     &     gchemSource,
     &     algae_ambient

#endif /* GCHEM_ADD2TR_TENDENCY */
CEOP
#endif /* ALLOW_GCHEM */

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
