#ifndef GCHEM_OPTIONS_H
#define GCHEM_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

#ifdef ALLOW_GCHEM

CBOP
C    !ROUTINE: GCHEM_OPTIONS.h
C    !INTERFACE:

C    !DESCRIPTION:
C options for biogeochemistry package
CEOP

C o Allow separated update of Geo-Chemistry and Advect-Diff
C    (fractional time-stepping type) for some gchem tracers
#define GCHEM_SEPARATE_FORCING

C o Allow single update of some gchem tracers, adding Geo-Chemistry
C    tendency to Advect-Diff tendency
#undef GCHEM_ADD2TR_TENDENCY
#ifdef ALLOW_CFC
# define GCHEM_ADD2TR_TENDENCY
#endif
#ifdef ALLOW_SPOIL
# define GCHEM_ADD2TR_TENDENCY
#endif

C     ========================
C     Here starts the modification:
C     ======================== 
C     By Alan 8, OCT, 2023
C     Looks like its easier to use add2tr_tendency
C     so Ill define both sparate forcing and use add2tr_tendency
#define GCHEM_ADD2TR_TENDENCY

#endif /* ALLOW_GCHEM */
#endif /* GCHEM_OPTIONS_H */
