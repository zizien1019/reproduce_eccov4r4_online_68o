      INTEGER FUNCTION CAL_BIWEEK_FROM_MYTIME(myTime, myIter, myThid)
      IMPLICIT NONE
#include "CPP_EEMACROS.h"
#include "CPP_OPTIONS.h"
#include "EEPARAMS.h"
#include "cal.h"
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
      INTEGER mydate(4)
      INTEGER yyyy, mm, dd
      INTEGER i, leap, yearday

      EXTERNAL CAL_GETDATE

C     call subroutine cal_getdate
C     and convert  myTime  into  mydate array
      CALL CAL_GETDATE(myIter, myTime, mydate, myThid)

C     extract yyyy, mm, dd from mydate(1) which is in form yyyymmdd
      yyyy = mydate(1) / 10000
      mm   = MOD(mydate(1)/100, 100)
      dd   = MOD(mydate(1), 100)

C     leap year check
      IF (MOD(yyyy,4).EQ.0 .AND. (MOD(yyyy,100).NE.0
     &      .OR. MOD(yyyy,400).EQ.0)) THEN
        leap = 2 
      ELSE
        leap = 1
      ENDIF

C     find which number of day in a year this is
      yearday = dd
      DO i = 1, mm-1
        yearday = yearday + nDayMonth(i, leap)
      ENDDO

C     Find which bi-week this is
C           since the ECCO data is bi-weekly averaged
C           there is little need to precisely retain
C           the variation of declination angle.
C     The following code is equivalent to
C           CAL_WEEK_FROM_MYDATE = CEILING ( yearday/14 )
C           because we are working with integers,
C           and the extra days exceeding 14
C           will be truncated by integer division.
      CAL_BIWEEK_FROM_MYTIME = (yearday + 13) / 14

      RETURN
      END
