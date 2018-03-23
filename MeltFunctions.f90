!
! MeltFunctions.f90
! Copyright (C) 2016 dlilien <dlilien90@gmail.com>
!
! Distributed under terms of the MIT license.
!
      Module MeltFunctions
          implicit None
          public
          contains

        FUNCTION RescaleByTime(t) RESULT(scalef)
            USE Types
            IMPLICIT NONE
            REAL(KIND=dp) :: t, scalef, frac
            INTEGER :: i
            INTEGER :: numyears
            REAL(KIND=dp), DIMENSION(8) :: years, melts
            numyears = 8
            years = (/1996.0, 2006.0, 2007.0, 2008.0, 2009.0, &
                                2010.0, 2011.7, 2014.5 /)
            melts = (/ 1.0, 1.3749287268,  1.4538309038, 1.5116462754 &
         , 1.54864336104, 1.59925859866, 1.56241501754, 1.65287227503 /)
            IF (t <= years(1)) THEN
                scalef = melts(1)
            ELSE IF (t >= years(numyears)) THEN
                scalef = melts(numyears)
            ELSE
                DO i = 2,(numyears)
                    IF (years(i) >= t) THEN
                        frac = (t - years(i - 1)) / (years(i) - years(i - 1))
                        scalef = melts(i - 1) * (1 - frac) + melts(i) * frac
                        EXIT
                    END IF 
                END DO
            END IF
            RETURN
        END

        FUNCTION BasalMeltFavier(z) RESULT(melt)
            USE Types
            IMPLICIT NONE
            REAL(KIND=dp) :: melt
            REAL(KIND=dp) :: z
            IF (z <= -600_dp) THEN
                melt = 200_dp
            ELSE IF (z >= -400) THEN
                melt = 0_dp
            ELSE
                melt = - 0.5_dp * z - 200_dp
            END IF
            RETURN
        END

        FUNCTION BasalMeltJoughin(z) RESULT(melt)
            Use Types
            IMPLICIT NONE
            REAL(KIND=dp) :: melt
            REAL(KIND=dp) :: z
            IF (z <= -600_dp) THEN
                melt = -0.273_dp * z - 154_dp
            ELSE IF (z >= -375_dp) THEN
                melt = 0_dp
            ELSE
                melt = -0.04_dp * z - 15_dp
            END IF
            RETURN
        END

        FUNCTION BasalMeltShean(z) RESULT(melt)
            Use Types
            IMPLICIT NONE
            REAL(KIND=dp) :: melt
            REAL(KIND=dp) :: z
            IF (z <= -400_dp) THEN
                melt = -0.182_dp * (z + 400.0_dp)
            ELSE
                melt = 0.0_dp
            END IF
            RETURN
        END
      END MODULE MeltFunctions
