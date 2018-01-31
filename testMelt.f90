!
! testMelt.f90
! Copyright (C) 2017 dlilien <dlilien@pfe22>
!
! Distributed under terms of the MIT license.
!

      PROGRAM testmf
          USE MeltFunctions
          USE Types
          implicit none
          INTEGER :: numyears, i
          Logical :: failed
          REAL(KIND=dp), DIMENSION(10) :: years, melts
          REAL rst
          numyears = 10
          years = (/1900.0, 1996.0, 2006.0, 2007.0, 2008.0, 2009.0, &
                            2010.0, 2011.7, 2014.5, 3000.0/)
          melts = (/ 1.0, 1.0, 1.37492872682,  1.45383090385, 1.51164627543, & 
                1.54864336104, 1.59925859866, 1.56241501754, 1.65287227503, 1.65287227503/)

          failed = .FALSE.
          do i = 1,numyears
            rst = RescaleByTime(years(i))
            print *, years(i), rst
            if ((abs(melts(i) - rst) / melts(i)) .GE. 1.0e-1) then
                failed = .TRUE.
            end if
          end do
          if (failed) then
            print *, 'Testing melt failed'
          else
            print *, 'Testing melt passed'
          end if
      END PROGRAM testmf
