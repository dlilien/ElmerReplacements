!
! Cost_Functions.f90
! Copyright (C) 2016 dlilien <dlilien@berens>
!
! Distributed under terms of the MIT license.
!
        FUNCTION AdjCostSquares( Model, nodenumber, vel) RESULT(cost)
        USE types
        implicit none
        TYPE(Model_t) :: Model
        Real(kind=dp) :: cost
        Real(kind=dp), dimension (1:4) :: vel
        INTEGER :: nodenumber
        cost=0.5*((vel(1)-vel(2))**2.0_dp + (vel(3)-vel(4))**2.0_dp)
        Return
        END 

        FUNCTION AdjCostSquares_der_x( Model, nodenumber, vel) RESULT(c)
        USE types
        implicit none
        TYPE(Model_t) :: Model
        Real(kind=dp) :: c
        Real(kind=dp), dimension (1:4) :: vel
        INTEGER :: nodenumber
        c=(vel(1)-vel(2))
        Return
        END

        FUNCTION AdjCostSquares_der_y( Model, nodenumber, vel) RESULT(c)
        USE types
        implicit none
        TYPE(Model_t) :: Model
        Real(kind=dp) :: c
        Real(kind=dp), dimension (1:4) :: vel
        INTEGER :: nodenumber
        c=(vel(3)-vel(4))
        Return
        END

        FUNCTION AdjCostNormedTotal( Model, nodenumber, vel) &
            RESULT(cost)
        USE types
        implicit none
        TYPE(Model_t) :: Model
        Real(kind=dp) :: cost
        Real(kind=dp), dimension (1:4) :: vel
        INTEGER :: nodenumber
        cost=0.5*((vel(1)-vel(2))**2.0_dp + &
            (vel(3)-vel(4))**2.0_dp)/(vel(2)**2.0_dp + vel(4)**2.0_dp)
        Return
        END

        FUNCTION AdjCostNormedTotal_der_x( Model, nodenumber, vel) RESULT(c)
        USE types
        implicit none
        TYPE(Model_t) :: Model
        Real(kind=dp) :: c
        Real(kind=dp), dimension (1:4) :: vel
        INTEGER :: nodenumber
        c=(vel(1)-vel(2))/(vel(2)**2.0_dp + vel(4)**2.0_dp)
        Return
        END

        FUNCTION AdjCostNormedTotal_der_y( Model, nodenumber, vel) RESULT(c)
        USE types
        implicit none
        TYPE(Model_t) :: Model
        Real(kind=dp) :: c
        Real(kind=dp), dimension (1:4) :: vel
        INTEGER :: nodenumber
        c=(vel(3)-vel(4))/(vel(2)**2.0_dp + vel(4)**2.0_dp)
        Return
        END
