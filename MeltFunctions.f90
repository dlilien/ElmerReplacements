        FUNCTION BMF_fsbad( model, nodenumber, z) & 
            RESULT(melt)
            Use Types
            implicit none
            Real(kind=dp) :: melt
            Real(kind=dp) :: z(2)
            Type(Model_t) :: model
            Integer :: nodenumber
            IF (Z(2) > 0.5_dp) THEN
                melt = 0_dp
            ELSE
                IF (z(1)<=-600_dp) THEN
                    melt = 200_dp
                ELSE IF (z(1)>=400) THEN
                    melt = 0_dp
                ELSE
                    melt = - 0.5_dp * z(1) - 200_dp
                END IF
            END IF
            return
        END

        FUNCTION BMJ_fsbad( model, nodenumber, z) RESULT(melt)
            Use Types
            implicit none
            Real(kind=dp) :: melt
            Real(kind=dp) :: z(2)
            Type(Model_t) :: model
            Integer :: nodenumber
            IF (Z(2) > 0.5_dp) THEN
                melt = 0_dp
            ELSE
                IF (z(1)<=-600_dp) THEN
                    melt = -0.273_dp * z(1) - 154_dp
                ELSE IF (z(1)>=375_dp) THEN
                    melt = 0_dp
                ELSE
                    melt = -0.04_dp * z(1) - 15_dp
                END IF
            END IF
            return
        END

        FUNCTION BMS_fsbad( model, nodenumber, z) RESULT(melt)
            Use Types
            implicit none
            Real(kind=dp) :: melt
            Real(kind=dp) :: z(2)
            Type(Model_t) :: model
            Integer :: nodenumber
            IF (Z(2) > 0.5_dp) THEN
                melt = 0_dp
            ELSE
                IF (z(1)<=-400_dp) THEN
                    melt = -0.182_dp * (z(1) + 400.0_dp)
                ELSE
                    melt = 0.0_dp
                END IF
            END IF
            return
        END

        FUNCTION BMZ_fsbad( model, nodenumber, z) RESULT(melt)
            Use Types
            implicit none
            Real(kind=dp) :: melt
            Real(kind=dp) :: z(2)
            Type(Model_t) :: model
            Integer :: nodenumber
            melt = 0_dp
            return
        END

        FUNCTION BMF_negbad( model, nodenumber, z) & 
            RESULT(melt)
            Use Types
            implicit none
            Real(kind=dp) :: melt
            Real(kind=dp) :: z(2)
            Type(Model_t) :: model
            Integer :: nodenumber
            IF (Z(2) > 0.5_dp) THEN
                melt = 0_dp
            ELSE
                IF (z(1)<=-600_dp) THEN
                    melt = -200_dp
                ELSE IF (z(1)>=400) THEN
                    melt = 0_dp
                ELSE
                    melt = 0.5 * z(1) + 200_dp
                END IF
            END IF
            return
        END

        FUNCTION BMJ_negbad( model, nodenumber, z) &
            RESULT(melt)
            Use Types
            implicit none
            Real(kind=dp) :: melt
            Real(kind=dp) :: z(2)
            Type(Model_t) :: model
            Integer :: nodenumber
            IF (Z(2) > 0.5_dp) THEN
                melt = 0_dp
            ELSE
                IF (z(1)<=-600_dp) THEN
                    melt = 0.273_dp * z(1) + 154_dp
                ELSE IF (z(1)>=375_dp) THEN
                    melt = 0_dp
                ELSE
                    melt = 0.04_dp * z(1) + 15_dp
                END IF
            END IF
            return
        END

        FUNCTION BMS_negbad( model, nodenumber, z) RESULT(melt)
            Use Types
            implicit none
            Real(kind=dp) :: melt
            Real(kind=dp) :: z(2)
            Type(Model_t) :: model
            Integer :: nodenumber
            IF (Z(2) > 0.5_dp) THEN
                melt = 0_dp
            ELSE
                IF (z(1)<=-400_dp) THEN
                    melt = 0.182_dp * (z(1) + 400.0_dp)
                ELSE
                    melt = 0.0_dp
                END IF
            END IF
            return
        END

        FUNCTION BMZ_negbad( model, nodenumber, z) &
            RESULT(melt)
            Use Types
            implicit none
            Real(kind=dp) :: melt
            Real(kind=dp) :: z(2)
            Type(Model_t) :: model
            Integer :: nodenumber
            melt = 0_dp
            return
        END
