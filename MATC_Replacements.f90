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

        FUNCTION BetaSquare( Model, nodenumber, beta ) Result(Bs)
        USE types
        implicit none
        TYPE(Model_t) :: Model
        Real(kind=dp) :: Bs
        Real(kind=dp) :: beta
        Real(kind=dp) :: yearinsec=365.25_dp*24.0_dp*60.0_dp*60.0_dp
        INTEGER :: nodenumber
        Bs=beta**2.0_dp/(1.0e6_dp*yearinsec)
        Return
        END

        FUNCTION MuSquare( Model, nodenumber, mu ) Result(Ms)
        USE types
        implicit none
        TYPE(Model_t) :: Model
        Real(kind=dp) :: Ms
        Real(kind=dp) :: mu
        Real(kind=dp) :: yearinsec=365.25_dp*24.0_dp*60.0_dp*60.0_dp
        INTEGER :: nodenumber
        Ms=mu**2.0_dp*1.0e-6_dp*(2.0_dp*yearinsec)**(1.0_dp/3.0_dp)
        Return
        END

        FUNCTION Mu( Model, nodenumber, T ) Result(m)
        USE types
        implicit none
        TYPE(Model_t) :: Model
        Integer :: nodenumber
        Real(kind=dp) :: T, m, glen
        Real(kind=dp) :: yearinsec=365.25_dp*24.0_dp*60.0_dp*60.0_dp
        m=sqrt(glen(T) * yearinsec ** (-1.0_dp/3.0_dp) * 1.0e-6_dp)
        return
        END

        FUNCTION Visc( Model, nodenumber, T ) Result(m)
        USE types
        implicit none
        TYPE(Model_t) :: Model
        Integer :: nodenumber
        Real(kind=dp) :: T, m, yearinsec=365.25_dp*24.0_dp*60.0_dp*60.0_dp
        Real(kind=dp) :: AF
        Real(kind=dp) :: UGC=8.314_dp, TT=273.15_dp
        IF (T<-10.0_dp) THEN
                AF = 3.985e-13_dp * exp( -60.0e3_dp/(UGC * (TT + T)))
        ELSE IF (T>0.0_dp) THEN
                AF = 1.916e03_dp * exp(-139.0e3_dp/(UGC * TT))
        ELSE
                AF = 1.916e03_dp * exp(-139.0e3_dp/(UGC *(TT + T)))
        END IF
        m=(2.0_dp * AF * yearinsec)** (-1.0_dp/3.0_dp) * 1.0e-6_dp
        return
        END

        FUNCTION glen( Model, nodenumber, T ) Result(g)
        USE types
        implicit none
        TYPE(Model_t) :: Model
        Real(kind=dp) :: T
        Real(kind=dp) :: g, ArrheniusFactor
        Real(kind=dp) :: yearinsec=365.25_dp*24.0_dp*60.0_dp*60.0_dp
        INTEGER :: nodenumber
        g=(2.0_dp*ArrheniusFactor(T))**(-1.0_dp/3.0_dp)
        Return
        END

        FUNCTION ArrheniusFactor(T) Result(AF)
        USE types
        implicit none
        Real(kind=dp) :: T
        Real(kind=dp) :: AF
        Real(kind=dp) :: UGC=8.314_dp, TT=273.15_dp
        IF (T<-10.0_dp) THEN
                AF = 3.985e-13_dp * exp( -60.0e3_dp/(UGC * (TT + T)))
        ELSE IF (T>0.0_dp) THEN
                AF = 1.916e03_dp * exp(-139.0e3_dp/(UGC * TT))
        ELSE
                AF = 1.916e03_dp * exp(-139.0e3_dp/(UGC *(TT + T)))
        END IF
        return
        END

        FUNCTION HeatCapacity(T) Result(HC)
        Use types
        implicit none
        Real(kind=dp) :: T
        Real(kind=dp) :: HC, yearinsec=365.25_dp*24.0_dp*60_dp*60_dp
        HC=(146.3_dp+7.253_dp * (T + 273.15_dp)) * yearinsec**2.0_dp
        return
        END

        FUNCTION ThermalConductivity(T) Result(TC)
        Use types
        implicit none
        ! I have played around with this function to see if the values
        ! here screw up the peclet number; there is weird behavior with
        ! low conductivity, some kind of numerical ringing i think, but
        ! everything looks good in the conductive limit
        Real(kind=dp) :: T
        Real(kind=dp) :: TC, yearinsec=365.25_dp*24.0_dp*60_dp*60_dp
        TC=9.828_dp * exp(-5.7e-3_dp * (T + 273.15_dp)) * yearinsec * 1.0e-6_dp
        return
        END

        FUNCTION Longitude( Model, nodenumber, coord) Result(Lon)
        Use types
        implicit none
        Real(kind=dp), dimension(1:2) :: coord
        Real(kind=dp) :: Lon
        TYPE(Model_t) :: Model
        Integer :: nodenumber
        Lon = -atan(coord(1)/coord(2))*180.0_dp/3.14159265359_dp
        return
        END

        FUNCTION Latitude( Model, nodenumber, coord) Result(Lat)
        Use types
        implicit none
        Real(kind=dp), dimension(1:2) :: coord
        Real(kind=dp) :: Lat
        TYPE(Model_t) :: Model
        Integer :: nodenumber
        Lat = (-3.14159265359_dp/2.0_dp + 2.0_dp * &
            atan(sqrt(coord(1)**2.0_dp+coord(2)**2.0_dp)/ &
            (2.0_dp * 6371225.0_dp * 0.97276_dp)))*180.0_dp/3.14159265359_dp
        return
        END

        FUNCTION SurfaceTemp( Model, nodenumber, LatX ) Result(T)
        Use types
        implicit none
        Real(kind=dp), dimension(1:2) :: LatX
        Real(kind=dp) :: T
        TYPE(Model_t) :: Model
        Integer :: nodenumber
        T=34.36_dp - 0.68775_dp * abs(LatX(1)) - 9.14e-3_dp * LatX(2)
        return
        END

        FUNCTION PressureMeltingPoint( Model, nodenumber, PIN ) &
            Result(PMP)
        Use types
        implicit none
        Real(kind=dp) :: PIN, PMP, P
        TYPE(Model_t) :: Model
        Integer :: nodenumber
        IF (PIN<0.0_dp) THEN
                P=0.0_dp
        ELSE
                P=PIN
        END IF
        PMP=-(9.8e-2_dp*P)
        return
        END

        FUNCTION FluxUnits( Model, nodenumber, fin ) RESULT(fout)
        Use Types
        implicit none
        Real(kind=dp) :: fin, fout
        TYPE(Model_t) :: Model
        Integer :: nodenumber
        fout=fin*365.25_dp*24.0_dp*60.0_dp*60.0_dp/1.0e6_dp
        return
        END

        FUNCTION BasalMeltFavier( model, nodenumber, z) & 
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

        FUNCTION BasalMeltJoughin( model, nodenumber, z) RESULT(melt)
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

        FUNCTION BasalMeltZero( model, nodenumber, z) RESULT(melt)
        Use Types
        implicit none
        Real(kind=dp) :: melt
        Real(kind=dp) :: z(2)
        Type(Model_t) :: model
        Integer :: nodenumber
        melt = 0_dp
        return
        END

        FUNCTION BasalMeltFavier_negative( model, nodenumber, z) & 
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

        FUNCTION BasalMeltJoughin_negative( model, nodenumber, z) RESULT(melt)
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

        FUNCTION BasalMeltZero_negative( model, nodenumber, z) &
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
