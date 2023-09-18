!
! sidewall_drag.f90
! Copyright (C) 2023 dlilien <dlilien@hozideh>
!
! Distributed under terms of the MIT license.
!


       Function BGlenTDum(Model, Nodenumber, Tc)

       USE DefUtils
       USE DefGrid
       
       Implicit None
       Real(kind=dp) :: BGlenTDum
       TYPE(Model_t),POINTER :: Model
       INTEGER :: NodeNumber
       Real(kind=dp), Intent(in) :: Tc                     ! Temperature en d Celsius
       Real(kind=dp), Dimension(7) :: Wn      ! Glen law parameters
       ! W(1)=B0 Glen's fluidity parameter
       ! W(2)=n Glen's law exponent
       ! W(3)=Q1 Energy activation for Tc<Tl
       ! W(4)=Q2 Energy activation for Tc>Tl
       ! W(5)=T0 temperature reference for B0
       ! W(6)=Tl temperature for Q1->Q2
       Real(kind=dp), parameter :: Tzero=273.15
       Real(kind=dp) :: Q, DT
       Logical :: GotIt, UnfoundFatal=.FALSE., FirstTime=.TRUE.
       TYPE(ValueList_t),POINTER :: Material

       SAVE :: FirstTime, Wn

       if ( FirstTime ) THEN
           FirstTime = .FALSE.

       Material => GetMaterial()
       Wn(7) = 8.314
      !------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------

      Wn(1) = ListGetConstReal( Material , 'Fluidity Parameter', GotIt,UnFoundFatal=UnFoundFatal)
      Wn(2) = ListGetConstReal( Material , 'Powerlaw Exponent', GotIt,UnFoundFatal=UnFoundFatal)
         !Previous default value: Wn(2) = 1.0
      WRITE(Message,'(A,F10.4)') 'Powerlaw Exponent = ',   Wn(2)
      CALL INFO('BGlenT', Message, Level = 20)

      Wn(3) = ListGetConstReal( Material, 'Activation Energy 1', GotIt,UnFoundFatal=UnFoundFatal)
         !Previous default value: Wn(3) = 1.0
      WRITE(Message,'(A,F10.4)') 'Activation Energy 1 = ',   Wn(3)
      CALL INFO('BGlenT', Message, Level = 20)

      Wn(4) = ListGetConstReal( Material, 'Activation Energy 2', GotIt,UnFoundFatal=UnFoundFatal)
         !Previous default value: Wn(4) = 1.0
      WRITE(Message,'(A,F10.4)') 'Activation Energy 2 = ',   Wn(4)
      CALL INFO('BGlenT', Message, Level = 20)

      Wn(5) = ListGetConstReal(Material, 'Reference Temperature', GotIt,UnFoundFatal=UnFoundFatal)
         !Previous default value: Wn(5) = -10.0 (Celsius)
      WRITE(Message,'(A,F10.4)') 'Reference Temperature = ',   Wn(5)
      CALL INFO('BGlenT', Message, Level = 20)

      Wn(6) = ListGetConstReal( Material, 'Limit Temperature', GotIt,UnFoundFatal=UnFoundFatal)
         !Previous default value: Wn(6) = -10.0 (Celsius)
      WRITE(Message,'(A,F10.4)') 'Limit Temperature = ',   Wn(6)
      CALL INFO('BGlenT', Message, Level = 20)
  END IF

       Q= Wn(3)
       If (Tc.GT.Wn(6)) Q=Wn(4)
       Q=Q/Wn(7)
       DT=-1./(Tc+Tzero)+1./(Wn(5)+Tzero)

       BGlenTDum=Wn(1) * Exp(Q*DT)

       End


       Function SidewallDrag(Model, Nodenumber, ins)
       USE DefUtils
       Implicit None

       INTERFACE
       Function BGlenTDum(Model, Nodenumber, Tc)
       USE DefUtils
       Real(kind=dp) :: BGlenTdum
       Real(kind=dp), INTENT(IN) :: Tc
       TYPE(Model_t),POINTER :: Model
       INTEGER :: NodeNumber
       END FUNCTION
       END INTERFACE
       Real(kind=dp) :: AVal, n, sidewalldrag, num, denom
       Real(kind=dp), Intent(in), Dimension(3) :: ins
       INTEGER :: NodeNumber
       TYPE(Model_t),POINTER :: Model

       Logical :: GotIt, UnfoundFatal=.FALSE.
       TYPE(ValueList_t),POINTER :: Material
       Material => GetMaterial()

       AVal = BGlenTDum(model, nodenumber, ins(3))
       n = ListGetConstReal( Material , 'Powerlaw Exponent', GotIt,UnFoundFatal=UnFoundFatal)


       num = -ins(1) * (n + 1.0) ** (1.0 / n) * abs(ins(1)) ** (1.0 / n - 1.0)
       denom = ins(2) ** ((n + 1.0) / n) * (2.0 * aval) ** (1.0 / n) 
       SidewallDrag = num / denom
       End
