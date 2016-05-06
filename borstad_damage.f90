!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! USF_Damage.f90
! Last modif : 12 December 2014
! Computes the evolution of damage and computes the Enhancement factor of the Glen's law
! as a function of damage evolution 
!
!
! (1) Enhancement Factor 
! Need some inputs in the sif file.
! Parameters: 
! Glen Exponent
!             
!
! (2) SourceDamage
! Need some inputs in the sif file.
! Parameters: 
! Damage Enhancement Factor
! Damage Parameter sigmath
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! EnhancementFactor
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION EnhancementFactor ( Model, nodenumber, D) RESULT(E)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   TYPE(ValueList_t), POINTER :: Material
   TYPE(Solver_t), TARGET :: Solver
   REAL(KIND=dp) :: D, E, n  
   INTEGER :: nodenumber
   LOGICAL :: FirstTime=.TRUE., GotIt

   SAVE FirstTime, n 

   IF (FirstTime) THEN
   FirstTime = .False.
    
      Material => GetMaterial()
      n = GetConstReal( Material, 'Glen Exponent', GotIt )
      IF (.NOT.GotIt) THEN
         WRITE(Message,'(A)') 'Variable Glen Exponent not found. &
              &Setting to 3.0'
         CALL INFO('Damage EnhancementFactor', Message, level=2)
         n = 3.0_dp
      ELSE
         WRITE(Message,'(A,F10.4)') 'n = ', n 
         CALL INFO('Damage EnhancementFactor', Message, level=2)
      END IF
   END IF

       E = (1.0 - D)**(-n) 
    
  ! write(*,*) D
  ! write(*,*)'E', E
END FUNCTION EnhancementFactor


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SourceDamage 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION SourceDamage (Model, nodenumber, D) RESULT(Source)


   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   USE GeneralUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   REAL (KIND=dp) :: D, Source          
   INTEGER :: nodenumber
 
   TYPE(Solver_t):: Solver 
   TYPE(ValueList_t), POINTER :: Material, Constants
   TYPE(Variable_t), POINTER :: StressVariable, FlowVariable, ChiVariable, PSeaDVariable
   TYPE(Variable_t), POINTER :: StrainRateVariable
   REAL(KIND=dp), POINTER :: StressValues(:), FlowValues(:), ChiValues(:), PSeaDValues(:)
   REAL(KIND=dp), POINTER :: StrainRateValues(:)
   INTEGER, POINTER :: StressPerm(:), FlowPerm(:), ChiPerm(:), PSeaDPerm(:)
   INTEGER, POINTER :: StrainRatePerm(:)

   INTEGER :: Ind(3,3), DIM, i, j, indice(3), infor
   REAL (KIND=dp) :: Sig(3,3), SigDev(3,3), tmp 
   REAL (KIND=dp) :: Eps(3,3)
   REAL (KIND=dp) :: Chi, B, pwater
   Real (KIND=dp) :: EffectiveStress, EffectiveStrainRate, Kappa
   REAL (KIND=DP) :: EI(3),Dumy(1),Work(24)
   REAL (KIND=DP) :: EpsilonZero, TauZero
   LOGICAL :: GotIt, FirstTime = .TRUE., Cauchy
   CHARACTER*20 :: USF_Name='SourceDamage'


   SAVE :: Ind, DIM, B
   SAVE :: FirstTime, Cauchy
   SAVE :: TauZero, EpsilonZero, Kappa
   

   IF (FirstTime) THEN
      FirstTime = .FALSE.  
      DIM = CoordinateSystemDimension()

      DO i=1, 3
         Ind(i,i) = i
      END DO
      Ind(1,2) = 4
      Ind(2,1) = 4
      Ind(2,3) = 5
      Ind(3,2) = 5
      Ind(3,1) = 6
      Ind(1,3) = 6

    Material => GetMaterial()
    !Read the coefficients B and sigmath  

    B = GetConstReal( Material, 'Damage Enhancement Factor', GotIt )
      IF (.NOT.GotIt) THEN
      CALL FATAL('Damage Source', 'Damage Enhancement Factor B not Found')
      ELSE
         WRITE(Message,'(A,F10.4)') 'Damage Enhancement Factor = ', B
         CALL INFO('Damage Source', Message, level=2)
      ENDIF

   ! Cauchy or deviatoric stresses ?
      Cauchy = ListGetLogical( Material , 'Cauchy', Gotit )
      WRITE(Message,'(A,L1)') 'Cauchy stress tensor computed ? ', Cauchy 
         CALL INFO('Damage Source', Message, level=2)

       ! Determination of the stress threshold   
       Constants => GetConstants()
       TauZero = GetConstReal( Constants, 'Max Stress', GotIt) 
        IF (.NOT.GotIt) THEN
          TauZero = 0.130_dp
          CALL INFO('USF_Damage','No "Max Stress" given, set &
          to 0.130', level=2)
        END IF
       Kappa = GetConstReal( Constants, 'Borstad Kappa', GotIt) 
        IF (.NOT.GotIt) THEN
          Kappa = 2.8_dp
          CALL INFO('USF_Damage','No "Borstad Kappa" given, set &
          to 0.130', level=2)
        END IF
   ! For now we just use pick a temperature value to calculate varepsilon_0
   ! Go with T=-10 from C&P, :: 3.5e-25 s^-1 Pa^-3
   ! So in a^-1 MPa^-3 this is
   EpsilonZero = 3.5e-25_dp * 1.0e18_dp * 365.25_dp * 24.0_dp * 60.0_dp * 60.0_dp * TauZero**3
   END IF ! FirstTime


   ! Get the Stress                     
   StressVariable => VariableGet( Model % Variables, 'Stress')
   StressPerm    => StressVariable % Perm
   StressValues  => StressVariable % Values

   ! Get the Stress                     
   StrainRateVariable => VariableGet( Model % Variables, 'StrainRate')
   StrainRatePerm    => StrainRateVariable % Perm
   StrainRateValues  => StrainRateVariable % Values

   ! Get Chi variable (positive where damage increases)
   ChiVariable => VariableGet( Model % Variables, 'Chi')
   ChiPerm    => ChiVariable % Perm
   ChiValues  => ChiVariable % Values
   
   ! Get the Sea Damage Pressure
   PSeaDVariable => VariableGet( Model % Variables, 'PSeaD' )
   IF ( ASSOCIATED( PSeaDVariable ) ) THEN
      PSeaDPerm   => PSeaDVariable % Perm
      PSeaDValues => PSeaDVariable % Values
   ELSE
      CALL WARN('Damage Source', 'PSeaD not associated, basal pressure not taken into account in damage formation' )
      CALL WARN('Damage Source', 'Taking default value PSeaD=0.0')
   END IF

   ! Get the variables to compute the hydrostatic pressure  
   FlowVariable => VariableGet( Model % Variables, 'Flow Solution')
   FlowPerm    => FlowVariable % Perm
   FlowValues  => FlowVariable % Values

   Sig = 0.0
   DO i=1, DIM
      DO j= 1, DIM
         Sig(i,j) =  &
              StressValues( 2*DIM *(StressPerm(Nodenumber)-1) + Ind(i,j) )
      END DO
   END DO
   IF (DIM==2) Sig(3,3) = StressValues( 2*DIM *(StressPerm(Nodenumber)-1) + Ind(3,3))

  ! S = sigma + p
  ! Need Cauchy Stress and Deviatoric Stress 
   IF (.NOT.Cauchy) THEN ! If Deviatoric Stress is computed, then, get the
                         ! Cauchy Stress
       SigDev = Sig
       DO i=1,3  
           Sig(i,i) = SigDev(i,i) - FlowValues((DIM+1)*FlowPerm(Nodenumber))
       END DO
   ELSE ! If the Cauchy Stress is computed, then get the Deviatoric Stress 
       DO i=1,3  
           SigDev(i,i) = Sig(i,i) + FlowValues((DIM+1)*FlowPerm(Nodenumber))
            !write(*,*)Sig(i,1), Sig(i,2), Sig(i,3)
       END DO
   END IF

   Eps = 0.0
   DO i=1, DIM
      DO j= 1, DIM
         Eps(i,j) =  &
              StrainRateValues( 2*DIM *(StrainRatePerm(Nodenumber)-1) + Ind(i,j) )
      END DO
   END DO
   IF (DIM==2) Eps(3,3) = StrainRateValues( 2*DIM *(StrainRatePerm(Nodenumber)-1) + Ind(3,3))

   EffectiveStress = SQRT(0.5_dp * (Sig(1,1)**2.0_dp + Sig(2,2)**2.0_dp + Sig(3,3)**2.0_dp + &
       Sig(3,1)**2.0_dp + Sig(2,1)**2.0_dp + Sig(3,2)**2.0_dp))
  
   EffectiveStrainRate = SQRT(0.5_dp * (Eps(1,1)**2.0_dp + Eps(2,2)**2.0_dp + Eps(3,3)**2.0_dp + &
       Eps(3,1)**2.0_dp + Eps(2,1)**2.0_dp + Eps(3,2)**2.0_dp))
   
   IF (EffectiveStress<TauZero) THEN
       Chi = 0.0_dp
   ELSE
       Chi = 1.0_dp - (EffectiveStrainRate / EpsilonZero)**(-1.0_dp/3.0_dp) * &
               EXP(-(EffectiveStrainRate - EpsilonZero) / (EpsilonZero * &
               (kappa - 1.0_dp)))
   END IF

   ChiValues(ChiPerm(nodenumber)) = Chi
 
   Source = Chi

END FUNCTION SourceDamage   
