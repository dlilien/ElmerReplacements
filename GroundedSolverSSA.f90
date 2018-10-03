!
! GroundedSolverSSA.f90
! Copyright (C) 2017 dlilien <dlilien@pfe22>
!
! Distributed under terms of the MIT license.
!

SUBROUTINE GroundedSolverSSAManual( Model,Solver,dt,TransientSimulation )
  USE DefUtils
  USE read_routines
  IMPLICIT NONE

  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation

  !------------------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------------------

  TYPE(Element_t),POINTER :: Element
  TYPE(ValueList_t), POINTER :: Material, SolverParams
  TYPE(Variable_t), POINTER :: PointerToVariable, bedrockVar, bedVar, TimeVar
  TYPE(Nodes_t), SAVE :: Nodes

  LOGICAL :: GotIt, stat,UnFoundFatal=.TRUE., firsttime=.TRUE.

  INTEGER :: nx, ny, TimesInd, Old_times_ind
  INTEGER :: i, mn, n, t, Nn, istat, DIM, MSum, ZSum, bedrockSource, Nnbr
  INTEGER, POINTER :: Permutation(:), bedrockPerm(:), bedPerm(:)

  REAL(KIND=dp), POINTER :: VariableValues(:), BedrockValues(:)
  REAL(KIND=dp) :: z, toler, x, y, time
  REAL(KIND=dp), ALLOCATABLE :: zb(:), zbs(:)
  REAL(KIND=dp) :: times(73)
  REAL,ALLOCATABLE :: dem(:,:),xx(:),yy(:)

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName = 'GroundedSolver', bedrockName
  CHARACTER(LEN=MAX_NAME_LEN) :: bedName, fmt_fname, filename

  INTEGER,PARAMETER :: MATERIAL_DEFAULT = 1, MATERIAL_NAMED = 2, VARIABLE = 3
       
  SAVE DIM, SolverName, zb, toler, xx, yy, nx, ny, dem, firsttime, times, old_times_ind, fmt_fname
  !------------------------------------------------------------------------------

!  NULLIFY(bedrockPerm,bedrockVar)

  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values

  TimeVar => VariableGet(Solver % Mesh % Variables, 'Time')
  Time = TimeVar % Values(1)

  IF (FirstTime) THEN
      FirstTime = .FALSE.
      DO i = 0,72
         times(i + 1) = 0.25_dp * i + 1996.0_dp
      END DO
      fmt_fname = "(A, F7.2, A)"
      old_times_ind = -20.0
  END IF

  CALL INFO(SolverName, 'Computing grounded mask from geometry', level=3)
   DIM = CoordinateSystemDimension()
   mn = Solver % Mesh % MaxElementNodes

   IF (ALLOCATED(zb)) DEALLOCATE(zb)     
      ALLOCATE(zb(mn), STAT=istat )
   IF ( istat /= 0 ) THEN
      CALL FATAL( SolverName, 'Memory allocation error.' )
   END IF

   IF (ALLOCATED(zbs)) DEALLOCATE(zbs)     
      ALLOCATE(zbs(mn), STAT=istat )
   IF ( istat /= 0 ) THEN
      CALL FATAL( SolverName, 'Memory allocation error.' )
   END IF
   CALL INFO( SolverName, 'Memory allocation done.',Level=1 )
  
  SolverParams => GetSolverParams()
  toler = GetConstReal(SolverParams, 'Toler', GotIt)
  IF (.NOT.GotIt) THEN
     CALL FATAL(SolverName, 'No tolerance given for the Grounded Mask.')
  END IF

  bedrockName = GetString(SolverParams, 'Bedrock VN', GotIt)
  IF (GotIt) THEN
     bedrockSource = VARIABLE
     CALL info(SolverName, 'Bedrock Variable name found', level=8)
  ELSE
     bedrockName = GetString(SolverParams, 'Bedrock Material', GotIt)
     IF (GotIt) THEN
        bedrockSource = MATERIAL_NAMED
        CALL info(SolverName, 'Bedrock Material name found', level=8)
     ELSE
        bedrockSource = MATERIAL_DEFAULT     
        CALL info(SolverName, 'No Bedrock Variable or Material; searching for material \"Min Zs Bottom\".', level=8)
     END IF
  END IF

  bedName = GetString(SolverParams, 'Bed VN', GotIt)
  IF (GotIt) THEN
     CALL info(SolverName, 'Bed Variable name found', level=8)
  ELSE
     CALL FATAL(SolverName, 'You need to input a \"Bed Variable\"')
 END IF
     
  Time = Time + 1996.0_dp
  IF (Time <= 2014.0_dp) THEN
      TimesInd = MINLOC(ABS(Times - time), 1)
      IF (old_times_ind .NE. TimesInd) THEN
         write(filename, fmt_fname) "/data/rd04/dal22/smith_inputs/masks/interp_gl/newmask_interp_", Times(TimesInd), ".xyz"
         call get_twod_grid(filename, xx, yy, dem)
         WRITE(Message,'(A)') 'Loaded new gmask'
         CALL INFO('GroundedSolver',Message,Level=3)
         nx = SIZE(xx)
         ny = SIZE(yy)
         old_times_ind = TimesInd
      END IF
  END IF

  !--------------------------------------------------------------
  ! Grounded/floating loop based on height of base above bedrock.
  !--------------------------------------------------------------
  DO t = 1, Solver % NumberOfActiveElements
     Element => GetActiveElement(t)
     n = GetElementNOFNodes()
     
     bedrockVar => VariableGet(Model % Mesh % Variables, bedrockName,UnFoundFatal=UnFoundFatal)
     bedrockPerm => bedrockVar % Perm
     zb(1:n) =  bedrockVar % values(bedrockPerm(Element % NodeIndexes)) + toler
     Material => GetMaterial( Element )


     bedVar => VariableGet(Model % Mesh % Variables, bedName, UnFoundFatal=UnFoundFatal)
     bedPerm => bedVar % Perm
     BedrockValues => bedrockVar % Values
     zbs(1:n) = bedVar % values(bedPerm(Element % NodeIndexes))
     NULLIFY(bedPerm)
     NULLIFY(bedVAR)

     CALL GetElementNodes( Nodes )
     
     DO i = 1, n
        Nn = Permutation(Element % NodeIndexes(i))
        Nnbr = bedrockPerm(Element % NodeIndexes(i))
        IF (Nn==0) CYCLE
        z = zbs(i)

        IF (Time <= 2014.0_dp) THEN
            x = Nodes % x( i )
            y = Nodes % y( i )
            VariableValues(Nn) = LID(dem, xx, yy, nx, ny, x, y, 1.0)
            IF (VariableValues(Nn) > 0.5_dp) THEN
                VariableValues(Nn) = -1.0_dp
                BedrockValues(Nnbr) = -9999.0_dp
            ELSE
                VariableValues(Nn) = 1.0_dp
            END IF
        ELSE
            ! Geometrical condition. Is the node is above the bedrock 
            ! (plus the tolerance)?  Note: zb includes tolerance.
            IF (z > zb(i)) THEN
               VariableValues(Nn) = -1.0_dp
            ELSE
               VariableValues(Nn) = 1.0_dp
            END IF
        END IF
     END DO
  END DO

  CALL info(SolverName, 'Finding Grounding Line', level=8)
  
  !--------------------------------------------------------------
  ! Grounding line loop to label grounded points at grounding Line.
  !--------------------------------------------------------------
  ! Loop over each element:
  !  if the sum of the element masks is lower than the element number 
  !  of nodes minus the number of zeros (i.e. if the element has at 
  !  least one floating node), then each mask equal to 1 is modified 
  !  to 0 (i.e. this node is on the grounding line).  
  DO t = 1, Solver % NumberOfActiveElements     
     Element => GetActiveElement(t)
     n = GetElementNOFNodes()
     CALL GetElementNodes( Nodes )
     MSum = 0
     ZSum = 0
     
     DO i = 1, n
        Nn = Permutation(Element % NodeIndexes(i))
        IF (Nn==0) CYCLE
        MSum = MSum + VariableValues(Nn)
        IF (VariableValues(Nn) == 0.0_dp) ZSum = ZSum + 1.0_dp
     END DO
     
     IF (MSum + ZSum < n) THEN
        DO i = 1, n
           Nn = Permutation(Element % NodeIndexes(i))
           IF (Nn==0) CYCLE
           IF (VariableValues(Nn) == 1.0_dp) THEN
              VariableValues(Nn) = 0.0_dp
           END IF
        END DO
     END IF
  END DO
  
  IF ( ParEnv % PEs>1 ) CALL ParallelSumVector( Solver % Matrix, VariableValues, 1 )
 
  CALL INFO( SolverName , 'Done', level=3)
END SUBROUTINE GroundedSolverSSAManual

SUBROUTINE GroundedSolverSSA( Model,Solver,dt,TransientSimulation )
  USE DefUtils
  IMPLICIT NONE

  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation

  !------------------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------------------

  TYPE(Element_t),POINTER :: Element
  TYPE(ValueList_t), POINTER :: Material, SolverParams
  TYPE(Variable_t), POINTER :: PointerToVariable, bedrockVar, bedVar
  TYPE(Nodes_t), SAVE :: Nodes

  LOGICAL :: GotIt, stat,UnFoundFatal=.TRUE.

  INTEGER :: i, mn, n, t, Nn, istat, DIM, MSum, ZSum, bedrockSource
  INTEGER, POINTER :: Permutation(:), bedrockPerm(:), bedPerm(:)

  REAL(KIND=dp), POINTER :: VariableValues(:)
  REAL(KIND=dp) :: z, toler
  REAL(KIND=dp), ALLOCATABLE :: zb(:), zbs(:)

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName = 'GroundedSolver', bedrockName
  CHARACTER(LEN=MAX_NAME_LEN) :: bedName

  INTEGER,PARAMETER :: MATERIAL_DEFAULT = 1, MATERIAL_NAMED = 2, VARIABLE = 3
       
  SAVE DIM, SolverName, zb, toler
  !------------------------------------------------------------------------------

!  NULLIFY(bedrockPerm,bedrockVar)

  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values

  CALL INFO(SolverName, 'Computing grounded mask from geometry', level=3)

  !--------------------------------------------------------------
  ! Allocate some permanent storage:
  !--------------------------------------------------------------
   DIM = CoordinateSystemDimension()
   mn = Solver % Mesh % MaxElementNodes

   IF (ALLOCATED(zb)) DEALLOCATE(zb)     
      ALLOCATE(zb(mn), STAT=istat )
   IF ( istat /= 0 ) THEN
      CALL FATAL( SolverName, 'Memory allocation error.' )
   END IF

   IF (ALLOCATED(zbs)) DEALLOCATE(zbs)     
      ALLOCATE(zbs(mn), STAT=istat )
   IF ( istat /= 0 ) THEN
      CALL FATAL( SolverName, 'Memory allocation error.' )
   END IF
   CALL INFO( SolverName, 'Memory allocation done.',Level=1 )
  
  SolverParams => GetSolverParams()
  toler = GetConstReal(SolverParams, 'Toler', GotIt)
  IF (.NOT.GotIt) THEN
     CALL FATAL(SolverName, 'No tolerance given for the Grounded Mask.')
  END IF

  bedrockName = GetString(SolverParams, 'Bedrock VN', GotIt)
  IF (GotIt) THEN
     bedrockSource = VARIABLE
     CALL info(SolverName, 'Bedrock Variable name found', level=8)
  ELSE
     bedrockName = GetString(SolverParams, 'Bedrock Material', GotIt)
     IF (GotIt) THEN
        bedrockSource = MATERIAL_NAMED
        CALL info(SolverName, 'Bedrock Material name found', level=8)
     ELSE
        bedrockSource = MATERIAL_DEFAULT     
        CALL info(SolverName, 'No Bedrock Variable or Material; searching for material \"Min Zs Bottom\".', level=8)
     END IF
  END IF

  bedName = GetString(SolverParams, 'Bed VN', GotIt)
  IF (GotIt) THEN
     CALL info(SolverName, 'Bed Variable name found', level=8)
  ELSE
     CALL FATAL(SolverName, 'You need to input a \"Bed Variable\"')
 END IF
     

    
  !--------------------------------------------------------------
  ! Grounded/floating loop based on height of base above bedrock.
  !--------------------------------------------------------------
  DO t = 1, Solver % NumberOfActiveElements
     Element => GetActiveElement(t)
     n = GetElementNOFNodes()
     
     SELECT CASE(bedrockSource)
     CASE (VARIABLE)
        bedrockVar => VariableGet(Model % Mesh % Variables, bedrockName,UnFoundFatal=UnFoundFatal)
        bedrockPerm => bedrockVar % Perm
        zb(1:n) =  bedrockVar % values(bedrockPerm(Element % NodeIndexes)) + toler
        NULLIFY(bedrockPerm)
        NULLIFY(bedrockVar)
     CASE (MATERIAL_NAMED)
        Material => GetMaterial( Element )
        zb(1:n) = ListGetReal( Material,bedrockName, n , & 
             Element % NodeIndexes, GotIt,UnFoundFatal=UnFoundFatal) + toler
     CASE (MATERIAL_DEFAULT)
        Material => GetMaterial( Element )
        zb(1:n) = ListGetReal( Material,'Min Zs Bottom',n , & 
             Element % NodeIndexes, GotIt,UnFoundFatal=UnFoundFatal) + toler
     END SELECT

     bedVar => VariableGet(Model % Mesh % Variables, bedname, UnFoundFatal=UnFoundFatal)
     bedPerm => bedVar % Perm
     zbs(1:n) = bedVar % values(bedPerm(Element % NodeIndexes))
     NULLIFY(bedPerm)

     NULLIFY(bedVAR)

     CALL GetElementNodes( Nodes )
     
     DO i = 1, n
        Nn = Permutation(Element % NodeIndexes(i))
        IF (Nn==0) CYCLE
        z = zbs(i)
        
        ! Geometrical condition. Is the node is above the bedrock 
        ! (plus the tolerance)?  Note: zb includes tolerance.
        IF (z > zb(i)) THEN
           VariableValues(Nn) = -1.0_dp
        ELSE
           VariableValues(Nn) = 1.0_dp
        END IF
     END DO
  END DO

  CALL info(SolverName, 'Finding Grounding Line', level=8)
  
  !--------------------------------------------------------------
  ! Grounding line loop to label grounded points at grounding Line.
  !--------------------------------------------------------------
  ! Loop over each element:
  !  if the sum of the element masks is lower than the element number 
  !  of nodes minus the number of zeros (i.e. if the element has at 
  !  least one floating node), then each mask equal to 1 is modified 
  !  to 0 (i.e. this node is on the grounding line).  
  DO t = 1, Solver % NumberOfActiveElements     
     Element => GetActiveElement(t)
     n = GetElementNOFNodes()
     CALL GetElementNodes( Nodes )
     MSum = 0
     ZSum = 0
     
     DO i = 1, n
        Nn = Permutation(Element % NodeIndexes(i))
        IF (Nn==0) CYCLE
        MSum = MSum + VariableValues(Nn)
        IF (VariableValues(Nn) == 0.0_dp) ZSum = ZSum + 1.0_dp
     END DO
     
     IF (MSum + ZSum < n) THEN
        DO i = 1, n
           Nn = Permutation(Element % NodeIndexes(i))
           IF (Nn==0) CYCLE
           IF (VariableValues(Nn) == 1.0_dp) THEN
              VariableValues(Nn) = 0.0_dp
           END IF
        END DO
     END IF
  END DO
  
  IF ( ParEnv % PEs>1 ) CALL ParallelSumVector( Solver % Matrix, VariableValues, 1 )
 
  CALL INFO( SolverName , 'Done', level=3)
 
END SUBROUTINE GroundedSolverSSA



SUBROUTINE UpdateZb( Model,Solver,dt,TransientSimulation )
  USE DefUtils

  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation

  !------------------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------------------

  TYPE(Element_t),POINTER :: Element
  TYPE(ValueList_t), POINTER :: Material, SolverParams
  TYPE(Variable_t), POINTER :: PointerToVariable, bedrockVar, bedVar
  TYPE(Variable_t), POINTER :: surfVar, groundedVar, thicknessVar
  TYPE(Nodes_t), SAVE :: Nodes

  LOGICAL :: GotIt, stat,UnFoundFatal=.TRUE.

  INTEGER :: i, mn, n, t, Ns, Nb, istat, DIM, MSum, ZSum, bedrockSource
  INTEGER, POINTER :: Permutation(:), bedrockPerm(:), bedPerm(:), surfPerm(:)
  INTEGER, POINTER :: groundedPerm(:), thicknessPerm(:)

  REAL(KIND=dp), POINTER :: VariableValues(:)
  REAL(KIND=dp), POINTER :: ThicknessValues(:)
  REAL(KIND=dp) :: z, toler=0., rhoi=910., rhow=1028., draft
  REAL(KIND=dp), ALLOCATABLE :: zb(:), grounded(:), thickness(:)
  REAL(KIND=dp), POINTER :: zbs(:), zs(:)

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName = 'GroundedSolver', bedrockName
  CHARACTER(LEN=MAX_NAME_LEN) :: bedName, surfName, groundedName, thicknessName

  INTEGER,PARAMETER :: MATERIAL_DEFAULT = 1, MATERIAL_NAMED = 2, VARIABLE = 3
       
  SAVE DIM, SolverName, zb, toler
  !------------------------------------------------------------------------------

!  NULLIFY(bedrockPerm,bedrockVar)

  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values

  CALL INFO(SolverName, 'Updating lower ice surface', level=3)

  !--------------------------------------------------------------
  ! Allocate some permanent storage:
  !--------------------------------------------------------------
   DIM = CoordinateSystemDimension()
   mn = Solver % Mesh % MaxElementNodes

   IF (ALLOCATED(zb)) DEALLOCATE(zb)     
      ALLOCATE(zb(mn), STAT=istat )
   IF ( istat /= 0 ) THEN
      CALL FATAL( SolverName, 'Memory allocation error.' )
   END IF

   IF (ALLOCATED(grounded)) DEALLOCATE(grounded)     
      ALLOCATE(grounded(mn), STAT=istat )
   IF ( istat /= 0 ) THEN
      CALL FATAL( SolverName, 'Memory allocation error.' )
   END IF

   IF (ALLOCATED(thickness)) DEALLOCATE(thickness)     
      ALLOCATE(thickness(mn), STAT=istat )
   IF ( istat /= 0 ) THEN
      CALL FATAL( SolverName, 'Memory allocation error.' )
   END IF

   CALL INFO( SolverName, 'Memory allocation done.',Level=1 )
  
  SolverParams => GetSolverParams()
  toler = GetConstReal(SolverParams, 'Toler', GotIt)
  IF (.NOT.GotIt) THEN
     CALL FATAL(SolverName, 'No tolerance given for the Grounded Mask.')
  END IF

  bedrockName = GetString(SolverParams, 'Bedrock VN', GotIt)
  IF (GotIt) THEN
     bedrockSource = VARIABLE
     CALL info(SolverName, 'Bedrock Variable name found', level=8)
  ELSE
     bedrockName = GetString(SolverParams, 'Bedrock Material', GotIt)
     IF (GotIt) THEN
        bedrockSource = MATERIAL_NAMED
        CALL info(SolverName, 'Bedrock Material name found', level=8)
     ELSE
        bedrockSource = MATERIAL_DEFAULT     
        CALL info(SolverName, 'No Bedrock Variable or Material; searching for material \"Min Zs Bottom\".', level=8)
     END IF
  END IF

  bedName = GetString(SolverParams, 'Bed VN', GotIt)
  IF (GotIt) THEN
     CALL info(SolverName, 'Bed Variable name found', level=8)
  ELSE
     CALL FATAL(SolverName, 'You need to input a \"Bed VN\"')
 END IF

  surfName = GetString(SolverParams, 'Surf VN', GotIt)
  IF (GotIt) THEN
     CALL info(SolverName, 'Surf Variable name found', level=8)
  ELSE
     CALL FATAL(SolverName, 'You need to input a \"Surf VN\"')
  END IF

  groundedName = GetString(SolverParams, 'Grounded VN', GotIt)
  IF (GotIt) THEN
     CALL info(SolverName, 'Grounded Variable name found', level=8)
  ELSE
     CALL FATAL(SolverName, 'You need to input a \"Grounded VN\"')
  END IF

  thicknessName = GetString(SolverParams, 'Thickness VN', GotIt)
  IF (GotIt) THEN
     CALL info(SolverName, 'Thickness Variable name found', level=8)
  ELSE
     CALL FATAL(SolverName, 'You need to input a \"Thickness VN\"')
  END IF
  !--------------------------------------------------------------
  ! Grounded/floating loop based on height of base above bedrock.
  !--------------------------------------------------------------
  DO t = 1, Solver % NumberOfActiveElements
     Element => GetActiveElement(t)
     n = GetElementNOFNodes()
     
     SELECT CASE(bedrockSource)
     CASE (VARIABLE)
        bedrockVar => VariableGet(Model % Mesh % Variables, bedrockName,UnFoundFatal=UnFoundFatal)
        bedrockPerm => bedrockVar % Perm
        zb(1:n) =  bedrockVar % values(bedrockPerm(Element % NodeIndexes)) + toler
        NULLIFY(bedrockPerm)
        NULLIFY(bedrockVar)
     CASE (MATERIAL_NAMED)
        Material => GetMaterial( Element )
        zb(1:n) = ListGetReal( Material,bedrockName, n , & 
             Element % NodeIndexes, GotIt,UnFoundFatal=UnFoundFatal) + toler
     CASE (MATERIAL_DEFAULT)
        Material => GetMaterial( Element )
        zb(1:n) = ListGetReal( Material,'Min Zs Bottom',n , & 
             Element % NodeIndexes, GotIt,UnFoundFatal=UnFoundFatal) + toler
     END SELECT

     bedVar => VariableGet(Model % Mesh % Variables, bedname, UnFoundFatal=UnFoundFatal)
     bedPerm => bedVar % Perm
     zbs => bedVar % values
     NULLIFY(bedVAR)

     surfVar => VariableGet(Model % Mesh % Variables, surfname, UnFoundFatal=UnFoundFatal)
     surfPerm => surfVar % Perm
     zs => surfVar % values
     NULLIFY(surfVAR)

     groundedVar => VariableGet(Model % Mesh % Variables, groundedname, UnFoundFatal=UnFoundFatal)
     groundedPerm => groundedVar % Perm
     grounded(1:n) = groundedVar % values(groundedPerm(Element % NodeIndexes))
     NULLIFY(groundedPerm)
     NULLIFY(groundedVAR)

     thicknessVar => VariableGet(Model % Mesh % Variables, thicknessname, UnFoundFatal=UnFoundFatal)
     thicknessPerm => thicknessVar % Perm
     thickness(1:n) = thicknessVar % values(thicknessPerm(Element % NodeIndexes))
     thicknessValues => thicknessVar % values
     DO i = 1, n
        IF (thickness(i) .LE. 0) THEN
            thicknessValues(thicknessPerm(Element % NodeIndexes(i))) = 0.1_dp
            thickness(i) = 0.1_dp
        END IF
     END DO
     NULLIFY(thicknessPerm)
     NULLIFY(thicknessVAR)
     NULLIFY(thicknessValues)

     CALL GetElementNodes( Nodes )
     
     DO i = 1, n
        Ns = surfPerm(Element % NodeIndexes(i))
        Nb = bedPerm(Element % NodeIndexes(i))

        IF (Ns==0) CYCLE
        IF (Nb==0) CYCLE
        
        IF (grounded(i) < 0.5) THEN
            draft = -thickness(i) * rhoi / rhow

            ! if the draft is larger than the bedrock, need to bump it up
            IF (draft .LE. zb(i)) THEN
                zbs(Nb) = zb(i)
            ELSE
                zbs(Nb) = draft
            END IF
        END IF

        zs(Ns) = zbs(Nb) + thickness(i)
     END DO
     NULLIFY(bedPerm)
     NULLIFY(surfPerm)
     NULLIFY(zbs)
     NULLIFY(zs)
  END DO
  CALL INFO( SolverName , 'Done', level=3)

  IF ( ParEnv % PEs>1 ) CALL ParallelSumVector( Solver % Matrix, VariableValues, 1 )
 
END SUBROUTINE UpdateZb
