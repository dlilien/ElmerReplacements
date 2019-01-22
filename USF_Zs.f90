FUNCTION ZsTopMzsIniLimited ( Model, nodenumber, Zs) RESULT(mu)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Variable_t), POINTER :: ZsSol
   INTEGER, POINTER :: ZsPerm(:)
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   INTEGER :: nodenumber,  NMAX, i, dim
   REAL(KIND=dp) :: mu,   Zs(2)       
   REAL(KIND=dp), ALLOCATABLE :: Zs0(:)       
   LOGICAL :: FirstTime=.True.,UnFoundFatal

   SAVE FirstTime
   SAVE Zs0 

   ZsSol => VariableGet( Model % Variables, 'Zs Top',UnFoundFatal=UnFoundFatal)
   ZsPerm => ZsSol % Perm

   IF (FirstTime) THEN
        FirstTime = .False.
        dim = CoordinateSystemDimension()
        NMAX = COUNT( ZsPerm > 0 ) 
        ALLOCATE(Zs0(NMAX))
        DO i = 1, Model % NumberOfNodes
          IF (ZsPerm(i)==0) CYCLE
          IF (dim==2) THEN
             Zs0(ZsPerm(i)) = Model % Nodes % y (i)
          ELSE
             Zs0(ZsPerm(i)) = Model % Nodes % z (i)
          END IF
        END DO
   END IF

      mu =  Zs(1) -  Zs0(ZsPerm(nodenumber)) 
      IF (mu > 100.0) THEN
          mu = 100.0
      END IF
      
      IF (Zs0(ZsPerm(nodenumber)) + mu < Zs(2) + 1.0_dp) THEN
          mu = Zs(2) + 1.0
      END IF

END FUNCTION ZsTopMZsIniLimited

FUNCTION ZsBotMzsIniLimited ( Model, nodenumber, Zs) RESULT(mu)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Variable_t), POINTER :: ZsSol
   INTEGER, POINTER :: ZsPerm(:)
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   INTEGER :: nodenumber,  NMAX, i, dim
   REAL(KIND=dp) :: mu,   Zs(2)       
   REAL(KIND=dp), ALLOCATABLE :: Zs0(:)       
   LOGICAL :: FirstTime=.True.,UnFoundFatal

   SAVE FirstTime
   SAVE Zs0 

   ZsSol => VariableGet( Model % Variables, 'Zs Bottom',UnFoundFatal=UnFoundFatal)
   ZsPerm => ZsSol % Perm

   IF (FirstTime) THEN
        FirstTime = .False.
        dim = CoordinateSystemDimension()
        NMAX = COUNT( ZsPerm > 0 ) 
        ALLOCATE(Zs0(NMAX))
        DO i = 1, Model % NumberOfNodes
          IF (ZsPerm(i)==0) CYCLE
          IF (dim==2) THEN
             Zs0(ZsPerm(i)) = Model % Nodes % y (i)
          ELSE
             Zs0(ZsPerm(i)) = Model % Nodes % z (i)
          END IF
        END DO
   END IF

   mu =  Zs(1) -  Zs0(ZsPerm(nodenumber)) 
   IF (Zs0(ZsPerm(nodenumber)) + mu > Zs(2) - 1.0_dp) THEN
       mu = Zs(2) - 1.0_dp
   END IF

   IF ((Zs(1) > 0.0_dp).AND.(Zs(2) <= 0.0_dp)) THEN
       mu = 0.0_dp
   END IF

END FUNCTION ZsBotMZsIniLimited
