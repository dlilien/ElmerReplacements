!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! * 
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
! ******************************************************************************
! *
! *  Authors: Olivier Gagliardini, Ga¨el Durand, Thomas Zwinger
! *  Email:   
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 
! *   2007/10/25. Gael Durand
! *   2008/04/06 OG 2D -> 3D
! *   2009/05/18 OG FirstTime in the SAVE !
! *****************************************************************************
!> USF_Sliding.f90
!> 
!> 
!>  Gives the basal drag for different sliding law
!> 
!>  (1) Sliding_Weertman
!>  Need some inputs in the sif file.
!>  Parameters: Weertman Friction Coefficient      -> C 
!>              Weertman Exponent         -> m
!>              Weertman Linear Velocity  -> ut0
!> 
!>  Compute the Bdrag coefficient such as tau_b = Bdrag ub
!>  for the non-linear Weertman law tau_b = C ub^m
!>  To linearize the Weertman law, we can distinguish 4 cases:
!>    1/ ut=0 , tau_b=0     => Bdrag = infinity (no sliding, first step)
!>    2/ ut=0 , tau_b =/0   => Bdrag = C^1/m tau_b^(1-1/m)
!>    3/ ut=/0 , tau_b=0    => Bdrag = Cub^(m-1)
!>    4/ ut=/0 , tau_b=/0   => case 3 
!>  For cases 3 and 4, if ut < ut0, Bdrag = C ut0^{m-1}
!> 
!> 
!>  (2) Friction_Coulomb Sliding Gag JGR 2007
!>  Need some inputs in the sif file.
!>  Parameters: Friction Law Sliding Coefficient      -> As 
!>              Friction Law Post-Peak Exponent         -> q >= 1
!>              Friction Law Maximum Value            -> C ~ max bed slope   
!>              Friction Law Linear Velocity          -> ut0
!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! * 
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
! ******************************************************************************
! *
! *  Authors: 
! *  Email:   
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 
! *  Date Modifications:2007/10/25. Gael Durand
! * 
! *****************************************************************************
!>  tests if resting ice becomes floating ice (GroundedMask from 0 or 1 to -1)
!>  
!>  Return a friction coefficient
!>  -from sliding weertman for resting ice (GroundedMask = 0 or 1)
!>  -of 0 for floating ice (GroundedMask = -1)
!>  2014 : Introduce 3 different ways of defining the grounding line (Mask = 0)
!>  Last Grounded ; First Floating ; Discontinuous

FUNCTION DummyCoef ( Model, nodenumber, y) RESULT(Bdrag)
  USE Types
  USE ElementDescription
  USE DefUtils

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  INTEGER :: nodenumber
  Real(KIND=dp) :: y, bdrag
  bdrag = y ** 2.0
END FUNCTION DummyCoef


FUNCTION SlidCoef_Contact ( Model, nodenumber, y) RESULT(Bdrag)

  USE Types
  USE ElementDescription
  USE DefUtils

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  TYPE(variable_t), POINTER :: TimeVar, NormalVar, VarSurfResidual, GroundedMaskVar, HydroVar, DistanceVar
  TYPE(ValueList_t), POINTER :: BC
  TYPE(Element_t), POINTER :: Element, CurElement, BoundaryElement
  TYPE(Nodes_t), SAVE :: Nodes

  REAL(KIND=dp), POINTER :: NormalValues(:), ResidValues(:), GroundedMask(:), Hydro(:), Distance(:)
  REAL(KIND=dp) :: Bdrag, t, told, thresh
  REAL(KIND=dp), ALLOCATABLE :: Normal(:), Fwater(:), Fbase(:)

  INTEGER, POINTER :: NormalPerm(:), ResidPerm(:), GroundedMaskPerm(:), HydroPerm(:), DistancePerm(:)
  INTEGER :: nodenumber, ii, DIM, GL_retreat, n, tt, Nn, jj, MSum, ZSum

  LOGICAL :: FirstTime = .TRUE., GotIt, Yeschange, GLmoves, Friction,UnFoundFatal=.TRUE.

  REAL (KIND=dp) ::  y, relChange, relChangeOld, Sliding_Budd, Sliding_Weertman, Friction_Coulomb, Sliding_lilien, DummyCoef

  REAL(KIND=dp) :: comp, cond, TestContact
  CHARACTER(LEN=MAX_NAME_LEN) :: USF_Name='SlidCoef_Contact', Sl_Law, GLtype

  SAVE FirstTime, yeschange, told, GLmoves, thresh, GLtype, TestContact
  SAVE DIM, USF_Name, Normal, Fwater, Fbase, relChangeOld, Sl_Law
  CALL Info(USF_name,'Called', Level=3)

!----------------------------------------------------------------------------

! Real time import
  Timevar => VariableGet( Model % Variables,'Time')
  t = TimeVar % Values(1)

! GroundedMask import
  GroundedMaskVar => VariableGet( Model % Mesh % Variables, 'GroundedMask',UnFoundFatal=UnFoundFatal)
  GroundedMask => GroundedMaskVar % Values
  GroundedMaskPerm => GroundedMaskVar % Perm
  
  relchange = Model % Solver % Variable % NonLinChange
  
  CALL Info(USF_name,'Got some basics', Level=3)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! First time step for the First time
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IF (FirstTime) THEN
     DIM = CoordinateSystemDimension()
     FirstTime = .FALSE.
     n = Model % MaxElementNodes
     told = t
     
! means have the possibility to change
     yesChange = .TRUE.
     ALLOCATE( Normal(DIM), Fwater(DIM), Fbase(DIM) )
     
     relChangeOld = relChange
    
     ! choice of the Sliding Law
     BoundaryElement => Model % CurrentElement
     BC => GetBC(BoundaryElement)
     
     Sl_Law = GetString( BC, 'Sliding Law', GotIt )
     IF (.NOT.Gotit) THEN
        CALL FATAL(USF_Name,'No "Sliding law" Name given')
     END IF
     
     GLtype = GetString( BC, 'Grounding Line Definition', GotIt )
     IF (.NOT.Gotit) THEN
        GLtype = 'last grounded'
        CALL Info(USF_Name, 'Grounded Line Defined as the last Grounded point', Level=3)
     ELSE
        WRITE(Message, '(A,A)') 'Grounding Line Defined as ', GLtype
        CALL Info(USF_Name, Message, Level=3)
     END IF
     
     ! Possiblity to fix the grounding line, default is a moving Grounding Line
     GLmoves = GetLogical( BC, 'Grounding line moves', GotIt )
     IF (.NOT.GotIt) THEN
        GLmoves = .TRUE.
     END IF
     IF (GLmoves) THEN
        CALL Info(USF_Name, 'GL may move by default', Level=3)
        CALL Info(USF_Name, 'If you want to fix the Grounding Line, put the keyword "Grounding line moves" to False', Level=3)
     ELSE
        CALL Info(USF_Name, 'GL will be fixed', Level=3)
     END IF
     
     TestContact = GetConstReal( BC, 'Test Contact Tolerance', GotIt )
     IF (.NOT.Gotit) THEN
        TestContact = 1.0e-3     
        CALL Info(USF_Name, 'Contact will be tested for a tolerance of 1.0e-3', Level=3)
     ELSE
        WRITE(Message, '(A,e14.8)') 'Contact tested for a tolerance of ', TestContact
        CALL Info(USF_Name, Message, Level=3)
     END IF
     
     ! Possibility to avoid detachement from nodes that are too far inland from the Grounding line
     ! Uses the DistanceSolver
     ! Default is non possible detachment
     thresh = GetConstReal( BC, 'non detachment inland distance', GotIt )
     IF (.NOT.GotIt) THEN
        thresh = -10000.0_dp
        CALL INFO( USF_Name, 'far inland nodes have the possibility to detach by default', Level=3)
        CALL INFO( USF_Name, 'to avoid detachment (when bedrock is well below sea level),', Level=3)
        CALL INFO( USF_Name, 'use the keyword "non detachment inland distance" to the distance you wish', Level=3)
        CALL INFO( USF_Name, 'This works with the DistanceSolver', Level=3)
     ELSE
        CALL INFO( USF_Name, 'far inland nodes will not detach', level=3)
     END IF
     CALL Info(USF_name,'First time complete', Level=3)
  ENDIF

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! First time step for a New time
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IF ( t > told ) THEN
     told = t
     yesChange = .TRUE.
     relChangeOld = relChange
  END IF
  
  ! to use the non detachment possibility when a grounded node is too far from the grounding line
  ! and positioned on a well below sea level bedrock
  IF (thresh.GT.0.0) THEN
     DistanceVar => VariableGet( Model % Mesh % Variables, 'Distance',UnFoundFatal=UnFoundFatal)
     Distance => DistanceVar % Values
     DistancePerm => DistanceVar % Perm
  END IF
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Look at the convergence of the FlowSolver.
  ! If relative change < TestContact, test if traction occurs. To apply one time
  !
  ! Only to release contact between bed and ice as hydrostatic pressure is higher than normal stress
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Normal = 0.0_dp
  Fwater = 0.0_dp
  Fbase = 0.0_dp

  IF ( (relChange.NE.relChangeOld) .AND. (relchange.GT.0.0_dp) .AND. & 
       &            (relchange.LT.TestContact) .AND. (yesChange) .AND. GLmoves ) THEN
     ! Change the basal condition just once per timestep
     yesChange = .FALSE.

     CALL Info(USF_name,'FLOW SOLVER HAS SLIGHTLY CONVERGED: look for new basal conditions', Level=3)

     VarSurfResidual => VariableGet( Model % Mesh % Variables, 'Flow Solution Loads',UnFoundFatal=UnFoundFatal)
     ResidPerm => VarSurfResidual  % Perm
     ResidValues => VarSurfResidual % Values

     NormalVar => VariableGet(Model % Variables,'Normal Vector',UnFoundFatal=UnFoundFatal)
     NormalPerm => NormalVar % Perm
     NormalValues => NormalVar % Values
     
     !Force exerted by the water, computed for each good boundary nodes (whatever on the bed or floating)
     !From GetHydrostaticLoads
     
     HydroVar => VariableGet( Model % Mesh % Variables, 'Fw',UnFoundFatal=UnFoundFatal)
     Hydro => HydroVar % Values
     HydroPerm => HydroVar % Perm
     
     ! Retreat of the Grounding line if Hydro loads higher than residual values
     GL_retreat = 0

     CurElement => Model % CurrentElement
     DO tt = 1, Model % NumberOfBoundaryElements

        Element => GetBoundaryElement(tt)
        IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
        n = GetElementNOFNodes(Element)

        CALL GetElementNodes(Nodes, Element)

        IF (ANY(GroundedMaskPerm(Element % NodeIndexes(1:n))==0)) CYCLE
        DO ii = 1,n

           Nn = GroundedMaskPerm(Element % NodeIndexes(ii))
           ! the grounded mask is not defined here
           IF (Nn==0) CYCLE
           IF (GroundedMask(Nn) < -0.5_dp) CYCLE
           
           jj = Element % NodeIndexes(ii)
           
           ! comparison between water load and reaction
           
           Normal = NormalValues(DIM*(NormalPerm(jj)-1)+1 : DIM*NormalPerm(jj))
           Fwater = Hydro(DIM*(HydroPerm(jj)-1)+1 : DIM*HydroPerm(jj))
           Fbase = ResidValues((DIM+1)*(ResidPerm(jj)-1)+1 : (DIM+1)*ResidPerm(jj)-1)

           ! comparison between water pressure and bed action
           comp = ABS( SUM( Fwater * Normal ) ) - ABS( SUM( Fbase * Normal ) )
           
           
           IF (comp .GE. 0.0_dp) THEN
              IF (thresh.LE.0.0_dp) THEN
                 GroundedMask(Nn) = -1.0_dp
                 GL_retreat = GL_retreat + 1
                 PRINT *, 'Retreat of the Grounding Line : '
                 PRINT *, Nodes % x(ii), Nodes % y(ii), Nodes % z(ii)
              ELSE
                 IF ( Distance(DistancePerm(Element % NodeIndexes(ii))).LE.thresh ) THEN
                    GroundedMask(Nn) = -1.0_dp
                    GL_retreat = GL_retreat + 1
                    PRINT *, 'Retreat of the Grounding Line : '
                    PRINT *, Nodes % x(ii), Nodes % y(ii), Nodes % z(ii)
                 END IF
              END IF
           END IF
        END DO
        
     END DO
     Model % CurrentElement => CurElement
     
     ! with the previous step
     ! Some 0 (Grounding line) may have been replaced by -1 (floating nodes)
     ! here replacement of some 1 by 0's
     
     IF (GL_retreat.GT.0) THEN
        CurElement => Model % CurrentElement
        DO tt = 1, Model % NumberOfBoundaryElements
           
           Element => GetBoundaryElement(tt)
           IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
           n = GetElementNOFNodes(Element)
           
           CALL GetElementNodes(Nodes, Element)
           MSum = 0
           ZSum = 0
           
           IF (ANY(GroundedMaskPerm(Element % NodeIndexes(1:n))==0)) CYCLE
           DO ii = 1,n
              
              Nn = GroundedMaskPerm(Element % NodeIndexes(ii))
              ! the grounded mask is not defined here
              IF (Nn==0) CYCLE
              MSum = MSum + INT(GroundedMask(Nn))
              IF (GroundedMask(Nn)==0.0_dp) ZSum = ZSum + 1
              
           END DO
           
           IF (MSum+ZSum .LT. n) THEN
              DO ii=1,n
                 Nn = GroundedMaskPerm(Element % NodeIndexes(ii))
                 IF (Nn==0) CYCLE
                 
                 IF (GroundedMask(Nn)==1.0_dp) THEN
                    GroundedMask(Nn)=0.0_dp
                 END IF

              END DO
           END IF
           
        END DO
        Model % CurrentElement => CurElement
     END IF
  END IF
  
  relChangeOld = relChange  
  
  
  IF (GroundedMaskPerm(nodenumber) > 0) THEN
  ! for the bottom surface, where the GroundedMask is defined
     cond = GroundedMask(GroundedMaskPerm(nodenumber))
     
  ! Definition of the Grounding line in term of friction 
  ! If GLtype = Last Grounded -> Bdrag = 0 if Nodal Mask < 0
  ! If GLtype = First Floating -> Bdrag = 0 if Nodal Mask <=0
  ! If GLtype = Discontinuous -> Bdrag = 0 if at one node of the element Mask<=0

     Friction = .FALSE.
     SELECT CASE(GLtype)
     CASE('last grounded')
        IF (cond > -0.5) Friction = .TRUE.  
     CASE('first floating')
        IF (cond > 0.5) Friction = .TRUE. 
     CASE('discontinuous')
        BoundaryElement => Model % CurrentElement
        IF (ALL(GroundedMask(GroundedMaskPerm(BoundaryElement % NodeIndexes))>-0.5)) Friction = .TRUE. 
     CASE DEFAULT
        WRITE(Message, '(A,A)') 'GL type not recognised ', GLtype 
        CALL FATAL( USF_Name, Message)
     END SELECT

     IF (Friction) THEN
        ! grounded node
        SELECT CASE(Sl_law)
        CASE ('weertman')
           Bdrag = Sliding_weertman(Model, nodenumber, y)
        CASE ('lilien')
           ! Bdrag = Sliding_lilien(Model, nodenumber, y)
           Bdrag = DummyCoef(Model, nodenumber, y)
        CASE ('budd')
           Bdrag = Sliding_Budd(Model, nodenumber, y)
        CASE ('coulomb')
           Bdrag = Friction_Coulomb(Model, nodenumber, y)
        CASE DEFAULT
           WRITE(Message, '(A,A)') 'Sliding law not recognised ',Sl_law
           CALL FATAL( USF_Name, Message)
        END SELECT
     ELSE
        ! floating node
        Bdrag = 0.0_dp
     END IF
  ELSE
     ! for other surfaces, typically for lateral surfaces within buttressing experiments
     Bdrag = Sliding_weertman(Model, nodenumber, y)
  END IF
END FUNCTION SlidCoef_Contact

FUNCTION Sliding_lilien (Model, nodenumber, x) RESULT(Bdrag)

  USE types
  USE CoordinateSystems
  USE SolverUtils
  USE ElementDescription
  USE DefUtils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  REAL (KIND=dp) :: y , x              
  INTEGER :: nodenumber
  
  TYPE(ValueList_t), POINTER :: BC
  TYPE(Variable_t), POINTER :: NormalVar, FlowVariable
  REAL(KIND=dp), POINTER :: NormalValues(:), FlowValues(:)
  INTEGER, POINTER :: NormalPerm(:), FlowPerm(:)
  INTEGER :: DIM, i, j, n
  REAL (KIND=dp) :: C, m, Bdrag 
  REAL (KIND=dp) :: ut, un, ut0
  REAL (KIND=dp), ALLOCATABLE :: normal(:), velo(:), AuxReal(:), mMat(:)
  LOGICAL :: GotIt, FirstTime = .TRUE., SSA = .FALSE., UnFoundFatal
  
  CHARACTER(LEN=MAX_NAME_LEN) :: FlowSolverName

  SAVE :: normal, velo, DIM, SSA
  SAVE :: FlowSolverName, FirstTime
   
  IF (FirstTime) THEN
     FirstTime = .FALSE.  
     DIM = CoordinateSystemDimension()
     n = Model % MaxElementNodes
     IF ((DIM == 2).OR.(DIM == 3))  THEN
        ALLOCATE(normal(DIM), velo(DIM))
     ELSE
        CALL FATAL('USF_sliding', 'Bad dimension of the problem')
     END IF
     
     !     BC => GetBC(Model % CurrentElement)  
     FlowSolverName = GetString( Model % Solver % Values , 'Flow Solver Name', GotIt )    
     IF (.NOT.Gotit) FlowSolverName = 'Flow Solution'
     SELECT CASE (FlowSolverName)
     CASE ('ssabasalflow') 
        SSA = .TRUE.
     END SELECT
  END IF
  
  !Read the coefficients C and m in the sif file
  BC => GetBC(Model % CurrentElement)
  IF (.NOT.ASSOCIATED(BC))THEN
     CALL Fatal('Sliding_Weertman', 'No BC Found')
  END IF
  
  n = GetElementNOFNodes()
  ALLOCATE (auxReal(n))
  auxReal(1:n) = GetReal( BC, 'Weertman Friction Coefficient', GotIt )
  IF (.NOT.GotIt) THEN
     CALL FATAL('USF_sliding', 'Need a Friction Coefficient for the Weertman sliding law')
  END IF
  DO i=1,n
     IF (nodenumber == Model % CurrentElement % NodeIndexes( i )) EXIT 
  END DO
  C = auxReal(i)
  DEALLOCATE(auxReal)
  
  ALLOCATE (mMat(n))
  mMat(1:n) = GetReal( BC, 'Weertman Exponent', GotIt )
  IF (.NOT.GotIt) THEN
     CALL FATAL('USF_sliding', 'Need an Exponent for the Weertman sliding law')
  END IF
  DO i=1,n
     IF (nodenumber == Model % CurrentElement % NodeIndexes( i )) EXIT 
  END DO
  m = mMat(i)
  DEALLOCATE(mMat)
  

  ut0 = GetConstReal( BC, 'Weertman Linear Velocity', GotIt )
  IF (.NOT.GotIt) THEN
     CALL FATAL('USF_sliding', 'Need a Linear Velocity for the Weertman sliding law')
  END IF
  ! Get the variables to compute ut
  FlowVariable => VariableGet( Model % Variables, FlowSolverName, UnFoundFatal)
  FlowPerm    => FlowVariable % Perm
  FlowValues  => FlowVariable % Values
  
  ! NS, AIFlow cases   
  IF (.NOT.SSA) THEN 
     ! Get the variable to compute the normal
     NormalVar =>  VariableGet(Model % Variables,'Normal Vector', UnFoundFatal)
     NormalPerm => NormalVar % Perm
     NormalValues => NormalVar % Values
     
     DO i=1, DIM
        normal(i) = -NormalValues(DIM*(NormalPerm(Nodenumber)-1) + i)      
        velo(i) = FlowValues( (DIM+1)*(FlowPerm(Nodenumber)-1) + i )
     END DO
     un = SUM(velo(1:DIM)*normal(1:DIM)) 
     ut = SQRT( SUM( (velo(1:DIM)-un*normal(1:DIM))**2.0 ) )
     ! SSA Flow case      
  ELSE
     DO i=1, DIM-1
        velo(i) = FlowValues( (DIM-1)*(FlowPerm(Nodenumber)-1) + i )
     END DO
     ut = SQRT(SUM( velo(1:DIM-1)**2.0 ))
  END IF
  
  ut = MAX(ut,ut0)
  Bdrag = MIN(C * ut**(m-1.0),1.0e20)
END FUNCTION Sliding_Lilien


FUNCTION Sliding_Weertman (Model, nodenumber, x) RESULT(Bdrag)

  USE types
  USE CoordinateSystems
  USE SolverUtils
  USE ElementDescription
  USE DefUtils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  REAL (KIND=dp) :: y , x              
  INTEGER :: nodenumber
  
  TYPE(ValueList_t), POINTER :: BC
  TYPE(Variable_t), POINTER :: NormalVar, FlowVariable
  REAL(KIND=dp), POINTER :: NormalValues(:), FlowValues(:)
  INTEGER, POINTER :: NormalPerm(:), FlowPerm(:)
  INTEGER :: DIM, i, j, n
  REAL (KIND=dp) :: C, m, Bdrag 
  REAL (KIND=dp) :: ut, un, ut0
  REAL (KIND=dp), ALLOCATABLE :: normal(:), velo(:), AuxReal(:), mMat(:)
  LOGICAL :: GotIt, FirstTime = .TRUE., SSA = .FALSE., UnFoundFatal
  
  CHARACTER(LEN=MAX_NAME_LEN) :: FlowSolverName

  SAVE :: normal, velo, DIM, SSA
  SAVE :: FlowSolverName, FirstTime
  WRITE(Message, '(A)') 'Called'
  CALL Info('USF_Weertman', Message, Level=19)
   
  IF (FirstTime) THEN
     FirstTime = .FALSE.  
     DIM = CoordinateSystemDimension()
     n = Model % MaxElementNodes
     IF ((DIM == 2).OR.(DIM == 3))  THEN
        ALLOCATE(normal(DIM), velo(DIM))
     ELSE
        CALL FATAL('USF_sliding', 'Bad dimension of the problem')
     END IF
     
     !     BC => GetBC(Model % CurrentElement)  
     FlowSolverName = GetString( Model % Solver % Values , 'Flow Solver Name', GotIt )    
     IF (.NOT.Gotit) FlowSolverName = 'Flow Solution'
     SELECT CASE (FlowSolverName)
     CASE ('ssabasalflow') 
        SSA = .TRUE.
     END SELECT

        WRITE(Message, '(A)') 'First Time Complete'
        CALL Info('USF_Weertman', Message, Level=3)
  END IF
  
  !Read the coefficients C and m in the sif file
  BC => GetBC(Model % CurrentElement)
  IF (.NOT.ASSOCIATED(BC))THEN
     CALL Fatal('Sliding_Weertman', 'No BC Found')
  END IF
  
  n = GetElementNOFNodes()
  ALLOCATE (auxReal(n))
  auxReal(1:n) = GetReal( BC, 'Weertman Friction Coefficient', GotIt )
  IF (.NOT.GotIt) THEN
     CALL FATAL('USF_sliding', 'Need a Friction Coefficient for the Weertman sliding law')
  END IF
  DO i=1,n
     IF (nodenumber == Model % CurrentElement % NodeIndexes( i )) EXIT 
  END DO
  C = auxReal(i)
  DEALLOCATE(auxReal)
  
  ALLOCATE (mMat(n))
  mMat(1:n) = GetReal( BC, 'Weertman Exponent', GotIt )
  IF (.NOT.GotIt) THEN
     CALL FATAL('USF_sliding', 'Need an Exponent for the Weertman sliding law')
  END IF
  DO i=1,n
     IF (nodenumber == Model % CurrentElement % NodeIndexes( i )) EXIT 
  END DO
  m = mMat(i)
  DEALLOCATE(mMat)
  

  ut0 = GetConstReal( BC, 'Weertman Linear Velocity', GotIt )
  IF (.NOT.GotIt) THEN
     CALL FATAL('USF_sliding', 'Need a Linear Velocity for the Weertman sliding law')
  END IF
  WRITE(Message, '(A)') 'Got variables'
  CALL Info('USF_Weertman', Message, Level=19)
  ! Get the variables to compute ut
  FlowVariable => VariableGet( Model % Variables, FlowSolverName, UnFoundFatal)
  FlowPerm    => FlowVariable % Perm
  FlowValues  => FlowVariable % Values
  WRITE(Message, '(A)') 'Got flow perm'
  CALL Info('USF_Weertman', Message, Level=19)
  
  ! NS, AIFlow cases   
  IF (.NOT.SSA) THEN 
     ! Get the variable to compute the normal
     NormalVar =>  VariableGet(Model % Variables,'Normal Vector', UnFoundFatal)
     NormalPerm => NormalVar % Perm
     NormalValues => NormalVar % Values
     WRITE(Message, '(A)') 'Got normal vector'
     CALL Info('USF_Weertman', Message, Level=19)
     
     DO i=1, DIM
        normal(i) = -NormalValues(DIM*(NormalPerm(Nodenumber)-1) + i)      
        velo(i) = FlowValues( (DIM+1)*(FlowPerm(Nodenumber)-1) + i )
     END DO
     un = SUM(velo(1:DIM)*normal(1:DIM)) 
     ut = SQRT( SUM( (velo(1:DIM)-un*normal(1:DIM))**2.0 ) )
     ! SSA Flow case      
  ELSE
     DO i=1, DIM-1
        velo(i) = FlowValues( (DIM-1)*(FlowPerm(Nodenumber)-1) + i )
     END DO
     ut = SQRT(SUM( velo(1:DIM-1)**2.0 ))
  END IF
  
  ut = MAX(ut,ut0)
  Bdrag = MIN(C * ut**(m-1.0),1.0e20)
END FUNCTION Sliding_Weertman

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>  (2) Sliding Gag JGR 2007
!> 
!>  Gagliardini, Cohen, Raback and Zwinger, 2007. Finite-Element Modelling of
!>  Subglacial Cavities and Related Friction Law. J. of Geophys. Res.,  Earth
!>  Surface, 112, F02027
!> 
!>  Need some inputs in the sif file.
!>  Parameters: Friction Law Sliding Coefficient      -> As 
!>              Friction Law Post-Peak Exponent         -> q >= 1
!>              Friction Law Maximum Value            -> C ~ max bed slope   
!>              Friction Law Linear Velocity          -> ut0
!>              Friction Law PowerLaw Exponent        -> m = (n Glen's law)
!> 
!>              Water Pressure (BC)    (Compressive - positive)
!> 
!>   tau_b = C.N.[ X . ub^-n / (1 + a.X^q) ]^1/n . ub
!>   with a = (q-1)^(q-1) / q^q and X = ub / (C^n N^n As)
!> 
!>  => Bdrag = C.N.[ X . ub^-n / (1 + a.X^q) ]^1/n 
FUNCTION Friction_Coulomb (Model, nodenumber, y) RESULT(Bdrag)

  USE types
  USE CoordinateSystems
  USE SolverUtils
  USE ElementDescription
  USE DefUtils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  REAL (KIND=dp) :: y , x              
  INTEGER :: nodenumber
  
  TYPE(ValueList_t), POINTER :: BC, Material
  TYPE(Variable_t), POINTER :: TimeVar, NVariable, StressVariable, NormalVar, FlowVariable
  TYPE(Element_t), POINTER ::  BoundaryElement, ParentElement
  REAL(KIND=dp), POINTER :: StressValues(:), NormalValues(:), FlowValues(:)
  REAL(KIND=dp), POINTER :: NValues(:)
  INTEGER, POINTER :: StressPerm(:), NormalPerm(:), FlowPerm(:), NPerm(:)
  INTEGER :: DIM, i, j, Ind(3,3), n, other_body_id
  REAL (KIND=dp) :: C, m, Bdrag, As, Ne, q, Xi, a, Pw 
  REAL (KIND=dp) :: Snt, Snn, ut, un, ut0, t, t0
  LOGICAL :: GotIt, FirstTime = .TRUE., Cauchy, UnFoundFatal
  REAL (KIND=dp), ALLOCATABLE :: Sig(:,:), normal(:), velo(:), Sn(:), AuxReal(:) 
  
  SAVE :: Sig, normal, velo, DIM, Ind, Sn 
  SAVE :: t0, FirstTime
  
  TimeVar => VariableGet( Model % Variables,'Time')
  t = TimeVar % Values(1)
  
  IF (FirstTime) THEN
     FirstTime = .FALSE.  
     t0 = t
     DIM = CoordinateSystemDimension()
     IF ((DIM == 2).OR.(DIM == 3))  THEN
        ALLOCATE(Sig(DIM,DIM),normal(DIM), velo(DIM), Sn(DIM))
     ELSE
        CALL FATAL('Friction_Coulomb', 'Bad dimension of the problem')
     END IF
     Do i=1, 3
        Ind(i,i) = i
     END DO
     Ind(1,2) = 4
     Ind(2,1) = 4
     Ind(2,3) = 5
     Ind(3,2) = 5
     Ind(3,1) = 6
     Ind(1,3) = 6
  END IF
  
  !Read the coefficients As, C, q, and m=1/n in the BC Section  
  BoundaryElement => Model % CurrentElement
  BC => GetBC(BoundaryElement)  
  n = GetElementNOFNodes()
  IF (.NOT.ASSOCIATED(BC))THEN
     CALL Fatal('Friction_Coulomb', 'No BC Found')
  END IF
  
  !  Friction Law Sliding Coefficient      -> As 
  ALLOCATE (auxreal(n))
  auxReal(1:n) = GetReal( BC, 'Friction Law Sliding Coefficient', GotIt )
  IF (.NOT.GotIt) THEN
     CALL FATAL('Friction_Coulomb', 'Need a Friction Law Sliding Coefficient for the Coulomb Friction Law')
  END IF
  DO i=1, n
     IF (NodeNumber== BoundaryElement % NodeIndexes( i )) EXIT 
  END DO
  As = auxReal(i)
  
  !  Friction Law Post-Peak Exponent         -> q >= 1
  auxReal(1:n) = GetReal( BC, 'Friction Law Post-Peak Exponent', GotIt )
  IF (.NOT.GotIt) THEN
     CALL FATAL('Friction_Coulomb', 'Need a Friction Law Post-Peak Exponent &
          &   (>= 1) for the Coulomb Friction Law')
  END IF
  DO i=1, n
     IF (NodeNumber== BoundaryElement % NodeIndexes( i )) EXIT 
  END DO
  q = auxReal(i)
  
  a = (q-1.0)**(q-1.0) / q**q
  
  !  Friction Law Maximum Value            -> C ~ max bed slope   
  auxReal(1:n) = GetReal( BC, 'Friction Law Maximum Value', GotIt )
  IF (.NOT.GotIt) THEN
     CALL FATAL('Friction_Coulomb', 'Need a Friction Law Maximum Value  &
          &   (~ Max Bed Slope) for the Coulomb Friction Law')
  END IF
  DO i=1, n
     IF (NodeNumber== BoundaryElement % NodeIndexes( i )) EXIT 
  END DO
  C = auxReal(i)
  
  
  !  Friction Law Linear Velocity          -> ut0
  ut0 = GetConstReal( BC, 'Friction Law Linear Velocity', GotIt )
  IF (.NOT.GotIt) THEN
     CALL FATAL('Friction_Coulomb', 'Need a Friction Law Linear Velocity for the Coulomb Friction Law ')
  END IF
  !    
  ! friction Law PowerLaw Exponent m
  m = GetConstReal( BC, 'Friction Law PowerLaw Exponent', GotIt )
  IF (.NOT.GotIt) THEN
     CALL FATAL('Friction_Coulomb', 'Need a Friction Law PowerLaw Exponent &
          &      (= n Glen law) for the Coulomb Friction Law')
  END IF
  !
  ! Effective Pressure is either given as a variable 
  ! or computed as N = -Snn - pw 
  ! Get the effective pressure         
  ! If NVariable does not exist, N will be computed as N = -Snn - pw         
  NVariable => VariableGet( Model % Variables, 'Effective Pressure' )
  IF ( ASSOCIATED( NVariable ) ) THEN
     NPerm    => NVariable % Perm
     NValues  => NVariable % Values

  ELSE 
  ! Get the water Pressure from the Stokes keyword 'External Pressure'
  ! Use the convention for the water pressure Pw > 0 => Compression
     auxReal(1:n) = GetReal( BC, 'External Pressure', GotIt )
     IF (.NOT.GotIt) THEN
        CALL FATAL('Friction_Coulomb', 'Need External Pressure &
          &      Or Variable Effective Pressure')
     END IF
     DO i=1, n
        IF (NodeNumber== BoundaryElement % NodeIndexes( i )) EXIT 
     END DO
     ! Because the convention is External Pressure < 0 => Compression,
     ! need to change the sign
     Pw = -auxReal(i)
  
     ! Get the variables to compute tau_b
     StressVariable => VariableGet( Model % Variables, 'Stress', UnFoundFatal)
     StressPerm    => StressVariable % Perm
     StressValues  => StressVariable % Values
     !
     ! Cauchy or deviatoric stresses ?
     !
     other_body_id = BoundaryElement % BoundaryInfo % outbody
     IF (other_body_id < 1) THEN ! only one body in calculation
        ParentElement => BoundaryElement % BoundaryInfo % Right
        IF ( .NOT. ASSOCIATED(ParentElement) ) ParentElement => BoundaryElement % BoundaryInfo % Left
     ELSE ! we are dealing with a body-body boundary and asume that the normal is pointing outwards
        ParentElement => BoundaryElement % BoundaryInfo % Right
        IF (ParentElement % BodyId == other_body_id) ParentElement => BoundaryElement % BoundaryInfo % Left
     END IF
  
     Material => GetMaterial(ParentElement)
     Cauchy = ListGetLogical( Material , 'Cauchy', Gotit )

  END IF
  DEALLOCATE(auxReal)
  
  ! Get the flow variables to compute ut
  FlowVariable => VariableGet( Model % Variables, 'Flow Solution', UnFoundFatal)
  FlowPerm    => FlowVariable % Perm
  FlowValues  => FlowVariable % Values
  
  ! Get the normal variable to compute the normal
  NormalVar =>  VariableGet(Model % Variables,'Normal Vector', UnFoundFatal)
  NormalPerm => NormalVar % Perm
  NormalValues => NormalVar % Values
  
  DO i=1, DIM
     normal(i) = -NormalValues(DIM*(NormalPerm(Nodenumber)-1) + i)      
     velo(i) = FlowValues( (DIM+1)*(FlowPerm(Nodenumber)-1) + i )
  END DO
  
  un = SUM(velo(1:DIM)*normal(1:DIM)) 
  ut = SQRT( SUM( (velo(1:DIM)-un*normal(1:DIM))**2.0 ) )
  
  ! Compute Effective Pressure Ne
  ! Effective pressure N >=0   
  IF ( ASSOCIATED( NVariable ) ) THEN
     Ne = NValues(NPerm(Nodenumber))
  ELSE
     DO i=1, DIM
        DO j= 1, DIM
           Sig(i,j) =  &
              StressValues( 2*DIM *(StressPerm(Nodenumber)-1) + Ind(i,j) )
        END DO
        IF (.NOT.Cauchy) THEN 
           Sig(i,i) = Sig(i,i) - FlowValues((DIM+1)*FlowPerm(Nodenumber))
        END IF
     END DO
     ! Stress vector Sn       
     DO i=1, DIM
        Sn(i) = SUM(Sig(i,1:DIM)*normal(1:DIM)) 
     END DO
     ! Convention is such that Snn should be negative (compressive)
     Snn = SUM( Sn(1:DIM) * normal(1:DIM) ) 
     Ne = -Snn -Pw
  ENDIF

  Bdrag = 0.0_dp
  IF ( Ne>0 ) THEN
     ut = MAX(ut,ut0)
     Xi = ut / (As * (C*Ne)**m ) 
     Xi = MIN(Xi,1.0e20_dp)
  ELSE
     Xi = 1.0e20_dp 
     write(*,*)'!!! Ne <=0, nodenumber',nodenumber, Ne
     Ne = 0.0       
  END IF
  
  Bdrag = C*Ne * ((Xi * ut**(-m)) / ( 1.0 + a * Xi**q))**(1.0/m)
  Bdrag = MIN(Bdrag,1.0e20_dp)
  
  ! Stress may be not known at first time / or first steady iteration  
  IF ((t==t0).AND.(.Not.ASSOCIATED( NVariable )).AND.(Snn.GE.0.0_dp)) Bdrag = 1.0e20
END FUNCTION Friction_Coulomb

! Sliding after Budd et al 1984, Annals of Glaciology 5, page 29-36.
!
! Magnitude of sliding is:
!
! tau_b = C.{u_b}^{m}*Zab^{q}
! 
! where Zab is height above bouyancy and C, m and q respectively are 
! given in the sif by:
!  Budd Friction Coefficient = Real 2.412579e-2        
!  Budd Velocity Exponent = Real $1.0/3.0
!  Budd Zab Exponent = Real 2.0
!
!  Budd Floatation = Logical False
! If this is set to true then the height above bouyancy will be based 
! on the floatation condition instead of inferred from the effective 
! pressure (i.e. depth is used instead of normal stress). Default is 
! false.
!
! Linearisation for small velocity is used, similar to Weertman 
! above:
!  Budd Linear Velocity = Real 0.00001
!
! Slip coefficient in Elmer is given as 
! C.{u_b}^{m - 1}*Zab^{q}
!
! 
! Pre-requisites are as for EffectivePressure below, plus:
! "Budd Ice Density" and "Budd Gravity" need to be defined in the 
! relevant (basal) boundary condition, in addition to the four 
! parameters above.
!
FUNCTION Sliding_Budd (Model, nodenumber, z) RESULT(Bdrag)

  USE types
  USE CoordinateSystems
  USE SolverUtils
  USE ElementDescription
  USE DefUtils

  IMPLICIT NONE

  TYPE(Model_t)  :: Model
  REAL (KIND=dp) :: z
  INTEGER        :: nodenumber

  REAL (KIND=dp) :: Bdrag 
 
  REAL (KIND=dp), ALLOCATABLE :: normal(:), velo(:)
  TYPE(ValueList_t), POINTER  :: BC, ParentMaterial
  TYPE(Variable_t), POINTER   :: NormalVar, FlowVariable, Hvar
  TYPE(Element_t), POINTER    :: parentElement, BoundaryElement
  REAL(KIND=dp), POINTER      :: NormalValues(:), FlowValues(:), HValues(:)
  INTEGER, POINTER :: NormalPerm(:), FlowPerm(:), HPerm(:)
  INTEGER          :: DIM, i, body_id, other_body_id, material_id
  REAL (KIND=dp)   :: C, m, q, g, rhoi, Zab, Zab_offset, ep, sl, H, rhow
  REAL (KIND=dp)   :: ut, un, ut0
  LOGICAL          :: GotIt, FirstTime = .TRUE., SSA = .FALSE., UseFloatation = .FALSE., H_scaling
  LOGICAL          :: UnFoundFatal
  CHARACTER(LEN=MAX_NAME_LEN) :: USF_name, FlowSolverName

  SAVE :: normal, velo, DIM, SSA, FirstTime, FlowSolverName, UseFloatation

  USF_name = "Sliding_Budd"

  IF (FirstTime) THEN
     FirstTime = .FALSE.  
     DIM = CoordinateSystemDimension()
     IF ((DIM == 2).OR.(DIM == 3))  THEN
        ALLOCATE(normal(DIM), velo(DIM))
     ELSE
        CALL FATAL(USF_name, 'Bad dimension of the problem')
     END IF
     
     FlowSolverName = GetString( Model % Solver % Values , 'Flow Solver Name', GotIt )    
     IF (.NOT.Gotit) FlowSolverName = 'Flow Solution'
     SELECT CASE (FlowSolverName)
     CASE ('ssabasalflow') 
        SSA = .TRUE.
     END SELECT
  END IF
  
  BC => GetBC(Model % CurrentElement)
  IF (.NOT.ASSOCIATED(BC))THEN
     CALL Fatal(USF_name, 'No BC Found')
  END IF

  !-----------------------------------------------------------------
  ! get some information upon active boundary element and its parent
  !-----------------------------------------------------------------
  BoundaryElement => Model % CurrentElement
  IF ( .NOT. ASSOCIATED(BoundaryElement) ) THEN
     CALL FATAL(USF_Name,'No boundary element found')
  END IF
  other_body_id = BoundaryElement % BoundaryInfo % outbody
  IF (other_body_id < 1) THEN ! only one body in calculation
     ParentElement => BoundaryElement % BoundaryInfo % Right
     IF ( .NOT. ASSOCIATED(ParentElement) ) ParentElement => BoundaryElement % BoundaryInfo % Left
  ELSE ! we are dealing with a body-body boundary and asume that the normal is pointing outwards
     ParentElement => BoundaryElement % BoundaryInfo % Right
     IF (ParentElement % BodyId == other_body_id) ParentElement => BoundaryElement % BoundaryInfo % Left
  END IF
  ! all the above was just so we can get the material properties of the parent element...
  body_id = ParentElement % BodyId
  material_id = ListGetInteger(Model % Bodies(body_id) % Values, 'Material', UnFoundFatal)
  ParentMaterial => Model % Materials(material_id) % Values
  IF (.NOT. ASSOCIATED(ParentMaterial)) THEN
     WRITE(Message,'(A,I10,A,I10)')&
          'No material values found for body no ', body_id,&
          ' under material id ', material_id
     CALL FATAL(USF_Name,Message)
  END IF

  rhoi = GetConstReal( ParentMaterial, 'Density', GotIt )
  IF (.NOT. GotIt) THEN
     CALL FATAL(USF_Name, 'Material property Density not found.')
  END IF
!  rhoi = GetConstReal( BC, 'Budd Ice Density', GotIt )
!  IF (.NOT.GotIt) THEN
!     CALL FATAL(USF_name, 'Need Ice Density for the Budd sliding law')
!  END IF

  C = GetConstReal( BC, 'Budd Friction Coefficient', GotIt )
  IF (.NOT. GotIt) THEN
     CALL FATAL(USF_name, 'Need a Friction Coefficient for the Budd sliding law')
  END IF
  
  m = GetConstReal( BC, 'Budd Velocity Exponent', GotIt )
  IF (.NOT. GotIt) THEN
     CALL FATAL(USF_name, 'Need a velocity Exponent for the Budd sliding law')
  END IF
  
  q = GetConstReal( BC, 'Budd Zab Exponent', GotIt )
  IF (.NOT. GotIt) THEN
     CALL FATAL(USF_name, 'Need a Zab Exponent for the Budd sliding law')
  END IF
  
  Zab_offset = GetConstReal( BC, 'Budd Zab Offset', GotIt )
  IF (.NOT. GotIt) THEN
     Zab_offset = 0.0_dp
  END IF
  
  ut0 = GetConstReal( BC, 'Budd Linear Velocity', GotIt )
  IF (.NOT. GotIt) THEN
     CALL FATAL(USF_name, 'Need a Linear Velocity for the Budd sliding law')
  END IF
  
  g = GetConstReal( BC, 'Budd Gravity', GotIt )
  IF (.NOT. GotIt) THEN
     CALL FATAL(USF_name, 'Need Gravity for the Budd sliding law')
  END IF
   
  UseFloatation = GetLogical( BC, 'Budd Floatation', GotIt )
  IF (.NOT. GotIt) THEN
     CALL FATAL(USF_name, 'Need Floatation for the Budd sliding law')
  END IF
 
  H_scaling = GetLogical( BC, 'Budd Thickness Scaling', GotIt )
  IF (.NOT. GotIt) THEN
     H_scaling = .FALSE.
  END IF
  
  FlowVariable => VariableGet( Model % Variables, FlowSolverName, UnFoundFatal)
  FlowPerm    => FlowVariable % Perm
  FlowValues  => FlowVariable % Values

  IF (.NOT.SSA) THEN 
     NormalVar =>  VariableGet(Model % Variables,'Normal Vector', UnFoundFatal)
     NormalPerm => NormalVar % Perm
     NormalValues => NormalVar % Values
     
     DO i=1, DIM
        normal(i) = -NormalValues(DIM*(NormalPerm(Nodenumber)-1) + i)      
        velo(i) = FlowValues( (DIM+1)*(FlowPerm(Nodenumber)-1) + i )
     END DO
     un = SUM(velo(1:DIM)*normal(1:DIM)) 
     ut = SQRT( SUM( (velo(1:DIM)-un*normal(1:DIM))**2.0 ) )
  ELSE
     DO i=1, DIM-1
        velo(i) = FlowValues( (DIM-1)*(FlowPerm(Nodenumber)-1) + i )
     END DO
     ut = SQRT(SUM( velo(1:DIM-1)**2.0 ))
  END IF

  ! Zab is height above bouyancy of the upper free surface.  This is 
  ! calculated based on the effective pressure at the bed.  The 
  ! effective pressure at the bed is calculated as the normal stress 
  ! at the lower boundary minus the External Pressure (which is set in 
  ! the boundary condition section of the sif).
  IF (UseFloatation) THEN

     Hvar => VariableGet( Model % Variables, "Depth", UnFoundFatal)
     HPerm    => Hvar % Perm
     HValues  => Hvar % Values
     H = HValues(HPerm(nodenumber))
     
     rhow = GetConstReal( BC, 'Budd Ocean Density', GotIt )
     IF (.NOT.GotIt) THEN
        CALL FATAL(USF_name, 'Need Ocean Density for the Budd sliding law')
     END IF
     
     sl = GetCReal( ParentMaterial, 'Sea level', GotIt )
     IF (.NOT.GotIt) THEN
        CALL FATAL(USF_Name, 'Material property Sea level not found.')
     END IF
     
     ! floatation condition
     ! (H - Zab) * rhoi = (sl - z) * rhow
     ! => Zab = H - (sl-z)*rhow/rhoi
     IF (z .LT. sl) THEN
        Zab = H - (sl-z)*rhow/rhoi
     ELSE
        Zab = H
     END IF
     
     ! this "offset" to the height above bouyancy is intended to provide a non-zero  
     ! basal drag due to contact with the bed, even when effective pressure is zero.
     ! Physically, this can be seen as a compromise between Elmer's "Weertman" 
     ! implementation and Elmer's "Budd" implementation.
     Zab = Zab + Zab_offset

  ELSE
     ep = effectivepressure (Model, nodenumber, z)
     Zab = - ep / (g * rhoi)
  END IF

  ut = MAX(ut,ut0) ! linearize for very low velocities

  IF (H_scaling) THEN
     Zab = Zab / H
  END IF

  Bdrag = C * ut**(m-1.0) * Zab**q
  
CONTAINS
  
  ! Effective Pressure is overburden pressure (or, in our case, normal 
  ! stress is more accurate) minus basal water pressure (or "external 
  ! pressure").
  !
  ! Pre-requisites:
  ! "External Pressure" needs to be defined in the relevant boundary 
  ! condition.
  ! The  "ComputeNormal" and "ComputeDevStressNS" solvers need to be 
  ! active.
  FUNCTION EffectivePressure (Model, nodenumber, y) RESULT(ep)
    
    USE types
    USE CoordinateSystems
    USE SolverUtils
    USE ElementDescription
    USE DefUtils
    IMPLICIT NONE
    
    TYPE(Model_t) :: Model
    REAL (KIND=dp) :: y , x              
    INTEGER :: nodenumber
    
    REAL (KIND=dp) :: ep
    
    TYPE(ValueList_t), POINTER :: BC, Material
    TYPE(Variable_t), POINTER :: TimeVar, StressVariable, NormalVar, FlowVariable
    TYPE(Element_t), POINTER ::  BoundaryElement, ParentElement
    REAL (KIND=dp), POINTER :: StressValues(:), NormalValues(:), FlowValues(:)
    INTEGER, POINTER :: StressPerm(:), NormalPerm(:), FlowPerm(:)
    INTEGER :: DIM, i, j, n, other_body_id, Ind(3,3)
    REAL (KIND=dp) :: Pext 
    REAL (KIND=dp) :: Snn, ut, un, t
    LOGICAL :: GotIt, FirstTime = .TRUE., Cauchy
    REAL (KIND=dp), ALLOCATABLE :: Sig(:,:), normal(:), velo(:), Sn(:), AuxReal(:) 
    CHARACTER(LEN=MAX_NAME_LEN) :: USF_name
    
    SAVE :: Sig, normal, velo, DIM, Ind, Sn, FirstTime
    
    USF_name = "EffectivePressure"

    IF (FirstTime) THEN
       FirstTime = .FALSE.  
       DIM = CoordinateSystemDimension()
       IF ((DIM == 2).OR.(DIM == 3))  THEN
          ALLOCATE(Sig(DIM,DIM),normal(DIM),Sn(DIM))
       ELSE
          CALL FATAL(USF_name, 'Bad dimension of the problem')
       END IF
       DO i=1, 3
          Ind(i,i) = i
       END DO
       Ind(1,2) = 4
       Ind(2,1) = 4
       Ind(2,3) = 5
       Ind(3,2) = 5
       Ind(3,1) = 6
       Ind(1,3) = 6
    END IF
    
    ! Check we have a boundary condition...
    BoundaryElement => Model % CurrentElement
    BC => GetBC(BoundaryElement)  
    IF (.NOT.ASSOCIATED(BC))THEN
       CALL Fatal(USF_name, 'No BC Found')
    END IF
    
    n = GetElementNOFNodes()
   ALLOCATE (auxReal(n))
    
    ! Get the external (probably water) pressure
    ! Use the convention Pext > 0 => Compression
    auxReal(1:n) = GetReal( BC, 'External Pressure', GotIt )
    DO i=1, n
       IF (NodeNumber== BoundaryElement % NodeIndexes( i )) EXIT 
    END DO
    Pext = auxReal(i)
    DEALLOCATE(auxReal)
    
    ! Get the variable to compute the normal
    NormalVar =>  VariableGet(Model % Variables,'Normal Vector', UnFoundFatal)
    NormalPerm => NormalVar % Perm
    NormalValues => NormalVar % Values
    
    ! Get the stress variable
    StressVariable => VariableGet( Model % Variables, 'Stress', UnFoundFatal)
    StressPerm    => StressVariable % Perm
    StressValues  => StressVariable % Values
    
    ! Cauchy or deviatoric stresses ?
    ! First, get parent element
    other_body_id = BoundaryElement % BoundaryInfo % outbody
    IF (other_body_id < 1) THEN ! only one body in calculation
       ParentElement => BoundaryElement % BoundaryInfo % Right
       IF ( .NOT. ASSOCIATED(ParentElement) ) ParentElement => BoundaryElement % BoundaryInfo % Left
    ELSE ! we are dealing with a body-body boundary and asume that the normal is pointing outwards
       ParentElement => BoundaryElement % BoundaryInfo % Right
       IF (ParentElement % BodyId == other_body_id) ParentElement => BoundaryElement % BoundaryInfo % Left
    END IF
    Material => GetMaterial(ParentElement)
    Cauchy = ListGetLogical( Material , 'Cauchy', Gotit )
  
    ! stress tensor
    DO i=1, DIM
       DO j= 1, DIM
          Sig(i,j) =  &
               StressValues( 2*DIM *(StressPerm(Nodenumber)-1) + Ind(i,j) )
       END DO
       IF (.NOT.Cauchy) THEN 
          Sig(i,i) = Sig(i,i) - FlowValues((DIM+1)*FlowPerm(Nodenumber))
       END IF
    END DO
    
    ! normal stress
    DO i=1, DIM
       normal(i) = -NormalValues(DIM*(NormalPerm(Nodenumber)-1) + i)      
    END DO
    DO i=1, DIM
       Sn(i) = SUM(Sig(i,1:DIM)*normal(1:DIM)) 
    END DO
    Snn = SUM( Sn(1:DIM) * normal(1:DIM) ) 
    
    ! effective pressure
    ep = -Snn -Pext
      
  END FUNCTION EffectivePressure
  
END FUNCTION Sliding_Budd
