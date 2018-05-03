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

FUNCTION BetaMzbZbFriction ( Model, nodenumber, y ) RESULT (BFric)
    USE Types
    IMPLICIT NONE
    TYPE(Model_t) :: Model
    TYPE(Solver_t) :: Solver
    INTEGER :: nodenumber
    Real(KIND=dp) :: BFric
    Real(KIND=dp), dimension (1:3) :: y

    if (y(3) > y(2) + 1.0) then
        bfric = 0.0_dp
    else
        bfric = y(1) ** 2.0
    end if
END FUNCTION BetaMzbZbFriction

FUNCTION IsGroundedSSA ( Model, nodenumber, y ) RESULT (grounded)
    USE Types
    IMPLICIT NONE
    TYPE(Model_t) :: Model
    TYPE(Solver_t) :: Solver
    INTEGER :: nodenumber
    Real(KIND=dp) :: grounded
    Real(KIND=dp), dimension (1:2) :: y

    if (y(2) > y(1) + 1.0) then
        grounded = 0.0_dp
    else
        grounded = 1.0_dp
    end if
END FUNCTION IsGroundedSSA

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


FUNCTION SlidCoef_Contact_lilien ( Model, nodenumber, y) RESULT(Bdrag)

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

  REAL (KIND=dp) ::  y, relChange, relChangeOld, Sliding_lilien, DummyCoef

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
        WRITE(Message, '(A)') 'Contact tolerance found'
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
        CASE ('lilien')
           ! Bdrag = Sliding_lilien(Model, nodenumber, y)
           Bdrag = DummyCoef(Model, nodenumber, y)
        CASE DEFAULT
           Bdrag = DummyCoef(Model, nodenumber, y)
        END SELECT
     ELSE
        ! floating node
        Bdrag = 0.0_dp
     END IF
  ELSE
     ! for other surfaces, typically for lateral surfaces within buttressing experiments
     ! Bdrag = Sliding_weertman(Model, nodenumber, y)
     Bdrag = 0
  END IF
END FUNCTION SlidCoef_Contact_lilien

FUNCTION Sliding_lilien_simple (Model, nodenumber, vel) RESULT(Bdrag)

  USE types
  USE DefUtils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  REAL (KIND=dp) :: y 
  Real(kind=dp), dimension (1:2) :: vel
  INTEGER :: nodenumber
  
  TYPE(ValueList_t), POINTER :: BC
  INTEGER, POINTER :: NormalPerm(:)
  INTEGER :: i, j, n
  REAL (KIND=dp) :: C, m, Bdrag 
  REAL (KIND=dp) :: ut, un, ut0
  REAL (KIND=dp), ALLOCATABLE :: velo(:), AuxReal(:), mMat(:)
  LOGICAL :: GotIt, FirstTime = .TRUE., SSA = .FALSE., UnFoundFatal

  
  CHARACTER(LEN=MAX_NAME_LEN) :: FlowSolverName, USF_name='Sliding_lilien_simple'

  SAVE :: velo, SSA
  SAVE :: FlowSolverName, FirstTime

  CALL Info(USF_name,'Called', Level=17)
   
  IF (FirstTime) THEN
     FirstTime = .FALSE.  
     ALLOCATE(velo(2))
     CALL Info(USF_name,'First time complete', Level=3)
  END IF
  
  !Read the coefficients C and m in the sif file
  BC => GetBC(Model % CurrentElement)
  IF (.NOT.ASSOCIATED(BC))THEN
     CALL Fatal('Sliding_Weertman', 'No BC Found')
  END IF
  CALL Info(USF_name,'Got BC', Level=17)
  
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
  

  ut0 = 1.0e-3
  ! Get the variables to compute ut

  ut = SQRT(vel(1)**2.0 + vel(2)**2 )
  
  ut = MAX(ut,ut0)
  Bdrag = MIN(C * ut**(m-1.0),1.0e20)
END FUNCTION Sliding_Lilien_simple
