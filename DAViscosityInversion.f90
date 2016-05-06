!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
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
!
!/******************************************************************************
! *
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 3.3.2008
! *  Modified Data: 27.9.2012
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!> Subroutine for projecting results in structured 3d mesh to a 2d surface.
!>  This solver assumes that the mesh is structural so that it could have 
!>  been obtained by extrusion in the direction of interest. For the given 
!>  direction the corresponding top and bottom node is computed for every node
!>  and this information is used to perform projection to the top or bottom
!>  plane, or alternatively to the whole body. 
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE MeanValue( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE CoordinateSystems
  USE MeshUtils
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL ::  Transient
  REAL(KIND=dp) :: dt
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
  TYPE(ValueList_t),POINTER :: Params
  TYPE(Solver_t), POINTER :: PSolver
  TYPE(Mesh_t), POINTER :: Mesh
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName, OldVarName, Name, &
      LevelsetName, TargetName
  INTEGER :: i,j,k,l,n,dim,Dofs,dof,itop,ibot,iup,jup,lup,ii,jj,Rounds,nsize,layer, &
      ActiveDirection,elem,TopNodes,NoVar,BotNodes
  INTEGER, POINTER :: MaskPerm(:),TopPointer(:),BotPointer(:),UpPointer(:),DownPointer(:),&
      NodeIndexes(:),TargetPointer(:),BotPerm(:),ThickPerm(:),&
      PermOut(:),PermIn(:),LevelsetPerm(:),TopPerm(:),UnitPerm(:)=>NULL()
  LOGICAL :: GotIt, Found, Visited = .FALSE., Initialized = .FALSE.,&
      Debug, MaskExist, GotVar, GotOldVar, &
      FirstTime = .TRUE.
  REAL(KIND=dp) :: dx,UnitVector(3),ElemVector(3),DotPro,Eps,Length,Level,val,q,depth,height,&
          thickness
  REAL(KIND=dp) :: at0,at1,at2
  REAL(KIND=dp), POINTER :: FieldOut(:), FieldIn(:), Levelset(:), Coord(:), TopField(:), ThickIn(:) 
  TYPE(Variable_t), POINTER :: Var, OldVar, ThickVar
  TYPE(Element_t), POINTER :: Element
  TYPE(Nodes_t) :: Nodes
  TYPE(ValueList_t),POINTER :: BC

  
  SAVE Visited,Nodes,Initialized,UnitVector,Coord,MaskExist,MaskPerm,TopPointer,BotPointer,&
      UpPointer,DownPointer,FieldOut,FieldIn,TopNodes,TopPerm, TopField, BotNodes, BotPerm, &
      nsize, UnitPerm, FirstTime
 
!------------------------------------------------------------------------------
!   Initialize the pointers to top and bottom nodes 
!------------------------------------------------------------------------------

  Debug = .TRUE.
  Params => GetSolverParams()
  Mesh => Solver % Mesh
  PSolver => Solver


  CALL Info('MeanValue','Starting MeanValue', level=9)

  IF( .NOT. Initialized ) THEN

    IF(Debug) CALL Info('MeanValue','start init')
    at0 = CPUTime()

    ! Choose active direction coordinate and set corresponding unit vector
    !---------------------------------------------------------------------
    PSolver => Solver
    CALL DetectExtrudedStructure( Mesh, PSolver, ExtVar = Var, &
        TopNodePointer = TopPointer, BotNodePointer = BotPointer, &
        UpNodePointer = UpPointer, DownNodePointer = DownPointer )
    MaskExist = ASSOCIATED( Var % Perm ) 
    IF( MaskExist ) MaskPerm => Var % Perm
    Coord => Var % Values
    nsize = SIZE( Coord )
    Initialized = .TRUE.

    TopNodes = 0
    ALLOCATE( TopPerm( Mesh % NumberOfNodes ) )
    TopPerm = 0
    DO i=1,Mesh % NumberOfNodes
      IF(TopPointer(i) == i) THEN
        TopNodes = TopNodes + 1
        TopPerm(i) = TopNodes
      END IF
    END DO
    IF( TopNodes > 0 ) THEN
      ALLOCATE( TopField( TopNodes ) ) 
      TopField = 0.0_dp
    END IF

    BotNodes = 0
    ALLOCATE( BotPerm( Mesh % NumberOfNodes ) )
    BotPerm = 0
    DO i=1,Mesh % NumberOfNodes
      IF(BotPointer(i) == i) THEN
        BotNodes = BotNodes + 1
        BotPerm(i) = BotNodes
      END IF
    END DO

  END IF
  at0 = CPUTime()
        

  !------------------------------------------------------------------------------
  ! Go through the variables and compute the desired projections
  !------------------------------------------------------------------------------
  GotVar  = .TRUE.
  GotOldVar = .FALSE.
  NULLIFY(OldVar)
  NoVar = 0
  debug = .FALSE.

  ! If it is our first time through or it is a transient simulation calculate the thickness
  IF(FirstTime .OR. Transient) THEN
    CALL Info('MeanValue','Initializing Thickness', level=3)
    TargetName = 'Thickness'
    Dofs = 1

    ! Create the projected variable if needed
    !-----------------------------------------------
    Var => VariableGet( Mesh % Variables, TRIM(TargetName) )
    IF ( .NOT. ASSOCIATED( Var ) )  THEN      
      IF( .NOT. ASSOCIATED( UnitPerm ) ) THEN
        ALLOCATE( UnitPerm( nsize ) ) 
        DO i=1,nsize
          UnitPerm(i) = i
        END DO
      END IF
     CALL VariableAddVector( Mesh % Variables, Solver % Mesh, PSolver, &
          TargetName, Dofs, Perm = UnitPerm )
      
      Var => VariableGet( Mesh % Variables, TRIM(TargetName) )
      IF( ASSOCIATED( Var ) ) THEN
        CALL Info('MeanValue','Created variable: '//TRIM(TargetName),Level=3)
      ELSE
        CALL Warn('MeanValue','Could not create variable: '//TRIM(TargetName))
      END IF 
    END IF
    IF( Var % Dofs /= Dofs ) THEN
      CALL Fatal('StructureProjectToPlane','Mismatch in the dofs in fields!')
    END IF

    FieldOut => Var % Values
    PermOut => Var % Perm    
    FieldOut = 0.0_dp
    TopField = 0.0_dp

    ! First get the thickness value at each top node 
    DO i=1,nsize
      itop = TopPointer(i)
      dx = 0.5*(Coord(UpPointer(i)) - Coord(DownPointer(i)))
      TopField(TopPerm(itop)) = TopField(TopPerm(itop)) + dx 
    END DO

    ! Now take these values and apply them to each node
    DO i=1,nsize
      IF( ASSOCIATED(PermOut)) THEN
        IF( PermOut(i)>0 ) THEN
          FieldOut(PermOut(i)) = TopField(TopPerm(TopPointer(i)))
        END IF
      ELSE
        FieldOut(i) = TopField(TopPerm(TopPointer(i)))
      END IF
    END DO
    FirstTime = .FALSE.
  END IF

  ! Now loop over the input variables we are depth-averaging
  DO WHILE(.TRUE.)

    NoVar = NoVar + 1    
    IF(Debug) PRINT *,'NoVar',NoVar

    WRITE (Name,'(A,I0)') 'Variable ',NoVar
    VarName = ListGetString( Params, TRIM(Name), GotVar )
    NULLIFY(Var)    
    IF(GotVar) THEN
      Var => VariableGet( Model % Variables, TRIM(VarName) )
      IF ( .NOT. ASSOCIATED( Var ) )  THEN
        CALL Fatal('MeanValue','Variable does not exist: '//TRIM(VarName))
      END IF
      FieldIn => Var % Values
      PermIn => Var % Perm
      Dofs = Var % Dofs
      GotOldVar = .TRUE.
      OldVarName = VarName
    ELSE
      Dofs = 1
    END IF

    ! Either new field or new operator is needed
    !-----------------------------------------------
    IF( .NOT. (GotVar ) ) THEN
      IF( NoVar == 1 ) THEN
        CALL Warn('MeanValue','Not even one field to treat?')
      END IF
      EXIT
    END IF

    ! Create the projected variable if needed
    !-----------------------------------------------
    WRITE (Name,'(A,I0)') 'Target Variable ',NoVar
    Call Info('MeanValue', Name, Level=2)
    TargetName = ListGetString( Params, TRIM(Name), GotIt )
    IF( .NOT. GotIt ) THEN
      WRITE (TargetName,'(A,A)') 'avg'//' '//TRIM(VarName)
    END IF
    Call Info('MeanValue', TargetName, Level=2)

    Var => VariableGet( Mesh % Variables, TRIM(TargetName) )
    IF ( .NOT. ASSOCIATED( Var ) )  THEN      
      IF(.NOT. ASSOCIATED( UnitPerm ) ) THEN
        ALLOCATE( UnitPerm( nsize ) ) 
        DO i=1,nsize
          UnitPerm(i) = i
        END DO
      END IF
     CALL VariableAddVector( Mesh % Variables, Solver % Mesh, PSolver, &
          TargetName, Dofs, Perm = UnitPerm )
      
      Var => VariableGet( Mesh % Variables, TRIM(TargetName) )
      IF( ASSOCIATED( Var ) ) THEN
        CALL Info('MeanValue','Created variable: '//TRIM(TargetName),Level=9)
      ELSE
        CALL Warn('MeanValue','Could not create variable: '//TRIM(TargetName))
      END IF 
    END IF
    IF( Var % Dofs /= Dofs ) THEN
      CALL Fatal('StructureProjectToPlane','Mismatch in the dofs in fields!')
    END IF

    FieldOut => Var % Values
    PermOut => Var % Perm    
    FieldOut = 0.0_dp


    ThickVar => VariableGet( Mesh % Variables, TRIM('Thickness') )
    IF( .NOT. ASSOCIATED( ThickVar )) THEN
            CALL Fatal('MeanValue', 'Cannot find the thickness variable')
    END IF
    ThickIn => ThickVar % Values
    ThickPerm => ThickVar % Perm

    ! Loop over components
    !------------------------------------------------
    DO dof = 1, Dofs
      ! First calculate the integrals
      TopField = 0.0_dp
      DO i=1,nsize
        itop = TopPointer(i)
        ! Note for top and bottom this will automatically reduce the distance to half
        !----------------------------------------------------------------------------
        dx = 0.5*(Coord(UpPointer(i)) - Coord(DownPointer(i)))
        j = i
        IF(ASSOCIATED(PermIn)) j = PermIn(i) 
        j = Dofs*(j-1) + dof
        TopField(TopPerm(itop)) = TopField(TopPerm(itop)) + dx * FieldIn(j)
      END DO

      ! Now put in the depth averaged values
      DO i=1,nsize
        j = i
        IF( ASSOCIATED(ThickPerm) ) j = ThickPerm(i)

        IF( ASSOCIATED(PermOut)) THEN
          IF( PermOut(i)>0 ) THEN
            FieldOut(PermOut(i)) = TopField(TopPerm(TopPointer(i))) / ThickIn(j)
          END IF
        ELSE
          FieldOut(i) = TopField(TopPerm(TopPointer(i))) / ThickIn(j)
        END IF
      END DO


    END DO


  END DO


  at1 = CPUTime()

  WRITE(Message,* ) 'Projection time: ',at1-at0
  CALL Info('MeanValue',Message)
  CALL Info( 'MeanValue','------------------------------------------')

!------------------------------------------------------------------------------
END SUBROUTINE MeanValue
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
SUBROUTINE ScaleViscosity( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE CoordinateSystems
  USE MeshUtils
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL ::  Transient
  REAL(KIND=dp) :: dt
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
  TYPE(ValueList_t),POINTER :: Params
  TYPE(Solver_t), POINTER :: PSolver
  TYPE(Mesh_t), POINTER :: Mesh
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName, OldVarName, Name, &
      LevelsetName, TargetName
  INTEGER :: i,j,k,l,n,dim,Dofs,dof,itop,ibot,iup,jup,lup,ii,jj,Rounds,nsize,layer, &
      ActiveDirection,elem,TopNodes,NoVar,BotNodes
  INTEGER, POINTER :: MaskPerm(:),TopPointer(:),BotPointer(:),UpPointer(:),DownPointer(:),&
      NodeIndexes(:),TargetPointer(:),BotPerm(:),ThickPerm(:),&
      PermOut(:),PermIn(:),LevelsetPerm(:),TopPerm(:),UnitPerm(:)=>NULL()
  LOGICAL :: GotIt, Found, Visited = .FALSE., Initialized = .FALSE.,&
      Debug, MaskExist, GotVar, GotOldVar, &
      FirstTime = .TRUE.
  REAL(KIND=dp) :: dx,UnitVector(3),ElemVector(3),DotPro,Eps,Length,Level,val,q,depth,height,&
          thickness
  REAL(KIND=dp) :: at0,at1,at2
  REAL(KIND=dp), POINTER :: FieldOut(:), FieldIn(:), Levelset(:), Coord(:), TopField(:), ThickIn(:) 
  TYPE(Variable_t), POINTER :: Var, OldVar, ThickVar
  TYPE(Element_t), POINTER :: Element
  TYPE(Nodes_t) :: Nodes
  TYPE(ValueList_t),POINTER :: BC

  
  SAVE Visited,Nodes,Initialized,UnitVector,Coord,MaskExist,MaskPerm,TopPointer,BotPointer,&
      UpPointer,DownPointer,FieldOut,FieldIn,TopNodes,TopPerm, TopField, BotNodes, BotPerm, &
      nsize, UnitPerm, FirstTime
  IF( .NOT. Initialized ) THEN

    IF(Debug) CALL Info('MeanValue','start init')
    at0 = CPUTime()

    ! Choose active direction coordinate and set corresponding unit vector
    !---------------------------------------------------------------------
    PSolver => Solver
    CALL DetectExtrudedStructure( Mesh, PSolver, ExtVar = Var, &
        TopNodePointer = TopPointer, BotNodePointer = BotPointer, &
        UpNodePointer = UpPointer, DownNodePointer = DownPointer )
    MaskExist = ASSOCIATED( Var % Perm ) 
    IF( MaskExist ) MaskPerm => Var % Perm
    Coord => Var % Values
    nsize = SIZE( Coord )
    Initialized = .TRUE.

    TopNodes = 0
    ALLOCATE( TopPerm( Mesh % NumberOfNodes ) )
    TopPerm = 0
    DO i=1,Mesh % NumberOfNodes
      IF(TopPointer(i) == i) THEN
        TopNodes = TopNodes + 1
        TopPerm(i) = TopNodes
      END IF
    END DO
    IF( TopNodes > 0 ) THEN
      ALLOCATE( TopField( TopNodes ) ) 
      TopField = 0.0_dp
    END IF

    BotNodes = 0
    ALLOCATE( BotPerm( Mesh % NumberOfNodes ) )
    BotPerm = 0
    DO i=1,Mesh % NumberOfNodes
      IF(BotPointer(i) == i) THEN
        BotNodes = BotNodes + 1
        BotPerm(i) = BotNodes
      END IF
    END DO

  END IF
        
!------------------------------------------------------------------------------
END SUBROUTINE ScaleViscosity
!------------------------------------------------------------------------------


