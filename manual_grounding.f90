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
! *  Authors: Olivier Gagliardini, Gael Durand
! *  Email:   
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 
! * 
! *****************************************************************************
!> Solver for creating a mask on whether the lower side of an ice sheet/shelf is
!>  grounded or not. +1=grounded,-1=detached, 0=grounding line (=last grounded node)
SUBROUTINE GroundedSolver( Model,Solver,dt,TransientSimulation )
  USE DefUtils
  USE read_routines

  IMPLICIT NONE
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation

  TYPE(Element_t),POINTER :: Element
  TYPE(ValueList_t), POINTER :: Material, SolverParams
  TYPE(Variable_t), POINTER :: PointerToVariable, bedrockVar, TimeVar
  TYPE(Nodes_t), SAVE :: Nodes

  LOGICAL :: AllocationsDone = .FALSE., GotIt, stat, firsttime=.TRUE.

  INTEGER :: nx,ny
  INTEGER :: i, mn, n, t, Nn, istat, DIM, MSum, ZSum, bedrockSource, Nnbr
  INTEGER :: TimesInd, old_times_ind
  INTEGER, POINTER :: Permutation(:), bedrockPerm(:)

  REAL(KIND=dp), POINTER :: VariableValues(:), BedrockValues(:)
  REAL(KIND=dp) :: x, y, z, toler, time
  REAL(KIND=dp), ALLOCATABLE :: zb(:)
  REAL(KIND=dp) :: times(73)
  REAL,ALLOCATABLE :: dem(:,:),xx(:),yy(:)


  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName = 'GroundedSolver', bedrockName, fmt_fname, filename

  INTEGER,PARAMETER :: MATERIAL_DEFAULT = 1, MATERIAL_NAMED = 2, VARIABLE = 3
       
  SAVE AllocationsDone, DIM, SolverName, zb, toler, fmt_fname, old_times_ind, times, firsttime, xx, yy, dem, ny, nx


  IF (FirstTime) THEN
      FirstTime = .FALSE.
      DO i = 0,72
         times(i + 1) = 0.25_dp * i + 1996.0_dp
      END DO
      fmt_fname = "(A, F7.2, A)"
      old_times_ind = -20.0
  END IF

!  NULLIFY(bedrockPerm,bedrockVar)

  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values

  CALL INFO(SolverName, 'Computing grounded mask from geometry', level=3)

  !--------------------------------------------------------------
  ! Allocate some permanent storage:
  !--------------------------------------------------------------
  IF ( (.NOT. AllocationsDone) .OR. Solver % Mesh % Changed ) THEN
     DIM = CoordinateSystemDimension()
     mn = Solver % Mesh % MaxElementNodes
     IF (AllocationsDone) DEALLOCATE(zb)     
     ALLOCATE(zb(mn), STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL FATAL( SolverName, 'Memory allocation error.' )
     END IF
     CALL INFO( SolverName, 'Memory allocation done.',Level=1 )
     AllocationsDone = .TRUE.
  END IF
  
  SolverParams => GetSolverParams()
  toler = GetConstReal(SolverParams, 'Toler', GotIt)
  IF (.NOT.GotIt) THEN
     CALL FATAL(SolverName, 'No tolerance given for the Grounded Mask.')
  END IF

  bedrockName = GetString(SolverParams, 'Bedrock Variable', GotIt)
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

  TimeVar => VariableGet(Solver % Mesh % Variables, 'Time')
  Time = TimeVar % Values(1)
    
  !--------------------------------------------------------------
  ! Grounded/floating loop based on height of base above bedrock.
  !--------------------------------------------------------------
  Time = Time + 1996.0_dp
  IF (Time <= 2014.0_dp) THEN
      TimesInd = MINLOC(ABS(Times - time), 1)
      IF (old_times_ind .NE. TimesInd) THEN
         write(filename, fmt_fname) "/nobackup/dlilien/smith_inputs/masks/interp_gl/newmask_interp_", Times(TimesInd), ".xyz"
         Firsttime=.False.
         call get_twod_grid(filename, xx, yy, dem)
         WRITE(Message,'(A)') 'Loaded new gmask'
         CALL INFO('GroundedSolver',Message,Level=3)
         nx = SIZE(xx)
         ny = SIZE(yy)
         old_times_ind = TimesInd
      END IF
  END IF


  DO t = 1, Solver % NumberOfActiveElements
     Element => GetActiveElement(t)
     n = GetElementNOFNodes()
     
     SELECT CASE(bedrockSource)
     CASE (VARIABLE)
        bedrockVar => VariableGet(Model % Mesh % Variables, bedrockName )
        IF (.NOT. ASSOCIATED(bedrockVar)) CALL FATAL(SolverName,"Could not find bedrock variable")
        bedrockPerm => bedrockVar % Perm
        zb(1:n) =  bedrockVar % values(bedrockPerm(Element % NodeIndexes)) + toler
     CASE (MATERIAL_NAMED)
        Material => GetMaterial( Element )
        zb(1:n) = ListGetReal( Material,bedrockName, n , & 
             Element % NodeIndexes, GotIt ) + toler
        IF (.NOT. GotIt) CALL FATAL(SolverName,"Could not find bedrock material")
     CASE (MATERIAL_DEFAULT)
        Material => GetMaterial( Element )
        zb(1:n) = ListGetReal( Material,'Min Zs Bottom',n , & 
             Element % NodeIndexes, GotIt ) + toler
        IF (.NOT. GotIt) CALL FATAL(SolverName,"Could not find bedrock material")
     END SELECT
     
     CALL GetElementNodes( Nodes )

     BedrockValues => bedrockVar % Values
     
     DO i = 1, n
        Nn = Permutation(Element % NodeIndexes(i))
        Nnbr = bedrockPerm(Element % NodeIndexes(i))
        IF (Nn==0) CYCLE
        IF (DIM == 2) THEN
           z = Nodes % y( i )
        ELSE IF (DIM == 3) THEN
           z = Nodes % z( i )
        END IF

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
  
  !--------------------------------------------------------------
  ! Grounding line loop to label grounded points at grounding line. Do
  ! regardless of time
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
              IF (DIM==2) PRINT *, 'Grounding Line, x', Nodes % x( i )
              IF (DIM==3) PRINT *, 'Grounding Line, (x,y)', Nodes % x( i ), Nodes % y( i )
           END IF
        END DO
     END IF
  END DO
  
  IF ( ParEnv % PEs>1 ) CALL ParallelSumVector( Solver % Matrix, VariableValues, 1 )
 
  CALL INFO( SolverName , 'Done')
 
END SUBROUTINE GroundedSolver 
