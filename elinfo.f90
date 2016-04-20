



      SUBROUTINE elinfo(Model,Solver,dt,TransientSimulation)
        USE DefUtils
        USE MaterialModels
          IMPLICIT NONE
        !------------------------------------------------------------------------------
          TYPE(Solver_t) :: Solver
          TYPE(Model_t) :: Model

          REAL(KIND=dp) :: dt
          LOGICAL :: TransientSimulation
        !  
        !!!!  Variables utiles pour les elements et les fonctions de base
          TYPE(Element_t),POINTER ::  Element
          TYPE(Nodes_t) :: ElementNodes
          TYPE(GaussIntegrationPoints_t) :: IntegStuff
          TYPE(ValueList_t), POINTER :: SolverParams,Material
          real(kind=dp),allocatable :: Basis(:),dBasisdx(:,:)
          real(kind=dp) :: u,v,w,SqrtElementMetric
          INTEGER, POINTER :: NodeIndexes(:)
          CHARACTER(LEN=MAX_NAME_LEN) :: SolverName

        !!!!! variables Elmer
          TYPE(Variable_t), POINTER :: GradVariable, Variable, VeloSolN,VeloSolD
          REAL(KIND=dp), POINTER ::  GradValues(:),VelocityN(:),VelocityD(:),Values(:)
          INTEGER, POINTER :: GradPerm(:), VeloNPerm(:),VeloDPerm(:),Perm(:)
          CHARACTER(LEN=MAX_NAME_LEN) ::GradSolName,NeumannSolName,DirichletSolName,VarSolName

        !! autres variables
          real(kind=dp),allocatable :: VisitedNode(:),db(:)
          real(kind=dp),allocatable,dimension(:) :: Ux,Uy,Uz
          real(kind=dp) :: Velo(3),dVelodx(3,3)
          real(kind=dp) :: s,ss,c2,c3
          real(kind=dp) :: mub,Viscosityb
          real(kind=dp),allocatable,dimension(:) :: c2n,c3n
          real(kind=dp),allocatable,dimension(:) :: NodalViscosityb


          integer :: i,j,t,n,NMAX,NpN,NActiveNodes,DIM,e,p,q

          CHARACTER(LEN=MAX_NAME_LEN) :: ViscosityFlag

          logical :: SquareFormulation
          Logical ::  Firsttime=.true.,Found,stat,gotit


          save Firsttime,DIM
          save ElementNodes
          save SolverName
          save NeumannSolName,DirichletSolName,VarSolName,GradSolName
          save SquareFormulation
          save VisitedNode,db,Basis,dBasisdx
          save Ux,Uy,Uz
          save c2n,c3n
          save NodalViscosityb

          !!!! Firsttime Do some allocation and initialisation
          If (Firsttime) then

              DIM = CoordinateSystemDimension()
              WRITE(SolverName, '(A)') 'DJDmu_Adjoint_Lilien'

              NMAX=Solver % Mesh % NumberOfNodes
              NpN=Model % MaxElementNodes

              allocate(VisitedNode(NMAX),db(NMAX), &
                       Basis(NpN),  &
                       dBasisdx(NpN,3), &
                       Ux(NpN),Uy(NpN),Uz(NpN),&
                       c2n(NpN),c3n(NpN),&
                       NodalViscosityb(NpN))

        !!!!!!!!!!! get Solver Variables
              SolverParams => GetSolverParams()

              NeumannSolName =  GetString( SolverParams,'Flow Solution Name', Found)
                  IF(.NOT.Found) THEN        
                       CALL WARN(SolverName,'Keyword >Flow Solution Name< not found in section >Solver<')
                       CALL WARN(SolverName,'Taking default value >Flow Solution<')
                       WRITE(NeumannSolName,'(A)') 'Flow Solution'
                  END IF
              DirichletSolName =  GetString( SolverParams,'Adjoint Solution Name', Found)
                  IF(.NOT.Found) THEN        
        [<32;5;0M               CALL WARN(SolverName,'Keyword >Adjoint Solution Name< not found in section >Solver<')
                       CALL WARN(SolverName,'Taking default value >VeloD<')
                       WRITE(DirichletSolName,'(A)') 'VeloD'
                  END IF
              VarSolName =  GetString( SolverParams,'Optimized Variable Name', Found)
                     IF(.NOT.Found) THEN
                            CALL WARN(SolverName,'Keyword >Optimized Variable Name< not found  in section >Solver<')
                            CALL WARN(SolverName,'Taking default value >Mu<')
                            WRITE(VarSolName,'(A)') 'Mu'
                      END IF
              GradSolName =  GetString( SolverParams,'Gradient Variable Name', Found)
                     IF(.NOT.Found) THEN
                            CALL WARN(SolverName,'Keyword >Gradient Variable Name< not found  in section >Solver<')
                            CALL WARN(SolverName,'Taking default value >DJDMu<')
                            WRITE(GradSolName,'(A)') 'DJDmu'
                     END IF
               SquareFormulation=GetLogical( SolverParams, 'SquareFormulation', Found)
                   IF(.NOT.Found) THEN
                           CALL WARN(SolverName,'Logical Keyword >SquareFormulation< not found  in section >Solver<')
                           CALL WARN(SolverName,'Taking default value >FALSE<')
                           SquareFormulation=.FALSE.
                   END IF
          
          !!! End of First visit
            Firsttime=.false.
          Endif


      END SUBROUTINE elinfo
