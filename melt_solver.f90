      FUNCTION cubicmeters2gt(m3) RESULT(gt)
          USE Types
          IMPLICIT NONE

          REAL(KIND=dp) :: m3
          REAL(KIND=dp) :: gt

          gt = m3 * 910.0_dp * 1.0e-12_dp
          RETURN
      END

      SUBROUTINE TimeVariableMeltrate( Model,Solver,dt,TransientSimulation )
        USE types
        USE DefUtils
        USE ParallelUtils
        USE MeltFunctions
        IMPLICIT NONE

        TYPE(Model_t) :: Model
        TYPE(Solver_t) :: Solver
        TYPE(Element_t),POINTER :: Element, CurrentElement
        TYPE(Variable_t), POINTER :: PointerToVariable, bmpointer, TimeVar
        TYPE(GaussIntegrationPoints_t) :: IP
        TYPE(Nodes_t) :: Nodes

        REAL(KIND=dp) :: dt, total_melt, element_area, detJ, global_melt
        REAL(KIND=dp) :: mr, mf, melt_ratio, cubicmeters2gt, time, element_melt
        REAL(KIND=dp) :: meltrate, max_melt_rate=3000.0
        REAL(KIND=dp) :: unpin_time
        LOGICAL :: TransientSimulation, stat, scalemelt, tv, UnfoundFatal
        LOGICAL :: getSecondDerivatives=.FALSE., found=.FALSE.
        LOGICAL :: firsttime=.TRUE.
        INTEGER :: t, n, k, i, p
        CHARACTER(LEN=MAX_NAME_LEN) :: melt_var,&
        solvername='TimeVariableMeltrate',&
        MeltFile='melt.dat'

        REAL(KIND=dp) :: Basis(3), dBasisdx(3,3), ddBasisddx(3,3,3), basemelt(6)
        INTEGER, POINTER :: Permutation(:), NodeIndexes(:), bmperm(:)
        REAL(KIND=dp), POINTER :: VariableValues(:), bmvals(:)

        SAVE tv, mr, mf, scalemelt, firsttime, melt_var, melt_ratio, global_melt, unpin_time, bmperm, bmvals

        PointerToVariable => Solver % Variable
        Permutation  => PointerToVariable % Perm
        VariableValues => PointerToVariable % Values

        IF ( firsttime ) THEN
            melt_var = GetString(GetSolverParams(), 'Melt Variable', found)

            IF (.NOT.Found) THEN
                melt_var = 'BaseMelt'
                IF (ParEnv % myPe == 0) THEN 
                    CALL WARN(SolverName,&
                              'Keyword Melt Variable not found')
                    CALL WARN(SolverName,'Using "BaseMelt"')
                END IF
            END IF 

            tv = GetLogical(GetSolverParams(), 'Time Variable', found)
            IF (.NOT.Found) THEN
                tv = .FALSE.
                IF (ParEnv % myPe == 0) THEN 
                    CALL INFO(SolverName,'Melt rate is time-independent', level=3)
                END IF
            END IF 

            IF (.NOT.tv) THEN
                mr = GetConstReal(GetSolverParams(), 'Melt Rate', found)
                IF (.NOT.Found) THEN
                    mr = HUGE(1.0_dp)
                    IF (ParEnv % myPe == 0) THEN 
                        CALL INFO(SolverName,&
                        'No MR or TV, you are only using MF', level=3)
                    END IF
                ELSE
                    IF (ParEnv % myPe == 0) THEN 
                        CALL INFO(SolverName,'Melt rate is pinned to MR', level=3)
                    END IF
                END IF
            END IF

            mf = GetConstReal(GetSolverParams(), 'Melt Factor', found)
            IF (.NOT.Found) THEN
                mf = 1.0_dp
                IF (ParEnv % myPe == 0) THEN 
                    CALL INFO(SolverName,'No mf found, using 1.0', level=3)
                END IF
            ELSE
                IF (ParEnv % myPe == 0) THEN 
                    CALL INFO(SolverName, 'MF Found', level=3)
                END IF
            END IF

            unpin_time = GetConstReal(GetSolverParams(), 'Unpin Time', found)
            IF (.NOT.Found) THEN
                unpin_time = HUGE(1.0_dp)
                IF (ParEnv % myPe == 0) THEN 
                    CALL INFO(SolverName,&
                              'Not unpinning from target time', level=3)
                END IF
            ELSE
                IF (ParEnv % myPe == 0) THEN 
                    CALL INFO(SolverName, 'Unpinning', level=3)
                END IF
            END IF

            bmpointer => VariableGet(Solver % Mesh % Variables, melt_var, UnfoundFatal=UnfoundFatal)

            IF (.NOT.ASSOCIATED( bmpointer).OR.(UnfoundFatal) ) THEN
                IF (ParEnv % myPe == 0) THEN 
                    CALL FATAL(SolverName, 'Base melt variable unfound')
                END IF
            ELSE
                bmperm => bmpointer % Perm
                bmvals => bmpointer % Values
                IF (ParEnv % myPe == 0) THEN 
                    CALL INFO(SolverName, 'Found BaseMelt Variable', level=3)
                END IF
            END IF

            OPEN (12, FILE=MeltFile)
                WRITE(12, *) 'Time, Target MR, Global MR, Melt Ratio'
            CLOSE(12)

            melt_ratio = 1.0_dp
        END IF

        ! We can skip most of the calculations if we are not time
        ! variable. Also i think we can just do this once per timestep
        IF (GetCoupledIter() == 1) THEN
            ! we start by integrating to find the total melt
            total_melt = 0.0_dp
            DO t=1,Solver % NumberOfActiveElements
                CurrentElement => GetActiveElement(t)
                CALL GetElementNodes(Nodes)

                n = GetElementNOFNodes(CurrentElement)
                IF (n <= 2) THEN
                    CYCLE
                END IF
                NodeIndexes => CurrentElement % NodeIndexes
                DO i=1,n
                    basemelt(i) = bmvals(bmperm(NodeIndexes(i)))
                END DO

                IP = GaussPoints(CurrentElement)
                stat = ElementInfo( CurrentElement, Nodes, IP % U(t), IP % V(t), &
                        IP % W(t), detJ, Basis, dBasisdx, ddBasisddx, &
                        getSecondDerivatives)
                ! The two in the next line is a little worrisome--i'm not
                ! sure why i need it
                element_melt = SUM( BaseMelt(1:IP % n) * Basis(1:IP % n) * detJ) / 2.0_dp
                meltrate = SUM( BaseMelt(1:IP % n) * Basis(1:IP % n))


                IF ((element_melt /= 0.0_dp).AND.(element_melt == element_melt).AND.(ABS(meltrate) < max_melt_rate)) THEN
                    total_melt = total_melt + element_melt
                END IF
            END DO

            ! this line merges parallel partitions to one. flag arg is sum/min/max
            global_melt = ParallelReduction(total_melt, 0)
        END IF

        ! Now if we are time variable we need to find our present melt ratio
        ! We need to avoid overriding a fixed rate from previous solvers
        TimeVar => VariableGet(Solver % Mesh % Variables, 'Time')
        Time = TimeVar % Values(1)

        IF (tv.AND.((Time <= unpin_time).OR.firsttime)) THEN
            mr = 41309491925.9_dp * RescaleByTime(Time + 1996.0_dp)
        END IF

        ! we only update the melt_ratio if we are less than the unpin time
        IF ((Time <= unpin_time).OR.firsttime) THEN
            IF (abs(mr) >= HUGE(1.0_dp) / 2) THEN
                melt_ratio = mf
            ELSE
                melt_ratio = abs(mr / global_melt * mf)
            END IF
            IF (melt_ratio /= melt_ratio) THEN
                melt_ratio = 1.0_dp
            END IF
            firsttime = .FALSE.
        END IF

        ! Always write stuff out
        IF ((ParEnv % myPe == 0).AND.(GetCoupledIter() == 1)) THEN 
            OPEN (12, FILE=MeltFile, POSITION='APPEND')
                WRITE(12, '(E14.7,2x,E14.7,2x,E14.7,2x,E14.7)') Time, mr, global_melt, melt_ratio
            CLOSE(12)
        END IF

        ! Now we need to update everything
        DO t=1,Solver % NumberOfActiveElements
            CurrentElement => GetActiveElement(t)
            n = GetElementNOFNodes(CurrentElement)
            NodeIndexes => CurrentElement % NodeIndexes
            DO i=1,n
                basemelt(i) = bmvals(bmperm(NodeIndexes(i)))
                k = Permutation(NodeIndexes(i))
                VariableValues(k) = basemelt(i) * melt_ratio
            END DO
        END DO
        RETURN
      END
