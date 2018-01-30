      FUNCTION cubicmeters2gt(m3) RESULT(gt)
          USE Types
          implicit NONE

          real(kind=dp) :: m3
          real(kind=dp) :: gt

          gt = m3 * 910.0_dp * 1.0e-12_dp
          RETURN
      END

      SUBROUTINE TimeVariableMeltrate( Model,Solver,dt,TransientSimulation )
        USE types
        USE DefUtils
        USE ParallelUtils
        USE MeltFunctions
        implicit none

        TYPE(Model_t) :: Model
        TYPE(Solver_t) :: Solver
        TYPE(Element_t),POINTER :: Element, CurrentElement
        TYPE(Variable_t), POINTER :: PointerToVariable, bmpointer, TimeVar
        TYPE(GaussIntegrationPoints_t) :: IP
        TYPE(Nodes_t) :: Nodes

        REAL(KIND=dp) :: dt, total_melt, element_area, detJ, global_melt
        REAL(KIND=dp) :: mr, mf, melt_ratio, cubicmeters2gt, time
        LOGICAL :: TransientSimulation, stat, scalemelt, tv
        LOGICAL :: getSecondDerivatives=.FALSE., found=.FALSE.
        LOGICAL :: firsttime=.TRUE.
        Integer :: t, n, k, i, p
        CHARACTER(LEN=MAX_NAME_LEN) :: melt_var,&
        solvername='TimeVariableMeltrate'

        REAL(KIND=dp) :: Basis(3), dBasisdx(3,3), ddBasisddx(3,3,3), basemelt(6)
        INTEGER, POINTER :: Permutation(:), NodeIndexes(:), bmperm(:)
        REAL(KIND=dp), POINTER :: VariableValues(:), bmvals(:)

        save tv, mr, mf, scalemelt, firsttime, melt_var, melt_ratio

        if ( firsttime ) THEN
            melt_var = GetString(GetSolverParams(), 'Melt Variable', found)
            IF (.NOT.Found) THEN
                melt_var = 'BaseMelt'
                CALL WARN(SolverName,'Keyword Melt Variable not found')
                CALL WARN(SolverName,'Using "BaseMelt"')
            END IF 

            tv = GetLogical(GetSolverParams(), 'Time Variable', found)
            IF (.NOT.Found) THEN
                tv = .FALSE.
                CALL WARN(SolverName,'Keyword Time Variable not found')
                CALL WARN(SolverName,'Melt rate is time-independent')
            END IF 

            IF (.NOT.tv) THEN
                mr = GetConstReal(GetSolverParams(), 'Melt Rate', found)
                IF (.NOT.Found) THEN
                    CALL FATAL(SolverName,'No MR or TV!!!!')
                END IF
            END IF

            mf = GetConstReal(GetSolverParams(), 'Melt Factor', found)
            IF (.NOT.Found) THEN
                mf = 1
                CALL WARN(SolverName,'No MR or MF, no scaling')
            END IF

            bmpointer => VariableGet(Solver % Mesh % Variables, melt_var)
            IF (.NOT.ASSOCIATED( bmpointer) ) THEN
                CALL FATAL(solvername, 'Base melt variable unfound')
            END IF
            firsttime = .FALSE.

        END IF


        bmpointer => VariableGet(Solver % Mesh % Variables, melt_var)
        bmperm => bmpointer % Perm
        bmvals => bmpointer % Values

        PointerToVariable => Solver % Variable
        Permutation  => PointerToVariable % Perm
        VariableValues => PointerToVariable % Values

        ! We can skip most of the calculations if we are not time
        ! variable
        IF (tv.OR.firsttime) THEN
            ! we start by integrating to find the total melt
            total_melt = 0.0_dp
            DO t=1,Solver % NumberOfActiveElements
                CurrentElement => GetActiveElement(t)
                CALL GetElementNodes(Nodes)

                n = GetElementNOFNodes(CurrentElement)
                NodeIndexes => CurrentElement % NodeIndexes
                do i=1,n
                    basemelt(i) = bmvals(bmperm(NodeIndexes(i)))
                end do

                IP = GaussPoints(CurrentElement)
                stat = ElementInfo( CurrentElement, Nodes, IP % U(t), IP % V(t), &
                        IP % W(t), detJ, Basis, dBasisdx, ddBasisddx, &
                        getSecondDerivatives)
                ! The two in the next line is a little worrisome--i'm not
                ! sure why i need it
                total_melt = total_melt + SUM( BaseMelt(1:IP % n) * Basis(1:IP % n) * detJ) / 2.0_dp
            END DO

            ! this line merges parallel partitions to one. flag arg is
            ! sum/min/max
            global_melt = ParallelReduction(total_melt, 0)
            global_melt = cubicmeters2gt(global_melt)
            
            ! Now if we are time variable we need to find our present
            ! melt ratio
            IF (tv) THEN
                TimeVar => VariableGet(Solver % Mesh % Variables, 'Time')
                Time = TimeVar % Values(1)
                mr = cubicmeters2gt(41309491925.9_dp) * RescaleByTime(Time + 1996.0_dp)
            END IF

            melt_ratio = mr / global_melt
        END IF

        ! Now we need to update everything
        DO t=1,Solver % NumberOfActiveElements
            CurrentElement => GetActiveElement(t)
            n = GetElementNOFNodes(CurrentElement)
            NodeIndexes => CurrentElement % NodeIndexes
            do i=1,n
                basemelt(i) = bmvals(bmperm(NodeIndexes(i)))
                k = Permutation(NodeIndexes(i))
                VariableValues(k) = basemelt(i) * melt_ratio
            end do
        END DO

        Return
      END
