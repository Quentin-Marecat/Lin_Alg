SUBROUTINE INV_MAT(DIM,MATRIX,VEC,SOL,CHECK)
    ! -------------------------------------------------------------- !
    ! --- THIS SUBROUTINE FIND LINEAR SOLUTION OF A.x = b SYSTEM --- !
    ! --- USING INVERSION OF MATRIX -------------------------------- !
    ! --- NOT RECOMMENDED TO SOLVE SUCH SYSTEMS -------------------- !
    ! -------------------------------------------------------------- !
    IMPLICIT NONE
    REAL*8, PARAMETER :: EPS0 = 1.D-15, CONV = 1.D-13, EPS = 1.D-12
    INTEGER, PARAMETER :: MAXSTEP = 1D2
    LOGICAL, INTENT(IN) :: CHECK
    LOGICAL :: TEST
    INTEGER, INTENT(IN) :: DIM
    REAL*8, INTENT(IN) :: MATRIX(DIM,DIM), VEC(DIM)
    REAL*8, INTENT(OUT) :: SOL(DIM)
    INTEGER :: I, J, ERR, STAT, COMPT
    REAL*8 :: INVM(DIM,DIM), VECTEST(DIM), XNORME
    REAL*8 :: VEC2(DIM), SOL2(DIM),VECTAMP(DIM)

    CALL INV_QR(DIM,MATRIX,INVM,.FALSE.)
!    CALL INV_LU(DIM,MATRIX,INVM,.FALSE.)

    VEC2 = VEC
    OPEN(UNIT = 97, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 2) THEN
        COMPT = 0
        XNORME = 1
        SOL2 = 0
        ! --- RESIDUAL METHOD IMPROVEMENT --- !
        DO WHILE(XNORME > CONV .AND. COMPT < MAXSTEP) 
            CALL MATAPPLI(DIM,INVM,VEC2,SOL)
            COMPT = COMPT + 1
            CALL NORMVEC(DIM,SOL,XNORME)
            CALL MATAPPLI(DIM,MATRIX,SOL,VECTAMP)
            ! --- RESIDUAL METHOD IMPROVEMENT --- !
            DO I = 1,DIM
                SOL2(I) = SOL2(I) + SOL(I)
                VEC2(I) = VEC2(I) - VECTAMP(I) 
            ENDDO
        ENDDO
    ENDIF
    SOL = SOL2

    ! --- VERIFICATION --- !
    IF (CHECK) THEN
        CALL MATAPPLI(DIM,MATRIX,SOL,VECTEST)
        TEST = .FALSE.
        DO I = 1,DIM
            IF (ABS(VECTEST(I) - VEC(I)) > EPS) TEST = .TRUE.
        ENDDO
        IF (TEST) THEN
            OPEN(UNIT = ERR,FILE = 'error')
            WRITE(ERR,'(A)') 'PROBLEM SOLUTION'
            WRITE(ERR,'(A)') 'A = '
            DO I = 1,DIM
                WRITE(ERR,'(100F14.5)') (MATRIX(I,J), J = 1,DIM)
            ENDDO
            WRITE(ERR,'(A)') '************'
            WRITE(ERR,'(A)') 'x = '
            WRITE(ERR,'(100F14.5)') (SOL(I), I = 1,DIM)
            WRITE(ERR,'(A)') '************'
            WRITE(ERR,'(A)') 'b = '
            WRITE(ERR,'(100F14.5)') (VEC(I), I = 1,DIM)
            WRITE(ERR,'(A)') '************'
            WRITE(ERR,'(A,4X,ES14.5)')'EPS',EPS
            WRITE(ERR,'(A)') 'ERROR = MAT*SOL - VEC '
            WRITE(ERR,'(100ES14.4)') (VECTEST(I)-VEC(I), I = 1,DIM)
            WRITE(ERR,'(A)') '************'
        ENDIF
    ENDIF

    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) WRITE(6,'(A)') 'ERROR LINEAR SOLUTION, SEE FILE error'

END SUBROUTINE