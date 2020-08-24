SUBROUTINE JACOBI(DIM,MATRIX,VEC,SOL,CHECK)
    ! --- FOR ANY REPORT OR SUGGESTION, PLEASE CONTACT quentin.marecat@etu.umontpellier.fr --- !
    ! ------------------------------------------------------------------------------------------------------ !
    ! --- THIS SUBROUTINE FIND THE SOLUTION OF THE SYSTEM MATRIX*SOL = VEC USING JACOBI ITERATIVE METHOD --- !
    ! --- MATRIX MUST BE SQUARE AND INVERSIBLE ------------------------------------------------------------- !
    ! --- RESIDUAL METHOD IS USED TO IMPROVE THE METHOD ---------------------------------------------------- !
    ! --- operation_mat.f90 ARE NECESSARY ----------------------------- !
    ! --- BAD AND OUTDATED METHOD, OFTEN WORKS WHEN MATRIX IS DIAGONAL DOMINANT ---------------------------- !
    ! ------------------------------------------------------------------------------------------------------ !
    IMPLICIT NONE
    LOGICAL,INTENT(IN) :: CHECK
    REAL*8, PARAMETER :: EPS0 = 1.D-15, CONV = 1.D-13, EPS = 1.D-10, CONV2 = 1.D-13
    INTEGER, PARAMETER :: MAXSTEP = 1D4, MAXSTEP2 = 1D2
    LOGICAL :: TEST
    INTEGER, INTENT(IN) :: DIM
    REAL*8, INTENT(IN) :: MATRIX(DIM,DIM),VEC(DIM)
    REAL*8,INTENT(OUT) :: SOL(DIM)
    REAL*8 ::  DM1(DIM,DIM),W(DIM,DIM) VECTEST(DIM),XNORME
    REAL*8 ::  DIFF, DM1W(DIM,DIM), DM1VEC(DIM), SOLP1(DIM)
    REAL*8 :: SOL2(DIM), VEC2(DIM)
    INTEGER :: I,J, COMPT,COMPT2, ERR, STAT

    ERR = 97
    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) CLOSE(ERR,STATUS = 'delete')

    DM1 = 0
    W = 0
    DO I = 1,DIM
        IF (ABS(MATRIX(I,I)) < EPS0) THEN
            OPEN(UNIT = ERR, FILE = 'error')
            WRITE(ERR,'(A)') 'ERROR 0 DIAGONAL'
            WRITE(6,'(A)') 'ERROR LINEAR SOLUTION, SEE FILE error'
            STOP
        ELSE
            DM1(I,I) = 1/MATRIX(I,I)
        ENDIF
        DO J = 1,DIM
            W(J,I) = MATRIX(J,I)
        ENDDO
        W(I,I) = 0
    ENDDO

    CALL PRODMAT(DIM,DM1,W,DM1W)
    VEC2 = VEC
    SOL2 = 0
    COMPT2 = 0
    XNORME = 1
    ! --- RESIDUAL METHOD IMPROVEMENT --- !
    DO WHILE(XNORME > CONV2 .AND. COMPT2 < MAXSTEP2) 
        CALL MATAPPLI(DIM,DM1,VEC2,DM1VEC)
        DIFF = 1
        COMPT = 0
        SOL = DM1VEC
        DO WHILE(DIFF > CONV .AND. COMPT < MAXSTEP)
            CALL MATAPPLI(DIM,DM1W,SOL,SOLP1)
            DO I = 1,DIM
                SOL(I) = SOL(I) - SOLP1(I)
            ENDDO
            ! --- CONVERGENCE --- !
            CALL MATAPPLI(DIM,MATRIX,SOL,VECTEST)
            VECTEST = VECTEST - VEC
            CALL NORMVEC(DIM,VECTEST,DIFF)
            ! ------------------- !
            COMPT = COMPT + 1
        ENDDO
        ! --- RESIDUAL METHOD IMPROVEMENT --- !
        COMPT2 = COMPT2 + 1
        CALL NORMVEC(DIM,SOL,XNORME)
        CALL MATAPPLI(DIM,MATRIX,SOL,VECTAMP)
        DO I = 1,DIM
            SOL2(I) = SOL2(I) + SOL(I)
            VEC2(I) = VEC2(I) - VECTAMP(I) 
        ENDDO
    ENDDO
SOL = SOL2

    ! --- VERIFICATION --- !
    IF(CHECK) THEN
        TEST = .FALSE.
        IF (COMPT == MAXSTEP) THEN
            TEST = .TRUE.
            OPEN(UNIT = ERR, FILE = 'error')
            WRITE(ERR,'(A,4X,A,4X,I15)') 'MAXSTEP REACHS','MAXSTEP =',MAXSTEP
        ENDIF
        CALL MATAPPLI(DIM,MATRIX,SOL,VECTEST)
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
        

