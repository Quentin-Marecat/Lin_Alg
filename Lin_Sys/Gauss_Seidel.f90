SUBROUTINE GAUSS_SEIDEL(DIM,MATRIX,VEC,SOL)
    ! --- FOR ANY REPORT OR SUGGESTION, PLEASE CONTACT quentin.marecat@etu.umontpellier.fr --- !
    ! ------------------------------------------------------------------------------------------------------ !
    ! --- THIS SUBROUTINE FIND THE SOLUTION OF THE SYSTEM MATRIX*SOL = VEC USING JACOBI ITERATIVE METHOD --- !
    ! --- MATRIX MUST BE SQUARE AND INVERSIBLE ------------------------------------------------------------- !
    ! --- RESIDUAL METHOD IS USED TO IMPROVE THE METHOD ---------------------------------------------------- !
    ! --- operation_mat.f90 InvQR.f90 QRhouse.f90 AND housestep.f90 ARE NECESSARY ----------------------------- !
    ! --- EXTREMELY BAD AND OUTDATED METHOD, OFTEN WORKS WHEN MATRIX IS DIAGONAL DOMINANT ------------------ !
    ! ------------------------------------------------------------------------------------------------------ !
    IMPLICIT NONE
    REAL*8, PARAMETER :: EPS0 = 1.D-15, CONV = 1.D-10
    INTEGER, PARAMETER :: MAXSTEP = 1D4
    LOGICAL :: TEST
    INTEGER, INTENT(IN) :: DIM
    REAL*8, INTENT(IN) :: MATRIX(DIM,DIM),VEC(DIM)
    REAL*8,INTENT(OUT) :: SOL(DIM)
    REAL*8 ::  VECTEST(DIM), AM1(DIM,DIM)
    REAL*8 :: A1(DIM,DIM), A2(DIM,DIM), DIFF, M(DIM,DIM), V(DIM),VECP1(DIM)
    INTEGER :: I,J, COMPT, ERR, STAT

    ERR = 97
    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) CLOSE(ERR,STATUS = 'delete')

    A1 = 0
    A2 = 0
    DO I = 1,DIM-1
        A1(I,I) = MATRIX(I,I)
        DO J = I+1,DIM
            A1(J,I) = MATRIX(J,I)
            A2(I,J) = MATRIX(I,J)
        ENDDO
    ENDDO
    A1(DIM,DIM) = MATRIX(DIM,DIM)
    CALL INVMATQR(DIM,A1,AM1)
    CALL PRODMAT(DIM,AM1,A2,M)
    DIFF = 1
    COMPT = 0
    DO I = 1,DIM
        SOL(I) = 1
    ENDDO
    M = -M
    CALL MATAPPLI(DIM,AM1,VEC,V)
    DO WHILE(DIFF > CONV .AND. COMPT < MAXSTEP)
        CALL MATAPPLI(DIM,M,SOL,VECP1)
        DO I = 1,DIM
            VECP1(I) = VECP1(I) + V(I)
        ENDDO
        CALL MATAPPLI(DIM,MATRIX,VECP1,VECTEST)
        VECTEST = VECTEST - VEC
        CALL NORMVEC(DIM,VECTEST,DIFF)
        SOL = VECP1
        COMPT = COMPT + 1
    ENDDO

    ! --- VERIFICATION --- !
    TEST = .FALSE.
    IF (COMPT == MAXSTEP) THEN
        TEST = .TRUE.
        OPEN(UNIT = ERR, FILE = 'error')
        WRITE(ERR,'(A,4X,A,4X,I15)') 'MAXSTEP REACHS','MAXSTEP =',MAXSTEP
    ENDIF
    CALL MATAPPLI(DIM,MATRIX,SOL,VECTEST)
    DO I = 1,DIM
        IF (ABS(VECTEST(I) - VEC(I)) > CONV) TEST = .TRUE.
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
        WRITE(ERR,'(A)') 'ERROR = '
        WRITE(ERR,'(100E14.4)') (VECTEST(I)-VEC(I), I = 1,DIM)
        WRITE(ERR,'(A)') '************'
    ENDIF

    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) WRITE(6,'(A)') 'ERROR LINEAR SOLUTION, SEE FILE error'

END SUBROUTINE
        

