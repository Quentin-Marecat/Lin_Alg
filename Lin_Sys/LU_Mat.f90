SUBROUTINE LU_MAT(DIM,MATRIX,VEC,SOL)
    ! ----------------------------------------------------------------------------------------------- !
    ! --- THIS SUBROUTINE FIND THE SOLUTION OF THE SYSTEM MATRIX*SOL = VEC USING LU DECOMPOS OF A --- !
    ! --- MATRIX MUST BE SQUARE AND INVERSIBLE ------------------------------------------------------ !
    ! ----------------------------------------------------------------------------------------------- !
    IMPLICIT NONE 
    REAL*8, PARAMETER :: EPS0 = 1.D-14
    INTEGER, INTENT(IN) :: DIM
    REAL*8, INTENT(IN) :: MATRIX(DIM,DIM),VEC(DIM)
    REAL*8, INTENT(OUT) :: SOL(DIM)
    REAL*8 :: L(DIM,DIM),U(DIM,DIM), VECTEST(DIM), SOLINT(DIM)
    LOGICAL :: TEST
    INTEGER :: I, J, ERR, STAT

    ERR = 97
    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) CLOSE(ERR,STATUS = 'delete')

    ! ---------------------------------------------------------------------------------------------------- !
    CALL LU(DIM,MATRIX,L,U)
    ! --- LU DECOMP SUBROUTINE FROM Quentin-Marecat PACKAGES IS USED, PLEASE FIND IT ON MY GITHUB PAGE --- !

    OPEN(UNIT = 97, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 2) THEN
        SOLINT(1) = VEC(1)/L(1,1)
        DO I = 2,DIM
            SOLINT(I) = VEC(I)
            DO J = 1,I-1
                SOLINT(I) = SOLINT(I) - SOLINT(J)*L(I,J)
            ENDDO
            SOLINT(I) = SOLINT(I)/L(I,I)
        ENDDO
        SOL(DIM) = SOLINT(DIM)/U(DIM,DIM)
        DO I = DIM-1,1,-1
            SOL(I) = SOLINT(I)
            DO J = DIM,I+1,-1
                SOL(I) = SOL(I) - SOL(J)*U(I,J)
            ENDDO
            SOL(I) = SOL(I)/U(I,I)
        ENDDO
    ENDIF

    ! --- VERIFICATION --- !
    CALL MATAPPLI(DIM,MATRIX,SOL,VECTEST)
    TEST = .FALSE.
    DO I = 1,DIM
        IF (ABS(VECTEST(I) - VEC(I)) > EPS0) TEST = .TRUE.
    ENDDO
    IF (TEST) THEN
        WRITE(6,*) SOL
        WRITE(6,*) VECTEST
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
    ENDIF

    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) WRITE(6,'(A)') 'ERROR LINEAR SOLUTION, SEE FILE error'

END SUBROUTINE
