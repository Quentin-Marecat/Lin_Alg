SUBROUTINE CHOLESKY_LU(DIM,MATRIX,L,CHECK)
    ! ------------------------------------------------------------------------------------ !
    ! --- THIS SUBROUTINE DO THE CHOLESKY DECOMPOSITION OF A POSITIVE SYMMETRIC MATRIX --- !
    ! --- USING LU FACTORISATION / LU SUBROUTINE IS NECESSARY ---------------------------- !
    ! ------------------------------------------------------------------------------------ !
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: CHECK
    REAL*8, PARAMETER :: EPS0 = 1.D-15, EPS = 1.D-9
    LOGICAL :: TEST
    INTEGER, INTENT(IN) :: DIM
    REAL*8, INTENT(IN) :: MATRIX(DIM,DIM)
    REAL*8, INTENT(OUT) :: L(DIM,DIM)
    REAL*8 :: U(DIM,DIM), DIAGMAT(DIM,DIM), LT(DIM,DIM),MATTEST(DIM,DIM)
    INTEGER :: I,J,K,ERR, STAT
    ERR = 97
    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) CLOSE(ERR,STATUS = 'delete')

    DIAGMAT = 0.
    CALL LU(DIM,MATRIX,L,U,.FALSE.)
    DO I = 1,DIM 
        IF (U(I,I) < EPS0 ) THEN
            OPEN(UNIT = ERR,FILE = 'error')
            WRITE(ERR,'(A,4X,F14.5)') 'ERROR NEGATIVE ELEMENT', U(I,I)
            WRITE(6,'(A)') 'ERROR, SEE FILE error'
            STOP
        ENDIF
        DIAGMAT(I,I) = SQRT(U(I,I))
    ENDDO
    CALL PRODMAT(DIM,L,DIAGMAT,L)

        ! --- VERIFICATION --- !
    IF (CHECK) THEN
        TEST = .FALSE.
        CALL TRANSPOSE(DIM,L,LT)
        CALL PRODMAT(DIM,L,LT,MATTEST)
        DO I = 1,DIM
            DO J = 1,DIM
                IF (ABS(MATTEST(I,J)-MATRIX(I,J)) > EPS) TEST = .TRUE.
            ENDDO
        ENDDO
        IF (TEST) THEN
            OPEN(UNIT = ERR,FILE = 'error')
            WRITE(ERR,'(A)') 'PROBLEM CHOLESKY DECOMPOSITION'
            WRITE(ERR,'(A)') 'L ='
            DO I = 1,DIM
                WRITE(ERR,'(100F14.5)') (L(I,J), J=1,DIM)
            ENDDO
            WRITE(ERR,'(A)')'***************'
            WRITE(ERR,'(A,4X,ES14.5)')'EPS',EPS
            WRITE(ERR,'(A)') 'MAT - L*LT ='
            DO I = 1,DIM
                WRITE(ERR,'(100ES14.5)') (MATRIX(I,J) - MATTEST(I,J), J=1,DIM)
            ENDDO
            WRITE(ERR,'(A)')'***************'
        ENDIF
    ENDIF

    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) WRITE(6,'(A)') 'ERROR, SEE FILE error'
END SUBROUTINE