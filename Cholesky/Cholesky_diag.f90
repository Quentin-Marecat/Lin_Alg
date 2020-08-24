SUBROUTINE CHOLESKY_DIAG(DIM,MATRIX,L,CHECK)
    ! ------------------------------------------------------------------------------------ !
    ! --- THIS SUBROUTINE DO THE CHOLESKY DECOMPOSITION OF A POSITIVE SYMMETRIC MATRIX --- !
    ! --- USING THE DIAGONALISATION OF M AND HIS HOUSEHOLDER DECOMPOSITION --------------- !
    ! --- DIAGONALISATION AND HOUSEHOLDER SUBROUTINE ARE NECESSARY ----------------------- !
    ! --- NOT NUMERICALLY INTERESTING AT ALL --------------------------------------------- !
    ! ------------------------------------------------------------------------------------ !
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: CHECK
    REAL*8, PARAMETER :: EPS0 = 1.D-15, EPS = 1.D-8
    LOGICAL :: TEST
    INTEGER, INTENT(IN) :: DIM
    REAL*8, INTENT(IN) :: MATRIX(DIM,DIM)
    REAL*8, INTENT(OUT) :: L(DIM,DIM)
    REAL*8 :: EIGENVAL(DIM), EIGENVEC(DIM,DIM)
    REAL*8 :: D(DIM,DIM), P(DIM,DIM), PT(DIM,DIM), DELTA(DIM,DIM)
    REAL*8 :: LT(DIM,DIM), MATTEST(DIM,DIM), ID(DIM,DIM)
    INTEGER :: I,J,K,ERR, STAT
    ERR = 97
    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) CLOSE(ERR,STATUS = 'delete')

    ID = 0.
    DO I = 1,DIM
        ID(I,I) = 1
    ENDDO

    ! --- STEP 1 : DIAGONALISATION --- !
!    CALL RQ_SHIFT_DIAG(DIM,MATRIX,EIGENVAL,EIGENVEC,.TRUE.,.TRUE.,.FALSE.)
    CALL JAC_DIAG(DIM,MATRIX,EIGENVAL,EIGENVEC)
    D = 0.
    DO I = 1,DIM
        IF (EIGENVAL(I) < EPS0) THEN
            OPEN(UNIT = ERR,FILE = 'error')
            WRITE(ERR,'(A,4X,F14.5)') 'ERROR INPUT MATRIX', EIGENVAL(I)
            WRITE(6,'(A)') 'ERROR, SEE FILE error'
            STOP
        ENDIF
        D(I,I) = SQRT(EIGENVAL(I))
    ENDDO

    ! --- STEP 2 : QR DECOMPOSITION --- !
    CALL PRODMAT(DIM,EIGENVEC,D,PT)
    CALL TRANSPOSE(DIM,PT,P)
    CALL QRHOUSE(DIM,P,DELTA,LT,.TRUE.)

    ! --- STEP 3 : CHOLESKY REDUCTION --- !
    CALL TRANSPOSE(DIM,LT,L)
    DO I = 1,DIM
        IF (L(I,I) < EPS0) THEN
            DO J = 1,DIM
                L(J,I) = -L(J,I)
            ENDDO
        ENDIF
    ENDDO

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