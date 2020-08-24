SUBROUTINE CHOLESKY_LDL(DIM,MATRIX,L,CHECK)
    ! ------------------------------------------------------------------------------------ !
    ! --- THIS SUBROUTINE DO THE CHOLESKY DECOMPOSITION OF A POSITIVE SYMMETRIC MATRIX --- !
    ! --- USING A DIAGONAL MATRIX, LIMITING SQUARE ROOTS CALCULATION --------------------- !
    ! ------------------------------------------------------------------------------------ !
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: CHECK
    REAL*8, PARAMETER :: EPS0 = 1.D-15, EPS = 1.D-9
    LOGICAL :: TEST
    INTEGER, INTENT(IN) :: DIM
    REAL*8, INTENT(IN) :: MATRIX(DIM,DIM)
    REAL*8, INTENT(OUT) :: L(DIM,DIM)
    REAL*8 :: FIRSTELEM, LT(DIM,DIM), MATTEST(DIM,DIM)
    REAL*8 :: D(DIM,DIM), ID(DIM,DIM)
    INTEGER :: I,J,K,ERR, STAT

    ERR = 97
    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) CLOSE(ERR,STATUS = 'delete')
    
    ID = 0.
    DO I = 1,DIM
        ID(I,I) = 1
    ENDDO
    L = ID
    D = 0.

    DO I = 1,DIM
        D(I,I) = MATRIX(I,I)
        IF (I > 1) THEN
            DO J = 1,I-1
                D(I,I) = D(I,I) - D(J,J)*(L(I,J)**2)
            ENDDO
        ENDIF
        IF (D(I,I) <= EPS0 ) THEN
            OPEN(UNIT = ERR,FILE = 'error')
            WRITE(ERR,'(A,4X,F14.5)') 'ERROR INPUT MATRIX', D(I,I)
            WRITE(6,'(A)') 'ERROR, SEE FILE error'
            STOP
        ENDIF
        IF (I /= DIM) THEN
            DO J = I+1,DIM
                L(J,I) = MATRIX(J,I)
                DO K = 1,I-1
                    L(J,I) = L(J,I) -  L(I,K)*L(J,K)*D(K,K)
                ENDDO
                L(J,I) = L(J,I)/D(I,I)
            ENDDO
        ENDIF
    ENDDO
    DO I = 1,DIM
        D(I,I) = SQRT(D(I,I))
    ENDDO
    CALL PRODMAT(DIM,L,D,L)
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
        ENDIF
        WRITE(ERR,'(A)')'***************'
        WRITE(ERR,'(A,4X,ES14.5)')'EPS',EPS
        WRITE(ERR,'(A)') 'MAT - L*LT ='
        DO I = 1,DIM
            WRITE(ERR,'(100ES14.5)') (MATRIX(I,J) - MATTEST(I,J), J=1,DIM)
        ENDDO
        WRITE(ERR,'(A)')'***************'
    ENDIF

    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) WRITE(6,'(A)') 'ERROR, SEE FILE error'
END SUBROUTINE