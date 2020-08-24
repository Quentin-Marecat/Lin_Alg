SUBROUTINE CHOLESKY(DIM,MATRIX,L,CHECK)
    ! ---------------------------------------------------------------------------------------- !
    ! --- FOR ANY REPORT OR SUGGESTION, PLEASE CONTACT quentin.marecat@etu.umontpellier.fr --- !
    ! --- THIS SUBROUTINE DO THE CHOLESKY DECOMPOSITION OF A POSITIVE SYMMETRIC MATRIX ------- !
    ! ---------------------------------------------------------------------------------------- !
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: CHECK
    REAL*8, PARAMETER :: EPS0 = 1.D-15, EPS = 1.D-9
    LOGICAL :: TEST
    INTEGER, INTENT(IN) :: DIM
    REAL*8, INTENT(IN) :: MATRIX(DIM,DIM)
    REAL*8, INTENT(OUT) :: L(DIM,DIM)
    REAL*8 :: FIRSTELEM, LT(DIM,DIM), MATTEST(DIM,DIM)
    INTEGER :: I,J,K,ERR,STAT

    ERR = 97
    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) CLOSE(ERR,STATUS = 'delete')
    
    IF (MATRIX(1,1) <= EPS0 ) THEN
        OPEN(UNIT = ERR,FILE = 'error')
        WRITE(ERR,'(A,4X,F14.5)') 'ERROR INPUT MATRIX', MATRIX(1,1)
        WRITE(6,'(A)') 'ERROR, SEE FILE error'
        STOP
    ENDIF
    L(1,1) = SQRT(MATRIX(1,1))
    DO I = 2,DIM
        L(I,1) = MATRIX(1,I)/L(1,1)
    ENDDO
    DO I = 2,DIM
        FIRSTELEM = MATRIX(I,I)
        DO J = 1,I-1
            FIRSTELEM = FIRSTELEM - L(I,J)**2
        ENDDO
        IF (FIRSTELEM <= EPS0 ) THEN
            OPEN(UNIT = ERR,FILE = 'error')
            WRITE(ERR,'(A,4X,F14.5)') 'ERROR INPUT MATRIX', FIRSTELEM
            WRITE(6,'(A)') 'ERROR, SEE FILE error'
            STOP
        ENDIF
        L(I,I) = SQRT(FIRSTELEM)
        IF (I /= DIM) THEN
            DO J = I+1,DIM
                L(J,I) = MATRIX(I,J)
                DO K = 1,I-1
                    L(J,I) = L(J,I) -  L(I,K)*L(J,K)
                ENDDO
                L(J,I) = L(J,I)/L(I,I)
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