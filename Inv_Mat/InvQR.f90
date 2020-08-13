SUBROUTINE INVMATQR(DIME,MATRIX,INVMATRIX)
    USE INVMATMOD
    ! ----------------------------------------------------------------------------- !
    ! --- THIS SUBROUTINE INVERSE THE MATRIX USING QR HOUSEHOLDER FACTORIZATION --- !
    ! ----------------------------------------------------------------------------- !
    IMPLICIT NONE
    INTEGER :: DIME
    REAL*8,INTENT(IN) :: MATRIX(DIME,DIME)
    REAL*8,INTENT(OUT) :: INVMATRIX(DIME,DIME)
    REAL*8 :: R(DIME,DIME), Q(DIME,DIME),RINV(DIME,DIME)
    REAL*8 :: QT(DIME,DIME), ID(DIME,DIME),ID2(DIME,DIME)
    INTEGER :: STAT
    DIM = DIME
    ID = 0.
    DO I = 1,DIM
        ID(I,I) = 1
    ENDDO

    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) CLOSE(ERR,STATUS = 'delete')
    

    CALL HOUSEHOLDER(MATRIX,R,Q)

    TESTLOG = .FALSE.
    DO I = 1,DIM
        IF(ABS(R(I,I)) < EPS0) TESTLOG = .TRUE.
    ENDDO
    IF (TESTLOG) THEN
        OPEN(UNIT = ERR,FILE = 'error')
        WRITE(ERR,'(A)') 'NON INVERSIBLE MATRIX'
        WRITE(ERR,'(A)') 'R-MATRIX :'
        DO I = 1,DIM
            WRITE(ERR,'(100F14.5)') (R(I,J), J = 1,DIM)
        ENDDO
        WRITE(ERR,'(A)') '********************'

    ELSE
        RINV = 0.
        DO I = 1,DIM
            RINV(I,I) = 1/R(I,I)
        ENDDO

        DO I = DIM,2,-1
            DO J = I-1,1,-1
                DO K = J+1,DIM
                    RINV(J,I) = RINV(J,I) - (RINV(K,I)*R(J,K)*(1/R(J,J)))
                ENDDO
            ENDDO
        ENDDO

        CALL TRANSPOSE(DIM,Q,QT)

        CALL PRODMAT(DIM,RINV,QT,INVMATRIX)
        CALL PRODMAT(DIM,MATRIX,INVMATRIX,ID2)
    ENDIF

    ! --- VERIFICATION --- !
    TESTLOG = .FALSE.
    DO I = 1,DIM
        DO J = 1,DIM
            IF (ABS(ID2(I,J)-ID(I,J)) > EPS0) TESTLOG = .TRUE.
        ENDDO
    ENDDO
    IF (TESTLOG) THEN
        OPEN(UNIT = ERR,FILE = 'error')
        WRITE(ERR,'(A)') 'PROBLEM MATRIX INVERSION'
        WRITE(ERR,'(A)') 'R = '
        DO I = 1,DIM
            WRITE(ERR,'(100F14.5)') (R(I,J), J = 1,DIM)
        ENDDO
        WRITE(ERR,'(A)') '************'
        WRITE(ERR,'(A)') 'M-1 = '
        DO I = 1,DIM
            WRITE(ERR,'(100F14.5)') (INVMATRIX(I,J), J = 1,DIM)
        ENDDO
        WRITE(ERR,'(A)') '************'
    ENDIF

    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) WRITE(6,'(A)') 'ERROR, SEE FILE error'
END SUBROUTINE
        
