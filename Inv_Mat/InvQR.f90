SUBROUTINE INVMATQR(DIM,MATRIX,INVMATRIX)
    ! ------------------------------------------------------------------------------------------------------------- !
    ! --- THIS SUBROUTINE INVERSE THE MATRIX USING QR HOUSEHOLDER FACTORIZATION ----------------------------------- !
    ! --- housestep.f90 AND QRhouse.f90 AND operation_mat.f90 ARE USEFULL FOR THIS ALGORTIHME AND MUS BE COMPIL --- !
    ! --- THEY CAN BE FOUND ON Quentin-Marecat GITHUB PAGE -------------------------------------------------------- !
    ! ------------------------------------------------------------------------------------------------------------- !
    IMPLICIT NONE
    INTEGER :: DIM
    LOGICAL :: TEST
    REAL*8, PARAMETER :: EPS0 = 1.D-14
    REAL*8,INTENT(IN) :: MATRIX(DIM,DIM)
    REAL*8,INTENT(OUT) :: INVMATRIX(DIM,DIM)
    REAL*8 :: R(DIM,DIM), Q(DIM,DIM),RINV(DIM,DIM)
    REAL*8 :: QT(DIM,DIM), ID(DIM,DIM),ID2(DIM,DIM)
    INTEGER :: I,J,K,ERR,STAT

    ID = 0.
    DO I = 1,DIM
        ID(I,I) = 1
    ENDDO

    ERR = 97
    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) CLOSE(ERR,STATUS = 'delete')
    

    CALL QRHOUSE(DIM,MATRIX,Q,R)

    TEST = .FALSE.
    DO I = 1,DIM
        IF(ABS(R(I,I)) < EPS0) TEST = .TRUE.
    ENDDO
    IF (TEST) THEN
        OPEN(UNIT = ERR,FILE = 'error')
        WRITE(6,'(A)') 'NON INVERSIBLE MATRIX'
        WRITE(ERR,'(A)') 'NON INVERSIBLE MATRIX'
        WRITE(ERR,'(A)') 'R-MATRIX :'
        DO I = 1,DIM
            WRITE(ERR,'(100F14.5)') (R(I,J), J = 1,DIM)
        ENDDO
        WRITE(ERR,'(A)') '********************'
        STOP

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
    TEST = .FALSE.
    DO I = 1,DIM
        DO J = 1,DIM
            IF (ABS(ID2(I,J)-ID(I,J)) > EPS0) TEST = .TRUE.
        ENDDO
    ENDDO
    IF (TEST) THEN
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
    IF (STAT == 0) WRITE(6,'(A)') 'ERROR INVERSE MATRIX, SEE FILE error'

END SUBROUTINE
        
