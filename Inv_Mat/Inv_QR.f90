!!! --- 3 SUBROUTINES : --------------------- !!!
!!! --- INV_QR (ANY MATRICES) --------------- !!!
!!! --- INV_QR_SYM (SYMMETRIC MATRICES) ----- !!!
!!! --- INV_QR_TRI (TRIDIAGONAL MATRICES) --- !!!


SUBROUTINE INV_QR(DIM,MATRIX,INVMATRIX,CHECK)
    ! -------------------------------------------------------------------------------------------- !
    ! --- THIS SUBROUTINE INVERSE THE MATRIX USING QR HOUSEHOLDER FACTORIZATION ------------------ !
    ! ---  QRhouse.f90 AND operation_mat.f90 ARE USEFULL FOR THIS ALGORTIHME AND MUS BE COMPIL --- !
    ! --- THEY CAN BE FOUND ON Quentin-Marecat GITHUB PAGE --------------------------------------- !
    ! -------------------------------------------------------------------------------------------- !
    IMPLICIT NONE
    INTEGER :: DIM
    LOGICAL :: TEST
    REAL*8, PARAMETER :: EPS0 = 1.D-15, EPS = 1.D-12
    LOGICAL, INTENT(IN) :: CHECK
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
    

    CALL QRHOUSE(DIM,MATRIX,Q,R,.FALSE.)

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
    ENDIF

    ! --- VERIFICATION --- !
    IF (CHECK) THEN
        CALL PRODMAT(DIM,MATRIX,INVMATRIX,ID2)
        TEST = .FALSE.
        DO I = 1,DIM
            DO J = 1,DIM
                IF (ABS(ID2(I,J)-ID(I,J)) > EPS) TEST = .TRUE.
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
            WRITE(ERR,'(A)') 'M*M-1 - Id = '
            DO I = 1,DIM
                WRITE(ERR,'(100ES14.5)') (ID2(I,J)-ID(I,J), J = 1,DIM)
            ENDDO
            WRITE(ERR,'(A)') '************'
        ENDIF
    ENDIF

    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) WRITE(6,'(A)') 'ERROR INVERSE MATRIX, SEE FILE error'

END SUBROUTINE









SUBROUTINE INV_QR_SYM(DIM,MATRIX,INVMATRIX,CHECK)
    ! -------------------------------------------------------------------------------------------- !
    ! --- THIS SUBROUTINE INVERSE THE MATRIX USING QR HOUSEHOLDER FACTORIZATION ------------------ !
    ! --- THE MATRIX MUST BE SYMMETRIC ----------------------------------------------------------- !
    ! ---  QRhouse.f90, tridiag_house.f90 AND operation_mat.f90 ARE USEFULL FOR THIS ALGORTIHME AND MUS BE COMPIL --- !
    ! --- THEY CAN BE FOUND ON Quentin-Marecat GITHUB PAGE --------------------------------------- !
    ! -------------------------------------------------------------------------------------------- !
    IMPLICIT NONE
    INTEGER :: DIM
    LOGICAL :: TEST
    REAL*8, PARAMETER :: EPS0 = 1.D-15, EPS = 1.D-12
    LOGICAL, INTENT(IN) :: CHECK
    REAL*8,INTENT(IN) :: MATRIX(DIM,DIM)
    REAL*8,INTENT(OUT) :: INVMATRIX(DIM,DIM)
    REAL*8 :: R(DIM,DIM), Q(DIM,DIM),RINV(DIM,DIM)
    REAL*8 :: QT(DIM,DIM), ID(DIM,DIM),ID2(DIM,DIM)
    REAL*8 :: TRIDIAG(DIM,DIM), ROTMAT(DIM,DIM), ROTMATT(DIM,DIM)
    INTEGER :: I,J,K,ERR,STAT

    ID = 0.
    DO I = 1,DIM
        ID(I,I) = 1
    ENDDO

    ERR = 97
    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) CLOSE(ERR,STATUS = 'delete')
    
    CALL TRIDIAG_HOUSE(DIM,MATRIX,TRIDIAG,ROTMAT,.FALSE.)

    CALL INV_QR_TRI(DIM,TRIDIAG,INVMATRIX,.FALSE.)

    CALL PRODMAT(DIM,ROTMAT,INVMATRIX,INVMATRIX)
    CALL TRANSPOSE(DIM,ROTMAT,ROTMATT)
    CALL PRODMAT(DIM,INVMATRIX,ROTMATT,INVMATRIX)

    ! --- VERIFICATION --- !
    IF (CHECK) THEN
        CALL PRODMAT(DIM,MATRIX,INVMATRIX,ID2)
        TEST = .FALSE.
        DO I = 1,DIM
            DO J = 1,DIM
                IF (ABS(ID2(I,J)-ID(I,J)) > EPS) TEST = .TRUE.
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
            WRITE(ERR,'(A)') 'M*M-1 - Id = '
            DO I = 1,DIM
                WRITE(ERR,'(100ES14.5)') (ID2(I,J)-ID(I,J), J = 1,DIM)
            ENDDO
            WRITE(ERR,'(A)') '************'
        ENDIF
    ENDIF

    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) WRITE(6,'(A)') 'ERROR INVERSE MATRIX, SEE FILE error'

END SUBROUTINE












SUBROUTINE INV_QR_TRI(DIM,MATRIX,INVMATRIX,CHECK)
    ! -------------------------------------------------------------------------------------------- !
    ! --- THIS SUBROUTINE INVERSE THE MATRIX USING QR HOUSEHOLDER FACTORIZATION ------------------ !
    ! --- THE MATRIX MUST BE TRIDIAGONAL --------------------------------------------------------- !
    ! ---  QRhouse.f90 AND operation_mat.f90 ARE USEFULL FOR THIS ALGORTIHME AND MUS BE COMPIL --- !
    ! --- THEY CAN BE FOUND ON Quentin-Marecat GITHUB PAGE --------------------------------------- !
    ! -------------------------------------------------------------------------------------------- !
    IMPLICIT NONE
    INTEGER :: DIM
    LOGICAL :: TEST
    REAL*8, PARAMETER :: EPS0 = 1.D-15, EPS = 1.D-12
    LOGICAL, INTENT(IN) :: CHECK
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
    

    CALL QRHOUSE_TRI(DIM,MATRIX,Q,R,.FALSE.)

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
    ENDIF

    ! --- VERIFICATION --- !
    IF (CHECK) THEN
        CALL PRODMAT(DIM,MATRIX,INVMATRIX,ID2)
        TEST = .FALSE.
        DO I = 1,DIM
            DO J = 1,DIM
                IF (ABS(ID2(I,J)-ID(I,J)) > EPS) TEST = .TRUE.
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
            WRITE(ERR,'(A)') 'M*M-1 - Id = '
            DO I = 1,DIM
                WRITE(ERR,'(100ES14.5)') (ID2(I,J)-ID(I,J), J = 1,DIM)
            ENDDO
            WRITE(ERR,'(A)') '************'
        ENDIF
    ENDIF

    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) WRITE(6,'(A)') 'ERROR INVERSE MATRIX, SEE FILE error'

END SUBROUTINE
        
