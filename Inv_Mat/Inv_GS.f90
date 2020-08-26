SUBROUTINE INV_GS(DIM,MATRIX,INVMATRIX,CHECK)
    ! ------------------------------------------------------------------------------------------------------------------ !
    ! --- THIS SUBROUTINE INVERSE THE MATRIX USING SOLUTION SUCCESIV LINEAR EQUATION SOLUTION WITH GAUSS_SEIGDEL ITERATION  ------------- !
    ! --- Jacobi.f90 AND operation_mat.f90 ARE USEFULL FOR THIS ALGORTIHME AND MUSt BE COMPIL --------------- !
    ! --- THEY CAN BE FOUND ON Quentin-Marecat GITHUB PAGE ------------------------------------------------------------- !
    ! ------------------------------------------------------------------------------------------------------------------ !
    IMPLICIT NONE
    REAL*8, PARAMETER :: EPS0 = 1.D-15, EPS = 1.D-12
    LOGICAL, INTENT(IN) :: CHECK
    INTEGER :: DIM
    REAL*8, INTENT(OUT) :: INVMATRIX(DIM,DIM)
    REAL*8 :: MATRIX(DIM,DIM), VEC(DIM),SOL(DIM), ID(DIM,DIM), ID2(DIM,DIM)
    LOGICAL :: TEST
    INTEGER :: I, J, STAT, ERR

    ERR = 97
    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) CLOSE(ERR,STATUS = 'delete')

    INVMATRIX = 0.
    ID = 0.
    DO I = 1,DIM
        ID(I,I) = 1
        VEC = 0
        VEC(I) = 1
        CALL GAUSS_SEIDEL(DIM,MATRIX,VEC,SOL,.FALSE.)
        DO J = 1,DIM
            INVMATRIX(J,I) = SOL(J)
        ENDDO
    ENDDO

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
            WRITE(ERR,'(A)') 'M-1 = '
            DO I = 1,DIM
                WRITE(ERR,'(100F14.5)') (INVMATRIX(I,J), J = 1,DIM)
            ENDDO
            WRITE(ERR,'(A)') '************'
            WRITE(ERR,'(A)') 'M*M-1 - ID = '
            DO I = 1,DIM
                WRITE(ERR,'(100ES14.5)') (ID2(I,J)-ID(I,J), J = 1,DIM)
            ENDDO
            WRITE(ERR,'(A)') '************'
        ENDIF
    ENDIF

    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) WRITE(6,'(A)') 'ERROR INVERSE MATRIX, SEE FILE error'

END SUBROUTINE

