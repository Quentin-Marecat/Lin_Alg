SUBROUTINE LU(DIM,MATRIX,L,U)
    ! --- FOR ANY REPORT OR SUGGESTION, PLEASE CONTACT quentin.marecat@etu.umontpellier.fr --- !
    ! ------------------------------------------------------------- !
    ! --- THIS ALGORITHME DO THE LU FACTORIZATION OF THE MATRIX --- !
    ! ------------------------------------------------------------- !
    REAL*8, PARAMETER :: EPS0 = 1.D-12
    INTEGER, INTENT(IN) :: DIM
    REAL*8, INTENT(IN) :: MATRIX(DIM,DIM)
    REAL*8, INTENT(OUT) :: L(DIM,DIM), U(DIM,DIM)
    REAL*8 :: ID(DIM,DIM), TAMP, LHIGHDIM(2*DIM,DIM,DIM), TAMPMAT(DIM,DIM) ! --- LHIGHDIM = ARRAY OF ELEMENTARY MATRIX --- !
    LOGICAL :: TEST
    INTEGER :: FIRSTNONZERO, COMPT, ERR, STAT
    INTEGER :: I,J,K
    ERR = 97
    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) CLOSE(ERR,STATUS = 'delete')

    U = MATRIX
    ID = 0.
    DO I = 1,DIM
        ID(I,I) = 1.
    ENDDO
    LHIGHDIM = 0
    L = ID

    COMPT = 0
    DO I = 1,DIM-1

        ! --- STEP 1 : EDIT THE SHAPE OF THE MATRIX --- !
        TEST = .FALSE.
        DO J = I,DIM
            IF (.NOT. TEST) FIRSTNONZERO = J
            IF (ABS(U(J,I)) > EPS0) TEST = .TRUE.
        ENDDO
        IF (TEST) THEN ! --- THE (I,I) ELEMENT IS EQUAL TO 0 --- ! 
            IF (FIRSTNONZERO /= I) THEN ! --- LINE PERMUTATION --- !
                LHIGHDIM(COMPT + 1,FIRSTNONZERO,FIRSTNONZERO) = 0
                LHIGHDIM(COMPT + 1,I,I) = 0
                LHIGHDIM(COMPT + 1,I,FIRSTNONZERO) = 1
                LHIGHDIM(COMPT + 1,FIRSTNONZERO,I) = 1
                COMPT = COMPT + 1    
            ENDIF
            DO J = 1,DIM  ! --- LINE PERMUTATION --- !
                TAMP = U(I,J)
                U(I,J) = U(FIRSTNONZERO,J)
                U(FIRSTNONZERO,J) = TAMP 
            ENDDO       

            DO J = 1,DIM
                LHIGHDIM(COMPT + 1,J,J) = 1 
            ENDDO   

            DO J = I+1,DIM
                TAMP = U(J,I)
                DO K = I,DIM
                    U(J,K) = U(J,K) - U(I,K)*TAMP/U(I,I) 
                ENDDO
                LHIGHDIM(COMPT + 1,J,I) = TAMP/U(I,I)
            ENDDO
            COMPT = COMPT + 1
        ENDIF
    ENDDO

    DO I = 1,COMPT
        DO J = 1,DIM
            DO K = 1,DIM
                TAMPMAT(J,K) = LHIGHDIM(COMPT - I + 1,J,K)
            ENDDO
        ENDDO
        CALL PRODMAT(DIM,TAMPMAT,L,L)
    ENDDO
    ! --- VERIFICATION --- !
    CALL PRODMAT(DIM,L,U,TAMPMAT)
    TEST = .FALSE.
    DO I = 1,DIM
        DO J = 1,DIM
            IF (ABS(TAMPMAT(I,J)-MATRIX(I,J)) > EPS0) TEST = .TRUE.
        ENDDO
    ENDDO 
    IF (TEST) THEN
        OPEN(UNIT = ERR,FILE = 'error')
        WRITE(ERR,*) 'ERROR LU FACTORIZATION'
        DO I = 1,DIM
            WRITE(ERR,'(100F14.5)') (L(I,J),J = 1,DIM)
        ENDDO
        WRITE(ERR,*) '*******************'
        DO I = 1,DIM
            WRITE(ERR,'(100F14.5)') (U(I,J),J = 1,DIM)
        ENDDO
    ENDIF
    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) WRITE(6,'(A)') 'ERROR LU DECOMPO, SEE FILE error'
END SUBROUTINE