SUBROUTINE DIAGMATSYM(DIM,MAT,EIGENVAL,EIGENVEC,SYM)
    ! --- FOR ANY REPORT OR SUGGESTION, PLEASE CONTACT quentin.marecat@etu.umontpellier.fr --- !
    ! --------------------------------------------------------------------------------- !
    ! --- THIS SUBROUTINE DIAGONALIZE THE MATRIX MAT ---------------------------------- !
    ! --- USING RQ ITERATIVE ALGORITHME ----------------------------------------------- !
    ! --- SYM = .TRUE. IF THE MATRIX IS SYMMETRIC ------------------------------------- !
    ! --- operation_mat.f90 tridiag.f90 QRhouse.f90 and housestep.f90 ARE NECESSARY --- !
    ! --------------------------------------------------------------------------------- !
    IMPLICIT NONE
    REAL*8, PARAMETER :: EPS = 1.D-10, EPS0 = 1.D-12
    INTEGER, PARAMETER :: MAXSTEP = 1D5
    LOGICAL :: TEST
    LOGICAL,INTENT(IN) :: SYM
    INTEGER, INTENT(IN) :: DIM
    REAL*8, INTENT(IN) :: MAT(DIM,DIM)
    REAL*8, INTENT(OUT) :: EIGENVEC(DIM,DIM), EIGENVAL(DIM)
    REAL*8 :: TRISUP(DIM,DIM),MATROT(DIM,DIM),MATTRIDIAG(DIM,DIM)
    REAL*8 :: ROTTRISUP(DIM,DIM), ID(DIM,DIM),MATTAMP(DIM,DIM)
    REAL*8 :: CONV
    INTEGER COMPT, VALPROPRE, STAT, ERR, I, J

    ERR = 97
    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) CLOSE(ERR,STATUS = 'delete')

    TRISUP = MAT
    MATTRIDIAG = MAT
    ID = 0.
    DO I = 1,DIM
        ID(I,I) = 1.
    ENDDO
    ! --- STEP 1 : TRIDIAGONALIZATION OF THE MATRIX --- !
    IF (SYM) CALL TRIDIAG(DIM,MAT,MATTRIDIAG,MATROT)
    ! --- STEP 2 : RQ ITERATION UNTIL CONVERGENCE --- !
    EIGENVEC = ID
    DO VALPROPRE = 1,DIM
        CONV = 1 
        COMPT = 0
        DO WHILE(CONV > EPS .AND. COMPT < MAXSTEP) 
            CALL QRHOUSE(DIM,MATTRIDIAG,ROTTRISUP,TRISUP)
            CALL PRODMAT(DIM,TRISUP,ROTTRISUP,MATTRIDIAG)
            CALL PRODMAT(DIM,EIGENVEC,ROTTRISUP,EIGENVEC)
            IF (COMPT > 0) CONV = ABS(EIGENVAL(VALPROPRE) - MATTRIDIAG(VALPROPRE,VALPROPRE))
            DO I = 1,DIM
                EIGENVAL(I) = MATTRIDIAG(I,I)
            ENDDO
            COMPT = COMPT + 1
        ! --- VERIFICATION --- !
            IF (SYM) THEN
                TEST = .FALSE.
                DO I = 1,DIM-2
                    DO J = I+2,DIM
                        IF (ABS(MATTRIDIAG(I,J)) > EPS .OR. ABS(MATTRIDIAG(J,I)) > EPS) TEST = .TRUE.
                    ENDDO
                ENDDO
                DO I = 1,DIM - 1
                    IF (ABS(MATTRIDIAG(I,I+1)-MATTRIDIAG(I+1,I)) > EPS0) TEST = .TRUE.
                ENDDO
                IF (TEST) THEN
                    OPEN(UNIT = ERR, FILE = 'error')
                    WRITE(ERR,'(A)')  'PROBLEM TRIDIAGONALISATION RQ ITERATIVE'
                    DO I = 1,DIM
                        WRITE(ERR,'(100F14.5)')(MATTRIDIAG(I,J),J=1,DIM)
                    ENDDO
                ENDIF
            ENDIF
        ENDDO
    ENDDO
    ! --- STEP 3 : PRODUCT WITH MATRIX ROTATION OF TRIDIAGONALE TRANFOR --- !
    IF (SYM) CALL PRODMAT(DIM,MATROT,EIGENVEC,EIGENVEC)
    ! --- STEP 4 : ORDERING EIGENVALS/EIGENVECTS --- !
    CALL ORDERING(DIM,EIGENVAL,EIGENVEC)
    ! --- VERIFICATION --- !
    TEST = .FALSE.
    IF (COMPT == MAXSTEP) THEN
        TEST = .TRUE.
        WRITE(ERR,'(A,4X,I5)') 'MAX STEP REACHS =',MAXSTEP
    ENDIF
    DO I = 1,DIM-1
        DO J = I+1,DIM
            IF (ABS(MATTRIDIAG(I,J)) > 100*EPS .OR. ABS(MATTRIDIAG(J,I)) > 100*EPS) TEST = .TRUE.
        ENDDO
    ENDDO
    MATTAMP = 0.
    DO I = 1,DIM
        MATTAMP(I,I) = EIGENVAL(I)
    ENDDO
    CALL PRODMAT(DIM,EIGENVEC,MATTAMP,MATTAMP)
    CALL TRANSPOSE(DIM,EIGENVEC,EIGENVEC)
    CALL PRODMAT(DIM,MATTAMP,EIGENVEC,MATTAMP)
    DO I = 1,DIM
        DO J = 1,DIM
            IF(ABS(MATTAMP(I,J)-MAT(I,J)) > 100*EPS) TEST = .TRUE.
        ENDDO
    ENDDO
    IF (TEST) THEN
        OPEN(UNIT = ERR, FILE = 'error')
        WRITE(ERR,'(A,10X,A,4X,ES14.1)') 'PROBLEM DIAGONALISATION','EPS =',EPS
        WRITE(ERR,'(A)') '***********************'
        WRITE(ERR,'(A)') 'DIAGONAL MATRIX OBTAINED'
        DO I = 1,DIM
            WRITE(ERR,*) (MATTRIDIAG(I,J),J=1,DIM)
        ENDDO
        WRITE(ERR,'(A)') '***********************'
        WRITE(ERR,'(A)') 'MAT_INI - P*DIAG*P-1'
        DO I = 1,DIM
            WRITE(ERR,*) (MAT(I,J) - MATTAMP(I,J),J=1,DIM)
        ENDDO
        WRITE(ERR,'(A)') '***********************'
    ENDIF

    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) WRITE(6,'(A)') 'ERROR DIAGONALIZATION, SEE FILE error'

END SUBROUTINE













SUBROUTINE ORDERING(DIM,EIGENVAL,EIGENVECT)
    ! -------------------------------------------------------------- !
    ! --- THIS SUBROUTINE ORDER THE EIGENVALUES AND EIGENVECTORS --- !
    ! -------------------------------------------------------------- !
    INTEGER :: DIM
    REAL*8 :: EIGENVAL(DIM),EIGENVECT(DIM,DIM)
    REAL*8 :: EIGENVALTAMP(DIM),EIGENVECTTAMP(DIM,DIM), MINI
    INTEGER :: ORD(DIM)
    INTEGER :: I,J
    EIGENVALTAMP = EIGENVAL
    EIGENVECTTAMP = EIGENVECT
    DO I = 1,DIM-1
        ORD(I) = I
        MINI = EIGENVALTAMP(I)
        DO J = I+1,DIM
            IF (EIGENVALTAMP(J) < MINI) THEN
                ORD(I) = J
                MINI = EIGENVALTAMP(J)
            ENDIF
        ENDDO
        EIGENVALTAMP(I) = MINI
        EIGENVALTAMP(ORD(I)) = EIGENVAL(I)
        DO J = 1,DIM
            EIGENVECTTAMP(J,I) = EIGENVECT(J,ORD(I))
            EIGENVECTTAMP(J,ORD(I)) = EIGENVECT(J,I)
        ENDDO
        EIGENVAL = EIGENVALTAMP
        EIGENVECT = EIGENVECTTAMP
    ENDDO
    RETURN
END SUBROUTINE