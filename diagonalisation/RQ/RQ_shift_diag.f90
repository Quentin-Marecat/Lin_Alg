SUBROUTINE RQ_shift_DIAG(DIM,MAT,EIGENVAL,EIGENVEC,SYM)
    ! --- FOR ANY REPORT OR SUGGESTION, PLEASE CONTACT quentin.marecat@etu.umontpellier.fr --- !
    ! --------------------------------------------------------------------------------- !
    ! --- THIS SUBROUTINE DIAGONALIZE THE MATRIX MAT ---------------------------------- !
    ! --- USING RQ ITERATIVE ALGORITHME WITH SHIFTING IMPROVEMENT --------------------- !
    ! --- NECESSARY IF IT INVOLVE SOME DEGENERANCIES WITH EIGENVALUES ----------------- !
    ! --- SYM = .TRUE. IF THE MATRIX IS SYMMETRIC ------------------------------------- !
    ! --- operation_mat.f90 tridiag.f90 QRhouse.f90 and housestep.f90 ARE NECESSARY --- !
    ! --------------------------------------------------------------------------------- !
    IMPLICIT NONE
    REAL*8, PARAMETER :: CONV = 1.D-12, EPS0 = 1.D-14
    INTEGER, PARAMETER :: MAXSTEP = 1D3
    LOGICAL :: TEST
    LOGICAL,INTENT(IN) :: SYM
    INTEGER, INTENT(IN) :: DIM
    REAL*8, INTENT(IN) :: MAT(DIM,DIM)
    REAL*8, INTENT(OUT) :: EIGENVEC(DIM,DIM), EIGENVAL(DIM)
    REAL*8 :: TRISUP(DIM,DIM),MATROT(DIM,DIM),MATTRIDIAG(DIM,DIM)
    REAL*8 :: ROTTRISUP(DIM,DIM), ID(DIM,DIM),MATTAMP(DIM,DIM)
    REAL*8 :: RESIDUAL, VEC(DIM), MATRES(DIM,DIM),MATROTT(DIM,DIM)
    REAL*8 ::  SHIFT, EIGENVECT(DIM,DIM)
    INTEGER COMPT, VALPROPRE, STAT, ERR, I, J, K, VALSHIFT

    TEST = .FALSE.
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
    MATROT = ID
    IF (SYM) CALL TRIDIAG_QR(DIM,MAT,MATTRIDIAG,MATROT)
    CALL TRANSPOSE(DIM,MATROT,MATROTT)
    CALL PRODMAT(DIM,MATROTT,MAT,MATRES)
    CALL PRODMAT(DIM,MATRES,MATROT)
    ! --- STEP 2 : RQ ITERATION WIT SHIFT UNTIL CONVERGENCE --- !
    EIGENVEC = ID
    VALSHIFT = DIM + 1
    SHIFT = 0
    DO VALPROPRE = 1,DIM
        RESIDUAL = 1 
        COMPT = 0
        DO WHILE(RESIDUAL > CONV .AND. COMPT < MAXSTEP) 
            IF (COMPT > 0 .OR. VALPROPRE /= 1) THEN
                VALSHIFT = VALSHIFT - 1
                IF (VALSHIFT < 1) VALSHIFT = DIM
                SHIFT = MATTRIDIAG(VALSHIFT,VALSHIFT)
            ENDIF
            DO I = 1,DIM
                 MATTRIDIAG(I,I) = MATTRIDIAG(I,I) - SHIFT
            ENDDO 
            CALL QRHOUSE(DIM,MATTRIDIAG,ROTTRISUP,TRISUP)
            CALL PRODMAT(DIM,TRISUP,ROTTRISUP,MATTRIDIAG)
            DO I = 1,DIM
                MATTRIDIAG(I,I) = MATTRIDIAG(I,I) + SHIFT
            ENDDO 
            CALL PRODMAT(DIM,EIGENVEC,ROTTRISUP,EIGENVEC)
            DO I = 1,DIM
                EIGENVAL(I) = MATTRIDIAG(I,I)
            ENDDO
            DO I = 1,DIM
                VEC(I) = EIGENVEC(I,VALPROPRE)
            ENDDO
            CALL MATAPPLI(DIM,MATRES,VEC,VEC)
            RESIDUAL = 0
            DO I = 1,DIM
                VEC(I) = VEC(I) - EIGENVAL(VALPROPRE)*EIGENVEC(I,VALPROPRE)
                RESIDUAL = RESIDUAL + VEC(I)**2
            ENDDO
            RESIDUAL = SQRT(RESIDUAL)
            COMPT = COMPT + 1
        ENDDO
        IF (COMPT == MAXSTEP) THEN
            TEST = .TRUE.
            OPEN(UNIT = ERR, FILE = 'error')
            WRITE(ERR,'(A,4X,E14.5)') 'EIGENVAL =',EIGENVAL(VALPROPRE)
            WRITE(ERR,'(A,4X,100E14.5)') 'EIGENVECTOR =',(EIGENVEC(I,VALPROPRE),I=1,DIM)
            WRITE(ERR,'(A,4X,I5,4X,A,E14.5)') 'MAX STEP REACHS =',MAXSTEP,'RESIDUAL =',CONV
            WRITE(ERR,'(A)') '****************'
        ENDIF
    ENDDO
    ! --- STEP 3 : PRODUCT WITH MATRIX ROTATION OF TRIDIAGONALE TRANFOR --- !
    IF (SYM) CALL PRODMAT(DIM,MATROT,EIGENVEC,EIGENVEC)
    ! --- STEP 4 : ORDERING EIGENVALS/EIGENVECTS --- !
    CALL ORDERING(DIM,EIGENVAL,EIGENVEC)

    ! --- VERIFICATION --- !
    MATTAMP = 0.
    DO I = 1,DIM
        MATTAMP(I,I) = EIGENVAL(I)
    ENDDO
    CALL PRODMAT(DIM,EIGENVEC,MATTAMP,MATTAMP)
    CALL TRANSPOSE(DIM,EIGENVEC,EIGENVECT)
    CALL PRODMAT(DIM,MATTAMP,EIGENVECT,MATTAMP)
    DO I = 1,DIM
        DO J = 1,DIM
            IF(ABS(MATTAMP(I,J)-MAT(I,J)) > CONV) TEST = .TRUE.
        ENDDO
    ENDDO
    IF (TEST) THEN
        OPEN(UNIT = ERR, FILE = 'error')
        WRITE(ERR,'(A,10X,A,4X,ES14.1)') 'PROBLEM DIAGONALISATION','CRIT CONV =',CONV
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


