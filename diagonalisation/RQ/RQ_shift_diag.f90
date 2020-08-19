SUBROUTINE RQ_SHIFT_DIAG(DIM,MAT,EIGENVAL,EIGENVEC,SYM,EGNVEC,CHECK)
    ! --- FOR ANY REPORT OR SUGGESTION, PLEASE CONTACT quentin.marecat@etu.umontpellier.fr --- !
    ! --------------------------------------------------------------------------------- !
    ! --- THIS SUBROUTINE DIAGONALIZE THE MATRIX MAT ---------------------------------- !
    ! --- USING RQ ITERATIVE ALGORITHME WITH SHIFTING IMPROVEMENT --------------------- !
    ! --- NECESSARY IF IT INVOLVE SOME DEGENERANCIES WITH EIGENVALUES ----------------- !
    ! --- SYM = .TRUE. IF THE MATRIX IS SYMMETRIC ------------------------------------- !
    ! --- operation_mat.f90 tridiag.f90 QRhouse.f90 and housestep.f90 ARE NECESSARY --- !
    ! --- CHECK = .TRUE. IF VERIFICATIONS HAS TO BE DONE ------------------------------ !
    ! --- EGNVEC = .TRUE. IF EIGENVECTORS MUST BE COMPUTE ----------------------------- !
    ! --- IF CHECK = .TRUE., EGNVEC = .TRUE. ------------------------------------------ !
    ! --------------------------------------------------------------------------------- !
    IMPLICIT NONE
    REAL*8, PARAMETER :: CONV = 1.D-6, EPS0 = 1.D-15
    INTEGER, PARAMETER :: MAXSTEP = 1D4
    LOGICAL,INTENT(IN) :: CHECK, EGNVEC
    LOGICAL :: TEST
    LOGICAL,INTENT(IN) :: SYM
    INTEGER, INTENT(IN) :: DIM
    REAL*8, INTENT(IN) :: MAT(DIM,DIM)
    REAL*8, INTENT(OUT) :: EIGENVEC(DIM,DIM), EIGENVAL(DIM)
    REAL*8 :: TRISUP(DIM,DIM),MATROT(DIM,DIM),MATTRI(DIM,DIM)
    REAL*8 :: ROTTRISUP(DIM,DIM), ID(DIM,DIM),MATTAMP(DIM,DIM)
    REAL*8 :: RES,RES2, VEC(DIM)
    REAL*8 ::  SHIFT, NOR(DIM)
    INTEGER COMPT, STAT, ERR, I, J, K, VALSHIFT

    TEST = .FALSE.
    ERR = 97
    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) CLOSE(ERR,STATUS = 'delete')

    ID = 0.
    DO I = 1,DIM
        ID(I,I) = 1.
    ENDDO
    ! --- STEP 1 : TRIDIAGONAL DECOMPO OF THE MATRIX --- !
    MATROT = ID
    MATTRI = MAT
    IF (SYM) CALL TRIDIAG_QR(DIM,MAT,MATTRI,MATROT,.TRUE.)
    MATTAMP = MATTRI
    ! --- STEP 2 : RQ ITERATION WITH SHIFT UNTIL CONVERGENCE --- !
    EIGENVEC = ID
    VALSHIFT = DIM + 1
    SHIFT = 0
    RES = 1 
    COMPT = 0
    DO WHILE(RES > CONV .AND. COMPT < MAXSTEP) 
        DO I = 1,DIM
            MATTRI(I,I) = MATTRI(I,I) - SHIFT
        ENDDO 
        CALL QRHOUSE_TRI(DIM,MATTRI,ROTTRISUP,TRISUP,.FALSE.)
        CALL PRODMAT_RQ(DIM,TRISUP,ROTTRISUP,MATTRI)
        ! --- SHIFT METHOD --- !
        DO I = 1,DIM
            MATTRI(I,I) = MATTRI(I,I) + SHIFT
        ENDDO 
        IF (EGNVEC) CALL PRODMAT(DIM,EIGENVEC,ROTTRISUP,EIGENVEC)
        ! --- 0 ELEMENT CALCULATION --- !
        RES = 0
        DO I = 1,DIM
            DO J = 1,DIM
                VEC(J) = EIGENVEC(J,I)
            ENDDO
            CALL MATAPPLI_TRI(DIM,MATTAMP,VEC,VEC)
            DO J = 1,DIM
                VEC(J) = VEC(J) - MATTRI(I,I)*EIGENVEC(J,I)
            ENDDO
            CALL NORMVEC(DIM,VEC,RES2)
            IF (RES2 > RES) RES = RES2
        ENDDO
        ! --- SHIFT PROCESS --- !
        VALSHIFT = VALSHIFT - 1
        IF (VALSHIFT < 1) VALSHIFT = DIM
        SHIFT = MATTRI(VALSHIFT,VALSHIFT)
        COMPT = COMPT + 1
    ENDDO
    DO I = 1,DIM
        EIGENVAL(I) = MATTRI(I,I)
    ENDDO
    ! --- STEP 3 : PRODUCT WITH MATRIX ROTATION OF TRIDIAG TRANFORMATION --- !
    IF (SYM) CALL PRODMAT(DIM,MATROT,EIGENVEC,EIGENVEC)
    ! --- STEP 4 : ORDERING EIGENVALS/EIGENVECTS --- !
    CALL ORDERING(DIM,EIGENVAL,EIGENVEC)

    ! --- VERIFICATION --- !
    IF (CHECK) THEN
        TEST = .FALSE.
        IF (COMPT == MAXSTEP) TEST = .TRUE.
        CALL PRODMAT(DIM,MAT,EIGENVEC,MATTAMP)
        DO I = 1,DIM
            DO J = 1,DIM
                VEC(J) = MATTAMP(J,I) - EIGENVAL(I)*EIGENVEC(J,I)
            ENDDO
            CALL NORMVEC(DIM,VEC,NOR(I))
        ENDDO
        IF (TEST) THEN
            OPEN(UNIT = ERR, FILE = 'error')
            WRITE(ERR,'(A,10X,A,4X,I14.5,10X,A,4X,I14.5)') 'PROBLEM DIAGONALISATION'&
            & ,'MAXSTEP =',MAXSTEP,'COMPT',COMPT
            WRITE(ERR,'(A)') '***********************'
            WRITE(ERR,'(A,4X,ES14.5,4X,A,4X,ES14.5)')'CRIT CONV =',CONV,'RESIDUAL =',RES
            WRITE(ERR,'(A)') '***********************'
            WRITE(ERR,'(A)') 'DIAGONAL MATRIX OBTAINED'
            DO I = 1,DIM
                WRITE(ERR,'(100ES14.5)') (MATTRI(I,J),J=1,DIM)
            ENDDO
            WRITE(ERR,'(A)') '***********************'
            WRITE(ERR,'(A)') 'MAT*EIGENVEC - EIGENVAL*EIGENVEC = 0'
            DO I = 1,DIM
                WRITE(ERR,'(100ES14.5)')(MATTAMP(J,I) - EIGENVAL(I)*EIGENVEC(J,I),J=1,DIM)
            ENDDO
            WRITE(ERR,'(A)') '***********************'
            WRITE(ERR,'(A)') 'NORME (MAT*EIGENVEC - EIGENVAL*EIGENVEC = 0)'
            WRITE(ERR,'(100ES14.5)') (NOR(I),I=1,DIM)
            WRITE(ERR,'(A)') '***********************'
        ENDIF
    ENDIF

    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) WRITE(6,'(A)') 'ERROR DIAGONALIZATION, SEE FILE error'

END SUBROUTINE


