SUBROUTINE RQ_DIAG(DIM,MAT,EIGENVAL,EIGENVEC,SYM,EGNVEC,CHECK)
    ! --- FOR ANY REPORT OR SUGGESTION, PLEASE CONTACT quentin.marecat@etu.umontpellier.fr --- !
    ! --------------------------------------------------------------------------------- !
    ! --- THIS SUBROUTINE DIAGONALIZE THE MATRIX MAT ---------------------------------- !
    ! --- USING RQ ITERATIVE ALGORITHME ----------------------------------------------- !
    ! --- SYM = .TRUE. IF THE MATRIX IS SYMMETRIC ------------------------------------- !
    ! --- operation_mat.f90 tridiag.f90 QRhouse.f90 and housestep.f90 ARE NECESSARY --- !
        ! --- CHECK = .TRUE. IF VERIFICATIONS HAS TO BE DONE -------------------------- !
    ! --- EGNVEC = .TRUE. IF EIGENVECTORS MUST BE COMPUTE ----------------------------- !
    ! --- IF CHECK = .TRUE., EGNVEC = .TRUE. ------------------------------------------ !
    ! --------------------------------------------------------------------------------- !
    IMPLICIT NONE
    REAL*8, PARAMETER :: EPS0 = 1.D-15, CONV = 1.D-12
    INTEGER, PARAMETER :: MAXSTEP = 1D4
    LOGICAL :: TEST
    LOGICAL,INTENT(IN) :: CHECK, EGNVEC, SYM
    INTEGER, INTENT(IN) :: DIM
    REAL*8, INTENT(IN) :: MAT(DIM,DIM)
    REAL*8, INTENT(OUT) :: EIGENVEC(DIM,DIM), EIGENVAL(DIM)
    REAL*8 :: TRISUP(DIM,DIM),MATROT(DIM,DIM),MATTRI(DIM,DIM)
    REAL*8 :: ROTTRISUP(DIM,DIM), ID(DIM,DIM),MATTAMP(DIM,DIM)
    REAL*8 :: RES,RES2, VEC(DIM),NOR(DIM), EIGENVECT(DIM,DIM)
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
    IF (SYM) CALL TRIDIAG_HOUSE(DIM,MAT,MATTRI,MATROT,.FALSE.)
    MATTAMP = MATTRI
    ! --- STEP 2 : RQ ITERATION WITH SHIFT UNTIL CONVERGENCE --- !
    EIGENVEC = ID
    RES = 1 
    COMPT = 0
    DO WHILE(RES > CONV .AND. COMPT < MAXSTEP) 
        IF (SYM) THEN
            CALL QRHOUSE_TRI(DIM,MATTRI,ROTTRISUP,TRISUP,.FALSE.)
            CALL PRODMAT_RQ(DIM,TRISUP,ROTTRISUP,MATTRI)
        ELSE
            CALL QRHOUSE(DIM,MATTRI,ROTTRISUP,TRISUP,.FALSE.)
            CALL PRODMAT(DIM,TRISUP,ROTTRISUP,MATTRI)
        ENDIF
        IF (EGNVEC) CALL PRODMAT(DIM,EIGENVEC,ROTTRISUP,EIGENVEC)
        ! --- 0 ELEMENT CALCULATION --- !
        RES = 0
        DO I = 2,DIM
            DO J = 1,I-1
                RES = RES + MATTRI(I,J)**2
            ENDDO
        ENDDO
        RES = SQRT(2*RES)
        ! ----------------------------- !
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
        IF (TEST) THEN
            OPEN(UNIT = ERR, FILE = 'error')
            WRITE(ERR,'(A,10X,A,4X,I14.5,10X,A,4X,I14.5)') 'PROBLEM DIAGONALISATION'&
            & ,'MAXSTEP =',MAXSTEP,'COMPT',COMPT
            WRITE(ERR,'(A)') '***********************'
            WRITE(ERR,'(A,4X,ES14.5,4X,A,4X,ES14.5)')'CRIT CONV =',CONV,'RESIDUAL =',RES
            WRITE(ERR,'(A)') '***********************'
            WRITE(ERR,'(A)') 'DIAGONAL MATRIX OBTAINED'
            DO I = 1,DIM
                WRITE(ERR,'(100F14.5)') (MATTRI(I,J),J=1,DIM)
            ENDDO
            WRITE(ERR,'(A)') '***********************'
            MATTRI = ID
            DO I = 1,DIM
                MATTRI(I,I) = EIGENVAL(I)
            ENDDO
            CALL PRODMAT(DIM,EIGENVEC,MATTRI,MATTAMP)
            CALL TRANSPOSE(DIM,EIGENVEC,EIGENVECT)
            CALL PRODMAT(DIM,MATTAMP,EIGENVECT,MATTAMP)
            WRITE(ERR,'(A)') ' Q*D*Q(T) - MAT'
            DO I = 1,DIM
                WRITE(ERR,'(100ES14.5)') (MATTAMP(I,J)-MAT(I,J), J=1,DIM)
            ENDDO
            WRITE(ERR,'(A)') '***********************'
            CALL PRODMAT(DIM,MAT,EIGENVEC,MATTAMP)
            DO I = 1,DIM
                DO J = 1,DIM
                    VEC(J) = MATTAMP(J,I) - EIGENVAL(I)*EIGENVEC(J,I)
                ENDDO
                CALL NORMVEC(DIM,VEC,NOR(I))
            ENDDO
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
