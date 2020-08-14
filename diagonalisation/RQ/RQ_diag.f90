SUBROUTINE RQ_DIAG(DIM,MAT,EIGENVAL,EIGENVEC,SYM)
    ! --- FOR ANY REPORT OR SUGGESTION, PLEASE CONTACT quentin.marecat@etu.umontpellier.fr --- !
    ! --------------------------------------------------------------------------------- !
    ! --- THIS SUBROUTINE DIAGONALIZE THE MATRIX MAT ---------------------------------- !
    ! --- USING RQ ITERATIVE ALGORITHME ----------------------------------------------- !
    ! --- SYM = .TRUE. IF THE MATRIX IS SYMMETRIC ------------------------------------- !
    ! --- operation_mat.f90 tridiag.f90 QRhouse.f90 and housestep.f90 ARE NECESSARY --- !
    ! --------------------------------------------------------------------------------- !
    IMPLICIT NONE
    REAL*8, PARAMETER :: CONV = 1.D-10, EPS0 = 1.D-15
    INTEGER, PARAMETER :: MAXSTEP = 1D4
    LOGICAL :: TEST
    LOGICAL,INTENT(IN) :: SYM
    INTEGER, INTENT(IN) :: DIM
    REAL*8, INTENT(IN) :: MAT(DIM,DIM)
    REAL*8, INTENT(OUT) :: EIGENVEC(DIM,DIM), EIGENVAL(DIM)
    REAL*8 :: TRISUP(DIM,DIM),MATROT(DIM,DIM),MATTRIDIAG(DIM,DIM)
    REAL*8 :: ROTTRISUP(DIM,DIM), ID(DIM,DIM),MATTAMP(DIM,DIM)
    REAL*8 :: RESIDUAL, VEC(DIM), MATRES(DIM,DIM),MATROTT(DIM,DIM)
    REAL*8 :: EIGENVECT(DIM,DIM)
    INTEGER COMPT, VALPROPRE, STAT, ERR, I, J, K

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
    ! --- THE TRIDIAGONALISATION METHOD CAN BE CHANGE --- !(DIM,MAT,MATTRIDIAG,MATROT)

! --- STEP 2 : RQ ITERATION WITh SHIFT UNTIL CONVERGENCE --- !
    EIGENVEC = ID
    RESIDUAL = 1 
    COMPT = 0
    DO WHILE(RESIDUAL > CONV .AND. COMPT < MAXSTEP) 
        CALL QRHOUSE(DIM,MATTRIDIAG,ROTTRISUP,TRISUP)
        CALL PRODMAT(DIM,TRISUP,ROTTRISUP,MATTRIDIAG)
        CALL PRODMAT(DIM,EIGENVEC,ROTTRISUP,EIGENVEC)

        DO I = 1,DIM
            EIGENVAL(I) = MATTRIDIAG(I,I)
        ENDDO

            ! --- FROEBIUS EXTRA-DIAG CALCULATION --- !
        RESIDUAL = 0
        DO I = 1,DIM-1
            DO J = I+1,DIM
                IF(ABS(MATTRIDIAG(I,J)) > EPS0) RESIDUAL = RESIDUAL + MATTRIDIAG(I,J)**2
            ENDDO
        ENDDO
        RESIDUAL = SQRT(RESIDUAL)
        COMPT = COMPT + 1
    ENDDO

    IF (COMPT == MAXSTEP) TEST = .TRUE.

    ! --- STEP 3 : PRODUCT WITH MATRIX ROTATION OF TRIDIAGONALE TRANFOR --- !
    IF (SYM) CALL PRODMAT(DIM,MATROT,EIGENVEC,EIGENVEC)
    ! --- STEP 4 : ORDERING EIGENVALS/EIGENVECTS --- !
    CALL ORDERING(DIM,EIGENVAL,EIGENVEC)

    ! --- VERIFICATION --- !
    IF (TEST) THEN
        OPEN(UNIT = ERR, FILE = 'error')
        WRITE(ERR,'(A,10X,A,4X,I14.5)') 'PROBLEM DIAGONALISATION','MAXSTEP =',MAXSTEP
        WRITE(ERR,'(A)') '***********************'
        WRITE(ERR,'(A,4X,ES14.5,4X,A,4X,ES14.5)')'CRIT CONV =',CONV,'RESIDUAL =',RESIDUAL
        WRITE(ERR,'(A)') '***********************'
        WRITE(ERR,'(A)') 'DIAGONAL MATRIX OBTAINED'
        DO I = 1,DIM
            WRITE(ERR,*) (MATTRIDIAG(I,J),J=1,DIM)
        ENDDO
        WRITE(ERR,'(A)') '***********************'
    ENDIF

    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) WRITE(6,'(A)') 'ERROR DIAGONALIZATION, SEE FILE error'


END SUBROUTINE

