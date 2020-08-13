SUBROUTINE QR_MAT(DIM,MATRIX,VEC,SOL)
    ! --- FOR ANY REPORT OR SUGGESTION, PLEASE CONTACT quentin.marecat@etu.umontpellier.fr --- !
    ! ----------------------------------------------------------------------------------------------- !
    ! --- THIS SUBROUTINE FIND THE SOLUTION OF THE SYSTEM MATRIX*SOL = VEC USING QR DECOMPOS OF A --- !
    ! --- MATRIX MUST BE SQUARE AND INVERSIBLE ------------------------------------------------------ !
    ! --- RESIDUAL METHOD IS USED TO IMPROVE THE METHOD --------------------------------------------- !
    ! --- operation_mat.f90 QRhouse.f90 AND housestep.f90 ARE NECESSARY ---------------------------------------------------- !
    ! ----------------------------------------------------------------------------------------------- !
    IMPLICIT NONE 
    REAL*8, PARAMETER :: EPS0 = 1.D-15, CONV = 1.D-13
    INTEGER, PARAMETER :: MAXSTEP = 1D5
    INTEGER, INTENT(IN) :: DIM
    REAL*8, INTENT(IN) :: MATRIX(DIM,DIM), VEC(DIM)
    REAL*8, INTENT(OUT) :: SOL(DIM)
    REAL*8 :: Q(DIM,DIM),R(DIM,DIM), QT(DIM,DIM), VECTEST(DIM), SOLINT(DIM), VEC2(DIM)
    REAL*8 :: SOL2(DIM), VECTAMP(DIM), XNORME
    LOGICAL :: TEST
    INTEGER :: I, J, ERR, STAT, COMPT

    ERR = 97
    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) CLOSE(ERR,STATUS = 'delete')

    ! ---------------------------------------------------------------------------------------------------- !
    CALL QRHOUSE(DIM,MATRIX,Q,R)
    CALL TRANSPOSE(DIM,Q,QT)
    ! --- QR DECOMP SUBROUTINE FROM Quentin-Marecat PACKAGES IS USED, PLEASE FIND IT ON MY GITHUB PAGE --- !

    VEC2 = VEC
    OPEN(UNIT = 97, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 2) THEN
        COMPT = 0
        XNORME = 1
        SOL2 = 0
        DO WHILE(XNORME > CONV .AND. COMPT < MAXSTEP) 
            CALL MATAPPLI(DIM,QT,VEC2,SOLINT)
            SOL(DIM) = SOLINT(DIM)/R(DIM,DIM)
            DO I = DIM-1,1,-1
                SOL(I) = SOLINT(I)
                DO J = DIM,I+1,-1
                    SOL(I) = SOL(I) - SOL(J)*R(I,J)
                ENDDO
                SOL(I) = SOL(I)/R(I,I)
            ENDDO
            COMPT = COMPT + 1
            CALL NORMVEC(DIM,SOL,XNORME)
            CALL MATAPPLI(DIM,MATRIX,SOL,VECTAMP)
            DO I = 1,DIM
                SOL2(I) = SOL2(I) + SOL(I)
                VEC2(I) = VEC2(I) - VECTAMP(I) 
            ENDDO
        ENDDO
    ENDIF
    SOL = SOL2
    ! --- VERIFICATION --- !
    CALL MATAPPLI(DIM,MATRIX,SOL,VECTEST)
    TEST = .FALSE.
    DO I = 1,DIM
        IF (ABS(VECTEST(I) - VEC(I)) > EPS0) TEST = .TRUE.
    ENDDO
    IF (TEST) THEN
        WRITE(6,*) SOL
        WRITE(6,*) VECTEST
        OPEN(UNIT = ERR,FILE = 'error')
        WRITE(ERR,'(A)') 'PROBLEM SOLUTION'
        WRITE(ERR,'(A)') 'A = '
        DO I = 1,DIM
            WRITE(ERR,'(100F14.5)') (MATRIX(I,J), J = 1,DIM)
        ENDDO
        WRITE(ERR,'(A)') '************'
        WRITE(ERR,'(A)') 'x = '
        WRITE(ERR,'(100F14.5)') (SOL(I), I = 1,DIM)
        WRITE(ERR,'(A)') '************'
        WRITE(ERR,'(A)') 'b = '
        WRITE(ERR,'(100F14.5)') (VEC(I), I = 1,DIM)
        WRITE(ERR,'(A)') '************'
    ENDIF

    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) WRITE(6,'(A)') 'ERROR LINEAR SOLUTION, SEE FILE error'

END SUBROUTINE
