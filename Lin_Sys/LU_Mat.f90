SUBROUTINE LU_MAT(DIM,MATRIX,VEC,SOL,CHECK)
    ! --- FOR ANY REPORT OR SUGGESTION, PLEASE CONTACT quentin.marecat@etu.umontpellier.fr --- !
    ! ----------------------------------------------------------------------------------------------- !
    ! --- THIS SUBROUTINE FIND THE SOLUTION OF THE SYSTEM MATRIX*SOL = VEC USING LU DECOMPOS OF A --- !
    ! --- MATRIX MUST BE SQUARE AND INVERSIBLE ------------------------------------------------------ !
    ! --- RESIDUAL METHOD IS USED TO IMPROVE THE METHOD --------------------------------------------- !
    ! --- operation_mat.f90 LU.f90 ARE NECESSARY ---------------------------------------------------- !
    ! ----------------------------------------------------------------------------------------------- !
    IMPLICIT NONE 
    LOGICAL, INTENT(IN) :: CHECK
    REAL*8, PARAMETER :: EPS0 = 1.D-15, CONV = 1.D-13, EPS = 1.D-10
    INTEGER, PARAMETER :: MAXSTEP = 1D2
    INTEGER, INTENT(IN) :: DIM
    REAL*8, INTENT(IN) :: MATRIX(DIM,DIM), VEC(DIM)
    REAL*8, INTENT(OUT) :: SOL(DIM)
    REAL*8 :: L(DIM,DIM),U(DIM,DIM), VECTEST(DIM), SOLINT(DIM), VEC2(DIM)
    REAL*8 :: SOL2(DIM), VECTAMP(DIM), XNORME
    LOGICAL :: TEST
    INTEGER :: I, J, ERR, STAT, COMPT

    ERR = 97
    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) CLOSE(ERR,STATUS = 'delete')

    ! ---------------------------------------------------------------------------------------------------- !
    CALL LU(DIM,MATRIX,L,U,.FALSE.)
    ! --- LU DECOMP SUBROUTINE FROM Quentin-Marecat PACKAGES IS USED, PLEASE FIND IT ON MY GITHUB PAGE --- !

    VEC2 = VEC
    OPEN(UNIT = 97, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 2) THEN
        COMPT = 0
        XNORME = 1
        SOL2 = 0
        ! --- RESIDUAL METHOD IMPROVEMENT --- !
        DO WHILE(XNORME > CONV .AND. COMPT < MAXSTEP) 
            SOLINT(1) = VEC2(1)/L(1,1)
            DO I = 2,DIM
                SOLINT(I) = VEC2(I)
                DO J = 1,I-1
                    SOLINT(I) = SOLINT(I) - SOLINT(J)*L(I,J)
                ENDDO
                SOLINT(I) = SOLINT(I)/L(I,I)
            ENDDO
            SOL(DIM) = SOLINT(DIM)/U(DIM,DIM)
            DO I = DIM-1,1,-1
                SOL(I) = SOLINT(I)
                DO J = DIM,I+1,-1
                    SOL(I) = SOL(I) - SOL(J)*U(I,J)
                ENDDO
                SOL(I) = SOL(I)/U(I,I)
            ENDDO
            COMPT = COMPT + 1
            CALL NORMVEC(DIM,SOL,XNORME)
            CALL MATAPPLI(DIM,MATRIX,SOL,VECTAMP)
            ! --- RESIDUAL METHOD IMPROVEMENT --- !
            DO I = 1,DIM
                SOL2(I) = SOL2(I) + SOL(I)
                VEC2(I) = VEC2(I) - VECTAMP(I) 
            ENDDO
        ENDDO
    ENDIF
    SOL = SOL2
    IF(COMPT == MAXSTEP) THEN
        OPEN(UNIT = ERR,FILE = 'error')
        WRITE(ERR,'(A,4X,A,4X,I5)') 'ERROR CONV RESIDUAL METHOD','MAXSTEP =',MAXSTEP
    ENDIF
    
    ! --- VERIFICATION --- !
    IF (CHECK) THEN
        CALL MATAPPLI(DIM,MATRIX,SOL,VECTEST)
        TEST = .FALSE.
        DO I = 1,DIM
            IF (ABS(VECTEST(I) - VEC(I)) > EPS) TEST = .TRUE.
        ENDDO
        IF (TEST) THEN
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
            WRITE(ERR,'(A,4X,ES14.5)')'EPS',EPS
            WRITE(ERR,'(A)') 'ERROR = MAT*SOL - VEC '
            WRITE(ERR,'(100ES14.4)') (VECTEST(I)-VEC(I), I = 1,DIM)
            WRITE(ERR,'(A)') '************'
        ENDIF
    ENDIF

    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) WRITE(6,'(A)') 'ERROR LINEAR SOLUTION, SEE FILE error'

END SUBROUTINE
