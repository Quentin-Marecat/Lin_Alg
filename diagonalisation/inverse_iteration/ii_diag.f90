SUBROUTINE II_DIAG(DIM,MAT,EVAL,EIGENVEC,SYM,CHECK)
    ! ------------------------------------------------------------------------------------------------------ !
    ! --- THIS SUBROUTINE FIND THE EIGENVECTOR ASSOCIATED TO THE EIGENVALUE USING POWER ITERATIVE METHOD --- !
    ! ------------------------------------------------------------------------------------------------------ !
    IMPLICIT NONE
    REAL*8, PARAMETER :: EPS0 = 1.D-15, CONV = 1.D-9, SHIFT = 1.D-10
    INTEGER, PARAMETER :: MAXSTEP = 1D5
    INTEGER, INTENT(IN) :: DIM
    LOGICAL, INTENT(IN) :: CHECK, SYM
    REAL*8, INTENT(IN) :: MAT(DIM,DIM), EVAL
    REAL*8, INTENT(OUT) :: EIGENVEC(DIM)
    LOGICAL :: TEST
    REAL*8 :: VEC(DIM), VECP1(DIM), MAT_SHIFT(DIM,DIM)
    REAL*8 :: NORM, MATTAMP(DIM,DIM), DIFF
    REAL*8 :: NOR
    INTEGER :: I, J, K, COMPT, ERR, STAT

    ERR = 97
    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) CLOSE(ERR,STATUS = 'delete')
    TEST = .FALSE.

    EIGENVEC = 0

    ! --- INVERSION OF THE MATRIX + SHIFT --- !
    MAT_SHIFT = MAT
    IF (ABS(SHIFT) > EPS0) THEN
        DO I = 1,DIM
            MAT_SHIFT(I,I) = MAT_SHIFT(I,I) - EVAL - SHIFT
        ENDDO
    ENDIF
    IF (SYM) THEN 
        CALL INV_QR_SYM(DIM,MAT_SHIFT,MATTAMP,.FALSE.)
    ELSE 
        CALL INV_QR(DIM,MAT_SHIFT,MATTAMP,.FALSE.)
    ENDIF
    ! -------------------------------- !

    COMPT = 0
    DIFF = 1
    CALL RANDOM_NUMBER(VEC)
    DO J = 1,DIM
        VEC(J) = 2*VEC(J) - 1
    ENDDO
    CALL NORMVECTOR(DIM,VEC,VEC)
    DO WHILE(DIFF > CONV .AND. COMPT < MAXSTEP)
        CALL EIGENV_II(DIM,MATTAMP,VEC,VECP1,NORM)
        DIFF = 0
        DO J = 1,DIM
            DIFF = DIFF + (VECP1(J) - NORM*VEC(J))**2
        ENDDO
        DIFF = SQRT(DIFF)
        COMPT = COMPT + 1
        CALL NORMVECTOR(DIM,VECP1,VECP1)
        VEC = VECP1
        IF (MOD(COMPT,10) == 0) THEN
            CALL MATAPPLI(DIM,MAT,VEC,VECP1)
            DO J = 1,DIM
                VECP1(J) = VECP1(J) - EVAL*VEC(J)
            ENDDO
            CALL NORMVEC(DIM,VECP1,DIFF)
        ENDIF
    ENDDO
    IF (CHECK .AND. COMPT == MAXSTEP) THEN
        TEST = .TRUE.
        OPEN(UNIT = ERR, FILE = 'error')
        WRITE(ERR,'(A)') 'ERROR CONVERGENCE EIGENVAL'
        WRITE(ERR,'(A,4X,I14.5,4X,A,4X,I14.5,4X,A,4X,ES14.5)') &
        & 'NB EIGENVAL',I,'MAXSTEP',MAXSTEP,'|MAT*VEC - EGNVAL*VEC|', DIFF
    ENDIF
    DO J = 1,DIM
        EIGENVEC(J) = VEC(J)
    ENDDO

    ! --- VERIFICATION --- !
IF (CHECK .AND. TEST) THEN
    OPEN(UNIT = ERR, FILE = 'error', STATUS = 'old')
    WRITE(ERR,'(A,10X,A,4X,ES14.5,4X,A,4X,I14.5)') 'PROBLEM DIAGONALISATION','CRIT CONV =',CONV,'MAXSTEP =',MAXSTEP
    WRITE(ERR,'(A)') '***********************'
    WRITE(ERR,'(A)') 'EIGENVALUE WANTED'
    WRITE(ERR,'(100F14.5)') EVAL
    WRITE(ERR,'(A)') '***********************'
    WRITE(ERR,'(A)') 'ROTATION MATRIX'
    DO I = 1,DIM
        WRITE(ERR,'(100F14.5)') EIGENVEC(I)
    ENDDO
    WRITE(ERR,*) '*******************'
    CALL MATAPPLI(DIM,MAT,EIGENVEC,VEC)
    DO J = 1,DIM
        VEC(J) = VEC(J) - EVAL*EIGENVEC(J)
    ENDDO
    CALL NORMVEC(DIM,VEC,NOR)
    WRITE(ERR,'(A)') 'MAT*EIGENVEC - EIGENVAL*EIGENVEC = 0'
    DO I = 1,DIM
        WRITE(ERR,'(100ES14.5)') VEC(I)
    ENDDO
    WRITE(ERR,'(A)') '***********************'
    WRITE(ERR,'(A)') 'NORME (MAT*EIGENVEC - EIGENVAL*EIGENVEC = 0)'
    WRITE(ERR,'(100ES14.5)') NOR
    WRITE(ERR,'(A)') '***********************'
ENDIF


OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
IF (STAT == 0) WRITE(6,'(A)') 'ERROR DIAGONALIZATION, SEE FILE error'

END SUBROUTINE








SUBROUTINE EIGENV_II(DIM,MAT,VEC,VECP1,EGNVAL)
    ! --------------------------------------------------------------------------------------- !
    ! --- THIS SUBROUTINE GIVE THE SUPPOSED EIGENVALUES OF THE MATRIX USING AVERAGE VALUE --- !
    ! --------------------------------------------------------------------------------------- !
    IMPLICIT NONE
    REAL*8, PARAMETER :: EPS0 = 1.D-15
    INTEGER, INTENT(IN) :: DIM
    REAL*8, INTENT(IN) :: MAT(DIM,DIM), VEC(DIM)
    REAL*8, INTENT(OUT) :: VECP1(DIM), EGNVAL
    REAL*8 :: VECTAMP(DIM)
    INTEGER :: COMPT, I, J
    VECTAMP = VEC
    VECP1 = 0
    EGNVAL = 0
    COMPT = 0
    CALL MATAPPLI(DIM,MAT,VECTAMP,VECP1)
    DO I = 1,DIM
        IF (ABS(VECTAMP(I)) > EPS0) THEN
            EGNVAL = EGNVAL + VECP1(I)/VECTAMP(I)
            COMPT = COMPT + 1
        ENDIF
    ENDDO  
    EGNVAL = EGNVAL/DBLE(COMPT)
END SUBROUTINE

        
