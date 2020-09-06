SUBROUTINE PI_DIAG(DIM,MAT,EIGENVAL,EIGENVEC,NB_EV,CHECK)
    ! ------------------------------------------------------------------------------------- !
    ! --- THIS SUBROUTINE FIND HIGHEST EIGENVALUES/VECTORS USING POWER ITERATIVE METHOD --- !
    ! ------------------------------------------------------------------------------------- !
    IMPLICIT NONE
    REAL*8, PARAMETER :: EPS0 = 1.D-15, CONV = 1.D-9
    INTEGER, PARAMETER :: MAXSTEP = 1D5
    INTEGER, INTENT(IN) :: DIM, NB_EV
    LOGICAL, INTENT(IN) :: CHECK
    REAL*8, INTENT(IN) :: MAT(DIM,DIM)
    REAL*8, INTENT(OUT) :: EIGENVAL(DIM), EIGENVEC(DIM,DIM)
    LOGICAL :: TEST
    REAL*8 :: VEC(DIM), VECP1(DIM)
    REAL*8 :: NORM, MATTAMP(DIM,DIM), DIFF, EIGENVECT(DIM,DIM)
    REAL*8 :: DIAG(DIM,DIM), NOR(DIM)
    INTEGER :: I, J, K, COMPT, ERR, STAT

    ERR = 97
    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) CLOSE(ERR,STATUS = 'delete')
    TEST = .FALSE.

    EIGENVAL = 0
    EIGENVEC = 0
    MATTAMP = MAT

    IF (NB_EV < 1 .OR. NB_EV > DIM) THEN
        OPEN(UNIT = ERR, FILE = 'error')
        WRITE(ERR,'(A,4X,I14.5)') 'ERROR NB EIGENVAL WANTED',NB_EV
        STOP
    ENDIF

    DO I = 1,NB_EV
        COMPT = 0
        DIFF = 1
        CALL RANDOM_NUMBER(VEC)
        DO J = 1,DIM
            VEC(J) = 2*VEC(J) - 1
        ENDDO
        CALL NORMVECTOR(DIM,VEC,VEC)
        DO WHILE(DIFF > CONV .AND. COMPT < MAXSTEP)
            CALL EIGENV(DIM,MATTAMP,VEC,VECP1,NORM)
            DIFF = 0
            DO J = 1,DIM
                DIFF = DIFF + (VECP1(J) - NORM*VEC(J))**2
            ENDDO
            DIFF = SQRT(DIFF)
            COMPT = COMPT + 1
            CALL NORMVECTOR(DIM,VECP1,VECP1)
            VEC = VECP1
        ENDDO
        IF (CHECK .AND. COMPT == MAXSTEP) THEN
            TEST = .TRUE.
            OPEN(UNIT = ERR, FILE = 'error')
            WRITE(ERR,'(A)') 'ERROR CONVERGENCE EIGENVAL'
            WRITE(ERR,'(A,4X,I14.5,4X,A,4X,I14.5,4X,A,4X,ES14.5)') &
            & 'NB EIGENVAL',I,'MAXSTEP',MAXSTEP,'|MAT*VEC - EGNVAL*VEC|', DIFF
        ENDIF
        EIGENVAL(I) = NORM
        DO J = 1,DIM
            EIGENVEC(J,I) = VEC(J)
        ENDDO
        DO J = 1,DIM
            DO K = 1,DIM
                MATTAMP(J,K) = MATTAMP(J,K) - EIGENVAL(I)*EIGENVEC(J,I)*EIGENVEC(K,I)
            ENDDO
        ENDDO
    ENDDO

    CALL ORDERING_0(DIM,EIGENVAL,EIGENVEC)

    ! --- VERIFICATION --- !
IF (CHECK .AND. TEST) THEN
    OPEN(UNIT = ERR, FILE = 'error', STATUS = 'old')
    WRITE(ERR,'(A,10X,A,4X,ES14.5,4X,A,4X,I14.5)') 'PROBLEM DIAGONALISATION','CRIT CONV =',CONV,'MAXSTEP =',MAXSTEP
    WRITE(ERR,'(A)') '***********************'
    WRITE(ERR,'(A)') 'EIGENVALUES OBTAINED'
    WRITE(ERR,'(100F14.5)') (EIGENVAL(J),J=1,NB_EV)
    WRITE(ERR,'(A)') '***********************'
    WRITE(ERR,'(A)') 'ROTATION MATRIX'
    DO I = 1,DIM
        WRITE(ERR,'(100F14.5)') (EIGENVEC(I,J),J = 1,NB_EV)
    ENDDO
    WRITE(ERR,*) '*******************'
    DIAG = 0
    DO I = 1,NB_EV
        DIAG(I,I) = EIGENVAL(I)
    ENDDO
    IF (NB_EV == DIM) THEN
        CALL PRODMAT(DIM,EIGENVEC,DIAG,MATTAMP)
        CALL TRANSPOSE(DIM,EIGENVEC,EIGENVECT)
        CALL PRODMAT(DIM,MATTAMP,EIGENVECT,MATTAMP)
        WRITE(ERR,'(A)') ' Q*D*Q(T) - MAT = 0'
        DO I = 1,DIM
            WRITE(ERR,'(100ES14.5)') (MATTAMP(I,J)-MAT(I,J), J=1,NB_EV)
        ENDDO
        WRITE(ERR,'(A)') '***********************'
    ENDIF
    CALL PRODMAT(DIM,MAT,EIGENVEC,MATTAMP)
    DO I = 1,NB_EV
        DO J = 1,DIM
            VEC(J) = MATTAMP(J,I) - EIGENVAL(I)*EIGENVEC(J,I)
        ENDDO
        CALL NORMVEC(DIM,VEC,NOR(I))
        ENDDO
    WRITE(ERR,'(A)') 'MAT*EIGENVEC - EIGENVAL*EIGENVEC = 0'
    DO I = 1,DIM
        WRITE(ERR,'(100ES14.5)')(MATTAMP(I,J) - EIGENVAL(J)*EIGENVEC(I,J),J=1,NB_EV)
    ENDDO
    WRITE(ERR,'(A)') '***********************'
    WRITE(ERR,'(A)') 'NORME (MAT*EIGENVEC - EIGENVAL*EIGENVEC = 0)'
    WRITE(ERR,'(100ES14.5)') (NOR(I),I=1,NB_EV)
    WRITE(ERR,'(A)') '***********************'
ENDIF


OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
IF (STAT == 0) WRITE(6,'(A)') 'ERROR DIAGONALIZATION, SEE FILE error'

END SUBROUTINE








SUBROUTINE EIGENV(DIM,MAT,VEC,VECP1,EGNVAL)
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

        
