SUBROUTINE II_DIAG(DIM,MAT,NB_EV,EIGENVAL,EIGENVEC,SYM,CHECK)
    ! --------------------------------------------------------------------------------------------------------- !
    ! --- THIS SUBROUTINE COMPUTE EIGENVECTORS OF MAT ASSOCIATED TO THE EIGENVALUES USING INVERSE ITERATION --- !
    ! --------------------------------------------------------------------------------------------------------- !
    IMPLICIT NONE
    REAL*8, PARAMETER :: EPS0 = 1.D-15, INV = 1.D-10, CONV = 1.D-10
    INTEGER, PARAMETER :: MAXSTEP = 1D4
    LOGICAL, INTENT(IN) :: SYM, CHECK
    INTEGER, INTENT(IN) :: DIM, NB_EV
    REAL*8, INTENT(IN) :: MAT(DIM,DIM), EIGENVAL(DIM)
    REAL*8, INTENT(OUT) :: EIGENVEC(DIM,DIM)
    LOGICAL :: TEST
    REAL*8 :: MATTAMP(DIM,DIM), EVAL, INVMAT(DIM,DIM), ROTMAT(DIM,DIM)
    REAL*8 :: VEC(DIM), DIFF, EIGENVECT(DIM,DIM), DIAG(DIM,DIM), NOR(DIM)
    REAL*8 :: VECP1(DIM), VECTAMP(DIM)
    INTEGER :: I,J, COMPT, ERR, STAT

    ERR = 97
    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) CLOSE(ERR,STATUS = 'delete')
    TEST = .FALSE.

    EIGENVEC = 0

    DO I = 1,NB_EV
        MATTAMP = MAT
        EVAL = EIGENVAL(I) + INV
        DO J = 1,DIM
            MATTAMP(J,J) = MATTAMP(J,J) - EVAL
        ENDDO

        IF (SYM) THEN
            CALL INV_QR_SYM(DIM,MATTAMP,INVMAT,ROTMAT,.FALSE.)
        ELSE
            CALL INV_QR(DIM,MATTAMP,INVMAT,.FALSE.)
        ENDIF

        CALL RANDOM_NUMBER(VEC)
        DO J = 1,DIM
            VEC(J) = 2*VEC(J) - 1
        ENDDO
        CALL NORMVECTOR(DIM,VEC,VEC)
        DIFF = 1
        COMPT = 0
        DO WHILE(DIFF > CONV .AND. COMPT < MAXSTEP)
            IF(SYM) THEN
                CALL MATAPPLI_TRI(DIM,INVMAT,VEC,VECP1)
            ELSE
                CALL MATAPPLI(DIM,INVMAT,VEC,VECP1)
            ENDIF
            CALL NORMVECTOR(DIM,VECP1,VECP1)
            DO J = 1,DIM
                IF (VECP1(1)*VEC(1) > EPS0) THEN
                    VECTAMP(J) = VECP1(J) - VEC(J)
                ELSE 
                    VECTAMP(J) = VECP1(J) + VEC(J)
                ENDIF
            ENDDO
            CALL NORMVEC(DIM,VECTAMP,DIFF)
            VEC = VECP1
            COMPT = COMPT + 1
        ENDDO

        IF (CHECK .AND. COMPT == MAXSTEP) THEN
            TEST = .TRUE.
            OPEN(UNIT = ERR, FILE = 'error')
            WRITE(ERR,'(A)') 'ERROR CONVERGENCE EIGENVAL'
            WRITE(ERR,'(A,4X,I14.5,4X,A,4X,I14.5,4X,A,4X,ES14.5)') &
            & 'NB EIGENVAL',I,'MAXSTEP',MAXSTEP,'NORME CONV VECTOR', DIFF
        ENDIF

        IF (SYM) THEN
             CALL MATAPPLI(DIM,ROTMAT,VEC,VEC)
             CALL NORMVECTOR(DIM,VEC,VEC)
        ENDIF

        DO J = 1,DIM
            EIGENVEC(J,I) = VEC(J)
        ENDDO

    ENDDO
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
    IF (STAT == 0) WRITE(6,'(A)') 'ERROR DIAGONALIZATION EIGENVECTOR, SEE FILE error'

END SUBROUTINE





