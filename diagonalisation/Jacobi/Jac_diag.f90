SUBROUTINE JAC_DIAG(DIM,MATRIX,EIGENVAL,EIGENVEC,CHECK)
    ! ----------------------------------------------------------------------------------- !
    ! --- THIS SUBROUTINE DIAGONALIZE THE MATRIX USING GIVENS_JACOBI ITERATIVE SCHEME --- !
    ! --- THE MATRIX MUST BE REAL AND SYMMETRIC ----------------------------------------- !
    ! ----------------------------------------------------------------------------------- !
    IMPLICIT NONE
    REAL*8, PARAMETER :: EPS0 = 1.D-14, CONV = 1.D-12, MAXSWEEP = 1D2, PI = 3.14159265359
    LOGICAL :: TEST
    LOGICAL, INTENT(IN) :: CHECK
    INTEGER, INTENT(IN) :: DIM
    REAL*8, INTENT(IN) :: MATRIX(DIM,DIM)
    REAL*8, INTENT(OUT) :: EIGENVAL(DIM), EIGENVEC(DIM,DIM)
    REAL*8 :: ID(DIM,DIM), MATTAMP(DIM,DIM),EIGENVECT(DIM,DIM)
    REAL*8 :: VEC(DIM), NOR(DIM)
    REAL*8 :: JAC(DIM,DIM), NORMF, THETA, T, C, S, TAMP
    INTEGER :: I, J, K, L, COMPT, ERR, STAT

    ERR = 97
    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) CLOSE(ERR,STATUS = 'delete')

    ID = 0
    DO I = 1,DIM
        ID(I,I) = 1
    ENDDO
    EIGENVEC = ID

    JAC = MATRIX
    NORMF = 0.
    COMPT = 0

    DO I = 1,DIM-1
        DO J = I+1,DIM
             NORMF = NORMF + JAC(J,I)**2
        ENDDO
    ENDDO
    NORMF = SQRT(NORMF)

    DO WHILE (NORMF > CONV .AND. COMPT < MAXSWEEP)
        DO I = 1,DIM-1
            DO J = I+1,DIM
                IF (ABS(JAC(J,I)) > EPS0) THEN
                    IF (ABS(JAC(I,I)-JAC(J,J)) > EPS0) THEN
                        THETA = DATAN(2*JAC(J,I)/(JAC(I,I)-JAC(J,J)))/2
                    ELSE 
                        THETA = PI/4
                    ENDIF
                    S = SIN(THETA)
                    C = COS(THETA)
                    TAMP = JAC(I,I)
                    JAC(I,I) = (S**2)*JAC(J,J) + 2*S*C*JAC(J,I) + (C**2)*JAC(I,I)
                    JAC(J,J) = (C**2)*JAC(J,J) - 2*S*C*JAC(J,I) + (S**2)*TAMP
                    DO K = 1,DIM
                        IF (K /= J .AND. K /= I) THEN
                            TAMP = JAC(K,J)
                            JAC(K,J) = C*JAC(K,J) - S*JAC(K,I)
                            JAC(K,I) = S*TAMP + C*JAC(K,I)
                            JAC(J,K) = JAC(K,J)
                            JAC(I,K) = JAC(K,I)
                        ENDIF
                        TAMP = EIGENVEC(K,J)
                        EIGENVEC(K,J) = C*EIGENVEC(K,J) - S*EIGENVEC(K,I)
                        EIGENVEC(K,I) = S*TAMP + C*EIGENVEC(K,I)
                    ENDDO
                    JAC(I,J) = 0
                    JAC(J,I) = 0
                ENDIF
            ENDDO
        ENDDO

        NORMF = 0
        DO I = 1,DIM-1
            DO J = I+1,DIM
                 NORMF = NORMF + JAC(J,I)**2
            ENDDO
        ENDDO
        NORMF = SQRT(NORMF)
        COMPT = COMPT + 1
    ENDDO

    DO I = 1,DIM
        EIGENVAL(I) = JAC(I,I)
    ENDDO

    CALL ORDERING(DIM,EIGENVAL,EIGENVEC)

        ! --- VERIFICATION --- !
    IF (CHECK) THEN
        TEST = .FALSE.
        IF (COMPT == MAXSWEEP) TEST = .TRUE.
        IF (TEST) THEN
            OPEN(UNIT = ERR, FILE = 'error')
            WRITE(ERR,'(A,10X,A,4X,ES14.1,4X,A,4X,ES14.1)') 'PROBLEM DIAGONALISATION','CRIT CONV =',CONV,'MAXSWEEP =',MAXSWEEP
            WRITE(ERR,'(A)') '***********************'
            WRITE(ERR,'(A)') 'DIAGONAL MATRIX OBTAINED'
            DO I = 1,DIM
                WRITE(ERR,'(100F14.5)') (JAC(I,J),J=1,DIM)
            ENDDO
            WRITE(ERR,'(A)') '***********************'
            WRITE(ERR,'(A)') 'ROTATION MATRIX'
            DO I = 1,DIM
                WRITE(ERR,'(100F14.5)') (EIGENVEC(I,J),J = 1,DIM)
            ENDDO
            WRITE(ERR,*) '*******************'
            JAC = ID
            DO I = 1,DIM
                JAC(I,I) = EIGENVAL(I)
            ENDDO
            CALL PRODMAT(DIM,EIGENVEC,JAC,MATTAMP)
            CALL TRANSPOSE(DIM,EIGENVEC,EIGENVECT)
            CALL PRODMAT(DIM,MATTAMP,EIGENVECT,MATTAMP)
            WRITE(ERR,'(A)') ' Q*D*Q(T) - MAT'
            DO I = 1,DIM
                WRITE(ERR,'(100ES14.5)') (MATTAMP(I,J)-MATRIX(I,J), J=1,DIM)
            ENDDO
            WRITE(ERR,'(A)') '***********************'
            CALL PRODMAT(DIM,MATRIX,EIGENVEC,MATTAMP)
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