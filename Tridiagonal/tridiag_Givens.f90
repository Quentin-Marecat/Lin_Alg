SUBROUTINE TRIDIAG_GIVENS(DIM,MATRIX,TRIDIAGMAT,ROTMAT,CHECK)
        ! --- FOR ANY REPORT OR SUGGESTION, PLEASE CONTACT quentin.marecat@etu.umontpellier.fr --- !
    ! --------------------------------------------------------------------------------------------- !
    ! --- THIS SUBROUTINE TRIDIAGONALIZE THE REAL SYMMETRIC MATRIX USING GIVENS ROTATION MATRIX --- !
    ! --- operation_mat.f90 IS NECESSARY ---------------------------------------------------------- !
    ! --- CHECK = .TRUE. IF VERIFICATIONS HAVE TO BE DONE ------------ !
    ! --------------------------------------------------------------------------------------------- !
    IMPLICIT NONE 
    REAL*8, PARAMETER :: EPS0 = 1.D-15, CONV = 1.D-12, MAXSWEEP = 1D2, PI = 3.14159265359
    LOGICAL,INTENT(IN) :: CHECK
    LOGICAL :: TEST
    INTEGER :: DIM
    REAL*8 :: MATRIX(DIM,DIM)
    REAL*8,INTENT(OUT) :: TRIDIAGMAT(DIM,DIM), ROTMAT(DIM,DIM)
    REAL*8 ::  ID(DIM,DIM), MATTAMP(DIM,DIM), ROTMATT(DIM,DIM)
    REAL*8 ::  NORMF, THETA, T, C, S, TAMP
    INTEGER :: I, J, K, L, COMPT, ERR, STAT

    ERR = 97
    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) CLOSE(ERR,STATUS = 'delete')

    ID = 0
    DO I = 1,DIM
        ID(I,I) = 1
    ENDDO
    ROTMAT = ID

    TRIDIAGMAT = MATRIX
    NORMF = 1.
    COMPT = 0

    DO WHILE (NORMF > CONV .AND. COMPT < MAXSWEEP)
        DO I = 1,DIM-2
            DO J = I+2,DIM
                IF (ABS(TRIDIAGMAT(J,I)) > EPS0) THEN
                    IF (ABS(TRIDIAGMAT(I,I)-TRIDIAGMAT(J,J)) > EPS0) THEN
                        THETA = DATAN(2*TRIDIAGMAT(J,I)/(TRIDIAGMAT(I,I)-TRIDIAGMAT(J,J)))/2
                    ELSE 
                        THETA = PI/4
                    ENDIF
                    S = SIN(THETA)
                    C = COS(THETA)
                    TAMP = TRIDIAGMAT(I,I)
                    TRIDIAGMAT(I,I) = (S**2)*TRIDIAGMAT(J,J) + 2*S*C*TRIDIAGMAT(J,I) + (C**2)*TRIDIAGMAT(I,I)
                    TRIDIAGMAT(J,J) = (C**2)*TRIDIAGMAT(J,J) - 2*S*C*TRIDIAGMAT(J,I) + (S**2)*TAMP
                    DO K = 1,DIM
                        IF (K /= J .AND. K /= I) THEN
                            TAMP = TRIDIAGMAT(K,J)
                            TRIDIAGMAT(K,J) = C*TRIDIAGMAT(K,J) - S*TRIDIAGMAT(K,I)
                            TRIDIAGMAT(K,I) = S*TAMP + C*TRIDIAGMAT(K,I)
                            TRIDIAGMAT(J,K) = TRIDIAGMAT(K,J)
                            TRIDIAGMAT(I,K) = TRIDIAGMAT(K,I)
                        ENDIF
                        TAMP = ROTMAT(K,J)
                        ROTMAT(K,J) = C*ROTMAT(K,J) - S*ROTMAT(K,I)
                        ROTMAT(K,I) = S*TAMP + C*ROTMAT(K,I)
                    ENDDO
                    TRIDIAGMAT(I,J) = 0
                    TRIDIAGMAT(J,I) = 0
                ENDIF
            ENDDO
        ENDDO

        NORMF = 0
        DO I = 1,DIM-2
            DO J = I+2,DIM
                NORMF = NORMF + TRIDIAGMAT(I,J)**2
            ENDDO
        ENDDO
        NORMF = SQRT(NORMF)
        COMPT = COMPT + 1
    ENDDO


        ! --- VERIFICATION --- !
    IF (CHECK) THEN
        TEST = .FALSE.
        IF (COMPT == MAXSWEEP) TEST = .TRUE.
        IF (TEST) THEN
            CALL PRODMAT(DIM,ROTMAT,TRIDIAGMAT,MATTAMP)
            CALL TRANSPOSE(DIM,ROTMAT,ROTMATT)
            CALL PRODMAT(DIM,MATTAMP,ROTMATT,MATTAMP)
            OPEN(UNIT = ERR, FILE = 'error')
            WRITE(ERR,'(A,10X,A,4X,ES14.1,4X,A,4X,ES14.1,4X,A,4X,ES14.5)') 'PROBLEM TRIGONALISATION','NORM FROEB EXTRA-TRIDIAG'&
            &,NORMF, 'CRIT CONV =',CONV,'MAXSWEEP =',MAXSWEEP
            WRITE(ERR,'(A)') '***********************'
            WRITE(ERR,'(A)') 'TRIDIAGONAL MATRIX OBTAINED'
            DO I = 1,DIM
                WRITE(ERR,'(100F14.5)') (TRIDIAGMAT(I,J),J=1,DIM)
            ENDDO
            WRITE(ERR,'(A)') '***********************'
            WRITE(ERR,'(A)') 'ROTATION MATRIX'
            DO I = 1,DIM
                WRITE(ERR,'(100F14.5)') (ROTMAT(I,J),J = 1,DIM)
            ENDDO
            WRITE(98,*) '*******************'
            WRITE(ERR,'(A)') ' Q*TRI*Q(T) - MAT'
            DO I = 1,DIM
                WRITE(ERR,'(100ES14.5)') (MATTAMP(I,J)-MATRIX(I,J), J=1,DIM)
            ENDDO
            WRITE(ERR,'(A)')  '****************'
        ENDIF
    ENDIF

    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) WRITE(6,'(A)') 'ERROR TRIGONALIZATION, SEE FILE error'

END SUBROUTINE