SUBROUTINE TRIDIAG_HOUSE(DIM,MATRIX,TRIDIAGMAT,ROTMAT,CHECK)
    ! --- FOR ANY REPORT OR SUGGESTION, PLEASE CONTACT quentin.marecat@etu.umontpellier.fr --- !
    ! ---------------------------------------------------------------- !
    ! --- THIS SUBROUTINE TRIDIAGONALIZE THE REAL SYMMETRIC MATRIX --- !
    ! --- CHECK = .TRUE. IF VERIFICATIONS HAVE TO BE DONE ------------ !
    ! ---------------------------------------------------------------- !
    IMPLICIT NONE 
    REAL*8, PARAMETER :: EPS0 = 1.D-15, EPS = 1.D-10
    LOGICAL,INTENT(IN) :: CHECK
    LOGICAL :: TEST
    INTEGER :: DIM
    REAL*8 :: MATRIX(DIM,DIM)
    REAL*8,INTENT(OUT) :: TRIDIAGMAT(DIM,DIM), ROTMAT(DIM,DIM)
    REAL*8 :: ROTMATT(DIM,DIM), ID(DIM,DIM)
    REAL*8 :: VEC(DIM), ROTMATTAMP(DIM,DIM)
    REAL*8 :: ROTMATTAMPT(DIM,DIM),MATTAMP(DIM,DIM)
    INTEGER :: LOWESTELEM, I,J,ERR, STAT

    ERR = 97
    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) CLOSE(ERR,STATUS = 'delete')

    TRIDIAGMAT = MATRIX

    ID = 0.
    DO I = 1,DIM
        ID(I,I) = 1.
    ENDDO
    ROTMAT = ID
    ROTMATT = ID
    DO LOWESTELEM = 1,DIM-2
        DO J = 1,DIM
            VEC(J) = TRIDIAGMAT(J,LOWESTELEM)
        ENDDO 
        CALL HOUSESTEP(DIM,VEC,LOWESTELEM+1,ROTMATTAMP)
        CALL PRODMAT(DIM,ROTMATTAMP,TRIDIAGMAT,TRIDIAGMAT)
        CALL TRANSPOSE(DIM,ROTMATTAMP, ROTMATTAMPT)
        CALL PRODMAT(DIM,TRIDIAGMAT,ROTMATTAMPT,TRIDIAGMAT)
        CALL PRODMAT(DIM,ROTMAT,ROTMATTAMP,ROTMAT)
    ENDDO
    ! --- VERIFICATION --- !
    IF (CHECK) THEN
        TEST = .FALSE.
        CALL PRODMAT(DIM,ROTMAT,TRIDIAGMAT,MATTAMP)
        CALL TRANSPOSE(DIM,ROTMAT,ROTMATT)
        CALL PRODMAT(DIM,MATTAMP,ROTMATT,MATTAMP)
        DO I = 1,DIM
            DO J = 1,DIM
                IF (ABS(MATTAMP(I,J)-MATRIX(I,J)) > EPS) TEST = .TRUE.
            ENDDO
        ENDDO
        IF (TEST) THEN
            OPEN(UNIT = ERR, FILE = 'error')
            WRITE(ERR,'(A)')  'PROBLEM TRIDIAGONALISATION'
            WRITE(ERR,'(A)') 'TRIDIAGONAL MATRIX OBTAINED'
            DO I = 1,DIM
                WRITE(ERR,'(100F14.5)') (TRIDIAGMAT(I,J), J=1,DIM)
            ENDDO
            WRITE(ERR,'(A)')  '****************'
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
    IF (STAT == 0) WRITE(6,'(A)') 'ERROR TRIDIAG HOUSEHOLDER, SEE FILE error'
END SUBROUTINE