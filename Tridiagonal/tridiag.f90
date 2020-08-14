SUBROUTINE TRIDIAG(DIM,MATRIX,TRIDIAGMAT,ROTMAT)
    ! --- FOR ANY REPORT OR SUGGESTION, PLEASE CONTACT quentin.marecat@etu.umontpellier.fr --- !
    ! ---------------------------------------------------------------- !
    ! --- THIS SUBROUTINE TRIDIAGONALIZE THE REAL SYMMETRIC MATRIX --- !
    ! --- housestep.f90 IS NECESSARY --------------------------------- !
    ! ---------------------------------------------------------------- !
    IMPLICIT NONE 
    REAL*8, PARAMETER :: EPS0 = 1.D-12
    LOGICAL :: TEST
    INTEGER :: DIM
    REAL*8 :: MATRIX(DIM,DIM)
    REAL*8,INTENT(OUT) :: TRIDIAGMAT(DIM,DIM), ROTMAT(DIM,DIM)
    REAL*8 :: ROTMATT(DIM,DIM), TRIANGSUP(DIM,DIM), ID(DIM,DIM)
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
                    ! --- VERIFICATION --- !
        TEST = .FALSE.
        DO I = LOWESTELEM+2,DIM
            IF (ABS(TRIDIAGMAT(I,LOWESTELEM)) > EPS0 .OR. ABS(TRIDIAGMAT(LOWESTELEM,I)) > EPS0 &
            & .OR. ABS(TRIDIAGMAT(LOWESTELEM,LOWESTELEM+1)-TRIDIAGMAT(LOWESTELEM+1,LOWESTELEM)) > EPS0) TEST = .TRUE.
        ENDDO
        IF (TEST) THEN
            OPEN(UNIT = ERR, FILE = 'error')
            WRITE(ERR,'(A)')  'PROBLEM HOUSESTEP'
            DO I = 1,DIM
                WRITE(ERR,'(100E14.5)') (TRIDIAGMAT(I,J), J=1,DIM)
            ENDDO
            WRITE(ERR,'(A)')  '****************'
        ENDIF
        CALL PRODMAT(DIM,ROTMAT,ROTMATTAMP,ROTMAT)
        CALL PRODMAT(DIM,ROTMATTAMPT,ROTMATT,ROTMATT)
    ENDDO
    ! --- VERIFICATION --- !
    TEST = .FALSE.
    DO I = 1,DIM-2
        DO J = I+2,DIM
            IF (ABS(TRIDIAGMAT(I,J)) > EPS0 .OR. ABS(TRIDIAGMAT(J,I)) > EPS0 &
            & .AND. ABS(TRIDIAGMAT(I,I+1)-TRIDIAGMAT(I+1,I)) < EPS0) TEST = .TRUE.
        ENDDO
    ENDDO
    CALL PRODMAT(DIM,ROTMAT,TRIDIAGMAT,MATTAMP)
    CALL PRODMAT(DIM,MATTAMP,ROTMATT,MATTAMP)
    DO I = 1,DIM
        DO J = 1,DIM
            IF(ABS(MATTAMP(I,J)-MATRIX(I,J)) > EPS0) TEST = .TRUE.
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
        WRITE(ERR,'(A)') ' DIFFERENCE'
        DO I = 1,DIM
            WRITE(ERR,'(100E14.5)') (MATTAMP(I,J)-MATRIX(I,J), J=1,DIM)
        ENDDO
        WRITE(ERR,'(A)')  '****************'
    ENDIF

    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) WRITE(6,'(A)') 'ERROR TRIDIAG HOUSEHOLDER, SEE FILE error'
END SUBROUTINE