SUBROUTINE HESSENBERG(MAT,HESSMATSUP,MATROT)
    USE DIAGMATMOD
    ! --- FOR ANY REPORT OR SUGGESTION, PLEASE CONTACT quentin.marecat@etu.umontpellier.fr --- !
    ! ---------------------------------------------------------------------------------------------- !
    ! --- THIS SUBROUTINE MAKE THE UPPER HESSENBERG TRANSFORMATION TO MATRIX AND GAVE HESSMATSUP --- !
    ! --- MATRIX AND HESSMATSUP CAN BE THE SAME ---------------------------------------------------- !
    ! ---------------------------------------------------------------------------------------------- !
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: MAT(DIM,DIM)
    REAL*8, INTENT(OUT) :: HESSMATSUP(DIM,DIM)
    REAL*8 :: MATTAMP(DIM:DIM), ID(DIM,DIM),MATROT(DIM,DIM), VEC(DIM)
    REAL*8 :: MATROTTAMP(DIM,DIM)
    INTEGER :: LOWESTELEM
    HESSMATSUP = MAT
    ID = 0.
    DO I = 1,DIM
        ID(I,I) = 1.
    ENDDO
    MATROT = ID
    DO LOWESTELEM = 1,DIM-2
        DO J = 1,DIM
            VEC(J) = HESSMATSUP(J,LOWESTELEM)
        ENDDO 
        CALL HOUSESTEP(VEC,LOWESTELEM+1,MATROTTAMP)
        CALL PRODMAT(DIM,MATROTTAMP,HESSMATSUP,HESSMATSUP)
                    ! --- VERIFICATION --- !
        TESTLOG = .FALSE.
        DO I = LOWESTELEM+2,DIM
            IF (ABS(HESSMATSUP(I,LOWESTELEM)) > EPS0 ) TESTLOG = .TRUE.
        ENDDO
        IF (TESTLOG) THEN
            OPEN(UNIT = ERR, FILE = 'error')
            WRITE(6,'(A)') 'PROBLEM HOUSESTEP'
            OPEN(UNIT = ERR,FILE='error',ACTION = 'WRITE')
            WRITE(ERR,'(A)')  'PROBLEM HOUSESTEP'
            CLOSE(ERR)
        ENDIF
        CALL PRODMAT(DIM,MATROT,MATROTTAMP,MATROT)
    ENDDO
            ! --- VERIFICATION --- !
    TESTLOG = .FALSE.
    DO I = 1,DIM-2
        DO J = I+2,DIM
            IF (ABS(HESSMATSUP(J,I)) > EPS0 ) TESTLOG = .TRUE.
        ENDDO
    ENDDO
    IF (TESTLOG) THEN
        WRITE(6,'(A)') 'PROBLEM HESSENBERG'
        OPEN(UNIT = ERR,FILE='error',ACTION = 'WRITE')
        WRITE(ERR,'(A)')  'PROBLEM HESSENBERG'
        CLOSE(ERR)
    ENDIF
END SUBROUTINE