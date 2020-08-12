SUBROUTINE HOUSEHOLDER(MAT,HOUSEMAT,MATROT)
    USE DIAGMATMOD
    ! --- FOR ANY REPORT OR SUGGESTION, PLEASE CONTACT quentin.marecat@etu.umontpellier.fr --- !
    ! ------------------------------------------------------------------------------------------ !
    ! --- THIS SUBROUTINE GENERATE THE HOUSEHOLDER TRANSFORMATION HOUSEMAT TO THE MATRIX MAT --- !
    ! --- HOUSEMAT AND MAT CAN BE THE SAME BUT MAT WILL BE ERASED ------------------------------ !
    ! ------------------------------------------------------------------------------------------ !
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: MAT(DIM,DIM)
    REAL*8, INTENT(OUT) :: HOUSEMAT(DIM,DIM)
    REAL*8 :: VEC(DIM),MATROT(DIM,DIM),MATROTTAMP(DIM,DIM)
    REAL*8 :: MATROTT(DIM,DIM),ID(DIM,DIM),MATTAMP(DIM,DIM)
    REAL*8 :: MATROTM1(DIM,DIM)
    INTEGER :: LOWESTELEM,COMPT
    HOUSEMAT = MAT 
    ID = 0.
    DO I = 1,DIM
        ID(I,I) = 1.
    ENDDO
    MATROT = ID
    MATROTM1 = ID
    DO LOWESTELEM = 1,DIM-1
            DO I = 1,DIM
                VEC(I) = HOUSEMAT(I,LOWESTELEM)
            ENDDO 
            CALL HOUSESTEP(VEC,LOWESTELEM,MATROTTAMP)
            CALL PRODMAT(DIM,MATROTTAMP,HOUSEMAT,HOUSEMAT)
            CALL PRODMAT(DIM,MATROT,MATROTTAMP,MATROT)
        ! --- VERIFICATION --- !
        TESTLOG = .FALSE.
        DO I = LOWESTELEM+1,DIM
            IF (ABS(HOUSEMAT(I,LOWESTELEM)) > EPS0) TESTLOG = .TRUE.
        ENDDO
        IF (TESTLOG) THEN
            OPEN(UNIT = ERR, FILE = 'error')
            WRITE(ERR,'(A,10X,A,I5)')  'PROBLEM HOUSEHOLDER STEP','LOWEST ELEMENT =', LOWESTELEM
            WRITE(ERR,*) (HOUSEMAT(I,LOWESTELEM),I=1,DIM)
        ENDIF
    ENDDO
END SUBROUTINE


