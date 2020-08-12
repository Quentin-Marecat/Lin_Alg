SUBROUTINE HOUSESTEP(VEC,LOWESTELEM,MATROT)
    USE INVMATMOD
    ! -------------------------------------------------------------------------------------- !
    ! --- THIS SUBROUTINE APPLIED A HOUSEHOLDER STEP ON THE VECTOR VEC AND GIVE VECHOUSE --- !
    ! --- VEC AND VECHOUSE CAN BE THE SAME BUT VEC WILL BE ERASED -------------------------- !
    ! --- INPUT  _ VEC : VECTOR TO PROCESS THE HOUSEHOLDER PROCESS
    ! ---        _ LOWESTELEM : LOWEST POSITION NON ZERO ELEMENT AFTER HOUSEHOLDER PROCESS, 2 <= LOWESTELEM <= DIM1 - 1 
    ! --- OUTPUT _ MATROT : ROTATION MATRIX : MATROT*VEC = VECHOUSE
    ! -------------------------------------------------------------------------------------- !
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: VEC(DIM)
    INTEGER, INTENT(IN) :: LOWESTELEM
    REAL*8, INTENT(OUT) :: MATROT(DIM,DIM)
    REAL*8 :: ID(DIM,DIM),TAMPVEC(DIM),VECTAMP2
    ID = 0.
    DO I = 1,DIM
        ID(I,I) = 1.
    ENDDO
    IF (LOWESTELEM < 1 .OR. LOWESTELEM > (DIM)) THEN
        WRITE(6,'(A)') 'ERROR HOUSESTEP DIMENSION'
        OPEN(UNIT = ERR,FILE='error',ACTION = 'WRITE')
        WRITE(ERR,'(A,10X,A,I2,10X,A,4X,I2)')  'ERROR HOUSESTEP DIMENSION','DIM =',DIM,'ELEM =',LOWESTELEM
        CLOSE(ERR)
    ELSEIF (LOWESTELEM == DIM) THEN
        MATROT = ID
    ELSE 
        TESTLOG = .TRUE.
        DO I = LOWESTELEM + 1, DIM  ! --- TEST IF THE VECTOR ALREADY HAVE THE GOOD SHAPE 
            IF (ABS(VEC(I)) > EPS0) TESTLOG = .FALSE. 
        ENDDO
        IF (TESTLOG) THEN
            MATROT = ID

        ELSE
            NORME = 0.
            DO I = LOWESTELEM,DIM 
                NORME = NORME + VEC(I)**2
            ENDDO
            TAMPVEC = 0.
            DO I = LOWESTELEM,DIM
                TAMPVEC(I) = SQRT(NORME)*ID(I,LOWESTELEM) - VEC(I)
            ENDDO
            NORME = 0.
            DO I = LOWESTELEM,DIM
                NORME = NORME + TAMPVEC(I)**2
            ENDDO
            MATROT = ID
            DO I = LOWESTELEM,DIM
                DO J = LOWESTELEM,DIM
                    MATROT(I,J) = MATROT(I,J) - (2/NORME)*TAMPVEC(I)*TAMPVEC(J) 
                ENDDO
            ENDDO
        ENDIF
    ENDIF
END SUBROUTINE





SUBROUTINE HOUSESTEPTRI(VEC,LOWESTELEM,MATROT)
    USE INVMATMOD
    ! -------------------------------------------------------------------------------------- !
    ! --- THIS SUBROUTINE APPLIED A HOUSEHOLDER STEP ON THE VECTOR VEC AND GIVE VECHOUSE --- !
    ! --- BUT VEC MUST COME FROM AN TRIDIAGONAL MATRIX ------------------------------------- !
    ! --- VEC AND VECHOUSE CAN BE THE SAME BUT VEC WILL BE ERASED -------------------------- !
    ! --- INPUT  _ VEC : VECTOR TO PROCESS THE HOUSEHOLDER PROCESS
    ! ---        _ LOWESTELEM : LOWEST POSITION NON ZERO ELEMENT AFTER HOUSEHOLDER PROCESS, 2 <= LOWESTELEM <= DIM1 - 1 
    ! --- OUTPUT _ MATROT : ROTATION MATRIX : MATROT*VEC = VECHOUSE
    ! -------------------------------------------------------------------------------------- !
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: VEC(DIM)
    INTEGER, INTENT(IN) :: LOWESTELEM
    REAL*8, INTENT(OUT) :: MATROT(DIM,DIM)
    REAL*8 :: ID(DIM,DIM),TAMPVEC(DIM),VECTAMP2
    INTEGER :: L
    L = LOWESTELEM
    ID = 0.
    DO I = 1,DIM
        ID(I,I) = 1.
    ENDDO
    IF (LOWESTELEM < 1 .OR. LOWESTELEM > (DIM)) THEN
        WRITE(6,'(A)') 'ERROR HOUSESTEP DIMENSION'
        WRITE(ERR,'(A,10X,A,I2,10X,A,4X,I2)')  'ERROR HOUSESTEP DIMENSION','DIM =',DIM,'ELEM =',LOWESTELEM
    ELSEIF (LOWESTELEM == DIM) THEN
        MATROT = ID
    ELSE 
        TESTLOG = .TRUE.
        DO I = L+1,DIM
            IF (ABS(VEC(I+1)) > EPS0) TESTLOG = .FALSE. ! --- TEST IF THE VECTOR ALREADY HAVE THE GOOD SHAPE 
        ENDDO
        IF (TESTLOG) THEN
            MATROT = ID
        ELSE
            NORME = 0.
            DO I = L,L+1
                NORME = NORME + VEC(I)**2
            ENDDO
            TAMPVEC = 0.
            DO I = L,L+1
                TAMPVEC(I) = SQRT(NORME)*ID(I,L) - VEC(I)
            ENDDO
            NORME = 0.
            DO I = L,L+1
                NORME = NORME + TAMPVEC(I)**2
            ENDDO
            MATROT = ID
            DO I = L,L+1
                DO J = L,L+1
                    MATROT(I,J) = MATROT(I,J) - (2/NORME)*TAMPVEC(I)*TAMPVEC(J) 
                ENDDO
            ENDDO
        ENDIF
    ENDIF
END SUBROUTINE

