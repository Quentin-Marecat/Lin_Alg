SUBROUTINE HOUSEHOLDER(MAT,HOUSEMAT,MATROT)
    USE INVMATMOD
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
!        TESTLOG = .TRUE.
!        COMPT = 0
!        DO WHILE (TESTLOG .AND. COMPT < MAXSTEPHOUSE)
            DO I = 1,DIM
                VEC(I) = HOUSEMAT(I,LOWESTELEM)
            ENDDO 
            CALL HOUSESTEP(VEC,LOWESTELEM,MATROTTAMP)
            CALL PRODMAT(MATROTTAMP,HOUSEMAT,HOUSEMAT)
            CALL PRODMAT(MATROT,MATROTTAMP,MATROT)
!            ! --- VERIFICATION --- !
!            TESTLOG = .FALSE.
!            DO J = LOWESTELEM+1,DIM
!                IF (ABS(HOUSEMAT(J,LOWESTELEM)) > EPS0) TESTLOG = .TRUE.
!            ENDDO
!            COMPT = COMPT + 1
!        ENDDO
        ! --- VERIFICATION --- !
        TESTLOG = .FALSE.
        DO I = LOWESTELEM+1,DIM
            IF (ABS(HOUSEMAT(I,LOWESTELEM)) > EPS0) TESTLOG = .TRUE.
        ENDDO
        IF (TESTLOG) THEN
!            WRITE(6,'(A)') 'PROBLEM HOUSEHOLDER STEP'
            WRITE(ERR,'(A,10X,A,I5)')  'PROBLEM HOUSEHOLDER STEP','LOWEST ELEMENT =', LOWESTELEM
            WRITE(ERR,*) (HOUSEMAT(I,LOWESTELEM),I=1,DIM)
        ENDIF
    ENDDO
END SUBROUTINE



SUBROUTINE HOUSEHOLDERTRI(MAT,HOUSEMAT,MATROT)
    USE INVMATMOD
    ! ------------------------------------------------------------------------------------------ !
    ! --- THIS SUBROUTINE GENERATE THE HOUSEHOLDER TRANSFORMATION HOUSEMAT TO THE MATRIX MAT --- !
    ! --- BUT MAT MUST BE AN TRIDIAGONAL MATRIX ------------------------------------------------ !
    ! --- HOUSEMAT AND MAT CAN BE THE SAME BUT MAT WILL BE ERASED ------------------------------ !
    ! ------------------------------------------------------------------------------------------ !
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: MAT(DIM,DIM)
    REAL*8, INTENT(OUT) :: HOUSEMAT(DIM,DIM)
    REAL*8 :: VEC(DIM),MATROT(DIM,DIM),MATROTTAMP(DIM,DIM)
    REAL*8 :: MATROTT(DIM,DIM),ID(DIM,DIM),MATTAMP(DIM,DIM)
    INTEGER :: LOWESTELEM
    HOUSEMAT = MAT 
    ID = 0.
    DO I = 1,DIM
        ID(I,I) = 1.
    ENDDO
    MATROT = ID
    DO LOWESTELEM = 1,DIM-1
        DO J = 1,DIM
            VEC(J) = HOUSEMAT(J,LOWESTELEM)
        ENDDO 
        CALL HOUSESTEPTRI(VEC,LOWESTELEM,MATROTTAMP)
        CALL PRODMAT(MATROTTAMP,HOUSEMAT,HOUSEMAT)
                ! --- VERIFICATION --- !
        TESTLOG = .FALSE.
        DO J = LOWESTELEM+1,DIM
            IF (ABS(HOUSEMAT(J,LOWESTELEM)) > EPS0) TESTLOG = .TRUE.
        ENDDO
        IF (TESTLOG) THEN
!            WRITE(6,'(A)') 'PROBLEM HOUSEHOLDER TRIDIAG STEP'
            WRITE(ERR,'(A,10X,A,I5)')  'PROBLEM HOUSEHOLDER STEP','LOWEST ELEMENT =', LOWESTELEM
            WRITE(ERR,*) (HOUSEMAT(I,LOWESTELEM),I=1,DIM)
        ENDIF
        CALL PRODMAT(MATROT,MATROTTAMP,MATROT)
    ENDDO
END SUBROUTINE