SUBROUTINE QRHOUSE(DIM,MATRIX,Q,R)
    ! --- FOR ANY REPORT OR SUGGESTION, PLEASE CONTACT quentin.marecat@etu.umontpellier.fr --- !
    ! -------------------------------------------------------------------------- !
    ! --- THIS SUBROUTINE GENERATE QR DECOMPOSITION USING HOUSEHOLDER SCHEME --- !
    ! --- HOUSEMAT AND MAT CAN BE THE SAME BUT MAT WILL BE ERASED -------------- !
    ! -------------------------------------------------------------------------- !
    IMPLICIT NONE
    REAL*8, PARAMETER :: EPS0 = 1.D-14
    LOGICAL :: TESTLOG
    INTEGER, INTENT(IN) :: DIM
    REAL*8, INTENT(IN) :: MATRIX(DIM,DIM)
    REAL*8, INTENT(OUT) :: R(DIM,DIM),Q(DIM,DIM)
    REAL*8 :: VEC(DIM),MATROTTAMP(DIM,DIM), PROD(DIM,DIM)
    REAL*8 :: MATROTT(DIM,DIM),ID(DIM,DIM),MATTAMP(DIM,DIM)
    REAL*8 :: MATROTM1(DIM,DIM)
    INTEGER :: I,J,K,LOWESTELEM,COMPT,ERR, STAT
    ERR = 97
    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) CLOSE(ERR,STATUS = 'delete')
    
    R = MATRIX 
    ID = 0.
    DO I = 1,DIM
        ID(I,I) = 1.
    ENDDO
    Q = ID
    MATROTM1 = ID
    DO LOWESTELEM = 1,DIM-1
            DO I = 1,DIM
                VEC(I) = R(I,LOWESTELEM)
            ENDDO 
            CALL HOUSESTEP(DIM,VEC,LOWESTELEM,MATROTTAMP)
            CALL PRODMAT(DIM,MATROTTAMP,R,R)
            CALL PRODMAT(DIM,Q,MATROTTAMP,Q)
        ! --- VERIFICATION --- !
        TESTLOG = .FALSE.
        DO I = LOWESTELEM+1,DIM
            IF (ABS(R(I,LOWESTELEM)) > EPS0) TESTLOG = .TRUE.
        ENDDO
        IF (TESTLOG) THEN
!            WRITE(6,'(A)') 'PROBLEM HOUSEHOLDER STEP'
            WRITE(ERR,'(A,10X,A,I5)')  'PROBLEM HOUSEHOLDER STEP','LOWEST ELEMENT =', LOWESTELEM
            WRITE(ERR,*) (R(I,LOWESTELEM),I=1,DIM)
        ENDIF
    ENDDO

        ! --- VERIFICATION --- !
    ERR = 97
    CALL PRODMAT(DIM,Q,R,PROD)
    TESTLOG = .FALSE.
    DO I = 1,DIM
        DO J = 1,DIM
            IF(ABS(PROD(I,J)-MATRIX(I,J)) > EPS0) TESTLOG = .TRUE.
        ENDDO
    ENDDO
    IF (TESTLOG) THEN
        OPEN(UNIT = ERR,FILE = 'error')
        WRITE(ERR,'(A)') 'ERROR QR DECOMPOSITION HOUSEHOLDER'
        WRITE(ERR,'(A)') 'Q ='
        DO I = 1,DIM
            WRITE(ERR,'(100F14.5)') (Q(I,J), J = 1,DIM)
        ENDDO 
        WRITE(ERR,'(A)') '*************'
        WRITE(ERR,'(A)') 'R ='
        DO I = 1,DIM
            WRITE(ERR,'(100F14.5)') (R(I,J), J = 1,DIM)
        ENDDO 
        WRITE(ERR,'(A)') '*************'
    ENDIF

    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) WRITE(6,'(A)') 'ERROR QR HOUSEHOLDER, SEE FILE error'
END SUBROUTINE

