SUBROUTINE QRHOUSE(DIM,MATRIX,Q,R,CHECK)
    ! --- FOR ANY REPORT OR SUGGESTION, PLEASE CONTACT quentin.marecat@etu.umontpellier.fr --- !
    ! -------------------------------------------------------------------------- !
    ! --- THIS SUBROUTINE GENERATE QR DECOMPOSITION USING HOUSEHOLDER SCHEME --- !
    ! --- HOUSEMAT AND MAT CAN BE THE SAME BUT MAT WILL BE ERASED -------------- !
    ! --- CHECK = .TRUE. IF VERIFICATIONS HAS TO BE DONE ----------------------- !
    ! -------------------------------------------------------------------------- !
    IMPLICIT NONE
    REAL*8, PARAMETER :: EPS0 = 1.D-15, EPS = 1.D-12
    LOGICAL :: TEST
    LOGICAL, INTENT(IN) :: CHECK
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
!        CALL PRODMAT(DIM,MATROTTAMP,R,R)
        CALL PRODMAT_BLOC_QR(DIM,LOWESTELEM-1,MATROTTAMP,R,R) ! --- IMPROVEMENT OF PRODMAT --- !
!        CALL PRODMAT(DIM,Q,MATROTTAMP,Q)
        CALL PRODMAT_BLOC_Q(DIM,LOWESTELEM-1,Q,MATROTTAMP,Q) ! --- IMPROVEMENT OF PRODMAT --- !
    ENDDO

        ! --- VERIFICATION --- !
    IF (CHECK) THEN
        ERR = 97
        TEST = .FALSE.
        CALL PRODMAT(DIM,Q,R,MATTAMP)
        DO I = 1,DIM
            DO J = 1,DIM
                IF (ABS(MATTAMP(J,I) - MATRIX(J,I)) > EPS) TEST = .TRUE.
            ENDDO
        ENDDO

        IF (TEST) THEN
            OPEN(UNIT = ERR,FILE = 'error')
            WRITE(ERR,'(A)') 'ERROR QR DECOMPOSITION HOUSEHOLDER'
            WRITE(ERR,'(A,4X,ES14.5)') 'CRIT ACPT',EPS
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
            WRITE(ERR,'(A)') 'QR - M ='
            DO I = 1,DIM
                WRITE(ERR,'(100ES14.5)') ((MATTAMP(J,I) - MATRIX(J,I)), J = 1,DIM)
            ENDDO 
            WRITE(ERR,'(A)') '*************'
        ENDIF
    ENDIF

    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) WRITE(6,'(A)') 'ERROR QR HOUSEHOLDER, SEE FILE error'
END SUBROUTINE




