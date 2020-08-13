SUBROUTINE QRGRAM(DIM,MATRIX,Q,R)
    ! --- FOR ANY REPORT OR SUGGESTION, PLEASE CONTACT quentin.marecat@etu.umontpellier.fr --- !
    ! --------------------------------------------------------------------------- !
    ! --- THIS SUBROUTINE GENERATE QR DECOMPOSITION USING GRAM_SCHMIDT SCHEME --- !
    ! --- HOUSEMAT AND MAT CAN BE THE SAME BUT MAT WILL BE ERASED --------------- !
    ! --- operation_mat.f90 AND GM.f90 ARE NECESSARY ---------------------------- !
    ! --------------------------------------------------------------------------- !
    IMPLICIT NONE
    REAL*8, PARAMETER :: EPS0 = 1.D-14
    LOGICAL :: TEST
    INTEGER, INTENT(IN) :: DIM
    REAL*8, INTENT(IN) :: MATRIX(DIM,DIM)
    REAL*8, INTENT(OUT) :: R(DIM,DIM),Q(DIM,DIM)
    INTEGER :: I,J,K, ERR, STAT
    REAL*8 :: TAMP, PROD(DIM,DIM), QT(DIM,DIM)
    ERR = 97
    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) CLOSE(ERR,STATUS = 'delete')

    R = 0.
    CALL GM(DIM,MATRIX,Q)
    CALL TRANSPOSE(DIM,Q,QT)
    CALL PRODMAT(DIM,QT,MATRIX,R)



    ! --- VERIFICATION --- !
    CALL PRODMAT(DIM,Q,R,PROD)
    TEST = .FALSE.
    DO I = 1,DIM
        DO J = 1,DIM
            IF(ABS(PROD(I,J)-MATRIX(I,J)) > EPS0) TEST = .TRUE.
        ENDDO
    ENDDO
    IF (TEST) THEN
        OPEN(UNIT = ERR,FILE = 'error')
        WRITE(ERR,'(A)') 'ERROR QR DECOMPOSITION GRAM_SCHMIDT'
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
    IF (STAT == 0) WRITE(6,'(A)') 'ERROR QR GRAM-SCHMIDT, SEE FILE error'
END SUBROUTINE
                        
