PROGRAM TEST
    IMPLICIT NONE
    REAL*8,ALLOCATABLE :: M(:,:), Q(:,:), R(:,:)
    INTEGER :: DIM,I,J, STAT
    OPEN(UNIT = 99,FILE = 'input',STATUS = 'OLD', ACTION= 'READ')
    OPEN(UNIT = 98,FILE = 'output')
!    READ(99,*) DIM
!    ALLOCATE(M(DIM,DIM),Q(DIM,DIM),R(DIM,DIM))
!    DO I = 1,DIM
!        READ(99,*) (M(I,J), J = 1,DIM)
!   ENDDO
    DIM = 10
    ALLOCATE(M(DIM,DIM),Q(DIM,DIM),R(DIM,DIM))
    M = 0
    CALL RANDOM_NUMBER(M)
    DO I = 1,DIM
        DO J = I,DIM
        M(I,J) = 200*M(I,J) - 100
        M(J,I) = M(I,J)
        ENDDO
    ENDDO
!    CALL QRGRAM(DIM,M,Q,R)
    CALL QRHOUSE(DIM,M,Q,R,.TRUE.)

    WRITE(98,*) 'QR DEC'
    WRITE(98,'(A)') 'Q UNITARY MATRIX'
    DO I = 1,DIM
        WRITE(98,'(100F14.5)') (Q(I,J),J = 1,DIM)
    ENDDO
    WRITE(98,*) '*******************'
    WRITE(98,'(A)') 'R UPPER MATRIX'
    DO I = 1,DIM
        WRITE(98,'(100F14.5)') (R(I,J),J = 1,DIM)
    ENDDO
    WRITE(98,*) '*******************'

    OPEN(UNIT = 97, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 2) WRITE(6,'(A)') 'RUN SUCCESFULLY, SEE output FILE' 
END