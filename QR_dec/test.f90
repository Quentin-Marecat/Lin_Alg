PROGRAM TEST
    IMPLICIT NONE
    REAL*8,ALLOCATABLE :: M(:,:), Q(:,:), R(:,:)
    INTEGER :: DIM,I,J, STAT
    LOGICAL :: CHECK
    OPEN(UNIT = 98,FILE = 'output')

!    OPEN(UNIT = 99,FILE = 'input',STATUS = 'OLD', ACTION= 'READ')
!    READ(99,*) DIM
!    ALLOCATE(M(DIM,DIM),Q(DIM,DIM),R(DIM,DIM))
!    DO I = 1,DIM
!        READ(99,*) (M(I,J), J = 1,DIM)
!   ENDDO

    DIM = 100
    ALLOCATE(M(DIM,DIM),Q(DIM,DIM),R(DIM,DIM))
    M = 0
    CALL RANDOM_NUMBER(M)
    DO I = 1,DIM
        DO J = I,DIM
        M(I,J) = 200*M(I,J) - 100
        M(J,I) = M(I,J)
        ENDDO
    ENDDO

    CHECK = .TRUE.
!    CALL QRGRAM(DIM,M,Q,R,CHECK)
    CALL QRHOUSE(DIM,M,Q,R,CHECK)

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