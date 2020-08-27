PROGRAM TEST
    IMPLICIT NONE
    REAL*8,ALLOCATABLE :: M(:,:), Q(:,:), R(:,:)
    INTEGER :: DIM,I,J, STAT
    OPEN(UNIT = 98,FILE = 'output')

!    OPEN(UNIT = 99,FILE = 'input',STATUS = 'OLD', ACTION= 'READ')
!    READ(99,*) DIM
!    ALLOCATE(M(DIM,DIM),Q(DIM,DIM),R(DIM,DIM))
!    DO I = 1,DIM
!        READ(99,*) (M(I,J), J = 1,DIM)
!   ENDDO

! --- M IS HESSENBERG MATRIX --- !
    DIM = 100
    ALLOCATE(M(DIM,DIM),Q(DIM,DIM),R(DIM,DIM))
    M = 0
    CALL RANDOM_NUMBER(M)
    DO J = 1,DIM
        M(1,J) = 200*M(1,J) - 100
    ENDDO
    DO I = 2,DIM
        DO J = I-1,DIM
        M(I,J) = 200*M(I,J) - 100
        ENDDO
    ENDDO
    DO I = 3,DIM
        DO J = 1,I-2
            M(I,J) = 0
        ENDDO
    ENDDO

    CALL QRHOUSE_HESS(DIM,M,Q,R,.TRUE.)

    WRITE(98,*) 'QR DEC'
    WRITE(98,'(A)') 'INITIAL HESSENBERG MATRIX'
    DO I = 1,DIM
        WRITE(98,'(100F14.5)') (M(I,J),J = 1,DIM)
    ENDDO
    WRITE(98,*) '*******************'
    WRITE(98,'(A)') 'Q UNITARY HESSENBERG MATRIX'
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