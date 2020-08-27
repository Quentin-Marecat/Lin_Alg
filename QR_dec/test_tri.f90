PROGRAM TEST
    IMPLICIT NONE
    REAL*8,ALLOCATABLE :: M(:,:), Q(:,:), R(:,:)
    INTEGER :: DIM,I,J, STAT
    REAL*8 :: X
    OPEN(UNIT = 98,FILE = 'output')

!    OPEN(UNIT = 99,FILE = 'input',STATUS = 'OLD', ACTION= 'READ')
!    READ(99,*) DIM
!    ALLOCATE(M(DIM,DIM),Q(DIM,DIM),R(DIM,DIM))
!    DO I = 1,DIM
!        READ(99,*) (M(I,J), J = 1,DIM)
!   ENDDO

! --- M IS TRIGONAL SYMETRIC --- !
    DIM = 100
    ALLOCATE(M(DIM,DIM),Q(DIM,DIM),R(DIM,DIM))
    M = 0
    DO I = 1,DIM-1
        CALL RANDOM_NUMBER(X)
        M(I,I) = 200*X - 100
        CALL RANDOM_NUMBER(X)
        M(I+1,I) = 200*X - 100
        M(I,I+1) = M(I+1,I)
    ENDDO
    CALL RANDOM_NUMBER(X)
    M(DIM,DIM) = 200*X - 100

    CALL QRHOUSE_TRI(DIM,M,Q,R,.TRUE.)

    WRITE(98,*) 'QR DEC'
    WRITE(98,'(A)') 'INITIAL TRIDIAG MATRIX'
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