PROGRAM TEST
    IMPLICIT NONE
    REAL*8,ALLOCATABLE :: M(:,:), Q(:,:), R(:,:)
    INTEGER :: DIM,I,J, STAT
    OPEN(UNIT = 99,FILE = 'input',STATUS = 'OLD', ACTION= 'READ')
    OPEN(UNIT = 98,FILE = 'output')
    READ(99,*) DIM
    ALLOCATE(M(DIM,DIM),Q(DIM,DIM),R(DIM,DIM))
    DO I = 1,DIM
        READ(99,*) (M(I,J), J = 1,DIM)
   ENDDO
    CALL QRGRAM(DIM,M,Q,R)
    WRITE(98,*) 'QR GRAM'
    DO I = 1,DIM
        WRITE(98,'(100F14.5)') (Q(I,J),J = 1,DIM)
    ENDDO
    WRITE(98,*) '*******************'
    DO I = 1,DIM
        WRITE(98,'(100F14.5)') (R(I,J),J = 1,DIM)
    ENDDO
    WRITE(98,*) '*******************'

    CALL QRHOUSE(DIM,M,Q,R)
    WRITE(98,*) 'QR HOUSE'
    DO I = 1,DIM
        WRITE(98,'(100F14.5)') (Q(I,J),J = 1,DIM)
    ENDDO
    WRITE(98,*) '*******************'
    DO I = 1,DIM
        WRITE(98,'(100F14.5)') (R(I,J),J = 1,DIM)
    ENDDO
    WRITE(98,*) '*******************'

    OPEN(UNIT = 97, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 2) WRITE(6,'(A)') 'RUN SUCCESFULLY, SEE output FILE' 
END