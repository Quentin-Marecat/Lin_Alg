PROGRAM TEST
    IMPLICIT NONE
    REAL*8,ALLOCATABLE :: M(:,:), ORTHOM(:,:)
    INTEGER :: DIM,I,J, STAT
    OPEN(UNIT = 99,FILE = 'input',STATUS = 'OLD', ACTION= 'READ')
    OPEN(UNIT = 98,FILE = 'output')
    READ(99,*) DIM
!    DIM = 10
    ALLOCATE(M(DIM,DIM),ORTHOM(DIM,DIM))
    DO I = 1,DIM
        READ(99,*) (M(I,J), J = 1,DIM)
   ENDDO
    CALL GM(DIM,M,ORTHOM)

    DO I = 1,DIM
        WRITE(98,*) (ORTHOM(I,J),J = 1,DIM)
    ENDDO
    WRITE(98,*) '*******************'
    OPEN(UNIT = 97, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 2) WRITE(6,'(A)') 'RUN SUCCESFULLY, SEE output FILE' 
END