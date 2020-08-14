PROGRAM TEST
    IMPLICIT NONE
    REAL*8,ALLOCATABLE :: M(:,:), INVM(:,:)
    INTEGER :: DIM,I,J, STAT
    CHARACTER :: TEXT
    OPEN(UNIT = 99,FILE = 'input',STATUS = 'OLD', ACTION= 'READ')
    OPEN(UNIT = 98,FILE = 'output')
    READ(99,*) DIM
    READ(99,'(A)') TEXT
    ALLOCATE(M(DIM,DIM),INVM(DIM,DIM))
    DO I = 1,DIM
        READ(99,*) (M(I,J), J = 1,DIM)
   ENDDO
!    CALL INVQR(DIM,M,INVM)
    CALL INVLU(DIM,M,INVM)
!    CALL INVJAC(DIM,M,INVM)

    WRITE(98,'(A)') 'INVERSE OF THE MATRIX ='
    DO I = 1,DIM
        WRITE(98,'(100F14.5)') (INVM(I,J),J = 1,DIM)
    ENDDO
    WRITE(98,*) '*******************'
    OPEN(UNIT = 97, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 2) WRITE(6,'(A)') 'RUN SUCCESFULLY, SEE output FILE' 
END