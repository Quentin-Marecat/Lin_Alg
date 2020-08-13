PROGRAM TEST
    IMPLICIT NONE
    REAL*8,ALLOCATABLE :: M(:,:), VEC(:), SOL(:)
    CHARACTER :: TEXT
    INTEGER :: DIM,I,J,STAT
    OPEN(UNIT = 99,FILE = 'input',STATUS = 'OLD', ACTION= 'READ')
    OPEN(UNIT = 98,FILE = 'output')
    READ(99,*) DIM
    READ(99,'(A)') TEXT
    ALLOCATE(M(DIM,DIM),VEC(DIM),SOL(DIM))
    DO I = 1,DIM
        READ(99,*) (M(I,J), J = 1,DIM)
    ENDDO
    READ(99,'(A)') TEXT
    READ(99,*) (VEC(I),I=1,DIM)
!    CALL INV_MAT(DIM,M,VEC,SOL)
    CALL LU_MAT(DIM,M,VEC,SOL)

    WRITE(98,'(A)') 'SOLUTION ='
    WRITE(98,'(100F14.5)') (SOL(I),I = 1,DIM)
    WRITE(98,*) '*******************'

    OPEN(UNIT = 97, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 2) WRITE(6,'(A)') 'RUN SUCCESFULLY, SEE output FILE' 
END