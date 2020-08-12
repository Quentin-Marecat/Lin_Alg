PROGRAM TEST
    IMPLICIT NONE
    REAL*8,ALLOCATABLE :: M(:,:), INVM(:,:)
    INTEGER :: DIM,I,J
    OPEN(UNIT = 99,FILE = 'input',STATUS = 'OLD', ACTION= 'READ')
    OPEN(UNIT = 98,FILE = 'output')
    READ(99,*) DIM
!    DIM = 10
    ALLOCATE(M(DIM,DIM),INVM(DIM,DIM))
    DO I = 1,DIM
        READ(99,*) (M(I,J), J = 1,DIM)
   ENDDO
    CALL INVMATQR(DIM,M,INVM)

    DO I = 1,DIM
        WRITE(6,*) (INVM(I,J),J = 1,DIM)
    ENDDO
    WRITE(6,*) '*******************'
END