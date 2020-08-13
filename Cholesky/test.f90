PROGRAM TEST
    IMPLICIT NONE
    REAL*8,ALLOCATABLE :: M(:,:), L(:,:)
    INTEGER :: DIM,I,J, STAT
    OPEN(UNIT = 99,FILE = 'input',STATUS = 'OLD', ACTION= 'READ')
    OPEN(UNIT = 98,FILE = 'output')
    READ(99,*) DIM
!    DIM = 10
    ALLOCATE(M(DIM,DIM),L(DIM,DIM))
    DO I = 1,DIM
        READ(99,*) (M(I,J), J = 1,DIM)
   ENDDO
    CALL CHOLESKY(DIM,M,L)
!    CALL CHOLESKY_LU(DIM,M,L)
!    CALL CHOLESKY_LDL(DIM,M,L)
!   CALL CHOLESKY_DIAGONALISATION(DIM,M,L)
    WRITE(98,'(A)') 'CHOLESKY DIAG'
    WRITE(98,'(A)') 'L MATRIX :'
    DO I = 1,DIM
        WRITE(98,'(100F14.5)') (L(I,J),J = 1,DIM)
    ENDDO
    WRITE(98,*) '*******************'
    OPEN(UNIT = 97, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 2) WRITE(6,'(A)') 'RUN SUCCESFULLY, SEE output FILE' 
END