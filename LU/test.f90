PROGRAM TEST
    IMPLICIT NONE 
    INTEGER :: DIM
    REAL*8, ALLOCATABLE :: M(:,:), L(:,:), U(:,:)
    INTEGER :: I,J, STAT
    OPEN(UNIT = 99,FILE = 'input',STATUS = 'OLD', ACTION= 'READ')
    OPEN(UNIT = 98,FILE = 'output')
    READ(99,*) DIM
!    DIM = 10
    ALLOCATE(M(DIM,DIM),L(DIM,DIM),U(DIM,DIM))
!    DO I = 1,DIM
!        DO J = 1,DIM
!            M(I,J) = I+J
!        ENDDO
!    ENDDO
    DO I = 1,DIM
        READ(99,*) (M(I,J), J = 1,DIM)
   ENDDO
    CALL LU(DIM,M,L,U)
    WRITE(98,'(A)') 'L ='
    DO I = 1,DIM
        WRITE(98,'(100F14.5)') (L(I,J),J = 1,DIM)
    ENDDO
    WRITE(98,*) '*******************'
    WRITE(98,'(A)') 'U ='
    DO I = 1,DIM
        WRITE(98,'(100F14.5)') (U(I,J),J = 1,DIM)
    ENDDO

    OPEN(UNIT = 97, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 2) WRITE(6,'(A)') 'RUN SUCCESFULLY, SEE output FILE' 
END