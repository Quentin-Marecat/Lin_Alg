PROGRAM TEST
    IMPLICIT NONE
    REAL*8,ALLOCATABLE :: M(:,:), INVM(:,:)
    INTEGER :: DIM,I,J, STAT
    LOGICAL :: CHECK
    CHARACTER :: TEXT
    OPEN(UNIT = 98,FILE = 'output')

!    OPEN(UNIT = 99,FILE = 'input',STATUS = 'OLD', ACTION= 'READ')
!    READ(99,*) DIM
!    READ(99,'(A)') TEXT
!    ALLOCATE(M(DIM,DIM),INVM(DIM,DIM))
!    DO I = 1,DIM
!        READ(99,*) (M(I,J), J = 1,DIM)
!   ENDDO

    DIM = 100
    ALLOCATE(M(DIM,DIM),INVM(DIM,DIM))
    M = 0
    CALL RANDOM_NUMBER(M)
    DO I = 1,DIM
        DO J = I,DIM
        M(I,J) = 200*M(I,J) - 100
        M(J,I) = M(I,J)
        ENDDO
    ENDDO

! --- GAUSS_JACOBI VS GAUSS_SEIDEL TEST ---!
!    DIM = 100
!    ALLOCATE(M(DIM,DIM),INVM(DIM,DIM))
!    DO I = 1,DIM
!        DO J = 1,DIM
!            M(J,I) = 1/SQRT(DBLE(1+I+J))
!        ENDDO
!        M(I,I) = I
!    ENDDO

    CHECK = .TRUE.
    CALL INV_QR(DIM,M,INVM,CHECK)
!    CALL INV_LU(DIM,M,INVM,CHECK)
!    CALL INV_JAC(DIM,M,INVM,CHECK)
!    CALL INV_GS(DIM,M,INVM,CHECK)

    WRITE(98,'(A)') 'INVERSE OF THE MATRIX ='
    DO I = 1,DIM
        WRITE(98,'(100F14.5)') (INVM(I,J),J = 1,DIM)
    ENDDO
    WRITE(98,*) '*******************'
    OPEN(UNIT = 97, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 2) WRITE(6,'(A)') 'RUN SUCCESFULLY, SEE output FILE' 
END