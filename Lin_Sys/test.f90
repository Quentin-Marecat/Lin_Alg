PROGRAM TEST
    IMPLICIT NONE
    REAL*8,ALLOCATABLE :: M(:,:), VEC(:), SOL(:)
    CHARACTER :: TEXT
    LOGICAL :: CHECK
    INTEGER :: DIM,I,J,STAT
    OPEN(UNIT = 98,FILE = 'output')

!    OPEN(UNIT = 99,FILE = 'input',STATUS = 'OLD', ACTION= 'READ')
!    READ(99,*) DIM
!    READ(99,'(A)') TEXT
!    ALLOCATE(M(DIM,DIM),VEC(DIM),SOL(DIM))
!    DO I = 1,DIM
!        READ(99,*) (M(I,J), J = 1,DIM)
!    ENDDO
!    READ(99,'(A)') TEXT
!    READ(99,*) (VEC(I),I=1,DIM)

    DIM = 100
    ALLOCATE(M(DIM,DIM),VEC(DIM),SOL(DIM))
    M = 0
    CALL RANDOM_NUMBER(M)
    CALL RANDOM_NUMBER(VEC)
    DO I = 1,DIM
        VEC(I) = 200*VEC(I) - 100
        DO J = I,DIM
        M(I,J) = 200*M(I,J) - 100
        M(J,I) = M(I,J)
        ENDDO
    ENDDO

    CHECK = .TRUE.
!    CALL INV_MAT(DIM,M,VEC,SOL,CHECK)
    CALL LU_MAT(DIM,M,VEC,SOL,CHECK)
!    CALL QR_MAT(DIM,M,VEC,SOL,CHECK)
!    CALL JACOBI(DIM,M,VEC,SOL,CHECK)
!    CALL GAUSS_SEIDEL(DIM,M,VEC,SOL,CHECK)

    WRITE(98,'(A)') 'SOLUTION ='
    WRITE(98,'(100F14.5)') (SOL(I),I = 1,DIM)
    WRITE(98,*) '*******************'

    OPEN(UNIT = 97, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 2) WRITE(6,'(A)') 'RUN SUCCESFULLY, SEE output FILE' 
END