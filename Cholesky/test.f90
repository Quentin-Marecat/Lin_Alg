PROGRAM TEST
    IMPLICIT NONE
    REAL*8,ALLOCATABLE :: M(:,:), L(:,:),MT(:,:)
    LOGICAL :: CHECK
    INTEGER :: DIM,I,J, STAT
    OPEN(UNIT = 99,FILE = 'input',STATUS = 'OLD', ACTION= 'READ')
    OPEN(UNIT = 98,FILE = 'output')
!    READ(99,*) DIM
!    ALLOCATE(M(DIM,DIM),L(DIM,DIM),MT(DIM,DIM))
!    DO I = 1,DIM
!        READ(99,*) (M(I,J), J = 1,DIM)
!   ENDDO
    DIM = 100
    ALLOCATE(M(DIM,DIM),L(DIM,DIM),MT(DIM,DIM))
    M = 0
    CALL RANDOM_NUMBER(M)
    DO I = 1,DIM
        DO J = I,DIM
        M(I,J) = 200*M(I,J)- 100
        ENDDO
    ENDDO

    CALL TRANSPOSE(DIM,M,MT)
    CALL PRODMAT(DIM,M,MT,M)

    CHECK = .TRUE.
    CALL CHOLESKY(DIM,M,L,CHECK)
!    CALL CHOLESKY_LU(DIM,M,L,CHECK)
!    CALL CHOLESKY_LDL(DIM,M,L,CHECK)
!    CALL CHOLESKY_DIAG(DIM,M,L,CHECK)
    WRITE(98,'(A)') 'CHOLESKY DIAG'
    WRITE(98,'(A)') 'L MATRIX :'
    DO I = 1,DIM
        WRITE(98,'(100F14.5)') (L(I,J),J = 1,DIM)
    ENDDO
    WRITE(98,*) '*******************'

    OPEN(UNIT = 97, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 2) WRITE(6,'(A)') 'RUN SUCCESFULLY, SEE output FILE' 
END