PROGRAM TEST
    IMPLICIT NONE
    CHARACTER :: TEXT
    REAL*8,ALLOCATABLE :: M(:,:), HESS(:,:), ROT(:,:)
    INTEGER :: DIM,I,J, STAT
    OPEN(UNIT = 98,FILE = 'output')

!    OPEN(UNIT = 99,FILE = 'input',STATUS = 'OLD', ACTION= 'READ')
!    READ(99,*) DIM
!    READ(99,'(A)') TEXT
!    ALLOCATE(M(DIM,DIM),TRI(DIM,DIM),ROT(DIM,DIM))
!    DO I = 1,DIM
!        READ(99,*) (M(I,J), J = 1,DIM)
!   ENDDO

    DIM = 100
    ALLOCATE(M(DIM,DIM),HESS(DIM,DIM),ROT(DIM,DIM))
    M = 0
    CALL RANDOM_NUMBER(M)
    DO I = 1,DIM
        DO J = I,DIM
        M(I,J) = 200*M(I,J) - 100
        M(J,I) = M(I,J)
        ENDDO
    ENDDO

    CALL HESSENBERG_HOUSE(DIM,M,ROT,HESS,.TRUE.)
    
    WRITE(98,*) 'HESSENBERG MATRIX'
    DO I = 1,DIM
        WRITE(98,'(100F14.5)') (HESS(I,J),J = 1,DIM)
    ENDDO
    WRITE(98,*) '*******************'

    OPEN(UNIT = 97, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 2) WRITE(6,'(A)') 'RUN SUCCESFULLY, SEE output FILE' 
END