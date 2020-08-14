PROGRAM TEST
    ! --- FOR ANY REPORT OR SUGGESTION, PLEASE CONTACT quentin.marecat@etu.umontpellier.fr --- !
    IMPLICIT NONE
    REAL*8,ALLOCATABLE :: M1(:,:)
    REAL*8, ALLOCATABLE :: EV(:),EF(:,:)
    INTEGER :: DIM, STAT
    INTEGER :: I,J,IN, OUT
    
    IN =99
    OUT = 98
    OPEN(UNIT = OUT,FILE = 'output')

    DIM = 100
    ALLOCATE(M1(DIM,DIM),EV(DIM),EF(DIM,DIM))
    M1 = 0
    CALL RANDOM_NUMBER(M1)
    DO I = 1,DIM
        DO J = I,DIM
        M1(I,J) = 200*M1(I,J) - 100
        M1(J,I) = M1(I,J)
        ENDDO
    ENDDO
!    OPEN(UNIT = IN,FILE = 'input', STATUS= 'OLD', ACTION = 'READ')
!    READ(IN,*) DIM
!    ALLOCATE(M1(DIM,DIM),EV(DIM),EF(DIM,DIM))
!    DO I = 1,DIM
!        READ(IN,*) (M1(I,J),J=1,DIM)
!    ENDDO
!    CLOSE(IN)
    CALL JAC_DIAG(DIM,M1,EV,EF)

    WRITE(OUT,'(A)') '***********************'
    WRITE(OUT,'(A)') 'EIGENVALUES :'
    WRITE(OUT,'(100F14.5)') EV
    WRITE(OUT,'(A)') '***********************'
    WRITE(OUT,'(A)') 'EIGENVECTORS : '
    DO J = 1,DIM
        WRITE(OUT,'(100F14.5)') (EF(J,I),I=1,DIM)
    ENDDO
    WRITE(OUT,'(A)') '***********************'

    OPEN(UNIT = 97, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 2) WRITE(6,'(A)') 'RUN SUCCESFULLY, SEE output FILE' 
END
