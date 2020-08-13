PROGRAM TEST
    ! --- FOR ANY REPORT OR SUGGESTION, PLEASE CONTACT quentin.marecat@etu.umontpellier.fr --- !
    IMPLICIT NONE
    REAL*8,ALLOCATABLE :: M1(:,:),M2(:,:),M3(:,:),MATROT(:,:)
    REAL*8, ALLOCATABLE :: EV(:),EF(:,:),ARRAY(:)
    INTEGER :: DIM, STAT
    INTEGER :: I,J,IN, OUT
    
    IN =99
    OUT = 98

!    DIM = 10
!    ALLOCATE(M1(DIM,DIM),M2(DIM,DIM),M3(DIM,DIM),MATROT(DIM,DIM),EV(DIM),EF(DIM,DIM),ARRAY(DIM))
!    M1 = 0
!    DO I = 1,DIM
!        DO J = 1,DIM
!        M1(I,J) = I*J
!        ENDDO
!    ENDDO
    OPEN(UNIT = OUT,FILE = 'output')
    OPEN(UNIT = IN,FILE = 'input', STATUS= 'OLD', ACTION = 'READ')
    READ(IN,*) DIM
    ALLOCATE(M1(DIM,DIM),M2(DIM,DIM),M3(DIM,DIM),MATROT(DIM,DIM),EV(DIM),EF(DIM,DIM))
    DO I = 1,DIM
        READ(IN,*) (M1(I,J),J=1,DIM)
    ENDDO
    CLOSE(IN)
!    CALL RQ_DIAG(DIM,M1,EV,EF,.TRUE.)
    CALL RQ_SHIFT_DIAG(DIM,M1,EV,EF,.TRUE.)
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
