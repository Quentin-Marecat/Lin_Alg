PROGRAM TEST
    ! --- FOR ANY REPORT OR SUGGESTION, PLEASE CONTACT quentin.marecat@etu.umontpellier.fr --- !
    IMPLICIT NONE
    LOGICAL :: EIGENPRINT
    REAL*8,ALLOCATABLE :: M1(:,:)
    REAL*8, ALLOCATABLE :: EV(:),EF(:,:),MATTAMP(:,:),VEC(:),NOR(:)
    INTEGER :: DIM, STAT
    INTEGER :: I,J,IN, OUT
    
    IN =99
    OUT = 98
    OPEN(UNIT = OUT,FILE = 'output')

    DIM = 20
    ALLOCATE(M1(DIM,DIM),EV(DIM),EF(DIM,DIM))
    ALLOCATE(MATTAMP(DIM,DIM),NOR(DIM),VEC(DIM))
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
    EIGENPRINT = .TRUE.

!    CALL RQ_DIAG(DIM,M1,EV,EF,.TRUE.,EIGENPRINT,.TRUE.)
    CALL RQ_SHIFT_DIAG(DIM,M1,EV,EF,.TRUE.,EIGENPRINT,.TRUE.)

    WRITE(OUT,'(A)') '***********************'
    WRITE(OUT,'(A)') 'EIGENVALUES :'
    WRITE(OUT,'(100F14.5)') EV
    WRITE(OUT,'(A)') '***********************'
    IF (EIGENPRINT) THEN
        CALL PRODMAT(DIM,M1,EF,MATTAMP)
        DO I = 1,DIM
            DO J = 1,DIM
                VEC(J) = MATTAMP(J,I) - EV(I)*EF(J,I)
            ENDDO
            CALL NORMVEC(DIM,VEC,NOR(I))
        ENDDO
        WRITE(OUT,'(A)') 'EIGENVECTORS : '
        DO J = 1,DIM
            WRITE(OUT,'(100F14.5)') (EF(J,I),I=1,DIM)
        ENDDO
        WRITE(OUT,'(A)') '***********************'
        WRITE(OUT,'(A)') 'NORME (MAT*EIGENVEC - EIGENVAL*EIGENVEC = 0)'
        WRITE(OUT,'(100ES14.5)') (NOR(I),I=1,DIM)
        WRITE(OUT,'(A)') '***********************'
    ENDIF

    OPEN(UNIT = 97, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 2) WRITE(6,'(A)') 'RUN SUCCESFULLY, SEE output FILE' 
END
