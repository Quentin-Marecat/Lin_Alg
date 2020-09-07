PROGRAM TEST
    ! --- FOR ANY REPORT OR SUGGESTION, PLEASE CONTACT quentin.marecat@etu.umontpellier.fr --- !
    IMPLICIT NONE
    REAL*8,ALLOCATABLE :: M1(:,:)
    REAL*8, ALLOCATABLE :: EV(:),EF(:,:)
    REAL*8 :: SHIFT
    LOGICAL :: SYM
    INTEGER :: DIM, STAT,NB_EV
    INTEGER :: I,J,IN, OUT
    
    IN =99
    OUT = 98
    OPEN(UNIT = OUT,FILE = 'output')

!    OPEN(UNIT = IN,FILE = 'input', STATUS= 'OLD', ACTION = 'READ')
!    READ(IN,*) DIM
!    ALLOCATE(M1(DIM,DIM),EV(DIM),EF(DIM,DIM))
!    DO I = 1,DIM
!        READ(IN,*) (M1(I,J),J=1,DIM)
!    ENDDO
!    CLOSE(IN)

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

    NB_EV = 10 !!! --- NBR OF ENGVAL TO COMPUTE --- !!!
    SYM = .TRUE.
    SHIFT = 100
!    CALL PI_DIAG(DIM,M1,EV,EF,NB_EV,.TRUE.)
    CALL PI_INV_DIAG(DIM,M1,EV,EF,NB_EV,SYM,.TRUE.)
!    CALL PI_INV_DIAG(DIM,M1,SHIFT,EV,EF,NB_EV,SYM,.TRUE.)

    WRITE(OUT,'(A)') '***********************'
    WRITE(OUT,'(A)') 'EIGENVALUES :'
    WRITE(OUT,'(100F14.5)') (EV(I),I=1,NB_EV)
    WRITE(OUT,'(A)') '***********************'
    WRITE(OUT,'(A)') 'EIGENVECTORS : '
    DO J = 1,DIM
        WRITE(OUT,'(100F14.5)') (EF(J,I),I=1,NB_EV)
    ENDDO
    WRITE(OUT,'(A)') '***********************'
!    CALL WRITE_DATA(DIM,M1,EF,EV,NB_EV,OUT)

    OPEN(UNIT = 97, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 2) WRITE(6,'(A)') 'RUN SUCCESFULLY, SEE output FILE' 
END










SUBROUTINE WRITE_DATA(DIM,MAT,EIGENVEC,EIGENVAL,NB_EV,OUT)
    ! ----------------------------------------------------------------------------- !
    ! --- THIS SUBROUTINE WRITE SUPPLEMENTARY DATA OF DIAGONALIZATION ALGORITHM --- !
    ! ----------------------------------------------------------------------------- !
    INTEGER, INTENT(IN) :: DIM, OUT, NB_EV
    REAL*8, INTENT(IN) :: MAT(DIM,DIM), EIGENVEC(DIM,DIM), EIGENVAL(DIM)
    REAL*8 :: DIAG(DIM,DIM), MATTAMP(DIM,DIM), EIGENVECT(DIM,DIM)
    REAL*8 :: NOR(DIM), VEC(DIM)
    INTEGER :: I,J
    DIAG = 0
    DO I = 1,NB_EV
        DIAG(I,I) = EIGENVAL(I)
    ENDDO
    IF (NB_EV == DIM) THEN
        CALL PRODMAT(DIM,EIGENVEC,DIAG,MATTAMP)
        CALL TRANSPOSE(DIM,EIGENVEC,EIGENVECT)
        CALL PRODMAT(DIM,MATTAMP,EIGENVECT,MATTAMP)
        WRITE(OUT,'(A)') ' Q*D*Q(T) - MAT = 0'
        DO I = 1,DIM
            WRITE(OUT,'(100ES14.5)') (MATTAMP(I,J)-MAT(I,J), J=1,NB_EV)
        ENDDO
        WRITE(OUT,'(A)') '***********************'
    ENDIF
    CALL PRODMAT(DIM,MAT,EIGENVEC,MATTAMP)
    DO I = 1,NB_EV
        DO J = 1,DIM
            VEC(J) = MATTAMP(J,I) - EIGENVAL(I)*EIGENVEC(J,I)
        ENDDO
        CALL NORMVEC(DIM,VEC,NOR(I))
        ENDDO
    WRITE(OUT,'(A)') 'MAT*EIGENVEC - EIGENVAL*EIGENVEC = 0'
    DO I = 1,DIM
        WRITE(OUT,'(100ES14.5)')(MATTAMP(I,J) - EIGENVAL(J)*EIGENVEC(I,J),J=1,NB_EV)
    ENDDO
    WRITE(OUT,'(A)') '***********************'
    WRITE(OUT,'(A)') 'NORME (MAT*EIGENVEC - EIGENVAL*EIGENVEC = 0)'
    WRITE(OUT,'(100ES14.5)') (NOR(I),I=1,NB_EV)
    WRITE(OUT,'(A)') '***********************'

END SUBROUTINE