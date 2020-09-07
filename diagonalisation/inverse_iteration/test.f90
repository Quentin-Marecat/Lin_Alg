PROGRAM TEST
    ! --- FOR ANY REPORT OR SUGGESTION, PLEASE CONTACT quentin.marecat@etu.umontpellier.fr --- !
    IMPLICIT NONE
    REAL*8 :: EPS0 = 1.D-15
    LOGICAL :: EIGENPRINT, SYM,CHECK
    REAL*8,ALLOCATABLE :: M1(:,:)
    REAL*8, ALLOCATABLE :: EV(:),EF(:,:), VEC(:)
    INTEGER :: DIM, STAT, NB_EV
    INTEGER :: I,J,IN, OUT
    IN =99
    OUT = 98
    OPEN(UNIT = OUT,FILE = 'output')

!    OPEN(UNIT = IN,FILE = 'input', STATUS= 'OLD', ACTION = 'READ')
!    READ(IN,*) DIM
!    ALLOCATE(M1(DIM,DIM),EV(DIM),EF(DIM,DIM))
!    ALLOCATE(MATTAMP(DIM,DIM),NOR(DIM),VEC(DIM))
!    DO I = 1,DIM
!        READ(IN,*) (M1(I,J),J=1,DIM)
!    ENDDO
!    CLOSE(IN)

    DIM = 100
    ALLOCATE(M1(DIM,DIM),EV(DIM),EF(DIM,DIM))
    ALLOCATE(VEC(DIM))
    M1 = 0
    CALL RANDOM_NUMBER(M1)
    DO I = 1,DIM
        DO J = I,DIM
        M1(I,J) = 200*M1(I,J) - 100
        M1(J,I) = M1(I,J)
        ENDDO
    ENDDO

    EIGENPRINT = .FALSE.
    CHECK = .TRUE.
    SYM = .TRUE.
    NB_EV = 1
!    CALL RQ_SHIFT_DIAG(DIM,M1,EV,EF,SYM,EIGENPRINT,CHECK)
    CALL PI_DIAG(DIM,M1,EV,EF,NB_EV,CHECK)
    CALL II_DIAG(DIM,M1,EV(NB_EV),VEC,SYM,CHECK)
    EIGENPRINT = .TRUE.

    WRITE(OUT,'(A)') '***********************'
    WRITE(OUT,'(A)') 'EIGENVALUE WANTED :'
    WRITE(OUT,'(100F14.5)') EV(1)
    WRITE(OUT,'(A)') '***********************'
    IF (EIGENPRINT) THEN
        IF (VEC(1)*EF(1,1) < EPS0) THEN
            DO J = 1,DIM
                VEC(J) = -VEC(J) 
            ENDDO
        ENDIF
        WRITE(OUT,'(A)') 'EIGENVECTOR FOUND / COMPUTED BEFORE: '
        DO J = 1,DIM
            WRITE(OUT,'(100F14.5)') VEC(J),EF(J,1)
        ENDDO
        WRITE(OUT,'(A)') '***********************'
    ENDIF
    CALL WRITE_DATA(DIM,M1,VEC,EV(NB_EV),OUT)

    OPEN(UNIT = 97, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 2) WRITE(6,'(A)') 'RUN SUCCESFULLY, SEE output FILE' 
END











SUBROUTINE WRITE_DATA(DIM,MAT,VEC,EVAL,OUT)
    ! ----------------------------------------------------------------------------- !
    ! --- THIS SUBROUTINE WRITE SUPPLEMENTARY DATA OF DIAGONALIZATION ALGORITHM --- !
    ! ----------------------------------------------------------------------------- !
    INTEGER, INTENT(IN) :: DIM, OUT
    REAL*8, INTENT(IN) :: MAT(DIM,DIM), VEC(DIM), EVAL
    REAL*8 :: NOR, VEC2(DIM)
    INTEGER :: I,J

    CALL MATAPPLI(DIM,MAT,VEC,VEC2)
    DO J = 1,DIM
        VEC2(J) = VEC2(J) - EVAL*VEC(J)
    ENDDO
    CALL NORMVEC(DIM,VEC2,NOR)
    WRITE(OUT,'(A)') 'MAT*EIGENVEC - EIGENVAL*EIGENVEC = 0'
    DO I = 1,DIM
        WRITE(OUT,'(100ES14.5)') (VEC2(I))
    ENDDO
    WRITE(OUT,'(A)') '***********************'
    WRITE(OUT,'(A)') 'NORME (MAT*EIGENVEC - EIGENVAL*EIGENVEC = 0)'
    WRITE(OUT,'(100ES14.5)') NOR
    WRITE(OUT,'(A)') '***********************'

END SUBROUTINE
