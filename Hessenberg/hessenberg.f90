SUBROUTINE HESSENBERG_HOUSE(DIM,MATRIX,ROTMAT,HESS,CHECK)
    ! ---------------------------------------------------------------------------------------- !
    ! --- THIS SUBROUTINE GIVE THE QR DECOMPOSITION OF MAT WHERE R IS AN HESSENBERG MATRIX --- !
    ! --- THIS ALGORTITHME IS SIMILARY TO TRIGONALISATION METHOD ----------------------------- !
    ! --- MAT MUST BE SYMMETRIC -------------------------------------------------------------- !
    ! --- CHECK = .TRUE. IF VERIFICATIONS HAVE TO BE DONE ------------------------------------ !
    ! ---------------------------------------------------------------------------------------- !
    LOGICAL,INTENT(IN) :: CHECK
    REAL*8, PARAMETER :: EPS0 = 1D-15, EPS = 1.D-12
    INTEGER, INTENT(IN) :: DIM
    LOGICAL :: TEST
    REAL*8, INTENT(IN) :: MATRIX(DIM,DIM)
    REAL*8, INTENT(OUT) :: ROTMAT(DIM,DIM),HESS(DIM,DIM)
    REAL*8 :: ID(DIM,DIM)
    REAL*8 :: VEC(DIM), ROTMATTAMP(DIM,DIM),MATTEST(DIM,DIM)
    INTEGER :: LOWESTELEM, I,J,ERR, STAT

    ERR = 97
    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) CLOSE(ERR,STATUS = 'delete')

    HESS = MATRIX

    ID = 0.
    DO I = 1,DIM
        ID(I,I) = 1.
    ENDDO
    ROTMAT = ID
    DO LOWESTELEM = 1,DIM-2
        DO J = 1,DIM
            VEC(J) = HESS(J,LOWESTELEM)
        ENDDO 
        CALL HOUSESTEP(DIM,VEC,LOWESTELEM+1,ROTMATTAMP)
        CALL PRODMAT_BLOC_QHESS(DIM,LOWESTELEM-1,ROTMATTAMP,HESS,HESS) ! --- MATRIX PRODUCT IMPROVEMENT --- !
        CALL PRODMAT_BLOC_Q(DIM,LOWESTELEM,ROTMAT,ROTMATTAMP,ROTMAT) ! --- MATRIX PRODUCT IMPROVEMENT --- !
    ENDDO
    ! --- VERIFICATION --- !
    IF (CHECK) THEN
        TEST = .FALSE.
        CALL PRODMAT(DIM,ROTMAT,HESS,MATTEST)
        DO I = 1,DIM
            DO J = 1,DIM
                IF (ABS(MATTEST(I,J)-MATRIX(I,J)) > EPS) TEST = .TRUE.
            ENDDO
        ENDDO

        IF (TEST) THEN
            OPEN(UNIT = ERR,FILE = 'error')
            WRITE(ERR,'(A)') 'PROBLEM HESSENBERG DECOMPOSITION'
            WRITE(ERR,'(A)') 'ROTATION MATRIX ='
            DO I = 1,DIM
                WRITE(ERR,'(100F14.5)') (ROTMAT(I,J), J=1,DIM)
            ENDDO
            WRITE(ERR,'(A)')'***************'
            WRITE(ERR,'(A)') 'HESSENBERG MATRIX ='
            DO I = 1,DIM
                WRITE(ERR,'(100F14.5)') (HESS(I,J), J=1,DIM)
            ENDDO
            WRITE(ERR,'(A)')'***************'
            WRITE(ERR,'(A,4X,ES14.5)')'EPS',EPS
            WRITE(ERR,'(A)') 'MATRIX - ROTMAT*HESSENBERG ='
            DO I = 1,DIM
                WRITE(ERR,'(100ES14.5)') (MATRIX(I,J) - MATTEST(I,J), J=1,DIM)
            ENDDO
            WRITE(ERR,'(A)')'***************'
        ENDIF
    ENDIF

    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) WRITE(6,'(A)') 'ERROR HESSENBERG HOUSEHOLDER, SEE FILE error'

END SUBROUTINE
