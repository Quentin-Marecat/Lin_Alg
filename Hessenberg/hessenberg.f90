SUBROUTINE HESSENBERG_HOUSE(DIM,MATRIX,ROTMAT,HESS,VERIF)
    ! ---------------------------------------------------------------------------------------- !
    ! --- THIS SUBROUTINE GIVE THE QR DECOMPOSITION OF MAT WHERE R IS AN HESSENBERG MATRIX --- !
    ! --- THIS ALGORTITHME IS SIMILARY TO TRIGONALISATION METHOD ----------------------------- !
    ! --- MAT MUST BE SYMMETRIC -------------------------------------------------------------- !
    ! --- VERIF = .TRUE. IF VERIFICATIONS HAVE TO BE DONE ------------------------------------ !
    ! ---------------------------------------------------------------------------------------- !
    REAL*8, PARAMETER :: EPS0 = 1D-15, CRITOK = 1D-11
    INTEGER, INTENT(IN) :: DIM
    LOGICAL,INTENT(IN) :: VERIF
    LOGICAL :: TEST
    REAL*8, INTENT(IN) :: MATRIX(DIM,DIM)
    REAL*8, INTENT(OUT) :: ROTMAT(DIM,DIM),HESS(DIM,DIM)
    REAL*8 :: ID(DIM,DIM)
    REAL*8 :: VEC(DIM), ROTMATTAMP(DIM,DIM), FROEB1, FROEB2
    REAL*8 :: ROTMATTAMPT(DIM,DIM),MATTAMP(DIM,DIM)
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
        CALL PRODMAT_BLOC_QHESS(DIM,LOWESTELEM-1,ROTMATTAMP,HESS,HESS)
        CALL PRODMAT_BLOC_Q(DIM,LOWESTELEM,ROTMAT,ROTMATTAMP,ROTMAT)
    ENDDO
    ! --- VERIFICATION --- !
    IF (VERIF) THEN
        TEST = .FALSE.
        CALL NORME_FROEB(DIM,MATRIX,FROEB1)
        CALL NORME_FROEB(DIM,HESS,FROEB2)
        IF (ABS(FROEB1 - FROEB2) > CRITOK) TEST = .TRUE.

        IF (TEST) THEN
            OPEN(UNIT = ERR, FILE = 'error')
            WRITE(ERR,'(A)')  'PROBLEM HESSENBERG'
            WRITE(ERR,'(A)') 'FROEBIUS NORM'
            WRITE(ERR,'(A,4X,ES14.5,10X,A,4X,ES14.5)') 'FROEBIUS MATRIX =',FROEB1,'FROEBIUS TRIDIAG =',FROEB2
            WRITE(ERR,'(A,4X,ES14.5,10X,A,4X,ES14.5)') 'DELTA FROEBIUS =',FROEB1-FROEB2,'CRITERE ACPT =',CRITOK
            WRITE(ERR,'(A)') 'HESSENBERG MATRIX OBTAINED'
            DO I = 1,DIM
                WRITE(ERR,'(100E14.5)') (HESS(I,J), J=1,DIM)
            ENDDO
            WRITE(ERR,'(A)')  '****************'
            CALL PRODMAT(DIM,ROTMAT,HESS,MATTAMP)
            WRITE(ERR,'(A)') ' DIFFERENCE'
            DO I = 1,DIM
                WRITE(ERR,'(100E14.5)') (MATTAMP(I,J)-MATRIX(I,J), J=1,DIM)
            ENDDO
            WRITE(ERR,'(A)')  '****************'
        ENDIF
    ENDIF

    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) WRITE(6,'(A)') 'ERROR HESSENBERG HOUSEHOLDER, SEE FILE error'

END SUBROUTINE
