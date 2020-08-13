SUBROUTINE GM(DIM,MATRIX,ORTHOMATRIX)
    ! --- FOR ANY REPORT OR SUGGESTION, PLEASE CONTACT quentin.marecat@etu.umontpellier.fr --- !
    ! ------------------------------------------------------------------------ !
    ! --- THIS SUBROUTINE ORTHOGONALIZE A MATRIX USING GRAM-SCHMIDT SCHEME --- !
    ! --- operation_mat.f90 is necessary ------------------------------------- !
    ! ------------------------------------------------------------------------ !
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: DIM
    REAL*8, INTENT(IN) :: MATRIX(DIM,DIM)
    REAL*8, INTENT(OUT) :: ORTHOMATRIX(DIM,DIM)
    REAL*8 :: ORTHOMATRIXT(DIM,DIM), MATTEST(DIM,DIM), ID(DIM,DIM)
    LOGICAL :: TEST
    REAL*8 :: PROD
    INTEGER :: I,J,K, STAT,ERR
    ERR = 97
    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) CLOSE(ERR,STATUS = 'delete')

    ORTHOMATRIX = MATRIX
    ! --- STEP 1 : DEFINITION OF THE FIRST VECTOR --- !
    CALL NORMVECMAT(DIM,ORTHOMATRIX,1)
    
    ! --- STEP 2 : OTHER VECTORS ARE ORTHOGONALE TO THE PREVIOUS ONE --- !
    DO I = 2,DIM
        DO J = 1,I-1
            CALL PRODVEC(DIM,MATRIX,ORTHOMATRIX,I,J,PROD)
            DO K = 1,DIM
                ORTHOMATRIX(K,I) = ORTHOMATRIX(K,I) - PROD*ORTHOMATRIX(K,J)
            ENDDO
        ENDDO
        CALL NORMVECMAT(DIM,ORTHOMATRIX,I)
    ENDDO

        ! --- VERIFICATION --- !
    TEST = .FALSE.
    CALL TRANSPOSE(DIM,ORTHOMATRIX,ORTHOMATRIXT)
    CALL PRODMAT(DIM,ORTHOMATRIX,ORTHOMATRIXT,MATTEST)
    ID = 0.
    DO I = 1,DIM
        ID(I,I) = 1
        DO J = 1,DIM
            IF (ABS(MATTEST(I,J)-ID(I,J)) > EPS0) TEST = .TRUE.
        ENDDO
    ENDDO
    IF (TEST) THEN
        OPEN(UNIT = ERR,FILE = 'error')
        WRITE(ERR,'(A)') 'PROBLEM GM DECOMPOSITION'
        WRITE(ERR,'(A)') 'L ='
        DO I = 1,DIM
            WRITE(ERR,'(100F14.5)') (ORTHOMATRIX(I,J), J=1,DIM)
        ENDDO
    ENDIF

    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) WRITE(6,'(A)') 'ERROR Gram-Schmidt, SEE FILE error'

END SUBROUTINE
    