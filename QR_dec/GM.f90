SUBROUTINE GM(DIM,MATRIX,ORTHOMATRIX)
    ! ------------------------------------------------------------------------ !
    ! --- THIS SUBROUTINE ORTHOGONALIZE A MATRIX USING GRAM-SCHMIDT SCHEME --- !
    ! ------------------------------------------------------------------------ !
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: DIM
    REAL*8, INTENT(IN) :: MATRIX(DIM,DIM)
    REAL*8, INTENT(OUT) :: ORTHOMATRIX(DIM,DIM)
    REAL*8 :: PROD
    INTEGER :: I,J,K

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

END SUBROUTINE
    