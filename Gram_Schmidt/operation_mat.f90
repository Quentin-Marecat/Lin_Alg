SUBROUTINE PRODMAT(DIM,M1,M2,MPROD)
    ! -------------------------------------------------------------------- !
    ! --- THIS SUBROUTINE COMPUTE MPROD, THE PRODUCT BETWEEN M1 AND M2 --- !
    ! --- M1/M2 CAN BE THE SAME AS MPROD BUT IT IS ERASED -----------------!
    ! -------------------------------------------------------------------- !
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: DIM
    REAL*8, INTENT(IN) :: M1(DIM,DIM), M2(DIM,DIM)
    REAL*8, INTENT(OUT) :: MPROD(DIM,DIM)
    REAL*8 :: M1TAMP(DIM,DIM), M2TAMP(DIM,DIM),MPRODTAMP(DIM,DIM)
    INTEGER :: INCR, LINE, COLUMN, I, J, K
    MPRODTAMP = 0.
    M1TAMP = M1
    M2TAMP = M2
    DO LINE = 1,DIM 
        DO COLUMN = 1,DIM
          DO K = 1,DIM
            MPRODTAMP(LINE,COLUMN) = MPRODTAMP(LINE,COLUMN) + M1TAMP(LINE,K)*M2TAMP(K,COLUMN)
          ENDDO
        ENDDO
    ENDDO
    MPROD = MPRODTAMP
END SUBROUTINE



SUBROUTINE MATAPPLI(DIM,ENDO,V0,V1)
  ! --------------------------------------------------------------------------- !
  ! --- THIS SUBROUTINE APPLY THE LINEAR APPLICATION ENDO TO V0 AND GIVE V1 --- !
  ! --- V0 AND V1 CAN BE THE SAME BUT V0 IS ERASED ---------------------------- !
  ! --------------------------------------------------------------------------- !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: DIM
  REAL*8,INTENT(IN) :: ENDO(DIM,DIM),V0(DIM)
  REAL*8,INTENT(OUT) :: V1(DIM)
  REAL*8 :: VTAMP(DIM)
  INTEGER :: I, J
  VTAMP = V0
  V1 = 0
  DO I = 1,DIM
      DO J = 1,DIM
          V1(I) = V1(I) + ENDO(I,J)*VTAMP(J)
      ENDDO
  ENDDO
END SUBROUTINE 




SUBROUTINE NORMVEC(DIM,VEC,NOR)
  ! ----------------------------------------------------------- !
  ! --- THIS SUBROUTINE COMPUTE THE NORME OF THE VECTOR VEC --- !
  ! ----------------------------------------------------------- !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: DIM
  REAL*8, INTENT(IN) :: VEC(DIM)
  REAL*8, INTENT(OUT) :: NOR
  INTEGER :: I

  NOR = 0.
  DO I = 1,DIM
    NOR = NOR + VEC(I)**2
  ENDDO
  NOR = SQRT(NOR)

END SUBROUTINE



SUBROUTINE NORMVECMAT(DIM,MAT,NUMBER)
  ! --------------------------------------------------------------- !
  ! --- THIS SUBROUTINE NORMALIZE THE I-TH VECTOR OF THE MATRIX --- !
  ! --------------------------------------------------------------- !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: DIM,NUMBER
  INTEGER :: I
  REAL*8 :: NOR, MAT(DIM,DIM)

  NOR = 0.
  DO I = 1,DIM
    NOR = NOR + MAT(I,NUMBER)**2
  ENDDO
  NOR = SQRT(NOR)
  DO I = 1,DIM
    MAT(I,NUMBER) = MAT(I,NUMBER)/NOR
  ENDDO
END SUBROUTINE


SUBROUTINE PRODVEC(DIM,MAT1,MAT2,POS1,POS2,PROD)
  ! ----------------------------------------------------------------------- !
  ! --- THIS SUBROUTINE COMPUTE THE VECTOR PRODUCT BETWEEN VEC1 ET VEC2 --- !
  ! ----------------------------------------------------------------------- !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: DIM,POS1, POS2
  REAL*8, INTENT(IN) :: MAT1(DIM,DIM), MAT2(DIM,DIM)
  REAL*8, INTENT(OUT) :: PROD
  INTEGER :: I

  PROD = 0.
  DO I = 1,DIM
    PROD = PROD + MAT1(I,POS1)*MAT2(I,POS2)
  ENDDO

END SUBROUTINE
