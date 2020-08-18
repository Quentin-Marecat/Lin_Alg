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




SUBROUTINE TRANSPOSE(DIM,M1,MT1)
  ! -------------------------------------------------------------- !
  ! --- THIS SUBROUTINE TRANSPOSE THE M1 MATRIX AND RETURN MT1 --- !
  ! --- M1 CAN BE THE SAME AS MT1 -------------------------------- !
  ! -------------------------------------------------------------- !
  IMPLICIT NONE 
  INTEGER, INTENT(IN) :: DIM
  REAL*8, INTENT(IN) :: M1(DIM,DIM)
  REAL*8, INTENT(OUT) :: MT1(DIM,DIM)
  REAL*8 :: MATTAMP(DIM,DIM)
  INTEGER :: LINE, COLUMN
  MATTAMP = M1
  MT1 = 0.
  DO LINE = 1,DIM
      DO COLUMN = 1,DIM
          MT1(LINE,COLUMN) = MATTAMP(COLUMN,LINE)
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




SUBROUTINE ORDERING(DIM,EIGENVAL,EIGENVECT)
  ! -------------------------------------------------------------- !
  ! --- THIS SUBROUTINE ORDER THE EIGENVALUES AND EIGENVECTORS --- !
  ! -------------------------------------------------------------- !
  INTEGER :: DIM
  REAL*8 :: EIGENVAL(DIM),EIGENVECT(DIM,DIM)
  REAL*8 :: EIGENVALTAMP(DIM),EIGENVECTTAMP(DIM,DIM), MINI
  INTEGER :: ORD(DIM)
  INTEGER :: I,J
  EIGENVALTAMP = EIGENVAL
  EIGENVECTTAMP = EIGENVECT
  DO I = 1,DIM-1
      ORD(I) = I
      MINI = EIGENVALTAMP(I)
      DO J = I+1,DIM
          IF (EIGENVALTAMP(J) < MINI) THEN
              ORD(I) = J
              MINI = EIGENVALTAMP(J)
          ENDIF
      ENDDO
      EIGENVALTAMP(I) = MINI
      EIGENVALTAMP(ORD(I)) = EIGENVAL(I)
      DO J = 1,DIM
          EIGENVECTTAMP(J,I) = EIGENVECT(J,ORD(I))
          EIGENVECTTAMP(J,ORD(I)) = EIGENVECT(J,I)
      ENDDO
      EIGENVAL = EIGENVALTAMP
      EIGENVECT = EIGENVECTTAMP
  ENDDO
  RETURN
END SUBROUTINE




SUBROUTINE NORME_FROEB(DIM,MATRIX,NORME)
  ! --------------------------------------------------------------- !
  ! --- THIS SUBROUTINE COMPUTE THE FROEBIUS NORM OF THE MATRIX --- !
  ! --------------------------------------------------------------- !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: DIM
  REAL*8, INTENT(IN) :: MATRIX(DIM,DIM)
  REAL*8, INTENT(OUT) :: NORME
  INTEGER :: I,J

  NORME = 0
  DO I = 1,DIM
    DO J = 1,DIM
      NORME = NORME + MATRIX(J,I)**2
    ENDDO
  ENDDO
  NORME = SQRT(NORME)
END SUBROUTINE






SUBROUTINE PRODMAT_BLOC_QR(DIM,BLOC,Q,RINT,R)
  ! --------------------------------------------------------------------------------- !
  ! --- THIS SUBROUTINE COMPUTE THE PRODUCT BETWEEN Q AND RINT ----------------------- !
  ! --- WHERE Q IS BLOCK DIAGONAL, Q = ((In,0),(0,FULL)) WITH DIM(IN)= BLOC*BLOC ---- !
  ! --- R IS BUILT TO BE UPPER MATRIX WITH 1 COLUMN OF ZERO ADDED TO RINT ----------- !
  ! --------------------------------------------------------------------------------- !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: DIM, BLOC
  REAL*8, INTENT(IN) :: Q(DIM,DIM), RINT(DIM,DIM)
  REAL*8, INTENT(OUT) :: R(DIM,DIM)
  INTEGER :: I,J,K
  REAL*8 :: RINTT(DIM,DIM)

  RINTT = RINT
  R = RINT

  IF (BLOC < DIM) THEN
    DO J = BLOC+1,DIM 
      DO I = BLOC+1,DIM 
        R(I,J) = 0
        IF(J /= BLOC+1 .OR. I == BLOC+1) THEN
          DO K = BLOC+1,DIM 
            R(I,J) = R(I,J) + Q(I,K)*RINTT(K,J)
          ENDDO
        ENDIF
      ENDDO
    ENDDO
  ENDIF
END SUBROUTINE


SUBROUTINE PRODMAT_BLOC_Q(DIM,BLOC,M1,Q,PROD)
  ! --------------------------------------------------------------------------------- !
  ! --- THIS SUBROUTINE COMPUTE THE PRODUCT BETWEEN M1 AND Q ----------------------- !
  ! --- WHERE Q IS BLOCK DIAGONAL, Q = ((In,0),(0,FULL)) WITH DIM(IN)= BLOC*BLOC ---- !
  ! --------------------------------------------------------------------------------- !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: DIM, BLOC
  REAL*8, INTENT(IN) :: M1(DIM,DIM), Q(DIM,DIM)
  REAL*8, INTENT(OUT) :: PROD(DIM,DIM)
  INTEGER :: I,J,K
  REAL*8 :: M1T(DIM,DIM)

  M1T = M1
  PROD = M1

  IF (BLOC < DIM) THEN
    DO J = BLOC+1,DIM 
      DO I = 1,DIM 
        PROD(I,J) = 0
          DO K = BLOC+1,DIM 
            PROD(I,J) = PROD(I,J) + M1T(I,K)*Q(K,J)
          ENDDO
      ENDDO
    ENDDO
  ENDIF
END SUBROUTINE
  






SUBROUTINE PRODMAT_BLOC_QR_TRI(DIM,BLOC,Q,RINT,R)
  ! --------------------------------------------------------------------------------- !
  ! --- THIS SUBROUTINE COMPUTE THE PRODUCT BETWEEN Q AND RINT ---------------------- !
  ! --- Q IS ID MATRIX WHERE ONLY 2*2 BLOC ELEMENT ON LOWESTELEM ARE MODIFIED ------- !
  ! --- R IS BUILT TO BE UPPER MATRIX WITH 1 COLUMN OF ZERO ADDED TO RINT ----------- !
  ! --------------------------------------------------------------------------------- !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: DIM, BLOC
  REAL*8, INTENT(IN) :: Q(DIM,DIM), RINT(DIM,DIM)
  REAL*8, INTENT(OUT) :: R(DIM,DIM)
  INTEGER :: I,J,K
  REAL*8 :: RINTT(DIM,DIM)

  RINTT = RINT
  R = RINT

  IF (BLOC < DIM) THEN
    DO J = BLOC+1,DIM
      DO I = BLOC+1,BLOC+2
        R(I,J) = 0
        IF(J /= BLOC+1 .OR. I == BLOC+1) THEN
          DO K = BLOC+1,BLOC+2 
            R(I,J) = R(I,J) + Q(I,K)*RINTT(K,J)
          ENDDO
        ENDIF
      ENDDO
    ENDDO
  ENDIF
END SUBROUTINE


SUBROUTINE PRODMAT_BLOC_Q_TRI(DIM,BLOC,M1,Q,PROD)
  ! ------------------------------------------------------------------------------------ !
  ! --- THIS SUBROUTINE COMPUTE THE PRODUCT BETWEEN M1 AND Q --------------------------- !
  ! --- M1 IS BLOCK DIAGONAL, M1 = ((HESS,0),(0,IN)) WHERE HESS IS AN HESSENBERG MAT --- !
  ! --- WITH DIM(HESS) = BLOC+2*BLOC+2
  ! --- WHERE Q IS BLOCK DIAGONAL, Q = ((In,0),(0,FULL)) WITH DIM(IN)= BLOC*BLOC ------- !
  ! ------------------------------------------------------------------------------------ !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: DIM, BLOC
  REAL*8, INTENT(IN) :: M1(DIM,DIM), Q(DIM,DIM)
  REAL*8, INTENT(OUT) :: PROD(DIM,DIM)
  INTEGER :: I,J,K
  REAL*8 :: M1T(DIM,DIM)

  M1T = M1
  PROD = M1

  IF (BLOC < DIM) THEN
    DO J = BLOC+1,BLOC+2
      DO I = 1,BLOC+2
        PROD(I,J) = 0
          DO K = BLOC+1,BLOC+2 
            PROD(I,J) = PROD(I,J) + M1T(I,K)*Q(K,J)
          ENDDO
      ENDDO
    ENDDO
  ENDIF
END SUBROUTINE
  
