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


SUBROUTINE NORMVECTOR(DIM,VEC,VECNORM)
  ! --- THIS SUBROUTINE NORMALIZE THE VECTOR VEC --- !
  IMPLICIT NONE
  INTEGER, INTENT(IN):: DIM
  REAL*8, INTENT(IN) :: VEC(DIM)
  REAL*8, INTENT(OUT) :: VECNORM(DIM)
  REAL*8 :: VECTAMP(DIM), NORM
  INTEGER :: I
  VECTAMP = VEC
  CALL NORMVEC(DIM,VECTAMP,NORM)
  DO I = 1,DIM
    VECNORM(I) = VECTAMP(I)/NORM
  ENDDO
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

SUBROUTINE MATAPPLI_TRI(DIM,ENDO,V0,V1)
  ! --------------------------------------------------------------------------- !
  ! --- THIS SUBROUTINE APPLY THE LINEAR APPLICATION ENDO TO V0 AND GIVE V1 --- !
  ! --- ENDO IS TRIDIAGONAL --------------------------------------------------- !
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
  DO J = 1,2
    V1(1) = V1(1) + ENDO(1,J)*VTAMP(J)
  ENDDO
  DO I = 2,DIM-1
      DO J = I-1,I+1
          V1(I) = V1(I) + ENDO(I,J)*VTAMP(J)
      ENDDO
  ENDDO
  DO J = DIM-1,DIM
    V1(DIM) = V1(DIM) + ENDO(DIM,J)*VTAMP(J)
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


SUBROUTINE PRODVEC2(DIM,VEC1,VEC2,PROD)
  ! ----------------------------------------------------------------------- !
  ! --- THIS SUBROUTINE COMPUTE THE VECTOR PRODUCT BETWEEN VEC1 ET VEC2 --- !
  ! ----------------------------------------------------------------------- !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: DIM
  REAL*8, INTENT(IN) :: VEC1(DIM), VEC2(DIM)
  REAL*8, INTENT(OUT) :: PROD
  INTEGER :: I

  PROD = 0.
  DO I = 1,DIM
    PROD = PROD + VEC1(I)*VEC2(I)
  ENDDO

END SUBROUTINE




SUBROUTINE ORDERING(DIM,EIGENVAL,EIGENVECT)
  ! -------------------------------------------------------------- !
  ! --- THIS SUBROUTINE ORDER THE EIGENVALUES AND EIGENVECTORS --- !
  ! -------------------------------------------------------------- !
  IMPLICIT NONE
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





SUBROUTINE ORDERING_0(DIM,EIGENVAL,EIGENVECT)
  ! -------------------------------------------------------------- !
  ! --- THIS SUBROUTINE ORDER THE EIGENVALUES AND EIGENVECTORS --- !
  ! --- 0 VALUE ARE SKIPED --------------------------------------- !
  ! -------------------------------------------------------------- !
  IMPLICIT NONE
  REAL*8, PARAMETER :: EPS0 = 1.D-15
  INTEGER, INTENT(IN) :: DIM
  REAL*8 :: EIGENVAL(DIM),EIGENVECT(DIM,DIM)
  REAL*8 :: EIGENVALTAMP(DIM),EIGENVECTTAMP(DIM,DIM), MINI
  INTEGER :: ORD(DIM)
  INTEGER :: I,J
  EIGENVALTAMP = EIGENVAL
  EIGENVECTTAMP = EIGENVECT
  DO I = 1,DIM-1
      ORD(I) = I
      MINI = EIGENVALTAMP(I)
      IF (ABS(MINI) < EPS0 .AND. ABS(EIGENVALTAMP(I+1)) < EPS0) EXIT
      DO J = I+1,DIM
        IF (ABS(EIGENVALTAMP(J)) < EPS0 ) EXIT
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



SUBROUTINE PRODMAT_BLOC_QHESS(DIM,BLOC,Q,RINT,R)
  ! --------------------------------------------------------------------------------- !
  ! --- THIS SUBROUTINE COMPUTE THE PRODUCT BETWEEN Q AND RINT ----------------------- !
  ! --- WHERE Q IS BLOCK DIAGONAL, Q = ((In,0),(0,FULL)) WITH DIM(IN)= BLOC*BLOC ---- !
  ! --- R IS BUILT TO BE HESSENBERG MATRIX WITH 1 COLUMN OF ZERO ADDED TO RINT ------ !
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
        IF(J /= BLOC+1 .OR. I == BLOC+1 .OR. I == BLOC+2) THEN
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







SUBROUTINE PRODMAT_RQ(DIM,R,Q,PROD)
  ! -------------------------------------------------------------------------------------- !
  ! --- THIS SUBROUTINE COMPUTE THE PRODUCT OF R*Q, WHERE R IS A UPPER DIAGONAL MATRIX --- !
  ! --- Q IS AN HESSENBERG MATRIX -------------------------------------------------------- !
  ! --- PROD IS TRIDIAGONAL -------------------------------------------------------------- !
  ! -------------------------------------------------------------------------------------- !
  INTEGER, INTENT(IN) :: DIM
  REAL*8, INTENT(IN) :: R(DIM,DIM), Q(DIM,DIM)
  REAL*8, INTENT(OUT) :: PROD(DIM,DIM)
  INTEGER :: I,J,K
  
  PROD = 0
  DO I = 1,DIM-1
    DO J = I,I+1
      DO K = I,J+1
        PROD(I,J) = PROD(I,J) + R(I,K)*Q(K,J)
      ENDDO
      PROD(J,I) = PROD(I,J)
    ENDDO
  ENDDO
  PROD(DIM,DIM) = R(DIM,DIM)*Q(DIM,DIM)

END SUBROUTINE




SUBROUTINE PRODMAT_RQ_HESS(DIM,R,Q,PROD)
  ! -------------------------------------------------------------------------------------- !
  ! --- THIS SUBROUTINE COMPUTE THE PRODUCT OF R*Q, WHERE R IS A UPPER DIAGONAL MATRIX --- !
  ! --- PROD AND Q ARE HESSENBERG -------------------------------------------------------- !
  ! -------------------------------------------------------------------------------------- !
  INTEGER, INTENT(IN) :: DIM
  REAL*8, INTENT(IN) :: R(DIM,DIM), Q(DIM,DIM)
  REAL*8, INTENT(OUT) :: PROD(DIM,DIM)
  INTEGER :: I,J,K
  
  PROD = 0
  DO J = 1,DIM-1
    DO K = 1,J+1
      PROD(1,J) = PROD(1,J) + R(1,K)*Q(K,J)
    ENDDO
  ENDDO
  DO K = 1,DIM
    PROD(1,DIM) = PROD(1,DIM) + R(1,K)*Q(K,DIM)
  ENDDO
  DO I = 2,DIM-1
    DO J = I-1,DIM
      DO K = I,DIM
        PROD(I,J) = PROD(I,J) + R(I,K)*Q(K,J)
      ENDDO
      PROD(J,I) = PROD(I,J)
    ENDDO
  ENDDO
  DO J = DIM-1,DIM
    PROD(DIM,J) = R(DIM,DIM)*Q(DIM,J)
  ENDDO

END SUBROUTINE






SUBROUTINE HOUSESTEP(DIM,VEC,LOWESTELEM,MATROT)
  ! -------------------------------------------------------------------------------------- !
  ! --- THIS SUBROUTINE APPLIED A HOUSEHOLDER STEP ON THE VECTOR VEC AND GIVE VECHOUSE --- !
  ! --- VEC AND VECHOUSE CAN BE THE SAME BUT VEC WILL BE ERASED -------------------------- !
  ! --- INPUT  _ VEC : VECTOR TO PROCESS THE HOUSEHOLDER PROCESS
  ! ---        _ LOWESTELEM : LOWEST POSITION NON ZERO ELEMENT AFTER HOUSEHOLDER PROCESS, 2 <= LOWESTELEM <= DIM1 - 1 
  ! --- OUTPUT _ MATROT : ROTATION MATRIX : MATROT*VEC = VECHOUSE
  ! -------------------------------------------------------------------------------------- !
  IMPLICIT NONE
  REAL*8, PARAMETER :: EPS0 = 1.D-14
  LOGICAL :: TESTLOG
  INTEGER, INTENT(IN) :: DIM
  REAL*8, INTENT(IN) :: VEC(DIM)
  INTEGER, INTENT(IN) :: LOWESTELEM
  REAL*8, INTENT(OUT) :: MATROT(DIM,DIM)
  REAL*8 :: ID(DIM,DIM),TAMPVEC(DIM),VECTAMP2, NORME
  INTEGER :: I,J,K, ERR
  ID = 0.
  DO I = 1,DIM
      ID(I,I) = 1.
  ENDDO
  IF (LOWESTELEM < 1 .OR. LOWESTELEM > (DIM)) THEN
      ERR = 97
      WRITE(6,'(A)') 'ERROR HOUSESTEP DIMENSION'
      OPEN(UNIT = ERR,FILE='error',ACTION = 'WRITE')
      WRITE(ERR,'(A,10X,A,I2,10X,A,4X,I2)')  'ERROR HOUSESTEP DIMENSION','DIM =',DIM,'ELEM =',LOWESTELEM
  ELSEIF (LOWESTELEM == DIM) THEN
      MATROT = ID
  ELSE 
      TESTLOG = .TRUE.
      DO I = LOWESTELEM + 1, DIM  ! --- TEST IF THE VECTOR ALREADY HAVE THE GOOD SHAPE 
          IF (ABS(VEC(I)) > EPS0) TESTLOG = .FALSE. 
      ENDDO
      IF (TESTLOG) THEN
          MATROT = ID

      ELSE
          NORME = 0.
          DO I = LOWESTELEM,DIM 
              NORME = NORME + VEC(I)**2
          ENDDO
          TAMPVEC = 0.
          DO I = LOWESTELEM,DIM
              TAMPVEC(I) = SQRT(NORME)*ID(I,LOWESTELEM) - VEC(I)
          ENDDO
          NORME = 0.
          DO I = LOWESTELEM,DIM
              NORME = NORME + TAMPVEC(I)**2
          ENDDO
          MATROT = ID
          DO I = LOWESTELEM,DIM
              DO J = LOWESTELEM,DIM
                  MATROT(I,J) = MATROT(I,J) - (2/NORME)*TAMPVEC(I)*TAMPVEC(J) 
              ENDDO
          ENDDO
      ENDIF
  ENDIF
END SUBROUTINE





SUBROUTINE HOUSESTEP_TRI(DIM,VEC,LOWESTELEM,MATROT)
  ! -------------------------------------------------------------------------------------- !
  ! --- THIS SUBROUTINE APPLIED A HOUSEHOLDER STEP ON THE VECTOR VEC AND GIVE VECHOUSE --- !
  ! --- THE MATRIX MUST BE TRIGONAL ------------------------------------------------------ !
  ! --- VEC AND VECHOUSE CAN BE THE SAME BUT VEC WILL BE ERASED -------------------------- !
  ! --- INPUT  _ VEC : VECTOR TO PROCESS THE HOUSEHOLDER PROCESS
  ! ---        _ LOWESTELEM : LOWEST POSITION NON ZERO ELEMENT AFTER HOUSEHOLDER PROCESS, 2 <= LOWESTELEM <= DIM1 - 1 
  ! --- OUTPUT _ MATROT : ROTATION MATRIX : MATROT*VEC = VECHOUSE
  ! -------------------------------------------------------------------------------------- !
  IMPLICIT NONE
  REAL*8, PARAMETER :: EPS0 = 1.D-15
  LOGICAL :: TESTLOG
  INTEGER, INTENT(IN) :: DIM
  REAL*8, INTENT(IN) :: VEC(DIM)
  INTEGER, INTENT(IN) :: LOWESTELEM
  REAL*8, INTENT(OUT) :: MATROT(DIM,DIM)
  REAL*8 :: ID(DIM,DIM),TAMPVEC(DIM),VECTAMP2, NORME
  INTEGER :: I,J,K, ERR
  ID = 0.
  DO I = 1,DIM
      ID(I,I) = 1.
  ENDDO
  IF (LOWESTELEM < 1 .OR. LOWESTELEM > (DIM)) THEN
      ERR = 97
      WRITE(6,'(A)') 'ERROR HOUSESTEP DIMENSION'
      OPEN(UNIT = ERR,FILE='error',ACTION = 'WRITE')
      WRITE(ERR,'(A,10X,A,I2,10X,A,4X,I2)')  'ERROR HOUSESTEP DIMENSION','DIM =',DIM,'ELEM =',LOWESTELEM
      CLOSE(ERR)
  ELSEIF (LOWESTELEM == DIM) THEN
      MATROT = ID
  ELSE 
      IF (ABS(VEC(LOWESTELEM + 1)) < EPS0) THEN
          MATROT = ID

      ELSE
          NORME = 0.
          DO I = LOWESTELEM,LOWESTELEM + 1
              NORME = NORME + VEC(I)**2
          ENDDO
          TAMPVEC = 0.
          DO I = LOWESTELEM,LOWESTELEM + 1
              TAMPVEC(I) = SQRT(NORME)*ID(I,LOWESTELEM) - VEC(I)
          ENDDO
          NORME = 0.
          DO I = LOWESTELEM,LOWESTELEM + 1
              NORME = NORME + TAMPVEC(I)**2
          ENDDO
          MATROT = ID
          DO I = LOWESTELEM,LOWESTELEM + 1
              DO J = LOWESTELEM,LOWESTELEM + 1
                  MATROT(I,J) = MATROT(I,J) - (2/NORME)*TAMPVEC(I)*TAMPVEC(J) 
              ENDDO
          ENDDO
      ENDIF
  ENDIF
END SUBROUTINE

  
