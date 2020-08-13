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