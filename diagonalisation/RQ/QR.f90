SUBROUTINE QR(MAT,Q,R)
    USE DIAGMATMOD
    ! ------------------------------------------------------------------- !
    ! --- THIS SUBROUTINE GIVE THE QR DECOMPOSITION OF THE MATRIX MAT --- !
    ! --- METHOD : HOUSEHOLDER DECOMPOSITION ---------------------------- !
    ! ------------------------------------------------------------------- !
    IMPLICIT NONE
    REAL*8,INTENT(IN) :: MAT(DIM,DIM)
    REAL*8, INTENT(OUT) :: Q(DIM,DIM),R(DIM,DIM)
    CALL HOUSEHOLDER(MAT,R,Q)
END SUBROUTINE