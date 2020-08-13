SUBROUTINE INV_MAT(DIM,MATRIX,VEC,SOL)
! --- FOR ANY REPORT OR SUGGESTION, PLEASE CONTACT quentin.marecat@etu.umontpellier.fr --- !
    ! ----------------------------------------------------------------------------------------------- !
    ! --- THIS SUBROUTINE FIND THE SOLUTION OF THE SYSTEM MATRIX*SOL = VEC USING THE INVERSE OF A --- !
    ! --- MATRIX MUST BE SQUARE AND INVERSIBLE ------------------------------------------------------ !
    ! --- THIS SOLUTION IS NOT THE BEST ONE BECAUSE IT INVOLVE SOME USELESS OPERATION --------------- !
    ! ----------------------------------------------------------------------------------------------- !
    IMPLICIT NONE 
    REAL*8, PARAMETER :: EPS0 = 1.D-14
    INTEGER, INTENT(IN) :: DIM
    REAL*8, INTENT(IN) :: MATRIX(DIM,DIM),VEC(DIM)
    REAL*8, INTENT(OUT) :: SOL(DIM)
    REAL*8 :: INVMATRIX(DIM,DIM), VECTEST(DIM)
    LOGICAL :: TEST
    INTEGER :: I,ERR, STAT

    ERR = 97
    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) CLOSE(ERR,STATUS = 'delete')

    ! ---------------------------------------------------------------------------------------------------- !
    CALL INVMATQR(DIM,MATRIX,INVMATRIX)
    ! --- INVERSION SUBROUTINE FROM Quentin-Marecat PACKAGES IS USED, PLEASE FIND IT ON MY GITHUB PAGE --- !
    OPEN(UNIT = 97, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 2) THEN CALL MATAPPLI(DIM,INVMATRIX,VEC,SOL)

    ! --- VERIFICATION --- !
    CALL MATAPPLI(DIM,MATRIX,SOL,VECTEST)
    TEST = .FALSE.
    DO I = 1,DIM
        IF (ABS(VECTEST(I) - VEC(I)) > EPS0) TEST = .TRUE.
    ENDDO
    IF (TEST) THEN
        OPEN(UNIT = ERR,FILE = 'error')
        WRITE(ERR,'(A)') 'PROBLEM SOLUTION'
        WRITE(ERR,'(A)') 'x = '
        WRITE(ERR,'(100F14.5)') (SOL(I), I = 1,DIM)
        WRITE(ERR,'(A)') '************'
    ENDIF

    OPEN(UNIT = ERR, FILE = 'error', IOSTAT = STAT, STATUS = 'old')
    IF (STAT == 0) WRITE(6,'(A)') 'ERROR LINEAR SOLUTION, SEE FILE error'

END SUBROUTINE
