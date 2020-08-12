MODULE INVMATMOD
    ! ------------------------------------------------------------------------------- !
    ! --- THIS MODULE CONTAIN ALL GLOBAL VARIABLE USEFULL FOR THE DIAGONALISATION --- !
    !-------------------------------------------------------------------------------- !
    IMPLICIT NONE 

    INTEGER, PARAMETER :: IN = 99, OUT = 98, ERR = 97 ! ---  USEFULL TO OPEN FILE --- !
    REAL*8,PARAMETER :: EPS0 = 1.D-14 ! --- NUMERICAL 0 --- !
    REAL*8, PARAMETER :: PI = 3.1415

    LOGICAL :: TESTLOG

    INTEGER :: DIM ! --- DIMENSION OF THE DIAGONALIZED MATRIX --- !
    INTEGER :: LINE, COLUMN ! --- USEFULL TO MANIPULATE LINE/COLUMN MATRIX ---!
    INTEGER :: I,J,K,BEGIN,END ! --- LOOP INDICES --- !

    REAL*8 :: NORME

END MODULE