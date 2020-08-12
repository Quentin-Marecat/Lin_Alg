MODULE DIAGMATMOD
    ! --- FOR ANY REPORT OR SUGGESTION, PLEASE CONTACT quentin.marecat@etu.umontpellier.fr --- !
    ! ------------------------------------------------------------------------------- !
    ! --- THIS MODULE CONTAIN ALL GLOBAL VARIABLE USEFULL FOR THE DIAGONALISATION --- !
    !-------------------------------------------------------------------------------- !
    IMPLICIT NONE 

    INTEGER, PARAMETER :: IN = 99, OUT = 98, ERR = 97 ! ---  USEFULL TO OPEN FILE --- !
    INTEGER, PARAMETER :: MAXSTEP = 1D7, MAXSTEPHOUSE = 1000
    REAL*8,PARAMETER :: EPS0 = 1.D-10 ! --- NUMERICAL 0 --- !
    REAL*8, PARAMETER :: PI = 3.1415
    REAL*8,PARAMETER :: EPS = 1.D-9 ! --- EIGENVAL CONVERGENCE CRITERIUM --- !

    LOGICAL :: TESTLOG

    INTEGER :: DIM ! --- DIMENSION OF THE DIAGONALIZED MATRIX --- !
    INTEGER :: LINE, COLUMN ! --- USEFULL TO MANIPULATE LINE/COLUMN MATRIX ---!
    INTEGER :: I,J,K,BEGIN,END ! --- LOOP INDICES --- !

    REAL*8, ALLOCATABLE :: MATRIX(:,:) ! --- MATRIX TO DIAGONALIZE --- !
    REAL*8 :: NORME

END MODULE