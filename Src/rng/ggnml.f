C   IMSL ROUTINE NAME   - GGNML                                         GGNH0010
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - DGC/SINGLE
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - NORMAL OR GAUSSIAN RANDOM DEVIATE GENERATOR
C
C   USAGE               - CALL GGNML (DSEED,NR,R)
C
C   ARGUMENTS    DSEED  - INPUT/OUTPUT DOUBLE PRECISION VARIABLE
C                           ASSIGNED AN INTEGER VALUE IN THE
C                           EXCLUSIVE RANGE (1.D0, 2147483647.D0).
C                           DSEED IS REPLACED BY A NEW VALUE TO BE
C                           USED IN A SUBSEQUENT CALL.
C                NR     - INPUT NUMBER OF DEVIATES TO BE GENERATED.
C                R      - OUTPUT VECTOR OF LENGTH NR CONTAINING THE
C                           NORMAL (0,1) RANDOM NUMBERS.
C
C   PRECISION/HARDWARE  - SINGLE/ALL
C
C   REQD. IMSL ROUTINES - GGUBS,MDNRIS,MERFI,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   COPYRIGHT           - 1980 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C

      SUBROUTINE GGNML  (DSEED,NR,R)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NR
      REAL               R(NR)
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER             IER
C                                  FIRST EXECUTABLE STATEMENT
C                                  GET NR RANDOM NUMBERS
      CALL GGUBS(DSEED,NR,R)
C                                  TRANSFORMS EACH UNIFORM DEVIATE
      DO 5 I=1,NR
         CALL MDNRIS(R(I),R(I),IER)
    5 CONTINUE
      RETURN
      END
      
