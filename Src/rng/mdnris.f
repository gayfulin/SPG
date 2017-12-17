C   IMSL ROUTINE NAME   - MDNRIS                                        MDRS0010
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - DGC/SINGLE
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - INVERSE STANDARD NORMAL (GAUSSIAN)
C                           PROBABILITY DISTRIBUTION FUNCTION
C
C   USAGE               - CALL MDNRIS (P,Y,IER)
C
C   ARGUMENTS    P      - INPUT VALUE IN THE EXCLUSIVE RANGE (0.0,1.0)
C                Y      - OUTPUT VALUE OF THE INVERSE NORMAL (0,1)
C                           PROBABILITY DISTRIBUTION FUNCTION
C                IER    - ERROR PARAMETER (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129 INDICATES P LIES OUTSIDE THE LEGAL
C                             RANGE. PLUS OR MINUS MACHINE INFINITY IS
C                             GIVEN AS THE RESULT (SIGN IS THE SIGN OF
C                             THE FUNCTION VALUE OF THE NEAREST LEGAL
C                             ARGUMENT).
C
C   PRECISION/HARDWARE  - SINGLE/ALL
C
C   REQD. IMSL ROUTINES - MERFI,UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MDNRIS (P,Y,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      REAL               P,Y
      INTEGER            IER
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL               EPS,G0,G1,G2,G3,H0,H1,H2,A,W,WI,SN,SD
      REAL               SIGMA,SQRT2,X,XINF
      DATA               XINF/0.7237E+76/
      DATA               SQRT2/1.414214/
      DATA               EPS/0.9537E-06/
      DATA               G0/.1851159E-3/,G1/-.2028152E-2/
      DATA               G2/-.1498384/,G3/.1078639E-1/
      DATA               H0/.9952975E-1/,H1/.5211733/
      DATA               H2/-.6888301E-1/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (P .GT. 0.0 .AND. P .LT. 1.0) GO TO 5
      IER = 129
      SIGMA = SIGN(1.0,P)
      Y = SIGMA * XINF
      GO TO 9000
    5 IF(P.LE.EPS) GO TO 10
      X = 1.0 -(P + P)
      CALL MERFI (X,Y,IER)
      Y = -SQRT2 * Y
      GO TO 9005
C                                  P TOO SMALL, COMPUTE Y DIRECTLY
   10 A = P+P
      W = SQRT(-ALOG(A+(A-A*A)))
C                                  USE A RATIONAL FUNCTION IN 1./W
      WI = 1./W
      SN = ((G3*WI+G2)*WI+G1)*WI
      SD = ((WI+H2)*WI+H1)*WI+H0
      Y = W + W*(G0+SN/SD)
      Y = -Y*SQRT2
      GO TO 9005
 9000 CONTINUE
      CALL UERTST(IER,'MDNRIS')
 9005 RETURN
      END
