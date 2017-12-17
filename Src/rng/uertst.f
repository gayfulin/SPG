C
C CALLING PROGR: UGETIO
C
C
      SUBROUTINE UERTST (IER,NAME)                                      
C-----------------------------------------------------------------------
C   PURPOSE             - PRINT A MESSAGE REFLECTING AN ERROR CONDITION 
C                                                                       
C   USAGE               - CALL UERTST (IER,NAME)                        
C   ARGUMENTS    IER    - ERROR PARAMETER. (INPUT)                      
C                           IER = I+J WHERE                             
C                             I = 128 IMPLIES TERMINAL ERROR,           
C                             I =  64 IMPLIES WARNING WITH FIX, AND     
C                             I =  32 IMPLIES WARNING.                  
C                             J = ERROR CODE RELEVANT TO CALLING        
C                                 ROUTINE.                              
C                NAME   - A SIX CHARACTER LITERAL STRING GIVING THE     
C                           NAME OF THE CALLING ROUTINE. (INPUT)        
C   PRECISION/HARDWARE  - SINGLE/ALL                                    
C   REQD. IMSL ROUTINES - UGETIO                                        
C   REMARKS      THE ERROR MESSAGE PRODUCED BY UERTST IS WRITTEN        
C                ONTO THE STANDARD OUTPUT UNIT. THE OUTPUT UNIT         
C                NUMBER CAN BE DETERMINED BY CALLING UGETIO AS          
C                FOLLOWS..   CALL UGETIO(1,NIN,NOUT).                   
C                THE OUTPUT UNIT NUMBER CAN BE CHANGED BY CALLING       
C                UGETIO AS FOLLOWS..                                    
C                                NIN = 0                                
C                                NOUT = NEW OUTPUT UNIT NUMBER          
C                                CALL UGETIO(3,NIN,NOUT)                
C                SEE THE UGETIO DOCUMENT FOR MORE DETAILS.              
C-----------------------------------------------------------------------
C                                  SPECIFICATIONS FOR ARGUMENTS         
      CHARACTER*6 NAME,NAMSET,NAMEQ
      CHARACTER*1 IEQ
      INTEGER  IER
C                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      DATA               NAMSET/'UERSET'/                         
      DATA               NAMEQ/'      '/                          
      DATA               LEVEL/4/,IEQDF/0/,IEQ/'='/                     
C                                  FIRST EXECUTABLE STATEMENT           
      IF (IER.GT.999) GO TO 25                                          
      IF (IER.LT.-32) GO TO 55                                          
      IF (IER.LE.128) GO TO 5                                           
      IF (LEVEL.LT.1) GO TO 30                                          
C                                  PRINT TERMINAL MESSAGE               
      CALL UGETIO(1,NIN,IOUNIT)                                         
      IF (IEQDF.EQ.1) WRITE(IOUNIT,35) IER,NAMEQ,IEQ,NAME               
      IF (IEQDF.EQ.0) WRITE(IOUNIT,35) IER,NAME                         
      GO TO 30                                                          
    5 IF (IER.LE.64) GO TO 10                                           
      IF (LEVEL.LT.2) GO TO 30                                          
C                                  PRINT WARNING WITH FIX MESSAGE       
      CALL UGETIO(1,NIN,IOUNIT)                                         
      IF (IEQDF.EQ.1) WRITE(IOUNIT,40) IER,NAMEQ,IEQ,NAME               
      IF (IEQDF.EQ.0) WRITE(IOUNIT,40) IER,NAME                         
      GO TO 30                                                          
   10 IF (IER.LE.32) GO TO 15                                           
C                                  PRINT WARNING MESSAGE                
      IF (LEVEL.LT.3) GO TO 30                                          
      CALL UGETIO(1,NIN,IOUNIT)                                         
      IF (IEQDF.EQ.1) WRITE(IOUNIT,45) IER,NAMEQ,IEQ,NAME               
      IF (IEQDF.EQ.0) WRITE(IOUNIT,45) IER,NAME                         
      GO TO 30                                                          
   15 CONTINUE                                                          
C                                  CHECK FOR UERSET CALL                
C
         IF (NAME.NE.NAMSET) GO TO 25                             
C   20 CONTINUE          
      LEVOLD = LEVEL                                                    
      LEVEL = IER                                                       
      IER = LEVOLD                                                      
      IF (LEVEL.LT.0) LEVEL = 4                                         
      IF (LEVEL.GT.4) LEVEL = 4                                         
      GO TO 30                                                          
   25 CONTINUE                                                          
      IF (LEVEL.LT.4) GO TO 30                                          
C                                  PRINT NON-DEFINED MESSAGE            
      CALL UGETIO(1,NIN,IOUNIT)                                         
CMMMMMMMMMMMMMMMMMMMMMMM
CC    IF (IEQDF.EQ.1) WRITE(IOUNIT,50) IER,NAMEQ,IEQ,NAME               
CC    IF (IEQDF.EQ.0) WRITE(IOUNIT,50) IER,NAME                         
CMMMMMMMMMMMMMMMMMMMMMMM
   30 IEQDF = 0                                                         
      RETURN                                                            
   35 FORMAT(' *** TERMINAL ERROR',10X,'(IER = ',I3,                   
     1       ') FROM IMSL ROUTINE ',A6,A1,A6)                        
   40 FORMAT(' *** WARNING WITH FIX ERROR  (IER = ',I3,                
     1       ') FROM IMSL ROUTINE ',A6,A1,A6)                        
   45 FORMAT(' *** WARNING ERROR',11X,'(IER = ',I3,                    
     1       ') FROM IMSL ROUTINE ',A6,A1,A6)                        
   50 FORMAT(' *** UNDEFINED ERROR',9X,'(IER = ',I5,                   
     1       ') FROM IMSL ROUTINE ',A6,A1,A6)                        
C                                  SAVE P FOR P = R CASE                
C                                    P IS THE PAGE NAME                 
C                                    R IS THE ROUTINE NAME              
   55 IEQDF = 1                                                         
C     DO 60 I=1,3                                                       
   60 NAMEQ = NAME                                             
   65 RETURN                                                            
      END
