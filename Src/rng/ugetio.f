C
C CALLING PROGR: NONE
C
C
C
      SUBROUTINE UGETIO(IOPT,NIN,NOUT)                                  
C-----------------------------------------------------------------------
C   PURPOSE             - TO RETRIEVE CURRENT VALUES AND TO SET NEW     
C                           VALUES FOR INPUT AND OUTPUT UNIT            
C                           IDENTIFIERS.                                
C   USAGE               - CALL UGETIO(IOPT,NIN,NOUT)                    
C   ARGUMENTS    IOPT   - OPTION PARAMETER. (INPUT)                     
C                           IF IOPT=1, THE CURRENT INPUT AND OUTPUT     
C                           UNIT IDENTIFIER VALUES ARE RETURNED IN NIN  
C                           AND NOUT, RESPECTIVELY.                     
C                           IF IOPT=2, THE INTERNAL VALUE OF NIN IS     
C                           RESET FOR SUBSEQUENT USE.                   
C                           IF IOPT=3, THE INTERNAL VALUE OF NOUT IS    
C                           RESET FOR SUBSEQUENT USE.                   
C                NIN    - INPUT UNIT IDENTIFIER.                        
C                           OUTPUT IF IOPT=1, INPUT IF IOPT=2.          
C                NOUT   - OUTPUT UNIT IDENTIFIER.                       
C                           OUTPUT IF IOPT=1, INPUT IF IOPT=3.          
C                                                                       
C   PRECISION/HARDWARE  - SINGLE/ALL                                    
C                                                                       
C   REQD. IMSL ROUTINES - NONE REQUIRED                                 
C   REMARKS      EACH IMSL ROUTINE THAT PERFORMS INPUT AND/OR OUTPUT    
C                OPERATIONS CALLS UGETIO TO OBTAIN THE CURRENT UNIT     
C                IDENTIFIER VALUES. IF UGETIO IS CALLED WITH IOPT=2 OR  
C                IOPT=3, NEW UNIT IDENTIFIER VALUES ARE ESTABLISHED.    
C                SUBSEQUENT INPUT/OUTPUT IS PERFORMED ON THE NEW UNITS. 
C-----------------------------------------------------------------------
C                                  SPECIFICATIONS FOR ARGUMENTS         
      INTEGER            IOPT,NIN,NOUT                                  
C                                  SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER            NIND,NOUTD                                     
CMMMMMMMMMMMMMMMMMMMMMMMMMMM !!!      !!! MMMMMMMMMMMMMMMMMMM
      DATA               NIND/5/,NOUTD/6/                               
C                                  FIRST EXECUTABLE STATEMENT           
      IF (IOPT.EQ.3) GO TO 10                                           
      IF (IOPT.EQ.2) GO TO 5                                            
      IF (IOPT.NE.1) GO TO 9005                                         
      NIN = NIND                                                        
      NOUT = NOUTD                                                      
      GO TO 9005                                                        
    5 NIND = NIN                                                        
      GO TO 9005                                                        
   10 NOUTD = NOUT                                                      
 9005 RETURN                                                            
      END                                                               
