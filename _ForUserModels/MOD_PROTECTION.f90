MODULE MOD_PROTECTION

	IMPLICIT NONE

CONTAINS
	
      SUBROUTINE SUB_OVDCPROT(istripvsc, counter, tinitial, udc, timenow, UDCMAX, TUDCMAX, busnumber, MODEL, LPDEV)

	      IMPLICIT NONE

	      INTEGER, INTENT(INOUT) :: istripvsc, counter
            REAL, INTENT(INOUT) :: tinitial

            INTEGER, INTENT(IN) :: LPDEV, busnumber
	      REAL, INTENT(IN) :: udc, timenow
	      REAL, INTENT(IN) :: UDCMAX, TUDCMAX
            CHARACTER(LEN=6), INTENT(IN) :: MODEL
	
	      IF (udc.GT.UDCMAX) THEN 

		      IF (counter .EQ. 0.0) THEN
			      tinitial = timenow   
			      counter = 1
                        WRITE (LPDEV,*) 'OVDCPROT of ',MODEL,' at bus ',busnumber,' armed at time = ', tinitial
		      END IF  

		      IF ((timenow - tinitial).GT.TUDCMAX) THEN
                        IF (istripvsc.EQ.0) THEN     
                              WRITE (LPDEV,*) MODEL,' at bus ',busnumber,' trips due to OVDCPROT at time = ', timenow
                        END IF     
			      istripvsc = 1
		      END IF 
	
	      ELSE
		      tinitial = timenow
		      counter = 0
	      END IF

      END SUBROUTINE SUB_OVDCPROT

      SUBROUTINE SUB_UVDCPROT(istripvsc, counter, tinitial, udc, timenow, UDCMIN, TUDCMIN, busnumber, MODEL, LPDEV)

	      IMPLICIT NONE

	      INTEGER, INTENT(INOUT) :: istripvsc, counter
            REAL, INTENT(INOUT) :: tinitial

            INTEGER, INTENT(IN) :: LPDEV, busnumber
	      REAL, INTENT(IN) :: udc, timenow
	      REAL, INTENT(IN) :: UDCMIN, TUDCMIN
            CHARACTER(LEN=6), INTENT(IN) :: MODEL
	
	      IF (udc.LT.UDCMIN) THEN 

		      IF (counter .EQ. 0.0) THEN
			      tinitial = timenow   
			      counter = 1
                        WRITE (LPDEV,*) 'UVDCPROT of ',MODEL,' at bus ',busnumber,' armed at time = ', tinitial
		      END IF  

		      IF ((timenow - tinitial).GT.TUDCMIN) THEN
                        IF (istripvsc.EQ.0) THEN     
                              WRITE (LPDEV,*) MODEL,' at bus ',busnumber,' trips due to UVDCPROT at time = ', timenow
                        END IF     
			      istripvsc = 1
		      END IF 
	
	      ELSE
		      tinitial = timenow
		      counter = 0
	      END IF

      END SUBROUTINE SUB_UVDCPROT
	
	! Two-stage undervoltage and -frequency protection     
      SUBROUTINE SUB_UACPROT(istripvsc,counter,tinitial,us,timenow,VUV2,VUV1,TUV2,TUV1,busnumber,MODEL,TYPE,LPDEV)
     
            REAL, INTENT(INOUT) :: tinitial
            INTEGER, INTENT(INOUT) :: istripvsc, counter
            
            INTEGER, INTENT(IN) :: busnumber, LPDEV
            REAL, INTENT(IN) :: VUV2, VUV1, TUV2, TUV1
            REAL, INTENT(IN) :: us, timenow    
            CHARACTER(LEN=6), INTENT(IN) :: MODEL    
            CHARACTER(LEN=1), INTENT(IN) :: TYPE   
                  
            IF (us.LT.VUV2) THEN  
      
                  IF (counter .EQ. 0.0) THEN
                        tinitial=timenow
                        counter=1
                        WRITE (LPDEV,*) TYPE,'-UACPROT of ',MODEL,' at bus ',busnumber,' armed at time = ', tinitial
                  END IF 
      

                  IF (((timenow-tinitial) .GT. TUV1) .AND. ((timenow-tinitial) .LE. TUV2)) THEN 
                        IF (us.LT.VUV1) THEN
                              IF (istripvsc.EQ.0) THEN
                                    WRITE (LPDEV,*) MODEL,' at bus ', busnumber, ' trips due to ',TYPE,'-UACPROT at time',timenow
                              END IF
                              istripvsc=1
                        END IF 
                  END IF ! t-tini1>TUV1 & t-tini1<=TUV2
      
                  IF ((timenow-tinitial) .GT. TUV2) THEN 
                        IF (istripvsc.EQ.0) THEN
                              WRITE (LPDEV,*) MODEL,' at bus ', busnumber, ' trips due to ',TYPE,'-UACPROT at time',timenow
                        END IF
                        istripvsc=1
                  END IF ! t-tini1>TUV2
      
            END IF ! v<VU2 	              
      
            RETURN

      END SUBROUTINE SUB_UACPROT
      
      ! Two-stage overvoltage and -frequency protection
      SUBROUTINE SUB_OACPROT(istripvsc,counter,tinitial,us,timenow,VOV2,VOV1,TOV2,TOV1,busnumber,MODEL,TYPE,LPDEV)
     
            REAL, INTENT(INOUT) :: tinitial
            INTEGER, INTENT(INOUT) :: istripvsc, counter
            
            INTEGER, INTENT(IN) :: busnumber, LPDEV
            REAL, INTENT(IN) :: VOV2, VOV1, TOV2, TOV1
            REAL, INTENT(IN) :: us, timenow    
            CHARACTER(LEN=6), INTENT(IN) :: MODEL  
            CHARACTER(LEN=1), INTENT(IN) :: TYPE     
                  
            IF (us.GT.VOV2) THEN  
      
                  IF (counter .EQ. 0.0) THEN
                        tinitial=timenow
                        counter=1
                        WRITE (LPDEV,*) TYPE,'-OACPROT of ',MODEL,' at bus ',busnumber,' armed at time = ', tinitial
                  END IF 
      

                  IF (((timenow-tinitial) .GT. TOV1) .AND. ((timenow-tinitial) .LE. TOV2)) THEN 
                        IF (us.GT.VOV1) THEN
                              IF (istripvsc.EQ.0) THEN
                                    WRITE (LPDEV,*) MODEL,' at bus ', busnumber, ' trips due to ',TYPE,'-OACPROT at time',timenow
                              END IF
                              istripvsc=1
                        END IF 
                  END IF ! t-tini1>TOV1 & t-tini1<=TOV2
      
                  IF ((timenow-tinitial) .GT. TOV2) THEN 
                        IF (istripvsc.EQ.0) THEN
                              WRITE (LPDEV,*) MODEL,' at bus ', busnumber, ' trips due to ',TYPE,'-OACPROT at time',timenow
                        END IF
                        istripvsc=1
                  END IF ! t-tini1>TOV2
      
            END IF ! v<VO2 	              
      
            RETURN

      END SUBROUTINE SUB_OACPROT

END MODULE MOD_PROTECTION