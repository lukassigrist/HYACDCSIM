MODULE MOD_TFBLOCKS

	IMPLICIT NONE
	
CONTAINS
	
	SUBROUTINE SUB_INTEGRATORWINDUP(y,x,dxdt,e,MODE,T,YMAX,YMIN)
	
		! T*x(s)*s = e(s) = u(s)
		! y(s) = min(max(x(s),XMIN),XMAX)
		! x0 = y0, u0 = 0
		
		INTEGER, INTENT(IN) :: MODE
		REAL, INTENT(IN) :: T, YMAX, YMIN, e
		REAL, INTENT(OUT) :: dxdt
		REAL, INTENT(INOUT) :: x, y
		
		dxdt = 0.0
		
		SELECT CASE (MODE)
			
			CASE (1)
				
				x = y
			
			CASE (2)
				
				IF ((y.GE.YMIN).or.(y.LE.YMAX)) THEN
					dxdt = e/T
				END IF
			
			CASE (3)
			
				y = MAX(MIN(x,YMAX),YMIN)
			
			CASE DEFAULT
		
		END SELECT
		
	END SUBROUTINE SUB_INTEGRATORWINDUP
	
	SUBROUTINE SUB_FIRSTORDERWINDUP(y,x,dxdt,u,MODE,INIDIR,DELTAT,T,YMAX,YMIN)
		
		! T*x(s)*s = -x(s) + u(s)
		! y(s) = min(max(x(s),XMIN),XMAX) 
		! x0 = y0 or x0 = K*u0
		
		INTEGER, INTENT(IN) :: MODE, INIDIR
		REAL, INTENT(IN) :: T, YMAX, YMIN, DELTAT, u
		REAL, INTENT(OUT) :: dxdt
		REAL, INTENT(INOUT) :: x, y
		
		dxdt = 0.0
		
		SELECT CASE (MODE)
			
			CASE (1)
				
				IF (INIDIR.GT.0) THEN
				
					! x = y
					CALL SUB_INTEGRATORWINDUP(y,x,dxdt,(u-x),1,T,YMAX,YMIN)
					
				ELSEIF (INIDIR.LT.0) THEN
				
					x = u
				
				END IF
			
			CASE (2)
				
				IF (T.LT.(2*DELTAT)) THEN
					x = u
				ELSE
					IF ((y.GE.YMIN).or.(y.LE.YMAX)) THEN
						! dxdt = (-x + u)/T
						CALL SUB_INTEGRATORWINDUP(y,x,dxdt,(u-x),2,T,YMAX,YMIN)
					END IF
				END IF
			
			CASE (3)
			
				! y = MAX(MIN(x,YMAX),YMIN)
				CALL SUB_INTEGRATORWINDUP(y,x,dxdt,(u-x),3,T,YMAX,YMIN)
			
			CASE DEFAULT
		
		END SELECT
		
	END SUBROUTINE SUB_FIRSTORDERWINDUP
	
	SUBROUTINE SUB_WASHOUTWINDUP(y,x,dxdt,u,MODE,DELTAT,T)
	
		! T*x(s)*s = y(s) = u(s) - x(s)
		! y(s) = -x(s) + u(s)  
		! x0 = u0, y0 = 0
		
		INTEGER, INTENT(IN) :: MODE
		REAL, INTENT(IN) :: T, DELTAT, u
		REAL :: YMAX, YMIN
		REAL, INTENT(OUT) :: dxdt
		REAL, INTENT(INOUT) :: y, x
		
		YMAX = 10e6
		YMIN = -10e6
		dxdt = 0.0
		
		SELECT CASE (MODE)
			
			CASE (1)
				
				! x = u
				CALL SUB_INTEGRATORWINDUP(y,x,dxdt,(u-x),1,T,YMAX,YMIN)
			
			CASE (2)
				
				! dxdt = y/T = (u-x)/T
				CALL SUB_INTEGRATORWINDUP(y,x,dxdt,(u-x),2,T,YMAX,YMIN)
			
			CASE (3)
			
				! y = u-x
				CALL SUB_INTEGRATORWINDUP(y,x,dxdt,(u-x),3,T,YMAX,YMIN)
				y = u - y
			
			CASE DEFAULT
		
		END SELECT
				
	END SUBROUTINE SUB_WASHOUTWINDUP
	
	SUBROUTINE SUB_LEADLAGWINDUP(y,x,dxdt,u,MODE,INIDIR,DELTAT,T1,T2)
		
		! T1*x(s)*s = -x(s) + u(s)
		! y(s) = (1-T2/T1)*x(s) + T2/T1*u(s)
		! x0 = y0 or x0 = u0
		
		INTEGER, INTENT(IN) :: MODE, INIDIR
		REAL, INTENT(IN) :: T1, T2, DELTAT, u
		REAL :: YMAX, YMIN
		REAL, INTENT(OUT) :: dxdt
		REAL, INTENT(INOUT) :: y, x
		
		YMAX = 10e6
		YMIN = -10e6
		dxdt = 0.0
		
		SELECT CASE (MODE)
			
			CASE (1)
				
				IF (INIDIR.GT.0) THEN
				
					! x = y
					CALL SUB_INTEGRATORWINDUP(y,x,dxdt,(u-x),1,T1,YMAX,YMIN)
					
				ELSEIF (INIDIR.LT.0) THEN
				
					x = u
				
				END IF
			
			CASE (2)
				
				! dxdt = (-x + u)/T1
				CALL SUB_INTEGRATORWINDUP(y,x,dxdt,(u-x),2,T1,YMAX,YMIN)
			
			CASE (3)
			
				! y = (1-T2/T1)*x + T2/T1*u
				CALL SUB_INTEGRATORWINDUP(y,x,dxdt,(u-x),3,T1,YMAX,YMIN)
				y = ((1-T2/T1)*y+T2/T1*u)
			
			CASE DEFAULT
		
		END SELECT
		
	END SUBROUTINE SUB_LEADLAGWINDUP
	
	SUBROUTINE SUB_PI(y,x,dxdt,e,MODE,DELTAT,KP,KI,YMAX,YMIN)
	
		! 1/Ki*x(s)*s = e(s)
		! y(s) = Kp*e(s) + x(s)  
		! x0 = y0, u0 = 0
		
		INTEGER, INTENT(IN) :: MODE
		REAL, INTENT(IN) :: KP, KI, YMAX, YMIN, DELTAT, e
		REAL, INTENT(OUT) :: dxdt
		REAL, INTENT(INOUT) :: y, x

		dxdt = 0.0
		
		SELECT CASE (MODE)
			
			CASE (1)
				
				! x = u
				CALL SUB_INTEGRATORWINDUP(y,x,dxdt,e,1,1/KI,YMAX,YMIN)
			
			CASE (2)
				
				! dxdt = e*KI
				IF (KI.NE.0.0) THEN
					CALL SUB_INTEGRATORWINDUP(y,x,dxdt,e,2,1/KI,YMAX,YMIN)
				END IF
				
			
			CASE (3)
			
				! y = KP*e+x
				CALL SUB_INTEGRATORWINDUP(y,x,dxdt,e,3,1/KI,YMAX,YMIN)
				y = y + KP*e
			
			CASE DEFAULT
		
		END SELECT
				
	END SUBROUTINE SUB_PI
	
END MODULE MOD_TFBLOCKS