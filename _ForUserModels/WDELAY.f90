! Dynamic model of communication delays between VSC converter stations for wide-area controls
! 
! The communication delays can be constant or they can vary randomly according to a triangular density function.
! 
! The supplementary active power control is implemented as a minimum excitation limiter model.  

SUBROUTINE WDELAY(I_MACH,I_SLOT)

	!DEC$ ATTRIBUTES DLLEXPORT, DECORATE, ALIAS: "WDELAY":: WDELAY
	!DEC$ ATTRIBUTES REFERENCE :: I_MACH,I_SLOT

	INCLUDE 'COMON4.INS'
	IMPLICIT none

	! Declaration
	! ----------- 
	INTEGER IB, I_SLOT, I_MACH, I_BUS
	INTEGER I_ICON, I_CON, I_STATE, I_VAR, I_VAR_EXC, I_VAR_STAB
	INTEGER i, IERR
	
	INTEGER MONOFF, NVSC

	REAL SIGMA, DELTAT
	REAL deltauref0, deltauref

	INTEGER, ALLOCATABLE :: I_BUSES(:), IB_SEQ(:), I_NZTAU(:)	
	REAL, ALLOCATABLE :: TAU(:), TAU_RANDOM(:), alpha1(:), alpha2(:), beta1(:), beta2(:) 
	REAL, ALLOCATABLE :: deltaFbus(:), Thetabus(:), deltau(:), deltaudelay(:), x(:), d_x(:), x1(:), d_x1(:), x2(:), d_x2(:) 
	REAL, ALLOCATABLE :: randvector1(:), randvector2(:)

	! Parameter assignment
	! --------------------
	IB = NUMTRM(I_MACH) 
	I_BUS = NUMBUS(IB) 

	! Index
	I_CON=STRTIN(1,I_SLOT)
	I_STATE=STRTIN(2,I_SLOT)
	I_VAR=STRTIN(3,I_SLOT)
	I_ICON=STRTIN(4,I_SLOT)

	CALL MDLIND(I_BUS, MACHID(I_MACH), 'EXC', 'VAR', I_VAR_EXC, IERR)
	CALL MDLIND(I_BUS, MACHID(I_MACH), 'STAB', 'VAR', I_VAR_STAB, IERR)
	
	! ICON assignment
	MONOFF = ICON(I_ICON) ! enable/disable the delays
	NVSC = ICON(I_ICON+1) ! number of converters participating in the coordinated control
	
	! CON assignment
	SIGMA = CON(I_CON + NVSC)

	! Allocate matrices
	ALLOCATE(I_BUSES(NVSC))
	ALLOCATE(I_NZTAU(NVSC))
	ALLOCATE(IB_SEQ(NVSC))
	ALLOCATE(TAU(NVSC))
	ALLOCATE(TAU_RANDOM(NVSC))
	ALLOCATE(alpha1(NVSC))
	ALLOCATE(alpha2(NVSC))
	ALLOCATE(beta1(NVSC))
	ALLOCATE(beta2(NVSC))
	
	ALLOCATE(randvector1(NVSC))
	ALLOCATE(randvector2(NVSC))

	ALLOCATE(deltaFbus(NVSC)) 
	ALLOCATE(Thetabus(NVSC))
	ALLOCATE(deltau(NVSC)) 
	ALLOCATE(deltaudelay(NVSC)) 
	ALLOCATE(x(2*NVSC)) 
	ALLOCATE(d_x(2*NVSC)) 
	ALLOCATE(x1(NVSC)) 
	ALLOCATE(d_x1(NVSC)) 
	ALLOCATE(x2(NVSC)) 
	ALLOCATE(d_x2(NVSC)) 

	! ICONs
	I_BUSES = ICON(I_ICON+2:I_ICON+NVSC+1)
	
	! CONs
	TAU = CON(I_CON:I_CON+NVSC-1)
	
	! Other parameters
	CALL RANDOM_NUMBER(randvector1)
	CALL RANDOM_NUMBER(randvector2)
	TAU_RANDOM = MAX(TAU + SIGMA*(randvector1-randvector2),0.0)	! triangular density function
	alpha1 = TAU_RANDOM/2 ! Padé with order 2x2
	alpha2 = TAU_RANDOM**2/12
	beta1 = -alpha1
	beta2 = alpha2

	! input	
	DO i=1,NVSC
			CALL BSSEQN(I_BUSES(i),IB_SEQ(i)) ! MTDC buses seq number
			deltaFbus(i) = 0.0 + BSFREQ(IB_SEQ(i))  
			Thetabus(i) = ATAN2(AIMAG(VOLT(IB_SEQ(i))),REAL(VOLT(IB_SEQ(i)))) 			
	END DO
	
	! VARs
	deltaudelay = VAR(I_VAR:I_VAR+NVSC-1)
	deltauref0 = VAR(I_VAR + NVSC)
	deltauref = VAR(I_VAR + NVSC+1)

	! STATEs
	x = STATE(I_STATE:I_STATE+2*NVSC-1) 
	d_x = DSTATE(I_STATE:I_STATE+2*NVSC-1) 
	x1 = x(1:NVSC)
	d_x1 = d_x(1:NVSC)
	x2 = x(NVSC+1:2*NVSC)
	d_x2 = d_x(NVSC+1:2*NVSC)

	CALL DSRVAL('DELT', 1, DELTAT, IERR) 
	I_NZTAU = PACK([ (i, i = 1, NVSC) ], (TAU>(2*DELTAT)))
	
	deltau = 0.0
	IF (MONOFF.EQ.1) THEN
		deltau = deltaFbus
	ELSE IF (MONOFF.EQ.4) THEN
		deltau = Thetabus
	END IF
	
	SELECT CASE (MODE)

	CASE (1) 

		! Inicialization 
		! ==============	
		
		x1 = deltau
		x2 = 0.0
		x(1:NVSC) = x1
		x(NVSC+1:2*NVSC) = x2
		deltaudelay = deltau
		deltauref0 = SUM(deltaudelay)/NVSC ! output is a variation of an average value around an initial one (f or theta)
								
	CASE (2) 
	   
		! Compute derivatives
		! ====================		
		d_x1 = 0.0
		d_x2 = 0.0
		IF (SIZE(I_NZTAU).GT.0) THEN
			d_x1(I_NZTAU) = x2(I_NZTAU)
			d_x2(I_NZTAU) = (1/alpha2(I_NZTAU))*(-x1(I_NZTAU) - alpha1(I_NZTAU)*x2(I_NZTAU) + deltau(I_NZTAU))
		END IF
		d_x(1:NVSC) = d_x1
		d_x(NVSC+1:2*NVSC) = d_x2

	CASE (3) 

		! Compute output
		! ==============	
		deltaudelay = deltau
		IF (SIZE(I_NZTAU).GT.0) THEN
			deltaudelay(I_NZTAU) = x1(I_NZTAU) + beta1(I_NZTAU)*x2(I_NZTAU) + beta2(I_NZTAU)*d_x2(I_NZTAU)
		END IF
		deltauref = SUM(deltaudelay)/NVSC - deltauref0
		
		IF ((MONOFF.EQ.1).OR.(MONOFF.EQ.4)) THEN	! active power control
			VAR(I_VAR_EXC+1) = deltauref 
		ELSE IF (MONOFF.EQ.2) THEN					! reactive power control
			VAR(I_VAR_STAB+1) = deltauref
		ELSE IF (MONOFF.EQ.3) THEN 					! both
			VAR(I_VAR_EXC+1) = deltauref 
			VAR(I_VAR_STAB+1) = deltauref
		END IF
		  
	CASE (4) 
		! Update number of STATEs
		! =======================
		NINTEG = MAX(NINTEG,I_STATE+2*NVSC-1)

	CASE (5)
		! Reporting mode
		! ==============

	CASE DEFAULT

	END SELECT

	! RE-ASSIGN VARIABLES
	! -------------------

	! VARs
	VAR(I_VAR:I_VAR+NVSC-1) = deltaudelay
	VAR(I_VAR + NVSC) = deltauref0
	VAR(I_VAR + NVSC+1) = deltauref
	
	! STATEs
	STATE(I_STATE:I_STATE+NVSC-1) = x
	DSTATE(I_STATE:I_STATE+NVSC-1) = d_x

END