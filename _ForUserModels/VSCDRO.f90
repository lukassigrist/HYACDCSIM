! Dynamic model of a voltage source with current limits
! 
! The voltage dynamics are represented by a simple first-order system. The control is formulated in a dq reference system aligned
! with the internal voltage. No PLL is used. The droop implementation controls the angle/frequency. 
! 
! ec,dq = ec,RI*exp(-j*delta) = ec + j*0.0
! ec,dq = [1/(1+s*Td) 0; 0 1/(1+s*Tq)]*ecref,dq 
! s*delta = Omegab*dw
! dw = 1/mp/(1+s*Tp)*(p0 - ps)
! ecref,dq = 1/mq/(1+s*Tq) 0]*(q0 - qs) + [ec0 0]
! 
! Supervisory controls such as active and reactive power controls, DC-voltage control, etc. are implemented in further
! user models. These user models are implemented as stabilizer (SQWDRD), excitation (SPWDRD), and governor type models 
! (DCGRID). Note that the calling sequence for the plant models is: 1. Generator models, 2. Current compensating models, 
! 3. Excitation stabilizer models, 4. Excitation system models, 5. Turbine-governor models. Turbine governor, stabilizer
! and excitation limiter models have no initialization duties other than STATEs and VARs.
! 
! This simple voltage source model, VOLSOU, is implemented as a coordinated-call generator model 
! (IC = 1 and IT = 01) since ZSORCE is either the transformer impedance or very small (0.001 pu) if a transformer is explicitely modelled. 
 
! ======================================================================================
! MODULE DECLARATION
! ======================================================================================
MODULE MOD_VSCDRO_INTERNAL
	! Internal variables
	INTEGER :: I_STATE_DCGRID_G = -1
END MODULE MOD_VSCDRO_INTERNAL      

! ==========================================================
! SOLUTION OF NETWORK EQUATIONS
! ==========================================================   
SUBROUTINE TSCDRO(I_MACH,I_SLOT)

	USE MOD_MISC 											! use module of miscellaneous functions
	INCLUDE 'COMON4.ins'
	IMPLICIT NONE

	! Declaration
	! ----------- 
	INTEGER :: I_MACH, I_SLOT, IB, I_STATE, I_CON, I_ICON		
	INTEGER :: ILIMITPRIORITY, ISHARDCLIMIT
	   
	REAL :: xdelta, xec, xecq
	REAL :: ectref, ICMAXIN, icmaxinpu

	COMPLEX :: ect_phasor_RI, ict_phasor_RI, us_phasor_RI, ss_phasor
	COMPLEX :: ect_phasor_dq
	COMPLEX :: zc

	COMPLEX	:: reference_transf
	
	! Parameter and variable assignment
	! ---------------------------------
	! Index
	IB = NUMTRM(I_MACH) 		! Bus sequence number
	I_CON = STRTIN(1,I_SLOT)	! initial CON index
	I_STATE = STRTIN(2,I_SLOT)	! initial STATE index
	I_ICON = STRTIN(4,I_SLOT)	! initial ICON index

	! STATES
	xec = STATE(I_STATE)         
	xdelta = STATE(I_STATE+1)

	! ICONs
	ISHARDCLIMIT = ICON(I_ICON+5) 		! Hard current limit
	ILIMITPRIORITY = ICON(I_ICON+6) 	! Current limit priority

	! CONs
	ICMAXIN = CON(I_CON+19) 						! maximum instantaneous current (pu)
	icmaxinpu = ICMAXIN*MBASE(I_MACH)/SBASE			! instantaneous current limit (120% of steady-state limit)

	! Common variables
	us_phasor_RI = VOLT(IB)	
	zc = ZSORCE(I_MACH)*SBASE/MBASE(I_MACH)

	
	SELECT CASE (MODE)
		CASE (1) 
			! Initialization
			! --------------
			! ec = ut + ZSORCE*ict, ISORCE = ec/ZSORCE -> ict = ISORCE - ut/ZSORCE
			! ISORCE: norton equivalent source current in pu of SBASE
			! ict: current at generator terminal bus in pu SBASE
			! ut: voltage at generator terminal bus in pu
			! ZSORCE: in pu of MBASE -> zc = ZSORCE*SBASE/MBASE
			ict_phasor_RI = ISORCE(I_MACH) - us_phasor_RI/zc
			ss_phasor = us_phasor_RI*CONJG(ict_phasor_RI)
			ect_phasor_RI = zc*ISORCE(I_MACH)	
			ectref = ABS(ect_phasor_RI)	
			xdelta = ATAN2(AIMAG(ect_phasor_RI),REAL(ect_phasor_RI))	
			xec = ectref		

		CASE (3) 
			! Current injections
			! ------------------         
			ect_phasor_dq = CMPLX(xec,xecq)	
	        reference_transf = CMPLX(COS(xdelta),SIN(xdelta))
			ect_phasor_RI = ect_phasor_dq*reference_transf
			
			! Hard instantaneous current limit (ict < icmaxinpu)
			IF (ISHARDCLIMIT.EQ.1) THEN
				ict_phasor_RI = (ect_phasor_RI - us_phasor_RI)/zc

				CALL SUB_ICMAXLIMITSECREF(ect_phasor_dq, (ict_phasor_RI/reference_transf), (us_phasor_RI/reference_transf), ILIMITPRIORITY, icmaxinpu, zc)	
				ect_phasor_RI = ect_phasor_dq*reference_transf
			END IF

			ict_phasor_RI = (ect_phasor_RI - us_phasor_RI)/zc
            ISORCE(I_MACH) = ect_phasor_RI/zc			 
			ss_phasor = us_phasor_RI*CONJG(ict_phasor_RI)

			!IF ((IFLAG.GT.0).AND.(TIME.GT.-1.5)) WRITE (LPDEV,*) 'TSCGFO - CASE 3: ',TIME,(ect_phasor_RI/reference_transf), (us_phasor_RI/reference_transf),(ict_phasor_RI/reference_transf),real(ect_phasor_RI*CONJG(ict_phasor_RI))
			
		CASE DEFAULT

	END SELECT

	! Re-assign algebraic variables
	! -----------------------------
	PELEC(I_MACH) = REAL(ss_phasor)
	QELEC(I_MACH) = AIMAG(ss_phasor)
	ETERM(I_MACH) = ABS(us_phasor_RI)
				
END SUBROUTINE TSCDRO

! ==========================================================
! SOLUTION OF DIFFERENTIAL EQUATIONS
! ==========================================================   
SUBROUTINE VSCDRO(I_MACH,I_SLOT)

	USE MOD_TFBLOCKS 										! use module of generic transfer functions
	USE MOD_READHYADCSIM 									! use module of reading HYADCSIM data
	USE MOD_PROTECTION 										! use module of protection functions
	USE MOD_MISC 											! use module of miscellaneous functions
	USE MOD_VSCDRO_INTERNAL									! use module for global VSCDRO-related variables
	INCLUDE 'COMON4.INS'									! common PSS/e variables and modules
	IMPLICIT NONE
	
	! Declaration
	! -----------
	INTEGER :: I_SLOT, I_MACH, I_VAR, I_CON, I_ICON, I_STATE	! PSS/e indices
	INTEGER :: IB												! Bus index
	INTEGER :: IDXCONVERTER, NDCBUS, NDCBUS_PREVIOUS, NDCLINES_PREVIOUS
	INTEGER, ALLOCATABLE :: m_VSCACDCBUS(:,:) ! ac and dc buses of each converter

    INTEGER :: ierr, I_STATE_DCGRID, ILIMITPRIORITY
		
	REAL :: xec, xdelta, xdc, xp, xq, xpmax1, xpmax2, xqmax1, xqmax2
	REAL :: d_xec, d_xdelta, d_xdc, d_xp, d_xq, d_xpmax1, d_xpmax2, d_xqmax1, d_xqmax2
	REAL :: ydc, ypmax1, yqmax1, ypmax2, yqmax2
	REAL :: deltap, deltaq, deltaec, dw
	
	REAL :: ps, qs, us, thetas
	REAL :: pct, qct, pctref, qctref, ect, ectref, deltapctref, deltaqctref
	REAL :: ictd, ictq, ict
	REAL :: pdc, udc, udcref
	REAL :: ploss
	REAL, ALLOCATABLE :: v_udc0(:,:)

	COMPLEX us_phasor_RI, ict_phasor_RI, ict_phasor_dq, ect_phasor_RI, ect_phasor_dq
    COMPLEX ss_phasor, sct_phasor
	COMPLEX	reference_transf

	REAL :: TAU, MP, TP, MQ, TQ, KPDC, KIDC, PCMAX, PCMIN, QCMAX, QCMIN, ALOSS_MW, BLOSS_kV, CLOSS_RECT_Ohm, CLOSS_INV_Ohm, ICMAXIN, KPPMAX, KPIMAX, KQPMAX, KQIMAX, ECMAX, ECMIN
	REAL :: mppu, mqpu, kpdcpu, kidcpu, pcmaxpu, pcminpu, qcmaxpu, qcminpu, aloss, bloss, c_inv, c_rect, icmaxinpu, kppmaxpu, kpimaxpu, kqpmaxpu, kqimaxpu

    REAL :: OMEGABASE, FBASE, PI, ZBASE
    PARAMETER (PI=3.14159265358979)
	REAL :: DELTAT
	COMPLEX :: zc

	CHARACTER(1) :: IDGRID
		
	! Parameter and variable assignment
	! ---------------------------------	  
	! Indexes
	IB = NUMTRM(I_MACH)
	I_CON = STRTIN(1,I_SLOT)
	I_STATE = STRTIN(2,I_SLOT) 
	I_VAR = STRTIN(3,I_SLOT)
	I_ICON = STRTIN(4,I_SLOT)
	
	! ICONs
	IDXCONVERTER = ICON(I_ICON) 		! absolute index number of all VSC
	NDCBUS = ICON(I_ICON+1) 			! Number of converters of the DC grid
	IDGRID = CHRICN(I_ICON+2)			! DC grid identifier
	NDCBUS_PREVIOUS = ICON(I_ICON+3) 	! accumulated number of DC buses
	NDCLINES_PREVIOUS = ICON(I_ICON+4) 	! accumulated numnber of DC lines
	ILIMITPRIORITY = ICON(I_ICON+6) 	! Current limit priority
			
	! CONs    
	TAU = CON(I_CON) 				! Inverter time constant (e.g., 0.01 s)
	MP = CON(I_CON+1) 			    ! Active power droop gain (e.g., 0.02 pu)
	TP = CON(I_CON+2) 			    ! Active power droop time constant (e.g., 0.1 s)
	MQ = CON(I_CON+3) 				! Reactive power droop control gain (e.g., 0.05 pu)
	TQ = CON(I_CON+4) 				! Reactive power droop time constant (e.g., 0.1 s)
	KPDC = CON(I_CON+5) 			! DC voltage control gain (e.g., 0.1 pu)
	KIDC = CON(I_CON+6) 			! DC voltage control integral gain (e.g., 0.1 pu)
    KPPMAX = CON(I_CON+7)           ! Active power limiter proportional gain (e.g., 0.01)
    KPIMAX = CON(I_CON+8)           ! Active power limiter integral gain (e.g., 0.1)
	PCMAX = CON(I_CON+9) 			! Maximum power (e.g., 1.0 pu)
	PCMIN = CON(I_CON+10) 			! Minimum power (e.g., -1.0 pu)
    KQPMAX = CON(I_CON+11)           ! Reactive power limiter proportional gain (e.g., 0.1)
    KQIMAX = CON(I_CON+12)          ! Reactive power limiter integral gain (e.g., 10.0)
	QCMAX = CON(I_CON+13) 			! Maximum reactive power (e.g., 0.1 pu)
	QCMIN = CON(I_CON+14) 			! Minimum reactive power (e.g., -0.1 pu)
	ALOSS_MW = CON(I_CON+15) 		! constant converter loss coefficient (MW): ploss = aloss + bloss*ict + c*ict^2
	BLOSS_kV = CON(I_CON+16) 		! linear converter loss coefficient (kV): ploss = aloss + bloss*ict + c*ict^2
	CLOSS_RECT_Ohm = CON(I_CON+17) 	! rectifier quadratic converter loss coefficient (ohm): ploss = aloss + bloss*ict + c*ict^2
	CLOSS_INV_Ohm = CON(I_CON+18) 	! constant converter loss coefficient (ohm): ploss = aloss + bloss*ict + c*ict^2
	ICMAXIN = CON(I_CON+19) 		! maximum instantaneous current (pu)
    ECMAX = CON(I_CON+20) 		    ! maximum voltage (pu)
    ECMIN = CON(I_CON+21) 		    ! minimum voltage (pu)
	
	! VARs
	pctref = VAR(I_VAR) 	
	qctref = VAR(I_VAR+1) 
	deltapctref = VAR(I_VAR+2)
	deltaqctref = VAR(I_VAR+3)
    ectref = VAR(I_VAR+4)
	udcref = VAR(I_VAR+5)

	ydc = VAR(I_VAR+6)
	ypmax1 = VAR(I_VAR+7)
	ypmax2 = VAR(I_VAR+8)
	yqmax1 = VAR(I_VAR+9)
	yqmax2 = VAR(I_VAR+10)

	ictd = VAR(I_VAR+11)
	ictq = VAR(I_VAR+12)
	
	! STATEs
	xec = STATE(I_STATE)         
	xdelta = STATE(I_STATE+1) 
	xdc = STATE(I_STATE+2) 
	xp = STATE(I_STATE+3) 
	xq = STATE(I_STATE+4) 
    xpmax1 = STATE(I_STATE+5)
    xpmax2 = STATE(I_STATE+6)
    xqmax1 = STATE(I_STATE+7)
    xqmax2 = STATE(I_STATE+8)

	! DSTATEs
	d_xec = DSTATE(I_STATE)         
	d_xdelta = DSTATE(I_STATE+1) 
	d_xdc = DSTATE(I_STATE+2) 
	d_xp = DSTATE(I_STATE+3) 
	d_xq = DSTATE(I_STATE+4) 
    d_xpmax1 = DSTATE(I_STATE+5)
    d_xpmax2 = DSTATE(I_STATE+6)
    d_xqmax1 = DSTATE(I_STATE+7)
    d_xqmax2 = DSTATE(I_STATE+8)

	! Base changes for parameters (in pu system rating)
    mppu = MP*SBASE/MBASE(I_MACH) 						
	mqpu = MQ*SBASE/MBASE(I_MACH)												
	kpdcpu = KPDC*MBASE(I_MACH)/SBASE			
	kidcpu = KIDC*MBASE(I_MACH)/SBASE
    kppmaxpu = KPPMAX*SBASE/MBASE(I_MACH)
    kpimaxpu = KPIMAX*SBASE/MBASE(I_MACH)
    kqpmaxpu = KQPMAX*SBASE/MBASE(I_MACH)
    kqimaxpu = KQIMAX*SBASE/MBASE(I_MACH)			
	pcmaxpu = PCMAX*MBASE(I_MACH)/SBASE			
	pcminpu = PCMIN*MBASE(I_MACH)/SBASE			
	qcmaxpu = QCMAX*MBASE(I_MACH)/SBASE			
	qcminpu = QCMIN*MBASE(I_MACH)/SBASE
	icmaxinpu = ICMAXIN*MBASE(I_MACH)/SBASE				
	ZBASE = BASVLT(IB)**2/SBASE             ! AC-side Zbase in ohms
	aloss = ALOSS_MW/SBASE                 
	bloss = BLOSS_kV/SQRT(3.0)/BASVLT(IB)
	c_rect = CLOSS_RECT_Ohm/ZBASE
	c_inv = CLOSS_INV_Ohm/ZBASE
	
	! Common variables
	ALLOCATE(v_udc0(NDCBUS,1)) ! allocate size
	ALLOCATE(m_VSCACDCBUS(NDCBUS,2))

	CALL DSRVAL('DELT', 1, DELTAT, ierr)
    CALL BASE_FREQUENCY_API(FBASE, 1, ierr)
    OMEGABASE = 2*PI*FBASE
	
	us_phasor_RI = VOLT(IB)
	us = ABS(us_phasor_RI)
	thetas = ATAN2(AIMAG(us_phasor_RI),REAL(us_phasor_RI))
	zc = ZSORCE(I_MACH)*SBASE/MBASE(I_MACH)

	I_STATE_DCGRID = I_STATE_DCGRID_G

    ! ec = ut + ZSORCE*ict, ISORCE = ec/ZSORCE -> ict = ISORCE - ut/ZSORCE
	! ISORCE: norton equivalent source current in pu of SBASE
	! ict: current at generator terminal bus in pu SBASE
	! ut: voltage at generator terminal bus in pu
	! ZSORCE: in pu of MBASE -> zc = ZSORCE*SBASE/MBASE
    ict_phasor_RI = ISORCE(I_MACH) - us_phasor_RI/zc
	ict = ABS(ict_phasor_RI)
	ss_phasor = us_phasor_RI*CONJG(ict_phasor_RI)
	ps = REAL(ss_phasor)
	qs = AIMAG(ss_phasor)
	ect_phasor_RI = zc*ISORCE(I_MACH)
	ect = ABS(ect_phasor_RI)
	sct_phasor = ect_phasor_RI*CONJG(ict_phasor_RI)
	pct = REAL(sct_phasor)
	qct = AIMAG(sct_phasor) 
    
	reference_transf = CMPLX(cos(xdelta),sin(xdelta))
	ict_phasor_dq = ict_phasor_RI/reference_transf
	ictd = REAL(ict_phasor_dq)
	ictq = AIMAG(ict_phasor_dq)

	! IF ((TIME.GT.-1.5)) WRITE (LPDEV,*) 'VSCGFO - CASE ',MODE ,': ',TIME,(ect_phasor_RI/reference_transf), (us_phasor_RI/reference_transf),(ict_phasor_dq),pct

	SELECT CASE (MODE)

		CASE (1) 
			   
			! Inicialization 
			! ==============

			! Get initial value of the DC voltage state from .txt files
			CALL SUB_READDCBUS(v_udc0, m_VSCACDCBUS, IDGRID, NDCBUS)
			
			! Get STATE index position of the dc grid model DCGRID (governor-type model)	  
			CALL MDLIND(m_VSCACDCBUS(1,2), MACHID(I_MACH), 'GOV', 'STATE', I_STATE_DCGRID, ierr)
			IF (ierr.NE.0) THEN ! No governor-type model
				WRITE (LPDEV,*) 'VSCDRO - CASE 1: No governor-type model found at bus ',NUMBUS(IB),' with id ',MACHID(I_MACH),'.'
			END IF
			IF (I_STATE_DCGRID_G.LT.0) THEN
				I_STATE_DCGRID_G = I_STATE_DCGRID
			END IF

			! Get variables of controlled AC bus (c)
			pctref = pct
			qctref = qct
			deltapctref = 0.0
			deltaqctref = 0.0
            ectref = ect	
			xdelta = ATAN2(AIMAG(ect_phasor_RI),REAL(ect_phasor_RI))
			
			reference_transf = CMPLX(cos(xdelta),sin(xdelta))
            
			ict_phasor_dq = ict_phasor_RI/reference_transf
			ictd = REAL(ict_phasor_dq)
			ictq = AIMAG(ict_phasor_dq)
			ict = ABS(ict_phasor_RI)
			IF (ict.GT.(MBASE(I_MACH)/SBASE)) THEN
				WRITE (LPDEV,*) 'VSCDRO - CASE 1: Current exceeds steady-state limit of ',1.0,' pu.'
			END IF		

			! DC-side voltage, current, power
			! Note that the VSC is behind a reactance not modelled in VSCGFO but in the power flow. bloss should incude its resistance, rc
			CALL SUB_COMPUTEPLOSS(ploss, pct, ict, aloss, bloss, c_inv, c_rect)
			pdc = -(pct + ploss)
			ydc = -pdc
			IF (I_STATE_DCGRID.GT.0) THEN
				udc = v_udc0(IDXCONVERTER,1) 	! extract initial dc voltage value if any DC grid model
			ELSE
				udc = 1.0
			END IF
			udcref = udc

            ypmax1 = 0.0
            ypmax2 = 0.0
            yqmax1 = 0.0
            yqmax2 = 0.0
					
			! State variables
			! ---------------			
			d_xec = 0.0
            d_xdelta = 0.0
			d_xdc = 0.0	
			d_xp = 0.0
			d_xq = 0.0
            d_xpmax1 = 0.0
            d_xpmax2 = 0.0
            d_xqmax1 = 0.0
            d_xqmax2 = 0.0

			! DC voltage control
			CALL SUB_PI(ydc,xdc,d_xdc,udc**2-udcref**2,1,DELTAT,kpdcpu,kidcpu,pcmaxpu,pcminpu)
			
			! droop
            deltap = MIN(MAX(ydc + deltapctref,pcminpu),pcmaxpu) - (pct + ploss)	
			CALL SUB_FIRSTORDERWINDUP(xp,xp,d_xp,deltap,1,-1,DELTAT,TP,10E9,-10E9)
			
            ! active power limiter
            CALL SUB_PI(ypmax1,xpmax1,d_xpmax1,(pcmaxpu-pct),1,DELTAT,kppmaxpu,kpimaxpu,0.0,-10E9)
            CALL SUB_PI(ypmax2,xpmax2,d_xpmax2,(pcminpu-pct),1,DELTAT,kppmaxpu,kpimaxpu,10E9,0.0)
            ypmax1 = MIN(ypmax1,0.0)
            ypmax2 = MAX(ypmax2,0.0)

            ! angle
			dw = mppu*xp + (ypmax1 + ypmax2)
            CALL SUB_INTEGRATORWINDUP(xdelta,xdelta,d_xdelta,dw,1,(1/OMEGABASE),10E9,-10E9)

			! AC voltage control
			deltaq = MIN(MAX(qctref + deltaqctref, qcminpu), qcmaxpu) - qct
			CALL SUB_FIRSTORDERWINDUP(xq,xq,d_xq,deltaq,1,-1,DELTAT,TQ,10E9,-10E9)
			
            ! reactive power limiter
            CALL SUB_PI(yqmax1,xqmax1,d_xqmax1,(qcmaxpu-qct),1,DELTAT,kqpmaxpu,kqimaxpu,0.0,-10E9)
            CALL SUB_PI(yqmax2,xqmax2,d_xqmax2,(qcminpu-qct),1,DELTAT,kqpmaxpu,kqimaxpu,10E9,0.0)
            yqmax1 = MIN(yqmax1,0.0)
            yqmax2 = MAX(yqmax2,0.0)

			! AC voltage dynamics
			deltaec = MIN(MAX(ectref + mqpu*xq + (yqmax1 + yqmax2),ECMIN),ECMAX)
			CALL SUB_FIRSTORDERWINDUP(xec,xec,d_xec,deltaec,1,-1,DELTAT,TAU,10e9,-10e9)
			
			WRITE (LPDEV,*) 'VSCDRO - CASE 1: Converter ',IDXCONVERTER,' at bus ',NUMBUS(IB),' with id ',MACHID(I_MACH),' initialized. Initial conditions of states K+5 to K+8 might be suspect.'  
			WRITE (LPDEV,*) 'VSCDRO - CASE 1: Converter ',IDXCONVERTER,' at bus ',NUMBUS(IB),' with id ',MACHID(I_MACH),': pct = ',pct,', qct = ',qct,', udc = ',udcref
			WRITE (LPDEV,*) 'VSCDRO - CASE 1: Converter ',IDXCONVERTER,' at bus ',NUMBUS(IB),' with id ',MACHID(I_MACH),': xec = ',xec,', xdelta = ',xdelta,', xdc = ',xdc
			WRITE (LPDEV,*) 'VSCDRO - CASE 1: Converter ',IDXCONVERTER,' at bus ',NUMBUS(IB),' with id ',MACHID(I_MACH),': deltap = ',deltap,', dw = ',dw,', deltaq = ',deltaq,', deltaec = ',deltaec
			EFD(I_MACH) = ectref
            SPEED(I_MACH) = dw
            ANGLE(I_MACH) = xdelta*180/PI
			PMECH(I_MACH) = pct + ploss
											
		CASE (2) 
		   
			! Compute derivatives
			! ===================	
			IF (I_STATE_DCGRID.GT.0) THEN
				udc = STATE(I_STATE_DCGRID + IDXCONVERTER + NDCLINES_PREVIOUS + NDCBUS_PREVIOUS - 1)
			ELSE
				udc = udcref
			END IF

			CALL SUB_COMPUTEPLOSS(ploss, pct, ict, aloss, bloss, c_inv, c_rect)

            deltap = MIN(MAX(ydc + deltapctref,pcminpu),pcmaxpu) - (pct + ploss)
            dw = mppu*xp + (ypmax1 + ypmax2)
            deltaq = MIN(MAX(qctref + deltaqctref, qcminpu), qcmaxpu) - qct
            deltaec = MIN(MAX(ectref + mqpu*xq + (yqmax1 + yqmax2),ECMIN),ECMAX)

			CALL SUB_PI(ydc,xdc,d_xdc,udc**2-udcref**2,2,DELTAT,kpdcpu,kidcpu,pcmaxpu,pcminpu)

			CALL SUB_FIRSTORDERWINDUP(xp,xp,d_xp,deltap,2,-1,DELTAT,TP,10E9,-10E9)
            CALL SUB_PI(ypmax1,xpmax1,d_xpmax1,(pcmaxpu-pct),2,DELTAT,kppmaxpu,kpimaxpu,0.0,-10E9)
            CALL SUB_PI(ypmax2,xpmax2,d_xpmax2,(pcminpu-pct),2,DELTAT,kppmaxpu,kpimaxpu,10E9,0.0)
            CALL SUB_INTEGRATORWINDUP(xdelta,xdelta,d_xdelta,dw,2,(1/OMEGABASE),10E9,-10E9)

			CALL SUB_FIRSTORDERWINDUP(xq,xq,d_xq,deltaq,2,-1,DELTAT,TQ,10E9,-10E9)
            CALL SUB_PI(yqmax1,xqmax1,d_xqmax1,(qcmaxpu-qct),2,DELTAT,kqpmaxpu,kqimaxpu,0.0,-10E9)
            CALL SUB_PI(yqmax2,xqmax2,d_xqmax2,(qcminpu-qct),2,DELTAT,kqpmaxpu,kqimaxpu,10E9,0.0)
			CALL SUB_FIRSTORDERWINDUP(xec,xec,d_xec,deltaec,2,-1,DELTAT,TAU,10e9,-10e9)

		CASE (3) 
		
			! Compute output
			! ==============
			IF (I_STATE_DCGRID.GT.0) THEN
				udc = STATE(I_STATE_DCGRID + IDXCONVERTER + NDCLINES_PREVIOUS + NDCBUS_PREVIOUS - 1)
			ELSE
				udc = udcref
			END IF

			CALL SUB_COMPUTEPLOSS(ploss, pct, ict, aloss, bloss, c_inv, c_rect)
			pdc = -(pct + ploss)
			
            deltap = MIN(MAX(ydc + deltapctref,pcminpu),pcmaxpu) - (pct + ploss)
            dw = mppu*xp + (ypmax1 + ypmax2)
            deltaq = MIN(MAX(qctref + deltaqctref, qcminpu), qcmaxpu) - qct
            deltaec = MIN(MAX(ectref + mqpu*xq + (yqmax1 + yqmax2),ECMIN),ECMAX)
			
			CALL SUB_FIRSTORDERWINDUP(xp,xp,d_xp,deltap,3,-1,DELTAT,TP,10E9,-10E9)
            CALL SUB_PI(ypmax1,xpmax1,d_xpmax1,(pcmaxpu-pct),3,DELTAT,kppmaxpu,kpimaxpu,0.0,-10E9)
            CALL SUB_PI(ypmax2,xpmax2,d_xpmax2,(pcminpu-pct),3,DELTAT,kppmaxpu,kpimaxpu,10E9,0.0)
            ypmax1 = MIN(ypmax1,0.0)
            ypmax2 = MAX(ypmax2,0.0)
            CALL SUB_INTEGRATORWINDUP(xdelta,xdelta,d_xdelta,dw,3,(1/OMEGABASE),10E9,-10E9)

			CALL SUB_FIRSTORDERWINDUP(xq,xq,d_xq,deltaq,3,-1,DELTAT,TQ,10E9,-10E9)
            CALL SUB_PI(yqmax1,xqmax1,d_xqmax1,(qcmaxpu-qct),3,DELTAT,kqpmaxpu,kqimaxpu,0.0,-10E9)
            CALL SUB_PI(yqmax2,xqmax2,d_xqmax2,(qcminpu-qct),3,DELTAT,kqpmaxpu,kqimaxpu,10E9,0.0)
            yqmax1 = MIN(yqmax1,0.0)
            yqmax2 = MAX(yqmax2,0.0)
			CALL SUB_FIRSTORDERWINDUP(xec,xec,d_xec,deltaec,3,-1,DELTAT,TAU,10e9,-10e9)

			! ect_phasor_dq = CMPLX(xec,0.0)

            ! RI -> dq reference system transformation: xRI = reference_transf*xdq		
	        ! reference_transf = CMPLX(COS(xdelta),SIN(xdelta))
			
			! Hard current limit approximation (should be actually placed in the network solution model due to the voltage dependency)	
			! ict_phasor_dq = (ect_phasor_dq - us_phasor_RI/reference_transf)/zc
			! CALL SUB_ICMAXLIMITSECREF(ect_phasor_dq, ict_phasor_dq, (us_phasor_RI/reference_transf), ILIMITPRIORITY, icmaxsspu, zc)		
			! ect_phasor_RI = ect_phasor_dq*reference_transf
            ! ISORCE(I_MACH) = ect_phasor_RI/zc       

			
            EFD(I_MACH) = xec
            SPEED(I_MACH) = dw
            ANGLE(I_MACH) = xdelta*180/PI
			PMECH(I_MACH) = pct + ploss
			
		CASE (4) 

			! Update number of STATEs.
			! ========================
			NINTEG = MAX(NINTEG,I_STATE+12)

		CASE (5)
			! reporting mode
			WRITE (LPDEV,*) 'VSCDRO model at bus ', NUMBUS(IB)
				
		CASE DEFAULT

	END SELECT

	! RE-ASSIGN VARIABLES
	! -------------------	  
	! VARs	  
	VAR(I_VAR) = pctref  	
	VAR(I_VAR+1) = qctref
	VAR(I_VAR+2) = deltapctref
	VAR(I_VAR+3) = deltaqctref
    VAR(I_VAR+4) = ectref	
	VAR(I_VAR+5) = udcref

	VAR(I_VAR+6) = ydc
	VAR(I_VAR+7) = ypmax1
	VAR(I_VAR+8) = ypmax2
	VAR(I_VAR+9) = yqmax1
	VAR(I_VAR+10) = yqmax2

	VAR(I_VAR+11) = ictd
	VAR(I_VAR+12) = ictq
	VAR(I_VAR+13) = deltap
	VAR(I_VAR+14) = (pct+ploss)

	VAR(I_VAR+25) = pdc					! To DCGRID
	VAR(I_VAR+22) = udc					! To SPWDRD
		
	! STATEs
	STATE(I_STATE) = xec        
	STATE(I_STATE+1) = xdelta
	STATE(I_STATE+2) = xdc 
	STATE(I_STATE+3) = xp 
	STATE(I_STATE+4) = xq 
    STATE(I_STATE+5) = xpmax1
    STATE(I_STATE+6) = xpmax2
    STATE(I_STATE+7) = xqmax1
    STATE(I_STATE+8) = xqmax2

	! DSTATEs
	DSTATE(I_STATE) = d_xec        
	DSTATE(I_STATE+1) = d_xdelta 
	DSTATE(I_STATE+2) = d_xdc 
	DSTATE(I_STATE+3) = d_xp 
	DSTATE(I_STATE+4) = d_xq 
    DSTATE(I_STATE+5) = d_xpmax1
    DSTATE(I_STATE+6) = d_xpmax2
    DSTATE(I_STATE+7) = d_xqmax1
    DSTATE(I_STATE+8) = d_xqmax2
      
END SUBROUTINE VSCDRO


