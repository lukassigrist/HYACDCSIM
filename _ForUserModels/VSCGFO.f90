! Dynamic model of a voltage source with current limits
! 
! The voltage dynamics are represented by a simple first-order system. The control is formulated in a dq reference system aligned
! with the internal voltage. No PLL is used. The VSM implementation controls the angle/frequency. 
! 
! ec,dq = ec,RI*exp(-j*delta) = ec + j*0.0
! ec,dq = [1/(1+s*Td) 0; 0 1/(1+s*Tq)]*ecref,dq 
! s*delta = Omegab*(omega - 1.0)
! s*omega = 1/2/H*(p0 - ps)
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
MODULE MOD_VSCGFO_INTERNAL
	! Internal variables
	INTEGER :: I_STATE_DCGRID_G = -1
END MODULE MOD_VSCGFO_INTERNAL      

! ==========================================================
! SOLUTION OF NETWORK EQUATIONS
! ==========================================================   
SUBROUTINE TSCGFO(I_MACH,I_SLOT)

	USE MOD_MISC 											! use module of miscellaneous functions
	INCLUDE 'COMON4.ins'
	IMPLICIT NONE

	! Declaration
	! ----------- 
	INTEGER :: I_MACH, I_SLOT, IB, I_STATE, I_CON, I_ICON		
	INTEGER :: ILIMITPRIORITY, ISHARDCLIMIT
	   
	REAL :: xdelta, xecd, xecq
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
	xecd = STATE(I_STATE)         
	xecq = STATE(I_STATE+1)  
	xdelta = STATE(I_STATE+2)

	! ICONs
	ISHARDCLIMIT = ICON(I_ICON+5) 		! Hard current limit
	ILIMITPRIORITY = ICON(I_ICON+6) 	! Current limit priority

	! CONs
	ICMAXIN = CON(I_CON+26) 						! maximum instantaneous current (pu)
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
			xecd = ectref
			xecq = 0.0			

		CASE (3) 
			! Current injections
			! ------------------         
			ect_phasor_dq = CMPLX(xecd,xecq)	
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
				
END SUBROUTINE TSCGFO

! ==========================================================
! SOLUTION OF DIFFERENTIAL EQUATIONS
! ==========================================================   
SUBROUTINE VSCGFO(I_MACH,I_SLOT)

	USE MOD_TFBLOCKS 										! use module of generic transfer functions
	USE MOD_READHYADCSIM 									! use module of reading HYADCSIM data
	USE MOD_PROTECTION 										! use module of protection functions
	USE MOD_MISC 											! use module of miscellaneous functions
	USE MOD_VSCGFO_INTERNAL									! use module for global VSCGFO-related variables
	INCLUDE 'COMON4.INS'									! common PSS/e variables and modules
	IMPLICIT NONE
	
	! Declaration
	! -----------
	INTEGER :: I_SLOT, I_MACH, I_VAR, I_CON, I_ICON, I_STATE	! PSS/e indices
	INTEGER :: IB												! Bus index
	INTEGER :: IDXCONVERTER, NDCBUS, NDCBUS_PREVIOUS, NDCLINES_PREVIOUS
	INTEGER, ALLOCATABLE :: m_VSCACDCBUS(:,:) ! ac and dc buses of each converter

    INTEGER :: ierr, I_STATE_DCGRID, ILIMITPRIORITY
		
	REAL :: xecd, xecq, xdelta, xomega, xdc, xpfc, xpod, xq, xitvid, xitviq, xudcpss, xdeltamax, xdeltamin
	REAL :: d_xecd, d_xecq, d_xdelta, d_xomega, d_xdc, d_xpod, d_xpfc, d_xq, d_xitvid, d_xitviq, d_xudcpss,  d_xdeltamax, d_xdeltamin
	REAL :: ydc, ypod, yitvid, yitviq, yudcpss, ydeltamax, ydeltamin
	REAL :: deltap, deltaq, deltaecd, deltaecq
	
	REAL :: ps, qs, us, thetas
	REAL :: pct, qct, pctref, qctref, ect, ectref, deltapctref, deltaqctref
	REAL :: ictd, ictq, ict, yecdviimax, yecqviimax
	REAL :: pdc, udc, udcref
	REAL :: ploss
	REAL, ALLOCATABLE :: v_udc0(:,:)

	COMPLEX us_phasor_RI, ict_phasor_RI, ict_phasor_dq, ect_phasor_RI, ect_phasor_dq
    COMPLEX ss_phasor, sct_phasor
	COMPLEX	reference_transf

	REAL :: TAU, H, D, KPFC, TPFC, KPOD, TPOD, KQ, TQ, RTVR, TTVR, KDCPSS, TUDCPSS, KD_P2, KD_I2, KPVI, SIGMAXR, PCMAX, PCMIN, QCMAX, QCMIN, ALOSS_MW, BLOSS_kV, CLOSS_RECT_Ohm, CLOSS_INV_Ohm, ICMAXSS, ICMAXIN, KD_P3, KD_I3, DELTAWMAX, KICPRIOR
	REAL :: hpu, dpu, kpfcpu, kpodpu, kqpu, rtvrpu, kdp2pu, kdi2pu, kpvipu, pcmaxpu, pcminpu, qcmaxpu, qcminpu, aloss, bloss, c_inv, c_rect, icmaxsspu, icmaxinpu, kdp3pu, kdi3pu

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
	H = CON(I_CON+1) 				! Virtual inertia (e.g., 4 s)
    D = CON(I_CON+2) 				! Virtual damping (e.g., 2 pu)
	KPFC = CON(I_CON+3) 			! Virtual power frequency control gain (e.g., 25 pu)
	TPFC = CON(I_CON+4) 			! Virtual power frequency control time constant (e.g., 0.1 s)
	KPOD = CON(I_CON+5) 			! Virtual power oscillation damping gain (e.g., 10 pu)
	TPOD = CON(I_CON+6) 			! Virtual power oscillation damping wash-out time constant (e.g., 5 s)
	KQ = CON(I_CON+7) 				! Virtual reactive power control gain (e.g., 20 pu)
	TQ = CON(I_CON+8) 				! Virtual reactive power control time constant (e.g., 0.1 s)
	RTVR = CON(I_CON+9) 			! Transient virtual resistance (e.g., 0.1 pu)
	TTVR = CON(I_CON+10) 			! Transient virtual resistance wash-out time constant (e.g., 0.02 s)
	KDCPSS = CON(I_CON+11) 			! DC voltage damping control gain (e.g., 0.1 pu)
	TUDCPSS = CON(I_CON+12) 		! DC voltage damping control wash-out time constant (e.g., 0.1 s)
	KD_P2 = CON(I_CON+13) 			! DC voltage control gain (e.g., 0.1 pu)
	KD_I2 = CON(I_CON+14) 			! DC voltage control integral gain (e.g., 0.1 pu)
	KPVI = CON(I_CON+15) 			! Virtual impedance current limiter gain (e.g., 0.1 pu)
	SIGMAXR = CON(I_CON+16) 		! Virtual impedance current limiter X/R relation (e.g., 3 pu)
	PCMAX = CON(I_CON+17) 			! Maximum power (e.g., 1.0 pu)
	PCMIN = CON(I_CON+18) 			! Minimum power (e.g., -1.0 pu)
	QCMAX = CON(I_CON+19) 			! Maximum reactive power (e.g., 0.1 pu)
	QCMIN = CON(I_CON+20) 			! Minimum reactive power (e.g., -0.1 pu)
	ALOSS_MW = CON(I_CON+21) 		! constant converter loss coefficient (MW): ploss = aloss + bloss*ict + c*ict^2
	BLOSS_kV = CON(I_CON+22) 		! linear converter loss coefficient (kV): ploss = aloss + bloss*ict + c*ict^2
	CLOSS_RECT_Ohm = CON(I_CON+23) 	! rectifier quadratic converter loss coefficient (ohm): ploss = aloss + bloss*ict + c*ict^2
	CLOSS_INV_Ohm = CON(I_CON+24) 	! constant converter loss coefficient (ohm): ploss = aloss + bloss*ict + c*ict^2
	ICMAXSS = CON(I_CON+25) 		! maximum steady-state current (pu)
	ICMAXIN = CON(I_CON+26) 		! maximum instantaneous current (pu)
	KICPRIOR = CON(I_CON+27) 		! fraction maximum stead-state current (e.g., 0.8)
	KD_P3 = CON(I_CON+28) 			! Max. power angle limiter control gain (e.g., 0.0 pu)
	KD_I3 = CON(I_CON+29) 			! Max. power angle limiter integral gain (e.g., 0.1 pu)
	DELTAWMAX = CON(I_CON+30) 		! Max. speed variation limit (pu) 
	
	! VARs
	pctref = VAR(I_VAR) 	
	qctref = VAR(I_VAR+1) 
	deltapctref = VAR(I_VAR+2)
	deltaqctref = VAR(I_VAR+3)
    ectref = VAR(I_VAR+4)
	udcref = VAR(I_VAR+5)

	ydc = VAR(I_VAR+6)
	ypod = VAR(I_VAR+7)
	yitvid = VAR(I_VAR+8)
	yitviq = VAR(I_VAR+9)
	yudcpss = VAR(I_VAR+10)
	ydeltamax = VAR(I_VAR+11)
	ydeltamin = VAR(I_VAR+12)

	ictd = VAR(I_VAR+13)
	ictq = VAR(I_VAR+14)
	
	! STATEs
	xecd = STATE(I_STATE)         
	xecq = STATE(I_STATE+1)       
	xdelta = STATE(I_STATE+2) 
	xomega = STATE(I_STATE+3) 
	xdc = STATE(I_STATE+4) 
	xpod = STATE(I_STATE+5) 
	xpfc = STATE(I_STATE+6) 
	xq = STATE(I_STATE+7) 
	xitvid = STATE(I_STATE+8) 
	xitviq = STATE(I_STATE+9) 
	xudcpss = STATE(I_STATE+10) 
	xdeltamax = STATE(I_STATE+11)
	xdeltamin = STATE(I_STATE+12)

	! DSTATEs
	d_xecd = DSTATE(I_STATE) 
	d_xecq = DSTATE(I_STATE+1)     
	d_xdelta = DSTATE(I_STATE+2)    
	d_xomega = DSTATE(I_STATE+3)   
	d_xdc = DSTATE(I_STATE+4) 
	d_xpod = DSTATE(I_STATE+5)
	d_xpfc = DSTATE(I_STATE+6) 
	d_xq = DSTATE(I_STATE+7) 
	d_xitvid = DSTATE(I_STATE+8) 
	d_xitviq = DSTATE(I_STATE+9) 
	d_xudcpss = DSTATE(I_STATE+10)   
	d_xdeltamax = DSTATE(I_STATE+11)
	d_xdeltamin = DSTATE(I_STATE+12)

	! Base changes for parameters (in pu system rating)
	hpu = H*MBASE(I_MACH)/SBASE 	
    dpu = D*MBASE(I_MACH)/SBASE 	
	kpfcpu = KPFC*MBASE(I_MACH)/SBASE 						
	kpodpu = KPOD*MBASE(I_MACH)/SBASE 						
	kqpu = KQ*MBASE(I_MACH)/SBASE								
	rtvrpu = RTVR*SBASE/MBASE(I_MACH) 						
	kdp2pu = KD_P2*MBASE(I_MACH)/SBASE			
	kdi2pu = KD_I2*MBASE(I_MACH)/SBASE			
	pcmaxpu = PCMAX*MBASE(I_MACH)/SBASE			
	pcminpu = PCMIN*MBASE(I_MACH)/SBASE			
	qcmaxpu = QCMAX*MBASE(I_MACH)/SBASE			
	qcminpu = QCMIN*MBASE(I_MACH)/SBASE
	icmaxsspu = ICMAXSS*MBASE(I_MACH)/SBASE
	icmaxinpu = ICMAXIN*MBASE(I_MACH)/SBASE				
	kpvipu = KPVI*SBASE/MBASE(I_MACH)
	ZBASE = BASVLT(IB)**2/SBASE             ! AC-side Zbase in ohms
	aloss = ALOSS_MW/SBASE                 
	bloss = BLOSS_kV/SQRT(3.0)/BASVLT(IB)
	c_rect = CLOSS_RECT_Ohm/ZBASE
	c_inv = CLOSS_INV_Ohm/ZBASE
	kdp3pu = KD_P3*MBASE(I_MACH)/SBASE			
	kdi3pu = KD_I3*MBASE(I_MACH)/SBASE

	
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
				WRITE (LPDEV,*) 'VSCGFO - CASE 1: No governor-type model found at bus ',NUMBUS(IB),' with id ',MACHID(I_MACH),'.'
				! I_STATE_DCGRID = -1
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
			IF (ict.GT.icmaxsspu) THEN
				WRITE (LPDEV,*) 'VSCGFO - CASE 1: Current exceeds steady-state limit of ',icmaxsspu,' pu.'
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
					
			! State variables
			! ---------------			
			xomega = 1.0
			d_xecd = 0.0
			d_xecq = 0.0
            d_xdelta = 0.0
            d_xomega = 0.0
			d_xdc = 0.0	
			d_xpod = 0.0
			d_xpfc = 0.0
			d_xq = 0.0
			d_xitvid = 0.0
			d_xitviq = 0.0
			d_xudcpss = 0.0
			ypod = 0.0
			yitvid = 0.0
			yitviq = 0.0
			yudcpss = 0.0
			ydeltamax = 10E9!asin(AIMAG(zc)*icmaxsspu)+thetas
			ydeltamin = -10E9!-asin(AIMAG(zc)*icmaxsspu)+thetas
			d_xdeltamax = 0.0
			d_xdeltamin = 0.0

			! DC voltage control
			CALL SUB_PI(ydc,xdc,d_xdc,udc**2-udcref**2,1,DELTAT,kdp2pu,kdi2pu,pcmaxpu,pcminpu)
			
			! virtual synchronous machine
			! To be done: add 
			! - virtual active power control (VAPC, Kanakesh, 2021)
			! - overload mitigation (Du, 2019)
			! - check stop inertia integration
			CALL SUB_WASHOUTWINDUP(ypod,xpod,d_xpod,-(xomega-1.0),1,DELTAT,TPOD,10e9,-10e9)
			CALL SUB_FIRSTORDERWINDUP(xpfc,xpfc,d_xpfc,-(xomega-1.0),1,-1,DELTAT,TPFC,(pcmaxpu-pctref)/kpfcpu,(pcminpu-pctref)/kpfcpu)
			
			!CALL SUB_PI(ydeltamax,xdeltamax,d_xdeltamax,(pcmaxpu/ect-ictd),1,DELTAT,0.0,kdi3pu,asin(AIMAG(zc*icmaxsspu))+thetas,thetas)
			!CALL SUB_PI(ydeltamin,xdeltamin,d_xdeltamin,(-pcmaxpu/ect-ictd),1,DELTAT,0.0,kdi3pu,thetas,-asin(AIMAG(zc*icmaxsspu))+thetas)
			
			deltap = min(max(ydc + deltapctref  + kpfcpu*xpfc + kpodpu*ypod - dpu*(xomega-1.0),pcminpu),pcmaxpu) - (pct + ploss)	
			CALL SUB_INTEGRATORWINDUP(xomega,xomega,d_xomega,deltap,1,hpu,1.0+DELTAWMAX,1.0-DELTAWMAX)
            CALL SUB_INTEGRATORWINDUP(xdelta,xdelta,d_xdelta,(xomega-1.0),1,(1/OMEGABASE),ydeltamax,ydeltamin)

			! AC voltage control
			deltaq = qctref + deltaqctref - qct
			CALL SUB_FIRSTORDERWINDUP(xq,xq,d_xq,deltaq,1,-1,DELTAT,TQ,(qcmaxpu-qctref)/kqpu,(qcminpu-qctref)/kqpu)
			
			! transient virtual resistance
			CALL SUB_WASHOUTWINDUP(yitvid,xitvid,d_xitvid,ictd,1,DELTAT,TTVR,10e9,-10e9)						
			CALL SUB_WASHOUTWINDUP(yitviq,xitviq,d_xitviq,ictq,1,DELTAT,TTVR,10e9,-10e9)						
			
			! UDC PSS
			CALL SUB_WASHOUTWINDUP(yudcpss,xudcpss,d_xudcpss,udcref,1,DELTAT,TUDCPSS,10e9,-10e9)				
			
			! AC voltage dynamics
			deltaecd = ectref - rtvrpu*yitvid - kqpu*xq - KDCPSS*yudcpss
			deltaecq = -rtvrpu*yitviq
			CALL SUB_FIRSTORDERWINDUP(xecd,xecd,d_xecd,deltaecd,1,-1,DELTAT,TAU,10e9,-10e9)
			CALL SUB_FIRSTORDERWINDUP(xecq,xecq,d_xecq,deltaecq,1,-1,DELTAT,TAU,10e9,-10e9)  	
			
			WRITE (LPDEV,*) 'VSCGFO - CASE 1: Converter ',IDXCONVERTER,' at bus ',NUMBUS(IB),' with id ',MACHID(I_MACH),' initialized. Initial conditions of states K+11/12 might be suspect.'  
			WRITE (LPDEV,*) 'VSCGFO - CASE 1: delta = ',xdelta,', deltamax = ',(asin(AIMAG(zc))+thetas),', thetas = ',thetas
			
			EFD(I_MACH) = ectref
            SPEED(I_MACH) = xomega-1.0
            ANGLE(I_MACH) = xdelta*180/PI
											
		CASE (2) 
		   
			! Compute derivatives
			! ===================	
			IF (I_STATE_DCGRID.GT.0) THEN
				udc = STATE(I_STATE_DCGRID + IDXCONVERTER + NDCLINES_PREVIOUS + NDCBUS_PREVIOUS - 1)
			ELSE
				udc = udcref
			END IF

			CALL SUB_COMPUTEPLOSS(ploss, pct, ict, aloss, bloss, c_inv, c_rect)
			! VI immpedance current limiter	modifying dq voltage references
			! If zvi = kpvipu*(1+jSIGMAXR) = zc*exp(j*phi), VI impedance and hard current limit are the same		
			CALL SUB_ZVILIMITSECREF(yecdviimax, yecqviimax, ict_phasor_dq, icmaxsspu, kpvipu, SIGMAXR)
			
			deltap = min(max(ydc + deltapctref  + kpfcpu*xpfc + kpodpu*ypod - dpu*(xomega-1.0),pcminpu),pcmaxpu) - (pct + ploss)
			deltaq = qctref + deltaqctref - qct
			deltaecd = ectref - rtvrpu*yitvid - kqpu*xq - KDCPSS*yudcpss - yecdviimax
			deltaecq = -rtvrpu*yitviq - yecqviimax

			CALL SUB_PI(ydc,xdc,d_xdc,udc**2-udcref**2,2,DELTAT,kdp2pu,kdi2pu,pcmaxpu,pcminpu)
			CALL SUB_WASHOUTWINDUP(ypod,xpod,d_xpod,-(xomega-1.0),2,DELTAT,TPOD,10e9,-10e9)
			CALL SUB_FIRSTORDERWINDUP(xpfc,xpfc,d_xpfc,-(xomega-1.0),2,0,DELTAT,TPFC,(pcmaxpu-pctref)/kpfcpu,(pcminpu-pctref)/kpfcpu)
			
			!CALL SUB_PI(ydeltamax,xdeltamax,d_xdeltamax,(pcmaxpu/ect-ictd),2,DELTAT,0.0,kdi3pu,asin(AIMAG(zc*icmaxsspu))+thetas,thetas)
			!CALL SUB_PI(ydeltamin,xdeltamin,d_xdeltamin,(-pcmaxpu/ect-ictd),2,DELTAT,0.0,kdi3pu,thetas,-asin(AIMAG(zc*icmaxsspu))+thetas)

			CALL SUB_INTEGRATORWINDUP(xomega,xomega,d_xomega,deltap,2,hpu,1.0+DELTAWMAX,1.0-DELTAWMAX)
            CALL SUB_INTEGRATORWINDUP(xdelta,xdelta,d_xdelta,(xomega-1.0),2,(1/OMEGABASE),ydeltamax,ydeltamin)
			
			CALL SUB_FIRSTORDERWINDUP(xq,xq,d_xq,deltaq,2,0,DELTAT,TQ,(qcmaxpu-qctref)/kqpu,(qcminpu-qctref)/kqpu)
			CALL SUB_WASHOUTWINDUP(yitvid,xitvid,d_xitvid,ictd,2,DELTAT,TTVR,10e9,-10e9)						
			CALL SUB_WASHOUTWINDUP(yitviq,xitviq,d_xitviq,ictq,2,DELTAT,TTVR,10e9,-10e9)						
			CALL SUB_WASHOUTWINDUP(yudcpss,xudcpss,d_xudcpss,udc,2,DELTAT,TUDCPSS,10e9,-10e9)
			
			CALL SUB_FIRSTORDERWINDUP(xecd,xecd,d_xecd,deltaecd,2,0,DELTAT,TAU,10e9,-10e9)
			CALL SUB_FIRSTORDERWINDUP(xecq,xecq,d_xecq,deltaecq,2,0,DELTAT,TAU,10e9,-10e9)

			! PELEC(I_MACH) = ps
            ! QELEC(I_MACH) = qs
			! ETERM(I_MACH) = us

		CASE (3) 
		
			! Compute output
			! ==============
			IF (I_STATE_DCGRID.GT.0) THEN
				udc = STATE(I_STATE_DCGRID + IDXCONVERTER + NDCLINES_PREVIOUS + NDCBUS_PREVIOUS - 1)
			ELSE
				udc = udcref
			END IF

			CALL SUB_COMPUTEPLOSS(ploss, pct, ict, aloss, bloss, c_inv, c_rect)
			! VI immpedance current limiter	modifying dq voltage references
			CALL SUB_ZVILIMITSECREF(yecdviimax, yecqviimax, ict_phasor_dq, icmaxsspu, kpvipu, SIGMAXR)
			
			deltap = min(max(ydc + deltapctref  + kpfcpu*xpfc + kpodpu*ypod - dpu*(xomega-1.0),pcminpu),pcmaxpu) - (pct + ploss)
			deltaq = qctref + deltaqctref - qct
			deltaecd = ectref - rtvrpu*yitvid - kqpu*xq - KDCPSS*yudcpss - yecdviimax
			deltaecq = -rtvrpu*yitviq - yecqviimax

			! Remember that the output is the ISORCE given in the RI refrence frame of the system
			CALL SUB_PI(ydc,xdc,d_xdc,udc**2-udcref**2,3,DELTAT,kdp2pu,kdi2pu,pcmaxpu,pcminpu)
			CALL SUB_WASHOUTWINDUP(ypod,xpod,d_xpod,-(xomega-1.0),3,DELTAT,TPOD,10e9,-10e9)
			CALL SUB_FIRSTORDERWINDUP(xpfc,xpfc,d_xpfc,-(xomega-1.0),3,0,DELTAT,TPFC,(pcmaxpu-pctref)/kpfcpu,(pcminpu-pctref)/kpfcpu)

			!CALL SUB_PI(ydeltamax,xdeltamax,d_xdeltamax,(pcmaxpu/ect-ictd),3,DELTAT,0.0,kdi3pu,asin(AIMAG(zc*icmaxsspu))+thetas,thetas)
			!CALL SUB_PI(ydeltamin,xdeltamin,d_xdeltamin,(-pcmaxpu/ect-ictd),3,DELTAT,0.0,kdi3pu,thetas,-asin(AIMAG(zc*icmaxsspu))+thetas)

			CALL SUB_INTEGRATORWINDUP(xomega,xomega,d_xomega,deltap,3,hpu,1.0+DELTAWMAX,1.0-DELTAWMAX)
            CALL SUB_INTEGRATORWINDUP(xdelta,xdelta,d_xdelta,(xomega-1.0),3,(1/OMEGABASE),ydeltamax,ydeltamin)
			
			CALL SUB_FIRSTORDERWINDUP(xq,xq,d_xq,deltaq,3,0,DELTAT,TQ,(qcmaxpu-qctref)/kqpu,(qcminpu-qctref)/kqpu)
			CALL SUB_WASHOUTWINDUP(yitvid,xitvid,d_xitvid,ictd,3,DELTAT,TTVR,10e9,-10e9)						
			CALL SUB_WASHOUTWINDUP(yitviq,xitviq,d_xitviq,ictq,3,DELTAT,TTVR,10e9,-10e9)						
			CALL SUB_WASHOUTWINDUP(yudcpss,xudcpss,d_xudcpss,udc,3,DELTAT,TUDCPSS,10e9,-10e9)
			
			CALL SUB_FIRSTORDERWINDUP(xecd,xecd,d_xecd,deltaecd,3,0,DELTAT,TAU,10e9,-10e9)
			CALL SUB_FIRSTORDERWINDUP(xecq,xecq,d_xecq,deltaecq,3,0,DELTAT,TAU,10e9,-10e9)
			ect_phasor_dq = CMPLX(xecd,xecq)

            ! RI -> dq reference system transformation: xRI = reference_transf*xdq		
	        reference_transf = CMPLX(COS(xdelta),SIN(xdelta))
			
			! Hard current limit approximation (should be actually placed in the network solution model due to the voltage dependency)	
			! ict_phasor_dq = (ect_phasor_dq - us_phasor_RI/reference_transf)/zc
			! CALL SUB_ICMAXLIMITSECREF(ect_phasor_dq, ict_phasor_dq, (us_phasor_RI/reference_transf), ILIMITPRIORITY, icmaxsspu, zc)		
			! ect_phasor_RI = ect_phasor_dq*reference_transf
            ! ISORCE(I_MACH) = ect_phasor_RI/zc       

			
            EFD(I_MACH) = ABS(ect_phasor_dq) ! xecd
            SPEED(I_MACH) = xomega-1.0
            ANGLE(I_MACH) = xdelta*180/PI
			
		CASE (4) 

			! Update number of STATEs.
			! ========================
			NINTEG = MAX(NINTEG,I_STATE+12)

		CASE (5)
			! reporting mode
			WRITE (LPDEV,*) 'VSCGFO model at bus ', NUMBUS(IB)
				
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
	VAR(I_VAR+7) = ypod
	VAR(I_VAR+8) = yitvid
	VAR(I_VAR+9) = yitviq
	VAR(I_VAR+10) = yudcpss
	VAR(I_VAR+11) = ydeltamax
	VAR(I_VAR+12) = ydeltamin

	VAR(I_VAR+13) = ictd
	VAR(I_VAR+14) = ictq
	VAR(I_VAR+15) = deltap
	VAR(I_VAR+16) = (pct+ploss)
		
	! STATEs
	STATE(I_STATE) = xecd         
	STATE(I_STATE+1) = xecq       
    STATE(I_STATE+2) = xdelta  
    STATE(I_STATE+3) = xomega  
	STATE(I_STATE+4) = xdc
	STATE(I_STATE+5) = xpod
	STATE(I_STATE+6) = xpfc
	STATE(I_STATE+7) = xq
	STATE(I_STATE+8) = xitvid
	STATE(I_STATE+9) = xitviq
	STATE(I_STATE+10) = xudcpss
	STATE(I_STATE+11) = xdeltamax
	STATE(I_STATE+12) = xdeltamin

	! DSTATEs
	DSTATE(I_STATE) = d_xecd      
	DSTATE(I_STATE+1) = d_xecq   
    DSTATE(I_STATE+2) = d_xdelta  
    DSTATE(I_STATE+3) = d_xomega 
	DSTATE(I_STATE+4) = d_xdc
	DSTATE(I_STATE+5) = d_xpod
	DSTATE(I_STATE+6) = d_xpfc
	DSTATE(I_STATE+7) = d_xq
	DSTATE(I_STATE+8) = d_xitvid
	DSTATE(I_STATE+9) = d_xitviq
	DSTATE(I_STATE+10) = d_xudcpss
	DSTATE(I_STATE+11) = d_xdeltamax
	DSTATE(I_STATE+12) = d_xdeltamin
      
END SUBROUTINE VSCGFO


