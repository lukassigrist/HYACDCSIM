! Dynamic model of a VSC converter of a VSC-based MTDC grid
! 
! The VSC converter computes the currents injected into the ac grid and the currents injected into the DC grid. The DC grid
! model, implemented as a governor type model, takes as input the DC current injected by VSC.
!
! The VSC converter model is based on a grid-following converter. The current control of the converter is simplified
! and reduced to a first-order transfer function. 
! 
! Supervisory controls such as active and reactive power controls, DC-voltage control, etc. are implemented in further
! user models. These user models are implemented as stabilizer (SQWDRD), excitation (SPWDRD), and governor type models 
! (DCGRID). Note that the calling sequence for the plant models is: 1. Generator models, 2. Current compensating models, 
! 3. Excitation stabilizer models, 4. Excitation system models, 5. Turbine-governor models. Turbine governor, stabilizer
! and excitation limiter models have no initialization duties other than STATEs and VARs.
! 
! The VSC converter model, VSCGFL, is implemented as a coordinated-call, current injecting generator model (IC = 1 
! and IT = 1). Coordinated-call implementation is required to cancel the effect of the admittance of the Norton
! equivalent applied to generators, which leads to a pure current source.  

! ======================================================================================
! MODULE DECLARATION
! ======================================================================================
MODULE MOD_VSCGFL_INTERNAL
	! Internal variables
	INTEGER :: I_STATE_DCGRID_G = -1
END MODULE MOD_VSCGFL_INTERNAL
   
! ======================================================================================
! CURRENT INJECTIONS FOR NETWORK SOLUTION
! ======================================================================================
SUBROUTINE TSCGFL(I_MACH,I_SLOT)
    
	INCLUDE 'COMON4.ins'					! common PSS/e variables and modules
	IMPLICIT NONE

	! Declaration
	! ----------- 
	INTEGER :: I_MACH, I_SLOT, IB, I_VAR		
	  
	REAL :: deltas, icd, icq 			 

	COMPLEX :: us_phasor_RI, ic_phasor_RI, ss_phasor 
	COMPLEX :: zc, reference_transf

	! Parameter assignment
	! --------------------
	! Index
	IB = NUMTRM(I_MACH) 					! Bus sequence number
	I_VAR = STRTIN(3,I_SLOT)				! initial VAR index

	! VARs
	icd=VAR(I_VAR+8) 
	icq=VAR(I_VAR+9) 

	! Common variables
	us_phasor_RI = VOLT(IB)	 
	zc = ZSORCE(I_MACH)*SBASE/MBASE(I_MACH)

	! RI -> dq reference system transformation: xRI = reference_transf*xdq
	deltas = ATAN2(AIMAG(us_phasor_RI),REAL(us_phasor_RI))			
	reference_transf = CMPLX(COS(deltas),SIN(deltas))

	SELECT CASE (MODE)
		CASE (1) 
			! Initialization
			! --------------
			! ec = ut + ZSORCE*ic, ISORCE = ec/ZSORCE -> ic = ISORCE - ut/ZSORCE
			! ISORCE: norton equivalent source current in pu of SBASE
			! ic: current at generator terminal bus in pu SBASE
			! ut: voltage at generator terminal bus in pu
			! ZSORCE: in pu of MBASE -> zc = ZSORCE*SBASE/MBASE
			ic_phasor_RI = ISORCE(I_MACH) - us_phasor_RI/zc 
			icd = REAL(ic_phasor_RI/reference_transf)
			icq = AIMAG(ic_phasor_RI/reference_transf) 
						 
			ss_phasor = us_phasor_RI*CONJG(ic_phasor_RI) ! in pu SBASE     

		CASE (3) 
			! Current injections
			! ------------------         
			! Set ISORCE such that the current through 1/ZSORCE is cancelled 
			! Note that the reference transformation should be done here (and not in mode 3 since the voltage 
			! bus angle depends on the TYSL iteration
			ic_phasor_RI = reference_transf*CMPLX(icd,icq)
			ISORCE(I_MACH) = ic_phasor_RI + us_phasor_RI/zc  

			ss_phasor = us_phasor_RI*CONJG(ic_phasor_RI)   
			       
		CASE DEFAULT

	END SELECT

	! Re-assign algebraic variables
	! -----------------------------
	PELEC(I_MACH) = REAL(ss_phasor)
	QELEC(I_MACH) = AIMAG(ss_phasor)
	ETERM(I_MACH) = ABS(us_phasor_RI)

	! VARs
	VAR(I_VAR+8) = icd
	VAR(I_VAR+9) = icq 
	 
END SUBROUTINE TSCGFL
      
! ==========================================================
! SOLUTION OF DIFFERENTIAL EQUATIONS
! ==========================================================      
SUBROUTINE VSCGFL(I_MACH,I_SLOT)

	USE MOD_TFBLOCKS 								! use module of generic transfer functions
	USE MOD_PROTECTION								! use module of protection functions
	USE MOD_READHYADCSIM							! use module for reading the HYACDCSIM text files
	USE MOD_MISC									! use module for miscellaneous functions
	USE MOD_VSCGFL_INTERNAL							! use module for global VSCGFL-related variables
	INCLUDE 'COMON4.INS'							! common PSS/e variables and modules
	IMPLICIT NONE
	
	! Declaration
	! -----------
	INTEGER IB, I_SLOT, I_MACH, I_VAR, I_CON, I_ICON, I_STATE       ! PSS/e indices
	INTEGER I_STATE_DCGRID      									! local index of the DCGRID STATE
	INTEGER IDXCONVERTER, DCNTRLTYPE, QCNTRLTYPE, ILIMITPRIORITY 
	
	INTEGER istripvsc                                             ! disconnects unit and sets current injections to 0
	INTEGER dcontroltype_aux, qcontroltype_aux     
	INTEGER NDCBUS 												  ! # of buses of the dc-grid
	INTEGER NDCBUS_PREVIOUS, NDCLINES_PREVIOUS
	
	INTEGER, ALLOCATABLE :: m_VSCACDCBUS(:,:) ! ac and dc buses of each converter
	INTEGER ierr
	REAL PI
	PARAMETER (PI=3.14159265358979)
	REAL DELTAT
		
	REAL us, deltas, deltas0
	REAL ps0, qs0 
	REAL icd, icq ! Currents in dq.
	REAL icd0, icq0 
	REAL icD_out, icQ_out ! Currents in DQ.
	REAL ps, qs

	REAL xd, xq 		! state variables: icd and icq
	REAL md, mq 		! state variables: integrals for active and reactive power control
	REAL ndc, nq 		! state variables: integrals for dc and ac voltage controls
	REAL eta_d 			! state variables: integrals angle control
	REAL d_xd, d_xq 	
	REAL d_md, d_mq 	
	REAL d_ndc, d_nq 	
	REAL d_eta_d 

	REAL aux_ilimit, antiwindupicd, antiwindupicq        

	REAL ic, ic_abs_puconv ! Current magnitude

	REAL psref, qsref, usref ! references
	REAL icdref, icqref
	REAL udcref

	REAL idc, idc0, udc, udc0
	REAL pdc

	REAL deltapsref 
	REAL deltaqsref       
	REAL ec, deltac      
	REAL pvsc, qvsc

	COMPLEX us_phasor_RI, ec_phasor_RI, ic_phasor_RI, sc_phasor ! bus voltage.
	COMPLEX reference_transf ! convert inverter reference (dq) to system reference (DQ).

	CHARACTER(1) IDGRID

	REAL, ALLOCATABLE :: v_udc0(:,:)

	REAL TAU	
	REAL KD_Pps, KD_Pudc, KD_Ips, KD_Iudc, KD_D2, KQ_Pqs, KQ_Iqs, KQ_Pus, KQ_Ius
	REAL ICMAX_PUconv
	REAL PS_MAX_MW, PS_MIN_MW, QS_MAX_Mvar, QS_MIN_Mvar
	REAL ALOSS_MW, BLOSS_kV, CLOSS_RECT_Ohm, CLOSS_INV_Ohm
	REAL UDC_MAX, UDC_MIN, TUDCMAX, TUDCMIN

	REAL icmax
	REAL psmax, psmin, qsmax, qsmin
	REAL aloss, bloss, c_rect, c_inv
	REAL CDC_uF, Cdc
	REAL UDC_NOMINAL_kV
	REAL ZBASE, ZDCBASE, SDCBASE_MVA
	COMPLEX zc
	REAL ploss
	
	INTEGER :: counter_uvdc, counter_ovdc
	REAL :: tinitial_uvdc, tinitial_ovdc

	REAL fmodulationpwm, FMODULATIONPWMMAX

	REAL :: yetad, yndc, ymd, ynq, ymq

	

	! COMMON /svsconinternal/ I_STATE_DCGRID_G, aux_var_GLOBAL ! This allows saving the value of I_STATE_DCGRID_G between calls. Note that I_STATE_DCGRID_G takes the value of the first DCGRID model read.

	! Parameter assignment
	! --------------------	  
	! Indexes
	IB=NUMTRM(I_MACH)
	I_CON = STRTIN(1,I_SLOT)
	I_STATE = STRTIN(2,I_SLOT) ! Model Calling Sequence Rules
	I_VAR = STRTIN(3,I_SLOT)
	I_ICON = STRTIN(4,I_SLOT)
		
	! CONs    
	TAU = CON(I_CON) 		! Inverter time constant (e.g., 0.1)
	KD_Pps = CON(I_CON+1) 	! Ps control
	KD_Ips = CON(I_CON+2)
	KD_Pudc = CON(I_CON+3) 	! Udc control 
	KD_Iudc = CON(I_CON+4)
	KD_D2 = CON(I_CON+5) ! gain of the differential control of Udc -> for the DC-voltage control, PID works much better than PI
	KQ_Pqs = CON(I_CON+6) 	! Qs control
	KQ_Iqs = CON(I_CON+7)
	KQ_Pus = CON(I_CON+8) 	! Us control
	KQ_Ius = CON(I_CON+9)
	ICMAX_PUconv = CON(I_CON+10) ! Maximum inverter current/susceptance in pu with respect to inverter rating (e.g., 1.1)
	PS_MAX_MW = CON(I_CON+11)
	PS_MIN_MW = CON(I_CON+12)
	QS_MAX_Mvar = CON(I_CON+13)
	QS_MIN_Mvar = CON(I_CON+14)
	UDC_MAX = CON(I_CON+15)
	UDC_MIN = CON(I_CON+16)
	FMODULATIONPWMMAX = CON(I_CON+17)
	ALOSS_MW = CON(I_CON+18) 		! constant converter loss coefficient (MW): ploss = aloss + bloss*ic + c*ic^2
	BLOSS_kV = CON(I_CON+19) 		! linear converter loss coefficient (kV): ploss = aloss + bloss*ic + c*ic^2
	CLOSS_RECT_Ohm = CON(I_CON+20) 	! rectifier quadratic converter loss coefficient (ohm): ploss = aloss + bloss*ic + c*ic^2
	CLOSS_INV_Ohm = CON(I_CON+21) 	! constant converter loss coefficient (ohm): ploss = aloss + bloss*ic + c*ic^2
	CDC_uF = CON(I_CON+22) 			! capacitor of the converter (micro-Faraday) (uF)
	UDC_NOMINAL_kV = CON(I_CON+23) 	! nominal dc-voltage of the converter (kV)
	TUDCMIN = CON(I_CON+24) 		! DC undervoltage protection delay
	TUDCMAX = CON(I_CON+25) 		! DC overvoltage protection delay

	! ICONs
	IDXCONVERTER = ICON(I_ICON) 		! absolute index number of all SVSCON
	DCNTRLTYPE = ICON(I_ICON+1) 		! Ps-control: 1 , Udc-control: 2, Deltas-control: 3
	QCNTRLTYPE = ICON(I_ICON+2) 		! Qs-control: 1 , Us-control: 2
	ILIMITPRIORITY = ICON(I_ICON+3) 	! Current limit: P-priority: 1, Q-priority: 2 and P-Q equal priority: 3 or any other integer
	NDCBUS = ICON(I_ICON+4) 			! Number of converters of the DC grid
	IDGRID = CHRICN(I_ICON+5)			! DC grid identifier
	NDCBUS_PREVIOUS = ICON(I_ICON+6) 	! accumulated number of DC buses
	NDCLINES_PREVIOUS = ICON(I_ICON+7) 	! accumulated numnber of DC lines

	IF (DCNTRLTYPE.EQ.3) THEN ! if feeding a passive grid, reactive power control is a voltage control
		QCNTRLTYPE = 2
	END IF
			
	! VARs
	ps0 = VAR(I_VAR) ! Reference active power, coincides with load flow values.
	qs0 = VAR(I_VAR+1) 
	icD_out = VAR(I_VAR+2) 
	icQ_out = VAR(I_VAR+3) 
	ps = VAR(I_VAR+4) ! Active model output power. 
	qs = VAR(I_VAR+5) 
	icd0 = VAR(I_VAR+6)
	icq0 = VAR(I_VAR+7)
	icd = VAR(I_VAR+8)
	icq = VAR(I_VAR+9)
	istripvsc = VAR(I_VAR+10)

	udcref = VAR(I_VAR+13)

	usref = VAR(I_VAR+15)
	aux_ilimit = VAR(I_VAR+16)
	antiwindupicd = VAR(I_VAR+17)
	antiwindupicq = VAR(I_VAR+18)
	ic_abs_puconv = VAR(I_VAR+19)
	dcontroltype_aux = VAR(I_VAR+20)
	qcontroltype_aux = VAR(I_VAR+21)
	udc = VAR(I_VAR+22)
	idc = VAR(I_VAR+23)
	idc0 = VAR(I_VAR+24)
	pdc = VAR(I_VAR+25)
	deltas0 = VAR(I_VAR+26)
	udc0 = VAR(I_VAR+27)
	icdref = VAR(I_VAR+28)
	deltapsref = VAR(I_VAR+29)
	icqref = VAR(I_VAR+30)
	deltaqsref = VAR(I_VAR+31)
	fmodulationpwm = VAR(I_VAR+32)
	counter_uvdc = VAR(I_VAR+33)
	tinitial_uvdc = VAR(I_VAR+34)
	counter_ovdc = VAR(I_VAR+35)
	tinitial_ovdc = VAR(I_VAR+36)

	! STATEs
	xd = STATE(I_STATE)      ! icd
	xq = STATE(I_STATE+1)    ! icq
	md = STATE(I_STATE+2)    ! d-integral state variable type1
	mq = STATE(I_STATE+3)    ! q-integral state variable type1
	ndc = STATE(I_STATE+4)   ! d-integral state variable type2 
	nq = STATE(I_STATE+5)    ! q-integral state variable type2
	eta_d = STATE(I_STATE+6) ! d-integral state variable type3 (passive grid)

	! DSTATEs
	d_xd = DSTATE(I_STATE)      
	d_xq = DSTATE(I_STATE+1)    
	d_md = DSTATE(I_STATE+2)    
	d_mq = DSTATE(I_STATE+3)    
	d_ndc = DSTATE(I_STATE+4)   
	d_nq = DSTATE(I_STATE+5)    
	d_eta_d = DSTATE(I_STATE+6) 

	! Common variables
	CALL DSRVAL('DELT', 1, DELTAT, IERR)
	
	ALLOCATE(v_udc0(NDCBUS,1)) ! allocate size
	ALLOCATE(m_VSCACDCBUS(NDCBUS,2))

	icmax = ICMAX_PUconv*MBASE(I_MACH)/SBASE ! current limit in system base
	psmax = PS_MAX_MW/SBASE               ! power limits in system base
	psmin = PS_MIN_MW/SBASE
	qsmax = QS_MAX_Mvar/SBASE
	qsmin = QS_MIN_Mvar/SBASE

	ZBASE = BASVLT(IB)**2/SBASE               ! AC-side Zbase in ohms
	SDCBASE_MVA = SBASE
	ZDCBASE = 2*UDC_NOMINAL_kV**2/SDCBASE_MVA ! DC-side Zbase in ohms

	aloss = ALOSS_MW/SBASE                  ! losses coefs in system base
	bloss = BLOSS_kV/SQRT(3.0)/BASVLT(IB)
	c_rect = CLOSS_RECT_Ohm/ZBASE
	c_inv = CLOSS_INV_Ohm/ZBASE
	zc = ZSORCE(I_MACH)*SBASE/MBASE(I_MACH)       ! conexion impedance in system base
	Cdc = (CDC_uF*1e-6)*ZDCBASE             ! [pu*s]

	I_STATE_DCGRID = I_STATE_DCGRID_G
	
	! ec = us + ZSORCE*ic, ISORCE = ec/ZSORCE -> ic = ISORCE - us/ZSORCE
	! ISORCE: norton equivalent source current in pu of SBASE
	! ic: current at generator terminal bus in pu SBASE
	! us: voltage at generator terminal bus in pu
	! ZSORCE: in pu of MBASE -> zc = ZSORCE*SBASE/MBASE
	us_phasor_RI = VOLT(IB)
	us = ABS(us_phasor_RI)
	ic_phasor_RI = ISORCE(I_MACH) - us_phasor_RI/zc
	

	! DQ -> dq reference system transformation: xDQ = reference_transf*xdq
	deltas = ATAN2(AIMAG(us_phasor_RI),REAL(us_phasor_RI))		
	IF (DCNTRLTYPE.EQ.3) THEN ! If feeding a passive grid, the VSC is a AC slack bus
		reference_transf = CMPLX(1.0, 0.0) 
	ELSE
		reference_transf = CMPLX(cos(deltas),sin(deltas))
	END IF

	icd = REAL(ic_phasor_RI/reference_transf) 
	icq = AIMAG(ic_phasor_RI/reference_transf) 

	! present power and current references
	psref = MIN(MAX(ps0 + deltapsref,psmin),psmax)
	qsref = MIN(MAX(qs0 + deltaqsref,qsmin),qsmax)

	SELECT CASE (MODE)

		CASE (1) 
			   
			! Inicialization 
			! ==============	  
			
			! Get initial value of the DC voltage state from .txt files
			CALL SUB_READDCBUS(v_udc0, m_VSCACDCBUS, IDGRID, NDCBUS)
			
			! Get STATE index position of the dc grid model DCGRID (governor-type model)	  
			CALL MDLIND(m_VSCACDCBUS(1,2), MACHID(I_MACH), 'GOV', 'STATE', I_STATE_DCGRID, ierr)
			IF (ierr.NE.0) THEN ! No governor-type model
				WRITE (LPDEV,*) 'VSCGFL - CASE 1: No governor-type model found at bus ',NUMBUS(IB),' with id ',MACHID(I_MACH),'.'
			END IF
			IF (I_STATE_DCGRID_G.LT.0) THEN
				I_STATE_DCGRID_G = I_STATE_DCGRID
			END IF

			! Algebraic variables (VARS)
			! --------------------------			
			deltas0 = deltas
			icd0 = icd
			icq0 = icq
			ps0 = us*icd
			qs0 = -us*icq
			
			! AC-side converter current out of range       
			ic = SQRT(icd**2+icq**2) 
			IF (ic.GT.icmax) THEN
				WRITE (LPDEV,*) 'VSCGFL::CASE1: Current exceeds limit of ',icmax,' pu.'
			END IF
			
			! AC-side converter voltages, power
			ec_phasor_RI = us_phasor_RI + zc*ic_phasor_RI
			ec = ABS(ec_phasor_RI)
			pvsc = REAL(ec_phasor_RI*CONJG(ic_phasor_RI))
			pvsc = ps0 + REAL(zc)*(ic**2) 

			! DC-side voltage, current, power
			CALL SUB_COMPUTEPLOSS(ploss, ps0, ic, aloss, bloss, c_inv, c_rect)
			pdc = -(pvsc + ploss)
			udc = v_udc0(IDXCONVERTER,1) ! extract initial dc voltage value
			udc0 = udc
			idc0 = pdc/udc0
			idc = idc0
		
			fmodulationpwm = ec/MAX(udc0,0.0001)
			IF (fmodulationpwm.GT.FMODULATIONPWMMAX) THEN
				WRITE (LPDEV,*) 'VSCGFL::CASE1: Modulation index exceeds limit of ',FMODULATIONPWMMAX,'.'
			END IF

			! references
			psref = ps0
			qsref = qs0
			usref = us
			udcref = udc0 
			deltapsref = 0.0 ! supplementary controller input (SPWDRD)
			deltaqsref = 0.0 ! supplementary controller input (SQWDRD)
			icdref = icd0
			icqref = icq0
			
			! indexes (stop simulation and anti-windup)
			istripvsc = 0
			antiwindupicd = 1.0
			antiwindupicq = 1.0

			! State variables (STATES)
			! ------------------------
			! Current control
			d_xd = 0.0
			d_xq = 0.0
			CALL SUB_FIRSTORDERWINDUP(icdref,xd,d_xd,icdref,1,-1,DELTAT,TAU,icmax,-icmax) 				! icd
			CALL SUB_FIRSTORDERWINDUP(icqref,xq,d_xq,icqref,1,-1,DELTAT,TAU,icmax,-icmax) 				! icq

			! Active power related control	
			d_eta_d = 0.0
			d_ndc = 0.0
			d_md = 0.0		
			yetad = 0.0
			yndc = 0.0
			ymd = 0.0
			CALL SUB_PI(yetad,eta_d,d_eta_d,(deltas0-deltas),2,DELTAT,100.0,1000.0,icmax,-icmax)			! Deltas control					
			CALL SUB_PI(yndc,ndc,d_ndc,(udcref-udc),2,DELTAT,KD_Pudc,KD_Iudc,icmax,-icmax)				! DC-voltage control													
			CALL SUB_PI(ymd,md,d_md,(psref-ps),2,DELTAT,KD_Pps,KD_Ips,icmax,-icmax)						! Ps control

			! Reactive power related control	
			d_nq = 0.0
			d_mq = 0.0
			ynq = 0.0
			ymq = 0.0
			CALL SUB_PI(ynq,nq,d_nq,(usref-us),2,DELTAT,KQ_Pus,KQ_Ius,icmax,-icmax)						! Us control													
			CALL SUB_PI(ymq,mq,d_mq,(qsref-qs),2,DELTAT,KQ_Pqs,KQ_Iqs,icmax,-icmax)						! Qs control

			PMECH(I_MACH) = -pdc				! pu system rating
			SPEED(I_MACH) = BSFREQ(IB)            		
			ETERM(I_MACH) = us
			PELEC(I_MACH) = ps0            		! pu system rating
			QELEC(I_MACH) = qs0            		! pu system rating
			
			WRITE (LPDEV,*) 'VSCGFL - CASE 1: Converter ',IDXCONVERTER,' at bus ',NUMBUS(IB),' with id ',MACHID(I_MACH),' initialized.'
			WRITE (LPDEV,*) 'VSCGFL - CASE 1: D-control: ',DCNTRLTYPE,' Q-control: ',QCNTRLTYPE,' I-limit priority: ',ILIMITPRIORITY
			WRITE (LPDEV,*) 'VSCGFL - CASE 1: us = ',us,' pu, ec = ',ec,' pu'
			WRITE (LPDEV,*) 'VSCGFL - CASE 1: ps = ',ps0,' pu, qs = ', qs0, ' pu'
			WRITE (LPDEV,*) 'VSCGFL - CASE 1: pc = ',pvsc,', pu ploss = ',ploss,' pu'
			WRITE (LPDEV,*) 'VSCGFL - CASE 1: udc = ',udc,' pu, pdc = ',-pdc,' pu'
			WRITE (LPDEV,*) 'VSCGFL - CASE 1: idc = ',idc, ' pu'
			WRITE (LPDEV,*) 'VSCGFL - CASE 1: icd = ',icd,' pu, icq = ',icq, ' pu'
							
	CASE (2) 
	   
		! Compute derivatives
		! ===================

		! Get DC voltage of the DC grid model (if any, else from VARs)
		IF (I_STATE_DCGRID.GT.0) THEN
			udc = STATE(I_STATE_DCGRID + IDXCONVERTER + NDCLINES_PREVIOUS + NDCBUS_PREVIOUS - 1)
		END IF

		! Current control
		CALL SUB_FIRSTORDERWINDUP(icd,xd,d_xd,icdref,2,0,DELTAT,TAU,icmax,-icmax) 						! icd
		CALL SUB_FIRSTORDERWINDUP(icq,xq,d_xq,icqref,2,0,DELTAT,TAU,icmax,-icmax) 						! icq

		! Active power related control
		IF(DCNTRLTYPE.EQ.3) THEN 							
			CALL SUB_PI(yetad,eta_d,d_eta_d,(deltas0-deltas),2,DELTAT,100.0,1000.0,icmax,-icmax)		! Deltas control
		ELSE IF(DCNTRLTYPE.EQ.2) THEN						
			CALL SUB_PI(yndc,ndc,d_ndc,(udcref-udc),2,DELTAT,KD_Pudc,KD_Iudc,icmax,-icmax)				! DC-voltage control
		ELSE 													
			CALL SUB_PI(ymd,md,d_md,(psref-ps),2,DELTAT,KD_Pps,KD_Ips,icmax,-icmax)						! Ps control
		END IF

		! Reactive power related control	
		IF(QCNTRLTYPE.EQ.2) THEN
			CALL SUB_PI(ynq,nq,d_nq,(usref-us),2,DELTAT,KQ_Pus,KQ_Ius,icmax,-icmax)						! Us control
		ELSE 													
			CALL SUB_PI(ymq,mq,d_mq,(qsref-qs),2,DELTAT,KQ_Pqs,KQ_Iqs,icmax,-icmax)						! Qs control
		END IF

		! Apply anti-wind up indicators of current control: if current reaches a limit, temporarily cancel corresponding active and reactive power control
		d_ndc = d_ndc*antiwindupicd
		d_md = d_md*antiwindupicd
		d_eta_d = d_eta_d*antiwindupicd
		d_nq = d_nq*antiwindupicq
		d_mq = d_mq*antiwindupicq
		   
	CASE (3) 

		! Compute output
		! ==============

		! Get DC voltage into the DC grid (if any, else from VARs)
		IF (I_STATE_DCGRID.GT.0) THEN
			udc = STATE(I_STATE_DCGRID + IDXCONVERTER + NDCLINES_PREVIOUS + NDCBUS_PREVIOUS - 1)
		END IF

		! Current control
		CALL SUB_FIRSTORDERWINDUP(icd,xd,d_xd,icdref,3,0,DELTAT,TAU,icmax,-icmax) 						! icd
		CALL SUB_FIRSTORDERWINDUP(icq,xq,d_xq,icqref,3,0,DELTAT,TAU,icmax,-icmax) 						! icq

		! Active power related control
		IF(DCNTRLTYPE.EQ.3) THEN 							
			CALL SUB_PI(yetad,eta_d,d_eta_d,(deltas0-deltas),3,DELTAT,100.0,1000.0,icmax,-icmax)		! Deltas control
			icdref = icd0 + yetad
		ELSE IF(DCNTRLTYPE.EQ.2) THEN						
			CALL SUB_PI(yndc,ndc,d_ndc,(udcref-udc),3,DELTAT,KD_Pudc,KD_Iudc,icmax,-icmax)				! DC-voltage control
			icdref = icd0 - yndc
		ELSE 													
			CALL SUB_PI(ymd,md,d_md,(psref-ps),3,DELTAT,KD_Pps,KD_Ips,icmax,-icmax)						! Ps control
			icdref = psref/us + ymd
		END IF

		! Reactive power related control	
		IF(QCNTRLTYPE.EQ.2) THEN
			CALL SUB_PI(ynq,nq,d_nq,(usref-us),3,DELTAT,KQ_Pus,KQ_Ius,icmax,-icmax)						! Us control
			icqref = -MIN(MAX(-us*(icq0 - ynq),qsmin),qsmax)/us
		ELSE 													
			CALL SUB_PI(ymq,mq,d_mq,(qsref-qs),3,DELTAT,KQ_Pqs,KQ_Iqs,icmax,-icmax)						! Qs control
			icqref = -qsref/us - ymq
		END IF
		
		! limit current references and set anti-wind up indicators
		! CALL SUB_LIMITIREF(icdref, icqref, antiwindupicd, antiwindupicq, us, udc, ILIMITPRIORITY, icmax, FMODULATIONPWMMAX, zc)		
		CALL SUB_IMAXLIMITSIREF(icdref, icqref, antiwindupicd, antiwindupicq, ILIMITPRIORITY, icmax)
		CALL SUB_ECMAXLIMITSIREF(icdref, icqref, us, udc, FMODULATIONPWMMAX, zc)


		! DC protection
		CALL SUB_OVDCPROT(istripvsc, counter_ovdc, tinitial_ovdc, udc, TIME, UDC_MAX, TUDCMAX, NUMBUS(IB),'VSCGFL',LPDEV)
		CALL SUB_UVDCPROT(istripvsc, counter_uvdc, tinitial_uvdc, udc, TIME, UDC_MIN, TUDCMIN, NUMBUS(IB),'VSCGFL',LPDEV)
				
		IF (istripvsc.GT.0) THEN		  
			icd = 0.0
			icq = 0.0
		END IF
		
		! AC-side converter current, voltage, power 
		ps = us*icd 
		qs = -us*icq ! qs = AIMAG(ss_vec)	
		ic_phasor_RI = reference_transf*CMPLX(icd,icq) ! dq -> DQ, pu in system rating	
		ec_phasor_RI = us_phasor_RI + zc*ic_phasor_RI ! converter voltage -> output		
		sc_phasor = ec_phasor_RI*CONJG(ic_phasor_RI) ! interesa tenerlas, pero de momento no las usamos
		ic = ABS(ic_phasor_RI) ! mag
		icD_out = REAL(ic_phasor_RI) 
		icQ_out = AIMAG(ic_phasor_RI)
		ec = ABS(ec_phasor_RI)
		deltac = ATAN2(AIMAG(ec_phasor_RI), REAL(ec_phasor_RI))
		pvsc = REAL(sc_phasor)
		qvsc = AIMAG(sc_phasor)
		fmodulationpwm = ec/MAX(udc,0.0001)

		! DC-side converter current current and power
		CALL SUB_COMPUTEPLOSS(ploss, ps, ic, aloss, bloss, c_inv, c_rect)	
		pdc = -(pvsc + ploss) 
		idc = pdc/udc ! output of the converter model; input of the DC-grid model
		
		PMECH(I_MACH) = -pdc 						! in pu system rating
		SPEED(I_MACH) = BSFREQ(IB)
		ETERM(I_MACH) = us
		PELEC(I_MACH) = ps 							! in pu system rating
		QELEC(I_MACH) = -us*icq 					! in pu system rating
		
				  
	CASE (4) 

		! Update number of STATEs.
		! ========================
		NINTEG = MAX(NINTEG,I_STATE+6)

	CASE (5)
		! reporting mode
		WRITE (LPDEV,*) 'Converter ', IDXCONVERTER, ' at bus ', NUMBUS(IB)
		WRITE (LPDEV,*) 'SVSCON - ', 'CON', I_CON
		WRITE (LPDEV,*) 'SVSCON - ', 'ICON', I_ICON
		WRITE (LPDEV,*) 'SVSCON - ', 'VAR', I_VAR
		WRITE (LPDEV,*) 'SVSCON - ', 'STATE', I_STATE
			
	CASE DEFAULT

	END SELECT

	! RE-ASSIGN VARIABLES
	! -------------------	  
	icd = xd ! s�lo era para identificar mejor los estados
	icq = xq

	! VARs	  
	VAR(I_VAR) = ps0					! To SPWDRD (not really needed) 
	VAR(I_VAR+1) = qs0 					! To SQWDRD (not really needed)
	VAR(I_VAR+26) = deltas0
	VAR(I_VAR+6) = icd0
	VAR(I_VAR+7) = icq0

	VAR(I_VAR+13) = udcref				! To SPWDRD
	VAR(I_VAR+15) = usref
	VAR(I_VAR+28) = icdref 
	VAR(I_VAR+30) = icqref
	VAR(I_VAR+29) = deltapsref 			! From SPWDRD
	VAR(I_VAR+31) = deltaqsref			! From SQWDRD

	VAR(I_VAR+2) = icD_out 
	VAR(I_VAR+3) = icQ_out 
	VAR(I_VAR+4) = ps 
	VAR(I_VAR+5) = qs 
	VAR(I_VAR+8) = icd
	VAR(I_VAR+9) = icq

	VAR(I_VAR+10) = istripvsc
	VAR(I_VAR+16) = aux_ilimit
	VAR(I_VAR+17) = antiwindupicd
	VAR(I_VAR+18) = antiwindupicq
	VAR(I_VAR+19) = ic_abs_puconv
	VAR(I_VAR+20) = dcontroltype_aux
	VAR(I_VAR+21) = qcontroltype_aux

	VAR(I_VAR+25) = pdc					! To DCGRID
	VAR(I_VAR+22) = udc					! To SPWDRD
	VAR(I_VAR+23) = idc

	VAR(I_VAR+27) = udc0
	VAR(I_VAR+24) = idc0

	


	VAR(I_VAR+32) = fmodulationpwm
	VAR(I_VAR+33) = counter_uvdc
	VAR(I_VAR+34) = tinitial_uvdc
	VAR(I_VAR+35) = counter_ovdc
	VAR(I_VAR+36) = tinitial_ovdc
		
	! STATEs
	STATE(I_STATE) = xd      ! icd
	STATE(I_STATE+1) = xq    ! icq
	STATE(I_STATE+2) = md    ! d-integral state variable type1
	STATE(I_STATE+3) = mq    ! q-integral state variable type1
	STATE(I_STATE+4) = ndc   ! d-integral state variable type2 
	STATE(I_STATE+5) = nq    ! q-integral state variable type2
	STATE(I_STATE+6) = eta_d ! d-integral state variable type3 (passive grid)

	! DSTATEs
	DSTATE(I_STATE) = d_xd      
	DSTATE(I_STATE+1) = d_xq    
	DSTATE(I_STATE+2) = d_md    
	DSTATE(I_STATE+3) = d_mq   
	DSTATE(I_STATE+4) = d_ndc   
	DSTATE(I_STATE+5) = d_nq    
	DSTATE(I_STATE+6) = d_eta_d 
      
END SUBROUTINE VSCGFL

! ======================================================================================
! OTHER SUBROUTINES
! ======================================================================================


