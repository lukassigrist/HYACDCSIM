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
! The VSC converter model, SVSCON, is implemented as a coordinated-call, current injecting generator model (IC = 1 
! and IT = 1). Coordinated-call implementation is needed due to the representation of the generatorby a Norton 
! equivalent.   
   
! ======================================================================================
! CURRENT INJECTIONS FOR NETWORK SOLUTION
! ======================================================================================
SUBROUTINE TVSCON(I_MACH,I_SLOT)
      
	!DEC$ ATTRIBUTES DLLEXPORT, DECORATE, ALIAS: "TVSCON":: TVSCON
	!DEC$ ATTRIBUTES REFERENCE :: I_MACH,I_SLOT
	
	INCLUDE 'COMON4.ins'
	IMPLICIT none

	! Declaration
	! ----------- 
	INTEGER I_MACH, I_SLOT, IB, I_VAR, I_STATE, I_ICON		
	INTEGER IDXCONVERTER, DCNTRLTYPE, QCNTRLTYPE, ILIMITPRIORITY
	INTEGER issimstop
	INTEGER NDCBUS, NDCBUS_PREVIOUS, NDCLINES_PREVIOUS
	  
	REAL ps_initial, qs_initial ! initial bus active and reactive powers  
	REAL ps, qs		 			! current bus active and reactive powers  
	REAL icD_out, icQ_out 		! D- and Q-axis current in the system reference (DQ)	 

	COMPLEX us_phasor ! bus voltage phasor
	COMPLEX ss_phasor ! apparent power 
	COMPLEX ic_phasor ! converter current phasor.

	CHARACTER(1) IDGRID
	
	! Parameter assignment
	! --------------------
	! Index
	IB=NUMTRM(I_MACH) ! Bus sequence number
	I_STATE=STRTIN(2,I_SLOT) ! initial STATE index
	I_VAR=STRTIN(3,I_SLOT)   ! initial VAR index
	I_ICON=STRTIN(4,I_SLOT)  ! initial ICON index

	! ICON	 
	IDXCONVERTER = ICON(I_ICON)     ! index of the converter
	DCNTRLTYPE = ICON(I_ICON+1)     ! Ps-control: 1 , Udc-control: 2, delta_s-control: 3 
	QCNTRLTYPE = ICON(I_ICON+2)     ! Qs-control: 1 , Us-control: 2
	ILIMITPRIORITY = ICON(I_ICON+3) ! Current limit: P-priority: 1, Q-priority: 2 and P-Q equal priority: 3 or any other integer
	NDCBUS = ICON(I_ICON+4)           ! number of converters of the DC grid
	IDGRID = CHRICN(I_ICON+5)       ! grid identifier
	NDCBUS_PREVIOUS = ICON(I_ICON+6)  ! accumulated number of previous DC grid converters
	NDCLINES_PREVIOUS = ICON(I_ICON+7)  ! accumulated number of previous DC grid converters
	
	
	IF (DCNTRLTYPE.EQ.3) THEN       ! if feeding a passive grid, reactive power control is a voltage control
		QCNTRLTYPE = 2
	END IF

	! VARs
	ps_initial=VAR(I_VAR) 
	qs_initial=VAR(I_VAR+1) 
	icD_out=VAR(I_VAR+2) 
	icQ_out=VAR(I_VAR+3) 
	ps=VAR(I_VAR+4) 
	qs=VAR(I_VAR+5) 
	issimstop=VAR(I_VAR+10)

	! Common variables
	us_phasor = VOLT(IB)	  

	SELECT CASE (MODE)
		CASE (1) 
			! Initialization
			! --------------
			! Eint = V + Zsorce*I, ISORCE = Eint/Zsorce -> I = ISORCE - V/ZSORCE
			! ISORCE: norton equivalent source current in pu of SBASE
			! I: current at generator terminal bus in pu SBASE
			! ZSORCE: in pu of MBASE
			ic_phasor = ISORCE(I_MACH) - us_phasor/(ZSORCE(I_MACH)*SBASE/MBASE(I_MACH)) ! 
						 
			ss_phasor = us_phasor*CONJG(ic_phasor) ! in pu SBASE     
			ps_initial = REAL(ss_phasor)
			qs_initial = AIMAG(ss_phasor)
			 
			issimstop = 0    

		CASE (3) 
			! Current injections
			! ------------------         
			IF ((issimstop .GT. 0)) THEN
				!WRITE (LPDEV,*) 'TVSCON::CASE3: Disconnect SVSCON at bus ', NUMTRM(IB), '.'
				icD_out = 0.0               
				icQ_out = 0.0
				ps = 0.0               
				qs = 0.0
			END IF
					  
			ic_phasor = CMPLX(icD_out,icQ_out)
			ss_phasor = us_phasor*CONJG(ic_phasor) ! in pu SBASE     
			!ps = REAL(ss_phasor)
			!qs = AIMAG(ss_phasor)

			! Add SVSCON currents       
			ISORCE(I_MACH) = ic_phasor + us_phasor/(ZSORCE(I_MACH)*SBASE/MBASE(I_MACH))         

		CASE DEFAULT

	END SELECT

	! Re-assign algebraic variables
	! -----------------------------
	! VARs
	VAR(I_VAR) = ps_initial 
	VAR(I_VAR+1) = qs_initial 
	VAR(I_VAR+2) = icD_out 
	VAR(I_VAR+3) = icQ_out 
	VAR(I_VAR+4) = ps 
	VAR(I_VAR+5) = qs 
	VAR(I_VAR+10) = issimstop
	 
END
      
      
! ==========================================================
! SOLUTION OF DIFFERENTIAL EQUATIONS
! ==========================================================      
SUBROUTINE SVSCON(I_MACH,I_SLOT)
      
	!DEC$ ATTRIBUTES DLLEXPORT, DECORATE, ALIAS: "SVSCON":: SVSCON
	!DEC$ ATTRIBUTES REFERENCE :: I_MACH,I_SLOT

	INCLUDE 'COMON4.INS'
	USE MOD_READHYADCSIM
	USE MOD_MISC
	USE MOD_PROTECTION

	IMPLICIT none
	
	! Declaration
	! -----------
	INTEGER IB, I_SLOT, I_MACH, I_VAR, I_CON, I_ICON, I_STATE         ! PSS/e indices
	INTEGER I_STATE_DCGRID, I_STATE_DCGRID_G, aux_var_GLOBAL      ! local and global varible array index of the DCGRID STATE
	INTEGER IDXCONVERTER, DCNTRLTYPE, QCNTRLTYPE, ILIMITPRIORITY 
	INTEGER issimstop                                             ! disconnects unit and sets current injections to 0
	INTEGER dcontroltype_aux, qcontroltype_aux     
	INTEGER NDCBUS 												  ! # of buses of the dc-grid
	INTEGER NDCBUS_PREVIOUS, NDCLINES_PREVIOUS
	
	INTEGER, ALLOCATABLE :: CONVERTER_ACDC_BUS(:,:) ! ac and dc buses of each converter
	INTEGER IERR
	REAL PI
	PARAMETER (PI=3.14159265358979)
	REAL DELTAT
		
	REAL us, delta_s, us_initial, delta_s_ini
	REAL ps_initial, qs_initial 
	REAL icd, icq ! Currents in dq.
	REAL icd_initial, icq_initial 
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
	REAL d_udc 			! DC voltage derivative (comes from DCGRID model)

	REAL aux_ilimit, awu_d, awu_q        

	REAL ic_abs, ic_abs_puconv ! Current magnitude

	REAL ps_ref, qs_ref, us_ref ! references
	REAL icd_ref, icq_ref
	REAL udc_ref
	REAL ps_ref_aux, udc_ref_aux, qs_ref_aux, us_ref_aux

	REAL idc, idc_ini, udc, udc_ini
	REAL pdc, icc_grid

	REAL ps_0, addps_ref 
	REAL qs_0, addqs_ref       
	REAL ec, delta_c      
	REAL pvsc, qvsc

	COMPLEX us_phasor ! bus voltage.
	COMPLEX reference_transf ! convert inverter reference (dq) to system reference (DQ).
	COMPLEX ss_phasor ! Apparent power consumed by the load.
	COMPLEX ic_phasor ! Currents in DQ.
	COMPLEX sc_phasor
	COMPLEX ec_phasor

	CHARACTER(1) IDGRID
     

	REAL, ALLOCATABLE :: Udc_ini_vector(:,:)

	REAL TAU	
	REAL KD_P1, KD_P2, KD_I1, KD_I2, KD_D2, KQ_P1, KQ_I1, KQ_P2, KQ_I2
	REAL ICMAX_PUconv
	REAL PS_MAX_MW, PS_MIN_MW, QS_MAX_Mvar, QS_MIN_Mvar
	REAL ALOSS_MW, BLOSS_kV, CLOSS_RECT_Ohm, CLOSS_INV_Ohm
	REAL UDC_MAX, UDC_MIN, TUDCMAX, TUDCMIN

	REAL ic_max
	REAL ps_max, ps_min, qs_max, qs_min
	REAL aloss, bloss, c_rect, c_inv
	REAL CDC_uF, Cdc
	REAL UDC_NOMINAL_kV
	REAL ZBASE, ZDCBASE, SDCBASE_MVA
	COMPLEX zc
	REAL ploss
	
	INTEGER :: counter_uvdc, counter_ovdc
	REAL :: tinitial_uvdc, tinitial_ovdc

	REAL m_modulation, m_MODULATION_MAX

	REAL icd_ref_p, icq_ref_p

	

	COMMON /svsconinternal/ I_STATE_DCGRID_G, aux_var_GLOBAL ! This allows saving the value of I_STATE_DCGRID_G between calls. Note that I_STATE_DCGRID_G takes the value of the first DCGRID model read.

	! Parameter assignment
	! --------------------	  
	! Indexes
	IB=NUMTRM(I_MACH)
	I_CON = STRTIN(1,I_SLOT)
	I_STATE = STRTIN(2,I_SLOT) ! Model Calling Sequence Rules
	I_VAR = STRTIN(3,I_SLOT)
	I_ICON = STRTIN(4,I_SLOT)
		
	! CONs    
	TAU = CON(I_CON) ! Inverter time constant (e.g., 0.1)
	KD_P1 = CON(I_CON+1) ! Ps control
	KD_I1 = CON(I_CON+2)
	KD_P2 = CON(I_CON+3) ! Udc control 
	KD_I2 = CON(I_CON+4)
	KD_D2 = CON(I_CON+5) ! gain of the differential control of Udc -> for the DC-voltage control, PID works much better than PI
	KQ_P1 = CON(I_CON+6) ! Qs control
	KQ_I1 = CON(I_CON+7)
	KQ_P2 = CON(I_CON+8) ! Us control
	KQ_I2 = CON(I_CON+9)
	ICMAX_PUconv = CON(I_CON+10) ! Maximum inverter current/susceptance in pu with respect to inverter rating (e.g., 1.1)
	PS_MAX_MW = CON(I_CON+11)
	PS_MIN_MW = CON(I_CON+12)
	QS_MAX_Mvar = CON(I_CON+13)
	QS_MIN_Mvar = CON(I_CON+14)
	UDC_MAX = CON(I_CON+15)
	UDC_MIN = CON(I_CON+16)
	m_MODULATION_MAX = CON(I_CON+17)
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
	ps_initial = VAR(I_VAR) ! Reference active power, coincides with load flow values.
	qs_initial = VAR(I_VAR+1) 
	icD_out = VAR(I_VAR+2) 
	icQ_out = VAR(I_VAR+3) 
	ps = VAR(I_VAR+4) ! Active model output power. 
	qs = VAR(I_VAR+5) 
	icd_initial = VAR(I_VAR+6)
	icq_initial = VAR(I_VAR+7)
	icd = VAR(I_VAR+8)
	icq = VAR(I_VAR+9)
	issimstop = VAR(I_VAR+10)
	us_initial = VAR(I_VAR+11)
	ps_ref = VAR(I_VAR+12)
	udc_ref = VAR(I_VAR+13)
	qs_ref = VAR(I_VAR+14)
	us_ref = VAR(I_VAR+15)
	aux_ilimit = VAR(I_VAR+16)
	awu_d = VAR(I_VAR+17)
	awu_q = VAR(I_VAR+18)
	ic_abs_puconv = VAR(I_VAR+19)
	dcontroltype_aux = VAR(I_VAR+20)
	qcontroltype_aux = VAR(I_VAR+21)
	udc = VAR(I_VAR+22)
	idc = VAR(I_VAR+23)
	idc_ini = VAR(I_VAR+24)
	pdc = VAR(I_VAR+25)
	delta_s_ini = VAR(I_VAR+26)
	udc_ini = VAR(I_VAR+27)
	icd_ref = VAR(I_VAR+28)
	addps_ref = VAR(I_VAR+29)
	icq_ref = VAR(I_VAR+30)
	addqs_ref = VAR(I_VAR+31)
	m_modulation = VAR(I_VAR+32)
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
	
	ALLOCATE(Udc_ini_vector(NDCBUS,1)) ! allocate size
	ALLOCATE(CONVERTER_ACDC_BUS(NDCBUS,2))

	ic_max = ICMAX_PUconv*MBASE(I_MACH)/SBASE ! current limit in system base
	ps_max = PS_MAX_MW/SBASE               ! power limits in system base
	ps_min = PS_MIN_MW/SBASE
	qs_max = QS_MAX_Mvar/SBASE
	qs_min = QS_MIN_Mvar/SBASE

	ZBASE = BASVLT(IB)**2/SBASE               ! AC-side Zbase in ohms
	SDCBASE_MVA = SBASE
	ZDCBASE = 2*UDC_NOMINAL_kV**2/SDCBASE_MVA ! DC-side Zbase in ohms

	aloss = ALOSS_MW/SBASE                  ! losses coefs in system base
	bloss = BLOSS_kV/SQRT(3.0)/BASVLT(IB)
	c_rect = CLOSS_RECT_Ohm/ZBASE
	c_inv = CLOSS_INV_Ohm/ZBASE
	zc = ZSORCE(I_MACH)*SBASE/MBASE(I_MACH)       ! conexion impedance in system base
	ic_abs = ic_abs_puconv*MBASE(I_MACH)/SBASE ! system base      
	Cdc = (CDC_uF*1e-6)*ZDCBASE             ! [pu*s]

	I_STATE_DCGRID = I_STATE_DCGRID_G
	
	us_phasor = VOLT(IB)
	delta_s = ATAN2(AIMAG(us_phasor),REAL(us_phasor))			
	! DQ -> dq reference system transformation: xDQ = reference_transf*xdq
	IF (DCNTRLTYPE.EQ.3) THEN ! If feeding a passive grid, the VSC is a AC slack bus
		reference_transf = CMPLX(1.0, 0.0) 
	ELSE
		reference_transf = CMPLX(cos(delta_s),sin(delta_s))
	END IF
	
	SELECT CASE (MODE)

		CASE (1) 
			   
			! Inicialization 
			! ==============	  
			
			! get initial value of the DC voltage state from .txt files
			CALL SUB_READDCBUS(Udc_ini_vector, CONVERTER_ACDC_BUS, IDGRID, NDCBUS)
			! CALL SUB_READFROMFILESVSCON(Udc_ini_vector, CONVERTER_ACDC_BUS, IDGRID, NDCBUS)
			
			! get STATE of dc grid model DCGRID (governor-type model)	  
			CALL MDLIND(CONVERTER_ACDC_BUS(1,2), MACHID(I_MACH), 'GOV', 'STATE', I_STATE_DCGRID, IERR)
			IF (IERR.NE.0) THEN ! No governor-type model
				WRITE (LPDEV,*) 'SVSCON - CASE 1: No DCGRID model at bus ',NUMBUS(IB),' with id ',MACHID(I_MACH)
			END IF
			IF (aux_var_GLOBAL.NE.7) THEN
				I_STATE_DCGRID_G = I_STATE_DCGRID
				aux_var_GLOBAL = 7
			END IF
			
			! Algebraic variables
			! -------------------			
			
			us_initial = ABS(us_phasor)
			delta_s_ini = delta_s

			! AC-side converter current (in the dynamic simulation the C-filter is omitted)
			ps = ps_initial
			qs = qs_initial
			ss_phasor = CMPLX(ps_initial,qs_initial)
			ic_phasor = CONJG(ss_phasor/us_phasor)   
			icd_initial = REAL(ic_phasor/reference_transf) ! in dq axes
			icq_initial = AIMAG(ic_phasor/reference_transf) ! in dq axes
			icd = icd_initial
			icq = icq_initial
			
			! AC-side converter current out of range       
			ic_abs = ABS(CMPLX(icd_initial,icq_initial)) ! mag
			IF (ic_abs.GT.ic_max) THEN
				WRITE (LPDEV,*) 'SVSCON::CASE1: Current above limit.'
				aux_ilimit = 0.0 ! out of limit -> disable integrals
				awu_d = 0.0
				awu_q = 0.0
			ELSE
				aux_ilimit = 1.0
				awu_d = 1.0
				awu_q = 1.0
			END IF

			! AC-side converter voltages
			ec_phasor = us_phasor + zc*ic_phasor
			ec = ABS(ec_phasor)
			delta_c = ATAN2(AIMAG(ec_phasor), REAL(ec_phasor))

			! AC-side converter power
			sc_phasor = ec_phasor*CONJG(ic_phasor) ! interesa tenerlas, pero de momento no las usamos
			pvsc = REAL(sc_phasor)
			qvsc = AIMAG(sc_phasor)
			
			! DC-side voltage, current, power
			udc = Udc_ini_vector(IDXCONVERTER,1) ! extract initial dc voltage value
			udc_ini = udc
			! CALL SUB_GETPLOSS(ploss, ps_initial, ic_abs, aloss, bloss, c_inv, c_rect)
			CALL SUB_COMPUTEPLOSS(ploss, ps_initial, ic_abs, aloss, bloss, c_inv, c_rect)
			pdc = -(pvsc+ploss)
			idc_ini = pdc/udc_ini
			idc = idc_ini
			 		
			m_modulation = ec/MAX(udc_ini,0.0001)

			! references
			ps_ref = ps_initial
			qs_ref = qs_initial
			us_ref = us_initial
			udc_ref = udc_ini 
			addps_ref = 0.0 ! supplementary controller input (SPWDRD)
			addqs_ref = 0.0 ! supplementary controller input (SQWDRD)
			icd_ref = icd_initial
			icq_ref = icq_initial
			
			issimstop = 0
			
			! State variables
			! ---------------
			! initial values of states
			xd = icd_initial 
			xq = icq_initial
			md = 0.0
			mq = 0.0
			ndc = 0.0
			nq = 0.0
			eta_d = 0.0

			ps_0 = 0.0
			qs_0 = 0.0
			icc_grid = idc
			d_udc = 0.0 
			
			! Initialize PSS/E arrays and print 
			! ---------------------------------
			PMECH(I_MACH) = -pdc*SBASE/MBASE(I_MACH) ! pu machine rating
			EFD(I_MACH) = ec
			SPEED(I_MACH) = BSFREQ(IB)
			ANGLE(I_MACH) = delta_c               ! in (deg), not given in rads! 
			ETERM(I_MACH) = us_initial
			PELEC(I_MACH) = ps_initial            ! pu system rating
			QELEC(I_MACH) = qs_initial            ! pu system rating
			
			WRITE (LPDEV,*) 'SVSCON - CASE 1: Converter ',IDXCONVERTER,' at bus ',NUMBUS(IB),' with id ',MACHID(I_MACH),' initialized.'
			!WRITE (LPDEV,*) 'SVSCON - CASE 1: us = ',us_initial,' pu, ec = ',ec,' pu'
			!WRITE (LPDEV,*) 'SVSCON - CASE 1: ps = ',ps,' pu, qs = ', qs, ' pu'
			!WRITE (LPDEV,*) 'SVSCON - CASE 1: pc = ',pvsc,', pu ploss = ',ploss,' pu'
			!WRITE (LPDEV,*) 'SVSCON - CASE 1: udc = ',udc_ini,' pu, pdc = ',-pdc,' pu'
			!WRITE (LPDEV,*) 'SVSCON - CASE 1: idc = ',idc, ' pu'
			!WRITE (LPDEV,*) 'SVSCON - CASE 1: icd = ',icd_initial,' pu, icq = ',icq_initial, ' pu'

								
	CASE (2) 
	   
		! Compute derivatives
		! ===================
		us = ABS(us_phasor) 		

		! get DC voltage of the DC grid model (if any, else from VARs)
		IF (I_STATE_DCGRID.GT.0) THEN
			udc = STATE(I_STATE_DCGRID + IDXCONVERTER + NDCLINES_PREVIOUS + NDCBUS_PREVIOUS - 1)
		END IF

		ps_ref_aux = ps_ref
		qs_ref_aux = qs_ref
		udc_ref_aux = udc_ref
		us_ref_aux = us_ref		
		
		icd_ref_p = icd_ref
		icq_ref_p = icq_ref
								   
		! Compute derivatives
		IF (TAU.LT.(2*DELTAT)) THEN								! current control
			d_xd = 0.0
			d_xq = 0.0
			xd = icd_ref_p
			xq = icq_ref_p
		ELSE
			d_xd = (icd_ref_p - xd)/TAU
			d_xq = (icq_ref_p - xq)/TAU
		END IF

		IF(dcontroltype_aux.EQ.3) THEN 							! delta_s control
			d_md = 0.0
			d_ndc = 0.0
			d_eta_d = (1000.0)*(delta_s_ini - delta_s)*awu_d
		ELSE IF(dcontroltype_aux.EQ.2) THEN						! DC-voltage control
			d_md = 0.0
			d_ndc = KD_I2*(udc_ref_aux - udc)*awu_d
			d_eta_d = 0.0
		ELSE 													! Ps control
			d_md = KD_I1*(ps_ref_aux - ps)*awu_d 
			d_ndc = 0.0
			d_eta_d = 0.0
		END IF

		IF(qcontroltype_aux.EQ.2) THEN							! Us control
			d_mq = 0.0
			d_nq = KQ_I2*(us_ref_aux - us)*awu_q 
		ELSE 													! Qs control
			d_mq = KQ_I1*(qs_ref_aux - qs)*awu_q 
			d_nq = 0.0
		END IF
		   
	CASE (3) 

		! Compute output
		! ==============
		us = ABS(us_phasor) 

		! get DC voltage into the DC grid (if any, else from VARs)
		IF (I_STATE_DCGRID.GT.0) THEN
			udc = STATE(I_STATE_DCGRID + IDXCONVERTER + NDCLINES_PREVIOUS + NDCBUS_PREVIOUS - 1)
		END IF

		! update and compute references	
		ps_ref = ps_initial + addps_ref ! P: operation reference + additional supplementary reference (for ex: droop)
		qs_ref = qs_initial + addqs_ref ! Q: operation reference + additional supplementary reference (for ex: droop)
		
		dcontroltype_aux = DCNTRLTYPE
		qcontroltype_aux = QCNTRLTYPE		
		CALL SUB_SETIREF(dcontroltype_aux, qcontroltype_aux, ps_ref, qs_ref, udc_ref, us_ref, icd_ref, icq_ref, & 
			md, mq, ndc, nq, eta_d, us, delta_s, ps, qs, udc, icd_initial, icq_initial, delta_s_ini, &
			ps_max, ps_min, qs_max, qs_min, UDC_MAX, UDC_MIN, KD_P1, KD_P2, KD_D2, KQ_P1, KQ_P2)
		
		! limit current references and set anti-wind up indicators
		m_modulation = ec/MAX(udc,0.0001)
		! CALL SUB_LIMITIREF(icd_ref, icq_ref, aux_ilimit, awu_d, awu_q, us_phasor, udc, ILIMITPRIORITY, ic_max, m_MODULATION_MAX, zc)		
		CALL SUB_IMAXLIMITSIREF(icd_ref, icq_ref, awu_d, awu_q, ILIMITPRIORITY, ic_max)
		CALL SUB_ECMAXLIMITSIREF(icd_ref, icq_ref, us, udc, m_MODULATION_MAX, zc)

		! DC protection
		CALL SUB_OVDCPROT(issimstop, counter_ovdc, tinitial_ovdc, udc, TIME, UDC_MAX, TUDCMAX, NUMBUS(IB),'SVSCON',LPDEV)
		CALL SUB_UVDCPROT(issimstop, counter_uvdc, tinitial_uvdc, udc, TIME, UDC_MIN, TUDCMIN, NUMBUS(IB),'SVSCON',LPDEV)
		
		icd = 0.0
		icq = 0.0
		IF (issimstop.LT.1) THEN		  
			icd = xd
			icq = xq
		END IF
		
		! AC-side converter current, voltage, power 
		ps = us*icd ! same as ps = REAL(ss_vec)
		qs = -us*icq ! qs = AIMAG(ss_vec)	
		ic_phasor = reference_transf*CMPLX(icd,icq) ! dq -> DQ, pu in system rating	
		ec_phasor = us_phasor + zc*ic_phasor ! converter voltage -> output		
		sc_phasor = ec_phasor*CONJG(ic_phasor) ! interesa tenerlas, pero de momento no las usamos
		ic_abs = ABS(ic_phasor) ! mag
		icD_out = REAL(ic_phasor) 
		icQ_out = AIMAG(ic_phasor)
		ec = ABS(ec_phasor)
		delta_c = ATAN2(AIMAG(ec_phasor), REAL(ec_phasor))
		pvsc = REAL(sc_phasor)
		qvsc = AIMAG(sc_phasor)

		! DC-side converter current current and power
		! CALL SUB_GETPLOSS(ploss, ps, ic_abs, aloss, bloss, c_inv, c_rect)	
		CALL SUB_COMPUTEPLOSS(ploss, ps, ic_abs, aloss, bloss, c_inv, c_rect)
		pdc = -(pvsc+ploss) ! the same... only necessary one of them output of the converter model; input of the DC-grid model
		idc = pdc/udc ! output of the converter model; input of the DC-grid model
		
		SPEED(I_MACH) = BSFREQ(IB)
		ANGLE(I_MACH) = delta_c ! in (deg), not given in rads!!
		ETERM(I_MACH) = us
		PELEC(I_MACH) = ps ! in pu system rating
		QELEC(I_MACH) = qs ! in pu system rating
		PMECH(I_MACH) = -pdc*SBASE/MBASE(I_MACH) ! in pu machine rating
		EFD(I_MACH) = ec
				  
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
	VAR(I_VAR) = ps_initial  
	VAR(I_VAR+1) = qs_initial 
	VAR(I_VAR+2) = icD_out 
	VAR(I_VAR+3) = icQ_out 
	VAR(I_VAR+4) = ps 
	VAR(I_VAR+5) = qs 
	VAR(I_VAR+6) = icd_initial
	VAR(I_VAR+7) = icq_initial
	VAR(I_VAR+8) = icd
	VAR(I_VAR+9) = icq
	VAR(I_VAR+10) = issimstop
	VAR(I_VAR+11) = us_initial
	VAR(I_VAR+12) = ps_ref
	VAR(I_VAR+13) = udc_ref
	VAR(I_VAR+14) = qs_ref
	VAR(I_VAR+15) = us_ref
	VAR(I_VAR+16) = aux_ilimit
	VAR(I_VAR+17) = awu_d
	VAR(I_VAR+18) = awu_q
	VAR(I_VAR+19) = ic_abs_puconv
	VAR(I_VAR+20) = dcontroltype_aux
	VAR(I_VAR+21) = qcontroltype_aux
	VAR(I_VAR+22) = udc
	VAR(I_VAR+23) = idc
	VAR(I_VAR+24) = idc_ini
	VAR(I_VAR+25) = pdc
	VAR(I_VAR+26) = delta_s_ini
	VAR(I_VAR+27) = udc_ini
	VAR(I_VAR+28) = icd_ref 
	VAR(I_VAR+29) = addps_ref 
	VAR(I_VAR+30) = icq_ref
	VAR(I_VAR+31) = addqs_ref
	VAR(I_VAR+32) = m_modulation
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
      
END 

SUBROUTINE SUB_READFROMFILESVSCON(Udc_ini_vector, CONVERTER_ACDC_BUS, IDGRID, NDCBUS)
	
	CHARACTER(1) :: IDGRID, line_header 
	INTEGER :: NDCBUS, CONVERTER_ACDC_BUS(NDCBUS,1), ivsc ! ac and dc buses of each converter
	REAL :: Udc_ini_vector(NDCBUS,1)
	
	OPEN(UNIT=20, FILE='.\data_acdcbus.txt') ! .txt files where .dll is located
	OPEN(UNIT=21, FILE='.\data_Udc_ini.txt')

	line_header = "+"
	DO WHILE (line_header.NE.IDGRID) ! read until header according to grid identifier IDGRID
		READ(20, *) line_header
	END DO

	line_header = "+"
	DO WHILE (line_header.NE.IDGRID) ! read until header according to grid identifier IDGRID
		READ(21, *) line_header
	END DO

	DO ivsc=1,NDCBUS
		READ(20,*)CONVERTER_ACDC_BUS(ivsc,1),CONVERTER_ACDC_BUS(ivsc,2)
		READ(21,*)Udc_ini_vector(ivsc,1)
	END DO
	CLOSE(20)
	CLOSE(21)

END 

SUBROUTINE SUB_GETPLOSS(ploss, ps, ic_abs, aloss, bloss, c_inv, c_rect)
	
	REAL :: ploss, ps, ic_abs, aloss, bloss, c_inv, c_rect

	IF (ps.GE.0.0) THEN 
		ploss = aloss + bloss*ic_abs + c_inv*ic_abs**2  ! inverter
	ELSE 
		ploss = aloss + bloss*ic_abs + c_rect*ic_abs**2 ! rectifier
	END IF

END

SUBROUTINE SUB_SETIREF(dcontroltype_aux, qcontroltype_aux, ps_ref, qs_ref, udc_ref, us_ref, icd_ref, icq_ref, & 
	md, mq, ndc, nq, eta_d, us, delta_s, ps, qs, udc, icd_initial, icq_initial, delta_s_ini, &
	ps_max, ps_min, qs_max, qs_min, UDC_MAX, UDC_MIN, KD_P1, KD_P2, KD_D2, KQ_P1, KQ_P2)

	INTEGER :: dcontroltype_aux, qcontroltype_aux
	REAL :: md, mq, ndc, nq, eta_d
	REAL :: us, delta_s, ps, qs, udc
	REAL :: ps_ref, qs_ref, udc_ref, idc_ref, us_ref, icd_ref, icq_ref, ps_ref_aux, udc_ref_aux, qs_ref_aux, us_ref_aux
	REAL :: icd_initial, icq_initial, delta_s_ini
	REAL :: ps_max, ps_min, qs_max, qs_min, UDC_MAX, UDC_MIN, KD_P1, KD_P2, KD_D2, KQ_P1, KQ_P2

	! limit total references (power and voltage)
	ps_ref_aux = ps_ref
	udc_ref_aux = udc_ref
	qs_ref_aux = qs_ref
	us_ref_aux = us_ref
	
	ps_ref_aux = MIN(MAX(ps_ref,ps_min),ps_max)
	qs_ref_aux = MIN(MAX(qs_ref,qs_min),qs_max)
	
	!IF ((udc_ref.GT.UDC_MAX).OR.udc_ref.LT.UDC_MIN) THEN	! DC voltage reference exceeding voltage limits. This could happen if DC-voltage control VSC is lost and the other VSC must take over, being so far under control of Ps (where udcref = udc). However, in MTDC the possibility for multiple PI DC voltage controller coud appear and by now aprubt current set point changes are resulting
	!	dcontroltype_aux = 2 								! control of Udc (DCNTRLTYPE = 2)
	!	udc_ref_aux = MIN(MAX(udc_ref,UDC_MIN),UDC_MAX)
	!END IF

	! current controller references: icd and icq
	IF(qcontroltype_aux.EQ.2) THEN 											! control of Us
		icq_ref = icq_initial - KQ_P2*(us_ref_aux - us) - nq     			! -(icq0 + Kq2*(usref  - us) + nq), icq = qs0/usref, usref = us0
		!qs_ref_aux = -us_ref_aux*icq_ref 						    		! reactive power needed for Us control		  
		!IF ((qs_ref_aux.GT.qs_max).OR.(qs_ref_aux.LT.qs_min)) THEN
		!	qcontroltype_aux = 1 											! control of Qs and set reactive power set point
		!	qs_ref_aux = MIN(MAX(qs_ref_aux,qs_min),qs_max)
		!	icq_ref = -qs_ref_aux/us - KQ_P1*(qs_ref_aux - qs) - mq
		!END IF		  
	ELSE 																	! control of Qs
		icq_ref = -qs_ref_aux/us - KQ_P1*(qs_ref_aux - qs) - mq 			! -(qsref/us + Kq1*(qsref - qs) + mq)		  
		!us_ref_aux = -qs_ref_aux/icq_ref ! for some calculus  
	END IF

	IF(dcontroltype_aux.EQ.3) THEN 									     	! Passive grid
		icd_ref = icd_initial + (100.0)*(delta_s_ini - delta_s) + eta_d      	! keep delta_s to 0 -> slack bus behavior	  
		!udc_ref_aux = udc
	ELSE IF (dcontroltype_aux.EQ.2) THEN							      	! DC voltage control
		!idc_ref = icc_grid - KD_D2*d_udc + KD_P2*(udc_ref_aux - udc) + ndc ! PID DC voltage control 
		!icd_ref = (-1/ec_d)*(ec_q*icq + udc*idc_ref + ploss) 				! comes from ' pvsc + pdc + ploss = 0 '
		icd_ref = icd_initial - KD_P2*(udc_ref_aux - udc) - ndc  				! icd0 - (Kd2*(udcref - udc) + ndc), icd0 = psref/usref
	ELSE 																    ! AC active power control
		icd_ref = ps_ref_aux/us + KD_P1*(ps_ref_aux - ps) + md 				! psref/us + (Kd1*(psref - ps) + mdc)		  
		!udc_ref_aux = udc												 	
	END IF
	
	ps_ref  = ps_ref_aux
	qs_ref = qs_ref_aux
	us_ref = us_ref_aux
	udc_ref = udc_ref_aux

END
      
SUBROUTINE SUB_LIMITIREF(icd_ref, icq_ref, aux_ilimit, awu_d, awu_q, us_phasor, udc, ILIMITPRIORITY, ic_max, m_MODULATION_MAX, zc)

	INTEGER :: ILIMITPRIORITY
	REAL :: icd_ref, icq_ref, aux_ilimit, awu_d, awu_q
	REAL :: icd_ref_aux, icq_ref_aux, ic_abs_ref, coseno_alpha_ref, seno_alpha_ref, ec_max, ec_ref, delta_c_ref, icd_ref_p, icq_ref_p
	COMPLEX :: us_phasor, ic_ref_phasor, ec_ref_phasor, ec_ref_phasor_p, ic_ref_phasor_p
	COMPLEX :: zc
	REAL :: udc, ic_max, m_MODULATION_MAX

	! Check current limits
	icd_ref_aux = icd_ref
	icq_ref_aux = icq_ref
	ic_abs_ref = SQRT(icd_ref**2 + icq_ref**2)
	coseno_alpha_ref = icd_ref/ic_abs_ref
	seno_alpha_ref = icq_ref/ic_abs_ref
	IF (ic_abs_ref.GT.ic_max) THEN
		aux_ilimit = 0.0 ! out of limit -> disable integrals
		IF (ILIMITPRIORITY.EQ.1) THEN ! P-priority
		  
			awu_q = 0.0 ! icq will be always on the limit in this case
			IF (ABS(icd_ref).GT.ic_max) THEN
				awu_d = 0.0 ! icd is on the limit
			ELSE
				awu_d = 1.0 ! icd<=ic_max
			END IF
			icd_ref_aux = MIN(ABS(icd_ref), ic_max)*SIGN(1.0,icd_ref) ! change the current reference: icq_ref = icq_ref_aux
			icq_ref_aux = SQRT(ic_max**2 - icd_ref_aux**2)*SIGN(1.0,icq_ref)
		  
		ELSE IF (ILIMITPRIORITY.EQ.2) THEN ! Q-priority
		  
			awu_d = 0.0 ! icq will be always on the limit in this case
			IF (ABS(icq_ref).GT.ic_max) THEN
				awu_q = 0.0 ! icq is on the limit
			ELSE
				awu_q = 1.0 ! icq<=ic_max
			END IF
			icq_ref_aux = MIN(ABS(icq_ref), ic_max)*SIGN(1.0,icq_ref) ! change the current reference: icq_ref = icq_ref_aux
			icd_ref_aux = SQRT(ic_max**2 - icq_ref_aux**2)*SIGN(1.0,icd_ref)
		  
		ELSE ! P-Q equal priority
		  
			awu_d = 0.0
			awu_q = 0.0
			icd_ref_aux = ic_max*coseno_alpha_ref
			icq_ref_aux = ic_max*seno_alpha_ref
		  
		END IF
	  
	ELSE
		aux_ilimit = 1.0
		awu_d = 1.0
		awu_q = 1.0
	END IF


	! Check Maximum modulation index (ec <= m_MODULATION_MAX*udc)
	ec_max = m_MODULATION_MAX*udc
	ic_ref_phasor = CMPLX(icd_ref_aux, icq_ref_aux) 
	ec_ref_phasor = us_phasor + zc*ic_ref_phasor
	ec_ref = ABS(ec_ref_phasor)
	delta_c_ref = ATAN2(AIMAG(ec_ref_phasor), REAL(ec_ref_phasor))
	

	IF (ec_ref.GT.ec_max) THEN	  
		ec_ref_phasor_p = CMPLX(ec_max*cos(delta_c_ref), ec_max*sin(delta_c_ref)) 	  
		ic_ref_phasor_p = (ec_ref_phasor_p - us_phasor)/zc
		icd_ref_p = REAL(ic_ref_phasor_p)
		icq_ref_p = AIMAG(ic_ref_phasor_p)
	ELSE
		icd_ref_p = icd_ref_aux
		icq_ref_p = icq_ref_aux  
	END IF

	icd_ref = icd_ref_p
	icq_ref = icq_ref_p
	
END
