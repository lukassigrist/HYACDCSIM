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
SUBROUTINE TSCGFO(I_MACH,I_SLOT)
      
	!DEC$ ATTRIBUTES DLLEXPORT, DECORATE, ALIAS: "TVSCON":: TVSCON
	!DEC$ ATTRIBUTES REFERENCE :: I_MACH,I_SLOT
	
	INCLUDE 'COMON4.ins'
	IMPLICIT none

	! Declaration
	! ----------- 
	INTEGER I_MACH, I_SLOT, IB, I_VAR, I_STATE, I_ICON		
	INTEGER issimstop
	   
	REAL ucR, ucI 

	COMPLEX us_phasor_RI, uc_phasor_RI, it_phasor_RI, sc_phasor, ss_phasor
	COMPLEX zc

	
	! Parameter assignment
	! --------------------
	! Index
	IB = NUMTRM(I_MACH) 		! Bus sequence number
	I_STATE = STRTIN(2,I_SLOT)	! initial STATE index
	I_VAR = STRTIN(3,I_SLOT)   	! initial VAR index

	! VARs
	ucR = VAR(I_VAR+5) 
	ucI = VAR(I_VAR+6) 
	issimstop = VAR(I_VAR+10)

	! Common variables
	us_phasor_RI = VOLT(IB)	
	zc = ZSORCE(I_MACH)*SBASE/MBASE(I_MACH)

	SELECT CASE (MODE)
		CASE (1) 
			! Initialization
			! --------------
			! Eint = V + Zsorce*I, ISORCE = Eint/Zsorce -> I = ISORCE - V/ZSORCE
			! ISORCE: norton equivalent source current in pu of SBASE
			! I: current at generator terminal bus in pu SBASE
			! ZSORCE: in pu of MBASE
			it_phasor_RI = ISORCE(I_MACH) - us_phasor_RI/zc 			
			uc_phasor_RI = ISORCE(I_MACH)*zc
						 
			ss_phasor = us_phasor_RI*CONJG(it_phasor_RI) ! in pu SBASE     
			sc_phasor = us_phasor_RI*CONJG(ISORCE(I_MACH))
			 
			issimstop = 0    

		CASE (3) 
			! Current injections
			! ------------------         
			IF ((issimstop .GT. 0)) THEN
				ucR = 0.0               
				ucI = 0.0
				zc = 10e6
				ZSORCE(I_MACH) = zc
			END IF
				
			! Add SVSCON currents       
			ISORCE(I_MACH) = CMPLX(ucR,ucI)/zc 
				
			it_phasor_RI = ISORCE(I_MACH) - us_phasor_RI/zc 					 
			sc_phasor = us_phasor_RI*CONJG(ISORCE(I_MACH))
		CASE DEFAULT

	END SELECT

	! Re-assign algebraic variables
	! -----------------------------
	PELEC(IMACH)=REAL(sc_phasor) 	! in pu system rating
    QELEC(IMACH)=AIMAG(sc_phasor) 	! in pu system rating
    ETERM(IMACH)=ABS(us_phasor_RI)
	EFD(IMACH)=ABS(uc_ph
	
	! VARs
	VAR(I_VAR+5) = ucR
	VAR(I_VAR+6) = ucQ
	VAR(I_VAR+10) = issimstop
	 
END
      
      
! ==========================================================
! SOLUTION OF DIFFERENTIAL EQUATIONS
! ==========================================================      
SUBROUTINE VSCGFO(I_MACH,I_SLOT)
      
	!DEC$ ATTRIBUTES DLLEXPORT, DECORATE, ALIAS: "SVSCON":: SVSCON
	!DEC$ ATTRIBUTES REFERENCE :: I_MACH,I_SLOT

	USE MOD_TFBLOCKS 		! use module of generic transfer functions
	INCLUDE 'COMON4.INS'
	IMPLICIT none
	
	! Declaration
	! -----------
	INTEGER IB, I_SLOT, I_MACH, I_VAR, I_CON, I_ICON, I_STATE       ! PSS/e indices
	INTEGER I_STATE_DCGRID, I_STATE_DCGRID_G, aux_var_GLOBAL      	! local and global varible array index of the DCGRID STATE
	INTEGER IDXCONVERTER, DCNTRLTYPE, QCNTRLTYPE, ILIMITPRIORITY 
	INTEGER issimstop                                             	! disconnects unit and sets current injections to 0
	INTEGER dcontroltype_aux, qcontroltype_aux     
	INTEGER NDCBUS 												  	! # of buses of the dc-grid
	INTEGER NDCBUS_PREVIOUS, NDCLINES_PREVIOUS
	
	INTEGER, ALLOCATABLE :: CONVERTER_ACDC_BUS(:,:) 				! ac and dc buses of each converter
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


	CALL BASE_FREQUENCY_API(FBASE, 1, IERR)
	WBASE = 2*PI*FBASE							! Check if OMEGAf exists	
	ZBASE = BASVLT(IB)**2/SBASE               	! AC-side Zbase in ohms
	SDCBASE_MVA = SBASE
	ZDCBASE = 2*UDC_NOMINAL_kV**2/SDCBASE_MVA 	! DC-side Zbase in ohms

	aloss = ALOSS_MW/SBASE                  ! losses coefs in system base
	bloss = BLOSS_kV/SQRT(3.0)/BASVLT(IB)
	c_rect = CLOSS_RECT_Ohm/ZBASE
	c_inv = CLOSS_INV_Ohm/ZBASE
	zc = ZSORCE(I_MACH)*SBASE/MBASE(I_MACH)       ! conexion impedance in system base
	ic_abs = ic_abs_puconv*MBASE(I_MACH)/SBASE ! system base      
	Cdc = (CDC_uF*1e-6)*ZDCBASE             ! [pu*s]

	I_STATE_DCGRID = I_STATE_DCGRID_G
	
	us_phasor_RI = VOLT(IB)
	delta_s = ATAN2(AIMAG(us_phasor_RI),REAL(us_phasor_RI))			
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
			CALL SUB_READFROMFILESVSCON(Udc_ini_vector, CONVERTER_ACDC_BUS, IDGRID, NDCBUS)
			
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
				
			! DC-side voltage, current, power
			udc = Udc_ini_vector(IDXCONVERTER,1) 			! extract initial dc voltage value
			udc_ini = udc
			
			! AC-side voltage, current, power
			ps_initial = ps
			qs_initial = qs
			ss_phasor = CMPLX(ps_initial,qs_initial)		
			it_phasor_RI = CONJG(ss_phasor/us_phasor_RI) 
			uc_phasor_RI = us_phasor_RI + zc*it_phasor_RI
			sc_phasor = uc_phasor_RI*CONJ(it_phasor_RI)
			pc_initial = REAL(sc_phasor)
			qc_initial = AIMAG(sc_phasor)
			
			it_phasor_dq = it_phasor_RI/reference_transf
			itd = REAL(it_phasor_dq)
			itq = AIMAG(it_phasor_dq)
			it_abs = ABS(it_phasor_dq)
			
			CALL SUB_GETPLOSS(ploss, pc_initial, it_abs, aloss, bloss, c_inv, c_rect)
			pdc = -(pc_initial+ploss)
			idc_ini = pdc/udc_ini
			idc = idc_ini

			! references
			uc_ref = uc_initial
			udc_ref = udc_ini 
			qc_ref = qc_initial
			addps_ref = 0.0 		! supplementary controller input (SPWDRD)
			addqs_ref = 0.0 		! supplementary controller input (SQWDRD)
			
			issimstop = 0
			
			! State variables
			! ---------------
			ypfc0 = 0.0
			ypod0 = 0.0
			yq0 = 0.0
			wvsm0 = 0.0
			ydcpss0 = 0.0
			deltavsm0 = ATAN2(ucd, ucq)
								
			! DC voltage control
			CALL SUB_PI(-pdc,xdc,dxdtdc,udc^2-udc_ref^2,1,DELTAT,KD_P2,KD_I2,PMAX,PMIN)
			
			! virtual synchronous machine
			deltap = -pdc - (pc_initial + ploss) - Kpfc*ypfc0 - Kd*ypod0		
			CALL SUB_INTEGRATORWINDUP(wvsm0,xwvsm,dxdtwvsm,deltap,1,H,10e6,-10e6)		
			CALL SUB_WASHOUTWINDUP(ypod0,xpod,dxdtpod,wvsm0,1,DELTAT,TPOD)
			CALL SUB_FIRSTORDERWINDUP(ypfc0,xpfc,dxdtpfc,wvsm0,1,1,DELTAT,TPFC,PMAX/Kpfd,PMIN/Kpfd)
			CALL SUB_INTEGRATORWINDUP(deltavsm0,xdvsm,dxdtdvsm,wvsm0,1,WBASE,10e6,-10e6)
			
			! AC voltage control
			deltaq = qcref - qc_initial
			CALL SUB_FIRSTORDERWINDUP(yq0,xq,dxdtq,deltaq,1,1,DELTAT,TQ,QMAX/Kq,QMIN/Kq)
			
			! transient virtual resistance
			CALL SUB_WASHOUTWINDUP(itvid0,xitvid,dxdtitvid,itd,1,DELTAT,TTVR)						! d-axis TVR current
			CALL SUB_WASHOUTWINDUP(itviq0,xitviq,dxdtitviq,itq,1,DELTAT,TTVR)						! q-axis TVR current
			
			! UDC PSS
			CALL SUB_WASHOUTWINDUP(ydcpss0,xudcpss,dxdtudcpss,udc,1,DELTAT,TUDCPSS)					! UDC PSS
			
						
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
			uc_phasor_dq = CMPLX(ucd,ucq) 										! from VARs (in steady state with current limits: ucq = 0)
			it_phasor_dq = (us_phasor_RI/reference_transf - uc_phasor_dq)/zc 	! dq axes
			itd = REAL(it_phasor_dq)
			itq = AIMAG(it_phasor_dq)
			it_abs = ABS(it_phasor_dq)
			pc = ucd*itd + ucq*itq
			qc = -ucd*itq + ucq*itd 
			
			CALL SUB_GETPLOSS(ploss, pc, it_abs, aloss, bloss, c_inv, c_rect)
			
			CALL SUB_WASHOUTWINDUP(ypod,xpod,dxdtpod,wvsm,3,DELTAT,TPOD)
			CALL SUB_FIRSTORDERWINDUP(ypfc,xpfc,dxdtpfc,wvsm,3,0,DELTAT,TPFC,PMAX,PMIN)
			CALL SUB_INTEGRATORWINDUP(wvsm,xvsm,dxdtvsm,deltap,3,H,-9999,9999)
			CALL SUB_PI(ydc,xdc,dxdtdc,udc^2-udc_ref^2,3,DELTAT,KD_P2,KD_I2,PMAX,PMIN)
			
			! virtual synchronous machine
			deltap = ydc + addps_ref - Kpfc*ypfc - Kd*ypod - (pc + ploss)		
			CALL SUB_INTEGRATORWINDUP(wvsm,xwvsm,dxdtwvsm,deltap,2,H,-9999,9999)					! virtual rotor			
			CALL SUB_WASHOUTWINDUP(ypod,xpod,dxdtpod,wvsm,2,DELTAT,TPOD)							! virtual POD
			CALL SUB_FIRSTORDERWINDUP(ypfc,xpfc,dxdtpfc,wvsm,2,0,DELTAT,TPFC,PMAX/Kpfc,PMIN/Kpfc)	! virtual PFC			
			CALL SUB_INTEGRATORWINDUP(deltavsm,xdvsm,dxdtdvsm,wvsm,2,WBASE,-9999,9999)			! virtual angle
			
			! AC voltage control
			deltaq = qc_ref - qc
			CALL SUB_FIRSTORDERWINDUP(yq,xq,dxdtq,deltaq,2,0,DELTAT,TQ,QMAX/Kq,QMIN/Kq)
			
			! DC voltage control
			CALL SUB_PI(ydc,xdc,dxdtdc,udc^2-udc_ref^2,2,DELTAT,KD_P2,KD_I2,PMAX,PMIN)
			
			! transient virtual resistance
			CALL SUB_WASHOUTWINDUP(itvid,xitvid,dxdtitvid,itd,2,DELTAT,TTVR)						! d-axis TVR current
			CALL SUB_WASHOUTWINDUP(itviq,xitviq,dxdtitviq,itq,2,DELTAT,TTVR)						! q-axis TVR current
			
			! UDC PSS
			CALL SUB_WASHOUTWINDUP(ydcpss,xudcpss,dxdtudcpss,udc,2,DELTAT,TUDCPSS)					! UDC PSS
			   
		CASE (3) 
		
			! get DC voltage into the DC grid (if any, else from VARs)
			IF (I_STATE_DCGRID.GT.0) THEN
				udc = STATE(I_STATE_DCGRID + IDXCONVERTER + NDCLINES_PREVIOUS + NDCBUS_PREVIOUS - 1)
			END IF	
			
			! DC protection
			CALL SUB_OVDCPROT(issimstop, udc, TIME, tinitial_ovdc, counter_ovdc, UDC_MAX, TUDCMAX)
			CALL SUB_UVDCPROT(issimstop, udc, TIME, tinitial_uvdc, counter_uvdc, UDC_MIN, TUDCMIN)
						
			! virtual synchronous machine
			CALL SUB_WASHOUTWINDUP(ypod,xpod,dxdtpod,wvsm,3,DELTAT,TPOD)					! virtual POD
			CALL SUB_FIRSTORDERWINDUP(ypfc,xpfc,dxdtpfc,wvsm,3,0,DELTAT,TPFC,PMAX,PMIN)		! virtual PFC	
			CALL SUB_INTEGRATORWINDUP(wvsm,xvsm,dxdtvsm,deltap,3,H,-9999,9999)				! virtual rotor
			CALL SUB_INTEGRATORWINDUP(deltavsm,xdvsm,dxdtdvsm,wvsm,3,WBASE,-9999,9999)	! virtual rotor angle
			reference_transf = CMPLX(cos(deltavsm),sin(deltavsm))
			
			! DC voltage control
			CALL SUB_PI(ydc,xdc,dxdtdc,udc^2-udc_ref^2,3,DELTAT,KD_P2,KD_I2,PMAX,PMIN)
			
			! AC voltage control
			CALL SUB_FIRSTORDERWINDUP(yq,xq,dxdtq,deltaq,3,0,DELTAT,TQ,QMAX/Kq,QMIN/Kq)			
			ucd = Kq*yq + uc_ref
			ucq = 0.0	
			
			! virtual impedance
			CALL SUB_WASHOUTWINDUP(itvid,xitvid,dxdtitvid,itd,2,DELTAT,TTVR)						! d-axis TVR current
			CALL SUB_WASHOUTWINDUP(itviq,xitviq,dxdtitviq,itq,2,DELTAT,TTVR)						! q-axis TVR current
			ucd = ucd - RTVR*itvid
			ucq = ucq - RTVR*itviq
			
			! UDC PSS
			CALL SUB_WASHOUTWINDUP(ydcpss,xudcpss,dxdtudcpss,udc,2,DELTAT,TUDCPSS)					! UDC PSS
			ucd = ucd + Kdcpss*ydcpss
			
			uc_phasor_dq = CMPLX(ucd, ucq)
			
			! virtual impedance-based current limit
			CALL SUB_VIILIMIT(uc_phasor_dq, us_phasor_dq, ic_max, zc, KPVI, SIGMAXR)
			
			! current limit
			us_phasor_dq = us_phasor_RI/reference_transf
			CALL SUB_LIMITIREF(uc_phasor_dq, us_phasor_dq, ic_max, zc)
			ucd = REAL(uc_phasor_dq)
			ucq = AIMAG(uc_phasor_dq)
			
			! AC-side converter current, voltage, power 
			it_phasor_dq = (us_phasor_dq - uc_phasor_dq)/zc
			it_abs = ABS(it_phasor_dq)
			pc = REAL(uc_phasor_dq*CONJG(it_phasor_dq))
			
			CALL SUB_GETPLOSS(ploss, pc, it_abs, aloss, bloss, c_inv, c_rect)
			pdc = -(pc + ploss) ! = ydc? - input of the DC-grid model
			
			! set voltage to zero if tripping
			IF (issimstop.EQ.1) THEN		  
				ucd = 0.0
				ucq = 0.0
				uc_phasor_dq = CMPLX(ucd,ucq)
				pdc = 0.0
				pc = 0.0
				qc = 0.0
			END IF
			
			uc_phasor_RI = uc_phasor_dq*reference_transf
			
			SPEED(I_MACH) = wvsm
			ANGLE(I_MACH) = deltavsm 					! in (deg), not given in rads!!
			ETERM(I_MACH) = ABS(us_phasor_dq)
			PELEC(I_MACH) = pc 							! in pu system rating
			QELEC(I_MACH) = qc 							! in pu system rating
			PMECH(I_MACH) = -pdc 						! in pu machine rating
			EFD(I_MACH) = ABS(uc_phasor_dq)
					  
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
	icd = xd ! sólo era para identificar mejor los estados
	icq = xq

	! VARs	  
	VAR(I_VAR) = addps_ref  				! references
	VAR(I_VAR+1) = addqs_ref 
	VAR(I_VAR+2) = udc_ref  
	VAR(I_VAR+3) = uc_ref 
	VAR(I_VAR+4) = qc_ref 
	VAR(I_VAR+5) = REAL(uc_phasor_RI) 		! output
	VAR(I_VAR+6) = AIMAG(uc_phasor_RI)
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

SUBROUTINE SUB_VIILIMIT(uc_phasor_dq, us_phasor_dq, ic_max, zc, KPVI, SIGMAXR)

	REAL :: it_abs, phi_it, phi_zc
	COMPLEX :: it_phasor_dq, uc_phasor_dq, uc_phasor_dq, zc
	REAL :: ic_max
	
	it_phasor_dq = (uc_phasor_dq - us_phasor_dq)/zc
	it_abs = abs(it_phasor_dq)
	deltait_abs = MAX(it_abs-ic_max,0.0)
	
	zvi = deltait_abs*KPVI*CMPLX(1,SIGMAXR)
	uc_phasor_dq = uc_phasor_dq - zvi*it_phasor_dq
	
END
      
SUBROUTINE SUB_LIMITIREF(uc_phasor_dq, us_phasor_dq, ic_max, zc)

	REAL :: it_abs, phi_it, phi_zc
	COMPLEX :: it_phasor_dq, uc_phasor_dq, uc_phasor_dq, zc
	REAL :: ic_max
	
	it_phasor_dq = (uc_phasor_dq - us_phasor_dq)/zc
	it_abs = abs(it_phasor_dq)
	phi_it = ATAN2(AIMAG(it_phasor_dq),REAL(it_phasor_dq))
	phi_zc = ATAN2(AIMAG(zc),REAL(zc))
	IF (it_abs.GT.ic_max) THEN
		phi_it = ASIN(AIMAG(us_phasor_dq)/ic_max/abs(zc))-phi_zc 	! ucq = usq + i*zc*sin(phiit+phizc) = 0 
		it_phasor_dq = it_abs*CMPLX(COS(phi_it),SIN(phi_it))
	END IF
	uc_phasor_dq = us_phasor_dq + it_phasor_dq*zc
	
	
END

SUBROUTINE SUB_OVDCPROT(issimstop, udc, timenow, tinitial1, counter1, UDC_MAX, TUDCMAX)

	IMPLICIT NONE
	INTEGER :: issimstop, counter1
	REAL :: udc, timenow, tinitial1
	REAL :: UDC_MAX, TUDCMAX
	
	IF (udc.GT.UDC_MAX) THEN 

		IF (counter1 .EQ. 0.0) THEN
			tinitial1 = timenow
			counter1 = 1
		END IF  

		IF ((timenow - tinitial1).GT.TUDCMAX) THEN
			issimstop = 1
		END IF 
	
	ELSE
		tinitial1 = timenow
		counter1 = 0
	END IF
END

SUBROUTINE SUB_UVDCPROT(issimstop, udc, timenow, tinitial1, counter1, UDC_MIN, TUDCMIN)

	IMPLICIT NONE
	INTEGER :: issimstop, counter1
	REAL :: udc, timenow, tinitial1
	REAL :: UDC_MIN, TUDCMIN
	
	IF (udc.LT.UDC_MIN) THEN 

		IF (counter1 .EQ. 0.0) THEN
			tinitial1 = timenow
			counter1 = 1
		END IF  

		IF ((timenow - tinitial1).GT.TUDCMIN) THEN
			issimstop = 1
		END IF 
	
	ELSE
		tinitial1 = timenow
		counter1 = 0
	END IF
END