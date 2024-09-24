! Dynamic model of a the supplementary reactive power control
! 
! The supplementary reactive power control includes:
! - AC voltage control
! - Frequency-based POD
! 
! The supplementary reactive power control is implemented as an stabilizer-type model.

SUBROUTINE SQWDRD(I_MACH,I_SLOT)

	!DEC$ ATTRIBUTES DLLEXPORT, DECORATE, ALIAS: "SQWDRD":: SQWDRD
	!DEC$ ATTRIBUTES REFERENCE :: I_MACH,I_SLOT

	INCLUDE 'COMON4.INS'
	IMPLICIT none
	
	! Declaration
	! ----------- 	
	INTEGER IB, I_SLOT, I_MACH, I_BUS
	INTEGER I_ICON, I_CON, I_STATE, I_VAR
	INTEGER I_VARCONV 
	INTEGER ISVACCNTRL, ISFCNTRL
	INTEGER delta_ac_aux, delta_f_aux
	INTEGER m_ac, m_f

	REAL KAC_NOM, Kac ! CONs
	REAL KQPOD_NOM, Kqpod ! CONs
	REAL TF, TW, UAC_THRES, W_THRES ! CONs
	REAL ADDQs_MAX_NOM, ADDQs_MIN_NOM, RAMP_QMAX_NOM ! CONs
	REAL addqs_max, addqs_min, ramp_qmax ! CONs
	REAL qs_ini, deltaw_ref, deltaw_bus, qs_ref, ew, us_ref, us! VARs
	REAL d_addqs_ref_supp, addqs_ref, addqs_ref_ac, addqs_ref_supp, addqs_ref_supp_previous
	REAL w_fil, d_w_fil, xw, d_xw ! STATEs and derivatives
	REAL DELTAT ! time step
	REAL uw
	REAL USUPPCNTRL_THRES

	INTEGER IERR

	! Parameter assignment
	! --------------------
	IB = NUMTRM(I_MACH) ! don't need this
	I_BUS = NUMBUS(IB) ! do not confuse with IB !!! ;)

	! Index
	I_CON=STRTIN(1,I_SLOT)
	I_STATE=STRTIN(2,I_SLOT)
	I_VAR=STRTIN(3,I_SLOT)
	I_ICON=STRTIN(4,I_SLOT)

	CALL MDLIND(I_BUS, MACHID(I_MACH), 'GEN', 'VAR', I_VARCONV, IERR) ! see API 7-73

	! ICON assignment
	ISVACCNTRL = ICON(I_ICON) 	! AC voltage control: 1, else: 0
	ISFCNTRL = ICON(I_ICON+1) 	! f-based POD: 1, else: 0

	! CON assignment
	TF = CON(I_CON) 						! Frequency filter time constant [s]
	TW = CON(I_CON+1) 						! Washout filter time constant [s]
	KAC_NOM = CON(I_CON+2) 					! AC-voltage droop gain
	KQPOD_NOM = CON(I_CON+3) 				! f-based POD gain (p.u-nom)
	UAC_THRES = CON(I_CON+4) 				! voltage threshold to enable voltage control [p.u] 
	W_THRES = CON(I_CON+5) 					! frequency threshold to enable POD [p.u] 
	ADDQs_MAX_NOM = CON(I_CON+6) 			! Qmax (p.u-nom)
	ADDQs_MIN_NOM = CON(I_CON+7) 			! Qmin (p.u-nom)
	RAMP_QMAX_NOM = CON(I_CON+8)			! maximum active power derivative (p.u/s)
	USUPPCNTRL_THRES = CON(I_CON+9) 		! Q-control strategy is activated only if V>USUPPCNTRL_THRES


	! VARs
	qs_ini = VAR(I_VAR) 
	deltaw_ref = VAR(I_VAR+1)
	deltaw_bus = VAR(I_VAR+2)
	qs_ref = VAR(I_VAR+3)
	addqs_ref = VAR(I_VAR+4)
	d_addqs_ref_supp = VAR(I_VAR+5)
	us_ref = VAR(I_VAR+6)
	us = VAR(I_VAR+7)
	m_ac = VAR(I_VAR+8)
	addqs_ref_ac = VAR(I_VAR+9)
	uw = VAR(I_VAR+10)
    
	! STATEs      
	w_fil = STATE(I_STATE) 	! filtered frequency
	xw = STATE(I_STATE+1) 	! state of the washout filter
    
	! DSTATEs     
	d_w_fil = DSTATE(I_STATE) 
	d_xw = DSTATE(I_STATE+1) 

	! Other variables and parameters	

	Kac = KAC_NOM*SBASE/MBASE(I_MACH)				! change to system rating
	addqs_max = ADDQs_MAX_NOM*MBASE(I_MACH)/SBASE
	addqs_min = ADDQs_MIN_NOM*MBASE(I_MACH)/SBASE
	ramp_qmax = RAMP_QMAX_NOM*MBASE(I_MACH)/SBASE
	Kqpod = KQPOD_NOM*MBASE(I_MACH)/SBASE

	CALL DSRVAL('DELT', 1, DELTAT, IERR)
	
	! AC voltage
	us = ABS(VOLT(IB)) ! variable for the converter model
	deltaw_bus = BSFREQ(IB)


	SELECT CASE (MODE)

		CASE (1) 
		   
			! Inicialization 
			! ========================
			   
			! Variables   
			us_ref = us

			deltaw_ref = 0.0
			deltaw_bus = deltaw_ref
			w_fil = deltaw_ref

			qs_ini = VAR(I_VARCONV+1) ! first variable form SVSCWI model
			qs_ref = qs_ini

			addqs_ref_ac = 0.0
			addqs_ref_supp = 0.0
			addqs_ref = addqs_ref_ac + addqs_ref_supp
			addqs_ref_supp_previous = addqs_ref_supp

			! States
			w_fil = deltaw_ref
			xw = w_fil 

		CASE (2) 
		   
			! Compute derivatives
			! ===================
			ew = deltaw_ref - deltaw_bus	
			IF (TF.LT.(2*DELTAT)) THEN
				d_w_fil = 0.0
				w_fil = ew
			ELSE
				d_w_fil = (-w_fil + ew)/TF
			END IF	

			IF (TW.EQ.0.0) THEN	
				d_xw = 0
				xw = 0.0
			ELSE
				d_xw  = (-xw + w_fil)/TW
			END IF

		CASE (3) 

			! Compute output
			! ==============

			addqs_ref_supp_previous = addqs_ref_supp
			ew = deltaw_ref - deltaw_bus
			uw = w_fil

			IF (ABS(us-us_ref).GE.UAC_THRES) THEN
			  delta_ac_aux = 1
			ELSE
			  delta_ac_aux = 0
			END IF

			IF (ABS(w_fil).GE.W_THRES) THEN ! before: IF (ABS(w_fil-w_ini_fil).GE.W_THRES) THEN
			  delta_f_aux = 1
			ELSE
			  delta_f_aux = 0
			END IF
			   

			m_ac = ISVACCNTRL*delta_ac_aux
			m_f = ISFCNTRL*delta_f_aux
			
			addqs_ref_ac = m_ac*(1.0/Kac)*(us_ref - us) ! AC-voltage droop
			addqs_ref_supp = -m_f*(Kqpod)*(w_fil - xw) 
			
			! reactive power limits
			addqs_ref_supp = MIN(MAX(addqs_ref_supp,addqs_min),addqs_max)

			! reactive power drivative limits
			d_addqs_ref_supp = (addqs_ref_supp-addqs_ref_supp_previous)/DELTAT ! derivative
			addqs_ref_supp = addqs_ref_supp_previous + MIN(MAX(d_addqs_ref_supp,-ramp_qmax),ramp_qmax)*DELTAT
			
			addqs_ref = addqs_ref_ac + addqs_ref_supp
			IF (us.LT.USUPPCNTRL_THRES) THEN
			  addqs_ref = 0.00
			END IF

			VAR(I_VARCONV+31) = addqs_ref ! put the REACTIVE power reference in the converter

		CASE (4) 

			! Update number of STATEs.
			! ========================
			NINTEG = MAX(NINTEG,I_STATE+1)

		CASE (5)
			! Reporting mode
			! ==============

		CASE DEFAULT

	END SELECT

	! RE-ASSIGN VARIABLES
	! -------------------

	! VARs
	VAR(I_VAR) = qs_ini 
	VAR(I_VAR+1) = deltaw_ref
	VAR(I_VAR+2) = deltaw_bus
	VAR(I_VAR+3) = qs_ref
	VAR(I_VAR+4) = addqs_ref 
	VAR(I_VAR+5) = d_addqs_ref_supp
	VAR(I_VAR+6) = us_ref
	VAR(I_VAR+7) = us
	VAR(I_VAR+8) = m_ac
	VAR(I_VAR+9) = addqs_ref_ac
	VAR(I_VAR+10) = uw

	! STATEs 
	STATE(I_STATE) = w_fil 
	STATE(I_STATE+1) = xw 

	! DSTATEs 
	DSTATE(I_STATE) = d_w_fil 
	DSTATE(I_STATE+1) = d_xw 

END