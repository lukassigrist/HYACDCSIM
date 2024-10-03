! Dynamic model of a the supplementary active power control
! 
! The supplementary active power control includes:
! - distributed DC voltage control
! - frequency controls including droop control, time optimal control Lyapunov function and synthetic inertia
! 
! The supplementary active power control is implemented as an excitation-type model.

SUBROUTINE SPWDRD(I_MACH,I_SLOT)

	!DEC$ ATTRIBUTES DLLEXPORT, DECORATE, ALIAS: "SPWDRD":: SPWDRD
	!DEC$ ATTRIBUTES REFERENCE :: I_MACH,I_SLOT

	INCLUDE 'COMON4.INS'
	IMPLICIT none

	! Declaration
	! -----------
	INTEGER IB, I_SLOT, I_MACH, I_BUS
	INTEGER I_ICON, I_CON, I_STATE, I_VAR, I_VARCONV 
	INTEGER ISDCCNTRL, ISFCNTRL, ISALPHACNTRL
	INTEGER delta_f, delta_toclf, delta_angle 
	INTEGER delta_dc_aux, delta_f_aux, delta_a_aux, delta_toclf_aux
	INTEGER m_dc, m_f, m_a, m_toclf
	INTEGER IERR
	
	REAL KDC_NOM, KF_NOM, K_ALPHA_NOM, KP_NOM 
	REAL TF, TW, UDC_THRSH, W_THRSH, ALPHA_THRSH 
	REAL ADDPs_MAX_NOM, ADDPs_MIN_NOM, RAMP_PMAX_NOM, BP_ERIK_NOM 
	
	REAL DELTAT ! time step
	
	REAL :: Kdc, Kf, K_alpha 
	REAL :: addps_max, addps_min, ramp_pmax, bp_toclf 
	REAL :: ps_ini, ps_ref, deltaw_ref, deltaw_bus, deltau_ini, theta_bus, ew, udc_ref, udc
	REAL :: addps_ref, addps_ref_dc, addps_ref_freq, addps_ref_inertia, addps_ref_supp, d_addps_ref_supp, addps_ref_supp_previous
	REAL :: w_fil, d_w_fil, xw, d_xw ! STATEs and derivatives
	
	REAL :: omega_est_i


	! Parameter assignment
	! --------------------
	IB = NUMTRM(I_MACH) 
	I_BUS = NUMBUS(IB) 

	! Index
	I_CON=STRTIN(1,I_SLOT)
	I_STATE=STRTIN(2,I_SLOT)
	I_VAR=STRTIN(3,I_SLOT)
	I_ICON=STRTIN(4,I_SLOT)

	CALL MDLIND(I_BUS, MACHID(I_MACH), 'GEN', 'VAR', I_VARCONV, IERR) ! see API 7-73

	! ICONs
	ISDCCNTRL = ICON(I_ICON) 		! enable dc-voltage droop
	ISFCNTRL = ICON(I_ICON+1) 		! frequency control: 1, time-optimal control lyapunov function [Eriksson]: 2, AC-line: 3, disabled: 0
	ISALPHACNTRL = ICON(I_ICON+2) 	! enable synthetic intertia
    
	! CONs
	TF = CON(I_CON) ! freq filter time constant [s]
	TW = CON(I_CON+1) ! Washout filter time constant [s]
	KDC_NOM = CON(I_CON+2) ! DC-voltage droop gain
	KP_NOM = CON(I_CON+3) ! prop. gain (p.u-nom)
	K_ALPHA_NOM = CON(I_CON+4) ! synthetic inertia gain
	UDC_THRSH = CON(I_CON+5) ! freq-deviation threshold to enable freq-droop [p.u] 
	W_THRSH = CON(I_CON+6) ! freq-deviation threshold to enable freq-droop [p.u] 
	ALPHA_THRSH = CON(I_CON+7) ! freq-derivative threshold to enable freq-droop [p.u/s]
	ADDPs_MAX_NOM = CON(I_CON+8) ! Pmax (p.u-nom)
	ADDPs_MIN_NOM = CON(I_CON+9) ! Pmin (p.u-nom)
	RAMP_PMAX_NOM = CON(I_CON+10) ! maximum active power derivative (p.u/s)
	BP_ERIK_NOM = CON(I_CON+11) ! maximum active power derivative (p.u/s)

	! VARs  
	ps_ini = VAR(I_VAR) ! Reference active power, coincides with load flow values.
	deltaw_ref = VAR(I_VAR+1)
	deltau_ini = VAR(I_VAR+2)
	ps_ref = VAR(I_VAR+3)
	addps_ref = VAR(I_VAR+4)
	d_addps_ref_supp = VAR(I_VAR+5)
	udc_ref = VAR(I_VAR+6)
	udc = VAR(I_VAR+7)
	m_dc = VAR(I_VAR+8)
	addps_ref_dc = VAR(I_VAR+9)
	addps_ref_freq = VAR(I_VAR+10)
	addps_ref_inertia = VAR(I_VAR+11)
	omega_est_i = VAR(I_VAR+14)
	
	! STATEs     
	w_fil = STATE(I_STATE) ! filtered frequency
	xw = STATE(I_STATE+1) ! state of the washout filter

	! DSTATEs 
	d_w_fil = DSTATE(I_STATE) ! filtered frequency
	d_xw = DSTATE(I_STATE+1) ! state of the washout filter


	! Change gain to 'inverse of the gain' (only conventions)
	KF_NOM = 1/KP_NOM

	! Change machine to system base
	Kdc = KDC_NOM*SBASE/MBASE(I_MACH)
	Kf = KF_NOM*SBASE/MBASE(I_MACH)
	K_alpha = K_ALPHA_NOM*MBASE(I_MACH)/SBASE
	addps_max = ADDPs_MAX_NOM*MBASE(I_MACH)/SBASE
	addps_min = ADDPs_MIN_NOM*MBASE(I_MACH)/SBASE
	ramp_pmax = RAMP_PMAX_NOM*MBASE(I_MACH)/SBASE
	bp_toclf = BP_ERIK_NOM*MBASE(I_MACH)/SBASE

	! Common variables
	CALL DSRVAL('DELT', 1, DELTAT, IERR) ! obtain time step (DELTAT) [s] (API 7.45)
	
	deltaw_bus = BSFREQ(IB) 							! deviation
	theta_bus = ATAN2(AIMAG(VOLT(IB)),REAL(VOLT(IB))) 	! angle


	! DC-voltage
	udc = VAR(I_VARCONV+22) ! variable for the converter model

	! frequency control and Eriksson's control 
	delta_f = 0
	delta_toclf = 0
	delta_angle = 0
	IF (ISFCNTRL.EQ.1) THEN
		delta_f = 1
	ELSE IF (ISFCNTRL.EQ.2) THEN
		delta_toclf = 1
	ELSE IF (ISFCNTRL.EQ.3) THEN
		delta_angle = 1
		delta_f = 1
		W_THRSH = 0.0
	END IF

	SELECT CASE (MODE)

	CASE (1) 
	   
		! Inicialization 
		! ========================	 

		! Algebraic variables
		udc_ref = VAR(I_VARCONV+27)		! DC voltage reference
		udc = udc_ref
		
		deltaw_ref = 0.0 				! frequency deviation reference
		
		ps_ini = VAR(I_VARCONV) 		! initial power SVSCON
		ps_ref = ps_ini

		addps_ref_dc = 0.0
		addps_ref_freq = 0.0
		addps_ref_inertia = 0.0
		addps_ref_supp = addps_ref_freq + addps_ref_inertia
		addps_ref = addps_ref_dc + addps_ref_supp
		addps_ref_supp_previous = addps_ref_supp
		
		omega_est_i = addps_ref*Kf + deltaw_bus

		! State variables	
		IF (delta_angle.EQ.0) THEN
			deltau_ini = deltaw_ref
		ELSE
			deltau_ini = theta_bus
		END IF	
		
		w_fil = 0.0
		xw = w_fil
		
		!WRITE (LPDEV,*) 'SPWDRD - MTDC at AC bus ',NUMBUS(NUMTRM(I_MACH)),' with id ',MACHID(I_MACH)
		!WRITE (LPDEV,*) 'SPWDRD - CASE 1: Kdc, Kf, K_alpha = ',Kdc, Kf, K_alpha
		!WRITE (LPDEV,*) 'SPWDRD - CASE 1: delta_f, delta_angle = ',delta_f, delta_angle
		!WRITE (LPDEV,*) 'SPWDRD - CASE 1: I_VARCONV,ps_initial = ',I_VARCONV,VAR(I_VARCONV)
					
	CASE (2) 
	   
		! Compute derivatives
		! ===================
		IF (delta_angle.EQ.0) THEN
			ew = deltaw_ref - deltaw_bus
		ELSE
			! thetai - thetaj = thetai0 + deltathetai - (thetaj0 + deltathetaj) =  thetai0 -  thetaj0 + deltathetai - deltathetaj
			! 				  = 2*deltathetai - 2*deltathetaverage
			! deltathetaaverage = (deltathetai + deltathetaj)/2 (from WDELAY)
			ew = 2*(theta_bus - deltau_ini) - 2*deltaw_ref
		END IF			
				
		
		
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
		addps_ref_supp_previous = addps_ref_freq + addps_ref_inertia

		! thresholds for activation
		IF (ABS(udc-udc_ref).GE.UDC_THRSH) THEN
			delta_dc_aux = 1
		ELSE
			delta_dc_aux = 0
		END IF

		IF (ABS(d_w_fil).GE.ALPHA_THRSH) THEN
			delta_a_aux = 1
		ELSE
			delta_a_aux = 0
		END IF

		IF (ABS(w_fil).GE.W_THRSH) THEN
			delta_f_aux = 1
			
			delta_toclf_aux = 1
			! This only for the controller similar to Eriksson proposal
			IF ((w_fil).GE.0) THEN
				delta_toclf_aux = 1
			ELSE IF ((w_fil).LT.0) THEN
				delta_toclf_aux = -1
			END IF
		ELSE
			delta_f_aux = 0
			delta_toclf_aux = 0
		END IF

		m_dc = ISDCCNTRL*delta_dc_aux
		m_f = delta_f*delta_f_aux
		m_a = ISALPHACNTRL*delta_a_aux
		m_toclf = delta_toclf*delta_toclf_aux

		addps_ref_dc = -m_dc*(1.0/Kdc)*(udc_ref - udc)
		
		addps_ref_freq = m_f*(1.0/Kf)*(w_fil - xw) + m_toclf*bp_toclf ! note that only one of them is activated
		addps_ref_inertia = m_a*(K_alpha)*d_w_fil
		addps_ref_supp = addps_ref_freq + addps_ref_inertia
		
		! power limits
		addps_ref_supp = MIN(MAX(addps_ref_supp,addps_min),addps_max)

		! power drivative limits
		d_addps_ref_supp = (addps_ref_supp-addps_ref_supp_previous)/DELTAT ! derivative
		addps_ref_supp = addps_ref_supp_previous + MIN(MAX(d_addps_ref_supp,-ramp_pmax),ramp_pmax)*DELTAT
		
		addps_ref = addps_ref_dc + addps_ref_supp
		
		omega_est_i = addps_ref*Kf + deltaw_bus

		VAR(I_VARCONV+29) = addps_ref ! put the power reference in the converter

	CASE (4) 

		! Update number of STATEs
		! ========================
		NINTEG = MAX(NINTEG,I_STATE+1)

	CASE (5)

	CASE DEFAULT

	END SELECT

	! VARs
	VAR(I_VAR) = ps_ini				! Reference active power, coincides with load flow values
	VAR(I_VAR+1) = deltaw_ref		! From WDELAY (else 0.0)
	VAR(I_VAR+2) = deltau_ini
	VAR(I_VAR+3) = ps_ref
	VAR(I_VAR+4) = addps_ref 
	VAR(I_VAR+5) = d_addps_ref_supp
	VAR(I_VAR+6) = udc_ref
	VAR(I_VAR+7) = udc
	VAR(I_VAR+8) = m_dc
	VAR(I_VAR+9) = addps_ref_dc
	VAR(I_VAR+10) = addps_ref_freq
	VAR(I_VAR+11) = addps_ref_inertia
	VAR(I_VAR+14) = omega_est_i
    
	! STATEs 
	STATE(I_STATE) = w_fil 			! filtered frequency deviation
	STATE(I_STATE+1) = xw 			! state of the washout filter
  
	! DSTATEs       
	DSTATE(I_STATE) = d_w_fil ! filtered frequency
	DSTATE(I_STATE+1) = d_xw ! state of the washout filter

 END