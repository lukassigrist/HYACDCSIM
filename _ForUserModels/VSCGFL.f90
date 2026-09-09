! Dynamic model of a grid-following VSC of a VSC-based MTDC grid
!
! The VSC converter computes the currents injected into the ac grid and the power injected into the DC grid. The DC grid
! model, implemented as a governor type model, uses the power injected into the DC grid as an input.
!
! The VSC converter model is based on a grid-following converter. The current dynamics are represented by a first-order system.
! The control is formulated in a dq reference system aligned with the terminal voltage. No PLL is used so far. The current
! references are computed from active and reactive power references.
! 
! us,dq = us,RI*exp(-j*deltas) = us + j*0.0
! icref,dq = (ssref/us,dq)* = psref/us - j*qsref/us
! ic,dq = [1/(1+s*Td) 0; 0 1/(1+s*Tq)]*icref,dq 
!
! Supervisory controls such as active and reactive power controls, DC-voltage control, etc., essentially modify the active and 
! reactive power references. These supervisory controls are implemented by means of different user models:
!
! - stabilizer (SQWDRD), modifying the reactive power reference,
! - excitation (SPWDRD), modifying the active power reference, and 
! - governor (DCGRID), modeling the DC grid dynamics.
!
! Note that the calling sequence for the plant models is: 1. Generator models, 2. Current compensating models, 
! 3. Excitation stabilizer models, 4. Excitation system models, 5. Turbine-governor models. Turbine governor, stabilizer
! and excitation limiter models have no initialization duties other than STATEs and VARs.
! 
! The grid-following VSC model, VSCGFL, is implemented as a coordinated-call, current injecting generator model (IC = 1 
! and IT = 1). Coordinated-call implementation is required to cancel the effect of the admittance of the Norton
! equivalent applied to generators, which leads to a pure current source.  
!
! ======================================================================================
! MODULE DECLARATION
! ======================================================================================
MODULE MOD_VSCGFL_INTERNAL
    ! Internal variables
    ! This stores the VAR/STATE index of the first DCGRID model found.
    INTEGER :: I_VAR_DCGRID_G = -1
    INTEGER :: I_STATE_DCGRID_G = -1
END MODULE MOD_VSCGFL_INTERNAL

! ======================================================================================
! CURRENT INJECTIONS FOR NETWORK SOLUTION
! ======================================================================================
SUBROUTINE TSCGFL(I_MACH,I_SLOT)

    INCLUDE 'COMON4.ins'
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
    IB = NUMTRM(I_MACH)
    I_VAR = STRTIN(3,I_SLOT)

    ! VARs
    icd = VAR(I_VAR+15)
    icq = VAR(I_VAR+16)

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
        ! ISORCE: Norton equivalent source current in pu of SBASE
        ! ic: current at generator terminal bus in pu SBASE
        ! ut: voltage at generator terminal bus in pu
        ! ZSORCE: in pu of MBASE -> zc = ZSORCE*SBASE/MBASE
        ic_phasor_RI = ISORCE(I_MACH) - us_phasor_RI/zc
        icd = REAL(ic_phasor_RI/reference_transf)
        icq = AIMAG(ic_phasor_RI/reference_transf)

        ss_phasor = us_phasor_RI*CONJG(ic_phasor_RI)

    CASE (3)
        ! Current injections
        ! ------------------
        ! Set ISORCE such that the current through 1/ZSORCE is cancelled.
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
    VAR(I_VAR+15) = icd
    VAR(I_VAR+16) = icq

END SUBROUTINE TSCGFL

! ==========================================================
! SOLUTION OF DIFFERENTIAL EQUATIONS
! ==========================================================
SUBROUTINE VSCGFL(I_MACH,I_SLOT)

    USE MOD_TFBLOCKS                                            ! use module of generic transfer functions
    USE MOD_PROTECTION                                          ! use module of protection functions
    USE MOD_READHYADCSIM                                        ! use module for reading the HYACDCSIM text files
    USE MOD_MISC                                                ! use module for miscellaneous functions
    USE MOD_VSCGFL_INTERNAL                                     ! use module for global VSCGFL-related variables
    USE MOD_DCGRID, ONLY: REGISTER_DCG               ! use automatic DC-grid registration for offset calculation
    INCLUDE 'COMON4.INS'                                        ! common PSS/e variables and modules
    IMPLICIT NONE

    ! Declaration
    ! -----------
    INTEGER IB, I_SLOT, I_MACH, I_VAR, I_CON, I_ICON, I_STATE        ! PSS/e indices
    INTEGER I_VAR_DCGRID, I_STATE_DCGRID                             ! local index of the DCGRID VAR/STATE
    INTEGER DCNTRLTYPE, QCNTRLTYPE, ILIMITPRIORITY
    INTEGER i, idx_converter, idx_DCGRID
    INTEGER, ALLOCATABLE :: temp_idx_converter(:), temp_idx_DCGRID(:)

    INTEGER istripvsc                                                ! disconnects unit and sets current injections to 0
    INTEGER dcontroltype_aux, qcontroltype_aux
    INTEGER NDCBUS, NDCLINES                                         ! # of buses and lines of the dc-grid
    INTEGER NDCBUS_PREVIOUS, NDCLINES_PREVIOUS
    INTEGER DCGRID_VAR_OFFSET

    INTEGER, ALLOCATABLE :: m_VSCACDCBUS(:,:)                  ! ac and dc buses of each converter
    INTEGER ierr
    REAL PI
    PARAMETER (PI=3.14159265358979)
    REAL DELTAT

    REAL us, deltas, deltas0
    REAL ps0, qs0
    REAL icd, icq
    REAL icd0, icq0
    REAL icD_out, icQ_out
    REAL ps, qs

    REAL xd, xq, md, mq, ndc, nq, eta_d                              ! state variables
    REAL d_xd, d_xq, d_md, d_mq, d_ndc, d_nq, d_eta_d                ! derivatives of state variables

    REAL antiwindupicd, antiwindupicq

    REAL ic, ic_abs_puconv

    REAL psref, qsref, usref
    REAL icdref, icqref
    REAL udcref

    REAL idc, idc0, udc, udc0
    REAL pdc

    REAL deltapsref
    REAL deltaqsref
    REAL ec, deltac
    REAL pvsc, qvsc

    COMPLEX us_phasor_RI, ec_phasor_RI, ic_phasor_RI, sc_phasor
    COMPLEX reference_transf

    CHARACTER(1) IDGRID

    REAL, ALLOCATABLE :: v_udc0(:,:)

    REAL TAU
    REAL KD_Pps, KD_Pudc, KD_Ips, KD_Iudc, KQ_Pqs, KQ_Iqs, KQ_Pus, KQ_Ius
    REAL ICMAX
    REAL PSMAX, PSMIN, QSMAX, QSMIN
    REAL ALOSS, BLOSS, CLOSSRECT, CLOSSINV
    REAL UDC_MAX, UDC_MIN, TUDCMAX, TUDCMIN

    REAL icmaxpu
    REAL psmaxpu, psminpu, qsmaxpu, qsminpu
    REAL alosspu, blosspu, clossrectpu, clossinvpu
    REAL ZBASE
    COMPLEX zc
    REAL ploss

    INTEGER :: counter_uvdc, counter_ovdc
    REAL :: tinitial_uvdc, tinitial_ovdc

    REAL fmodulationpwm, FMODULATIONPWMMAX

    REAL :: yetad, yndc, ymd, ynq, ymq

    ! Parameter assignment
    ! --------------------
    ! Indexes
    IB = NUMTRM(I_MACH)
    I_CON = STRTIN(1,I_SLOT)
    I_STATE = STRTIN(2,I_SLOT)
    I_VAR = STRTIN(3,I_SLOT)
    I_ICON = STRTIN(4,I_SLOT)

    ! CONs (if in pu, in pu of converter rating)
    TAU = CON(I_CON)                      ! Inverter time constant (e.g., 0.1)
    KD_Pps = CON(I_CON+1)                 ! Ps control
    KD_Ips = CON(I_CON+2)
    KD_Pudc = CON(I_CON+3)                ! Udc control 
    KD_Iudc = CON(I_CON+4)
    KQ_Pqs = CON(I_CON+6)                 ! Qs control
    KQ_Iqs = CON(I_CON+7)
    KQ_Pus = CON(I_CON+8)                 ! Us control
    KQ_Ius = CON(I_CON+9)
    ICMAX = CON(I_CON+10)                 ! Maximum inverter current (e.g., 1.1)
    PSMAX = CON(I_CON+11)
    PSMIN = CON(I_CON+12)
    QSMAX = CON(I_CON+13)
    QSMIN = CON(I_CON+14)
    UDC_MAX = CON(I_CON+15)
    UDC_MIN = CON(I_CON+16)
    FMODULATIONPWMMAX = CON(I_CON+17)
    ALOSS = CON(I_CON+18)             ! constant converter loss coefficient (MW): ploss = aloss + bloss*ic + c*ic^2
    BLOSS = CON(I_CON+19)             ! linear converter loss coefficient (kV): ploss = aloss + bloss*ic + c*ic^2
    CLOSSRECT = CON(I_CON+20)       ! rectifier quadratic converter loss coefficient (ohm): ploss = aloss + bloss*ic + c*ic^2
    CLOSSINV = CON(I_CON+21)        ! constant converter loss coefficient (ohm): ploss = aloss + bloss*ic + c*ic^2
    TUDCMIN = CON(I_CON+22)              ! DC undervoltage protection delay
    TUDCMAX = CON(I_CON+23)              ! DC overvoltage protection delay

    ! ICONs
    DCNTRLTYPE = ICON(I_ICON)          	 ! Ps-control: 1 , Udc-control: 2, Deltas-control: 3
    QCNTRLTYPE = ICON(I_ICON+1)          ! Qs-control: 1 , Us-control: 2
    ILIMITPRIORITY = ICON(I_ICON+2)      ! Current limit: P-priority: 1, Q-priority: 2 and P-Q equal priority: 3 or any other integer
    NDCBUS = ICON(I_ICON+3)              ! Number of converters of the DC grid
    IDGRID = CHRICN(I_ICON+5)            ! DC grid identifier
    NDCLINES = ICON(I_ICON+4)

    CALL REGISTER_DCG(IDGRID, NDCBUS, NDCLINES, NDCBUS_PREVIOUS, NDCLINES_PREVIOUS, ierr)

    IF (ierr .NE. 0) THEN
        WRITE (LPDEV,*) 'VSCGFL - ERROR: automatic DC-grid offset registration failed. IDGRID = ', IDGRID
        RETURN
    END IF

    ! DCGRID VAR block order is: Pdc(NDCBUS), Idc(NDCBUS), Udc(NDCBUS), Icc(NDCLINES)
    DCGRID_VAR_OFFSET = 3*NDCBUS_PREVIOUS + NDCLINES_PREVIOUS

    IF (DCNTRLTYPE.EQ.3) THEN           ! if feeding a passive grid, reactive power control is a voltage control
        QCNTRLTYPE = 2
    END IF

    ! VARs
    deltapsref = VAR(I_VAR)     ! From SPWDRD
    deltaqsref = VAR(I_VAR+1)   ! From SQWDRD    
    udcref = VAR(I_VAR+2)       ! To SPWDRD
    udc = VAR(I_VAR+3)          ! To SPWDRD
    pdc = VAR(I_VAR+4)          ! To DDCGRD/SDCGRD
    
    ps0 = VAR(I_VAR+5)          
    qs0 = VAR(I_VAR+6)          
    deltas0 = VAR(I_VAR+7)
    icd0 = VAR(I_VAR+8)
    icq0 = VAR(I_VAR+9)

    usref = VAR(I_VAR+10)
    icdref = VAR(I_VAR+11)
    icqref = VAR(I_VAR+12)
    
    ps = VAR(I_VAR+13)
    qs = VAR(I_VAR+14)
    icd = VAR(I_VAR+15)
    icq = VAR(I_VAR+16)

    istripvsc = VAR(I_VAR+17)
    antiwindupicd = VAR(I_VAR+18)
    antiwindupicq = VAR(I_VAR+19)
    ic_abs_puconv = VAR(I_VAR+20)
    dcontroltype_aux = VAR(I_VAR+21)
    qcontroltype_aux = VAR(I_VAR+22)

    fmodulationpwm = VAR(I_VAR+23)
    counter_uvdc = VAR(I_VAR+24)
    tinitial_uvdc = VAR(I_VAR+25)
    counter_ovdc = VAR(I_VAR+26)
    tinitial_ovdc = VAR(I_VAR+27)
    idx_converter = VAR(I_VAR+28)
    
    ! STATEs
    xd = STATE(I_STATE)                      ! icd
    xq = STATE(I_STATE+1)                    ! icq
    md = STATE(I_STATE+2)                    ! d-integral state variable type1
    mq = STATE(I_STATE+3)                    ! q-integral state variable type1
    ndc = STATE(I_STATE+4)                   ! d-integral state variable type2
    nq = STATE(I_STATE+5)                    ! q-integral state variable type2
    eta_d = STATE(I_STATE+6)                 ! d-integral state variable type3 (passive grid)

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

    ALLOCATE(v_udc0(NDCBUS,1))                ! allocate size
    ALLOCATE(m_VSCACDCBUS(NDCBUS,2))
    
    icmaxpu = ICMAX*MBASE(I_MACH)/SBASE       ! current limit in system base
    psmaxpu = PSMAX/SBASE                     ! power limits in system base
    psminpu = PSMIN/SBASE
    qsmaxpu = QSMAX/SBASE
    qsminpu = QSMIN/SBASE

    ZBASE = BASVLT(IB)**2/SBASE              ! AC-side Zbase in ohms
	
	alosspu = ALOSS*MBASE(I_MACH)/SBASE                  ! losses coefs in system base
	blosspu = BLOSS
	clossrectpu = CLOSSRECT*SBASE/MBASE(I_MACH) 
	clossinvpu = CLOSSINV*SBASE/MBASE(I_MACH) 

    zc = ZSORCE(I_MACH)*SBASE/MBASE(I_MACH)  ! conexion impedance in system base

    I_VAR_DCGRID = I_VAR_DCGRID_G            ! stored DCGRID VAR base index
    I_STATE_DCGRID = I_STATE_DCGRID_G        ! stored DCGRID STATE base index

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
    IF (DCNTRLTYPE.EQ.3) THEN      ! If feeding a passive grid, the VSC is a AC slack bus
        reference_transf = CMPLX(1.0, 0.0)
    ELSE
        reference_transf = CMPLX(COS(deltas),SIN(deltas))
    END IF

    icd = REAL(ic_phasor_RI/reference_transf)
    icq = AIMAG(ic_phasor_RI/reference_transf)

    ! Present power and current references
    psref = MIN(MAX(ps0 + deltapsref,psminpu),psmaxpu)
    qsref = MIN(MAX(qs0 + deltaqsref,qsminpu),qsmaxpu)

    SELECT CASE (MODE)

    CASE (1)

        ! Initialization
        ! ==============

        ! Get initial value of the DC voltage state from .txt files
        CALL SUB_READDCBUS(v_udc0, m_VSCACDCBUS, IDGRID, NDCBUS)
        
        ! Get VAR index of the DCGRID governor-type model, since DCGRID now updates DC variables through VARs.
		temp_idx_DCGRID = pack(m_VSCACDCBUS(:,2), m_VSCACDCBUS(:,2) /= -1)
		idx_DCGRID = temp_idx_DCGRID(1)

        ! Get DCGRID base index, supporting both STATE-based and VAR-based implementations.
        I_STATE_DCGRID = 0
        I_VAR_DCGRID   = 0
        
        CALL MDLIND(idx_DCGRID, MACHID(I_MACH), 'GOV', 'STATE', I_STATE_DCGRID, ierr)
        IF (ierr.NE.0) THEN
            I_STATE_DCGRID = 0

            CALL MDLIND(idx_DCGRID, MACHID(I_MACH), 'GOV', 'VAR', I_VAR_DCGRID, ierr)
            IF (ierr.NE.0) THEN
                I_VAR_DCGRID = 0
            END IF
        END IF    
        
        IF ((I_VAR_DCGRID.LE.0).AND.(I_STATE_DCGRID.LE.0)) THEN
            WRITE (LPDEV,*) 'VSCGFL - CASE 1: No DCGRID VAR or STATE at bus ', NUMBUS(IB), ' with id ', MACHID(I_MACH), '.'

        ELSE IF ((I_VAR_DCGRID_G.LT.0).AND.(I_STATE_DCGRID_G.LT.0)) THEN
            I_STATE_DCGRID_G = I_STATE_DCGRID   ! Stored DCGRID STATE base index.
            I_VAR_DCGRID_G   = I_VAR_DCGRID     ! Stored DCGRID VAR base index.      
        END IF

        ! Algebraic variables (VARS)
        ! --------------------------
        deltas0 = deltas
        icd0 = icd
        icq0 = icq
        ps0 = us*icd
        qs0 = -us*icq

        ! AC-side converter current out of range
        ic = SQRT(icd**2 + icq**2)
        IF (ic.GT.icmaxpu) THEN
            WRITE (LPDEV,*) 'VSCGFL - CASE 1: Current exceeds limit of ', icmaxpu, ' pu.'
        END IF

        ! AC-side converter voltages, power
        ec_phasor_RI = us_phasor_RI + zc*ic_phasor_RI
        ec = ABS(ec_phasor_RI)
        deltac = ATAN2(AIMAG(ec_phasor_RI), REAL(ec_phasor_RI))
        pvsc = REAL(ec_phasor_RI*CONJG(ic_phasor_RI))
        pvsc = ps0 + REAL(zc)*(ic**2)

        temp_idx_converter = pack([(i, i = 1, size(m_VSCACDCBUS(:,2)))], m_VSCACDCBUS(:,2) == NUMBUS(IB))
		idx_converter = temp_idx_converter(1)
           
        ! DC-side voltage, current, power
        CALL SUB_COMPUTEPLOSS(ploss, ps0, ic, alosspu, blosspu, clossinvpu, clossrectpu)
        pdc = -(pvsc + ploss)
        udc = v_udc0(idx_converter,1) ! extract initial dc voltage value
        udc0 = udc
        idc0 = pdc/udc0
        idc = idc0

        fmodulationpwm = ec/MAX(udc0,0.0001)
        IF (fmodulationpwm.GT.FMODULATIONPWMMAX) THEN
            WRITE (LPDEV,*) 'VSCGFL - CASE 1: Modulation index exceeds limit of ', FMODULATIONPWMMAX, '.'
        END IF

        ! References
        psref = ps0
        qsref = qs0
        usref = us
        udcref = udc0
        deltapsref = 0.0            ! supplementary controller input (SPWDRD)
        deltaqsref = 0.0            ! supplementary controller input (SQWDRD)
        icdref = icd0
        icqref = icq0

        ! Indexes (stop simulation and anti-windup)
        istripvsc = 0
        antiwindupicd = 1.0
        antiwindupicq = 1.0

        ! State variables (STATES)
        ! ------------------------
        ! Current control
        d_xd = 0.0
        d_xq = 0.0
        CALL SUB_FIRSTORDERWINDUP(icdref,xd,d_xd,icdref,1,-1,DELTAT,TAU,icmaxpu,-icmaxpu)        ! icd
        CALL SUB_FIRSTORDERWINDUP(icqref,xq,d_xq,icqref,1,-1,DELTAT,TAU,icmaxpu,-icmaxpu)        ! icq

        ! Active power related control
        d_eta_d = 0.0
        d_ndc = 0.0
        d_md = 0.0
        yetad = 0.0
        yndc = 0.0
        ymd = 0.0
        CALL SUB_PI(yetad,eta_d,d_eta_d,(deltas0-deltas),2,DELTAT,100.0,1000.0,icmaxpu,-icmaxpu)   ! Deltas control

        CALL SUB_PI(yndc,ndc,d_ndc,(udcref-udc),2,DELTAT,KD_Pudc,KD_Iudc,icmaxpu,-icmaxpu)         ! DC-voltage control

        CALL SUB_PI(ymd,md,d_md,(psref-ps),2,DELTAT,KD_Pps,KD_Ips,icmaxpu,-icmaxpu)                ! Ps control

        ! Reactive power related control
        d_nq = 0.0
        d_mq = 0.0
        ynq = 0.0
        ymq = 0.0
        CALL SUB_PI(ynq,nq,d_nq,(usref-us),2,DELTAT,KQ_Pus,KQ_Ius,icmaxpu,-icmaxpu)                ! Us control

        CALL SUB_PI(ymq,mq,d_mq,(qsref-qs),2,DELTAT,KQ_Pqs,KQ_Iqs,icmaxpu,-icmaxpu)                ! Qs control

        PMECH(I_MACH) = -pdc          ! pu system rating
        SPEED(I_MACH) = BSFREQ(IB)
		ANGLE(I_MACH) = deltac
        ETERM(I_MACH) = us
        PELEC(I_MACH) = ps0           ! pu system rating
        QELEC(I_MACH) = qs0           ! pu system rating

        WRITE (LPDEV,*) 'VSCGFL - CASE 1: Converter ', idx_converter, ' at bus ', NUMBUS(IB), ' with id ', MACHID(I_MACH), ' initialized.'
        WRITE (LPDEV,*) 'VSCGFL - CASE 1: D-control: ', DCNTRLTYPE, ' Q-control: ', QCNTRLTYPE, ' I-limit priority: ', ILIMITPRIORITY


    CASE (2)

        ! Compute derivatives
        ! ===================

        ! Get DC voltage from DCGRID: use STATE if available, otherwise use VAR.
        IF (I_STATE_DCGRID.GT.0) THEN
            udc = STATE(I_STATE_DCGRID + idx_converter + NDCLINES_PREVIOUS + NDCBUS_PREVIOUS - 1)
        ELSE IF (I_VAR_DCGRID.GT.0) THEN
            udc = VAR(I_VAR_DCGRID + DCGRID_VAR_OFFSET + 2*NDCBUS + idx_converter - 1)
        END IF

        ! Current control
        CALL SUB_FIRSTORDERWINDUP(icd,xd,d_xd,icdref,2,0,DELTAT,TAU,icmaxpu,-icmaxpu)          ! icd
        CALL SUB_FIRSTORDERWINDUP(icq,xq,d_xq,icqref,2,0,DELTAT,TAU,icmaxpu,-icmaxpu)          ! icq

        ! Active power related control
        IF (DCNTRLTYPE.EQ.3) THEN
            CALL SUB_PI(yetad,eta_d,d_eta_d,(deltas0-deltas),2,DELTAT,100.0,1000.0,icmaxpu,-icmaxpu)    ! Deltas control
        ELSE IF (DCNTRLTYPE.EQ.2) THEN
            CALL SUB_PI(yndc,ndc,d_ndc,(udcref-udc),2,DELTAT,KD_Pudc,KD_Iudc,icmaxpu,-icmaxpu)          ! DC-voltage control
        ELSE
            CALL SUB_PI(ymd,md,d_md,(psref-ps),2,DELTAT,KD_Pps,KD_Ips,icmaxpu,-icmaxpu)                 ! Ps control
        END IF

        ! Reactive power related control
        IF (QCNTRLTYPE.EQ.2) THEN
            CALL SUB_PI(ynq,nq,d_nq,(usref-us),2,DELTAT,KQ_Pus,KQ_Ius,icmaxpu,-icmaxpu)                 ! Us control
        ELSE
            CALL SUB_PI(ymq,mq,d_mq,(qsref-qs),2,DELTAT,KQ_Pqs,KQ_Iqs,icmaxpu,-icmaxpu)                 ! Qs control
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

        ! Get DC voltage from DCGRID: use STATE if available, otherwise use VAR.
        IF (I_STATE_DCGRID.GT.0) THEN
            udc = STATE(I_STATE_DCGRID + idx_converter + NDCLINES_PREVIOUS + NDCBUS_PREVIOUS - 1)
        ELSE IF (I_VAR_DCGRID.GT.0) THEN
            udc = VAR(I_VAR_DCGRID + DCGRID_VAR_OFFSET + 2*NDCBUS + idx_converter - 1)
        END IF

        ! Current control
        CALL SUB_FIRSTORDERWINDUP(icd,xd,d_xd,icdref,3,0,DELTAT,TAU,icmaxpu,-icmaxpu)          ! icd
        CALL SUB_FIRSTORDERWINDUP(icq,xq,d_xq,icqref,3,0,DELTAT,TAU,icmaxpu,-icmaxpu)          ! icq

        ! Active power related control
        IF (DCNTRLTYPE.EQ.3) THEN
            CALL SUB_PI(yetad,eta_d,d_eta_d,(deltas0-deltas),3,DELTAT,100.0,1000.0,icmaxpu,-icmaxpu)       ! Deltas control
            icdref = icd0 + yetad
        ELSE IF (DCNTRLTYPE.EQ.2) THEN
            CALL SUB_PI(yndc,ndc,d_ndc,(udcref-udc),3,DELTAT,KD_Pudc,KD_Iudc,icmaxpu,-icmaxpu)             ! DC-voltage control
            icdref = icd0 - yndc
        ELSE
            CALL SUB_PI(ymd,md,d_md,(psref-ps),3,DELTAT,KD_Pps,KD_Ips,icmaxpu,-icmaxpu)                    ! Ps control
            icdref = psref/us + ymd
        END IF

        ! Reactive power related control
        IF (QCNTRLTYPE.EQ.2) THEN
            CALL SUB_PI(ynq,nq,d_nq,(usref-us),3,DELTAT,KQ_Pus,KQ_Ius,icmaxpu,-icmaxpu)                    ! Us control
            icqref = -MIN(MAX(-us*(icq0 - ynq),qsminpu),qsmaxpu)/us
        ELSE
            CALL SUB_PI(ymq,mq,d_mq,(qsref-qs),3,DELTAT,KQ_Pqs,KQ_Iqs,icmaxpu,-icmaxpu)                    ! Qs control
            icqref = -qsref/us - ymq
        END IF

        ! Limit current references and set anti-wind up indicators.
        CALL SUB_IMAXLIMITSIREF(icdref, icqref, antiwindupicd, antiwindupicq, ILIMITPRIORITY, icmaxpu)
        CALL SUB_ECMAXLIMITSIREF(icdref, icqref, us, udc, FMODULATIONPWMMAX, zc)

        ! DC protection
        CALL SUB_OVDCPROT(istripvsc, counter_ovdc, tinitial_ovdc, udc, TIME, UDC_MAX, TUDCMAX, NUMBUS(IB), 'VSCGFL', LPDEV)
        CALL SUB_UVDCPROT(istripvsc, counter_uvdc, tinitial_uvdc, udc, TIME, UDC_MIN, TUDCMIN, NUMBUS(IB), 'VSCGFL', LPDEV)

        IF (istripvsc.GT.0) THEN
            icd = 0.0
            icq = 0.0
        END IF

        ! AC-side converter current, voltage, power
        ps = us*icd
        qs = -us*icq
        ic_phasor_RI = reference_transf*CMPLX(icd,icq)         ! dq -> DQ, pu in system rating	
        ec_phasor_RI = us_phasor_RI + zc*ic_phasor_RI          ! converter voltage -> output
        sc_phasor = ec_phasor_RI*CONJG(ic_phasor_RI)
        ic = ABS(ic_phasor_RI)
        icD_out = REAL(ic_phasor_RI)
        icQ_out = AIMAG(ic_phasor_RI)
        ec = ABS(ec_phasor_RI)
        deltac = ATAN2(AIMAG(ec_phasor_RI), REAL(ec_phasor_RI))
        pvsc = REAL(sc_phasor)
        qvsc = AIMAG(sc_phasor)
        fmodulationpwm = ec/MAX(udc,0.0001)

        ! DC-side converter current and power
        CALL SUB_COMPUTEPLOSS(ploss, ps, ic, alosspu, blosspu, clossinvpu, clossrectpu)
        pdc = -(pvsc + ploss)
        idc = pdc/udc             ! output of the converter model; input of the DC-grid model

        PMECH(I_MACH) = -pdc                      ! in pu system rating
        SPEED(I_MACH) = BSFREQ(IB)
		ANGLE(I_MACH) = deltac
        ETERM(I_MACH) = us
        PELEC(I_MACH) = ps                        ! in pu system rating
        QELEC(I_MACH) = -qs                   ! in pu system rating

    CASE (4)

        ! Update number of STATEs.
        ! ========================
        NINTEG = MAX(NINTEG,I_STATE+6)

    CASE (5)
        ! Reporting mode
        WRITE (LPDEV,*) 'Converter ', idx_converter, ' at bus ', NUMBUS(IB)
        WRITE (LPDEV,*) 'VSCGFL - ', 'CON', I_CON
        WRITE (LPDEV,*) 'VSCGFL - ', 'ICON', I_ICON
        WRITE (LPDEV,*) 'VSCGFL - ', 'VAR', I_VAR
        WRITE (LPDEV,*) 'VSCGFL - ', 'STATE', I_STATE

    CASE DEFAULT

    END SELECT

    ! RE-ASSIGN VARIABLES
    ! -------------------
    icd = xd
    icq = xq

    ! VARs
    VAR(I_VAR) = deltapsref      ! From SPWDRD
    VAR(I_VAR+1) = deltaqsref    ! From SQWDRD    
    VAR(I_VAR+2) = udcref        ! To SPWDRD
    VAR(I_VAR+3) = udc           ! To SPWDRD
    VAR(I_VAR+4) = pdc           ! To DDCGRD/SDCGRD
    
    VAR(I_VAR+5) = ps0           
    VAR(I_VAR+6) = qs0          
    VAR(I_VAR+7) = deltas0
    VAR(I_VAR+8) = icd0
    VAR(I_VAR+9) = icq0

    VAR(I_VAR+10) = usref
    VAR(I_VAR+11) = icdref
    VAR(I_VAR+12) = icqref
    
    VAR(I_VAR+13) = ps
    VAR(I_VAR+14) = qs
    VAR(I_VAR+15) = icd
    VAR(I_VAR+16) = icq

    VAR(I_VAR+17) = istripvsc
    VAR(I_VAR+18) = antiwindupicd
    VAR(I_VAR+19) = antiwindupicq
    VAR(I_VAR+20) = ic_abs_puconv
    VAR(I_VAR+21) = dcontroltype_aux
    VAR(I_VAR+22) = qcontroltype_aux

    VAR(I_VAR+23) = fmodulationpwm
    VAR(I_VAR+24) = counter_uvdc
    VAR(I_VAR+25) = tinitial_uvdc
    VAR(I_VAR+26) = counter_ovdc
    VAR(I_VAR+27) = tinitial_ovdc
    VAR(I_VAR+28) = idx_converter               

    ! STATEs
    STATE(I_STATE) = xd              ! icd
    STATE(I_STATE+1) = xq            ! icq
    STATE(I_STATE+2) = md            ! d-integral state variable type1
    STATE(I_STATE+3) = mq            ! q-integral state variable type1
    STATE(I_STATE+4) = ndc           ! d-integral state variable type2
    STATE(I_STATE+5) = nq            ! q-integral state variable type2
    STATE(I_STATE+6) = eta_d         ! d-integral state variable type3 (passive grid)

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
