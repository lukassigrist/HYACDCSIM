MODULE MOD_MISC

	IMPLICIT NONE

CONTAINS

    ! Subroutine to compute the power losses of the VSC operating as inverter or rectifier
    SUBROUTINE SUB_COMPUTEPLOSS(ploss, ps, ic, aloss, bloss, cinv, crect)
        
        REAL, INTENT(OUT) :: ploss
        REAL, INTENT(IN) :: ps, ic, aloss, bloss, cinv, crect

        IF (ps.GE.0.0) THEN 
            ploss = aloss + bloss*ic + cinv*ic**2  ! inverter
        ELSE 
            ploss = aloss + bloss*ic + crect*ic**2 ! rectifier
        END IF

    END SUBROUTINE SUB_COMPUTEPLOSS

    ! Subroutine to limit the current references according to the maximum current and maximum modulation index
    SUBROUTINE SUB_IMAXLIMITSIREF(icdref, icqref, antiwindupicd, antiwindupicq, ILIMITPRIORITY, icmax)

        INTEGER, INTENT(IN) :: ILIMITPRIORITY
        REAL, INTENT(IN) :: icmax
        
        REAL, INTENT(INOUT) :: icdref, icqref, antiwindupicd, antiwindupicq
    
        REAL :: icref      
    
        ! Check current limits (icref <= icmax)
        icref = SQRT(icdref**2 + icqref**2)
        IF (icref.GT.icmax) THEN
            IF (ILIMITPRIORITY.EQ.1) THEN 								! P-priority		  
                antiwindupicq = 0.0 									! icq will be always on the limit in this case
                IF (ABS(icdref).GT.icmax) THEN
                    antiwindupicd = 0.0 								! icd is on the limit
                ELSE
                    antiwindupicd = 1.0 										
                END IF
                icdref = MIN(ABS(icdref), icmax)*SIGN(1.0,icdref) 		
                icqref = SQRT(icmax**2 - icdref**2)*SIGN(1.0,icqref)		  
            ELSE IF (ILIMITPRIORITY.EQ.2) THEN 							! Q-priority
                antiwindupicd = 0.0 									! icq will be always on the limit in this case
                IF (ABS(icqref).GT.icmax) THEN
                    antiwindupicq = 0.0 								! icq is on the limit
                ELSE
                    antiwindupicq = 1.0 										
                END IF
                icqref = MIN(ABS(icqref), icmax)*SIGN(1.0,icqref) 		! change the current reference: icqref = icq_ref_aux
                icdref = SQRT(icmax**2 - icqref**2)*SIGN(1.0,icdref)		  
            ELSE 														! P-Q equal priority 
                antiwindupicd = 0.0
                antiwindupicq = 0.0
                icdref = icmax*icdref/icref
                icqref = icmax*icqref/icref	  
            END IF  
        ELSE
            antiwindupicd = 1.0
            antiwindupicq = 1.0
        END IF
        
    END SUBROUTINE SUB_IMAXLIMITSIREF

    ! Subroutine to limit the current references according to the maximum current and maximum modulation index
    SUBROUTINE SUB_ECMAXLIMITSIREF(icdref, icqref, us, udc, FMODULATIONPWMMAX, zc)

        REAL, INTENT(IN) :: us, udc, FMODULATIONPWMMAX
        COMPLEX, INTENT(IN) :: zc
        
        REAL, INTENT(INOUT) :: icdref, icqref
    
        REAL :: ecmax, ecref, deltacref
        COMPLEX :: icref_phasor_dq, ecref_phasor_dq, us_phasor_dq
            
        ! Check maximum modulation index (ec <= FMODULATIONPWMMAX*udc)
        ecmax = FMODULATIONPWMMAX*udc
        us_phasor_dq = CMPLX(us, 0.0)
        icref_phasor_dq = CMPLX(icdref, icqref) 
        ecref_phasor_dq = us_phasor_dq + zc*icref_phasor_dq
        ecref = ABS(ecref_phasor_dq)
        deltacref = ATAN2(AIMAG(ecref_phasor_dq), REAL(ecref_phasor_dq))
        IF (ecref.GT.ecmax) THEN	  
            ecref_phasor_dq = CMPLX(ecmax*cos(deltacref), ecmax*sin(deltacref)) 	  
            icref_phasor_dq = (ecref_phasor_dq - us_phasor_dq)/zc
            icdref = REAL(icref_phasor_dq)
            icqref = AIMAG(icref_phasor_dq) 
        END IF
        
    END SUBROUTINE SUB_ECMAXLIMITSIREF

    SUBROUTINE SUB_ZVILIMITSECREF(yecdviimax, yecqviimax, ic_phasor_dq, icmax, kpvi, sigmaxr)

        
        REAL, INTENT(OUT) :: yecdviimax, yecqviimax

        COMPLEX, INTENT(IN) :: ic_phasor_dq
        REAL, INTENT(IN) :: icmax, kpvi, sigmaxr

        REAL :: ic, deltaic
        COMPLEX :: zvi
        
        ic = abs(ic_phasor_dq)
        deltaic = MAX(ic-icmax,0.0)
        
        zvi = deltaic*kpvi*CMPLX(1,sigmaxr)
        yecdviimax = REAL(zvi*ic_phasor_dq)
        yecqviimax = AIMAG(zvi*ic_phasor_dq) 
        
    END SUBROUTINE SUB_ZVILIMITSECREF

    SUBROUTINE SUB_ICMAXLIMITSECREF(ec_phasor_dq, ic_phasor_dq, us_phasor_dq, ILIMITPRIORITY, icmax, zc)

        
        COMPLEX, INTENT(OUT) :: ec_phasor_dq

        COMPLEX, INTENT(IN) :: ic_phasor_dq, us_phasor_dq, zc
        REAL, INTENT(IN) :: icmax
        INTEGER, INTENT(IN) :: ILIMITPRIORITY

        REAL :: ic, phic, icd, icq, us, thetas, deltac
        
        us = ABS(us_phasor_dq)
        thetas = ATAN2(AIMAG(us_phasor_dq), REAL(us_phasor_dq))
        deltac = ATAN2(AIMAG(ec_phasor_dq), REAL(ec_phasor_dq))
        ic = MIN(abs(ic_phasor_dq),icmax)
        phic = ATAN2(AIMAG(ic_phasor_dq), REAL(ic_phasor_dq))   ! P-Q equal priority 
        icd = ic*cos(phic)
        icq = ic*sin(phic)
        IF (ILIMITPRIORITY.EQ.1) THEN                           ! P-priority
            icd = MIN(ABS(REAL(ic_phasor_dq)),icmax)*SIGN(1.0,icd)
            icq = min(SQRT(icmax**2 - icd**2),abs(icq))*SIGN(1.0,icq)
            phic = ATAN2(icq, icd)
        ELSE IF (ILIMITPRIORITY.EQ.2) THEN                      ! Q-priority
            icq = MIN(ABS(AIMAG(ic_phasor_dq)),icmax)*SIGN(1.0,icq)
            icd = min(SQRT(icmax**2 - icq**2),abs(icd))*SIGN(1.0,icd)
            phic = ATAN2(icq, icd)
        ELSE IF (ILIMITPRIORITY.EQ.3) THEN                      ! Modify |ec| with deltac = const
            phic = ASIN(-us*sin(thetas-deltac)/ABS(zc)/ic) + deltac - ATAN2(AIMAG(zc), REAL(zc))   
            icd = ic*cos(phic)
            icq = ic*sin(phic)                                                              
        END IF
       
                
        ec_phasor_dq = us_phasor_dq + zc*ic*CMPLX(cos(phic), sin(phic))
        
    END SUBROUTINE SUB_ICMAXLIMITSECREF

    SUBROUTINE SUB_IDQMAXPRIORITY(icdmax, icqmax, ic_phasor_dq, ILIMITPRIORITY, KICPRIORITY, icmax)

        
        REAL, INTENT(OUT) :: icdmax, icqmax
        COMPLEX, INTENT(IN) :: ic_phasor_dq
        REAL, INTENT(IN) :: KICPRIORITY, icmax
        INTEGER, INTENT(IN) :: ILIMITPRIORITY

        REAL :: icd, icq, phic
        
        icd = REAL(ic_phasor_dq)
        icq = AIMAG(ic_phasor_dq)
        IF (ILIMITPRIORITY.EQ.1) THEN                           ! P-priority
            icdmax = icmax*KICPRIORITY
            IF (ABS(icd).GT.icmax) THEN
                icqmax = 0.0
            ELSE
                icqmax = SQRT(icmax**2 - icd**2)
            END IF
        ELSE IF (ILIMITPRIORITY.EQ.2) THEN                      ! Q-priority
            icdmax = icmax*KICPRIORITY
            IF (ABS(icd).GT.icmax) THEN
                icqmax = 0.0
            ELSE
                icqmax = SQRT(icmax**2 - icd**2)
            END IF   
        ELSE                                                     ! P-Q equal priority
            phic = ATAN2(AIMAG(ic_phasor_dq), REAL(ic_phasor_dq)) 
            icdmax = icmax*cos(phic)
            icqmax = icmax*sin(phic)                                                             
        END IF

    END SUBROUTINE SUB_IDQMAXPRIORITY

END MODULE MOD_MISC