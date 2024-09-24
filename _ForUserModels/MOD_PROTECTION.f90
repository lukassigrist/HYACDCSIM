MODULE MOD_PROTECTION

	IMPLICIT NONE
	CONTAINS
	
	SUBROUTINE FUN_VDCPROT()
		CHARACTER(*), INTENT(IN) :: name
		PRINT *, 'Hello, ', name
	END SUBROUTINE VDCPROT
	
	! LVRT Capability    
    SUBROUTINE FUN_LVRT(busnumber,fin_sim,voltage_abs,VMN,VLVRT2,VLVRT1,timenow,&
		TMN,TLVRT2,TLVRT1,tinitial1,counter1,LPDEV) 
      
		REAL :: VMN,VLVRT2,VLVRT1,TMN,TLVRT2,TLVRT1
		REAL voltage_abs,timenow,tinitial1,aux_LVRT1,aux_LVRT2
		INTEGER busnumber,fin_sim,counter1,LPDEV

		IF (voltage_abs.LT.VLVRT2) THEN 

			IF (counter1 .EQ. 0.0) THEN
				tinitial1=timenow
				counter1=1
				WRITE (LPDEV,*) 'PVIMBL::LVRT: LVRT counter activated for PV inverter at ', busnumber
			END IF  

			IF (voltage_abs.LT.VMN) THEN 
				WRITE (LPDEV,*) 'PVIMBL::LVRT: PV inverter at ', busnumber, ' trips due to LVRT VMN: ', 'TIME ', timenow, 'TINIT ', tinitial1
				fin_sim=1
			END IF ! v<VMN

			IF (((timenow-tinitial1).GT.TMN) .AND. ((timenow-tinitial1).LE.TLVRT1)) THEN
				aux_LVRT1 = VMN+(VLVRT1-VMN)/(TLVRT1-TMN)*(timenow-tinitial1-TMN)
				IF (voltage_abs.LT.aux_LVRT1) THEN 
					WRITE (LPDEV,*) 'PVIMBL::LVRT: PV inverter at ', busnumber, ' trips due to LVRT VMN-VLVRT1: ', 'TIME ', timenow, 'TINIT ', tinitial1
				fin_sim=1
				END IF
			END IF ! t-tini1>=TMN & t-tini1<TLVRT1

			IF (((timenow-tinitial1) .GT. TLVRT1) .AND. ((timenow-tinitial1) .LE. TLVRT2)) THEN
				aux_LVRT2=VLVRT1
				IF (voltage_abs.LT.aux_LVRT2)  THEN 
					WRITE (LPDEV,*) 'PVIMBL::LVRT: PV inverter at ', busnumber, ' trips due to LVRT VLVRT1-VLVRT2: ', 'TIME ', timenow, 'TINIT ', tinitial1
					fin_sim=1
				END IF
			END IF ! t-tini1>=TLVRT1+TMN & t-tini1<TLVRT2

			IF ((timenow-tinitial1) .GT. TLVRT2) THEN
				WRITE (LPDEV,*) 'PVIMBL::LVRT: PV inverter at ', busnumber, ' trips due to LVRT VLVRT2: ', 'TIME ', timenow, 'TINIT ', tinitial1
				fin_sim=1
			END IF ! t-tini1>TLVRT2

		END IF ! v<VLVRT2           

		RETURN
    END SUBROUTINE FUN_LVRT
 
	! Undervoltage Protection     
    SUBROUTINE FUN_UVPROT(busnumber,fin_sim,voltage_abs,VUV2,VUV1,&
		timenow,TUV2,TUV1,tinitial1,counter1,LPDEV)
     
      REAL VUV2,VUV1,TUV2,TUV1
      REAL voltage_abs,timenow,tinitial1
      INTEGER busnumber,fin_sim,counter1,LPDEV
                  
      IF (voltage_abs.LT.VUV2) THEN  
      
      IF (counter1 .EQ. 0.0) THEN
      tinitial1=timenow
      counter1=1
      WRITE (LPDEV,*) 'PVIMBL::UVPROT: UVPROT counter activated for PV 
     *inverter at ', busnumber
      END IF 
      
c      IF (voltage_abs.LT.VMN) THEN
c      WRITE (LPDEV,*) 'PVIMBL::UVPROT: PV inverter at ', busnumber, ' trips
c     * due to undervoltage: ', 'TIME ', timenow, 'TINIT ', tinitial1
c      fin_sim=1
c      END IF ! v<VMN

      IF (((timenow-tinitial1) .GT. TUV1) .AND. ((timenow-tinitial1) .LE. TUV2)) THEN 
      IF (voltage_abs.LT.VUV1) THEN
      WRITE (LPDEV,*) 'PVIMBL::UVPROT: PV inverter at ', busnumber, ' trips
     * due to undervoltage: ', 'TIME ', timenow, 'TINIT ', tinitial1
      fin_sim=1
      END IF 
      END IF ! t-tini1>TUV1 & t-tini1<=TUV2
      
      IF ((timenow-tinitial1) .GT. TUV2) THEN 
      WRITE (LPDEV,*) 'PVIMBL::UVPROT: PV inverter at ', busnumber, ' trips
     * due to undervoltage: ', 'TIME ', timenow, 'TINIT ', tinitial1
      fin_sim=1
      END IF ! t-tini1>TUV2
      
      END IF ! v<VU2
             	              
      
      RETURN
      END SUBROUTINE FUN_UVPROT
      
C     Undervoltage Protection
C     =======================   
      SUBROUTINE FUN_OVPROT(busnumber,fin_sim,voltage_abs,VOV2,VOV1,
     *timenow,TOV2,TOV1,tinitial1,counter1,LPDEV)
     
      REAL VOV2,VOV1,TOV2,TOV1
      REAL voltage_abs,timenow,tinitial1
      INTEGER busnumber,fin_sim,counter1,LPDEV
      
      IF (voltage_abs.GT.VOV2) THEN 
      
      IF (counter1 .EQ. 0.0) THEN
      tinitial1=timenow
      counter1=1
      WRITE (LPDEV,*) 'PVIMBL::OVPROT: OVPROT counter activated for PV 
     *inverter at ', busnumber
      END IF    
          	
c      IF (voltage_abs.GT.VMX) THEN
c      WRITE (LPDEV,*) 'PVIMBL::CASE3: Voltage is higher than limit
c     * and PV inverter trips  ', busnumber, ': V= ', voltage_abs, 'VMAX= ', VMX
c      fin_sim=1
c      END IF
      
      IF (((timenow-tinitial1) .GT. TOV1) .AND. ((timenow-tinitial1) .LE. TOV2)) THEN
      IF (voltage_abs.GT.VOV1) THEN  
      WRITE (LPDEV,*) 'PVIMBL::OVPROT: PV inverter at ', busnumber, ' trips
     * due to overvoltage: ', 'TIME ', timenow, 'TINIT ', tinitial1
      fin_sim=1
      END IF 
      END IF ! t-tini1>TOV1 & t-tini1<=TOV2
         	              
      IF ((timenow-tinitial1) .GT. TOV2) THEN 
      WRITE (LPDEV,*) 'PVIMBL::OVPROT: PV inverter at ', busnumber, ' trips
     * due to overvoltage: ', 'TIME ', timenow, 'TINIT ', tinitial1
      fin_sim=1
      END IF ! t-tini1>TOV2
      
      END IF ! v>VO2
            	              

      
      RETURN
      END SUBROUTINE FUN_OVPROT
      
C     Underfrequency Protection
C     =========================       
      SUBROUTINE FUN_UFPROT(busnumber,fin_sim,freq,FUF2,FUF1,timenow,TUF2,
     *TUF1,tinitial1,counter1,LPDEV)
      
      REAL FUF2,FUF1,TUF2,TUF1
      REAL freq,timenow,tinitial1
      INTEGER busnumber,fin_sim,counter1,LPDEV
      
      IF (freq.LT.FUF2) THEN 
      
      IF (counter1 .EQ. 0.0) THEN
      tinitial1=timenow
      counter1=1
      WRITE (LPDEV,*) 'PVIMBL::UFPROT: UFPROT counter activated for PV 
     *inverter at ', busnumber
      END IF 
             	                   
      IF (((timenow-tinitial1) .GT. TUF1) .AND. ((timenow-tinitial1) .LE. TUF2)) THEN 
      IF (freq.LT.FUF1) THEN  
      WRITE (LPDEV,*) 'PVIMBL::UFPROT: PV inverter at ', busnumber, ' trips
     * due to underfrequency: ', 'TIME ', timenow, 'TINIT ', tinitial1
      fin_sim=1
      END IF 
      END IF ! t-tini1>TUF1 & t-tini1<=TUF2
      
      IF ((timenow-tinitial1) .GT. TUF2) THEN 
      WRITE (LPDEV,*) 'PVIMBL::UFPROT: PV inverter at ', busnumber, ' trips
     * due to underfrequency: ', 'TIME ', timenow, 'TINIT ', tinitial1
      fin_sim=1
      END IF ! t-tini1>TUF2
      
      END IF ! f<FUF2
      
      RETURN
      END SUBROUTINE FUN_UFPROT

C     Overfrequency Protection
C     =========================            
      SUBROUTINE FUN_OFPROT(busnumber,fin_sim,freq,FOF2,FOF1,timenow,TOF2,
     *TOF1,tinitial1,counter1,LPDEV)
      
      REAL FOF2,FOF1,TOF2,TOF1
      REAL freq,timenow,tinitial1
      INTEGER busnumber,fin_sim,counter1,LPDEV
      
      IF (freq.GT.FOF2) THEN 
      
      IF (counter1 .EQ. 0.0) THEN
      tinitial1=timenow
      counter1=1
      WRITE (LPDEV,*) 'PVIMBL::OFPROT: OFPROT counter activated for PV 
     *inverter at ', busnumber
      END IF     
         	              	              
      IF (((timenow-tinitial1).GT. TOF1) .AND. ((timenow-tinitial1) .LE. TOF2)) THEN 
      IF (freq.GT.FOF1) THEN 
      WRITE (LPDEV,*) 'PVIMBL::OFPROT: PV inverter at ', busnumber, ' trips
     * due to overfrequency: ', 'TIME ', timenow, 'TINIT ', tinitial1
      fin_sim=1
      END IF 
      END IF ! t-tini1>TOF1 & t-tini1<=TOF2 
      
      IF ((timenow-tinitial1) .GT. TOF2) THEN 
      WRITE (LPDEV,*) 'PVIMBL::OFPROT: PV inverter at ', busnumber, ' trips
     * due to overfrequency: ', 'TIME ', timenow, 'TINIT ', tinitial1
      fin_sim=1
      END IF ! t-tini1>TOF2
      
      END IF ! f>FOF2
      
      RETURN
      END SUBROUTINE FUN_OFPROT
	
END MODULE MOD_PROTECTION