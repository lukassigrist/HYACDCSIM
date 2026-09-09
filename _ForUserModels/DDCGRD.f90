! Dynamic model of the DC grid
!
! The DC grid makes use of the Pi-model of lines. The bus-end susceptances are combined into a single bus susceptance. The
! RL element forms the series branch. 
!
! A DC grid consists of NDCBUS DC buses and NDCLINES DC lines. Congruently, there are (NDCBUS + NDCLINES) of state variables,
! i.e., NDCBUS DC bus voltages (udc) and NDCLINES DC line currents (icc).
!
! x = A*x + B*u
! y = C*x
! 
! where
!
! x = [udc, icc], u = idc, y = x
! idc = pdc/udc, where pdc comes from the converter model
!
! The DC grid of a MTDC is modelled as a governor-type model. In case of several MTDCs, the DC grid model must account for
! the fact that a model instance cannot be called with different numbers of state and algebraic variables. Therefore, the 
! number of states variables are set when building the AC/DC model.
!
! To reflect this dynamic DC-grid representation and distinguish it from the static DC-grid model, the model is named
! DDCGRD.
	
MODULE MOD_DCGRID

    IMPLICIT NONE

	REAL, ALLOCATABLE :: A_GLOBAL(:,:), B_GLOBAL(:,:)
	REAL, ALLOCATABLE :: C_GLOBAL(:,:), D_GLOBAL(:,:)


	INTEGER, ALLOCATABLE :: I_VARCONV_GLOBAL(:,:)

	INTEGER :: aux_var_GLOBAL = 0

	INTEGER, ALLOCATABLE :: CONVERTER_ACDC_BUS_GLOBAL(:,:)           ! Global buffer to store converter-AC/DC bus data for all DC networks

    INTEGER, PARAMETER :: MAX_DCG = 9

    INTEGER :: DCGRID_REGISTERED(MAX_DCG) = 0

    INTEGER :: NDCBUS_OFFSET_GLOBAL(MAX_DCG) = 0
    INTEGER :: NDCLINES_OFFSET_GLOBAL(MAX_DCG) = 0

    INTEGER :: NDCBUS_SIZE_GLOBAL(MAX_DCG) = 0
    INTEGER :: NDCLINES_SIZE_GLOBAL(MAX_DCG) = 0

    INTEGER :: NEXT_NDCBUS_OFFSET_GLOBAL = 0
    INTEGER :: NEXT_NDCLINES_OFFSET_GLOBAL = 0

CONTAINS

    SUBROUTINE REGISTER_DCG(IDGRID, NDCBUS_IN, NDCLINES_IN, NDCBUS_PREVIOUS, NDCLINES_PREVIOUS, IERR_REG)

        CHARACTER(1), INTENT(IN) :: IDGRID
        INTEGER, INTENT(IN) :: NDCBUS_IN, NDCLINES_IN
        INTEGER, INTENT(OUT) :: NDCBUS_PREVIOUS, NDCLINES_PREVIOUS
        INTEGER, INTENT(OUT) :: IERR_REG

        INTEGER :: IDGRID_num

        IERR_REG = 0
        NDCBUS_PREVIOUS = 0
        NDCLINES_PREVIOUS = 0

        READ(IDGRID, *, ERR=900) IDGRID_num

        IF (IDGRID_num .LT. 1 .OR. IDGRID_num .GT. MAX_DCG) THEN
            IERR_REG = 1
            RETURN
        END IF

        IF (NDCBUS_IN .LT. 1 .OR. NDCLINES_IN .LT. 0) THEN
            IERR_REG = 4
            RETURN
        END IF

        IF (DCGRID_REGISTERED(IDGRID_num) .EQ. 0) THEN

            NDCBUS_OFFSET_GLOBAL(IDGRID_num) = NEXT_NDCBUS_OFFSET_GLOBAL
            NDCLINES_OFFSET_GLOBAL(IDGRID_num) = NEXT_NDCLINES_OFFSET_GLOBAL

            NDCBUS_SIZE_GLOBAL(IDGRID_num) = NDCBUS_IN
            NDCLINES_SIZE_GLOBAL(IDGRID_num) = NDCLINES_IN

            NEXT_NDCBUS_OFFSET_GLOBAL = NEXT_NDCBUS_OFFSET_GLOBAL + NDCBUS_IN
            NEXT_NDCLINES_OFFSET_GLOBAL = NEXT_NDCLINES_OFFSET_GLOBAL + NDCLINES_IN

            DCGRID_REGISTERED(IDGRID_num) = 1

        ELSE

            IF (NDCBUS_SIZE_GLOBAL(IDGRID_num) .NE. NDCBUS_IN .OR. &
                NDCLINES_SIZE_GLOBAL(IDGRID_num) .NE. NDCLINES_IN) THEN

                IERR_REG = 2
                RETURN

            END IF

        END IF

        NDCBUS_PREVIOUS = NDCBUS_OFFSET_GLOBAL(IDGRID_num)
        NDCLINES_PREVIOUS = NDCLINES_OFFSET_GLOBAL(IDGRID_num)

        RETURN

900     CONTINUE
        IERR_REG = 3
        RETURN

    END SUBROUTINE REGISTER_DCG


    SUBROUTINE ENSURE_GLOBAL_STORAGE(IERR_MEM)

        INTEGER, INTENT(OUT) :: IERR_MEM

        INTEGER :: NTOTDCBUS_AUTO
        INTEGER :: NTOTDCLINES_AUTO
        INTEGER :: NTOTSTATES_AUTO

        IERR_MEM = 0

        NTOTDCBUS_AUTO   = NEXT_NDCBUS_OFFSET_GLOBAL
        NTOTDCLINES_AUTO = NEXT_NDCLINES_OFFSET_GLOBAL
        NTOTSTATES_AUTO  = NTOTDCBUS_AUTO + NTOTDCLINES_AUTO

        IF (NTOTDCBUS_AUTO .LE. 0 .OR. NTOTSTATES_AUTO .LE. 0) THEN
            IERR_MEM = 10
            RETURN
        END IF

        CALL GROW_REAL2D(A_GLOBAL, NTOTSTATES_AUTO, NTOTSTATES_AUTO, IERR_MEM)
        IF (IERR_MEM .NE. 0) RETURN

        CALL GROW_REAL2D(B_GLOBAL, NTOTSTATES_AUTO, NTOTDCBUS_AUTO, IERR_MEM)
        IF (IERR_MEM .NE. 0) RETURN

        CALL GROW_REAL2D(C_GLOBAL, NTOTSTATES_AUTO, NTOTSTATES_AUTO, IERR_MEM)
        IF (IERR_MEM .NE. 0) RETURN

        CALL GROW_REAL2D(D_GLOBAL, NTOTSTATES_AUTO, NTOTDCBUS_AUTO, IERR_MEM)
        IF (IERR_MEM .NE. 0) RETURN

        CALL GROW_INT2D(I_VARCONV_GLOBAL, NTOTDCBUS_AUTO, 1, IERR_MEM)
        IF (IERR_MEM .NE. 0) RETURN

        CALL GROW_INT2D(CONVERTER_ACDC_BUS_GLOBAL, NTOTDCBUS_AUTO, 2, IERR_MEM)
        IF (IERR_MEM .NE. 0) RETURN

        aux_var_GLOBAL = 7

    END SUBROUTINE ENSURE_GLOBAL_STORAGE


    SUBROUTINE GROW_REAL2D(X, NROW_NEW, NCOL_NEW, IERR_MEM)

        REAL, ALLOCATABLE, INTENT(INOUT) :: X(:,:)
        INTEGER, INTENT(IN) :: NROW_NEW, NCOL_NEW
        INTEGER, INTENT(OUT) :: IERR_MEM

        REAL, ALLOCATABLE :: TMP(:,:)
        INTEGER :: NROW_OLD, NCOL_OLD
        INTEGER :: NROW_ALLOC, NCOL_ALLOC
        INTEGER :: STAT_ALLOC

        IERR_MEM = 0

        IF (NROW_NEW .LE. 0 .OR. NCOL_NEW .LE. 0) THEN
            IERR_MEM = 20
            RETURN
        END IF

        IF (.NOT. ALLOCATED(X)) THEN

            ALLOCATE(X(NROW_NEW, NCOL_NEW), STAT=STAT_ALLOC)

            IF (STAT_ALLOC .NE. 0) THEN
                IERR_MEM = 21
                RETURN
            END IF

            X = 0.0
            RETURN

        END IF

        NROW_OLD = SIZE(X,1)
        NCOL_OLD = SIZE(X,2)

        IF (NROW_NEW .LE. NROW_OLD .AND. NCOL_NEW .LE. NCOL_OLD) THEN
            RETURN
        END IF

        NROW_ALLOC = MAX(NROW_NEW, NROW_OLD)
        NCOL_ALLOC = MAX(NCOL_NEW, NCOL_OLD)

        ALLOCATE(TMP(NROW_ALLOC, NCOL_ALLOC), STAT=STAT_ALLOC)

        IF (STAT_ALLOC .NE. 0) THEN
            IERR_MEM = 22
            RETURN
        END IF

        TMP = 0.0
        TMP(1:NROW_OLD, 1:NCOL_OLD) = X(1:NROW_OLD, 1:NCOL_OLD)

        DEALLOCATE(X)

        ALLOCATE(X(NROW_ALLOC, NCOL_ALLOC), STAT=STAT_ALLOC)

        IF (STAT_ALLOC .NE. 0) THEN
            IERR_MEM = 23
            DEALLOCATE(TMP)
            RETURN
        END IF

        X = TMP

        DEALLOCATE(TMP)

    END SUBROUTINE GROW_REAL2D


    SUBROUTINE GROW_INT2D(X, NROW_NEW, NCOL_NEW, IERR_MEM)

        INTEGER, ALLOCATABLE, INTENT(INOUT) :: X(:,:)
        INTEGER, INTENT(IN) :: NROW_NEW, NCOL_NEW
        INTEGER, INTENT(OUT) :: IERR_MEM

        INTEGER, ALLOCATABLE :: TMP(:,:)
        INTEGER :: NROW_OLD, NCOL_OLD
        INTEGER :: NROW_ALLOC, NCOL_ALLOC
        INTEGER :: STAT_ALLOC

        IERR_MEM = 0

        IF (NROW_NEW .LE. 0 .OR. NCOL_NEW .LE. 0) THEN
            IERR_MEM = 30
            RETURN
        END IF

        IF (.NOT. ALLOCATED(X)) THEN

            ALLOCATE(X(NROW_NEW, NCOL_NEW), STAT=STAT_ALLOC)

            IF (STAT_ALLOC .NE. 0) THEN
                IERR_MEM = 31
                RETURN
            END IF

            X = 0
            RETURN

        END IF

        NROW_OLD = SIZE(X,1)
        NCOL_OLD = SIZE(X,2)

        IF (NROW_NEW .LE. NROW_OLD .AND. NCOL_NEW .LE. NCOL_OLD) THEN
            RETURN
        END IF

        NROW_ALLOC = MAX(NROW_NEW, NROW_OLD)
        NCOL_ALLOC = MAX(NCOL_NEW, NCOL_OLD)

        ALLOCATE(TMP(NROW_ALLOC, NCOL_ALLOC), STAT=STAT_ALLOC)

        IF (STAT_ALLOC .NE. 0) THEN
            IERR_MEM = 32
            RETURN
        END IF

        TMP = 0
        TMP(1:NROW_OLD, 1:NCOL_OLD) = X(1:NROW_OLD, 1:NCOL_OLD)

        DEALLOCATE(X)

        ALLOCATE(X(NROW_ALLOC, NCOL_ALLOC), STAT=STAT_ALLOC)

        IF (STAT_ALLOC .NE. 0) THEN
            IERR_MEM = 33
            DEALLOCATE(TMP)
            RETURN
        END IF

        X = TMP

        DEALLOCATE(TMP)

    END SUBROUTINE GROW_INT2D

END MODULE MOD_DCGRID


SUBROUTINE DDCGRD(I_MACH,I_SLOT)

	! for the global variables (matrices)
	USE MOD_DCGRID
	USE MOD_READHYADCSIM
	INCLUDE 'COMON4.INS'
	IMPLICIT none


	INTEGER I_SLOT, I_MACH
	INTEGER I_ICON, I_CON, I_STATE, I_VAR
	INTEGER NDCBUS, NDCLINES, NPOLES									! number of DC buses, lines, poles
	INTEGER NDCBUS_PREVIOUS, NDCLINES_PREVIOUS	! total number of DC buses and lines, and previously used number of DC buses and lines
	INTEGER, ALLOCATABLE :: I_VARCONV(:,:)							! array position of the converters
	INTEGER, ALLOCATABLE :: CONVERTER_ACDC_BUS(:,:)					! ac and dc buses of each converter									! vector with the VARs position index of each converter
	INTEGER i,j
    INTEGER IERR_REG, IERR_MEM

	CHARACTER(1) IDGRID
    REAL :: IDGRID_num
	
	REAL, ALLOCATABLE :: Ydc(:,:)										! admittance matrix
	REAL, ALLOCATABLE :: Zdc(:,:)										! Impedance matrix
	REAL, ALLOCATABLE :: Ac(:,:), Ac_T(:,:)								! incidence matrix
	REAL, ALLOCATABLE :: Gdc(:,:), Cdc(:,:)								! bus conductance and capacitor matrix (diagonal)
	REAL, ALLOCATABLE :: Rsdc(:,:), Ldc(:,:)							! line resistance and inductance matrix (diagonal)
	REAL, ALLOCATABLE :: Cdc_inv(:,:), Ldc_inv(:,:), Rsdc_inv(:,:)
	REAL, ALLOCATABLE :: A11(:,:), A12(:,:), A21(:,:), A22(:,:)
	REAL, ALLOCATABLE :: B1(:,:), B2(:,:)
	REAL, ALLOCATABLE :: A(:,:), B(:,:), C(:,:), D(:,:)
	REAL, ALLOCATABLE :: I_sts(:,:)

	REAL, ALLOCATABLE :: Idc_ini(:,:), Udc_ini(:,:), Pdc_ini(:,:)		! Initial DC-side currents, voltages, power of converters (NDCBUS x 1)
	REAL, ALLOCATABLE :: Idc(:,:), Udc(:,:), Pdc(:,:)					! DC-side currents, voltages, power of converters (NDCBUS x 1)
	REAL, ALLOCATABLE :: Icc_ini(:,:), Icc(:,:)						! Line currents (NDCLINES x 1)

	REAL, ALLOCATABLE :: x(:,:), d_x(:,:), u(:,:), y(:,:)				! x = [udc, icc], u = idc, y = x

	! Indexes
	I_CON=STRTIN(1,I_SLOT)
	I_STATE=STRTIN(2,I_SLOT)
	I_VAR=STRTIN(3,I_SLOT)
	I_ICON=STRTIN(4,I_SLOT)
	
	! ICONs assignment
	NDCBUS = ICON(I_ICON) 				! number of DC buses of the ith MTDC
	NDCLINES = ICON(I_ICON+1) 			! number of DC lines of the ith MTDC
	NPOLES = ICON(I_ICON+2) 			! number of poles of each converter
	IDGRID = CHRICN(I_ICON+3)
    
    CALL REGISTER_DCG(IDGRID, NDCBUS, NDCLINES, NDCBUS_PREVIOUS, NDCLINES_PREVIOUS, IERR_REG)

	IF (IERR_REG .NE. 0) THEN
		WRITE (LPDEV,*) 'DDCGRD - ERROR: automatic DC-grid offset registration failed. IDGRID = ', IDGRID, ' IERR_REG = ', IERR_REG
		RETURN
	END IF
    
	! Allocate matrices
	ALLOCATE(I_VARCONV(NDCBUS,1))
	ALLOCATE(CONVERTER_ACDC_BUS(NDCBUS,2))
	ALLOCATE(Ydc(NDCBUS,NDCBUS))
	ALLOCATE(Zdc(NDCBUS,NDCBUS))
	ALLOCATE(Ac(NDCBUS,NDCLINES))
	ALLOCATE(Ac_T(NDCLINES,NDCBUS))
	ALLOCATE(Gdc(NDCBUS,NDCBUS))
	ALLOCATE(Cdc(NDCBUS,NDCBUS))
	ALLOCATE(Rsdc(NDCLINES,NDCLINES))
	ALLOCATE(Ldc(NDCLINES,NDCLINES))
	ALLOCATE(Cdc_inv(NDCBUS,NDCBUS))
	ALLOCATE(Ldc_inv(NDCLINES,NDCLINES))
	ALLOCATE(Rsdc_inv(NDCLINES,NDCLINES))
	ALLOCATE(A11(NDCBUS,NDCBUS))
	ALLOCATE(A12(NDCBUS,NDCLINES))
	ALLOCATE(A21(NDCLINES,NDCBUS))
	ALLOCATE(A22(NDCLINES,NDCLINES))
	ALLOCATE(B1(NDCBUS,NDCBUS))
	ALLOCATE(B2(NDCLINES,NDCBUS))
	ALLOCATE(A(NDCBUS+NDCLINES,NDCBUS+NDCLINES))
	ALLOCATE(B(NDCBUS+NDCLINES,NDCBUS))
	ALLOCATE(C(NDCBUS+NDCLINES,NDCBUS+NDCLINES))
	ALLOCATE(D(NDCBUS+NDCLINES,NDCBUS))
	ALLOCATE(I_sts(NDCBUS+NDCLINES,NDCBUS+NDCLINES))

	ALLOCATE(Idc_ini(NDCBUS,1))
	ALLOCATE(Udc_ini(NDCBUS,1))
	ALLOCATE(Pdc_ini(NDCBUS,1))
	ALLOCATE(Icc_ini(NDCLINES,1))

	ALLOCATE(Idc(NDCBUS,1))
	ALLOCATE(Udc(NDCBUS,1))
	ALLOCATE(Pdc(NDCBUS,1))
	ALLOCATE(Icc(NDCLINES,1))

	ALLOCATE(x(NDCBUS+NDCLINES,1))
	ALLOCATE(d_x(NDCBUS+NDCLINES,1))
	ALLOCATE(u(NDCBUS,1))
	ALLOCATE(y(NDCBUS+NDCLINES,1))

	! STATEs, DSTATEs and VARs assignment
	DO i=1,(NDCBUS+NDCLINES) 
		x(i,1) = STATE(I_STATE + i - 1)
		d_x(i,1) = DSTATE(I_STATE + i - 1)
	  
		IF (i.LE.NDCBUS) THEN 
			Udc(i,1) = x(i,1) 						! DC voltage
			Pdc(i,1) = VAR(I_VAR + i - 1) 			! DC power of converters (stored as VARs)
			Idc(i,1) = VAR(I_VAR + NDCBUS - 1 + i ) ! DC current of converters (stored as VARs)
		ELSE
			Icc(i-NDCBUS,1) = x(i,1) 				! DC branch currents
		END IF
	END DO

	! set the matrices to zero first
	CALL SUB_JZEROS(Ydc, NDCBUS, NDCBUS)
	CALL SUB_JZEROS(Zdc, NDCBUS, NDCBUS)
	CALL SUB_JZEROS(Gdc, NDCBUS, NDCBUS)
	CALL SUB_JZEROS(Cdc, NDCBUS, NDCBUS)
	CALL SUB_JZEROS(Cdc_inv, NDCBUS, NDCBUS)
	CALL SUB_JZEROS(Rsdc, NDCLINES, NDCLINES)
	CALL SUB_JZEROS(Rsdc_inv, NDCLINES, NDCLINES)
	CALL SUB_JZEROS(Ldc, NDCLINES, NDCLINES)
	CALL SUB_JZEROS(Ldc_inv, NDCLINES, NDCLINES)
	CALL SUB_JZEROS(Ac, NDCBUS, NDCLINES)
    
	SELECT CASE (MODE)
	  
	CASE (1) ! INITIAL STATE

		CALL ENSURE_GLOBAL_STORAGE(IERR_MEM)

		IF (IERR_MEM .NE. 0) THEN
			WRITE (LPDEV,*) 'DDCGRD - ERROR: global storage allocation failed. IDGRID = ', IDGRID, ' IERR_MEM = ', IERR_MEM
			RETURN
		END IF
	  
		! Read DC-grid data from .txt file
		CALL SUB_READDCGRID(I_VARCONV, Udc_ini, Ydc, Zdc, Ac, Gdc, Cdc_inv, Rsdc, Rsdc_inv, Ldc_inv, I_MACH, CONVERTER_ACDC_BUS, IDGRID, NDCBUS, NDCLINES)
        DO i = 1, NDCBUS
    		CONVERTER_ACDC_BUS_GLOBAL(NDCBUS_PREVIOUS + i, 1) = CONVERTER_ACDC_BUS(i,1)
    		CONVERTER_ACDC_BUS_GLOBAL(NDCBUS_PREVIOUS + i, 2) = CONVERTER_ACDC_BUS(i,2)
		END DO
       
		Ac_T = TRANSPOSE(Ac)

		! identity matrix: (NDCBUS+NDCLINES)x(NDCBUS+NDCLINES)
		CALL SUB_JZEROS(I_sts, NDCBUS+NDCLINES, NDCBUS+NDCLINES)
		DO i=1,(NDCBUS+NDCLINES)
			I_sts(i,i) = 1.0   
		END DO

		! Build state-space matrices
		! --------------------------	
		! 	dx/dt = Ax + Bu
		! 	y = Cx+Du
		! where A = [A11 A12; A21 A22] and B = [B1;B2]
		A11 = -MATMUL(Cdc_inv, Gdc) ! nxn
		A12 = -MATMUL(Cdc_inv, Ac) ! nxnL
		A21 = MATMUL(Ldc_inv, Ac_T) ! nLxn
		A22 = -MATMUL(Ldc_inv, Rsdc) ! nLxnL

		DO i=1,(NDCBUS+NDCLINES)
			DO j=1,(NDCBUS+NDCLINES)
				IF (i.LE.NDCBUS) THEN
					IF (j.LE.NDCBUS) THEN
						A(i,j) = A11(i,j)
					ELSE ! j>NDCBUS
						A(i,j) = A12(i,j-NDCBUS)
					END IF
				ELSE ! i>NDCBUS
					IF (j.LE.NDCBUS) THEN
						A(i,j) = A21(i-NDCBUS,j)
					ELSE ! j>NDCBUS
						A(i,j) = A22(i-NDCBUS,j-NDCBUS)
					END IF
				END IF      
			END DO
		END DO

		B1 = Cdc_inv ! nxn
		CALL SUB_JZEROS(B2, NDCLINES, NDCBUS)

		DO i=1,(NDCBUS+NDCLINES)
			DO j=1,NDCBUS
				IF (i.LE.NDCBUS) THEN
					B(i,j) = B1(i,j) 
				ELSE ! i>NDCBUS
					B(i,j) = B2(i-NDCBUS,j) 
				END IF
			END DO
		END DO

		C = I_sts ! identity: (NDCBUS+NDCLINES)x(NDCBUS+NDCLINES) 
		CALL SUB_JZEROS(D, NDCBUS+NDCLINES, NDCBUS)

		DO i=1,NDCBUS
		    IF (CONVERTER_ACDC_BUS(i,2) .NE. -1.0) THEN
		       Pdc_ini(i,1) = VAR(I_VARCONV(i,1) + 4) 			! DC-power from each converter model 
		    ELSE
		       Pdc_ini(i,1) = 0.0   ! No converter on this DC bus
		    END IF   
		END DO
		
		CALL SUB_DIVTBT(Pdc_ini, Udc_ini, Idc_ini, NDCBUS, 1) 
		Icc_ini = MATMUL(Rsdc_inv, MATMUL(Ac_T,Udc_ini))  		! Icc_ini = Rsdc_inv*(Ac'*Udc_ini))
		
		! store global matrices
		CALL SUB_LOCALTOGLOBAL(A, A_GLOBAL, NDCBUS+NDCLINES, NDCBUS+NDCLINES, SIZE(A_GLOBAL,1), SIZE(A_GLOBAL,2), NDCBUS_PREVIOUS + NDCLINES_PREVIOUS, NDCBUS_PREVIOUS + NDCLINES_PREVIOUS)
		CALL SUB_LOCALTOGLOBAL(B, B_GLOBAL, NDCBUS+NDCLINES, NDCBUS, SIZE(B_GLOBAL,1), SIZE(B_GLOBAL,2), NDCBUS_PREVIOUS + NDCLINES_PREVIOUS, NDCBUS_PREVIOUS)
		CALL SUB_LOCALTOGLOBAL(C, C_GLOBAL, NDCBUS+NDCLINES, NDCBUS+NDCLINES, SIZE(C_GLOBAL,1), SIZE(C_GLOBAL,2), NDCBUS_PREVIOUS + NDCLINES_PREVIOUS, NDCBUS_PREVIOUS + NDCLINES_PREVIOUS)
        
		DO i = 1,NDCBUS
			I_VARCONV_GLOBAL(NDCBUS_PREVIOUS+i,1) = I_VARCONV(i,1)
		END DO
		
		! Initialization
		! --------------
		Udc = Udc_ini
		Idc = Idc_ini
		Pdc = Pdc_ini
		Icc = Icc_ini
		  
		! states
		DO i=1,NDCBUS
			x(i,1) = Udc(i,1)
		END DO
		  
		DO i=1,NDCLINES
			x(NDCBUS+i,1) = Icc(i,1)
		END DO
		  
		! inputs
		u = Idc
		  
		! outputs
		y = MATMUL(C,x) + MATMUL(D,u) ! y = x
		  
		! report mode 
		WRITE (LPDEV,*) 'DDCGRD - CASE 1: DC grid at AC bus ',NUMBUS(NUMTRM(I_MACH)),' with id ',MACHID(I_MACH),' initialized.'
 
	CASE (2)

		! Compute derivatives
		! ===================
		! u = idc = pdc/uc (where pdc comes from converter) 
		! x = [udc;icc]: states
		! y = x

		CALL SUB_GLOBALTOLOCAL(A, A_GLOBAL, NDCBUS+NDCLINES, NDCBUS+NDCLINES, SIZE(A_GLOBAL,1), SIZE(A_GLOBAL,2), NDCBUS_PREVIOUS + NDCLINES_PREVIOUS, NDCBUS_PREVIOUS + NDCLINES_PREVIOUS)
		CALL SUB_GLOBALTOLOCAL(B, B_GLOBAL, NDCBUS+NDCLINES, NDCBUS, SIZE(B_GLOBAL,1), SIZE(B_GLOBAL,2), NDCBUS_PREVIOUS + NDCLINES_PREVIOUS, NDCBUS_PREVIOUS)
		CALL SUB_JZEROS(D, NDCBUS+NDCLINES, NDCBUS)

		u = Idc  
		
		d_x = MATMUL(A,x) + MATMUL(B,u) 		! derivatives

  	
	CASE (3) 

		! Outputs
		! =======
		CALL SUB_GLOBALTOLOCAL(C, C_GLOBAL, NDCBUS+NDCLINES, NDCBUS+NDCLINES, SIZE(C_GLOBAL,1), SIZE(C_GLOBAL,2), NDCBUS_PREVIOUS + NDCLINES_PREVIOUS, NDCBUS_PREVIOUS + NDCLINES_PREVIOUS)
		CALL SUB_JZEROS(D, NDCBUS+NDCLINES, NDCBUS)
	   
		DO i=1,NDCBUS
			I_VARCONV(i,1) = I_VARCONV_GLOBAL(NDCBUS_PREVIOUS+i,1)
		END DO

		DO i = 1, NDCBUS
			CONVERTER_ACDC_BUS(i,1) = CONVERTER_ACDC_BUS_GLOBAL(NDCBUS_PREVIOUS + i, 1)
			CONVERTER_ACDC_BUS(i,2) = CONVERTER_ACDC_BUS_GLOBAL(NDCBUS_PREVIOUS + i, 2)
		END DO
		DO i=1,NDCBUS
		   IF (CONVERTER_ACDC_BUS(i,2) .NE. -1.0) THEN
			   Pdc(i,1) = VAR(I_VARCONV(i,1) + 4)
			ELSE
			   Pdc(i,1) = 0.0
			END IF      
		END DO
		
		CALL SUB_DIVTBT(Pdc, Udc, Idc, NDCBUS, 1) ! Idc = Pdc/Udc
		u = Idc
		
		y = MATMUL(C,x) + MATMUL(D,u) 			! outputs y = x

	CASE (4)
		NINTEG = MAX(NINTEG, I_STATE+NDCBUS+NDCLINES-1)
	  
	CASE (5)
		! reporting mode

		WRITE (LPDEV,*) 'DDCGRD - ', 'CON', I_CON
		WRITE (LPDEV,*) 'DDCGRD - ', 'ICON', I_ICON
		WRITE (LPDEV,*) 'DDCGRD - ', 'VAR', I_VAR
		WRITE (LPDEV,*) 'DDCGRD - ', 'STATE', I_STATE

	CASE DEFAULT
		  
	END SELECT

	! STATEs and DSTATEs assignment
	DO i=1,(NDCBUS+NDCLINES)
		STATE(I_STATE+i-1) = x(i,1)
		DSTATE(I_STATE+i-1) = d_x(i,1)
	END DO

	! VARs assignment
	DO i=1,NDCBUS
		VAR(I_VAR + i - 1) = Pdc(i,1) 
		VAR(I_VAR + NDCBUS - 1 + i ) = Idc(i,1)
	END DO

END

SUBROUTINE SUB_JZEROS(A, NROW, NCOL)
	! subroutine to create a zero matrix of size mxn 
	! A = zeros(NROW,NCOL)
	INTEGER :: NROW,NCOL
	REAL :: A(NROW,NCOL)
	INTEGER :: i,j

	DO i=1,NROW
		DO j=1,NCOL
			A(i,j) = 0.0
		END DO
	END DO
END

SUBROUTINE SUB_JONES(A, NROW, NCOL)
	! subroutine to create a zero matrix of size mxn 
	! A = ones(NROW,NCOL)
	INTEGER :: NROW,NCOL
	REAL :: A(NROW,NCOL)
	INTEGER :: i,j

	DO i=1,NROW
		DO j=1,NCOL
			A(i,j) = 1.0
		END DO
	END DO
END

SUBROUTINE SUB_MULTBT(A, B, C, NROW, NCOL)
	! subroutine to obtain the term-by-term
	! multiplication of two matrices mxn: C = A.*B
	INTEGER :: NROW,NCOL
	REAL :: A(NROW,NCOL), B(NROW,NCOL), C(NROW,NCOL)
	INTEGER :: i,j
	
	DO i=1,NROW
		DO j=1,NCOL
			C(i,j) = A(i,j)*B(i,j)
		END DO
	END DO
END

SUBROUTINE SUB_DIVTBT(A, B, C, NROW, NCOL)
	! subroutine to obtain the term-by-term
	! division of two matrices mxn: C = A./B
	INTEGER :: NROW,NCOL
	REAL :: A(NROW,NCOL), B(NROW,NCOL), C(NROW,NCOL)
	INTEGER :: i,j
	
	DO i=1,NROW
		DO j=1,NCOL
			C(i,j) = A(i,j)/B(i,j)
		END DO
	END DO
END

SUBROUTINE SUB_GLOBALTOLOCAL(X, X_GLOBAL, NROW_LOCAL, NCOL_LOCAL, NROW_GLOBAL, NCOL_GLOBAL, NROW_PREVIOUS, NCOL_PREVIOUS)
	! subroutine to obtain the current grid A,B,C,D matrix from the global matrices
	! (NROW_LOCAL,NCOL_LOCAL) Local X dimensions, (NROW_GLOBAL, NCOL_GLOBAL) Global X dimensions, (NROW_PREVIOUS,NCOL_PREVIOUS) pointer of Global matrix where Local matrix is saved
	INTEGER NROW_LOCAL, NCOL_LOCAL, NROW_GLOBAL, NCOL_GLOBAL, NROW_PREVIOUS, NCOL_PREVIOUS
	REAL X(NROW_LOCAL,NCOL_LOCAL), X_GLOBAL(NROW_GLOBAL,NCOL_GLOBAL)
	INTEGER  i, j

	DO i = 1, NROW_LOCAL
		DO j = 1, NCOL_LOCAL
			X(i,j) = X_GLOBAL(NROW_PREVIOUS + i, NCOL_PREVIOUS + j)
		END DO
	END DO
END      

SUBROUTINE SUB_LOCALTOGLOBAL(X, X_GLOBAL, NROW_LOCAL, NCOL_LOCAL, NROW_GLOBAL, NCOL_GLOBAL, NROW_PREVIOUS, NCOL_PREVIOUS)
	! subroutine to save the current grid A,B,C,D matrix into the global matrices
	! (NROW_LOCAL,NCOL_LOCAL) Local X dimensions, (NROW_GLOBAL, NCOL_GLOBAL) Global X dimensions, (NROW_PREVIOUS,NCOL_PREVIOUS) pointer of Global matrix where Local matrix is saved
	INTEGER :: NROW_LOCAL, NCOL_LOCAL, NROW_GLOBAL, NCOL_GLOBAL, NROW_PREVIOUS, NCOL_PREVIOUS
	REAL :: X(NROW_LOCAL,NCOL_LOCAL), X_GLOBAL(NROW_GLOBAL,NCOL_GLOBAL)
	INTEGER :: i, j

	DO i = 1, NROW_LOCAL
		DO j = 1, NCOL_LOCAL
			X_GLOBAL(NROW_PREVIOUS + i, NCOL_PREVIOUS + j) = X(i,j)
		END DO
	END DO
END
