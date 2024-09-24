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
	
MODULE module_globalmatrices
	REAL, ALLOCATABLE :: A_GLOBAL(:,:), B_GLOBAL(:,:)
	REAL, ALLOCATABLE :: C_GLOBAL(:,:), D_GLOBAL(:,:)
	REAL, ALLOCATABLE :: I_VAR_SVSCON_GLOBAL(:,:)
	INTEGER aux_var_GLOBAL
END MODULE module_globalmatrices


SUBROUTINE DCGRID(I_MACH,I_SLOT)

	! for the global variables (matrices)
	USE module_globalmatrices
	INCLUDE 'COMON4.INS'
	IMPLICIT none


	INTEGER I_SLOT, I_MACH
	INTEGER I_ICON, I_CON, I_STATE, I_VAR
	INTEGER NDCBUS, NDCLINES, NPOLES									! number of DC buses, lines, poles
	INTEGER NTOTDCBUS, NTOTDCLINES, NDCBUS_PREVIOUS, NDCLINES_PREVIOUS	! total number of DC buses and lines, and previously used number of DC buses and lines
	INTEGER, ALLOCATABLE :: I_VAR_SVSCON(:,:) 							! array position of the converters
	INTEGER, ALLOCATABLE :: CONVERTER_ACDC_BUS(:,:) 					! ac and dc buses of each converter													! vector with the VARs position index of each converter
	INTEGER i,j

	CHARACTER(1) IDGRID

	REAL UDC_BASE_kV, PDC_BASE_MW
	
	REAL, ALLOCATABLE :: Ydc(:,:) 										! admittance matrix
	REAL, ALLOCATABLE :: Ac(:,:), Ac_T(:,:)								! incidence matrix
	REAL, ALLOCATABLE :: Gdc(:,:), Cdc(:,:)								! bus conductance and capacitor matrix (diagonal)
	REAL, ALLOCATABLE :: Rsdc(:,:), Ldc(:,:)							! line resistance and inductance matrix (diagonal)
	REAL, ALLOCATABLE :: Cdc_inv(:,:), Ldc_inv(:,:), Rsdc_inv(:,:)
	REAL, ALLOCATABLE :: A11(:,:), A12(:,:), A21(:,:), A22(:,:)
	REAL, ALLOCATABLE :: B1(:,:), B2(:,:)
	REAL, ALLOCATABLE :: A(:,:), B(:,:), C(:,:), D(:,:)
	REAL, ALLOCATABLE :: I_sts(:,:)

	REAL, ALLOCATABLE :: Idc_ini(:,:), Udc_ini(:,:), Pdc_ini(:,:) 		! Initial DC-side currents, voltages, power of converters (NDCBUS x 1)
	REAL, ALLOCATABLE :: Idc(:,:), Udc(:,:), Pdc(:,:) 					! DC-side currents, voltages, power of converters (NDCBUS x 1)
	REAL, ALLOCATABLE :: Icc_ini(:,:), Icc(:,:) 						! Line currents (NDCLINES x 1)

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
	NTOTDCBUS = ICON(I_ICON+4)			! total number of DC buses of all MTDCs
	NTOTDCLINES = ICON(I_ICON+5)		! total number of DC lines of all MTDCs
	NDCBUS_PREVIOUS = ICON(I_ICON+6)	! number of DC buses of the previous (ith-1) MTDC if any
	NDCLINES_PREVIOUS = ICON(I_ICON+7)	! number of DC lines of the previous (ith-1) MTDC if any

	! CONs assignment
	UDC_BASE_kV = CON(I_CON) 			! DC base voltage (kV)
	PDC_BASE_MW = CON(I_CON+1) 			! DC base power (MW)

	! Allocate matrices
	ALLOCATE(I_VAR_SVSCON(NDCBUS,1))
	ALLOCATE(CONVERTER_ACDC_BUS(NDCBUS,2))
	ALLOCATE(Ydc(NDCBUS,NDCBUS))
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

		IF (aux_var_global.NE.7) THEN
			ALLOCATE(A_GLOBAL(NTOTDCBUS+NTOTDCLINES,NTOTDCBUS+NTOTDCLINES)) ! size of total number of state variables
			ALLOCATE(B_GLOBAL(NTOTDCBUS+NTOTDCLINES,NTOTDCBUS))
			ALLOCATE(C_GLOBAL(NTOTDCBUS+NTOTDCLINES,NTOTDCBUS+NTOTDCLINES))
			ALLOCATE(D_GLOBAL(NDCBUS+NDCLINES,NDCBUS))
			ALLOCATE(I_VAR_SVSCON_GLOBAL(NTOTDCBUS,1))
			aux_var_global = 7	
		END IF
	  
		! Read DC-grid data from .txt file
		CALL SUB_READFROMFILEDCGRID(I_VAR_SVSCON, Udc_ini, Ydc, Ac, Gdc, Cdc_inv, Rsdc, Rsdc_inv, Ldc_inv, I_MACH, CONVERTER_ACDC_BUS, IDGRID, NDCBUS, NDCLINES)	

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

		!Idc_ini = MATMUL(Ydc, Udc_ini) 						! Idc_ini = Ydc*Udc_ini where Udc_ini from file
		!CALL SUB_MULTBT(Udc_ini, Idc_ini, Pdc_ini, NDCBUS, 1) 	! Pdc_ini = Udc_ini.*Idc_ini
		DO i=1,NDCBUS
			Pdc_ini(i,1) = VAR(I_VAR_SVSCON(i,1) + 25) 			! DC-power from each converter model 
		END DO
		CALL SUB_DIVTBT(Pdc_ini, Udc_ini, Idc_ini, NDCBUS, 1) 
		Icc_ini = MATMUL(Rsdc_inv, MATMUL(Ac_T,Udc_ini))  		! Icc_ini = Rsdc_inv*(Ac'*Udc_ini))
		
		! store global matrices
		CALL SUB_LOCALTOGLOBAL(A, A_GLOBAL, NDCBUS+NDCLINES, NDCBUS+NDCLINES, NTOTDCBUS + NTOTDCLINES, NTOTDCBUS + NTOTDCLINES, NDCBUS_PREVIOUS + NDCLINES_PREVIOUS, NDCBUS_PREVIOUS + NDCLINES_PREVIOUS)
		CALL SUB_LOCALTOGLOBAL(B, B_GLOBAL, NDCBUS+NDCLINES, NDCBUS, NTOTDCBUS + NTOTDCLINES, NTOTDCBUS, NDCBUS_PREVIOUS + NDCLINES_PREVIOUS, NDCBUS_PREVIOUS)
		CALL SUB_LOCALTOGLOBAL(C, C_GLOBAL, NDCBUS+NDCLINES, NDCBUS+NDCLINES, NTOTDCBUS + NTOTDCLINES, NTOTDCBUS + NTOTDCLINES, NDCBUS_PREVIOUS + NDCLINES_PREVIOUS, NDCBUS_PREVIOUS + NDCLINES_PREVIOUS)

		DO i = 1,NDCBUS
			I_VAR_SVSCON_GLOBAL(NDCBUS_PREVIOUS+i,1) = I_VAR_SVSCON(i,1)
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
		WRITE (LPDEV,*) 'DCGRID - CASE1 1: DC grid at AC bus ',NUMBUS(NUMTRM(I_MACH)),' with id ',MACHID(I_MACH),' initialized.'
		! WRITE (LPDEV,*) 'DCGRID - CASE 1: Pdc_ini = ',Pdc
		! WRITE (LPDEV,*) 'DCGRID - CASE 1: Udc_ini = ',Udc
		! WRITE (LPDEV,*) 'DCGRID - CASE 1: Idc_ini = ',Idc
		! WRITE (LPDEV,*) 'DCGRID - CASE 1: Icc_ini = ',Icc

		  
	CASE (2)

		! Compute derivatives
		! ===================
		! u = idc = pdc/uc (where pdc comes from converter) 
		! x = [udc;icc]: states
		! y = x

		CALL SUB_GLOBALTOLOCAL(A, A_GLOBAL, NDCBUS+NDCLINES, NDCBUS+NDCLINES, NTOTDCBUS + NTOTDCLINES, NTOTDCBUS + NTOTDCLINES, NDCBUS_PREVIOUS + NDCLINES_PREVIOUS, NDCBUS_PREVIOUS + NDCLINES_PREVIOUS)
		CALL SUB_GLOBALTOLOCAL(B, B_GLOBAL, NDCBUS+NDCLINES, NDCBUS, NTOTDCBUS + NTOTDCLINES, NTOTDCBUS, NDCBUS_PREVIOUS + NDCLINES_PREVIOUS, NDCBUS_PREVIOUS)
		CALL SUB_JZEROS(D, NDCBUS+NDCLINES, NDCBUS)

		u = Idc  
		d_x = MATMUL(A,x) + MATMUL(B,u) 		! derivatives
	  	  
	CASE (3) 
		! Outputs
		! =======
		CALL SUB_GLOBALTOLOCAL(C, C_GLOBAL, NDCBUS+NDCLINES, NDCBUS+NDCLINES, NTOTDCBUS + NTOTDCLINES, NTOTDCBUS + NTOTDCLINES, NDCBUS_PREVIOUS + NDCLINES_PREVIOUS, NDCBUS_PREVIOUS + NDCLINES_PREVIOUS)
		CALL SUB_JZEROS(D, NDCBUS+NDCLINES, NDCBUS)
	  
		DO i=1,NDCBUS
			I_VAR_SVSCON(i,1) = I_VAR_SVSCON_GLOBAL(NDCBUS_PREVIOUS+i,1)
		END DO
		DO i=1,NDCBUS
			Pdc(i,1) = VAR(I_VAR_SVSCON(i,1) + 25) 
		END DO
		CALL SUB_DIVTBT(Pdc, Udc, Idc, NDCBUS, 1) ! Idc = Pdc/Udc
		u = Idc
		
		y = MATMUL(C,x) + MATMUL(D,u) 			! outputs y = x
	  	  
	CASE (4)
		NINTEG = MAX(NINTEG, I_STATE+NDCBUS+NDCLINES-1)
	  
	CASE (5)
		! reporting mode

		WRITE (LPDEV,*) 'DCGRID - ', 'CON', I_CON
		WRITE (LPDEV,*) 'DCGRID - ', 'ICON', I_ICON
		WRITE (LPDEV,*) 'DCGRID - ', 'VAR', I_VAR
		WRITE (LPDEV,*) 'DCGRID - ', 'STATE', I_STATE

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

SUBROUTINE SUB_READFROMFILEDCGRID(I_VAR_SVSCON, Udc_ini, Ydc, Ac, Gdc, Cdc_inv, Rsdc, Rsdc_inv, Ldc_inv, I_MACH, CONVERTER_ACDC_BUS, IDGRID, NDCBUS, NDCLINES)

	INCLUDE 'COMON4.INS'
	IMPLICIT none
	
	CHARACTER(1) :: IDGRID, line_header 
	INTEGER :: I_VAR_SVSCON(NDCBUS,1), CONVERTER_ACDC_BUS(NDCBUS,2), i, j, NDCBUS, NDCLINES
	INTEGER :: IERR, IVAL, I_MACH
	REAL :: Ydc(NDCBUS,NDCBUS), Ac(NDCBUS,NDCLINES), Gdc(NDCBUS,NDCBUS), Cdc(NDCBUS,NDCBUS), Rsdc(NDCLINES,NDCLINES)
	REAL :: Ldc(NDCLINES,NDCLINES), Cdc_inv(NDCBUS,NDCBUS), Ldc_inv(NDCLINES,NDCLINES), Rsdc_inv(NDCLINES,NDCLINES)
	REAL :: Udc_ini(NDCBUS,1)

	OPEN(UNIT=10, FILE='.\data_Udc_ini.txt')
	OPEN(UNIT=11, FILE='.\data_Ydc.txt')
	OPEN(UNIT=12, FILE='.\data_Ac.txt')
	OPEN(UNIT=13, FILE='.\data_Gdc.txt')
	OPEN(UNIT=14, FILE='.\data_Cdc.txt')
	OPEN(UNIT=15, FILE='.\data_Rsdc.txt')
	OPEN(UNIT=16, FILE='.\data_Ldc.txt')
	OPEN(UNIT=20, FILE='.\data_acdcbus.txt')


	! Salta las líneas hasta encontrar la letra de su RED [header]
	line_header = "+"
	DO WHILE (line_header.NE.IDGRID)
		READ(10, *) line_header
	END DO
	
	line_header = "+"
	DO WHILE (line_header.NE.IDGRID)
		READ(11, *) line_header
	END DO
	
	line_header = "+"
	DO WHILE (line_header.NE.IDGRID)
		READ(12, *) line_header
	END DO
	
	line_header = "+"
	DO WHILE (line_header.NE.IDGRID)
		READ(13, *) line_header
	END DO
	
	line_header = "+"
	DO WHILE (line_header.NE.IDGRID)
		READ(14, *) line_header
	END DO
	
	line_header = "+"
	DO WHILE (line_header.NE.IDGRID)
		READ(15, *) line_header
	END DO
	
	line_header = "+"
	DO WHILE (line_header.NE.IDGRID)
		READ(16, *) line_header
	END DO
	
	line_header = "+"
	DO WHILE (line_header.NE.IDGRID)
		READ(20, *) line_header
	END DO

	DO i=1,NDCBUS
		READ(10,*)Udc_ini(i,1)
		READ(11,*)(Ydc(i,j),j=1,NDCBUS)
		READ(12,*)(Ac(i,j),j=1,NDCLINES)
		READ(13,*)Gdc(i,i)
		READ(14,*)Cdc(i,i)
		READ(20,*)CONVERTER_ACDC_BUS(i,1),CONVERTER_ACDC_BUS(i,2)
		Cdc_inv(i,i) = 1/Cdc(i,i)
		CALL MDLIND(CONVERTER_ACDC_BUS(i,2), MACHID(I_MACH), 'GEN', 'VAR', IVAL, IERR)
		I_VAR_SVSCON(i,1) = IVAL ! starting array index of each converter model
	END DO

	DO i=1,NDCLINES
		READ(15,*)Rsdc(i,i)
		READ(16,*)Ldc(i,i)
		Rsdc_inv(i,i) = 1/Rsdc(i,i)
		Ldc_inv(i,i) = 1/Ldc(i,i)
	END DO

	CLOSE(10)
	CLOSE(11)
	CLOSE(12)
	CLOSE(13)
	CLOSE(14)
	CLOSE(15)
	CLOSE(16)
	CLOSE(20)

END