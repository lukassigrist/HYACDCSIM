! Dynamic model of the DC grid
!
! The DC grid makes use of the Pi-model of lines. The bus-end susceptances are combined into a single bus susceptance. The
! RL element forms the series branch. 
!
! In this implementation, the DC grid is solved through an algebraic DC power-flow formulation. The DC-bus voltages and
! DC-line currents are calculated in CASE 3 from the converter DC injections, stored as VARs, and provided to the converter
! models for subsequent control calculations.
!
! To reflect this static DC-grid representation, the model is named SDCGRD.
	
MODULE MOD_DCGRID

    IMPLICIT NONE

    REAL, ALLOCATABLE :: A_GLOBAL(:,:), B_GLOBAL(:,:)
    REAL, ALLOCATABLE :: C_GLOBAL(:,:), D_GLOBAL(:,:)

    ! Stores the VAR index of each converter. This must be INTEGER because it is used as an array index.
    INTEGER, ALLOCATABLE :: I_VARCONV_GLOBAL(:,:)

    INTEGER :: aux_var_GLOBAL = 0

    INTEGER, ALLOCATABLE :: CONVERTER_ACDC_BUS_GLOBAL(:,:)
    REAL, ALLOCATABLE :: Ydc_GLOBAL(:,:)
    REAL, ALLOCATABLE :: Ac_GLOBAL(:,:)

    ! Automatic offset and total-size calculation for multiple DC grids.
    ! The user does not need to enter: NDCBUS_PREVIOUS, NDCLINES_PREVIOUS, NTOTDCBUS, NTOTDCLINES
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

        CALL GROW_REAL2D(Ydc_GLOBAL, NTOTDCBUS_AUTO, NTOTDCBUS_AUTO, IERR_MEM)
        IF (IERR_MEM .NE. 0) RETURN

        CALL GROW_REAL2D(Ac_GLOBAL, NTOTDCBUS_AUTO, NTOTDCLINES_AUTO, IERR_MEM)
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


SUBROUTINE SDCGRD(I_MACH,I_SLOT)

	! for the global variables (matrices)
	USE MOD_DCGRID
	USE MOD_READHYADCSIM
	INCLUDE 'COMON4.INS'
	IMPLICIT none


	INTEGER I_SLOT, I_MACH
	INTEGER I_ICON, I_CON, I_STATE, I_VAR
	INTEGER NDCBUS, NDCLINES, NPOLES									! number of DC buses, lines, poles
	INTEGER NDCBUS_PREVIOUS, NDCLINES_PREVIOUS	                        ! automatically calculated offsets for DC buses and lines
	INTEGER, ALLOCATABLE :: I_VARCONV(:,:) 							! array position of the converters
	INTEGER, ALLOCATABLE :: CONVERTER_ACDC_BUS(:,:) 	                ! ac and dc buses of each converter
	INTEGER i,j
    INTEGER IERR_REG, IERR_MEM

	CHARACTER(1) IDGRID
    INTEGER :: IDGRID_num
	
	REAL, ALLOCATABLE :: Ydc(:,:) 										! admittance matrix
	REAL, ALLOCATABLE :: Zdc(:,:) 										! Impedance matrix
	REAL, ALLOCATABLE :: Ac(:,:), Ac_T(:,:)								! incidence matrix
	REAL, ALLOCATABLE :: Gdc(:,:), Cdc(:,:)								! bus conductance and capacitor matrix (diagonal)
	REAL, ALLOCATABLE :: Rsdc(:,:), Ldc(:,:)							! line resistance and inductance matrix (diagonal)
	REAL, ALLOCATABLE :: Cdc_inv(:,:), Ldc_inv(:,:), Rsdc_inv(:,:)
	REAL, ALLOCATABLE :: A11(:,:), A12(:,:), A21(:,:), A22(:,:)
	REAL, ALLOCATABLE :: B1(:,:), B2(:,:)
	REAL, ALLOCATABLE :: A(:,:), B(:,:), C(:,:), D(:,:)
	REAL, ALLOCATABLE :: I_sts(:,:)

	REAL, ALLOCATABLE :: Idc_ini(:,:), Udc_ini(:,:), Pdc_ini(:,:)
	REAL, ALLOCATABLE :: Idc(:,:), Udc(:,:), Pdc(:,:) 					! DC-side currents, voltages, power of converters (NDCBUS x 1)
	REAL, ALLOCATABLE :: Icc_ini(:,:), Icc(:,:) 						! Line currents (NDCLINES x 1)

	REAL, ALLOCATABLE :: x(:,:), d_x(:,:), u(:,:), y(:,:)				! x = [udc, icc], u = idc, y = x

    !Related to DC power flow (in case 3) 
    INTEGER :: MAX_ITER
    Real, PARAMETER :: tol_dc=1.0d-6
    INTEGER :: ii, jj, ll, iter_dc, from_bus, to_bus
    Real :: sumYdc
	INTEGER :: info
	Real :: Gline 
    INTEGER :: i_slack, n_active, idx, ii2, jj2, kk
    INTEGER, ALLOCATABLE :: active_bus(:)
    REAL, ALLOCATABLE :: Udc_active(:,:), Pdc_active(:,:)
    REAL, ALLOCATABLE :: Ydc_active(:,:)
    REAL, ALLOCATABLE :: F_active(:,:), delta_Udc_active(:,:), Jac_active(:,:)
    INTEGER, ALLOCATABLE :: ipiv_active(:)
    REAL :: sumYdc_active
    INTRINSIC :: ABS, MAX, MAXVAL
    REAL, PARAMETER :: epsU=0.0000000001

	! Indexes
	I_CON=STRTIN(1,I_SLOT)
	I_VAR=STRTIN(3,I_SLOT)
	I_ICON=STRTIN(4,I_SLOT)
   
	! ICONs assignment
	NDCBUS = ICON(I_ICON) 				! number of DC buses of the ith MTDC
	NDCLINES = ICON(I_ICON+1) 			! number of DC lines of the ith MTDC
	NPOLES = ICON(I_ICON+2) 			! number of poles of each converter
	IDGRID = CHRICN(I_ICON+3)
	MAX_ITER = ICON(I_ICON+4)			! max number of DC load flow iteration

    ! Automatic calculation of NDCBUS_PREVIOUS and NDCLINES_PREVIOUS
	
    CALL REGISTER_DCG(IDGRID, NDCBUS, NDCLINES, NDCBUS_PREVIOUS, NDCLINES_PREVIOUS, IERR_REG)

	IF (IERR_REG .NE. 0) THEN
		WRITE (LPDEV,*) 'SDCGRD - ERROR: automatic DC-grid offset registration failed. IDGRID = ', IDGRID, ' IERR_REG = ', IERR_REG
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

	ALLOCATE(u(NDCBUS,1))
	ALLOCATE(y(NDCBUS+NDCLINES,1))
    
	! VARs assignment
	DO i=1,(NDCBUS+NDCLINES)   
		IF (i.LE.NDCBUS) THEN 
			Pdc(i,1) = VAR(I_VAR + i - 1) 			                      ! DC power of converters (stored as VARs)
			Idc(i,1) = VAR(I_VAR + NDCBUS - 1 + i )                       ! DC current of converters (stored as VARs)

			Udc(i,1) = VAR(I_VAR + 2*NDCBUS - 1 + i ) 				      ! DC voltage
		ELSE
			Icc(i-NDCBUS,1) = VAR(I_VAR + 3*NDCBUS + (i-NDCBUS) - 1) 	  ! DC branch currents
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
			WRITE (LPDEV,*) 'SDCGRD - ERROR: global storage allocation failed. IDGRID = ', IDGRID, ' IERR_MEM = ', IERR_MEM
			RETURN
		END IF
	  
		! Read DC-grid data from .txt file
		CALL SUB_READDCGRID(I_VARCONV, Udc_ini, Ydc, Zdc, Ac, Gdc, Cdc_inv, Rsdc, Rsdc_inv, Ldc_inv, I_MACH, CONVERTER_ACDC_BUS, IDGRID, NDCBUS, NDCLINES)

        READ(IDGRID, *) IDGRID_num ! Convert char to number
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
		Icc_ini = MATMUL(Rsdc_inv, MATMUL(Ac_T,Udc_ini))  		    ! Icc_ini = Rsdc_inv*(Ac'*Udc_ini))
		
		! store global matrices
		CALL SUB_LOCALTOGLOBAL(A, A_GLOBAL, NDCBUS+NDCLINES, NDCBUS+NDCLINES, SIZE(A_GLOBAL,1), SIZE(A_GLOBAL,2), NDCBUS_PREVIOUS + NDCLINES_PREVIOUS, NDCBUS_PREVIOUS + NDCLINES_PREVIOUS)

		CALL SUB_LOCALTOGLOBAL(B, B_GLOBAL, NDCBUS+NDCLINES, NDCBUS, SIZE(B_GLOBAL,1), SIZE(B_GLOBAL,2), NDCBUS_PREVIOUS + NDCLINES_PREVIOUS, NDCBUS_PREVIOUS)

		CALL SUB_LOCALTOGLOBAL(C, C_GLOBAL, NDCBUS+NDCLINES, NDCBUS+NDCLINES, SIZE(C_GLOBAL,1), SIZE(C_GLOBAL,2), NDCBUS_PREVIOUS + NDCLINES_PREVIOUS, NDCBUS_PREVIOUS + NDCLINES_PREVIOUS)

        CALL SUB_LOCALTOGLOBAL(Ydc, Ydc_GLOBAL, NDCBUS, NDCBUS, SIZE(Ydc_GLOBAL,1), SIZE(Ydc_GLOBAL,2), NDCBUS_PREVIOUS, NDCBUS_PREVIOUS)

        CALL SUB_LOCALTOGLOBAL(Ac, Ac_GLOBAL, NDCBUS, NDCLINES, SIZE(Ac_GLOBAL,1), SIZE(Ac_GLOBAL,2), NDCBUS_PREVIOUS, NDCLINES_PREVIOUS)

		DO i = 1,NDCBUS
			I_VARCONV_GLOBAL(NDCBUS_PREVIOUS+i,1) = I_VARCONV(i,1)
		END DO
		
		! Initialization
		! --------------
		Udc = Udc_ini
		Idc = Idc_ini
		Pdc = Pdc_ini
		Icc = Icc_ini

		  
		! report mode 
		WRITE (LPDEV,*) 'SDCGRD - CASE 1: DC grid at AC bus ',NUMBUS(NUMTRM(I_MACH)),' with id ',MACHID(I_MACH),' initialized.'
 
	CASE (2)

        ! CASE 2 is kept inactive; DC grid variables are updated using the DC power-flow formulation in CASE 3 instead of state derivatives.
	  	
	CASE (3) 
        
		CALL SUB_GLOBALTOLOCAL(Ydc, Ydc_GLOBAL, NDCBUS, NDCBUS, SIZE(Ydc_GLOBAL,1), SIZE(Ydc_GLOBAL,2), NDCBUS_PREVIOUS, NDCBUS_PREVIOUS)

        CALL SUB_GLOBALTOLOCAL(Ac, Ac_GLOBAL, NDCBUS, NDCLINES, SIZE(Ac_GLOBAL,1), SIZE(Ac_GLOBAL,2), NDCBUS_PREVIOUS, NDCLINES_PREVIOUS)

        DO i = 1, NDCBUS
		    I_VARCONV(i,1) = I_VARCONV_GLOBAL(NDCBUS_PREVIOUS + i, 1)
        END DO
        
		READ(IDGRID, *) IDGRID_num ! Convert char to number
        DO i = 1, NDCBUS
    		CONVERTER_ACDC_BUS(i,1) = CONVERTER_ACDC_BUS_GLOBAL(NDCBUS_PREVIOUS + i, 1)
    		CONVERTER_ACDC_BUS(i,2) = CONVERTER_ACDC_BUS_GLOBAL(NDCBUS_PREVIOUS + i, 2)
		END DO

        DO i = 1, NDCBUS
            IF (CONVERTER_ACDC_BUS(i,2) .NE. -1.0) THEN
                Pdc(i,1) = VAR(I_VARCONV(i,1) + 4)
            ELSE
                Pdc(i,1) = 0.0
            END IF
        END DO
        
		! DC power-flow solution: converter powers are used as inputs, and DC bus voltages and line currents are updated algebraically.
        i_slack = 0

        DO i = 1, NDCBUS
            IF (CONVERTER_ACDC_BUS(i,2) .NE. -1.0 .AND. ABS(Pdc(i,1)) .GE. 0.0000001) THEN
                i_slack = i
                GO TO 100
            END IF
        END DO

		100 CONTINUE

		IF (i_slack .EQ. 0) THEN
            WRITE (LPDEV,*) 'SDCGRD - ERROR: No slack bus found (no converter bus).'
            RETURN
        END IF
		
        ALLOCATE(active_bus(NDCBUS))
        n_active = 0
        DO i = 1, NDCBUS
            IF (i .NE. i_slack) THEN
                n_active = n_active + 1
                active_bus(n_active) = i
            END IF
        END DO

        ALLOCATE(Udc_active(n_active,1))
        ALLOCATE(Pdc_active(n_active,1))
        ALLOCATE(Ydc_active(n_active,n_active))
        ALLOCATE(F_active(n_active,1))
        ALLOCATE(delta_Udc_active(n_active,1))
        ALLOCATE(Jac_active(n_active,n_active))
        ALLOCATE(ipiv_active(n_active))

        DO i = 1, n_active
            ii = active_bus(i)
            Udc_active(i,1) = Udc(ii,1)
            Pdc_active(i,1) = Pdc(ii,1)
            DO j = 1, n_active
            	jj = active_bus(j)
            	Ydc_active(i,j) = Ydc(ii,jj)
        	END DO
    	END DO

        DO iter_dc = 1, MAX_ITER
        	DO i = 1, n_active
			    ii = active_bus(i)
            	F_active(i,1) = 0.0
           		DO j = 1, n_active
				    jj = active_bus(j)
                	F_active(i,1) = F_active(i,1) + Ydc_active(i,j) * Udc_active(j,1)
            	END DO
                
				F_active(i,1) = F_active(i,1) + Ydc(ii, i_slack) * Udc(i_slack,1) 

                IF (ABS(Udc_active(i,1)) .LT. epsU) THEN
				    WRITE(LPDEV,*) 'SDCGRD - ERROR: Udc too small at bus ', active_bus(i), '  Udc=', Udc_active(i,1)
                    RETURN
				END IF

                F_active(i,1) = F_active(i,1) - Pdc_active(i,1) / Udc_active(i,1)

        	END DO

        	IF (MAXVAL(ABS(F_active)) .LT. tol_dc) GO TO 200
            
            DO i = 1, n_active

			    IF (ABS(Udc_active(i,1)) .LT. epsU) THEN
				    WRITE(LPDEV,*) 'SDCGRD - ERROR: Udc too small at bus ', active_bus(i), '  Udc=', Udc_active(i,1)
				    RETURN
				END IF	

            	DO j = 1, n_active
                	IF (i .EQ. j) THEN

                    	Jac_active(i,i) = Ydc_active(i,i) + Pdc_active(i,1) / (Udc_active(i,1) * Udc_active(i,1))

                	ELSE
                    	Jac_active(i,j) = +Ydc_active(i,j)
                	END IF
            	END DO
        	END DO

            delta_Udc_active(:,1) = -F_active(:,1)
        	CALL SUB_DGESV(n_active, 1, Jac_active, n_active, ipiv_active, delta_Udc_active, n_active, info)

        	IF (info .NE. 0) THEN
            	WRITE (LPDEV,*) 'DC power flow failed at bus ', active_bus(info)
            	RETURN
        	END IF

        	Udc_active(:,1) = Udc_active(:,1) + delta_Udc_active(:,1)

    	END DO
       
		GO TO 300

		200 CONTINUE

		300 CONTINUE

		DO i = 1, n_active
        	ii = active_bus(i)
        	Udc(ii,1) = Udc_active(i,1)
    	END DO
       
	    DEALLOCATE(active_bus)
    	DEALLOCATE(Udc_active, Pdc_active, Ydc_active)
    	DEALLOCATE(F_active, delta_Udc_active, Jac_active, ipiv_active)
        
		DO i = 1, NDCBUS
    		IF (Udc(i,1) .NE. 0.0) THEN
				!CALL SUB_DIVTBT(Pdc, Udc, Idc, NDCBUS, 1) ! Idc = Pdc/Udc
				Idc(i,1) = Pdc(i,1) / Udc(i,1)
    		ELSE
        		Idc(i,1) = 0.0
    		END IF
		END DO

        !!!! DC line current calculation using Ac and Ydc !!!!
		DO ll = 1, NDCLINES
        	from_bus = -1
        	to_bus   = -1

			! Find "from" and "to" bus of line ll
        	DO ii = 1, NDCBUS
            	IF (Ac(ii,ll) ==  1.0) from_bus = ii   
            	IF (Ac(ii,ll) == -1.0) to_bus   = ii
        	END DO

			! Safety check
        	IF (from_bus .LT. 1 .OR. to_bus .LT. 1) THEN
            	WRITE (LPDEV,*) 'Error: Line ', ll, ' bus not defined in Ac matrix!'
            	RETURN
        	END IF
        	Gline = ABS(Ydc(from_bus, to_bus))                                       ! Calculate line conductance from Ydc
        	Icc(ll,1) = Gline * (Udc(from_bus,1) - Udc(to_bus,1))                    ! Calculate DC line current
    	END DO
        
		! Outputs
		! =======

        ! Fill y: first voltages, then line currents
		DO ii = 1, NDCBUS
        	y(ii,1) = Udc(ii,1)
    	END DO
    	DO ii = 1, NDCLINES
        	y(NDCBUS + ii, 1) = Icc(ii,1)
    	END DO

	CASE (4)
	  
	CASE (5)
		! reporting mode

		WRITE (LPDEV,*) 'SDCGRD - ', 'CON', I_CON
		WRITE (LPDEV,*) 'SDCGRD - ', 'ICON', I_ICON
		WRITE (LPDEV,*) 'SDCGRD - ', 'VAR', I_VAR
		WRITE (LPDEV,*) 'SDCGRD - ', 'STATE', I_STATE

	CASE DEFAULT
		  
	END SELECT

	! VARs assignment
	DO i=1,NDCBUS
        VAR(I_VAR + i - 1) = Pdc(i,1)       
        VAR(I_VAR + NDCBUS + i - 1) = Idc(i,1)  
		VAR(I_VAR + 2*NDCBUS + i - 1) = Udc(i,1)
	END DO

    DO ii = 1, NDCLINES
        VAR(I_VAR + 3*NDCBUS + ii - 1)   = Icc(ii,1) 
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
	! Local X dimensions, Global X dimensions, and previous offsets are passed for compatibility.
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
	! Local X dimensions, Global X dimensions, and previous offsets are passed for compatibility.
	INTEGER :: NROW_LOCAL, NCOL_LOCAL, NROW_GLOBAL, NCOL_GLOBAL, NROW_PREVIOUS, NCOL_PREVIOUS
	REAL :: X(NROW_LOCAL,NCOL_LOCAL), X_GLOBAL(NROW_GLOBAL,NCOL_GLOBAL)
	INTEGER :: i, j

	DO i = 1, NROW_LOCAL
		DO j = 1, NCOL_LOCAL
			X_GLOBAL(NROW_PREVIOUS + i, NCOL_PREVIOUS + j) = X(i,j)
		END DO
	END DO
END


!**********************************************************
!   Simple replacement for DGESV - solves A * X = B
!   Gaussian elimination with partial pivoting
!   Matrix A and vector B are overwritten (like DGESV itself)
!**********************************************************
SUBROUTINE SUB_DGESV(N, NRHS, A, LDA, IPIV, B, LDB, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, NRHS, LDA, LDB
    REAL, INTENT(INOUT) :: A(LDA,*), B(LDB,*) 
    INTEGER, INTENT(OUT) :: IPIV(*), INFO
    
    INTEGER :: I, J, K, P
    REAL :: FACTOR, TEMP, MAX_VAL 
    
    INFO = 0
    IF (N == 0) RETURN
    
    DO K = 1, N-1

        P = K
        MAX_VAL = ABS(A(K,K))
        DO I = K+1, N
            IF (ABS(A(I,K)) > MAX_VAL) THEN
                MAX_VAL = ABS(A(I,K))
                P = I
            END IF
        END DO
        
        IPIV(K) = P
        
        IF (MAX_VAL == 0.0D0) THEN
            INFO = K
            RETURN
        END IF
        
        IF (P /= K) THEN
            DO J = 1, N
                TEMP = A(K,J)
                A(K,J) = A(P,J)
                A(P,J) = TEMP
            END DO

            DO J = 1, NRHS
                TEMP = B(K,J)
                B(K,J) = B(P,J)
                B(P,J) = TEMP
            END DO
        END IF
        
        DO I = K+1, N
            FACTOR = A(I,K) / A(K,K)
            A(I,K) = FACTOR  
            DO J = K+1, N
                A(I,J) = A(I,J) - FACTOR * A(K,J)
            END DO

            DO J = 1, NRHS
                B(I,J) = B(I,J) - FACTOR * B(K,J)
            END DO
        END DO
    END DO
    
    IPIV(N) = N
    IF (A(N,N) == 0.0D0) THEN
        INFO = N
        RETURN
    END IF
    
    DO J = 1, NRHS
        B(N,J) = B(N,J) / A(N,N)
        DO I = N-1, 1, -1
            TEMP = B(I,J)
            DO K = I+1, N
                TEMP = TEMP - A(I,K) * B(K,J)
            END DO
            B(I,J) = TEMP / A(I,I)
        END DO
    END DO
    
END SUBROUTINE SUB_DGESV
