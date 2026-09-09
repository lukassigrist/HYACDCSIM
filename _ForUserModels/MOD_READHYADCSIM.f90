MODULE MOD_READHYADCSIM

	IMPLICIT NONE

CONTAINS

    SUBROUTINE SUB_READDCBUS(v_udc0, m_VSCACDCBUS, IDGRID, NDCBUS)
        
        CHARACTER(1), INTENT(IN) :: IDGRID
        INTEGER, INTENT(IN) :: NDCBUS

        INTEGER, INTENT(INOUT) :: m_VSCACDCBUS(NDCBUS,2)
        REAL,INTENT(INOUT) :: v_udc0(NDCBUS,1)

        REAL :: IDGRID_num 
        INTEGER :: network_id
        INTEGER :: i 

        REAL :: m1, m2, v1

        READ(IDGRID, *) IDGRID_num ! Convert char to number

        OPEN(UNIT=20, FILE='.\data_Buses_base.txt') ! .txt files where .dll is located

        i = 1
        DO WHILE (.TRUE.)
            READ(20, *, END=100) network_id, m1, m2, v1

            IF (network_id .NE. IDGRID_num) CYCLE 

            m_VSCACDCBUS(i,1) = m1
            m_VSCACDCBUS(i,2) = m2
            v_udc0(i,1)        = v1

            i = i + 1
            IF (i > NDCBUS) EXIT
        END DO
        100 CONTINUE

        CLOSE(20)

    END SUBROUTINE SUB_READDCBUS

    SUBROUTINE SUB_READDCGRID(v_IVARVSC, v_udc0, m_Ydc, m_Zdc, m_Ac, m_Gdc, m_Cdcinv, m_Rsdc, m_Rsdcinv, m_Ldcinv, I_MACH, m_VSCACDCBUS, IDGRID, NDCBUS, NDCLINES)

        INCLUDE 'COMON4.INS'
        IMPLICIT NONE
        
        CHARACTER(1), INTENT(IN) :: IDGRID 
        INTEGER, INTENT(IN) :: I_MACH, NDCBUS, NDCLINES

        INTEGER, INTENT(INOUT) ::  v_IVARVSC(NDCBUS,1), m_VSCACDCBUS(NDCBUS,2)
        REAL, INTENT(INOUT) :: m_Ydc(NDCBUS,NDCBUS), m_Zdc(NDCBUS,NDCBUS), m_Ac(NDCBUS,NDCLINES), m_Gdc(NDCBUS,NDCBUS), m_Rsdc(NDCLINES,NDCLINES)
        REAL, INTENT(INOUT)  ::  m_Cdcinv(NDCBUS,NDCBUS), m_Ldcinv(NDCLINES,NDCLINES), m_Rsdcinv(NDCLINES,NDCLINES)
        REAL, INTENT(INOUT)  :: v_udc0(NDCBUS,1)

        REAL :: m_Cdc(NDCBUS,NDCBUS), m_Ldc(NDCLINES,NDCLINES)
        INTEGER :: ierr, IVAL, i, j
        CHARACTER(1) :: line_header

        REAL :: IDGRID_num 
        INTEGER :: network_id

        REAL :: m1, m2, v1, m3, m4, m5, m6

        READ(IDGRID, *) IDGRID_num ! Convert char to number

        OPEN(UNIT=10, FILE='.\data_Buses_base.txt')

        OPEN(UNIT=11, FILE='.\data_Ydc.txt')
        OPEN(UNIT=12, FILE='.\data_Zdc.txt')
        OPEN(UNIT=13, FILE='.\data_Ac.txt')
    
        OPEN(UNIT=14, FILE='.\data_Lines_base.txt')
        
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
          
        i = 1
        DO WHILE (.TRUE.)
            READ(10, *, END=100) network_id, m1, m2, v1, m3, m4

            IF (network_id .NE. IDGRID_num) CYCLE 

            m_VSCACDCBUS(i,1) = m1
            m_VSCACDCBUS(i,2) = m2
            v_udc0(i,1) = v1
            m_Gdc(i,i) = m3
            m_Cdc(i,i) = m4

            m_Cdcinv(i,i) = 1/m_Cdc(i,i)
            
            IF (m2 .NE. -1) THEN

                CALL MDLIND(m_VSCACDCBUS(i,2), MACHID(I_MACH), 'GEN', 'VAR', IVAL, ierr)
                v_IVARVSC(i,1) = IVAL ! starting array index of each converter model

            ELSE

                v_IVARVSC(i,1) = 0  ! For DC buses without VSC

            END IF
            
            i = i + 1    
            IF (i > NDCBUS) EXIT
        END DO
        100 CONTINUE
        
        DO i=1,NDCBUS

            READ(11,*)(m_Ydc(i,j),j=1,NDCBUS)
            READ(12,*)(m_Zdc(i,j),j=1,NDCBUS)
            READ(13,*)(m_Ac(i,j),j=1,NDCLINES)
            
        END DO

        i = 1
        DO WHILE (.TRUE.)
            READ(14, *, END=200) network_id, m5, m6

            IF (network_id .NE. IDGRID_num) CYCLE 

            m_Rsdc(i,i) = m5
            m_Ldc(i,i) = m6
            
            m_Rsdcinv(i,i) = 1/m_Rsdc(i,i)
            m_Ldcinv(i,i) = 1/m_Ldc(i,i)

            i = i + 1
            IF (i > NDCLINES) EXIT
        END DO
        200 CONTINUE
    
        CLOSE(10)
        CLOSE(11)
        CLOSE(12)
        CLOSE(13)
        CLOSE(14)

    END SUBROUTINE SUB_READDCGRID

END MODULE MOD_READHYADCSIM