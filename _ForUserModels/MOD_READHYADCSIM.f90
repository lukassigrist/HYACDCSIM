MODULE MOD_READHYADCSIM

	IMPLICIT NONE

CONTAINS

    SUBROUTINE SUB_READDCBUS(v_udc0, m_VSCACDCBUS, IDGRID, NDCBUS)
        
        CHARACTER(1), INTENT(IN) :: IDGRID
        INTEGER, INTENT(IN) :: NDCBUS

        INTEGER, INTENT(INOUT) :: m_VSCACDCBUS(NDCBUS,2)
        REAL,INTENT(INOUT) :: v_udc0(NDCBUS,1)

        CHARACTER(1) :: line_header 
        INTEGER :: ivsc ! ac and dc buses of each converter
        
        
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
            READ(20,*)m_VSCACDCBUS(ivsc,1),m_VSCACDCBUS(ivsc,2)
            READ(21,*)v_udc0(ivsc,1)
        END DO
        CLOSE(20)
        CLOSE(21)

    END SUBROUTINE SUB_READDCBUS

    SUBROUTINE SUB_READDCGRID(v_IVARVSC, v_udc0, m_Ydc, m_Ac, m_Gdc, m_Cdcinv, m_Rsdc, m_Rsdcinv, m_Ldcinv, I_MACH, m_VSCACDCBUS, IDGRID, NDCBUS, NDCLINES)

        INCLUDE 'COMON4.INS'
        IMPLICIT NONE
        
        CHARACTER(1), INTENT(IN) :: IDGRID 
        INTEGER, INTENT(IN) :: I_MACH, NDCBUS, NDCLINES

        INTEGER, INTENT(INOUT) ::  v_IVARVSC(NDCBUS,1), m_VSCACDCBUS(NDCBUS,2)
        REAL, INTENT(INOUT) :: m_Ydc(NDCBUS,NDCBUS), m_Ac(NDCBUS,NDCLINES), m_Gdc(NDCBUS,NDCBUS), m_Rsdc(NDCLINES,NDCLINES)
        REAL, INTENT(INOUT)  ::  m_Cdcinv(NDCBUS,NDCBUS), m_Ldcinv(NDCLINES,NDCLINES), m_Rsdcinv(NDCLINES,NDCLINES)
        REAL, INTENT(INOUT)  :: v_udc0(NDCBUS,1)

        REAL :: m_Cdc(NDCBUS,NDCBUS), m_Ldc(NDCLINES,NDCLINES)
        INTEGER :: ierr, IVAL, i, j
        CHARACTER(1) :: line_header
    
        OPEN(UNIT=10, FILE='.\data_Udc_ini.txt')
        OPEN(UNIT=11, FILE='.\data_Ydc.txt')
        OPEN(UNIT=12, FILE='.\data_Ac.txt')
        OPEN(UNIT=13, FILE='.\data_Gdc.txt')
        OPEN(UNIT=14, FILE='.\data_Cdc.txt')
        OPEN(UNIT=15, FILE='.\data_Rsdc.txt')
        OPEN(UNIT=16, FILE='.\data_Ldc.txt')
        OPEN(UNIT=20, FILE='.\data_acdcbus.txt')
    
    
        ! Salta las l�neas hasta encontrar la letra de su RED [header]
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
            READ(10,*)v_udc0(i,1)
            READ(11,*)(m_Ydc(i,j),j=1,NDCBUS)
            READ(12,*)(m_Ac(i,j),j=1,NDCLINES)
            READ(13,*)m_Gdc(i,i)
            READ(14,*)m_Cdc(i,i)
            READ(20,*)m_VSCACDCBUS(i,1),m_VSCACDCBUS(i,2)
            m_Cdcinv(i,i) = 1/m_Cdc(i,i)
            CALL MDLIND(m_VSCACDCBUS(i,2), MACHID(I_MACH), 'GEN', 'VAR', IVAL, ierr)
            v_IVARVSC(i,1) = IVAL ! starting array index of each converter model
        END DO
    
        DO i=1,NDCLINES
            READ(15,*)m_Rsdc(i,i)
            READ(16,*)m_Ldc(i,i)
            m_Rsdcinv(i,i) = 1/m_Rsdc(i,i)
            m_Ldcinv(i,i) = 1/m_Ldc(i,i)
        END DO
    
        CLOSE(10)
        CLOSE(11)
        CLOSE(12)
        CLOSE(13)
        CLOSE(14)
        CLOSE(15)
        CLOSE(16)
        CLOSE(20)
    
    END SUBROUTINE SUB_READDCGRID

END MODULE MOD_READHYADCSIM