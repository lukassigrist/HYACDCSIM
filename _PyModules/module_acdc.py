"""
This code sets up and solves the sequential AC/DC power flow algorithm for multi-terminal VSC-HVDC systems

PSS/E solves the AC grid, whereas a Python module solves the DC grid and the copupling between both the AC and DC grids. 

This codes enables:

    1. solving the AC/DC power flow
    2. showing the solution of the AC/DC power flow
    3. writing data files for dynamic simulations

The module contains the following functions:

- main_runacdclf: principal function to prepare and run the AC/DC power flow
- fun_readMTDCdata: extract the MTDC system data from a predefined Excel file
- fun_getT2Piequivalent: build Pi equivalent of T model of transformer and LC filter
- fun_setinitialvaluesacdclf: initialize sequential AC/DC power flow 
- fun_setACsetpoints: set AC-side generator set points
- fun_backwardcomputation: do a backward computation by explicitly representing the transformer and shunt branch, and then recompute the voltage and power set points at point f
- fun_addnewelements: Add an artificial bus, a shunt filter, and a transformer filter, which are connected between the AC bus and the artificial bus
- fun_relatedArtificialbuses: retrieve the artificial bus corresponding to a specific AC bus
- fun_getACresults: extract results from the AC power flow run
- fun_calculateVSCtransformerlosses: calculate the active and reactive power losses of VSC-connected transformers
- fun_initializeDCfromAC: initialize DC-side variables from AC results
- fun_calldclf: call the DC power flow (DC power flow and DC slack bus iteration)
- fun_getacdcoutput: get AC/DC power flow results

- fun_showacdclfoutput: show results

- fun_setacdcdynamicdata: create network and element dynamic data
- fun_updatedynamicfile: update .dyr file with VSC models and DC grid model
- fun_setDCnetworkdata: set DC network data for dynamic simulation
- fun_setMTDCdyrdata: set VSC model and DC grid model parameters
- fun_write2txtfile: write DC network data to .txt file
- fun_setDCnetworkdata: organize and write MTDC grid data, including buses, branches, and grid parameters, to text files

Authors:
- Aurelio Garcia Cerrada
- Javier Renedo
- Carlos Prieto
- Lukas Sigrist
- Saeed Rezaeian-Marjani
"""

import os, sys, shutil

import psse34 # PSSE
# PSS/e
import psspy, redirect
redirect.psse2py()
from psspy import _i, _f, _s, _o
# import numpy and math
import numpy as np
import pandas as pd                      # Added a new library due to Excel file reading
from math import pi
# import hybrid sequential ACDC power flow
from _PyModules.module_dclf import fun_getDCconductance, fun_getDCincidence, fun_mainDCloadflow, fun_slackiterationDC 

# indexes for PYFORMAT data 
BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)
F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, RATIO, ANGLE, STATUS, ANGMIN, ANGMAX, PF, QF, PT, QT = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)
GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, PC1, PC2, QC1MIN, QC1MAX, QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)
FIL_DCBUS, FIL_ACBUS, FIL_RT, FIL_XT, FIL_BF, FIL_RC, FIL_XC, FIL_A, FIL_B, FIL_CRECT, FIL_CINV, FIL_PMAX, FIL_PMIN, FIL_QMAX, FIL_QMIN, FIL_RATE_A, FIL_RATE_B, FIL_RATE_C, FIL_RATIO, FIL_ANGLE, FIL_STATUS, FIL_ANGMIN, FIL_ANGMAX = (0, 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19, 20, 21, 22)
DC_BUS, DC_TYPE, DC_TYPE2, DC_US, DC_DELTAS, DC_PS, DC_QS, DC_UDC, DC_PDC, DC_IDC, DC_GDC, DC_CDC, DC_AREA, DC_BASE_KV, DC_ZONE, DC_VMAX, DC_VMIN = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)
DC_F_BUS, DC_T_BUS, DC_BR_R, DC_BR_L, DC_BR_C, DC_RATE_A, DC_RATE_B, DC_RATE_C, DC_RATIO, DC_ANGLE, DC_STATUS, DC_PIJ, DC_PJI, DC_ICCIJ = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)

DC_PF = DC_PS
DC_UF = DC_US
DC_DELTAF = DC_DELTAS
DC_QF = DC_QS

# Constants
j = 1j

def main_runacdclf(d_acdcoptions, str_savfileorig, str_pathdclfresults, str_MTDCdatafile):
    """
    Sequential AC/DC Power Flow algorithm for multi-terminal VSC-HVDC systems   
    """
    # ===================================================================
    # Preparing sequential ACDC power flow
    # ===================================================================    
     
    # get acdc options
    isprogress = d_acdcoptions.get('isprogress')
    issavedclfresults = d_acdcoptions.get('issavedclfresults')
    issetdynamicfiles = d_acdcoptions.get('issetdynamicfiles')
    idowner = d_acdcoptions.get('idowner')           
    idconv = d_acdcoptions.get('idconv')  
    lftol = d_acdcoptions.get('lftol')
    lfmaxiter = d_acdcoptions.get('lfmaxiter')      
    
    # progress options
    if isprogress==2:
        str_progressfile = os.path.join(str_pathdclfresults,r"progress.dat")
        psspy.progress_output(isprogress,str_progressfile,[0,0])
    else:
        psspy.progress_output(isprogress,r"",[0,0]) 
    
    # modified file name
    str_savfilemtdc = str_savfileorig[:-4] + r"_MTDCg.sav" # update AC power flow (VSC as generators)
    str_rawfilemtdc = str_savfilemtdc[:-4] + r".raw"
    
    # open PSS/e file
    psspy.case(str_savfileorig)

    # system VA base
    ACSbase = psspy.sysmva()

    # System AC bus numbers and the maximum count of AC buses
    ierr, l_acbus = psspy.abusint(-1, 2, 'NUMBER') # For determining AC buses and its bus No. (SID<0 & FLAG = 2: All buses (in-Service and out-of-service) 
    nACbus = len(l_acbus[0])
    maxnACbus = max(l_acbus[0])
    
    # get static invariant DC data
    [DCSbase, nDCbus, m_converter, m_dcbus, m_dcbranch, l_MTDCgrids, l_MTDCnetworkid, nVSC] = fun_readMTDCdata(str_MTDCdatafile, ACSbase)  
    
    # get converter losses parameters
    v_aloss = m_converter[:, FIL_A]
    v_bloss = m_converter[:, FIL_B]
    v_crect = m_converter[:, FIL_CRECT]
    v_cinv = m_converter[:, FIL_CINV]
    
    # To consider the explicit model of the low-pass filter and transformer filter, only the phase reactor between points f and c is included
    [v_Am, v_Bm, v_Cm, v_Dm, v_Y1, v_Z2, v_Y3] = fun_getT2Piequivalent(m_converter)
   
    # ===================================================================
    # Sequential AC/DC power flow iterations
    # =================================================================== 

    # initialize powers and voltages of sequential acdc power flow
    [v_Us, v_Theta_s, v_Ps_pu, v_Qs_pu, v_Pdc_pu, v_PlossVSC] = fun_setinitialvaluesacdclf(DCSbase, nDCbus, m_dcbus, nVSC)
    v_Us0 = v_Us
    v_Qs0_pu = v_Qs_pu
    
    # identify the first running
    isprerun = True    
    
    # Assign Ps, Qs and Us set points to AC-side converters (with the converters still being connected to the point s)   
    l_indexDCbusVSC = fun_setACsetpoints(v_Us, v_Theta_s, v_Ps_pu, v_Us0, v_Qs0_pu, DCSbase, m_converter, m_dcbus, idconv, idowner, isprerun, None)
    isprerun = False
    
    # Solve AC power flow
    psspy.fnsl(
        options1=0, # disable tap stepping adjustment.
        options5=0, # disable switched shunt adjustment.
      # options6=0, # flat start (first time).
        )
    
    [v_Uf, v_Theta_f, v_Uf_phasor, v_Po_pu, v_Qo_pu, v_Qf_pu, v_QLosses_tf_pu, v_Pdc_pu, m_dcbus_f] = fun_backwardcomputation(m_converter, m_dcbus, nVSC, v_Ps_pu, v_Qs_pu, DCSbase, l_indexDCbusVSC)
    v_Uf0 = v_Uf
    v_Theta_f0 = v_Theta_f
    v_Po0_pu = v_Po_pu
    v_Qo0_pu = v_Qo_pu
    v_QLosses_tf0_pu = v_QLosses_tf_pu  # To control the variations in reactive power losses of the VSCs transformer based on potential changes in the current passing through it
    
    l_artificialACbus = fun_addnewelements(nVSC, maxnACbus, DCSbase, idconv, m_converter, v_Uf0, v_Theta_f0, v_Qf_pu)      # Add new elements
    
    # DC slack indexes
    v_indexDCslack = np.where(m_dcbus[l_indexDCbusVSC, DC_TYPE] == 2)[0].tolist()
    
    isconverged = 0 # convergence loop
    it = 0
    v_k_int_ac = []
    while not isconverged: # BEGIN WHILE

        it = it + 1
        v_Po_no_previous = v_Po_pu[v_indexDCslack,0].astype(float)
       
        # ===================================================================
        # AC grid power flow
        # ===================================================================
        
        # 1. Assign Po, Qo and Uf set points to AC-side converters (In this section, the converters are connected to the point f)     
        fun_setACsetpoints(v_Uf, v_Theta_f, v_Po_pu, v_Uf0, v_Qo0_pu, DCSbase, m_converter, m_dcbus, idconv, idowner, isprerun, l_artificialACbus)
         
        # 2. Solve AC power flow
        psspy.fnsl(
            options1=0, # disable tap stepping adjustment.
            options5=0, # disable switched shunt adjustment.
          # options6=0, # flat start (first time).
        )

        v_k_int_ac.append(psspy.iterat())
        
        # 3. Extract solutions of AC power flow 
        [v_Uf, v_Theta_f, v_Po_pu, v_Qo_pu] = fun_getACresults(v_Uf, v_Theta_f, v_Po_pu, v_Qo_pu, DCSbase, m_converter, m_dcbus_f, idconv, l_artificialACbus)      
        
        [_, v_QLosses_tf_pu] = fun_calculateVSCtransformerlosses(l_artificialACbus, DCSbase) # Extraction of VSCs transformer losses after running the power flow in PSS/E.
        v_QLosses_corr_pu = v_QLosses_tf0_pu - v_QLosses_tf_pu  # Correction value for VSCs transformer reactive power losses compensation  
        v_Qo_pu = v_Qo_pu - v_QLosses_corr_pu                   # Updated value of Qo after addition correction value
        
        # 4. Compute the DC-side active power for the subsequent DC-side power flow
        [v_If_phasor, v_Ic_phasor, v_Uc_phasor, v_Pc_pu, v_Qc_pu, v_Pdc_pu, v_PlossVSC] = fun_initializeDCfromAC(v_Uf,
         v_Theta_f, v_Po_pu, v_Qo_pu, v_Pdc_pu, v_PlossVSC, v_Am, v_Bm, v_Cm, v_Dm, v_aloss, v_bloss, v_crect, v_cinv)
        
        m_dcbus_f[l_indexDCbusVSC,DC_PDC] = np.array(v_Pdc_pu).reshape(-1)*DCSbase        
        
        # ===================================================================
        # DC grid power flow
        # ===================================================================
        
        [v_Po_no_current, k_int_dc, v_mismatch_dc, k_int_dcslack, j_int_dcslack, v_Po_pu, v_Pc_pu, v_Qc_pu, v_Pdc_pu, v_Uc_phasor,
         m_dcbus_f, m_dcbranch] = fun_calldclf(v_Uf, v_Theta_f, v_Uc_phasor, v_Po_pu, v_Qo_pu, v_Pc_pu, v_Qc_pu, v_Pdc_pu, DCSbase,
         m_dcbus_f, m_dcbranch, v_Y1, v_Z2, v_Y3, v_aloss, v_bloss, v_crect, v_cinv, l_MTDCgrids, v_indexDCslack, lftol, lfmaxiter,
        l_MTDCnetworkid, l_indexDCbusVSC)
        
        # ===================================================================
        # Convergence 
        # ===================================================================

        error = v_Po_no_current - v_Po_no_previous
        if abs(error).max()<=lftol:
            isconverged = 1
            print('\nConverged after '+str(it)+' iterations.')
            
        if it>=lfmaxiter:
            isconverged = 2
            print('\nMaximum iterations reached ('+str(lfmaxiter)+').')
        
    k_int_ac = max(v_k_int_ac)
    
    # ========================================================================= 
    # Output
    # =========================================================================
    d_acdcoutput = fun_getacdcoutput(v_If_phasor, v_Uf, v_Theta_f, v_Uc_phasor, v_Ic_phasor, v_Pc_pu, v_Qc_pu, v_PlossVSC,
    DCSbase, m_dcbus_f, m_dcbranch, m_converter, l_artificialACbus, l_indexDCbusVSC)

    psspy.fnsl(
        options1=0, # disable tap stepping adjustment.
        options5=0, # disable switched shunt adjustment.
    ) # we have finished, but we run a power flow again to avoid strange things with the tricky thing of PQ generators

    # write AC-side to file
    psspy.rawd_2(0,1,[1,1,1,0,0,0,0],0,str_rawfilemtdc) # save .raw
    psspy.save(str_savfilemtdc) # save .sav   
         
    return d_acdcoutput, isconverged, it, k_int_ac, k_int_dc, k_int_dcslack, j_int_dcslack, v_mismatch_dc, l_MTDCgrids, l_indexDCbusVSC, l_artificialACbus

def fun_readMTDCdata(str_MTDCdatafile, ACSbase):

    ACSbase = float(ACSbase)

    # Read data from different sheets of an Excel file
    df_baseMVA = pd.read_excel(str_MTDCdatafile, sheet_name='DCbase', header=0)
    df_converter = pd.read_excel(str_MTDCdatafile, sheet_name='converter', header=0)
    df_dcbus = pd.read_excel(str_MTDCdatafile, sheet_name='dcbus', header=0)
    df_dcbranch = pd.read_excel(str_MTDCdatafile, sheet_name='dcbranch', header=0)

    l_MTDCnetworkid = df_dcbus.iloc[:, 0].values          # Dterminig the ID of MTDC networks

    # Filter data for each unique MTDC network ID and append the relevant data (baseMVA, dcbus, converter, and dcbranch) into a list of dictionaries
    l_MTDCgrids = []
    for network_id in np.unique(l_MTDCnetworkid):
        
        baseMVA_data = df_baseMVA[df_baseMVA.iloc[:, 0] == network_id].iloc[:, 1:].values
        converter_data = df_converter[df_converter.iloc[:, 0] == network_id].iloc[:, 1:].values
        dcbus_data = df_dcbus[df_dcbus.iloc[:, 0] == network_id].iloc[:, 1:].values
        dcbranch_data = df_dcbranch[df_dcbranch.iloc[:, 0] == network_id].iloc[:, 1:].values

        DCSbasei = float(baseMVA_data[0, 0]) if baseMVA_data.size > 0 else np.nan
        
        converter_data = 1.0 * np.array(converter_data)
        dcbus_data = 1.0 * np.array(dcbus_data)
        dcbranch_data = 1.0 * np.array(dcbranch_data)

        # -------------------------------
        # Converter data conversion
        # -------------------------------

        # Transformer and converter impedance parameters
        v_MbaseVSC = converter_data[:, FIL_RATE_A]
        converter_data[:, FIL_RT] = np.divide(converter_data[:, FIL_RT],v_MbaseVSC) * (ACSbase)
        converter_data[:, FIL_XT] = np.divide(converter_data[:, FIL_XT],v_MbaseVSC) * (ACSbase)
        converter_data[:, FIL_RC] = np.divide(converter_data[:, FIL_RC],v_MbaseVSC) * (ACSbase)
        converter_data[:, FIL_XC] = np.divide(converter_data[:, FIL_XC],v_MbaseVSC) * (ACSbase)

        # Filter susceptance
        converter_data[:, FIL_BF] = np.multiply(converter_data[:, FIL_BF],v_MbaseVSC) / ACSbase

        # Loss coefficients: still in pu, but converted to ACSbase
        # ploss_pu = a + b*ic + c*ic^2
        converter_data[:, FIL_A] = np.multiply(converter_data[:, FIL_A],v_MbaseVSC) / (ACSbase)
        converter_data[:, FIL_B] = converter_data[:, FIL_B]  # unchanged
        converter_data[:, FIL_CRECT] = np.divide(converter_data[:, FIL_CRECT],v_MbaseVSC) * (ACSbase)
        converter_data[:, FIL_CINV] = np.divide(converter_data[:, FIL_CINV],v_MbaseVSC) * (ACSbase)

        # -------------------------------
        # DC bus conversion
        # -------------------------------

        # Shunt conductance and capacitance
        dcbus_data[:, DC_GDC] = dcbus_data[:, DC_GDC] * (DCSbasei / ACSbase)
        dcbus_data[:, DC_CDC] = dcbus_data[:, DC_CDC] * (DCSbasei / ACSbase)

        # -------------------------------
        # DC branch conversion
        # -------------------------------

        # DC branch resistance and inductance
        dcbranch_data[:, DC_BR_R] = dcbranch_data[:, DC_BR_R] * (ACSbase / DCSbasei)
        dcbranch_data[:, DC_BR_L] = dcbranch_data[:, DC_BR_L] * (ACSbase / DCSbasei)

        # DC branch capacitance
        dcbranch_data[:, DC_BR_C] = dcbranch_data[:, DC_BR_C] * (DCSbasei / ACSbase)
        
        
        l_MTDCgrids.append({
            "baseMVA": ACSbase,
            "original_baseMVA": DCSbasei,
            "dcbus": dcbus_data,
            "converter": converter_data,
            "dcbranch": dcbranch_data
        })  

    # After all conversions, the DC system base is the same as the AC system base
    DCSbase = ACSbase
    
    l_converter = []
    l_dcbus = []
    l_dcbranch = []
    l_nDClines = []
    l_nDCbus = []
    for i in range(0,len(l_MTDCgrids)):
        l_converter.append(l_MTDCgrids[i]["converter"])
        a_dcbus = l_MTDCgrids[i]["dcbus"]
        v_sortedDCbus = np.argsort(a_dcbus[:,DC_BUS]) # sort DC buses
        a_dcbus = a_dcbus[v_sortedDCbus,:]
        l_dcbus.append(a_dcbus)
        l_dcbranch.append(l_MTDCgrids[i]["dcbranch"])
        l_nDClines.append(len(l_MTDCgrids[i]["dcbranch"]))
        l_nDCbus.append(len(a_dcbus))

    # all matrices are numpy matrices with float numbers
    m_converter = 1.0*np.matrix(np.vstack(l_converter))
    m_dcbus = 1.0*np.matrix(np.vstack(l_dcbus))
    m_dcbranch = 1.0*np.matrix(np.vstack(l_dcbranch))
    nDCbus = np.shape(m_dcbus[:,0])[0] # # of nodes of the DC grid
    nVSC = np.shape(m_converter[:,0])[0] # # of converters of the DC grid
        
    return DCSbase, nDCbus, m_converter, m_dcbus, m_dcbranch, l_MTDCgrids, l_MTDCnetworkid, nVSC

def fun_getT2Piequivalent(m_converter):

    # Phase reactor (p.u)
    v_Rc = m_converter[:, FIL_RC]
    v_Xc = m_converter[:, FIL_XC]
    v_Zc = v_Rc + 1j * v_Xc  

    # T model of phase reactor to equivalent Pi model (Simplified Steinmetz)
    v_Y1 = np.matrix(np.zeros(v_Zc.shape, dtype=v_Zc.dtype))  
    v_Z2 = v_Zc  
    v_Y3 = np.matrix(np.zeros(v_Zc.shape, dtype=v_Zc.dtype)) 

    v_Am = np.matrix(np.ones(v_Zc.shape, dtype=v_Zc.dtype))  
    v_Bm = v_Z2
    v_Cm = np.matrix(np.zeros(v_Zc.shape, dtype=v_Zc.dtype))  
    v_Dm = np.matrix(np.ones(v_Zc.shape, dtype=v_Zc.dtype)) 
        
    return v_Am, v_Bm, v_Cm, v_Dm, v_Y1, v_Z2, v_Y3

def fun_setinitialvaluesacdclf(DCSbase, nDCbus, m_dcbus, nVSC):

    v_Ps_MW = m_dcbus[:, DC_PS]
    v_Qs_Mvar = m_dcbus[:, DC_QS]
    v_Ps_pu = v_Ps_MW/DCSbase # initial values
    v_Qs_pu = v_Qs_Mvar/DCSbase # initial values
    v_Ss_pu = v_Ps_pu + j*v_Qs_pu

    v_Pdc_pu = -v_Ps_pu # only to create the vector
    
    v_Us = m_dcbus[:,DC_US] # initial value
    v_Theta_s = m_dcbus[:,DC_DELTAS]
    v_Us_phasor = np.multiply(v_Us, np.exp(j*v_Theta_s))
 
    # Pre-alocate memory
    v_PlossVSC = np.matrix(np.ones(shape=(nVSC, 1)))
    
    return v_Us, v_Theta_s, v_Ps_pu, v_Qs_pu, v_Pdc_pu, v_PlossVSC
   
def fun_setACsetpoints(v_Uref, v_Theta_ref, v_Pref_pu, v_Uref0, v_Qref0_pu, DCSbase, m_converter, m_dcbus, idconv, idowner, isprerun, l_artificialACbus):

    # get list of DC buses that have converters
    l_DCbusVSC = m_converter[:, DC_F_BUS].astype(int)
    
    # indexes of DC buses with VSC
    l_indexDCbusVSC = np.in1d(m_dcbus[:, DC_BUS].astype(int), l_DCbusVSC)
    
    m_dcbusVSC = m_dcbus[l_indexDCbusVSC, :]  # filter data of DC buses with VSC (keep only buses with converters)
    
    # get indices of filtered DC buses according to ac-side gen bus type (qs: type2=1, us: type2=2, slack: type2=3)
    l_indexDCtype1_filtered = list(np.array((np.nonzero(m_dcbusVSC[:,DC_TYPE2]==1)[0]).astype(int)).reshape(-1,)) # make indices a integer-type list
    l_indexDCtype2_filtered = list(np.array((np.nonzero(m_dcbusVSC[:,DC_TYPE2]==2)[0]).astype(int)).reshape(-1,))
    l_indexDCtype3_filtered = list(np.array((np.nonzero(m_dcbusVSC[:,DC_TYPE2]==3)[0]).astype(int)).reshape(-1,))

    # map filtered DC bus indices to converter indices for each type
    l_indexDCtype1 = [np.where(l_DCbusVSC == m_dcbusVSC[i, DC_BUS].astype(int))[0][0] for i in l_indexDCtype1_filtered]
    l_indexDCtype2 = [np.where(l_DCbusVSC == m_dcbusVSC[i, DC_BUS].astype(int))[0][0] for i in l_indexDCtype2_filtered]
    l_indexDCtype3 = [np.where(l_DCbusVSC == m_dcbusVSC[i, DC_BUS].astype(int))[0][0] for i in l_indexDCtype3_filtered]

    # get corresponding AC-side bus 
    l_busACtype1 = list(np.array((m_converter[l_indexDCtype1, FIL_ACBUS]).astype(int)).reshape(-1,))
    l_busACtype2 = list(np.array((m_converter[l_indexDCtype2, FIL_ACBUS]).astype(int)).reshape(-1,))
    l_busACtype3 = list(np.array((m_converter[l_indexDCtype3, FIL_ACBUS]).astype(int)).reshape(-1,))
    
    nVSCslack = len(l_busACtype3) # number of slack-type AC-side converters
    nVSCpv = len(l_busACtype2) # number of PV-type AC-side converters
    nVSCpq = len(l_busACtype1) # number of PQ-type AC-side converters
    
    # Transformer filter (p.u)
    v_Rt = m_converter[:, FIL_RT]
    v_Xt = m_converter[:, FIL_XT]

    # Phase reactor (p.u)
    v_Rc = m_converter[:, FIL_RC]
    v_Xc = m_converter[:, FIL_XC]

    # Power rating of the converters (MW)
    v_MbaseVSC = m_converter[:, FIL_RATE_A]
    
    # PQ Limits of the converters (p.u)
    v_Pcmax = np.matrix(1.0*m_converter[:, FIL_PMAX])
    v_Pcmin = np.matrix(1.0*m_converter[:, FIL_PMIN])
    v_Qcmax = np.matrix(1.0*m_converter[:, FIL_QMAX])
    v_Qcmin = np.matrix(1.0*m_converter[:, FIL_QMIN])
    
    # Converters that control Qs or Qf         
    for i in range(0, nVSCpq):

        # Zsorce is only used in the dynamics  
        Rsorce_pu_conv = 0.0
        Xsorce_pu_conv = 0.1 # any value since not used later on
            
        if not isprerun:

            psspy.purgplnt(l_busACtype1[i]) # Removal of the converters from point s

            # Pf and Qf injection before the transformer (Qmin = Qmax = Qf)
            Rsorce_pu_conv = (v_Rc[l_indexDCtype1[i],0])*v_MbaseVSC[l_indexDCtype1[i],0]/DCSbase # p.u-machine
            Xsorce_pu_conv = (v_Xc[l_indexDCtype1[i],0])*v_MbaseVSC[l_indexDCtype1[i],0]/DCSbase # p.u-machine
            # Zsorce is only used in the dynamics
            
            l_busACtype1[i] = fun_relatedArtificialbuses(l_artificialACbus, l_busACtype1[i])      # Replacement of the main bus with the virtual bus installed at point f
        
        # Include a new generator with Ps (or Pf) and Qs (or Qf)   
        psspy.plant_data(l_busACtype1[i],_i,[ v_Uref[l_indexDCtype1[i],0],_f]) # Control Us or Uf -- it is not necessary           
        psspy.machine_data_2(l_busACtype1[i], idconv,[_i,idowner,_i,_i,_i,_i],[ v_Pref_pu[l_indexDCtype1[i],0]*DCSbase, v_Qref0_pu[l_indexDCtype1[i],0]*DCSbase, v_Qcmax[l_indexDCtype1[i],0],v_Qcmin[l_indexDCtype1[i],0], v_Pcmax[l_indexDCtype1[i],0], v_Pcmin[l_indexDCtype1[i],0],_f,_f, 0.0001,_f,_f,0.0,_f,0.0,0.0,0.0,_f])
        psspy.machine_data_2(l_busACtype1[i], idconv,[_i,_i,_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f, v_MbaseVSC[l_indexDCtype1[i],0],_f,_f,_f,_f,_f,_f,_f,_f,_f,_f]) # de momento pongo potencia nominal del conv = a 100MVA
        psspy.machine_data_2(l_busACtype1[i], idconv,[_i,_i,_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f,_f, Rsorce_pu_conv,_f,_f,_f,_f,_f,_f,_f,_f,_f]) # de momento sin filtro paralelo
        psspy.machine_data_2(l_busACtype1[i], idconv,[_i,_i,_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f,_f,_f, Xsorce_pu_conv,_f,_f,_f,_f,_f,_f,_f,_f])
        psspy.machine_data_2(l_busACtype1[i], idconv,[_i,_i,_i,_i,_i,_i],[_f, _f, v_Qref0_pu[l_indexDCtype1[i],0]*DCSbase, v_Qref0_pu[l_indexDCtype1[i],0]*DCSbase,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f])
        psspy.bus_data_2(l_busACtype1[i],[2,_i,_i,_i],[_f,_f,_f],_s)

    # Converters that control Us or Uf
    for i in range(0, nVSCpv):
    
        # Zsorce is only used in the dynamics  
        Rsorce_pu_conv = 0.0
        Xsorce_pu_conv = 0.1

        if not isprerun:

            psspy.purgplnt(l_busACtype2[i]) # Removal of the converters from point s

            # Pf injection after the transformer and control voltage at point f
            Rsorce_pu_conv = (v_Rc[l_indexDCtype2[i],0])*v_MbaseVSC[l_indexDCtype2[i],0]/DCSbase # p.u-machine
            Xsorce_pu_conv = (v_Xc[l_indexDCtype2[i],0])*v_MbaseVSC[l_indexDCtype2[i],0]/DCSbase # p.u-machine
            # Zsorce is only used in the dynamics
            
            l_busACtype2[i] = fun_relatedArtificialbuses(l_artificialACbus, l_busACtype2[i])      # Replacement of the main bus with the virtual bus installed at point f

        # Include a new generator with Ps (or Pf) and Us (or Uf)
        psspy.plant_data(l_busACtype2[i],_i,[ v_Uref0[l_indexDCtype2[i],0],_f]) # Control Us or Uf
        psspy.machine_data_2(l_busACtype2[i], idconv,[_i,idowner,_i,_i,_i,_i],[ v_Pref_pu[l_indexDCtype2[i],0]*DCSbase, _f, v_Qcmax[l_indexDCtype2[i],0], v_Qcmin[l_indexDCtype2[i],0], v_Pcmax[l_indexDCtype2[i],0], v_Pcmin[l_indexDCtype2[i],0],_f,_f, 0.0001,_f,_f,0.0,_f,0.0,0.0,0.0,_f])
        psspy.machine_data_2(l_busACtype2[i], idconv,[_i,_i,_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f, v_MbaseVSC[l_indexDCtype2[i],0],_f,_f,_f,_f,_f,_f,_f,_f,_f,_f]) # de momento pongo potencia nominal del conv = a 100MVA
        psspy.machine_data_2(l_busACtype2[i], idconv,[_i,_i,_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f,_f, Rsorce_pu_conv,_f,_f,_f,_f,_f,_f,_f,_f,_f]) # de momento sin filtro paralelo
        psspy.machine_data_2(l_busACtype2[i], idconv,[_i,_i,_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f,_f,_f, Xsorce_pu_conv,_f,_f,_f,_f,_f,_f,_f,_f])
        psspy.bus_data_2(l_busACtype2[i],[2,_i,_i,_i],[_f,_f,_f],_s) # say to psse that the bus is a PV node (code = 2)

    # Converters that control Us or Uf (only if the converter is feeding a passive grid)
    for i in range(0, nVSCslack):

        # Zsorce is only used in the dynamics  
        Rsorce_pu_conv = 0.0
        Xsorce_pu_conv = 0.1
        
        if not isprerun:

            psspy.purgplnt(l_busACtype3[i]) # Removal of the converters from point s
            
            # Pf injection after the transformer and control voltage at point f
            Rsorce_pu_conv = (v_Rc[l_indexDCtype3[i],0])*v_MbaseVSC[l_indexDCtype3[i],0]/DCSbase # p.u-machine
            Xsorce_pu_conv = (v_Xc[l_indexDCtype3[i],0])*v_MbaseVSC[l_indexDCtype3[i],0]/DCSbase # p.u-machine
            # Zsorce is only used in the dynamics

            l_busACtype3[i] = fun_relatedArtificialbuses(l_artificialACbus, l_busACtype3[i])      # Replacement of the main bus with the virtual bus installed at point f

        # Include a new generator with Ps (or Pf) and Us (or Uf)
        psspy.plant_data(l_busACtype3[i],_i,[ v_Uref0[l_indexDCtype3[i],0],_f]) # Control Us or Uf
        psspy.machine_data_2(l_busACtype3[i], idconv,[_i,idowner,_i,_i,_i,_i],[ v_Pref_pu[l_indexDCtype3[i],0]*DCSbase, _f, v_Qcmax[l_indexDCtype3[i],0], v_Qcmin[l_indexDCtype3[i],0], v_Pcmax[l_indexDCtype3[i],0], v_Pcmin[l_indexDCtype3[i],0],_f,_f, 0.0001,_f,_f,0.0,_f,0.0,0.0,0.0,_f])
        psspy.machine_data_2(l_busACtype3[i], idconv,[_i,_i,_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f, v_MbaseVSC[l_indexDCtype3[i],0],_f,_f,_f,_f,_f,_f,_f,_f,_f,_f]) # de momento pongo potencia nominal del conv = a 100MVA
        psspy.machine_data_2(l_busACtype3[i], idconv,[_i,_i,_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f,_f, Rsorce_pu_conv,_f,_f,_f,_f,_f,_f,_f,_f,_f]) # de momento sin filtro paralelo
        psspy.machine_data_2(l_busACtype3[i], idconv,[_i,_i,_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f,_f,_f, Xsorce_pu_conv,_f,_f,_f,_f,_f,_f,_f,_f])
        psspy.bus_data_2(l_busACtype3[i],[3,_i,_i,_i],[_f,_f,0.0],_s) # say to psse that the bus is a ac-slack node (code = 2)

    return l_indexDCbusVSC

def fun_backwardcomputation(m_converter, m_dcbus, nVSC, v_Ps_pu, v_Qs_pu, DCSbase, l_indexDCbusVSC):

    m_dcbus_f = m_dcbus  # Resetting the new settings due to solving the AC/DC power flow at point f

    # Transformer filter (p.u)
    v_Rt = m_converter[:, FIL_RT]
    v_Xt = m_converter[:, FIL_XT]
    v_Zt = v_Rt + j*v_Xt

    # Low pass filter (p.u)
    v_Bf = m_converter[:, FIL_BF]
    
    # Zf = -j/v_Bf # division term by term
    v_Yf = j*v_Bf

    t_AC_bus = m_converter[:, FIL_ACBUS]
 
    # Backward Converter calculations 
    # converter reactor voltages and power
    v_Us = []
    v_Theta_s = []
    for i in range(0, nVSC):
       
        ierr, voltage_s = psspy.busdat(t_AC_bus[i], "PU")     # Extracting the voltage magnitude at point s [p.u.] (API)
        ierr, angle_s = psspy.busdat(t_AC_bus[i], "ANGLE")    # Extracting the voltage phase at point s [rad] (API)
        
        v_Us.append(voltage_s)
        v_Theta_s.append(angle_s)
    
    v_Us = np.reshape(np.matrix(v_Us), (-1, 1))
    v_Theta_s = np.reshape(np.matrix(v_Theta_s), (-1, 1))
    v_Us_phasor = np.multiply(v_Us, np.exp(j*v_Theta_s))      # Calculating the voltage phasor at point s

    v_Itf = np.conj((v_Ps_pu[l_indexDCbusVSC] +j*v_Qs_pu[l_indexDCbusVSC])/ v_Us_phasor)        # Calculating transformer current
  
    v_Uf_phasor = v_Us_phasor + np.multiply(v_Itf, v_Zt)      # Calculating the voltage phasor at point f (filter side) 
    v_Uf = abs(v_Uf_phasor)
    v_Theta_f = np.angle(v_Uf_phasor)

    v_Qf_pu = np.multiply(v_Bf, np.multiply(np.abs(v_Uf_phasor), np.abs(v_Uf_phasor)))    # Calculating filter generated reactive power (PSS/e sign convention for shunts)

    # converter apparent, active and reactive powers
    v_So_pu = np.multiply(v_Uf_phasor, np.conj(v_Itf))      # Calculating filter side transformer complex power (Converter Output Power, Gen & Zc)
    v_Po_pu = v_So_pu.real
    v_Qo_pu = v_So_pu.imag - v_Qf_pu

    v_Pdc_pu = -v_Po_pu # only to create the vector ////////////// In general, Not true because: Pdc = Po - Ploss
    
    # Update the dcbus matrix with new settings at point f and reflect changes in the dcbus_f matrix
    m_dcbus_f[l_indexDCbusVSC,DC_PF] = np.array(v_Po_pu).reshape(-1)*DCSbase
    m_dcbus_f[l_indexDCbusVSC,DC_UF] = np.array(v_Uf).reshape(-1)
    m_dcbus_f[l_indexDCbusVSC,DC_DELTAF] = np.array(v_Theta_f*180/pi).reshape(-1)
    m_dcbus_f[l_indexDCbusVSC,DC_QF] = np.array(v_Qo_pu*DCSbase).reshape(-1)
    
    v_QLosses_tf_pu = np.multiply(v_Xt, np.multiply(np.abs(v_Itf), np.abs(v_Itf))) # Calculating reactive power losses in the transformer (p.u.)
    
    return v_Uf, v_Theta_f, v_Uf_phasor, v_Po_pu, v_Qo_pu, v_Qf_pu, v_QLosses_tf_pu, v_Pdc_pu, m_dcbus_f

def fun_addnewelements(nVSC, maxnACbus, DCSbase, idconv, m_converter, v_Uf0, v_Theta_f0, v_Qf_pu):

    v_Rt = m_converter[:, FIL_RT]
    v_Xt = m_converter[:, FIL_XT]

    # Add an artificial bus f, a shunt filter, and a transformer connected between the AC bus and the artificial bus f
    l_artificialACbus = []
    for i in range(0, nVSC):
    
        artificalACbus = maxnACbus + i + 1

        ierr, VBASE_KV = psspy.busdat(int(m_converter[i, FIL_ACBUS]), 'BASE')   # Extracting the voltage base of buses at point s to be used in the artificial bus at point f
        psspy.bus_data_4(artificalACbus, 0,[_i,_i,_i,_i],[VBASE_KV, v_Uf0[i], v_Theta_f0[i]*(180/pi),_f,_f,_f,_f],"BUS" + str(artificalACbus))
        psspy.shunt_data(artificalACbus, idconv,_i,[_f, v_Qf_pu[i]*DCSbase])
        psspy.two_winding_data_6(int(m_converter[i, FIL_ACBUS]),artificalACbus, idconv,[_i,_i,_i,_i,_i,_i,_i,_i, int(m_converter[i, FIL_ACBUS]),_i,_i,_i, 0,_i,_i,_i],[v_Rt[i], v_Xt[i],_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f],[_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f],"","")
         
        l_artificialACbus.append([int(m_converter[i, FIL_ACBUS]), artificalACbus])

    return l_artificialACbus

def fun_relatedArtificialbuses(l_artificialACbus, l_busACtype):
    
    value = next((pair[1] for pair in l_artificialACbus if pair[0] == l_busACtype), None)

    return value

def fun_getACresults(v_Uf, v_Theta_f, v_Po_pu, v_Qo_pu, DCSbase, m_converter, m_dcbus_f, idconv, l_artificialACbus):

    # get list of DC buses that have converters
    l_DCbusVSC = m_converter[:, DC_F_BUS].astype(int)

    # create mask to filter only DC buses that have converters
    l_indexDCbusVSC = np.in1d(m_dcbus_f[:, DC_BUS].astype(int), l_DCbusVSC)
    
    filtered_m_dcbus_f = m_dcbus_f[l_indexDCbusVSC, :]  # filter DC bus data by mask (keep only buses with converters)
    
    # get indices of filtered DC buses according to ac-side gen bus type (qs: type2=1, us: type2=2, slack: type2=3)
    l_indexDCtype1_filtered = list(np.array((np.nonzero(filtered_m_dcbus_f[:,DC_TYPE2]==1)[0]).astype(int)).reshape(-1,)) # make indices a integer-type list
    l_indexDCtype2_filtered = list(np.array((np.nonzero(filtered_m_dcbus_f[:,DC_TYPE2]==2)[0]).astype(int)).reshape(-1,))
    l_indexDCtype3_filtered = list(np.array((np.nonzero(filtered_m_dcbus_f[:,DC_TYPE2]==3)[0]).astype(int)).reshape(-1,))

    # map filtered DC bus indices to converter indices for each type
    l_indexDCtype1 = [np.where(l_DCbusVSC == filtered_m_dcbus_f[i, DC_BUS].astype(int))[0][0] for i in l_indexDCtype1_filtered]
    l_indexDCtype2 = [np.where(l_DCbusVSC == filtered_m_dcbus_f[i, DC_BUS].astype(int))[0][0] for i in l_indexDCtype2_filtered]
    l_indexDCtype3 = [np.where(l_DCbusVSC == filtered_m_dcbus_f[i, DC_BUS].astype(int))[0][0] for i in l_indexDCtype3_filtered]

    # get corresponding AC-side bus 
    l_busACtype1 = list(np.array((m_converter[l_indexDCtype1, FIL_ACBUS]).astype(int)).reshape(-1,))
    l_busACtype2 = list(np.array((m_converter[l_indexDCtype2, FIL_ACBUS]).astype(int)).reshape(-1,))
    l_busACtype3 = list(np.array((m_converter[l_indexDCtype3, FIL_ACBUS]).astype(int)).reshape(-1,))
       
    nVSCslack = len(l_busACtype3) # number of slack-type AC-side converters
    nVSCpv = len(l_busACtype2) # number of PV-type AC-side converters
    nVSCpq = len(l_busACtype1) # number of PQ-type AC-side converters
            
    # Converters that control Qf
    for i in range(0, nVSCpq):
        
        l_busACtype1[i] = fun_relatedArtificialbuses(l_artificialACbus, l_busACtype1[i])    # Replacement of the main bus with the virtual bus installed at point f
        
        ierr, v_Uf[l_indexDCtype1[i],0] = psspy.busdat(l_busACtype1[i], 'PU')         # Voltage magnitude at ponit f [p.u.] (API)
        ierr, v_Theta_f[l_indexDCtype1[i],0] = psspy.busdat(l_busACtype1[i], 'ANGLE') # Voltage phase at ponit f [rad] (API)
        
        ierr, rval = psspy.macdat(l_busACtype1[i], idconv, 'P') # Pgen [MW] (API)
        v_Po_pu[l_indexDCtype1[i],0] = rval/DCSbase
        ierr, rval = psspy.macdat(l_busACtype1[i], idconv, 'Q') # Qgen [MVAr] (API)
        v_Qo_pu[l_indexDCtype1[i],0] = rval/DCSbase
        
    # Converters that control Uf
    for i in range(0, nVSCpv):

        l_busACtype2[i] = fun_relatedArtificialbuses(l_artificialACbus, l_busACtype2[i])    # Replacement of the main bus with the virtual bus installed at point f

        ierr, v_Uf[l_indexDCtype2[i],0] = psspy.busdat(l_busACtype2[i], 'PU')          # Voltage magnitude at ponit f [p.u.] (API)
        ierr, v_Theta_f[l_indexDCtype2[i],0] = psspy.busdat(l_busACtype2[i], 'ANGLE')  # Voltage phase at ponit f [rad] (API)
        
        ierr, rval = psspy.macdat(l_busACtype2[i], idconv, 'P') # Pgen [MW] (API)
        v_Po_pu[l_indexDCtype2[i],0] = rval/DCSbase
        ierr, rval = psspy.macdat(l_busACtype2[i], idconv, 'Q') # Qgen [MVAr] (API)
        v_Qo_pu[l_indexDCtype2[i],0] = rval/DCSbase

    # ac-slack converter
    for i in range(0, nVSCslack):

        l_busACtype3[i] = fun_relatedArtificialbuses(l_artificialACbus, l_busACtype3[i])    # Replacement of the main bus with the virtual bus installed at point f

        ierr, v_Uf[l_indexDCtype3[i],0] = psspy.busdat(l_busACtype3[i], 'PU')           # Voltage magnitude at ponit f [p.u.] (API)
        ierr, v_Theta_f[l_indexDCtype3[i],0] = psspy.busdat(l_busACtype3[i], 'ANGLE')   # Voltage phase at ponit f [rad] (API)
        
        ierr, rval = psspy.macdat(l_busACtype3[i], idconv, 'P') # Pgen [MW] (API)
        v_Po_pu[l_indexDCtype3[i],0] = rval/DCSbase
        ierr, rval = psspy.macdat(l_busACtype3[i], idconv, 'Q') # Qgen [MVAr] (API)
        v_Qo_pu[l_indexDCtype3[i],0] = rval/DCSbase 
    
    return v_Uf, v_Theta_f, v_Po_pu, v_Qo_pu

def fun_calculateVSCtransformerlosses(l_artificialACbus, DCSbase):

    ACbuses_Artificialbuses_array = np.array(l_artificialACbus)
    all_branches_strings = ['fromnumber', 'tonumber']                    
    ierr, idata = psspy.aflowint(-1, 1, 5, 1,all_branches_strings)     # Extraction of the list of sending (From) and receiving (To) buses
    
    from_buses = np.array(idata[0])  
    to_buses = np.array(idata[1])    
    
    VSC_tr_branch_indices = np.concatenate([np.where((from_buses == branch[0]) & (to_buses == branch[1]))[0] for branch in ACbuses_Artificialbuses_array]) # Identification of the VSCs transformer branch number
    
    [ierr, rarray_PLosses_tf] = psspy.aflowreal(-1, 1, 1, 1, 'PLOSS')  # Active power losses of all branches [MW] (API)  
    [ierr, rarray_QLosses_tf] = psspy.aflowreal(-1, 1, 1, 1, 'QLOSS')  # Reactive power losses of all branches [MVAr] (API)  
    
    v_PLosses_tf_MW = np.array([rarray_PLosses_tf[0][i] for i in VSC_tr_branch_indices]).reshape(-1, 1)  # Active power losses of VSCs transformer [MW]
    v_PLosses_tf_pu = v_PLosses_tf_MW/DCSbase
    
    v_QLosses_tf_MVAr = np.array([rarray_QLosses_tf[0][i] for i in VSC_tr_branch_indices]).reshape(-1, 1)  # Reactive power losses of VSCs transformer [MVAr]
    v_QLosses_tf_pu = v_QLosses_tf_MVAr/DCSbase
    
    return v_PLosses_tf_pu, v_QLosses_tf_pu

def fun_initializeDCfromAC(v_Uf, v_Theta_f, v_Po_pu, v_Qo_pu, v_Pdc_pu, v_PlossVSC, v_Am, v_Bm, v_Cm, v_Dm,
    v_aloss, v_bloss, v_crect, v_cinv):

    v_Uf_phasor = np.multiply(v_Uf, np.exp(j*v_Theta_f))
    v_So_pu = v_Po_pu + j*v_Qo_pu
    v_If_phasor = np.conj(v_So_pu/v_Uf_phasor)
    v_Uc_phasor = np.multiply(v_Am, v_Uf_phasor) + np.multiply(v_Bm, v_If_phasor) # v_Uc_phasor = v_Uf_phasor + np.multiply(v_Zt,v_If_phasor)
    v_Ic_phasor = np.multiply(v_Cm, v_Uf_phasor) + np.multiply(v_Dm, v_If_phasor) # v_Ic_phasor = v_Is_phasor + np.multiply(v_Yf,v_Uf_phasor)
    v_Sc_pu = np.multiply( v_Uc_phasor, np.conj(v_Ic_phasor) )
    v_Pc_pu = np.real(v_Sc_pu)
    v_Qc_pu = np.imag(v_Sc_pu)

    # Losses
    ind_rect = list(np.array((np.nonzero(v_Pc_pu<0)[0]).astype(int)).reshape(-1,)) # integer-type list
    ind_inv = list(np.array((np.nonzero(v_Pc_pu>=0)[0]).astype(int)).reshape(-1,))
    if len(ind_rect)>0:
        v_PlossVSC[ind_rect] = v_aloss[ind_rect] + np.multiply( v_bloss[ind_rect], abs(v_Ic_phasor[ind_rect]) ) \
                          + np.multiply( v_crect[ind_rect], np.power(abs(v_Ic_phasor[ind_rect]),2) )
        v_Pdc_pu[ind_rect] = (abs(v_Pc_pu[ind_rect]) - v_PlossVSC[ind_rect])      
    if len(ind_inv)>0:
        v_PlossVSC[ind_inv] = v_aloss[ind_inv] + np.multiply( v_bloss[ind_inv], abs(v_Ic_phasor[ind_inv]) ) \
                         + np.multiply( v_crect[ind_inv], np.power(abs(v_Ic_phasor[ind_inv]),2) )
        v_Pdc_pu[ind_inv] = -(abs(v_Pc_pu[ind_inv]) + v_PlossVSC[ind_inv])
    
    return v_If_phasor, v_Ic_phasor, v_Uc_phasor, v_Pc_pu, v_Qc_pu, v_Pdc_pu, v_PlossVSC  
   
def fun_calldclf(v_Uf, v_Theta_f, v_Uc_phasor, v_Po_pu, v_Qo_pu, v_Pc_pu, v_Qc_pu, v_Pdc_pu, DCSbase, 
    m_dcbus_f, m_dcbranch, v_Y1, v_Z2, v_Y3, v_aloss, v_bloss, v_crect, v_cinv, l_MTDCgrids, v_indexDCslack, 
    lftol, lfmaxiter, l_MTDCnetworkid, l_indexDCbusVSC):
    
    # iteration and mismatch parameters
    k_int_dc = []
    v_mismatch_dc = []
    k_int_dcslack = []
    j_int_dcslack = []
    v_Po_no_current = []
    
    nMTDCgrids = len(l_MTDCgrids)

    ntotDCbus = sum([len(grid["dcbus"]) for grid in l_MTDCgrids])
    ntotDCbranch = sum([len(grid["dcbranch"]) for grid in l_MTDCgrids])
    ntotVSC = sum([len(grid["converter"]) for grid in l_MTDCgrids])

    # 2. Run DC-side power flow
    m_dcbuses = m_dcbus_f[:ntotDCbus, :]
    
    m_dcbranches = m_dcbranch[:ntotDCbranch, :]
    
    [m_dcbuses, m_dcbranches, k_int_dc_dcpf, mismatch_dci] = fun_mainDCloadflow(DCSbase, lftol, lfmaxiter, m_dcbuses, m_dcbranches, l_MTDCnetworkid)
    Pdc_dcpf = m_dcbuses[l_indexDCbusVSC, DC_PDC] / DCSbase
    
    # 3. Run slack-bus iteration
    Uf_dcpf = v_Uf
    Th_f_dcpf = v_Theta_f
    Po_dcpf = v_Po_pu
    Qo_dcpf = v_Qo_pu
    Pc_dcpf = v_Pc_pu
    Qc_dcpf = v_Qc_pu
    Uc_phasor_dcpf = v_Uc_phasor
    
    # get DC slack    
    n_dcslack_dcpf = v_indexDCslack
    
    po_n = Po_dcpf[n_dcslack_dcpf]
    po_n = po_n.astype(float)
    pc_n = Pc_dcpf[n_dcslack_dcpf]
    pc_n = pc_n.astype(float)
    qc_n = Qc_dcpf[n_dcslack_dcpf]
    qc_n = qc_n.astype(float)
    uf_n = Uf_dcpf[n_dcslack_dcpf]
    th_n = Th_f_dcpf[n_dcslack_dcpf]
    uc_n = abs(Uc_phasor_dcpf[n_dcslack_dcpf])
    thc_n = np.angle(Uc_phasor_dcpf[n_dcslack_dcpf])
    qo_n = Qo_dcpf[n_dcslack_dcpf]
    pdc_n = Pdc_dcpf[n_dcslack_dcpf]
    po_n = Po_dcpf[n_dcslack_dcpf]
    
    y1_n = v_Y1[n_dcslack_dcpf]
    z2_n = v_Z2[n_dcslack_dcpf]
    y3_n = v_Y3[n_dcslack_dcpf]
    a_n = v_aloss[n_dcslack_dcpf]
    b_n = v_bloss[n_dcslack_dcpf]
    c_inv_n = v_cinv[n_dcslack_dcpf]
    c_rect_n = v_crect[n_dcslack_dcpf]
    
    # DC slack iteration
    [Po_dcpf[n_dcslack_dcpf], Pc_dcpf[n_dcslack_dcpf], Qc_dcpf[n_dcslack_dcpf], 
    Uc_phasor_dcpf[n_dcslack_dcpf], k_int_dcslack_dcpf, j_int_dcslack_dcpf, convergence] = fun_slackiterationDC(uf_n, 
    th_n, qo_n, pdc_n, po_n, y1_n, z2_n, y3_n, a_n, b_n, c_inv_n, c_rect_n, lftol, lfmaxiter)
    
    # update of dcbus matrix with the AC-side information of the converter
    m_dcbuses[l_indexDCbusVSC,DC_PF] = np.array(Po_dcpf).reshape(-1)*DCSbase # put the new Ps in the dcbus matrix 
    m_dcbuses[l_indexDCbusVSC,DC_UF] = np.array(Uf_dcpf).reshape(-1)
    m_dcbuses[l_indexDCbusVSC,DC_PDC] = np.array(Pdc_dcpf).reshape(-1)
    m_dcbuses[l_indexDCbusVSC,DC_DELTAF] = np.array(Th_f_dcpf*180/pi).reshape(-1)
    m_dcbuses[l_indexDCbusVSC,DC_QF] = np.array(Qo_dcpf).reshape(-1)*DCSbase

    # 4. Update variables for outer AC power flow        
    m_dcbus_f = m_dcbuses
    v_Pdc_pu = Pdc_dcpf
    v_Po_pu = Po_dcpf
    v_Pc_pu = Pc_dcpf
    v_Qc_pu = Qc_dcpf
    v_Uc_phasor = Uc_phasor_dcpf
    m_dcbranch = m_dcbranches
    
    k_int_dc.append(k_int_dc_dcpf)
    v_mismatch_dc.append(mismatch_dci)
    k_int_dcslack.append(k_int_dcslack_dcpf)
    j_int_dcslack.append(j_int_dcslack_dcpf)

    v_Po_no_current = Po_dcpf[n_dcslack_dcpf, 0].astype(float)
    
    return v_Po_no_current, k_int_dc, v_mismatch_dc, k_int_dcslack, j_int_dcslack, v_Po_pu, v_Pc_pu, v_Qc_pu, v_Pdc_pu, v_Uc_phasor, m_dcbus_f, m_dcbranch
  
def fun_getacdcoutput(v_If_phasor, v_Uf, v_Theta_f, v_Uc_phasor, v_Ic_phasor, v_Pc_pu, v_Qc_pu, v_PlossVSC, DCSbase, m_dcbus_f, m_dcbranch, m_converter, l_artificialACbus, l_indexDCbusVSC):

    # Low pass filter current
    v_Bf = m_converter[:, FIL_BF]
    v_Yf = j*v_Bf
    v_Uf_phasor = np.multiply(v_Uf, np.exp(j*v_Theta_f))

    v_Is_phasor = v_If_phasor + np.conj(np.dot(v_Yf.T, v_Uf_phasor))  # Calculation of VSCs transformer current: Is = If + Yf*Uf
  
    # get losses   >>>>>>>>>>>>>>>>>>>>>> The following presents transformer losses calculations using the API for improved comparison
    v_Rt = m_converter[:, FIL_RT]
    v_Rc = m_converter[:, FIL_RC]
    v_Ploss_tf = np.multiply(v_Rt, abs(np.power(v_Is_phasor, 2)))  
    v_Ploss_filter = np.multiply(v_Rc, abs(np.power(v_Ic_phasor, 2)));   
    
    # get DC variables
    m_converter = np.array(m_converter)
    m_dcbus_f = np.array(m_dcbus_f)
    m_dcbranch = np.array(m_dcbranch)
    
    # get AC variables
    sid = -1 # all buses
    flag = 1 # for only in-service buses
    [ierr, v_Um_bus] = psspy.abusreal(sid, flag, 'PU')        # List of voltage magnitudes for the in-service buses in the system [p.u.] (API)
    [ierr, v_Theta_bus] = psspy.abusreal(sid, flag, 'ANGLE')  # List of voltage angles for the in-service buses in the system [rad] (API)
    v_Um_bus = np.transpose(np.matrix(v_Um_bus))              # column vectors (type matrix)
    v_Theta_bus = np.transpose(np.matrix(v_Theta_bus)) 
    [ierr, rarray_P] = psspy.aflowreal(-1, 1, 1, 1, 'P')      # List of active power flow for the branches [MW] (API)  
    Pf_ac_MW = np.matrix(rarray_P)
    Pf_ac_MW = np.transpose(np.matrix(Pf_ac_MW))
    [ierr, rarray_Q] = psspy.aflowreal(-1, 1, 1, 1, 'Q')      # List of reactive power flow for the branches [MVAr] (API)
    Qf_ac_MW = np.matrix(rarray_Q)
    Qf_ac_MW = np.transpose(np.matrix(Qf_ac_MW))
    
    # Extracting S-point values
    # Power setpoints
    v_Po_MW = m_dcbus_f[l_indexDCbusVSC,DC_PF]
    v_Qo_MVAr = m_dcbus_f[l_indexDCbusVSC,DC_QF]
    
    [v_PLosses_tf_pu, v_QLosses_tf_pu] = fun_calculateVSCtransformerlosses(l_artificialACbus, DCSbase) # Extraction of VSCs transformer losses after running the power flow in PSS/E.

    v_PLosses_tf_MW = (v_PLosses_tf_pu * DCSbase).reshape(-1)
    v_QLosses_tf_MVAr = (v_QLosses_tf_pu * DCSbase).reshape(-1)
    
    v_Ps_MW = (v_Po_MW - v_PLosses_tf_MW).reshape(-1)
    v_Qs_MVAr = (v_Qo_MVAr - v_QLosses_tf_MVAr).reshape(-1)

    # Voltage setpoints
    ACbuses_Artificialbuses_array = np.array(l_artificialACbus)
    ierr, l_acbus = psspy.abusint(sid, flag, 'NUMBER')       # For determining AC buses and its bus No. (SID<0 & FLAG = 1: only in-service buses
    l_acbus = np.array(l_acbus).flatten()
    t_AC_bus_indices = np.searchsorted(l_acbus, ACbuses_Artificialbuses_array[:, 0])   # Identification of the indices of AC buses connected to the VSCs
    
    v_Us = [v_Um_bus[i, 0] for i in t_AC_bus_indices]                          # Voltage magnitude at point s [p.u.]
    v_Theta_s = [(v_Theta_bus[i, 0] * (180/np.pi)) for i in t_AC_bus_indices]  # Voltage phase at point s [deg]
    
    # If we want to see f point results, we can deactive this section and change m_dcbus to m_dcbus_f in line (809 or A++)
    ##################################################################
    m_dcbus = m_dcbus_f
    m_dcbus[l_indexDCbusVSC,DC_PS] = v_Ps_MW.reshape(-1)
    m_dcbus[l_indexDCbusVSC,DC_US] = v_Us
    m_dcbus[l_indexDCbusVSC,DC_DELTAS] = v_Theta_s
    m_dcbus[l_indexDCbusVSC,DC_QS] = v_Qs_MVAr.reshape(-1)
    ##################################################################
    d_acdcoutput = {"version": '2'}
    d_acdcoutput["DCSbase"] = DCSbase
    d_acdcoutput["acbus"] = np.concatenate((v_Um_bus, v_Theta_bus*180/pi),1)
    d_acdcoutput["acbranch"] = np.concatenate((Pf_ac_MW, Qf_ac_MW),1)
    d_acdcoutput["vsc"] = np.concatenate((abs(v_Uc_phasor), np.angle(v_Uc_phasor), v_Pc_pu, v_Qc_pu),1)
    d_acdcoutput["converter"] = m_converter
    d_acdcoutput["dcbus"] = m_dcbus            # based on s or f set points // Command line A++;
    d_acdcoutput["dcbranch"] = m_dcbranch
    d_acdcoutput["losses"] = np.concatenate((v_Ploss_tf, v_Ploss_filter, v_PlossVSC),1)
    
    return d_acdcoutput
    
def fun_showacdclfoutput(isprogress, issavedclfresults, str_pathdclfresults, d_acdcoutput):
    
    if isprogress!=6:
        DCSbase = d_acdcoutput.get("DCSbase")
        m_dcbus = d_acdcoutput.get("dcbus")
        m_dcbranch = d_acdcoutput.get("dcbranch")
        m_losses = d_acdcoutput.get("losses")
        m_vsc = d_acdcoutput.get("vsc")
        v_Ploss_tf = m_losses[:,0]
        v_Ploss_filter = m_losses[:,1]
        v_PlossVSC = m_losses[:,2]
        v_Uc = m_vsc[:,0]
        v_Theta_c = m_vsc[:,1]
        v_Pc_pu = m_vsc[:,2]
        v_Qc_pu = m_vsc[:,3]
        Ploss_coup = v_Ploss_tf + v_Ploss_filter + v_PlossVSC; # losses transformer + reactor + converter
        
        np.set_printoptions(precision=8, linewidth=200, suppress=True)
        print ('===================================================================================================')
        print ('DC grid - dcbus')
        print ('===================================================================================================')
        print ('---------------------------------------------------------------------------------------------------')
        print ('| dcbus | dctype1 | dctype2 | us (pu) | deltas (deg) | Ps (MW) | Qs (Mvar) |  udc (pu) | Pdc (MW) |')
        print ('---------------------------------------------------------------------------------------------------')
        print (m_dcbus[:,DC_BUS:DC_PDC+1]) # [DC_BUS DC_TYPE DC_TYPE2 DC_US DC_DELTAS DC_PS DC_QS DC_UDC DC_PDC])
        print ('===================================================================================================')
        print ('DC grid - dcbranch')
        print ('===================================================================================================')
        print ('---------------------------------------------------------------------------------------------------')
        print ('| dc From bus | dc To bus  | Pdc,ij (MW)| Pdc,ji (MW) | Icc,ij (pu) ')
        print ('---------------------------------------------------------------------------------------------------')
        print (m_dcbranch[:,(DC_F_BUS, DC_T_BUS, DC_PIJ, DC_PJI, DC_ICCIJ)])
        print ('===================================================================================================')
        print ('VSC converters - losses (MW) and AC/DC coupling')
        print ('===================================================================================================')
        print ('---------------------------------------------------------------------------------------------------')
        print ('| transformer | reactor  | converter| total -Coupling AC/DC| ')
        print ('---------------------------------------------------------------------------------------------------')
        print (np.concatenate((v_Ploss_tf*DCSbase, v_Ploss_filter*DCSbase, v_PlossVSC*DCSbase, Ploss_coup*DCSbase),axis=1))
        print ('===================================================================================================')
        print ('VSC converters - Internal voltage')
        print ('===================================================================================================')
        print ('---------------------------------------------------------------------------------------------------')
        print ('| ec (p.u) | deltac (deg) | Pc (MW) | Qc (Mvar) | ')
        print ('---------------------------------------------------------------------------------------------------')
        print (np.concatenate((v_Uc, v_Theta_c*180/pi, v_Pc_pu*DCSbase, v_Qc_pu*DCSbase),axis=1))
        print ('===================================================================================================')
        print ('MisMatches')
        print ('===================================================================================================')
    
    if issavedclfresults:
        with open(os.path.join(str_pathdclfresults, 'result_dcbus.csv'), 'w') as f_Rdcbus:
            np.savetxt(f_Rdcbus, m_dcbus[:,DC_BUS:DC_PDC+1], delimiter=',')
        f_Rdcbus.close()
        with open(os.path.join(str_pathdclfresults, 'result_dcbranch.csv'), 'w') as f_Rdcbranch:
            np.savetxt(f_Rdcbranch, m_dcbranch[:,(DC_F_BUS, DC_T_BUS, DC_PIJ, DC_PJI, DC_ICCIJ)], delimiter=',')
        f_Rdcbranch.close()
        with open(os.path.join(str_pathdclfresults, 'result_losses.csv'), 'w') as f_Rlosses:
            np.savetxt(f_Rlosses, v_PlossVSC*DCSbase, delimiter=',')
        f_Rlosses.close()
        with open(os.path.join(str_pathdclfresults, 'result_coupling.csv'), 'w') as f_Rcoupling:
            np.savetxt(f_Rcoupling, np.concatenate((v_Ploss_tf*DCSbase, v_Ploss_filter*DCSbase, v_PlossVSC*DCSbase, Ploss_coup*DCSbase),axis=1), delimiter=',')
        f_Rcoupling.close()
        with open(os.path.join(str_pathdclfresults, 'result_voltage.csv'), 'w') as f_Rvoltage:
            np.savetxt(f_Rvoltage, np.concatenate((v_Uc, v_Theta_c*180/pi, v_Pc_pu*DCSbase, v_Qc_pu*DCSbase),axis=1), delimiter=',')
        f_Rvoltage.close()

def main_setacdcdynamicdata(d_acdcoutput, str_dyrfileorig, str_path4dynamics, d_acdcoptions, l_MTDCgrids, l_indexDCbusVSC, l_artificialACbus):

    if d_acdcoptions["issetdynamicfiles"]==True:
    
        idconv = d_acdcoptions.get('idconv') 
        
        # write to text files DC grid data
        fun_setDCnetworkdata(str_path4dynamics, d_acdcoutput, l_MTDCgrids, l_indexDCbusVSC, l_artificialACbus)
       
        # update .dyr file with VSC and DC grid models
        fun_updatedynamicfile(str_dyrfileorig, idconv, d_acdcoutput, l_MTDCgrids, l_indexDCbusVSC, l_artificialACbus)
   
def fun_updatedynamicfile(str_dyrfileorig, idconv, d_acdcoutput, l_MTDCgrids, l_indexDCbusVSC, l_artificialACbus):
    """
    This function updates the existing .dyr file by adding MTDC network and converter and control models with pre-set parameters
    """
        
    # modified file name and copy original file into it
    str_dyrfilemtdc = str_dyrfileorig[:-4] + r"_MTDCg.dyr" # update AC power flow (VSC as generators)
    shutil.copy2(str_dyrfileorig, str_dyrfilemtdc)
    
    DCSbase = d_acdcoutput.get("DCSbase")
    m_dcbus = d_acdcoutput.get("dcbus")
    m_dcbranch = d_acdcoutput.get("dcbranch") 
    m_converter = d_acdcoutput["converter"]
           
    nDClines = np.shape(m_dcbranch[:,0])[0] # # of branches of the DC grid
    nDCbus = np.shape(m_dcbus[:,0])[0] # # of nodes of the DC grid
    nVSC = np.shape(m_converter[:,0])[0] # # of VSCs of the DC grid
    l_nDClines = []
    l_nDCbus = []
    l_nVSC = []
    for i in range(0,len(l_MTDCgrids)):
        a_dcbus = l_MTDCgrids[i]["dcbus"]
        l_nDClines.append(len(l_MTDCgrids[i]["dcbranch"]))
        l_nDCbus.append(len(a_dcbus))
        l_nVSC.append(len(l_MTDCgrids[i]["converter"]))
        
    # max states and variables for DCGRID model (the same model must be always called with the same NS and NV values although not all of them are used) 
    v_nVSC, v_nDCbus, v_nDClines = np.array(l_nVSC), np.array(l_nDCbus), np.array(l_nDClines)
  
    max_NS = max(v_nDCbus+v_nDClines)
    max_NV = max(v_nDCbus*2)

    ntotDCbus = sum(l_nDCbus)
    ntotDClines = sum(l_nDClines)
    ntotVSC = sum(l_nVSC)
    
    # Conversion of converter parameters into converter rating# Loss coefficients: still in pu, but converted to ACSbase
    # Convert loss coefficients from pu on ACSbase to physical units for Fortran
    ACSbase = float(psspy.sysmva())
    v_MbaseVSC = m_converter[:, FIL_RATE_A]
    m_converter[:, FIL_A] = np.divide(m_converter[:, FIL_A],v_MbaseVSC) * (ACSbase)
    m_converter[:, FIL_B] = m_converter[:, FIL_B]  # unchanged
    m_converter[:, FIL_CRECT] = np.multiply(m_converter[:, FIL_CRECT],v_MbaseVSC) / (ACSbase)
    m_converter[:, FIL_CINV] = np.multiply(m_converter[:, FIL_CINV],v_MbaseVSC) / (ACSbase)
    
    # AC bus, idconv, DC bus type
    v_indexACbus = np.searchsorted(m_dcbus[l_indexDCbusVSC,DC_BUS], m_converter[:, FIL_DCBUS])
    m_vscdynmodelbusparam = np.concatenate((m_converter[v_indexACbus, FIL_ACBUS].reshape(-1,1),
        np.ones((ntotVSC,1))*float(idconv),m_dcbus[l_indexDCbusVSC,BUS_TYPE].reshape(-1,1), m_converter[v_indexACbus, FIL_A].reshape(-1,1),
        m_converter[v_indexACbus, FIL_B].reshape(-1,1), m_converter[v_indexACbus, FIL_CRECT].reshape(-1,1), m_converter[v_indexACbus, FIL_CINV].reshape(-1,1)),axis=1)

    m_artificialACbus = np.array(l_artificialACbus)
    
    # adding dynamic data iteratively
    iMTDC = 0
    nDCbus_previous = 0
    nDClines_previous = 0
    v_sumnDCbus = np.cumsum(l_nDCbus)
    
    ivsc = -1
    for idcbus in range(0,ntotDCbus):
    
        # set 
        if (idcbus >= v_sumnDCbus[iMTDC]):
            nDCbus_previous += l_nDCbus[iMTDC]
            nDClines_previous += l_nDClines[iMTDC]
            iMTDC += 1
            
        nDCbusi = l_nDCbus[iMTDC]
        nDClinesi = l_nDClines[iMTDC]   
        
        idDCgrid = "'" + str(iMTDC + 1) + "'"
        
        if m_dcbus[idcbus, DC_BUS] in m_converter[:, FIL_DCBUS]:

            ivsc += 1
 
            with open(str_dyrfilemtdc,'a') as f_dyr:
                str_dynmodelinstance = fun_setMTDCdyrdata(m_artificialACbus[ivsc, 1] if m_artificialACbus[ivsc, 0] == m_vscdynmodelbusparam[ivsc, 0] else None,
                                                      m_vscdynmodelbusparam[ivsc,1], m_vscdynmodelbusparam[ivsc,2], m_vscdynmodelbusparam[ivsc,3], m_vscdynmodelbusparam[ivsc,4],
                                                      m_vscdynmodelbusparam[ivsc,5], m_vscdynmodelbusparam[ivsc,6], nDClinesi, nDCbusi, (idcbus-nDCbus_previous+1), idDCgrid, max_NS, max_NV,
                                                      ntotDCbus, ntotDClines, nDCbus_previous, nDClines_previous) # Updated to consider the bus number at point 'f'
                f_dyr.write(str_dynmodelinstance)
    f_dyr.close()
    
def fun_setMTDCdyrdata(ACbusi, ACbusidi, dctypei, aloss, bloss, crect, cinv, nDClinesi, nDCbusi, ivsci, idDCgrid, max_NS, max_NV, ntotDCbus, ntotDClines, nDCbus_previous, nDClines_previous):

    """ 
    This function sets the pre-set dynamic data which is added to the existing .dyr file. It essentially creates a string with the dynamic data information (bus number, model names, id, parameters).
    """

    #ierr, Vacbase = psspy.busdat(int(ACbusi), 'BASE')
    #Vacbase = float(Vacbase)

    #aloss = aloss * ACSbase
    #bloss = bloss * np.sqrt(3.0) * Vacbase
    #crect = crect * (Vacbase**2 / ACSbase)
    #cinv  = cinv  * (Vacbase**2 / ACSbase)

    # Convert inputs to strings
    ACbusi = str(int(ACbusi))
    ACbusidi = str(int(ACbusidi))
    dctypei = int(dctypei)
    aloss = str(aloss)
    bloss = str(bloss)
    crect = str(crect)
    cinv = str(cinv)
    ivsci = str(ivsci)
    idDCgrid = str(idDCgrid)
    nDCbusi = str(nDCbusi)
    nDClinesi = str(nDClinesi)
    ntotDCbus = str(int(ntotDCbus))
    ntotDClines = str(int(ntotDClines))
    nDCbus_previous = str(int(nDCbus_previous))
    nDClines_previous = str(int(nDClines_previous))
    
    # max number of states and variables
    NS = str(max_NS)
    NV = str(max_NV)
    
    # standard SPWDRD and SQWDRD options
    dyr_text_ln1 = '1 0 0'
    dyr_text_ln2 = '0 0 0'
 
    # in case of DC_TYPE = 2, add DC grid model
    if dctypei == 2:
        #   Modificar 320 U base, 100 MVA base                                           Data list --> 
        slack_dyr = '\n'+ ACbusi +' \'USRMDL\' '+ ACbusidi +' \'DDCGRD\'  5  0  4  0  '+ NS +'  '+ NV +'   '+ nDCbusi +' '+ nDClinesi +' 2 ' + idDCgrid + '/'
    else:
        slack_dyr = ''
        
    # single string containing dynamic data information    
    str_dynmodelinstance = '\n/ bus '+ ACbusi +' - converter\
        \n'+ ACbusi +' \'USRMDL\' '+ ACbusidi +' \'VSCGFL\'  1  1  6  24  7  29   1 1 1 '+ nDCbusi +' ' + nDClinesi + ' '+ idDCgrid +' 0.005 0.00 0.0 20.0 1.2 0.0 0.0 0.0 150.0 3000.0 1.00 500 -500 200 -200 1.10 0.90 1.31 '+ aloss +' ' + bloss + ' '+ crect +' '+ cinv +' 0.1 0.1/\
        \n'+ ACbusi +' \'USRMDL\' '+ ACbusidi +' \'SPWDRD\' 4 0 3 12 2 15 	'+ dyr_text_ln1 +' 	0.100 10.0 0.1 200.0 10.0 0.00001 0.00001 0.0001 1.0 -1.0 999.45 0.2    /\
        \n'+ ACbusi +' \'USRMDL\' '+ ACbusidi +' \'SQWDRD\' 3 0 3 13 2 14 	'+ dyr_text_ln2 +' 	0.100 10.0 0.1 200.0 10.0 0.00001 0.00001 0.0001 0.4 -0.4 999.45 1.5 9999.0   /\
        '+ slack_dyr # +'\\n'+ ACbusi +' \'USRMDL\' '+ ACbusidi +' \'WDELAY\' 9 0   6 5 8 6 	1 4     5 6 11 10  0.000 0.000 0.000 0.000   0.000 /\n'

    return str_dynmodelinstance
    
def fun_write2txtfile(str_openformat, m_data, str_path, str_filename, str_idDCgridi, str_writeformat):
    
    with open(os.path.join(str_path,str_filename) ,str_openformat) as f_file:
        np.savetxt(f_file, m_data, fmt = str_writeformat, delimiter=',', header = str_idDCgridi, comments='')
    f_file.close()

def fun_setDCnetworkdata(str_path4dynamics, d_acdcoutput, l_MTDCgrids, l_indexDCbusVSC, l_artificialACbus):

    # delete existing .txt files (if any)
    for file in os.listdir(str_path4dynamics):
        if file.endswith(".txt"):
            os.remove(os.path.join(str_path4dynamics, file))

    # get DC grid information (number of lines, etc.)
    DCSbase = d_acdcoutput.get("DCSbase")
    m_dcbus = d_acdcoutput.get("dcbus")
    m_dcbranch = d_acdcoutput.get("dcbranch") 
    m_converter = d_acdcoutput.get("converter")
    
    nDClines = np.shape(m_dcbranch[:,0])[0] # # of branches of the DC grid
    nDCbus = np.shape(m_dcbus[:,0])[0] # # of nodes of the DC grid
    #nVSC = np.shape(m_converter[:,0])[0] # # of VSCs of the DC grid
    l_nDClines = []
    l_nDCbus = []
    #l_nVSC = []
    nMTDCgrids = len(l_MTDCgrids)
    m_artificialACbus = np.array(l_artificialACbus)
    for i in range(0,nMTDCgrids):
        a_dcbus = l_MTDCgrids[i]["dcbus"]
        l_nDClines.append(len(l_MTDCgrids[i]["dcbranch"]))
        l_nDCbus.append(len(a_dcbus))
        #l_nVSC.append(len(l_MTDCgrids[i]["converter"]))
    
    v_Udc_ini = m_dcbus[:,DC_UDC]    
    m_Ydc = fun_getDCconductance(DCSbase, nDCbus, nDClines, m_dcbus, m_dcbranch)
    m_Ac = fun_getDCincidence(nDCbus, nDClines, m_dcbus, m_dcbranch)
    v_Gdc = m_dcbus[:,DC_GDC].reshape(nDCbus)
    v_Cdc = m_dcbus[:,DC_CDC].reshape(nDCbus)
    v_Rsdc = m_dcbranch[:,DC_BR_R].reshape(nDClines)
    v_Ldc = m_dcbranch[:,DC_BR_L].reshape(nDClines)
    m_acdc_buses = np.column_stack((
        m_dcbus[:,DC_BUS].astype(int),[m_artificialACbus[m_converter[:,FIL_DCBUS].astype(int) == int(dc_bus_id), 1][0]
                                       if np.any(m_converter[:,FIL_DCBUS].astype(int) == int(dc_bus_id)) else -1
                                       for dc_bus_id in m_dcbus[:,DC_BUS]])) # Changed by Saeed / For adding DC buses without a converter to create the required bus information file for dynamic simulation
   
    # set data for each DC grid with a given header
    str_openformat = 'a'

    nDCbus_previous = 0
    nDClines_previous = 0
    #nVSC_previous = 0    # No need. I will remove it after validation
    for iMTDC in range(0,nMTDCgrids):
        
        str_idDCgrid = str(iMTDC + 1)
       
        # extract data from one grid
        Ac_new = m_Ac[nDCbus_previous:nDCbus_previous+l_nDCbus[iMTDC],nDClines_previous:nDClines_previous+l_nDClines[iMTDC]]        
        Ydc_new = m_Ydc[nDCbus_previous:nDCbus_previous+l_nDCbus[iMTDC],nDCbus_previous:nDCbus_previous+l_nDCbus[iMTDC]]
        Zdc_new = np.linalg.inv(Ydc_new)
        Udc_ini_new = v_Udc_ini[nDCbus_previous:nDCbus_previous+l_nDCbus[iMTDC]]
        Gdc_new = v_Gdc[nDCbus_previous:nDCbus_previous+l_nDCbus[iMTDC]]
        Cdc_new = v_Cdc[nDCbus_previous:nDCbus_previous+l_nDCbus[iMTDC]]
        Rsdc_new = v_Rsdc[nDClines_previous:nDClines_previous+l_nDClines[iMTDC]]
        Ldc_new = v_Ldc[nDClines_previous:nDClines_previous+l_nDClines[iMTDC]]
        acdcbus_new = m_acdc_buses[nDCbus_previous:nDCbus_previous+l_nDCbus[iMTDC]]
        ####################################   Newly added   ###############################################################
        MTDC_numbering_bus_new = np.full(Udc_ini_new.size, iMTDC +1, dtype=int)
        Buses_base_new = np.column_stack((MTDC_numbering_bus_new, acdcbus_new, Udc_ini_new, Gdc_new, Cdc_new))
        MTDC_numbering_line_new = np.full(Rsdc_new.size, iMTDC +1, dtype=int)
        Lines_base_new = np.column_stack((MTDC_numbering_line_new, Rsdc_new, Ldc_new))
        ####################################   %%%%%%%%%%%   ###############################################################

        # write the extracted data below the header str_idDCgrid
        fun_write2txtfile(str_openformat, Ac_new, str_path4dynamics, 'data_Ac.txt', str_idDCgrid, '%-7.4f')
        fun_write2txtfile(str_openformat, Ydc_new, str_path4dynamics, 'data_Ydc.txt', str_idDCgrid, '%-7.8f')
        fun_write2txtfile(str_openformat, Zdc_new, str_path4dynamics, 'data_Zdc.txt', str_idDCgrid, '%-7.8f')
        ####################################   Newly added   ###############################################################
        fun_write2txtfile(str_openformat, Buses_base_new, str_path4dynamics, 'data_Buses_base.txt', [], '%-7.8f')
        fun_write2txtfile(str_openformat, Lines_base_new, str_path4dynamics, 'data_Lines_base.txt', [], '%-7.8f')
        ####################################   %%%%%%%%%%%   ###############################################################

        nDCbus_previous += l_nDCbus[iMTDC]
        nDClines_previous += l_nDClines[iMTDC]
        #nVSC_previous += l_nVSC[iMTDC] 

    return
