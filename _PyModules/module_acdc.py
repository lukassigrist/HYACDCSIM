"""
This code sets up and solves the sequential AC/DC load flow algorithm for multi-terminal VSC-HVDC systems

PSS/E solves the AC grid, whereas a Python module solves the DC grid and the copupling between both the AC and DC grids. 

This codes enables:

    1. solving the AC/DC load flow
    2. showing the solution of the AC/DC load flow
    3. writing data files for dynamic simulations

The module contains the following functions:

- main_runacdclf: principal function to prepare and run the AC/DC load flow
- fun_getT2Piequivalent: build Pi equivalent of T model of transformer and LC filter
- fun_getDCdata: extract DC data
- fun_setinitialvaluesacdclf: initialize sequential AC/DC load flow 
- fun_setACsetpoints: set AC-side generator set points
- fun_getACresults: extract results from the AC load flow run
- fun_initializeDCfromAC: initialize DC-side variables from AC results
- fun_calldclf: call the DC load flow (DC load flow and DC slack bus iteration)
- fun_getacdcoutput: get AC/DC load flow results

- fun_showacdclfoutput: show results

- fun_setacdcdynamicdata: create network and element dynamic data
- fun_updatedynamicfile: update .dyr file with VSC models and DC grid model
- fun_setDCnetworkdata: set DC network data for dynamic simulation
- fun_setMTDCdyrdata: set VSC model and DC grid model parameters
- fun_write2txtfile: write DC network data to .txt file

Authors:
- Aurelio Garcia Cerrada
- Javier Renedo
- Carlos Prieto
- Lukas Sigrist
"""

import os, sys, shutil

import psse34 # PSSE
# PSS/e
import psspy, redirect
redirect.psse2py()
from psspy import _i, _f, _s, _o
# import numpy and math
import numpy as np
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

# Constants
j = 1j

def main_runacdclf(d_acdcoptions, str_savfileorig, str_pathdclfresults, l_MTDCgrids):
    """
    Sequential AC/DC Power Flow algorithm for multi-terminal VSC-HVDC systems   
    """
    # ===================================================================
    # Preparing sequential ACDC load flow
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
    str_savfilemtdc = str_savfileorig[:-4] + r"_MTDCg.sav" # update AC load flow (VSC as generators)
    str_rawfilemtdc = str_savfilemtdc[:-4] + r".raw"
    
    # open PSS/e file
    psspy.case(str_savfileorig)

    # system VA base
    ACSbase = psspy.sysmva()    

    # get static invariant DC data
    [DCSbase, nDCbus, v_indexDCslack, m_converter, m_dcbus, m_dcbranch] = fun_getDCdata(l_MTDCgrids)

    # get converter losses parameters
    v_Sn_converter_MVA = m_converter[:, FIL_RATE_A]
    v_aloss = m_converter[:, FIL_A]#/v_Sn_converter_MVA*DCSbase
    v_bloss = m_converter[:, FIL_B]#/v_Sn_converter_MVA*DCSbase
    v_crect = m_converter[:, FIL_CRECT]#/v_Sn_converter_MVA*DCSbase
    v_cinv = m_converter[:, FIL_CINV]#/v_Sn_converter_MVA*DCSbase
    
    # Transformer and LC filter
    [v_Am, v_Bm, v_Cm, v_Dm, v_Y1, v_Z2, v_Y3] = fun_getT2Piequivalent(m_converter)
    
    # ===================================================================
    # Sequential AC/DC load flow iterations
    # =================================================================== 

    # initialize powers and voltages of sequential acdc load flow
    [v_Us, v_Theta_s, v_Ps_pu, v_Qs_pu, v_Pdc_pu, v_PlossVSC] = fun_setinitialvaluesacdclf(DCSbase, v_indexDCslack, nDCbus, m_dcbus)
    v_Us0 = v_Us
    v_Qs0 = v_Qs_pu  
     
    isconverged = 0 # convergence loop
    it = 0
    v_k_int_ac = []
    while not isconverged: # BEGIN WHILE
        
        it = it + 1
        v_Ps_ns_previous = v_Ps_pu[v_indexDCslack,0].astype(float)

        # ===================================================================
        # AC grid power flow
        # ===================================================================
        
        # 1. Assign Ps, Qs and Us set points to AC-side converters     
        fun_setACsetpoints(v_Us, v_Theta_s, v_Ps_pu, v_Us0, v_Qs0, DCSbase, m_converter, m_dcbus, idconv, idowner)
             
        # 2. Solve AC load flow
        psspy.fnsl(
            options1=0, # disable tap stepping adjustment.
            options5=0, # disable switched shunt adjustment.
            # options6=0, # flat start (first time).
        )
        v_k_int_ac.append(psspy.iterat())

        # 3. Extract solutions of AC load flow 
        [v_Us, v_Theta_s, v_Ps_pu, v_Qs_pu] = fun_getACresults(v_Us, v_Theta_s, v_Ps_pu, v_Qs_pu, 
            DCSbase, m_converter, m_dcbus, idconv)      
                    
        # 4. Compute the DC-side active power for the subsequent DC-side power flow
        [ v_Uc_phasor, v_Is_phasor, v_Ic_phasor, v_Pc_pu, v_Qc_pu, v_Pdc_pu, v_PlossVSC] = fun_initializeDCfromAC(v_Us, 
            v_Theta_s, v_Ps_pu, v_Qs_pu, v_Pdc_pu, v_PlossVSC, v_Am, v_Bm, v_Cm, v_Dm, v_aloss, v_bloss, v_crect, v_cinv)

        m_dcbus[:,DC_PDC] = v_Pdc_pu*DCSbase        
            
        # ===================================================================
        # DC grid power flow
        # ===================================================================

        [v_Ps_ns_current, k_int_dc, v_mismatch_dc, k_int_dcslack, j_int_dcslack, v_Ps_pu, v_Pc_pu, v_Qc_pu, v_Pdc_pu, v_Uc_phasor, 
            m_dcbus, m_dcbranch] = fun_calldclf(v_Us, v_Theta_s, v_Uc_phasor, v_Ps_pu, v_Qs_pu, v_Pc_pu, v_Qc_pu, v_Pdc_pu, DCSbase, 
            m_dcbus, m_dcbranch, v_Y1, v_Z2, v_Y3, v_aloss, v_bloss, v_crect, v_cinv, l_MTDCgrids, v_indexDCslack, lftol, lfmaxiter)

        # ===================================================================
        # Convergence 
        # ===================================================================
        v_Ps_ns_current = np.transpose(np.matrix(v_Ps_ns_current))
        error = v_Ps_ns_current - v_Ps_ns_previous
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
    d_acdcoutput = fun_getacdcoutput(v_Is_phasor, v_Uc_phasor, v_Ic_phasor, v_Pc_pu, v_Qc_pu, v_PlossVSC, 
        DCSbase, m_dcbus, m_dcbranch, m_converter)

    psspy.fnsl(
        options1=0, # disable tap stepping adjustment.
        options5=0, # disable switched shunt adjustment.
    ) # we have finished, but we run a power flow again to avoid strange things with the tricky thing of PQ generators

    # write AC-side to file
    psspy.rawd_2(0,1,[1,1,1,0,0,0,0],0,str_rawfilemtdc) # save .raw
    psspy.save(str_savfilemtdc) # save .sav   
         
    return d_acdcoutput, isconverged, it, k_int_ac, k_int_dc, k_int_dcslack, j_int_dcslack, v_mismatch_dc

def fun_getT2Piequivalent(m_converter):
    # Transformer filter (p.u)
    v_Rt = m_converter[:, FIL_RT]
    v_Xt = m_converter[:, FIL_XT]
    v_Zt = v_Rt + j*v_Xt

    # Low pass filter (p.u)
    v_Bf = m_converter[:, FIL_BF]
    # Zf = -j/v_Bf # division term by term
    v_Yf = j*v_Bf

    # Phase reactor (p.u)
    v_Rc = m_converter[:, FIL_RC]
    v_Xc = m_converter[:, FIL_XC]
    v_Zc = v_Rc + j*v_Xc

    # T model of filter and transformer to equivalent Pi model (Steinmetz)
    v_Zctf = np.multiply( v_Zc, np.multiply(v_Zt, v_Yf) ) 
    v_Y1 = np.multiply( v_Yf, v_Zc/(v_Zc + v_Zt + v_Zctf) )
    v_Z2 = v_Zc + v_Zt + v_Zctf
    v_Y3 = np.multiply( v_Yf, v_Zt/(v_Zc + v_Zt + v_Zctf) ) 
    v_Am = 1.0 + np.multiply(v_Z2, v_Y1)
    v_Bm = v_Z2
    v_Cm = v_Y1 + v_Y3 + np.multiply( v_Z2, np.multiply(v_Y1, v_Y3) )
    v_Dm = 1.0 + np.multiply(v_Z2, v_Y3)
    
    return v_Am, v_Bm, v_Cm, v_Dm, v_Y1, v_Z2, v_Y3

def fun_getDCdata(l_MTDCgrids):
    DCSbase = l_MTDCgrids[0]["baseMVA"] # DC base power is the one of the first MTDC grid 
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
    
    '''
    # common bases - to be completed       
    m_converter[:, FIL_A] = np.divide(m_converter[:, FIL_A],v_Sn_converter_MVA)*DCSbase
    m_converter[:, FIL_B] = np.divide(m_converter[:, FIL_B],v_Sn_converter_MVA)*DCSbase
    m_converter[:, FIL_CRECT] = np.divide(m_converter[:, FIL_CRECT],v_Sn_converter_MVA)*DCSbase
    m_converter[:, FIL_CINV] = np.divide(m_converter[:, FIL_CINV],v_Sn_converter_MVA)*DCSbase
    m_converter[:, FIL_RT] = np.multiply(m_converter[:, FIL_RT],v_Sn_converter_MVA)/DCSbase
    m_converter[:, FIL_XT] = np.multiply(m_converter[:, FIL_XT],v_Sn_converter_MVA)/DCSbase
    m_converter[:, FIL_BF] = np.divide(m_converter[:, FIL_BF],v_Sn_converter_MVA)*DCSbase
    m_converter[:, FIL_RC] = np.multiply(m_converter[:, FIL_RC],v_Sn_converter_MVA)/DCSbase
    m_converter[:, FIL_XC] = np.multiply(m_converter[:, FIL_XC],v_Sn_converter_MVA)/DCSbase
    '''
    
    # DC slack indexes
    v_indexDCslack = [i for i, x in enumerate(m_dcbus[:,DC_TYPE]) if x == 2]
    
    return DCSbase, nDCbus, v_indexDCslack, m_converter, m_dcbus, m_dcbranch

def fun_setinitialvaluesacdclf(DCSbase, v_indexDCslack, nDCbus, m_dcbus):
    v_Ps_MW = m_dcbus[:, DC_PS]
    v_Qs_Mvar = m_dcbus[:, DC_QS]
    v_Ps_pu = v_Ps_MW/DCSbase # initial values
    v_Qs_pu = v_Qs_Mvar/DCSbase # initial values
    v_Ss_pu = v_Ps_pu + j*v_Qs_pu

    v_Pdc_pu = -v_Ps_pu # only to create the vector

    v_Us = m_dcbus[:,DC_US] # initial value
    v_Theta_s = m_dcbus[:,DC_DELTAS]
    v_Us_phasor = np.multiply( v_Us, np.exp(j*v_Theta_s) )
        
    # Pre-alocate memory
    v_PlossVSC = np.matrix(np.ones(shape=(nDCbus, 1)))
    
    return v_Us, v_Theta_s, v_Ps_pu, v_Qs_pu, v_Pdc_pu, v_PlossVSC

def fun_setACsetpoints(v_Us, v_Theta_s, v_Ps_pu, v_Us0, v_Qs0, DCSbase, m_converter, m_dcbus, idconv, idowner):
    
    # get DC-side indices according to ac-side gen bus type (qs: type2=1, us: type2=2, slack: type2=3)
    l_indexDCtype1 = list(np.array((np.nonzero(m_dcbus[:,DC_TYPE2]==1)[0]).astype(int)).reshape(-1,)) # make indices a integer-type list
    l_indexDCtype2 = list(np.array((np.nonzero(m_dcbus[:,DC_TYPE2]==2)[0]).astype(int)).reshape(-1,)) 
    l_indexDCtype3 = list(np.array((np.nonzero(m_dcbus[:,DC_TYPE2]==3)[0]).astype(int)).reshape(-1,))
    
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
    v_Sn_converter_MVA = m_converter[:, FIL_RATE_A]

    # PQ Limits of the converters (p.u)
    v_Pcmax = np.matrix(1.0*m_converter[:, FIL_PMAX])
    v_Pcmin = np.matrix(1.0*m_converter[:, FIL_PMIN])
    v_Qcmax = np.matrix(1.0*m_converter[:, FIL_QMAX])
    v_Qcmin = np.matrix(1.0*m_converter[:, FIL_QMIN])
    
    # Converters that control Qs          
    for i in range(0, nVSCpq):
        # Ps and Qs injection after the transformer (Qmin = Qmax = Qs)
        Rsorce_pu_conv = (v_Rt[l_indexDCtype1[i],0]+v_Rc[l_indexDCtype1[i],0])*v_Sn_converter_MVA[l_indexDCtype1[i],0]/DCSbase # p.u-machine
        Xsorce_pu_conv = (v_Xt[l_indexDCtype1[i],0]+v_Xc[l_indexDCtype1[i],0])*v_Sn_converter_MVA[l_indexDCtype1[i],0]/DCSbase # p.u-machine
        
        # Zsorce is only used in the dynamics.. note that we are neglecting the shunt filter. we should decide what to do with that
        # print(Rsorce_pu_conv,Xsorce_pu_conv) -> do a backward computation by representing the transformer and shunt explicitely and the recompute voltage and power set points...
        
        psspy.plant_data(l_busACtype1[i],_i,[ v_Us[l_indexDCtype1[i],0],_f]) # Control Us -- it is not necessary           
        psspy.machine_data_2(l_busACtype1[i], idconv,[_i,idowner,_i,_i,_i,_i],[ v_Ps_pu[l_indexDCtype1[i],0]*DCSbase, v_Qs0[l_indexDCtype1[i],0]*DCSbase, v_Qcmax[l_indexDCtype1[i],0],v_Qcmin[l_indexDCtype1[i],0], v_Pcmax[l_indexDCtype1[i],0], v_Pcmin[l_indexDCtype1[i],0],_f,_f, 0.0001,_f,_f,0.0,_f,0.0,0.0,0.0,_f])
        psspy.machine_data_2(l_busACtype1[i], idconv,[_i,_i,_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f, v_Sn_converter_MVA[l_indexDCtype1[i],0],_f,_f,_f,_f,_f,_f,_f,_f,_f,_f]) # de momento pongo potencia nominal del conv = a 100MVA
        psspy.machine_data_2(l_busACtype1[i], idconv,[_i,_i,_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f,_f, Rsorce_pu_conv,_f,_f,_f,_f,_f,_f,_f,_f,_f]) # de momento sin filtro paralelo
        psspy.machine_data_2(l_busACtype1[i], idconv,[_i,_i,_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f,_f,_f, Xsorce_pu_conv,_f,_f,_f,_f,_f,_f,_f,_f])
        psspy.machine_data_2(l_busACtype1[i], idconv,[_i,_i,_i,_i,_i,_i],[_f, _f, v_Qs0[l_indexDCtype1[i],0]*DCSbase, v_Qs0[l_indexDCtype1[i],0]*DCSbase,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f])
        psspy.bus_data_2(l_busACtype1[i],[2,_i,_i,_i],[_f,_f,_f],_s)

    # Converters that control Us
    for i in range(0, nVSCpv): 
        # include a new generator with Ps and Us
        Rsorce_pu_conv = (v_Rt[l_indexDCtype2[i],0]+v_Rc[l_indexDCtype2[i],0])*v_Sn_converter_MVA[l_indexDCtype2[i],0]/DCSbase # p.u-machine
        Xsorce_pu_conv = (v_Xt[l_indexDCtype2[i],0]+v_Xc[l_indexDCtype2[i],0])*v_Sn_converter_MVA[l_indexDCtype2[i],0]/DCSbase # p.u-machine
        # Zsorce is only used in the dynamics.. note that we are neglecting the shunt filter. we should decide what to do with that
        psspy.plant_data(l_busACtype2[i],_i,[ v_Us0[l_indexDCtype2[i],0],_f]) # Control Us
        psspy.machine_data_2(l_busACtype2[i], idconv,[_i,idowner,_i,_i,_i,_i],[ v_Ps_pu[l_indexDCtype2[i],0]*DCSbase, _f, v_Qcmax[l_indexDCtype2[i],0], v_Qcmin[l_indexDCtype2[i],0], v_Pcmax[l_indexDCtype2[i],0], v_Pcmin[l_indexDCtype2[i],0],_f,_f, 0.0001,_f,_f,0.0,_f,0.0,0.0,0.0,_f])
        psspy.machine_data_2(l_busACtype2[i], idconv,[_i,_i,_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f, v_Sn_converter_MVA[l_indexDCtype2[i],0],_f,_f,_f,_f,_f,_f,_f,_f,_f,_f]) # de momento pongo potencia nominal del conv = a 100MVA
        psspy.machine_data_2(l_busACtype2[i], idconv,[_i,_i,_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f,_f, Rsorce_pu_conv,_f,_f,_f,_f,_f,_f,_f,_f,_f]) # de momento sin filtro paralelo
        psspy.machine_data_2(l_busACtype2[i], idconv,[_i,_i,_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f,_f,_f, Xsorce_pu_conv,_f,_f,_f,_f,_f,_f,_f,_f])
        psspy.bus_data_2(l_busACtype2[i],[2,_i,_i,_i],[_f,_f,_f],_s) # say to psse that the bus is a PV node (code = 2)

    # Converters that control Us (only if the converter is feeding a passive grid)
    for i in range(0, nVSCslack): 
        # include a new generator with Ps and Us
        Rsorce_pu_conv = (v_Rt[l_indexDCtype3[i],0]+v_Rc[l_indexDCtype3[i],0])*v_Sn_converter_MVA[l_indexDCtype3[i],0]/DCSbase # p.u-machine
        Xsorce_pu_conv = (v_Xt[l_indexDCtype3[i],0]+v_Xc[l_indexDCtype2[i],0])*v_Sn_converter_MVA[l_indexDCtype3[i],0]/DCSbase # p.u-machine
        # Zsorce is only used in the dynamics.. note that we are neglecting the shunt filter. we should decide what to do with that
        psspy.plant_data(l_busACtype3[i],_i,[ v_Us0[l_indexDCtype3[i],0],_f]) # Control Us
        psspy.machine_data_2(l_busACtype3[i], idconv,[_i,idowner,_i,_i,_i,_i],[ v_Ps_pu[l_indexDCtype3[i],0]*DCSbase, _f, v_Qcmax[l_indexDCtype3[i],0], v_Qcmin[l_indexDCtype3[i],0], v_Pcmax[l_indexDCtype3[i],0], v_Pcmin[l_indexDCtype3[i],0],_f,_f, 0.0001,_f,_f,0.0,_f,0.0,0.0,0.0,_f])
        psspy.machine_data_2(l_busACtype3[i], idconv,[_i,_i,_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f, v_Sn_converter_MVA[l_indexDCtype3[i],0],_f,_f,_f,_f,_f,_f,_f,_f,_f,_f]) # de momento pongo potencia nominal del conv = a 100MVA
        psspy.machine_data_2(l_busACtype3[i], idconv,[_i,_i,_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f,_f, Rsorce_pu_conv,_f,_f,_f,_f,_f,_f,_f,_f,_f]) # de momento sin filtro paralelo
        psspy.machine_data_2(l_busACtype3[i], idconv,[_i,_i,_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f,_f,_f, Xsorce_pu_conv,_f,_f,_f,_f,_f,_f,_f,_f])
        psspy.bus_data_2(l_busACtype3[i],[3,_i,_i,_i],[_f,_f,0.0],_s) # say to psse that the bus is a ac-slack node (code = 2)

def fun_getACresults(v_Us, v_Theta_s, v_Ps_pu, v_Qs_pu, DCSbase, m_converter, m_dcbus, idconv):
       
    # get DC-side indices according to ac-side gen bus type (qs: type2=1, us: type2=2, slack: type2=3)
    l_indexDCtype1 = list(np.array((np.nonzero(m_dcbus[:,DC_TYPE2]==1)[0]).astype(int)).reshape(-1,)) # make indices a integer-type list
    l_indexDCtype2 = list(np.array((np.nonzero(m_dcbus[:,DC_TYPE2]==2)[0]).astype(int)).reshape(-1,)) 
    l_indexDCtype3 = list(np.array((np.nonzero(m_dcbus[:,DC_TYPE2]==3)[0]).astype(int)).reshape(-1,))
    
    # get corresponding AC-side bus 
    l_busACtype1 = list(np.array((m_converter[l_indexDCtype1, FIL_ACBUS]).astype(int)).reshape(-1,))
    l_busACtype2 = list(np.array((m_converter[l_indexDCtype2, FIL_ACBUS]).astype(int)).reshape(-1,))
    l_busACtype3 = list(np.array((m_converter[l_indexDCtype3, FIL_ACBUS]).astype(int)).reshape(-1,))
       
    nVSCslack = len(l_busACtype3) # number of slack-type AC-side converters
    nVSCpv = len(l_busACtype2) # number of PV-type AC-side converters
    nVSCpq = len(l_busACtype1) # number of PQ-type AC-side converters
            
    # Converters that control Qs
    for i in range(0, nVSCpq):
        ierr, v_Us[l_indexDCtype1[i],0] = psspy.busdat(l_busACtype1[i], 'PU')
        ierr, v_Theta_s[l_indexDCtype1[i],0] = psspy.busdat(l_busACtype1[i], 'ANGLE')

        ierr, rval = psspy.macdat(l_busACtype1[i], idconv, 'P') # Pgen [MW] (API)
        v_Ps_pu[l_indexDCtype1[i],0] = rval/DCSbase
        ierr, rval = psspy.macdat(l_busACtype1[i], idconv, 'Q') # Pgen [MW] (API)
        v_Qs_pu[l_indexDCtype1[i],0] = rval/DCSbase
        
    # Vac control
    for i in range(0, nVSCpv):
        ierr, v_Us[l_indexDCtype2[i],0] = psspy.busdat(l_busACtype2[i], 'PU')
        ierr, v_Theta_s[l_indexDCtype2[i],0] = psspy.busdat(l_busACtype2[i], 'ANGLE')
        
        ierr, rval = psspy.macdat(l_busACtype2[i], idconv, 'P') # Pgen [MW] (API)
        v_Ps_pu[l_indexDCtype2[i],0] = rval/DCSbase
        ierr, rval = psspy.macdat(l_busACtype2[i], idconv, 'Q') # Pgen [MW] (API)
        v_Qs_pu[l_indexDCtype2[i],0] = rval/DCSbase

    # ac-slack converter
    for i in range(0, nVSCslack):
        ierr, v_Us[l_indexDCtype3[i],0] = psspy.busdat(l_busACtype3[i], 'PU')
        ierr, v_Theta_s[l_indexDCtype3[i],0] = psspy.busdat(l_busACtype3[i], 'ANGLE')
        
        ierr, rval = psspy.macdat(l_busACtype3[i], idconv, 'P') # Pgen [MW] (API)
        v_Ps_pu[l_indexDCtype3[i],0] = rval/DCSbase
        ierr, rval = psspy.macdat(l_busACtype3[i], idconv, 'Q') # Pgen [MW] (API)
        v_Qs_pu[l_indexDCtype3[i],0] = rval/DCSbase 

    return v_Us, v_Theta_s, v_Ps_pu, v_Qs_pu
    
def fun_initializeDCfromAC(v_Us, v_Theta_s, v_Ps_pu, v_Qs_pu, v_Pdc_pu, v_PlossVSC, v_Am, v_Bm, v_Cm, v_Dm,
    v_aloss, v_bloss, v_crect, v_cinv):
       
    v_Us_phasor = np.multiply(v_Us, np.exp(j*v_Theta_s))
    v_Ss_pu = v_Ps_pu + j*v_Qs_pu
    v_Is_phasor = np.conj(v_Ss_pu/v_Us_phasor)
    v_Uc_phasor = np.multiply(v_Am, v_Us_phasor) + np.multiply(v_Bm, v_Is_phasor) # v_Uc_phasor = v_Uf_phasor + np.multiply(v_Zt,v_If_phasor)
    v_Ic_phasor = np.multiply(v_Cm, v_Us_phasor) + np.multiply(v_Dm, v_Is_phasor) # v_Ic_phasor = v_Is_phasor + np.multiply(v_Yf,v_Uf_phasor)
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
    
    return v_Uc_phasor, v_Is_phasor, v_Ic_phasor, v_Pc_pu, v_Qc_pu, v_Pdc_pu, v_PlossVSC  

def fun_calldclf(v_Us, v_Theta_s, v_Uc_phasor, v_Ps_pu, v_Qs_pu, v_Pc_pu, v_Qc_pu, v_Pdc_pu, DCSbase, 
    m_dcbus, m_dcbranch, v_Y1, v_Z2, v_Y3, v_aloss, v_bloss, v_crect, v_cinv, l_MTDCgrids, v_indexDCslack, 
    lftol, lfmaxiter):
    
    # initial index of each MTDC grid data entries
    index_ini_bus = 0
    index_ini_conv = 0
    index_ini_branch = 0

    # iteration and mismatch parameters
    k_int_dc = []
    v_mismatch_dc = []
    k_int_dcslack = []
    j_int_dcslack = []
    v_Ps_ns_current = []
    
    nMTDCgrids = len(l_MTDCgrids)
    for i in range (0, nMTDCgrids):
        
        nDCbusi = len(l_MTDCgrids[i]["dcbus"])
        nDCbranchi = len(l_MTDCgrids[i]["dcbranch"])
        nVSCi = len(l_MTDCgrids[i]["converter"])
    
        # 2. Run DC-side load flow
        m_dcbusi = m_dcbus[index_ini_bus:nDCbusi+index_ini_bus]
        m_dcbranchi = m_dcbranch[index_ini_branch:nDCbranchi+index_ini_branch]
               
        [m_dcbusi, m_dcbranchi, k_int_dc_dcpf, mismatch_dci] = fun_mainDCloadflow(DCSbase, lftol, lfmaxiter, 
            m_dcbusi, m_dcbranchi);                      
        Pdc_dcpf = m_dcbusi[:,DC_PDC]/DCSbase
        
        # 3. Run slack-bus iteration
        Us_dcpf = v_Us[index_ini_bus:nDCbusi+index_ini_bus]
        Th_s_dcpf = v_Theta_s[index_ini_bus:nDCbusi+index_ini_bus]
        Ps_dcpf = v_Ps_pu[index_ini_bus:nDCbusi+index_ini_bus]
        Qs_dcpf = v_Qs_pu[index_ini_bus:nDCbusi+index_ini_bus]
        Pc_dcpf = v_Pc_pu[index_ini_bus:nDCbusi+index_ini_bus]
        Qc_dcpf = v_Qc_pu[index_ini_bus:nDCbusi+index_ini_bus]
        Uc_phasor_dcpf = v_Uc_phasor[index_ini_bus:nDCbusi+index_ini_bus]
        
        # get DC slack
        n_dcslack_dcpf = v_indexDCslack[i]-index_ini_bus          
        ps_n = float(Ps_dcpf[n_dcslack_dcpf, 0])
        pc_n = float(Pc_dcpf[n_dcslack_dcpf, 0])
        qc_n = float(Qc_dcpf[n_dcslack_dcpf, 0])
        us_n = Us_dcpf[n_dcslack_dcpf, 0]
        th_n = Th_s_dcpf[n_dcslack_dcpf]
        uc_n = abs(Uc_phasor_dcpf[n_dcslack_dcpf, 0])
        thc_n = np.angle(Uc_phasor_dcpf[n_dcslack_dcpf, 0])
        qs_n = Qs_dcpf[n_dcslack_dcpf, 0]
        pdc_n = Pdc_dcpf[n_dcslack_dcpf, 0]
        ps_n = Ps_dcpf[n_dcslack_dcpf, 0]
        
        y1_n = v_Y1[index_ini_conv:nVSCi+index_ini_conv][n_dcslack_dcpf, 0]
        z2_n = v_Z2[index_ini_conv:nVSCi+index_ini_conv][n_dcslack_dcpf, 0]
        y3_n = v_Y3[index_ini_conv:nVSCi+index_ini_conv][n_dcslack_dcpf, 0]
        a_n = v_aloss[index_ini_conv:nVSCi+index_ini_conv][n_dcslack_dcpf, 0]
        b_n = v_bloss[index_ini_conv:nVSCi+index_ini_conv][n_dcslack_dcpf, 0]
        c_inv_n = v_cinv[index_ini_conv:nVSCi+index_ini_conv][n_dcslack_dcpf, 0]
        c_rect_n = v_crect[index_ini_conv:nVSCi+index_ini_conv][n_dcslack_dcpf, 0]

        # DC slack iteration
        [Ps_dcpf[n_dcslack_dcpf], Pc_dcpf[n_dcslack_dcpf], Qc_dcpf[n_dcslack_dcpf], 
        Uc_phasor_dcpf[n_dcslack_dcpf], k_int_dcslack_dcpf, j_int_dcslack_dcpf, convergence] = fun_slackiterationDC(us_n, 
        th_n, qs_n, pdc_n, ps_n, y1_n, z2_n, y3_n, a_n, b_n, c_inv_n, c_rect_n, lftol, lfmaxiter)
        
        # update of dcbus matrix with the AC-side information of the converter
        m_dcbusi[:,DC_PS] = Ps_dcpf*DCSbase # put the new Ps in the dcbus matrix 
        m_dcbusi[:,DC_US] = Us_dcpf
        m_dcbusi[:,DC_PDC] = Pdc_dcpf
        m_dcbusi[:,DC_DELTAS] = Th_s_dcpf*180/pi
        m_dcbusi[:,DC_QS] = Qs_dcpf*DCSbase 
        

        # 4. Update variables for outer AC power flow        
        m_dcbus[index_ini_bus:nDCbusi+index_ini_bus] = m_dcbusi
        v_Pdc_pu[index_ini_bus:nDCbusi+index_ini_bus] = Pdc_dcpf
        v_Ps_pu[index_ini_bus:nDCbusi+index_ini_bus] = Ps_dcpf
        v_Pc_pu[index_ini_bus:nDCbusi+index_ini_bus] = Pc_dcpf
        v_Qc_pu[index_ini_bus:nDCbusi+index_ini_bus] = Qc_dcpf
        v_Uc_phasor[index_ini_bus:nDCbusi+index_ini_bus] = Uc_phasor_dcpf
        m_dcbranch[index_ini_branch:nDCbranchi+index_ini_branch] = m_dcbranchi

        k_int_dc.append(k_int_dc_dcpf)
        v_mismatch_dc.append(mismatch_dci)
        k_int_dcslack.append(k_int_dcslack_dcpf)
        j_int_dcslack.append(j_int_dcslack_dcpf)
        v_Ps_ns_current.append(float(Ps_dcpf[n_dcslack_dcpf, 0]))
        
        # next MTDC grid indexes            
        index_ini_bus = len(l_MTDCgrids[i]["dcbus"])+index_ini_bus
        index_ini_conv = len(l_MTDCgrids[i]["converter"])+index_ini_conv
        index_ini_branch = len(l_MTDCgrids[i]["dcbranch"])+index_ini_branch
        
    return v_Ps_ns_current, k_int_dc, v_mismatch_dc, k_int_dcslack, j_int_dcslack, v_Ps_pu, v_Pc_pu, v_Qc_pu, v_Pdc_pu, v_Uc_phasor, m_dcbus, m_dcbranch
    
def fun_getacdcoutput(v_Is_phasor, v_Uc_phasor, v_Ic_phasor, v_Pc_pu, v_Qc_pu, v_PlossVSC, DCSbase, m_dcbus, m_dcbranch, m_converter):
    # get losses
    v_Rt = m_converter[:, FIL_RT]
    v_Rc = m_converter[:, FIL_RC]
    v_Ploss_tf = np.multiply(v_Rt, abs(np.power(v_Is_phasor, 2)))
    v_Ploss_filter = np.multiply(v_Rc, abs(np.power(v_Ic_phasor, 2)));

    # get DC variables
    m_converter = np.array(m_converter)
    m_dcbus = np.array(m_dcbus)
    m_dcbranch = np.array(m_dcbranch)
    
    # get AC variables
    sid = -1 # all buses
    flag = 1
    [ierr, v_Um_bus] = psspy.abusreal(sid, flag, 'PU') # the result is a list (v_Um_bus and v_Theta_bus)
    [ierr, v_Theta_bus] = psspy.abusreal(sid, flag, 'ANGLE')
    v_Um_bus = np.transpose(np.matrix(v_Um_bus)) # column vectors (type matrix)
    v_Theta_bus = np.transpose(np.matrix(v_Theta_bus))
    [ierr, rarray] = psspy.aflowreal(-1, 1, 1, 1, 'P')
    Pf_ac_MW = np.matrix(rarray)
    Pf_ac_MW = np.transpose(np.matrix(Pf_ac_MW))
    [ierr, rarray] = psspy.aflowreal(-1, 1, 1, 1, 'Q')
    Qf_ac_MW = np.matrix(rarray)
    Qf_ac_MW = np.transpose(np.matrix(Qf_ac_MW))
    
    d_acdcoutput = {"version": '2'}
    d_acdcoutput["DCSbase"] = DCSbase
    d_acdcoutput["acbus"] = np.concatenate((v_Um_bus, v_Theta_bus*180/pi),1)
    d_acdcoutput["acbranch"] = np.concatenate((Pf_ac_MW, Qf_ac_MW),1)
    d_acdcoutput["vsc"] = np.concatenate((abs(v_Uc_phasor), np.angle(v_Uc_phasor), v_Pc_pu, v_Qc_pu),1)
    d_acdcoutput["converter"] = m_converter
    d_acdcoutput["dcbus"] = m_dcbus
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

def main_setacdcdynamicdata(d_acdcoutput, str_dyrfileorig, str_path4dynamics, d_acdcoptions, l_MTDCgrids):


    if d_acdcoptions["issetdynamicfiles"]==True:
    
        idconv = d_acdcoptions.get('idconv') 
        
        # write to text files DC grid data
        fun_setDCnetworkdata(str_path4dynamics, d_acdcoutput, l_MTDCgrids)
       
        # update .dyr file with VSC and DC grid models
        fun_updatedynamicfile(str_dyrfileorig, idconv, d_acdcoutput, l_MTDCgrids)
   
def fun_updatedynamicfile(str_dyrfileorig, idconv, d_acdcoutput, l_MTDCgrids):
    """
    This function updates the existing .dyr file by adding MTDC network and converter and control models with pre-set parameters
    """
    # modified file name and copy original file into it
    str_dyrfilemtdc = str_dyrfileorig[:-4] + r"_MTDCg.dyr" # update AC load flow (VSC as generators)
    shutil.copy2(str_dyrfileorig, str_dyrfilemtdc)
    
    DCSbase = d_acdcoutput.get("DCSbase")
    m_dcbus = d_acdcoutput.get("dcbus")
    m_dcbranch = d_acdcoutput.get("dcbranch") 
    m_converter = d_acdcoutput["converter"]
           
    nDClines = np.shape(m_dcbranch[:,0])[0] # # of branches of the DC grid
    nDCbus = np.shape(m_dcbus[:,0])[0] # # of nodes of the DC grid
    l_nDClines = []
    l_nDCbus = []        
    for i in range(0,len(l_MTDCgrids)):
        a_dcbus = l_MTDCgrids[i]["dcbus"]
        l_nDClines.append(len(l_MTDCgrids[i]["dcbranch"]))
        l_nDCbus.append(len(a_dcbus))
        
    # max states and variables for DCGRID model (the same model must be always called with the same NS and NV values although not all of them are used) 
    v_nDCbus, v_nDClines = (np.array(l_nDCbus), np.array(l_nDClines))
    max_NS = max(v_nDCbus+v_nDClines)
    max_NV = max(v_nDCbus*2)

    ntotDCbus = sum(l_nDCbus)
    ntotDClines = sum(l_nDClines)

    # AC bus, idconv, DC bus type
    v_indexACbus = np.searchsorted(m_dcbus[:,DC_BUS], m_converter[:, FIL_DCBUS])
    m_vscdynmodelbusparam = np.concatenate((m_converter[v_indexACbus, FIL_ACBUS].reshape(-1,1),
        np.ones((ntotDCbus,1))*float(idconv),m_dcbus[:,BUS_TYPE].reshape(-1,1)),axis=1)
   
    # adding dynamic data iteratively
    iMTDC = 0
    nDCbus_previous = 0
    nDClines_previous = 0
    v_sumnDCbus = np.cumsum(l_nDCbus)
    for ivsc in range(0,ntotDCbus):
    
        # set 
        if (ivsc >= v_sumnDCbus[iMTDC]):
            nDCbus_previous += l_nDCbus[iMTDC]
            nDClines_previous += l_nDClines[iMTDC]
            iMTDC += 1
            
        nDCbusi = l_nDCbus[iMTDC]
        nDClinesi = l_nDClines[iMTDC]   
        
        idDCgrid = "'" + chr(ord('A') + iMTDC % 26) + "'"  
        
        with open(str_dyrfilemtdc,'a') as f_dyr:
            str_dynmodelinstance = fun_setMTDCdyrdata(m_vscdynmodelbusparam[ivsc,0], m_vscdynmodelbusparam[ivsc,1], m_vscdynmodelbusparam[ivsc,2], nDClinesi, nDCbusi, (ivsc-nDCbus_previous+1), idDCgrid, max_NS, max_NV, ntotDCbus, ntotDClines, nDCbus_previous, nDClines_previous)
            f_dyr.write(str_dynmodelinstance)
    f_dyr.close()
    
def fun_setMTDCdyrdata(ACbusi, ACbusidi, dctypei, nDClinesi, nDCbusi, ivsci, idDCgrid, max_NS, max_NV, ntotDCbus, ntotDClines, nDCbus_previous, nDClines_previous):

    """ 
    This function sets the pre-set dynamic data which is added to the existing .dyr file. It essentially creates a string with the dynamic data information (bus number, model names, id, parameters).
    """  

    # Convert inputs to strings
    ACbusi = str(int(ACbusi))
    ACbusidi = str(int(ACbusidi))
    dctypei = int(dctypei)
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
        slack_dyr = '\n'+ ACbusi +' \'USRMDL\' '+ ACbusidi +' \'DCGRID\'  5  0  8  2  '+ NS +'  '+ NV +'   '+ nDCbusi +' '+ nDClinesi +' 2 ' + idDCgrid + ' '+ ntotDCbus +' '+ ntotDClines +' '+ nDCbus_previous +' '+ nDClines_previous +' 220 100/'
    else:
        slack_dyr = ''
        
    # single string containing dynamic data information    
    str_dynmodelinstance = '\n/ bus '+ ACbusi +' - converter\
        \n'+ ACbusi +' \'USRMDL\' '+ ACbusidi +' \'SVSCON\'  1  1  8  26  7  37   '+ ivsci +'  1 1 1 '+ nDCbusi +' ' + idDCgrid + ' '+ nDCbus_previous +' '+ nDClines_previous +' 0.005 0.00 0.0 20.0 1.2 0.0 0.0 0.0 150.0 3000.0 1.00 500 -500 200 -200 1.10 0.90 1.31 2.625 0.6287 0.2033 0.3040 195.00 220 0.1 0.1/\
        \n'+ ACbusi +' \'USRMDL\' '+ ACbusidi +' \'SPWDRD\' 4 0 3 12 2 15 	'+ dyr_text_ln1 +' 	0.100 10.0 0.1 200.0 10.0 0.00001 0.00001 0.0001 1.0 -1.0 999.45 0.2    /\
        \n'+ ACbusi +' \'USRMDL\' '+ ACbusidi +' \'SQWDRD\' 3 0 3 13 2 14 	'+ dyr_text_ln2 +' 	0.100 10.0 0.1 200.0 10.0 0.00001 0.00001 0.0001 0.4 -0.4 999.45 1.5 9999.0   /\
        '+ slack_dyr # +'\\n'+ ACbusi +' \'USRMDL\' '+ ACbusidi +' \'WDELAY\' 9 0   6 5 8 6 	1 4     5 6 11 10  0.000 0.000 0.000 0.000   0.000 /\n'

    return str_dynmodelinstance
    
def fun_write2txtfile(str_openformat, m_data, str_path, str_filename, str_idDCgridi, str_writeformat):
    
    with open(os.path.join(str_path,str_filename) ,str_openformat) as f_file:
        np.savetxt(f_file, m_data, fmt = str_writeformat, delimiter=',', header = str_idDCgridi, comments='')
    f_file.close()

def fun_setDCnetworkdata(str_path4dynamics, d_acdcoutput, l_MTDCgrids):

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
    l_nDClines = []
    l_nDCbus = []
    nMTDCgrids = len(l_MTDCgrids)   
    for i in range(0,nMTDCgrids):
        a_dcbus = l_MTDCgrids[i]["dcbus"]
        l_nDClines.append(len(l_MTDCgrids[i]["dcbranch"]))
        l_nDCbus.append(len(a_dcbus))

    v_Udc_ini = m_dcbus[:,DC_UDC]
    m_Ydc = fun_getDCconductance(DCSbase, nDCbus, nDClines, m_dcbus, m_dcbranch)
    m_Ac = fun_getDCincidence(nDCbus, nDClines, m_dcbus, m_dcbranch)
    v_Gdc = m_dcbus[:,DC_GDC].reshape(nDCbus)
    v_Cdc = m_dcbus[:,DC_CDC].reshape(nDCbus)
    v_Rsdc = m_dcbranch[:,DC_BR_R].reshape(nDClines)
    v_Ldc = m_dcbranch[:,DC_BR_L].reshape(nDClines)
    m_acdc_buses = m_converter[:,(FIL_DCBUS, FIL_ACBUS)]

    # set data for each DC grid with a given header
    str_openformat = 'a'

    nDCbus_previous = 0
    nDClines_previous = 0 
    for iMTDC in range(0,nMTDCgrids):
        
        str_idDCgrid = chr(ord('A') + iMTDC % 26)
       
        # extract data from one grid
        Ac_new = m_Ac[nDCbus_previous:nDCbus_previous+l_nDCbus[iMTDC],nDClines_previous:nDClines_previous+l_nDClines[iMTDC]]        
        Ydc_new = m_Ydc[nDCbus_previous:nDCbus_previous+l_nDCbus[iMTDC],nDCbus_previous:nDCbus_previous+l_nDCbus[iMTDC]]       
        Udc_ini_new = v_Udc_ini[nDCbus_previous:nDCbus_previous+l_nDCbus[iMTDC]]
        Gdc_new = v_Gdc[nDCbus_previous:nDCbus_previous+l_nDCbus[iMTDC]]
        Cdc_new = v_Cdc[nDCbus_previous:nDCbus_previous+l_nDCbus[iMTDC]]
        Rsdc_new = v_Rsdc[nDClines_previous:nDClines_previous+l_nDClines[iMTDC]]
        Ldc_new = v_Ldc[nDClines_previous:nDClines_previous+l_nDClines[iMTDC]]
        acdcbus_new = m_acdc_buses[nDCbus_previous:nDCbus_previous+l_nDCbus[iMTDC]]      
        
        # write the extracted data below the header str_idDCgrid
        fun_write2txtfile(str_openformat, Ac_new, str_path4dynamics, 'data_Ac.txt', str_idDCgrid, '%-7.4f')
        fun_write2txtfile(str_openformat, Ydc_new, str_path4dynamics, 'data_Ydc.txt', str_idDCgrid, '%-7.4f')
        fun_write2txtfile(str_openformat, Udc_ini_new, str_path4dynamics, 'data_Udc_ini.txt', str_idDCgrid, '%-7.8f') 
        fun_write2txtfile(str_openformat, Gdc_new, str_path4dynamics, 'data_Gdc.txt', str_idDCgrid, '%-7.8f') 
        fun_write2txtfile(str_openformat, Cdc_new, str_path4dynamics, 'data_Cdc.txt', str_idDCgrid, '%-7.8f') 
        fun_write2txtfile(str_openformat, Rsdc_new, str_path4dynamics, 'data_Rsdc.txt', str_idDCgrid, '%-7.8f') 
        fun_write2txtfile(str_openformat, Ldc_new, str_path4dynamics, 'data_Ldc.txt', str_idDCgrid, '%-7.8f') 
        fun_write2txtfile(str_openformat, acdcbus_new, str_path4dynamics, 'data_acdcbus.txt', str_idDCgrid, '%d')
        
        nDCbus_previous += l_nDCbus[iMTDC]
        nDClines_previous += l_nDClines[iMTDC]        

    
    return