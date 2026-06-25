"""
This module contains those functions used in the AC/DC power flow to prepare and solve the DC grid power flow

The module contains the following functions:

- fun_mainDCloadflow: principal function
- fun_getDCjacobian : builds the DC grid Jacobian matrix
- fun_getDCconductance: builds the DC grid admittance matrix
- fun_getDCincidence: builds the DC grid incidence matrix
- fun_solveDCloadflow: solves the DC power flow

- fun_slackiterationDC: iteration of the DC slack bus converter


Authors:
- Aurelio Garcia Cerrada
- Javier Renedo
- Lukas Sigrist
- Saeed Rezaeian-Marjani
"""

import numpy as np
from math import pi

DC_BUS, DC_TYPE, DC_TYPE2, DC_US, DC_DELTAS, DC_PS, DC_QS, DC_UDC, DC_PDC, DC_IDC, DC_GS, DC_BS, DC_AREA, DC_BASE_KV, DC_ZONE, DC_VMAX, DC_VMIN = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)
DC_F_BUS, DC_T_BUS, DC_BR_R, DC_BR_X, DC_BR_B, DC_RATE_A, DC_RATE_B, DC_RATE_C, DC_RATIO, DC_ANGLE, DC_STATUS, DC_PIJ, DC_PJI, DC_ICCIJ = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,13)

# DC_BUS, DC_TYPE, DC_UDC, DC_PDC = (0, 1, 7, 8)
# DC_F_BUS, DC_T_BUS, DC_BR_R, DC_PIJ, DC_PJI, DC_ICCIJ = (0, 1, 2, 11, 12,13)

# Constants
j = 1j

def fun_getDCjacobian(v_udc, v_pdcbus, m_Ydc, last_indices):
    """ 
    This function builds the Jacobian matrix of a DC grid. It is assumed that each vector and matrix are numpy matrix type
     
    Inputs:
    v_udc: Voltages vector (all the buses of the grid)
    v_pdcbus: Powers vector (all the buses of the grid)
    m_Ydc: Admittance bus matrix of the DC grid 
     
    Output:
    m_Jdc: DC Jacobian matrix
    """

    v_udc = np.matrix(v_udc) # i want matrices, not arrays
    v_pdcbus = np.matrix(v_pdcbus)
    m_Ydc = np.matrix(m_Ydc)   
    v_udc = 1.0*v_udc # ensure that the numbers are float type and not int
    v_pdcbus = 1.0*v_pdcbus
    m_Ydc = 1.0*m_Ydc

    m_Ydc_reduced = np.delete(m_Ydc, last_indices, axis=0)  
    m_Ydc_reduced = np.delete(m_Ydc_reduced, last_indices, axis=1)  

    nDC = np.shape(v_udc)[0] # number of DC buses

    Udc_x = np.delete(v_udc, last_indices, axis=0) 
    Pdc_x = np.delete(v_pdcbus, last_indices, axis=0)
    
    m_Jdc = np.matrix(np.zeros(shape=(nDC-len(last_indices),nDC-len(last_indices))))
     
    for i in range(0,nDC-len(last_indices)):
        for j in range(i,nDC-len(last_indices)): # symmetric matrix
            if i!=j:
                m_Jdc[i,j] = Udc_x[i]*Udc_x[j]*m_Ydc_reduced[i,j]
                m_Jdc[j,i] = m_Jdc[i,j]
            else:
                m_Jdc[i,i] = Pdc_x[i] + m_Ydc_reduced[i,i]*Udc_x[i]**2

    return m_Jdc

def fun_getDCconductance(baseMVA, nDCbus, nDClines, m_dcbus, m_dcbranch):
    """ 
    This functions builds the bus admittance matrix of a DC grid, m_Ydc. 
    
    Inputs:
    baseMVA: Basis of the system
    m_dcbus: structure of the DC bus data in PYPOWER format
    m_dcbranch: structure of the DC branch data in PYPOWER format
    
    Output:
    m_Ydc: DC conductance matrix
    """

    v_From_bus = np.array(m_dcbranch[:,0]).reshape(-1,)
    v_To_bus = np.array(m_dcbranch[:,1]).reshape(-1,)
    v_Rdc = m_dcbranch[:,2]
        
    v_DCbus = np.array(m_dcbus[:,0]).reshape(-1,) # DC bus number
    nDCbus = np.shape(m_dcbus[:,0])[0]
    nDClines = np.shape(v_From_bus)[0]
    m_Ydc = np.zeros(shape=(nDCbus,nDCbus))

    v_Gshunt = m_dcbus[:, DC_GS]
    
    for k in range(0,nDClines):
                       
        ix = np.searchsorted(v_DCbus, v_From_bus[k])
        jx = np.searchsorted(v_DCbus, v_To_bus[k])

        m_Ydc[ix, jx] = -1/v_Rdc[k]
        m_Ydc[jx, ix] = m_Ydc[ix, jx]

    for k in range(0,nDCbus):
        m_Ydc[k, k] = -m_Ydc[k,:].sum(0) + v_Gshunt[k] # sum the terms of row k
    
    return m_Ydc

def fun_getDCincidence(nDCbus, nDClines, m_dcbus, m_dcbranch):
    """ 
    This function builds the incidence matrix of a DC grid, m_Ac. The matrix is given by:
        m_Ac(i,j) = +1 if i and j are connected and the line is Lij
        m_Ac(i,j) = -1 if i and j are connected and the line is Lji
        m_Ac(i,j) = 0 if i and j are not connected 
    
    Inputs:
    m_dcbus: structure of the DC bus data in PYPOWER format
    m_dcbranch: structure of the DC branch data in PYPOWER format
    
    Output: 
    m_Ac
    """

    v_From_bus = np.array(m_dcbranch[:,0]).reshape(-1,)
    v_To_bus = np.array(m_dcbranch[:,1]).reshape(-1,) 
    v_DCbus = np.array(m_dcbus[:,0]).reshape(-1,) # DC bus number
        
    m_Ac = np.zeros(shape=(nDCbus,nDClines))

    for k in range(0,nDClines):
        
        ix = np.searchsorted(v_DCbus, v_From_bus[k])
        jx = np.searchsorted(v_DCbus, v_To_bus[k])

        m_Ac[ix,k] = 1.0
        m_Ac[jx,k] = -1.0

    return m_Ac
    
def fun_solveDCloadflow(v_Udc_0, v_Pdc_0, lftol, lfmaxiter, nDCbus, m_Ydc, last_indices):
    """ This function obtains the power flows of a DC grid. The initial values are:
    v_Udc_0: Voltages (p.u) (all buses) v_Pdc_0: Injected powers  (p.u) (all buses) m_Ydc: Ybus of the DC grid (p.u)
    it is assumed that the dc-slack bus is the last node """
    
    v_mismatch_DC = 0.0*v_Pdc_0

    # Boundary conditions: Pref and Uref
    v_Pdc_ref = np.delete(v_Pdc_0, last_indices, axis=0) # pu
    v_Udc_bus = v_Udc_0 # vector with all the voltages
    
    # DC power flow 
    convergence = 0 # convergence bucle
    it = 0
    while not convergence:
        it = it + 1
        
        v_Pdc_bus = np.multiply(v_Udc_bus, m_Ydc*v_Udc_bus)
        v_delta_Pdc_x = v_Pdc_ref - np.delete(v_Pdc_bus, last_indices, axis=0) # mismatch
        m_Jdc = fun_getDCjacobian (v_Udc_bus, v_Pdc_bus, m_Ydc, last_indices)    
        v_delta_U_Ux = np.dot(np.linalg.inv(m_Jdc), v_delta_Pdc_x)
                
        v_Udc_bus[np.delete(np.arange(nDCbus), last_indices)] = np.multiply(v_Udc_bus[np.delete(np.arange(nDCbus), last_indices)], np.ones(shape=(nDCbus-len(last_indices),1)) + v_delta_U_Ux )
        v_mismatch_DC[np.delete(np.arange(nDCbus), last_indices)] = v_delta_Pdc_x
        
        if np.max(np.abs(v_delta_Pdc_x))<=lftol:
            convergence = 1
        if it>=lfmaxiter:
            convergence = 2
    
    return v_Udc_bus, v_Pdc_bus, it, convergence, v_mismatch_DC

def fun_mainDCloadflow(baseMVA, lftol, lfmaxiter, m_dcbus, m_dcbranch, MTDC_network_ids):
    """ 
    This function solves the power flow equations of a DC grid.
    
    To do so: 
        (1) the dc slack is put in the last position
        (2) the DC power flow is solved, 
        (3) undo the position change of the dc slack
    """
    v_DCbus = np.array(m_dcbus[:,DC_BUS]).reshape(-1,)
    nDCbus = np.shape(v_DCbus)[0]
    nDClines = np.shape(m_dcbranch[:,0])[0]
    
    v_Rdc = np.matrix(m_dcbranch[:, DC_BR_R]) # p.u

    # Build the admittance matrices
    m_Ydc_original = fun_getDCconductance(baseMVA, nDCbus, nDClines, m_dcbus, m_dcbranch)
    
    # initial states as pero original positions
    v_Udc_0_original = m_dcbus[:, DC_UDC]
    v_Pdc_0_original = m_dcbus[:, DC_PDC]/baseMVA
    
    # change dc slack position  
    A = np.zeros((nDCbus, nDCbus))
    start_idx = 0
    for network_id in np.unique(MTDC_network_ids):
        
        DC_network_buses = m_dcbus[MTDC_network_ids == network_id]
        nDCbus_networks = DC_network_buses.shape[0]
        indexdcslack = np.nonzero(DC_network_buses[:, DC_TYPE] == 2)[0]
        Aa = np.matrix(np.eye(nDCbus_networks))
        if indexdcslack != (nDCbus_networks - 1):
            Aa[indexdcslack, nDCbus_networks-1] = 1
            Aa[indexdcslack, indexdcslack] = 0
            Aa[nDCbus_networks-1, indexdcslack] = 1
            Aa[nDCbus_networks-1, nDCbus_networks-1] = 0

        A[start_idx:start_idx + nDCbus_networks, start_idx:start_idx + nDCbus_networks] = Aa
        start_idx += nDCbus_networks

    # initial states according to modified positions
    v_Udc_0 = A*v_Udc_0_original
    v_Pdc_0 = A*v_Pdc_0_original

    m_Ydc = np.dot(np.dot(A, m_Ydc_original), np.linalg.inv(A))

    unique_ids, indices = np.unique(MTDC_network_ids, return_inverse=True)
    last_indices = np.array([np.max(np.where(indices == i)) for i in range(len(unique_ids))])
    
    # DC power flow
    [v_Udc_bus, v_Pdc_bus, it, convergence, v_mismatch_DC] = fun_solveDCloadflow(v_Udc_0, v_Pdc_0, lftol, lfmaxiter, nDCbus, m_Ydc, last_indices)
    
    # Final states: original numeration again
    v_Udc_original = A*v_Udc_bus
    v_Pdc_original = A*v_Pdc_bus
    v_mismatch_DC_original = A*v_mismatch_DC
    
    # compute branch flow and currents
    v_indexfrom = np.searchsorted(v_DCbus, np.array(m_dcbranch[:, DC_F_BUS]).reshape(-1,))
    v_indexto = np.searchsorted(v_DCbus, np.array(m_dcbranch[:, DC_T_BUS]).reshape(-1,))

    Pdc_original_ij = np.multiply( 1/v_Rdc, np.multiply( v_Udc_original[v_indexfrom], v_Udc_original[v_indexfrom] - v_Udc_original[v_indexto] ) )
    Pdc_original_ji = np.multiply( 1/v_Rdc, np.multiply( v_Udc_original[v_indexto], v_Udc_original[v_indexto] - v_Udc_original[v_indexfrom] ) )
    Icc_original_ij = np.multiply( 1/v_Rdc,  v_Udc_original[v_indexfrom] - v_Udc_original[v_indexto]  )
    
    # put the results in the PYPOWER format data structure
    m_dcbus[:, DC_UDC] = v_Udc_original
    m_dcbus[:, DC_PDC] = v_Pdc_original*baseMVA
    m_dcbranch[:, DC_PIJ] = Pdc_original_ij*baseMVA
    m_dcbranch[:, DC_PJI] = Pdc_original_ji*baseMVA
    m_dcbranch[:, DC_ICCIJ] = Icc_original_ij # p.u-system
    
    return m_dcbus, m_dcbranch, it, v_mismatch_DC_original

def fun_slackiterationDC(u_f, delta_f_real, q_o, p_dc, p_o0, y1, z2, y3, a, b, c_inv, c_rect, lftol, lfmaxiter):
    """ 
    This function carries out the dc lack iteration.

    The AC bus active power injection of the DC slack bus is calculated from its DC power by accounting for the converter losses. The converter lossses are however dependent of the yet unknown converter current, requiring an additional iteration. The idea is to obtain p_s such that is compatible with pdc.
    
    This is solved by an external iteration that updates the p_c according to p_dc and ploss
    
    Ploss is computed by solving an internal mini power flow that computes u_c and delta_c such that q_o (from AC power flow) and p_c (from external iteration) are maintained. Equation (57) in Beerten et al. (2012) is adopted but by making use of the PI model of the transformer and LC filter not approximation is made,

    Tick of the nodes
    f: AC side node (The corresponding powers  are indicated by the index 'o', for example: p_o and q_o).
    c: converter node (AC side)
    Same notation as in J. Beerten et. al. (2012)
    
    INPUTS
    u_f, delta_f: Voltage (modulus and argument) of the 'f' node at that iteration
    q_o: reactive power injected from the 'f' node to the AC grid at that iteration
    p_dc: Power of the DC-grid of the dc-slack bus. 
    p_o0: Initial guess of the active power injected from the 'o' node to the AC grid at that iteration
    y1, z2, y3: admitances and impedances of the 'pi' scheme of the converter filter 
    a, b, c_inv, c_rect: converter losses parameters
    
    OUTPUTS
    p_o_out: active power injected from the 'f' node to the AC grid at that iteration
    p_c_out: active power of the node 'c'
    q_c_out: reactive power of the node 'c'
    u_c_fasor_out: complex voltage of the node 'c'
    it_ext: # of external iterations required to convergence
    vit_int: vector: # of internal iterations required to convergence at each external iteration
    convergence: YES (1); NO (0) 
    """
    
    # to improve convergence, refer all angles to bus s (delta_s), which is at then corrected at the very end
    delta_f = 0
    
    # build the pi equivalent of the T transformer and filter model
    # uc = Am*us + Cm*is
    # ic = Cm*us + Dm*is
    Am = 1 + np.multiply(z2, y1)
    Bm = z2
    Cm = y1 + y3 + np.multiply(np.multiply(z2, y1), y3)
    Dm = 1 + np.multiply(z2, y3)

    # complex addmitance matrix [is, ic]' = Ym*[us, uc]' (from c to s)
    N = len(y1)
    Y_mbus = np.zeros((2*N, 2*N), dtype=complex)

    for i in range(N):
        y_mbus = np.matrix([[-(y1[i,0] + 1/z2[i,0]), 1/z2[i,0]],[-1/z2[i,0], (1/z2[i,0] + y3[i,0])]])

        Y_mbus[2*i:2*i+2, 2*i:2*i+2] = y_mbus

    G_mbus, B_mbus = (Y_mbus.real, Y_mbus.imag)
   
    # declare DC slack jacobian
    Jm = np.matrix(np.zeros(shape=(2,2)))

    # initialize dc slack iteration
    [u_f_phasor, s_o, u_c_phasor, s_c] = fun_initializeDCslackiteration(u_f, delta_f, p_o0, q_o, p_dc, Am, Bm, Cm, Dm, a, b, c_inv, c_rect)
    p_c = s_c.real
    
    #################################################################################
    # External iteration
    #################################################################################
    
    convergence = 0; # convergence loop
    k_ext = 0;
    while not convergence:

        k_ext = k_ext + 1
        p_c_previous = p_c
      
        u_f, delta_f = (abs(u_f_phasor), np.angle(u_f_phasor))        
        u_c, delta_c = (abs(u_c_phasor), np.angle(u_c_phasor))

        p_c = s_c.real
        q_o = s_o.imag
        
        num_buses = u_f_phasor.shape[0] 

        U_mbus = np.zeros((2 * num_buses, 1), dtype=complex)

        U_mbus[0:2 * num_buses:2] = u_f_phasor  
        U_mbus[1:2 * num_buses:2] = u_c_phasor  

        S_mbus = np.multiply(U_mbus, np.conj(Y_mbus.dot(U_mbus)))
     
        #################################################################################
        # Mini power flow 
        #################################################################################
        
        convergenceint = 0
        j_int = 0
        while not convergenceint:

            j_int = j_int + 1

            S_mbus = np.multiply(U_mbus, np.conj(Y_mbus.dot(U_mbus)))
            p_o_calc, q_o_calc = (S_mbus[0:2*num_buses:2, 0].real.astype(float),S_mbus[0:2*num_buses:2, 0].imag.astype(float))
            p_c_calc, q_c_calc = (S_mbus[1:2*num_buses:2, 0].real.astype(float),S_mbus[1:2*num_buses:2, 0].imag.astype(float))

            Delta_p_c = p_c - (np.matrix(p_c_calc)).T
            Delta_q_o = q_o - (np.matrix(q_o_calc)).T
            
            v_mismatchint = np.zeros((2 * num_buses, 1))
            v_mismatchint[0:2 * num_buses:2] = Delta_p_c  
            v_mismatchint[1:2 * num_buses:2] = Delta_q_o

            u_f = np.asarray(u_f)  
            u_c = np.asarray(u_c)
            delta_f = np.array(delta_f)
            delta_c = np.array(delta_c)

            Jm = np.zeros((num_buses * 2, num_buses * 2))

            idx = np.arange(num_buses) * 2
            Jm[idx, idx] = -q_c_calc - B_mbus[idx + 1, idx + 1] * u_c.flatten()**2
            Jm[idx, idx + 1] = p_c_calc + G_mbus[idx + 1, idx + 1] * u_c.flatten()**2
            Jm[idx + 1, idx] = u_f.flatten() * u_c.flatten() * (-G_mbus[idx, idx + 1] * np.cos(delta_f.flatten() - delta_c.flatten()) + B_mbus[idx, idx + 1] * np.sin(delta_f.flatten() - delta_c.flatten()))
            Jm[idx + 1, idx + 1] = -u_f.flatten() * u_c.flatten() * (G_mbus[idx, idx + 1] * np.sin(delta_f.flatten() - delta_c.flatten()) + B_mbus[idx, idx + 1] * np.cos(delta_f.flatten() - delta_c.flatten()))

            # update delta2 and u_c
            Delta_X =  np.dot(np.linalg.inv(Jm), v_mismatchint)
            delta_c = delta_c + Delta_X[0::2] # update the angle
            u_c = u_c*(1.0 + Delta_X[1::2]) # update the voltage
            u_c_phasor = u_c * np.exp(1j * delta_c)
            
            U_mbus[1::2] = u_c_phasor

            if np.max(np.abs(v_mismatchint))<=lftol:
                convergenceint = 1
            if j_int>=lfmaxiter:
                convergenceint = 2

        # re-compute current and powers
        I_mbus = np.dot(Y_mbus, U_mbus)
        S_mbus = np.multiply(U_mbus, np.conj(Y_mbus.dot(U_mbus)))
        
        # switch to normal notation (node s follows load convention and the flow is in the same direction as c)
        s_o = S_mbus[0::2]        
        s_c = S_mbus[1::2]
        p_c = s_c.real
        i_c_phasor = I_mbus[1::2]
        
        p_loss = a + np.multiply(b, np.abs(i_c_phasor)) + np.where(
        p_dc <= 0, np.multiply(c_inv, np.multiply(np.abs(i_c_phasor), np.abs(i_c_phasor))),   # inverter
        np.multiply(c_rect, np.multiply(np.abs(i_c_phasor), np.abs(i_c_phasor))))             # rectifier

        p_c = np.where(
        p_dc <= 0,
        np.abs(p_dc) - p_loss,         # inverter
        - (np.abs(p_dc) + p_loss))     # rectifier
         
        p_c_current = p_c
        
        if np.max(abs(p_c_current - p_c_previous) <= lftol):
            convergence = 1
        if k_ext>=lfmaxiter:
            convergence = 2

    # correct angles acknowledging that the reference node s has a angle delta_s_real
    delta_c = np.angle(u_c_phasor)
    delta_c_real = delta_c + delta_f_real
    u_c_phasor = np.multiply( abs(u_c_phasor), np.exp(j*delta_c_real) )
    
    p_o = s_o.real.astype(float)
    p_c = p_c.astype(float)
    q_c = s_c.imag.astype(float)
    
    return p_o, p_c, q_c, u_c_phasor, k_ext, j_int, convergence

def fun_initializeDCslackiteration(u_f, delta_f, p_o0, q_o, p_dc, Am, Bm, Cm, Dm, a, b, c_inv, c_rect):

    ## Bus s
    p_o = p_o0
    u_f_phasor = u_f*np.exp(j*delta_f)
    s_o = p_o + j*q_o # same as: complex(p_o, q_o)
    i_f_phasor = (s_o/u_f_phasor).conj() # same to put conj(s_s/u_s_phasor) both are valids

    ## Bus c
    u_c_phasor = np.multiply(Am, u_f_phasor) + np.multiply(Bm, i_f_phasor)
    i_c_phasor = np.multiply(Cm, u_f_phasor) + np.multiply(Dm, i_f_phasor)
    s_c = np.multiply(u_c_phasor, np.conj(i_c_phasor))
    p_c = s_c.real
    q_c = s_c.imag
   
    p_loss = a + np.multiply(b, np.abs(i_c_phasor)) + np.where(
        p_dc <= 0, np.multiply(c_inv, np.multiply(np.abs(i_c_phasor), np.abs(i_c_phasor))),   # inverter
        np.multiply(c_rect, np.multiply(np.abs(i_c_phasor), np.abs(i_c_phasor))))             # rectifier

    p_c = np.where(
        p_dc <= 0,
        np.abs(p_dc) - p_loss,         # inverter
        - (np.abs(p_dc) + p_loss))     # rectifier

    s_c = p_c + j*q_c
   
    return u_f_phasor, s_o, u_c_phasor, s_c
