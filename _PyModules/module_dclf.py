"""
This module contains those functions used in the AC/DC power flow to prepare and solve the DC grid load flow

The module contains the following functions:

- fun_mainDCloadflow: principal function
- fun_getDCjacobian : builds the DC grid Jacobian matrix
- fun_getDCconductance: builds the DC grid admittance matrix
- fun_getDCincidence: builds the DC grid incidence matrix
- fun_solveDCloadflow: solves the DC load flow

- fun_slackiterationDC: iteration of the DC slack bus converter


Authors:
- Aurelio Garcia Cerrada
- Javier Renedo
- Lukas Sigrist
"""

import numpy as np
from math import pi

DC_BUS, DC_TYPE, DC_TYPE2, DC_US, DC_DELTAS, DC_PS, DC_QS, DC_UDC, DC_PDC, DC_IDC, DC_GS, DC_BS, DC_AREA, DC_BASE_KV, DC_ZONE, DC_VMAX, DC_VMIN = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)
DC_F_BUS, DC_T_BUS, DC_BR_R, DC_BR_X, DC_BR_B, DC_RATE_A, DC_RATE_B, DC_RATE_C, DC_RATIO, DC_ANGLE, DC_STATUS, DC_PIJ, DC_PJI, DC_ICCIJ = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,13)

# DC_BUS, DC_TYPE, DC_UDC, DC_PDC = (0, 1, 7, 8)
# DC_F_BUS, DC_T_BUS, DC_BR_R, DC_PIJ, DC_PJI, DC_ICCIJ = (0, 1, 2, 11, 12,13)

# Constants
j = 1j

def fun_getDCjacobian(v_udc, v_pdcbus, m_Ydc):
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
    

    nDC = np.shape(v_udc)[0] # number of DC buses

    Udc_x = v_udc[0:nDC-1]
    Pdc_x = v_pdcbus[0:nDC-1]

    m_Jdc = np.matrix(np.zeros(shape=(nDC-1,nDC-1)))

    for i in range(0,nDC-1):
        for j in range(i,nDC-1): # symmetric matrix
            if i!=j:
                m_Jdc[i,j] = Udc_x[i]*Udc_x[j]*m_Ydc[i,j]
                m_Jdc[j,i] = m_Jdc[i,j]
            else:
                m_Jdc[i,i] = Pdc_x[i] + m_Ydc[i,i]*Udc_x[i]**2
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

    for k in range(0,nDClines):
                       
        ix = np.searchsorted(v_DCbus, v_From_bus[k])
        jx = np.searchsorted(v_DCbus, v_To_bus[k])

        m_Ydc[ix, jx] = -1/v_Rdc[k]
        m_Ydc[jx, ix] = m_Ydc[ix, jx]

    for k in range(0,nDCbus):
        m_Ydc[k, k] = -m_Ydc[k,:].sum(0) # sum the terms of row k

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
    
def fun_solveDCloadflow(v_Udc_0, v_Pdc_0, lftol, lfmaxiter, nDCbus, m_Ydc):
    """ This function obtains the power flows of a DC grid. The initial values are:
    v_Udc_0: Voltages (p.u) (all buses) v_Pdc_0: Injected powers  (p.u) (all buses) m_Ydc: Ybus of the DC grid (p.u)
    it is assumed that the dc-slack bus is the last node """

    v_mismatch_DC = 0.0*v_Pdc_0

    # Boundary conditions: Pref and Uref
    v_Pdc_ref = v_Pdc_0[0:nDCbus-1] # pu
    v_Udc_bus = v_Udc_0 # vector with all the voltages

    # DC load flow 
    convergence = 0 # convergence bucle
    it = 0
    while not convergence:
        it = it + 1

        v_Pdc_bus = np.multiply(v_Udc_bus, m_Ydc*v_Udc_bus)
        v_delta_Pdc_x = v_Pdc_ref - v_Pdc_bus[0:nDCbus-1] # mismatch
        m_Jdc = fun_getDCjacobian (v_Udc_bus, v_Pdc_bus, m_Ydc)
        v_delta_U_Ux = np.linalg.inv(m_Jdc)*v_delta_Pdc_x

        v_Udc_bus[0:nDCbus-1] = np.multiply(v_Udc_bus[0:nDCbus-1], np.ones(shape=(nDCbus-1,1)) + v_delta_U_Ux ) 
        v_mismatch_DC[0:nDCbus-1] = v_delta_Pdc_x
        
        if np.max(np.abs(v_delta_Pdc_x))<=lftol:
            convergence = 1
        if it>=lfmaxiter:
            convergence = 2

    return v_Udc_bus, v_Pdc_bus, it, convergence, v_mismatch_DC

def fun_mainDCloadflow(baseMVA, lftol, lfmaxiter, m_dcbus, m_dcbranch):
    """ 
    This function solves the load flow equations of a DC grid.
    
    To do so: 
        (1) the dc slack is put in the last position
        (2) the DC load flow is solved, 
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
    indexdcslack = np.nonzero(m_dcbus[:,DC_TYPE]==2)[0] # identify in which position is the dc slack
    indexdcslack = indexdcslack[0] # if the user has typed more than one dcslack, only the first one will be considered
    
    A = np.matrix(np.eye(nDCbus))
    if indexdcslack != (nDCbus-1):
        A[indexdcslack, nDCbus-1] = 1
        A[indexdcslack, indexdcslack] = 0
        A[nDCbus-1, indexdcslack] = 1
        A[nDCbus-1, nDCbus-1] = 0

    # initial states according to modified positions
    v_Udc_0 = A*v_Udc_0_original
    v_Pdc_0 = A*v_Pdc_0_original
    m_Ydc = A*m_Ydc_original*np.linalg.inv(A)

    # DC load flow
    [v_Udc_bus, v_Pdc_bus, it, convergence, v_mismatch_DC] = fun_solveDCloadflow(v_Udc_0, v_Pdc_0, lftol, lfmaxiter, nDCbus, m_Ydc)

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

def fun_slackiterationDC(u_s, delta_s_real, q_s, p_dc, p_s0, y1, z2, y3, a, b, c_inv, c_rect, lftol, lfmaxiter):
    """ 
    This function carries out the dc lack iteration.

    The AC bus active power injection of the DC slack bus is calculated from its DC power by accounting for the converter losses. The converter lossses are however dependent of the yet unknown converter current, requiring an additional iteration. The idea is to obtain p_s such that is compatible with pdc.
    
    This is solved by an external iteration that updates the p_c according to p_dc and ploss
    
    Ploss is computed by solving an internal mini load flow that computes u_c and delta_c such that q_s (from AC load flow) and p_c (from external iteration) are maintained. Equation (57) in Beerten et al. (2012) is adopted but by making use of the PI model of the transformer and LC filter not approximation is made,

    Tick of the nodes
    s: AC side node
    c: converter node (AC side)
    Same notation as in J. Beerten et. al. (2012)
    
    INPUTS
    u_s, delta_s: Voltage (modulus and argument) of the 's' node at that iteration
    q_s: reactive power injected from the 's' node to the AC grid at that iteration
    p_dc: Power of the DC-grid of the dc-slack bus. 
    p_s0: Initial guess of the active power injected from the 's' node to the AC grid at that iteration
    y1, z2, y3: admitances and impedances of the 'pi' scheme of the converter filter 
    a, b, c_inv, c_rect: converter losses parameters
    
    OUTPUTS
    p_s_out: active power injected from the 's' node to the AC grid at that iteration
    p_c_out: active power of the node 'c'
    q_c_out: reactive power of the node 'c'
    u_c_fasor_out: complex voltage of the node 'c'
    it_ext: # of external iterations required to convergence
    vit_int: vector: # of internal iterations required to convergence at each external iteration
    convergence: YES (1); NO (0) 
    """

    # to improve convergence, refer all angles to bus s (delta_s), which is at then corrected at the very end
    delta_s = 0
    
    # build the pi equivalent of the T transformer and filter model
    # uc = Am*us + Cm*is
    # ic = Cm*us + Dm*is
    Am = 1+ z2*y1
    Bm = z2
    Cm = y1 + y3 + z2*y1*y3
    Dm = 1 + z2*y3
    
    # complex addmitance matrix [is, ic]' = Ym*[us, uc]' (from c to s)
    Y_mbus = np.matrix([[-(y1+1/z2),1/z2],[-1/z2,(1/z2+y3)]]) 

    G_mbus, B_mbus = (Y_mbus.real, Y_mbus.imag)
    
    # declare DC slack jacobian
    Jm = np.matrix(np.zeros(shape=(2,2)))

    # initialize dc slack iteration
    [u_s_phasor, s_s, u_c_phasor, s_c] = fun_initializeDCslackiteration(u_s, delta_s, p_s0, q_s, p_dc, Am, Bm, Cm, Dm, a, b, c_inv, c_rect)
    p_c = s_c.real

    #################################################################################
    # External iteration
    #################################################################################
    
    convergence = 0; # convergence loop
    k_ext = 0;
    while not convergence:

        k_ext = k_ext + 1
        p_c_previous = p_c
    
        u_s, delta_s = (abs(u_s_phasor), np.angle(u_s_phasor))        
        u_c, delta_c = (abs(u_c_phasor), np.angle(u_c_phasor))

        p_c = s_c.real
        q_s = s_s.imag

        U_mbus = j*np.matrix([[1.0],[1.0]])
        U_mbus[0] = u_s_phasor # parece que asigando asi da menos problemas 
        U_mbus[1] = u_c_phasor
        S_mbus = np.multiply(U_mbus, np.conj(Y_mbus*U_mbus))

        #################################################################################
        # Mini power flow 
        #################################################################################
        
        convergenceint = 0
        j_int = 0
        while not convergenceint:

            j_int = j_int + 1

            S_mbus = np.multiply(U_mbus, np.conj(Y_mbus*U_mbus)) # S_mbus = [(p_s + j*q_s); p_c + j*p_c]
            p_s_calc, q_s_calc = ( float(S_mbus[0].real), float(S_mbus[0].imag) )
            p_c_calc, q_c_calc = ( float(S_mbus[1].real), float(S_mbus[1].imag) )

            Delta_p_c = p_c - p_c_calc
            Delta_q_s = q_s - q_s_calc 

            v_mismatchint = np.matrix([[1.0],[1.0]]) 
            v_mismatchint[0] = Delta_p_c
            v_mismatchint[1] = Delta_q_s
            
            # set Jacobian matrix Jm = [dpc/ddeltac u_c*dpc/duc; dqs/ddeltac u_c*dqs/duc] with u_s_phasor = u_s*exp(j*0)
            Jm[0,0] = -q_c_calc - B_mbus[1,1]*u_c**2 # J_p2_delta2 "M=M(i,j)"
            Jm[0,1] = p_c_calc + G_mbus[1,1]*u_c**2 # J_p2_u2 "N=N(i,i)"
            Jm[1,0] = u_s*u_c*(-G_mbus[0,1]*np.cos(delta_s - delta_c) + B_mbus[0,1]*np.sin(delta_s - delta_c)) # J_q1_delta2 = "M=M(i,j)"
            Jm[1,1] = -u_s*u_c*(G_mbus[0,1]*np.sin(delta_s - delta_c) + B_mbus[0,1]*np.cos(delta_s - delta_c)) # J_q1_u2 "L=L(i,j)"

            # update delta2 and u_c
            Delta_X = np.linalg.inv(Jm)*v_mismatchint
            delta_c = delta_c + Delta_X[0] # update the angle
            u_c = u_c*(1.0 + Delta_X[1]) # update the voltage
            u_c_phasor = u_c*np.exp( complex(0, delta_c) )

            U_mbus[1] = u_c_phasor

            if np.max(np.abs(v_mismatchint))<=lftol:
                convergenceint = 1
            if j_int>=lfmaxiter:
                convergenceint = 2

        # re-compute current and powers    
        I_mbus = Y_mbus*U_mbus # I_mbus = [i_s, i_c]'
        S_mbus = np.multiply(U_mbus, np.conj(Y_mbus*U_mbus))

        # switch to normal notation (node s follows load convention and the flow is in the same direction as c)
        s_s = S_mbus[0]        
        s_c = S_mbus[1]
        p_c = s_c.real
        i_c_phasor = I_mbus[1]

        if p_dc<=0: # inverter
            p_loss = a + b*abs(i_c_phasor) + c_inv*abs(i_c_phasor)**2
            p_c = abs(p_dc) - p_loss
        else: # rectifier
            p_loss = a + b*abs(i_c_phasor) + c_rect*abs(i_c_phasor)**2
            p_c = -(abs(p_dc) + p_loss)

        p_c_current = p_c

        # Convergence of p_c
        if abs(p_c_current - p_c_previous)<=lftol:
            convergence = 1
        if k_ext>=lfmaxiter:
            convergence = 2
    
    # correct angles acknowledging that the reference node s has a angle delta_s_real
    delta_c = np.angle(u_c_phasor)
    delta_c_real = delta_c + delta_s_real
    u_c_phasor = np.multiply( abs(u_c_phasor), np.exp(j*delta_c_real) )

    p_s = float(s_s.real)
    p_c = float(p_c) # float instead of matrix 1x1, but is the same
    q_c = float(s_c.imag) # U_c_phasor ya esta calculado

    
    return p_s, p_c, q_c, u_c_phasor, k_ext, j_int, convergence

def fun_initializeDCslackiteration(u_s, delta_s, p_s0, q_s, p_dc, Am, Bm, Cm, Dm, a, b, c_inv, c_rect):

    ## Bus s
    p_s = p_s0
    u_s_phasor = u_s*np.exp(j*delta_s)
    s_s = p_s + j*q_s # same as: complex(p_s, q_s)
    i_s_phasor = (s_s/u_s_phasor).conj() # same to put conj(s_s/u_s_phasor) both are valids

    ## Bus c
    u_c_phasor = Am*u_s_phasor + Bm*i_s_phasor
    i_c_phasor = Cm*u_s_phasor + Dm*i_s_phasor
    s_c = u_c_phasor*np.conj(i_c_phasor)
    p_c = s_c.real
    q_c = s_c.imag

    if p_dc<=0: # inverter
        p_loss = a + b*abs(i_c_phasor) + c_inv*abs(i_c_phasor)**2
        p_c = abs(p_dc) - p_loss
    else: # rectifier
        p_loss = a + b*abs(i_c_phasor) + c_rect*abs(i_c_phasor)**2
        p_c = -(abs(p_dc) + p_loss)


    s_c = p_c + j*q_c

    return u_s_phasor, s_s, u_c_phasor, s_c

        
        

        

        
        
        
    
    
    
    






