'''
This code provides the MTDC power flow data in PYWOER format. 

As many grids can be defined as needed by using a def statement and following the PYPOWER format.

So far
- the bus numbering of DC grids is: 1, 2, 3, 4, etc.
- only DC buses with a converter are considered. Although DC buses without converter work well for the static analysis, the DCbus matrix should be re-organized (DC buses with VSC, then others) for the dynamic analysis

Author: 
- Javier Renedo

To do:
- Confirm generalized DC bus numbering
- Convert pu input data on invidual ratings to pu on common rating
- Implement DC buses without VSC


Copyright (C) 1996-2011 Power System Engineering Research Center
Copyright (C) 2010-2011 Richard Lincoln

PYPOWER is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License,
or (at your option) any later version.

PYPOWER is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with PYPOWER. If not, see <http://www.gnu.org/licenses/>.
'''

from numpy import array

def MTDC_1():
    """
    Kundur case with an MTDC connected to buses 5, 7, and 11
    Typical values for DC grids (in SI units, DOI: 10.1109/PowerTech46648.2021.9494939)
        Ldcij = 1 mH/km = 1e-3 H/km
        Rdcij = 0.02 Ohm/km
        Cdcij = 0.02 uF/km = 0.02e-6 F/km
        Cdci = 195 uF  = 195e-6 F (Cdci = Cvsci + sum(j,Cdcij/2))
        
    Typical values for transformer and filter (here in pu of the converter, DOI: 10.1109/PowerTech46648.2021.9494939))
        rt = 0.01
        xt = 0.1
        bf = 0
        rc = 0.005
        xc = 0.1
     
    Typical values for losses (here in pu of the converters, DOI: 10.1109/PowerTech46648.2021.9494939))
        a = 10e-3
        b = 3e-3
        crect = 4.5e-3
        cinv = 6.5e-3
    
    DC per unit (PhD Beerten, 2013)
        Sdc,b = Sac,b
        Zdc,b = Udc,b^2/(Sdc,b/p)
        Idc,b = (Sdc,b/p)/Udc,b
        
        where Udc,b = Uptg,dc,N (the nominal pole-to-ground voltage of the DC grid) and where p = 1 for a monopolar scheme 
        and p = 2 for a bipolar or symmetrically grounded monopolar scheme
        
        Note that the Udc,b is not independent of the Uac,b. The nominal DC voltage should give rise to nominal fundamental AC
        voltage without overmodulation for instance (m = 1):
        
        Uac,b = m*0.866/sqrt(2)*Udc,b = sqrt(3/2)/2*Udc,b if m = 1
        
        rdc = Rdc/Zdc,b [pu]
        ldc = Ldc/Zdc,b [pu*s]
        cdc = Cdc*Zdc,b [pu*s]
        
    """
    ppc = {"version": '1'}

    ##-----  Power Flow Data  -----##
    ## system MVA base - all pu values (rt, xt, a, b, Ldc, Ccc, Cdc, etc.) are expressed in this rating
    ppc["baseMVA"] = 100.0

    ## converters and AC/DC coupling data
    # f_DC_bus, t_AC_bus, rt, xt, bf, rc, xc, a, b, crect, cinv, Pmax, Pmin, Qmax, Qmin, rateA, rateB, rateC, ratio, angle, status, angmin, angmax
    ppc["converter"] = array([
        [1, 7,  0.01*100/500, 0.1*100/500, 0.0000, 0.005*100/500, 0.1*100/500, 0.0, 0.0, 0.0, 0.0, 9999, -9999, 9999, -9999, 500, 0, 0, 0, 0, 1, -360, 360],
        [2, 9,  0.01*100/500, 0.1*100/500, 0.0000, 0.005*100/500, 0.1*100/500, 0.0, 0.0, 0.0, 0.0, 9999, -9999, 9999, -9999, 500, 0, 0, 0, 0, 1, -360, 360],
#        [3, 11, 0.01*100/500, 0.1*100/500, 0.0000, 0.005*100/500, 0.1*100/500, 10e-3*500/100, 3e-3*500/100, 4.5e-3*500/100, 6.5e-3*500/100, 9999, -9999, 9999, -9999, 500, 0, 0, 0, 0, 1, -360, 360],
        ])

    ## DC grid data: buses and branches
    ## DC bus data
    # bus_i	type1	type2    Us     delta_s     Ps      Qs      Udc     Pdc_iny    Idc  Gdc      Cdc     area          baseKV      zone        Vmax        Vmin
    # type1: 1: node P, 2: dc-slack
    # type2: 1: control of Qs, 2: control of u_s
    ppc["dcbus"] = array([
        [1, 2, 1, 1.00, 0, -400, 0.0, 1.00, 0, 0, 0, 0.39936, 1, 320, 1, 1.06, 0.94],
        [2, 1, 1, 1.00, 0,  400, 0.0, 1.00, 0, 0, 0, 0.39936, 1, 320, 1, 1.06, 0.94],
#        [3, 1, 1, 1.00, 0,  45, 0.0, 1.00, 0, 0, 0, 0.39936, 1, 320, 1, 1.06, 0.94],
        ])

    ## DC branch data
  
    linelength = 200 # km
    #	fbus	tbus	r       Ldc       Ccc       rateA	rateB	rateC	ratio	angle	status   
    ppc["dcbranch"] = array([
	[1,       2,       9.77*1e-6*linelength,   4.88*1e-7*linelength,      4.1*1e-5*linelength,      1000,     0,      0,      0,       0,       1,         0,      0,       0],
#	[1,       3,       9.77*1e-6*linelength,   4.88*1e-7*linelength,      4.1*1e-5*linelength,      1000,     0,      0,      0,       0,       1,         0,      0,       0],
#    [2,       3,       9.77*1e-6*linelength,   4.88*1e-7*linelength,      4.1*1e-5*linelength,      1000,     0,      0,      0,       0,       1,         0,      0,       0],
	])


    return ppc

########################################################################

def MTDC_2():
    """
    Caso Kundur
    Con MTDC conectado a (6), (8) Y (10)
    """
    ppc = {"version": '1'}

    ##-----  Power Flow Data  -----##
    ## system MVA base
    ppc["baseMVA"] = 100.0

    ## converters and AC/DC coupling data
    #   f_DC_bus, t_AC_bus, rt, xt, bf, rc, xc, a, b, crect, cinv, Pmax, Pmin, Qmax, Qmin, rateA, rateB, rateC, ratio, angle, status, angmin, angmax
    ppc["converter"] = array([
        [4, 6,  0.0001, 0.0001, 0.0000, 0.004, 0.04, 26.25e-3, 1.65e-3, 4.2e-4, 6.28e-4, 9999, -9999, 999, -999, 500, 0, 0, 0, 0, 1, -360, 360],
        [5, 8,  0.0001, 0.0001, 0.0000, 0.004, 0.04, 26.25e-3, 1.65e-3, 4.2e-4, 6.28e-4, 9999, -9999, 999, -999, 500, 0, 0, 0, 0, 1, -360, 360],
        [6, 10, 0.0001, 0.0001, 0.0000, 0.004, 0.04, 26.25e-3, 1.65e-3, 4.2e-4, 6.28e-4, 9999, -9999, 999, -999, 500, 0, 0, 0, 0, 1, -360, 360],
        ])

    ## DC grid data: buses and branches
    ## DC bus data
    #	bus_i	type1	type2    Us     delta_s     Ps      Qs      Udc     Pdc_iny    Idc  Gdc      Cdc     area          baseKV      zone        Vmax        Vmin
    # type1: 1: node P, 2: dc-slack
    # type2: 1: control of Qs, 2: control of u_s
    ppc["dcbus"] = array([
        [4, 2, 1, 1.00, 0, -90, 0.0, 1.00, 0, 0, 0, 0.39936, 1, 220, 1, 1.06, 0.94],
        [5, 1, 1, 1.00, 0, 45, 0.0, 1.00, 0, 0, 0, 0.39936, 1, 220, 1, 1.06, 0.94],
        [6, 1, 1, 1.00, 0, 45, 0.0, 1.00, 0, 0, 0, 0.39936, 1, 220, 1, 1.06, 0.94]
        ])

    ## DC branch data
    #	fbus	tbus	r       Ldc       Ccc       rateA	rateB	rateC	ratio	angle	status   
    ppc["dcbranch"] = array([
	[4,        5,       2.06612*1e-4*120.0,    1.2981788*1e-3*120.0,       3.34517*1e-3*120.0,       1000,     0,      0,      0,       0,       1,         0,      0,       0],
	[4,        6,       2.06612*1e-4*240.0,    1.2981788*1e-3*240.0,       3.34517*1e-3*240.0,       1000,     0,      0,      0,       0,       1,         0,      0,       0],
    [5,        6,       2.06612*1e-4*120.0,    1.2981788*1e-3*120.0,       3.34517*1e-3*120.0,       1000,     0,      0,      0,       0,       1,         0,      0,       0]
	])


    return ppc

########################################################################

def MTDC_3():
    """
    Caso Kundur
    Con MTDC conectado a (6), (8) Y (10) IDEAL 
    Con MTDC conectado a (5), (7), (11)
    """
    ppc = {"version": '1'}

    ##-----  Power Flow Data  -----##
    ## system MVA base
    ppc["baseMVA"] = 100.0

    ## bus data
    # bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin


    ## generator data
    # bus, Pg, Qg, Qmax, Qmin, Vg, mBase, status, Pmax, Pmin, Pc1, Pc2,
    # Qc1min, Qc1max, Qc2min, Qc2max, ramp_agc, ramp_10, ramp_30, ramp_q, apf


    ## branch data
    # fbus, tbus, r, x, b, rateA, rateB, rateC, ratio, angle, status, angmin, angmax

    ## converters and AC/DC coupling data
    #   f_DC_bus, t_AC_bus, rt, xt, bf, rc, xc, a, b, crect, cinv, Pmax, Pmin, Qmax, Qmin, rateA, rateB, rateC, ratio, angle, status, angmin, angmax
    ppc["converter"] = array([
        [1, 6, 0.0001, 0.0001, 0.0000, 0.004, 0.04, 26.25e-3, 1.65e-3, 4.2e-4, 6.28e-4, 9999, -9999, 999, -999, 500, 0, 0, 0, 0, 1, -360, 360],
        [2, 8, 0.0001, 0.0001, 0.0000, 0.004, 0.04, 26.25e-3, 1.65e-3, 4.2e-4, 6.28e-4, 9999, -9999, 999, -999, 500, 0, 0, 0, 0, 1, -360, 360],
        [3, 10, 0.0001, 0.0001, 0.0000, 0.004, 0.04, 26.25e-3, 1.65e-3, 4.2e-4, 6.28e-4, 9999, -9999, 999, -999, 500, 0, 0, 0, 0, 1, -360, 360],
        ])

    ## DC grid data: buses and branches
    ## DC bus data
    #	bus_i	type1	type2    Us     delta_s     Ps      Qs      Udc     Pdc_iny    Idc  Gdc      Cdc     area          baseKV      zone        Vmax        Vmin
    # type1: 1: node P, 2: dc-slack
    # type2: 1: control of Qs, 2: control of u_s
    ppc["dcbus"] = array([
        [1, 2, 1, 1.00, 0, -30, 0.0, 1.00, 0, 0, 0, 0.39936, 1, 320, 1, 1.06, 0.94],
        [2, 1, 1, 1.00, 0, -30, 0.0, 1.00, 0, 0, 0, 0.39936, 1, 320, 1, 1.06, 0.94],
        [3, 1, 1, 1.00, 0,  60, 0.0, 1.00, 0, 0, 0, 0.39936, 1, 320, 1, 1.06, 0.94]
        ])

    ## DC branch data
    #	fbus	tbus	r       Ldc       Ccc       rateA	rateB	rateC	ratio	angle	status   
    ppc["dcbranch"] = array([
	[1,        2,       0.001*125.0/100.0,    6.84*1e-5*125.0/100.0,       0,       1000,     0,      0,      0,       0,       1,         0,      0,       0],
	[1,        3,       0.001*250.0/100.0,    6.84*1e-5*250.0/100.0,       0,       1000,     0,      0,      0,       0,       1,         0,      0,       0],
        [2,        3,       0.001*125.0/100.0,    6.84*1e-5*125.0/100.0,       0,       1000,     0,      0,      0,       0,       1,         0,      0,       0]
	])


    return ppc

if __name__ == "__main__":
    print(r"Module.")