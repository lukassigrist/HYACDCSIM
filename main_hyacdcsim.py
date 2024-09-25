"""
This code implements the building of the static and dynamic MTDC grids and the solutions of the resulting AC/DC load flow.

From a static point of view, MTDC grids are modelled as:

- PV/PQ generators at the AC side
- Transformers and LC filters at the AC side included within ZSORCE
- Vdc/Pdc sources at the DC side
- RL branches at the DC side

This code first calls the subroutines of a sequential AC/DC load flow. Note that the corresponding call can substitute the standard PSS/E power flow solution API to simulat AC/DC load flows. The results of the load flow can be optionally shown.

From the dynamic point of view, the dynamics of the multiple MTDC grids are modelled by now through a single dynamic model with the number of states updated according to the number of MTDC grids and its elements. Converters and their controls have individual models.

This model finally creates the data files necessary for the dynamic simulations. In particular the following classe of files are generated:

- .txt files containing the DC grid data (topology, R, L, C, etc.)
- .dyr file, updated to include the models of the converters (SVSCON), their supplementary controls (SPWDRD o SQWDRD), and the grid (DCGRID)

Inputs:
- Static MTDC data (data in PYPOWER format - .py file)
- Static load flow file (in PSS/e format - .sav file)

Authors: 
- Aurelio Garcia Cerrada
- Javier Renedo
- Carlos Prieto
- Lukas Sigrist
"""

# ---------
# Libraries
# ---------

import os, sys
# PSSE
import psse34 
import psspy, redirect
redirect.psse2py()
# Hybrid AC/DC load flow solution and dynamic simulation preparation
from _PyModules import module_acdc


# ------------------
# User-defined input
# ------------------

# file name and paths
str_lffile = r"kundur32_noAC.sav" # initial AC load flow"
str_dyrfile = r"kundur.dyr" # initial AC load flow"
str_pathlffile = r"C:\Users\lsigrist\OneDrive - Universidad Pontificia Comillas\PSSE\Tools\HYACDCSIM\Input\KundurDC" # Path
str_path4dynamics = r"C:\Users\lsigrist\OneDrive - Universidad Pontificia Comillas\PSSE\Tools\HYACDCSIM\Simulation\KundurDC" # Path

sys.path.append(str_pathlffile)
import define_grids_mtdc # MTDC grid definition (.py) - USER DEFINED

# saving and showing options
issetdynamicfiles = True # build dynamic files
issavedclfresults = False
isprogress = 2 # 6 no progress, 1 standard destination, 2 file

# AC and DC load flow
lftol = 1e-4 # lftolerance
lfmaxiter = 20 # Max number of iterations 
 
# ----------------------
# Not user-defined input
# ----------------------
 
# converter identifier (all converters have an ID and belong to the same owner
idconv = r"""14""" # use an identifier not used priorly to ease subsystem definition
idowner = 8 # idem

d_acdcoptions = dict(issetdynamicfiles=issetdynamicfiles,issavedclfresults=issavedclfresults,isprogress=isprogress,
    lftol=lftol,lfmaxiter=lfmaxiter,idconv=idconv,idowner=idowner)

str_savfileorig = os.path.join(str_pathlffile,str_lffile)  # initial AC load flow
str_dyrfileorig = os.path.join(str_pathlffile,str_dyrfile)  # initial AC load flow

str_pathdclfresults = str_pathlffile

if __name__ == "__main__":

    # add multiple MTDC dictionaries to the MTDCgrid list
    l_MTDCgrids = []
    l_MTDCgrids.append(define_grids_mtdc.MTDC_1())
    # l_MTDCgrids.append(define_grids_mtdc.MTDC_2())
    
    # initialize PSS/e
    psspy.psseinit(2000)

    # run ACDC load flow
    [d_acdcoutput, success, it, k_int_ac, k_int_dc, k_int_dcslack, j_int_dcslack, MM_Pdc_bus] = module_acdc.main_runacdclf(d_acdcoptions, 
    str_savfileorig, str_pathdclfresults, l_MTDCgrids)
    
    # show load flow results if needed
    module_acdc.fun_showacdclfoutput(isprogress, issavedclfresults, str_pathdclfresults, d_acdcoutput)  
    
    # write the files needed for dynamic simulations (.txt for DC grid) and update .dyr file
    module_acdc.main_setacdcdynamicdata(d_acdcoutput, str_dyrfileorig, str_path4dynamics, d_acdcoptions, l_MTDCgrids)  
