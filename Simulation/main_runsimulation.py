# ---------
# Libraries
# ---------
import os, sys
# PSS/E
import psse34
import psspy, redirect
from psspy import _i, _f, _s, _o
import dyntools
# Plotting
from matplotlib import pyplot as plt
from matplotlib.font_manager import FontProperties
# Math
import numpy as np

# ------------------
# User-defined input
# ------------------

# file name and paths
str_pathinputfiles = r"C:\Users\lsigrist\OneDrive - Universidad Pontificia Comillas\PSSE\Tools\HYACDCSIM\Input\KundurDC" # Path
str_pathsimfiles = r"C:\Users\lsigrist\OneDrive - Universidad Pontificia Comillas\PSSE\Tools\HYACDCSIM\Simulation\KundurDC"  # Path
str_lffile = r"kundur32_noAC_MTDCg.sav" # initial AC load flow"
str_dyrfile = r"kundur_MTDCg.dyr" # initial AC load flow"
str_dllfile = r"MTDCDyn.dll"

# Solver parameters
TSTEP = 0.001
TYSLACCFAC = 0.9
TYSLTOL = 0.0001
TYSLMAXITER = 75
TINI = 1.1      
TFINAL = 12

# ------------
# Preparations
# ------------
str_pathlffile = os.path.join(str_pathinputfiles, str_lffile)
str_pathdyrfile = os.path.join(str_pathinputfiles, str_dyrfile)
str_pathdllfile = os.path.join(str_pathsimfiles, str_dllfile)
str_pathoutfile = os.path.join(str_pathsimfiles, r"output.out")
str_pathconecfile = os.path.join(str_pathsimfiles, r"conec")
str_pathconetfile = os.path.join(str_pathsimfiles, r"conet")
str_pathcompilefile = os.path.join(str_pathsimfiles, r"compile")

# Set working directory to the simulation folder,
# where the DLL, Fortran-related files, and auxiliary .txt files are located
print("Current Working Directory:", os.getcwd())

# ----------
# Simulation
# ----------
if __name__ == "__main__":
    redirect.psse2py()
   
    psspy.psseinit(2000)
    psspy.progress_output(2,"report",[0,0])

    # Open file and run load flow
    psspy.case(str_pathlffile)
    psspy.fnsl([0,0,0,1,0,0,99,0])
    
    # Preparing dynamic simulation
    # ----------------------------
    # Convert generator and loads
    psspy.cong(0)
    psspy.conl(0,1,1,[0,0],[ 100.0,0.0,0.0, 100.0])
    psspy.conl(0,1,2,[0,0],[ 100.0,0.0,0.0, 100.0])
    psspy.conl(0,1,3,[0,0],[ 100.0,0.0,0.0, 100.0])

    # Pre-run current-based load flow (TYSL, like SITER)
    psspy.ordr(0)
    psspy.fact()
    psspy.tysl(0)

    # Add dynamic file and set up conec, conet, and compile files
    psspy.dyre_new([1,1,1,1],str_pathdyrfile,str_pathconecfile,str_pathconetfile,str_pathcompilefile)

    # Frequency-dependent reactances and susceptances of branches
    psspy.set_netfrq(1)

    # Define solver parameters
    psspy.dynamics_solution_param_2([TYSLMAXITER,_i,_i,_i,_i,_i,_i,_i],[TYSLACCFAC,TYSLTOL,TSTEP,_f,_f,_f,_f,_f])
    
    # Choose pre-defined output channels
    psspy.change_channel_out_file(str_pathoutfile)
    psspy.chsb(0,1,[-1,-1,-1,1,1,0]) # angle
    psspy.chsb(0,1,[-1,-1,-1,1,2,0]) # pelec
    psspy.chsb(0,1,[-1,-1,-1,1,3,0]) # qelec
    psspy.chsb(0,1,[-1,-1,-1,1,4,0]) # eterm
    psspy.chsb(0,1,[-1,-1,-1,1,5,0]) # efd
    psspy.chsb(0,1,[-1,-1,-1,1,6,0]) # pmech
    psspy.chsb(0,1,[-1,-1,-1,1,7,0]) # speed

    # Add user-defined output channels (e.g., state variables of a user-defined model)
    ierr, I_STATE_VSC = psspy.mdlind(12,"""14""",'GEN','STATE')
    ierr, I_VAR_VSC = psspy.mdlind(12,"""14""",'GEN','VAR')
    ierr, I_STATE_DCGRID = psspy.mdlind(12,"""14""",'GOV','STATE') 
    # ierr, I_STATE_VSC = psspy.mdlind(9,"""14""",'GEN','STATE') 
    # ierr, I_STATE_DCGRID = psspy.mdlind(7,"""14""",'GOV','STATE') 
    psspy.state_channel([-1, I_STATE_VSC], r"xecd")
    psspy.state_channel([-1, I_STATE_VSC+1], r"xecq")
    # psspy.state_channel([-1, I_STATE_DCGRID], r"udc1")
    # psspy.state_channel([-1, I_STATE_DCGRID+1], r"udc2")
    psspy.var_channel([-1, I_VAR_VSC+22], r"udc")
    psspy.var_channel([-1, I_VAR_VSC+23], r"idc")
    psspy.var_channel([-1, I_VAR_VSC+25], r"pdc")

    # Add user-defined dynamic model
    psspy.addmodellibrary(str_pathdllfile)
    # os.system(str_pathcompilefile)

    # Running dynamic simulation
    # --------------------------
    # Initialization
    psspy.strt_2([0,1],str_pathoutfile)
    psspy.run(0, TINI,0,0,0)
    
    # Apply a disturbance
    # psspy.change_plmod_var(7, r"14", r"VSCGFL", 1, -3.5)
    psspy.change_gref(3, r"1", 0.65)

    # Simulate
    psspy.run(0, TFINAL,0,0,0)
    
    # Plotting
    # --------
    # Get indices
    ierr,NMACHINES = psspy.amachcount(-1,1)
    v_idxangle= range(0,NMACHINES)
    v_idxpelec= range(NMACHINES,2*NMACHINES)
    v_idxqelec= range(2*NMACHINES,3*NMACHINES)
    v_idxeterm = range(3*NMACHINES,4*NMACHINES)
    v_idxefd = range(4*NMACHINES,5*NMACHINES)
    v_idxpmech = range(5*NMACHINES,6*NMACHINES)
    v_idxspeed = range(6*NMACHINES,7*NMACHINES)
    v_idxvscstates = range(7*NMACHINES,7*NMACHINES+2) 
    # v_idxdcgridstates = range(7*NMACHINES+2,7*NMACHINES+4) 
    v_idxvscvars = range(7*NMACHINES+2,7*NMACHINES+5) 

    fontP = FontProperties()
    fontP.set_size('small')

    # Read channel file
    chnfobj = dyntools.CHNF(str_pathoutfile)
    short_title, chanid, chandata = chnfobj.get_data()

    # Extract time vector
    v_t = chandata['time']
    
    # Figure with 3 subplots: PELEC and QELEC and ETERM
    fig, axs = plt.subplots(3,1)
    
    for imach in v_idxpelec:
        axs[0].plot(v_t,chandata[imach+1],linewidth=2,label=chanid.values()[imach+1]) 
    axs[0].set_ylabel("Active power (pu)")
    axs[0].set_ylim(3,8)
    axs[0].legend(loc='upper right',fontsize=6,bbox_to_anchor=(1,1))

    for imach in v_idxqelec:
        axs[1].plot(v_t,chandata[imach+1],linewidth=2,label=chanid.values()[imach+1]) 
    axs[1].set_ylabel("Reactive power (pu)")
    axs[1].legend(loc='upper right',fontsize=6,bbox_to_anchor=(1,1))

    for imach in v_idxeterm:
        axs[2].plot(v_t,chandata[imach+1],linewidth=2,label=chanid.values()[imach+1]) 
    axs[2].set_ylabel("Terminal voltage (pu)")
    axs[2].set_xlabel("Time (s)")
    axs[2].legend(loc='upper right',fontsize=6,bbox_to_anchor=(1,1))

    
    plt.subplots_adjust(wspace=0.5)
    plt.savefig(os.path.join(str_pathsimfiles,"MachinesPQV.png"),bbox_inches="tight")
    plt.show()
    fig.clf()
    plt.close(fig)

    # Figure with 2 subplots: ANGLE and SPEED
    fig, axs = plt.subplots(2,1)
    
    for imach in v_idxspeed:
        axs[0].plot(v_t,chandata[imach+1],linewidth=2,label=chanid.values()[imach+1]) 
    axs[0].set_ylabel("Speed (pu)")
    axs[0].legend(loc='upper right',fontsize=6,bbox_to_anchor=(1,1))

    for imach in v_idxangle:
        axs[1].plot(v_t,chandata[imach+1],linewidth=2,label=chanid.values()[imach+1]) 
    axs[1].set_ylabel("Angle (pu)")
    axs[1].set_xlabel("Time (s)")
    axs[1].legend(loc='upper right',fontsize=6,bbox_to_anchor=(1,1))

    plt.subplots_adjust(wspace=0.5)
    plt.savefig(os.path.join(str_pathsimfiles,"MachinesDW.png"),bbox_inches="tight")
    plt.show()
    fig.clf()
    plt.close(fig)
    
    # Figure with two subplots: VSC states and currents
    fig, axs = plt.subplots(2,1) 
    for ix in v_idxvscstates:
        axs[0].plot(v_t,chandata[ix+1],linewidth=2,label=chanid.values()[ix+1]) 
    axs[0].set_ylabel("VSC states (pu)")
    # axs[0].set_ylim(-0.1,1.2)
    axs[0].legend(loc='upper right',fontsize=6,bbox_to_anchor=(1,1))

    for ix in v_idxvscvars:
        axs[1].plot(v_t,chandata[ix+1],linewidth=2,label=chanid.values()[ix+1]) 
    axs[1].set_ylabel("VSC variable (pu)")
    # axs[1].set_ylim(0.974,1.006)
    axs[1].set_xlabel("Time (s)")    
    axs[1].legend(loc='upper right',fontsize=6,bbox_to_anchor=(1,1))

    plt.savefig(os.path.join(str_pathsimfiles,"VSC.png"),bbox_inches="tight")
    plt.show()
    fig.clf()
    plt.close(fig)

