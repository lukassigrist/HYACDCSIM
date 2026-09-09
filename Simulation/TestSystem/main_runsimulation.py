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
str_lffile = r"kundur32_noAC_MTDCg.sav" # initial AC load flow"
str_dyrfile = r"kundur_MTDCg.dyr" # initial AC load flow"
str_dllfile = r"USRSOU.dll"
str_pathfiles = r"C:\Users\lsigrist\OneDrive - Universidad Pontificia Comillas\PSSE\Tools\HYACDCSIM\Simulation\TestSystem" # Path

# Solver parameters
TSTEP = 0.002
TYSLACCFAC = 0.8
TYSLTOL = 0.0001
TYSLMAXITER = 125
TINI = 1
TFINAL = 20

# ------------
# Preparations
# ------------
str_pathlffile = os.path.join(str_pathfiles,str_lffile)
str_pathdyrfile = os.path.join(str_pathfiles,str_dyrfile)
str_pathdllfile = os.path.join(str_pathfiles,str_dllfile)
str_pathoutfile = os.path.join(str_pathfiles,r"output.out")
str_pathconecfile = os.path.join(str_pathfiles,r"conec")
str_pathconetfile = os.path.join(str_pathfiles,r"conet")
str_pathcompilefile = os.path.join(str_pathfiles,r"compile")

# Guaranty that working directory is where Fortran code is located, which is where the input files (.txt) are located
print("Current Working Directory:", os.getcwd())
os.chdir(str_pathfiles)
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
    psspy.conl(0,1,1,[0,0],[ 100.0,00.0,0.0, 100.0])
    psspy.conl(0,1,2,[0,0],[ 100.0,00.0,0.0, 100.0])
    psspy.conl(0,1,3,[0,0],[ 100.0,00.0,0.0, 100.0])

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
    # ierr, I_STATE_VSC = psspy.mdlind(4,"""1""",'GEN','STATE')
    # ierr, I_VAR_VSC = psspy.mdlind(4,"""1""",'GEN','VAR')
    ierr, I_STATE_VSC = psspy.mdlind(9,"""14""",'GEN','STATE') 
    ierr, I_VAR_VSC = psspy.mdlind(9,"""14""",'GEN','VAR')
    ierr, I_STATE_DCGRID = psspy.mdlind(7,"""14""",'GOV','STATE') 
    # VSCGFO
    '''psspy.state_channel([-1, I_STATE_VSC], r"xecd")
    psspy.state_channel([-1, I_STATE_VSC+1], r"xecq")
    psspy.state_channel([-1, I_STATE_VSC+2], r"delta")
    psspy.state_channel([-1, I_STATE_VSC+3], r"omega")
    # psspy.state_channel([-1, I_STATE_DCGRID], r"udc1")
    # psspy.state_channel([-1, I_STATE_DCGRID+1], r"udc2")
    psspy.var_channel([-1, I_VAR_VSC+13], r"id")
    psspy.var_channel([-1, I_VAR_VSC+14], r"iq")
    psspy.var_channel([-1, I_VAR_VSC+15], r"deltap")
    psspy.var_channel([-1, I_VAR_VSC+16], r"pdc")'''
    # VSCODRO
    psspy.state_channel([-1, I_STATE_VSC], r"xec")
    psspy.state_channel([-1, I_STATE_VSC+1], r"delta")
    psspy.state_channel([-1, I_STATE_DCGRID], r"udc1")
    psspy.state_channel([-1, I_STATE_DCGRID+1], r"udc2")
    psspy.var_channel([-1, I_VAR_VSC+11], r"id")
    psspy.var_channel([-1, I_VAR_VSC+12], r"iq")
    psspy.var_channel([-1, I_VAR_VSC+13], r"deltap")
    psspy.var_channel([-1, I_VAR_VSC+14], r"pc+ploss")

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
    # psspy.change_plmod_var(7, r"14", r"SVSCON", 1, -3.5)
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
    # v_idxvscstates = range(7*NMACHINES+2,7*NMACHINES+4) # VSCGFO
    v_idxvscstates = range(7*NMACHINES,7*NMACHINES+2) # VSCDRO
    v_idxdcgridstates = range(7*NMACHINES+2,7*NMACHINES+4) # DCGRID
    # v_idxvscvarsi = range(7*NMACHINES+4,7*NMACHINES+6) # VSCGFO
    # v_idxvscvarsp = range(7*NMACHINES+7,7*NMACHINES+8) # VSCGFO
    v_idxvscvarsi = range(7*NMACHINES+4,7*NMACHINES+6) # VSCDRO
    v_idxvscvarsp = range(7*NMACHINES+7,7*NMACHINES+8) # VSCDRO

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
    plt.savefig(os.path.join(str_pathfiles,"MachinesPQV.png"),bbox_inches="tight")
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
    plt.savefig(os.path.join(str_pathfiles,"MachinesDW.png"),bbox_inches="tight")
    plt.show()
    fig.clf()
    plt.close(fig)
    

    # Figure with two subplots: VSC states and currents
    fig, axs = plt.subplots(3,1,sharex=True) 
    '''for ix in v_idxvscstates:
        v_x = chandata[ix+1]
        axs[0].plot(v_t,v_x,linewidth=2)
        l_legend1.append(chanid.values()[ix+1])'''
    axs[0].plot(v_t,chandata[v_idxvscstates[0]+1],linewidth=2,color='blue',label=chanid.values()[v_idxvscstates[0]+1])
    ax2 = axs[0].twinx()
    ax2.plot(v_t,chandata[v_idxvscstates[1]+1],linewidth=2,color='orange',label=chanid.values()[v_idxvscstates[1]+1])
    axs[0].set_ylabel("VSC states (pu)")
    # axs[0].set_ylim(-0.1,1.2)
    handles1, labels1 = axs[0].get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    axs[0].legend(handles1 + handles2, labels1 + labels2, loc="upper right",fontsize=6,bbox_to_anchor=(1,1),frameon=True)

    '''for ix in v_idxvscvarsi:
        v_i = chandata[ix+1]
        axs[1].plot(v_t,v_i,linewidth=2)
        l_legend2.append(chanid.values()[ix+1])'''
    axs[1].plot(v_t,chandata[v_idxvscvarsi[0]+1],linewidth=2,color='blue',label=chanid.values()[v_idxvscvarsi[0]+1])
    ax2 = axs[1].twinx()
    ax2.plot(v_t,-np.array(chandata[v_idxvscvarsi[1]+1]),linewidth=2,color='orange',label="-"+chanid.values()[v_idxvscvarsi[1]+1])
    v_i = np.sqrt(np.square(chandata[v_idxvscvarsi[0]+1])+np.square(chandata[v_idxvscvarsi[1]+1]))
    axs[1].plot(v_t,v_i,linewidth=2,color='green',label="Iabs")
    axs[1].set_ylabel("VSC currents (pu)")
    # axs[1].set_ylim(0.974,1.006)
    axs[1].set_xlabel("Time (s)")  
    handles1, labels1 = axs[1].get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()  
    axs[1].legend(handles1 + handles2, labels1 + labels2, loc="upper right",fontsize=6,bbox_to_anchor=(1,1),frameon=True)

    for ix in v_idxvscvarsp:
        v_p = chandata[ix+1]
        axs[2].plot(v_t,v_p,linewidth=2,label=chanid.values()[ix+1])
    axs[2].plot(v_t,chandata[v_idxpelec[3]+1],linewidth=2,label=chanid.values()[v_idxpelec[3]+1])
    axs[2].plot(v_t,np.array(chandata[v_idxqelec[3]+1]),linewidth=2,label=chanid.values()[v_idxqelec[3]+1])
    axs[2].set_ylabel("Powers (pu)")
    # axs[1].set_ylim(0.974,1.006)
    axs[2].set_xlabel("Time (s)")    
    axs[2].legend(loc='upper right',fontsize=6,bbox_to_anchor=(1,1))

    

    plt.savefig(os.path.join(str_pathfiles,"VSC.png"),bbox_inches="tight")
    plt.show()
    fig.clf()
    plt.close(fig)

    # Figure with two plots
    l_legend4 = []
    fig, axs = plt.subplots(2,1) 
    for ix in v_idxvscstates:
        v_x = chandata[ix+1]
        axs[0].plot(v_t,v_x,linewidth=2,label=chanid.values()[ix+1]) 
    axs[0].set_ylabel("VSC states (pu)")
    # axs[0].set_ylim(-0.5,4.5)

    l_legend5 = []
    for ix in v_idxdcgridstates:
        v_udc = chandata[ix+1]
        axs[1].plot(v_t,v_udc,linewidth=2,label=chanid.values()[ix+1])
    axs[1].set_ylabel("DC grid states (pu)")
    # axs[1].set_ylim(0.974,1.006)
    axs[1].set_xlabel("Time (s)")

    axs[0].legend(loc='upper right',fontsize=6,bbox_to_anchor=(1,1))
    axs[1].legend(loc='upper right',fontsize=6,bbox_to_anchor=(1,1))

    plt.savefig(os.path.join(str_pathfiles,"States.png"),bbox_inches="tight")
    plt.show()
    fig.clf()
    plt.close(fig)
