# File:"C:\Users\lsigrist\OneDrive - Universidad Pontificia Comillas\PSSE\Tools\HY_ACDC_SIM2\Version2023\Dynamic\TestRun.py", generated on TUE, JUN 11 2024  11:48, PSS(R)E release 34.06.01
import psse34
import psspy, redirect
redirect.psse2py()
from psspy import _i, _f, _s, _o
import dyntools
from matplotlib import pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np

psspy.psseinit(2000)

psspy.case(r"""C:\Users\lsigrist\OneDrive - Universidad Pontificia Comillas\PSSE\Tools\HY_ACDC_SIM2\Version2023\Input\KundurACDC\kundur32_AC_MTDCg.sav""")
psspy.fnsl([0,0,0,1,0,0,99,0])
psspy.cong(0)
psspy.conl(0,1,1,[0,0],[ 100.0,0.0,0.0, 100.0])
psspy.conl(0,1,2,[0,0],[ 100.0,0.0,0.0, 100.0])
psspy.conl(0,1,3,[0,0],[ 100.0,0.0,0.0, 100.0])
psspy.ordr(0)
psspy.fact()
psspy.tysl(0)
psspy.dyre_new([1,1,1,1],r"""C:\Users\lsigrist\OneDrive - Universidad Pontificia Comillas\PSSE\Tools\HY_ACDC_SIM2\Version2023\Input\KundurACDC\kundur_MTDCg.dyr""",
r"""C:\Users\lsigrist\OneDrive - Universidad Pontificia Comillas\PSSE\Tools\HY_ACDC_SIM2\Version2023\Simulation\KundurACDC\conec""",
r"""C:\Users\lsigrist\OneDrive - Universidad Pontificia Comillas\PSSE\Tools\HY_ACDC_SIM2\Version2023\Simulation\KundurACDC\conet""",
r"""C:\Users\lsigrist\OneDrive - Universidad Pontificia Comillas\PSSE\Tools\HY_ACDC_SIM2\Version2023\Simulation\KundurACDC\compile""")
psspy.set_netfrq(1)
psspy.dynamics_solution_param_2([_i,_i,_i,_i,_i,_i,_i,_i],[_f,_f, 0.001,_f,_f,_f,_f,_f])
psspy.change_channel_out_file(r"""C:\Users\lsigrist\OneDrive - Universidad Pontificia Comillas\PSSE\Tools\HY_ACDC_SIM2\Version2023\Simulation\KundurACDC\output""")
psspy.chsb(0,1,[-1,-1,-1,1,1,0])
psspy.chsb(0,1,[-1,-1,-1,1,2,0])
psspy.chsb(0,1,[-1,-1,-1,1,3,0])
psspy.chsb(0,1,[-1,-1,-1,1,4,0])
psspy.chsb(0,1,[-1,-1,-1,1,5,0])
psspy.chsb(0,1,[-1,-1,-1,1,6,0])
psspy.chsb(0,1,[-1,-1,-1,1,7,0])

ierr, ival = psspy.mdlind(5,"""14""",'GEN','VAR') 
I_VAR_SVSCON_1 = ival
ierr, ival = psspy.mdlind(5,"""14""",'GEN','GOV') 
I_VAR_DCGRID_1 = ival

psspy.state_channel([-1, 92], r"VDC1_DC1")
psspy.state_channel([-1, 93], r"VDC2_DC1")
# psspy.state_channel([-1, 145], r"VDC3_DC1")
# psspy.state_channel([-1, 149], r"VDC1_DC2")
# psspy.state_channel([-1, 150], r"VDC2_DC2")
# psspy.state_channel([-1, 151], r"VDC3_DC2")
psspy.var_channel([-1, 145], r"PDC1_DC1")
psspy.var_channel([-1, 146], r"PDC2_DC1")
# psspy.var_channel([-1, 412], r"PDC3_DC1")
# psspy.var_channel([-1, 416], r"PDC1_DC2")
# psspy.var_channel([-1, 417], r"PDC2_DC2")
# psspy.var_channel([-1, 418], r"PDC3_DC2")
psspy.var_channel([-1, 126], r"ADDPREF_DC1_7")
psspy.var_channel([-1, 106], r"MEANDELTA")

psspy.addmodellibrary(r"""C:\Users\lsigrist\OneDrive - Universidad Pontificia Comillas\PSSE\Tools\HY_ACDC_SIM2\Version2023\Simulation\KundurACDC\MTDC.dll""")
psspy.strt_2([0,1],r"""C:\Users\lsigrist\OneDrive - Universidad Pontificia Comillas\PSSE\Tools\HY_ACDC_SIM2\Version2023\Simulation\KundurACDC\output.out""")
psspy.run(0, 1,0,0,0)
psspy.change_plmod_var(7, r"14", r"SVSCON", 1, 0.8)

psspy.run(0, 10,0,0,0)

v_idxdeltavsc = range(4,6)
v_idxpemach = range(6,10)
v_idxpevsc = range(10,12)
v_idxqevsc = range(16,18)
v_idxvdcmtdc1 = range(42,44)
v_idxvdcmtdc2 = range(42,44)
v_idxpdcmtdc1 = range(44,46)
v_idxpdcmtdc2 = range(44,46)
v_idxaddprefdc17 = 46
v_idxmeandelta = 47

fontP = FontProperties()
fontP.set_size('small')

chnfobj = dyntools.CHNF(r"C:\Users\lsigrist\OneDrive - Universidad Pontificia Comillas\PSSE\Tools\HY_ACDC_SIM2\Version2023\Simulation\KundurACDC\output.out")
short_title, chanid, chandata = chnfobj.get_data()
v_t = chandata['time']

# equivalent
fig, axs = plt.subplots(3,1)
l_legend1 = []
for imach in v_idxpemach:
    v_pe = chandata[imach+1]
    axs[0].plot(v_t,v_pe,linewidth=2)
    l_legend1.append(chanid.values()[imach+1]) 
axs[0].set_ylabel("Active power (pu)")
#axs[0].set_ylim(-0.1,0.1)
l_legend2 = []
for imach in v_idxpevsc:
    v_etm = chandata[imach+1]
    axs[1].plot(v_t,v_etm,linewidth=2)
    l_legend2.append(chanid.values()[imach+1]) 
axs[1].set_ylabel("Active power (pu)")
# axs[1].set_ylim(0,1.5)
l_legend3 = []
for imach in v_idxqevsc:
    v_etm = chandata[imach+1]
    axs[2].plot(v_t,v_etm,linewidth=2)
    l_legend3.append(chanid.values()[imach+1]) 
axs[2].set_ylabel("Reactive power (pu)")
# axs[1].set_ylim(0,1.5)
axs[0].legend(l_legend1,loc='upper right',fontsize=6,bbox_to_anchor=(1,1))
axs[1].legend(l_legend2,loc='upper right',fontsize=6,bbox_to_anchor=(1,1))
plt.subplots_adjust(wspace=0.5)
plt.savefig("Machines.png",bbox_inches="tight")
fig.clf()
plt.close(fig)

fig, axs = plt.subplots(2,1) 
for imach in v_idxdeltavsc:
    v_vdc = chandata[imach+1]
    axs[0].plot(v_t,v_vdc,linewidth=2)
    l_legend1.append(chanid.values()[imach+1]) 
axs[0].plot(v_t,chandata[v_idxmeandelta+1],linewidth=2)
axs[0].set_ylabel("Delta (pu)")
v_vdc = np.array(chandata[v_idxdeltavsc[0]+1]) - np.array(chandata[v_idxdeltavsc[1]+1])
axs[1].plot(v_t,v_vdc,linewidth=2)
axs[1].plot(v_t,chandata[v_idxaddprefdc17+1],linewidth=2)
axs[1].set_ylabel("DiffDelta (pu)")
plt.savefig("MTDCAngle.png",bbox_inches="tight")
fig.clf()
plt.close(fig)

fig, axs = plt.subplots(2,1) 
for imach in v_idxvdcmtdc1:
    v_vdc = chandata[imach+1]
    axs[0].plot(v_t,v_vdc,linewidth=2)
    l_legend1.append(chanid.values()[imach+1]) 
axs[0].set_ylabel("Vdc (pu)")
for imach in v_idxvdcmtdc2:
    v_vdc = chandata[imach+1]
    axs[1].plot(v_t,v_vdc,linewidth=2)
    l_legend1.append(chanid.values()[imach+1]) 
axs[1].set_ylabel("Vdc (pu)")
plt.savefig("MTDCVdc.png",bbox_inches="tight")
fig.clf()
plt.close(fig)

fig, axs = plt.subplots(2,1) 
for imach in v_idxpdcmtdc1:
    v_vdc = chandata[imach+1]
    axs[0].plot(v_t,v_vdc,linewidth=2)
    l_legend1.append(chanid.values()[imach+1]) 
axs[0].set_ylabel("Pdc (pu)")
for imach in v_idxpdcmtdc2:
    v_vdc = chandata[imach+1]
    axs[1].plot(v_t,v_vdc,linewidth=2)
    l_legend1.append(chanid.values()[imach+1]) 
axs[1].set_ylabel("Pdc (pu)")
plt.savefig("MTDCPdc.png",bbox_inches="tight")
fig.clf()
plt.close(fig)