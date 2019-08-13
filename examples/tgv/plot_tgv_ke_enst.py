"""
        FILE: plot_tgv_ke_enst.py
      AUTHOR: Paul Bartholomew
 DESCRIPTION: Plots kinetic energy and enstrophy over time.
              Expects a file containing:
              < t > < KE > < IGNORED > < ENSTROPHY >
              which you provide by setting FILEPATH.
"""
import numpy as np
import matplotlib.pyplot as plt
plt.rc("text", usetex=True)
plt.rc("font", family="serif")
plt.rc("font", size=11)

FILEPATH="/media/yorgos/HardDrive2/TGV/iSVV/time_evol.dat"

def read_stats_first(filename):
  t = []
  enst = []
  ke = []

  with open(filename, "r") as data:
    print("Reading " + filename)
    for row in data:
      if not (row[0]=="#"):
        words = row.split()
        t.append(float(words[0]))
        enst.append(float(words[2]))
        ke.append(float(words[1]))

  #enst=-np.gradient(ke,t)
  return t, enst, ke

def read_stats_second(filename):
  t = []
  enst = []
  ke = []

  with open(filename, "r") as data:
    print("Reading " + filename)
    for row in data:
      if not (row[0]=="#"):
        words = row.split()
        t.append(float(words[0]))
        enst.append(float(words[4]))
        ke.append(float(words[1]))

  return t, enst, ke

def plot_stats(x3d_t, x3d_dat, x3d_lab, e3d_t, e3d_dat, e3d_lab,
               xlab, ylab, outfile, figsize=(5.0, 3.5)):

  plt.figure(figsize=figsize)

  plt.plot(x3d_t, x3d_dat, label=x3d_lab)
  plt.plot(e3d_t, e3d_dat, label=e3d_lab)

  plt.xlabel(xlab)
  plt.ylabel(ylab)
  plt.legend(prop={"family":"serif",
                   "size":11})

  plt.savefig(outfile, bbox_inches="tight")
  plt.close()

FILEPATH_iSVV="/media/yorgos/HardDrive2/TGV/iSVV/time_evol.dat"
iSVV_t, iSVV_enst, iSVV_ke = read_stats_second(FILEPATH_iSVV)
FILEPATH_SS="/media/yorgos/HardDrive2/TGV/SS/time_evol.dat"

SS_t, SS_enst, SS_ke = read_stats_first(FILEPATH_SS)
SS_enst2=-np.gradient(SS_ke,SS_t)

FILEPATH_Reference="/media/yorgos/HardDrive2/TGV/DNS_reference/TGV_Re10000.dat"
A=np.genfromtxt(FILEPATH_Reference,skip_header=47,delimiter='')
Ref_t=A[:,0]
Ref_ke=A[:,1]
Ref_enst=A[:,2]

fig = plt.figure(1,figsize=(10,4),edgecolor='none')
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.25, hspace=None)
ax1 = fig.add_subplot(121)
ax1.plot(Ref_t, Ref_enst, 'r',linewidth=2.0,label="DNS")
ax1.plot(iSVV_t, iSVV_enst, 'b', linewidth=2.0,label="iSVV")
ax1.plot(SS_t, SS_enst, 'g',linewidth=2.0, label="SS-SD")
ax1.plot(SS_t, SS_enst2, 'c',linewidth=2.0, label="SS-TD")
ax1.set_xlabel(r"$t$",fontsize=18)
ax1.set_ylabel(r"$\varepsilon_t$",fontsize=18)
ax1.set_yticks([0,0.005,0.01,0.015,0.02])
ax1.set_xticks([0,5,10,15,20])
ax1.tick_params(labelsize=16)
ax1.set_xlim(0,20)
ax1.set_ylim(0,0.022)
ax1.legend(prop={"family":"serif",
                 "size":14})
ax2 = fig.add_subplot(122)
ax2.plot(iSVV_t, iSVV_ke, 'r',linewidth=2.0 ,label="DNS")
ax2.plot(iSVV_t, iSVV_ke, 'b',linewidth=2.0,label="iSVV")
ax2.plot(SS_t, SS_ke, 'g',linewidth=2.0,label="SS")
ax2.set_xlabel(r"$t$",fontsize=18)
ax2.set_ylabel(r"$E_k$",fontsize=18)
ax2.set_yticks([0,0.02,0.04,0.06,0.08,0.1,0.12,0.14])
ax2.set_xticks([0,5,10,15,20])
ax2.tick_params(labelsize=16)
ax2.set_xlim(0,20)
ax2.set_ylim(0,0.14)
ax2.legend(prop={"family":"serif",
                 "size":14})

plt.savefig('TKE_Dissipation_LES_TGV.pdf',format='pdf',dpi=100,bbox_inches="tight")

plt.show()
plt.close()
