"""
        FILE: plot_tgv_ke_enst.py
      AUTHOR: Paul Bartholomew
 DESCRIPTION: Plots kinetic energy and enstrophy over time.
              Expects a file containing:
              < t > < KE > < IGNORED > < ENSTROPHY >
              which you provide by setting FILEPATH.
"""

import matplotlib.pyplot as plt
plt.rc("text", usetex=True)
plt.rc("font", family="serif")
plt.rc("font", size=11)

FILEPATH="./ke_enst.dat"

def read_stats(filename):
  t = []
  enst = []
  ke = []

  with open(filename, "r") as data:
    print "Reading " + filename
    for row in data:
      if not (row[0]=="#"):
        words = row.split()
        t.append(float(words[0]))
        enst.append(float(words[3]))
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

x3d_t, x3d_enst, x3d_ke = read_stats(FILEPATH)

plt.figure(figsize=(5.0, 3.5))
plt.plot(x3d_t, x3d_enst, label="X3D")
plt.xlabel(r"$t$")
plt.ylabel(r"$\varepsilon$")
plt.legend(prop={"family":"serif",
                 "size":11})
plt.savefig("tgv_enstrophy.eps", bbox_inches="tight")
plt.close()

plt.figure(figsize=(5.0, 3.5))
plt.plot(x3d_t, x3d_ke, label="X3D")
plt.xlabel(r"$t$")
plt.ylabel(r"$k$")
plt.legend(prop={"family":"serif",
                 "size":11})
plt.savefig("tgv_ke.eps", bbox_inches="tight")
plt.close()
