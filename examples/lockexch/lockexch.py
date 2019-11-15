"""
"""

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.rc("text", usetex=True)
plt.rc("font", family="serif")
plt.rc("font", size=11)
plt.rc("axes", axisbelow=True)

import matplotlib.pyplot as plt

#####################################################
## Read, normalise and plot energy budget
def main():

    ## Read energies and dissipation rates
    # t, ek, dek, ep, dep = read_budget("/home/paul/DATA/Incompact3d/lock-exchange/non-boussinesq-r07/out/budget")
    t, ek, dek, ep, dep = read_budget("/home/paul/DATA/Incompact3d/lock-exchange/boussinesq-uset02/out/budget")

    # ## Subtract background potential energy
    # L = 18.0
    # H = 1.0
    # D = 2.0
    # rho0 = 0.998
    # Fr2 = 1 - rho0
    # Fr = Fr2**0.5
    # print min(ep), max(ep)
    # for i in range(len(ep)):
    #     ep[i] -= rho0 * L * (2 * H)**2 * D / (Fr**2) / 2
    # # print min(ep), max(ep)

    ## Integrate dissipation rates
    et = integrate(dek, t)
    em = integrate(dep, t)
    for i in range(len(em)):
        em[i] *= -1

    ## Compute total energy budget (should be 1 for all time)
    etot = calc_etot(ek, et, ep, em)

    ## Normalise wrt initial total energy budget
    etot0 = etot[0]
    ek = normalise(ek, etot0)
    et = normalise(et, etot0)
    ep = normalise(ep, etot0)
    em = normalise(em, etot0)
    etot = normalise(etot, etot0)

    ref_data = read_digitised_data("/home/paul/DATA/benchmarking/lockexch/necker2002-ebudg.csv")

    ## Plot and save
    plt.plot(t, etot, color="black", label=r"$E_{tot}$")
    plt.plot(t, ek, color="red", label=r"$K$")
    plt.plot(t, ep, color="green", label=r"$E_p$")
    plt.plot(t, et, color="blue", label=r"$E_d$")
    plt.plot(t, em, color="purple", label=r"$E_{\nabla^2 \rho}$")
    #plt.axhline(0, ls=":", color="black")
    #plt.axhline(1, ls=":", color="black")
    plt.grid()

    # Ref data
    # plt.plot(ref_data["Etot"][0], ref_data["Etot"][1], color="black",
    #          ls="--")
    plt.plot(ref_data["KE"][0], ref_data["KE"][1], color="red",
             ls="--")
    plt.plot(ref_data["Ep"][0], ref_data["Ep"][1], color="green",
             ls="--")
    plt.plot(ref_data["Ed"][0], ref_data["Ed"][1], color="blue",
             ls="--")
    plt.plot(ref_data["Es"][0], ref_data["Es"][1], color="purple",
             ls="--")

    plt.xlim(xmin=0, xmax=30)
    plt.ylim((-0.2, 1.2))

    plt.xlabel(r"$t$")
    plt.ylabel(r"$E / E^0_{tot}$")

    plt.legend()
    plt.savefig("lockexch-budget.eps", bbox_inches="tight")

def read_budget(budgfile):

    t = []
    ek = []
    dek = []
    ep = []
    dep = []

    with open(budgfile, "r") as data:
        #next(data) # skip header
        for row in data:
            words = row.split()

            t.append(float(words[0]))
            ek.append(float(words[1]))
            dek.append(float(words[2]))
            ep.append(float(words[3]))
            dep.append(float(words[4]))

    return t, ek, dek, ep, dep

def normalise(e, eref):

    enorm = []
    for i in range(len(e)):
        enorm.append(e[i] / eref)

    return enorm

def integrate(diss, t):

    e = [0]

    for i in range(1, len(diss)):

        dedt = 0.5 * (diss[i] + diss[i - 1])
        dt = t[i] - t[i - 1]
        e.append(e[i - 1] + dedt * dt)

    return e

def calc_etot(ek, dek, ep, dep):

    etot = []
    for i in range(len(ek)):
        etot.append(ek[i] + dek[i] + ep[i] + dep[i])

    return etot

def read_digitised_data(datafile, sep=","):

    datadict = {}
    with open(datafile, "r") as data:
        header = data.readline()
        words = header.split(sep)
        col = 0
        coldict = {}
        for key in words:
            if (len(key) and (key != "\n")):
                coldict[key] = col
                datadict[key] = [[], []]
                col += 1

        data.next() # Skip "x" "y" row

        for row in data:
            words = row.split(sep)
            print words

            for key in datadict.keys():
                idx = 2 * coldict[key]
                if (len(words[idx + 0])):
                    print words[idx + 1], words[idx + 1]
                    datadict[key][0].append(float(words[idx + 0]))
                    datadict[key][1].append(float(words[idx + 1]))

    return datadict

#####################################################
## Run as a script
if __name__ == "__main__":
    main()
