from REFC_Lee import REFC
from new_dataframe import new
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib import rc, cm
import matplotlib.patches as patches

mi_name = "Samalas_MIs.csv"

df = pd.read_csv(mi_name)
df_melts = pd.read_csv("MagmaSat_model.csv")
# inputs
# dM_re = [0.00001, 0.00002, 0.00004, 0.00008]

rM_reM_x = [2, 1.5, 1.2, 1, 0.8]
dM_e = -0.0
dM_cc = 0
overturn = 3
c0_A = 1.3 # MgO
c0_i = 8 # Th
c0_ferric = 1.75 #ferric iron
c0_al = 16.5
c0_com = 1
c0_s = 500
kd_s = 75

L = len(rM_reM_x)
c = 4 #number of columns
r = 4 #number of rows
sp = 0
my_list_resutls = []
# re = REFC(N=N, dM_re=0.002, dM_cc=dM_cc, dM_x=-0.002, dM_e=dM_e)
# c_water, cfluid_x, mass_fluid = re.water_solubility(solubility=4.5, c_re=4, c_cc=3, c0=4)
# plt.figure(4)
# plt.plot(re.M_re, mass_fluid)
# plt.show()
i = 0
for j in range(0, L):
    if rM_reM_x[j] == 0.5:
        dM_re = 0.00001
        N = int(overturn/dM_re)
    elif rM_reM_x[j] == 0.25:
        dM_re = 0.000005
        N = int(overturn/dM_re)
    else:
        dM_re = 0.0001
        N = int(overturn/dM_re)

    new_d = new(N)
    results = new_d.df
    re = REFC(N=N, dM_re=dM_re, dM_cc=dM_cc, dM_x=-dM_re / rM_reM_x[j], dM_e=dM_e)
    Mt = re.M_ch
    ch_A = re.incompatible(c_re=c0_A, c_cc=0.2, D=2.2, c0=c0_A)  # compatible element, MgO
    ch_i = re.incompatible(c_re=c0_i, c_cc=10, D=0.002, c0=c0_i)  # incompatible element, Th
    ch_al = re.incompatible(c_re=c0_al, c_cc=20, D=1.1, c0=c0_al)  # compatible element, Al
    kd_ferric = (-48.378*(ch_A)+80.488)
    cs_mt = (ch_A*(-4)+7.5)/100
    print((kd_ferric))

    ch_ferric = re.ferric_iron(c_re=c0_ferric, c_cc=20, DFe3=kd_ferric, c0=c0_ferric, cs_mt=cs_mt)  # compatible element, Al
    print(len(ch_ferric[0]))
    # ch_com = re.incompatible(c_re=c0_com, c_cc=0.2, D=3, c0=c0_com)
    # ch_S, ch_S_x = re.solubility_control(solubility=1600, c_re=c0_s, c_cc=500, c0=c0_s)
    # ch_cu, ch_cu_x = re.chalcophile(Dsf=800, c_re=1, c_cc=0.4, c0=1, cs_x=ch_S_x)
    c_water, cfluid_x, mass_fluid = re.water_solubility(solubility=4.3, c_re=4.7, c_cc=3, c0=4.7)
    C_S, C_Sfluid, S_fluid = re.volatile_partition(kd=kd_s, c_re=c0_s, c_cc=100, c0=c0_s, cfluid_x=cfluid_x,
                                                   m_fluid=mass_fluid)
    mass_cumulate = re.M_x
    mass_recharge = re.M_re
    mass_chamber = re.M_ch
    results["M_re"] = mass_recharge
    results["M_x"] = mass_cumulate
    results["M_ch"] = mass_chamber
    results["MgO"] = ch_A
    results["Al2O3"] = ch_al
    results["Th"] = ch_i
    results["Fe3+"] = ch_ferric[0]
    results["S_m"] = C_S
    results["S_f"] = S_fluid
    results["fluid"] = mass_fluid
    filepath = Path(f'outputs/dM_re={dM_re}, dM_x = {dM_re / rM_reM_x[j]}.txt')
    filepath.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(filepath)
    my_list_resutls.append(results)

ferric = 0.5 * (71.844*2+15.99)*df["Fe3+/total Fe"]*df["FeO"]/(71.844*100)

rc('font',**{'size': 20})
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams["xtick.major.size"] = 4 # Sets length of ticks
plt.rcParams["ytick.major.size"] = 4 # Sets length of ticks
plt.rcParams["xtick.labelsize"] = 8 # Sets size of numbers on tick marks
plt.rcParams["ytick.labelsize"] = 8 # Sets size of numbers on tick marks
plt.rcParams["axes.titlesize"] = 10
plt.rcParams["axes.labelsize"] = 10 # Axes labels

sz_sm = 80
sz = 150
fig8, ax8 = plt.subplots(3, 2, figsize=(9.7, 6))
ax8 = ax8.flatten()
fig8.tight_layout()
plt.subplots_adjust(left=0.1, top=0.9, bottom=0.1)

ax8[0].plot(my_list_resutls[i * L]["M_re"], my_list_resutls[i * L]["Al2O3"])
ax8[0].plot(my_list_resutls[i * L + 1]["M_re"], my_list_resutls[i * L + 1]["Al2O3"])
ax8[0].plot(my_list_resutls[i * L + 2]["M_re"], my_list_resutls[i * L + 2]["Al2O3"])
ax8[0].plot(my_list_resutls[i * L + 3]["M_re"], my_list_resutls[i * L + 3]["Al2O3"])
ax8[0].plot(my_list_resutls[i * L + 4]["M_re"], my_list_resutls[i * L + 4]["Al2O3"])
ax8[0].set_ylabel("$\mathregular{Al_2O_3}$ (wt.%)")
ax8[0].set_xlim([0, 3])
ax8[0].set_ylim([14, 18])
ax8[0].annotate("a", xy=(-0.1, 1.05), xycoords="axes fraction", fontsize=12, weight="bold")

ax8[1].plot(my_list_resutls[i * L]["M_re"], my_list_resutls[i * L]["Th"])
ax8[1].plot(my_list_resutls[i * L + 1]["M_re"], my_list_resutls[i * L + 1]["Th"])
ax8[1].plot(my_list_resutls[i * L + 2]["M_re"], my_list_resutls[i * L + 2]["Th"])
ax8[1].plot(my_list_resutls[i * L + 3]["M_re"], my_list_resutls[i * L + 3]["Th"])
ax8[1].plot(my_list_resutls[i * L + 4]["M_re"], my_list_resutls[i * L + 4]["Th"])
ax8[1].set_ylabel("Th (ppm)")
ax8[1].set_xlim([0, 3])
ax8[1].set_ylim([0, 50])
ax8[1].annotate("b", xy=(-0.1, 1.05), xycoords="axes fraction", fontsize=12, weight="bold")

ax8[2].plot(my_list_resutls[i * L]["M_re"], my_list_resutls[i * L]["MgO"])
ax8[2].plot(my_list_resutls[i * L + 1]["M_re"], my_list_resutls[i * L + 1]["MgO"])
ax8[2].plot(my_list_resutls[i * L + 2]["M_re"], my_list_resutls[i * L + 2]["MgO"])
ax8[2].plot(my_list_resutls[i * L + 3]["M_re"], my_list_resutls[i * L + 3]["MgO"])
ax8[2].plot(my_list_resutls[i * L + 4]["M_re"], my_list_resutls[i * L + 4]["MgO"])
ax8[2].set_ylim([0.4, 1.4])
ax8[2].set_xlim([0, 3])
ax8[2].set_ylabel("MgO (wt.%)")
ax8[2].annotate("c", xy=(-0.1, 1.05), xycoords="axes fraction", fontsize=12, weight="bold")

ax8[3].plot(my_list_resutls[i * L]["M_re"], my_list_resutls[i * L]["S_m"])
ax8[3].plot(my_list_resutls[i * L + 1]["M_re"], my_list_resutls[i * L + 1]["S_m"])
ax8[3].plot(my_list_resutls[i * L + 2]["M_re"], my_list_resutls[i * L + 2]["S_m"])
ax8[3].plot(my_list_resutls[i * L + 3]["M_re"], my_list_resutls[i * L + 3]["S_m"])
ax8[3].plot(my_list_resutls[i * L + 4]["M_re"], my_list_resutls[i * L + 4]["S_m"])
ax8[3].set_ylabel("$\mathregular{S_{melt}}$ (ppm)")
ax8[3].set_xlabel("Magma reservoir overturns")
ax8[3].set_xlim([0, 3])
ax8[3].set_ylim([0, 500])
ax8[3].annotate("d", xy=(-0.1, 1.05), xycoords="axes fraction", fontsize=12, weight="bold")

ax8[4].plot(my_list_resutls[i * L]["M_re"], my_list_resutls[i * L]["fluid"])
ax8[4].plot(my_list_resutls[i * L + 1]["M_re"], -my_list_resutls[i * L + 1]["fluid"])
ax8[4].plot(my_list_resutls[i * L + 2]["M_re"], -my_list_resutls[i * L + 2]["fluid"])
ax8[4].plot(my_list_resutls[i * L + 3]["M_re"], -my_list_resutls[i * L + 3]["fluid"])
ax8[4].plot(my_list_resutls[i * L + 4]["M_re"], -my_list_resutls[i * L + 4]["fluid"])
ax8[4].set_ylim([-0.1, 0.2])
ax8[4].set_xlim([0, 3])
ax8[4].set_xlabel("Magma reservoir overturns")
ax8[4].set_ylabel("$\mathregular{M_{fluid}}$")
ax8[4].annotate("e", xy=(-0.1, 1.05), xycoords="axes fraction", fontsize=12, weight="bold")

ax8[5].plot(-my_list_resutls[i * L]["fluid"], my_list_resutls[i * L]["S_f"])
ax8[5].plot(-my_list_resutls[i * L + 1]["fluid"], my_list_resutls[i * L + 1]["S_f"])
ax8[5].plot(-my_list_resutls[i * L + 2]["fluid"], my_list_resutls[i * L + 2]["S_f"])
ax8[5].plot(-my_list_resutls[i * L + 3]["fluid"], my_list_resutls[i * L + 3]["S_f"])
ax8[5].plot(-my_list_resutls[i * L + 4]["fluid"], my_list_resutls[i * L + 4]["S_f"])
ax8[5].set_xlabel("$\mathregular{M_{fluid}}$")
ax8[5].set_ylabel("$\mathregular{S_{fluid}}$ (wt.%)")
ax8[5].set_xlim([0, 0.2])
ax8[5].set_ylim([1, 4])
ax8[5].legend(labels=[f"{rM_reM_x[0]}", f"{rM_reM_x[1]}", f"{rM_reM_x[2]}", f"{rM_reM_x[3]}", f"{rM_reM_x[4]}"],
              loc =(0.50,0.5), labelspacing = 0.2, handletextpad = 0.5,handlelength = 0.5, prop={'size': 10}, frameon=False, ncol=3)
ax8[5].annotate("Recharge/Crystallization", xy=(0.5, 0.85), xycoords="axes fraction", fontsize=12)
ax8[5].annotate("f", xy=(-0.1, 1.05), xycoords="axes fraction", fontsize=12, weight="bold")
rect = patches.Rectangle((0.095, 2.5), 0.1, 1.4, linewidth=0.8, edgecolor='black', facecolor='none')
ax8[5].add_patch(rect)

# plt.figure(2)
# plt.plot(my_list_resutls[i * L]["MgO"], my_list_resutls[i * L]["S_m"])
# plt.plot(my_list_resutls[i * L+1]["MgO"], my_list_resutls[i * L+1]["S_m"])
# plt.plot(my_list_resutls[i * L + 2]["MgO"], my_list_resutls[i * L + 2]["S_m"])
# plt.plot(my_list_resutls[i * L + 3]["MgO"], my_list_resutls[i * L + 3]["S_m"])
# plt.plot(my_list_resutls[i * L + 4]["MgO"], my_list_resutls[i * L + 4]["S_m"])
# plt.plot(df["MgO"][0:41], df["S"][0:41], "o", markerfacecolor="white", markeredgecolor="grey", markersize=6)
# plt.plot(df["MgO"][41:61], df["S"][41:61], "o", markerfacecolor="black", markeredgecolor="white", markersize=10)
# plt.plot(df["MgO"][61:65], df["S"][61:65], "d", markerfacecolor="black", markeredgecolor="black", markersize=10)
# plt.annotate("(a) S ", xy=(0.01, 0.90), xycoords="axes fraction", fontsize=12)
# # plt.xlabel("MgO (wt.%)")
# plt.ylabel("$\mathregular{S_{melt}}$ (ppm)")
# plt.xlim([0.5, 1.4])
# plt.ylim([0, 500])
# plt.legend(labels=[f"{rM_reM_x[0]}", f"{rM_reM_x[1]}", f"{rM_reM_x[2]}",
#             f"{rM_reM_x[3]}", f"{rM_reM_x[4]}"], loc =(0.800,0.05))
#
# plt.subplot(3, 1, 2)
# plt.plot(my_list_resutls[i * L]["MgO"], my_list_resutls[i * L]["Al2O3"])
# plt.plot(my_list_resutls[i * L+1]["MgO"], my_list_resutls[i * L+1]["Al2O3"])
# plt.plot(my_list_resutls[i * L + 2]["MgO"], my_list_resutls[i * L + 2]["Al2O3"])
# plt.plot(my_list_resutls[i * L + 3]["MgO"], my_list_resutls[i * L + 3]["Al2O3"])
# plt.plot(my_list_resutls[i * L + 4]["MgO"], my_list_resutls[i * L + 4]["Al2O3"])
# plt.plot(df["MgO"][0:41], df["Al2O3"][0:41], "o", markerfacecolor="white", markeredgecolor="grey", markersize=6)
# plt.plot(df["MgO"][41:61], df["Al2O3"][41:61], "o", markerfacecolor="black", markeredgecolor="white", markersize=10)
# plt.plot(df["MgO"][61:65], df["Al2O3"][61:65], "d", markerfacecolor="black", markeredgecolor="black", markersize=10)
# # plt.plot(df_melts["MgO"][2:], df_melts["Al2O3"][2:])
# # plt.xlabel("MgO (wt.%)")
# plt.ylabel("$\mathregular{Al_2O_3}$ (wt.%)")
# plt.xlim([0.5, 1.4])
# plt.ylim([15, 17])
# plt.annotate("(b) $\mathregular{Al_2O_3}$", xy=(0.01, 0.90), xycoords="axes fraction", fontsize=12)
# plt.subplot(3, 1, 3)
# plt.plot(my_list_resutls[i * L]["MgO"], my_list_resutls[i * L]["Th"])
# plt.plot(my_list_resutls[i * L + 1]["MgO"], my_list_resutls[i * L + 1]["Th"])
# plt.plot(my_list_resutls[i * L+2]["MgO"], my_list_resutls[i * L+2]["Th"])
# plt.plot(my_list_resutls[i * L + 3]["MgO"], my_list_resutls[i * L + 3]["Th"])
# plt.plot(my_list_resutls[i * L + 4]["MgO"], my_list_resutls[i * L + 4]["Th"])
# plt.plot(df["MgO"][0:40], df["Th"][0:40], "o", markerfacecolor="white", markeredgecolor="grey", markersize=6)
# # plt.plot(df["MgO"][41:60], df["Th"][41:60], "o", markerfacecolor="black", markeredgecolor="white", markersize=10)
# plt.xlabel("MgO (wt.%)")
# plt.ylabel("Th (ppm)")
# plt.xlim([0.5, 1.4])
# plt.ylim([6, 20])
# plt.annotate("(c) Th", xy=(0.01, 0.90), xycoords="axes fraction", fontsize=12)
# plt.subplot(2, 2, 4)
# plt.plot(my_list_resutls[i * L]["MgO"], my_list_resutls[i * L]["Fe3+"])
# plt.plot(my_list_resutls[i * L + 1]["MgO"], my_list_resutls[i * L + 1]["Fe3+"])
# plt.plot(my_list_resutls[i * L+2]["MgO"], my_list_resutls[i * L+2]["Fe3+"])
# plt.plot(my_list_resutls[i * L + 3]["MgO"], my_list_resutls[i * L + 3]["Fe3+"])
# plt.plot(my_list_resutls[i * L + 4]["MgO"], my_list_resutls[i * L + 4]["Fe3+"])
# plt.plot(df_melts["MgO"][2:], df_melts["Fe2O3"][2:])
# plt.plot(df["MgO"][41:61], ferric[41:61], "o", markerfacecolor="black", markeredgecolor="white", markersize=10)
# plt.plot(df["MgO"][61:65], ferric[61:65], "d", markerfacecolor="black", markeredgecolor="black", markersize=10)
# plt.xlabel("MgO (wt.%)")
# plt.ylabel("$\mathregular{Fe_2O_3}$ (wt.%)")
# plt.xlim([0.5, 1.4])
# plt.ylim([0.5, 2.5])

plt.figure(3, figsize=(6, 8))
plt.subplot(3,1,1)
plt.plot(my_list_resutls[i * L]["M_re"], my_list_resutls[i * L]["S_m"])
plt.plot(my_list_resutls[i * L + 1]["M_re"], my_list_resutls[i * L + 1]["S_m"])
plt.plot(my_list_resutls[i * L + 2]["M_re"], my_list_resutls[i * L + 2]["S_m"])
plt.plot(my_list_resutls[i * L + 3]["M_re"], my_list_resutls[i * L + 3]["S_m"])
plt.plot(my_list_resutls[i * L + 4]["M_re"], my_list_resutls[i * L + 4]["S_m"])
# plt.plot(df["M_re"], df["S"], "o")
plt.ylabel("$\mathregular{S_melt}$ (ppm)")
plt.xlim([0, 3])
plt.ylim([0, 500])
plt.legend([f"{rM_reM_x[0]}", f"{rM_reM_x[1]}", f"{rM_reM_x[2]}",
            f"{rM_reM_x[3]}", f"{rM_reM_x[4]}"])
plt.subplot(3, 1, 2)
plt.plot(my_list_resutls[i * L]["M_re"], my_list_resutls[i * L]["fluid"])
plt.plot(my_list_resutls[i * L + 1]["M_re"], -my_list_resutls[i * L + 1]["fluid"])
plt.plot(my_list_resutls[i * L + 2]["M_re"], -my_list_resutls[i * L + 2]["fluid"])
plt.plot(my_list_resutls[i * L + 3]["M_re"], -my_list_resutls[i * L + 3]["fluid"])
plt.plot(my_list_resutls[i * L + 4]["M_re"], -my_list_resutls[i * L + 4]["fluid"])
plt.ylim([-0.1, 0.2])
plt.xlim([0, 3])
plt.xlabel("Overtun")
plt.ylabel("$\mathregular{M_fluid}$")

plt.subplot(3,1,3)
plt.plot(-my_list_resutls[i * L]["fluid"], my_list_resutls[i * L]["S_f"])
plt.plot(-my_list_resutls[i * L + 1]["fluid"], my_list_resutls[i * L + 1]["S_f"])
plt.plot(-my_list_resutls[i * L + 2]["fluid"], my_list_resutls[i * L + 2]["S_f"])
plt.plot(-my_list_resutls[i * L + 3]["fluid"], my_list_resutls[i * L + 3]["S_f"])
plt.plot(-my_list_resutls[i * L + 4]["fluid"], my_list_resutls[i * L + 4]["S_f"])
plt.xlabel("$\mathregular{M_fluid}$")
plt.ylabel("$\mathregular{S_fluid}$ (wt.%)")
plt.xlim([0, 0.2])
plt.ylim([1, 4])

plt.show()








# plt.figure(1)
# plt.plot(ch_A, ch_i)
# plt.plot(ch_A, ch_com)
#
# plt.figure(2)
# plt.plot(ch_A, ch_S)
# plt.plot(ch_A, ch_S_x)
#
# plt.figure(3)
# plt.plot(mass_cumulate, ch_cu)
# plt.plot(mass_cumulate, ch_cu_x)
#

#
# plt.figure(5)
# plt.plot(mass_fluid, S_fluid)
# # plt.plot(mass_fluid, C_Sfluid)
#
# plt.figure(6)
# plt.subplot(2,2,1)
# plt.plot(ch_A, C_S)
# plt.plot(df["MgO"], df["S"], "o")
#
# plt.subplot(2,2,2)
# plt.plot(ch_i, C_S)
# plt.plot(df["Th"], df["S"], "o")
#
# plt.subplot(2,2,3)
# plt.plot(ch_A, ch_al)
# plt.plot(df["MgO"], df["Al2O3"], "o")
# plt.show()
