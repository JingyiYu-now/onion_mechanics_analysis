# script: Study progression of elasticity and plasticity with compliance method
# output value of T & E & P compliance and their incremental compliance;
# plot: line-plot for exact value, stack bar chart for incremental value

# import necessary package
import pandas as pd
import numpy as np
import os
import glob
from matplotlib import pyplot as plt
import matplotlib.lines as mlines
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Avenir'
plt.rcParams['font.size'] = 15
plt.rcParams['axes.linewidth'] = 2
mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['ytick.major.size'] = 5
mpl.rcParams['ytick.major.width'] = 2
mpl.use('macosx')

# User input
folder2 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion2/5mmpermin_9.22.20'
Onion3 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion3/5mmpermin'
Onion5 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion5'
Onion6 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion6'

group = [folder2, Onion3, Onion5, Onion6]

thickness = 7
width = 3
# linear fitting percantage (%)
fitpc = 5

# for load sequence of 4-8-12-16-20
loadx = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20]

# define a function that calls the index of certain value
def ind(df, value):
    for i in range(len(df)):
        if df.position[i] == value:
            return i


# define a function that find the index at certain loads
def ten_ind(pull, load):
    for i in range(len(pull) - 5):
        if (pull.load[i] > load) & (pull.load[i + 1] > load) & (pull.load[i + 2] > load) & \
                (pull.load[i + 3] > load) & (pull.load[i + 4] > load):
            index = i
            return index


def pos_p(pull, p_pull):
    for i in range(len(pull) - 5):
        n = np.max(p_pull.strain)
        if (pull.strain[i] > n) & (pull.strain[i + 1] > n) & (pull.strain[i + 2] > n) & \
                (pull.strain[i + 3] > n) & (pull.strain[i + 4] > n):
            index = i
            return index
        elif i == (len(pull) - 6):
            return i

# define a function that find the position where peels are under tension ( constant > 0.1g )
def ten_p(pull):
    for i in range(len(pull) - 5):
        if (pull.load[i] > 0.1) & (pull.load[i + 1] > 0.1) & (pull.load[i + 2] > 0.1) & \
                (pull.load[i + 3] > 0.1) & (pull.load[i + 4] > 0.1):
            point = pull.position[i]
            return point


# define a function that calculate the stress and strain of the curve
def norm(pull):
    # strain calculated based on very original length (before plastic deformation)
    ori_p = ten_p(first_pull)

    pull['force_N'] = pull.load * 0.0098

    # engineering stress & strain
    pull['stress'] = pull.force_N / (thickness * width * 0.001)
    pull['strain'] = (pull.position - ori_p) / 5000

def purify_p(pull, target):
    for i in range(len(pull)):
        if i < (len(pull)-1) :
            if pull.load[i] > target:
                return i
                break
        else:
            return len(pull) -1

def purify_revrs(pull):
    b = np.array(pull.position)
    maxp_index = np.argmax(b)
    threshload = pull.load[maxp_index]
    for i in range(len(pull)):
        if pull.load[i] > threshload:
            rmindex.append(i)

# extract the fitting part of the data and calculate the compliance according to the fitting part
# (option: fit using load/stress)
def fit(pull, percentage, target_load):
    if target_load == 'N':
        fitting_part = pull[int(len(pull) * (100 - percentage) / 100):].reset_index(drop=True)
        # fit using stress (MPa)
        z = np.polyfit(fitting_part.strain, fitting_part.stress, 1)
        pull['fit_e'] = pull.strain * z[0] + z[1]
    else:
        cutted_data = pull[:ten_ind(pull, target_load)].reset_index(drop=True)
        fitting_part = cutted_data[int(len(cutted_data) * (100 - percentage) / 100):].reset_index(drop=True)
        # fit using stress (MPa)
        z = np.polyfit(fitting_part.strain, fitting_part.stress, 1)
        pull['fit_e'] = pull.strain * z[0] + z[1]
    modulus = z[0]
    return modulus

def fit_p(pull, p_pull, percentage):
    cutted_data = pull[pull.strain < np.max(p_pull.strain)].reset_index(drop=True)
    fitting_part = cutted_data[int(len(cutted_data) * (100 - percentage) / 100):].reset_index(drop=True)
    z = np.polyfit(fitting_part.strain, fitting_part.stress, 1)
    pull['fit_e'] = pull.strain * z[0] + z[1]

    modulus = z[0]
    return modulus


def tar_ind(pull, load):
    if load == 'N':
        return len(pull)
    else:
        stress = load * 0.0098/ (thickness * width * 0.001)
        for i in range(len(pull) - 5):
            if (pull.stress[i] > stress) & (pull.stress[i + 1] > stress) & (pull.stress[i + 2] > stress) & \
                    (pull.stress[i + 3] > stress) & (pull.stress[i + 4] > stress):
                return i

# extract data from each pull
# Use if multiple files need analysis
g = 0
mark = '8'
total = pd.DataFrame(columns=loadx)
ela = pd.DataFrame(columns=loadx)
retra = pd.DataFrame(columns=loadx)
pla = pd.DataFrame(columns=loadx)

for folder in group:


    n = 0
    for file in sorted(glob.glob(os.path.join(folder, '*SHORT.csv'))):
        if n % 2 == 0:
            print(n)
            whole_df = pd.read_table(file,
                                     delimiter=',',
                                     header=None,
                                     names=range(0, 2),
                                     )
            whole_df.columns = ['position', 'load']

            first_pull = whole_df[ind(whole_df, 'FIRST PULL') + 1:ind(whole_df, 'FIRST RETRACT')].astype(
                float).reset_index(drop=True)
            first_retract = whole_df[ind(whole_df, 'FIRST RETRACT') + 1:ind(whole_df, 'SECOND PULL')].astype(
                float).iloc[::-1].reset_index(drop=True)
            second_pull = whole_df[ind(whole_df, 'SECOND PULL') + 1:ind(whole_df, 'SECOND RETRACT')].astype(
                float).reset_index(drop=True)
            second_retract = whole_df[ind(whole_df, 'SECOND RETRACT') + 1:ind(whole_df, 'THIRD PULL')].astype(
                float).iloc[::-1].reset_index(drop=True)
            third_pull = whole_df[ind(whole_df, 'THIRD PULL') + 1:ind(whole_df, 'THIRD RETRACT')].astype(
                float).reset_index(drop=True)
            third_retract = whole_df[ind(whole_df, 'THIRD RETRACT') + 1:ind(whole_df, 'FOURTH PULL')].astype(
                float).iloc[::-1].reset_index(drop=True)
            fourth_pull = whole_df[ind(whole_df, 'FOURTH PULL') + 1:ind(whole_df, 'FOURTH RETRACT')].astype(
                float).reset_index(drop=True)
            fourth_retract = whole_df[ind(whole_df, 'FOURTH RETRACT') + 1:ind(whole_df, 'FIFTH PULL')].astype(
                float).iloc[::-1].reset_index(drop=True)
            fifth_pull = whole_df[ind(whole_df, 'FIFTH PULL') + 1:ind(whole_df, 'FIFTH RETRACT')].astype(
                float).reset_index(drop=True)
            fifth_retract = whole_df[ind(whole_df, 'FIFTH RETRACT') + 1:ind(whole_df, 'SIXTH PULL')].astype(float).iloc[
                            ::-1].reset_index(drop=True)
            sixth_pull = whole_df[ind(whole_df, 'SIXTH PULL') + 1:ind(whole_df, 'SIXTH RETRACT')].astype(
                float).reset_index(drop=True)
            sixth_retract = whole_df[ind(whole_df, 'SIXTH RETRACT') + 1:].astype(float).iloc[::-1].reset_index(
                drop=True)

            ori_p = ten_p(first_pull)
            n += 1
            continue

        if n % 2 == 1:
            print(n)
            whole_df = pd.read_table(file,
                                     delimiter=',',
                                     header=None,
                                     names=range(0, 2),
                                     )
            whole_df.columns = ['position', 'load']

            seventh_pull = whole_df[ind(whole_df, 'FIRST PULL') + 1:ind(whole_df, 'FIRST RETRACT')].astype(
                float).reset_index(drop=True)
            seventh_retract = whole_df[ind(whole_df, 'FIRST RETRACT') + 1:ind(whole_df, 'SECOND PULL')].astype(
                float).iloc[::-1].reset_index(drop=True)
            eighth_pull = whole_df[ind(whole_df, 'SECOND PULL') + 1:ind(whole_df, 'SECOND RETRACT')].astype(
                float).reset_index(drop=True)
            eighth_retract = whole_df[ind(whole_df, 'SECOND RETRACT') + 1:ind(whole_df, 'THIRD PULL')].astype(
                float).iloc[::-1].reset_index(drop=True)
            ninth_pull = whole_df[ind(whole_df, 'THIRD PULL') + 1:ind(whole_df, 'THIRD RETRACT')].astype(
                float).reset_index(drop=True)
            ninth_retract = whole_df[ind(whole_df, 'THIRD RETRACT') + 1:ind(whole_df, 'FOURTH PULL')].astype(
                float).iloc[::-1].reset_index(drop=True)
            tenth_pull = whole_df[ind(whole_df, 'FOURTH PULL') + 1:ind(whole_df, 'FOURTH RETRACT')].astype(
                float).reset_index(drop=True)
            tenth_retract = whole_df[ind(whole_df, 'FOURTH RETRACT') + 1:ind(whole_df, 'FIFTH PULL')].astype(
                float).iloc[::-1].reset_index(drop=True)
            eleventh_pull = whole_df[ind(whole_df, 'FIFTH PULL') + 1:ind(whole_df, 'FIFTH RETRACT')].astype(
                float).reset_index(drop=True)
            eleventh_retract = whole_df[ind(whole_df, 'FIFTH RETRACT') + 1:ind(whole_df, 'SIXTH PULL')].astype(
                float).iloc[::-1].reset_index(drop=True)
            twlveth_pull = whole_df[ind(whole_df, 'SIXTH PULL') + 1:ind(whole_df, 'SIXTH RETRACT')].astype(
                float).reset_index(drop=True)
            twlveth_retract = whole_df[ind(whole_df, 'SIXTH RETRACT') + 1:].astype(float).iloc[::-1].reset_index(
                drop=True)

            n += 1

            pull = [first_pull, second_pull, third_pull, fourth_pull, fifth_pull, seventh_pull, eighth_pull,
                    ninth_pull, tenth_pull, eleventh_pull, twlveth_pull]
            retract = [first_retract, second_retract, third_retract, fourth_retract, fifth_retract,
                       seventh_retract, eighth_retract, ninth_retract, tenth_retract, eleventh_retract, twlveth_retract]

        # calculate the stress and strain of the curve
        for i in range(len(pull)):
            norm(pull[i])
            norm(retract[i])

        target = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 20]
        # loop for seuqntial Instron
        for i in range(len(pull)):
            retract[i].drop(retract[i].loc[purify_p(retract[i], target[i]): len(retract[i])].index, inplace=True)

            retract[i].reset_index(inplace=True)
            rmindex = []
            purify_revrs(retract[i])
            retract[i].drop(retract[i].loc[rmindex].index, inplace=True)

            retract[i] = retract[i].reset_index(drop=True)
            pull[i] = pull[i].reset_index(drop=True)

    # for sequential Instron - total compliance
        total_modulus = []
        for i in range(len(pull)-1):
            total_modulus.append(fit(pull[i], fitpc, 'N'))

    #     t_modulus.loc[len(t_modulus)] = total_modulus
        # plot average onion
        total.loc[len(total)] = total_modulus

    # for sequential Instron - elastic compliance
        ela_modulus = []

        ## fitting using the position
        for i in range(1, len(pull)):
            ela_modulus.append(fit_p(pull[i], pull[i - 1], fitpc))

        # plot average onion
        ela.loc[len(ela)] = ela_modulus

        pla_modulus = []
        for i in range(len(loadx)):
            pla_modulus.append(1/(1/total_modulus[i] - 1/ela_modulus[i]))
        # pla_forplot = [ -x for x in pla_comp]
        pla.loc[len(pla)] = pla_modulus

        # for sequential Instron - retract compliance
        retra_modulus = []
        for i in range(len(pull) - 1):
            retra_modulus.append(fit(retract[i], fitpc, 'N'))

        # plot average onion
        retra.loc[len(retra)] = retra_modulus
        print(retra)

    ####### plot the stress-train curve and fitting curve
        ela_target = ['N', 2, 4, 6, 8, 10, 12, 14, 16, 18, 'N']

        # Line-plot with error bar
ax = plt.subplot(111)
#plot as stress
x = [a * 0.098 /3/7/0.01 for a in loadx]

C_t = []
C_e = []
C_p = []
C_r = []

er_t = []
er_e = []
er_p = []
er_r = []
for i in range(len(loadx)):
    C_t.append(np.mean(total.iloc[:, i]))
    C_e.append(np.mean(ela.iloc[:, i]))
    C_r.append(np.mean(retra.iloc[:, i]))
    er_t.append(np.std(total.iloc[:, i]))
    er_e.append(np.std(ela.iloc[:, i]))
    er_r.append(np.std(retra.iloc[:, i]))

print(retra)
# plot with error bar
ms = 4
lw = 2
plt.errorbar(x, C_t, yerr = er_t, color= 'C7', alpha = 0.6, marker = mark, markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)
plt.errorbar(x, C_e, yerr = er_e, color='C2', alpha = 0.6, marker = mark, markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)
plt.errorbar(x, C_r, yerr = er_r, color='C1', alpha = 0.6, marker = mark, markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)

xlabel = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
ax.set_xticks(xlabel)
ax.set_yticks(range(0,400, 50))
ax.set_xticklabels(xlabel)
ax.set(ylim = [-50,400], xlabel='Stress (MPa)', ylabel='Modulus (MPa)')
t_line = mlines.Line2D([], [], color = 'C7', label = 'Total loading modulus')
fr_line = mlines.Line2D([], [], color = 'C2', label = 'Elastic loading modulus')
blue_line = mlines.Line2D([], [], color = 'C1', label = 'Elastic unloading modulus')
ax.legend(handles = [t_line, fr_line, blue_line], loc = 2, fontsize=10, frameon = False)
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().set_size_inches(6.4, 5.5)

# plt.savefig(''f'{output}/Strain_modulus_seq.pdf', transparent = True)
plt.show()
