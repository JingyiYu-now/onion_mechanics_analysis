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
strainT = [5, 10, 15, 20, 25]


#Onion experiment with rehydration
Con1 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Sequential_Instron_With_Retracts_Tstrain_Data/40percent_PEG8k/Open cell onions/open_onion3/buffer'
PEG1 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Sequential_Instron_With_Retracts_Tstrain_Data/40percent_PEG8k/Open cell onions/open_onion3/PEG_2h'
Reh1 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Sequential_Instron_With_Retracts_Tstrain_Data/40percent_PEG8k/Open cell onions/open_onion3/PEG re-hydrated'

Con2 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Sequential_Instron_With_Retracts_Tstrain_Data/40percent_PEG8k/Open cell onions/open_onion4/Control'
PEG2 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Sequential_Instron_With_Retracts_Tstrain_Data/40percent_PEG8k/Open cell onions/open_onion4/PEG2h'
Reh2 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Sequential_Instron_With_Retracts_Tstrain_Data/40percent_PEG8k/Open cell onions/open_onion4/PEG_2h_rehydrated'

Con3 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Sequential_Instron_With_Retracts_Tstrain_Data/40percent_PEG8k/Open cell onions/open_onion5/Control'
PEG3 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Sequential_Instron_With_Retracts_Tstrain_Data/40percent_PEG8k/Open cell onions/open_onion5/PEG_2h'
Reh3 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Sequential_Instron_With_Retracts_Tstrain_Data/40percent_PEG8k/Open cell onions/open_onion5/PEG_2h_rehydration'

Ucon = [Con1, Con2, Con3]
UPEG = [PEG1, PEG2, PEG3]
UReH = [Reh1, Reh2, Reh3]
group = [Ucon, UPEG]

thickness = 7
width = 3
# linear fitting percantage (%)
fitpc = 5

color = ['green', 'orange', 'blue', 'black']


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
    # Xm = np.median(fitting_part.position)
    # compliance = 100 / (4500 + Xm) / z[0]
    modulus = z[0]
    return modulus

def fit_p(pull, p_pull, percentage):
    cutted_data = pull[pull.strain < np.max(p_pull.strain)].reset_index(drop=True)
    fitting_part = cutted_data[int(len(cutted_data) * (100 - percentage) / 100):].reset_index(drop=True)
    # fit using stress (MPa)
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

for unit in group:

    total = pd.DataFrame(columns= strainT)
    ela = pd.DataFrame(columns= strainT)
    retra = pd.DataFrame(columns= strainT)

    ## for data outpur name
    if g == 0:
        name = 'Con'
    elif g == 1:
        name = 'PEG'

    if g == 0:
        ls = '-'
    elif g == 1:
        ls = '--'
    elif g == 2:
        mark = 'p'
    elif g == 3:
        mark = 'P'
    for folder in unit:

        n = 0
        for file in sorted(glob.glob(os.path.join(folder, '*SHORT.csv'))):

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

            pull = [first_pull, second_pull, third_pull, fourth_pull, fifth_pull, sixth_pull]
            retract = [first_retract, second_retract, third_retract, fourth_retract, fifth_retract, sixth_retract]

            # calculate the stress and strain of the curve
            for i in range(len(pull)):
                norm(pull[i])
                norm(retract[i])

            # loop for seuqntial Instron
            for i in range(len(pull)):
                rmindex = []
                purify_revrs(retract[i])
                retract[i].drop(retract[i].loc[rmindex].index, inplace=True)

                retract[i] = retract[i].reset_index(drop=True)
                pull[i] = pull[i].reset_index(drop=True)

        # for sequential Instron - total compliance
            total_modulus = []
            for i in range(len(pull)-1):
                total_modulus.append(fit(pull[i], fitpc, 'N'))

            total.loc[len(total)] = total_modulus

        # for sequential Instron - elastic compliance
            ela_modulus = []

            ## fitting using the position
            for i in range(1, len(pull)):
                ela_modulus.append(fit_p(pull[i], pull[i - 1], fitpc))

            # for sequential Instron - retract compliance
            retra_modulus = []
            for i in range(len(pull) - 1):
                retra_modulus.append(fit(retract[i], fitpc, 'N'))

            # plot average onion
            retra.loc[len(retra)] = retra_modulus
            print(retra)

            # Line-plot with error bar
    ax = plt.subplot(111)
    x = strainT

    C_t = []
    C_e = []
    C_p = []
    C_r = []
    # C_h = []
    er_t = []
    er_e = []
    er_p = []
    er_r = []
    # er_h = []
    for i in range(len(strainT)):
        C_t.append(np.mean(total.iloc[:, i]))
        C_e.append(np.mean(ela.iloc[:, i]))
        C_r.append(np.mean(retra.iloc[:, i]))
        er_t.append(np.std(total.iloc[:, i]))
        er_e.append(np.std(ela.iloc[:, i]))
        er_r.append(np.std(retra.iloc[:, i]))

    total_chart = pd.DataFrame(columns=strainT)
    total_chart = ela
    total_chart.loc[len(total_chart)] = C_e
    total_chart.loc[len(total_chart)] = er_e
    # total_chart.to_csv(''f'{output}/' + name + '_Modulus_loading.csv')

    total_chart = pd.DataFrame(columns=strainT)
    total_chart = retra
    total_chart.loc[len(total_chart)] = C_r
    total_chart.loc[len(total_chart)] = er_r
    # total_chart.to_csv(''f'{output}/' + name + '_Modulus_unloading.csv')

    # plot with error bar
    ms = 6
    lw = 2
    # plt.errorbar(x, C_t, yerr = er_t, color= 'C7', linestyle = ls, marker = '8', markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)
    # plt.errorbar(x, C_t, yerr = er_t, color= color[g], linestyle = '-', marker = '8', markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)

    plt.errorbar(x, C_e, yerr = er_e, color='C2', linestyle = ls, marker = '8', markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)
    # plt.errorbar(x, C_p, yerr = er_p, color='C9', alpha = 0.6, marker = mark, markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)
    plt.errorbar(x, C_r, yerr = er_r, color='C1', linestyle = ls, marker = '8', markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)

    g += 1

x = [0, 5, 10, 15, 20, 25]
###### plot the total modulus for different treatments (Fig. 5B)
ax.set_xticks(x)
ax.set_yticks(range(0,700, 100))
ax.set_xticklabels(x)
ax.set(xlim = [0, 28], ylim = (0, 600), xlabel='Strain (%)', ylabel='Modulus (MPa)')
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().set_size_inches(6.4, 5.5)
green_line = mlines.Line2D([], [], color='green', label='Elastic loading modulus')
orange_line = mlines.Line2D([], [], color='orange', label='Elastic unloading modulus')
ax.legend(handles=[green_line, orange_line], loc=2, fontsize=10, frameon=False)

plt.show()



