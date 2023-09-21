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
from numpy import trapz
import statsmodels.api as sm
strainT = [5, 10, 15, 20, 25]

# User input
# folder = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion2/10mmpermin_9.22.20'
# folder2 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion2/5mmpermin_9.22.20'
# folder3 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion2/5mmpermin_9.17.20'
# folder4 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion1'
# Onion3 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion3/5mmpermin'
# Onion4 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion4'
# Onion5 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion5'
# Onion6 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion6'
# # Con = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/Onion_PEG/For_code/Con_pH7.0'
# # PEG = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/Onion_PEG/For_code/40percentPEG8k_pH7.0'
# output = '/Users/jingyiyu/Documents/Cosgrovelab/manuscript/CW_mechanics_Instron/Data_collection/Python_plot/For new Fig4'


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


output = '/Users/jingyiyu/Documents/Cosgrovelab/manuscript/CW_mechanics_Instron/Data_collection/Python_plot/PEG'
# group = [ folder2, Onion3, Onion5]
# group = [Con4, PEG4]
# group = [ folder2]
Ucon = [Con1, Con2, Con3]
UPEG = [PEG1, PEG2, PEG3]
UReH = [Reh1, Reh2, Reh3]
group = [Ucon, UPEG]
group2 = [Ucon, UPEG, UReH]
# total = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_seq_instron/4,8,12,16,20_sequence/Total_3onion'
# folder4 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_seq_instron/4,8,12,16,20_sequence/seq_Onion_10'
# folder5 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_seq_instron/4,8,12,16,20_sequence/Onion_6.23.20_9'
# group = [folder1, folder2, folder3, folder4, folder5]

# 8,12,16,20,24 sequence
seq1 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_seq_instron/8,12,16,20,24_sequence/Onion1'
seq2 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_seq_instron/8,12,16,20,24_sequence/onion2_10'

name = 'Onion4'
# file = ''f'{folder}/Sequential_Instron__02_15_2020__14_30_58_SHORT.csv'
# output = ''f'{folder1}/summary1.csv'

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
    # strain calculated based on length at the beginning of each pull
    # ori_p = ten_p(pull)

    pull['force_N'] = pull.load * 0.0098

    # engineering stress & strain
    pull['stress'] = pull.force_N / (thickness * width * 0.001)
    pull['strain'] = (pull.position - ori_p) / 5000
    # true stress & strain
    # pull['strain'] = np.log(1 + (pull.position - ori_p) / (ori_p + 4500) )
    # pull['stress'] = pull.force_N/(thickness * width * (ori_p + 4500) / (pull.position + 4500) * 0.001)

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
        # fit using gram as unit
        # z = np.polyfit(fitting_part.position, fitting_part.load, 1)

        # fit using stress (MPa)
        z = np.polyfit(fitting_part.strain, fitting_part.stress, 1)
        pull['fit_e'] = pull.strain * z[0] + z[1]
    else:
        cutted_data = pull[:ten_ind(pull, target_load)].reset_index(drop=True)
        fitting_part = cutted_data[int(len(cutted_data) * (100 - percentage) / 100):].reset_index(drop=True)
        # fit using gram as unit
        # z = np.polyfit(fitting_part.position, fitting_part.load, 1)

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
    # fit using gram as unit
    # z = np.polyfit(fitting_part.position, fitting_part.load, 1)

    # fit using stress (MPa)
    z = np.polyfit(fitting_part.strain, fitting_part.stress, 1)
    pull['fit_e'] = pull.strain * z[0] + z[1]

    # Xm = np.median(fitting_part.position)
    # modulus = 100 / (4500 + Xm) / z[0]
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
# total = pd.DataFrame(columns=loadx)
# ela = pd.DataFrame(columns=loadx)
# retra = pd.DataFrame(columns=loadx)
# pla = pd.DataFrame(columns=loadx)
# t_modulus = pd.DataFrame(columns=loadx)
# e_modulus = pd.DataFrame(columns=loadx)
# r_modulus = pd.DataFrame(columns=loadx)
for unit in group2:

    total = pd.DataFrame(columns= strainT)
    ela = pd.DataFrame(columns= strainT)
    retra = pd.DataFrame(columns= strainT)
    pla = pd.DataFrame(columns= strainT)

    if g == 0:
        ls = '-'
    elif g == 1:
        ls = '--'
    elif g == 2:
        mark = 'p'
    elif g == 3:
        mark = 'P'
    for folder in unit:

        # create data frame for results from each property
        # t_modulus = pd.DataFrame(columns=loadx)
        # e_modulus = pd.DataFrame(columns=loadx)
        # r_modulus = pd.DataFrame(columns=loadx)
        # t_inc_comp = pd.DataFrame(columns=loadx)
        # e_inc_comp = pd.DataFrame(columns=loadx)
        # r_inc_comp = pd.DataFrame(columns=loadx)



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
                # retract[i].drop(retract[i].loc[purify_p(retract[i], target[i]): len(retract[i])].index, inplace=True)
                # retract[i].drop(retract[i].loc[0:ten_ind(retract[i], 0.1)].index, inplace=True)
                # pull[i].drop(pull[i].loc[0:ten_ind(pull[i], 0.1)].index, inplace=True)

                # retract[i].reset_index(inplace=True)
                rmindex = []
                purify_revrs(retract[i])
                retract[i].drop(retract[i].loc[rmindex].index, inplace=True)

                retract[i] = retract[i].reset_index(drop=True)
                pull[i] = pull[i].reset_index(drop=True)

        # for sequential Instron - total compliance
            total_modulus = []
            for i in range(len(pull)-1):
                total_modulus.append(fit(pull[i], fitpc, 'N'))

        # add data of this file to the results data frame
            # plot the individual onion
        #     t_modulus.loc[len(t_modulus)] = total_modulus
            # plot average onion
            total.loc[len(total)] = total_modulus

        # incremental total compliance
        #     inc_t = [fit(first_pull, fitpc, 'N')]
        #     for i in range(len(total_comp)):
        #         if i > 0:
        #             inc_t.append(total_comp[i] - total_comp[i - 1])
        #     t_inc_comp.loc[len(t_inc_comp)] = inc_t

        # to save total fit in a different column
        #     for i in range(len(pull)-1):
        #         pull[i]['fit'] = pull[i].fit_e


        # for sequential Instron - elastic compliance
            ela_modulus = []
            # for i in range(1, len(pull)-1):
            #     ela_comp.append(fit(pull[i], fitpc, loadx[i-1]))
            # ela_comp.append(fit(pull[-1], fitpc, 'N'))

            ## fitting using the position
            for i in range(1, len(pull)):
                ela_modulus.append(fit_p(pull[i], pull[i - 1], fitpc))

            # plot individual onion
            # e_modulus.loc[len(e_modulus)] = ela_modulus

            # plot average onion
            ela.loc[len(ela)] = ela_modulus

            # incremental elastic compliance
            # inc_e = [ela_comp[0]]
            # for i in range(len(ela_comp)):
            #     if i > 0:
            #         inc_e.append(ela_comp[i] - ela_comp[i - 1])
            # e_inc_comp.loc[len(e_inc_comp)] = inc_e


        # for sequential Instron - plastic compliance
        #     pla_comp = []
        #     for i in range(len(pull)-2):
        #         pla_comp.append(fit(pull[i], fitpc, 'N') - fit(pull[i+1], fitpc, loadx[0]))
        #     pla_comp.append(fit(pull[-2], fitpc, 'N') - fit(pull[-1], fitpc, 'N'))
        # # pla_forplot = [ -x for x in pla_comp]
        #     p_comp.loc[len(p_comp)] = pla_comp

            # pla_modulus = []
            # for i in range(len(loadx)):
            #     pla_modulus.append(1/(1/total_modulus[i] - 1/ela_modulus[i]))
            # pla_forplot = [ -x for x in pla_comp]
            # pla.loc[len(pla)] = pla_modulus

        # incremental plastic compliance
        #     inc_p = [fit(first_pull, fitpc, 'N') - fit(second_pull, fitpc, loadx[0])]
        #     for i in range(len(pla_comp)):
        #         if i > 0:
        #             inc_p.append(pla_comp[i] - pla_comp[i - 1])
        #     p_inc_comp.loc[len(p_inc_comp)] = inc_p

            # for sequential Instron - retract compliance
            retra_modulus = []
            for i in range(len(pull) - 1):
                retra_modulus.append(fit(retract[i], fitpc, 'N'))
            # pla_forplot = [ -x for x in pla_comp]
            #plot individual onion
            # r_modulus.loc[len(r_modulus)] = retra_modulus

            # plot average onion
            retra.loc[len(retra)] = retra_modulus
            print(retra)
            # fit(retract[-1], fitpc, 'N')

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
        # C_p.append(np.mean(pla.iloc[:, i]))
        er_t.append(np.std(total.iloc[:, i]))
        er_e.append(np.std(ela.iloc[:, i]))
        er_r.append(np.std(retra.iloc[:, i]))
        # er_p.append(np.std(pla.iloc[:, i]))

    print(retra)
    # er_t = [np.std(t_comp[loadx[0]]), np.std(t_comp[loadx[1]]), np.std(t_comp[loadx[2]]), np.std(t_comp[loadx[3]]), np.std(t_comp[loadx[4]])]
    # er_e = [np.std(e_comp[loadx[0]]), np.std(e_comp[loadx[1]]), np.std(e_comp[loadx[2]]), np.std(e_comp[loadx[3]]), np.std(e_comp[loadx[4]])]
    # er_p = [np.std(p_comp[loadx[0]]), np.std(p_comp[loadx[1]]), np.std(p_comp[loadx[2]]), np.std(p_comp[loadx[3]]), np.std(p_comp[loadx[4]])]

    # plot with error bar
    ms = 6
    lw = 2
    # plt.errorbar(x, C_t, yerr = er_t, color= 'C7', linestyle = ls, marker = '8', markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)
    plt.errorbar(x, C_t, yerr = er_t, color= color[g], linestyle = '-', marker = '8', markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)

    # plt.errorbar(x, C_e, yerr = er_e, color='C2', linestyle = ls, marker = '8', markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)
    # plt.errorbar(x, C_p, yerr = er_p, color='C9', alpha = 0.6, marker = mark, markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)
    # plt.errorbar(x, C_r, yerr = er_r, color='C1', linestyle = ls, marker = '8', markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)

    g += 1


#### plot the elastic loading and unloading curves (Fig. 5C)
# ax.set_xticks(strainT)
# ax.set_yticks(range(0,650, 100))
# ax.set_xticklabels(strainT)
# ax.set(ylim = (0, 650), xlabel='Strain (%)', ylabel='Modulus (MPa)')
# # t_line = mlines.Line2D([], [], color = 'C7', label = 'Total loading modulus')
# # green_line = mlines.Line2D([], [], color = 'C9', label = 'Plastic modulus')
# # orange_line = mlines.Line2D([], [], color = 'C3', label = 'Hysteresis compliance', linestyle = '--')
# fr_line = mlines.Line2D([], [], color = 'C2', label = 'Elastic loading modulus')
# blue_line = mlines.Line2D([], [], color = 'C1', label = 'Elastic unloading modulus')
# # ax.legend(handles = [fr_line, blue_line], loc = 2, fontsize=10, frameon = False)
# ax.legend(handles = [fr_line, blue_line], loc = 2, fontsize=10, frameon = False)
# plt.gcf().subplots_adjust(bottom=0.15)
# plt.gcf().set_size_inches(6.4, 5.5)
# ax.grid(alpha=0.4, linestyle='--', linewidth = 2)
x = [0, 5, 10, 15, 20, 25]
###### plot the total modulus for different treatments (Fig. 5B)
ax.set_xticks(x)
# ax.set_yticks(range(0,700, 100))
ax.set_xticklabels(x)
ax.set(xlim = [0, 28], ylim = (0, 100), xlabel='Strain (%)', ylabel='Modulus (MPa)')
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().set_size_inches(6.4, 5.5)
green_line = mlines.Line2D([], [], color='green', label='Elastic loading modulus')
orange_line = mlines.Line2D([], [], color='orange', label='Elastic unloading modulus')
# blue_line = mlines.Line2D([], [], color='blue', label='40% PEG8k rehydrated')
# ax.legend(handles=[green_line, orange_line], loc=2, fontsize=10, frameon=False)
# ax.legend(handles=[green_line, orange_line], loc=2, fontsize=10, frameon=False)


plt.savefig(''f'{output}/Total_Modulus.pdf', transparent = True)
plt.show()



