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

# User input
folder = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion2/10mmpermin_9.22.20'
folder2 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion2/5mmpermin_9.22.20'
folder3 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion2/5mmpermin_9.17.20'
folder4 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion1'
Onion3 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion3/5mmpermin'
Onion4 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion4'
Onion5 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion5'
Onion6 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion6'
# Con = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/Onion_PEG/For_code/Con_pH7.0'
# PEG = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/Onion_PEG/For_code/40percentPEG8k_pH7.0'
output = '/Users/jingyiyu/Documents/Cosgrovelab/manuscript/CW_mechanics_Instron/Data_collection/Python_plot/For new Fig4'


group = [folder2, Onion3, Onion5, Onion6]


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

# for load sequence of 4-8-12-16-20
loadx = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
# for load sequence of 8-12-16-20-24
# loadx = [8, 12, 16, 20, 24]


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
total = pd.DataFrame(columns=loadx)
ela = pd.DataFrame(columns=loadx)
retra = pd.DataFrame(columns=loadx)
pla = pd.DataFrame(columns=loadx)
# t_modulus = pd.DataFrame(columns=loadx)
# e_modulus = pd.DataFrame(columns=loadx)
# r_modulus = pd.DataFrame(columns=loadx)
ax = plt.subplot(111)
for folder in group:

    # create data frame for results from each property
    # t_modulus = pd.DataFrame(columns=loadx)
    # e_modulus = pd.DataFrame(columns=loadx)
    # r_modulus = pd.DataFrame(columns=loadx)
    # t_inc_comp = pd.DataFrame(columns=loadx)
    # e_inc_comp = pd.DataFrame(columns=loadx)
    # r_inc_comp = pd.DataFrame(columns=loadx)

    # if g == 0:
    #     mark = '8'
    # elif g == 1:
    #     mark = 's'
    # elif g == 2:
    #     mark = 'p'
    # elif g == 3:
    #     mark = 'P'

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

        # if n % 2 == 1:
        #     print(n)
        #     whole_df = pd.read_table(file,
        #                              delimiter=',',
        #                              header=None,
        #                              names=range(0, 2),
        #                              )
        #     whole_df.columns = ['position', 'load']
        #
        #     seventh_pull = whole_df[ind(whole_df, 'FIRST PULL') + 1:ind(whole_df, 'FIRST RETRACT')].astype(
        #         float).reset_index(drop=True)
        #     seventh_retract = whole_df[ind(whole_df, 'FIRST RETRACT') + 1:ind(whole_df, 'SECOND PULL')].astype(
        #         float).iloc[::-1].reset_index(drop=True)
        #     eighth_pull = whole_df[ind(whole_df, 'SECOND PULL') + 1:ind(whole_df, 'SECOND RETRACT')].astype(
        #         float).reset_index(drop=True)
        #     eighth_retract = whole_df[ind(whole_df, 'SECOND RETRACT') + 1:ind(whole_df, 'THIRD PULL')].astype(
        #         float).iloc[::-1].reset_index(drop=True)
        #     ninth_pull = whole_df[ind(whole_df, 'THIRD PULL') + 1:ind(whole_df, 'THIRD RETRACT')].astype(
        #         float).reset_index(drop=True)
        #     ninth_retract = whole_df[ind(whole_df, 'THIRD RETRACT') + 1:ind(whole_df, 'FOURTH PULL')].astype(
        #         float).iloc[::-1].reset_index(drop=True)
        #     tenth_pull = whole_df[ind(whole_df, 'FOURTH PULL') + 1:ind(whole_df, 'FOURTH RETRACT')].astype(
        #         float).reset_index(drop=True)
        #     tenth_retract = whole_df[ind(whole_df, 'FOURTH RETRACT') + 1:ind(whole_df, 'FIFTH PULL')].astype(
        #         float).iloc[::-1].reset_index(drop=True)
        #     eleventh_pull = whole_df[ind(whole_df, 'FIFTH PULL') + 1:ind(whole_df, 'FIFTH RETRACT')].astype(
        #         float).reset_index(drop=True)
        #     eleventh_retract = whole_df[ind(whole_df, 'FIFTH RETRACT') + 1:ind(whole_df, 'SIXTH PULL')].astype(
        #         float).iloc[::-1].reset_index(drop=True)
        #     twlveth_pull = whole_df[ind(whole_df, 'SIXTH PULL') + 1:ind(whole_df, 'SIXTH RETRACT')].astype(
        #         float).reset_index(drop=True)
        #     twlveth_retract = whole_df[ind(whole_df, 'SIXTH RETRACT') + 1:].astype(float).iloc[::-1].reset_index(
        #         drop=True)


        pull = [first_pull, second_pull, third_pull, fourth_pull, fifth_pull]
        retract = [first_retract, second_retract, third_retract, fourth_retract, fifth_retract]

        # calculate the stress and strain of the curve
        for i in range(len(pull)):
            norm(pull[i])
            norm(retract[i])

        # target = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 20]
        # loop for seuqntial Instron
        for i in range(len(pull)):
            retract[i].drop(retract[i].loc[purify_p(retract[i], target[i]): len(retract[i])].index, inplace=True)
            # retract[i].drop(retract[i].loc[0:ten_ind(retract[i], 0.1)].index, inplace=True)
            # pull[i].drop(pull[i].loc[0:ten_ind(pull[i], 0.1)].index, inplace=True)

            retract[i].reset_index(inplace=True)
            rmindex = []
            purify_revrs(retract[i])
            retract[i].drop(retract[i].loc[rmindex].index, inplace=True)

            retract[i] = retract[i].reset_index(drop=True)
            pull[i] = pull[i].reset_index(drop=True)



# xlabel = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
# ax.set_xticks(xlabel)
# ax.set_yticks(range(0,400, 50))
# ax.set_xticklabels(xlabel)
# ax.set(ylim = [-50,400], xlabel='Stress (MPa)', ylabel='Modulus (MPa)')
# t_line = mlines.Line2D([], [], color = 'C7', label = 'Total loading modulus')
# # green_line = mlines.Line2D([], [], color = 'C9', label = 'Plastic modulus')
# # orange_line = mlines.Line2D([], [], color = 'C3', label = 'Hysteresis compliance', linestyle = '--')
# fr_line = mlines.Line2D([], [], color = 'C2', label = 'Elastic loading modulus')
# blue_line = mlines.Line2D([], [], color = 'C1', label = 'Elastic unloading modulus')
# ax.legend(handles = [t_line, fr_line, blue_line], loc = 2, fontsize=10, frameon = False)
# plt.gcf().subplots_adjust(bottom=0.15)
# plt.gcf().set_size_inches(6.4, 5.5)
# ax.grid(alpha=0.4, linestyle='--', linewidth = 2)

        # plt.show()

        # bar-plot with error bar
        # bx = plt.subplot(122)
        # x = [loadx[1], loadx[2], loadx[3], loadx[4]]
        # inc_t = [np.mean(t_inc_comp[loadx[1]]), np.mean(t_inc_comp[loadx[2]]), np.mean(t_inc_comp[loadx[3]]), np.mean(t_inc_comp[loadx[4]])]
        # inc_e = [np.mean(e_inc_comp[loadx[1]]), np.mean(e_inc_comp[loadx[2]]), np.mean(e_inc_comp[loadx[3]]), np.mean(e_inc_comp[loadx[4]])]
        # inc_p = [np.mean(p_inc_comp[loadx[1]]), np.mean(p_inc_comp[loadx[2]]), np.mean(p_inc_comp[loadx[3]]), np.mean(p_inc_comp[loadx[4]])]
        # er_t = [np.std(t_inc_comp[loadx[1]]), np.std(t_inc_comp[loadx[2]]), np.std(t_inc_comp[loadx[3]]), np.std(t_inc_comp[loadx[4]])]
        # er_e = [np.std(e_inc_comp[loadx[1]]), np.std(e_inc_comp[loadx[2]]), np.std(e_inc_comp[loadx[3]]), np.std(e_inc_comp[loadx[4]])]
        # er_p = [np.std(p_inc_comp[loadx[1]]), np.std(p_inc_comp[loadx[2]]), np.std(p_inc_comp[loadx[3]]), np.std(p_inc_comp[loadx[4]])]
        #
        # bx.bar(xt, inc_t, yerr = er_t, color = 'blue')
        # bx.bar(xe, inc_e, yerr = er_e, color = 'green')
        # bx.bar(xp, inc_p, yerr = er_p, color = 'orange')
        # bx.legend(['total', 'elastic', 'plastic'], fontsize = 8)
        # bx.set_xticks(xe)
        # bx.set_xticklabels([loadx[1], loadx[2], loadx[3], loadx[4]])
        # bx.set(xlabel='loading (g)', ylabel='Incremental Compliance(∆L%/MPa)', title='Incremental Compliance_' + name)

# plt.savefig(''f'{output}/Strain_modulus_seq.pdf', transparent = True)
plt.show()




