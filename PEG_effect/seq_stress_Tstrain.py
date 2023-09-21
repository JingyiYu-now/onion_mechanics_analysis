# script: Study progression of elasticity and plasticity with deformation method
# output value of T & E & P deformation and their incremental deformation;
# plot: line-plot for exact value, stack bar chart for incremental value

# import necessary package
import pandas as pd
import numpy as np
import os
import glob
from matplotlib import pyplot as plt
import matplotlib.lines as mlines
import matplotlib as mpl
import matplotlib
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Avenir'
plt.rcParams['font.size'] = 15
plt.rcParams['axes.linewidth'] = 2
mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['ytick.major.size'] = 5
mpl.rcParams['ytick.major.width'] = 2
import matplotlib as mpl
mpl.use('macosx')
from numpy import trapz
import statsmodels.api as sm

# User input
# folder1 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion2/10mmpermin_9.22.20'
# folder2 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion2/5mmpermin_9.22.20'
# folder3 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion2/5mmpermin_9.17.20'
# folder4 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion1'
# Onion3 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion3/5mmpermin'
# Onion4 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion4'
# Onion5 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion5'
# Con = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Sequential_Instron_With_Retracts_Tstrain_Data/40percent_PEG8k/40perctPEG8k/Control'
# PEG = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Sequential_Instron_With_Retracts_Tstrain_Data/40percent_PEG8k/40perctPEG8k/PEG'

# PEG Ca2+ sample
# ConE = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Sequential_Instron_With_Retracts_Tstrain_Data/40percent_PEG8k/Ca2_40perctPEG8k/buffer_con'
# Con_Ca = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Sequential_Instron_With_Retracts_Tstrain_Data/40percent_PEG8k/Ca2_40perctPEG8k/Ca2+_con'
# PEG_conE = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Sequential_Instron_With_Retracts_Tstrain_Data/40percent_PEG8k/Ca2_40perctPEG8k/buffer+PEG'
# PEG_Ca = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Sequential_Instron_With_Retracts_Tstrain_Data/40percent_PEG8k/Ca2_40perctPEG8k/Ca2+_PEG'

# open cell onions
# Con2 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Sequential_Instron_With_Retracts_Tstrain_Data/40percent_PEG8k/Open cell onions/onion1/Control'
# PEG2 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Sequential_Instron_With_Retracts_Tstrain_Data/40percent_PEG8k/Open cell onions/onion1/40pc_PEG_8k'
# Con3 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Sequential_Instron_With_Retracts_Tstrain_Data/40percent_PEG8k/Open cell onions/onion2/pH7.0'
# PEG3 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Sequential_Instron_With_Retracts_Tstrain_Data/40percent_PEG8k/Open cell onions/onion2/40pc_PEG8k'

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
group = [Ucon, UPEG, UReH]


# file = ''f'{folder}/Sequential_Instron__02_15_2020__14_30_58_SHORT.csv'
# output = ''f'{folder1}/summary1.csv'

thickness = 7
width = 3

strainT = [5,10, 15, 20, 25]


# define a function that calls the index of certain value
def ind(df, value):
    for i in range(len(df)):
        if df.position[i] == value:
            return i


# define a function that find the index where peels are under tension ( constant > 0.1g )
def ten_ind(pull, load):
    for i in range(len(pull) - 5):
        if (pull.load[i] > load) & (pull.load[i + 1] > load) & (pull.load[i + 2] > load) & \
                (pull.load[i + 3] > load) & (pull.load[i + 4] > load):
            index = i
            return index

def norm(pull):
    # strain calculated based on very original length (before plastic deformation)
    ori_p = ten_p(first_pull)
    # strain calculated based on length at the beginning of each pull
    # ori_p = ten_p(pull)

    pull['force_N'] = pull.load * 0.0098

    # engineering stress & strain
    pull['stress'] = pull.force_N / (thickness * width * 0.001)
    pull['strain'] = (pull.position - ori_p) / (5000 + ori_p)
    # true stress & strain
    # pull['strain'] = np.log(1 + (pull.position - ori_p) / (ori_p + 4500) )
    # pull['stress'] = pull.force_N/(thickness * width * (ori_p + 4500) / (pull.position + 4500) * 0.001)


# define a function that find the position where peels are under tension ( constant > 0.1g )
def ten_p(pull):
    for i in range(len(pull) - 5):
        if (pull.load[i] > 0.1) & (pull.load[i + 1] > 0.1) & (pull.load[i + 2] > 0.1) & \
                (pull.load[i + 3] > 0.1) & (pull.load[i + 4] > 0.1):
            point = pull.position[i]
            return point



# create data frame for results from each property
# t_defo = pd.DataFrame(columns=loadx)
# e_defo = pd.DataFrame(columns=loadx)
# pr_defo = pd.DataFrame(columns=loadx)
# p_defo = pd.DataFrame(columns=loadx)
# real_e_defo = pd.DataFrame(columns=loadx)
# t_inc_defo = pd.DataFrame(columns=loadx)
# e_inc_defo = pd.DataFrame(columns=loadx)
# p_inc_defo = pd.DataFrame(columns=loadx)


# extract data from each pull
# Use if multiple files need analysis
n = 0
g = 0
mark = '8'
color = ['green', 'orange','blue',  'grey']
for unit in group:
    Stress = pd.DataFrame(columns=strainT)

    if g == 0:
        mark = '8'
    elif g == 1:
        mark = 's'
    elif g == 2:
        mark = 'p'
    elif g == 3:
        mark = 'P'
    print(g)
    for folder in unit:

        for file in sorted(glob.glob(os.path.join(folder, '*SHORT.csv'))):
            whole_df = pd.read_table(file,
                                     delimiter=',',
                                     header=None,
                                     names=range(0, 2),
                                     )
            whole_df.columns = ['position', 'load']

            first_pull = whole_df[ind(whole_df, 'FIRST PULL') + 1:ind(whole_df, 'FIRST RETRACT')].astype(float).reset_index(
                drop=True)
            first_retract = whole_df[ind(whole_df, 'FIRST RETRACT') + 1:ind(whole_df, 'SECOND PULL')].astype(float).iloc[
                            ::-1].reset_index(drop=True)
            second_pull = whole_df[ind(whole_df, 'SECOND PULL') + 1:ind(whole_df, 'SECOND RETRACT')].astype(
                float).reset_index(drop=True)
            second_retract = whole_df[ind(whole_df, 'SECOND RETRACT') + 1:ind(whole_df, 'THIRD PULL')].astype(float).iloc[
                             ::-1].reset_index(drop=True)
            third_pull = whole_df[ind(whole_df, 'THIRD PULL') + 1:ind(whole_df, 'THIRD RETRACT')].astype(float).reset_index(
                drop=True)
            third_retract = whole_df[ind(whole_df, 'THIRD RETRACT') + 1:ind(whole_df, 'FOURTH PULL')].astype(float).iloc[
                            ::-1].reset_index(drop=True)
            fourth_pull = whole_df[ind(whole_df, 'FOURTH PULL') + 1:ind(whole_df, 'FOURTH RETRACT')].astype(
                float).reset_index(drop=True)
            fourth_retract = whole_df[ind(whole_df, 'FOURTH RETRACT') + 1:ind(whole_df, 'FIFTH PULL')].astype(float).iloc[
                             ::-1].reset_index(drop=True)
            fifth_pull = whole_df[ind(whole_df, 'FIFTH PULL') + 1:ind(whole_df, 'FIFTH RETRACT')].astype(float).reset_index(
                drop=True)
            fifth_retract = whole_df[ind(whole_df, 'FIFTH RETRACT') + 1:ind(whole_df, 'SIXTH PULL')].astype(float).iloc[
                            ::-1].reset_index(drop=True)
            sixth_pull = whole_df[ind(whole_df, 'SIXTH PULL') + 1:ind(whole_df, 'SIXTH RETRACT')].astype(float).reset_index(
                drop=True)
            sixth_retract = whole_df[ind(whole_df, 'SIXTH RETRACT') + 1:].astype(float).iloc[::-1].reset_index(drop=True)

            ori_p = ten_p(first_pull)


            curve = [first_pull, second_pull, third_pull, fourth_pull, fifth_pull, sixth_pull]
            retract = [first_retract, second_retract, third_retract, fourth_retract, fifth_retract,
                       sixth_retract]
            n += 1

            # calculate the stress and strain of the curve
            for i in range(len(curve)):
                norm(curve[i])

        # for sequential Instron - total strain
            Stress_l = []
            Stress_l.clear()
            for i in range(len(curve) - 1):
                Stress_l.append((curve[i].stress[len(curve[i]) - 1]))
        # add data of this file to the results data frame
            Stress.loc[len(Stress)] = Stress_l


    ## plot total value curve
    # ax = plt.subplot(121)
    # ax.plot(range(5), total_defo)
    # ax.plot(range(5), ela_defo)
    # ax.plot(range(5), pla_defo)
    # ax.legend(['total', 'elastic', 'plastic'])
    # ax.set_xticks(range(5))
    # ax.set_xticklabels(loadx)
    # ax.set(xlabel='loading (g)', ylabel='deformation(∆L/100gram)', title='Sequential deformation')

    ## plot incremental data
    # bx = plt.subplot(122)
    # bx = plt.subplot()
    # inc_t = [np.mean(t_inc_defo[8]), np.mean(t_inc_defo[12]), np.mean(t_inc_defo[16]), np.mean(t_inc_defo[20])]
    # inc_e = [np.mean(e_inc_defo[8]), np.mean(e_inc_defo[12]), np.mean(e_inc_defo[16]), np.mean(e_inc_defo[20])]
    # inc_p = [np.mean(p_inc_defo[8]), np.mean(p_inc_defo[12]), np.mean(p_inc_defo[16]), np.mean(p_inc_defo[20])]
    # er_t = [np.std(t_inc_defo[8]), np.std(t_inc_defo[12]), np.std(t_inc_defo[16]), np.std(t_inc_defo[20])]
    # er_e = [np.std(e_inc_defo[8]), np.std(e_inc_defo[12]), np.std(e_inc_defo[16]), np.std(e_inc_defo[20])]
    # er_p = [np.std(p_inc_defo[8]), np.std(p_inc_defo[12]), np.std(p_inc_defo[16]), np.std(p_inc_defo[20])]
    #
    # bx.bar(xt, inc_t, yerr = er_t)
    # bx.bar(xe, inc_e, yerr = er_e)
    # bx.bar(xp, inc_p, yerr = er_p)
    # bx.legend(['total', 'elastic', 'plastic'])
    # bx.set_xticks(xe)
    # bx.set_xticklabels(['8', '12', '16', '20'])
    # bx.set(xlabel='loading (g)', ylabel='Incremental deformation(∆L/100gram)', title='Incremental deformation')
    # plt.show()


# output the results in a file
# defo = pd.concat([t_defo, e_defo, p_defo], axis = 1)
# inc_defo = pd.concat([t_inc_defo, e_inc_defo, p_inc_defo], axis = 1)
# summary = pd.concat([defo, inc_defo], axis = 1)
# print(summary)
# summary.to_csv(output)

# Line-plot with error bar
    ax = plt.subplot(111)
    x = strainT
# C_t = [np.mean(t_defo[4]), np.mean(t_defo[8]), np.mean(t_defo[12]), np.mean(t_defo[16]), np.mean(t_defo[20])]
# C_e = [np.mean(e_defo[4]), np.mean(e_defo[8]), np.mean(e_defo[12]), np.mean(e_defo[16]), np.mean(e_defo[20])]
# C_p = [np.mean(p_defo[4]), np.mean(p_defo[8]), np.mean(p_defo[12]), np.mean(p_defo[16]), np.mean(p_defo[20])]
# er_t = [np.std(t_defo[4]), np.std(t_defo[8]), np.std(t_defo[12]), np.std(t_defo[16]), np.std(t_defo[20])]
# er_e = [np.std(e_defo[4]), np.std(e_defo[8]), np.std(e_defo[12]), np.std(e_defo[16]), np.std(e_defo[20])]
# er_p = [np.std(p_defo[4]), np.std(p_defo[8]), np.std(p_defo[12]), np.std(p_defo[16]), np.std(p_defo[20])]

    ave_stress = []
    std_stress = []
    for i in range(len(strainT)):
        ave_stress.append(np.mean(Stress.iloc[:, i]))
        std_stress.append(np.std(Stress.iloc[:, i]))

# jx = plt.subplot(111)
# jx.bar(target, mean_pull)

# plot with error bar
#     plt.errorbar(x, C_e, yerr = er_e, color='green')
#     plt.errorbar(x, C_p, yerr = er_p, color='orange')

# plot without error bar, without total value
    ms = 6
    # plt.plot(x, ave_stress, color='grey', linestyle = '--', marker = mark, markersize = ms)
    plt.errorbar(x, ave_stress, yerr = std_stress, capsize = 5, color= color[g], marker = '8', markersize = ms)
    # plt.plot(x, ave_stress, color='grey', linestyle = '--', marker = mark, markersize = ms)

    g += 1

    # plt.scatter(x, C_t, color='grey', s=12)
    #     plt.plot(x, C_e, color='C2')
    # plt.scatter(x, C_re, color='C2', alpha=0.6, s=12)
    # plt.scatter(x, C_p, color='C9', alpha=0.6, s=12)
    # plt.scatter(x, C_er, color='C3', alpha=0.6, s=12)
# ay = ax.twinx()
# ay.set_ylabel('Load (grams)')
# B = -1
# T = 40
# ay.set(ylim=[B, T])
# SB = B * 0.0098 / (thickness * width * 0.001)
# ST = T * 0.0098 / (thickness * width * 0.001)
# ax.set(ylim=[SB, ST])
x = [0, 5, 10, 15, 20, 25]
ax.set_xticks(x)
ax.set_xticklabels(x)
ax.set(xlim = [0, 30], ylim = [0, 17], xlabel='Target strain(%)', ylabel='Stress (MPa)')
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().set_size_inches(6.4, 5.5)
# ax.grid(alpha = 0.4, linestyle = '--')
blue_line = mlines.Line2D([], [], color = color[0], label = 'Control',  marker = '8')
# green_line = mlines.Line2D([], [], color = 'C6', label = 'elastic strain' )
C7_line = mlines.Line2D([], [], color = color[1], label = '40% PEG8k', marker = '8' )
orange_line = mlines.Line2D([], [], color = color[2], label = '40% PEG8k_rehydrated', marker = '8')
# C12_line = mlines.Line2D([], [], color = color[3], label = '40% PEG8k & 1mM Ca2+' , linestyle = '--', marker = 'P')
ax.legend(handles = [blue_line,C7_line, orange_line], loc = 2)


    # plt.show()

    # bar-plot with error bar
    # bx = plt.subplot(122)
    # x = [1,2,3,4]
    # inc_t = [np.mean(t_inc_defo[8]), np.mean(t_inc_defo[12]), np.mean(t_inc_defo[16]), np.mean(t_inc_defo[20])]
    # inc_e = [np.mean(e_inc_defo[8]), np.mean(e_inc_defo[12]), np.mean(e_inc_defo[16]), np.mean(e_inc_defo[20])]
    # inc_p = [np.mean(p_inc_defo[8]), np.mean(p_inc_defo[12]), np.mean(p_inc_defo[16]), np.mean(p_inc_defo[20])]
    # er_t = [np.std(t_inc_defo[8]), np.std(t_inc_defo[12]), np.std(t_inc_defo[16]), np.std(t_inc_defo[20])]
    # er_e = [np.std(e_inc_defo[8]), np.std(e_inc_defo[12]), np.std(e_inc_defo[16]), np.std(e_inc_defo[20])]
    # er_p = [np.std(p_inc_defo[8]), np.std(p_inc_defo[12]), np.std(p_inc_defo[16]), np.std(p_inc_defo[20])]

    # separate normal bar chart
    # bx.bar(xt, inc_t, yerr = er_t, color = 'blue')
    # bx.bar(xe, inc_e, yerr = er_e, color = 'green')
    # bx.bar(xp, inc_p, yerr = er_p, color = 'orange')
    # bx.legend(['total', 'elastic', 'plastic'])
    # bx.set_xticks(xe)
    # bx.set_xticklabels(['8', '12', '16', '20'])
    # bx.set(xlabel='loading (g)', ylabel='Incremental deformation(∆L%/MPa)', title='Incremental deformation_onion1')
    # plt.show()

    # stack bar chart

    # bx.bar(x, inc_e, yerr = er_e)
    # bx.bar(x, inc_p, yerr = er_p, bottom = inc_e)
    # bx.legend(['elastic', 'plastic'])
    # bx.set_xticks(x)
    # bx.set_xticklabels(['8', '12', '16', '20'])
    # bx.set(xlabel='loading (g)', ylabel='Incremental deformation(%)', title='Incremental deformation_3onion')
# plt.axhline(y = 0, color = 'black')

plt.savefig(''f'{output}/Ave_stress.pdf', transparent = True)
plt.show()