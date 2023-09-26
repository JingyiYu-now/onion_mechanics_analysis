# script: Study progression of elasticity and plasticity with deformation method
# output value of T & E & P deformation
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

loadx = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20]



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


# define a function that find the position where peels are under tension ( constant > 0.1g )
def ten_p(pull):
    for i in range(len(pull) - 5):
        if (pull.load[i] > 0.1) & (pull.load[i + 1] > 0.1) & (pull.load[i + 2] > 0.1) & \
                (pull.load[i + 3] > 0.1) & (pull.load[i + 4] > 0.1):
            point = pull.position[i]
            return point


# extract data from each pull
# Use if multiple files need analysis
n = 0
g = 0
mark = '8'
t_defo = pd.DataFrame(columns=loadx)
rec_defo = pd.DataFrame(columns=loadx)
hys_defo = pd.DataFrame(columns=loadx)
p_defo = pd.DataFrame(columns=loadx)
e_defo = pd.DataFrame(columns=loadx)
t_inc_defo = pd.DataFrame(columns=loadx)
e_inc_defo = pd.DataFrame(columns=loadx)
p_inc_defo = pd.DataFrame(columns=loadx)
for folder in group:
    # if g == 0:
    #     mark = '8'
    # else:
    #     mark = 's'
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

            curve = [first_pull, second_pull, third_pull, fourth_pull, fifth_pull, seventh_pull, eighth_pull,
                    ninth_pull, tenth_pull, eleventh_pull, twlveth_pull]
            retract = [first_retract, second_retract, third_retract, fourth_retract, fifth_retract,
                       seventh_retract, eighth_retract, ninth_retract, tenth_retract, eleventh_retract, twlveth_retract]
            n += 1

        ori_p = ten_p(first_pull)

    # for sequential Instron - total strain
        total_defo = []
        total_defo.clear()
        for i in range(len(curve) - 1):
            total_defo.append((curve[i].position[len(curve[i]) - 1] - ori_p) / 5000 * 100)
    # add data of this file to the results data frame
        t_defo.loc[len(t_defo)] = total_defo


    # for sequential Instron - elastic deformation
        recoverable_defo = []
        for i in range(len(curve) - 1):
            recoverable_defo.append((curve[i].position[len(curve[i]) - 1] - ten_p(curve[i + 1])) / 5000 * 100)
        rec_defo.loc[len(rec_defo)] = recoverable_defo


    # for sequential Instron - plastic deformation
        pla_defo = []
        for i in range(len(curve) - 1):
            pla_defo.append((ten_p(curve[i + 1]) - ori_p) / 5000 * 100)
        # pla_forplot = [ -x for x in pla_defo]
        p_defo.loc[len(p_defo)] = pla_defo

        # date reporter
        print(p_defo)

    # for sequential Instron - hysteresis deformation
        hysteresis_defo = []
        for i in range(len(curve) - 1):
            hysteresis_defo.append((ten_p(retract[i]) - ten_p(curve[i + 1])) / 5000 * 100)
        hys_defo.loc[len(rec_defo)] = hysteresis_defo

    # elastic deformation
        ela_defo = []
        for i in range(len(curve) - 1):
            ela_defo.append(recoverable_defo[i] - hysteresis_defo[i])
        e_defo.loc[len(rec_defo)] = ela_defo


# Line-plot with error bar
ax = plt.subplot(111)
# plot as stress
x = [a * 0.098 /3/7/0.01 for a in loadx]

C_t = []
C_rec = []
C_e = []
C_p = []
C_hys = []
er_t = []
er_e = []
er_p = []
er_hys = []
for i in range(len(loadx)):
    C_t.append(np.mean(t_defo.iloc[:, i]))
    C_e.append(np.mean(e_defo.iloc[:, i]))
    C_p.append(np.mean(p_defo.iloc[:, i]))
    C_hys.append(np.mean(hys_defo.iloc[:, i]))
    er_t.append(np.std(t_defo.iloc[:, i]))
    er_e.append(np.std(e_defo.iloc[:, i]))
    er_p.append(np.std(p_defo.iloc[:, i]))
    er_hys.append(np.std(hys_defo.iloc[:, i]))


# plot with error bar
ms = 4
lw = 2
mark = '8'
plt.errorbar(x, C_t, yerr = er_t, color= 'grey', alpha = 0.6, marker = mark, markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)
plt.errorbar(x, C_e, yerr = er_e, color='C2', alpha = 0.6, marker = mark, markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)
plt.errorbar(x, C_p, yerr = er_p, color='C9', alpha = 0.6, marker = mark, markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)
plt.errorbar(x, C_hys, yerr = er_hys, color='C3', alpha = 0.6, marker = mark, markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)
xlabel = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
ax.set_xticks(xlabel)
ax.set_xticklabels(xlabel)
ax.set_yticks(range(0, 35, 5))
ax.set(ylim = [-5, 35], xlabel='Stress (MPa)', ylabel='Deformation (%)')
blue_line = mlines.Line2D([], [], color = 'grey', label = 'Total strain')
C7_line = mlines.Line2D([], [], color = 'C2', label = 'Elastic strain')
orange_line = mlines.Line2D([], [], color = 'C9', label = 'Plastic strain')
C12_line = mlines.Line2D([], [], color = 'C3', label = 'Delayed elastic strain')
ax.legend(handles = [blue_line,C7_line, orange_line, C12_line], loc = 2, fontsize=10, frameon = False)
plt.gcf().subplots_adjust(bottom=0.15)
fig = plt.gcf()
fig.set_size_inches(6.4, 5.5)
size = fig.get_size_inches()
print(size)

# plt.savefig(''f'{output}/Strain_evolution.pdf', transparent = True)
plt.show()