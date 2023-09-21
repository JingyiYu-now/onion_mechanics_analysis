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
output = '/Users/jingyiyu/Documents/Cosgrovelab/manuscript/CW_mechanics_Instron/Data_collection/Python_plot/For new Fig4'

# User input
folder1 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion2/10mmpermin_9.22.20'
folder2 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion2/5mmpermin_9.22.20'
folder3 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion2/5mmpermin_9.17.20'
folder4 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion1'
Onion3 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion3/5mmpermin'
Onion4 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion4'
Onion5 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion5'
Onion6 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load/Onion6'

RT = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/contineous_rate_change/deformation'

# group = [folder2, Onion3, Onion5, Onion6]
group = [RT]


# file = ''f'{folder}/Sequential_Instron__02_15_2020__14_30_58_SHORT.csv'

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
e_defo = pd.DataFrame(columns=loadx)
pr_defo = pd.DataFrame(columns=loadx)
p_defo = pd.DataFrame(columns=loadx)
real_e_defo = pd.DataFrame(columns=loadx)
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

    # incremental total deformation
        inc_t = [total_defo[0]]
        for i in range(len(total_defo)):
            if i > 0:
                inc_t.append(total_defo[i] - total_defo[i - 1])
        t_inc_defo.loc[len(t_inc_defo)] = inc_t

    # for sequential Instron - elastic deformation
        ela_defo = []
        for i in range(len(curve) - 1):
            ela_defo.append((curve[i].position[len(curve[i]) - 1] - ten_p(curve[i + 1])) / 5000 * 100)
        e_defo.loc[len(e_defo)] = ela_defo


    # incremental elastic deformation
        inc_e = []
        for i in range(len(curve) - 1):
            if i == 0:
                inc_e.append((curve[i].position[len(curve[i]) - 1] - ten_p(curve[i + 1])) / 5000 * 100)
            else:
                inc_e.append(((curve[i].position[len(curve[i]) - 1] - ten_p(curve[i + 1])) - (
                            curve[i - 1].position[len(curve[i - 1]) - 1] - ten_p(curve[i]))) / 5000 * 100)
        e_inc_defo.loc[len(e_inc_defo)] = inc_e

    # for sequential Instron - plastic deformation
        pla_defo = []
        for i in range(len(curve) - 1):
            pla_defo.append((ten_p(curve[i + 1]) - ori_p) / 5000 * 100)
        # pla_forplot = [ -x for x in pla_defo]
        p_defo.loc[len(p_defo)] = pla_defo

        # date reporter
        print(p_defo)

    # incremental plastic deformation
        inc_p = []
        for i in range(len(curve) - 1):
            inc_p.append((ten_p(curve[i + 1]) - ten_p(curve[i])) / 5000 * 100)
        p_inc_defo.loc[len(p_inc_defo)] = inc_p

    # for sequential Instron - recoverable elastic deformation
        pla_rc_defo = []
        for i in range(len(curve)):
            pla_rc_defo.append((ten_p(retract[i]) - ten_p(curve[i])) / 5000 * 100)
        print(pla_rc_defo)
        pr_defo.loc[len(e_defo)] = pla_rc_defo

    # real elastic deformation
        real_ela_defo = []
        for i in range(len(curve) - 1):
            real_ela_defo.append(ela_defo[i] - pla_rc_defo[i])
        real_e_defo.loc[len(e_defo)] = real_ela_defo
    ###### results plotting
    ##### set x
    # n = 1
    # t = 3
    # d = 4
    # w = 0.8
    # xt = [t * x + w * n for x in range(d)]
    #
    # n = 2
    # t = 3
    # d = 4
    # w = 0.8
    # xe = [t * x + w * n for x in range(d)]
    #
    # n = 3
    # t = 3
    # d = 4
    # w = 0.8
    # xp = [t * x + w * n for x in range(d)]

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
# x = loadx
# plot as stress
x = [a * 0.098 /3/7/0.01 for a in loadx]

# C_t = [np.mean(t_defo[4]), np.mean(t_defo[8]), np.mean(t_defo[12]), np.mean(t_defo[16]), np.mean(t_defo[20])]
# C_e = [np.mean(e_defo[4]), np.mean(e_defo[8]), np.mean(e_defo[12]), np.mean(e_defo[16]), np.mean(e_defo[20])]
# C_p = [np.mean(p_defo[4]), np.mean(p_defo[8]), np.mean(p_defo[12]), np.mean(p_defo[16]), np.mean(p_defo[20])]
# er_t = [np.std(t_defo[4]), np.std(t_defo[8]), np.std(t_defo[12]), np.std(t_defo[16]), np.std(t_defo[20])]
# er_e = [np.std(e_defo[4]), np.std(e_defo[8]), np.std(e_defo[12]), np.std(e_defo[16]), np.std(e_defo[20])]
# er_p = [np.std(p_defo[4]), np.std(p_defo[8]), np.std(p_defo[12]), np.std(p_defo[16]), np.std(p_defo[20])]

C_t = []
C_e = []
C_re = []
C_p = []
C_er = []
er_t = []
er_re = []
er_p = []
er_er = []
for i in range(len(loadx)):
    C_t.append(np.mean(t_defo.iloc[:, i]))
    # C_e.append(np.mean(e_defo.iloc[:, i]))
    C_re.append(np.mean(real_e_defo.iloc[:, i]))
    C_p.append(np.mean(p_defo.iloc[:, i]))
    C_er.append(np.mean(pr_defo.iloc[:, i]))
    er_t.append(np.std(t_defo.iloc[:, i]))
    er_re.append(np.std(real_e_defo.iloc[:, i]))
    er_p.append(np.std(p_defo.iloc[:, i]))
    er_er.append(np.std(pr_defo.iloc[:, i]))

# jx = plt.subplot(111)
# jx.bar(target, mean_pull)

# plot with error bar
ms = 4
lw = 2
mark = '8'
plt.errorbar(x, C_t, yerr = er_t, color= 'grey', alpha = 0.6, marker = mark, markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)
plt.errorbar(x, C_re, yerr = er_re, color='C2', alpha = 0.6, marker = mark, markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)
plt.errorbar(x, C_p, yerr = er_p, color='C9', alpha = 0.6, marker = mark, markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)
plt.errorbar(x, C_er, yerr = er_er, color='C3', alpha = 0.6, marker = mark, markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)
# plot without error bar, without total value
# ms = 4
# mark = '8'
# plt.plot(x, C_t, color='grey', linestyle = '--', marker = mark, markersize = ms)
#     plt.plot(x, C_e, color='C2')
# plt.plot(x, C_re, color='C2', alpha = 0.6, linestyle = '--', marker = mark, markersize = ms)
# plt.plot(x, C_p, color='C9', alpha = 0.6, linestyle = '--', marker = mark, markersize = ms)
# plt.plot(x, C_er, color='C3', alpha = 0.6, linestyle = '--', marker = mark, markersize = ms)

# plt.errorbar(x, C_t, color='grey', linestyle = '--', marker = mark, markersize = ms)
#     plt.plot(x, C_e, color='C2')
# plt.errorbar(x, C_re, color='C2', alpha = 0.6, linestyle = '--', marker = mark, markersize = ms)
# plt.errorbar(x, C_p, color='C9', alpha = 0.6, linestyle = '--', marker = mark, markersize = ms)
# plt.errorbar(x, C_er, color='C3', alpha = 0.6, linestyle = '--', marker = mark, markersize = ms)


# g =+ 1

    # plt.scatter(x, C_t, color='grey', s=12)
    #     plt.plot(x, C_e, color='C2')
    # plt.scatter(x, C_re, color='C2', alpha=0.6, s=12)
    # plt.scatter(x, C_p, color='C9', alpha=0.6, s=12)
    # plt.scatter(x, C_er, color='C3', alpha=0.6, s=12)
xlabel = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
ax.set_xticks(xlabel)
ax.set_xticklabels(xlabel)
ax.set_yticks(range(0, 35, 5))
ax.set(ylim = [-5, 35], xlabel='Stress (MPa)', ylabel='Deformation (%)')
# ax.grid(alpha = 0.4, linestyle = '--')
blue_line = mlines.Line2D([], [], color = 'grey', label = 'Total strain')
# green_line = mlines.Line2D([], [], color = 'C6', label = 'elastic strain' )
C7_line = mlines.Line2D([], [], color = 'C2', label = 'Elastic strain')
orange_line = mlines.Line2D([], [], color = 'C9', label = 'Plastic strain')
C12_line = mlines.Line2D([], [], color = 'C3', label = 'Hysteresis strain')
ax.legend(handles = [blue_line,C7_line, orange_line, C12_line], loc = 2, fontsize=10, frameon = False)
plt.gcf().subplots_adjust(bottom=0.15)
fig = plt.gcf()
fig.set_size_inches(6.4, 5.5)
size = fig.get_size_inches()
print(size)
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

# plt.savefig(''f'{output}/Strain_evolution.pdf', transparent = True)
plt.show()