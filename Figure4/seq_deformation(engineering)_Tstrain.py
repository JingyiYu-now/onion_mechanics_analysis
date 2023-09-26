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


# file = ''f'{folder}/Sequential_Instron__02_15_2020__14_30_58_SHORT.csv'
# output = ''f'{folder1}/summary1.csv'

thickness = 7
width = 3

strainT = [5, 10, 15, 20, 25]


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
ls = ['-', '--']
for unit in group:
    t_defo = pd.DataFrame(columns=strainT)
    rec_defo = pd.DataFrame(columns=strainT)
    hys = pd.DataFrame(columns=strainT)
    p_defo = pd.DataFrame(columns=strainT)
    e_defo = pd.DataFrame(columns=strainT)
    t_inc_defo = pd.DataFrame(columns=strainT)
    e_inc_defo = pd.DataFrame(columns=strainT)
    p_inc_defo = pd.DataFrame(columns=strainT)

    ## for data outpur name
    if g == 0:
        name = 'Con'
    elif g == 1:
        name = 'PEG'


    for folder in unit:

        for file in sorted(glob.glob(os.path.join(folder, '*SHORT.csv'))):
            print(n)
            print(file)
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

        # for sequential Instron - total strain
            total_defo = []
            total_defo.clear()
            for i in range(len(curve) - 1):
                total_defo.append((curve[i].position[len(curve[i]) - 1] - ori_p) / 5000 * 100)
        # add data of this file to the results data frame
            t_defo.loc[len(t_defo)] = total_defo

        # for sequential Instron - recoverable deformation
            recoverable_defo = []
            for i in range(len(curve) - 1):
                recoverable_defo.append((curve[i].position[len(curve[i]) - 1] - ten_p(curve[i + 1])) / 5000 * 100)
            rec_defo.loc[len(rec_defo)] = recoverable_defo
            print(rec_defo)

        # for sequential Instron - plastic deformation
            pla_defo = []
            for i in range(len(curve) - 1):
                pla_defo.append((ten_p(curve[i + 1]) - ori_p) / 5000 * 100)
            # pla_forplot = [ -x for x in pla_defo]
            p_defo.loc[len(p_defo)] = pla_defo
            print(p_defo)

        # for sequential Instron - hysteresis strain
            hysteresis_defo = []
            for i in range(len(curve) - 1):
                hysteresis_defo.append((ten_p(retract[i]) - ten_p(curve[i + 1])) / 5000 * 100)
            hys.loc[len(hys)] = hysteresis_defo
            print(hys)

        # elastic deformation
            ela_defo = []
            for i in range(len(curve) - 1):
                ela_defo.append(recoverable_defo[i] - hysteresis_defo[i])
            e_defo.loc[len(e_defo)] = ela_defo
            print(e_defo)

    # Line-plot with error bar
    ax = plt.subplot(111)
    x = strainT

    C_rec = []
    C_ela = []
    C_p = []
    C_hys = []
    C_e_sd = []
    C_p_sd = []
    C_hys_sd = []
    for i in range(len(strainT)):
        C_rec.append(np.mean(rec_defo.iloc[:, i]))

        C_ela.append(np.mean(e_defo.iloc[:, i]))
        C_p.append(np.mean(p_defo.iloc[:, i]))
        C_hys.append(np.mean(hys.iloc[:, i]))

        C_e_sd.append(np.std(e_defo.iloc[:, i]))
        C_p_sd.append(np.std(p_defo.iloc[:, i]))
        C_hys_sd.append(np.std(hys.iloc[:, i]))

    total_chart = pd.DataFrame(columns=strainT)
    total_chart = e_defo
    total_chart.loc[len(total_chart)] = C_ela
    total_chart.loc[len(total_chart)] = C_e_sd
    # total_chart.to_csv(''f'{output}/' + name + '_deformation_elastic.csv')

    total_chart = pd.DataFrame(columns=strainT)
    total_chart = p_defo
    total_chart.loc[len(total_chart)] = C_p
    total_chart.loc[len(total_chart)] = C_p_sd
    # total_chart.to_csv(''f'{output}/' + name + '_deformation_plastic.csv')

    total_chart = pd.DataFrame(columns=strainT)
    total_chart = hys
    total_chart.loc[len(total_chart)] = C_hys
    total_chart.loc[len(total_chart)] = C_hys_sd
    # total_chart.to_csv(''f'{output}/' + name + '_deformation_hysteresis.csv')

# plot without error bar, without total value
    ms = 8
    lw = 2

    plt.errorbar(x, C_ela, yerr = C_e_sd, color='C2', linestyle = ls[g], marker = mark, markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)
    plt.errorbar(x, C_p, yerr = C_p_sd, color='C9', linestyle = ls[g], marker = mark, markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)
    plt.errorbar(x, C_hys, yerr = C_hys_sd, color='C3', linestyle = ls[g], marker = mark, markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)

    g += 1

x = [0, 5, 10, 15, 20, 25]
ax.set_xticks(x)
ax.set_xticklabels(x)
ax.set(xlim = [0, 28], ylim = [0, 14], xlabel='Total strain (%)', ylabel='Strain (%)')
C2_line = mlines.Line2D([], [], color = 'C2', label = 'Elastic strain')
C9_line = mlines.Line2D([], [], color ='C9', label ='Plastic strain')
C3_line = mlines.Line2D([], [], color ='C3', label ='hysteresis strain')
ax.legend(handles = [C9_line, C2_line, C3_line], loc = 2, frameon = False)


plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().set_size_inches(6.4, 5.5)
# plt.savefig(''f'{output}/Strain_hysteresis.pdf', transparent = True)
plt.show()