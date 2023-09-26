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

# total = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_seq_instron/4,8,12,16,20_sequence/Total_3onion'
# folder4 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_seq_instron/4,8,12,16,20_sequence/seq_Onion_10'
# folder5 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_seq_instron/4,8,12,16,20_sequence/Onion_6.23.20_9'
# group = [folder1, folder2, folder3, folder4, folder5]

# 8,12,16,20,24 sequence
# seq1 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_seq_instron/8,12,16,20,24_sequence/Onion1'
# seq2 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_seq_instron/8,12,16,20,24_sequence/onion2_10'

name = 'Onion4'
# file = ''f'{folder}/Sequential_Instron__02_15_2020__14_30_58_SHORT.csv'
# output = ''f'{folder1}/summary1.csv'

thickness = 7
width = 3
# linear fitting percantage (%)
fitpc = 5

# for load sequence of 4-8-12-16-20
strainT = [5,10, 15, 20, 25]
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

def strain_p(pull, strain):
    for i in range(len(pull) - 5):
        if (pull.strain[i] > strain/100) & (pull.strain[i + 1] > strain/100) & (pull.strain[i + 2] > strain/100) & \
                (pull.strain[i + 3] > strain/100) & (pull.strain[i + 4] > strain/100):
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
    pull['strain'] = (pull.position - ori_p) / (5000 + ori_p)
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
def fit(pull, percentage, target_strain):
    if target_strain == 'N':
        fitting_part = pull[int(len(pull) * (100 - percentage) / 100):].reset_index(drop=True)
        # fit using gram as unit
        # z = np.polyfit(fitting_part.position, fitting_part.load, 1)

        # fit using stress (MPa)
        z = np.polyfit(fitting_part.strain, fitting_part.stress, 1)
        pull['fit_e'] = pull.strain * z[0] + z[1]
    else:
        cutted_data = pull[:strain_p(pull, target_strain)].reset_index(drop=True)
        fitting_part = cutted_data[int(len(cutted_data) * (100 - percentage) / 100):].reset_index(drop=True)
        # fit using gram as unit
        # z = np.polyfit(fitting_part.position, fitting_part.load, 1)

        # fit using stress (MPa)
        z = np.polyfit(fitting_part.strain, fitting_part.stress, 1)
        pull['fit_e'] = pull.strain * z[0] + z[1]
    # Xm = np.median(fitting_part.position)
    # compliance = 100 / (4500 + Xm) / z[0]
    compliance = 1/z[0]
    return compliance

def fit_p(pull, p_pull, percentage):
    cutted_data = pull[pull.strain < np.max(p_pull.strain)].reset_index(drop=True)
    fitting_part = cutted_data[int(len(cutted_data) * (100 - percentage) / 100):].reset_index(drop=True)
    # fit using gram as unit
    # z = np.polyfit(fitting_part.position, fitting_part.load, 1)

    # fit using stress (MPa)
    z = np.polyfit(fitting_part.strain, fitting_part.stress, 1)
    pull['fit_e'] = pull.strain * z[0] + z[1]

    # Xm = np.median(fitting_part.position)
    # compliance = 100 / (4500 + Xm) / z[0]
    compliance = 1/z[0]
    return compliance


def tar_ind(pull, load):
    if load == 'N':
        return len(pull)
    else:
        stress = load * 0.0098/ (thickness * width * 0.001)
        for i in range(len(pull) - 5):
            if (pull.stress[i] > stress) & (pull.stress[i + 1] > stress) & (pull.stress[i + 2] > stress) & \
                    (pull.stress[i + 3] > stress) & (pull.stress[i + 4] > stress):
                return i


# create data frame for results from each property


# extract data from each pull
# Use if multiple files need analysis
g = 0
mark = '8'
ls = ['-', '--']
for unit in group:

    # create data frame for results from each property
    t_comp = pd.DataFrame(columns=strainT)
    e_comp = pd.DataFrame(columns=strainT)
    p_comp = pd.DataFrame(columns=strainT)
    # r_comp = pd.DataFrame(columns=loadx)
    h_comp = pd.DataFrame(columns=strainT)
    # t_inc_comp = pd.DataFrame(columns=loadx)
    # e_inc_comp = pd.DataFrame(columns=loadx)
    # p_inc_comp = pd.DataFrame(columns=loadx)
    # r_inc_comp = pd.DataFrame(columns=loadx)
    # h_inc_comp = pd.DataFrame(columns=loadx)

    # if g == 0:
    #     ls
    # elif g == 1:
    #     mark = 's'
    # elif g == 2:
    #     mark = 'p'
    # elif g == 3:
    #     mark = 'P'

    for folder in unit:



        n = 0
        for file in sorted(glob.glob(os.path.join(folder, '*SHORT.csv'))):

            whole_df = pd.read_table(file,
                                     delimiter=',',
                                     header=None,
                                     names=range(0, 2),
                                     )
            whole_df.columns = ['position', 'load']

            first_pull = whole_df[ind(whole_df, 'FIRST PULL') + 1:ind(whole_df, 'FIRST RETRACT')].astype(
                float).reset_index(
                drop=True)
            first_retract = whole_df[ind(whole_df, 'FIRST RETRACT') + 1:ind(whole_df, 'SECOND PULL')].astype(
                float).iloc[
                            ::-1].reset_index(drop=True)
            second_pull = whole_df[ind(whole_df, 'SECOND PULL') + 1:ind(whole_df, 'SECOND RETRACT')].astype(
                float).reset_index(drop=True)
            second_retract = whole_df[ind(whole_df, 'SECOND RETRACT') + 1:ind(whole_df, 'THIRD PULL')].astype(
                float).iloc[
                             ::-1].reset_index(drop=True)
            third_pull = whole_df[ind(whole_df, 'THIRD PULL') + 1:ind(whole_df, 'THIRD RETRACT')].astype(
                float).reset_index(
                drop=True)
            third_retract = whole_df[ind(whole_df, 'THIRD RETRACT') + 1:ind(whole_df, 'FOURTH PULL')].astype(
                float).iloc[
                            ::-1].reset_index(drop=True)
            fourth_pull = whole_df[ind(whole_df, 'FOURTH PULL') + 1:ind(whole_df, 'FOURTH RETRACT')].astype(
                float).reset_index(drop=True)
            fourth_retract = whole_df[ind(whole_df, 'FOURTH RETRACT') + 1:ind(whole_df, 'FIFTH PULL')].astype(
                float).iloc[
                             ::-1].reset_index(drop=True)
            fifth_pull = whole_df[ind(whole_df, 'FIFTH PULL') + 1:ind(whole_df, 'FIFTH RETRACT')].astype(
                float).reset_index(
                drop=True)
            fifth_retract = whole_df[ind(whole_df, 'FIFTH RETRACT') + 1:ind(whole_df, 'SIXTH PULL')].astype(float).iloc[
                            ::-1].reset_index(drop=True)
            sixth_pull = whole_df[ind(whole_df, 'SIXTH PULL') + 1:ind(whole_df, 'SIXTH RETRACT')].astype(
                float).reset_index(
                drop=True)
            sixth_retract = whole_df[ind(whole_df, 'SIXTH RETRACT') + 1:].astype(float).iloc[::-1].reset_index(
                drop=True)

            ori_p = ten_p(first_pull)
            n += 1


            pull = [first_pull, second_pull, third_pull, fourth_pull, fifth_pull, sixth_pull]
            retract = [first_retract, second_retract, third_retract, fourth_retract, fifth_retract,sixth_retract]

            # calculate the stress and strain of the curve
            for i in range(len(pull)):
                norm(pull[i])
                norm(retract[i])

            target = [5, 10, 15, 20, 25, 'N']
            # loop for seuqntial Instron
            for i in range(len(pull)):
                # retract[i].drop(retract[i].loc[purify_p(retract[i], target[i]): len(retract[i])].index, inplace=True)
                # retract[i].drop(retract[i].loc[0:ten_ind(retract[i], 0.1)].index, inplace=True)
                # pull[i].drop(pull[i].loc[0:ten_ind(pull[i], 0.1)].index, inplace=True)

                retract[i].reset_index(inplace=True)
                rmindex = []
                purify_revrs(retract[i])
                retract[i].drop(retract[i].loc[rmindex].index, inplace=True)

                retract[i] = retract[i].reset_index(drop=True)
                pull[i] = pull[i].reset_index(drop=True)

        # for sequential Instron - total compliance
            total_comp = []
            for i in range(len(pull)-1):
                total_comp.append(fit(pull[i], fitpc, 'N'))

        # add data of this file to the results data frame
            t_comp.loc[len(t_comp)] = total_comp
            print(t_comp)
        # incremental total compliance
        #     inc_t = [fit(first_pull, fitpc, 'N')]
        #     for i in range(len(total_comp)):
        #         if i > 0:
        #             inc_t.append(total_comp[i] - total_comp[i - 1])
        #     t_inc_comp.loc[len(t_inc_comp)] = inc_t
        #
        # to save total fit in a different column
        #     for i in range(len(pull)-1):
        #         pull[i]['fit'] = pull[i].fit_e


        # for sequential Instron - elastic compliance
            ela_comp = []
            for i in range(1, len(pull)-1):
                ela_comp.append(fit(pull[i], fitpc, strainT[i-1]))
            ela_comp.append(fit(pull[-1], fitpc, 'N'))

            ## fitting using the position
            # for i in range(1, len(pull)):
            #     ela_comp.append(fit_p(pull[i], pull[i - 1], fitpc))
            e_comp.loc[len(e_comp)] = ela_comp
            #
            # incremental elastic compliance
            # inc_e = [ela_comp[0]]
            # for i in range(len(ela_comp)):
            #     if i > 0:
            #         inc_e.append(ela_comp[i] - ela_comp[i - 1])
            # e_inc_comp.loc[len(e_inc_comp)] = inc_e
            #
            #
        # for sequential Instron - plastic compliance
        #     pla_comp = []
        #     for i in range(len(pull)-2):
        #         pla_comp.append(fit(pull[i], fitpc, 'N') - fit(pull[i+1], fitpc, loadx[0]))
        #     pla_comp.append(fit(pull[-2], fitpc, 'N') - fit(pull[-1], fitpc, 'N'))
        # # pla_forplot = [ -x for x in pla_comp]
        #     p_comp.loc[len(p_comp)] = pla_comp

            pla_comp = []
            for i in range(len(strainT)):
                pla_comp.append(total_comp[i] - ela_comp[i])
            # pla_forplot = [ -x for x in pla_comp]
            p_comp.loc[len(p_comp)] = pla_comp

        # incremental plastic compliance
        #     inc_p = [fit(first_pull, fitpc, 'N') - fit(second_pull, fitpc, loadx[0])]
        #     for i in range(len(pla_comp)):
        #         if i > 0:
        #             inc_p.append(pla_comp[i] - pla_comp[i - 1])
        #     p_inc_comp.loc[len(p_inc_comp)] = inc_p

            # for sequential Instron - retract compliance
            # retra_comp = []
            # for i in range(len(pull) - 1):
            #     retra_comp.append(fit(retract[i], fitpc, 'N'))
            # # pla_forplot = [ -x for x in pla_comp]
            # r_comp.loc[len(r_comp)] = retra_comp
            # fit(retract[-1], fitpc, 'N')

            # incremental retract compliance
            # inc_r = [retra_comp[0]]
            # for i in range(len(retra_comp)):
            #     if i > 0:
            #         inc_r.append(retra_comp[i] - retra_comp[i - 1])
            # r_inc_comp.loc[len(r_inc_comp)] = inc_r
            #
            # for sequential Instron - hysteresis compliance
            # hys_comp = []
            # for i in range(len(pull) - 1):
            #     hys_comp.append(ela_comp[i] - retra_comp[i])
            # # pla_forplot = [ -x for x in pla_comp]
            # h_comp.loc[len(h_comp)] = hys_comp

            # incremental hysteresis compliance
            # inc_h = [hys_comp[0]]
            # for i in range(len(hys_comp)):
            #     if i > 0:
            #         inc_h.append(hys_comp[i] - hys_comp[i - 1])
            # h_inc_comp.loc[len(h_inc_comp)] = inc_h

            # if n > 1:
            #     break
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
            # ax.plot(range(5), total_comp)
            # ax.plot(range(5), ela_comp)
            # ax.plot(range(5), pla_comp)
            # ax.legend(['total', 'elastic', 'plastic'])
            # ax.set_xticks(range(5))
            # ax.set_xticklabels(loadx)
            # ax.set(xlabel='loading (g)', ylabel='Compliance(∆L/100gram)', title='Sequential Compliance')

            ## plot incremental data
            # bx = plt.subplot(122)
            # bx = plt.subplot()
            # inc_t = [np.mean(t_inc_comp[loadx[1]]), np.mean(t_inc_comp[loadx[2]]), np.mean(t_inc_comp[loadx[3]]), np.mean(t_inc_comp[loadx[4]])]
            # inc_e = [np.mean(e_inc_comp[loadx[1]]), np.mean(e_inc_comp[loadx[2]]), np.mean(e_inc_comp[loadx[3]]), np.mean(e_inc_comp[loadx[4]])]
            # inc_p = [np.mean(p_inc_comp[loadx[1]]), np.mean(p_inc_comp[loadx[2]]), np.mean(p_inc_comp[loadx[3]]), np.mean(p_inc_comp[loadx[4]])]
            # er_t = [np.std(t_inc_comp[loadx[1]]), np.std(t_inc_comp[loadx[2]]), np.std(t_inc_comp[loadx[3]]), np.std(t_inc_comp[loadx[4]])]
            # er_e = [np.std(e_inc_comp[loadx[1]]), np.std(e_inc_comp[loadx[2]]), np.std(e_inc_comp[loadx[3]]), np.std(e_inc_comp[loadx[4]])]
            # er_p = [np.std(p_inc_comp[loadx[1]]), np.std(p_inc_comp[loadx[2]]), np.std(p_inc_comp[loadx[3]]), np.std(p_inc_comp[loadx[4]])]
            #
            # bx.bar(xt, inc_t, yerr = er_t)
            # bx.bar(xe, inc_e, yerr = er_e)
            # bx.bar(xp, inc_p, yerr = er_p)
            # bx.legend(['total', 'elastic', 'plastic'])
            # bx.set_xticks(xe)
            # bx.set_xticklabels(['8', '12', '16', '20'])
            # bx.set(xlabel='loading (g)', ylabel='Incremental Compliance(∆L/100gram)', title='Incremental Compliance')
            # plt.show()


            # output the results in a file
            # comp = pd.concat([t_comp, e_comp, p_comp], axis = 1)
            # inc_comp = pd.concat([t_inc_comp, e_inc_comp, p_inc_comp], axis = 1)
            # summary = pd.concat([comp, inc_comp], axis = 1)
            # print(summary)
            # summary.to_csv(output)


        ####### plot the stress-train curve and fitting curve
            # ela_target = ['N', 2, 4, 6, 8, 10, 12, 14, 16, 18, 'N']
            # a = plt.subplot(111)
            # for j in range(len(pull)):
            #     # for target loading fitting plot
            #     # PL = int(tar_ind(pull[j], ela_target[j]) * (100 - fitpc) / 100)
            #     # PR = tar_ind(pull[j], ela_target[j])
            #     # for target position fitting plot
            #     if j > 0:
            #         PR = pos_p(pull[j], pull[j - 1])
            #
            #     else:
            #         PR = pos_p(pull[j], pull[j])
            #     PL = int(PR * (100 - fitpc) / 100)
            #     # for total comliance fitting plot border
            #     L = int((100-fitpc) * 0.01 * len(pull[j]))
            #     # a.plot(smpull[j].strain , smpull[j].stress, color = color[j])
            #     a.plot(pull[j].strain * 100, pull[j].stress)
            #     a.plot(retract[j].strain * 100, retract[j].stress, linestyle = '--')
            #     # a.scatter(smretract[j].strain *100, smretract[j].stress, c=derv(smretract[j], intv), cmap = map, s = 4)
            #     # a.plot(smpull[j].strain[int(len(smpull[j]) * (100 - per) / 100):] * 100 , smpull[j].fit_e[int(len(smpull[j]) * (100 - per) / 100):], color = 'black')
            #     # a.plot(smretract[j].strain[int(len(smretract[j]) * (100 - per) / 100):] * 100, smretract[j].fit_e[int(len(smretract[j]) * (100 - per) / 100):], color = 'black')
            #     if j < len(pull)-1:
            #         a.plot(pull[j].strain[L:] * 100, pull[j].fit[L:], color='black', linewidth=2)
            #     if j > 0:
            #         a.plot(pull[j].strain[PL:PR] * 100, pull[j].fit_e[PL:PR], color='grey', linewidth=2)
            #     if j < len(pull)-1:
            #         a.plot(retract[j].strain[int(len(retract[j]) * (100 - fitpc) / 100):] * 100,
            #            retract[j].fit_e[int(len(retract[j]) * (100 - fitpc) / 100):], linestyle = '--', color='black')
            # a.set(xlabel = 'Strain(%)', ylabel = 'Stress(MPa)', title = 'Sequential cyclic loading')
            # # setting for gram y axis
            # ax0y = a.twinx()
            # ax0y.set_ylabel('Load (grams)')
            # B = -1
            # T = 21
            # ax0y.set(ylim=[B, T])
            # SB = B * 0.0098 / (thickness * width * 0.001)
            # ST = T * 0.0098 / (thickness * width * 0.001)
            # a.set(ylim=[SB, ST])
            # a.grid(alpha=0.4, linestyle='--')
            # plt.show()

            # Line-plot with error bar
    ax = plt.subplot(111)
    x = strainT

    C_t = []
    C_e = []
    C_p = []
    C_r = []
    C_h = []
    er_t = []
    er_e = []
    er_p = []
    er_r = []
    er_h = []
    for i in range(len(strainT)):
        C_t.append(np.mean(t_comp.iloc[:, i]) * 100)
        C_e.append(np.mean(e_comp.iloc[:, i]) * 100)
        C_p.append(np.mean(p_comp.iloc[:, i]) * 100)
        # C_r.append(np.mean(r_comp.iloc[:, i]) * 100)
        # C_h.append(np.mean(h_comp.iloc[:, i]) * 100)
        er_t.append(np.std(t_comp.iloc[:, i]) * 100)
        er_e.append(np.std(e_comp.iloc[:, i]) * 100)
        er_p.append(np.std(p_comp.iloc[:, i]) * 100)
        # er_r.append(np.std(r_comp.iloc[:, i]) * 100)
        # er_h.append(np.std(h_comp.iloc[:, i]) * 100)
        # er_t = [np.std(t_comp[loadx[0]]), np.std(t_comp[loadx[1]]), np.std(t_comp[loadx[2]]), np.std(t_comp[loadx[3]]), np.std(t_comp[loadx[4]])]
        # er_e = [np.std(e_comp[loadx[0]]), np.std(e_comp[loadx[1]]), np.std(e_comp[loadx[2]]), np.std(e_comp[loadx[3]]), np.std(e_comp[loadx[4]])]
        # er_p = [np.std(p_comp[loadx[0]]), np.std(p_comp[loadx[1]]), np.std(p_comp[loadx[2]]), np.std(p_comp[loadx[3]]), np.std(p_comp[loadx[4]])]

    # plot with error bar
    # plt.errorbar(x, C_t, yerr = er_t, color= 'C7', alpha = 0.6)
    # plt.errorbar(x, C_e, yerr = er_e, color='C1', alpha = 0.6)
    # plt.errorbar(x, C_p, yerr = er_p, color='C9', alpha = 0.6)
    # plt.errorbar(x, C_r, yerr = er_e, color='C2', alpha = 0.6)
    # plt.errorbar(x, C_h, yerr = er_p, color='C3', alpha = 0.6)

    # plot without error bar, without total value
    ms = 6
    lw = 2
    mark = '8'
    # plt.plot(x,C_t, color = 'C7', alpha = 0.6, linestyle = '--', marker = mark, markersize = ms)
    # plt.plot(x,C_p, color = 'C9', alpha = 0.6, linestyle = '--', marker = mark, markersize =ms)
    # plt.errorbar(x, C_t, yerr = er_t, color= 'C7', linestyle = ls[g], marker = mark, markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)
    plt.errorbar(x, C_p, yerr = er_p, color= 'C9', linestyle = ls[g], marker = mark, markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)
    plt.errorbar(x, C_e, yerr = er_e, color= 'C2', linestyle = ls[g], marker = mark, markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)

    # plt.plot(x,C_h, color = 'C3', alpha = 0.6, linestyle = '--', marker = mark, markersize = ms)
    # plt.plot(x,C_e, color = 'C2', alpha = 0.6, linestyle = '--', marker = mark, markersize = ms)
    # plt.plot(x,C_r, color = 'C1', alpha = 0.6, linestyle = '--', marker = mark, markersize = ms)
    g += 1

        # plt.scatter(x, C_t, color='C7', alpha=0.6, s=12)
        # plt.scatter(x, C_p, color='C9', alpha=0.6, s=12)
        # plt.scatter(x, C_h, color='C3', alpha=0.6, s=12)
        # plt.scatter(x, C_e, color='C2', alpha=0.6, s=12)
        # plt.scatter(x, C_r, color='C1', alpha=0.6, s=12)
ax.set_xticks(strainT)
ax.set_xticklabels(strainT)
ax.set(xlabel='Strain (%)', ylabel='Compliance(∆L%/MPa)')
# t_line = mlines.Line2D([], [], color = 'C7', label = 'Total compliance')
green_line = mlines.Line2D([], [], color = 'C9', label = 'Plastic compliance')
# orange_line = mlines.Line2D([], [], color = 'C3', label = 'Hysteresis compliance', linestyle = '--')
fr_line = mlines.Line2D([], [], color = 'C2', label = 'Elastic loading compliance')
# blue_line = mlines.Line2D([], [], color = 'C1', label = 'Retract compliance', linestyle = '--')
ax.legend(handles = [green_line, fr_line], loc = 2, fontsize = 8, frameon = False)
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().set_size_inches(6.4, 5.5)
# ax.grid(alpha=0.4, linestyle='--')

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

# plt.savefig(''f'{output}/Compliance_withE.pdf', transparent = True)
plt.show()

# stack bar chart
# bx.bar(x, inc_e, yerr = er_e, alpha = 0.4)
# bx.bar(x, inc_p, yerr = er_p, bottom = inc_e, alpha = 0.4)
# plt.axhline(y = 0, color = 'black')
# plt.show()






##### find plotting part: end-fitting part
# fire = first_pull[int(len(first_pull) * 0.9):]
# se = second_pull[int(len(second_pull) * 0.9):]
# te = third_pull[int(len(third_pull) * 0.9):]
# fe = fourth_pull[int(len(fourth_pull) * 0.9):]
# fie = fifth_pull[int(len(fifth_pull) * 0.9):]
# sie = sixth_pull[int(len(sixth_pull) * 0.9):]


### find plotting part: mid-fitting part
# sef = second_pull[int(ten_ind(second_pull,loadx[0]) * 0.9):int(ten_ind(second_pull,loadx[0]))].reset_index(drop=True)
# tf = third_pull[int(ten_ind(third_pull,loadx[1]) * 0.9):int(ten_ind(third_pull,loadx[1]))].reset_index(drop=True)
# fof = fourth_pull[int(ten_ind(fourth_pull,loadx[2]) * 0.9):int(ten_ind(fourth_pull,loadx[2]))].reset_index(drop=True)
# fif = fifth_pull[int(ten_ind(fifth_pull,loadx[3]) * 0.9):int(ten_ind(fifth_pull,loadx[3]))].reset_index(drop=True)
# sif = sixth_pull[int(len(sixth_pull) * (100 - 10) / 100):].reset_index(drop=True)

### plot position - load curve
# for i in range(len(curve)):
#     plt.plot(curve[i].position, curve[i].load, color = 'blue')

#### plot the end-fit linear curve
# plt.plot(fire.position, fire.fit_e, color = 'red')
# plt.plot(se.position, se.fit_e, color = 'red')
# plt.plot(te.position, te.fit_e, color = 'red')
# plt.plot(fe.position, fe.fit_e, color = 'red')
# plt.plot(fie.position, fie.fit_e, color = 'red')
# plt.plot(sie.position, sie.fit_e, color = 'red')

#### plot the mid-fit linear curve
# plt.plot(sef.position, sef.fit, color = 'green', linewidth = 2)
# plt.plot(tf.position, tf.fit, color = 'green', linewidth = 2)
# plt.plot(fof.position, fof.fit, color = 'green', linewidth = 2)
# plt.plot(fif.position, fif.fit, color = 'green', linewidth = 2)

#### plot setting
# ax = plt.subplot(111)
# blue_line = mlines.Line2D([], [], color = 'blue', label = 'load-position curve' )
# green_line = mlines.Line2D([], [], color = 'green', label = 'Mid-curve_fit' )
# red_line = mlines.Line2D([], [], color = 'red', label = 'End_curve_fit' )
# ax.legend(handles = [blue_line, red_line, green_line])
# ax.set(xlabel = 'position', ylabel = 'load(g)')
# plt.show()

############## paralle comparison
