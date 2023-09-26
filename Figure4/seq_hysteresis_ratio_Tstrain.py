# script: Study progression of elasticity and plasticity with strain energy density method
# output value of T & E & P strain energy density and their incremental strain energy density;
# plot: line-plot for exact value, stack bar chart for incremental value

# import necessary package
import pandas as pd
import numpy as np
import os
import glob
from matplotlib import pyplot as plt
import matplotlib.lines as mlines
import statsmodels.api as sm
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

color = ['black', 'blue', 'green', 'orange']
thickness = 7
width = 3
smf = 0.04

# for load sequence of 4-8-12-16-20
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
        elif i > len(pull) - 7:
            print('no fit')

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
    pull['strain'] = (pull.position - ori_p) / (5000 + ori_p)


def smooth(pull, smf):
    lowess = sm.nonparametric.lowess
    z = lowess(pull.stress, pull.strain,frac= smf)
    return z

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


# define a function to calculate the area for total pull
def area_t(pull):
    area = 0
    for i in range(len(pull) - 1):
        area += (pull.stress[i] + pull.stress[i+1]) * abs(pull.strain[i + 1] - pull.strain[i]) / 2
    # print(str(area) + 'total')
    return area


# calculate the area according to the previous pull max strain
def area_s(pull, pre_pull):
    area = 0
    cutted_data = pull[pull.strain < np.max(pre_pull.strain)].reset_index(drop=True)
    for i in range(len(cutted_data)):
        area += (pull.stress[i] + pull.stress[i+1]) * abs(pull.strain[i + 1] - pull.strain[i]) / 2
    # print(str(area) + 'strain')
    return area

def FindZeroStrain(pull):
    for i in range(len(pull) - 5):
        if (pull.position[i] > 0) & (pull.position[i + 1] > 0) & (pull.position[i + 2] > 0) & \
                (pull.position[i + 3] > 0) & (pull.position[i + 4] > 0):
            return i

def max(pull):
    max = 0
    for i in range(len(pull)):
        if pull.stress[i] > max:
            max = pull.stress[i]
    return max


n = 0
# extract data from each pull
# Use if multiple files need analysis
g = 0
h = ['use default', '/']
ls = ['-', '--']
for unit in group:
    t_energ = pd.DataFrame(columns=strainT)
    e_energ = pd.DataFrame(columns=strainT)
    p_energ = pd.DataFrame(columns=strainT)
    store_energ = pd.DataFrame(columns=strainT)
    hyst_energy = pd.DataFrame(columns=strainT)
    t_inc_energ = pd.DataFrame(columns=strainT)
    e_inc_energ = pd.DataFrame(columns=strainT)
    p_inc_energ = pd.DataFrame(columns=strainT)
    h_inc_energ = pd.DataFrame(columns=strainT)
    s_inc_energ = pd.DataFrame(columns=strainT)

    t_hys = pd.DataFrame(columns=strainT)
    p_hys = pd.DataFrame(columns=strainT)
    e_hys = pd.DataFrame(columns=strainT)
    T_stress = pd.DataFrame(columns=strainT)
    sto_ratio_DF = pd.DataFrame(columns=strainT)

## for data outpur name
    if g == 0:
        name = 'Con'
    elif g == 1:
        name = 'PEG'

    for folder in unit:

        for file in sorted(glob.glob(os.path.join(folder, '*SHORT.csv'))):
        # if n%2 == 0:
            print(n)
            whole_df = pd.read_table(file,
                                 delimiter=',',
                                 header=None,
                                 names=range(0, 2),
                                 )
            whole_df.columns = ['position', 'load']


            first_pull = whole_df[ind(whole_df, 'FIRST PULL') + 1:ind(whole_df, 'FIRST RETRACT')].astype(float).reset_index(drop=True)
            first_retract = whole_df[ind(whole_df, 'FIRST RETRACT') + 1:ind(whole_df, 'SECOND PULL')].astype(float).iloc[::-1].reset_index(drop=True)
            second_pull = whole_df[ind(whole_df, 'SECOND PULL') + 1:ind(whole_df, 'SECOND RETRACT')].astype(float).reset_index(drop=True)
            second_retract = whole_df[ind(whole_df, 'SECOND RETRACT') + 1:ind(whole_df, 'THIRD PULL')].astype(float).iloc[::-1].reset_index(drop=True)
            third_pull = whole_df[ind(whole_df, 'THIRD PULL') + 1:ind(whole_df, 'THIRD RETRACT')].astype(float).reset_index(drop=True)
            third_retract = whole_df[ind(whole_df, 'THIRD RETRACT') + 1:ind(whole_df, 'FOURTH PULL')].astype(float).iloc[::-1].reset_index(drop=True)
            fourth_pull = whole_df[ind(whole_df, 'FOURTH PULL') + 1:ind(whole_df, 'FOURTH RETRACT')].astype(float).reset_index(drop=True)
            fourth_retract = whole_df[ind(whole_df, 'FOURTH RETRACT') + 1:ind(whole_df, 'FIFTH PULL')].astype(float).iloc[::-1].reset_index(drop=True)
            fifth_pull = whole_df[ind(whole_df, 'FIFTH PULL') + 1:ind(whole_df, 'FIFTH RETRACT')].astype(float).reset_index(drop=True)
            fifth_retract = whole_df[ind(whole_df, 'FIFTH RETRACT') + 1:ind(whole_df, 'SIXTH PULL')].astype(float).iloc[::-1].reset_index(drop=True)
            sixth_pull = whole_df[ind(whole_df, 'SIXTH PULL') + 1:ind(whole_df, 'SIXTH RETRACT')].astype(float).reset_index(drop=True)
            sixth_retract = whole_df[ind(whole_df, 'SIXTH RETRACT') + 1:].astype(float).iloc[::-1].reset_index(drop=True)

            ori_p = ten_p(first_pull)
            print(ori_p)
            n += 1

            pull = [first_pull, second_pull, third_pull, fourth_pull, fifth_pull, sixth_pull]
            retract = [first_retract, second_retract, third_retract, fourth_retract, fifth_retract, sixth_retract]


            for i in range(len(pull)):
                norm(pull[i])
                norm(retract[i])
            # loop for seuqntial Instron
            for i in range(len(pull)):
                pull[i].drop(pull[i].loc[0:FindZeroStrain(pull[i])].index, inplace=True)

                # retract[i].reset_index(inplace=True)
                rmindex = []
                purify_revrs(retract[i])

                retract[i].drop(retract[i].loc[rmindex].index, inplace=True)
                retract[i] = retract[i].reset_index(drop=True)
                pull[i] = pull[i].reset_index(drop=True)

            sm1pull = pd.DataFrame(data=smooth(first_pull, smf), columns=['strain', 'stress'])
            sm2pull = pd.DataFrame(data=smooth(second_pull, smf), columns=['strain', 'stress'])
            sm3pull = pd.DataFrame(data=smooth(third_pull, smf), columns=['strain', 'stress'])
            sm4pull = pd.DataFrame(data=smooth(fourth_pull, smf), columns=['strain', 'stress'])
            sm5pull = pd.DataFrame(data=smooth(fifth_pull, smf), columns=['strain', 'stress'])
            sm6pull = pd.DataFrame(data=smooth(sixth_pull, smf), columns=['strain', 'stress'])

            sm1retract = pd.DataFrame(data=smooth(first_retract, smf), columns=['strain', 'stress'])
            sm2retract = pd.DataFrame(data=smooth(second_retract, smf), columns=['strain', 'stress'])
            sm3retract = pd.DataFrame(data=smooth(third_retract, smf), columns=['strain', 'stress'])
            sm4retract = pd.DataFrame(data=smooth(fourth_retract, smf), columns=['strain', 'stress'])
            sm5retract = pd.DataFrame(data=smooth(fifth_retract, smf), columns=['strain', 'stress'])
            sm6retract = pd.DataFrame(data=smooth(sixth_retract, smf), columns=['strain', 'stress'])

            smpull = [sm1pull, sm2pull, sm3pull, sm4pull, sm5pull, sm6pull]
            smretract = [sm1retract, sm2retract, sm3retract, sm4retract, sm5retract, sm6retract]



            for i in range(len(smpull)):
                smpull[i]['load'] = pull[i].load
                smretract[i]['load'] = retract[i].load

            # incremental plastic strain energy density
            inc_p = []
            for i in range(len(smpull)-1):
                if i == len(smpull)-2:
                    inc_p.append(area_t(smpull[i]) - area_t(smpull[i+1]))
                else:
                    inc_p.append(area_t(smpull[i]) - area_s(smpull[i+1], smpull[i]))

            p_inc_energ.loc[len(p_inc_energ)] = inc_p

            # for sequential Instron - plastic strain energy density
            pla_energ = []
            ener = 0
            for i in range(len(inc_p)):
                ener += inc_p[i]
                pla_energ.append(ener)
            p_energ.loc[len(p_energ)] = pla_energ

            # for sequential Instron - elastic strain energy density
            ela_energ = []
            for i in range(len(smpull)-1):
                if i == len(smpull)-2:
                    ela_energ.append(area_t(smpull[i+1]))
                else:
                    ela_energ.append(area_s(smpull[i+1], smpull[i]))
            e_energ.loc[len(e_energ)] = ela_energ

            # for sequential Instron - total strain energy density
            total_energ = []
            for i in range(len(ela_energ)):
                total_energ.append(ela_energ[i] + pla_energ[i])

            # stored energy
            retra = []
            for i in range(len(smretract)-1):
                retra.append(area_t(smretract[i]))
            store_energ.loc[len(store_energ)] = retra

            # hysteresis dissipated energy
            hys = []
            for i in range(len(strainT)):
                hys.append(ela_energ[i] - retra[i])
            hyst_energy.loc[len(hyst_energy)] = hys

        # add data of this file to the results data frame
            t_energ.loc[len(t_energ)] = total_energ

            # plastic hysteresis ratio calculation
            p_hys_list = []
            for i in range(len(strainT)):
                p_hys_list.append(pla_energ[i]/total_energ[i] * 100)
            p_hys.loc[len(p_hys)] = p_hys_list

            # plastic hysteresis ratio calculation
            e_hys_list = []
            for i in range(len(strainT)):
                e_hys_list.append(hys[i]/total_energ[i] * 100)
            e_hys.loc[len(e_hys)] = e_hys_list

            sto_ratio = []
            for i in range(len(strainT)):
                sto_ratio.append(retra[i]/total_energ[i] * 100)
            sto_ratio_DF.loc[len(sto_ratio_DF)] = sto_ratio

    C_t = []
    C_p = []
    C_e = []
    C_s = []
    C_stress = []
    Ce_t = []
    Ce_p = []
    Ce_e = []
    Ce_s = []
    Ce_stress = []
    for i in range(len(strainT)):
        C_p.append(np.mean(p_hys.iloc[:, i]))
        C_e.append(np.mean(e_hys.iloc[:, i]))
        C_s.append(np.mean(sto_ratio_DF.iloc[:, i]))
        Ce_p.append(np.std(p_hys.iloc[:, i]))
        Ce_e.append(np.std(e_hys.iloc[:, i]))
        Ce_s.append(np.std(sto_ratio_DF.iloc[:, i]))


    total_chart = pd.DataFrame(columns=strainT)
    total_chart = p_hys
    total_chart.loc[len(total_chart)] = C_p
    total_chart.loc[len(total_chart)] = Ce_p
    total_chart.to_csv(''f'{output}/' + name + '_Elastic_Energy_ratio_plasticity.csv')

    total_chart = pd.DataFrame(columns=strainT)
    total_chart = e_hys
    total_chart.loc[len(total_chart)] = C_e
    total_chart.loc[len(total_chart)] = Ce_e
    total_chart.to_csv(''f'{output}/' + name + '_Elastic_Energy_ratio_elastic hysteresis.csv')

    total_chart = pd.DataFrame(columns=strainT)
    total_chart= sto_ratio_DF
    total_chart.loc[len(total_chart)] = C_s
    total_chart.loc[len(total_chart)] = Ce_s
    total_chart.to_csv(''f'{output}/' + name + '_Elastic_Energy_ratio_stored.csv')

    lw = 2
    ms = 6
# plot with error bar
    allx = plt.subplot(111)
    allx.set(xlim = [0, 28], ylim = [0, 80], ylabel = 'Percentage (%)', xlabel = 'Strain (%)')
    x = [0, 5, 10, 15, 20, 25]
    allx.set_xticks(x)
    allx.set_xticklabels(x)

    # ax.errorbar(strainT, C_stress, yerr = Ce_stress, marker = '8',color= color[g], alpha = 0.6, markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)
    allx.errorbar(strainT, C_e, yerr = Ce_e, linestyle = ls[g], marker = '8', color= 'C3', markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)
    allx.errorbar(strainT, C_p, yerr = Ce_p,  linestyle = ls[g], marker = '8', color='C9', markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)
    allx.errorbar(strainT, C_s, yerr = Ce_s,  linestyle = ls[g], marker = '8', color='C2', markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw)



    blue_line = mlines.Line2D([], [], color='C2', label='Stored energy')
    green_line = mlines.Line2D([], [], color='C3', label='Dissipated energy (elastic hysteresis)')
    orange_line = mlines.Line2D([], [], color='C9', label='Dissipated energy(plasticity)')
    allx.legend(handles=[blue_line, orange_line, green_line], loc=2, fontsize=10, frameon=False)

    g += 1

plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().set_size_inches(6.4, 5.5)
plt.show()

