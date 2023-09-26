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
from math import e
from numpy import trapz
import statsmodels.api as sm
import copy
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
fig = mpl.pyplot.gcf()
output = '/Users/jingyiyu/Documents/Cosgrovelab/manuscript/CW_mechanics_Instron/Data_collection/Python_plot/'

# User input
smf = 0.04

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

# for load sequence of 4-8-12-16-20
loadx = [5, 10, 15, 20, 25]
# for load sequence of 8-12-16-20-24
# loadx = [8, 12, 16, 20, 24]


# define a function that calls the index of certain value
def ind(df, value):
    for i in range(len(df)):
        if df.position[i] == value:
            return i


# define a function that find the index at certain loads
def ten_ind(pull, load):
    if load == 'N':
        return len(pull)
    else:
        for i in range(len(pull) - 5):
            if (pull.load[i] > load) & (pull.load[i + 1] > load) & (pull.load[i + 2] > load) & \
                    (pull.load[i + 3] > load) & (pull.load[i + 4] > load):
                index = i
                return index

def strain_ind(pull, strain):
    if strain == 'N':
        return len(pull)
    else:
        for i in range(len(pull) - 5):
            if (pull.strain[i] > strain) & (pull.strain[i + 1] > strain) & (pull.strain[i + 2] > strain) & \
                    (pull.strain[i + 3] > strain) & (pull.strain[i + 4] > strain):
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

# define a function that smooth the curve
def smooth(pull,f):
    lowess = sm.nonparametric.lowess
    z = lowess(pull.stress, pull.strain,frac= f)
    return z

# define a function that calculate the stress and strain of the curve
def norm(pull):
    # strain calculated based on very original length (before plastic deformation)
    ori_p = ten_p(first_pull)
    pull['force_N'] = pull.load * 0.0098

    # engineering stress & strain
    pull['stress'] = pull.force_N / (thickness * width * 0.001)
    pull['strain'] = (pull.position - ori_p) / (5000 + ori_p)

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
    # fit using stress (MPa)
    z = np.polyfit(fitting_part.strain, fitting_part.stress, 1)
    pull['fit_e'] = pull.strain * z[0] + z[1]
    modulus = z[0]
    return modulus

def fit_expo(pull, weight):
    ## fit accroding to original data
    z = np.polyfit(pull.strain, np.log(pull.stress), 1, w = weight)
    pull['fit_e'] = e**z[1] * e**(z[0]*pull.strain)
    equation = [z[0], z[1]]
    return equation

def fit_expo_l(pull, weight):
    ## fit accroding to original data
    z = np.polyfit(pull.strain, np.log(pull.stress), 1, w = weight)
    pull['fit_e'] = e**z[1] * e**(z[0]*pull.strain)

    equation = [z[0], z[1]]
    return equation


def fit_linear(pull, weight):
    z = np.polyfit(pull.strain, pull.stress, 1, w = weight)
    pull['fit_e'] = z[0]* pull.strain + z[1]

    equation = [z[0], z[1]]
    return equation

def tar_ind(pull, load):
    if load == 'N':
        return len(pull)
    else:
        stress = load * 0.0098/ (thickness * width * 0.001)
        for i in range(len(pull) - 5):
            if (pull.stress[i] > stress) & (pull.stress[i + 1] > stress) & (pull.stress[i + 2] > stress) & \
                    (pull.stress[i + 3] > stress) & (pull.stress[i + 4] > stress):
                return i

def Rsquare(curve):
    mean = curve.stress.mean()
    sum_total = 0
    sum_residue = 0
    r = 0
    for i in range(len(curve)):
        sum_total += (curve.stress[i] - mean)**2
        sum_residue += (curve.stress[i] - curve.fit_e[i])**2
    r = 1 - sum_residue/sum_total
    return r



# extract data from each pull
# Use if multiple files need analysis
g = 0
mark = '8'
total = pd.DataFrame(columns=loadx)
ela = pd.DataFrame(columns=loadx)
retra = pd.DataFrame(columns=loadx)

loadx_lexpo = [5, 10, 15, 20, 25]
xlabel = [0, 5, 10, 15, 20, 25]

## plot the elasic loading fitting curve
sa = plt.subplot(111)
sa.set(xlim = [0, 30], ylim = [-0.5, 80],xlabel='Strain (%)', ylabel='Linear slope')
sa.set_xticks(xlabel)
sa.set_xticklabels(xlabel)
say2 = sa.twinx()
say2.set(ylim = [0, 85])
say2.set_ylabel('Loading exponent', rotation=270, va='bottom')

## plot the elasic unloading fitting curve
# sa2 = plt.subplot(111)
# sa2.set(xlim = [0, 30], ylim = [0, 180])
# sa2.set_xticks(xlabel)
# sa2.set_xticklabels(xlabel)

for unit in group:
    if g == 0:
        style = '-'
    elif g == 1:
        style = '--'
#     elif g == 2:
#         mark = 'p'
#     elif g == 3:
#         mark = 'P'
#
    n = 0
    color = ['C0','C1','C2', 'C3', 'C4', 'C5','C6','C8','C9','C10','C11', 'C12']
    retract_exp_0 = pd.DataFrame(columns=loadx)
    retract_exp_1 = pd.DataFrame(columns=loadx)
    load_exp_0 = pd.DataFrame(columns=loadx)
    load_exp_1 = pd.DataFrame(columns=loadx)
    load_lnr_0 = pd.DataFrame(columns=loadx)
    load_lnr_1 = pd.DataFrame(columns=loadx)
    f = 0
    markj = 0
    for folder in unit:
        for file in sorted(glob.glob(os.path.join(folder, '*SHORT.csv'))):
            # print(n)
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

            pull = [first_pull, second_pull, third_pull, fourth_pull, fifth_pull, sixth_pull]
            retract = [first_retract, second_retract, third_retract, fourth_retract, fifth_retract, sixth_retract]

            # calculate the stress and strain of the curve
            for i in range(len(pull)):
                norm(pull[i])
                norm(retract[i])

            target = [0.05, 0.1, 0.15, 0.20, 0.25, 0.25]
            ela_target = ['N', 0.05, 0.10, 0.15, 0.20, 'N']
            ela1_epon = first_pull
            ela2_epon = second_pull
            ela3_epon = third_pull
            ela4_epon = fourth_pull
            ela5_epon = fifth_pull
            ela6_epon = sixth_pull

            ela1_lnr = first_pull
            ela2_lnr = second_pull
            ela3_lnr = third_pull
            ela4_lnr = fourth_pull
            ela5_lnr = fifth_pull
            ela6_lnr = sixth_pull

            re_fit1 = first_retract
            re_fit2 = second_retract
            re_fit3 = third_retract
            re_fit4 = fourth_retract
            re_fit5 = fifth_retract
            re_fit6 = sixth_retract

            orip1 = copy.deepcopy(first_pull)
            orip2 = copy.deepcopy(second_pull)
            orip3 = copy.deepcopy(third_pull)
            orip4 = copy.deepcopy(fourth_pull)
            orip5 = copy.deepcopy(fifth_pull)
            orip6 = copy.deepcopy(sixth_pull)


            orir1 = copy.deepcopy(first_retract)
            orir2 = copy.deepcopy(second_retract)
            orir3 = copy.deepcopy(third_retract)
            orir4 = copy.deepcopy(fourth_retract)
            orir5 = copy.deepcopy(fifth_retract)
            orir6 = copy.deepcopy(sixth_retract)

            elaexpo = [ela1_epon, ela2_epon,ela3_epon, ela4_epon, ela5_epon, ela5_epon]
            elalnr = [ela1_lnr, ela2_lnr, ela3_lnr, ela4_lnr, ela5_lnr, ela6_lnr]
            retract_fit = [re_fit1, re_fit2, re_fit3, re_fit4, re_fit5, re_fit6]
            orip = [orip1, orip2, orip3, orip4, orip5, orip6]
            orir = [orir1, orir2, orir3, orir4, orir5, orir6]

        # loop for seuqntial Instron
            for i in range(len(pull)):

                retract[i].drop(retract[i].loc[0:ten_ind(retract[i], 0.1)].index, inplace=True)
                retract[i].reset_index(inplace=True)

                pull[i].drop(pull[i].loc[0:ten_ind(pull[i], 0.1)].index, inplace = True)
                pull[i].reset_index(inplace=True)

                #extract elastic part of the curve
                pull[i].drop(pull[i].loc[strain_ind(pull[i],ela_target[i]):].index, inplace = True)
                pull[i].reset_index(inplace=True)

                elalnr[i] = pull[i].loc[int(len(pull[i])*0.1):int(len(pull[i])*0.35)]

                elaexpo[i] = pull[i].loc[int(len(pull[i])*0.35):int(len(pull[i])*0.88)]

                retract[i].reset_index(inplace=True)
                rmindex = []
                retract_fit[i] = retract[i].loc[int(len(retract[i])*0.2):int(len(retract[i])*0.95)]

                retract_fit[i] = retract_fit[i].reset_index(drop=True)
                elaexpo[i] = elaexpo[i].reset_index(drop=True)
                elalnr[i] = elalnr[i].reset_index(drop=True)
                pull[i] = pull[i].reset_index(drop=True)

            sm1retract = pd.DataFrame(data=smooth(first_retract,smf), columns=['strain','stress'])
            sm2retract = pd.DataFrame(data=smooth(second_retract,smf), columns=['strain','stress'])
            sm3retract = pd.DataFrame(data = smooth(third_retract,smf),columns=['strain','stress'])
            sm4retract = pd.DataFrame(data = smooth(fourth_retract,smf),columns=['strain','stress'])
            sm5retract = pd.DataFrame(data = smooth(fifth_retract,smf),columns=['strain','stress'])



            smretract = [sm1retract, sm2retract, sm3retract, sm4retract, sm5retract]

        ####### plot the stress-train curve and fitting curve

            for j in range(len(retract_fit)-1):
                weight_re = []
                for i in range(len(retract_fit[j])):
                    if i <= int(len(retract_fit[j]) * 0.3):
                        weight_re.append(1)
                    elif int(len(retract_fit[j]) * 0.3) < i <= int(len(retract_fit[j]) * 0.95):
                        weight_re.append(1)
                    else:
                        weight_re.append(1)

                fitNmnow = fit_expo(retract_fit[j], weight_re)
                retract_fit_line = mlines.Line2D([], [], color = 'grey', label = 'Exponential fitting curve')
                retract_exp_0.loc[f,loadx[j]] = fitNmnow[0]
                retract_exp_1.loc[f,loadx[j]] = fitNmnow[1]



            for j in range(1, len(loadx)+1):

                weight_re = []
                for i in range(len(elalnr[j])):
                    if i < int(len(elalnr[j]) * 0.05):
                        weight_re.append(1)
                    else:
                        weight_re.append(1)

                fitl = (fit_linear(elalnr[j], weight_re))
                load_lnr_0.loc[f,loadx[j-1]] = fitl[0]
                load_lnr_1.loc[f,loadx[j-1]] = fitl[1]
                load_expo_line = mlines.Line2D([], [], color = 'grey', label = 'Exponential fitting curve')
                load_lnr_line = mlines.Line2D([], [], color = 'black', label = 'Linear fitting curve')

            for j in range(1, len(loadx)+1):
                weight_re = []
                markj = j
                for i in range(len(elaexpo[j])):
                    if i > int(len(elaexpo[j]) * 0.9):
                        weight_re.append(1)
                    else:
                        weight_re.append(1)

                fitNmnow = (fit_expo_l(elaexpo[j], weight_re))
                load_exp_0.loc[f,loadx[j-1]] = fitNmnow[0]
                load_exp_1.loc[f,loadx[j-1]] = fitNmnow[1]
                print(load_exp_0)

            # plt.show()
            fig.set_size_inches(8, 6)

            master = plt.subplot(111)
            lw = 2
            for i in range(len(retract)):
                master.plot(orir[i].strain *100, orir[i].stress, color = color[i],linestyle = '--', linewidth = lw)
                master.plot(orip[i].strain *100, orip[i].stress, color = color[i], linewidth = lw)
            for i in range(len(retract_fit)-1):
                master.plot(retract_fit[i].strain *100, retract_fit[i].fit_e, color = 'grey', linewidth = lw)
            for i in range(1,len(elalnr)):
                master.plot(elalnr[i].strain *100, elalnr[i].fit_e, color = 'black', linewidth = lw)
            for i in range(2,len(elaexpo)):
                master.plot(elaexpo[i].strain *100, elaexpo[i].fit_e, color = 'grey', linewidth = lw)
            # master.set(xlabel = 'Strain(%)', ylabel = 'Stress(MPa)', title = 'Fitting of elastic curves')
            # master.grid(alpha=0.4, linestyle='--')
            load_expo_line = mlines.Line2D([], [], color = 'grey', label = 'Exponential fitting curve')
            load_lnr_line = mlines.Line2D([], [], color = 'black', label = 'Linear fitting curve')
            # master.legend(handles = [load_lnr_line, load_expo_line], loc = 2, fontsize=10, frameon = False)

            print(f)
            f += 1
            # plt.gcf().subplots_adjust(bottom=0.12, right = 0.8)
            # plt.savefig(''f'{output}/Seq_cyclic_fitting_master.pdf', transparent = True)

            # plt.show()

    R_0 = []
    L_l0 = []
    L_e0 = []
    R_0ero = []
    L_l0ero = []
    L_e0ero = []
    for i in range(len(loadx)):
        R_0.append(np.mean(retract_exp_0.iloc[:, i]))
        L_l0.append(np.mean(load_lnr_0.iloc[:, i]))

        R_0ero.append(np.std(retract_exp_0.iloc[:, i]))
        L_l0ero.append(np.std(load_lnr_0.iloc[:, i]))

        L_e0.append(np.mean(load_exp_0.iloc[:, i]))
        L_e0ero.append(np.std(load_exp_0.iloc[:, i]))


    ms = 4
    lw = 2
    mark = '8'
    # sa2.errorbar(loadx, R_0, yerr = R_0ero, color= 'C7', alpha = 0.6, marker = mark, markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw, linestyle = style)
    sa.errorbar(loadx, L_l0, yerr = L_l0ero, color='black', marker = mark, markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw, linestyle = style)
    say2.errorbar(loadx, L_e0, yerr = L_e0ero, color='gray', marker = mark, markersize = ms, linewidth = lw, capsize=5, markeredgewidth= lw, linestyle = style, alpha = 0.8)

    t_line = mlines.Line2D([], [], color = 'C7', label = 'Retract expo base')
    orange_line = mlines.Line2D([], [], color = 'C9', label = 'Load linear slope')
    blue_line = mlines.Line2D([], [], color = 'C3', label = 'Load expo base')

    g += 1

    # plt.show()

plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().set_size_inches(6, 5.5)
# plt.savefig(''f'{output}/PEG_loading_fit.pdf', transparent = True)

plt.show()

