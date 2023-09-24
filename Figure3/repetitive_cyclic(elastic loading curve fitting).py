import pandas as pd
import glob
import os
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.lines as mlines
import statsmodels.api as sm
from math import e
import matplotlib as mpl
mpl.use('macosx')
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Avenir'
plt.rcParams['font.size'] = 15
plt.rcParams['axes.linewidth'] = 2
mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['ytick.major.size'] = 5
mpl.rcParams['ytick.major.width'] = 2
fig = mpl.pyplot.gcf()
# fig.set_size_inches(6.4, 9.6)
output = '/Users/jingyiyu/Documents/Cosgrovelab/manuscript/CW_mechanics_Instron/Data_collection/Python_plot/'

lw = 2
thickness = 7
width = 3
smf = 0.04
intv = 5
l = -0.5
r = 45

repetitive10 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/repetitive pulling_8.1.20/10g'
repetitive10_2 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/repetitive_pulling_10g_3mpm_10.6.20'

# define a function that calls the index of certain value
def ind(df, value):
    for i in range(len(df)):
        if df.position[i] == value:
            return i

def ten_p(pull):
    for i in range(len(pull) - 5):
        if (pull.load[i] > 0.1) & (pull.load[i + 1] > 0.1) & (pull.load[i + 2] > 0.1) & \
                (pull.load[i + 3] > 0.1) & (pull.load[i + 4] > 0.1):
            point = pull.position[i]
            return point

def purify_p(pull, target):
    for i in range(len(pull)):
        if i < (len(pull)-1) :
            if pull.load[i] > target:
                return i
                break
        else:
            return (len(pull) -1)


# define a function that calculate the stress and strain of the curve
def norm(pull):
    pull['force_N'] = pull.load * 0.0098

    # engineering stress & strain
    pull['stress'] = pull.force_N / (thickness * width * 0.001)
    pull['strain'] = (pull.position - ori_p) / (5000 + ori_p)

def purify_revrs(pull):
    b = np.array(pull.position)
    maxp_index = np.argmax(b)
    threshload = pull.load[maxp_index]
    for i in range(len(pull)):
        if pull.load[i] > threshload:
            rmindex.append(i)


# define a function that smooth the curve
def smooth(pull,f):
    lowess = sm.nonparametric.lowess
    z = lowess(pull.stress, pull.strain,frac= f)
    return z

def derv(pull, inter):
    deriv = []
    for i in range(inter):
        deriv.append((pull.stress[i + inter] - pull.stress[i]) / (pull.strain[i + inter] - pull.strain[i]))
    for i in range(inter, len(pull) - inter):
        deriv.append(
            (pull.stress[i + inter] - pull.stress[i - inter]) / (pull.strain[i + inter] - pull.strain[i - inter]))
    for i in range(len(pull) - inter, len(pull)):
        deriv.append((pull.stress[i] - pull.stress[i - inter]) / (pull.strain[i] - pull.strain[i - inter]))
    return deriv

def dervuni(y, x, inter):
    deriv = []
    for i in range(inter):
        deriv.append((y[i + inter] - y[i]) / (x[i + inter] - x[i]))
    for i in range(inter, len(y) - inter):
        deriv.append(
            (y[i + inter] - y[i - inter]) / (x[i + inter] - x[i - inter]))
    for i in range(len(y) - inter, len(y)):
        deriv.append((y[i] - y[i - inter]) / (x[i] - x[i - inter]))
    return deriv

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

def fit_expo(pull, weight):
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

color = ['blue','green', 'orange', 'red', 'purple', 'black']
n = 0
markj = 0
fitNm = ['z[0]', 'z[1]']
equationlist = pd.DataFrame(columns= fitNm)
Lln0_lst = []
Lexpo0_lst = []
Rexpo0_lst = []
Lln1_lst = []
Lexpo1_lst = []
Rexpo1_lst = []

R2_Lln_lst = []
R2_Lexpo_lst = []
R2_Rexpo_lst = []

for file in sorted(glob.glob(os.path.join(repetitive10_2, '*SHORT.csv'))):
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


    pull = [second_pull, third_pull, fourth_pull, fifth_pull, sixth_pull]
    retract = [second_retract, third_retract, fourth_retract, fifth_retract, sixth_retract]
    pull = [second_pull]
    retract = [second_retract]

    ori_p = ten_p(first_pull)

    ### normalize curve
    for i in range(len(pull)):
        norm(pull[i])
        norm(retract[i])

    # ela1_epon = first_pull
    ela2_epon = second_pull
    ela3_epon = third_pull
    ela4_epon = fourth_pull
    ela5_epon = fifth_pull
    ela6_epon = sixth_pull

    # ela1_lnr = first_pull
    ela2_lnr = second_pull
    ela3_lnr = third_pull
    ela4_lnr = fourth_pull
    ela5_lnr = fifth_pull
    ela6_lnr = sixth_pull

    re_1 = first_retract
    re_2 = second_retract
    re_3 = third_retract
    re_4 = fourth_retract
    re_5 = fifth_retract

    elaexpo = [ela2_epon]
    elalnr = [ela2_lnr]
    elaret = [re_2]



    for i in range(len(elaexpo)):
        pull[i].drop(pull[i].loc[0:ten_ind(pull[i], 0.1)].index, inplace = True)
        pull[i].reset_index(inplace=True)

        print(retract[i])
        retract[i].drop(retract[i].loc[0:ten_ind(retract[i], 0.1)].index, inplace=True)

        retract[i].reset_index(inplace=True)
        rmindex = []
        purify_revrs(retract[i])
        retract[i].drop(retract[i].loc[rmindex].index, inplace=True)

        retract[i].reset_index(inplace=True)
        elaret[i] = retract[i].loc[int(len(retract[i])*0.2):int(len(retract[i])*0.95)]

#extract elastic part of the curve

        elaexpo[i] = pull[i].loc[int(len(pull[i])*0.35):int(len(pull[i])*0.95)]
        # elalnr[i] = pull[i].drop(pull[i].loc[:].index)
        elalnr[i] = pull[i].loc[int(len(pull[i])*0.1):int(len(pull[i])*0.35)]

        elaexpo[i] = elaexpo[i].reset_index(drop=True)
        elalnr[i] = elalnr[i].reset_index(drop=True)
        elaret[i] = elaret[i].reset_index(drop=True)
        pull[i] = pull[i].reset_index(drop=True)


# smpull = [sm1pull, sm2pull]
    # smretract = [sm1retract, sm2retract]

    #new_data = pd.concat([first_pull, first_retract, second_pull, second_retract], axis=1)

    #summary = pd.concat([summary, new_data],axis = 1)
    map = 'Paired'
    ax = plt.subplot(111)

    # bx = plt.subplot(312)
    # cx = plt.subplot(313)

    # ax.grid(alpha=0.4, linestyle='--')
    # bx.grid(alpha=0.4, linestyle='--')
    # cx.grid(alpha=0.4, linestyle='--')

# for i in range(len(pull)):
# for only first two cycle
    load_lnr0 = []
    load_lnr1 = []

    for j in range(len(elaexpo)):
        markj = j
        weight_re = []
        for i in range(len(elalnr[j])):
            if i < int(len(elalnr[j]) * 0.05):
                weight_re.append(1)
            else:
                weight_re.append(1)

        fitNmnow = (fit_linear(elalnr[j], weight_re))
        print('elanr0 = ' + str(fitNmnow[0]))
        print('elanr1 = ' + str(fitNmnow[1]))
        elanr0 = fitNmnow[0]
        elanr1 = fitNmnow[1]
        load_lnr0.append(fitNmnow[0])
        load_lnr1.append(fitNmnow[1])
# print(fitNmnow)
        ax.plot(pull[j].strain * 100, pull[j].stress, color = 'red', linewidth = lw)
        ax.plot(elalnr[j].strain * 100, elalnr[j].fit_e, color = 'black', linewidth = lw)
        # bx.scatter(elalnr[j].strain * 100, elalnr[j].stress - elalnr[j].fit_e,color = color[n], s=2)
        # bx.set(title = 'Residue of linear fitting', xlabel = 'Strain(%)', ylabel = 'Residue')
        load_expo_line = mlines.Line2D([], [], color = 'grey', label = 'Exponential fitting curve', linewidth = lw)
        load_lnr_line = mlines.Line2D([], [], color = 'black', label = 'Linear fitting curve', linewidth = lw)
        # ax.legend(handles = [load_lnr_line, load_expo_line], loc = 2, fontsize = 8)

        # unsmoothed data
        weight_re = []
        for i in range(len(elaexpo[j])):
            if i > int(len(elaexpo[j])):
                weight_re.append(0)
            else:
                weight_re.append(1)

        fitNmnow = (fit_expo(elaexpo[j], weight_re))
        print('elaexpo0 = ' + str(fitNmnow[0]))
        print('elaexpo1 = ' + str(fitNmnow[1]))
        elaexpo0 = fitNmnow[0]
        elaexpo1 = fitNmnow[1]
        ax.set(xlabel = 'Strain(%)', ylabel = 'Stress(MPa)')
        ax.plot(elaexpo[j].strain * 100, elaexpo[j].fit_e, color = 'grey', linewidth = lw)

        weight_re = []
        for i in range(len(elaret[j])):
            if i <= int(len(elaret[j]) * 0.2):
                weight_re.append(1)
            elif int(len(elaret[j]) * 0.2) < i < int(len(elaret[j]) * 0.95):
                weight_re.append(1)
            else:
                weight_re.append(1)

        fitNmnow = fit_expo(elaret[j], weight_re)
        ax.plot(retract[j].strain * 100, retract[j].stress, color = 'red', linestyle = '--', linewidth = lw)
        expt_line = mlines.Line2D([], [], color = 'red', label = 'Elastic cyclic loading')
        linear_fit_line = mlines.Line2D([], [], color = 'black', label = 'Linear fitting')
        retract_fit_line = mlines.Line2D([], [], color = 'grey', label = 'Exponential fitting', linewidth = lw)
        ax.legend(handles = [expt_line, linear_fit_line, retract_fit_line], loc = 2, fontsize = 8, frameon = False)
        ax.plot(elaret[j].strain * 100, elaret[j].fit_e, color = 'grey', linewidth = lw)
        ax.set(xlim = [0,11])
        plt.gcf().subplots_adjust(bottom=0.2, wspace= 0.4)
        print('retrexpo0 = ' + str(fitNmnow[0]))
        print('retrexpo1 = ' + str(fitNmnow[1]))
        retrexpo0 = fitNmnow[0]
        retrexpo1 = fitNmnow[1]


    Lln0_lst.append(elanr0)
    Lexpo0_lst.append(elaexpo0)
    Rexpo0_lst.append(retrexpo0)
    Lln1_lst.append(elanr1)
    Lexpo1_lst.append(elaexpo1)
    Rexpo1_lst.append(retrexpo1)
    R2_Lln_lst.append(Rsquare(elalnr[0]))
    R2_Lexpo_lst.append(Rsquare(elaexpo[0]))
    R2_Rexpo_lst.append(Rsquare(elaret[0]))

    fig.set_size_inches(6.4, 5.5)
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.savefig(''f'{output}/Repetitive_cyclic_fitting_partI.pdf', transparent=True)
    plt.show()

    bx = plt.subplot(211)
    cx = plt.subplot(212)
    sm2pull = pd.DataFrame(data=smooth(second_pull, smf), columns=['strain', 'stress'])
    sm2retract = pd.DataFrame(data=smooth(second_retract, smf), columns=['strain', 'stress'])
    bx.plot(sm2pull.strain * 100, derv(sm2pull, 10))
    bx.plot(sm2retract.strain * 100, derv(sm2retract, 10))

    cx.plot(sm2pull.strain * 100, derv(sm2pull, 10))
    cx.plot(sm2retract.strain * 100, derv(sm2retract, 10))
    bx.set(xlabel='Strain (%)', ylabel='Modulus (MPa)')
    cx.set(xlabel='Strain (%)', ylabel='Modulus (MPa)')
    cx.set_yscale('log')

    fig.set_size_inches(6.4, 4.8)
    plt.gcf().subplots_adjust(bottom=0.2, hspace = 0.5, right = 0.8)
    # plt.savefig(''f'{output}/Repetitive_cyclic_fitting_modulus_map.pdf', transparent=True)
    plt.show()

    nx = plt.subplot(111)
    nx.plot(sm2retract.stress, derv(sm2retract, 10))
    plt.show()

R2_Lln_mean = np.mean(R2_Lln_lst)
R2_Lexpo_mean = np.mean(R2_Lexpo_lst)
R2_Rexpo_mean = np.mean(R2_Rexpo_lst)
Lln0_mean = np.mean(Lln0_lst)
Lexpo0_mean = np.mean(Lexpo0_lst)
Rexpo0_mean = np.mean(Rexpo0_lst)
Lln1_mean = np.mean(Lln1_lst)
Lexpo1_mean = np.mean(Lexpo1_lst)
Rexpo1_mean = np.mean(Rexpo1_lst)

R2_Lln_std = np.std(R2_Lln_lst)
R2_Lexpo_std = np.std(R2_Lexpo_lst)
R2_Rexpo_std = np.std(R2_Rexpo_lst)
Lln0_std = np.std(Lln0_lst)
Lexpo0_std = np.std(Lexpo0_lst)
Rexpo0_std = np.std(Rexpo0_lst)
Lln1_std = np.std(Lln1_lst)
Lexpo1_std = np.std(Lexpo1_lst)
Rexpo1_std = np.std(Rexpo1_lst)

print(Lln0_mean, Lln1_mean, Lexpo0_mean, Lexpo1_mean, Rexpo0_mean, Rexpo1_mean)
print(R2_Lln_mean, R2_Lln_std, R2_Lexpo_mean, R2_Lexpo_std, R2_Rexpo_mean, R2_Rexpo_std)




