import pandas as pd
import glob
import os
import numpy as np
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
fig = mpl.pyplot.gcf()
fig.set_size_inches(6.4, 4.8)

output = '/Users/jingyiyu/Documents/Cosgrovelab/manuscript/CW_mechanics_Instron/Data_collection/Python_plot/'

thickness = 7
width = 3
target = 15
inter = 5
Title = ''

# onion1
# R10 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/Onion1_0.1_1_10/10mm_min'
# R1 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/Onion1_0.1_1_10/1mm_min'
# R01 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/Onion1_0.1_1_10/0.1mm_min'

#onion2
# R10 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/Onion2_0.1_1_10/10mm_min'
# R1 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/Onion2_0.1_1_10/1mm_min'
# R01 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/Onion2_0.1_1_10/0.1mm_min'

#all onion
#26C
R10 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/ForAnalysis/26C/10'
R1 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/ForAnalysis/26C/1'
R01 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/ForAnalysis/26C/0.1'
#4C
# R10 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/ForAnalysis/4C/10'
# R1 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/ForAnalysis/4C/1'
# R01 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/ForAnalysis/4C/0.1'
## RT control
# R01_4c = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/2nd_set/26C/Onion3/2d_inFridge_0.1'

# fresh = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Dif scale/Onion1/5th'
# nextday = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Dif scale/Onion1/5th/the next day test'


### file output directory
# fR10 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/ForAnalysis/26C-10.csv'
# fR1 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/ForAnalysis/26C-1.csv'
# fR01 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/ForAnalysis/26C-0.1.csv'
#4C
# fR10_4c = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/ForAnalysis/4C-10.csv'
# fR1_4c = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/ForAnalysis/4C-1.csv'
# fR01_4c = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/ForAnalysis/4C-0.1.csv'

# fnextday = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Dif scale/Onion1/5th/5th.csv'

#1 onion
# R10 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/2nd_set/26C/Onion1/10'
# R1 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/2nd_set/26C/Onion1/1'
# R01 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/2nd_set/26C/Onion1/0.1'
#2 onion
# R10 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/2nd_set/26C/Onion2/10'
# R1 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/2nd_set/26C/Onion2/1'
# R01 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/2nd_set/26C/Onion2/0.1'
#3 onion
# R10 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/2nd_set/26C/Onion3/10'
# R1 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/2nd_set/26C/Onion3/1'
# R01 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/2nd_set/26C/Onion3/0.1'

# 4C experiment
# R10_4c = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/4C_onion3/10mm_min'
# R1_4c = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/4C_onion3/1mm_min'
# R01_4c = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/4C_onion3/0.1mm_min'


group = [R10, R1, R01]
# group = [R10, R1, R01, R10_4c, R1_4c, R01_4c]

summary = pd.DataFrame()


# define a function that calls the index of certain value
def ind(df, value):
    for i in range(len(df)):
        if df.position[i] == value:
            return i

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

def ten_p(pull):
    for i in range(len(pull) - 5):
        if (pull.load[i] > 0.1) & (pull.load[i + 1] > 0.1) & (pull.load[i + 2] > 0.1) & \
                (pull.load[i + 3] > 0.1) & (pull.load[i + 4] > 0.1):
            point = pull.position[i]
            return point

# define a function that calculate the stress and strain of the curve
def norm(pull):
    ori_p = ten_p(first_pull)
    ## stress output
    pull['force'] = pull.load * 0.0098
    pull['stress'] = pull.force / (thickness * width * 0.001)
    ## load output
    # pull['stress'] = pull.load
    # for data output as strain (0 as 5mm)
    # pull['strain'] = (pull.position - ori_p) / (5000 + ori_p)
    # for data output as strain (500 as 5mm)
    pull['strain'] = (pull.position - ori_p) / (5000 + ori_p)
    # for data output as strain
    # pull['strain'] = pull.position * 4500/5000

    # pull['T_strain'] = np.log(1 + pull['strain'])
    # pull['T_stress'] = pull.force/(thickness * width * (ori_p + 4500) / (pull.position + 4500) * 0.001)

# define a function that smooth the curve
def smooth(pull):
    lowess = sm.nonparametric.lowess
    z = lowess(pull.stress, pull.strain, frac=0.05)
    return z

def derv(pull, inter):
    deriv = []
    for i in range(inter):
        deriv.append((pull.stress[i + inter] - pull.stress[i]) / (pull.strain[i + inter] - pull.strain[i]))
    for i in range(inter, len(pull) - inter):
        deriv.append(
            (pull.stress[i + inter] - pull.stress[i - inter]) / (pull.strain[i + inter] - pull.strain[i - inter]))
            # interval.append(pull.strain[i + inter] - pull.strain[i - inter])
    for i in range(len(pull) - inter, len(pull)):
        deriv.append((pull.stress[i] - pull.stress[i - inter]) / (pull.strain[i] - pull.strain[i - inter]))
    return deriv

def strain(pull, stress_OI):
    for i in range(len(pull)):
        if pull.stress[i] > stress_OI:
            strainT = pull.strain[i]
            return strainT

# path = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Nitro_glove_dif_rate_20g/Instron_with_Retract__03_18_2020__15_06_56_SHORT.csv'

#for file in sorted(glob.glob(os.path.join(path, '*SHORT.csv'))):

color = ['orange', 'red', 'purple']
# style = ['solid', 'solid', 'solid', 'dashed', 'dashed', 'dashed']
# marker = ['', '', '', 'v', 'v', 'v']
n = 1
g = 0
a = 0.8
# pla_outputlist = pd.DataFrame(columns= group_name)

for folder in group:
    p_defo = pd.DataFrame(columns= ['defo'])
    fig.set_size_inches(6.4, 4.8)

    for file in sorted(glob.glob(os.path.join(folder, '*SHORT.csv'))):
    # if n <= 3:
    #     n += 1
    #     continue
        whole_df = pd.read_table(file,
                             delimiter=',',
                             header=None,
                             names=range(0, 2),
                             )
        whole_df.columns = ['position', 'load']


        first_pull = whole_df[ind(whole_df, 'FIRST PULL') + 1:ind(whole_df, 'FIRST RETRACT')].astype(float).reset_index(drop=True)
        first_retract = whole_df[ind(whole_df, 'FIRST RETRACT') + 1:ind(whole_df, 'SECOND PULL')].astype(float).reset_index(drop=True).iloc[::-1]
        second_pull = whole_df[ind(whole_df, 'SECOND PULL') + 1:ind(whole_df, 'SECOND RETRACT')].astype(float).reset_index(drop=True)
        second_retract = whole_df[ind(whole_df, 'SECOND RETRACT') + 1:].astype(float).reset_index(drop=True).iloc[::-1]

        # first_retract.drop(first_retract.index[[purify_p(first_retract, target), len(first_retract)-1]], inplace = True)
        # second_retract.drop(second_retract.index[[purify_p(second_retract, target), len(second_retract)-1]], inplace = True)
        ori_p = ten_p(first_pull)
        norm(first_pull)
        # norm(first_retract)
        norm(second_pull)
        # norm(second_retract)


        retract = [first_retract, second_retract]
        pull = [first_pull, second_pull]

        for i in range(len(pull)):
            # retract[i].drop(retract[i].loc[purify_p(retract[i], 10): len(retract[i])].index, inplace=True)
            # retract[i].drop(retract[i].loc[0:ten_ind(retract[i], 0.1)].index, inplace=True)
            # pull[i].drop(pull[i].loc[0:ten_ind(pull[i], 0.1)].index, inplace=True)

            # retract[i].reset_index(inplace=True)
            rmindex = []
            # purify_revrs(retract[i])
            # retract[i].drop(retract[i].loc[rmindex].index, inplace=True)

            # retract[i] = retract[i].reset_index(drop=True)
            pull[i] = pull[i].reset_index(drop=True)

        sm1pull = pd.DataFrame(data=smooth(first_pull), columns=['strain', 'stress'])
        # sm1retract = pd.DataFrame(data=smooth(first_retract), columns=['strain', 'stress'])
        sm2pull = pd.DataFrame(data=smooth(second_pull), columns=['strain', 'stress'])
        # sm2retract = pd.DataFrame(data=smooth(second_retract), columns=['strain', 'stress'])

        #summary = pd.concat([summary, new_data],axis = 1)
        ax = plt.subplot(111)
        ax.plot(first_pull.strain *100, first_pull.stress, color = color[g], alpha = a)
        # ax.plot(second_pull.strain *100, second_pull.stress, color = color[g], alpha = a)


        n += 1

    ax.set(xlim = [-1, 22], ylabel='Stress (MPa)', xlabel = 'Strain', ylim = [-0.5, 5.5])

    plt.gcf().subplots_adjust(bottom=0.15, right=0.83)
    if g == 0:
        plt.savefig(''f'{output}/Strain_rate_dependency(plasticity_26)_10_firstpull.pdf', transparent=True)
    elif g == 1:
        plt.savefig(''f'{output}/Strain_rate_dependency(plasticity_26)_1_firstpull.pdf', transparent=True)
    elif g == 2:
        plt.savefig(''f'{output}/Strain_rate_dependency(plasticity_26)_01_firstpull.pdf', transparent=True)

    plt.show()

    g += 1

    # ax.set(ylabel='Stress (MPa)', xlabel = 'Strain', ylim = [-0.5, 5.5])

# plt.gcf().subplots_adjust(bottom=0.15, right= 0.83)
# plt.savefig(''f'{output}/Strain_rate_dependency(plasticity_4)_secondpull.pdf', transparent=True)
# plt.show()

