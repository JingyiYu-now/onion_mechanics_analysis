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

thickness = 7
width = 3
target = 15
inter = 5
Title = ''

#all onion
#26C
R10 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/ForAnalysis/26C/10'
R1 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/ForAnalysis/26C/1'
R01 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/ForAnalysis/26C/0.1'
#4C
# R10 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/ForAnalysis/4C/10'
# R1 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/ForAnalysis/4C/1'
# R01 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/Instron with retraction/Plasticity_strain_rate/ForAnalysis/4C/0.1'

group = [R10, R1, R01]

summary = pd.DataFrame()

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


color = ['orange', 'red', 'purple']
n = 1
g = 0
a = 0.8

for folder in group:
    p_defo = pd.DataFrame(columns= ['defo'])
    fig.set_size_inches(6.4, 4.8)

    for file in sorted(glob.glob(os.path.join(folder, '*SHORT.csv'))):
        whole_df = pd.read_table(file,
                             delimiter=',',
                             header=None,
                             names=range(0, 2),
                             )
        whole_df.columns = ['position', 'load']


        first_pull = whole_df[ind(whole_df, 'FIRST PULL') + 1:ind(whole_df, 'FIRST RETRACT')].astype(float).reset_index(drop=True)
        second_pull = whole_df[ind(whole_df, 'SECOND PULL') + 1:ind(whole_df, 'SECOND RETRACT')].astype(float).reset_index(drop=True)

        ori_p = ten_p(first_pull)
        norm(first_pull)
        norm(second_pull)

        pull = [first_pull, second_pull]

        for i in range(len(pull)):
            pull[i] = pull[i].reset_index(drop=True)

        ax = plt.subplot(111)
        ax.plot(first_pull.strain *100, first_pull.stress, color = color[g], alpha = a)

        n += 1

    ax.set(xlim = [-1, 22], ylabel='Stress (MPa)', xlabel = 'Strain', ylim = [-0.5, 5.5])

    plt.gcf().subplots_adjust(bottom=0.15, right=0.83)
    # if g == 0:
        # plt.savefig(''f'{output}/Strain_rate_dependency(plasticity_26)_10_firstpull.pdf', transparent=True)
    # elif g == 1:
        # plt.savefig(''f'{output}/Strain_rate_dependency(plasticity_26)_1_firstpull.pdf', transparent=True)
    # elif g == 2:
        # plt.savefig(''f'{output}/Strain_rate_dependency(plasticity_26)_01_firstpull.pdf', transparent=True)

    plt.show()

    g += 1

