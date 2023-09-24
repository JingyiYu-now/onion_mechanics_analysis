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

# Make dataframes that include all curves for average
Slist = np.arange(0, 4.5, 0.05)
# first loading curve
R10c = pd.DataFrame(columns=Slist)
R1c = pd.DataFrame(columns=Slist)
R01c = pd.DataFrame(columns=Slist)

# second loading curve
R10c2 = pd.DataFrame(columns=Slist)
R1c2 = pd.DataFrame(columns=Slist)
R01c2 = pd.DataFrame(columns=Slist)


S_group = [R10c, R1c, R01c]
S_group2 = [R10c2, R1c2, R01c2]
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
    pull['strain'] = (pull.position - ori_p) / (5000 + ori_p)

# define a function that smooth the curve
def smooth(pull):
    lowess = sm.nonparametric.lowess
    z = lowess(pull.stress, pull.strain, frac=0.05)
    return z

def strain(pull, stress_OI):
    for i in range(len(pull)):
        if pull.stress[i] > stress_OI:
            strainT = pull.strain[i]
            return strainT

color = ['blue', 'orange', 'red']
n = 1
g = 0
a = 0.8
# pla_outputlist = pd.DataFrame(columns= group_name)

for folder in group:
    p_defo = pd.DataFrame(columns= ['defo'])

    for file in sorted(glob.glob(os.path.join(folder, '*SHORT.csv'))):
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

        ori_p = ten_p(first_pull)
        norm(first_pull)
        norm(second_pull)

        retract = [first_retract, second_retract]
        pull = [first_pull, second_pull]

        for i in range(len(pull)):
            pull[i] = pull[i].reset_index(drop=True)

        sm1pull = pd.DataFrame(data=smooth(first_pull), columns=['strain', 'stress'])
        sm2pull = pd.DataFrame(data=smooth(second_pull), columns=['strain', 'stress'])
        strain_num = []
        strain_num2 = []

        for i in range(len(Slist)):
            strain_num.append(strain(sm1pull, Slist[i]))
            strain_num2.append(strain(sm2pull, Slist[i]))
        S_group[g].loc[len(S_group[g])] = strain_num
        S_group2[g].loc[len(S_group2[g])] = strain_num2

        pla_defo = []

        pla_defo.append((ten_p(second_pull) - ori_p) / (5000 + ori_p) * 100)
        p_defo.loc[len(p_defo)] = pla_defo

        n += 1
        print(n)


    g += 1

ax = plt.subplot(111)

y = Slist

s10 = []
s1 = []
s01 = []
se10 = []
se1 = []
se01 = []

#second pull
s10s = []
s1s = []
s01s = []
se10s = []
se1s = []
se01s = []

for i in range(len(Slist)):
    ### average curve
    # first curve
    s10.append(np.mean(R10c.iloc[:, i]) *100)
    s1.append(np.mean(R1c.iloc[:, i]) *100)
    s01.append(np.mean(R01c.iloc[:, i]) *100)
    #second curve
    s10s.append(np.mean(R10c2.iloc[:, i]) *100)
    s1s.append(np.mean(R1c2.iloc[:, i]) *100)
    s01s.append(np.mean(R01c2.iloc[:, i]) *100)

    ### stand error of the mean
    # first curve
    se10.append((np.std(R10c.iloc[:, i]) *100)/(len(R10c[0]) ** 1/2))
    se1.append((np.std(R1c.iloc[:, i]) *100)/(len(R1c[0]) ** 1/2))
    se01.append((np.std(R01c.iloc[:, i]) *100)/(len(R01c[0]) ** 1/2))
    # second curve
    se10s.append((np.std(R10c2.iloc[:, i]) *100)/(len(R10c2[0]) ** 1/2))
    se1s.append((np.std(R1c2.iloc[:, i]) *100)/(len(R1c2[0]) ** 1/2))
    se01s.append((np.std(R01c2.iloc[:, i]) *100)/(len(R01c2[0]) ** 1/2))

# plot the average first loading curve
ax.errorbar(s10 , y, color= 'orange', alpha = 0.6)
ax.errorbar(s1 , y, color='red', alpha = 0.6)
ax.errorbar(s01 , y, color='purple', alpha = 0.6)

#plot the average second loading curve
ax.errorbar(s10s , y, color= 'orange', alpha = 0.6)
ax.errorbar(s1s , y, color='red', alpha = 0.6)
ax.errorbar(s01s , y, color='purple', alpha = 0.6)
#

ax.set(ylabel='Stress (MPa)', xlabel = 'Strain', ylim = [-0.5, 5.5])

orange_line = mlines.Line2D([], [], color = 'orange', label = '10 mm/min')
blue_line = mlines.Line2D([], [], color = 'red', label = '1 mm/min')
green_line = mlines.Line2D([], [], color = 'purple', label = '0.1 mm/min')
ax.legend(handles = [orange_line, blue_line, green_line], loc = 2, fontsize=10, frameon=False)


plt.gcf().subplots_adjust(bottom=0.15, right= 0.83)
plt.show()

