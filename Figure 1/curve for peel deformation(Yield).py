import matplotlib.lines as mlines

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
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
output = '/Users/jingyiyu/Documents/Cosgrovelab/manuscript/CW_mechanics_Instron/Data_collection/Python_plot/'

# sample dimension
thickness = 7
width = 3
#plot boundry
left_lim = -1
right_lim = 51
color = ['blue', 'green', 'orange', 'black']

# peel curve used in manuscript v1 - v6
V1 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/Video_peel_deformation/Onion2_4.22.21_+-coverslide/Yield__04_22_2021__14_46_39_SHORT.csv'


whole_df = pd.read_table(V1,
                         delimiter=',',
                         header=None,
                         names=range(0, 2),
                         )

whole_df.columns = ['position', 'load']

def ten_p(pull):
    for i in range(len(pull) - 5):
        if (pull.load[i] > 0.1) & (pull.load[i + 1] > 0.1) & (pull.load[i + 2] > 0.1) & \
                (pull.load[i + 3] > 0.1) & (pull.load[i + 4] > 0.1):
            point = pull.position[i]
            return point

# define a function to normalize the data
def norm(pull):
    ori_p = ten_p(pull)
    pull['force'] = pull.load * 0.0098
    pull['stress'] = pull.force/(thickness * width * 0.001)
    pull['strain'] = (pull.position - ori_p) / (5000 + ori_p)

# define a function that calls the index of certain value
def ind(df, value):
    for i in range(len(df)):
        if df.position[i] == value:
            return i

# define a function that smooth the curve
def smooth(pull, x, y):
    lowess = sm.nonparametric.lowess
    z = lowess(pull[y], pull[x], frac=0.03)
    return z

# calculate derivative of the curve using ∆y/∆x
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

pull = whole_df[ind(whole_df, 'STRAIN [MICRONS]') + 1:].astype(float).reset_index(drop=True)

norm(pull)

# plot setting for curve
ax = plt.subplot(111)
ax.set(xlim = [left_lim, right_lim], xlabel ='Strain(%)', ylabel ='Stress(MPa)')

#setting for gram y axis
axy = ax.twinx()
axy.set_ylabel('Force (mN)', rotation = 270, va = 'bottom')
B = -1 * 9.8
T = 36 * 9.8
axy.set(ylim = [B,T])
SB = B * 0.001/(thickness * width * 0.001)
ST = T * 0.001/(thickness * width * 0.001)
ax.set(ylim = [SB,ST])

ax.plot(pull.strain * 100, pull.stress, color = 'blue')

fig = mpl.pyplot.gcf()
fig.set_size_inches(8, 5)
plt.gcf().subplots_adjust(bottom = 0.5)
# plt.savefig(''f'{output}/Yield_curve_deformation.pdf', transparent = True)

plt.show()
