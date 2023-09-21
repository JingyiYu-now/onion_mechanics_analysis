# script: calculate derivative of stress-strain curve

import pandas as pd
import os
import glob
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.lines as mlines
import scipy.stats as stats
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
fig.set_size_inches(6.4, 8)
# for onion
thickness = 7
width = 3

output = '/Users/jingyiyu/Documents/Cosgrovelab/manuscript/CW_mechanics_Instron/Data_collection/Python_plot/'

#Anisotropical
Londi = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/transverse_longi_compare_11.10.20/Longitudinal'
Trans = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/transverse_longi_compare_11.10.20/transverse'
Londi2 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/longi_vs_trans/1.6.21/Onion1/longi'
Trans2 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/longi_vs_trans/1.6.21/Onion1/Trans'
Londi3 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/longi_vs_trans/1.7.21/longi'
Trans3 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/longi_vs_trans/1.7.21/trans'

represent_longi = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/longi_vs_trans/represent/longi'
represent_trans = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/longi_vs_trans/represent/trans'
# Lon = [Londi, Londi2, Londi3]
# Tran = [Trans, Trans2, Trans3]

Lon = [represent_longi]
Tran = [represent_trans]


group = [Lon, Tran]

orientation = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/Strecthed from 3 direction'

# file = ''f'{folder}/Instron__02_15_2020__13_45_10_SHORT.csv'
# output = ''f'{folder}/derivative.csv'
# group = [folder2, folder3, folder4, folder5, folder6, folder7]
# group = [Londi, Trans]
# group = [Londi2, Trans2]
# group = [Londi3, Trans3]
# group = [orientation]
# group = [mm3, mm5]
# group = [Ara_3mm, Ara_5mm]
# group = [onion_3_10, onion_2_10]
# L = [3000, 5000]
L = [5000, 5000]
# color = ['blue', 'green', 'red', 'purple', 'black', 'orange']
color = ['blue', 'red', 'green', 'purple', 'black', 'orange']

# define a function that calls the index of certain value
def ind(df, value):
    for i in range(len(df)):
        if df.position[i] == value:
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
    ori_p = ten_p(pull)
    pull['force'] = pull.load * 0.0098
    ## engineering stress/strain
    pull['stress'] = pull.force / (thickness * width * 0.001)
    pull['strain'] = (pull.position - ori_p) / (5000 + ori_p)
    ## Convert to true stress,strain
    # pull['stress'] = pull.force/(thickness * width * 5000 / (pull.position + 5000) * 0.001)
    # pull['strain'] = np.log(1 + (pull.position - ori_p) / 5000)

    # pull['T_strain'] = np.log(1 + pull['strain'])
    # pull['T_stress'] = pull.force/(thickness * width * (ori_p + 4500) / (pull.position + 4500) * 0.001)


# define a function that smooth the curve
def smooth(pull):
    lowess = sm.nonparametric.lowess
    # stress
    z = lowess(pull.stress, pull.strain, frac=0.025)
    # load
    # z = lowess(pull.load, pull.strain, frac=0.025)
    return z


interval = []
# define a function that calculate the derivative of s-s curve using ∆y/∆x
def derv(pull, inter):
    deriv = []
    for i in range(inter):
        # stress
        deriv.append((pull.stress[i + inter] - pull.stress[i]) / (pull.strain[i + inter] - pull.strain[i]))
        # load
        # deriv.append((pull.load[i + inter] - pull.load[i]) / (pull.strain[i + inter] - pull.strain[i]))
    for i in range(inter, len(pull) - inter):
        # stress
        deriv.append(
            (pull.stress[i + inter] - pull.stress[i - inter]) / (pull.strain[i + inter] - pull.strain[i - inter]))
        # interval.append(pull.strain[i + inter] - pull.strain[i - inter])
        # load
        # deriv.append(
        #     (pull.load[i + inter] - pull.load[i - inter]) / (pull.strain[i + inter] - pull.strain[i - inter]))
    for i in range(len(pull) - inter, len(pull)):
        # stress
        deriv.append((pull.stress[i] - pull.stress[i - inter]) / (pull.strain[i] - pull.strain[i - inter]))
        # load
        # deriv.append((pull.load[i] - pull.load[i - inter]) / (pull.strain[i] - pull.strain[i - inter]))
    return deriv

ax0 = plt.subplot(211)
ax0.set(xlabel = 'Strain (%)', ylabel = 'Stress (MPa)')
# ax0.grid(alpha = 0.4, linestyle = '--', linewidth = 2)
# blue_line = mlines.Line2D([], [], color = 'blue', label = 'Onion 1' , alpha = 0.6)
# green_line = mlines.Line2D([], [], color = 'green', label = 'Onion 2' , alpha = 0.6)
# red_line = mlines.Line2D([], [], color = 'red', label = 'Onion 3' , alpha = 0.6)
# purple_line = mlines.Line2D([], [], color = 'purple', label = 'Onion 4' , alpha = 0.6)
# black_line = mlines.Line2D([], [], color = 'black', label = 'Onion 5' , alpha = 0.6)
# orange_line = mlines.Line2D([], [], color = 'orange', label = 'Onion 6' , alpha = 0.6)
# longitudinal and transverse comparison
blue_line = mlines.Line2D([], [], color = 'blue', label = 'Longitudinal' , alpha = 0.6)
green_line = mlines.Line2D([], [], color = 'green', label = 'Transverse' , alpha = 0.6)


# ax0.legend(handles = [blue_line, green_line, red_line, purple_line, black_line, orange_line], loc = 2, fontsize = 8)
# ax0.legend(handles = [blue_line, green_line], loc = 2, fontsize = 8)




ax = plt.subplot(212)
ax.set(ylim = [0, 150], xlabel = 'Strain (%)', ylabel = 'Modulus (MPa)')
# ax.grid(alpha = 0.4, linestyle = '--', linewidth = 2)
# ax.legend(handles = [blue_line, green_line], loc = 2, fontsize = 8)

# # ax.legend(handles = [blue_line, green_line, red_line, purple_line], loc = 1, fontsize = 8)
#
# ax2 = plt.subplot(212)
# ax2.set(ylim = [0, 120], xlabel = 'Stress(MPa)', ylabel = 'Derivative(MPa)', title = 'Derivative(∆y/∆x) as function of stress')
# ax2.grid(alpha = 0.4, linestyle = '--')
# ax2.legend(handles = [blue_line, green_line, red_line, purple_line], loc = 2, fontsize = 8)

deri_list = pd.DataFrame()
n = 0
count = 0
lw = 2
for folder in group:
    for unit in folder:
        c = 0
        for file in sorted(glob.glob(os.path.join(unit, '*SHORT.csv'))):
            if c < 5:
                whole_df = pd.read_csv(file, header=None,names = range(2))
                whole_df.columns = ['position', 'load']

                pull = whole_df[ind(whole_df, 'STRAIN [MICRONS]') + 1:].astype(float).reset_index(drop=True)

                norm(pull)

                # stress
                # smpull = pd.DataFrame(data=smooth(pull), columns=['strain', 'stress'])
                # load
                smpull = pd.DataFrame(data=smooth(pull), columns=['strain', 'stress'])


                #### plot the stress-strain curve
                        # for in Onion 1&2 curve plot
                # ax0.plot(pull.strain * 100 , pull.stress, color = 'blue')
                       # for group wise single pull plot
                ax0.plot(pull.strain * 100 , pull.stress, color = color[n], alpha = 0.6,  linewidth = lw)
                # ax0.plot(pull.strain * 100 , pull.load, color = color[n], alpha = 0.6)
                # ax0y.plot(pull.strain * 100 , pull.load, color = color[n], alpha = 0)

                ### derivative from ∆y/∆x with window
                smpull['deriv'] = derv(smpull, 10)


                        # plot derivative over strain
                        # for in Onion 1&2 curve plot
                # ax.plot(pull.strain * 100, smpull.deriv, color='blue')

                        # for group wise single pull plot
                # ax.plot(smpull.strain * 100, smpull.deriv, color = color[n], alpha = 0.6)


                        # plot derivative over stress
                        # for in Onion 1&2 curve plot
                # ax2.plot(pull.stress, smpull.deriv, color=color[n], alpha = 0.6)

                        # for group wise single pull plot
                ax.plot(smpull.strain *100, smpull.deriv, color = color[n], alpha = 0.6, linewidth = lw)
                ax.set(ylim = [0, 120])

                c += 1
                count += 1


    n += 1
# setting for gram y axis
ax0y = ax0.twinx()
ax0y.set_ylabel('Loading (mN)', rotation = 270, va = 'bottom')
ax0.set_yticks([0, 5, 10, 15, 20, 25])
ax0y.set_yticks([0, 100, 200, 300, 400, 500])
B = -1 * 9.8
T = 55 * 9.8
ax0y.set(ylim = [B,T])
SB = B * 0.001/(thickness * width * 0.001)
ST = T * 0.001/(thickness * width * 0.001)
ax0.set(ylim = [SB,ST])
print(count)
ax.set_yticks([0, 25, 50, 75, 100])
# ax0.xaxis.set_tick_params(which='major', size=10, width=2, direction='out')
# ax0.yaxis.set_tick_params(which='major', size=10, width=2, direction='out')
# ax0y.yaxis.set_tick_params(which='major', size=10, width=2, direction='out')
# ax.xaxis.set_tick_params(which='major', size=10, width=2, direction='out')
# ax.yaxis.set_tick_params(which='major', size=10, width=2, direction='out')
blue_line = mlines.Line2D([], [], color = 'blue', label = 'Longitudinal' , alpha = 0.6)
green_line = mlines.Line2D([], [], color = 'red', label = 'Transverse' , alpha = 0.6)
ax.legend(handles = [blue_line, green_line], loc = 2, fontsize=10, frameon = False)
ax0.legend(handles = [blue_line, green_line], loc = 2, fontsize=10, frameon = False)

plt.gcf().subplots_adjust(bottom=0.15, right = 0.8, left = 0.25, hspace= 0.3)
# plt.subplots_adjust(hspace= 0.6)
plt.savefig(''f'{output}/Yield_anisotropy_3rep.pdf', transparent = True)
size = fig.get_size_inches()
print(size)
plt.show()
