# script: calculate derivative of stress-strain curve

test
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
fig.set_size_inches(12, 4)
# for onion
thickness = 7
width = 3

output = '/Users/jingyiyu/Documents/Cosgrovelab/manuscript/CW_mechanics_Instron/Data_collection/Python_plot/'
output = '/Users/jingyiyu/Downloads'

# trail data
# folder1 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/trial_6.4.20/Onion1/5mm'
# folder2 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/trial_6.4.20/Onion2/5mm'
# folder3 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/trial_6.5.20/Onion3_peellevel_0.8'
mm5 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/trial_6.4.20/Onion2/5mm'
mm3 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/trial_6.4.20/Onion2/3mm'
# mm5 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/trial_6.4.20/Onion1/5mm'
# mm3 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/trial_6.4.20/Onion1/3mm'

# onion data
folder1 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/Yield_6.8.20/Onion1_8'
folder2 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/Yield_6.8.20/Onion2_10'
folder3 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/Yield_6.8.20/Onion3_10+'
# folder4 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/Yield.6.18.20/Onion1_10(last2_mixed_cell)'
folder4 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/Yield_6.23.20/Onion1_9_sample_is_good_quality'
folder5 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/Yield_7.22.20/onion1_10_allopen'
folder6 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/Yield_7.22.20/onion2_10+'
folder7 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/Yield_9.18.20'

O1 = '/Users/jingyiyu/Documents/Cosgrovelab/manuscript/CW_mechanics_Instron/Data_collection/Experimental_data/Yield/representative curves/O1'
O2 = '/Users/jingyiyu/Documents/Cosgrovelab/manuscript/CW_mechanics_Instron/Data_collection/Experimental_data/Yield/representative curves/O2'
O3 = '/Users/jingyiyu/Documents/Cosgrovelab/manuscript/CW_mechanics_Instron/Data_collection/Experimental_data/Yield/representative curves/O3'
O4 = '/Users/jingyiyu/Documents/Cosgrovelab/manuscript/CW_mechanics_Instron/Data_collection/Experimental_data/Yield/representative curves/O4'
O5 = '/Users/jingyiyu/Documents/Cosgrovelab/manuscript/CW_mechanics_Instron/Data_collection/Experimental_data/Yield/representative curves/O5'
O6 = '/Users/jingyiyu/Documents/Cosgrovelab/manuscript/CW_mechanics_Instron/Data_collection/Experimental_data/Yield/representative curves/O6'

# file = ''f'{folder}/Instron__02_15_2020__13_45_10_SHORT.csv'
# output = ''f'{folder}/derivative.csv'
group = [folder2, folder3, folder4, folder5, folder6, folder7]
group = [O1, O2, O3, O4, O5, O6]

# L = [3000, 5000]
L = [5000, 5000]
color = ['blue', 'green', 'red', 'purple', 'black', 'orange']

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

ax0 = plt.subplot(121)
ax0.set(title = 'Yield test', xlabel = 'Strain (%)', ylabel = 'Stress (MPa)')
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




ax = plt.subplot(122)
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
    c = 0
    for file in sorted(glob.glob(os.path.join(folder, '*SHORT.csv'))):
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
            ax0.plot(pull.strain * 100 , pull.stress, color = color[n],  linewidth = lw)
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
            ax.plot(smpull.strain *100, smpull.deriv, color = color[n], linewidth = lw)
            ax.set(ylim = [0, 100])

            c += 1
            count += 1


    n += 1
# setting for gram y axis
ax0y = ax0.twinx()
ax0y.set_ylabel('Force (mN)', rotation = 270, va = 'bottom')
ax0.set_yticks([0, 4, 8, 12, 16, 20])
ax0y.set_yticks([0, 80, 160, 240, 320])
B = -1 * 9.8
T = 40 * 9.8
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
plt.gcf().subplots_adjust(left = 0.1, right = 0.95, wspace= 0.4, bottom = 0.2)
# plt.subplots_adjust(hspace= 0.6)
fig = plt.gcf()
size = fig.get_size_inches()
print(size)
# plt.savefig(''f'{output}/Yield_rep.pdf', transparent = True)
plt.show()
