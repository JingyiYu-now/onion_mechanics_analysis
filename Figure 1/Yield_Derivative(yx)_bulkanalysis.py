# script: Plot the stress-strain curve and modulus-strain curve for representative s-s curve and all curves
# Fig 1 A, B use group2 for input
# Fig S1 use group 1 for input

import pandas as pd
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
fig = mpl.pyplot.gcf()
fig.set_size_inches(12, 4)
# for onion
thickness = 7
width = 3

output = '/Users/jingyiyu/Documents/Cosgrovelab/manuscript/CW_mechanics_Instron/Data_collection/Python_plot/'
output = '/Users/jingyiyu/Downloads'

# all onion data
folder2 = '/Users/jingyiyu/Documents/Cosgrovelab/manuscript/CW_mechanics_Instron/Data_collection/Experimental_data/Yield/Onion2_10'
folder3 = '/Users/jingyiyu/Documents/Cosgrovelab/manuscript/CW_mechanics_Instron/Data_collection/Experimental_data/Yield/Onion3_10+'
folder4 = '/Users/jingyiyu/Documents/Cosgrovelab/manuscript/CW_mechanics_Instron/Data_collection/Experimental_data/Yield/Onion1_9_sample_is_good_quality'
folder5 = '/Users/jingyiyu/Documents/Cosgrovelab/manuscript/CW_mechanics_Instron/Data_collection/Experimental_data/Yield/onion1_10_allopen'
folder6 = '/Users/jingyiyu/Documents/Cosgrovelab/manuscript/CW_mechanics_Instron/Data_collection/Experimental_data/Yield/onion2_10+'
folder7 = '/Users/jingyiyu/Documents/Cosgrovelab/manuscript/CW_mechanics_Instron/Data_collection/Experimental_data/Yield/Yield_9.18.20'
group1 = [folder2, folder3, folder4, folder5, folder6, folder7]

# representative curves for each onion
O1 = '/Users/jingyiyu/Documents/Cosgrovelab/manuscript/CW_mechanics_Instron/Data_collection/Experimental_data/Yield/representative curves/O1'
O2 = '/Users/jingyiyu/Documents/Cosgrovelab/manuscript/CW_mechanics_Instron/Data_collection/Experimental_data/Yield/representative curves/O2'
O3 = '/Users/jingyiyu/Documents/Cosgrovelab/manuscript/CW_mechanics_Instron/Data_collection/Experimental_data/Yield/representative curves/O3'
O4 = '/Users/jingyiyu/Documents/Cosgrovelab/manuscript/CW_mechanics_Instron/Data_collection/Experimental_data/Yield/representative curves/O4'
O5 = '/Users/jingyiyu/Documents/Cosgrovelab/manuscript/CW_mechanics_Instron/Data_collection/Experimental_data/Yield/representative curves/O5'
O6 = '/Users/jingyiyu/Documents/Cosgrovelab/manuscript/CW_mechanics_Instron/Data_collection/Experimental_data/Yield/representative curves/O6'
group2 = [O1, O2, O3, O4, O5, O6]

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
    for i in range(inter, len(pull) - inter):
        # stress
        deriv.append(
            (pull.stress[i + inter] - pull.stress[i - inter]) / (pull.strain[i + inter] - pull.strain[i - inter]))
    for i in range(len(pull) - inter, len(pull)):
        # stress
        deriv.append((pull.stress[i] - pull.stress[i - inter]) / (pull.strain[i] - pull.strain[i - inter]))
    return deriv

ax0 = plt.subplot(121)
ax0.set(title = 'Yield test', xlabel = 'Strain (%)', ylabel = 'Stress (MPa)')
blue_line = mlines.Line2D([], [], color = 'blue', label = 'Longitudinal' , alpha = 0.6)
green_line = mlines.Line2D([], [], color = 'green', label = 'Transverse' , alpha = 0.6)

ax = plt.subplot(122)
ax.set(ylim = [0, 150], xlabel = 'Strain (%)', ylabel = 'Modulus (MPa)')
deri_list = pd.DataFrame()
n = 0
count = 0
lw = 2
for folder in group1:
    c = 0
    for file in sorted(glob.glob(os.path.join(folder, '*SHORT.csv'))):
        if c < 5:
            whole_df = pd.read_csv(file, header=None,names = range(2))
            whole_df.columns = ['position', 'load']

            pull = whole_df[ind(whole_df, 'STRAIN [MICRONS]') + 1:].astype(float).reset_index(drop=True)

            norm(pull)

            # smooth
            smpull = pd.DataFrame(data=smooth(pull), columns=['strain', 'stress'])

            #### plot the stress-strain curve
            ax0.plot(pull.strain * 100 , pull.stress, color = color[n],  linewidth = lw)

            ### derivative from ∆y/∆x with window
            smpull['deriv'] = derv(smpull, 10)

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
plt.gcf().subplots_adjust(left = 0.1, right = 0.95, wspace= 0.4, bottom = 0.2)
fig = plt.gcf()
size = fig.get_size_inches()
print(size)
# plt.savefig(''f'{output}/Yield_rep.pdf', transparent = True)
plt.show()
