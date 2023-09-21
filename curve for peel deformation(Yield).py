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


thickness = 7
width = 3
l = -1
r = 51
color = ['blue', 'green', 'orange', 'black']

Y1 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/Yield_6.8.20/Onion3_10+/Yield__06_08_2020__14_16_28_SHORT.csv'
Y2 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/Yield_6.8.20/Onion3_10+/Yield__06_08_2020__14_23_07_SHORT.csv'

S11 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/slipage_check_curve7.7.20/peel1_peel2/Yield__07_07_2020__14_13_28_SHORT.csv'
S12 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/slipage_check_curve7.7.20/peel1_peel2/Yield__07_07_2020__14_14_41_SHORT.csv'
S21 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/slipage_check_curve7.7.20/peel1_peel2/Yield__07_07_2020__14_19_55_SHORT.csv'
S22 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/slipage_check_curve7.7.20/peel1_peel2/Yield__07_07_2020__14_20_56_SHORT.csv'
S31 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/slipage_check_curve7.7.20/pre-test_sample_from_same_onion/Yield__07_07_2020__13_45_03_SHORT.csv'
S32 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/slipage_check_curve7.7.20/pre-test_sample_from_same_onion/Yield__07_07_2020__13_46_12_SHORT.csv'

#for vedio data
V11 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/Yield_7.22.20/onion2_10+/Yield__07_22_2020__10_31_22_SHORT.csv'
V12 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/Yield_7.22.20/onion2_10+/Yield__07_22_2020__10_32_44_SHORT.csv'

# peel curve used in manuscript v1 - v6
V13 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/Video_peel_deformation/Onion2_4.22.21_+-coverslide/Yield__04_22_2021__14_46_39_SHORT.csv'


# video from 5.3.21 (v7)
V2 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/5th_yield/Video_peel_deformation/Onion4_5.3.21_+-coer/Yield__05_03_2021__15_32_01_SHORT.csv'

# whole_df = pd.read_table(V13,
#                          delimiter=',',
#                          header=None,
#                          names=range(0, 2),
#                          )
#
# whole_df.columns = ['position', 'load']

whole_df2 = pd.read_table(V2,
                         delimiter=',',
                         header=None,
                         names=range(0, 2),
                         )

whole_df2.columns = ['position', 'load']


def ten_p(pull):
    for i in range(len(pull) - 5):
        if (pull.load[i] > 0.1) & (pull.load[i + 1] > 0.1) & (pull.load[i + 2] > 0.1) & \
                (pull.load[i + 3] > 0.1) & (pull.load[i + 4] > 0.1):
            point = pull.position[i]
            return point


# define a function to normalize the data
def norm(pull):
    ori_p = ten_p(ori_pull)
    pull['force'] = pull.load * 0.0098
    pull['stress'] = pull.force/(thickness * width * 0.001)
    pull['strain'] = (pull.position - ori_p) / 5000
    pull['T_strain'] = np.log(1 + pull['strain'])
    ##### poisson's ratio = 0.5
    pull['T_stress'] = pull.force/(thickness * width * 5000 / (pull.position + 5000) * 0.001)

    ##### poisson's ratio = 0.3
    # pull['T_stress_3'] = pull.force/(thickness * width * 0.001 * (1 - 0.5 * pull.T_strain) ** 2 )

    ##### poisson's ratio = 0.5, if no comparison in true and engi
    # pull['strain'] = np.log(1 + (pull.position - ori_p) / (ori_p+4500))
    # pull['stress'] = pull.force/(thickness * width * (ori_p + 4500) / (pull.position + 4500) * 0.001)


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

def fsmooth(x, y):
    lowess = sm.nonparametric.lowess
    z = lowess(y, x, frac=0.05)
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

ori_pull = whole_df2[ind(whole_df2, 'STRAIN [MICRONS]') + 1:].astype(float).reset_index(drop=True)
# pull = whole_df[ind(whole_df, 'STRAIN [MICRONS]') + 1:].astype(float).reset_index(drop=True)

# norm(pull)
norm(ori_pull)
# print(ten_p(pull)/5000 *100)

osmee = pd.DataFrame(data=smooth(ori_pull, 'strain','stress'), columns=['strain', 'stress'])
# smee = pd.DataFrame(data=smooth(pull, 'strain','stress'), columns=['strain', 'stress'])
# smet = pd.DataFrame(data=smooth(pull, 'T_strain','stress'), columns=['strain', 'stress'])
# smte = pd.DataFrame(data=smooth(pull, 'strain','T_stress'), columns=['strain', 'stress'])
# smtt = pd.DataFrame(data=smooth(ori_pull, 'T_strain','T_stress'), columns=['strain', 'stress'])

# plot setting for curve
ax = plt.subplot(111)
ax.set(xlim = [l,r], xlabel = 'Strain(%)', ylabel = 'Stress(MPa)')
# ax.grid(alpha = 0.4, linestyle = '--')
# for time indication in video
# ax2 = ax.twiny()
# ax2.set_xlim([l*0.01*5000/400* 1.08, r*0.01*5000/400*1.08])
# ax2.set_xlabel('recording time(s)')

#setting for gram y axis
axy = ax.twinx()
axy.set_ylabel('Load (mN)', rotation = 270, va = 'bottom')
B = -1 * 9.8
T = 36 * 9.8
axy.set(ylim = [B,T])
SB = B * 0.001/(thickness * width * 0.001)
ST = T * 0.001/(thickness * width * 0.001)
ax.set(ylim = [SB,ST])

# bx = plt.subplot(212)
# bx.set(xlim = [l,r], xlabel = 'Strain(%)', ylabel = 'Modulus(MPa)')
# bx.grid(alpha = 0.4, linestyle = '--')
# for time indication in video
# bx2 = bx.twiny()
# bx2.set_xlim([l*0.01*5000/400, r*0.01*5000/400])
# bx2.set_xlabel('recording time(s)')

# cx = plt.subplot(413)
# cx.set(xlabel = 'Stress(MPa)', ylabel = 'Derivative(MPa)')
# cx.grid(alpha = 0.4, linestyle = '--')
# dx = plt.subplot(414)
# dx.set(xlabel = 'Stain(%)', ylabel = 'Stress(MPa)')
# dx.grid(alpha = 0.4, linestyle = '--')
# dx.set_yscale('log')

ax.plot(osmee.strain * 100, osmee.stress, color = 'blue')
# ax.plot(smee.strain * 100, smee.stress, color = 'green')
# ax.plot(smet.strain * 100, smet.stress, color = 'green')
# ax.plot(smte.strain * 100, smte.stress, color = 'orange')
# ax.plot(smtt.strain * 100, smtt.stress, color = 'green')

# bx.plot(osmee.strain * 100, derv(osmee, 5), color = 'blue')
# bx.plot(smee.strain * 100, derv(smee, 5), color = 'green')
# bx.plot(smet.strain * 100, derv(smet, 5), color = 'green')
# bx.plot(smte.strain * 100, derv(smte, 5), color = 'orange')
# bx.plot(smtt.strain * 100, derv(smtt, 5), color = 'green')

# cx.plot(osmee.stress, derv(osmee, 5), color = 'blue')
# cx.plot(smee.stress, derv(smee, 5), color = 'green')

# dx.plot(osmee.strain * 100, osmee.stress, color = 'blue')
# dx.plot(smee.strain * 100, smee.stress, color = 'green')

# true curve comparison
# blue_line = mlines.Line2D([], [], color = 'blue', label = 'Engineering' )
# green_line = mlines.Line2D([], [], color = 'green', label = 'True(v=0.5)' )
# orange_line = mlines.Line2D([], [], color = 'Orange', label = 'True(v=0.3)' )
# ax.legend(handles = [blue_line, green_line])

# plot width and thickness
# ax = plt.subplot(211)
# ax.set(xlabel = 'Engineering Strain(%)', ylabel = 'Width(mm)', title = 'Changes of width along stretching (v=0.5)')
# bx = plt.subplot(212)
# bx.set(xlabel = 'Engineering Strain(%)', ylabel = 'Thickness(µm)', title = 'Changes of thickness along stretching (v=0.5)')
#
# ax.plot(pull.strain * 100, width * (1 - 0.5 * pull.T_strain))
# bx.plot(pull.strain * 100, thickness * (1 - 0.5 * pull.T_strain))



# reason for tilt comparison
# blue_line = mlines.Line2D([], [], color = 'blue', label = 'E stress & E strain' )
# green_line = mlines.Line2D([], [], color = 'green', label = 'E stress & T strain' )
# orange_line = mlines.Line2D([], [], color = 'Orange', label = 'T stress & E strain' )
# black_line = mlines.Line2D([], [], color = 'black', label = 'T stress & T strain' )

# ax.legend(handles = [blue_line, green_line, orange_line, black_line])
# leg = ax.legend(handles = [black_line, orange_line, green_line, blue_line],fontsize=7)
# for line, text in zip(leg.get_lines(), leg.get_texts()):
#     text.set_color(line.get_color())
# leg = bx.legend(handles = [black_line, orange_line, green_line, blue_line],fontsize=7)
# for line, text in zip(leg.get_lines(), leg.get_texts()):
#     text.set_color(line.get_color())
fig = mpl.pyplot.gcf()
fig.set_size_inches(8, 5)
plt.gcf().subplots_adjust(bottom = 0.5)
# plt.savefig(''f'{output}/Yield_curve_deformation.pdf', transparent = True)

plt.show()
