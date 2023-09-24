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
# fig.set_size_inches(5.5, 4.8)
thickness = 7
width = 3
smf = 0.05
intv = 10
per = 10
num = 40

output = '/Users/jingyiyu/Documents/Cosgrovelab/manuscript/CW_mechanics_Instron/Data_collection/Python_plot/'

Title = 'Cyclic loading(10 loads)_20g_'
### room temp fruit strips
strip1 = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/room temperature/strip1'
strip2 = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/room temperature/strip2'
strip3 = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/room temperature/strip3'
strip4 = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/room temperature/strip4'
### strips area
# area_list = [3552762.637, 2927490.149, 3650148.381, 2944460.14]
# group = [strip1, strip2, strip3, strip4]

#dif target force test
# strip_s1 = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/room temperature/strips-1_20g'
# strip_s2 = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/room temperature/strips-2_20g'
# strip_s3 = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/room temperature/strips-3_10g'
# strip_s4 = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/room temperature/strips-4_30g'
# strip_s5 = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/room temperature/strips-5_10g'
### standard strips area
# area_list = [3389297.281, 2679108.159, 3082287.14, 3597204.676, 2974816.086]
# group = [strip_s1, strip_s2, strip_s3, strip_s4, strip_s5]

#representative curve (room temp)
Strip_26C = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/room temperature/strips-3_10g'
area_list = [3082287.14]
group = [Strip_26C]

### 4 degree fruit strips
# strip1 = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/4C/strip1_10g'
# strip2 = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/4C/strip2_10g'
# strip3 = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/4C/strip3_10g'
# strip4 = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/4C/strip4_10g'
# area_list = [4373323.381, 4215275.754, 3970144.4, 3568115.645]
# group = [strip1, strip2, strip3, strip4]

#dif target foce test
strip_s1 = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/4C/strips-1_20g'
strip_s2 = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/4C/strips-2_20g'
strip_s3 = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/4C/strips-3_10g'
strip_s4 = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/4C/strips-4_30g'
strip_s5 = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/4C/strips-5_10g'
### standard strips area
# area_list = [3251820.35, 2591324.295, 3610690.651, 3374927.8, 3102991.195]
# group = [strip_s1, strip_s2, strip_s3, strip_s4, strip_s5]

# representative 4 degree strip
# strip_4C = '/Users/jingyiyu/Documents/Cosgrovelab/Fruit strip/fruit strip seqWretract/strain rate dependency/4C/strips-3_10g'
# area_list = [3610690.651]


# target for seuqntial Instron
target = [2,4,6,8,10,10,12,14,16,18,20,20]

Pull = pd.DataFrame(columns= target)
Retract = pd.DataFrame(columns= target)

Strain_4C = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/4C_experiment/Onion_4c'
Strain_4C_2 = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/4C_experiment/Onion_4c/Onion3'
Strain_RT_Con = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/4C_experiment/Onion_RT_control'

# 20g
Strain_4c_20g = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/4C_experiment/Onion_4c_20g'
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

def ten_ind(pull):
    for i in range(len(pull) - 5):
        if (pull.load[i] > 0.1) & (pull.load[i + 1] > 0.1) & (pull.load[i + 2] > 0.1) & \
                (pull.load[i + 3] > 0.1) & (pull.load[i + 4] > 0.1):
            return i


# define a function that calculate the stress and strain of the curve
def norm(pull, area_um2):
    # strain calculated based on very original length (before plastic deformation)
    ori_p = ten_p(first_pull)

    pull['force_N'] = pull.load * 0.0098

    # engineering stress & strain
    pull['stress'] = pull.force_N / area_um2 * 10 ** 9
    pull['strain'] = (pull.position - ori_p) / (5000 + ori_p)


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
            # interval.append(pull.strain[i + inter] - pull.strain[i - inter])
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
            # interval.append(pull.strain[i + inter] - pull.strain[i - inter])
    for i in range(len(y) - inter, len(y)):
        deriv.append((y[i] - y[i - inter]) / (x[i] - x[i - inter]))
    return deriv

color = ['C0','C1','C2', 'C3', 'C4', 'C5','C6','C7','C8','C9','C10','C11']
n = 0
cc = 0
for folder in group:
    area = area_list[cc]
    n = 0
    for file in sorted(glob.glob(os.path.join(folder, '*SHORT.csv'))):

        if n%2 == 0:
            print(n)
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

            ori_p = ten_p(first_pull)
            n += 1
            continue

        if n%2 == 1:
            print(n)
            whole_df = pd.read_table(file,
                                 delimiter=',',
                                 header=None,
                                 names=range(0, 2),
                                 )
            whole_df.columns = ['position', 'load']


            seventh_pull = whole_df[ind(whole_df, 'FIRST PULL') + 1:ind(whole_df, 'FIRST RETRACT')].astype(float).reset_index(drop=True)
            seventh_retract = whole_df[ind(whole_df, 'FIRST RETRACT') + 1:ind(whole_df, 'SECOND PULL')].astype(float).iloc[::-1].reset_index(drop=True)
            eighth_pull = whole_df[ind(whole_df, 'SECOND PULL') + 1:ind(whole_df, 'SECOND RETRACT')].astype(float).reset_index(drop=True)
            eighth_retract = whole_df[ind(whole_df, 'SECOND RETRACT') + 1:ind(whole_df, 'THIRD PULL')].astype(float).iloc[::-1].reset_index(drop=True)
            ninth_pull = whole_df[ind(whole_df, 'THIRD PULL') + 1:ind(whole_df, 'THIRD RETRACT')].astype(float).reset_index(drop=True)
            ninth_retract = whole_df[ind(whole_df, 'THIRD RETRACT') + 1:ind(whole_df, 'FOURTH PULL')].astype(float).iloc[::-1].reset_index(drop=True)
            tenth_pull = whole_df[ind(whole_df, 'FOURTH PULL') + 1:ind(whole_df, 'FOURTH RETRACT')].astype(float).reset_index(drop=True)
            tenth_retract = whole_df[ind(whole_df, 'FOURTH RETRACT') + 1:ind(whole_df, 'FIFTH PULL')].astype(float).iloc[::-1].reset_index(drop=True)
            eleventh_pull = whole_df[ind(whole_df, 'FIFTH PULL') + 1:ind(whole_df, 'FIFTH RETRACT')].astype(float).reset_index(drop=True)
            eleventh_retract = whole_df[ind(whole_df, 'FIFTH RETRACT') + 1:ind(whole_df, 'SIXTH PULL')].astype(float).iloc[::-1].reset_index(drop=True)
            twlveth_pull = whole_df[ind(whole_df, 'SIXTH PULL') + 1:ind(whole_df, 'SIXTH RETRACT')].astype(float).reset_index(drop=True)
            twlveth_retract = whole_df[ind(whole_df, 'SIXTH RETRACT') + 1:].astype(float).iloc[::-1].reset_index(drop=True)


            pre_p = [first_pull, second_pull, third_pull, fourth_pull, fifth_pull, sixth_pull, seventh_pull]
            pre_r = [first_retract, second_retract, third_retract, fourth_retract, fifth_retract, sixth_retract, seventh_retract]
            rate_p = [eighth_pull, ninth_pull, tenth_pull, eleventh_pull, twlveth_pull]
            rate_r = [eighth_retract, ninth_retract, tenth_retract, eleventh_retract, twlveth_retract]

            ### normalize curve
            for i in range(len(pre_p)):
                norm(pre_p[i], area)
                norm(pre_r[i], area)
                # print((ten_p(pre_p[i]) - ori_p)/5000)
            for i in range(len(rate_p)):
                norm(rate_p[i], area)
                norm(rate_r[i], area)
                # print(ten_p(rate_p[i])/5000)

            for i in range(len(pre_r)):
                pre_r[i].drop(pre_r[i].loc[purify_p(pre_r[i], 10): len(pre_r[i])].index, inplace=True)
            for i in range(len(rate_r)):
                rate_r[i].drop(rate_r[i].loc[purify_p(rate_r[i], 10): len(rate_r[i])].index, inplace=True)

            map = 'Paired'
            a = plt.subplot(111)
            for j in range(len(pre_p)):
                a.plot(pre_p[j].strain *100, pre_p[j].stress, color='grey', alpha = 0.4)
                a.plot(pre_r[j].strain *100, pre_r[j].stress, color='grey', linestyle = '--', alpha = 0.4)
            for j in range(len(rate_p)):
                a.plot(rate_p[j].strain *100, rate_p[j].stress, color= color[j])
                a.plot(rate_r[j].strain *100, rate_r[j].stress, color=color[j], linestyle='--')


            a.set(xlim = [-0.4, 12], ylabel = 'Stress(kPa)', xlabel = 'Strain(%)')

            ### setting for gram y axis
            ay = a.twinx()
            ay.set_ylabel('Load (mN)', rotation = 270, va = 'bottom')
            B = -1 * 9.8
            T = 13 * 9.8
            ay.set(ylim=[B, T])
            SB = B * 0.001 / area * 10 ** 9
            ST = T * 0.001 / area * 10 ** 9
            a.set(ylim=[SB, ST])

            grey_line = mlines.Line2D([], [], color = 'grey', label = 'Pre-stretch (10mm/min)' )
            C0_line = mlines.Line2D([], [], color = 'C0', label = '10 mm/min' )
            C1_line = mlines.Line2D([], [], color = 'C1', label = '5 mm/min' )
            C2_line = mlines.Line2D([], [], color = 'C2', label = '1 mm/min' )
            C3_line = mlines.Line2D([], [], color = 'C3', label = '0.5 mm/min' )
            C4_line = mlines.Line2D([], [], color = 'C4', label = '0.1 mm/min' )
            a.set(xlabel = 'Strain (%)', ylabel = 'Stress (kPa)')
            a.legend(handles = [grey_line, C0_line, C1_line, C2_line, C3_line, C4_line], loc = 2, fontsize=10, frameon = False)



            n += 1
            cc += 1


            plt.gcf().subplots_adjust(bottom=0.15, right= 0.815)
            plt.show()





