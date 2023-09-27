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
fig.set_size_inches(12.8, 4.8)

thickness = 7
width = 3
smf = 0.05
intv = 10
per = 10
num = 40

output = '/Users/jingyiyu/Documents/Cosgrovelab/manuscript/CW_mechanics_Instron/Data_collection/Python_plot/'
RT = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/contineous_rate_change/Plot'
cold = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/4C_experiment/Onion_4c/plot'


# define a function that calls the index of certain value
def ind(df, value):
    for i in range(len(df)):
        if df.position[i] == value:
            return i

def ten_p(pull):
    for i in range(len(pull) - 5):
        if (pull.load[i] > 0.1) & (pull.load[i + 1] > 0.1) & (pull.load[i + 2] > 0.1) & \
                (pull.load[i + 3] > 0.1) & (pull.load[i + 4] > 0.1) :
            point = pull.position[i]
            return point

def pos_p(pull, p_pull):
    for i in range(len(pull) - 5):
        n = np.max(p_pull.strain)
        if (pull.strain[i] > n) & (pull.strain[i + 1] > n) & (pull.strain[i + 2] > n) & \
                (pull.strain[i + 3] > n) & (pull.strain[i + 4] > n):
            index = i
            return index
        elif i == (len(pull) - 6):
            return i

def ten_ind(pull):
    for i in range(len(pull) - 5):
        if (pull.load[i] > 0.1) & (pull.load[i + 1] > 0.1) & (pull.load[i + 2] > 0.1) & \
                (pull.load[i + 3] > 0.1) & (pull.load[i + 4] > 0.1):
            return i

# define a function that calculate the stress and strain of the curve
def norm(pull):
    # strain calculated based on very original length (before plastic deformation)
    ori_p = ten_p(first_pull)

    pull['force_N'] = pull.load * 0.0098
    # engineering stress & strain
    pull['stress'] = pull.force_N / (thickness * width * 0.001)
    pull['strain'] = (pull.position - ori_p) / (5000 +ori_p)

# define a function that smooth the curve
def smooth(pull,f):
    lowess = sm.nonparametric.lowess
    z = lowess(pull.stress, pull.strain,frac= f)
    return z

color = ['C0','C1','C2', 'C3', 'C4', 'C5','C6','C7','C8','C9','C10','C11']
n = 0
map = 'Paired'
a = plt.subplot(121)
b = plt.subplot(122)
lw = 2
for file in sorted(glob.glob(os.path.join(RT, '*SHORT.csv'))):
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

        print(ten_p(first_pull), ten_p(second_pull), ten_p(third_pull), ten_p(fourth_pull), ten_p(fifth_pull), ten_p(sixth_pull), ten_p(seventh_pull))
        print(ten_p(eighth_pull), ten_p(ninth_pull), ten_p(tenth_pull), ten_p(eleventh_pull), ten_p(twlveth_pull))

        ### normalize curve
        for i in range(len(pre_p)):
            norm(pre_p[i])
            norm(pre_r[i])
        for i in range(len(rate_p)):
            norm(rate_p[i])
            norm(rate_r[i])

        sm1pull = pd.DataFrame(data = smooth(first_pull,smf),columns=['strain','stress'])
        sm2pull = pd.DataFrame(data = smooth(second_pull,smf),columns=['strain','stress'])
        sm3pull = pd.DataFrame(data = smooth(third_pull,smf),columns=['strain','stress'])
        sm4pull = pd.DataFrame(data = smooth(fourth_pull,smf),columns=['strain','stress'])
        sm5pull = pd.DataFrame(data = smooth(fifth_pull,smf),columns=['strain','stress'])
        sm6pull = pd.DataFrame(data = smooth(sixth_pull,smf),columns=['strain','stress'])
        sm7pull = pd.DataFrame(data = smooth(seventh_pull,smf),columns=['strain','stress'])
        sm8pull = pd.DataFrame(data = smooth(eighth_pull,smf),columns=['strain','stress'])
        sm9pull = pd.DataFrame(data = smooth(ninth_pull,smf),columns=['strain','stress'])
        sm10pull = pd.DataFrame(data = smooth(tenth_pull,smf),columns=['strain','stress'])
        sm11pull = pd.DataFrame(data = smooth(eleventh_pull,smf),columns=['strain','stress'])
        sm12pull = pd.DataFrame(data = smooth(twlveth_pull,smf),columns=['strain','stress'])

        sm1retract = pd.DataFrame(data=smooth(first_retract, smf), columns=['strain', 'stress'])
        sm2retract = pd.DataFrame(data=smooth(second_retract, smf), columns=['strain', 'stress'])
        sm3retract = pd.DataFrame(data=smooth(third_retract, smf), columns=['strain', 'stress'])
        sm4retract = pd.DataFrame(data=smooth(fourth_retract, smf), columns=['strain', 'stress'])
        sm5retract = pd.DataFrame(data=smooth(fifth_retract, smf), columns=['strain', 'stress'])
        sm6retract = pd.DataFrame(data = smooth(sixth_retract,smf),columns=['strain','stress'])
        sm7retract = pd.DataFrame(data=smooth(seventh_retract, smf), columns=['strain', 'stress'])
        sm8retract = pd.DataFrame(data=smooth(eighth_retract, smf), columns=['strain', 'stress'])
        sm9retract = pd.DataFrame(data=smooth(ninth_retract, smf), columns=['strain', 'stress'])
        sm10retract = pd.DataFrame(data=smooth(tenth_retract, smf), columns=['strain', 'stress'])
        sm11retract = pd.DataFrame(data=smooth(eleventh_retract, smf), columns=['strain', 'stress'])
        sm12retract = pd.DataFrame(data=smooth(twlveth_retract, smf), columns=['strain', 'stress'])

        smpre_p = [sm1pull, sm2pull, sm3pull, sm4pull, sm5pull, sm6pull, sm7pull]
        smpre_r = [sm1retract, sm2retract, sm3retract, sm4retract, sm5retract, sm6retract,sm7retract]
        smrate_p = [sm8pull, sm9pull, sm10pull, sm11pull, sm12pull]
        smrate_r = [sm8retract, sm9retract, sm10retract, sm11retract, sm12retract]

        for j in range(len(pre_p)):
            a.plot(pre_p[j].strain *100, pre_p[j].stress, color='grey', alpha = 0.6, linewidth = lw)
            a.plot(pre_r[j].strain *100, pre_r[j].stress, color='grey', linestyle = '--', alpha = 0.6, linewidth = lw)
        for j in range(len(rate_p)):
            a.plot(rate_p[j].strain *100, rate_p[j].stress, color= color[j], linewidth = lw)
            a.plot(rate_r[j].strain *100, rate_r[j].stress, color=color[j], linestyle='--', linewidth = lw)

        a.set(ylim = [-1,10], ylabel = 'Stress(MPa)', xlabel = 'Strain(%)')
        # setting for gram y axis
        ay = a.twinx()
        ay.set_ylabel('Load (mN)', rotation = 270, va = 'bottom')
        B = -1 * 9.8
        T = 12 * 9.8
        ay.set(ylim=[B, T])
        SB = B * 0.001 / (thickness * width * 0.001)
        ST = T * 0.001 / (thickness * width * 0.001)
        a.set(ylim=[SB, ST])
        grey_line = mlines.Line2D([], [], color = 'grey', label = 'Pre-stretch (10mm/min)' )
        C0_line = mlines.Line2D([], [], color = 'C0', label = '10 mm/min' )
        C1_line = mlines.Line2D([], [], color = 'C1', label = '5 mm/min' )
        C2_line = mlines.Line2D([], [], color = 'C2', label = '1 mm/min' )
        C3_line = mlines.Line2D([], [], color = 'C3', label = '0.5 mm/min' )
        C4_line = mlines.Line2D([], [], color = 'C4', label = '0.1 mm/min' )
        a.set(xlabel = 'Strain (%)', ylabel = 'Stress (MPa)')
        a.legend(handles = [grey_line, C0_line, C1_line, C2_line, C3_line, C4_line], loc = 2, fontsize=10, frameon = False)

        n += 1

n = 0

for file in sorted(glob.glob(os.path.join(cold, '*SHORT.csv'))):
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

        print(ten_p(first_pull), ten_p(second_pull), ten_p(third_pull), ten_p(fourth_pull), ten_p(fifth_pull), ten_p(sixth_pull), ten_p(seventh_pull))
        print(ten_p(eighth_pull), ten_p(ninth_pull), ten_p(tenth_pull), ten_p(eleventh_pull), ten_p(twlveth_pull))

        ### normalize curve
        for i in range(len(pre_p)):
            norm(pre_p[i])
            norm(pre_r[i])
        for i in range(len(rate_p)):
            norm(rate_p[i])
            norm(rate_r[i])

        for j in range(len(pre_p)):
            b.plot(pre_p[j].strain *100, pre_p[j].stress, color='grey', alpha = 0.6, linewidth = lw)
            b.plot(pre_r[j].strain *100, pre_r[j].stress, color='grey', linestyle = '--', alpha = 0.6, linewidth = lw)
        for j in range(len(rate_p)):
            b.plot(rate_p[j].strain *100, rate_p[j].stress, color= color[j], linewidth = lw)
            b.plot(rate_r[j].strain *100, rate_r[j].stress, color=color[j], linestyle='--', linewidth = lw)


        b.set(ylim = [-1,10], ylabel = 'Stress(MPa)', xlabel = 'Strain(%)')
        by = b.twinx()
        by.set_ylabel('Load (mN)', rotation = 270, va = 'bottom')
        B = -1 * 9.8
        T = 12 * 9.8
        by.set(ylim=[B, T])
        SB = B * 0.001 / (thickness * width * 0.001)
        ST = T * 0.001 / (thickness * width * 0.001)
        b.set(ylim=[SB, ST])
        grey_line = mlines.Line2D([], [], color = 'grey', label = 'Pre-stretch (10mm/min)' )
        C0_line = mlines.Line2D([], [], color = 'C0', label = '10 mm/min' )
        C1_line = mlines.Line2D([], [], color = 'C1', label = '5 mm/min' )
        C2_line = mlines.Line2D([], [], color = 'C2', label = '1 mm/min' )
        C3_line = mlines.Line2D([], [], color = 'C3', label = '0.5 mm/min' )
        C4_line = mlines.Line2D([], [], color = 'C4', label = '0.1 mm/min' )
        b.set(xlabel = 'Strain (%)', ylabel = 'Stress (MPa)')
        b.legend(handles = [grey_line, C0_line, C1_line, C2_line, C3_line, C4_line], loc = 2, fontsize=10, frameon = False)



        n += 1


plt.gcf().subplots_adjust(bottom=0.2, wspace= 0.4)
plt.show()





