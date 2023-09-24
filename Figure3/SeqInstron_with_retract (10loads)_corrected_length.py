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
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 2
mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['ytick.major.size'] = 5
mpl.rcParams['ytick.major.width'] = 2
mpl.use('macosx')
fig = mpl.pyplot.gcf()
fig.set_size_inches(8, 6)

output = '/Users/jingyiyu/Documents/Cosgrovelab/manuscript/CW_mechanics_Instron/Data_collection/Python_plot/'


thickness = 7
width = 3
smf = 0.05
intv = 10
per = 10
num = 40

example = '/Users/jingyiyu/Documents/Cosgrovelab/Onion_mechanics/SeqInstronWRetract/10_seq_load'

# target for seuqntial Instron
target = [2,4,6,8,10,10,12,14,16,18,20,20]

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

def tar_ind(pull, load):
    if load == 'N':
        return len(pull)
    else:
        stress = load * 0.0098/ (thickness * width * 0.001)
        for i in range(len(pull) - 5):
            if (pull.stress[i] > stress) & (pull.stress[i + 1] > stress) & (pull.stress[i + 2] > stress) & \
                    (pull.stress[i + 3] > stress) & (pull.stress[i + 4] > stress):
                return i


def purify_p(pull, target):
    for i in range(len(pull)):
        if i < (len(pull)-1) :
            if pull.load[i] > target:
                return i
                break
        else:
            return (len(pull) -1)

def purify_revrs(pull):
    b = np.array(pull.position)
    maxp_index = np.argmax(b)
    threshload = pull.load[maxp_index]
    for i in range(len(pull)):
        if pull.load[i] > threshload:
            rmindex.append(i)


# define a function that calculate the stress and strain of the curve
def norm(pul):
    # strain calculated based on very original length (before plastic deformation)
    # ori_p = ten_p(first_pull)
    # strain calculated based on length at the beginning of each pull
    ori_p = ten_p(pull[ini])

    pul['force_N'] = pul.load * 0.0098

    # engineering stress & strain
    pul['stress'] = pul.force_N / (thickness * width * 0.001)
    pul['strain'] = (pul.position - ori_p) / (5000 + ori_p)


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
ini = 0
for file in sorted(glob.glob(os.path.join(example, '*SHORT.csv'))):
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


        # pull = [first_pull, second_pull, third_pull, fourth_pull, fifth_pull, sixth_pull, seventh_pull, eighth_pull, ninth_pull, tenth_pull, eleventh_pull, twlveth_pull]
        # retract = [first_retract, second_retract, third_retract, fourth_retract, fifth_retract, sixth_retract, seventh_retract, eighth_retract, ninth_retract, tenth_retract, eleventh_retract, twlveth_retract]
        pull = [first_pull, second_pull, third_pull, fourth_pull, fifth_pull, sixth_pull,seventh_pull , eighth_pull, ninth_pull, tenth_pull, eleventh_pull, twlveth_pull]
        retract = [first_retract, second_retract, third_retract, fourth_retract, fifth_retract, sixth_retract,seventh_retract, eighth_retract, ninth_retract, tenth_retract, eleventh_retract, twlveth_retract]
        lw = 2
        ### normalize curve
        for i in range(len(pull)):
            ini = i
            norm(pull[i])
            norm(retract[i])

        pull = [first_pull, second_pull, third_pull, fourth_pull, fifth_pull, seventh_pull , eighth_pull, ninth_pull, tenth_pull, eleventh_pull, twlveth_pull]
        retract = [first_retract, second_retract, third_retract, fourth_retract, fifth_retract,seventh_retract, eighth_retract, ninth_retract, tenth_retract, eleventh_retract, twlveth_retract]

        map = 'Paired'
        a = plt.subplot(211)
        for j in range(len(pull)):
            a.plot(pull[j].strain *100, pull[j].stress, color=color[j], linewidth = lw)
            a.plot(retract[j].strain *100, retract[j].stress, color=color[j], linestyle='--', linewidth = lw)
            # a0.scatter(retract[j].strain *100, retract[j].stress, color=color[j], s=2)
        a.set_xlim(left = -1)
        # a.grid(alpha = 0.4, linestyle = '--')


        # loop for seuqntial Instron
        for i in range(len(pull)):

            retract[i].reset_index(inplace = True)
            rmindex = []
            purify_revrs(retract[i])
            retract[i].drop(retract[i].loc[rmindex].index, inplace = True)

            retract[i].drop(retract[i].loc[0:ten_ind(retract[i])].index, inplace = True)
            pull[i].drop(pull[i].loc[0:ten_ind(pull[i])].index, inplace = True)
            retract[i].reset_index(inplace = True)
            pull[i].reset_index(inplace=True)


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

        sm1retract = pd.DataFrame(data=smooth(first_retract,smf), columns=['strain','stress'])
        sm2retract = pd.DataFrame(data=smooth(second_retract,smf), columns=['strain','stress'])
        sm3retract = pd.DataFrame(data = smooth(third_retract,smf),columns=['strain','stress'])
        sm4retract = pd.DataFrame(data = smooth(fourth_retract,smf),columns=['strain','stress'])
        sm5retract = pd.DataFrame(data = smooth(fifth_retract,smf),columns=['strain','stress'])
        sm6retract = pd.DataFrame(data = smooth(sixth_retract,smf),columns=['strain','stress'])
        sm7retract = pd.DataFrame(data=smooth(seventh_retract,smf), columns=['strain','stress'])
        sm8retract = pd.DataFrame(data=smooth(eighth_retract,smf), columns=['strain','stress'])
        sm9retract = pd.DataFrame(data = smooth(ninth_retract,smf),columns=['strain','stress'])
        sm10retract = pd.DataFrame(data = smooth(tenth_retract,smf),columns=['strain','stress'])
        sm11retract = pd.DataFrame(data = smooth(eleventh_retract,smf),columns=['strain','stress'])
        sm12retract = pd.DataFrame(data = smooth(twlveth_retract,smf),columns=['strain','stress'])

        smpull = [sm1pull, sm2pull, sm3pull, sm4pull, sm5pull, sm7pull, sm8pull, sm9pull, sm10pull, sm11pull, sm12pull]
        smretract = [sm1retract, sm2retract, sm3retract, sm4retract, sm5retract, sm7retract, sm8retract, sm9retract, sm10retract, sm11retract, sm12retract]

        ela_target = ['N', 2, 4, 6, 8, 'N', 10, 12, 14, 16, 18, 'N']
        mo_pull = []
        mo_re = []

        cx = plt.subplot(427)
        c0x = plt.subplot(425)
        for j in range(len(smpull)):
            # cx.plot(smpull[j].strain, np.log(smpull[j].stress), color=color[j])
            c0x.plot(smpull[j].strain *100, derv(smpull[j], intv), color=color[j], linewidth = lw)
            cx.plot(smpull[j].stress , derv(smpull[j], intv), color=color[j], linewidth = lw)


        cx.set(ylabel='Modulus (MPa)', xlabel='Stress (MPa)')
        c0x.set(ylabel='Modulus (MPa)', xlabel='Strain (%)')

        dx = plt.subplot(428)
        d0x = plt.subplot(426)
        for j in range(len(smpull)):
            dx.plot(smretract[j].stress, derv(smretract[j], intv), color=color[j], linewidth = lw)
            d0x.plot(smretract[j].strain *100, derv(smretract[j], intv), color=color[j], linewidth = lw)

        #     dx.plot(smretract[j].strain, np.log(smretract[j].stress), color = color[j], linestyle = '--')
        dx.set(ylabel='Modulus (MPa)', xlabel='Stress (MPa)')
        d0x.set(ylabel='Modulus (MPa)', xlabel='Strain (%)')

        cx.set_xticks([0, 2.5, 5.0, 7.5])
        dx.set_xticks([0, 2.5, 5.0, 7.5])
        a.set(ylabel = 'Stress (MPa)', xlabel = 'Strain (%)')

        # setting for gram y axis
        ay = a.twinx()
        ay.set_ylabel('Loading (mN)', rotation = 270, va = 'bottom')
        B = -1 * 9.8
        T = 22 * 9.8
        ay.set(ylim=[B, T])
        SB = B * 0.001 / (thickness * width * 0.001)
        ST = T * 0.001 / (thickness * width * 0.001)
        a.set(ylim=[SB, ST])

        n += 1

        fig = mpl.pyplot.gcf()
        fig.set_size_inches(8, 8)
        plt.gcf().subplots_adjust(wspace= 0.3, hspace= 0.4)
        # plt.savefig(''f'{output}/Seq_cyclic_corrected.pdf', transparent = True)

        plt.show()





